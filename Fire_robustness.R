# Fire_robustness.R
#
# Tests whether the boundary results (Fig. 4) are robust to how fire is
# implemented, including the scale choices of the percolation-like spread
# algorithm.
#
# Runs the boundary simulation (same model as Boundary.R) for the three
# Fig. 4 parameter sets under:
#   - the baseline seeded fire (one ignition point in the left half with
#     Pignite >= 0.5, 5-m lattice, 8-cell connectivity)
#   i)   no spread: every potential ignition burns (independent ignition)
#   ii)  2.5 x 2.5 and 10 x 10 m fire lattices
#   iii) 4-cell (von Neumann) connectivity
#   iv)  a lower ignition-point threshold (Pignite >= 0.3)
#   v)   fire entering along the whole left edge of the domain
# It also records, for every fire year, the fraction of trees with a
# potential ignition and the fraction that actually burned, so the strength
# of the spread constraint can be reported (table `burn_tab`).
# It reports the basal-area contrast between the initially woody and grassy
# halves over time (Fig. S8; one panel per scenario, saved as FigS8.png) and
# at the end of each run (table `fin`, printed to the console). Raw results
# are saved to Fire_robustness_output/res.rds.

library(ggplot2)
library(cowplot)
library(dplyr)
library(parallel)
library(Rcpp)

sourceCpp("spatial_ops.cpp")
source("fire_functions.R")

Tmax  <- 5000 # Duration (years), as in Boundary.R / Simulation_examples.R
scale <- 3    # Scale of simulation domain (x100 m)
N     <- 400  # Stem density (ha-1)
NREP  <- 10   # Replicate runs per scenario x fire variant
FIRE  <- FIRE_BOUNDARY  # not used directly; kept so fire_exports resolves

# Fire variants -------------------------------------------------------------
# Baseline = FIRE_BOUNDARY (one seed in the left half with Pignite >= 0.5,
# 5-m lattice, 8-cell connectivity, wrapping torus).
fire_variants <- list(
  "Seeded (baseline)"      = FIRE_BOUNDARY,
  "No spread"              = default_fire(mode = "independent"),
  "Seeded, 2.5 m lattice"  = default_fire("seeded", cell = 2.5, pign_min = 0.5,
                                          seed_xfrac = c(0, 0.5)),
  "Seeded, 10 m lattice"   = default_fire("seeded", cell = 10, pign_min = 0.5,
                                          seed_xfrac = c(0, 0.5)),
  "Seeded, 4-cell"         = default_fire("seeded", connectivity = 4,
                                          pign_min = 0.5, seed_xfrac = c(0, 0.5)),
  "Seeded, Pignite >= 0.3" = default_fire("seeded", pign_min = 0.3,
                                          seed_xfrac = c(0, 0.5)),
  "Left-edge entry"        = default_fire(mode = "edge")
)

# Fig. 4 parameter sets -------------------------------------------------------
scenarios <- data.frame(
  Scenario  = c("a) Pfire = 0.2, z = 10", "b) Pfire = 0.8, z = 10",
                "c) Pfire = 0.8, z = 20, no growth variation"),
  Pfire     = c(0.2, 0.8, 0.8),
  z         = c(10, 10, 20),
  sigma.ind = c(0.15, 0.15, 0.0),
  sigma.res = c(0.15, 0.15, 0.0),
  mu        = c(0.01, 0.01, 0.005)
)

# Model: identical to tree_mod() in Boundary.R except that the fire
# implementation is passed as an argument (fire_cfg).
tree_mod_fire <- function(Pfire, z, sigma.ind, sigma.res, mu, fire_cfg){
  
  # Parameters -----------------------------------------
  K     <- N * scale ^ 2     # Number of trees
  gmax  <- 800               # Max grass biomass in g m-2
  alpha <- 0.18              # Grass biomass reduction due to tree shading
  lambda    <- 0.23          # Mean tree annual growth rate (cm)
  lambdar <- 0.61            # Resprout growth rate
  shs       <- 0.8           # Resprout growth half-saturation constant
  rho       <- 0.05          # Tree crowding parameter
  gamma     <- 0.69          # Tree growth reduction due to grass
  a0 <- -0.1                 # Height allometry intercept
  a1 <-  0.77                # Height allometry slope
  RHref <- 20                # Reference RH
  r     <-  0.41             # RH increase per unit basal area
  Wref  <-  5                # Reference wind speed
  w     <-  0.047            # Wind speed decay with basal area
  Mref  <- 40                # Fuel moisture
  i0 <- 0.15                 # Ignition parameters
  i1 <- 0.23
  t0 <- 1.48                 # Topkill parameters (Higgins et al. 2000)
  t1 <- 3.31
  t2 <- 0.026
  m0 <- 2729                 # Fire intensity parameters (Higgins et al. 2000)
  m1 <- 8.684
  m2 <- 530
  m3 <- 0.907
  m4 <- 596

  # Spatial indexing setup -----------------------------------------
  domain    <- 100 * scale
  cell_size <- min(z, domain / 2)  # ba grid: cell = z so 3x3 block covers z
  nn_cell_size <- 5                 # NN grid: 5 m cells, 3x3 covers 15 m
  
  # Initialize model ------------------------------------------
  # Set baseline growth rate for each tree (deltas.nres, i.e. not a resprout)
  Sdev <- rnorm(K, mean = 0, sd = sigma.ind)
  Sdev <- ifelse(lambda + Sdev < 0, -lambda, Sdev)
  # Preserve Sdev to apply to future topkilled trees
  # Flag to indicate if a tree is resprouting following topkill or not
  res  <- rep(FALSE, K)
  # Initialize with mature trees and allow to reach steady state
  S    <- rep(20, K)
  Smax <- S
  x    <- runif(K, min = 0, max = domain)
  y    <- runif(K, min = 0, max = domain)
  
  # Output
  out <- data.frame(B.op  = numeric(Tmax), B.cl  = numeric(Tmax),
                    D.break = numeric(Tmax), B.diff.dyn = numeric(Tmax),
                    fire = FALSE, pot = NA_real_, burn = NA_real_,
                    pot.cl = NA_real_, burn.cl = NA_real_)
  
  for(t in 1:Tmax){
    
    if (t == 500){
      S[x < domain / 2] <- 0.5 # Reset size of trees on lhs (open)
      Smax[x < domain / 2] <- 0.5
    }
    
    # ------------------------------------------------------------------
    # Neighborhood basal area (C++)
    # compute_ba() uses grid-based spatial indexing with toroidal distances.
    # Each tree's own basal area is included in its neighborhood total.
    # ------------------------------------------------------------------
    s2_scaled <- S ^ 2 / (4 * z ^ 2)
    B <- compute_ba(x, y, s2_scaled, z, domain, cell_size)

    # Grass biomass and competitive effect
    G <- gmax * exp(-alpha * B)
    f <- (1 - G / gmax * (1 - gamma))
    
    # Growth
    # If a tree is resprouting and less than 5 cm, use faster topkill growth
    # If not, keep current growth rate
    deltas <- ifelse(res,
      lambdar * Smax / (shs + Smax) + Sdev + rnorm(K, 0, sigma.res),
      lambda + Sdev + rnorm(K, 0, sigma.res))
    deltas <- ifelse(deltas < 0, 0, deltas)
    S <- S + deltas * f * exp(- rho * B)
    
    # Resprout status
    res <- ifelse(res, ifelse(S < 2 & S < Smax, TRUE, FALSE), FALSE)
    Smax <- ifelse(S > Smax, S, Smax)
    
    # Fire
    fire <- ifelse(t > 500, runif(1) < Pfire, FALSE)
    # Local ignition probability if fire has spread
    Pignite <- 1 / (1 + i0 * exp(i1 * B))
    ignite <- runif(K) < Pignite
    pot    <- mean(ignite)              # potential ignitions
    pot.cl <- mean(ignite[x > domain / 2])   # ... in the initially woody half
    # Update ignite with connectivity map
    ignite <- fire_ignition(x, y, ignite, Pignite, domain, fire_cfg)
    if (fire) {                         # burn diagnostics
      out$fire[t] <- TRUE; out$pot[t] <- pot; out$burn[t] <- mean(ignite)
      out$pot.cl[t] <- pot.cl; out$burn.cl[t] <- mean(ignite[x > domain / 2])
    }
    RH <- RHref + r * B
    W <- Wref * exp(-w * B)
    Q <- m0 + m1 * G - m2 * sqrt(Mref) -
      m3 * RH ^ 2 - m4 / W
    Q <- ifelse(Q < 0, 0, Q)
    
    if (fire){
      # Convert diameter to height
      logH <- a1 * log(S) + a0
      eta <- t0 - t1 * logH + t2 * sqrt(Q)
      Ptopkill <- ignite * exp(eta) / (1 + exp(eta))
      topkill <- runif(K) < Ptopkill
      # Reset size of topkilled trees
      S[topkill] <- 0.5
      res[topkill] <- TRUE
    }
    
    # Background mortality
    mort <- runif(K) < mu
    S[mort] <- 0.5
    x[mort] <- runif(sum(mort), min = 0, max = domain)
    y[mort] <- runif(sum(mort), min = 0, max = domain)
    Sdev[mort] <- rnorm(sum(mort), mean = 0, sd = sigma.ind)
    Sdev[mort] <- ifelse(lambda + Sdev[mort] < 0, -lambda, Sdev[mort])
    res[mort] <- FALSE

    # Express basal area in m2 / ha
    out$B.op[t] <- sum((S[x < domain / 2] / 2) ^ 2 * pi) / (1e4 * 0.5 * scale ^ 2)
    out$B.cl[t] <- sum((S[x > domain / 2] / 2) ^ 2 * pi) / (1e4 * 0.5 * scale ^ 2)
    
    # Find the x location (in 5-m increments) that maximises the ba
    # difference between left and right segments. Candidate breakpoints are
    # the unique x coordinates of all trees, so the search is over at most
    # K values.
    
    x_candidates <- seq(20, domain - 20, by = 5)
    b_diff <- sapply(x_candidates, function(xb) {
      b.left  <- sum((S[x <= xb] / 2) ^ 2 * pi) / (1e4 * scale ^ 2 * xb / domain)
      b.right <- sum((S[x >  xb] / 2) ^ 2 * pi) /
        (1e4 * scale ^ 2 * (domain - xb) / domain)
      abs(b.left - b.right)
    })
    best <- which.max(b_diff)
    out$D.break[t] <- x_candidates[best]
    out$B.diff.dyn[t] <- b_diff[best]
  }
  # Output all but first 500 years
  return(out[501:Tmax,])
}

# Run all scenario x variant x replicate combinations ------------------------
df <- expand.grid(Rep = 1:NREP, Variant = names(fire_variants),
                  Scen = seq_len(nrow(scenarios)), stringsAsFactors = FALSE)
df$Iter <- seq_len(nrow(df))

cl <- makeCluster(detectCores())
clusterExport(cl, varlist = c(fire_exports, "tree_mod_fire", "fire_variants",
                              "scenarios", "scale", "Tmax", "N", "df"))
cpp_path <- file.path(getwd(), "spatial_ops.cpp")
clusterExport(cl, varlist = "cpp_path")
clusterEvalQ(cl, { library(Rcpp); sourceCpp(cpp_path) })

df.sim <- parLapply(cl, df$Iter, function(i) {
  sc  <- scenarios[df$Scen[i], ]
  out <- tree_mod_fire(sc$Pfire, sc$z, sc$sigma.ind, sc$sigma.res, sc$mu,
                       fire_variants[[df$Variant[i]]])
  out$Iter <- i
  out$Year <- 501:Tmax
  out
})
stopCluster(cl)

res <- bind_rows(df.sim) %>%
  inner_join(df, by = "Iter") %>%
  mutate(Scenario = scenarios$Scenario[Scen],
         Variant  = factor(Variant, levels = names(fire_variants)),
         B.diff   = B.cl - B.op)

# Fig. S8: basal-area contrast at the initial boundary location over time
# (one panel per scenario; scenario details go in the figure caption)
ts <- res %>% filter(Year %% 50 == 0) %>%
  group_by(Scenario, Variant, Year) %>%
  summarize(B.diff = mean(B.diff), .groups = "drop")
figS8 <- ggplot(ts, aes(x = (Year - 500) / 1e3, y = B.diff,
                            colour = Variant)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Scenario, ncol = 1,
             labeller = as_labeller(function(x) sub("\\).*$", ")", x))) +
  labs(x = 'Time since fire switched on (1000 y)',
       y = expression(paste(italic('B'),' difference (m'^2,' ha'^-1,')')),
       colour = 'Fire implementation') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off",
        strip.text = element_text(face = "bold", size = rel(1.5), hjust = 0,
                                  margin = margin(t = 4, b = 4, l = -25)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.5)),
        axis.title.x = element_text(vjust = 0, size = rel(1.5)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = "right",
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

# Final contrast (last 100 y), fixed and moving breakpoint (table `fin`,
# printed below)
fin <- res %>% filter(Year > Tmax - 100) %>%
  group_by(Scenario, Variant, Rep) %>%
  summarize(B.diff = mean(B.diff), B.diff.dyn = mean(B.diff.dyn),
            .groups = "drop") %>%
  group_by(Scenario, Variant) %>%
  summarize(across(c(B.diff, B.diff.dyn), list(mn = mean, sd = sd)),
            .groups = "drop")

# Strength of the spread constraint -------------------------------------------
# For fire years only: fraction of trees with a potential ignition (pot),
# fraction that burned (burn), burned / potential, the same for the initially
# woody half (x > 150 m; .cl), where the spread rule should matter most, and
# the share of fire years
# in which nothing burned (e.g. no qualifying ignition point in the left half).
burn_tab <- res %>% filter(fire) %>%
  group_by(Scenario, Variant) %>%
  summarize(pot = mean(pot), burn = mean(burn),
            burn_over_pot = mean(burn) / mean(pot),
            pot.cl = mean(pot.cl), burn.cl = mean(burn.cl),
            burn_over_pot.cl = mean(burn.cl) / mean(pot.cl),
            no_fire_years = mean(burn == 0), .groups = "drop")

# Results tables (printed to the console)
print(as.data.frame(fin), digits = 3)
print(as.data.frame(burn_tab), digits = 3)

# Save raw results so they can be inspected without re-running
dir.create("Fire_robustness_output", showWarnings = FALSE)
saveRDS(res, "Fire_robustness_output/res.rds")

# Save figure
ggsave("FigS8.png", figS8, width = 7, height = 9, dpi = 300)