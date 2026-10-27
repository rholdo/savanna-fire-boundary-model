# Spatial_structure.R
#
# Characterizes the spatial structure of the clustering runs beyond average
# nearest-neighbor distance.
#
# Runs the no-boundary clustering model (same as Boundary_NN.R) for:
#   lo  -- Fig. 6a (Pfire = 0.2, z = 20)
#   hi  -- Fig. 6b and Fig. S11 (Pfire = 0.6, z = 10)
#   gap -- Fig. S12 (Pfire = 0.2, z = 5); map correlations only, not in Fig. S9
# saving tree maps every 500 y over the second half of each run, and computes
# the pair correlation function g(r) of large trees (S > 10 cm) on the torus
# (no edge correction needed), with a pointwise 95% CSR envelope from 99
# simulations (Fig. S9, saved as FigS9.png).
#
# Map-correlation check (reported in the text, no figure): gridded basal-area
# maps are correlated between snapshots at increasing time lags. Clusters
# that stay in place would keep maps positively correlated; clusters that
# turn over drive the correlation towards zero (for two unrelated 30 x 30
# maps, about 0 +/- 0.03). The summary table is printed to the console.

library(ggplot2)
library(dplyr)
library(parallel)
library(Rcpp)

sourceCpp("spatial_ops.cpp")
source("fire_functions.R")

Tmax  <- 5000 # Duration (years)
scale <- 3    # Scale of simulation domain (x100 m)
N     <- 400  # Stem density (ha-1)
NREP  <- 5    # Replicate runs per scenario
pcf_years  <- seq(Tmax / 2, Tmax, by = 500)  # snapshots plotted in Fig. S9
lag_base   <- Tmax - 500                     # base year for the short lags
short_lags <- c(10, 25, 50, 100, 250)        # short lags for the map correlation
snap_years <- sort(unique(c(pcf_years, lag_base + short_lags)))
FIG_SCEN   <- c("lo", "hi")  # scenarios shown in Fig. S9
MAP_CELL   <- 10  # grid cell (m) for the basal-area maps (= z in Fig. 6b)
S_LARGE <- 10 # size threshold (cm) for "large" trees, as in Figs. 5 and 6
# Same fire implementation as Boundary_NN.R (open-area seeding)
FIRE <- FIRE_CLUSTER

scenarios <- data.frame(
  Scenario  = c("lo", "hi", "gap"),
  Label     = c("Pfire = 0.2, z = 20 (Fig. 6a)",
                "Pfire = 0.6, z = 10 (Fig. 6b)",
                "Pfire = 0.2, z = 5 (Fig. S12)"),
  Pfire     = c(0.2, 0.6, 0.2),
  z         = c(20, 10, 5),
  sigma.ind = 0.15, sigma.res = 0.0, mu = 0.0075
)

# Model: identical to tree_mod() in Boundary_NN.R except that it returns
# snapshots of tree locations and sizes instead of nearest-neighbor summaries.
tree_mod_snap <- function(Pfire, z, sigma.ind, sigma.res, mu, snap_years){
  
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
  nn_cell_size <- 5                 # NN grid cell (m); speed only, search is exact

  # Initialize model ------------------------------------------
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
  
  # Output: snapshots of all trees at the requested years
  snaps <- list()
  if (0 %in% snap_years)                        # initial state
    snaps[[1]] <- data.frame(Year = 0, x = x, y = y, S = S)

  for (t in 1:Tmax){
    
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
    deltas <- ifelse(res,
                     lambdar * Smax / (shs + Smax) + Sdev + rnorm(K, 0, sigma.res),
                     lambda + Sdev + rnorm(K, 0, sigma.res))
    deltas <- ifelse(deltas < 0, 0, deltas)
    S <- S + deltas * f * exp(-rho * B)

    # Resprout status
    res  <- ifelse(res, ifelse(S < 2 & S < Smax, TRUE, FALSE), FALSE)
    Smax <- ifelse(S > Smax, S, Smax)

    # Fire
    fire    <- runif(1) < Pfire
    Pignite <- 1 / (1 + i0 * exp(i1 * B))
    ignite  <- runif(K) < Pignite
    ignite  <- fire_ignition(x, y, ignite, Pignite, domain, FIRE)
    RH <- RHref + r * B
    W  <- Wref * exp(-w * B)
    Q  <- m0 + m1 * G - m2 * sqrt(Mref) - m3 * RH ^ 2 - m4 / W
    Q  <- ifelse(Q < 0, 0, Q)

    if (fire){
      logH     <- a1 * log(S) + a0
      eta      <- t0 - t1 * logH + t2 * sqrt(Q)
      Ptopkill <- ignite * exp(eta) / (1 + exp(eta))
      topkill  <- runif(K) < Ptopkill
      S[topkill]   <- 0.5
      res[topkill] <- TRUE
    }

    # Background mortality
    mort <- runif(K) < mu
    S[mort]    <- 0.5
    x[mort]    <- runif(sum(mort), min = 0, max = domain)
    y[mort]    <- runif(sum(mort), min = 0, max = domain)
    Sdev[mort] <- rnorm(sum(mort), mean = 0, sd = sigma.ind)
    Sdev[mort] <- ifelse(lambda + Sdev[mort] < 0, -lambda, Sdev[mort])
    res[mort]  <- FALSE
    
    if (t %in% snap_years) {
      snaps[[length(snaps) + 1]] <- data.frame(Year = t, x = x, y = y, S = S)
    }
    
  } # end time loop
  
  return(do.call(rbind, snaps))
}

# Run simulations ------------------------------------------------------------
df <- expand.grid(Rep = 1:NREP, Scen = seq_len(nrow(scenarios)))
df$Iter <- seq_len(nrow(df))

cl <- makeCluster(min(detectCores(), nrow(df)))
clusterExport(cl, varlist = c(fire_exports, "tree_mod_snap", "scenarios",
                              "snap_years", "scale", "Tmax", "N", "df"))
cpp_path <- file.path(getwd(), "spatial_ops.cpp")
clusterExport(cl, varlist = "cpp_path")
clusterEvalQ(cl, { library(Rcpp); sourceCpp(cpp_path) })
snaps <- parLapply(cl, df$Iter, function(i) {
  sc  <- scenarios[df$Scen[i], ]
  out <- tree_mod_snap(sc$Pfire, sc$z, sc$sigma.ind, sc$sigma.res, sc$mu,
                       snap_years)
  out$Scenario <- sc$Scenario
  out$Rep <- df$Rep[i]
  out
})
stopCluster(cl)
snaps  <- bind_rows(snaps)
domain <- 100 * scale
lab    <- setNames(scenarios$Label, scenarios$Scenario)
snaps$Scenario <- factor(snaps$Scenario, levels = scenarios$Scenario)  # panel order: 6a, 6b

# 1. Pair correlation function on the torus ----------------------------------
pcf_torus <- function(x, y, domain, rmax = 50, dr = 1) {
  n   <- length(x)
  cnt <- toroidal_pcf_counts(x, y, domain, rmax, dr)
  r0  <- seq(0, by = dr, length.out = length(cnt))
  data.frame(r = r0 + dr / 2,
             g = 2 * cnt * domain ^ 2 /
               (n * (n - 1) * pi * ((r0 + dr) ^ 2 - r0 ^ 2)))
}
csr_envelope <- function(n, domain, nsim = 99, rmax = 50, dr = 1) {
  gs <- sapply(1:nsim, function(k)
    pcf_torus(runif(n, 0, domain), runif(n, 0, domain), domain, rmax, dr)$g)
  data.frame(r  = seq(dr / 2, by = dr, length.out = nrow(gs)),
             lo = apply(gs, 1, quantile, 0.025),
             hi = apply(gs, 1, quantile, 0.975))
}

pcf <- snaps %>% filter(Scenario %in% FIG_SCEN, Year %in% pcf_years, S > S_LARGE) %>%
  group_by(Scenario, Rep, Year) %>% filter(n() >= 20) %>%
  group_modify(~ cbind(pcf_torus(.x$x, .x$y, domain), n = nrow(.x))) %>%
  ungroup()
# CSR envelope for the median number of large trees in each scenario
env <- pcf %>% distinct(Scenario, Rep, Year, n) %>%
  group_by(Scenario) %>% summarize(n = round(median(n)), .groups = "drop") %>%
  group_by(Scenario) %>% group_modify(~ csr_envelope(.x$n, domain)) %>%
  ungroup()

pcf_plot <- pcf %>% filter(Year >= Tmax / 2) %>%
  group_by(Scenario, Year, r) %>% summarize(g = mean(g), .groups = "drop")
figS9 <- ggplot() +
  geom_ribbon(data = env, aes(x = r, ymin = lo, ymax = hi), fill = 'grey85') +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(data = pcf_plot, aes(x = r, y = g, colour = factor(Year))) +
  facet_wrap(~ Scenario, ncol = 1, labeller = as_labeller(lab), drop = TRUE) +
  labs(x = 'Distance r (m)', y = 'Pair correlation g(r), trees > 10 cm',
       colour = 'Year') +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = rel(1)))

# 2. Map correlation through time (cluster turnover) -------------------------
ba_grid <- function(x, y, S, domain, cell = MAP_CELL) {
  n  <- as.integer(round(domain / cell))
  ix <- pmin(floor(x / cell) + 1, n)
  iy <- pmin(floor(y / cell) + 1, n)
  ba <- pi * (S / 200) ^ 2                       # m2 per tree
  m  <- tapply(ba, list(factor(iy, levels = 1:n), factor(ix, levels = 1:n)),
               sum)
  m[is.na(m)] <- 0
  as.vector(m) / (cell ^ 2 / 1e4)                # m2 ha-1
}
maps <- snaps %>% group_by(Scenario, Rep, Year) %>%
  summarize(ba = list(ba_grid(x, y, S, domain)), .groups = "drop")

# Pairs of snapshot years: short lags from lag_base, plus all pairs of the
# 500-y snapshots (lags 500 to 2500 y)
pairs <- rbind(data.frame(t1 = lag_base, t2 = lag_base + c(short_lags, 500)),
               subset(expand.grid(t1 = pcf_years, t2 = pcf_years), t2 > t1))
pairs <- unique(pairs)
pairs$lag <- pairs$t2 - pairs$t1

mapcor <- maps %>% group_by(Scenario, Rep) %>%
  group_modify(function(d, key) {
    get <- function(t) d$ba[[match(t, d$Year)]]
    data.frame(pairs,
               r = mapply(function(a, b) cor(get(a), get(b)),
                          pairs$t1, pairs$t2))
  }) %>% ungroup()
mapcor_summary <- mapcor %>% group_by(Scenario, lag) %>%
  summarize(r_mean = mean(r), r_sd = sd(r), n = n(), .groups = "drop")
print(as.data.frame(mapcor_summary), digits = 2)

# Save figure ---------------------------------------------------------------
# Map correlations are printed above (mapcor_summary); no files are written.
ggsave("FigS9.png", figS9, width = 7, height = 9, dpi = 300)