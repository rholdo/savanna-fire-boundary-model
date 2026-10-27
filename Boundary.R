# Tree-grass-fire spatial simulation
# Run on a 300 x 300 torus with N trees/ha
# Examination of woody-grassy boundary dynamics

# This script generates figs. 2 & 3
# Set N = 1000 on line 23 to generate figs. S3 & S7

library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)

# Compile the C++ spatial operations functions
sourceCpp("spatial_ops.cpp")

Tmax  <- 5000 # Duration (years)
scale <- 3    # Scale of simulation domain (x100 m)
N <- 400      # Stem density (ha-1); 400 is default, 1000 also tested

# Shared fire-ignition helpers (see fire_functions.R)
source("fire_functions.R")
# Fire implementation: in each fire year one ignition point
# is drawn at random among trees in the LEFT HALF of the domain with a
# potential ignition and Pignite >= 0.5 (i.e. in open, grassy areas); fire
# then spreads through the full-domain 5-m lattice (see fire_functions.R).
# Edit pign_min / n_seeds / cell / connectivity here to change it.
FIRE <- FIRE_BOUNDARY

# tree_mod(): simulates the persistence or erosion of an imposed
# woody/grassy boundary. Trees are initialized as a mature (S = 20 cm)
# population and left to reach a fire-free steady state for 500 y; at
# t = 500 all trees on the left half of the arena (x < domain/2) are reset
# to saplings (S = 0.5 cm), creating a sharp boundary at x = domain/2. Fire
# is then switched on and the simulation runs a further Tmax - 500 years.
#   Pfire     -- exogenous fire probability per year (P_fire)
#   z         -- neighborhood radius (m) for competition and fire effects
#   sigma.ind -- SD of growth rate across individuals (cm y-1)
#   sigma.res -- residual (within-individual) growth SD (cm y-1)
#   mu        -- background (adult) mortality rate (y-1)
# Returns a data frame (years 501:Tmax only) with:
#   B.op       -- basal area (m2 ha-1) on the initially "open"/grassy
#                 (left, x < domain/2) side of the arena
#   B.cl       -- basal area on the initially "closed"/woody (right) side
#   D.break    -- x location (5-m increments) that maximizes the basal-area
#                 difference between left and right segments in a given year
#   B.diff.dyn -- the basal-area difference at that maximizing location
# B.op and B.cl together give the boundary contrast at its initial,
# fixed location (x = domain/2); D.break/B.diff.dyn instead track a moving
# breakpoint, letting a persistent-but-shifted boundary still register as a
# large difference (see manuscript Methods 2.2 and Fig. 3).
tree_mod <- function(Pfire, z, sigma.ind, sigma.res, mu){
  
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
                    D.break = numeric(Tmax), B.diff.dyn = numeric(Tmax))
  
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
    # Update ignite with connectivity map
    ignite <- fire_ignition(x, y, ignite, Pignite, domain, FIRE)
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

# Run model through scenarios
Pfire <- round(seq(from = 0.2, to = 0.8, by = 0.1), 1) # round() so Pfire == 0.6 etc. match exactly
z <- c(5, 10, 20)
sigma.ind <- seq(from = 0.0, to = 0.3, by = 0.15)
sigma.res <- c(0, 0.15)
mu <- seq(from = 0.005, to = 0.01, by = 0.0025)

# Run model through driving variable scenarios
ITER <- length(Pfire) * length(z) * length(sigma.ind) * length(sigma.res) * 
  length(mu) * 10

df <- data.frame(Iter = 1:ITER)
df$Pfire <- rep(Pfire, each = length(z) * length(sigma.ind) *
                  length(sigma.res) * length(mu) * 10)
df$z <- rep(z, each = length(sigma.ind) * length(sigma.res) * length(mu) * 10,
           times = length(Pfire))
df$sigma.ind <- rep(sigma.ind, each = length(sigma.res) * length(mu) * 10,
                    times = length(Pfire) * length(z))
df$sigma.res <- rep(sigma.res, each = length(mu) * 10,
                    times = length(Pfire) * length(z) * length(sigma.ind))
df$mu <- rep(mu, each = 10,
            times = length(Pfire) * length(z) * length(sigma.ind) *
              length(sigma.res))
df$Rep <- rep(1:10,
             times = length(Pfire) * length(z) * length(sigma.ind) *
               length(sigma.res) * length(mu))

# Set up cluster
numCores <- detectCores()
cl <- makeCluster(numCores)

# Export everything each worker needs
clusterExport(cl, varlist = c(fire_exports, "tree_mod", "scale", "Tmax", "N", "df"))

# Compile C++ on every worker (path must be absolute)
cpp_path <- file.path(getwd(), "spatial_ops.cpp")
clusterExport(cl, varlist = "cpp_path")
clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp(cpp_path)   # each worker compiles once and caches
})

df.sim <- parLapply(cl, 1:ITER, function(i) {
  tree_mod(df$Pfire[i], df$z[i], df$sigma.ind[i], df$sigma.res[i], df$mu[i])
})

stopCluster(cl)

df.sim <- do.call(rbind, df.sim)
Trun <- Tmax - 500
df.sim$Iter <- rep(1:ITER, each = Trun)
df.sim$Year <- rep(501:Tmax, times = ITER)
boun <- inner_join(df, df.sim, by = "Iter")

# Plot results
# Reference line: the initial (year 501, i.e. right after the boundary
# reset + burn-in) basal-area contrast between the two sides, for comparison
# against the persisting/eroded contrast at the end of each run.
boun.init <- dplyr::filter(boun, Year == 501)
boun.init<- boun.init %>% mutate(B.diff = B.cl - B.op)
Bdiff.init <- mean(boun.init$B.diff)

# Fig. 2: mean +/- SD basal-area difference at the fixed boundary location
# (x = domain/2 = 150 m) over the last 100 y of each run, across parameter
# combinations. Values well above the dotted initial-contrast line indicate
# a persistent boundary; values near/below zero indicate erosion.
boun.fin <- dplyr::filter(boun, Year > 4900)
# Calculate BA difference at fixed location (150 m)
boun.fin <- boun.fin %>% mutate(B.diff = B.cl - B.op)
boun.fin <- boun.fin %>% group_by(Pfire, z, sigma.ind, sigma.res, mu, Rep) %>% 
  summarize(B.diff = mean(B.diff), B.diff.dyn = mean(B.diff.dyn))
boun.fin <- boun.fin %>% group_by(Pfire, z, sigma.ind, sigma.res, mu) %>% 
  summarize(B.diff.mn = mean(B.diff),
            B.diff.sd = sd(B.diff),
            B.diff.dyn.mn = mean(B.diff.dyn),
            B.diff.dyn.sd = sd(B.diff.dyn))
boun.fin <- boun.fin %>% dplyr::filter(Pfire != 0.3 &
                                         Pfire != 0.5 & Pfire != 0.7)
boun.fin <- boun.fin %>% mutate(mu = as.factor(mu),
                                Pfire = as.factor(Pfire),
                                z = as.factor(z))
boun.fin <- boun.fin %>% dplyr::filter(sigma.res == 0.15)
boun.fin <- boun.fin %>% mutate(
  z = case_match(
    z,
    '5' ~ 'z = 5',
    '10' ~ 'z = 10',
    '20' ~ 'z = 20'
  )
)
boun.fin <- boun.fin %>% mutate(
  mu = case_match(
    mu,
    '0.005' ~ '0.0050',
    '0.0075' ~ '0.0075',
    '0.01' ~ '0.0100'
  )
)
boun.fin <- boun.fin %>%
  mutate(z = factor(z, levels = c('z = 5', 'z = 10', 'z = 20')))
boun.fin <- boun.fin %>% mutate(
  Pfire = case_match(
    Pfire,
    '0.2' ~ 'Pfire = 0.2',
    '0.4' ~ 'Pfire = 0.4',
    '0.6' ~ 'Pfire = 0.6',
    '0.8' ~ 'Pfire = 0.8'
  )
)

fig2 <- ggplot(boun.fin, aes(x = sigma.ind, y = B.diff.mn, 
                                 group = mu, col = mu)) + 
  geom_point(size = 3, position = position_dodge(width = 0.05)) +
  geom_errorbar(aes(x = sigma.ind, ymin = B.diff.mn - B.diff.sd,
                    ymax = B.diff.mn + B.diff.sd),
                position = position_dodge(width = 0.05), width = 0.05, size = 1) +
  geom_hline(yintercept = Bdiff.init, linetype = 3) +
  geom_hline(yintercept = 0) +
  facet_grid(rows = vars(Pfire), cols = vars(z)) +
  labs(x = expression(paste(sigma[ind],' (cm y'^-1,')')),
       y = expression(paste(italic('B'),' difference (m'^2,' ha'^-1,')'))) +
  scale_colour_manual(name = expression(mu), values = c('red','blue','green')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.5)),
        axis.title.x = element_text(vjust = 0, size = rel(1.5)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2), hjust = 0.5),
        legend.text = element_text(size = rel(1.2)))

# Fig. 3: boundary trajectory figure. Tracks the *moving* breakpoint
# (D.break, B.diff.dyn) over time for a single parameter combination
# (sigma.res = 0, sigma.ind = 0, mu = 0.005), plotting basal-area contrast
# against its x location every 20 y. A path that stays near x = 150 m
# indicates a stationary boundary; drift away from 150 m indicates a
# traveling front. Start (green) and end (red) points mark the trajectory
# endpoints.
boun.traj <- dplyr::filter(boun, sigma.res == 0.0)
boun.traj <- dplyr::filter(boun.traj, Pfire != 0.3 &
                             Pfire != 0.5 & Pfire != 0.7)
boun.traj <- dplyr::filter(boun.traj, mu == 0.005)
boun.traj <- dplyr::filter(boun.traj, sigma.ind == 0.0)

# Average across Rep for each combination of grouping variables and Year
boun.traj <- boun.traj %>%
  group_by(Pfire, z, sigma.ind, sigma.res, mu, Year) %>%
  summarize(D.break = mean(D.break), B.diff.dyn = mean(B.diff.dyn),
            .groups = "drop")

# Retain only every 10th year
boun.traj <- dplyr::filter(boun.traj, Year %% 20 == 0)

boun.traj <- boun.traj %>%
  mutate(sigma.ind = as.factor(sigma.ind),
         Pfire = as.factor(Pfire),
         z = as.factor(z))

boun.traj <- boun.traj %>% mutate(
  z = case_match(z,
                 '5'  ~ 'z = 5',
                 '10' ~ 'z = 10',
                 '20' ~ 'z = 20'),
  z = factor(z, levels = c('z = 5', 'z = 10', 'z = 20')))

boun.traj <- boun.traj %>% mutate(
  Pfire = case_match(Pfire,
                     '0.2' ~ 'Pfire = 0.2',
                     '0.4' ~ 'Pfire = 0.4',
                     '0.6' ~ 'Pfire = 0.6',
                     '0.8' ~ 'Pfire = 0.8'))

traj.start <- boun.traj %>% group_by(Pfire, z) %>% slice_min(Year, n = 1)
traj.end   <- boun.traj %>% group_by(Pfire, z) %>% slice_max(Year, n = 1)

fig3 <- ggplot(boun.traj, aes(x = D.break, y = B.diff.dyn)) +
  geom_vline(xintercept = 150, linetype = 'dashed', colour = 'grey50') +
  geom_path(linewidth = 0.7, col = 'blue') +
  geom_point(data = traj.start, aes(x = D.break, y = B.diff.dyn),
             col = 'green', size = 3) +
  geom_point(data = traj.end, aes(x = D.break, y = B.diff.dyn),
             col = 'red', size = 3) +
  facet_grid(rows = vars(Pfire), cols = vars(z)) +
  labs(x = 'Breakpoint location (m)',
       y = expression(paste(italic('B'),' difference (m'^2,' ha'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.5)),
        axis.title.x = element_text(vjust = 0, size = rel(1.5)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"))

# Save figures: Figs. 2 and 3 (N = 400) or Figs. S3 and S7 (N = 1000)
fig_names <- if (N == 400) c("Fig2.png", "Fig3.png") else c("FigS3.png", "FigS7.png")
ggsave(fig_names[1], fig2, width = 10, height = 10, dpi = 300)
ggsave(fig_names[2], fig3, width = 10, height = 10, dpi = 300)