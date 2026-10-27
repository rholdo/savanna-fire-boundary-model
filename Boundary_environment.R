# Tree-grass-fire spatial simulation
# Run on a 300 x 300 torus with N trees/ha
# Examination of woody-grassy boundary dynamics with
# an environmental transition affecting tree growth rate

# This script generates fig. S13

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

# tree_mod_env(): boundary-persistence simulation with an underlying
# environmental gradient in growth rate, instead of the homogeneous
# template used in Boundary.R. Trees on the right-hand ("woody") side of
# the domain (x >= domain / 2) grow at `env` times the default rate for both
# intact stems (lambda) and resprouts (lambdar); trees on the left keep
# the default rate. As in Boundary.R, all trees on the left half are reset
# to saplings at t = 500 to create the initial grassy/woody boundary, and
# fire is on for the remaining Tmax - 500 y. z, sigma.ind, sigma.res and mu
# are fixed core parameters here rather than function arguments (see
# manuscript Methods 2.4/2.5).
#   Pfire -- exogenous fire probability per year (P_fire)
#   env   -- right-to-left growth-rate ratio (1 = homogeneous, as in
#            Boundary.R; up to 2 = right side grows twice as fast)
# Returns a data frame (years 501:Tmax) with the same B.op/B.cl/D.break/
# B.diff.dyn columns as Boundary.R's tree_mod().
tree_mod_env <- function(Pfire, env){

  # Core parameters (fixed here)
  z <- 10 # Neighborhood size in m
  sigma.ind <- 0.2 # SD in growth rate attributed to individual
  sigma.res <- 0.2 # Residual SD in growth rate
  mu <- 0.01 # Background mortality rate
  
  # Parameters -----------------------------------------
  K     <- N * scale ^ 2     # Number of trees
  gmax  <- 800               # Max grass biomass in g m-2
  alpha <- 0.18              # Grass biomass reduction due to tree shading
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
  # Flag to indicate if a tree is resprouting following topkill or not
  res  <- rep(FALSE, K)
  # Initialize with mature trees and allow to reach steady state
  S    <- rep(20, K)
  Smax <- S
  x    <- runif(K, min = 0, max = domain)
  y    <- runif(K, min = 0, max = domain)
  # Growth is location-dependent, unlike the default case
  lambda <- rep(0.23, K)
  lambda[x >= domain / 2] <- 0.23 * env
  lambdar <- rep(0.61, K)
  lambdar[x >= domain / 2] <- 0.61 * env
  # Set baseline growth rate for each tree (deltas.nres, i.e. not a resprout)
  Sdev <- rnorm(K, mean = 0, sd = sigma.ind)
  Sdev <- ifelse(lambda + Sdev < 0, -lambda, Sdev)
  # Preserve Sdev to apply to future topkilled trees
  
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
    # Reset growth rates, since some array positions are now in new locations
    lambda <- rep(0.23, K)
    lambda[x >= domain / 2] <- 0.23 * env
    lambdar <- rep(0.61, K)
    lambdar[x >= domain / 2] <- 0.61 * env
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
Pfire <- round(seq(from = 0.0, to = 0.8, by = 0.2), 1) # round() so Pfire == 0.6 etc. match exactly
env <- seq(from = 1.0, to = 2.0, by = 0.1)

# Run model through driving variable scenarios
ITER <- length(Pfire) * length(env) * 10

df <- data.frame(Iter = 1:ITER)
df$Pfire <- rep(Pfire, each = length(env) * 10)
df$env <- rep(env, each = 10, times = length(Pfire))
df$Rep <- rep(1:10, times = length(Pfire) * length(env))

# Set up cluster
numCores <- detectCores()
cl <- makeCluster(numCores)

# Export everything each worker needs
clusterExport(cl, varlist = c(fire_exports, "tree_mod_env", "scale", "Tmax", "N", "df"))

# Compile C++ on every worker (path must be absolute)
cpp_path <- file.path(getwd(), "spatial_ops.cpp")
clusterExport(cl, varlist = "cpp_path")
clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp(cpp_path)   # each worker compiles once and caches
})

df.sim <- parLapply(cl, 1:ITER, function(i) {
  tree_mod_env(df$Pfire[i], df$env[i])
})

stopCluster(cl)

df.sim <- do.call(rbind, df.sim)
Trun <- Tmax - 500
df.sim$Iter <- rep(1:ITER, each = Trun)
df.sim$Year <- rep(501:Tmax, times = ITER)
boun.env <- inner_join(df, df.sim, by = "Iter")

# Fig. S13: mean +/- SD basal-area difference between the two sides of the
# arena over the last 1000 y of each run, as a function of the growth-rate
# ratio (env), for three representative fire frequencies. binit.mean (the
# dotted reference line) is the initial right-side ("woody") basal area at
# year 501, for comparison against the persisting/eroded contrast later on.
boun.init <- dplyr::filter(boun.env, Year == 501)
binit.mean <- mean(boun.init$B.cl)
# Analyze mean difference between left and right in last 100 years
boun.fin <- dplyr::filter(boun.env, Year > 4500)
boun.fin <- boun.fin %>% mutate(B.diff = B.cl - B.op)
boun.fin <- boun.fin %>% group_by(Pfire, env, Rep) %>% 
  summarize(B.diff = mean(B.diff))
boun.fin <- boun.fin %>% group_by(Pfire, env) %>% 
  summarize(B.diff.mn = mean(B.diff), B.diff.sd = sd(B.diff))

boun.fin <- boun.fin %>% dplyr::filter(Pfire == 0.0 |
                                         Pfire == 0.4 | Pfire == 0.8)
boun.fin <- boun.fin %>% mutate(Pfire = as.factor(Pfire))

# Plot env on x-axis with different lines for Pfire
figS13 <- ggplot(boun.fin, aes(x = env, y = B.diff.mn,
                                 group = Pfire, col = Pfire)) +
  geom_point(size = 3, position = position_dodge(width = 0.05)) +
  geom_errorbar(aes(x = env, ymin = B.diff.mn - B.diff.sd,
                    ymax = B.diff.mn + B.diff.sd),
                position = position_dodge(width = 0.05), width = 0.05, size = 1) +
  geom_hline(yintercept = binit.mean, linetype = 3) +
  geom_hline(yintercept = 0) +
  labs(x = 'Growth ratio',
       y = expression(paste(italic('B'),' difference (m'^2,' ha'^-1,')'))) +
  xlim(0.95,2.05) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.5)),
        axis.title.x = element_text(vjust = 0, size = rel(1.5)),
        axis.text = element_text(size = rel(1.2)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2), hjust = 0.5),
        legend.text = element_text(size = rel(1.2)))

# Save figure
ggsave("FigS13.png", figS13, width = 8, height = 5, dpi = 300)