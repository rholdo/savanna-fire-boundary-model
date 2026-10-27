# Tree-grass-fire spatial simulation
# Run on a 300 x 300 torus with N trees/ha
# Examination of woody-grassy boundary dynamics
# Across a wider range of parameter values
# Parameter uncertainty is informed by variability in
# the empirical data used for parameterization

# This script produces fig. S4

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
NP <- 15 # Number of parameters in uncertainty analysis

# Input bootstrapped parameter values
p1 <- read.csv('Boot_pars_1.csv')
p2 <- read.csv('Boot_pars_2.csv')
p3 <- read.csv('Boot_pars_3.csv')
p4 <- read.csv('Boot_pars_4.csv')

# Shared fire-ignition helpers (see fire_functions.R)
source("fire_functions.R")
# Fire implementation: in each fire year one ignition point
# is drawn at random among trees in the LEFT HALF of the domain with a
# potential ignition and Pignite >= 0.5 (i.e. in open, grassy areas); fire
# then spreads through the full-domain 5-m lattice (see fire_functions.R).
# Edit pign_min / n_seeds / cell / connectivity here to change it.
FIRE <- FIRE_BOUNDARY

# getpars(): draws one random parameter set (length NP = 15) for a single
# call to tree_mod(), to propagate empirical parameter uncertainty through
# the model (manuscript Methods 2.3). Groups of parameters that were fitted
# together (and so are statistically correlated) are resampled as a single
# bootstrapped row from the corresponding Boot_pars_*.csv, keeping those
# correlations intact; lambda and gamma, which are means across studies,
# are instead drawn from independent normal distributions using their
# reported mean +/- SD.
#
# pvals index -> tree_mod() parameter (see that function's Parameters block):
#   1, 2      -- gmax, alpha            (bootstrapped together, p1)
#   3         -- lambda                 (independent normal draw)
#   4, 5      -- lambdar, shs           (bootstrapped together, p3)
#   6         -- gamma                  (independent normal draw)
#   7, 8      -- a0, a1                 (bootstrapped together, p2)
#   9, 10     -- r, w                   (bootstrapped together, p1)
#   11, 12    -- i0, i1                 (bootstrapped together, p1)
#   13, 14, 15 -- t0, t1, t2            (bootstrapped together, p4)
getpars <- function(){
  pvals <- numeric(NP)
  r1 <- p1[sample(nrow(p1),1),]
  r2 <- p2[sample(nrow(p2),1),]
  r3 <- p3[sample(nrow(p3),1),]
  r4 <- p4[sample(nrow(p4),1),]
  pvals[c(1,2,9:12)] <- as.numeric(r1)
  pvals[7:8] <- as.numeric(r2)
  pvals[4:5] <- as.numeric(r3)
  pvals[13:15] <- as.numeric(r4)
  pvals[3] <- rnorm(1, 0.23, 0.06)
  pvals[6] <- rnorm(1, 0.72, 0.18)
  return(pvals)
}

# tree_mod(): identical boundary-persistence simulation to Boundary.R's
# tree_mod(), except that the parameters listed under "Parameters" below
# are supplied via `pars` (see getpars(), above) rather than fixed at their
# default values -- this is what propagates parameter uncertainty through
# to the boundary-persistence outcome. Output columns (B.op, B.cl, D.break,
# B.diff.dyn) are as described for Boundary.R's tree_mod().
#   Pfire, z, sigma.ind, sigma.res, mu -- as in Boundary.R
#   pars -- length-15 vector of parameter values from getpars()
tree_mod <- function(Pfire, z, sigma.ind, sigma.res, mu, pars){

  # Parameters -----------------------------------------
  K     <- N * scale ^ 2      # Number of trees
  gmax  <- pars[1]            # Max grass biomass in g m-2
  alpha <- pars[2]            # Grass biomass reduction due to tree shading
  lambda    <- pars[3]        # Mean tree annual growth rate (cm)
  lambdar <- pars[4]          # Resprout growth rate
  shs       <- pars[5]        # Resprout growth half-saturation constant
  rho       <- 0.05           # Tree crowding parameter
  gamma     <- pars[6]        # Tree growth reduction due to grass
  a0 <- pars[7]               # Height allometry intercept
  a1 <- pars[8]               # Height allometry slope
  RHref <- 20                 # Reference RH
  r     <-  pars[9]           # RH increase per unit basal area
  Wref  <-  5                 # Reference wind speed
  w     <-  pars[10]          # Wind speed decay with basal area
  Mref  <- 40                 # Fuel moisture
  i0 <- pars[11]              # Ignition parameters
  i1 <- pars[12]
  t0 <- pars[13]              # Topkill parameters
  t1 <- pars[14]
  t2 <- pars[15]
  m0 <- 2729                  # Fire intensity parameters (Higgins et al. 2000)
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
  length(mu) * 100

# Generate parameter set
df <- data.frame(Iter = 1:ITER)
df$Pfire <- rep(Pfire, each = length(z) * length(sigma.ind) *
                  length(sigma.res) * length(mu) * 100)
df$z <- rep(z, each = length(sigma.ind) * length(sigma.res) * length(mu) * 100,
           times = length(Pfire))
df$sigma.ind <- rep(sigma.ind, each = length(sigma.res) * length(mu) * 100,
                    times = length(Pfire) * length(z))
df$sigma.res <- rep(sigma.res, each = length(mu) * 100,
                    times = length(Pfire) * length(z) * length(sigma.ind))
df$mu <- rep(mu, each = 100,
            times = length(Pfire) * length(z) * length(sigma.ind) *
              length(sigma.res))
df$Rep <- rep(1:100,
             times = length(Pfire) * length(z) * length(sigma.ind) *
               length(sigma.res) * length(mu))

# Set up cluster
numCores <- detectCores()
cl <- makeCluster(numCores)

# Export everything each worker needs
clusterExport(cl, varlist = c(fire_exports, "tree_mod", "scale", "Tmax", "N", "df",
                              "getpars", "NP", "p1", "p2", "p3", "p4"))

# Compile C++ on every worker (path must be absolute)
cpp_path <- file.path(getwd(), "spatial_ops.cpp")
clusterExport(cl, varlist = "cpp_path")
clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp(cpp_path)   # each worker compiles once and caches
})

df.sim <- parLapply(cl, 1:ITER, function(i) {
  tree_mod(df$Pfire[i], df$z[i], df$sigma.ind[i], df$sigma.res[i],
           df$mu[i], getpars())
})

stopCluster(cl)

df.sim <- do.call(rbind, df.sim)
Trun <- Tmax - 500
df.sim$Iter <- rep(1:ITER, each = Trun)
df.sim$Year <- rep(501:Tmax, times = ITER)
boun <- inner_join(df, df.sim, by = "Iter")

# Plot results (Fig. S4)
# Same approach as Boundary.R: reference the initial (year 501) basal-area
# contrast, then summarize the mean +/- SD contrast at the fixed boundary
# location over the last 100 y of each run, across parameter combinations
# and the (here, 100 vs. 10) replicate runs -- but here each replicate also
# draws its own random parameter set via getpars(), so the spread across
# replicates reflects both stochastic and parameter uncertainty.
boun.init <- dplyr::filter(boun, Year == 501)
boun.init<- boun.init %>% mutate(Bdiff = B.cl - B.op)
Bdiff.init <- mean(boun.init$Bdiff) # Initial basal area difference (for reference)
# Analyze mean difference between left and right in last 100 years
boun.fin <- dplyr::filter(boun, Year > 4900)
boun.fin <- boun.fin %>% mutate(B.diff = B.cl - B.op)
boun.fin <- boun.fin %>% group_by(Pfire, z, sigma.ind, sigma.res, mu, Rep) %>%
  summarize(B.diff = mean(B.diff))
boun.fin <- boun.fin %>% group_by(Pfire, z, sigma.ind, sigma.res, mu) %>%
  summarize(B.diff.mn = mean(B.diff), B.diff.sd = sd(B.diff))
boun.fin <- boun.fin %>% dplyr::filter(Pfire != 0.3 &
                                         Pfire != 0.5 & Pfire != 0.7)
boun.fin<- boun.fin %>% mutate(mu = as.factor(mu),
                               Pfire = as.factor(Pfire),
                               z = as.factor(z))
boun.fin <- boun.fin %>% dplyr::filter(sigma.res == 0.0)
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

figS4 <- ggplot(boun.fin, aes(x = sigma.ind, y = B.diff.mn,
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

# Save figure
ggsave("FigS4.png", figS4, width = 10, height = 10, dpi = 300)