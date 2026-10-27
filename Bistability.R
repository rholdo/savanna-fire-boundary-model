# Tree-grass-fire spatial simulation
# Run on a (100 * scale) x (100 * scale) m torus with N trees/ha
# This script runs the simulation over contrasting fire frequencies
# for a range of initial size conditions
# The first case (function tree_mod_mf) assumes an infinite neighborhood size
# The second (function tree_mod) tests a both range of initial size conditions
# and neighborhood sizes
# Fire spreads in both cases, using the same open-area seeding as the
# clustering runs (FIRE_CLUSTER; see below).

# This script generates fig. 1
# For fig. S2, set SIGMA_IND <- 0.15 below

library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)

# Compiled helpers and shared fire functions
sourceCpp("spatial_ops.cpp")
source("fire_functions.R")
# Fire implementation: open-area seeding, as in the clustering runs. In each
# fire year one ignition point is drawn among trees with a potential ignition
# anywhere in the domain (chance proportional to Pignite), and fire spreads
# through the full-domain 5-m lattice.
FIRE <- FIRE_CLUSTER

Tmax  <- 5000 # Duration (years)
scale <- 3    # Scale of simulation domain (x100 m)
# Individual growth variation (sigma.ind, cm/y) used by both model
# functions: 0.0 reproduces Fig. 1; 0.15 reproduces Fig. S2.
SIGMA_IND <- 0.0
N <- 400      # Stem density (ha-1); 400 is default, 1000 also tested

# tree_mod_mf(): nonspatial ("mean-field") version of the simulation, i.e. z = Inf.
# Growth, ignition, and fire intensity all depend on arena-wide total basal
# area (B.tot) rather than a local neighborhood, so there is no spatial
# structure to this version -- it is the baseline against which the
# neighborhood-limited tree_mod() below is compared.
#   sinit -- initial diameter (cm) assigned to every tree at t = 1
#   Pfire -- exogenous fire probability per year (P_fire)
# Returns a data frame with one row per simulated year giving total basal
# area (B), used to classify the long-run "Grassy" vs "Woody" end state.
tree_mod_mf <- function(sinit, Pfire){

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
  # These parameters are kept fixed in these simulations
  sigma.ind <- SIGMA_IND   # set at top of script
  sigma.res <- 0.0
  mu <- 0.005
  
  # Spatial indexing setup -----------------------------------------
  domain    <- 100 * scale
  
  # Initialize model ------------------------------------------
  # Set baseline growth rate for each tree (deltas.nres, i.e. not a resprout)
  Sdev <- rnorm(K, mean = 0, sd = sigma.ind)
  Sdev <- ifelse(lambda + Sdev < 0, -lambda, Sdev)
  # Preserve Sdev to apply to future topkilled trees
  # Flag to indicate if a tree is resprouting following topkill or not
  res  <- rep(FALSE, K)
  S    <- rep(sinit, K)
  Smax <- S
  x    <- runif(K, min = 0, max = domain)
  y    <- runif(K, min = 0, max = domain)


  # Output: one row per simulated year
  Time <- 1:Tmax
  out <- data.frame(Time)
  out$sinit <- sinit
  out$pfire <- Pfire
  out$B <- NA # Total basal area (m2 ha-1) each year
  
  for(t in 1:Tmax){
    
    # Calculate total basal area
    B.tot <- sum(pi * (S / 2) ^ 2 / 1e4) / scale ^ 2 # per ha
    # Grass biomass and competitive effect
    G <- gmax * exp(-alpha * B.tot)
    f <- (1 - G / gmax * (1 - gamma))

    # Growth
    # If a tree is resprouting and less than 5 cm, use faster topkill growth
    # If not, keep current growth rate
    deltas <- ifelse(res,
                     lambdar * Smax / (shs + Smax) + Sdev + rnorm(K, 0, sigma.res),
                     lambda + Sdev + rnorm(K, 0, sigma.res))
    deltas <- ifelse(deltas < 0, 0, deltas)
    S <- S + deltas * f * exp(- rho * B.tot)

    # Resprout status
    res <- ifelse(res, ifelse(S < 2 & S < Smax, TRUE, FALSE), FALSE)
    Smax <- ifelse(S > Smax, S, Smax)

    # Fire
    fire <- runif(1) < Pfire
    # Local ignition probability if fire has spread
    Pignite <- 1 / (1 + i0 * exp(i1 * B.tot))
    ignite <- runif(K) < Pignite
    # Fire spread (Pignite is the same for every tree here)
    ignite <- fire_ignition(x, y, ignite, rep(Pignite, K), domain, FIRE)
    RH <- RHref + r * B.tot
    W <- Wref * exp(-w * B.tot)
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

    # Store output
    out$Time[t] <- t
    out$sinit <- sinit
    out$pfire <- Pfire
    out$B[t] <- B.tot
  }
  return(out)
}

# tree_mod(): spatially-explicit version of the simulation.
# Growth, ignition probability and fire intensity depend on neighborhood
# basal area B, computed for each tree within radius z via compute_ba()
# (see spatial_ops.cpp), rather than on the arena-wide total used by
# tree_mod_mf() above.
#   sinit -- initial diameter (cm) assigned to every tree at t = 1
#   Pfire -- exogenous fire probability per year (P_fire)
#   z     -- neighborhood radius (m) for competition and fire effects
# Returns a data frame with one row per simulated year giving total basal
# area (B), used to classify the long-run "Grassy" vs "Woody" end state.
tree_mod <- function(sinit, Pfire, z){
  
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
  # These parameters are kept fixed in these simulations
  sigma.ind <- SIGMA_IND   # set at top of script
  sigma.res <- 0.0
  mu <- 0.005
  
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
  S    <- rep(sinit, K)
  Smax <- S
  x    <- runif(K, min = 0, max = domain)
  y    <- runif(K, min = 0, max = domain)


  # Output: one row per simulated year
  Time <- 1:Tmax
  out <- data.frame(Time)
  out$sinit <- sinit
  out$pfire <- Pfire
  out$z <- z
  out$B <- NA # Total basal area (m2 ha-1) each year
  
  for(t in 1:Tmax){
    
    # ------------------------------------------------------------------
    # Neighborhood basal area (C++)
    # compute_ba() uses grid-based spatial indexing with toroidal distances.
    # Each tree's own basal area is included in its neighborhood total.
    # ------------------------------------------------------------------
    s2_scaled <- S ^ 2 / (4 * z ^ 2)
    B <- compute_ba(x, y, s2_scaled, z, domain, cell_size)
    # Calculate total basal area
    B.tot <- sum(pi * (S / 2) ^ 2 / 1e4) / scale ^ 2 # per ha
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
    fire <- runif(1) < Pfire
    # Local ignition probability if fire has spread
    Pignite <- 1 / (1 + i0 * exp(i1 * B))
    ignite <- runif(K) < Pignite
    ignite <- fire_ignition(x, y, ignite, Pignite, domain, FIRE) # fire spread
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
    
    # Store output
    out$Time[t] <- t
    out$sinit <- sinit
    out$pfire <- Pfire
    out$z <- z
    out$B[t] <- B.tot
  }
  return(out)
}

# Run model through scenarios (mean-field first)
# Cross every initial size (sinit) with every fire frequency (pfire) and
# run tree_mod_mf() once per combination, in parallel on a PSOCK cluster
# (parLapply), because the fire spread step needs the compiled
# burn_connected() on every worker.
sinit <- seq(from = 1, to = 20, by = 0.1)
pfire<- round(seq(from = 0.2, to = 0.8, by = 0.2), 1) # round() so Pfire == 0.6 etc. match exactly

ITER <- length(sinit) * length(pfire)
df <- data.frame(Iter = 1:ITER)
df$sinit <- rep(sinit, each = length(pfire))
df$pfire <- rep(pfire, times = length(sinit))

numCores <- detectCores()
cl <- makeCluster(numCores)
clusterExport(cl, varlist = c(fire_exports, "tree_mod_mf", "scale", "Tmax",
                              "N", "df", "SIGMA_IND"))
cpp_path <- file.path(getwd(), "spatial_ops.cpp")
clusterExport(cl, varlist = "cpp_path")
clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp(cpp_path)   # each worker compiles once and caches
})
df.sim <- parLapply(cl, 1:ITER, function(i) {
  tree_mod_mf(df$sinit[i], df$pfire[i])
})
stopCluster(cl)
df.sim <- do.call(rbind, df.sim)

# Classify each run's long-run outcome from its final-year basal area,
# using a 5 m2 ha-1 threshold to separate "Grassy" from "Woody" states,
# then attach that classification back onto the full time series for plotting.
df.sim.fin <- df.sim %>% dplyr::filter(Time == max(Time))
df.sim.fin <- df.sim.fin %>%
  mutate(Fstate = ifelse(B < 5, 'Grassy', 'Woody'))
df.sim.fin <- df.sim.fin %>% dplyr::select(sinit,pfire,Fstate)
df.sim.mf <- left_join(df.sim,df.sim.fin)

# Now the spatial version
# Add range of neighborhood sizes
z <- c(5, 10, 20)

# Run model through scenarios
# Same design as above, but crossed with neighborhood radius z; tree_mod()
# also needs the compiled compute_ba() Rcpp function on each worker.
ITER <- length(sinit) * length(pfire) * length(z)
df <- data.frame(Iter = 1:ITER)
df$sinit <- rep(sinit, each = length(pfire) * length(z))
df$pfire <- rep(pfire, each = length(z), times = length(sinit))
df$z <- rep(z, times = length(sinit) * length(pfire))

# Set up cluster
numCores <- detectCores()
cl <- makeCluster(numCores)

# Export everything each worker needs
clusterExport(cl, varlist = c(fire_exports, "tree_mod", "scale", "Tmax", "N", "df",
                              "SIGMA_IND"))

# Compile C++ on every worker (path must be absolute)
cpp_path <- file.path(getwd(), "spatial_ops.cpp")
clusterExport(cl, varlist = "cpp_path")
clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp(cpp_path)   # each worker compiles once and caches
})

df.sim <- parLapply(cl, 1:ITER, function(i) {
  tree_mod(df$sinit[i], df$pfire[i], df$z[i])
})

stopCluster(cl)

df.sim <- do.call(rbind, df.sim)

df.sim.fin <- df.sim %>% dplyr::filter(Time == max(Time))
df.sim.fin <- df.sim.fin %>%
  # Same 5 m2 ha-1 threshold as above
  mutate(Fstate = ifelse(B < 5, 'Grassy', 'Woody'))
df.sim.fin <- df.sim.fin %>% dplyr::select(sinit,pfire,z,Fstate)
df.sim <- left_join(df.sim,df.sim.fin)

# Combine results from the mean-field and spatial versions
# Prepare for plotting: label the mean-field runs as "z = Inf" so they can
# share a facet axis with the finite-z spatial runs, restrict to the four
# representative fire frequencies shown in Fig. 1, and format factor labels.
df.sim.mf$z <- "z = Inf"
df.sim.mf <- df.sim.mf[,c(1:3,6,4,5)]
df.sim  <- df.sim  %>% mutate(
  z = case_match(
    z,
    5 ~ "z = 5 m",
    10 ~ "z = 10 m",
    20 ~ "z = 20 m",
  )
)
df.sim <- rbind(df.sim.mf,df.sim)
df.sim <- df.sim %>% dplyr::filter(pfire == 0.2 | pfire == 0.4 |
                                     pfire == 0.6 | pfire == 0.8)
df.sim$z <- factor(df.sim$z, 
                   levels = c("z = 5 m","z = 10 m","z = 20 m","z = Inf"))
df.sim$pfire <- factor(df.sim$pfire)
df.sim <- df.sim %>% mutate(
  pfire = case_match(
    pfire,
    '0.2' ~ 'Pfire = 0.2',
    '0.4' ~ 'Pfire = 0.4',
    '0.6' ~ 'Pfire = 0.6',
    '0.8' ~ 'Pfire = 0.8'
  )
)

# Plot: basal-area trajectories over time, one panel per fire frequency (rows)
# x neighborhood size (columns), colored by each run's final Grassy/Woody
# state. Bimodal coloring within a panel indicates bistability -- the
# steady state depends on initial tree size (sinit) rather than being
# uniquely determined by Pfire and z. Reproduces Fig. 1.
fig1 <- ggplot(df.sim, aes(x = Time/1e3, y = B, group = sinit, color = Fstate)) +
  geom_line(alpha = 0.3) +
  facet_grid(rows = vars(pfire), cols = vars(z)) +
  labs(x = 'Time (1000 y)', 
       y = expression(paste('Basal area (m'^2,' ha'^-1,')'))) +
  scale_x_continuous(breaks = seq(0, Tmax / 1e3, by = 1)) + # follows Tmax
  scale_colour_manual(name = "Final state", values = c('red','blue')) +
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
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)))

# Save figure: Fig. 1 (SIGMA_IND = 0) or Fig. S2 (SIGMA_IND = 0.15)
fig_name <- if (SIGMA_IND == 0) "Fig1.png" else "FigS2.png"
ggsave(fig_name, fig1, width = 10, height = 7, dpi = 300)