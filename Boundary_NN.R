# Tree-grass-fire spatial simulation
# Run on a 300 x 300 torus with N trees/ha
# Generates nearest-neighbor outcomes for two tree size classes

# This script produces figures 5, 6 and S10

library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)
library(tidyr)
library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)

# Compile the C++ spatial operations (do this once at startup)
# spatial_ops.cpp must be in the same directory as this script,
# or provide the full path below.
sourceCpp("spatial_ops.cpp")

Tmax  <- 5000 # Duration (years)
scale <- 3    # Scale of simulation domain (x100 m)
N <- 400      # Stem density (ha-1); 400 is default, 1000 also tested

# Shared fire-ignition helpers (see fire_functions.R)
source("fire_functions.R")
# Fire implementation: open-area seeding. In each fire year
# one ignition point is drawn at random among trees with a potential
# ignition anywhere in the domain, so the chance of a location being the
# ignition point is proportional to its Pignite (fires start mostly in open
# areas). Fire then spreads through the full-domain 5-m lattice, so no
# left-right burning gradient is imposed on the homogeneous template. See
# fire_functions.R.
FIRE <- FIRE_CLUSTER

# tree_mod(): tracks transient spatial clustering rather than an imposed
# boundary. Trees are initialized as a mature (S = 20 cm), randomly-located
# population and simulated for Tmax years with fire on throughout; in each
# year, the mean nearest-neighbor distance is recorded for two size classes
# (S > 5 cm and S > 10 cm) via toroidal_nn().
#   Pfire     -- exogenous fire probability per year (P_fire)
#   z         -- neighborhood radius (m) for competition and fire effects
#   sigma.ind -- SD of growth rate across individuals (cm y-1)
#   sigma.res -- residual (within-individual) growth SD (cm y-1)
#   mu        -- background (adult) mortality rate (y-1)
# Returns a data frame with one row per year:
#   D5, D10   -- mean nearest-neighbor distance (m) among trees with
#                S > 5 cm / S > 10 cm, respectively (NA if fewer than 2
#                qualifying trees)
#   N5, N10   -- number of qualifying trees in each size class that year
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
  
  # Output
  out <- data.frame(D5  = numeric(Tmax), N5  = numeric(Tmax),
                    D10 = numeric(Tmax), N10 = numeric(Tmax))
  
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
    
    # ------------------------------------------------------------------
    # Nearest-neighbor output (C++)
    # toroidal_nn() uses grid-based spatial indexing with toroidal distances.
    # toroidal_nn() returns the exact nearest-neighbor distance for every
    # tree, however sparse the size class.
    # ------------------------------------------------------------------
    
    # --- S > 5 size class ---
    idx5 <- which(S > 5)
    n5   <- length(idx5)
    if (n5 > 1) {
      dmin5       <- toroidal_nn(x[idx5], y[idx5], domain, nn_cell_size)
      out$D5[t]   <- mean(dmin5[is.finite(dmin5)])
      out$N5[t]   <- n5
    } else {
      out$D5[t] <- NA
      out$N5[t] <- n5
    }

    # --- S > 10 size class ---
    idx10 <- which(S > 10)
    n10   <- length(idx10)
    if (n10 > 1) {
      dmin10      <- toroidal_nn(x[idx10], y[idx10], domain, nn_cell_size)
      out$D10[t]  <- mean(dmin10[is.finite(dmin10)])
      out$N10[t]  <- n10
    } else {
      out$D10[t] <- NA
      out$N10[t] <- n10
    }
    
  } # end time loop
  
  return(out)
}

# Run model through scenarios
Pfire     <- round(seq(from = 0.2, to = 0.8, by = 0.1), 1) # round() so Pfire == 0.6 etc. match exactly
z         <- c(5, 10, 20)
sigma.ind <- seq(from = 0.0, to = 0.3, by = 0.15)
sigma.res <- c(0, 0.15)
mu        <- seq(from = 0.005, to = 0.01, by = 0.0025)

ITER <- length(Pfire) * length(z) * length(sigma.ind) * length(sigma.res) *
  length(mu) * 10

df <- data.frame(Iter = 1:ITER)
df$Pfire     <- rep(Pfire, each = length(z) * length(sigma.ind) *
                      length(sigma.res) * length(mu) * 10)
df$z         <- rep(z, each = length(sigma.ind) * length(sigma.res) *
                      length(mu) * 10,
                    times = length(Pfire))
df$sigma.ind <- rep(sigma.ind, each = length(sigma.res) * length(mu) * 10,
                    times = length(Pfire) * length(z))
df$sigma.res <- rep(sigma.res, each = length(mu) * 10,
                    times = length(Pfire) * length(z) * length(sigma.ind))
df$mu        <- rep(mu, each = 10,
                    times = length(Pfire) * length(z) * length(sigma.ind) *
                      length(sigma.res))
df$Rep       <- rep(1:10,
                    times = length(Pfire) * length(z) * length(sigma.ind) *
                      length(sigma.res) * length(mu))

# Set up cluster
numCores <- detectCores()
cl <- makeCluster(numCores)

# Export everything each worker needs
clusterExport(cl, varlist = c(fire_exports, "tree_mod",
                              "scale", "Tmax", "N", "df"))

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

df.sim$Iter <- rep(1:ITER, each = Tmax)
df.sim$Year <- rep(1:Tmax, times = ITER)
df <- inner_join(df, df.sim, by = "Iter")

# Plot results

# Expected nearest-neighbor distance under complete spatial randomness (CSR),
# i.e. the benchmark against which observed clustering is judged, following
# Clark & Evans (1954): for N points scattered uniformly at random over area
# A, the expected nearest-neighbor distance is 1 / (2*sqrt(N/A)) with
# standard error sqrt((4-pi) / (4*pi*N*(N/A))). D.med/D.lo/D.hi give the
# expected distance and its 95% CI for every observed value of N5/N10, so
# each simulated year can be compared against the CSR expectation for that
# year's actual tree count.
# N.obs is the observed tree count (distinct from the stem-density setting
# N), and A is computed from scale so the benchmark matches the domain.
N.obs <- unique(c(unique(df$N5),unique(df$N10)))
bound <- data.frame(N = N.obs)

A <- (100 * scale) ^ 2  # domain area in m²
bound$D.med <- 1 / (2 * sqrt(bound$N / A))
D.se  <- sqrt((4 - pi) / (4 * pi * bound$N * (bound$N / A)))
bound$D.lo  <- bound$D.med - 1.96 * D.se
bound$D.hi  <- bound$D.med + 1.96 * D.se

# Attach the CSR benchmark to each simulated year (matched on that year's
# tree count), then flag years as significantly clustered (*.sig = TRUE)
# when the observed mean nearest-neighbor distance falls below the lower
# 95% CI of the CSR expectation, i.e. trees are closer together than
# expected under random placement. *.diff is the gap between the CSR
# expectation and the observed distance (used in fig6 to show clustering
# strength over time).
df <- left_join(df, bound, by = c("N5" = "N"))
df <- left_join(df, bound, by = c("N10" = "N"))
names(df)[13:18] <- c('D5.med','D5.lo','D5.hi','D10.med','D10.lo','D10.hi')
df <- df %>% mutate(
  D5.diff = D5.med - D5,
  D5.sig = ifelse(D5 < D5.lo, TRUE, FALSE),
  D10.diff = D10.med - D10,
  D10.sig = ifelse(D10 < D10.lo, TRUE, FALSE)
)
# Restrict datasets to cases with 20 or more trees and last 500 years
# (sufficient trees for adequate power to detect clustering, and burned-in
# so early transients don't bias the estimate; see manuscript Methods 2.4)
df5 <- dplyr::filter(df, N5 >= 20 & Year > Tmax - 500)
df10 <- dplyr::filter(df, N10 >= 20 & Year > Tmax - 500)

# Find prevalence of large tree clustering across iterations
# Agg = % of the last 500 y that were flagged as significantly clustered
# (D.sig == TRUE), first averaged within each of the 10 replicate runs
# (Agg), then across replicates for each parameter combination (Agg.mn/sd).
clus5 <- df5 %>% group_by(Pfire, z, sigma.ind, sigma.res, mu, Rep) %>%
  summarize(Agg = mean(D5.sig, na.rm = TRUE) * 100)
clus5 <- clus5 %>% group_by(Pfire, z, sigma.ind, sigma.res, mu) %>%
  summarize(Agg.mn = mean(Agg), Agg.sd = sd(Agg))
# Reserve for later: mean clustering prevalence by fire frequency alone,
# averaged across all other parameter combinations (used in figS10)
clus5.ag <- clus5 %>% group_by(Pfire) %>%
  summarize(Agg.fire = mean(Agg.mn, na.rm = TRUE))

clus5 <- clus5 %>% dplyr::filter(Pfire == 0.2 |
                                   Pfire == 0.4 | Pfire == 0.6 | Pfire == 0.8)
clus5<- clus5 %>% mutate(mu = as.factor(mu),
                         Pfire = as.factor(Pfire),
                         z = as.factor(z))
clus5 <- clus5 %>% dplyr::filter(sigma.res == 0.0)
clus5 <- clus5 %>% mutate(
  z = case_match(
    as.character(z), # as.character() avoids a factor-cast error when a level has no data
    '5' ~ 'z = 5',
    '10' ~ 'z = 10',
    '20' ~ 'z = 20'
  )
)
clus5 <- clus5 %>% mutate(
  mu = case_match(
    as.character(mu), # as.character() avoids a factor-cast error when a level has no data
    '0.005' ~ '0.0050',
    '0.0075' ~ '0.0075',
    '0.01' ~ '0.0100'
  )
)
clus5 <- clus5 %>%
  mutate(z = factor(z, levels = c('z = 5', 'z = 10', 'z = 20')))
clus5 <- clus5 %>% mutate(
  Pfire = case_match(
    as.character(Pfire), # as.character() avoids a factor-cast error when a level has no data
    '0.2' ~ 'Pfire = 0.2',
    '0.4' ~ 'Pfire = 0.4',
    '0.6' ~ 'Pfire = 0.6',
    '0.8' ~ 'Pfire = 0.8'
  )
)

# Clustering of trees > 5 cm could be plotted here (not shown)

# Same clustering-prevalence calculation as clus5 above, for the > 10 cm
# size class -- this is the class actually plotted in fig5.
clus10 <- df10 %>% group_by(Pfire, z, sigma.ind, sigma.res, mu, Rep) %>%
  summarize(Agg = mean(D10.sig, na.rm = TRUE) * 100)
clus10 <- clus10 %>% group_by(Pfire, z, sigma.ind, sigma.res, mu) %>% 
  summarize(Agg.mn = mean(Agg), Agg.sd = sd(Agg))
# Reserve overall fig for later
clus10.ag <- clus10 %>% group_by(Pfire) %>% 
  summarize(Agg.fire = mean(Agg.mn, na.rm = TRUE))

clus10 <- clus10 %>% dplyr::filter(Pfire == 0.2 |
                                     Pfire == 0.4 | Pfire == 0.6 | Pfire == 0.8)
clus10<- clus10 %>% mutate(mu = as.factor(mu),
                           Pfire = as.factor(Pfire),
                           z = as.factor(z))
clus10 <- clus10 %>% dplyr::filter(sigma.res == 0.0)
clus10 <- clus10 %>% mutate(
  z = case_match(
    as.character(z), # as.character() avoids a factor-cast error when a level has no data
    '5' ~ 'z = 5',
    '10' ~ 'z = 10',
    '20' ~ 'z = 20'
  )
)
clus10 <- clus10 %>% mutate(
  mu = case_match(
    as.character(mu), # as.character() avoids a factor-cast error when a level has no data
    '0.005' ~ '0.0050',
    '0.0075' ~ '0.0075',
    '0.01' ~ '0.0100'
  )
)
clus10 <- clus10 %>%
  mutate(z = factor(z, levels = c('z = 5', 'z = 10', 'z = 20')))
clus10 <- clus10 %>% mutate(
  Pfire = case_match(
    as.character(Pfire), # as.character() avoids a factor-cast error when a level has no data
    '0.2' ~ 'Pfire = 0.2',
    '0.4' ~ 'Pfire = 0.4',
    '0.6' ~ 'Pfire = 0.6',
    '0.8' ~ 'Pfire = 0.8'
  )
)

# Fig. 5: prevalence of significant spatial clustering (% of years, > 10 cm
# trees) as a function of individual growth variability (sigma.ind), fire
# frequency (rows) and neighborhood size (columns), colored by mortality
# rate.
fig5 <- ggplot(clus10, aes(x = sigma.ind, y = Agg.mn,
                                 group = mu, col = mu)) + 
  geom_point(size = 3, position = position_dodge(width = 0.05)) +
  geom_errorbar(aes(x = sigma.ind, ymin = Agg.mn - Agg.sd, ymax = Agg.mn + Agg.sd),
                position = position_dodge2(width = 0.05, preserve = "single"),
                width = 0.05, size = 1) +
  facet_grid(rows = vars(Pfire), cols = vars(z)) +
  labs(x = expression(paste(sigma[ind],' (cm y'^-1,')')), 
       y = 'Spatial clustering (%)') +
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

# Fig. S10: clustering prevalence as a function of fire frequency alone
# (averaged across all other parameter combinations), for both size classes.
fire.ag <- left_join(clus5.ag, clus10.ag, by = 'Pfire')
names(fire.ag)[2:3] <- c('5 cm', '10 cm')
fire.ag <-  fire.ag %>% pivot_longer(cols = 2:3,
                                     names_to = "Size",
                                     values_to = "Agg")

figS10 <- ggplot(fire.ag, aes(x = Pfire, y = Agg, group = Size, col = Size)) +
  geom_line(size = 1.5) +
  labs(x = expression(paste('P'[fire],' (y'^-1,')')), 
       y = 'Spatial clustering (%)') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

# Fig. 6: show transient dynamics of clustering for two specific cases --
# a low-fire/large-neighborhood scenario (dfclus.lo, using the > 10 cm
# clustering signal D10.diff) and a high-fire/small-neighborhood scenario
# (dfclus.hi), both using the > 10 cm class (D10.diff). Each of the 10
# replicate runs is plotted as a separate line, to illustrate that
# individual clusters form and dissolve over time even though clustering
# prevalence overall is stable.
df<- df %>%
  mutate(Rep = factor(Rep))
dfclus.lo <- dplyr::filter(df, Pfire == 0.2 & mu == 0.0075 & sigma.ind == 0.15 &
                             z == 20 & sigma.res == 0)
dfclus.hi <- dplyr::filter(df, Pfire == 0.6 & mu == 0.0075 & sigma.ind == 0.15 &
                             z == 10 & sigma.res == 0)

fig6a <- ggplot(dfclus.lo, aes(x = Year, y = D10.diff, group = Rep)) + 
  geom_line(aes(color = Rep)) +
  labs(x = 'Year', 
       y = expression(paste('ANN'[exp],' - ANN'[obs],' (m)'))) +
  ylim(-1,12) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = "none")
fig6b <- ggplot(dfclus.hi, aes(x = Year, y = D10.diff, group = Rep)) + 
  geom_line(aes(color = Rep)) +
  labs(x = 'Year', 
       y = expression(paste('ANN'[exp],' - ANN'[obs],' (m)'))) +
  ylim(-1,12) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = "none")
fig6 <- plot_grid(
  fig6a,
  fig6b,
  align = 'vh',
  labels = c('a)','b)'),
  hjust = 0,
  nrow = 1
)

# Save figures
ggsave("Fig5.png", fig5, width = 10, height = 10, dpi = 300)
ggsave("Fig6.png", fig6, width = 12, height = 5, dpi = 300)
ggsave("FigS10.png", figS10, width = 8, height = 4, dpi = 300)