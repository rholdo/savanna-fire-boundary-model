# Tree-grass-fire spatial simulation
# Run on a 300 x 300 torus with N trees/ha
# This script outputs location/size data for individual trees
# for single runs

# This script generates figures 4, S5, S6, S11 & S12

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

# tree_mod_example(): same boundary-persistence simulation as Boundary.R's
# tree_mod() (mature trees for 500 y, then the left half is reset to
# saplings and fire switched on for the remaining Tmax - 500 y), but here
# the full per-tree state (location, size, and neighborhood grass biomass/
# fire intensity) is recorded every year for a single run, rather than
# just arena-wide summary statistics. This is what lets the figures below
# show individual tree locations/sizes directly instead of only aggregate
# basal area.
#   Pfire     -- exogenous fire probability per year (P_fire)
#   z         -- neighborhood radius (m) for competition and fire effects
#   sigma.ind -- SD of growth rate across individuals (cm y-1)
#   sigma.res -- residual (within-individual) growth SD (cm y-1)
#   mu        -- background (adult) mortality rate (y-1)
#   RESET     -- TRUE (default; Figs. 4, S5, S6):
#                boundary simulation as described above. FALSE (Figs. S11,
#                S12): no boundary, i.e. the same no-boundary model as
#                Boundary_NN.R / Spatial_structure.R: the left half is NOT
#                reset at t = 500, fire is on from year 1, and the ignition
#                point can fall anywhere in the domain (FIRE_CLUSTER, open-
#                area seeding). FIRE_BOUNDARY cannot be used without the
#                reset: it only seeds fire in the left half where
#                Pignite >= 0.5, which no tree satisfies in an initially
#                closed stand, so fire would never start.
# Returns a data frame with one row per tree per year (all Tmax years,
# including the pre-reset burn-in), giving Size, x, y, G (neighborhood
# grass biomass) and Q (neighborhood fire intensity).
tree_mod_example <- function(Pfire, z, sigma.ind, sigma.res, mu, RESET = TRUE){
  fire_cfg <- if (RESET) FIRE else FIRE_CLUSTER
  
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
  sim <- data.frame(Size = S)
  sim$Year <- 0
  sim$x <- x
  sim$y <- y
  sim$G <- NA
  sim$Q <- NA

  for(t in 1:Tmax){

    if (RESET && t == 500){
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
    # Without the reset, fire is on from year 1
    fire <- ifelse(t > 500 || !RESET, runif(1) < Pfire, FALSE)
    # Local ignition probability if fire has spread
    Pignite <- 1 / (1 + i0 * exp(i1 * B))
    ignite <- runif(K) < Pignite
    # Update ignite with connectivity map
    ignite <- fire_ignition(x, y, ignite, Pignite, domain, fire_cfg)
    RH <- RHref + r * B
    W <- Wref * exp(-w * B)
    Q <- m0 + m1 * G - m2 * sqrt(Mref) -
      m3 * RH ^ 2 - m4 / W
    Q <- ifelse(Q < 0, 0, Q)
    # Convert diameter to height
    logH <- a1 * log(S) + a0
    eta <- t0 - t1 * logH + t2 * sqrt(Q)
    Ptopkill <- ignite * exp(eta) / (1 + exp(eta))

    if (fire){
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

    # Collect output
    sim.year <- data.frame(Size = S)
    sim.year$Year <- t
    sim.year$x <- x
    sim.year$y <- y
    sim.year$G <- G
    sim.year$Q <- Q
    sim <- rbind(sim, sim.year)
  }
  return(sim)
}

# Run scenarios
# Five single-run examples, each illustrating a different qualitative
# outcome discussed in the manuscript. Scenarios 1-3 feed the tree-map
# panels (fig4a-c), size-distribution histograms (figS6a-c), and
# basal-area/grass/fire-intensity transects (figS5a-i) below; scenarios 4
# and 5 feed the multi-year cluster-tracking maps (figS11, figS12).
#
# Scenarios 4 and 5 (Figs. S11, S12) run with RESET = FALSE, i.e. the
# no-boundary model (no reset at year 500, fire on from year 1, open-area
# seeding), matching Boundary_NN.R and Fig. 6.

# Scenario 1: low fire -- illustrates boundary erosion under Pfire = 0.2
sim1 <- tree_mod_example(0.2,10,0.15,0.15,0.01)

# Scenario 2: high fire -- illustrates faster boundary erosion/homogenization
# under Pfire = 0.8
sim2 <- tree_mod_example(0.8,10,0.15,0.15,0.01)

# Scenario 3: possible stable trajectory (traveling wave) -- deterministic
# growth (sigma.ind = sigma.res = 0) and a large neighborhood (z = 20),
# conditions expected to best preserve boundary sharpness (if not location)
sim3 <- tree_mod_example(0.8,20,0.0,0.0,0.005)

# Scenario 4: example shown in fig. S11 -- illustrates transient clustering
# of large trees within an otherwise grassy matrix
sim4 <- tree_mod_example(0.6,10,0.15,0.0,0.0075, RESET = FALSE)

# Scenario 5: example shown in fig. S12 -- as above, but illustrating
# transient grassy patches within an otherwise woody matrix
sim5 <- tree_mod_example(0.2,5,0.15,0.0,0.0075, RESET = FALSE)

# Create plots
# --- Scenario 1 (sim1) processing -- this same pattern (relabel to 3
# snapshot years relative to the fire-on reset at year 500, cap tree size
# for the histograms, zoom to the central 100 x 100 m to avoid edge/wrap
# artifacts, then bin x into 10 transects) is repeated below for sim2 and
# sim3. ---
# Keep 3 snapshots: right after the boundary reset (Year 501 = "Year 1"
# post-reset), partway through (Year 700 = "Year 200"), and at the end of
# the run (Year 5000 = "Year 4500").
sim1 <- dplyr::filter(sim1, Year == 501 | Year == 700 | Year == 5000)
sim1 <- sim1 %>% mutate(
  Year2 = case_match(
    Year,
    501 ~ "Year 1",
    700 ~ "Year 200",
    5000 ~ "Year 4500"
  )
)

# Cap displayed size at 30 cm so the size-distribution histograms (figS6)
# aren't dominated by a long tail of a few very large trees.
sim.size <- sim1 %>% mutate(Size_cap = pmin(Size, 30))
# Zoom in on central 100 x 100 m sector (avoids the toroidal wraparound
# edge showing up as a visual artifact in the tree-map panels)
sim1 <- dplyr::filter(sim1, (x > 100 & x < 200) & (y > 100 & y < 200))
# Divide the zoomed window into 10 x-direction transects for the
# basal-area/grass/fire-intensity profile figures (figS5) below
sim1 <- sim1 %>% mutate(x_bin = ntile(x, n = 10))

# Low fire example
fig4a <- ggplot(sim1, aes(x = x, y = y, size = Size)) +
  geom_point(pch = 21, alpha = 0.8, color = "black", fill = "green") +
  geom_vline(xintercept = 150, linetype = "dashed", linewidth = 1) +
  facet_wrap(~ factor(Year2, c('Year 1', 'Year 200', 'Year 4500')),
             nrow = 1) +
  scale_size_area(max_size = 4, 
                  breaks = c(1, 10, 25, 50, 100),
                  name = "Size (cm)") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

figS6a <- ggplot(sim.size, aes(x = Size_cap)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  facet_wrap(~ factor(Year2), nrow = 1) +
  labs(x = "Size (cm)", y = "Count") +
  coord_cartesian(xlim = c(0, 30)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

# Plot change in variables over time across x transects
# Average grass biomass and fire intensity (both already neighborhood-level
# quantities per tree), and sum tree basal area, within each 10-m-wide x
# transect and snapshot year -- these become the x-axis profiles in figS5,
# used to visualize how sharply the woody/grassy contrast (and its drivers)
# vary across the boundary at each time point.
sim.ag <- sim1 %>% group_by(Year, x_bin) %>%
  summarise(Grass = mean(G), Energy = mean(Q),
            Basal = sum((Size / 2) ^ 2 * pi) / 1e4 * 10)
sim.ag$x <- sim.ag$x_bin * 10 + 95
sim.ag <- sim.ag %>% mutate(
  Year = case_match(
    Year,
    501 ~ "1",
    700 ~ "200",
    5000 ~ "4500"
  )
)
sim.ag <- sim.ag %>%
  mutate(Year = factor(Year,
                       levels = c("1", "200", "4500")))

figS5a <- ggplot(sim.ag, aes(x = x, y = Basal, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Basal area (m'^2,'ha'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))
figS5b <- ggplot(sim.ag, aes(x = x, y = Grass, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Grass biomass (g m'^-2,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))
figS5c <- ggplot(sim.ag, aes(x = x, y = Energy, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Fire intensity (kJ s'^-1,' m'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))

# --- Scenario 2 (sim2): same processing pattern as sim1 above ---
sim2 <- dplyr::filter(sim2, Year == 501 | Year == 700 | Year == 5000)
sim2 <- sim2 %>% mutate(
  Year2 = case_match(
    Year,
    501 ~ "Year 1",
    700 ~ "Year 200",
    5000 ~ "Year 4500"
  )
)

sim.size <- sim2 %>% mutate(Size_cap = pmin(Size, 30))
# Zoom in on central 100 x 100 m sector
sim2 <- dplyr::filter(sim2, (x > 100 & x < 200) & (y > 100 & y < 200))
sim2 <- sim2 %>% mutate(x_bin = ntile(x, n = 10))

# High fire example
fig4b <- ggplot(sim2, aes(x = x, y = y, size = Size)) +
  geom_point(pch = 21, alpha = 0.8, color = "black", fill = "green") +
  geom_vline(xintercept = 150, linetype = "dashed", linewidth = 1) +
  facet_wrap(~ factor(Year2, c('Year 1', 'Year 200', 'Year 4500')),
             nrow = 1) +
  scale_size_area(max_size = 4, 
                  breaks = c(1, 10, 25, 50, 100),
                  name = "Size (cm)") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

figS6b <- ggplot(sim.size, aes(x = Size_cap)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  facet_wrap(~ factor(Year2), nrow = 1) +
  labs(x = "Size (cm)", y = "Count") +
  coord_cartesian(xlim = c(0, 30)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

# Plot change in variables over time across x transects
sim.ag <- sim2 %>% group_by(Year, x_bin) %>% 
  summarise(Grass = mean(G), Energy = mean(Q),
            Basal = sum((Size / 2) ^ 2 * pi) / 1e4 * 10)
sim.ag$x <- sim.ag$x_bin * 10 + 95
sim.ag <- sim.ag %>% mutate(
  Year = case_match(
    Year,
    501 ~ "1",
    700 ~ "200",
    5000 ~ "4500"
  )
)
sim.ag <- sim.ag %>%
  mutate(Year = factor(Year,
                       levels = c("1", "200", "4500")))

figS5d <- ggplot(sim.ag, aes(x = x, y = Basal, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Basal area (m'^2,'ha'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))
figS5e <- ggplot(sim.ag, aes(x = x, y = Grass, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Grass biomass (g m'^-2,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))
figS5f <- ggplot(sim.ag, aes(x = x, y = Energy, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Fire intensity (kJ s'^-1,' m'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))

# --- Scenario 3 (sim3): same processing pattern as sim1 above ---
sim3 <- dplyr::filter(sim3, Year == 501 | Year == 700 | Year == 5000)
sim3 <- sim3 %>% mutate(
  Year2 = case_match(
    Year,
    501 ~ "Year 1",
    700 ~ "Year 200",
    5000 ~ "Year 4500"
  )
)

sim.size <- sim3 %>% mutate(Size_cap = pmin(Size, 30))
# Zoom in on central 100 x 100 m sector
sim3 <- dplyr::filter(sim3, (x > 100 & x < 200) & (y > 100 & y < 200))
sim3 <- sim3 %>% mutate(x_bin = ntile(x, n = 10))

# Boundary persistence example
fig4c <- ggplot(sim3, aes(x = x, y = y, size = Size)) +
  geom_point(pch = 21, alpha = 0.8, color = "black", fill = "green") +
  geom_vline(xintercept = 150, linetype = "dashed", linewidth = 1) +
  facet_wrap(~ factor(Year2, c('Year 1', 'Year 200', 'Year 4500')),
             nrow = 1) +
  scale_size_area(max_size = 4, 
                  breaks = c(1, 10, 25, 50, 100),
                  name = "Size (cm)") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

figS6c <- ggplot(sim.size, aes(x = Size_cap)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  facet_wrap(~ factor(Year2), nrow = 1) +
  labs(x = "Size (cm)", y = "Count") +
  coord_cartesian(xlim = c(0, 30)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

# Plot change in variables over time across x transects
sim.ag <- sim3 %>% group_by(Year, x_bin) %>% 
  summarise(Grass = mean(G), Energy = mean(Q),
            Basal = sum((Size / 2) ^ 2 * pi) / 1e4 * 10)
sim.ag$x <- sim.ag$x_bin * 10 + 95
sim.ag <- sim.ag %>% mutate(
  Year = case_match(
    Year,
    501 ~ "1",
    700 ~ "200",
    5000 ~ "4500"
  )
)
sim.ag <- sim.ag %>%
  mutate(Year = factor(Year,
                       levels = c("1", "200", "4500")))

figS5g <- ggplot(sim.ag, aes(x = x, y = Basal, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Basal area (m'^2,'ha'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))
figS5h <- ggplot(sim.ag, aes(x = x, y = Grass, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Grass biomass (g m'^-2,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))
figS5i <- ggplot(sim.ag, aes(x = x, y = Energy, group = Year, col = Year)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 150, linetype = "dashed") +
  labs(x = 'x', 
       y = expression(paste('Fire intensity (kJ s'^-1,' m'^-1,')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)))

# Fig. 4: combine the three tree-map panels (fig4a = low fire, fig4b =
# high fire, fig4c = deterministic-growth/large-neighborhood) into one
# stacked figure, each showing tree size/location at 3 snapshot years
# relative to the boundary reset.
fig4 <- plot_grid(
  fig4a,
  fig4b,
  fig4c,
  align = 'vh',
  labels = c('a)','b)','c)'),
  hjust = 0,
  nrow = 3
)

prow <- plot_grid(
  figS5a + theme(legend.position = "none",
                 plot.margin = margin(t = 20, r = 5, b = 5, l = 5)),
  figS5b + theme(legend.position = "none",
                 plot.margin = margin(t = 20, r = 5, b = 5, l = 5)),
  figS5c + theme(legend.position = "none",
                 plot.margin = margin(t = 20, r = 5, b = 5, l = 5)),
  figS5d + theme(legend.position = "none"),
  figS5e + theme(legend.position = "none"),
  figS5f + theme(legend.position = "none"),
  figS5g + theme(legend.position = "none"),
  figS5h + theme(legend.position = "none"),
  figS5i + theme(legend.position = "none"),
  align = 'vh',
  labels = c('a)','b)','c)','d)','e)','f)','g)','h)','i)'),
  hjust = -1,
  label_y = 1.03,
  nrow = 3
)
legend <- get_legend(
  # create some space to the left of the legend
  figS5a + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# Fig. S5: 3x3 grid of basal-area/grass-biomass/fire-intensity transects
# (rows a-c = scenario 1, d-f = scenario 2, g-i = scenario 3), each row
# showing all 3 snapshot years overlaid, giving the profile-view counterpart
# to the tree maps in fig4.
figS5 <- plot_grid(prow, legend, rel_widths = c(3, .55))

# Fig. S6: tree size-distribution histograms (a = scenario 1, b = scenario
# 2, c = scenario 3), one row per scenario with all 3 snapshot years faceted
figS6 <- plot_grid(
  figS6a,
  figS6b,
  figS6c,
  align = 'vh',
  labels = c('a)','b)','c)'),
  hjust = 0,
  nrow = 3
)

sim4 <- dplyr::filter(sim4, (x > 100 & x < 200) & (y > 100 & y < 200))

# Fig. S11: track the clustering dynamics of a single site over 500-y
# intervals (years 0, 500, ..., 4000), from scenario 4 (moderate fire).
# Just show the central 100 x 100 m. Illustrates transient clusters of
# large trees forming and dissolving within an otherwise grassy matrix.
sim4 <- sim4 %>% dplyr::filter(Year %% 500 == 0 & Year <= 4000)
figS11 <- ggplot(sim4, aes(x = x, y = y, size = Size)) +
  geom_point(pch = 21, alpha = 0.8, color = "black", fill = "green") +
  facet_wrap(~ factor(Year), nrow = 3) +
  scale_size_area(max_size = 4, 
                  breaks = c(1, 10, 25, 50, 100, 140),
                  name = "Size (cm)") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

# Fig. S12: same approach as figS11, but from scenario 5 (low fire, small
# neighborhood), illustrating transient grassy patches forming and being
# filled in within an otherwise woody matrix.
sim5 <- dplyr::filter(sim5, (x > 100 & x < 200) & (y > 100 & y < 200))
sim5 <- sim5 %>% dplyr::filter(Year %% 500 == 0 & Year <= 4000)
figS12 <- ggplot(sim5, aes(x = x, y = y, size = Size)) +
  geom_point(pch = 21, alpha = 0.8, color = "black", fill = "green") +
  facet_wrap(~ factor(Year), nrow = 3) +
  scale_size_area(max_size = 4, 
                  breaks = c(1, 10, 25, 50, 100, 140),
                  name = "Size (cm)") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

# Save figures
ggsave("Fig4.png", fig4, width = 10, height = 10, dpi = 300)
ggsave("FigS5.png", figS5, width = 12, height = 10, dpi = 300)
ggsave("FigS6.png", figS6, width = 10, height = 10, dpi = 300)
ggsave("FigS11.png", figS11, width = 10, height = 10, dpi = 300)
ggsave("FigS12.png", figS12, width = 10, height = 10, dpi = 300)