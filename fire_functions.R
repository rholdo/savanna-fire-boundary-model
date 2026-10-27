# fire_functions.R
#
# Shared fire-ignition helpers, sourced by every simulation script.
# Requires spatial_ops.cpp to be compiled (burn_connected()).
#
# Fire-spread lattice
# -------------------
# The domain is divided into square lattice cells of side `cell` (5 m by
# default) over its FULL extent: a 300 x 300 m domain gives a 60 x 60
# lattice, a 100 x 100 m domain a 20 x 20 lattice.
# A cell can carry fire if it contains at least one tree with a potential
# ignition (drawn with probability Pignite), or contains no trees at all;
# cells whose trees all failed to ignite block the fire.
#
# Ignition modes (set with default_fire(mode = ...))
# --------------------------------------------------
# "independent"  Every potential ignition burns; no spread (used only in
#                Fire_robustness.R).
# "seeded"       In each fire year, n_seeds ignition point(s) are drawn at
#                random from the trees that (i) have a potential ignition,
#                (ii) have Pignite >= pign_min and (iii) lie within the band
#                seed_xfrac[1] * domain <= x < seed_xfrac[2] * domain.
#                Fire spreads from the seed cell(s) through connected
#                passable cells; only trees in reached cells burn. With
#                pign_min = 0 the chance of a tree being chosen as a seed is
#                proportional to its Pignite, so fires start preferentially
#                in open areas ("open-area seeding"). If no tree qualifies,
#                nothing burns that year.
#                  Boundary runs:   FIRE_BOUNDARY (left half, Pignite >= 0.5)
#                  Clustering runs: FIRE_CLUSTER  (whole domain, no threshold)
# "edge"         Fire enters along the whole left edge (x = 0) and spreads
#                rightward. Used only for comparison in Fire_robustness.R.
#
# The lattice wraps in x and y for the seeded mode, matching the torus used
# for every other model process; the edge mode never wraps in x.

# Option empty_pass:
#   "always"  (default; used by every main-text/SI simulation): a cell with
#             no trees always carries fire.
#   "pignite" (optional sensitivity test): a cell with no trees carries fire
#             with probability equal to the mean Pignite of the trees in the
#             surrounding 3 x 3 block of cells (i.e. it is as flammable as
#             its surroundings). This gives the spread rule real force in
#             woody areas and tests whether the weak constraint of the
#             default rule matters for the boundary results.
default_fire <- function(mode = "seeded", cell = 5, connectivity = 8,
                         n_seeds = 1, pign_min = 0, seed_xfrac = c(0, 1),
                         wrap_x = (mode != "edge"), wrap_y = TRUE,
                         empty_pass = "always") {
  stopifnot(mode %in% c("independent", "seeded", "edge"),
            connectivity %in% c(4, 8), cell > 0, n_seeds >= 1,
            length(seed_xfrac) == 2, seed_xfrac[1] < seed_xfrac[2],
            empty_pass %in% c("always", "pignite"))
  list(mode = mode, cell = cell, connectivity = connectivity,
       n_seeds = n_seeds, pign_min = pign_min, seed_xfrac = seed_xfrac,
       wrap_x = wrap_x, wrap_y = wrap_y, empty_pass = empty_pass)
}

# Presets used by the scripts. pign_min = 0.5 corresponds to a neighborhood
# basal area of about 8.2 m2 ha-1 (Pignite = 1 / (1 + 0.15 exp(0.23 B))).
FIRE_BOUNDARY <- default_fire(mode = "seeded", pign_min = 0.5,
                              seed_xfrac = c(0, 0.5))
FIRE_CLUSTER  <- default_fire(mode = "seeded", pign_min = 0,
                              seed_xfrac = c(0, 1))

# Lattice helpers ------------------------------------------------------------
lattice_cells <- function(x, y, domain, cell) {
  n <- as.integer(ceiling(domain / cell - 1e-9))   # cells per side
  list(n = n, xc = pmin(floor(x / cell) + 1, n),
       yc = pmin(floor(y / cell) + 1, n))
}
passable_lattice <- function(lc, ig, Pignite = NULL, empty_pass = "always") {
  occupied <- matrix(FALSE, lc$n, lc$n)
  occupied[cbind(lc$yc, lc$xc)] <- TRUE
  has_ig <- matrix(FALSE, lc$n, lc$n)
  if (any(ig)) has_ig[cbind(lc$yc[ig], lc$xc[ig])] <- TRUE
  empty_ok <- !occupied
  if (identical(empty_pass, "pignite")) {
    n <- lc$n
    idx <- (lc$xc - 1) * n + lc$yc               # column-major cell index
    rs   <- rowsum(Pignite, idx)                 # summed Pignite per cell
    psum <- numeric(n * n)
    psum[as.integer(rownames(rs))] <- rs[, 1]
    psum <- matrix(psum, n, n)
    pcnt <- matrix(tabulate(idx, n * n), n, n)
    sh <- function(m, dy, dx) {                  # periodic shift
      m[((seq_len(n) - 1 + dy) %% n) + 1, ((seq_len(n) - 1 + dx) %% n) + 1]
    }
    S3 <- C3 <- matrix(0, n, n)
    for (dy in -1:1) for (dx in -1:1) {
      S3 <- S3 + sh(psum, dy, dx); C3 <- C3 + sh(pcnt, dy, dx)
    }
    # fallback if the whole 3 x 3 block is empty: Pignite at B = 0 (i0 = 0.15)
    p_empty <- ifelse(C3 > 0, S3 / pmax(C3, 1), 1 / (1 + 0.15))
    empty_ok <- empty_ok & (matrix(runif(n * n), n, n) < p_empty)
  }
  m <- empty_ok | has_ig
  storage.mode(m) <- "integer"
  m
}

# spread_fire(): converts potential ignitions (ig) into actual ignitions.
spread_fire <- function(x, y, ig, Pignite, domain, fire_cfg) {
  lc <- lattice_cells(x, y, domain, fire_cfg$cell)
  seeds <- matrix(0L, lc$n, lc$n)
  if (fire_cfg$mode == "edge") {
    seeds[, 1] <- 1L
  } else {                                   # seeded
    cand <- which(ig & Pignite >= fire_cfg$pign_min &
                  x >= fire_cfg$seed_xfrac[1] * domain &
                  x <  fire_cfg$seed_xfrac[2] * domain)
    if (length(cand) == 0) return(rep(FALSE, length(x)))
    pick <- cand[sample.int(length(cand), min(fire_cfg$n_seeds, length(cand)))]
    seeds[cbind(lc$yc[pick], lc$xc[pick])] <- 1L
  }
  ep <- if (is.null(fire_cfg$empty_pass)) "always" else fire_cfg$empty_pass
  burn <- burn_connected(passable_lattice(lc, ig, Pignite, ep), seeds,
                         as.integer(fire_cfg$connectivity),
                         fire_cfg$wrap_x, fire_cfg$wrap_y)
  ig & (burn[cbind(lc$yc, lc$xc)] == 1L)
}

# fire_ignition(): single entry point used inside the model loops.
fire_ignition <- function(x, y, ig, Pignite, domain, fire_cfg) {
  if (fire_cfg$mode == "independent") return(ig)
  spread_fire(x, y, ig, Pignite, domain, fire_cfg)
}

# Objects that must be exported to parallel workers
fire_exports <- c("default_fire", "lattice_cells", "passable_lattice",
                  "spread_fire", "fire_ignition", "FIRE")