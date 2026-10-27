# Neighborhood size and variability in tree growth modulate the effects of fire on woody cover spatial structure

This repository contains the R and C++ code used to run the simulations and generate the figures in the manuscript by Ricardo M. Holdo (Odum School of Ecology, University of Georgia) and William A. Hoffmann (Department of Plant and Microbial Biology, North Carolina State University).

## Associated manuscript

Holdo, R.M. & Hoffmann, W.A. (in review). Neighborhood size and variability in tree growth modulate the effects of fire on woody cover spatial structure. *Ecology Letters*.

This reference (and a DOI) will be updated on publication. If you use this code, please cite the manuscript.

The code implements an individual-based, spatially explicit tree-grass-fire model on a 300 x 300 m torus (periodic boundary) domain (`scale <- 3` in each script). Individual trees are tracked by location and diameter, with growth constrained by local tree-tree competition and grass competition within a neighborhood of radius *z*. Fire ignition and intensity depend on neighborhood basal area, grass biomass, relative humidity, and wind speed, and determine the probability of topkill for each tree. Full model equations are given in Appendix S1 of the manuscript's Supporting Information; parameter values and their sources are summarized below (see [Model parameters and sources](#model-parameters-and-sources)). Variable names in the code follow the same notation used there (e.g. `S` = size/diameter, `B` = basal area, `G` = grass biomass, `Q` = fire intensity), so the code and the manuscript equations can be read side by side.

## Repository contents

| File | Description | Figures produced |
|---|---|---|
| `Bistability.R` | Tests for bistability in steady-state tree cover as a function of initial tree size and fire frequency, comparing a spatially-explicit model (`tree_mod`) against a non-spatial, mean-field version (`tree_mod_mf`) across a range of neighborhood sizes *z*. | Fig. 1 (set `SIGMA_IND <- 0.15` at the top of the script for Fig. S2) |
| `Boundary.R` | Simulates the persistence or erosion of an imposed woody/grassy boundary over 4500 y across combinations of fire frequency, neighborhood size, growth variability, and mortality. | Figs. 2 & 3 (set `N <- 400` to `N <- 1000` for Figs. S3 & S7) |
| `Boundary_NN.R` | Tracks nearest-neighbor distances for two tree size classes (> 5 cm, > 10 cm) to quantify transient spatial clustering, without an imposed boundary. | Figs. 5, 6 & S10 |
| `Boundary_uncertainty.R` | Repeats the boundary-persistence analysis while drawing parameter values from empirical error distributions (bootstrapped from `Boot_pars_1.csv`-`Boot_pars_4.csv`), to propagate parameter uncertainty through the model. | Fig. S4 |
| `Boundary_environment.R` | Repeats the boundary-persistence analysis with an underlying environmental gradient in tree growth rate (1-2x) across the arena, instead of a homogeneous template. | Fig. S13 |
| `Simulation_examples.R` | Runs and visualizes individual example simulations (tree locations, size distributions, and spatial profiles of basal area/grass/fire intensity over time) for several representative parameter combinations. | Figs. 4, S5, S6, S11 & S12 |
| `Fire_robustness.R` | Runs the Fig. 4 boundary scenarios under the baseline seeded fire, independent ignition, 2.5 and 10 m lattices, 4-cell connectivity, a lower seeding threshold, and left-edge entry. | Fig. S8 |
| `Spatial_structure.R` | Pair correlation functions of large trees for the two Fig. 6 scenarios, and correlations between basal-area maps at increasing time lags (printed to the console; reported in the text) for the Fig. 6 and Fig. S12 scenarios. | Fig. S9 |
| `fire_functions.R` | Shared fire-ignition helpers and presets (`FIRE_BOUNDARY`, `FIRE_CLUSTER`), sourced by every simulation script. |  |
| `spatial_ops.cpp` | Rcpp/C++ helper functions: `compute_ba()` (neighborhood basal area), `toroidal_nn()` (exact nearest-neighbor distances), `toroidal_pcf_counts()` (pair counts for g(r)) and `burn_connected()` (fire-spread flood fill from seed cells). Compiled at the start of each script via `sourceCpp()`. |  |
| `Boot_pars_1.csv` … `Boot_pars_4.csv` | Bootstrapped draws of empirically estimated parameters, used only by `Boundary_uncertainty.R` (contents described [below](#bootstrapped-parameter-files)). |  |

Every simulation script sources `spatial_ops.cpp` and `fire_functions.R` and can otherwise be run independently in a fresh R session.

### Fire implementation

Fire spread uses a lattice of 5 × 5 m cells over the full domain (60 × 60 cells for 300 × 300 m). A cell carries fire if it contains at least one tree with a potential ignition (probability `Pignite`) or no trees. Every simulation uses fire spread. Each script sets its fire implementation near the top (`FIRE <- ...`; see `fire_functions.R`):

- `FIRE_BOUNDARY` (boundary scripts): in each fire year one ignition point is drawn at random among trees in the left (initially grassy) half of the domain with a potential ignition and `Pignite >= 0.5`; fire spreads from it through connected cells.
- `FIRE_CLUSTER` (`Bistability.R`, `Boundary_NN.R`, `Spatial_structure.R`): open-area seeding. One ignition point is drawn among trees with a potential ignition anywhere in the domain, so locations are chosen in proportion to `Pignite`.
- `default_fire(mode = "independent")`: every potential ignition burns, no spread; used only for comparison in `Fire_robustness.R`.
- `default_fire(mode = "edge")`: fire enters along the whole left edge; used only for comparison in `Fire_robustness.R`.

Options: `cell` (lattice cell size), `connectivity` (8 or 4), `n_seeds`, `pign_min`, `seed_xfrac` (x band for ignition points, as fractions of the domain) and `wrap_x`/`wrap_y` (the lattice wraps like the torus).

## Requirements

- R (>= 4.0)
- A working C++ compiler toolchain for `Rcpp` (e.g. Rtools on Windows, Xcode command line tools on macOS, or `build-essential` on Linux)
- R packages: `Rcpp`, `ggplot2`, `grid`, `gridExtra`, `cowplot`, `dplyr`, `tidyr`, `parallel`, `foreach`, `doParallel`

```r
install.packages(c("Rcpp", "ggplot2", "gridExtra", "cowplot", "dplyr",
                    "tidyr", "foreach", "doParallel"))
```

## Running the code

All scripts assume the working directory is this repository (they source `spatial_ops.cpp` and, for `Boundary_uncertainty.R`, read the `Boot_pars_*.csv` files using relative paths). In RStudio, open the project/folder and either source a script directly or run:

```r
setwd("path/to/this/repository")
source("Bistability.R")
```

Each script registers a parallel cluster (via `parallel`/`doParallel`) sized to the number of available cores and runs the full set of parameter combinations described in the manuscript Methods. Several scripts (particularly `Boundary.R`, `Boundary_NN.R`, `Boundary_uncertainty.R`, and `Boundary_environment.R`) simulate thousands of parameter combinations over runs of up to 5000 simulated years, and can take a long time to complete depending on the number of cores available. Each script produces its output figures as named `ggplot`/`cowplot` objects in the R environment (e.g. `fig1`, `fig2`, `figS4`) and saves each one as a PNG (e.g. `Fig1.png`, `FigS4.png`) in the working directory.

## Model parameters and sources

Code variable names, manuscript notation (Appendix S1 of the Supporting Information), default values and sources. Parameters marked † are varied in `Boundary_uncertainty.R` (see [Bootstrapped parameter files](#bootstrapped-parameter-files)). Values are hard-coded in the `# Parameters` block at the top of each model function (e.g. `tree_mod()`), so they can be changed there for new experiments.

### Core parameters (varied across scenarios)

| Code | Manuscript | Description | Values used |
|---|---|---|---|
| `N` | — | Stem density (trees ha⁻¹) | 400 (default); 1000 |
| `z` | *z* | Neighborhood radius for competition and fire effects (m) | 5, 10, 20 (∞ = mean field in `Bistability.R`) |
| `Pfire` | *P*<sub>fire</sub> | Exogenous fire probability (y⁻¹) | 0–0.8, in steps of 0.1 or 0.2 |
| `sigma.ind` | σ<sub>ind</sub> | SD of growth rate across individuals (cm y⁻¹) | 0, 0.15, 0.3 (0.2 in `Boundary_environment.R`) |
| `sigma.res` | σ<sub>res</sub> | Residual (within-individual) growth SD (cm y⁻¹) | 0, 0.15 (0.2 in `Boundary_environment.R`) |
| `mu` | μ | Adult (background) mortality rate (y⁻¹) | 0.005, 0.0075, 0.01 |

### Empirically estimated and fixed parameters

| Code | Manuscript | Description | Equation (Appendix S1) | Default | Source |
|---|---|---|---|---|---|
| `gmax` † | *g*<sub>max</sub> | Maximum grass biomass (g m⁻²) | 2, 3 | 800 | Newberry *et al.* (2020) |
| `alpha` † | α | Grass biomass reduction per unit tree basal area (ha m⁻²) | 3 | 0.18 | Newberry *et al.* (2020) |
| `lambda` † | λ | Tree diameter growth rate (cm y⁻¹) | 1 | 0.23 (SD 0.06) | Holdo (2006); Holdo, unpublished data; Rossatto *et al.* (2009) |
| `lambdar` † | λ<sub>r</sub> | Maximum resprout growth rate (cm y⁻¹) | 4 | 0.61 | NeSmith *et al.* (2021) |
| `shs` † | *s*<sub>½</sub> | Resprout growth half-saturation size (cm) | 4 | 0.8 | NeSmith *et al.* (2021) |
| `gamma` † | γ | Tree growth with vs. without grass competition | 2 | 0.69 (SD 0.18) | Riginos (2009); February *et al.* (2013); Laris *et al.* (2021) |
| `rho` | ρ | Tree crowding (tree–tree competition) parameter (ha m⁻²) | 1 | 0.05 | Calibrated to give realistic steady-state basal area |
| `a0`, `a1` † | *a*<sub>0</sub>, *a*<sub>1</sub> | Height–diameter allometry (log–log) intercept and slope | — | −0.10, 0.77 | Holdo (2007) |
| `RHref` | *RH*<sub>ref</sub> | Reference relative humidity in the open (%) | 7a | 20 | Author estimate |
| `r` † | *r* | Increase in relative humidity per unit basal area (% per m² ha⁻¹) | 7a | 0.41 | Newberry *et al.* (2020) |
| `Wref` | *W*<sub>ref</sub> | Reference wind speed in the open (m s⁻¹) | 7b | 5 | Author estimate |
| `w` † | *w* | Decay in wind speed per unit basal area (ha m⁻²) | 7b | 0.047 | Newberry *et al.* (2020) |
| `Mref` | *M*<sub>ref</sub> | Fuel moisture (%) | 6 | 40 | Higgins *et al.* (2000) |
| `i0`, `i1` † | *i*<sub>0</sub>, *i*<sub>1</sub> | Ignition probability constant and basal-area effect | 5 | 0.15, 0.23 | Newberry *et al.* (2020) |
| `m0`–`m4` | *m*<sub>0</sub>–*m*<sub>4</sub> | Fire intensity (*Q*) constant and effects of grass biomass, fuel moisture, relative humidity and wind speed | 6 | 2729, 8.684, 530, 0.907, 596 | Higgins *et al.* (2000) |
| `t0`–`t2` † | *t*<sub>0</sub>–*t*<sub>2</sub> | Topkill probability constant and effects of tree height and fire intensity | 8 | 1.48, 3.31, 0.026 | Higgins *et al.* (2000) |

Initial size at recruitment (0.5 cm) and the size at which resprouts revert to normal growth (2 cm) are set in the model loops. The fire-spread lattice settings are described under [Fire implementation](#fire-implementation).

State variables computed within the simulation loop follow the same convention: `S` (diameter, cm), `Smax` (maximum diameter reached, used for resprout growth), `B` (neighborhood basal area, m² ha⁻¹), `G` (grass biomass, g m⁻²), `Q` (fire intensity, kJ s⁻¹ m⁻¹), `RH` (relative humidity, %), `W` (wind speed, m s⁻¹) and `Pignite` (local ignition probability).

### Bootstrapped parameter files

`Boundary_uncertainty.R` propagates parameter uncertainty through the boundary simulations. Parameters estimated together from the same dataset are correlated, so they are drawn jointly: each simulation run samples one row at random from each of the four files below. Each file holds 1000 parameter sets obtained by refitting the relevant model to the original dataset resampled with replacement. The two parameters estimated as means across studies (`lambda` and `gamma`) are instead drawn from normal distributions with the means and SDs given above. The mapping from file columns to model parameters is documented in `getpars()` in `Boundary_uncertainty.R`.

| File | Columns (code names) | Model component | Data source |
|---|---|---|---|
| `Boot_pars_1.csv` | `gmax`, `alpha`, `r`, `w`, `i0`, `i1` | Grass biomass, relative humidity, wind speed and ignition probability as functions of neighborhood basal area (eqs. 3, 5, 7a, 7b) | Newberry *et al.* (2020) |
| `Boot_pars_2.csv` | `a0`, `a1` | Height–diameter allometry | Holdo (2007) |
| `Boot_pars_3.csv` | `lambdares` (= `lambdar`), `shs` | Resprout growth as a function of pre-topkill size (eq. 4) | NeSmith *et al.* (2021) |
| `Boot_pars_4.csv` | `t0`, `t1`, `t2` | Topkill probability (eq. 8) | Higgins *et al.* (2000) |

### Parameter source references

- February, E.C., Higgins, S.I., Bond, W.J. & Swemmer, L. (2013). Influence of competition and rainfall manipulation on the growth responses of savanna trees and grasses. *Ecology*, 94, 1155–1164.
- Higgins, S.I., Bond, W.J. & Trollope, W.S.W. (2000). Fire, resprouting and variability: a recipe for grass–tree coexistence in savanna. *Journal of Ecology*, 88, 213–229.
- Holdo, R.M. (2006). Tree growth in an African woodland savanna affected by disturbance. *Journal of Vegetation Science*, 17, 369–378.
- Holdo, R.M. (2007). Elephants, fire, and frost can determine community structure and composition in Kalahari woodlands. *Ecological Applications*, 17, 558–568.
- Laris, P., Yang, L., Dembele, F. & Rodrigue, C.M. (2021). Fire and water: the role of grass competition on juvenile tree growth and survival rates in a mesic savanna. *Plant Ecology*, 222, 861–875.
- NeSmith, J.E., Twine, W. & Holdo, R.M. (2021). Interspecific variation in post-disturbance growth responses of a savanna tree community and its implications for escaping the fire trap. *Biotropica*, 53, 896–905.
- Newberry, B.M., Power, C.R., Abreu, R.C.R., Durigan, G., Rossatto, D.R. & Hoffmann, W.A. (2020). Flammability thresholds or flammability gradients? Determinants of fire across savanna–forest transitions. *New Phytologist*, 228, 910–921.
- Riginos, C. (2009). Grass competition suppresses savanna tree growth across multiple demographic stages. *Ecology*, 90, 335–340.
- Rossatto, D.R., Hoffmann, W.A. & Franco, A.C. (2009). Differences in growth patterns between co-occurring forest and savanna trees affect the forest–savanna boundary. *Functional Ecology*, 23, 689–698.

## License

This code is released under the MIT License (see `LICENSE`).

## Contact

Ricardo M. Holdo — rholdo@uga.edu — Odum School of Ecology, University of Georgia
