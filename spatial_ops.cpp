// spatial_ops.cpp
// Rcpp implementations of the spatial hotspots in the tree-grass IBM:
//   1. compute_ba()          -- neighborhood basal area for all K trees
//   2. toroidal_nn()         -- nearest-neighbor distances for a size class
//                               (exact search; see below)
//   3. toroidal_pcf_counts() -- pair counts for the pair correlation
//                               function
//   4. burn_connected()      -- flood fill for the fire-spread lattice
//
// Compile from R with: Rcpp::sourceCpp("spatial_ops.cpp")
// Requires: Rcpp package and a C++ compiler (Rtools on Windows, Xcode on Mac)

#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------
// compute_ba
//
// For each of K trees, sums the basal-area contributions of all neighbors
// within radius z, using a grid to avoid the full K x K distance matrix.
// Toroidal (wraparound) distances are used throughout.
//
// Arguments:
//   x, y       -- tree coordinates (length K)
//   s2_scaled  -- s^2 / (4*z^2) for each tree (per-tree ba contribution)
//   z          -- neighborhood radius (m)
//   domain     -- side length of the square domain (100 * scale)
//   cell_size  -- grid cell size (should equal z so 3x3 block suffices)
//
// Returns: numeric vector of length K with neighborhood ba for each tree.
// -----------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector compute_ba(NumericVector x, NumericVector y,
                         NumericVector s2_scaled,
                         double z, double domain, double cell_size) {
  int K = x.size();
  double half_domain = domain / 2.0;
  double z2 = z * z;

  int nx = (int)std::ceil(domain / cell_size);
  int ny = (int)std::ceil(domain / cell_size);

  // Assign each tree to a grid cell (0-indexed internally)
  IntegerVector cx(K), cy(K);
  for (int i = 0; i < K; i++) {
    cx[i] = std::min((int)std::ceil(x[i] / cell_size), nx) - 1;
    cy[i] = std::min((int)std::ceil(y[i] / cell_size), ny) - 1;
  }

  // Build a list of tree indices per cell for fast lookup
  // cell index = cx * ny + cy  (0-indexed)
  int n_cells = nx * ny;
  std::vector<std::vector<int>> cell_members(n_cells);
  for (int i = 0; i < K; i++) {
    cell_members[cx[i] * ny + cy[i]].push_back(i);
  }

  NumericVector ba(K, 0.0);

  for (int i = 0; i < K; i++) {
    // Iterate over 3x3 block of neighboring cells (toroidal wraparound)
    for (int dxi = -1; dxi <= 1; dxi++) {
      int ncx = ((cx[i] + dxi) % nx + nx) % nx;  // wrap
      for (int dyi = -1; dyi <= 1; dyi++) {
        int ncy = ((cy[i] + dyi) % ny + ny) % ny;
        int cell_idx = ncx * ny + ncy;

        // Check every tree in this candidate cell
        for (int j : cell_members[cell_idx]) {
          double dx = std::abs(x[i] - x[j]);
          if (dx > half_domain) dx = domain - dx;
          double dy = std::abs(y[i] - y[j]);
          if (dy > half_domain) dy = domain - dy;
          // Include if within radius z (the focal tree itself is included)
          if (dx * dx + dy * dy <= z2) {
            ba[i] += s2_scaled[j];
          }
        }
      }
    }
  }

  return ba;
}

// -----------------------------------------------------------------------
// toroidal_nn
//
// For each point in a set of n points, finds the distance to its nearest
// neighbor (excluding self), using a grid and toroidal distances.
//
// Exact search: the search expands ring by ring around each point's cell and
// stops only when no unvisited cell can contain a closer point (any point in
// ring r+1 is at least r * cell_size away), so the result is the true
// toroidal nearest neighbor, however sparse the points.
//
// Arguments:
//   px, py     -- coordinates of the point set (length n)
//   domain     -- side length of the square domain (100 * scale)
//   cell_size  -- grid cell size for indexing (affects speed only)
//
// Returns: numeric vector of length n with nearest-neighbor distances
//          (R_PosInf only if n == 1).
// -----------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector toroidal_nn(NumericVector px, NumericVector py,
                          double domain, double cell_size) {
  int n = px.size();
  double half_domain = domain / 2.0;

  int nx = (int)std::ceil(domain / cell_size);
  int ny = (int)std::ceil(domain / cell_size);

  IntegerVector cx(n), cy(n);
  for (int i = 0; i < n; i++) {
    cx[i] = std::min((int)std::ceil(px[i] / cell_size), nx) - 1;
    cy[i] = std::min((int)std::ceil(py[i] / cell_size), ny) - 1;
    if (cx[i] < 0) cx[i] = 0;
    if (cy[i] < 0) cy[i] = 0;
  }

  int n_cells = nx * ny;
  std::vector<std::vector<int>> cell_members(n_cells);
  for (int i = 0; i < n; i++) {
    cell_members[cx[i] * ny + cy[i]].push_back(i);
  }

  // Beyond this ring every cell has been visited at least once
  int max_ring = std::max(nx, ny) / 2 + 1;

  NumericVector dmin(n, R_PosInf);

  for (int i = 0; i < n; i++) {
    for (int ring = 0; ring <= max_ring; ring++) {
      for (int dxi = -ring; dxi <= ring; dxi++) {
        int ncx = ((cx[i] + dxi) % nx + nx) % nx;
        for (int dyi = -ring; dyi <= ring; dyi++) {
          // Only cells on the perimeter of the current ring
          if (std::abs(dxi) != ring && std::abs(dyi) != ring) continue;
          int ncy = ((cy[i] + dyi) % ny + ny) % ny;
          int cell_idx = ncx * ny + ncy;
          for (int j : cell_members[cell_idx]) {
            if (j == i) continue;
            double dx = std::abs(px[i] - px[j]);
            if (dx > half_domain) dx = domain - dx;
            double dy = std::abs(py[i] - py[j]);
            if (dy > half_domain) dy = domain - dy;
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist < dmin[i]) dmin[i] = dist;
          }
        }
      }
      // Points in ring + 1 are at least ring * cell_size away
      if (dmin[i] <= ring * cell_size) break;
    }
  }

  return dmin;
}

// -----------------------------------------------------------------------
// toroidal_pcf_counts
//
// Counts unordered point pairs in distance bins [0, dr), [dr, 2 dr), ...
// up to rmax, using toroidal distances. On a torus no edge correction is
// needed, so the pair correlation function is obtained in R as
//   g(r) = 2 * count(r) * A / (n * (n - 1) * annulus_area(r)).
// Valid for rmax < domain / 2.
// -----------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector toroidal_pcf_counts(NumericVector px, NumericVector py,
                                  double domain, double rmax, double dr) {
  int n = px.size();
  int nb = (int)std::ceil(rmax / dr);
  double half_domain = domain / 2.0;
  NumericVector counts(nb, 0.0);
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      double dx = std::abs(px[i] - px[j]);
      if (dx > half_domain) dx = domain - dx;
      double dy = std::abs(py[i] - py[j]);
      if (dy > half_domain) dy = domain - dy;
      double d = std::sqrt(dx * dx + dy * dy);
      if (d < rmax) {
        int b = (int)(d / dr);
        if (b < nb) counts[b] += 1.0;
      }
    }
  }
  return counts;
}

// -----------------------------------------------------------------------
// burn_connected
//
// Flood fill used by the fire-spread algorithm (see fire_functions.R).
//   mat    0/1 lattice (rows = y cells, columns = x cells); 1 = can carry fire
//   seeds  0/1 lattice of the same size; fire starts in every cell where
//          seeds == 1 and mat == 1
//   connectivity  4 (von Neumann) or 8 (Moore) neighbors
//   wrap_x, wrap_y  TRUE = the lattice wraps in that direction (torus)
// Returns an integer matrix (1 = reached by fire, i.e. burns).
// -----------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix burn_connected(IntegerMatrix mat, IntegerMatrix seeds,
                             int connectivity, bool wrap_x, bool wrap_y) {
  int nr = mat.nrow();
  int nc = mat.ncol();
  IntegerMatrix out(nr, nc);
  std::vector<int> qr, qc;
  qr.reserve(nr * nc);
  qc.reserve(nr * nc);
  for (int r = 0; r < nr; r++) {
    for (int c = 0; c < nc; c++) {
      if (seeds(r, c) == 1 && mat(r, c) == 1) {
        out(r, c) = 1;
        qr.push_back(r);
        qc.push_back(c);
      }
    }
  }
  size_t head = 0;
  while (head < qr.size()) {
    int r = qr[head];
    int c = qc[head];
    head++;
    for (int dr = -1; dr <= 1; dr++) {
      for (int dc = -1; dc <= 1; dc++) {
        if (dr == 0 && dc == 0) continue;
        if (connectivity == 4 && dr != 0 && dc != 0) continue;
        int rr = r + dr;
        int cc = c + dc;
        if (cc < 0 || cc >= nc) {
          if (!wrap_x) continue;
          cc = (cc % nc + nc) % nc;
        }
        if (rr < 0 || rr >= nr) {
          if (!wrap_y) continue;
          rr = (rr % nr + nr) % nr;
        }
        if (mat(rr, cc) == 1 && out(rr, cc) == 0) {
          out(rr, cc) = 1;
          qr.push_back(rr);
          qc.push_back(cc);
        }
      }
    }
  }
  return out;
}
