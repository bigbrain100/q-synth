# How the synthetic data are generated (methods documentation)

This note documents, in full, how Q-SYNTH generates its synthetic Q-sorts and how
the results are reproduced. It is written so a reviewer or co-author can follow
the data generation without reading the source, and so the study can be
reproduced and checked end to end.

## What this is (and what it is not)

The data are produced by a **deterministic Monte Carlo simulation**: a seeded
pseudo-random-number generator drives a linear factor model. There is no
generative AI model anywhere in the pipeline. Given the same seed, parameters,
and code, the simulation reproduces byte-identical data.

## Implementation and provenance

- **Language / runtime:** JavaScript, run on Node.js (v24). The files are ES
  modules (`.mjs`). Python is used only to build the reference-verification
  worksheet, not for the simulation.
- **Analysis engine:** self-contained. `engine.mjs` holds open implementations of
  the standard routines (Pearson correlation, a Jacobi symmetric eigensolver for
  principal components, Kaiser varimax, Brown-weighted factor scores, forced-grid
  rendering); `centroid.mjs` adds centroid factoring. These are textbook
  algorithms, not tied to any closed product, and their output is verified against
  a real KADE analysis (see `validation/`): correlation bit-exact, unrotated PCA
  loadings and eigenvalues to floating point, varimax to a fourth-decimal rounding
  level. The harness no longer depends on the closed product engine.
- **Design provenance:** the plant-and-recover loop is the standard
  factor-recovery simulation design used to study sample size in factor analysis
  (MacCallum, Widaman, Zhang, & Hong, 1999), adapted to Q's by-person matrix.
  Its components each follow an established source: Tucker's congruence for
  matching recovered to planted factors (Tucker, 1951), parallel analysis for
  retention (Horn, 1965), and the forced quasi-normal distribution for rendering
  (Brown, 1980; Watts & Stenner, 2012). The seven-step statement in the
  manuscript is our articulation of this design, not a procedure attributed to a
  single author.

## The generative model, step by step

Notation: K planted viewpoints, M statements, c communality of a defining sort,
rho distinctness (correlation among viewpoints), and a forced grid whose column
counts sum to M.

**1. Seeded randomness (`makeRng(seed)`), seeded per cell.** A mulberry32-style
PRNG. `.u()` returns a uniform in [0,1); `.n()` returns a standard normal via
Box-Muller. Each design cell seeds its own generator deterministically from its
parameters, and each parallel-analysis threshold seeds from its (N, M, rendering),
via a fixed base (20260711 for the main sweep, 20260712 for robustness) hashed
with the cell coordinates. A cell's result is therefore independent of how many
workers run or in what order cells are processed, so the parallel sweep is
bit-identical to a serial one and fully reproducible.

**2. Planted viewpoints (`makeIdeals`).** Draw a K by M matrix Z of independent
standard normals. Build the lower-triangular Cholesky factor L of the K by K
correlation matrix with 1 on the diagonal and rho off it (`cholCS`). For
viewpoint k, form y[m] = sum over j <= k of L[k][j] * Z[j][m]; this makes the K
viewpoints mutually correlated at rho. Render y into the forced grid and
z-score. The result is the k-th **factor array**, the ground-truth ordering a
person perfectly saturated on factor k would give.

**3. Defining sorts (`makeParticipants`).** For each viewpoint k and each of its
defining participants, form a continuous score

    score[m] = sqrt(c) * idealZ[k][m] + sqrt(1 - c) * noise[m]

with noise a fresh standard-normal vector. The weight sqrt(c) makes the expected
squared loading of the sort on its factor equal to c, so **c is the communality
of a defining sort**. Render `score` into the forced grid to get an integer
Q-sort. Optionally add non-defining sorts that are pure noise rendered into the
grid, belonging to no viewpoint.

**4. Rendering to the grid (`render` / `generateCompositeSort`).** The continuous
score is turned into a Q-sort by ranking its M values and assigning grid
positions according to the forced distribution counts (ipsative). A free,
non-ipsative rendering (round and clamp) is available for contrast.

**5. One raw dataset.** Steps 2 to 4 produce an N by M integer matrix (N sorts,
M statements). This matrix is the "raw data" of one replicate; it is analysed
exactly as a real Q study would be (correlation, extraction, rotation, factor
arrays).

**6. Monte Carlo.** Steps 2 to 5 are repeated R times (R = 1000 in the submission
sweep) for each design cell, the RNG advancing each time, and the recovery rate
is the fraction of replicates in which every planted viewpoint is matched by a
distinct recovered factor at Tucker congruence at or above the threshold.

## Where the "raw data" lives

A Monte Carlo study of this kind does not archive a single data file. Its raw
data is the full set of generated matrices, hundreds of thousands of them across
all cells and replicates. They are not stored; they are **regenerated on demand**
from (seed + parameters + code), which reproduces them exactly. The reproducible
object is therefore the code plus the fixed seeds, not a static data dump. A
sample matrix can be exported for illustration on request.

## Reproducing every result

From `q-synth/` on Node.js v24 (self-contained; the only external library is
`xlsx`, used solely to read the local KADE file in the cross-check):

| Result | Command | Output |
|---|---|---|
| Main PCA sweep + figures | `node run.mjs 1000` (parallel, worker_threads) | `results.csv`, `figures/*.svg` |
| PCA vs centroid robustness | `node robustness.mjs 1000` | `robustness.csv`, `figures/fig_robust_pca_vs_centroid.svg` |
| Real-KADE cross-check | `node validation/kade-crosscheck.mjs <KADE_results.xlsx>` (run from repo root) | console report |
| Lipset ecological calibration | `cd seeding && node seed-calibrate.mjs` | console report |

Fixed seeds: `run.mjs` and `robustness.mjs` seed the PRNG per cell from a fixed
base (20260711 and 20260712 respectively) hashed with the cell parameters, so
each table is reproducible to the digit regardless of worker count.

## Public release

The harness is self-contained: `engine.mjs` and `centroid.mjs` provide open
implementations of every routine, so `q-synth/` can be released and re-run by
anyone under a permissive license, and the KADE cross-check validates those open
implementations. The closed product engine is not part of this release. Regenerating the
sweep with the open engine reproduces the published tables to within Monte Carlo
noise; the tiny residual differences come from the order in which the Jacobi and
the reference eigensolvers return near-equal eigenvalues, which is invariant for
recovery because factors are matched by congruence.
