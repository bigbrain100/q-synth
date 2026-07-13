# KADE cross-check — engine fidelity on real data

Validates q-synth's self-contained engine (the correlation / PCA / varimax
routines in `engine.mjs`, and centroid in `centroid.mjs`) against **real KADE
analysis exports** — not synthetic data.

## IRB firewall

`kade-crosscheck.mjs` reads **only numeric matrices** and reports **only
engine-vs-KADE agreement metrics**. It never surfaces statement text or any
participant viewpoint. The KADE `.xlsx` files are human-subjects data and are
**not committed** (kept local / gitignored); pass their paths as arguments.
Eigenvalues below are aggregate variance statistics (standard published Q
output), not participant-level content.

## Run

```
node kade-crosscheck.mjs <KADE_results.xlsx> [more.xlsx ...]
```

Geometry (N sorts, M statements, extracted/retained factor counts) and the
extraction method (PCA vs centroid) are auto-detected, so one script covers KADE
exports of different sizes and methods.

## Result (2026-07-12)

Two KADE exports, one per extraction method:

| Export | N × M | Extraction | Correlation matrix | Unrotated loadings | Eigenvalues | Varimax loadings |
|--------|-------|-----------|--------------------|--------------------|-------------|------------------|
| A | 7 × 22 | **PCA** | `maxAbsDiff = 0` (bit-exact ×100) | `maxAbsDiff = 0` @ 4dp | `4.8e-9` | `≤ 2e-4` (all 28 cells; rounding-level) |
| B | 24 × 40 | **centroid** | `maxAbsDiff = 0` (bit-exact ×100) | per-factor col-corr `0.9995–0.9999` | within ~2% (see note) | — |

**PCA path (export A): faithful to KADE.** The engine reproduces KADE's
correlation matrix bit-exactly, its unrotated PCA loadings and eigenvalues to
floating-point precision, and its varimax loadings to a uniform 4th-decimal
rounding level (no cell off by more than `2e-4`; varimax convergence tolerance
differs slightly between implementations, never enough to change a defining-sort
flag).

**Method detection.** Export B's correlation matrix also matches bit-exactly, but
its unrotated loadings do not match PCA — because it was extracted with
**centroid**, not PCA. PCA's first component maximizes single-factor variance, so
PCA factor-1 SS is an upper bound; KADE's `11.64 < 12.08` (our PCA) proves the
extraction was not PCA. The script uses this bound as its automatic method
detector.

**Centroid path (export B): structure faithful to KADE.** `q-synth/centroid.mjs`
(classic centroid, iterated communalities; public math, not in the closed engine)
reproduces KADE's centroid factor **structure** to a per-factor column
correlation of `0.9995–0.9999`, loadings `maxAbsDiff ≈ 0.024`. Extracted SS
("eigenvalues") land ~2% above KADE's:

```
our centroid eig : 11.8690, 4.3030, 2.1077   (converged in 11 iters)
KADE centroid eig: 11.6357, 4.2051, 2.0404
```

This ~2% is a communality-convergence convention difference (KADE's final
communalities average slightly lower). It is **scale-invariant for factor
recovery**: which sorts define which factor depends on the loading pattern, which
matches KADE at ≥0.9995 — so centroid recovery in the sweep is equivalent to
KADE's.

## Status / next

- **Both extraction paths are validated against real KADE output.** PCA (the path
  q-synth uses) is bit-exact; centroid (the traditional Q method) reproduces
  KADE's structure to ≥0.9995 column correlation.
- **Next:** run the recovery sweep under centroid as well as PCA, to show the M×N
  design standard is robust to the choice of extraction method (answers the
  traditional-Q reviewer). Centroid is a robustness check only; the tool's
  default path is PCA.
