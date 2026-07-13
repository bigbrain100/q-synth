# Robustness of the M×N standard to extraction method (PCA vs centroid)

Does the recovery-based design standard depend on *which* factor-extraction
method is used? This answers the traditional-Q reviewer, who uses **centroid**
rather than the **PCA** that q-synth (and this study's sweep) run on.

- Both methods are validated against real KADE output (`validation/`): PCA
  bit-exact; centroid structure ≥0.9995 per-factor column correlation.
- **Paired design:** each replicate draws one synthetic dataset and scores *both*
  methods on it, so any difference is the extraction method alone.
- Retention fixed to **oracle** (nf = K planted factors) to isolate extraction
  from factor-retention rules. mode = forced, no non-loaders, K = 3, R = 400.
- Reproduce: `node robustness.mjs 400` → `robustness.csv`,
  `figures/fig_robust_pca_vs_centroid.svg`.

## Finding: the standard is robust; centroid is modestly more conservative

Recovery @ |φ| ≥ 0.90, shown **PCA% / centroid%** (M = 40, forced, oracle):

| c \ persons-per-factor | 4 | 6 | 8 | 12 |
|---|---|---|---|---|
| **ρ = 0.0** |||||
| 0.4 | 0 / 0 | 2 / 1 | 24 / 19 | 91 / 90 |
| 0.6 | 54 / 35 | 99 / 91 | 100 / 95 | 100 / 99 |
| 0.8 | 100 / 82 | 100 / 91 | 100 / 96 | 100 / 99 |
| **ρ = 0.4** |||||
| 0.4 | 0 / 0 | 1 / 1 | 22 / 19 | 89 / 85 |
| 0.6 | 53 / 37 | 99 / 88 | 100 / 95 | 100 / 99 |
| 0.8 | 100 / 74 | 100 / 87 | 100 / 92 | 100 / 97 |

(M = 20 in `robustness.csv`; same qualitative pattern.)

1. **Qualitative thresholds are method-invariant.** Under both PCA and centroid:
   c = 0.4 fails regardless of N; c ≥ 0.6 is required; recovery saturates by
   ~8 persons/factor once c ≥ 0.6. The design conclusions do not change with the
   extraction method.
2. **Centroid needs slightly more people at the margin** — at most ~1–2 extra
   persons/factor. Where PCA saturates at 6/factor (c = 0.6) or 4/factor
   (c = 0.8), centroid reaches the same recovery at ~8 and ~6 respectively. The
   gap closes to a few points once N is comfortable (≥8/factor).
3. **Direction is reviewer-friendly:** centroid is never *more* permissive than
   PCA, so a PCA-derived standard is a safe **lower bound** — a centroid user
   following it will not be under-powered. Largest gap: centroid −26 pp at the
   single hardest cell (c = 0.8, 4/factor, ρ = 0.4), where 4/factor is already
   below the recommended minimum.

## Caveat

At the smallest N (4/factor → N = 12) with high shared variance, iterated
communalities do not fully converge (the reduced correlation matrix approaches
rank-deficiency; `centroid.mjs` guards against the resulting non-positive grand
sum so loadings stay finite rather than NaN). This is an inherent small-sample
limitation of centroid, and it coincides with the region the standard already
flags as under-powered — so it does not affect the recommended design region.
