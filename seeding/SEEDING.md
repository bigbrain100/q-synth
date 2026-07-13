# Ecological calibration — do the swept c/ρ ranges match real Q data?

The sweep varies shared variance **c** ∈ {0.4, 0.6, 0.8} and factor distinctness
**ρ** ∈ {0.0, 0.4}. Are these realistic, or arbitrary? We calibrate against a
**published** Q dataset (firewall: published factor arrays only).

## Anchor: the Lipset dataset

The classic study of democratic values — Lipset (1963), analysed in Brown (1980,
*Political Subjectivity*) — distributed as the `lipset` dataset in the `qmethod`
R package (GPL, public teaching data). **9 Q-sorts × 33 statements**, forced
−4..+4, standard **3-factor** solution. Raw sorts are committed here
(`lipset-sorts.json`) so the calibration is self-contained and reproducible.

Reproduce: `node seed-calibrate.mjs`.

## How c and ρ map to the generator

- The generator builds each defining sort as `√c · ideal + √(1−c) · noise`, so a
  defining sort's loading on its own factor has expectation `√c` and **loading²
  has expectation c**. Therefore **mean defining-sort loading² in real data
  directly estimates c.**
- ρ is the correlation among the planted ideals; for recovered factors this
  equals the correlation among the factor arrays. So **|φ| between real factor
  arrays estimates ρ.**

## Result — the ranges bracket reality, with c = 0.6 as the realistic centre

| Quantity | Lipset (PCA) | Lipset (centroid) | Swept range |
|---|---|---|---|
| c = mean defining loading² | **0.624** (0.44–0.78) | 0.572 (0.35–1.02*) | 0.4 / **0.6** / 0.8 |
| mean communality (K=3) | 0.638 | 0.513 | — |
| ρ = mean \|φ\| between arrays | **0.185** (pairs 0.03–0.33) | 0.091 (0.07–0.12) | 0.0 / 0.4 |
| defining sorts / factor | 3, 3, 2 | 2, 3, 2 | — |

\* the 1.02 is a mild Heywood case on this very small N = 9 study; PCA gives a
clean 0.44–0.78 spread. Reported for transparency.

**Takeaways**

1. **c ≈ 0.6 is the realistic centre.** Real defining sorts share ~60% of their
   variance with their factor — exactly the sweep's middle level. So the headline
   result (at c = 0.6, Brown's 4–5/factor recovers only ~49–63%; ~8/factor is
   needed) speaks to **real** Q studies, not an artificial regime.
2. **ρ ≈ 0.1–0.3 sits between the swept ρ = 0 and ρ = 0.4.** Real viewpoints are
   fairly distinct but not perfectly orthogonal; the sweep brackets this.
3. The sweep's grid is therefore **ecologically anchored**, not arbitrary — a
   point to state explicitly in Paper B when defending the parameter choices.

## Caveat / next

One canonical anchor (N = 9) fixes the *centre* of the realistic c/ρ region well.
For Paper B, optionally add 1–2 open-access applied studies (larger N, more
factors) to show the range of c/ρ across topics; the method here transfers
directly (drop their raw sorts into `seed-calibrate.mjs`).
