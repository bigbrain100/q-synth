# Prior-art positioning — Q-SYNTH recovery standard

Prior-art scan (2026-07-12) to situate Q-SYNTH and avoid collision, especially
with Akhtar-Danesh and with bootstrap-based Q stability work. Bottom line: **no
existing work builds simulation-based (known-ground-truth) recovery standards for
Q design geometry (N Q-sorts × M statements).** Current guidance is either
rule-of-thumb or post-hoc stability on already-collected data.

## The gap

- **Sample size (P-set):** rule-of-thumb only. Brown (1980, 1993): retain a
  factor with ≥2 significant loadings, ~4–5 defining sorts per factor, never a
  single-person factor; composite reliability rises with the number of defining
  sorts. Common practice 12–40 participants. None derived from a recovery study.
- **Q-set size (statements):** rule-of-thumb only. Watts & Stenner (2012): ~40–80
  statements; others 30–60. No statistical basis tying M to recoverability.
- Q's own literature frames M as "breadth of the concourse," never as the
  paired-observation count that governs each by-person correlation — the reframing
  Q-SYNTH makes explicit (M = observations per correlation; N = variables).

## Nearest neighbours — and how Q-SYNTH differs

1. **Zabala & Pascual (2016), "Bootstrapping Q Methodology" (PLoS ONE
   11(2):e0148087); `qmethod` R package `qindtest`/bootstrap.** Closest work on
   *confidence* in Q. It **resamples an existing dataset** to put CIs on loadings
   and assess factor stability, handling PCA bootstrap sign/order indeterminacy.
   → **Post-hoc** stability of one collected study. Q-SYNTH is **a-priori design**:
   with a *known* planted structure it asks how large N and M must be to recover
   it — a planning tool, not a diagnostic on data you already have. Complementary,
   not overlapping.

2. **Akhtar-Danesh (2017), "A Comparison between Major Factor Extraction and
   Factor Rotation Techniques in Q-Methodology" (Open J. Applied Sciences
   7:147–156).** Compares PCA/PAF extraction and varimax/oblimin rotation on
   **three real datasets**; concludes rotation effects are substantial. → Uses
   real data with **no ground truth** and yields **no sample-size standard**. Our
   PCA-vs-centroid robustness sweep (ROBUSTNESS.md) answers the same
   which-method-matters question but with **recovery against known truth**, so
   cite as the qualitative real-data complement to our recovery result.

3. **"Impact of factor rotation on Q-methodology analysis" (PLoS ONE 2023,
   e0290728).** Real-data rotation-sensitivity study; same class as (2). Cite to
   show rotation/extraction choices are an active concern → motivates our
   method-robustness check.

4. **R-methodology sample-size simulation tradition — MacCallum, Widaman, Zhang &
   Hong (1999), "Sample size in factor analysis" (Psychological Methods).** Shows
   required N depends on communality and overdetermination, not a fixed rule. →
   This is Q-SYNTH's **theoretical anchor**: we transpose it onto Q's by-person
   matrix (roles reversed: M = observations, N = variables) and show the same
   determinants (communality c, distinctness ρ, persons/factor) govern Q recovery.
   Our finding that Brown's 4–5/factor holds only at high c is the Q analogue of
   MacCallum's communality result.

5. **Akhtar-Danesh & Wingreen (2022), `qpair` (Stata Journal); Akhtar-Danesh
   (2018), `qfactor`.** Software commands for Q analysis (paired sorts / factor
   extraction). → Tooling, not design standards; relevant because our product
   line (QPAIR/Q-TWO) is adjacent. The Q-SYNTH tool (JOSS) is a *validation /
   design-planning* utility, a different category from analysis commands.

## Residual scoop risk — low

Searched: Q sample-size simulation / Monte Carlo / factor recovery / power
analysis / synthetic Q-sorts (2016–2025). No dedicated Q recovery-simulation for
design geometry surfaced. Adjacent-but-different: a Delphi latent-structure
recovery benchmark (arXiv, ordinal data, not Q); "seven-step Q-sample design"
(PMC12875056, concourse→statement *content* design, not statistical recovery);
"Foundations of Q Methodology" (2025 Springer chapter — check for any updated
size conventions to cite). Recommend a final targeted check of these two before
submission, but none appears to pre-empt the recovery-standard contribution.

## How to frame the contribution (for Paper B / OS)

- **Not** "first ever" (per red-team ruling). Frame: *"Q's design conventions rest
  on rules of thumb (Brown; Watts & Stenner) and post-hoc bootstrap stability
  (Zabala & Pascual); we complement these with a recovery-simulation account that
  ties recoverability to communality, distinctness and persons/factor, extending
  the R-methodology tradition (MacCallum et al.) to Q's transposed matrix, and
  show the standard is robust across PCA and the traditional centroid method."*
- Cite Akhtar-Danesh (2017) + PLoS 2023 as the extraction/rotation-sensitivity
  motivation for our robustness check; cite Zabala & Pascual (2016) as the
  post-hoc counterpart to our a-priori design standard.

*Citations to verify page/year at write-up time; all located via 2026-07-12
web scan and consistent with known sources.*
