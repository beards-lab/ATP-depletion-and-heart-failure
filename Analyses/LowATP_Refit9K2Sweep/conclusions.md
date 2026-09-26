# Low ATP on the current incumbent: a rough k2 (A2-detachment) sweep

2026-09-26. Core `robust1_seed` = `refit9_opt` parameters + the robust feature reader
(8 mM cost 0.424047, reproduced bit-for-bit before and after the floored-campaign merge;
`refit9_opt` 0.207910 likewise — `results/baseline8_{premerge,postmerge}.log`).

## Question
With **only** A2 detachment feeling ATP (hydrolysis not gated — the physiological choice),
how far does the 8 mM incumbent get towards the pooled 2 mM data?

## Design
* [MgATP] = 8 mM (reference) and 2 mM. At 2 mM, R2 is multiplied by g = (2/8)^p via the merged
  power gate `pAtp_k2` (exactly 1 at 8 mM). g ∈ {1, 0.75, 0.54, 0.44, 0.35, 0.30}.
* Scored **absolutely** against the pooled 2 mM target (`pooledTarget2mM`, regenerated with the
  current feature set; `pooledKtrTarget(2)`, new), incumbent objective moved onto the 2 mM spread
  by `rescaleFnToCondition`. Dropped at 2 mM: `FV_normpowerAvg` (runFVExperiment always simulates
  8 mM) and the `steady[1]` 52.5–87.7 kPa guard (an 8 mM window; 2 mM steady ≈ 87–91 kPa).
* Added: FV2, the in-preparation 2 ML/s force–velocity point, pooled over the two preparations that
  have it (04/03, 04/10): P/P0 0.284 ± 0.062 (8 mM), 0.203 ± 0.041 (2 mM), ratio 0.715 ± 0.012.
* **score_rel**: each χ² term × min(z,3)/3, z = the data's ATP-effect z-score |r−1|/√(s²+sr²)
  (`t0_crossing` takes t0's z); guards/bounds weight 1. **score_plain**: all terms.

## Result (`results/k2sweep.txt`)

| g | 1.00 | 0.75 | **0.54** | 0.44 | 0.35 | 0.30 |
|---|---|---|---|---|---|---|
| score_rel | 22.71 | 7.58 | **4.52** | 6.46 | 11.02 | 15.17 |
| score_plain | 27.13 | 9.80 | **7.86** | 12.91 | 22.76 | 31.61 |

2/8 ratios, model at g = 0.54 (0.44) vs pooled data:

| feature | model | data |
|---|---|---|
| steady force | 1.24 (1.32) | 1.32 ± 0.09 |
| re-stretch peak | 1.39 (1.54) | 1.40 ± 0.12 |
| slack ktr (ktrfit) | 0.75 (0.68) | 0.57 ± 0.10 |
| ktr protocol k80 | 0.68 (n/a) | 0.74 ± 0.14 |
| rsK80 | 0.49 (0.39) | 0.38 ± 0.07 |
| t0 | 1.54 (1.77) | 1.43 ± 0.17 |
| rsTfitE | 1.04 (0.98) | 1.68 ± 0.23 |
| valley time | 1.19 (1.26) | 0.82 ± 0.11 |
| FV2 P/P0 | 0.60 (0.47) | 0.715 ± 0.012 |

## Conclusions
1. **One physiological lever gets most of the way.** Slowing A2 detachment ~2× (g ≈ 0.5) takes the
   relevance-weighted 2 mM χ² from 22.7 to 4.5 and moves force, re-stretch peak, rsK80 and t0 in the
   right direction with roughly the right size. Both scores pick the same g = 0.54 — close to the
   poster's k2 exponent (0.44 ⇒ g = 0.54), found here on a different core and target.
2. **But k2 alone cannot match force and kinetics together.** Force ×1.32 needs g ≈ 0.44, where
   slack ktr is still ×0.68 (data 0.57) while t0 (×1.77) and FV2 (×0.47) are already overshooting.
   This is the rigid F0–ktr coupling through k2, now on the incumbent. (The floored-campaign lean5
   core matched force and ktr together at g ≈ 0.35; refit9 is less ktr-sensitive per unit of force.)
3. **Two features move the wrong way**: recovery time rsTfitE barely changes (data ×1.68) and the
   valley time lengthens (data shortens ×0.82). Both are re-stretch timing features — the same block
   that tops the 8 mM misfit (OverstretchedHeads). Consistent with the floored-campaign verdict:
   the remaining low-ATP error is largely the 8 mM re-stretch defect, not missing ATP chemistry.
4. **F–V is too ATP-sensitive** with a k2-only gate (FV2 ×0.60 vs 0.715 at best g) — the model's
   F–V shoulder collapses when detachment slows; the in-preparation data say the depression is ~29 %.
5. Largest remaining 2 mM terms at g = 0.54: ktrfit_shape (ladder slope), vall_y, vall_t,
   vall2_dy, peak2, rsTfitE — re-stretch shape/timing and the ktr ladder.

## Caveats
* Rough 1-D sweep, no refit; 6 points. The 2 mM target is n = 3 (FV2 n = 2).
* ktr-protocol k80 is undefined for g ≤ 0.44: the 2 mM protocol holds only 58.5 ms after the
  release, and slow redevelopment never reaches 80 % in the window (k63 stays defined:
  0.76 at g = 0.54, 0.69 at 0.44 vs data 0.73). Those NaN terms may be scored differently by
  evalFeatureCost — treat score_plain at g ≤ 0.44 with care.
* The trace figure shows single-prep 03/27 data in its own kPa; the objective scores the pooled
  features (Fref 69.35 kPa), hence the visible level offset.

## Files
`RunRefit9K2Sweep.m` (sweep), `SummarizeK2Sweep.m` (table/JSON/figures), `fv2Target.m`,
`fv2Read.m`, `fv2Model.m`, `RegenPooledTargets.m`, `RunBaseline8mM.m`, `PlotRefit9Transitions.m`;
results in `results/`.
