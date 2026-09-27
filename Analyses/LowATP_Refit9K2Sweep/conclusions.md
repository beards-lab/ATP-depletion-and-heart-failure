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
* ktr-protocol k80 was undefined for g ≤ 0.44 (the 2 mM hold is 58.5 ms; slow redevelopment
  never reaches 80 %). `extractKtrAttributes` now fits A(1−e^{−k(t−t0)}) from the 20 % crossing and
  reads k80 off the fit in that case (flag `ktrP_k80_extrap`); the sweep was re-scored
  (`PatchK2SweepK80.m`): ratios 0.56 / 0.43 / 0.35 at g = 0.44 / 0.35 / 0.30, best g unchanged.
  A floating A on a truncated span biases k low (~86 %), so extrapolated values are a lower bound.
* The trace figure (`PlotPooledTrace.m`) shows the POOLED band: each prep rundown-lifted, divided by
  its own 8 mM isometric and multiplied by Fref = 69.35 kPa — the scale of the scored features. Late
  plateau: 8 mM data 69 vs model 71.5 kPa (+3 %), 2 mM 90 vs 88 kPa (−2.5 %). An earlier version drew
  prep 03/27 in raw kPa (its isometric is 15 % above the pool), which made the model look ~15 % low.

## Second lever on top of k2 (`RunAtpLeverSensitivity.m`, `results/levers.txt`)
One lever at a time, ×0.9 / ×1.1 at 2 mM only, on top of R2 ×0.54; same score. Candidates ordered
by chemical plausibility (tiers from `Analyses/AtpEffectHypotheses/RunPlausibleSets.m`):

| lever (tier) | Δ score_rel ×0.9 | ×1.1 | better side moves toward data |
|---|---|---|---|
| ksrd2sr (A) | +0.009 | −0.005 | nothing |
| **ksr0 (A)** | **−0.084** | +0.090 | ktrP_k80, t0, FV2 (ktrfit worse) |
| kmsrd (AB) | −0.026 | +0.034 | ktrfit (ktrP, t0, FV2 worse) |
| A0 (B) | +0.303 | −0.091 | force, ktrfit, rsK80, t0; five worse |
| tau_reg (B) | +0.002 | −0.003 | nothing |
| k_1 (C, Pi) | +0.107 | −0.058 | ktrP, rsTfitE; force/rates worse |
| k1 (C, Pi) | −0.098 | +0.143 | ktrP, rsTfitE; force/rates worse |
| kah (X) | +0.182 | **−0.148** | force, peak, rsK80, t0, FV2 |
| ka (X) | +0.422 | −0.195 | force, ktrfit, rsK80, t0; five worse |
| kd (X) | **−0.204** | +0.265 | force, ktrfit, rsK80, t0, rsTfitE |

* The two tier-A levers are weak; the largest effects come from the controls (kd, ka, kah).
* kah helps when **raised** at low ATP — the opposite sign to the H1 gate (×0.66). Together with
  the literature check (`kah_atp_dependence.md`: hydrolysis acts on bound ATP; binding is saturated
  far below 2 mM) this says the H1 kah gate was compensating, not describing chemistry.

### Follow-up: k2 × SRX-entry grid (`RunKsr0Grid.m`, `results/ksr0grid.txt`)
ksr0 is the most plausible lever that moved the right way. Grid g ∈ {0.54, 0.44, 0.35} × ksr0
factor ∈ {1, 0.7, 0.4, 0.2}: best **g = 0.54, ksr0 ×0.4, score_rel 4.13** (k2 only 4.52, −9 %).
It corrects the F–V over-sensitivity (FV2 ×0.60 → 0.74, data 0.72), the ktr protocol (0.68 → 0.72,
data 0.74) and t0 (1.54 → 1.35, data 1.43), but slack ktr gets worse (0.75 → 0.80, data 0.57) and
force, rsK80, rsTfitE and valley time do not move. It does not allow a stronger k2 gate (g = 0.44
stays worse at every ksr0). **A modest, chemically sensible partner — not the forward-limb lever
that would decouple force from slack ktr.** The residual is still dominated by re-stretch timing.

### Extension on top of B (`RunAtpLeverSensitivity.m` with BASE = ksr0 ×0.4; `extB*.log`)
Tested the two levers that moved the features B leaves wrong: kmsrd (improved slack ktr) and k1 (Pi,
the only lever that moved rsTfitE the right way).

| on top of B | score_rel | plain |
|---|---|---|
| B (R2 ×0.54, ksr0 ×0.4) | 4.13 | 6.73 |
| + kmsrd ×0.8 / ×0.6 / **×0.4** / ×0.25 | 3.88 / 3.62 / **3.33** / 3.16 | 6.50 / 6.29 / **6.13** / 6.16 |
| + k1 ×0.8 / ×0.6 | 4.10 / 5.22 | 6.65 / 8.40 |

**Set C = R2 ×0.54 + ksr0 ×0.4 + kmsrd ×0.4** (best plain score; diminishing returns beyond). Slower
return from the SRX-ADP pool slows slack ktr without touching force: ktrfit ×0.58 (data 0.57) at
steady ×1.24. But it also slows the ktr protocol (×0.59 vs data 0.74) and shortens t0 (×1.19 vs 1.43):
the data's slack-vs-protocol ktr gap (0.57 vs 0.74) cannot be closed by one pool lever — consistent
with the reviewers' compliance explanation of slack ktr. **Plausibility is weak** (tier AB; the
AtpGateHillForm analysis rejected an ATP-binding reading of kmsrd on its implied Km), so C is a
hypothesis to test (2 vs 8 mM SRX/mant-ATP data), not a result. rsTfitE and valley time remain wrong
in A, B and C — the re-stretch defect.

Full-slide figures of the three sets (slack + ktr protocol + F–V, pooled bands):
`results/bestfit_{k2only,k2ksr0,k2ksr0kmsrd}.png` (`RunBestFitFigures.m`).

## Feature table, 8 mM vs 2 mM sets A / B / C (`BuildFeatureTable.m`, `results/feature_table.txt`)
Every scored or read-out feature, from the saved runs (no new simulation):
* Force level, re-stretch peak, redeveloped force: within 1 SD in all three sets.
* **Valley after the re-stretch peak and second peak: the model collapses twice as deep at 2 mM**
  (valley 55 vs 82 kPa; vall2_dy −31 vs −14.5) in A, B and C — the re-stretch defect amplified at low ATP.
* Recovery time (rsTfitE ×1.0 vs 1.68) and valley time (×1.19 vs 0.82) wrong in every set.
* Slack ktr: only C matches (25.3 vs 25.6 s⁻¹); ktr protocol k63 too fast in A/B.
* F–V: peak power 0.44–0.48 vs Baker 2 mM 0.52 (not scored); FV2 in range only in B/C.
* Read-outs: XTOR 13–15 s⁻¹ (bound 3–15, i.e. at the ceiling at 8 mM); attached fraction 0.29 (8 mM)
  → 0.37 (2 mM), inside the objective's 0.10–0.45 bound but above the literature 0.12–0.15.

## Steady force with its spread (`PlotSteadyStateSD.m`, `results/steady_state_sd.png`)
Trace plots use prep 03/27, whose isometric (79.8 kPa) is 15 % above the pooled level (69.3 kPa),
so the model plateau looks low there. Against the pooled per-segment values the model is within
+0.9…+1.2 raw SD at 8 mM and −0.2…−0.35 raw SD at 2 mM. The raw between-prep SD of steady (relative
to each prep's own isometric) is 4 % (8 mM) / 8 % (2 mM); the moderated SD used by the objective is
17 % / 20 %, i.e. steady's weight is ~17× below its own precision — the argument for raising it
(caveat: steady and the anchor are the same-length isometric, so part of the 4 % is
self-normalisation).

## Files
`RunRefit9K2Sweep.m` (sweep), `SummarizeK2Sweep.m` (table/JSON/figures), `PatchK2SweepK80.m`,
`RunAtpLeverSensitivity.m` + `SummarizeLevers.m`, `RunKsr0Grid.m` + `SummarizeKsr0Grid.m`,
`PlotSteadyStateSD.m`, `kah_atp_dependence.md` (literature), `fv2Target.m`,
`fv2Read.m`, `fv2Model.m`, `RegenPooledTargets.m`, `RunBaseline8mM.m`, `PlotRefit9Transitions.m`;
results in `results/`.
