# ROADMAP — consolidated next steps

The project's TODOs, pulled out of the lab diary, the restretch reports, and the meeting
abstract into one prioritized backlog grouped by theme. Source docs keep their history; this is
the live list.

_Last consolidated: 2026-06-16._

## Attachment kinetics (the recurring lever)

- [ ] **Reach a ~30% attached pool (p1+p2) at steady state** — adjust rate constants / `xrate`
  toward the physiological attached fraction. Pursue via the kinetics the sensitivity analysis
  flags (`k2`, `kah`, transition geometry), **not** via the mis-grounded per-bin occupancy flag.
  (source: `notes/labdiary.md`; see `Analyses/AttachedPool_vs_ATPase`, `Analyses/BindingSiteOccupancy`)
- [ ] **Test a velocity-dependent vernier effect** —
  `f_vernier(|v_hs|) = 1 + alpha * |v|/(|v| + v_ref)` on the attachment rate, as a
  complementary test of the vernier hypothesis. (source: `notes/labdiary.md`)
- [ ] **Fix the throttled-hydrolysis gap** identified in AttachedPool_vs_ATPase so turnover and
  attached pool become mutually consistent at baseline. (source: `Analyses/AttachedPool_vs_ATPase`)

## Restretch shape (last-slack double-peak)

- [ ] **Raise `kSE` (try 3–5×)**, let `UseA2AttachmentShift` fire at threshold strain, then tune
  `slope` and `s_threshold_R` for valley/second-peak sharpness. Use `TuneRestretch.m` for the
  sweep. (source: `Analyses/RestretchMechanisms` discrepancy appendix)
- [ ] **Fill in the sweep table** from `AnalyseRestretchMechanisms.m` output. (source: RestretchMechanisms)
- [ ] **Phase-2 fit:** fminsearch on the best variant, 10 free params, ~500 evals; then **verify the
  FV curve is unaffected**, and **write the resulting g values to a named params file**.
  (source: RestretchMechanisms)
- [ ] **Fit the time-course of a single slack** starting from `ModelParamsSlackKtrOpt`, matching the
  full transient shape. (source: `notes/labdiary.md`)
- [ ] **Fix `vall2_dy` (second-peak rise rate):** best config so far is `Catch k=0.05 + kSE=3000`
  (feature cost −3%, E_slack −55%); remaining issue is the second-peak rise (`vall2_dy` 3.4 vs
  baseline 1.0). Investigate DRX pool dynamics and `ka` for reattachment speed after the valley.
  (source: `notes/labdiary.md`)

## ATP-dependence (the headline science)

- [ ] **Mechanistic identification of the ATP-dependent contractile differences** — explain the
  rundown-corrected force *increase* at low ATP. Foundation for coupling to mitochondrial ATP
  production models for heart-failure interventions. (source: meeting abstract;
  `Analyses/LowATP_ForceEnhancement`)

## Identifiability / methodology

- [ ] **Understand why steady-state force is so sensitive to the right side of the rate-transition
  curve** — investigate piecewise strain-dependent detachment (`PieceWiseStrainDepParams/X`) and
  why right-side breakpoints disproportionately affect isometric force. (source: `notes/labdiary.md`)
- [ ] **Respect the 11-DOF ceiling** — concentrate fitting on `ekSE`, `estiff`, `dr2`/`dr`, `k2`,
  `kah`; adding parameters does not add constraint. (source: `Analyses/SensitivityAnalysis`)

## Housekeeping (from the reorganization)

- [ ] **Fill the `LastSlackIdentification` conclusions stub.**
- [ ] **Decide the fate of `_archive/`** (purge, or move to git-LFS / out-of-repo store; ~400 MB).
- [ ] **Optional in-file refactors** tracked separately in `Docs/refactoring-plan.md`.
- [ ] **Prune `PassiveTitin/`** internal clutter (out of scope for the layout reorg).
