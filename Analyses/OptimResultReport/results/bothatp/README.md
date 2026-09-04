# Both-[ATP] slide figures — 03/27, 8 mM vs 2 mM

Built 2026-09-03. Extends the single-condition slide set one level up: the same
three time-series panels plus one panel per feature, with the 2 mM condition
drawn alongside the 8 mM one — **once per rung of the ATP hypothesis ladder**,
in `H1/`, `H2/`, `H3/`, `H4/`.

## The models

The 8 mM side is the same in every subfolder: `params/spreadw_opt.m`
(= `spreadw_iter_018`), from `results/fit_report_spreadw_opt.mat`. Every ATP gate
is exactly 1 at `MgATP_ref = 8`, so **all four rungs share one 8 mM model** —
the ladder's frozen core. Only the red 2 mM curve moves between subfolders.

The 2 mM side is `params/atphyp_<rung>_opt.m` from the **`ATP-wt-lowatp`
worktree**, run at `MgATP 2 / MgADP 0 / Pi 0.8` as the ladder scored it.

| rung | free gates | what is free | **fitted non-zero** | 2 mM cost |
|---|---|---|---|---|
| H1 |  3 | chemistry: `k2 kah kamh` | `k2 0.44`, `kah 0.30` | 8.2013 |
| H2 |  7 | + SRX: `kmsrd ksr0 ksr2srd ksrd2sr` | `k2 0.44`, `kah 0.30` | 8.2013 |
| H3 | 10 | + mechanical: `eta_M kSE mu_neg` | + `eta_M −0.15` | 7.6429 |
| H4 | 21 | + `kstiff2 kSE_M dr dr2 P_bound_max ka kd k1 k_1 tau_reg A0` | + `kSE 0.07`, `kSE_M −0.15` | 6.4333 |

**H1 and H2 are the same point.** Freeing the four SRX nucleotide-state gates
bought exactly nothing — all four fitted to zero and the cost did not move. The
`H1/` and `H2/` figures are therefore identical by construction, and that
identity is the result, not a bug.

Read the middle column as *free*, not *changed*: 16 of H4's 21 gates fit to
exactly zero. What actually separates the four models is five numbers.

**The gate code lives only in `ATP-wt-lowatp`.** `params0.pAtp_*` and
`MgATP_ref` are read by that worktree's `Model/getParams.m` and
`Model/evaluateModel.m`; here they would be set and never used, and a "2 mM" run
would come back silently as an 8 mM run. Hence `RunBothAtpSlideRuns.m` runs
there and exports one `lowatp_runs.mat` per rung; nothing in this worktree
re-simulates 2 mM.

## Scripts

| script | worktree | does |
|---|---|---|
| `Analyses/AtpEffectHypotheses/RunBothAtpSlideRuns.m` | ATP-wt-lowatp | simulates 2 mM slack + staircase for every rung, writes `<rung>/lowatp_runs.mat` here |
| `Analyses/AtpEffectHypotheses/RunFVFineAllRungs.m` | ATP-wt-lowatp | fine F-V sweep (0 : 0.25 : 7 ML/s), writes `<rung>/fv_fine.mat` here |
| `RunEnvelopeBothATP.m` | this one | 3-preparation feature envelope at both `[ATP]` → `envelope_bothatp.mat` (model-free, shared by all rungs) |
| `MakeSlideFiguresBothATP.m` | this one | `<rung>/slide2_restretch_all/det/stairs.png` |
| `MakeFVFiguresBothATP.m` | this one | `<rung>/slide2_fv.png` (force+power, 2 cells) + `slide2_fv_force.png` / `slide2_fv_power.png` (1 cell each) |
| `MakeFeatureFiguresBothATP.m` | this one | `<rung>/feat_<name>.png` ×20 + `feat_contact.png` |
| `MakeAllRungsBothATP.m` | this one | drives the three figure scripts over all four rungs |

All figure scripts take `RUNG` (default `'H4'`). Re-run order: the two
simulations only after a refit; then `MakeAllRungsBothATP.m`, seconds per rung.

`RunFVFineAllRungs.m` runs the **8 mM leg once** and reuses it for all four
rungs — every gate is 1 at `MgATP_ref = 8`, so that curve is the frozen core and
running it per rung would produce four copies of one result.

## Force-velocity notes

`normalised power = (F/F_iso) · v`, so its unit is **ML/s = s⁻¹**, not watts.
The un-normalised product `F[kPa] · v[s⁻¹]` is stress × strain rate =
**kW/m³**, power per unit muscle volume; dividing by `F_iso` expresses that in
units of the preparation's own isometric stress. Multiply back by `F_iso` to
recover kW/m³.

The band on these two panels is the spread between the **two** source FV curves
(Baker2022, IsovelocityForFilip2021), each over its own isometric — *not* the
3-preparation slack envelope the feature panels use.

Two things carried over from `RunFVFineBothATP.m`:
- **0.25 ML/s is excluded from the drawn curve.** `runFVExperiment` takes the
  isometric reference at SL 2.0 but starts every shortening run at SL 2.2, so at
  the slowest speed the 0.4 s ramp spends most of its time at longer, stronger
  lengths and the point lands *above* isometric. It stays in the `.mat`.
- **Each model curve ends at its own V_max** rather than running flat along
  zero, which would imply the model still predicts something past it.

| | 8 mM | H1/H2/H3 2 mM | H4 2 mM |
|---|---|---|---|
| F_iso | 64.68 kPa | 78.68 (2/8 = 1.216) | 78.21 (1.209) |
| peak power | 0.6546 at 1.25 ML/s | 0.4317 at 1.00 | 0.4341 at 1.00 |
| V_max | > 7 ML/s | 5.50 | 5.50 |

The F-V leg is essentially rung-independent: H1→H4 moves peak power by 0.6 %
and V_max not at all. Whatever the mechanical gates buy, they do not buy it here.

## Two choices worth knowing about

**Normalisation.** Force-dimension features divide by that preparation's mean
**8 mM** `steady` — *both* conditions by the 8 mM one. `RunDataEnvelope.m`
divides each condition by its own steady, which is right for judging shape and
wrong here: it makes 2 mM force 1.0 by construction and deletes the ATP effect.
The slack panels themselves are in raw kPa for the same reason.

**SD, not SE.** The band is mean ± SD over the three 2026 preparations. SE would
be the precision of a mean that is not on the figure (the marker is 03/27, one
preparation), and — decisively — the objective's weights are λ/s² built from
this same between-preparation SD, so cost 1.0 *means* one SD. Banding with
SE = SD/√3 would draw a tolerance 42 % tighter than the optimiser enforces. Set
`BAND = 'se'` to see the alternative; the label follows.

## Caveats that limit what these figures can claim

- **No 2 mM staircase exists.** All three protocol days drive the staircase from
  a `*_Merged_8mM_Active.txt` file. That panel is `F/F_iso`, each trace over its
  own isometric, so the 2 mM curve claims only *shape*. It is a prediction
  twice over: never fitted, at an `[ATP]` never measured.
- **The 2/8 force ratio carries rundown.** 1.180 / 1.296 / 1.651 across
  03/27 / 04/03 / 04/10 — and 04/10 deliberately ran 2 mM *first*. That spread
  lands in the 2 mM band width, which is why the 2 mM bands are the wide ones on
  every force feature.

- **The 2 mM bands are not all n = 3.** 04/03 has no cycle-5 at 2 mM, so the
  last point of every 2 mM band is a TWO-preparation mean while cycles 1–4 are
  three. `MakeFeatureFiguresBothATP.m` prints a NOTE for each affected feature
  (`steady`, `vall2_dy`); the panels draw a 2-prep and a 3-prep SD identically,
  so the console is the only place this shows.

## Quoting the ATP force effect — pick a number and say which

The `steady force` panel draws **03/27 markers** against a **3-preparation
band**. Those are different quantities and they differ a lot here, so a caption
that names one number while pointing at the other will be wrong whichever way
it is read:

| | 2 mM / 8 mM |
|---|---|
| **03/27 alone** — the markers, and the preparation both models are fitted to | **×1.18** |
| 3-preparation mean — the band centre | ×1.37 |
| per preparation (cycles 1–4) | 1.18 / 1.22 / 1.64 |
| **model** | ×1.186 (H1, H3) · ×1.179 (H4) |

So the model reproduces **03/27's ×1.18 to within 0.5 %**, and does *not*
reproduce the 3-prep mean — it sits 0.7 SD below it (`steady`, `z 2 mM ≈ −0.7`).

And ×1.37 is a weak population estimate: it is carried by 04/10's ×1.64, the
day that ran 2 mM **first**, so its rundown pushes the ratio the opposite way
from the other two. The three values are 1.18 / 1.22 / **1.64** — one deviant
protocol order, not scatter. Do not quote ×1.37 as "the" ATP effect without
that caveat, and do not quote anything near ×1.30: it matches neither the
panel nor the file.
- **`SLslack` and `v_restretch` are protocol panels, not tests.** They document
  the ladder that the x label used to carry. Between-preparation SD is a few
  parts in 10⁴ of a commanded quantity, so their z is meaningless and is
  reported as `NaN` rather than as the ~6 an unguarded ratio produces. Their
  model lines are dashed so the coincident measured trace shows through — and so
  a model that failed to reach the commanded slack would still be visible.

## Naming

`ktr` is the rate constant of the exponential force redevelopment after a slack
release and nothing else, so the fit's asymptote `A` (= `ae.A + F0`, the force
it redevelops *to*) is **redeveloped force**, not a "k_tr amplitude". `Am` is
the same quantity read off the trace at the end of the recovery window instead
of from the fit — **redeveloped force (end)**. The rate itself is **slack k_tr**,
to keep it distinct from the k_tr *protocol*, a separate experiment here.

`ktr` and `rsK` share y limits: both are first-order recovery rates in s⁻¹, one
after the slack and one after the re-stretch, and the point of reading them
together is that low ATP slows them by different factors (×0.54 vs ×0.34).

## What the figures say

`z` = mean over cycles of (model − preparation mean) / SD; |z| > 1 puts the
model outside the spread of the three preparations.

Inside the envelope at both `[ATP]`: `steady`, `A`, `vall_y`, `peak2`,
`restretchSlopeStart`, `t0`, `vall2_t`.

Outside, and in the same direction at both:

| feature | z 8 mM | z 2 mM | reading |
|---|---|---|---|
| `vall_t` | +2.9 | +7.8 | valley too late |
| `peak1_dSL` | +2.8 | +6.4 | re-stretch peak too far into the stretch |
| `peak1_t` | +2.8 | +5.6 | and too late |
| `rsR2` | +3.6 | +5.3 | model recovery is a *cleaner* exponential than the data |
| `rsT0` | −4.1 | −0.4 | recovery delay too short (unscored) |

The lengthening transient is mistimed the same way at both `[ATP]`, and the
error is roughly twice as large in SD units at 2 mM — but the 2 mM bands are
wider, so that is a real widening of the gap, not a narrower tolerance. This is
the same defect the 8 mM report named; low ATP does not introduce a new one.

Where the ATP transfer itself is weakest: **`rsK` z = −0.2 at 8 mM but +4.6 at
2 mM** — the model does not slow the re-stretch recovery nearly enough at low
ATP (model 29 vs data 15 s⁻¹, against a measured ×0.34 effect). `vall2_dy`
(−1.6) and `rsA` (+1.9) go the same way: the 2 mM second-valley/recovery
amplitudes are overshot. The gate set reproduces the force elevation and most of
the `k_tr` slowing, and under-delivers on the re-stretch recovery slowing.
