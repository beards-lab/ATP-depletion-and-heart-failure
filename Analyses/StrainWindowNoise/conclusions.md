# Where the objective's numerical noise comes from, and the fix

**Status:** resolved, 2026-09-24. The full 8 mM objective at `refit8_opt` goes from **19.5 % noise
to 0.3 %** under 1e-6 parameter perturbations, and its answer converges as precision increases.
The fix is four opt-in parameters; with them off the model reproduces the old runs bit-for-bit
(55/55 slack features identical). All numbers below are measured, on `params/refit8_opt.m`.

## Verdict

| Question | Answer |
|---|---|
| Do the strain-window shifts delete force-bearing bridges? | **No.** Over the 11 shifts of a slack run, the deleted bins hold at most **1e-5 of the attached mass (0.01 %)** and **0.03 kPa**. Heads occupy −34…+66 nm; the deleted bins sit at ±74–112 nm. The earlier diagnosis ("the truncation is load-bearing") was wrong. |
| Why did widening the grid "break" the model? | Two separate bugs, neither physics. (1) The piecewise strain-dependent rates are evaluated past their last knot by **extrapolating the end cubic**, which turns **negative** (R2 = −460 s⁻¹ at +150 nm), i.e. creates mass ("Leaking probs"). (2) `MaxStrainArraySize = 128` did not widen the grid at all: `getParams` only centres the grid when the cap is *strictly less* than the `Slim_l:Slim_r` length (128 here), so it built a lopsided −250…+4 nm grid. With the extrapolation clamped and the grid centred, 60, 90 and 200 bins give **identical** features. |
| So where did the noise come from? | **Reading features off the solver's own time steps.** `evaluateModel` returned the trace only at `ode15s` steps (at `RelTol 1e-2`), and `rsK80`, `rsT20`, `ktr80`, `peak1_dSL` are crossing-time features read from that trace. A 1e-6 parameter change moves the step sequence, hence the samples, hence the features. A dense output grid alone takes the slack objective from **8.4 % to 1.0 %**. |
| Why not just tighten the tolerance? | It stalled. At `RelTol 1e-3` a slack run never left the initial **isometric hold** (1.02 M RHS calls, stuck at t = 16.9 s). Cause: the tension-only titin Maxwell element charges as `kSE_M·max(0, dSLc)` — a **corner exactly at isometric rest** (`dSLc = 0`) inside a force feedback loop, so the solver chatters there. A C∞ rectifier (SiLU, exactly 0 at 0) removes it: 30 s isometric at `RelTol 1e-3` goes from a stall to **5 000 calls**, same force (71.028 kPa). |
| With the fix, is the answer converged? | **Yes.** Slack objective noise 0.2 % at `RelTol 1e-3`, 0.03 % at `1e-4`; unperturbed cost 0.3557 vs 0.3567 (0.3 %). Full objective 0.3 % at `1e-3`. |
| Is the 2 nm strain bin fine enough? | **Nearly.** Converged `rsK80` segment 4 ≈ 58.0; the current grid (2 nm, even bin count) gives 56.9 (−1.8 %), every other scored rate within 0.5 %, isometric force within 0.01 %. An **odd** bin count at 2 nm is much worse (60.5, slack cost +26 %): keep the count even. |

## The fix (all default off)

| Parameter | Set to | What it does |
|---|---|---|
| `SmoothRectifierEps` | `0.01` µm/s | SiLU instead of `max(0,·)` in the titin Maxwell element ([dPUdT_CombinedTransitions.m](../../Model/dPUdT_CombinedTransitions.m)). Makes tight tolerances affordable. |
| `odeRelTol`, `odeAbsTol` | `1e-3`, `1e-6` | Integrator tolerance, formerly hardcoded at `1e-2`/`1e-4` in [evaluateModel.m](../../Model/evaluateModel.m). `1e-4`/`1e-7` costs ~1.3× more for 0.03 % noise. |
| `odeOutputDt` | `5e-4` s | Fixed output grid, so features read a trace that does not depend on where the steps landed. **Non-negotiable**: at `RelTol 1e-3` *without* it, per-segment ktr read 40.0 and 36.2 s⁻¹ against a true flat ~45.3. |
| `ClampStrainDepExtrapolation` | `true` | Holds each piecewise rate at its edge-knot value beyond the knots. Inert on the 60-bin grid (cost 0.3565 vs 0.3564); prevents mass creation on any wider grid. |

Also kept, as a diagnostic: `params.LogWindowShifts` → [logWindowShift.m](../../Auxiliary/diagnostics/logWindowShift.m)
records what each shift deletes (serial runs, global `WINDOW_SHIFT_LOG`).

## Measurements

**Slack objective, 9 points 1e-6 apart** ([RunSlackNoise.m](RunSlackNoise.m)):

| Setting | Slack cost SD | `ktr80` relSD | `rsK80` | `rsT20` | s/run (5 parallel) |
|---|---|---|---|---|---|
| baseline (`RelTol 1e-2`, solver-step output) | **8.4 %** | 0.76 % | 0.80 % | 2.78 % | 94 |
| dense output | 1.0 % | 0.00 % | 0.26 % | 0.98 % | 135 |
| SiLU + dense | 3.5 % | 0.01 % | 0.46 % | 1.92 % | 156 |
| SiLU + `RelTol 1e-3` + dense | **0.2 %** | 0.00 % | 0.03 % | 0.12 % | 143 |
| SiLU + `RelTol 1e-4` + dense | **0.03 %** | 0.00 % | 0.01 % | 0.01 % | 181 |
| `RelTol 1e-3` without SiLU | did not finish in 900 s | | | | |

**Full objective at `refit8_opt`** ([RunObjectiveNoise.m](RunObjectiveNoise.m)), SiLU + `RelTol 1e-3` +
dense + clamp: base 0.3557, mean 0.3564, **SD 0.0009 (0.3 %)**, range 0.3551–0.3580. Largest per-term
scatter `rsK80` 0.8 %, `rsT20` 0.5 %; every other term ≤ 0.1 %. (The ktr thread measured 19.5 %
under the old settings.)

**Where the stall was** ([ProbeIsometricStall.m](ProbeIsometricStall.m), isometric 0–30 s, `RelTol 1e-3`):
as fitted — stalled at 8.6 s after 166 k calls; `mu_neg = mu` — 5.6 k calls; no `A_reg` — 167 k; smooth
dashpot blend (w = 0.1–2 kPa) — still stalled; **SiLU rectifier, ε = 0.01 or 0.001 — 5.1–6.9 k calls**,
also at `RelTol 1e-4` (7.5 k).

**Strain grid** ([ProbeGridConvergence.m](ProbeGridConvergence.m); smooth rectifier, `RelTol 1e-4`, dense, clamp):

| dS | phase | F(70 s) | slack cost | `rsK80` 1–4 |
|---|---|---|---|---|
| 2 nm | even (current) | 71.028 | 0.357 | 45.5 48.0 51.7 56.9 |
| 2 nm | odd | 70.937 | 0.449 | 45.8 48.4 52.9 60.5 |
| 1 nm | even | 71.072 | 0.359 | 45.5 47.8 51.1 57.1 |
| 1 nm | odd | 70.977 | 0.397 | 45.5 48.0 52.2 58.6 |
| 0.5 nm | even | 71.017 | 0.384 | 45.3 47.7 51.7 58.1 |
| 0.5 nm | odd | 71.022 | 0.378 | 45.5 47.9 51.7 57.8 |

Width at 2 nm, centred: 60, 90 and 200 bins identical to four digits (F 71.0276, same `ktr80`/`rsK80`).

## Consequences for fitting

- **Old costs are not comparable.** `refit8_opt` recorded ~0.14; under the fix it is **0.356**. Part of what
  refit7/refit8 "fitted" was sampling artefact (e.g. the ladder's segment-5 value moved 0.974 → 0.998
  with dense output alone). Re-score every candidate seed under the new settings before choosing one.
- **Floors built from numerical noise should be rebuilt.** refit8's deadbands add a measured `s_num²`
  (1–9 % per feature under the old settings). Numerical spread is now ≤ 0.8 % of any term; drop or
  re-measure `s_num`, and noise-aware acceptance (`CONFIRM_N`) is no longer needed.
- **Keep the bin count even** and do not "widen" by raising `MaxStrainArraySize` to the `Slim` length —
  widen `Slim_l`/`Slim_r` and keep the cap below it, so the grid stays centred. Verify a final optimum at
  dS = 0.5 nm (≈ 1.8 % on `rsK80` segment 4).

## Physiology, for a separate decision

Not numerical, and not changed here:
- On re-stretch, attached heads reach **+66 nm** of strain; on release, −34 nm. Beyond ~15 nm a head
  should be gone. The fitted rates there are modest: at +15 nm A1 detaches at ~120 s⁻¹, A2 at ~2 100 s⁻¹;
  R2 only exceeds 3 500 s⁻¹ past +30 nm. `UseStrictDetachmentAt` already exists to enforce a strain bound.
- The series elastic element is `Force ∝ LSE^0.845` — an exponent **below 1**, i.e. infinitely stiff at
  zero load and softening with load. Tendon-like series compliance stiffens (toe region). It is also
  non-Lipschitz exactly at slack, though with the rectifier fix `RelTol 1e-4` runs regardless.
- Beyond their knots the rates were arbitrary cubics (now clampable); the knot ranges are
  R1D −51…+50 nm, R12 −1000…+100 nm, R21 −71…+33 nm, R2 −31…+79 nm.

## Files

[ProbeStrainWindow.m](ProbeStrainWindow.m) rates vs strain + shift log ·
[RunSlackNoise.m](RunSlackNoise.m) noise per setting · [ProbeStiffness.m](ProbeStiffness.m) /
[dPUdT_Counted.m](dPUdT_Counted.m) where `ode15s` spends its calls · [ProbeIsometricStall.m](ProbeIsometricStall.m) ·
[ProbeGridConvergence.m](ProbeGridConvergence.m) · [RunObjectiveNoise.m](RunObjectiveNoise.m) ·
`results/` holds every `.mat` and `.log`.
