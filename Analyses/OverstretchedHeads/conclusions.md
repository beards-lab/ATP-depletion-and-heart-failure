# OverstretchedHeads — are the re-stretch heads held too long, and what does fixing it cost?

**Question (fjezek, 2026-09-24).** The current best fit (`refit9_opt`) keeps heads attached
far into the re-stretch. Is that physiological? If not, can earlier detachment be added,
isolated to the re-stretch, and what does it do at high and low ATP? Can a torn-head
state D* (no nucleotide, must bind ATP to reach DT/PT) help?

Base: `params/refit9_opt.m` (8 mM, pooled 3-prep target, objective 0.2079).
Reproduce: `ProbeRestretchDwell.m` → `PlotRestretchDwell.m` → `RunYieldGateVariants.m` →
`BuildYieldSeeds.m` → `RunYieldRefit.m` (TAG `yieldA`, `yieldB`, `yieldC`) →
`EvaluateYieldCampaigns.m` (runs `RunStairsCheck.m`, the held-out staircase).
`ProbeOverstretch.m` / `RunOverstretchVariants.m` are the earlier probe of the same question
(symmetric cut-off, exp-rip D*) and are superseded.

## 0. Verdict (2026-09-25)

* **Not physiological; modify rather than justify.** The long-held heads are *weak* (A1/p1)
  heads dragged to 20–30 nm for ~4 ms, only while lengthening. The data's own re-stretch
  timing is biased the way that predicts (§1–2). The objective could not see it because the
  timing features were guards (§3).
* **No existing parameter fixes it.** A 3 h control re-fit with every attachment and
  detachment lever compulsory left the timing terms unchanged to 3 digits (§5, A).
* **A stretched-side yield gate fixes the physiology and the valley timing**
  (A1 detaches at 6–8 nm instead of 12–24; heads beyond +20 nm fall from 11–31 % to 0–3 %;
  `vall_t` solved; `peak2` onto the data). It works **only with fast rebinding** of the torn
  weak heads. Plain A1 → PD (campaign C) leaves the valley empty and costs 15× more.
* **The gated fit still costs more than the control: 0.763 vs 0.529.** The bill is the valley
  level (too low), a slightly high first peak, and a too-fast post-restretch recovery. On the
  held-out staircase it fixes the nadir and worsens the step amplitude (§7).
  The residual points to one missing mechanism: **lengthening-driven recruitment**. A no-refit
  probe moves every one of those terms the right way but needs a re-fit (§8).
* **D\* carries nothing, at 8 or 2 mM** (`rsK80` 2/8 ×0.98–1.00 against ×0.38 in the data).
  ATP binding is sub-ms at both concentrations, and very little mass reaches the rupture (§6).
* **Adopt:** `UseYieldGate` with `YieldHop1 = 1`, plus a lengthening-recruitment re-fit (§8)
  before this replaces refit9 as the incumbent. **Do not adopt:** D\*, the plain A1 → PD route.
* **Follow-up (2026-09-25, §10):** `peak1_dSL` was a step function because the model is read at
  1 kHz (28 % steps). A sub-sample vertex fixes it and halves the data-side noise. **A1 slip to the
  neighbouring actin site** (the A2 attachment-shift logic, applied to A1) fills the valley:
  slip 5.5 nm with s_yield 12 nm scores **0.501 with no re-fit**, better than the 0.529 control.
  A three-arm campaign for a 16+ core machine is built and pre-flighted (§10c, §11).

---

## 1. The diagnosis: it is the weak heads, and only while lengthening

`ProbeRestretchDwell.m` rebuilds strain, masses and every detachment flux per bin from the
saved state with the RHS's own rate functions (the rebuilt XB force matches the logged one to
5e-16). Slack protocol, re-stretches 1–5:

| seg | p1 detach at (flux-weighted) | p2 detach at | p1 life past +15 nm | p2 life past +15 nm | XB force from >15 nm heads, max in 40 ms |
|---|---|---|---|---|---|
| 1 | 16.5 nm | 10.2 nm | 4.3 ms | 0.27 ms | 72 % |
| 2 | 19.5 | 11.1 | 4.0 | 0.26 | 77 % |
| 3 | 22.0 | 12.0 | 3.7 | 0.26 | 79 % |
| 4 | 23.9 | 12.7 | 3.4 | 0.25 | 81 % |
| 5 | 11.8 | 8.9 | 5.2 | 0.28 | 57 % |

* **The strong (p2) heads are fine.** They let go at 9–13 nm and survive 0.26 ms beyond 15 nm:
  `R2(s)` is already 1400–2100 s⁻¹ at 10–15 nm.
* **The weak (p1) heads are the problem.** Strained p2 heads reverse into p1
  (`R21` 700–1700 s⁻¹ beyond 15 nm), and p1 then detaches at only 120–300 s⁻¹ (`R1D`)
  while `R12` sends part of it back: a p1↔p2 shuttle that holds heads for ~4 ms and drags
  them to 20–30 nm. p1 outnumbers p2 at the peak (p1:p2 1.0–1.3).
* **They carry the valley and the second rise, not the first peak.** Heads beyond +15 nm
  carry 0–1 % of XB force *at* the first peak but 57–81 % within the next 40 ms
  (`results/restretch_dwell_refit9_opt.png`).
* **It is a lengthening-only phenomenon.** Share of attached heads beyond +15 nm on the
  *stretched* side: isometric plateau 0.0 %, slack shortening 0.0 %, re-stretch ≤ 66 %,
  ktr overstretch ≤ 82 %. So a mechanism acting on positive strain beyond ~12 nm is
  isolated to re-stretch protocols (slack re-stretch, ktr overstretch, staircase) by
  construction. It cannot touch the isometric operating point or F–V.

## 2. Is it physiological? No, and the data already say so

* Literature: when highly strained bridges pass a critical extension they detach fast
  (Lombardi & Piazzesi 1990; Campbell & Lakie 1998; the SRS elastic-limit literature puts
  the extension at ~8–10 nm per half-sarcomere, from a review summary I did not verify
  against the primary papers). In steady lengthening the force reaches a plateau after
  ~20 nm/hs, held up by rapid detach–reattach cycling. A *weak* bond in particular is in
  rapid equilibrium with actin (Brenner: k₊ ≫ k₋, dissociation 10²–10³ s⁻¹ even in
  force-generating states). A weak bond that holds 15–30 nm for 4 ms is outside all of this.
* **The data's own timing points the same way.** The three re-stretch timing features are
  biased in the direction late yield predicts, on every segment:

  | feature | model (seg 1–5) | pooled data | mean relative bias |
  |---|---|---|---|
  | `peak1_dSL` (µm SL to first peak) | 0.024–0.042 | 0.023–0.031 | **+20 %** (peaks later) |
  | `vall_t` | 10–14 ms | 6–11 ms | **+34 %** (valley later) |
  | `restretchSlopeStart` | 915–1056 | 1107–1571 | **−28 %** (rises softer) |

  The data rise stiffly, peak at 11.5–15.5 nm/hs of stretch (series compliance included,
  so the XB extension is less), then drop within a few ms to a flat plateau. That is the
  fast-yield / fast-reattach picture above. The model rises softly, yields late, and holds
  force up with stretched weak heads.

## 3. Why the fit is allowed to do this: the timing trio is a guard, not a term

`BuildRefit9Seed` scores `peak1_dSL`, `vall_t`, `restretchSlopeStart` (and the force shape) as
**guards**, floored at 1.2× the *seed's own* residual. At `refit9_opt` their raw χ² are
0.43 / 0.60 / 0.45 against floors of 0.52 / 0.72 / 0.54, so all three pay **exactly zero** while
sitting an RMS 1.5–1.8 preparation-widths off per segment, all with the same sign. The design froze whatever timing the seed had.
On the honest rule every ATP term already uses (fl = nEl·s²/3, pooled s 0.19–0.20) the
floors are ~0.06, and refit9_opt scores **0.558 instead of 0.208**. The timing trio alone
is 0.35 of that.

## 4. Earlier detachment without a re-fit: the valley collapses

Wired as `UseYieldGate` (default off, inert: `as fitted` = 0.2079 exactly). A bounded
logistic in s, stretched side only: `R1D += k_yield1·g(s)` (p1 → PD),
`XB_Ripped += k_yield2·g(s)` (p2 → PT, or → D0 when `UseD0State`),
g = 1/(1+e^{−(s−s_yield)/w_yield}). `YieldHop1` (added below) lets the p1 share rebind at zero
strain instead. `RunYieldGateVariants.m`, full objective at 8 mM unless stated, nothing re-fitted:

| variant | cost | peak1_y seg1 | vall_y seg1 | vall_t | peak2 seg1 | rsK80 | trio raw χ² |
|---|---|---|---|---|---|---|---|
| **data** | | **80.1** | **63.4** | 6–11 ms | **69.5** | 46–54 | |
| as fitted | 0.208 | 82.4 | 56.7 | 10–14 | 65.3 | 44–57 | 0.43 / 0.60 / 0.45 |
| p1 gate 15 nm, 1000/s | 9.87 | 81.8 | 28.7 | 6–10 | 48.7 | 43–53 | 0.27 / 0.19 / 0.45 |
| p1 gate 12 nm, 1000/s | 13.0 | 80.3 | 25.6 | 6–9 | 45.1 | 44–52 | 0.13 / 0.31 / 0.47 |
| p1 gate 10 nm, 2000/s | 20.3 | 77.6 | 17.9 | 5–7 | 39.2 | 45–52 | 0.25 / 0.78 / 0.58 |
| p2 gate 12 nm, 3000/s → PT | **0.78** | 79.6 | 51.3 | 9–13 | 61.1 | 47–58 | 0.13 / 0.29 / 0.50 |
| p1+p2 gate 12 nm → PT | 19.2 | 78.4 | 17.8 | 4–7 | 40.6 | 45–52 | 0.26 / 0.93 / 0.56 |
| p1 gate 12 nm 1000/s, `A_def` off | 9.86 | 80.5 | 34.4 | 6–8 | 52.8 | **71–76** | |
| **p1 HOP 12 nm, 2000/s** (rebind at 0) | **2.41** | 80.6 | 46.9 | 7–10 | 57.0 | 49–60 | |

1. **The timing moves the right way immediately.** `peak1_dSL` and `vall_t` drop into the data's range.
2. **The first peak barely moves** (82.4 → 78–81). The premise "earlier detachment brings the
   peaks down" holds only weakly, because the first peak is carried by heads *below* 15 nm (§1).
3. **The valley collapses** (56.7 → 18–29 kPa; data 63). Heads that let go are not replaced:
   attachment while lengthening runs at `ka·A_reg` ≈ 292 × 0.21 ≈ **61 s⁻¹**
   (`RegAvailShorteningOnly` pins availability at A0 during lengthening), and the load-realign
   deficit `A_def` adds to it (more mass stripped under load → less attachment).
   Switching `A_def` off only half-rescues the valley (25.6 → 34.4) and breaks the
   post-restretch rate (`rsK80` ×1.5). So `A_def` is not the main brake; slow attachment is.
4. **Rebinding is what the physiology says and what the fit needs.** Letting stretched weak
   heads rebind near zero strain (rapid-equilibrium weak binding, mass-conserving, so no
   `A_def`) holds the valley at 46.9. At 2000 s⁻¹ it costs 2.4 against 13–22 for plain
   detachment at the same 12 nm (1000–5000 s⁻¹), with nothing else touched.
5. **The stretched strong heads are nearly free to cut** (p2 gate: 0.21 → 0.78).

**Read:** the long-attached weak heads are not a free-standing error. They are how this fit
holds force up after the yield, because re-attachment during lengthening is too slow.
Earlier detachment is only viable together with fast rebinding. Whether the rest of the
parameters can absorb the change needs a re-fit (§5).

## 5. The controlled re-fit

Three campaigns, one objective (`yieldA_seed`'s fn: the trio floored at nEl·s²/3), the same
optimiser and draw sizes as refit9. The compulsory set in every round is the re-stretch
machinery: `ka, kd, k_1, kstiff1, kstiff2, k_cr, tau_cr, A0` (+ `s_yield` in B, C). Gate
rates are fixed at 2000 / 3000 s⁻¹, so a fit cannot switch the premise off by shrinking them.
All four fits re-scored in one process by `EvaluateYieldCampaigns.m` (`results/evaluate.log`):

| fit | stretched A1 heads | seed | final | rounds | wall |
|---|---|---|---|---|---|
| refit9_opt | (incumbent, no gate) | | 0.558 | | |
| **A** `yieldA_opt` | no gate (control) | 0.558 | **0.529** (stalled) | 2 | 3 h |
| **B** `yieldB_opt` | yield gate, A1 **rebinds** at zero strain | 4.287 | **0.763** (s_yield 14.1 nm) | 9 | 15 h |
| **C** `yieldC_opt` | yield gate, A1 → **PD** (no rebinding) | 19.59 | **11.54** | 3 | 3 h |

B's trajectory 4.29 → 2.29 → 1.48 → 0.91 → 0.86 → 0.83 → 0.79 → 0.78 → 0.76 → 0.76 has
flattened (last round −0.4 %). Overnight throughput was ~3× below expectation because another
MATLAB session (not part of this analysis) held a core all night.

**A: the existing parameters cannot move the yield timing.** After 3 h the trio is where it
started: `vall_t` 0.146 → 0.146, `peak1_dSL` 0.103 → 0.103, `restretchSlopeStart` 0.101 →
0.099. The whole 0.030 gain is in ktr / rsK80 / rsT20 / steady. So the late yield is structural
to this parametrisation, not a search failure.

**B: the gate buys the timing and the physiology, and pays in the valley level:**

| term | A | B | B, the feature |
|---|---|---|---|
| `vall_t` | 0.146 | **0** | 11/10/8/7/10 ms, data 11.3/10.1/7.9/6.0/10.3 |
| `restretchSlopeStart` | 0.099 | 0.077 | 974–1094, data 1107–1571 (A 921–1059) |
| `peak2` | 0 | 0 | 64–70, data 70 (A overshoots seg 3/4 by +6 / +11 %) |
| `peak1_dSL` | 0.103 | 0.090 | unchanged on seg 1–4 (caveat below) |
| `vall_y` | 0 | **0.191** | 51–54 on seg 1–4, data 57–63 |
| `peak1_y` | 0 | **0.081** | 74–86 on seg 1–4, data 74–80 |
| `rsK80` + `rsT20` | 0.078 | **0.242** | recovery too fast: rsT20 4.7–6.6 ms, data 8.2–9.2 |
| total | 0.529 | 0.763 | |

Physiology, `ProbeRestretchDwell` on the final fits (`results/dwell_*.mat`):

| | refit9_opt | A | **B** | C |
|---|---|---|---|---|
| A1 detach strain (flux-weighted) | 12–24 nm | 12–24 nm | **6–8 nm** | 11–15 nm |
| A2 detach strain | 9–13 nm | 9–13 nm | 7–10 nm | 7–10 nm |
| A1 life past +15 nm (detach *or* rebind) | 3.4–5.2 ms | 3.5–5.2 ms | **0.53–0.57 ms** | 0.51–0.55 ms |
| attached beyond +20 nm, re-stretch | 11–32 % | 11–32 % | **0–3 %** | 3–14 % |
| XB force from > 15 nm heads, max in 40 ms | 57–81 % | 58–82 % | **21–46 %** | 32–78 % |
| ktr overstretch, beyond +15 nm | 82 % | 82 % | 67 % | 70 % |
| isometric / shortening, stretched side | 0 % | 0 % | 0 % | 0 % |

`results/restretch_dwell_refit9_vs_yieldB.png` shows it in the trace. B drops right after
the peak, when the data does, not in refit9's slow sag, and its stretched-head load is a
fraction of refit9's. But through the rest of the ramp the data holds a high plateau (~80 kPa
on seg 1) that neither model reaches. That is the missing lengthening recruitment of §8.

**C: the literal A1 → PD route does not work.** The valley cannot be refilled through the
rate-limited `ka·PD` attachment (`A_reg` at A0 during lengthening, `A_def` on top): `vall_y`
22–24 against 57–63, cost 11.5. The fit's answer was to keep C's heads longer (A1 detach
strain back up to 11–15 nm, force from > 15 nm up to 78 %). Rebinding (rapid-equilibrium weak
binding) is what makes earlier detachment viable, not a detail.

**Caveat on `peak1_dSL` (corrected 2026-09-25, see §10a).** The model trace is read at 1 kHz
by `extractSlackAttributes` (fixed output grid → resampled to TARGET_FS), so the sample peak
moves in 1 ms steps = 0.0078 µm of the 7.79 µm/s ramp, **28 %** of the feature's mean (not the
13 % first written here). It reads identically to 4 digits in refit9, A and B on segments 1–4:
the peak lands on the same sample. The cost is a step function the optimiser cannot descend,
and the same fixed tax in every arm. §10a replaces it with a sub-sample version.

## 6. Low ATP and D\*

The data's 2 mM / 8 mM ratios on the re-stretch are **`rsK80` ×0.38** (paired t −26),
`ktr` ×0.53–0.57, `peak1_dSL` ×1.17, `vall_t` ×0.79 (`KtrSegmentSlope/results/atp_relevance.mat`).
So at low ATP the yield comes *later*, and the post-restretch recovery slows about a third more
than the post-slack one. That is the restretch-specific ATP signal D\* was meant to carry.

The p2 rupture flux routed through D0 (`UseD0State`, `k_D0T` 500 s⁻¹, exit ∝ [MgATP]/(K+[MgATP])):

| | K_D0T 0.15 mM (physiological) | K_D0T 5 mM (phenomenological) | data |
|---|---|---|---|
| no-refit sweep, gate 12 nm: `rsK80` 2/8 | ×1.00 | ×0.99 | **×0.38** |
| **fitted B**: `rsK80` 2/8 | ×0.999 | ×0.980 | **×0.38** |
| fitted B: `ktrfit`, `peak1_dSL`, `vall_t`, `steady` 2/8 | ×1.000 | ×1.000 | ×0.55 / ×1.17 / ×0.79 / ×1.34 |

**D\* contributes no ATP effect**, even with a K_D0T 30× above the measured ATP-binding K_m.
Nothing else in these fits depends on [MgATP], so every ratio here is D\*'s alone. Two reasons:
* **Kinetic.** ATP binding to a nucleotide-free head is > 10⁴ s⁻¹ at 8 mM and ~3–4·10³ s⁻¹ at
  2 mM (α-MHC, K₁k₊₂ 1.6–2.0 µM⁻¹s⁻¹, `parameterBounds.m`). That is a dwell of 0.1–0.3 ms,
  invisible on the ms re-stretch timescale at either concentration. `ApoPoolDetachment`
  reached the same bound.
* **Structural.** Very little mass goes through it. Strained A2 heads leave by `R2` (> 1400 s⁻¹
  past 10 nm) or reverse into A1 long before the rupture strain, so D0 is fed by a trickle.
  Also chemically, an A2 head torn off still holds ADP; only the rigor (post-ADP) fraction
  would be apo. That fraction is lumped into p2 in the 2-state tree and cannot be separated.

**Where the low-ATP re-stretch effect has to come from instead.** The later yield at 2 mM
(`peak1_dSL` ×1.17) means the yield itself must be ATP-sensitive, i.e. slower
chemical detachment of the stretched heads. That is `UseAtpDetachR2` territory, acting on
the heads *before* they reach `s_yield`, not a refractory pool after it. The fixed-strain gate
is ATP-blind by construction, so it will not produce that ratio on its own. That is correct: the
mechanical limit should not depend on [ATP]. Test it on the gated fit with the low-ATP gate
vector, not with D\*.

## 7. Held-out staircase (never in any objective)

`RunStairsCheck.m`, 03/27 8 mM staircase, window-mean features (F₀ units):

| | E_stairs | step amplitude `st_peak` (1 %) | nadir `st_nadir` | end plateau |
|---|---|---|---|---|
| **data** | | 0.073–0.107 | −0.019 … −0.034 | 0.191 |
| refit9_opt | 0.0372 | 0.130–0.157 | −0.085 … −0.112 | 0.330 |
| A | 0.0359 | 0.133–0.159 | −0.082 … −0.108 | 0.329 |
| **B** | 0.0445 | 0.167–0.188 | **−0.014 … −0.025** | 0.329 |
| C | 0.0672 | 0.152–0.177 | −0.050 … −0.058 | 0.357 |

**Mixed, and informative.** The gate fixes the phase-2 undershoot, which was 3× too deep in
every ungated fit and is now at the data (slightly shallow: −0.014…−0.025 vs −0.019…−0.034). That is an independent confirmation, on
a protocol the fit never saw, that the stretched heads were the cause of the undershoot. The
step amplitude gets *worse* (1.6× → 2.0× the data), because B's first peak is higher (§5).
The same shortfall shows in slack `peak1_y`, so both protocols point to the same fix. The end
plateau (`0.33` vs `0.19`) is untouched by anything here; it is a length–tension /
passive property, not a re-stretch one.

## 8. What B is missing: lengthening recruitment (probe, no re-fit)

B's whole residual (valley too low, recovery too fast, initial slope too soft) is what you get
when heads let go at the right strain but **nothing extra binds during lengthening**: attachment
runs at `ka·A0` ≈ 61 s⁻¹ whenever the fibre is not shortening (`RegAvailShorteningOnly = 1`).
Physiologically, lengthening recruits the second myosin head (Brunello et al. 2007, PNAS) and
the steady-lengthening force is carried by fast detach–reattach cycling (Lombardi & Piazzesi
1990). `RegAvailLengtheningGain` (default 0 = shortening-only exactly) lets lengthening count at a
fraction of its speed in the availability target. On `yieldB_opt`, no re-fit:

| gain | cost | `peak1_y` | `vall_y` | `restretchSlopeStart` | `rsK80` | `rsT20` |
|---|---|---|---|---|---|---|
| **data** | | 74–80 | 57–63 | 1107–1571 | 46–54 | 8.2–9.2 ms |
| 0 (= B) | 0.763 | 74–86 | 51–54 | 974–1094 | 41–56 | 4.7–6.6 ms |
| 0.03 | 8.98 | 89–100 | 75–79 | 1143–1253 | 36–43 | 7.2–8.4 ms |
| 0.08 | 26.2 | 103–117 | 89–96 | 1247–1389 | 31–34 | 8.8–9.8 ms |
| 1 (symmetric) | 72.6 | 116–147 | 95–115 | 1389–1667 | 29–30 | |

Every one of B's residual terms moves the right way (valley up, slope up, recovery slower), and
at gain 0.03 `rsT20` is nearly on the data. The lever is strong, because availability sits at A0 = 0.21,
so the no-refit dose overshoots the peaks. Interpolated, the valley lands on the data near gain ≈ 0.01.
F–V is untouched (shortening-only path, FV cost 0 at gain 1).

## 9. Recommendation

1. **Next campaign, `yieldD`:** seed `yieldB_opt` + `RegAvailLengtheningGain` ≈ 0.01, with
   the gain in the pool and compulsory, same objective, ~6–8 h on a quiet machine. Expected: `vall_y`,
   `peak1_y`, `rsK80/rsT20` and `restretchSlopeStart` all pull toward the data together. This is
   the one mechanism the residual names. Success criterion: at or below A's 0.529 with the §5
   physiology intact, and the staircase step amplitude not worse than refit9's.
2. **Objective hygiene before any re-fit.** Keep the trio on the nEl·s²/3 floor (not guards).
   Make `peak1_dSL` resolvable (finer `odeOutputDt` in the re-stretch window, or sub-sample peak
   location), otherwise it is a flat 0.09–0.10 tax.
3. **Low ATP:** evaluate the gated fit at 2 mM with the ATP-gate vector used in
   `StairsPrediction` (the gate itself stays ATP-blind). Do not pursue D\*.
4. **Model changes (uncommitted, all opt-in, exact no-ops when off):** `UseYieldGate`,
   `s_yield`, `w_yield`, `k_yield1`, `k_yield2`, `YieldHop1`, `RegAvailLengtheningGain` in
   `Model/dPUdT_CombinedTransitions.m` + `Model/getParams.m`. Regression: refit9_opt with all
   off = 0.2079, unchanged. Not in the MEX path (`UseCompiledMex = 0` everywhere here).

**Pitfall recorded for anyone scripting around `RunBakersExp`:** it runs in the caller's
workspace and overwrites `E`, `v`, `s`, `b`, `h`, `out`, `params`, `features_*` and more.
A results struct named `E` is silently replaced (this crashed the first evaluation run).

## 10. Follow-up (2026-09-25): a resolvable `peak1_dSL`, A1 slip-rebinding, and the campaign for a big machine

### 10a. `peak1_dSL` is `peak1_t`, and the hard cost was the read, not the metric

`peak1_dSL = SL(t_peak) − SLslack`, and each re-stretch is one constant-velocity ramp, so it is
`peak1_t` × 7.79 µm/s. It is still the right observable, the stretch at which the heads let go
plus the series compliance. The problem was how it was read. The model trace goes through
`extractSlackAttributes` at **1 kHz** (fixed output grid → resampled), so the sample peak moves in
1 ms = 0.0078 µm steps, **28 % of the mean**. The data are read raw at 4 kHz and take one noisy
sample as the peak.

**Fix (cheap, no solver change):** `peak1_tq`, `peak1_dSLq`, `vall_tq` = the vertex of a
least-squares quadratic over ±1 ms (peak) / ±2 ms (valley) around the sample extremum, on data
and model alike (`vertexTime` in `extractSlackAttributes.m`; the old features are kept).
`TestPeakSubsample.m`:

| test | sample-based | vertex |
|---|---|---|
| **noise**: smooth model trace + the data's white noise (σ 1.3 kPa, 4 kHz, 60 reps), `peak1_dSL` jitter / bias | 10.0 % / +2.1 % | **5.2 % / +0.3 %** |
| same, `vall_t` | 6.8 % / +0.5 % | **3.5 % / −0.3 %** |
| **smoothness**: `kstiff2` ×0.98…×1.02 (9 points), `peak1_dSL` | constant (35.60 at all 9) | **monotone 34.56 → 33.97** |
| **data**, between-prep rel. SD, `peak1_dSL` | 0.179 | 0.201 (real prep differences: 04/10 peaks later) |
| same, `vall_t` | 0.191 | **0.119** (a third of the spread was estimator noise) |

The pooled target now carries the three vertex features (`UpdatePooledTargetVertex.m`: re-extract,
stage, diff). **Every existing value is unchanged** (max relative change 0). Moderated spreads
drift ≤ 1.5 % through the shared variance prior, which nothing in the objective reads. Backups:
`Analyses/DataPooling/results/*_20260925_1306.bak`. The slip campaign scores
`peak1_dSLq|fl:0.0649|0.273` and `vall_tq|fl:0.0523|0.339` in place of the sample-based pair.

### 10b. A1 slip to the neighbouring actin site (the A2 attachment-shift logic, applied to A1)

`YieldHopShift` > 0: the yield-gated A1 heads that rebind (share `YieldHop1`) land
`YieldHopShift` back along the filament, still strained, split over the two neighbouring
bins as in `UseA2AttachmentShift`. The rest (`1 − YieldHop1`) unbinds to PD. Mass-conserving, and
**neither path uses ATP**: an A1 head holds ADP·Pi, so slipping keeps it and unbinding returns it to
PD, the pre-hydrolysed state. The extra re-stretch binding therefore costs no turnover, in line with
lengthening being energetically cheap. `RunA1SlipVariants.m`, from `yieldB_opt`, no re-fit, the
campaigns' common objective (sample-based timing terms):

| variant | cost | peak1_y | vall_y | peak2 | peak1_dSLq | rsK80 |
|---|---|---|---|---|---|---|
| **data** | | 74–80 | 57–63 | ~70 | 0.025–0.031 | 46–54 |
| yieldB_opt (land at 0) | 0.763 | 74–86 | 51–54 | 64–70 | 0.033–0.036 | 41–56 |
| slip 5.5 nm, all rebind | 0.785 | 76–86 | **60–62** | 68–79 | 0.035–0.040 | 39–50 |
| slip 5.5 nm, ¾ rebind ¼ PD | 3.08 | 75–86 | 40–42 | 59–62 | 0.034–0.037 | 38–46 |
| slip 5.5 nm, ½ rebind ½ PD | 7.20 | 74–85 | 30–32 | 53–55 | 0.033–0.036 | 37–45 |
| slip 11 nm, all rebind | 0.630 | 75–86 | 55–56 | 65–74 | 0.033–0.038 | 40–54 |
| **slip 5.5 nm, s_yield 12 nm** | **0.501** | 72–84 | 55–57 | 64–74 | 0.032–0.035 | 41–56 |

* **Slipping fills the valley** that zero-strain rebinding left empty, because the head keeps
  part of its extension. With a slightly lower yield strain it beats the no-gate control's 0.529
  **with no re-fit** (0.501).
* **Unbinding to PD collapses the valley** at any appreciable share. Even 2 % moves `vall_y`
  from < 0.08 (all rebind, where it is not among the top three terms) to 0.17 (slipA seed, below). PD heads wait on the slow lengthening attachment
  (`ka·A_reg` at A0). So A1 → PD raises re-stretch attachment only *together with* faster
  attachment while lengthening, which is arm B.
* What slip does not fix: the initial slope and the too-fast post-restretch recovery
  (`rsK80`, `rsT20`). That is the lengthening-recruitment residual of §8, and the reason arm B exists.

### 10c. The campaign (for a 16+ core machine)

`RunSlipCampaign.m`, one MATLAB session per arm, 5 workers each (a session saturates at 5–6: the
parallelism is over the five slack segments; 16 cores = 3 arms). One objective, `BuildSlipSeeds.m`
(the trio floored + the vertex pair), pre-flight re-scores each seed against its recorded cost:

| arm | seed | seed cost | free on top of the usual pool | compulsory adds |
|---|---|---|---|---|
| **a** | yieldB_opt + A1 slip 5.5 nm, s_yield 12 nm, YieldHop1 0.98 | 0.585 | YieldHop1, YieldHopShift | s_yield, YieldHop1 |
| **b** | a + `RegAvailLengtheningGain` 0.005 | 0.768 | + the gain | + the gain |
| **c** | yieldA_opt, no gate (control) | 0.546 | — | — |

Every arm also has `ka, kd, k_1, kstiff1, kstiff2, k_cr, tau_cr, A0` compulsory. Gate rates stay
fixed (2000 / 3000 s⁻¹), so the premise cannot be switched off. The per-solve guard is 600 s
(gated arms ran 210–380 s per evaluation under load here; 300 s tripped the serial coverage check).

```
# Linux / node          bash Analyses/OverstretchedHeads/launch_slip.sh              # pre-flight
                        GO=1 HOURS=20 bash Analyses/OverstretchedHeads/launch_slip.sh
# SLURM                 sbatch --account=<acct> Analyses/OverstretchedHeads/slip.sbatch
# Windows               .\Analyses\OverstretchedHeads\launch_slip.ps1 -Hours 20 -Go
# afterwards            matlab -batch "SNAPS={'yieldA_opt','slip_a_opt','slip_b_opt','slip_c_opt'}; run('Analyses/OverstretchedHeads/EvaluateYieldCampaigns.m')"
```
**Tested here (2026-09-25):** all three arms pre-flight OK through `launch_slip.ps1`, each seed
re-scoring to its recorded cost to 4 digits (0.5848 / 0.7678 / 0.5457); `launch_slip.sh`
exercised with a stub MATLAB (argument passing, modes, logs, summary); **smoke round of arm a**
end to end (state file, snapshot, clean-up). The smoke round already moved **0.5848 → 0.4939**
in its initial simplex, below the control seed (0.5457) on the same objective. (A smoke round is
~12 evaluations, not 4: `fminsearch` scores the whole initial simplex before it honours the cap.)
Budget guide: an evaluation is ~70–100 s on a quiet 5-worker session, so a 30-evaluation round is
~45 min, and 20 h is ~25 rounds per arm.

Resume is automatic (re-run the same command). The other machine needs the working-tree files
listed in §11. They are uncommitted here, and `data/` (not in git) needs the 03/27 slack, ktr
and passive protocol files.

## 11. Moving the campaign to another machine

The campaign runs on the **current working tree**, which carries uncommitted changes, both mine
and ones that predate this analysis (e.g. `Model/evaluateModel.m`, `Model/experiments/*`). Sync
the whole tree (commit + pull, or rsync), not a hand-picked subset. This analysis added or changed:

| file | what |
|---|---|
| `Model/dPUdT_CombinedTransitions.m` | `UseYieldGate` (+ `YieldHop1` split, clamped to [0,1]; `YieldHopShift` slip landing), `RegAvailLengtheningGain`; all exact no-ops at their defaults |
| `Model/getParams.m` | defaults for the above |
| `Model/extractSlackAttributes.m` | `peak1_tq`, `peak1_dSLq`, `vall_tq` (`vertexTime`); existing features unchanged |
| `Model/pooledTarget8mM.m` | regenerated with the three vertex features; existing values identical |
| `Analyses/DataPooling/buildPooledTarget.m` | FEAT list + the three features |
| `params/slip{A,B,C}_seed.m` | the arm seeds; `params/yield{A,B,C}_opt.m`, `yieldA_seed.m` for the evaluation |
| `Analyses/OverstretchedHeads/*.m, *.sh, *.ps1, *.sbatch` | the analysis, the runner and launchers |

`data/` (gitignored) must hold `protocol_03_27_2026_8mM_slack.mat`,
`protocol_03_27_2026_velocitytable_ktr.mat`, `protocol_03_27_2026_ActivePNBMava_slack.mat`
(the pre-flight checks all three and stops if any is missing). The pre-flight also re-scores each
seed and refuses to start if it is more than 5 % off the cost recorded here: a different MATLAB
release shifts this objective by ~1–2 %, and a larger gap means code or data differ.
