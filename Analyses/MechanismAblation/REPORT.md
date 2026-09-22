# The lean model — methods, mechanisms, and what it does and does not predict

**2026-09-15.** Built on the 35-arm mechanism ablation (`conclusions.md`).
Snapshot **`params/lean_opt.m`**, 2 h refit, 14 rounds.
Figures: `Analyses/OptimResultReport/results/report_lean_opt.html` (12 inlined).

---

## 1. Headline

| | |
|---|---|
| lean objective, seed → refit | **0.8224 → 0.5288** (−35.7 %), still descending when the clock stopped |
| same objective, incumbent vs refit | **0.7598 → 0.5311** (−30.1 %) |
| staircase (never scored) | **clearly better** — peak heights now on the data, where the incumbent overshot every step |
| kʈr (never scored) | **worse** — recovery too fast, re-stretch peak 1.17 vs data 0.72 |

The refit bought a large, real gain on everything it was shown and on one
protocol it was not shown — and paid for part of it on another protocol it was
not shown. That trade is the most important thing in this report and it is
discussed in §5.

---

## 2. Methods

**Objective.** Feature-based (`OptimizeOn = 'Feats'`), 24 entries, each a
relative chi-square against the **pooled three-preparation 8 mM target** with a
weight `λ/s²` from the measured between-preparation spread, and a **one-
preparation-width deadband**: a term inside its own noise pays nothing but still
defends itself against regression.

Three changes from the incumbent objective, all earned by the ablation:

1. **`coolDownLS` removed.** An unbounded raw MSE (`extractSlackAttributes.m:456`),
   not a chi-2, it was 16.3 % of the incumbent cost and dominated the movers of
   nearly every ablation arm (−706 on one arm while data terms moved by 0.01).
   A term that can swing by 700 in a 0.9 objective is a guard rail, not a
   measurement.
2. **`rsK80` weight doubled** to inherit the re-stretch-recovery role. *In
   practice this changed nothing*: `rsK80` is already inside its deadband and
   pays 0 in both fits. The honest description is that `coolDownLS` was removed
   and nothing needed to replace it.
3. **F-V weight left alone.** Reducing the 0.40 pin was requested and is a
   no-op: under the deadband F-V sits at chi-2 0.0515 against a floor of 0.0515
   and **already pays exactly nothing**. Its weight only bites if F-V regresses
   out of its own noise band, which is the deadband's intended behaviour.

**Search.** `optimizeFeatures`, 2 h wall-clock, rounds of 30 simplex evaluations
on a 12-parameter draw (7 compulsory: `k2, kstiff2, ka, v_ref_reg, k_cr, eta_M,
kmsrd`). Short rounds by design — the budget is only checked *between* rounds,
so a single long round cannot be stopped on time.

**Pool: 70 → 59 entries.** Removed as inert (`sigma1`, `sigma_srd2`,
`LXBpivot`), as belonging to deleted mechanisms (`kA2re`, `kA2hop`, `sA2hop`,
`x_M_slack`), and as redundant (`d10_ref`, `SL_ref_lattice`, `R_thick`,
`R_thin`, `L_hbare` — see §4).

**Reproducibility.** The objective is bit-reproducible on an idle machine
(4 evaluations, spread exactly 0). Its only non-determinism is the `MaxRunTime`
watchdog, a **wall-clock-valued ODE event** that terminates the integration and
therefore fires on machine load rather than on model behaviour. Separately,
snapshots are written with `%g` — six significant digits — which moves the cost
by ~0.004. Both are visible here: the optimiser reports 0.5288, the re-evaluated
snapshot scores 0.5311.

---

## 3. Mechanisms in the lean model

**Active, and demonstrated load-bearing by ablation** (cost of removing, after a
budget-matched re-optimisation):

| mechanism | parameters | cost to remove |
|---|---|---|
| registration availability (`UseRegistrationAvailability`, shortening-only) | `A0`, `v_ref_reg`, `tau_reg` | **+793** |
| free-form R1D strain law (`UsePieceWiseStrainDep`) | 9 knots | **+82** |
| super-relaxed pool (`UseSuperRelaxed` + `ADP`) | 8 | **+10.1** |
| passive Maxwell element (`UseMaxwellDashpot`) | `kSE_M`, `eta_M` | +0.97 |
| global occupancy saturation | `P_bound_max` | +1.22 |
| load realignment (`UseLoadRealign`) | `k_cr`, `tau_cr`, `F_cr` | +0.73 |
| filament overlap (`UseOverlap`) | `L_thick`, `L_thin` | +0.46 |
| negative-strain stiffness | `kstiff1_n`, `kstiff2_n` | +0.25 |
| lattice spacing | `d_optimal`, `sigma_lattice` (§4) | +0.12 |
| serial elastic element | `kSE`, `ekSE` | always on — its flag is dead code |

**Deleted, free:** `UseA2Reattaching` (3 parameters, +0.0001),
`UseMaxwellTensionOnly` (1 parameter, +0.037 and `rsK80` 0.06 *better*),
`UseDynamicPassive` (exactly inert).

**Tested and rejected** (adding them costs): stretch activation +19.0,
velocity-Gaussian attachment as a substitute for availability +746, `F_passive`
as the SRX driver +3.2, A1 attachment kernel +0.40, catch bonds +0.15,
R2 ceiling +0.09, direct SRX routing +0.09.

**Thick-filament topology, as it actually is.** `kmsr = 0` and `ksrd = 0`, so
the model has no direct force-dependent SRX→DRX return and no PD→SRD route:
`PT →(ksr0)→ SRX ⇄ SRD →(kmsrd)→ PD`, a one-way cascade. Reviving either route
(R1, R2) lands inside the noise, so this is not a defect — but it was never a
choice either, and it is worth stating in any write-up of the model.

---

## 4. The lattice mechanism, reduced

`f_lattice` is written with six parameters and has **two**:

```
f_lattice = exp( −(1/√SL − a)² / (2b²) )
     a = (R_thick + R_thin + d_optimal) / (d10_ref·√SL_ref_lattice)
     b = sigma_lattice / (d10_ref·√SL_ref_lattice)
```

`d10_ref·√SL_ref_lattice` is a product and `R_thick + R_thin + d_optimal` a sum;
neither is separable. What the mechanism says is **attachment is attenuated
below SL ≈ 2.09 µm**, roll-off width 0.047 in 1/√SL — inert at the 2.2 µm
operating point, `f ≈ 0.95` at 2.0, `≈ 0.52` at 1.8.

The four redundant parameters are exactly the four that are *observations* —
X-ray `d10`, the length it was measured at, and two EM filament radii. They are
now fixed constants. `d_optimal` and `sigma_lattice`, the two genuinely
unmeasured ones, remain free. The mechanism therefore costs **2 parameters for
+0.12**, not 6, and should be kept.

---

## 5. Model performance

### What the model now gets right

- **Absolute force** moved to 72 kPa against a three-preparation mean of 66
  (the incumbent sat at 81, i.e. on the single preparation it was fitted to).
- **`t0`, the force-onset timing** — on the pooled mean across all five cycles.
- **`rsK`, post-re-stretch recovery rate** — inside the envelope, *and rising
  across segments as the data does*. This was the open structural defect for
  months (the model's amplitude slope had the wrong sign); it no longer does.
- **`A`, `Am`, `peak2`, `vall2_t`, `vall2_dy`, `restretchSlopeStart`** — inside
  the three-preparation envelope.
- **The staircase, unscored:** step peaks now sit on the data. The incumbent
  overshot every one.

### What it still gets wrong

| feature | model | data (3-prep mean) | |
|---|---|---|---|
| `peak1_dSL` | 0.034–0.040 | 0.008–0.030 | outside envelope — **31 % of the remaining cost** |
| `peak1_t` | 4.4e-3 | 1.1–2.7e-3 | first peak arrives too late |
| `vall_t` | 0.010–0.014 | 0.004–0.008 | valley too slow |
| `peak1_y` | 1.27 | 1.02 | re-stretch overshoot, cycle 1 |
| `vall_y` | 0.79 | 0.93 | valley too deep — *worse* than the incumbent |
| `ktr` | 56 | 51 | recovery too fast |
| `rsT0` | 0.001–0.002 | 0.004–0.007 | outside envelope |
| `rsR2` | 1.00 | 0.59–0.80 | the model's recovery is a perfect exponential; real ones are not |

The remaining cost is now **peak1_dSL (31 %)** plus the **passive terms
(36 % combined:** `PS_holdDecayRMSE` 19 %, `PS_rampupRMSE` 17 %**)**, then
`vall_t` 16 %. The re-stretch *transient shape* and passive mechanics are what
is left; force level, recovery rate and F-V are done.

### The trade that has to be reported

**The kʈr prediction got worse.** Force redevelopment now reaches plateau by
~1.625 s where the data is still climbing at 1.66 s, and the re-stretch peak
went from 1.02 (incumbent) to 1.17 against a measured 0.72. `ktr` was inside its
deadband throughout — it paid nothing, so the search was free to spend it, and
it did.

This is the documented `ktr ↔ FV` weld showing up again: both are set by the
same effective detachment rate, and nothing in the scored objective was
defending kʈr. It is the clearest argument in this campaign that **a deadband on
a solved term is not the same as a constraint**: `ktr` sat inside one
preparation-width and drifted straight out of it unpunished.

---

## 6. What I would do next, in order

1. **Score kʈr.** It is a measured protocol, it is currently a free variable,
   and it just moved the wrong way by ~10 %. Either add it to `fn` or accept
   that every future refit may spend it.
2. **`peak1_dSL` is now the single largest term** (31 %) and it is a *timing*
   feature the model has never fitted. Worth its own analysis rather than
   another refit.
3. **The passive block is 36 %** of what is left and has never had a campaign of
   its own. `PS_holdDecayRMSE` and `PS_rampupRMSE` both went the wrong way
   slightly in this fit.
4. **Fix `writeParamsToMFile` to `%.17g`** (one line, two branches). Every
   snapshot in `params/` is currently a 6-significant-digit quantisation of the
   point that was fitted.
5. **Wire in the retired-knob registry** (`conclusions.md` §9) so `sigma1` and
   friends stop being drawn.
6. **Re-run B4, C2 and a full-length B9** — all three were refused by the
   snapshot guard or under-budgeted, not by the data.
