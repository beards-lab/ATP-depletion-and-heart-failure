# Methods — correcting skinned-fibre force for rundown on effective activated time

*Appendix. Self-contained: notation, procedure, calibration, validation, limitations,
and a fully worked example. All figures are embedded below and live in `results/`;
all scripts are in this folder.*

---

## M.1 Rationale

A permeabilised cardiac trabecula loses force progressively over a session. If two
experimental conditions (here 2 mM and 8 mM MgATP) are measured sequentially on the
same preparation, the second is measured on a weaker fibre, and the raw ratio between
them confounds the treatment effect with that decay. Because the decay is monotonic
and the treatment effect is not, the confound does not average out: it biases every
sequential comparison in a direction set purely by the running order.

Two observations motivate the specific form of the correction developed here.

**(i) The fibre degrades while it generates force, not while it waits.** Between
activations the preparation sits in relaxing solution. The natural independent
variable is therefore accumulated *activated* time, not wall-clock time.

**(ii) The decay is linear in that variable, not exponential.** On 03/27 the same
8 mM protocol was run twice, 694 s apart. The two runs show within-run force slopes
of **−0.4592** and **−0.4618 kPa/s** — indistinguishable in *absolute* terms despite
the second run being 17.5 % weaker. A force-proportional (exponential) decay would
have predicted −0.379 kPa/s for the second run. Rundown is thus a constant
**kPa per activated second**, and the separation between two runs is properly
expressed as an *effective time*:

```
tau_eff  =  dF / |slope|  =  11.98 / 0.4592  =  26.1 s
```

against 694 s of wall clock. (This is the quantity `T0` computed by
`DataCuration/FitRundownCorrection.m`.)

---

## M.2 The bracket design

The correction requires at least one pair of measurements made under **identical
conditions at different points in the session**. On 03/27 the 8 mM protocol was
repeated at the end of the session, so the 2 mM run is *bracketed in time* by two
8 mM runs.

| run | session time | plateau force at SL 2.0 |
|---|---|---|
| 8 mM (fresh) | 0 s | 68.43 kPa |
| 2 mM | 132 s | 83.59 kPa |
| PNB + mavacamten | 348 s | 0.60 kPa |
| 8 mM (repeat) | 694 s | 56.45 kPa |

> **"Bracket" refers to TIME, not to force.** Force never rises: the 8 mM condition
> falls monotonically 68.43 → 56.45 kPa (−17.5 %). The 2 mM run reads higher than both
> only because low ATP increases force. The two 8 mM measurements are the only pair in
> the session recorded under identical conditions, so the trajectory between them is
> **pure rundown, uncontaminated by the treatment**, and everything else is measured
> against it.

*File-naming trap.* The 03/27 repeat is stored as `04_..._repeat` but was recorded
**last** (15:29:40, after the PNB+Mava run at 15:23:54). File numbering is not
chronological; session times must be taken from the `Created:` header of each `.txt`.

![The bracket, explained from the raw timecourses](results/bracket_explained.png)

**Panel (a)** the four raw force timecourses on the protocol clock — 2 mM (orange)
highest, 8 mM fresh (light blue), 8 mM repeat (dark blue) lowest, PNB+Mava (grey) near
zero. **(b)** the same four placed on a single session clock. **(c)** the damage
staircase (M.5). **(d)** the effect of the correction on this one preparation.

---

## M.3 Measured quantities

For each run *i*, from the merged `[t, L, F]` trace:

| symbol | definition | how measured |
|---|---|---|
| `F_i` | plateau force | mean over the isometric window at ML 1.0, `t ∈ [71.50, 72.35] s` (a second window at `[77.70, 78.30] s` is used for the slope) |
| `s_i` | within-run decay rate (kPa/s) | slope of a straight line fitted jointly to force in **both** isometric windows |
| `T_i` | activation duration (s) | interval over which smoothed force exceeds 20 % of its 99th percentile |
| `a_i` | force-generating? | false for blebbistatin / mavacamten runs |

Both plateau windows sit at the same commanded length (ML 1.0) and at the same
protocol times in every run, so within-run decay cancels when runs are ratioed.

Measured values in this dataset: `T_i = 15.3 ± 0.4 s` in every force-generating run;
the PNB+Mava runs have `|s| ≈ 0.016 kPa/s` and are assigned `T_i = 0`.

---

## M.4 The permanent fraction, phi

Naively integrating each run's own slope over its activation predicts a loss of

```
dF_pred  =  SUM_i |s_i| * T_i  =  26.52 kPa
```

whereas the bracket measures only **11.98 kPa**. Most of the within-run force decline
is therefore **reversible** — consistent with accumulation of hydrolysis products
during activation that washes out during the relaxed interval — and only a fraction is
permanent damage. Define

```
                dF_measured            11.98
    phi  =  ------------------   =  ---------  =  0.452
             SUM_i |s_i| * T_i         26.52
```

(0.342 if evaluated at SL 2.2 instead of SL 2.0). About **55 % of the within-run
decline recovers** between activations; ~45 % is permanent.

---

## M.5 The correction

**Damage accrued between two runs.** For runs separated by a set `A` of intervening
force-generating activations:

```
    dmg  =  phi * SUM_{i in A} |s_i| * T_i        [kPa]
```

For two consecutive activations this reduces to `dmg = phi * |s_e| * T_e`, where *e*
denotes the **earlier** run — damage accrues during *its* activation.

**Which run to correct.** The damage is done by the time the *later* run is recorded,
so the later run is the degraded one. Referencing the loss to that run's own force:

```
    f  =  dmg / F_later
```

and for the ratio of the low-ATP to the high-ATP condition:

```
    order 8 -> 2   (the 2 mM run is degraded):   R_corr  =  R_raw * (1 + f)

    order 2 -> 8   (the 8 mM run is degraded):   R_corr  =  R_raw / (1 + f)
```

Normalising by `F_later` rather than `F_earlier` matters: for a reversed-order
preparation the two differ substantially, and using the wrong one under-corrects.

**Interpretation.** Multiplying the degraded run upward is an estimate of what that
run *would have measured on an undamaged fibre*. It is **not** a claim that force
recovered.

**The damage staircase** — panel (c) of the figure in M.2 — is the visual form of this
model. The fibre ages only while activated, so the trajectory is flat in relaxing
solution and drops only during the shaded activation windows:

```
  68.43  --(8 mM activation, -3.18)-->  65.25  ----flat, 117 s----
         --(2 mM activation, -8.81)-->  56.44  ----flat, 547 s---->  repeat = 56.45
```

The ATP effect is the vertical gap at t = 132 s: the measured 2 mM force (83.59)
against where the 8 mM condition would have been at that instant (65.25).

---

## M.5b Refinement — reference each run to its own activation

The correction in M.5 compares two runs at a **fixed protocol window**. That is only
valid if both conditions are at the same point of their own activation, and they are
not.

**Force takes 3.3–5.8 s to develop** after activation begins, and rundown runs
throughout that rise, so the observed peak is not the undamaged force. The two ATP
conditions develop force at different rates:

| | time to peak after onset | state at W1 (68.7 s) |
|---|---|---|
| 8 mM | 5.3 – 5.8 s | still rising |
| 2 mM | 3.6 – 4.6 s | already falling |

At the W2 window the 2 mM run has therefore been decaying ~3.0 s against ~0.7 s for
the 8 mM run, **and** its slope is ~2× steeper (M.5c). On 03/27 that is 3.65 kPa lost
from the 2 mM peak versus 0.34 kPa from the 8 mM peak — the fixed-window comparison
understates the ATP effect by ~4 %, entirely *inside* the runs, where no between-run
correction can reach it.

**Procedure.** For each run, take the decay line from W2–W3 and evaluate it at that
run's own peak time and its own relaxation time:

```
    t_pk    = argmax of the ML=1.0 force envelope
    t_off   = last sample above 20 % of peak
    F_pk    = v2 + slope * (t_pk  - t_W2)      "force at full activation"
    F_end   = v2 + slope * (t_off - t_W2)      "force at relaxation"
    T_dec   = t_off - t_pk                     the run's own decay phase
```

Then form ratios from `F_pk` rather than from the window value, and compute damage
over `T_dec` rather than over a common activation duration:

| referencing | 03/27 | 04/03 | 04/10 | mean | CV |
|---|---|---|---|---|---|
| fixed window W2 | 1.225 | 1.196 | 1.678 | 1.366 | 20 % |
| own-peak | 1.269 | 1.233 | 1.649 | 1.384 | 17 % |
| own-peak + damage over `T_dec` | 1.306 | 1.313 | 1.486 | **1.369** | **7 %** |

The own-peak step costs nothing — no phi, no damage model — and should be used
regardless.

### M.5c Decay rate depends on ATP

| | 8 mM | 2 mM | ratio |
|---|---|---|---|
| absolute (kPa/s) | 0.560 ± 0.142 | 1.092 ± 0.113 | **×1.95** |
| relative (%/s) | 1.06 ± 0.37 | 1.50 ± 0.12 | ×1.42 |
| paired within preparation | | | ×2.66, ×1.37, ×1.77 |

All three preparations agree in direction. The 2 mM *relative* rate is tight across
preparations (CV 8 %) while the 8 mM is scattered (CV 35 %). Note that 03/27's 8 mM
run has the shallowest relative slope in the entire set (0.67 %/s): the preparation
the correction is calibrated on is the most robust one available.

### M.5d Does decay occur between activations?

Extrapolating each run to its own `F_pk` and `F_end` (above) makes the session a
chain. With `a` the ATP effect and `r` the between-activation recovery factor
(`r > 1` = force partially recovers while relaxed), the 03/27 chain gives two
equations in two unknowns:

```
    F_pk(2mM #2)  / F_end(8mM #1)  =  1.355  =  r * a
    F_pk(8mM #4)  / F_end(2mM #2)  =  0.798  =  r^n / a
```

| rest intervals assumed | `a` | `r` |
|---|---|---|
| one each (n = 1) | **1.303** | +4.0 % per rest |
| two after the 2 mM run (n = 2) | **1.320** | +2.6 % per rest |

**The "no decay between activations" assumption essentially holds**: the residual is
a small *recovery* of 2.6–4.0 % per rest interval, not further decay. This route
estimates the ATP effect as **×1.30–1.32 from one preparation, with no phi and no
damage model** — an independent confirmation of the value obtained in M.7.

![Decay geometry](results/decay_geometry.png)

*(a) activation envelopes normalised to their own peak — 2 mM (red) peaks earlier
than 8 mM (blue). (b) 03/27 decay lines extrapolated to each run's own peak (stars).
(c) paired decay rates. (d) the continuity chain. (e) CV by referencing scheme.
(f) decay-phase duration.*

---

## M.6 Worked example (03/27)

| step | quantity | value |
|---|---|---|
| damage during the fresh 8 mM run | `0.452 × 0.4592 × 15.3` | 3.18 kPa |
| damage during the 2 mM run | `0.452 × 1.1940 × 16.3` | 8.81 kPa |
| total predicted across the bracket | | 11.99 kPa |
| **measured** across the bracket | `68.43 − 56.45` | **11.98 kPa** |
| 8 mM interpolated to the 2 mM moment | `68.43 − 3.18` | 65.25 kPa |
| ATP ratio, raw | `83.59 / 68.43` | ×1.222 |
| **ATP ratio, corrected** | `83.59 / 65.25` | **×1.281** |

The 2 mM run decays roughly 2.6× faster within-run than the 8 mM run
(−1.194 vs −0.459 kPa/s), so it contributes most of the session's damage. A
preparation that runs 2 mM *first* therefore incurs a substantially larger correction
than one that runs 8 mM first — the reversed order is not merely a sign flip.

---

## M.7 Validation

The staircase terminating on the measured repeat (M.6) is **not** an independent
check: phi was calibrated from precisely that difference. Two genuine validations
exist.

**(1) Two unrelated calibrations of phi agree.** phi can be obtained either from the
bracket (mechanical, no reference to ATP) or by choosing the value that minimises the
between-preparation spread of the corrected ATP effect. These share no information:

| route | phi |
|---|---|
| 03/27 bracket | **0.452** |
| minimum cross-preparation CV | **0.625** |

Over that entire range the corrected ATP force effect is stable at ≈ ×1.35.

**(2) The correction reduces between-preparation scatter.** Across three
preparations, two of which ran 8→2 and one 2→8:

| phi | 03/27 | 04/03 | 04/10 | mean | CV |
|---|---|---|---|---|---|
| 0 (uncorrected) | 1.222 | 1.194 | 1.685 | 1.367 | **20.2 %** |
| 0.452 (bracket) | 1.281 | 1.325 | 1.441 | 1.349 | **6.1 %** |
| 0.625 (optimum) | 1.305 | 1.384 | 1.365 | 1.351 | **3.0 %** |

A correction derived from one preparation, applied blind to two others including one
of opposite running order, reduces scatter more than threefold. Applied at the feature
level across six force features, CV falls from 22 % to 11 %.

![Calibration and validation of the correction](results/rundown_correction.png)

**(a)** the linearity test of M.1. **(b)** integrated within-run slope versus the loss
the bracket actually measures — the basis of phi (M.4). **(c)** the two independent
calibrations of phi. **(d)** the three preparations converging as phi rises.
**(e)** raw versus corrected. **(f)** every run on one effective-time axis.

**Sensitivity to the damage model.** An alternative in which every activation does the
same damage per second (a universal rate `r = dF_meas / SUM T_i = 0.379 kPa/s`) rather
than a slope-proportional one gives 1.335 / 1.341 / 1.486 (CV 6.2 %). The result is
robust to this choice.

The downstream consequence, across all four preparations:

![The reconciled ATP effect](../ATPEffectReconciliation/results/atp_reconciliation.png)

---

## M.8 What must NOT be corrected

The same bracket shows the rate of force redevelopment falls by ×0.884 over the
session. It is nonetheless **wrong to apply a rundown correction to `ktr`**.

This was tested properly rather than asserted (`RunJointCorrection.m`). The coupling
between the two corrections was left free and scanned:

```
    ktr_frac(gap)  =  gamma * force_frac(gap)
```

| gamma | meaning | CV of the ktr ratio |
|---|---|---|
| 0 | no correction | **9.4 %** |
| +0.68 | implied by the bracket (force ×0.829, ktr ×0.884) | 15.8 % |
| +0.71 | predicted by the series-creep lesion (M.9) | 16.1 % |

Both physically motivated couplings make consistency **worse**. `CV(gamma)` is
monotone over the whole scanned range, so there is no interior optimum; the apparent
improvement at negative gamma is not a mechanism (it would mean rundown makes `ktr`
faster) and works only by lifting one preparation toward the other two — one free
parameter fitted to three points.

There is also no dose-response to follow: `ktr` ratio against force-damage dose gives
**r = −0.63, p = 0.57 (n = 3)**, non-monotone.

**The reason is mechanistic.** The lesion (M.9) reaches force through the
length–tension relation, a strong channel, and `ktr` through added series compliance,
a weak one — ×0.829 versus ×0.884 over the same bracket. A scaled correction assumes
precisely the proportionality the mechanism denies.

**A within-run `ktr` decay cannot currently be measured.** Comparing `ktr` from the
ktr burst (t ≈ 71.6 s) with slack segment 2 (t ≈ 75.5 s — same SL 2.0, same
activation, ~4 s apart) gives ratios of 1.09 / 0.86 / 1.45 / 1.22 / 0.94 / 0.73 across
the six runs: scatter far larger than any plausible 4-second decay. The two are
mechanically different perturbations, and the burst fit window spans only ~2.8 time
constants. A hi-res repeat run would be needed.

> **Recommendation:** correct force; report `ktr` raw; never compare `ktr` across long
> session gaps. Over the 132–308 s gaps actually used the implied `ktr` bias is ≈ 2 %,
> smaller than the between-preparation scatter.

![Joint correction](results/joint_correction.png)

> **Recommendation:** correct force; report `ktr` raw; never compare `ktr` across long
> session gaps.

Measuring the `ktr` change requires care: the 03/27 repeat is log-rate (10 ms) while
the time constant is ≈ 20 ms. The fresh trace must be resampled onto the repeat's own
sample times and both fitted identically, making the sampling bias common-mode.
Sampling alone accounts for ×0.956 of the apparent change.

![Length-tension and ktr evidence from the data](results/rundown_mechanism.png)

---

## M.9 Mechanistic basis, and why the correction is length-dependent

The correction above is empirical. Its form is nonetheless supported by an
identification of the underlying lesion (`RunMechanismSimulation.m`). Each candidate
lesion was imposed on the cross-bridge model, tuned to reproduce the observed force
loss exactly, so that `ktr` and the shape of the length–tension relation became free
discriminators:

| lesion | model lever | `ktr` at matched force | L–T shape |
|---|---|---|---|
| fewer attached heads | `kstiff1,2` | 0.970 | vertical scale |
| lost sarcomere length | `SL0` | 0.982 | **horizontal shift** |
| serial creep | `kSE` | strong `ktr` lever, negligible force lever | — |
| uniform kinetic slowdown | `xrate` | pure `ktr` lever, force unchanged | — |
| reduced attachment | `ka` | 1.045 (wrong sign) | vertical scale |
| **observed** | | **0.884** | **horizontal shift** |

No single lesion suffices. The `ktr` loss requires **added series compliance** — no
alternative reaches it. The force half is settled by a second, stronger test below.

![Model perturbation study](results/mechanism_simulation.png)

### M.9b Which force lever? The five-segment profile decides

The bracket yields not one number but a **profile**: `ktr` and amplitude at each of
the five slack segments, fresh and repeat (`RunModelRundownFit.m`). Observed
(repeat / fresh, sampling-matched):

| | s1 (2.04 µm) | s2 | s3 | s4 | s5 (1.88 µm) | mean |
|---|---|---|---|---|---|---|
| `ktr` ratio | 0.947 | 0.906 | 0.847 | 0.829 | 0.893 | 0.884 |
| amplitude ratio | 0.823 | 0.805 | 0.808 | 0.802 | 0.792 | **0.806** |

The amplitude profile is essentially **flat** (gradient −0.031 across the ladder):

| candidate | RMSE `ktr` | RMSE amplitude | total |
|---|---|---|---|
| SL0 only (−0.098 µm) | 0.116 | 0.058 | 0.130 |
| kSE only (×0.65) | 0.069 | 0.187 | 0.199 |
| SL0 + kSE | 0.066 | 0.052 | 0.084 |
| kstiff only (×0.84) | 0.108 | 0.023 | 0.111 |
| **kstiff + kSE (×0.84, ×0.65)** | **0.058** | **0.017** | **0.060** |

A length shift of the required size predicts a **steep** amplitude gradient
(0.910 → 0.726) because the deep segments lose more length–tension; the data show
−0.031, and a force scale predicts −0.001. **On five points the force loss is a
scale, not a shift.**

The earlier length-shift verdict came from the 6-point length–tension curve, in which
the **ML 1.10** pre-slack plateau (ratio 0.887) is the only point breaking an
otherwise flat 0.79–0.83 pattern — and it lies outside the slack ladder. Passive
force is largest there (~5 kPa of 80, and passive does not run down), which explains
part of it but not all. **The length-shift evidence rests on one measurement and
should be treated as provisional.**

**Revised mechanism:** rundown = loss of force-generating material in parallel (≈16 %)
**plus** added series compliance (≈35 %) — consistent with myofibrils rupturing near
the attachments, which removes parallel force-generating material and leaves a
compliant damaged region in series.

![Rundown as a decaying model parameter](results/model_rundown_fit.png)

### M.9c Rundown as a fitted nuisance parameter

The practical alternative to correcting data: give each run a **single rundown
coordinate** mapped jointly onto (`kstiff`, `kSE`), fit it alongside the physiology,
and the treatment parameters are estimated on a common footing automatically. This is
one nuisance dimension per run, not two free knobs — the two model parameters move
together because they are two consequences of one lesion. It also removes the need
for phi, `T_dec` and the own-peak referencing entirely, since the model then carries
the fibre state explicitly.

Two consequences follow for practice:

1. A **scalar** force correction is inadequate; the correction must be
   length-dependent. This is why `DataCuration/FitRundownCorrection.m` fits a surface
   `r(F, SL)` — it is absorbing a length shift in force coordinates. Its hand-tuned
   parameters (`r0 = 1.214`, `k = −0.6`) reproduce the bracket to within 2 %.
2. When representing rundown *inside* a model, prefer a per-run **length offset** to a
   per-run force scale. A corollary is that the nominal sarcomere-length axis is
   systematically too long for later activations (≈ 0.017 µm over a 132 s gap,
   ≈ 0.09 µm over a full session).

---

## M.10 Assumptions and limitations

1. **phi is calibrated on a single bracket** and assumed transferable between
   preparations. This is the principal limitation. A second end-of-session repeat on
   another preparation is the highest-value additional measurement.
2. **Damage is assumed proportional to each run's own within-run slope.** The
   universal-rate alternative gives a comparable answer (M.7), but the two are not
   formally distinguishable with one bracket.
3. **Blebbistatin / mavacamten runs are assigned zero damage.** Justified by their
   near-zero force and slope (`|s| ≈ 0.016 kPa/s`), but not independently verified.
4. **Linearity is established over one 17.5 % decrement.** Behaviour under larger
   decrements is untested.
5. **Recovery between activations is treated as complete and instantaneous.** The data
   cannot resolve a recovery time course; only the net permanent fraction is
   identified.
6. **The correction is applied as a fractional change to every force feature**, using
   the factor measured at the isometric plateau. Since the underlying lesion is a
   length shift (M.9), the true correction is mildly length-dependent; residual error
   from this is smaller than the between-preparation scatter that remains.
7. **One preparation (04/10) remains high** after correction (×1.54 vs ×1.24 / ×1.29).
   Its 8 mM run is the weakest 8 mM recording in the dataset. Either that run was
   additionally compromised, or a first activation at 2 mM is more damaging than the
   model allows. The reported mean is, if anything, an over-estimate; the two
   same-order preparations alone give ×1.27.

---

## M.11 Reproducing

| script | produces |
|---|---|
| `RunBracketExplainer.m` | `results/bracket_explained.png` — the bracket and the damage staircase |
| `RunDecayGeometry.m` | `results/decay_geometry.png`, `.mat` — activation rise vs decay, ATP-dependence of the decay rate, continuity/recovery test, own-peak referencing |
| `RunJointCorrection.m` | `results/joint_correction.png`, `.mat` — phi on the own-peak/`T_dec` footing; free scan of the force→`ktr` correction coupling |
| `RunModelRundownFit.m` | `results/model_rundown_fit.png`, `.mat` — rundown as a decaying model parameter, tested against the 5-segment `ktr` and amplitude profiles |
| `RunRundownCorrection.m` | `results/rundown_correction.png`, `.mat` — linearity test, phi calibration, phi scan, sensitivity |
| `RunRundownMechanisms.m` | `results/rundown_mechanism.png` — length–tension and `ktr` evidence from the data |
| `RunMechanismSimulation.m` | `results/mechanism_simulation.png` — model perturbation study |
| `../ATPEffectReconciliation/RunATPReconciliation.m` | the corrected ATP effect across all preparations |

Each runs standalone from a clean MATLAB session after
`cd(<repo>); addpath(genpath('.'))`.

---

## M.12 Summary of the procedure

1. Measure `F_i`, `s_i`, `T_i` for every run (M.3).
2. **Reference each run to its own activation**: find `t_pk`, extrapolate the decay
   line to it, and use `F_pk` in place of the window value; use each run's own decay
   phase `T_dec = t_off - t_pk` in place of a common activation duration (M.5b).
   This step is free and should always be applied.
3. Obtain phi from a same-condition bracket: `phi = dF_measured / SUM |s_i| T_i` (M.4).
4. Compute `dmg = phi * SUM |s_i| T_dec_i` over the intervening force-generating
   activations (M.5).
5. Form `f = dmg / F_later` and correct the **later** run upward (M.5).
6. Report corrected force ratios; report `ktr` **uncorrected** (M.8).
7. Validate by between-preparation scatter, not by the bracket the correction was
   calibrated on (M.7); cross-check against the continuity route (M.5d), which is
   independent of phi.
