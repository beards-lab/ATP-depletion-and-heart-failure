# Methods — correcting skinned-fibre force for rundown on effective activated time

*Appendix. Self-contained: notation, procedure, calibration, validation, limitations,
and a fully worked example. All figures referenced are in `results/`; all scripts are
in this folder.*

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
expressed as an *effective time*,

$$\tau_{\text{eff}} = \frac{\Delta F}{|\mathrm{slope}|} = \frac{11.98}{0.4592} = 26.1\ \text{s},$$

against 694 s of wall clock. (This is the quantity `T0` computed by
`DataCuration/FitRundownCorrection.m`.)

> **Figure `results/rundown_correction.png`, panel (a)** — measured slopes of the two
> 8 mM runs against the linear and exponential predictions.

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

*Note.* "Bracket" refers to **time**, not to force. Force never rises: the 8 mM
condition falls monotonically 68.43 → 56.45 kPa (−17.5 %). The 2 mM run reads higher
than both only because low ATP increases force. The two 8 mM measurements are the only
pair in the session recorded under identical conditions, so the trajectory between
them is **pure rundown, uncontaminated by the treatment**, and everything else is
measured against it.

*Note on file naming.* The 03/27 repeat is stored as `04_..._repeat` but was recorded
**last** (15:29:40, after the PNB+Mava run at 15:23:54). File numbering is not
chronological; session times must be taken from the `Created:` header of each `.txt`.

> **Figure `results/bracket_explained.png`, panels (a)–(b)** — the four raw force
> timecourses, and the same four placed on a single session clock.

---

## M.3 Measured quantities

For each run *i*, from the merged `[t, L, F]` trace:

| symbol | definition | how measured |
|---|---|---|
| $F_i$ | plateau force | mean over the isometric window at ML 1.0, `t ∈ [71.50, 72.35] s` (a second window at `[77.70, 78.30] s` is used for the slope) |
| $s_i$ | within-run decay rate (kPa s⁻¹) | slope of a straight line fitted jointly to force in **both** isometric windows |
| $T_i$ | activation duration (s) | interval over which smoothed force exceeds 20 % of its 99th percentile |
| $a_i$ | force-generating? | false for blebbistatin/mavacamten runs |

Both plateau windows sit at the same commanded length (ML 1.0) and at the same
protocol times in every run, so within-run decay cancels when runs are ratioed.

Measured values (this dataset): $T_i = 15.3 \pm 0.4$ s in every force-generating run;
the PNB+Mava runs have $|s| \approx 0.016$ kPa s⁻¹ and are assigned $T_i = 0$.

---

## M.4 The permanent fraction φ

Naively integrating each run's own slope over its activation predicts the loss

$$\Delta F_{\text{pred}} = \sum_i |s_i|\,T_i = 26.52\ \text{kPa},$$

whereas the bracket measures only **11.98 kPa**. Most of the within-run force decline
is therefore **reversible** — consistent with accumulation of hydrolysis products
during activation that washes out during the relaxed interval — and only a fraction is
permanent damage. Define

$$\boxed{\ \phi \;=\; \frac{\Delta F_{\text{measured}}}{\sum_i |s_i| T_i}\ }
\qquad \phi = \frac{11.98}{26.52} = \mathbf{0.452}$$

(0.342 if evaluated at SL 2.2 instead of SL 2.0). About **55 % of the within-run
decline recovers** between activations; ~45 % is permanent.

> **Figure `results/rundown_correction.png`, panel (b)** — integrated within-run slope
> versus the loss the bracket actually measures.

---

## M.5 The correction

**Damage accrued between two runs.** For runs separated by a set $\mathcal{A}$ of
intervening force-generating activations,

$$\boxed{\ \Delta_{\text{dmg}} \;=\; \phi \sum_{i\,\in\,\mathcal{A}} |s_i|\,T_i\ \ [\text{kPa}]\ }$$

For two consecutive activations this reduces to $\Delta_{\text{dmg}} = \phi\,|s_e|\,T_e$,
where *e* denotes the **earlier** run — damage accrues during *its* activation.

**Which run to correct.** The damage is done by the time the *later* run is recorded,
so the later run is the degraded one. Referencing the loss to that run's own force,

$$f \;=\; \frac{\Delta_{\text{dmg}}}{F_{\text{later}}},$$

and for the ratio of the low-ATP to the high-ATP condition,

$$\boxed{\ R_{\text{corr}} = \begin{cases}
R_{\text{raw}}\,(1+f), & \text{order } 8\!\rightarrow\!2 \ \ (\text{the 2 mM run is degraded}) \\[4pt]
R_{\text{raw}}\,/\,(1+f), & \text{order } 2\!\rightarrow\!8 \ \ (\text{the 8 mM run is degraded})
\end{cases}}$$

Normalising by $F_{\text{later}}$ rather than $F_{\text{earlier}}$ matters: for a
reversed-order preparation the two differ substantially, and using the wrong one
under-corrects.

**Interpretation.** Multiplying the degraded run upward is an estimate of what that
run *would have measured on an undamaged fibre*. It is not a claim that force recovered.

> **Figure `results/bracket_explained.png`, panel (c)** — the damage **staircase**:
> flat while the fibre sits in relaxing solution, dropping only during the shaded
> activation windows. It starts at the measured fresh 8 mM force, steps down by
> $\phi|s|T$ during each activation, and terminates on the measured repeat. The
> vertical bar at 132 s is the ATP effect: the gap between the measured 2 mM force
> and where the 8 mM condition would have been at that instant.

---

## M.6 Worked example (03/27)

| step | quantity | value |
|---|---|---|
| damage during the fresh 8 mM run | $\phi\,|s_1|T_1 = 0.452 \times 0.4592 \times 15.3$ | 3.18 kPa |
| damage during the 2 mM run | $\phi\,|s_2|T_2 = 0.452 \times 1.1940 \times 16.3$ | 8.81 kPa |
| total predicted across the bracket | | 11.99 kPa |
| **measured** across the bracket | $68.43 - 56.45$ | **11.98 kPa** |
| 8 mM interpolated to the 2 mM moment | $68.43 - 3.18$ | 65.25 kPa |
| ATP ratio, raw | $83.59 / 68.43$ | ×1.222 |
| **ATP ratio, corrected** | $83.59 / 65.25$ | **×1.281** |

The 2 mM run decays roughly 2.6× faster within-run than the 8 mM run
(−1.194 vs −0.459 kPa s⁻¹), so it contributes most of the session's damage. A
preparation that runs 2 mM *first* therefore incurs a substantially larger correction
than one that runs 8 mM first — the reversed order is not merely a sign flip.

---

## M.7 Validation

The staircase terminating on the measured repeat (M.6) is **not** an independent
check: φ was calibrated from precisely that difference. Two genuine validations exist.

**(1) Two unrelated calibrations of φ agree.** φ can be obtained either from the
bracket (mechanical, no reference to ATP) or by choosing the value that minimises the
between-preparation spread of the corrected ATP effect. These share no information:

| route | φ |
|---|---|
| 03/27 bracket | **0.452** |
| minimum cross-preparation CV | **0.625** |

Over that entire range the corrected ATP force effect is stable at ≈ ×1.35.

**(2) The correction reduces between-preparation scatter.** Across three
preparations, two of which ran 8→2 and one 2→8:

| φ | 03/27 | 04/03 | 04/10 | mean | CV |
|---|---|---|---|---|---|
| 0 (uncorrected) | 1.222 | 1.194 | 1.685 | 1.367 | **20.2 %** |
| 0.452 (bracket) | 1.281 | 1.325 | 1.441 | 1.349 | **6.1 %** |
| 0.625 (optimum) | 1.305 | 1.384 | 1.365 | 1.351 | **3.0 %** |

A correction derived from one preparation, applied blind to two others including one
of opposite running order, reduces scatter more than threefold. Applied at the
feature level across six force features, CV falls from 22 % to 11 %.

> **Figures** `results/rundown_correction.png` panels (c)–(f);
> `../ATPEffectReconciliation/results/atp_reconciliation.png`.

**Sensitivity to the damage model.** An alternative in which every activation does
the same damage per second (a universal rate $r = \Delta F_{\text{meas}}/\sum T_i
= 0.379$ kPa s⁻¹) rather than a slope-proportional one gives 1.335 / 1.341 / 1.486
(CV 6.2 %). The result is robust to this choice.

---

## M.8 What must **not** be corrected

The same bracket shows the rate of force redevelopment falls by ×0.884 over the
session. It is nonetheless **wrong to apply a rundown correction to $k_{tr}$**: doing
so *increases* between-preparation scatter (CV 9 % → 13–16 %), under both wall-clock
and effective-time models. The bracket's $k_{tr}$ loss does not transfer to the much
shorter 132–308 s gaps. Over those gaps the implied bias is ≈ 2 %.

**Recommendation:** correct force; report $k_{tr}$ raw; never compare $k_{tr}$ across
long session gaps.

*(Measuring the $k_{tr}$ change requires care: the 03/27 repeat is log-rate, 10 ms,
while τ ≈ 20 ms. The fresh trace must be resampled onto the repeat's own sample times
and both fitted identically, making the sampling bias common-mode. Sampling alone
accounts for ×0.956 of the apparent change.)*

---

## M.9 Mechanistic basis, and why the correction is length-dependent

The correction above is empirical. Its form is nonetheless supported by an
identification of the underlying lesion (`RunMechanismSimulation.m`). Each candidate
lesion was imposed on the cross-bridge model, tuned to reproduce the observed force
loss exactly, so that $k_{tr}$ and the shape of the length–tension relation became
free discriminators:

| lesion | $k_{tr}$ at matched force | L–T shape |
|---|---|---|
| fewer attached heads (`kstiff`) | 0.970 | vertical scale |
| lost sarcomere length (`SL0`) | 0.982 | **horizontal shift** |
| serial creep (`kSE`) | strong $k_{tr}$ lever, negligible force lever | — |
| uniform kinetic slowdown (`xrate`) | pure $k_{tr}$ lever, force unchanged | — |
| reduced attachment (`ka`) | 1.045 (wrong sign) | vertical scale |
| **observed** | **0.884** | **horizontal shift** |

No single lesion suffices. Force loss together with the observed *bend* in the
length–tension relation requires a length shift; the $k_{tr}$ loss requires added
series compliance. **These are the same physical lesion** — a series elastic element
that creeps longer and softer simultaneously. Imposing both (`kSE` ×0.65, SL −0.098 µm)
reproduces force ×0.827, $k_{tr}$ ×0.877 and the bend, against observed 0.829 / 0.884 /
bend. An alternative pairing reaching the same $(F, k_{tr})$ point by fewer heads plus
slower motors predicts a pure vertical scale and is rejected by the shape.

Two consequences follow for practice:

1. A **scalar** force correction is inadequate; the correction must be
   length-dependent. This is why `DataCuration/FitRundownCorrection.m` fits a surface
   $r(F, \mathrm{SL})$ — it is absorbing a length shift in force coordinates. Its
   hand-tuned parameters ($r_0 = 1.214$, $k = -0.6$) reproduce the bracket to within 2 %.
2. When representing rundown *inside* a model, prefer a per-run **length offset** to a
   per-run force scale. A corollary is that the nominal sarcomere-length axis is
   systematically too long for later activations (≈ 0.017 µm over a 132 s gap,
   ≈ 0.09 µm over a full session).

> **Figure `results/mechanism_simulation.png`.**

---

## M.10 Assumptions and limitations

1. **φ is calibrated on a single bracket** and assumed transferable between
   preparations. This is the principal limitation. A second end-of-session repeat on
   another preparation is the highest-value additional measurement.
2. **Damage is assumed proportional to each run's own within-run slope.** The
   universal-rate alternative gives a comparable answer (M.7), but they are not
   formally distinguishable with one bracket.
3. **Blebbistatin/mavacamten runs are assigned zero damage.** Justified by their
   near-zero force and slope (|s| ≈ 0.016 kPa s⁻¹), but not independently verified.
4. **Linearity is established over one 17.5 % decrement.** Behaviour under larger
   decrements is untested.
5. **Recovery between activations is treated as complete and instantaneous.** The
   data cannot resolve a recovery time course; only the net permanent fraction is
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
| `RunRundownCorrection.m` | `results/rundown_correction.png`, `.mat` — linearity test, φ calibration, φ scan, sensitivity |
| `RunRundownMechanisms.m` | `results/rundown_mechanism.png` — length–tension and $k_{tr}$ evidence from the data |
| `RunMechanismSimulation.m` | `results/mechanism_simulation.png` — model perturbation study |
| `../ATPEffectReconciliation/RunATPReconciliation.m` | the corrected ATP effect across all preparations |

Each runs standalone from a clean MATLAB session after
`cd(<repo>); addpath(genpath('.'))`.

---

## M.12 Summary of the procedure

1. Measure $F_i$, $s_i$, $T_i$ for every run (M.3).
2. Obtain φ from a same-condition bracket: $\phi = \Delta F_{\text{meas}} / \sum |s_i| T_i$ (M.4).
3. Compute $\Delta_{\text{dmg}} = \phi \sum |s_i| T_i$ over the intervening
   force-generating activations (M.5).
4. Form $f = \Delta_{\text{dmg}} / F_{\text{later}}$ and correct the **later** run
   upward (M.5).
5. Report corrected force ratios; report $k_{tr}$ **uncorrected** (M.8).
6. Validate by between-preparation scatter, not by the bracket the correction was
   calibrated on (M.7).
