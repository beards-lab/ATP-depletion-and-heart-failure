# ATP gates: power law → saturable Hill/Michaelis

**Question.** The low-ATP fit expresses every ATP effect as a power law,
`rate_eff = rate·(ATP/ATP_ref)^p`. Can it use the classic saturable form instead —
`rate_eff = rate·Corr_ref·ATP/(K + ATP)`, with `Corr_ref` chosen so 8 mM is untouched —
what do the coefficients become, and does the model still hold at both 8 and 2 mM?

**Verdict.** Yes, exactly, and with no refit. All 9 active gates convert in closed
form to a classic Michaelis term (Hill coefficient ±1). 8 mM is **bit-identical**
(full battery, `Et = 6.50884430` under every configuration) and 2 mM reproduces the
power-law fit to within the model's own numerical noise floor. The two forms are
**indistinguishable on the data by construction** — but they are not the same model,
and the saturable one is the better object: it is bounded, and it converts
uninterpretable exponents into per-step Km values that can be checked against
chemistry. Doing that check finds one gate that does not survive it — `kmsrd`,
which carries the largest ATP effect in the fit.

**Read §3a even if you skip the rest.** Verifying the 2 mM equivalence turned up
something with wider consequences than the re-parameterisation itself: the 2 mM
objective carries **±0.3–0.6 of numerical scatter** at the production ODE tolerance.
Writing the incumbent's exponents to a snapshot at `%.6g` — a 3.8×10⁻⁷ parameter
change — moves its cost by **0.60**. The model is perfectly reproducible; the
*objective* is rough. Snapshot headers are therefore not comparable to fresh
evaluations, and optimiser accept-thresholds of `1e-4` are fitting the solver.

---

## 1. Why it converts exactly, and why no refit is needed

The data live at exactly **two** [ATP], 8 and 2 mM. Both forms are pinned to 1 at
8 mM by construction, so *a gate is completely described by one number*: the factor
it applies at 2 mM,

    g₂ = (2/8)^p = 0.25^p

`p` is not identified by the data — `g₂` is. `p` is one parameterisation of it. Write
the renormalised Hill gate with coefficient `n` and half-saturating concentration `K`
(`A = 2`, `A_ref = 8`, `r = A/A_ref`):

    f(A) = [Aⁿ/(Kⁿ + Aⁿ)]·[(Kⁿ + A_refⁿ)/A_refⁿ] = rⁿ·(Kⁿ + A_refⁿ)/(Kⁿ + Aⁿ)

Two limits matter, and they are what make the whole thing work:

| limit | meaning | value at 2 mM |
|---|---|---|
| `Kⁿ → 0` | step is **saturated** at both concentrations | `f = 1` — no ATP effect at all |
| `Kⁿ → ∞` | step is far from saturation | `f = rⁿ` — **the power law** |

So **the power law is the unsaturated limit of the Hill gate**, and every value in
between is reachable. Setting `κ = Kⁿ` and solving `f = g₂`:

    κ = Aⁿ(1 − g₂)/(g₂ − rⁿ),        K = κ^(1/n)

with a positive solution iff `g₂` lies strictly between 1 and `rⁿ` — i.e. `n > p` for
a suppressive gate, `n < p` (both negative) for an enhancing one. `n = ±1`, the classic
Michaelis case the question asks for, therefore works for **any gate with 0 < p < 1**,
which is all of them here. `solveHillK.m` is that formula; `atpGateFactor.m` evaluates
`f` in log space.

Two consequences worth stating plainly:

- **No refit is possible or needed.** The conversion is a bijection on the only
  quantity the data constrain. Re-optimising under the Hill form searches the same
  cost surface in different coordinates and must land on the same optimum.
- **The Michaelis form is not a constraint here.** With `n = ±1` the reachable range
  at 2 mM is `g₂ ∈ (0.25, 4)`. The strongest incumbent gate is `kmsrd` at 0.414 —
  comfortably interior, so nothing is being clipped. (A gate needing suppression
  beyond ×0.25 or enhancement beyond ×4 *would* need `|n| > 1`; none does.)

## 2. The recalculated coefficients

From `params/lowATP_powergate2_opt.m` (gates with `|p| < 10⁻³` are optimiser dust and
are dropped — they perturb their rate by < 10⁻⁴ relative):

| rate | tier | `p` (power) | `g₂` at 2 mM | `n` | **`K` (mM)** | rate `k` (s⁻¹) | implied `k/K` (M⁻¹s⁻¹) |
|---|---|---|---|---|---|---|---|
| `k2` strong-bridge detachment | A | +0.3167 | 0.6446 | +1 | **1.801** | 100.5 | 5.6×10⁴ |
| `kmsrd` SRX-ADP return | A | +0.6362 | 0.4140 | +1 | **7.148** | 19.8 | 2.8×10³ |
| `kah` hydrolysis / re-priming | A | +0.0997 | 0.8709 | +1 | **0.4159** | 101.6 | 2.4×10⁵ |
| `ksr0` SRX entry | A | +0.1000 | 0.8705 | +1 | **0.4174** | 328.5 | 7.9×10⁵ |
| `ksr2srd` SRX → SRX-ADP | A | +0.0519 | 0.9306 | +1 | **0.2040** | 110.4 | 5.4×10⁵ |
| `ksrd2sr` SRX-ADP → SRX | A | +0.1000 | 0.8705 | +1 | **0.4174** | 8.68 | 2.1×10⁴ |
| `kamh` reverse priming | A | −0.0500 | 1.0718 | −1 | **81.56** | 8.20 | ATP-*inhibited* |
| `kSE` series elastic | B | +0.0700 | 0.9075 | +1 | **0.2814** | — | (not a rate) |
| `mu_neg` negative-strain viscosity | B | +0.1500 | 0.8122 | +1 | **0.6680** | — | (not a rate) |

Written to **`params/lowATP_hillgate_opt.m`** (`AtpGateForm = 'hill'`, `pAtp_* = ±1`,
`KAtp_*` as above). Converted factors match the power-law factors to **3.4×10⁻¹⁶**.

**Independent corroboration.** The model already contained two hand-built Michaelis
gates from earlier work, with Km chosen on physiological grounds and never fitted:
`K_T2 = 2.5 mM` for the same detachment step and `K_srx = 9.3 mM` for the same SRX-ADP
return. The conversion arrives at **1.80** and **7.15 mM** from a completely different
route — a power-law fit that knew nothing about them. Agreement within ~30 % on both.

## 3. Consistency at 8 mM and 2 mM

**8 mM — bit-identical, full battery** (`RunFullBatteryCheck.m`, force-velocity +
slack + passive, scored on the complete `fn` list):

| configuration | `Et` |
|---|---|
| (a) `rskR2_w025_opt`, no ATP gates | 6.50884430 |
| (b) + power-law gates @ 8 mM | 6.50884430 |
| (c) + saturable Hill gates @ 8 mM | 6.50884430 |

Zero deviation in every feature. This is structural, not lucky: `evaluateModel`
short-circuits the whole gate block when `MgATP == MgATP_ref`, and both forms equal 1
there for any coefficients anyway.

**2 mM — equal to the noise floor, and the noise floor is the interesting part.**

| comparison | 2 mM feature cost | Δ |
|---|---|---|
| power law (incumbent) | 16.586911 | — |
| power law, re-run identically | 16.586911 | **0** |
| Hill, all gates converted (rates equal to 3×10⁻¹⁶) | 16.570962 | −0.0159 |
| Hill, `|p| < 10⁻³` gates dropped (rates equal to 8×10⁻⁵) | 16.519905 | −0.0670 |

The pipeline is **bit-reproducible** (row 2), so the residual is not run-to-run
variability. It is not the reparameterisation either — the rate constants going into
row 3 are identical to the last bit. What it is: the objective is a **rough function**
of the parameters at a scale ~10⁻⁹, and any change to the arithmetic re-draws from that
roughness. The cost change is **non-monotone in the size of the perturbation** —
enlarging the input by 11 orders of magnitude (row 4 vs row 3) changes the cost *less*
than a 4×10⁻⁷ nudge does (§3a). Sensitivity does not behave like that; scatter does.

**The ODE tolerance is a contributor but NOT the main source** — measured, not assumed
(`RunCostNoiseAmplitude.m`, 11-point delta ladder at each tolerance):

| RelTol | peak-to-peak | σ | mean | wall time |
|---|---|---|---|---|
| 1e-2 (production) | 0.2495 | **0.0761** | 16.5186 | 45 s/pair |
| 1e-3 | 0.1987 | **0.0574** | 16.4548 | 61 s/pair |
| 1e-4 | — | — | — | **killed at 1200 s** |

A full decade of tolerance buys only **×1.3 in σ**, for ~1.4× runtime, and the next
decade is unaffordable. So **tightening the solver is not the fix.** The dominant source
has to be something tolerance-insensitive, i.e. **discrete branch decisions in the
solve/extract chain** — event-triggered strain-grid shifts (`UseSpaceExtension` moves
the grid by an integer `nds`), force-onset event crossings, slack-onset detection, and
the windows the `rsK`/`ktr`/`A` fits are run over. One of those flipping side changes a
feature discontinuously. That is corroborated by where the wobble lands: amplitudes are
stable to 0.002–0.03 % while exactly the *fitted-rate* features move 0.1–1.7 %.

> **Retracted:** an earlier reading of this section claimed the difference "collapses
> ×24 when the tolerance is tightened", from the single power-vs-hill pair
> (0.067 → 0.0028 at 1e-3). The ladder shows that pair was two draws landing close.
> The honest number is ×1.3.

Where it lands, in absolute units (identical model, both forms):

| feature | relative wobble | verdict |
|---|---|---|
| `steady`, `peak1_y`, `vall_y`, `peak2` (amplitudes) | 0.002 – 0.03 % | solid |
| `A`, `ktr`, `rsK` (rates fitted to a transient) | 0.1 – 1.7 % | noise-limited |
| `ovrsht_dy` | one element is 6×10⁻⁵ vs exactly 0 | already retired by [SlackDataAnalysis](../SlackDataAnalysis/conclusions.md) |

### 3a. How big the floor really is — the incumbent snapshot measures it for us

The incumbent's header records `2 mM absolute cost 15.9849` (2026-08-10), but the
committed snapshot scores **16.587** today. `RunIncumbentAudit.m` splits that gap, and
the split is the most useful number in this document:

| | 2 mM cost |
|---|---|
| optimiser's recorded `best.c` | 15.984853 |
| `best.x` re-evaluated at full precision, today | **15.984853** |
| the committed snapshot | **16.586911** |
| | |
| drift since the fit (code + data) | **+0.0000** |
| `%.6g` serialisation of the exponents | **+0.6021** |

The model is **perfectly reproducible over two weeks** — zero drift. The entire 0.60
gap comes from writing the exponents to the snapshot at 6 significant figures. The
largest parameter difference that causes it is **3.79×10⁻⁷**.

So a perturbation ~10⁹ times smaller than any physically meaningful change moves the
objective by **0.60 — 3.8 % of it**. `RunCostNoiseAmplitude.m` sizes this properly by
nudging **one** exponent (`pAtp_kmsrd = 0.636193`) along a ladder of deltas far too
small to mean anything, at the production tolerance:

| δ | 0 | +1e−9 | −1e−9 | +1e−8 | −1e−8 | +1e−7 | −1e−7 | +1e−6 | −1e−6 | +1e−5 | −1e−5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| cost | 16.587 | 16.520 | 16.493 | 16.534 | 16.458 | **16.663** | **16.413** | 16.579 | 16.481 | 16.417 | 16.559 |

**Peak-to-peak 0.250, σ = 0.076** over 11 samples — and the values have no relationship
to either the magnitude or the sign of δ. A **10⁻⁹** nudge moves the cost by 0.067; the
largest excursions come from the middle of the ladder. This is not a cost surface with
a small floor, it is a smooth surface plus **numerical scatter of σ ≈ 0.08 per
perturbed parameter**. The snapshot's 0.60 is the same effect with all 21 exponents
rounded at once — independent contributions adding to ≈ √21 × 0.076 ≈ 0.35 σ, which
puts 0.60 at a routine ~1.7 σ.

Three consequences, in increasing order of importance:

1. **The Hill-vs-power comparison is unaffected.** 16.520 vs 16.587 is a 0.067
   difference against σ ≈ 0.076 and a peak-to-peak of 0.25 — an unremarkable draw.
   (The ladder's `δ = +1e−9` entry is 16.5199; the Hill conversion's is 16.5199 to four
   decimals. Same distribution, different draw.) The equivalence rests on the
   **algebra** — rate factors identical to 3×10⁻¹⁶ (§1, §2) and 8 mM bit-identical —
   not on the 2 mM costs matching, which they cannot do better than the scatter allows.
2. **Never compare against a snapshot header.** A reader seeing 15.98 in the header
   and 16.52 from the Hill conversion would record a 0.54 "regression" that is entirely
   `%.6g`. Snapshots should serialise at `%.17g`, or record the cost *of the snapshot
   as written* rather than of the optimiser's internal vector.
3. **`RunPowerGateOptim2.m` accepts a coordinate move when it improves the cost by
   `1e-4`** — roughly **three orders of magnitude below σ**. Much of the late
   fine-tuning in that run (and in any driver using a comparable threshold) was
   fitting the scatter, not the data. The gate exponents sitting at the 10⁻⁶ level in
   the incumbent snapshot are exactly that: noise the optimiser chased and recorded.
   This does not overturn the fit's headline improvement (18.58 → 15.98 is 3.4× the
   peak-to-peak), but the last few tenths of it are not real.

**What to do about it**, in order of cost:

- **Raise accept thresholds to ~0.25** (the measured peak-to-peak) in any driver that
  fine-tunes against the 2 mM cost. Free, and it stops the overfitting immediately.
- **Serialise snapshots at `%.17g`** in `writeParamsToMFile`/the driver footers, so a
  saved optimum *is* the optimum that was found.
- **Do not reach for the tolerance** — ×1.3 for 1.4× runtime, and 1e-4 does not finish.
- **The real fix, if it is ever worth it**, is to remove the discrete decisions from
  feature extraction: anchor the `rsK`/`ktr`/`A` fit windows to protocol time rather
  than to a detected crossing. That is a separate piece of work and this analysis does
  not attempt it — it only establishes that the roughness is there and roughly how big.

`params.odeRelTol` / `params.odeAbsTol` are now overridable (defaults unchanged) so
this stays measurable. **All comparisons in this document are like-for-like**, measured
within one session against a freshly re-evaluated incumbent.

## 4. What the change actually buys

The two forms cannot be told apart by any fit to 8 and 2 mM data. So the reason to
prefer one is not fit quality — it is these three:

**(a) Boundedness, and a different depletion prediction.** The power law is scale-free:
as ATP → 0 every gated rate → 0 (or → ∞ for the enhancing gate), and it keeps changing
without limit above 8 mM. The Hill form saturates above and goes *linear* in [ATP]
below Km — which is the physically correct behaviour, and it is not a small difference
where this project actually points. Gate factor, power / hill (`RunAtpDoseResponse.m`):

| gate | 4 mM | **2 mM** | 1 mM | 0.5 mM | 0.25 mM |
|---|---|---|---|---|---|
| `k2` | 0.803 / 0.845 | **0.645 / 0.645** | 0.518 / 0.437 | 0.416 / 0.266 | 0.334 / **0.149** |
| `kmsrd` | 0.643 / 0.679 | **0.414 / 0.414** | 0.266 / 0.232 | 0.171 / 0.124 | 0.110 / **0.064** |
| `kah`, `ksr0`, `ksrd2sr` | 0.933 / 0.953 | **0.871 / 0.871** | 0.813 / 0.742 | 0.758 / 0.573 | 0.708 / **0.394** |
| `kamh` (enhancing) | 1.035 / 1.047 | **1.072 / 1.072** | 1.110 / 1.085 | 1.149 / 1.091 | 1.189 / **1.095** |

Welded at 2 mM, crossing just above it, and **~2× apart by 0.25 mM** — the saturable
form predicts a far sharper collapse of detachment and SRX return under depletion,
while the enhancing gate stays bounded instead of growing without limit.

That propagates to the model's observables (`RunAtpDoseResponse.m`, each form
normalised to its own 8 mM value, ratio hill/power in bold):

| [ATP] mM | 0.25 | 0.5 | 1 | **2** | 4 | **8** | 12 | 20 |
|---|---|---|---|---|---|---|---|---|
| `ktr`, power | 0.118 | 0.195 | 0.326 | 0.490 | 0.709 | 1.000 | 1.191 | 1.485 |
| `ktr`, hill | 0.068 | 0.106 | 0.265 | 0.489 | 0.750 | 1.000 | 1.110 | 1.215 |
| **hill/power** | **0.58** | **0.54** | 0.81 | **1.00** | 1.06 | **1.00** | 0.93 | **0.82** |
| `rsK`, hill/power | **0.60** | 0.66 | 0.79 | **1.00** | 1.04 | **1.00** | 0.93 | **0.76** |
| `steady`, hill/power | 1.05 | 1.08 | 1.04 | **1.00** | 0.99 | **1.00** | 1.02 | 1.07 |
| `A`, hill/power | 0.95 | 1.23 | 1.07 | **1.00** | 0.98 | **1.00** | 1.03 | 1.08 |

Three things to take from it:

- **The weld holds at the level of full model observables**, not just gate factors:
  every ratio is 1.00 at 2 and 8 mM. That is the strongest form of the consistency
  claim — two different gate laws, one model.
- **Kinetics are form-sensitive, amplitudes are not.** By 0.25 mM the two forms differ
  by **40 % on `ktr` and `rsK`** but only ~5 % on force. Sensible: force is set by where
  the gated rates balance, kinetics by the rates themselves.
- **The power law's unboundedness bites above 8 mM too** — at 20 mM it predicts `ktr`
  ×1.49 and `rsK` ×1.43 and still climbing, versus ×1.22/×1.09 saturating. 20 mM MgATP
  is an achievable bath condition, so this is a **cheap discriminating experiment** —
  and it does not even require going into depletion.

Two parameterisations no existing measurement can separate thus make materially
different predictions exactly in the regime the project is named after. The conclusion
is to **measure a third [ATP]**, not to prefer one fit.

*(Caveat: the 0.25–0.5 mM column is far outside the fitted range and the model is
marginal there — `A` and `peak1_y` are non-monotone between 0.5 and 0.25 mM under the
Hill form. Read the ×0.4 kinetic divergence as the robust signal, not those points.)*

**(b) Falsifiability — the big one.** `p` cannot be checked against anything. `K` can.
If a gated step is ATP-binding-limited, its second-order association rate is `k/K`,
and the textbook band for ATP-induced actomyosin dissociation is ~10⁵–10⁶ M⁻¹s⁻¹
(lower in a fibre lattice than in solution). From the table in §2:

- `kah`, `ksr0`, `ksr2srd` → **2.4–7.9×10⁵ M⁻¹s⁻¹**. Squarely in the band. These gates
  are consistent with real ATP binding at the step they are attached to.
- `k2` → 5.6×10⁴ and `ksrd2sr` → 2.1×10⁴ M⁻¹s⁻¹. One to two orders of magnitude slow;
  defensible for a strained, lattice-constrained head, but on the edge.
- **`kmsrd` → 2.8×10³ M⁻¹s⁻¹. Two to three orders of magnitude too slow to be ATP
  binding.** Read as a rapid-equilibrium dissociation constant instead, `K = 7.1 mM`
  says the SRX-ADP return is only half-saturated at 7 mM ATP — i.e. barely saturated at
  physiological [ATP], which no measured nucleotide affinity supports.

  `kmsrd` is the **ktr lever** and carries the **largest ATP effect in the fit**
  (×0.414). So the model's dominant low-ATP mechanism is the one that fails the
  chemistry check. The power law could never have said this.

**(c) A testable structural hypothesis falls out.** Pin Km to something physiological
and ask what Hill coefficient each gate then needs (§2 of `RunHillGateConversion`
output; at `K = 2.5 mM`):

| gate | required `n` |
|---|---|
| `kmsrd` | **2.15** |
| `k2` | 0.77 |
| every other gate | 0.10 – 0.33 |

At a plausible Km, every gate but one needs a *shallow* (sub-unity) dependence — the
signature of weak phenomenological corrections — while the SRX-ADP return demands
`n ≈ 2`. That is not an arbitrary number: SRX is the **two-headed** interacting-heads
motif, so an exit step requiring **two** ATP molecules is mechanistically natural.
Worth testing directly: `pAtp_kmsrd = 2`, `KAtp_kmsrd ≈ 2.5` is a *physiologically
constrained* version of the same fit, and it is the first ATP-gate hypothesis in this
project that predicts a specific curve shape rather than absorbing one.
It is **not** what the current snapshot does — the snapshot is the faithful `n = 1`
conversion of the incumbent fit.

Note the flip side, from the same table: with `K = 1 mM` most gates have **no
solution at any `n`**. A step whose Km is below both measured concentrations is
saturated at both and cannot carry an ATP effect at all. So the fit is not merely
compatible with mM-range Km values — it *requires* them.

## 5. Caveats

1. **Km is not identifiable from these data, and this analysis does not claim it is.**
   Two concentrations fix one number per gate. Every `(n, K)` pair on the curve in §1
   fits identically. The values in §2 are "what Km must be **if** the Hill coefficient
   is 1" — a choice, made explicit. Identifying Km needs a third [ATP]; a 4–5 point
   ATP series over 0.5–8 mM would separate the forms outright (§4a).
2. Tier-B gates (`kSE`, `mu_neg`) are stiffness/viscosity, not rates. Their "Km" is a
   phenomenological interpolation coefficient — the `k/K` reading does not apply, and
   the weaker justification inherited from [LowATP_AbsoluteFit](../LowATP_AbsoluteFit/conclusions.md) stands.
3. `kamh`'s inhibitory branch (`n = −1`, `K = 81.6 mM`) is a 7 % effect on a reverse
   step, fitted at the edge of the optimiser's resolution. Do not read physiology into it.
4. The residual low-ATP problem is untouched by any of this. `rsK` is still the
   dominant term and is still ~×3 too fast; re-parameterising a gate cannot fix a
   defect that is shared with the 8 mM condition.
5. `KAtp_*` is not in `parameterBounds`/`assertActiveParams` (neither were `pAtp_*`,
   `K_T2` or `K_srx`). If a future fit frees Km, add bounds there first.
6. **Artifact state:** `results/tolerancefloor.mat` holds only the 1e-2 and 1e-3 rows —
   the original run was killed by the unaffordable 1e-4 leg, and `RunToleranceFloor.m`
   was then hardened to record such a row instead of aborting but was **not re-run**
   (it would cost the same 20-minute timeout for nothing new). The tolerance numbers
   quoted in §3 come from `results/costnoise.mat`, which is the better measurement
   anyway — a spread over 11 samples rather than a single pair.

## 6. Files

| file | what it does |
|---|---|
| `RunHillGateConversion.m` | converts `p → (n, K)`, writes `params/lowATP_hillgate_opt.m`, verifies 2 mM equality + 8 mM invariance. `SKIP_SIM = true` for the arithmetic only |
| `RunFullBatteryCheck.m` | the full 8 mM battery under no gates / power gates / Hill gates |
| `RunEquivalenceDiagnostic.m` | separates "conversion is inexact" from "pipeline is noisy"; runs the same parameters twice as the control |
| `RunEquivalenceBreakdown.m` | post-processes the above per feature (no simulation) |
| `RunToleranceFloor.m` | sweeps `odeRelTol` on the power/hill pair; records 1e-4 as unaffordable rather than aborting |
| `RunIncumbentAudit.m` | reconciles the incumbent's recorded 2 mM cost with its value today |
| `RunCostNoiseAmplitude.m` | measures the objective's numerical scatter by nudging one exponent |
| `RunAtpDoseResponse.m` | where the two forms differ: gate factors and model observables over 0.25–20 mM |
| `solveHillK.m`, `atpGateFactor.m` | the conversion formula and a numerically safe gate evaluator |

Model changes: `params0.AtpGateForm` (`'power'` default | `'hill'`) and `params0.KAtp_<rate>`
(default `Inf` = the power law) in `getParams`; the gate block in `evaluateModel`;
`params0.odeRelTol` / `odeAbsTol` made overridable at their existing values.
Every previous snapshot is byte-identical under the defaults.

    matlab -batch "run('Analyses/AtpGateHillForm/RunHillGateConversion.m')"
