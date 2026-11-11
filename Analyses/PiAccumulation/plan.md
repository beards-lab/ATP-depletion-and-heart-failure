# Pi accumulation in the permeabilised trabecula — analysis design (2026-10-07)

*Design only; nothing here has been run. Feeds manuscript §3.1 (decision D2 in
`Docs/manuscript/REVIEW-corrections.md`).*

## Question

During each activation, how much inorganic phosphate builds up inside the trabecula, at
8 mM and at 2 mM bath MgATP? And how much of the measured 2 mM/8 mM effect could that
build-up explain — force ×1.32, `ktr` ×0.57 (slack) / ×0.74 (protocol), FV2 ×0.715?

## Why it matters

- **The preparation.** It is a permeabilised mouse RV trabecula, roughly 140 × 100 µm in
  cross-section (Beard 2022: 137 ± 15 × 97 ± 7 µm). Every metabolite reaches the
  myofilaments only by diffusion from the bath.
- **What the buffer clamps.** CrP/CK (12 mM CrP, 10 IU/mL CK) regenerates ATP from ADP
  locally, so it clamps ADP. **Nothing clamps Pi**: Pi has to diffuse out.
- **What Pi does.** Pi depresses isometric force roughly linearly in log[Pi], with little
  effect on Vmax (Cooke & Pate 1985; Ebus & Stienen 2001 for cardiac).
- **Why the confound runs in our direction.** Pi production equals ATPase, and ATPase is
  lower at 2 mM. The low-ATP condition therefore carries *less* Pi. Part of its higher force
  could be **Pi relief**, not an ATP effect. The confound has the same sign as the effect we
  report.
- **A second confound.** The same geometry could deplete ATP / CrP in the core at a 2 mM bath.
  Ogut & Brozovich 2003: lowering CrP at constant MgATP halves `ktr` in trabeculae, but not
  in myofibrils.

## What already exists

- **Lab-meeting deck, 2026-10-06 (slide `o5-mava`, not committed).** Steady radial
  estimate, D_Pi ≈ 3×10⁻¹⁰ m²/s, 0.15 mM heads, turnover 9.4 s⁻¹ (8 mM) and 5.9 s⁻¹ (2 mM):

  | diameter | mean excess Pi, 8 mM | mean excess Pi, 2 mM |
  |---|---|---|
  | 100 µm | 1.5 mM | 0.9 mM |
  | 150 µm | 3.3 mM | 2.1 mM |
  | 200 µm | 5.9 mM | 3.7 mM |

- **Model code.** The Pi and ADP terms (Beard 2022 forms, `K_Pi = 4.007`, `K_D = 0.194`
  mM) exist in `Model/dPUdT_CombinedTransitions.m` but are commented out.
- **Literature.** `Analyses/HumanizeXB/literature_log.md` gives a measured ATPase for rat
  trabeculae at 20 °C of **0.43–0.49 mM/s** (Ebus & Stienen).
- **The "paper on Pi diffusion".** This is **not yet identified**. Neither the repo nor the
  bib holds it. Candidates are Ogut & Brozovich 2003 (CrP/diffusion) or a skinned-fibre Pi
  gradient paper. Ask the user.

## Order-of-magnitude bounds (computed 2026-10-07)

The model is a cylinder of radius R with uniform production q (mM/s).
- Steady centre excess is qR²/4D; steady mean excess is qR²/8D.
- The slowest diffusive mode relaxes with τ₁ = R²/(5.78 D).

| diameter | R²/8D | τ₁ | mean excess, q = 0.45 mM/s (measured ATPase) | mean excess, q = 2.1 mM/s (model XTOR 14 s⁻¹ × 0.15 mM) |
|---|---|---|---|---|
| 100 µm | 1.0 s | 1.4 s | 0.5 mM | 2.2 mM |
| 150 µm | 2.3 s | 3.2 s | 1.1 mM | 4.9 mM |
| 200 µm | 4.2 s | 5.8 s | 1.9 mM | 8.8 mM |

Reading the table:
1. **τ₁ ≈ 1–6 s**, comparable to how long the trabecula is activated before and during the
   slack train (force peaks 3.6–5.8 s after activation; the train spans 2.7 s). So Pi is
   still rising during the protocols, and the **transient**, not the steady state, is the
   relevant calculation. Pi build-up may also contribute to the within-record force decay
   that the rundown correction treats as damage.
2. **The answer is uncertain by ~4× in q.** The model's gross XTOR is at the top of its
   physiological bound (`project_xtor_is_gross`); the measured ATPase is ~4× lower. That
   uncertainty matters more than the uncertainty in D.
3. **Even the low estimate gives ~1 mM excess Pi** in a 150 µm preparation. Against an
   unknown bath Pi, that is not negligible: the force–log[Pi] slope is steepest at low Pi.

## Design

### Stage 0 — closed-form transient (cheap, do first)

- **Physics.** Solve ∂c/∂t = D∇²c + q(t) on a cylinder, with an elliptical-equivalent R
  taken from each preparation's width and thickness. Use the Bessel-series solution.
- **Inputs.**
  - q(t) from the model's own ATPase trace for each protocol, at 8 and 2 mM.
  - The measured ATPase as the lower bound.
- **Outputs.** Mean and centre [Pi](t) over each protocol window, per preparation and per
  ATP level.
- **Sweep.** D ∈ {2, 3, 6}×10⁻¹⁰ m²/s; bath Pi ∈ {0, 0.5, 1} mM, which must come from the
  solution recipe (TODO, Baker lab).

### Stage 1 — give the model a Pi sensitivity

- **Re-enable the Pi term.** It reverses the force-generating step, as in Beard 2022.
- **Normalise it to a reference [Pi]**, the same way the ATP gates are normalised, so the
  8 mM fit is unchanged at its own mean Pi.
- **Check its strength against the literature.** Simulate isometric force versus added Pi at
  8 mM. Compare the slope d(F/F₀)/d log[Pi] with Ebus & Stienen 2001 (cardiac). Then:
  - adjust only the Pi-binding constant;
  - record the slope;
  - do not refit anything else.

### Stage 2 — mean-field attribution

- **Run.** Run the 2 mM sweep core (`robust1_seed` + the k2/ksr0 gates) with [Pi] = the
  Stage 0 mean at each ATP level, instead of 0.
- **Read out.** The 2/8 ratios of the §1.6 features with and without Pi accumulation.
  Their difference is the **Pi-relief share** of each ATP effect.
- **Refit check.** Re-fit the single ATP lever (k2 exponent) with Pi accumulation included,
  and report how much the gate strength drops.

### Stage 3 — radially resolved (only if Stage 2 says it matters)

- **Discretise into shells.** Use N = 5–8 concentric shells. Each shell runs the
  cross-bridge model at its own local [Pi] (and [ATP] / [CrP] if Stage 0 shows ATP-side
  gradients).
- **Couple the shells.** Shells exchange Pi by finite-volume diffusion. Measured force is
  the area-weighted sum.
- **Cost.** N× a normal evaluation, so run it on the parpool. Stage 3 matters because a
  rigor-loaded core plus a well-supplied shell is a different *shape* signature
  (force up, rate down, Vmax unchanged), not only a different level.

### ATP / CrP side (in parallel with Stage 0)

- **Compute CrP consumption and diffusion** in the same geometry. CK-buffered ATP stays
  near bath until local CrP falls, so the question is whether core CrP runs down within an
  activation at 2 mM bath ATP.
- **This is the Ogut & Brozovich mechanism.**

## Decision rule

| Pi-relief share of the force effect | consequence |
|---|---|
| < 20 % (ratio shift < 0.07) | report as a quantified limitation in §3.1; headline stands |
| 20–50 % | report the Pi-corrected ATP effect alongside the raw one; refit the gate with Pi included |
| > 50 % | the headline is about metabolite gradients, not ATP; rewrite §1 |

## Experiments that decide it (for the Baker lab)

1. **Width *and* thickness of every preparation.** This is free and needed for Stage 0.
   It also enables test 3.
2. **Added Pi (0 / 5 / 10 mM) at 8 mM ATP** on the same preparation. This measures the
   force–Pi slope in *these* preparations, which Stage 1 currently borrows from the
   literature.
3. **ATP effect versus thickness across preparations.** A transport artefact scales with R².
   A molecular effect does not.
4. **Pi scavenger.** Add the MESG / purine-nucleoside-phosphorylase Pi mop to both ATP
   solutions. If the force gain at 2 mM persists with Pi held near zero, it is not Pi
   relief. This is the decisive single experiment.
5. **Optional:** omit CrP/CK in one run to measure the ADP side directly.

## Files to create when run

- `RunPiTransient.m` (Stage 0)
- `CalibratePiSlope.m` (Stage 1)
- `RunPiAttribution.m` (Stage 2)
- `results/`
- `conclusions.md`
