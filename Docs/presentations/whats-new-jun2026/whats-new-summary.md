# Cross-Bridge / ATP-Depletion Model — What's New (late May – June 2026)

*Draft summary for the combined project presentation. Audience: PI update (mixed technical). Two slides' worth of material below — I'll trim once you're happy with it.*

---

## One-line version

- I re-grounded the cardiac cross-bridge model in α-myosin physiology and tightened an already-working fit (output error from roughly **~20 down to ~9**).
- The real win: I **identified and added a missing mechanism** that finally lets the model produce the force–velocity *shoulder* near isometric.
- The new mechanism **beats the optimizer's best fit** and **cuts the force–velocity error by ~64%** — but the shoulder is improved, **not yet fully closed**.

---

## The headline result — the force–velocity "shoulder"

The phenomenon:

- The experimental force–velocity curve has a **flat shoulder**: when the muscle starts to shorten slowly, force barely drops (only ~8% at v = 0.5 ML/s).
- Every version of my model fell off a cliff instead (~40% drop).
- This gap was stuck across *every* parameter I tried — it looked like a structural limit.

What I found and did:

- It isn't structural. I identified the likely missing physics as **target-zone registration availability**: at a fixed (isometric) position only a fraction of myosin heads line up with an available actin binding site; *sliding* continuously refreshes that alignment, so more heads bind and force is *maintained* as shortening begins.
- I implemented it as a single extra slow state variable that gates attachment, driven by the imposed sliding velocity — not the noisy instantaneous velocity that broke my earlier attempts.

The result — it helped, but isn't finished:

- At v = 0.5 ML/s the model goes from 0.62 to **0.73** (data 0.92).
- The force–velocity sub-error drops **~64%**, and the overall fit **beats the optimizer's best** (9.9 → 8.9).
- It does this *without* breaking force redevelopment (ktr) or the restretch transients that earlier velocity-gated attempts destroyed.
- It's **backward-compatible** (off by default; identical to the old model when disabled) and makes a **testable prediction**: the shoulder should fade at low ATP.
- A clear gap to the data remains — closing it fully is ongoing (deeper availability, residual ktr/restretch retuning).

![Force–velocity shoulder — improved, not yet closed](figures/fig_fv_shoulder.png)

---

## The fit campaign — a working fit, refined and grounded

- My starting point was an already-working fit with output error in the **~20 range**. (An early reseed briefly spiked the error while I moved parameters around, but that was a transient artifact, not the real baseline.)
- From there I ran a systematic, physiologically-bounded campaign down to **~9**, with each improvement traceable to a specific mechanism rather than blind fitting:
  - **XTOR fix** — pinned ATP turnover via the hydrolysis rate (kah) and detachment (k2), bringing ATPase into its physiological window.
  - **Manual rate tuning** — found attachment rate (ka) is the master force/duty knob; nailed all slack-force magnitudes and removed a redevelopment oscillation artifact.
  - **Length-dependent activation via lattice spacing** — got the force–length slope right, sourced from *attachment* rather than ad-hoc feedback.
  - **Binding-site occupancy** (below) — capped the bound fraction at the physical ceiling and fixed an over-binding problem.
  - **Force-bearing pre-stroke state + registration availability** — decoupled ktr from force–velocity and reopened the previously-stuck shoulder.

The bar chart below is the cleanest single comparison: **the same optimized parameter set with the new mechanism off vs on** (same base, mechanism toggled). It lowers total error (9.9 → 8.9) and cuts the force–velocity sub-error by ~64%.

![New mechanism vs optimizer's best fit (same base parameters, mechanism off vs on)](figures/fig_fit_cost.png)

---

## Mechanisms I switched on and tuned

Several model switches that were previously off (or mis-used) had to be turned on, grounded, and tuned together to reach this fit. Each carries distinct physics rather than being a free fitting knob:

- **Length-dependent activation via lattice spacing** (`UseLatticeSpacing`, d_optimal ≈ 0.0242, σ ≈ 0.0025) — makes attachment depend on filament lattice spacing, producing the cardiac force–length slope from physiology instead of ad-hoc feedback.
- **Binding-site occupancy** (`UseGlobalOccupancySaturation`, linear form, P_bound_max ≈ 0.3–1.0) — a global bound-fraction cap (~12–15% rat-cardiac ceiling) that fixes the over-binding and the ATPase (XTOR) over-shoot.
- **Force-bearing pre-stroke state** ("Route B", dr1 ≈ 0.0047, kstiff1/kstiff2 retuned) — lets the weakly-bound p1 state carry force, which decouples force redevelopment (ktr) from the force–velocity curve.
- **Registration availability** (`UseRegistrationAvailability`, A0 ≈ 0.6, v_ref_reg ≈ 0.8, tau_reg ≈ 0.02 ≈ 1/ktr) — the new sliding-gated attachment state behind the FV shoulder.
- **Base cycle rates re-tuned and grounded** — attachment (ka), power-stroke (k1), detachment (k2) and hydrolysis (kah, fixed at the α-myosin floor of 80), plus the R2 detachment-vs-strain curve, set ATP turnover, duty ratio and absolute force.

---

## Binding-site occupancy — a conceptual correction

- While auditing how the model limits attachment, I found the pre-existing per-strain-bin "occupancy" term was a **category error**: the model's strain axis is *not* a binding-site axis, so a per-bin saturation term penalizes heads that merely share a strain value, not heads competing for a site.
- I replaced it with the only faithful variable in this model class — a **global bound-fraction limit** — grounded in the rat-cardiac physical ceiling (~12–15% of motors bound, from the same lab that produced my data).
- The global form is grid-independent (identical to 5 digits across a 2× grid refinement) where the old per-bin form drifted ~45%.
- With stiffness recompensation it cuts fit error and fixes the ATPase over-binding.

---

## Infrastructure & grounding

- **Parameter bounds re-grounded to α-myosin** (rat/mouse cardiac), so every fitted rate now sits in a defensible physiological range rather than floating.
- **Repository reorganized** into a clear knowledge map: a single `Model/` for core code, one self-contained folder per analysis (each with its own `conclusions.md`), and a structured `Docs/` knowledge base. Dead code archived; two clear "front doors" (`Docs/README.md`, `Analyses/README.md`).
- **Diagnostic probes and reporting tools** added (SRX, ktr, mechanics, isometric) plus a sorted cost-breakdown reporter that makes each fit residual traceable.

---

## Where it stands / open threads

- **Solid:** force–length relation, force redevelopment (ktr, clean single-order), ATP turnover, slack-force magnitudes.
- **Improved but open:** the FV shoulder — the new mechanism helped substantially but a gap to the data remains.
- **Next:** fully closing the shoulder, residual ktr-vs-FV tuning, the restretch transient (catch-bond territory), and the **low-ATP shoulder-fade test** — a direct check of the new mechanism that ties straight back to the project's ATP-depletion / heart-failure motivation.

---

*Figures generated 2026-06-22 from the logged fit results (`Analyses/FV_Shoulder/conclusions.md` §11, `Docs/notes/labdiary_boundedfit.md`). The fit-error bars compare `params/params_reseeded_opt.m` (mechanism off, cost 9.90) against `params/params_reseeded_regavail.m` (same base, mechanism on + kstiff×1.28, cost 8.89). The before/after FV curves use the logged @v=0.5 values (0.62 → 0.73; data 0.92) with representative shapes for the other velocities — I can regenerate them from actual MATLAB output if you want them exact.*
