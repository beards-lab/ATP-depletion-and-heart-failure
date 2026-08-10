# Can a rundown-driven passive difference serve as a second handle on high-vs-low ATP? (2026-08-10)

**Hypothesis (user).** The PNB+Mava passive run was recorded **last**, so it measures passive in
the **most** rundown state. If rundown removes force-generating material in parallel, the torn
myofibrils also stop carrying **passive** (titin) force, so

    passive(8 mM run)  >  passive(2 mM run)  >  passive(PNB+Mava, as fitted)

The model fits passive to the PNB+Mava level, therefore under-assigns passive to *both* active
runs, and more so at 8 mM — which would bias the ATP comparison and offer a second, independent
handle on the high-vs-low ATP difference.

Reproduce: `RunPassiveRundownTest.m`.

## Step 1 — consistency with [`../RundownCorrection`](../RundownCorrection/conclusions.md)

| premise | status |
|---|---|
| `kSE` changed (series compliance) | ✅ **confirmed** — rundown adds ≈35 % series compliance |
| `SL0` changed | ⚠️ **revised**: "provisional and resting on a single measurement" (§4). On five points the loss is a **scale**, not a shift. |
| PNB+Mava ran last / reads the most-rundown state | ✅ confirmed — it "is charged zero damage" (§3) |
| passive falls with rundown | ⚠️ **contested inside the source** — see below |

**A genuine tension in the existing analysis, and it favours the hypothesis.** §4 asserts in
passing that *"passive does not run down"*. But the same section's **revised** lesion is *"loss of
force-generating material in parallel (≈16 %) plus added series compliance (≈35 %)"* — the "tearing
the cell" picture. Torn-away myofibrils cannot carry titin force either, so the revised mechanism
**logically implies passive falls with rundown**, by roughly the same parallel fraction. The
hypothesis is consistent with the revised mechanism and conflicts only with that older aside.

## Step 2 — quantitative prediction, made before touching the model

Passive ≈ 5 kPa (the PNB+Mava level, §4); parallel loss ≈ 16 % at the **full** bracket dose;
03/27's dose is **λ = 0.17**. So between the 8 mM and 2 mM runs:

    material lost ≈ 0.17 × 16 % = 2.7 %   →   passive difference ≈ 5 kPa × 0.027 ≈ 0.14 kPa

Against a total force of ~90 kPa and an ATP force effect of ×1.36 (≈25 kPa), that is **0.15 % of
force and ~180× smaller than the signal it would correct.**

> **Prediction: real in direction, negligible in magnitude.**

## Step 3 — the test

Free the passive scale at the 2 mM condition only and ask what the data wants.

| `k_pas` × | 2 mM cost | peak passive | steady | **peak1_y** | vall_y |
|---|---|---|---|---|---|
| 0.400 | 24.12 | 17.01 | 85.64 | **113.52** | 72.74 |
| 0.800 | 17.97 | 18.97 | 87.16 | **113.53** | 73.14 |
| **0.973** ← rundown predicts | 16.94 | 19.83 | 87.76 | **113.53** | 73.31 |
| 1.000 (baseline) | 16.79 | 19.96 | 87.85 | **113.53** | 73.34 |
| **1.100** ← data prefers | **16.58** | 20.45 | 88.19 | **113.54** | 73.43 |
| 1.300 | 16.79 | 21.44 | 88.85 | **113.54** | 73.62 |
| 3.000 | 65.09 | 29.52 | 95.27 | **113.59** | 75.12 |
| *data* | | | *91.25* | *111.93* | *85.19* |

### Three findings, all negative for the hypothesis as a usable handle

1. **The direction is wrong.** The 2 mM data prefers **more** passive (×1.10), not the ×0.973 that
   rundown predicts.
2. **The magnitude is unresolvable.** The whole passive axis buys **0.209 of 16.79 = 1.2 %** of the
   cost, and the curve is flat: over ×0.9–1.3 the cost moves only 0.63. The predicted 2.7 % change
   sits far inside that flat region — **this data cannot resolve it even in principle.**
3. **The proposed probe has no sensitivity.** The hypothesis pointed at the restretch **peaks**.
   But `peak1_y` moves **113.52 → 113.59 (0.06 %) while passive changes 7.5×**. The restretch peak
   is set by cross-bridge stiffness, not titin. **Peaks cannot detect a passive difference.**
   The passive-sensitive observable is **`steady`** (85.6 → 95.3 over the same range).

## Verdict

The hypothesis is **mechanistically coherent** — and it exposes a real internal inconsistency in
the rundown write-up worth fixing. But as a *second handle on the ATP difference for 03/27* it
fails on all three counts: the predicted effect is ~180× smaller than the ATP signal, the observable
it proposed (peaks) has essentially zero passive sensitivity, and what little preference the data
does express points the other way.

**Do not use passive as an ATP discriminator on 03/27.**

## Where it could still be worth testing

**04/10** is the one prep with leverage:
- dose **λ = 0.83** — 5× 03/27's — so the predicted passive difference is ≈13 %, ≈0.65 kPa;
- it ran **2 mM first**, so the sign **flips** (there the 8 mM run should show the *lower* passive).

A sign reversal between preps is a much stronger test than a magnitude on one prep, because no
competing model error would flip with acquisition order. The data (`protocol_04_10_2026_*`) is
already curated and present in `data/`.

Two caveats to carry into that test:
- The model's settled passive is **1.6–2.1** in force units where steady ≈ 88, i.e. ~2 %, whereas
  the rundown analysis quotes ~5 kPa of ~80, i.e. ~6 %. **The model may under-represent passive by
  ~3×**, which would need resolving before passive is used quantitatively for anything.
- Peak passive differs slightly between the two runs here (20.28 vs 19.96) purely because the two
  protocols have different length trajectories — not an ATP or rundown effect. Any future
  comparison must be made at matched length.
