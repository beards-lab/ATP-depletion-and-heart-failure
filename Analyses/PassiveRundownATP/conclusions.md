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

---

# CORRECTION + confirmation (2026-08-10, second pass)

The first pass got two things wrong. Both are fixed here.
Script: `RunPassivePostRestretch.m`.

## 1. Direction — I used the wrong reference point; the hypothesis is RIGHT

The model's passive is fitted to **PNB+Mava**, which ran **last** and is therefore the **lowest**
of the three. So relative to the model, **both** active runs should want **more** passive:

    model (= PNB+Mava)  x1.00   |   2 mM  x~1.04   |   8 mM  x~1.08

(03/27 loses 3.8 % per activation interval; PNB+Mava sits after both activations, so ~7.6 % below
the 8 mM state.) My first pass compared the 2 mM optimum against the *2-vs-8 difference* (x0.973)
and called the observed x1.10 "wrong direction". **That reference was wrong** — against the
PNB+Mava-fitted baseline, x1.10 is the *expected* sign, and close to the predicted x1.04.
**Direction confirmed.**

## 2. Locus — passive acts on the SECOND peak and the decay, not the first

The first pass reported only `peak1_y` and `steady`, which is why it saw nothing. Scanning `k_pas`
at 2 mM and reading the post-restretch observables:

| `k_pas` × | peak2 | vall2_dy | ovrsht_dy | rsK | **peak1_y** | steady |
|---|---|---|---|---|---|---|
| 0.5 | 80.41 | −22.07 | 0.036 | 58.9 | **113.52** | 86.05 |
| 1.0 | 82.29 | −21.51 | 0.212 | 62.0 | **113.53** | 87.85 |
| 1.5 | 84.11 | −20.99 | 0.321 | 65.6 | **113.55** | 89.50 |
| 2.0 | 85.88 | −20.83 | 0.058 | 67.5 | **113.56** | 91.44 |
| 3.0 | 89.25 | −20.40 | 0.151 | 65.0 | **113.59** | 95.27 |
| *data 2 mM* | *92.43* | *−14.58* | *0.445* | *15.9* | *111.93* | *91.25* |

**`peak2` moves +11 % across the scan while `peak1_y` moves 0.06 %.** The hypothesis that passive
shapes the *second* peak and the post-restretch decay is **confirmed**, and it corrects the first
pass's "the probe is blind" verdict — the probe was blind only because I read the wrong peak.

## 3. Does the rundown differential help? Yes, correctly signed — but small

| case | rsK8 | rsK2 | **rsK ratio** | peak2 (2 mM) | 2 mM cost |
|---|---|---|---|---|---|
| model as fitted (1.00 / 1.00) | 102.5 | 62.0 | 0.605 | 82.29 | 16.789 |
| **rundown differential (1.08 / 1.04)** | 110.1 | 62.3 | **0.565** | 82.44 | **16.680** |
| exaggerated 1.5 / 1.2 | 108.5 | 64.7 | 0.596 | 83.03 | 16.800 |
| exaggerated 2.0 / 1.3 | 109.1 | 65.1 | 0.597 | 83.39 | 16.792 |
| *data* | *43.7* | *15.9* | ***0.363*** | *92.43* | |

The predicted differential moves the `rsK` ratio 0.605 → 0.565 (~13 % of the way to 0.363) and
improves the cost. **Real and correctly signed, but it is not the `rsK` answer** — and note that
*more* passive makes `rsK` itself *faster* (58.9 → 67.5), i.e. away from the data.

## 4. The bigger finding that fell out: the model's passive looks ~3× too small

`peak2` is a 10-unit residual (model 82.29 vs data 92.43) and **more passive closes it** — x3.0
reaches 89.25. Independently, the model's settled passive is **~2 % of steady force where the
rundown analysis measures ~6 %** (5 kPa of 80). **Two independent routes both say the model's
passive is roughly 3× under-scaled.**

That is far beyond anything rundown can explain (~4–8 %), so it is a *separate and larger* issue
than the ATP question — but it is plausibly why `peak2`, `vall2_dy` and `ovrsht_dy` are all
persistent residuals, and those are three of the top terms in the 2 mM cost.

## Revised verdict

- The rundown→passive→ATP-bias chain is **real, correctly signed, and now correctly located**
  (second peak + decay). Its magnitude on 03/27 is small (cost 16.79 → 16.68).
- It is **not** the explanation for `rsK` — passive moves `rsK` the wrong way.
- **The actionable finding is the passive scale itself**: fix the ~3× under-scaling (and check it
  against the PNB+Mava `PS_*` residuals) before drawing any further conclusion from `peak2`,
  `vall2_dy` or `ovrsht_dy`. That single correction touches three of the largest 2 mM residuals.
- 04/10 remains the strongest test of the rundown-differential itself, because its acquisition
  order is reversed and the sign should **flip**.

---

# Third pass — CAN the passive be scaled up, and by how much? (2026-08-10)

Script: `RunPassiveScaleFeasibility.m` · figure: `results/passive_feasibility.png`

## The acquisition order — the user was right, and it is worse than assumed

From the raw folder `data/03 27 2026 M/`:

    01_Relax → 02_8mM_Active → 03_2mM_Active → 04_8mM_Active_repeat → [05 MISSING] → 06_8mM_PNB_Mava

The passive protocol sits after **four** activations, not two. (`05` is absent from the series —
something was run and not exported.) The bracket **02 vs 04** is the −17.5 % force measurement that
[RundownCorrection](../RundownCorrection/conclusions.md) §4b is built on, and its parallel-material
component is **~16 % at full dose**. So passive at the time of PNB+Mava should sit **~16–20 % below**
the first 8 mM run — more than the 7.6 % used in the second pass.

## Do we have SL? **No.**

The raw ASI600A export carries three channels only:

    Time (ms) | L in (Lo) | F in (kPa)

That is **motor length in Lo units — there is no measured sarcomere length**. The `SL` column in the
curated `.mat` is *derived* from motor length. Consequence: the passive–SL relation cannot be
established independently of the model's own series-compliance assumption, and since rundown
*adds ~35 % series compliance*, the ML→SL map is itself changing through the session. Any passive
claim resting on SL is therefore model-internal, not measured.

## The feasibility test

Three constraints, one axis (`k_pas`). (A) is the **direct** measurement — what passive was fitted
to. (B) is the **indirect** demand from the 2 mM slack run.

| `k_pas` × | (A) PNB+Mava passive cost | (B) peak2 (2 mM) | (B) 2 mM total cost |
|---|---|---|---|
| 0.50 | 2.383 | 80.41 | 22.00 |
| 0.75 | 1.433 | 81.36 | 18.51 |
| **1.00** | **1.115** ← (A) optimum | 82.29 | 16.79 |
| **1.25** | 1.398 | 83.21 | **16.65** ← (B) optimum |
| 1.50 | 2.240 | 84.11 | 17.16 |
| 2.00 | 5.475 | 85.88 | 23.03 |
| 3.00 | 17.262 | 89.25 | 65.09 |
| **4.00** | **34.711** (31× worse) | **92.43** = data | 201.45 |

![passive feasibility](results/passive_feasibility.png)

## Answers

**1. My "~3× under-scaled" claim (second pass) was WRONG — retracted.**
The direct passive measurement prefers **×1.00**. The model's passive is correctly scaled *against
the passive data it was fitted to*. The earlier inference came from comparing the model's settled
passive (~2 % of steady) with the rundown write-up's "~5 kPa of 80" (~6 %) — those are different
length points and different definitions, so the comparison was invalid.

**2. Can we reasonably scale the passive up? YES — to ≈×1.2, and three independent lines agree.**

| line of evidence | says |
|---|---|
| rundown, 4 activations × 16 % parallel loss | ×1.16 – 1.20 |
| what the 2 mM slack total cost prefers | **×1.25** |
| what the direct passive measurement tolerates | ×1.20 costs only **+0.23** (1.115 → 1.341) |

A ~20 % upward passive correction for the **active** runs is therefore **defensible on three
independent grounds**, cheap on the direct fit, and it is exactly the user's hypothesis quantified.
It should be applied as a *per-run* correction (8 mM more than 2 mM), not as a change to `k_pas`
itself — the fitted `k_pas` is right *for the PNB+Mava condition*.

**3. Must we postulate a scaling "because it fits"? NO — and the data forbids the big one.**
Closing `peak2` on passive alone needs **×4.0**, where the direct passive cost is **34.7 vs 1.1 —
31× worse**. That is a hard experimental rejection, on data already in the cost function.
**`peak2` is therefore NOT a passive-scaling problem**; something else drives it.

## Consequences

- Apply the ≈×1.2 passive correction to the active runs (differential: 8 mM > 2 mM). Expected
  effect is modest — the second pass measured `rsK` ratio 0.605 → 0.565 and cost 16.79 → 16.68 for
  the 1.08/1.04 differential.
- Stop attributing `peak2` to passive. Note also that
  [SlackDataAnalysis](../SlackDataAnalysis/conclusions.md) already retires **`ovrsht_*`** as
  unreliable across preps — and `ovrsht_dy` is the *second largest* term in the current 2 mM cost
  (2.77). Two of the top four 2 mM residuals (`ovrsht_dy`, and `rsK` via slack-compliance
  contamination) are therefore of questionable evidential value, which materially changes what
  "fitting the 2 mM data" should even mean.
- Getting further on the passive question needs **measured SL** (striation/laser diffraction),
  which this rig did not record. Without it the ML→SL map is model-internal and moves with rundown.
