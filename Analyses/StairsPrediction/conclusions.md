# Staircase prediction, and the `estiff` nonlinearity

**Date:** 2026-09-04, revised 2026-09-05 · **Incumbent:** `params/spreadw_opt.m`
(e = 0.9653, `Et` = 0.9786) · **Challenger (adopt this):** `params/estiff125_opt.m` (e = 1.25, `Et` = 0.9559). Also
tested: `params/ec080_opt.m`, `params/ec090_opt.m`, `params/ec097_opt.m` (the control),
`params/ec110_opt.m`, `params/ec125_opt.m`, `params/ec140_opt.m`, `params/ec160_opt.m`

---

## 0. The answer, in one place

**Does moving `estiff` help the staircase, and does it cost anything already paid for?**

> **Yes — `estiff` = 1.25 is better on everything measured, and impairs nothing.** It
> beats the incumbent on the scored objective, on all three held-out staircases, and on
> the 2 mM ATP response. The incumbent value 0.9653 is not a fitted optimum; it is where
> a 12 h search got stuck, for a reason that is now measured (§4.2–4.3).

**Everything at a matched 1.5 h budget, plus the control:**

| | `Et` (scored) | stairs 03/27 | 04/03 | 04/10 | budget |
|---|---|---|---|---|---|
| incumbent, e = 0.9653 | 0.9786 | 5.415 | 5.230 | 2.484 | 12 h |
| e = 0.965 refit *(control)* | 0.9572 | 5.157 | — | — | 0.5 h |
| e = 0.90 | 0.9580 | 4.000 | — | — | 1.5 h |
| **e = 1.25** | **0.9559** | **1.595** | **2.875** | **1.369** | 1.5 h |

Three things make this stand up where the previous two verdicts did not:

1. **The control (§5.1).** Refitting at the *incumbent* exponent, same optimiser, same
   budget, **fully recovers the objective** (0.9572 — better than the incumbent itself)
   and **leaves the staircase at 5.16**. The staircase gain is therefore not "another
   half hour of optimiser". It is the exponent.
2. **The matched budget (§5.3).** At 0.5 h the campaign says *lower* `estiff` is better —
   the old, retracted answer, reproduced for the old reason (at that budget the refit has
   barely left its seed, and the seed staircase is monotone in `estiff` by construction).
   Given the same 1.5 h as e = 1.25, e = 0.90 improves its objective to 0.958 but its
   staircase gets **worse**, 3.42 → 4.00. Head to head at equal budget, e = 1.25 wins on
   the staircase by 2.5× at the same objective.
3. **Low ATP agrees, independently (§6.2).** e = 1.25 is closest to the measured 2 mM
   ratios on all three ATP-sensitive quantities (force +12 %, `ktr` +34 %, `rsK` +21 %
   closer than the incumbent), while e = 0.80 and e = 0.90 are *worse* than the incumbent
   on all three. This comes from neither the objective nor the staircase.

**What it costs.** Inside the objective the only term that pays is FV (+0.085), while
`rsK`, `ktr`, `doublePeak`, `vall2_dy` and `PS_holdDecayRMSE` all improve. Outside it,
the quasi-static staircase envelope gets worse (block 0.023 → 0.065) and the two coarse
days' raw trace LSQ worsens (04/03 0.100 → 0.151, 04/10 0.0153 → 0.0280) even as their
feature costs improve. That is §2's point made concrete — the trace MSE is ~90 % hold
plateaux — and the measured weights back the transient (envelope 35 % between-prep
spread; 1 % peak 21 %; recovery 17 %). But it *is* a trade, and it is the one real cost.

**What changes structurally.** The blocks that move are the ones that were *wrong*, not
merely mis-scaled: phase-2 nadir 2.93 → 0.20, delayed recovery 1.22 → 0.26, elastic peak
1.25 → 1.07. The model finally produces rip-off-and-recover instead of a slow decay.

> ### Note on the history of this document
> An earlier version of this report concluded the exact opposite — that e > 1 was
> "the wrong sign, monotonically, on all three preparations". That verdict came from
> sweeping `estiff` with `kstiff2` set by an iso-force rule, and it was an artefact of
> the sweep, not a property of the exponent. §4 explains the failure in full, because
> the failure mode is more generally useful than the result: **a clean, smooth,
> internally-consistent low-dimensional profile can still be the wrong experiment.**
> The retracted material is summarised in §4.4 rather than deleted.

---

## 1. Why this is a prediction and not a fit

`params/spreadw_opt.m` carries `RunStairs = 0`. The staircase has never been in any
objective, so everything below is genuine out-of-sample. There are three staircase
records, and two of them use a protocol the model has never been shown:

| day | protocol | steps |
|---|---|---|
| 03/27 M | fine | 7 × 1 % ML, then 1 × 3 % |
| 04/03 F | coarse | 1 × 1 %, then 3 × 3 % |
| 04/10 M | coarse | 1 × 1 %, then 3 × 3 % |

Each step is a ~20 ms ramp followed by an ~80 ms hold, on a 4 kHz record.

## 2. How to score it: features, not least squares

Every step gives the classic three-phase response

```
base --(ramp)--> PEAK --(fast decay)--> NADIR --(slow rise)--> PLATEAU
```

**The raw MSE is mostly noise.** Detrended residuals inside the isometric holds give a
per-sample σ of 0.0245 (03/27), 0.0333 (04/03), 0.0279 (04/10) in units of F₀. The
`runStairsExperiment` cost is `mean(residual²)·10`, so a *perfect* model still scores
`10σ² = 0.0060 / 0.0111 / 0.0078`. The incumbent scores 0.0267 on 03/27 — only 4.4× the
floor, i.e. roughly a quarter of the number is irreducible scatter. Worse, ~90 % of the
samples sit in the holds, so the MSE spends ~90 % of its weight there and almost none on
the ramp transient, which is where the model is actually wrong.

**Window means are ~20× more precise than the trace.** Each feature is a mean over a
fixed window, never an `argmax` (a max over a ramp window is biased upward by ~1σ, about
15 % of the signal). With n ≈ 24–72 samples per window the standard errors are

| feature | SE (03/27) |
|---|---|
| peak − base | 0.0059 |
| nadir − base | 0.0048 |
| plateau − F₀ | 0.0041 |

against a peak amplitude of 0.091 — so the per-step features resolve the response at
~1.5 % of F₀. `extractStairsAttributes.m` implements this; `stairsCost.m` sums four
block costs, each `mean_k (m_k−d_k)² / rms_k(d_k)²`, giving the four *phases* equal say
instead of letting sample count decide.

**Between-preparation spread** (step-size matched across the three days, the
`Analyses/FeatureWeighting` χ² convention `w = λ/s²`):

| quantity | n | mean | rel. SD | w ∝ |
|---|---|---|---|---|
| peak @1 % | 9 | 0.0850 | 21 % | 21.8 λ |
| peak @3 % | 7 | 0.0971 | 46 % | 4.7 λ |
| nadir @1 % | 9 | −0.0256 | 26 % | 14.4 λ |
| nadir @3 % | 7 | −0.0154 | 217 % | 0.2 λ |
| recovery @1 % | 9 | 0.0498 | 17 % | 35.1 λ |
| recovery @3 % | 7 | 0.0955 | 27 % | 14.2 λ |
| envelope 2.0→2.2 µm | 3 | 0.2593 | 35 % | 8.2 λ |

The **1 % step peak, nadir and recovery** are the well-determined observables. The 3 %
nadir is not measurable at all (217 % spread) and should never be scored.

## 3. What the incumbent predicts

| block | 03/27 | 04/03 | 04/10 |
|---|---|---|---|
| elastic peak | **1.25** | **4.11** | 0.37 |
| phase-2 nadir | **2.93** | 0.65 | 1.58 |
| envelope | 0.02 | 0.17 | 0.01 |
| delayed recovery | **1.22** | 0.30 | 0.52 |
| total | 5.42 | 5.23 | 2.48 |
| LSQ (floor) | 0.0267 (0.0060) | 0.1001 (0.0111) | 0.0153 (0.0078) |

**Right:** the quasi-static staircase envelope. Plateau-vs-step slope is off by 6.5 %,
block cost 0.02. The model's steady force–length behaviour over 2.0 → 2.2 µm is good.

**Wrong, and by a lot:**

1. **The elastic peak is 2.0–4.1× too large.** Data 0.073–0.107 F₀ per 1 % step; model
   0.184–0.198. That is 15–20 SE.
2. **The phase-2 undershoot has the wrong sign.** Data dips to −0.028 F₀ *below* the
   pre-step baseline; the model sits at +0.025 *above* it for the first six steps.
3. **The delayed recovery is missing.** Data climbs +0.05 F₀ back during each hold; the
   model climbs ~0 and on steps 3–7 drifts slightly down.

So the model gets the *equilibrium* right and the *transient* wrong: too stiff on the
way up, then no rip-off and no re-recruitment.

---

## 4. Why every `estiff` slice lied — and why the optimiser never moved it

`estiff` is the exponent e in the post-stroke force sum
`F₂ = kstiff2 · Σ sign(s+dr)|s+dr|^e p₂(s) dS`.

### 4.1 The coupling, analytically

At the isometric operating point the whole attached population is a narrow spike:
⟨x⟩ = 0.0102 µm, SD 0.0023 µm (this needs the *dynamic* `out.LXBPivot`, not
`params.LXBpivot`). So

```
d(log F₂)/d(estiff) = log⟨x⟩ = −4.59   per unit e
```

A 1 % change in `estiff` is a 4.4 % change in active force unless `kstiff2` follows it
exactly. Between e = 0.965 and e = 1.5 the required `kstiff2` correction is a factor of
12. **Any experiment that changes `estiff` without re-optimising `kstiff2` to high
precision is measuring `kstiff2`, not `estiff`.**

### 4.2 The wall and the valley floor — measured

`RunEstiffWall.m`. Two curves about the incumbent: `estiff` moved alone, and `estiff`
moved with the iso-force `kstiff2` correction.

| Δe | −0.05 | −0.02 | −0.01 | −0.005 | **0** | +0.005 | +0.01 | +0.02 | +0.05 |
|---|---|---|---|---|---|---|---|---|---|
| F_iso error if uncompensated | +25.7 % | +9.6 % | +4.7 % | +2.3 % | 0 | −2.3 % | −4.5 % | −8.7 % | −20.4 % |
| **`estiff` alone** | 197.9 | 5.450 | 1.517 | 1.088 | **0.979** | 1.124 | 1.770 | 58.75 | 108.97 |
| **`estiff` + iso `kstiff2`** | 1.220 | 1.029 | 0.995 | 0.988 | **0.979** | **0.966** | 0.970 | 0.979 | 1.091 |

and the same wall in the other coordinate — `kstiff2` moved alone:

| `kstiff2` | ×0.90 | ×0.95 | ×0.98 | ×0.99 | **×1** | ×1.01 | ×1.02 | ×1.05 | ×1.10 |
|---|---|---|---|---|---|---|---|---|---|
| `Et` | 59.20 | 57.12 | 1.084 | 0.999 | **0.979** | 0.990 | 1.044 | 1.560 | 5.970 |

Three things follow, and together they are the whole story:

1. **The valley is ~60× deep at a 2 % excursion.** Moving `estiff` by 2 % with `kstiff2`
   held costs `Et` 0.979 → 58.75.
2. **The compensating tube is ~1 % wide.** `kstiff2` ×0.98 already costs +11 %; ×0.95
   costs ×58.
3. **The valley floor is flat.** Along the compensated ridge `Et` stays inside
   0.966–1.09 over the whole range e ∈ [0.915, 1.015] — and is *minimised off the
   incumbent*, at e ≈ 0.970 (0.9660 vs 0.9786). `estiff` is locally **unidentifiable**.

### 4.3 So why did a 12 h optimiser leave it alone?

Not because it couldn't. `estiff` **was** in the spreadw pool
([`RunOptimSpreadW.m:76`](../FeatureWeighting/RunOptimSpreadW.m)), the pool holds 88
parameters, `N_DRAW = 10` of which 6 are compulsory — and `kstiff2` is one of the
compulsory ones. So `estiff` was drawn with probability 4/82 = **4.9 % per round**, and
on every one of those rounds `kstiff2` was drawn *with* it. A 12 h run drew the pair
together something like nine times.

What it did with them, from the `params/spreadw_iter/` gallery:

| iter | 001 | 002 | 003 | 004 | 005 … 018 |
|---|---|---|---|---|---|
| `estiff` | 0.965542 | 0.965542 | 0.964708 | 0.965322 | 0.965322 (unchanged) |

A total excursion of **0.09 %**, then nothing across fourteen further *improved*
iterations. §4.2 says why: a flat floor gives the search no reason to move along the
ridge, and 60×-tall walls a few percent to either side punish every step that leaves it.
Reaching e = 1.25 is not a 2-D move at all — §4.4 shows it needs the FV strain-shape
terms to travel with the exponent, which is a coordinated ~12-parameter step that no
random 10-parameter draw is going to stumble onto.

**This is a statement about the search, not about the model.** "The optimiser did not
move it" is not evidence that the incumbent value is right.

### 4.4 The retracted slices — what was done and why it was wrong

Four experiments concluded e > 1 was the wrong sign. All four fixed `kstiff2` by the
iso-force rule and are retracted:

- a 1-D `estiff` scan on all three preparations (`RunEstiffScan`, `RunEstiffAllDays`),
- the same on the objective (`Et` 0.979 → 5.48 → 40.2 for e = 0.965 → 1.25 → 2.0),
- a 2-D `estiff` × `kSE` grid (`RunGrid2D`), in which no e > 1 node reached even the
  incumbent's staircase cost,
- an `estiff` **profile** (`RunEstiffProfile`) that re-optimised `kstiff2` against the
  paid-for objective at every exponent and *still* put both minima at e = 0.965.

The profile is the instructive one, because it is internally clean: it validated at the
incumbent (e = 0.965 re-found `kstiff2` ×1.001, `Et` 0.9765) and the correction to the
iso-force rule was smooth and small (×0.84 → ×1.05 across e = 0.85 … 1.60). It was still
wrong, because the coupling is not two-way. A full field-by-field diff of the repair fit
shows round 1 drew ten parameters and moved **twelve**: `kstiff2` ×3.98, plus three
strain-dependence shape terms, `d10_ref`, `gamma`, `kSE`, `k2` and `v_ref_reg` — the FV
levers, and FV was 84 % of the seed's cost. Decomposing that point:

| held back to the incumbent | `Et` | staircase |
|---|---|---|
| nothing (the repair result) | **1.163** | **1.674** |
| the FV shape terms | 6.420 | 16.24 |
| `estiff` + `kstiff2` | 39.0–39.2 | 3.48–3.62 |
| `k2` + `v_ref_reg` only moved | 1.280 | 5.314 |

**Removing either half destroys it.** One coupled solution, not a sum of separable
effects.

### 4.5 Two basins, and why one of them hides

Testing the exponent from both sides, each time re-optimising `kstiff2` so it is never
judged against a de-optimised stiffness. Rows are the exponent, columns are whose other
parameters you keep (`Et` / staircase):

| | incumbent's other params | refit's other params |
|---|---|---|
| **e = 0.965** | **0.975 / 5.415** | 5.669 / 7.638 |
| **e = 1.25** | 5.030 / 18.13 | **0.956 / 1.595** |

Every 1-D path between them goes uphill. The exponent is **load-bearing**: reverting it
to 0.965 at the refit point, with the best `kstiff2` available, costs 6× in objective and
4.8× in staircase.

The two basins also differ in *shape*. Displacing `kstiff2` by the same ×0.898 from each
optimum costs

| | `Et` at optimum | after ×0.898 on `kstiff2` |
|---|---|---|
| incumbent, e = 0.965 | 0.975 | **59.43** |
| refit, e = 1.25 | 0.956 | **5.48** |

an 11-fold difference in brittleness. The incumbent sits where a 10 % stiffness error
costs a factor of 60; the strain-stiffening solution tolerates the same error for a
factor of 5.6. That is an identifiability result in its own right, and it is the
mechanical reason a 9 %-off compensation could produce a confident, four-way, completely
wrong verdict.

---

## 5. The matched-budget campaign

`RunEstiffCampaign.m` → `CollectCampaign.m` → `estiff_campaign.png`. Seven frozen
exponents, each given

- an **iso-force-compensated seed** (`kstiff2` and `kstiff2_n` × `pw(e₀)/pw(e)`, so
  every seed starts on the analytic ridge of §4.1), and
- the **same optimiser, the same pool minus `estiff`** (88 → 87 parameters, `kSE` and
  `ekSE` deliberately kept in), with the **same 0.5 h budget**,

then scored on the held-out staircase. Equal budget from an equally-constructed seed is
the only comparison that is fair across exponents. e = 0.965 is included on purpose: its
iso-rule seed *is* the incumbent, so that row is the **matched control** — "can the same
staircase gain be bought at the incumbent exponent?"

| e | `kstiff2` × | seed `Et` | **opt `Et`** | seed stairs | **opt stairs** | final `kstiff2` |
|---|---|---|---|---|---|---|
| 0.80 | 0.470 | 54.515 | 1.288 | 1.811 | **1.983** | 10 178 |
| 0.90 | 0.742 | 1.376 | **0.982** | 3.653 | **3.424** | 16 160 |
| **0.965** *(control)* | 0.999 | 0.974 | **0.957** | 5.523 | **5.157** | 22 126 |
| 1.10 | 1.849 | 2.003 | 1.004 | 10.620 | 10.702 | 42 655 |
| 1.25 | 3.662 | 5.480 | 1.329 | 18.236 | 14.486 | 84 559 |
| 1.40 | 7.247 | 10.571 | 1.406 | 29.372 | 15.723 | 195 480 |
| 1.60 | 17.971 | 18.905 | 2.015 | 43.297 | 23.574 | 510 920 |

*(incumbent for reference: `Et` 0.9786, stairs 5.415)*

### 5.1 The control: the gain is the exponent, not the refit

**This is the experiment the first version of this report was missing.** Refitting *at*
the incumbent exponent for the same half hour recovers the objective completely —
`Et` 0.9744 → **0.9572**, as good as anything else in the table, and better than the
incumbent itself — while the staircase barely moves, 5.523 → **5.157**.

So "half an hour of the same optimiser" is not what fixes the staircase. At e = 0.965
the staircase stays at ~5.2 no matter what the other 87 parameters do with that budget.
The exponent is doing real work.

### 5.2 At a SHORT budget, the useful direction looks like DOWN

The two headline points:

- **e = 0.90 is close to free.** `Et` 0.982 against the incumbent's 0.979 — inside the
  objective's own numerical roughness — with the staircase at **3.42 vs 5.42**, a 1.6×
  improvement, for 30 minutes of compute.
- **e = 0.80 buys more**, staircase **1.98** (2.7× better than the incumbent), for a 33 %
  objective cost (`Et` 1.288).

And the cost composition says this is not robbery. At e = 0.90, `FV_normpowerAvg` — the
single largest term, and the one the e = 1.25 solution had to *pay* — actually
**improves**, 0.1803 → 0.1623, as does `rsK` (0.0635 → 0.0582). At e = 0.80 FV improves
further (0.1421) and the bill is paid instead by `vall2_dy` (0.2008), `coolDownLS`
(0.1431) and `rsK` (0.1118).

**`doublePeak` is not the barrier §7 previously claimed.** The old note recorded that
`estiff` below ~0.95 trips a spurious second slack-restretch peak at cost 50–102. Under a
refit it does not: `doublePeak` scores **0.0215 at e = 0.80** and 0.0164 at e = 0.90,
against 0.0226 for the incumbent itself. That figure of 50–102 was measured on the same
uncompensated slices §4.4 retracts — the artefact, one more time. The low band is
admissible; it simply is not where the best solution turns out to be.

> **Read §5.3 before using any of this.** Everything in §5.2 is the 0.5 h column, and
> the low-`estiff` advantage in it does **not** survive a budget-matched comparison. It
> is retained because it is the exact shape of the trap this whole report is about: a
> clean, monotone, internally-consistent result that measures the search rather than the
> model.

### 5.3 …but at a MATCHED budget, e = 1.25 wins outright

§5.2 was written from the 0.5 h column. Extending e = 0.90 to **1.5 h**, the same budget
the e = 1.25 point had, reverses it:

| e = 0.90 | `Et` | stairs |
|---|---|---|
| seed (iso-force) | 1.3757 | 3.653 |
| after 0.5 h | 0.9821 | **3.424** |
| after 1.5 h | **0.9580** | **4.000** |

The extra hour bought a better objective (0.982 → 0.958) and a **worse** staircase
(3.42 → 4.00). At e = 0.90 the two criteria pull against each other: as the search
converges on the paid-for objective it drifts away from the staircase. Head to head at
1.5 h:

| | `Et` | stairs 03/27 | stairs 04/03 | stairs 04/10 |
|---|---|---|---|---|
| **e = 0.90** | 0.9580 | 4.000 | — | — |
| **e = 1.25** | **0.9559** | **1.595** | **2.875** | **1.369** |

Same objective to within the noise, and a **2.5× better staircase** for e = 1.25. The
low-`estiff` advantage in §5.2 was a short-budget artefact — the mirror image of the one
§5.4 warns about, caught by running the control.

### 5.4 The budget is load-bearing, and this cuts both ways

At 0.5 h the staircase column is **monotone increasing in `estiff`** — which is exactly
the retracted §4.1 verdict. The reason is visible in the table: at this budget the
refit barely moves the staircase off its seed (1.81→1.98, 10.62→10.70, 3.65→3.42), and
the *seed* staircase is monotone in `estiff` as a mechanical consequence of the iso-force
ridge. **A short-budget campaign therefore reproduces the old, wrong answer, for the old,
wrong reason: it measures the seed.**

§5.3 is the counter-example in the other direction. And e = 1.25, which at 0.5 h reaches only `Et` 1.329 / stairs 14.486,
but at **1.5 h** reaches `Et` 0.956 / stairs **1.595** — the best point on both criteria
anyone has found. Tripling the budget changed that row by a factor of nine.

So the honest statement is not "lower `estiff` is better". It is:

> The staircase is set by **which basin you land in**, and `estiff` selects which basins
> are *available*. Going down is the cheap direction — e = 0.90 gets most of the way in
> half an hour. Going up reaches a better place, but only if you pay for the search.

A caveat on the budgets themselves: the seven 0.5 h runs were executed three-at-a-time on
a 10-core box, so each got fewer optimiser rounds than an uncontended 0.5 h would give,
whereas the 1.5 h e = 1.25 run had the machine to itself. The comparison is conservative
towards the long run, not towards the short ones.

### 5.5 An honest failure, recorded

The e = 1.40 shard **died**. `optimizeFeatures` ends by re-evaluating `state.best_params`
and refusing to write a snapshot if the cost has moved by more than 2 %; here it moved
2.4 % (1.4321 → 1.3980, i.e. *better*). That is the objective's known numerical roughness
tripping a guard that is otherwise correct. The improving round had already written
`params/ec140_opt.m`, so `CollectCampaign.m` recovers the point by re-scoring it from
disk (`Et` 1.4056 — reproducing the 1.398 the guard saw, not the 1.4321 it recorded). The
campaign script now catches that throw so one exponent cannot abort a whole shard.

## 6. The low-ATP leg

Everything above is `MgATP = 8`. The main tree has no ATP gate code at all, so nothing in
it *could* have tested 2 mM. This leg runs in the `ATP-wt-lowatp` worktree
(`Analyses/StairsPredictionLowAtp/`).

**Why it is a clean experiment.** The gates there are pure power laws
`rate_eff = rate · (MgATP/MgATP_ref)^p` with `MgATP_ref = 8`, so they equal exactly 1 at
8 mM for any `p`. Two consequences: switching them on cannot change either snapshot's
8 mM behaviour, and transplanting **one** gate vector onto **both** snapshots compares
nothing except how each parameter set *responds* to lowering ATP. The vector used is the
H4 hypothesis (`params/atphyp_H4_opt.m`): `pAtp_k2` 0.44 (slower strong-bridge
detachment), `pAtp_kah` 0.30 (slower re-priming), plus small Maxwell/series terms; 2 mM
additionally carries Pi = 0.8 mM as in every other head-to-head there.

### 6.1 The branch-drift guard — and the bug it caught in itself

The worktree is a different branch, so the snapshots are re-scored at 8 mM before
anything is believed. The **first** version of this guard reported catastrophic drift:
`Et` ≈ 12 531 for *both* snapshots. That was a bug in the guard, not the model — it ran
under the reduced `iCfg` configuration, which switches FV and the passive protocol off
while the `fn` list still demands `FV_*` and `PS_*` features, and missing features are
penalised.

Run under the configuration the snapshots were actually fitted with
(`RunDriftGuard.m`), the worktree reproduces the main tree **exactly**:

| snapshot | e | `Et` main tree | `Et` worktree | drift |
|---|---|---|---|---|
| `spreadw_opt` | 0.9653 | 0.9786 | 0.9787 | **+0.0 %** |
| `ec080_opt` | 0.80 | 1.2878 | 1.2881 | **+0.0 %** |
| `ec090_opt` | 0.90 | 0.9821 | 0.9823 | **+0.0 %** |
| `ec097_opt` | 0.9653 | 0.9572 | 0.9573 | **+0.0 %** |
| `estiff125_opt` | 1.25 | 0.9559 | 0.9548 | **−0.1 %** |

and every 8 mM slack feature matches to five significant figures (`steady` 77.243,
`ktr` 52.929, `rsK` 51.324 … identical in both trees). The ratios below are trustworthy.

### 6.2 The flip does not impair low ATP — and 2 mM breaks the tie

2 mM / 8 mM ratios, one transplanted gate vector, five snapshots. Force-class features are
normalised by each condition's own `steady` before ratioing, so `steady` carries the
amplitude and the rest carry shape.

| 2 mM / 8 mM | e = 0.80 | e = 0.90 | e = 0.965 *(refit)* | e = 0.9653 *(incumbent)* | **e = 1.25** | data |
|---|---|---|---|---|---|---|
| `steady` (force) | 1.171 | 1.176 | 1.178 | 1.179 | **1.202** | 1.376 |
| `ktr` | 0.656 | 0.652 | 0.648 | 0.650 | **0.613** | 0.540 |
| `rsK` | 0.613 | 0.629 | 0.601 | 0.566 | **0.520** | 0.345 |
| `A` | 1.025 | 1.026 | 1.025 | 1.020 | 1.016 | 0.993 |
| `peak1_y` | 1.109 | 1.111 | 1.108 | 1.108 | 1.110 | 1.039 |
| `vall_y` | 0.962 | 0.956 | 0.955 | 0.952 | 0.944 | 1.027 |
| `peak2` | 0.928 | 0.923 | 0.922 | 0.919 | 0.915 | 1.011 |

Distance from the measured ratio on the three ATP-sensitive quantities:

| |gap| | e = 0.80 | e = 0.90 | incumbent | **e = 1.25** |
|---|---|---|---|---|
| `steady` | 0.205 | 0.200 | 0.197 | **0.174** |
| `ktr` | 0.116 | 0.112 | 0.110 | **0.073** |
| `rsK` | 0.268 | 0.284 | 0.221 | **0.175** |

Two things fall out, and the second was not expected:

1. **The exponent barely touches the ATP response at all.** The four shape features are
   identical to the third digit across the whole range e = 0.80 … 1.25, and `steady` and
   `ktr` move by ~2 % and ~6 %. Low ATP does not veto any of these solutions. **The
   question this leg was run to settle is answered: the flip costs nothing at 2 mM.**
2. **But where it does discriminate, it favours e > 1 — and penalises e < 1.** On all
   three ATP-sensitive quantities e = 1.25 is the closest to the data (12 %, 34 % and
   21 % closer than the incumbent), while **e = 0.80 and e = 0.90 are *worse* than the
   incumbent on all three**, most clearly on `rsK` (0.613 / 0.629 vs 0.566, against data
   0.345). `rsK` is a ×0.3 ATP effect and the model already under-produces it; lowering
   the exponent makes that worse, raising it helps.

This is the first evidence that separates the two directions, and it does **not** come
from the staircase or from the scored objective — both of which are nearly indifferent.
§5 says e = 0.90 is the cheap way to a better staircase; §6 says e = 1.25 is the one that
also moves low ATP the right way. **Where the two disagree, 2 mM is the independent
witness.**

The residual low-ATP gap (force gain 1.17–1.20 vs 1.38, both rate slowdowns too shallow)
persists at every exponent, so it remains a separate, pre-existing defect.

---

## 7. What the data themselves say about the exponent

Model-free reads, no simulation involved. These are unaffected by everything above.

**Against e > 1 — the amplitude does not scale.** Tripling the step size should triple
the elastic peak for a Hookean bridge far from its rest strain:

- pooled across all three preparations: **peak(3 %)/peak(1 %) = 1.14**
- paired within 03/27 (step 8 vs step 7, same prep, adjacent): **1.40**

against 3.00 for a linear spring. The apparent amplitude exponent is 0.12–0.30. Forced
detachment during the 3× faster ramp pushes this ratio down too, so these are *lower*
bounds — but they are nowhere near 1.

**For e > 1 — the step-to-step growth.** On the coarse days the same ratio is data
1.48 / 1.34 versus model 0.90 / 0.88: the model's second, larger step produces a
*smaller* peak than its first. Raising e moves that ratio the right way (0.90 → 1.07 at
e = 1.15 → 1.28 at e = 1.40). Likewise the data's peak grows +5.5 % per step along the
03/27 staircase while the model grows +1.2 %.

These two are not actually in contradiction: the first is about absolute gain, the second
about shape, and `estiff` multiplies both. **Fix the gain first; the ratio only becomes
identifiable afterwards** — which is what the e = 1.25 refit does by moving `kSE` and the
shape terms along with the exponent (peak block 1.25 → 1.07 with `kSE` barely moved,
5324 → 4775).

**One dead end and one retraction:**

- **`dr` cannot substitute for the exponent.** peak/base ≈ e·Δx/x₀, so raising the
  resting strain x₀ should work; but `dr` is pinned by the isometric force (dr ×0.8 →
  46.9 kPa, ×1.4 → 96.8 kPa) and `Et` explodes to 627 / 4689. This one stands — it was
  measured against the isometric force directly, not through a fixed-`kstiff2` slice.
- **RETRACTED: "`estiff` below ~0.95 trips `doublePeak`."** The earlier note recorded a
  spurious second slack-restretch peak at cost 50–102 below e ≈ 0.95, and used it to rule
  the low band out. Under a refit it simply does not happen: `doublePeak` scores **0.0215
  at e = 0.80** and 0.0164 at e = 0.90, against **0.0226 for the incumbent itself**. The
  50–102 came from the same uncompensated slices §4.4 retracts. The low band is
  admissible, and §5 shows it is the cheap direction.

---

## 8. What to do next

**The right experiment is no longer another `estiff` study.** §4 shows why every
low-dimensional probe of this parameter misleads, §5.1 shows the exponent is real, and
§5.3 shows that ranking exponents by a fixed search budget mostly ranks seeds. The way
out is not a better sweep:

1. **Put the staircase in the objective.** It is the only thing here that actually
   discriminates between these solutions, and it has never been scored. §2 gives the
   features and the measured weights; `extractStairsAttributes.m` is written to the same
   contract as `Model/extractSlackAttributes.m` and is ready to be promoted to `Model/`.
   Score the **1 % peak, nadir and recovery** (rel. SD 21 / 26 / 17 %); never score the
   3 % nadir (217 %).
2. **Then let `estiff` float**, with `kstiff2` compulsory. A search scored on the thing
   that distinguishes the basins will move it; a search that is not, provably will not.

**Open, in order of how much they matter:**

- **Is there a third basin at e = 0.965 that no search has found?** §5.1 shows half an
  hour does not find one, and the older 1.5 h control (`estiff965`, seeded off-optimum)
  reached only `Et` 4.9955 / stairs 5.780 without re-entering the incumbent basin at all.
  Neither is proof of absence. This is the one result that would overturn §0.
- ~~How far does e = 0.90 go with a real budget?~~ **Answered in §5.3:** budget-matched to
  1.5 h it reaches `Et` 0.9580 with staircase **4.000** — a better objective than at 0.5 h
  but a *worse* staircase, and 2.5× behind e = 1.25.
- **Which of e ≈ 0.90 and e = 1.25 is physically right?** They are different mechanisms —
  softening versus stiffening — and the evidence is genuinely split: §7's model-free
  amplitude ratio points *down* (apparent exponent 0.12–0.30), the step-to-step growth
  and **§6.2's low-ATP ratios** point *up*, and both the objective and the staircase are
  nearly indifferent. §6.2 is currently the strongest tiebreaker and favours e = 1.25. A
  protocol that varies step size systematically would settle the mechanism directly.
- **The residual low-ATP gap** (force gain 1.17–1.20 vs 1.38; `ktr` 0.61–0.66 vs 0.54;
  `rsK` 0.52–0.63 vs 0.35) persists at every exponent and is therefore a separate
  problem — see `Analyses/AtpEffectHypotheses` in the `ATP-wt-lowatp` worktree.
- **Tighten the optimiser's final guard.** §5.4: a 2.4 % benign re-evaluation drift
  killed a shard. The tolerance should be set from the objective's measured roughness
  rather than a flat 2 %.

## 9. Files

| file | what |
|---|---|
| `extractStairsAttributes.m` | per-step peak / nadir / plateau / recovery, window means, one instrument for data and model |
| `stairsCost.m` | four-block, amplitude-normalised staircase cost |
| `stairsIsoRef.m` | isometric operating point + iso-force compensation handle |
| `scanLever.m`, `showBlocks.m` | one-parameter sweeps and their tables |
| `RunStairsBaseline.m` | the prediction test + noise-floor diagnostics → `baseline.png` |
| `RunDataSpread.m` | measurement precision, between-prep spread, model-free exponent |
| `RunCrossProtocol.m` | prediction of the two coarse staircases + the `dr` lever |
| **`RunEstiffWall.m`** | **§4.2/§4.3 — the wall, the valley floor, why the optimiser stalled** |
| **`RunEstiffCampaign.m`** | **§5 — matched-budget refits at seven frozen exponents** |
| **`estiffCampaignPool.m`** | the spreadw pool minus `estiff` (the frozen treatment) |
| **`CollectCampaign.m`** | re-scores every `ec*_{seed,opt}.m` from disk → `campaign_table.mat`; the authoritative §5 table, immune to shard clobbering and to shard death |
| **`MakeCampaignFigure.m`** | pure plotting: `estiff_wall.png`, `estiff_campaign.png` |
| **`../../../ATP-wt-lowatp/…/RunDriftGuard.m`** | §6.1 — does the low-ATP worktree reproduce the main tree's 8 mM model |
| `RunRepairFit.m` | the original frozen-e = 1.25 repairability test (§4.4) |
| `RunEstiffScan.m`, `RunLevers.m`, `RunGrid2D.m`, `RunEstiffProfile.m` | **retracted slices** (§4.4); kept as the record of the artefact |
| `MakeFigures.m` | `estiff_verdict.png`, `traces.png` (both show the retracted verdict) |
| `../../../ATP-wt-lowatp/Analyses/StairsPredictionLowAtp/RunLowAtpLeg.m` | §6 — the 2 mM leg |

`extractStairsAttributes.m` is written to the same contract as
`Model/extractSlackAttributes.m` and is ready to be promoted to `Model/` if the staircase
is ever wired into `RunBakersExp` as a scored experiment. The weights in §2 are what it
should be given.

**MATLAB gotcha found the hard way:** MATLAB caches parsed script files, so `loadParams`
on a snapshot a background optimiser has just rewritten silently returns the OLD one.
`clear functions; rehash` before re-reading, or you will score a stale point and believe
it.
