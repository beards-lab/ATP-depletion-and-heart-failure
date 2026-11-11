# Fitting Strategy: Getting a Plausible Fit

> **Working plan, iterative review.** Pushback inline as `> [!CAUTION]`
> for open items. Resolved comments are silently folded into the body.
>
> **Status: revision 4.** Sections reordered by implementation phase.
> Open: §1f dashboard placement; §1c fit-RMSE feature handling; §4 basin
> informativeness; §3 outcome row 4; §5 sensitivity stability; §7 decisions.

---

## 0. Diagnosis (context)

### Symptom
~16 parametrisations, none fits well. Different protocols fail on different runs (A, t0, ktr alternate). Even "good" fits give physiologically implausible parameters.

### Hypothesis
Search problem, not model problem:
1. Local-minimum trapping
2. Parameter unidentifiability (flat ridges)
3. No physiological anchoring

**Path forward**: constrain the search before enriching the model. Implementation phases below.

---

## 1. Phase 1 — Bounded fminsearch infrastructure

**Primary deliverable.** Everything else depends on this.

### 1a. Bounds source

Each parameter has a literature-justified range. Cardiac biophysics, not generic ±X%. Full table in `parameterBounds.m` (§1e).

### 1b. PCHIP convention

PCHIP central value (at strain s ≈ 0) is fixed at ~1 by convention. The base rate (`ka`, `kd`, `k1`, `k2`, ...) then equals the rate at zero strain — what literature reports — so bounds apply directly to the base rate. PCHIP shape values away from centre stay unconstrained (Tier D in §1e).

### 1c. Cost decomposition

Existing cost structure already covers feature-based and time-course terms. **Add a single new term: physiology penalty.**

```
E_total(g) = w_feat * E_features(g)        % current
           + w_tc   * E_timecourse(g)      % current
           + w_phys * E_physiology(g)      % NEW
```

If `w_phys = 0`, the term is not computed (no overhead).

**Fit-quality RMSE** is *not* a separate cost term — it's already accessible as a feature (e.g. `ktr_rmse`, `slack_fit_rmse` from `extractSlackAttributes`). Include it in `params0.fn` with target value 0 and an appropriate weight, e.g. `'ktr_rmse|0|5'`. This pushes the optimiser away from regions where exponential fitting fails.

> [!CAUTION]
> **§1c — confirm fit-RMSE features exist & their names.**
> `extractSlackAttributes` returns the features struct used in §1c. Need to verify which `*_rmse` fields are emitted and their typical magnitudes. If they don't exist, add to the extractor (one-line addition: store `gof.rmse` from `fit()`). FJ to confirm or specify the field names.

### 1d. Penalty form

```matlab
function p = softbound(x, lb, ub)
    % Quadratic exterior penalty: 0 inside [lb, ub], rises as squared
    % fractional violation outside.
    p = max(0, (lb - x) / (ub - lb))^2 ...
      + max(0, (x - ub) / (ub - lb))^2;
end
```

### 1e. `parameterBounds.m` — structure and content

Single file, organised by penalty weight. Penalty weight per parameter, no tier abstraction at the algorithm level (the `%% TIER` headings are organisational only).

Citations live as **code comments** next to each entry, never as struct fields — keeps the runtime struct lean.

```matlab
function bounds = parameterBounds()
% PARAMETERBOUNDS  Physiological bounds for ALL params in dPUdT_CombinedTransitions.
%
% Returns: bounds.(paramname).{lb, ub, weight}
% Applied universally — to every parameter, optimised or not. Wrong frozen
% values still register on the dashboard (§1f).
%
% Note on SRX nomenclature: ksr0, kmsr, sigma1, sigma2 implement the IHM
% mechanosensing ON/OFF transition (Linari 2015, Brunello 2020) — NOT the
% biochemical SRX of McNamara 2015 (~0.004 s⁻¹). Bounds reflect the
% mechanosensing interpretation. No rename — comments only.

%% TIER A — penalty weight 100
% Well-anchored from direct measurement; bound width ~ literature uncertainty.

% Power stroke distance dr [µm]. Single-molecule optical-trap measurements:
%   Guilford 1997  https://doi.org/10.1016/S0006-3495(97)78753-8  PMID:9138566
%   Mehta 1999     https://doi.org/10.1038/22146                  PMID:10440376
%   Capitanio 2006 https://doi.org/10.1073/pnas.0506830103        PMID:16371472
bounds.dr.lb = 0.008; bounds.dr.ub = 0.011; bounds.dr.weight = 100;

% Cardiac titin passive-force exponent gamma [dimensionless].
%   Granzier 1995  https://doi.org/10.1016/S0006-3495(95)80278-X  PMID:7756567
%   Wu 2000        https://doi.org/10.1161/01.RES.86.6.682        PMID:10747009
%   Terui 2008     https://doi.org/10.1085/jgp.200709864          PMID:18316560
bounds.gamma.lb = 7; bounds.gamma.ub = 10; bounds.gamma.weight = 100;

% ... L_thick, L_thin, L_hbare, d_actin

%% TIER B — penalty weight 10
% Order-of-magnitude from cardiac literature.

% Attachment rate ka [s⁻¹] (cardiac).
%   Little 2012 https://doi.org/10.1074/jbc.M112.380048 PMID:22685295
bounds.ka.lb = 50; bounds.ka.ub = 200; bounds.ka.weight = 10;

% Detachment rate kd [s⁻¹] (cardiac). Little 2012 (same ref).
bounds.kd.lb = 40; bounds.kd.ub = 150; bounds.kd.weight = 10;

% ADP-release rate k2 [s⁻¹] (cardiac).
%   Malmqvist 2004    https://doi.org/10.1074/jbc.M310974200 PMID:14976210
%   Siemankowski 1985 https://doi.org/10.1073/pnas.82.3.658  PMID:3156150
bounds.k2.lb = 100; bounds.k2.ub = 400; bounds.k2.weight = 10;

% ... kah, kamh, k1, k_1, kSE, kstiff1, kstiff2-ratio, sigma1/2

%% TIER C — penalty weight 1
% Defensible-but-uncertain (mechanosensing IHM rates, viscosity).

% IHM mechanosensing OFF rate ksr0 [s⁻¹] (NOT biochemical SRX).
%   Linari 2015   https://doi.org/10.1038/nature15727 PMID:26560032
%   Brunello 2020 https://doi.org/10.1073/pnas.1920632117 PMID:32193335
bounds.ksr0.lb = 0.5; bounds.ksr0.ub = 30; bounds.ksr0.weight = 1;

% ... kmsr, ksr2srd, ksrd2sr, mu, mu_neg

%% TIER D — no penalty (weight 0)
% PCHIP shape values away from centre; A2-hopping rates — no cardiac
% measurement available. Listed for coverage; weight = 0 means no penalty.
% bounds.PieceWiseStrainDepParams__2.weight = 0;
% bounds.slope.weight = 0;

end
```

**Coverage requirement**: every parameter that appears in `dPUdT_CombinedTransitions.m` must be listed. Use `Docs/parameter-reference.md` as the master checklist. Tier-D entries are intentional acknowledgements of "no prior available", not omissions.

### 1f. Bound dashboard — integration with `plotFeatures`

> [!CAUTION]
> **§1f open — dashboard placement.**
>
> | Option | Pros | Cons |
> |---|---|---|
> | **A: Extra panel in `figure(80085)`**, `plotFeatures` calls dashboard with bounds struct | Single visual artefact per fit; existing `BatchRunAllParams` workflow captures it | Couples plotting to bounds-loading; extra arg to plotFeatures |
> | **B: Bound violations stored as fields in `features_model.physiology`** | Clean architectural fit; `evalFeatureCost` could compute physiology penalty as a "feature" | Mixes apples and oranges; normalisation logic doesn't apply |
> | **C: Separate `plotPhysiologyDashboard`** function called next to `plotFeatures` | Clean separation; doesn't disturb feature workflow | Extra plot to track; not auto-included in batch |
>
> **Recommendation: Option A.** Bounds struct loaded at the top of `plotFeatures` via `parameterBounds()`. If absent, panel is omitted (legacy callers unbroken). FJ to confirm.

### 1g. Universal application

Penalties apply to every parameter in `parameterBounds.m`, regardless of `mods`. Reasoning: the rules are universal. A wrong frozen value (e.g. `kstiff1` left from a stale param file) should still light up the dashboard — the optimiser can't change it, but the user gets alerted that the baseline is unphysical.

This also resolves the "data/physiology balance" issue — physiology cost is a fixed scalar that doesn't depend on `numel(mods)`.

### 1h. Implementation hooks

- **New file**: `Model/parameterBounds.m`
- **New helper**: `Model/evalPhysiologyCost.m` — takes `params` struct + bounds, returns scalar `E_phys`
- **Modify** `Model/RunBakersExp.m` (or `runSlackExperiment.m`):
  - At end of cost calculation, if `params0.w_phys > 0`, call `evalPhysiologyCost` and add to `E`
- **Modify** `Auxiliary/plotFeatures.m`: optionally call dashboard helper if bounds available
- **No change** to `OptimizeMechanismEvaluation.m` — cost picks up the new term automatically

---

## 2. Phase 2 — Decisive multi-start of current best

**Goal**: distinguish search problem from model problem before sinking time into full multi-start.

### 2a. Procedure
Once Phase 1 is in place:
1. Take the current best parametrisation (e.g. `ModelOptParams_SRXTD_iter_9`)
2. Run 10× fminsearch restarts from small random perturbations (~5% of each param)
3. Evaluate residual spread per feature

### 2b. Outcome interpretation

| Spread of feature residuals | Diagnosis | Next |
|---|---|---|
| Consistent across restarts (low variance) | Mechanism issue — search isn't the problem | Pause Phases 3-5; revisit mechanisms (§4 mechanism A/B) |
| Scattered (high variance) | Search issue confirmed | Proceed to Phase 3 |

### 2c. Compute (FJ's measured numbers)
- `RunBakersExp` parallel ≈ 25 s (slack 20 s + FV 5 s)
- One restart ≈ 2000 evals × 25 s = ~14 hr parallel
- 10 restarts × 14 hr = ~140 hr ≈ 6 days (overnight × 6, or weekend run)
- With MEX (3.6× speedup): ~40 hr ≈ 1.7 days

---

## 3. Phase 3 — Full multi-start (LHS)

**Goal**: explore the full physiology-bounded parameter region.

### 3a. Procedure
1. Build LHS sample of N=50 points within the Tier A/B bounded region
2. Cheap prefilter: evaluate each at one cost call (~25 s × 50 = ~21 min)
3. Top 10 by cost → full fminsearch (~14 hr each parallel; ~140 hr total)
4. Cluster results in parameter space; save top 3 basins as separate param files

### 3b. LHS coverage caveat
50 samples in 20-D ≈ 1.2 points/dim — sparse. Mitigations:
- Run multi-start over **active mods only** (10–15 dims, not 20)
- Bounds-aware LHS (samples land inside Tier A/B bounded region by construction)

### 3c. Outcome interpretation

> [!CAUTION]
> **§3c open — fourth row needed.**
>
> | Outcome | Diagnosis | Action |
> |---|---|---|
> | All cluster, low cost, **physiology OK** | Search problem solved | Lock and verify |
> | 3–5 distinct basins, moderate cost | Genuine multimodal | Pick most plausible, do Phase 4 |
> | All scattered, high cost | Model structurally wrong | Mechanism A/B (Phase 4) |
> | **All cluster, low cost, but in a region violating Tier A bounds** | Model has strong unphysical attractor | Mechanism investigation, NOT more search |
>
> **FJ to confirm row 4 makes sense.**

### 3d. Compute
- Without MEX: ~6 days
- With MEX: ~1.7 days

---

## 4. Phase 4 — Identifiability + mechanism A/B

**Prerequisite**: Phase 3 has produced ≥ 1 defensible basin (low cost AND most Tier A bounds satisfied).

### 4a. Identifiability analysis
At each candidate basin:
1. Compute Jacobian `J = ∂features/∂params` via `Auxiliary/ResidualAndJacobian.m`
2. SVD: `J = U Σ V'`
3. Sort singular values; cut at elbow
4. Right-singular-vectors of small `σ_i` reveal unidentifiable parameter combinations

**Hard ceiling**: with ~17 features, rank(J) ≤ 17. Active `mods > 15` guarantees structural unidentifiability. Use this as upper limit on optimisation set size.

> [!CAUTION]
> **§4 open — basin informativeness criterion.**
> Need a falsifiable rule for "is this basin good enough that the local Hessian is informative?". Proposal: `E_features < threshold` AND `≥ 50% of Tier A bounds satisfied`. Otherwise the Hessian reflects the bad basin, not the model structure. **FJ to confirm thresholds.**

### 4b. Mechanism A/B
For each questionable mechanism (negative-kstiff, A2-hopping, SRD, target zone saturation, vernier):
```
E_with_M     = current cost with M enabled
Disable M, re-optimise (50-iter cap, fminsearch local)
E_without_M  = best cost without M

if E_without_M  ≤ E_with_M + ε   → M is doing nothing, permanently disable
elif E_without_M >= 1.5 * E_with_M → M is load-bearing, keep
else                              → M is load-shifting (compensating elsewhere); flag
```

### 4c. Output
- Drop unused mechanisms
- Lock unidentifiable params at literature defaults (remove from `mods`)
- Reduce `mods` to ≤ 15 informative parameters

---

## 5. Phase 5 — Feature-pinning (final polish)

**Goal**: fine-tune residual misfits without disturbing the bulk fit.

### 5a. Feature ↔ parameter map (revised per FJ)

| Feature | Constrains | Notes |
|---|---|---|
| `t0` | One-sided XB transition rate; defines Vmax | Tightly coupled with FV; calibrate together with FV, not slack |
| `ktr` rate | Sum of forward + backward XB rates | Single ktr → 1 constraint, not enough alone |
| `A` (steady amplitude) | Max attached fraction × stiffness × overlap | Magnitude product, not individual factors |
| `peak1_dSL` (slack overshoot) | Detachment kinetics + serial elasticity | k2, R2 PCHIP, kSE |
| `restretchSlopeStart` | **Serial elastic stiffness `kSE`** | SE behaviour during fast restretch |
| `vall_y`, `vall2_dy` | Slack zero-load force | XB detachment kinetics |

### 5b. How to lock
Drop the param from `mods` — it stays at its current value. No tight bounds, no hardcoding.

### 5c. Sensitivity-stability problem

> [!CAUTION]
> **§5c open — sensitivities aren't stable across parameter regions.**
> A sensitivity matrix at point A may rank features differently than at point B. **FJ asks: how to battle this?**
>
> **Proposed approach** (needs FJ buy-in):
> 1. Compute Jacobian at the top-3 multi-start basins (Phase 3 output)
> 2. Lock only params whose dominant feature *agrees across all 3 basins*
> 3. Disagreements → leave the param free
>
> If multi-start collapsed to one basin, the question is moot — sensitivities at that single point are the answer. **FJ to confirm.**

---

## 6. Recommended sequence

| Phase | Action | Compute |
|---|---|---|
| 1 | Build bounds infrastructure (§1) | 3–5 days dev, no compute |
| 2 | 10× restart of current best (§2) | 1.7 days (MEX) / 6 days |
| **Gate** | Diagnose: search vs model problem | — |
| 3 | Full LHS multi-start (§3) | 1.7 days (MEX) / 6 days |
| 4 | Identifiability + mechanism A/B (§4) | ~2 days |
| 5 | Feature-pinning final polish (§5) | open-ended |

**Total**: ~2 weeks (MEX), ~5 weeks (no MEX). Phase 2 is a hard gate — if it shows model issues, abort search-side phases and revisit mechanisms.

---

## 7. Open decisions for FJ

1. **§1c** — confirm fit-RMSE feature names from `extractSlackAttributes`. Do `*_rmse` fields exist? What weight to use in `params0.fn`?
2. **§1f** — Option A (panel in 80085) for dashboard placement?
3. **§3c** — fourth outcome row (cluster but bounds-violating) makes sense?
4. **§4** — basin informativeness criterion (Tier A satisfaction threshold)?
5. **§5c** — sensitivity-stability strategy (lock only consensus across basins)?
6. **MEX freshness** — is current MEX consistent with Maxwell-dashpot logic? Multi-start budget assumes yes.
7. **Coverage scope** — confirm `parameterBounds.m` covers every param in `dPUdT_CombinedTransitions.m`, even Tier D entries listed for documentation.

---

## 8. Things explicitly NOT to do

- Don't add new mechanisms before Phase 2 (search vs model gate)
- Don't switch to fmincon
- Don't refactor parameter naming (SRX → mechanosensing rename) — comment in `parameterBounds.m` is sufficient
- Don't run identifiability (Phase 4) until Phase 2 or Phase 3 has produced a defensible basin
- Don't trust Tier-A-violating fits even if `E_total` is low — read the dashboard
