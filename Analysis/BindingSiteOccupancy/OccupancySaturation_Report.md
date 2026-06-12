# Binding-Site Occupancy / Attachment Saturation in the Cardiac Cross-Bridge Model

**Audit, reframe, implementation, and fit-potential evaluation**
Repo: `ATP-depletion-and-heart-failure` · Model: `dPUdT_CombinedTransitions` · Date: 2026-06-11

---

## 1. Motivation

Physically, the probability of a *new* myosin head attaching falls as more heads
are already attached and more actin binding sites are occupied. Three effects
contribute: a finite number of actin binding sites per unit length (target zones
of ~2–4 actin subunits between troponins), steric crowding around bound heads,
and competition for the same site. The question posed: **does the model represent
this, and if not, how should it — via optional `params0` flags — and would it
improve the general fit?**

The short answer changed twice during investigation:

1. A per-strain-bin "occupancy" mechanism (`UseTargetZoneSaturation`) **already
   existed** but is mis-grounded.
2. Crucially, **the model's strain axis is not a binding-site axis**, so any
   *local* (per-bin) occupancy term is a category error. In this model class the
   only faithful occupancy variable is the **global** bound fraction.

This report documents the physics, the code audit, the implemented fix, the
verification, and an honest assessment of fit-improvement potential.

---

## 2. Physical mechanism and literature

### 2.1 Why attachment saturates
- **Finite sites / target zones.** Strong-binding attachments are confined to
  discrete target zones (~2–4 actin subunits midway between troponins). The
  myosin (42.9 nm) / actin (~37 nm) periodicity mismatch limits how many heads
  can reach a free site.
- **Steric crowding & competition.** Bound heads sterically hinder neighbours;
  multiple crowns compete for a small pool of exposed sites. In vitro motility
  velocity peaks when *actin sites saturate* before stall — i.e. site
  availability becomes rate-limiting (**Stewart, Murthy, Dugan & Baker 2021,
  JBC**; K_M ≈ 16–26 myosin heads per µm actin). *This paper is from the same
  laboratory that produced the data this model is fit to.*
- **Cooperativity caveat (opposing sign).** Thin-filament activation is
  *positively* cooperative (bound heads open neighbouring tropomyosin), Hill
  n≈2; site saturation is *negatively* cooperative. The two partly offset. The
  present model has **no Ca/thin-filament activation term** (`UseCa=false`), so
  there is no positive-cooperativity term to double-count against (see §5).

### 2.2 How model classes handle it
- **Mean-field / strain-distribution models** (Huxley 1957; Razumova 1999; Rice
  et al. 2008; Tewari/Beard 2016; ToR-Land 2017): the attachment "saturation" is
  only the `(1−n)` probability-normalisation of Huxley — **not** a biological
  site term. These models **cannot** represent local/neighbour site competition:
  the mean-field assumption *is* actin positional symmetry, which erases site
  identity (**Smith & Mijailovich 2008**, Ann Biomed Eng). Rice et al. (2008,
  PSB) further warn that *global feedback onto Ca-affinity* produces spurious
  hysteresis — a hazard we avoid (§5).
- **Spatially-explicit models** (Daniel/Chase 1998; **Tanner/Daniel/Regnier
  2007**, PLoS CB — "ODE models cannot capture XB-induced XB recruitment";
  **Mijailovich MUSICO 2016**, JGP — "mass-action models neglect discrete
  positions and site-number constraints"; **FiberSim 2022**): these track each
  actin site and enforce per-site steric exclusion. **True local site
  competition lives here, and is out-of-class for the present model.**

### 2.3 How occupied does cardiac muscle actually get?
The commonly cited **~30–43% bound fraction is frog/rabbit *skeletal*** muscle
(**Linari et al. 1998**, stiffness) — not applicable here. The data in this
project are **rat RV cardiac trabeculae** (Baker lab, UCSF;
`Docs/2026 Myofilament-meeting-abstract.md`), and rat ventricle is **α-MHC
dominant** (shifts toward β-MHC under failure). For **rat cardiac at maximal-Ca
isometric, direct trabecula measurement gives ≈12% of motors bound** (range
~8–20%; α-/β-MHC duty ratios ~10–16% concur). SRX sequestration and submaximal
Ca push the realised value lower.

> **Physical ceiling adopted: P_bound_max ≈ 0.12–0.15** (a hard upper bound on
> the total bound fraction at maximal-Ca isometric).

---

## 3. Code audit

### 3.1 Where attachment happens
`RD1 = ka_eff*PD*N_overlap*f_lattice` ([dPUdT_CombinedTransitions.m:309](../../Model/dPUdT_CombinedTransitions.m)).
`RD1` is the total attachment flux, then spread spatially into `p1` near zero
strain. The detached pool `PD`/`PT` is conserved
(`PT = 1-NP-(p1_0+p2_0+p3_0+PD+P_SR+P_SRD)`, :180).

### 3.2 The pre-existing per-bin mechanism is mis-specified
`UseTargetZoneSaturation` multiplied the attachment flux *per strain bin* by
`f_sat = max(0, 1 - (p1+p2+p3)*dS / max_attached_per_bin)` (default off, zero
optimisation weight). Three defects:

- **D1 — heads vs sites, and both matter.** It counts attached *heads*; the
  physical limiter is occupied *sites*. The head pool is *already* self-limited
  in the model — but only as a **mean-field ensemble** (`PD`/`PT` depletion), with
  no strain- or neighbour-specificity. The *site* pool is genuinely missing.
- **D2 — `max_attached_per_bin` ungrounded & wrong-tissue ceiling.** It floated
  with a placeholder default; and the literature number usually reached for
  (~30–43%) is skeletal, not the rat-cardiac ~12–15% (§2.3).
- **D3 (the decisive one) — the strain axis is not a site axis.** `s` is
  cross-bridge **strain** (displacement from the attachment point; see
  :146–154 *"strain — above the myosin heads is zero; negative = shorter,
  positive = longer"*). There is **no spatial coordinate along the thin filament
  and no site index anywhere**. Two heads in the same `s`-bin merely share a
  *strain value*; they are **not** at the same actin site and do not compete.
  Heads competing for one site generally have *different* strains. So a
  per-strain-bin factor penalises "heads sharing a strain," which is **meaningless
  as site competition** — a category error, deeper than the dS issue.
- **D3b — and it is grid-dependent.** Because `p1,p2,p3` are densities, the old
  `(p1+p2+p3)*dS` is per-bin *mass*, which shrinks as `dS` refines. Even after
  rewriting to a density-vs-density form, point attachment makes the injected
  density diverge as `dS→0`, so **any pointwise occupancy is ill-posed** (shown
  empirically in §6).

### 3.3 The correct variable here is the global integral
In a strain-distribution model the only faithful occupancy proxy is the total
bound fraction `P_bound = p1_0+p2_0+p3_0 = ∫(p1+p2+p3)dx`. A global factor
`(1 − P_bound/P_bound_max)` on the attachment flux is a recognised mean-field
approximation, with direct precedent in **Stewart/Baker 2021** (§2.1).

---

## 4. Implementation (optional flags, all default off)

All changes are behaviour-preserving when the flags are off, so existing fits and
saved param snapshots are unaffected.

**`Model/getParams.m`** — new defaults:

| Param | Default | Meaning |
|---|---|---|
| `UseGlobalOccupancySaturation` | `false` | enable global occupancy attenuation of attachment |
| `P_bound_max` | `0.15` | ceiling on total bound fraction (rat-cardiac isometric) |
| `OccupancyForm` | `'linear'` | `'linear'` → `1−P_bound/P_bound_max`; `'langmuir'` → `1/(1+P_bound/P_bound_max)` |
| `UseTargetZoneSaturation` | `false` | **[DEPRECATED as occupancy]** per-strain-bin modifier |
| `rho_attach_max` | `16.0` | density cap `[1/µm]` for the (deprecated) dS-invariant per-bin form |

**`Model/dPUdT_CombinedTransitions.m`** — the global term, applied to `RD1`
immediately after it is computed (before spatial spreading):

```matlab
if params.UseGlobalOccupancySaturation
    P_bound = p1_0 + p2_0 + p3_0;
    if strcmpi(params.OccupancyForm, 'langmuir')
        f_sat_global = 1 / (1 + P_bound / params.P_bound_max);
    else
        f_sat_global = max(0, 1 - P_bound / params.P_bound_max);
    end
    RD1 = RD1 * f_sat_global;
end
```

The deprecated per-bin form was relabelled and made dS-invariant (density vs
density: `f_sat = max(0, 1 - (p1+p2+p3)/rho_attach_max)`), but see §6 — it remains
ill-posed and should not be used as a site-occupancy mechanism.

**`Model/parameterBounds.m`** — added `bounds.P_bound_max` (lb 0.02, ub 0.20 =
rat-cardiac ceiling) and `bounds.rho_attach_max`; the legacy
`max_attached_per_bin` bound is marked deprecated.

---

## 5. Avoiding double-counting

The known hazard (Rice et al. 2008 PSB: feeding total bound fraction back onto
**Ca-affinity** causes spurious hysteresis) does not arise here: the model runs
`UseCa=false` with no permissive-fraction / thin-filament activation term. The
global occupancy factor is therefore the clean, sole home for site saturation and
does not collide with cooperativity. If a Ca-activation term is added later, keep
the two orthogonal: permissive-fraction = *chemically available* sites
(cooperative, Ca-driven); `(1−P_bound/P_bound_max)` = *physically unoccupied*
sites (anti-cooperative, occupancy-driven).

---

## 6. Verification (empirical)

Run `Analysis/BindingSiteOccupancy/TestOccupancySaturation.m` (isometric,
maximal-Ca, unfitted defaults; for *mechanism* verification, not absolute force).

**A/B — the global term clamps the bound fraction to the ceiling:**

| Condition | Force | P_bound |
|---|---:|---:|
| baseline (off) | 6.16 | **0.534** |
| global **linear**, cap 0.15 | 1.61 | **0.139** (≈ ceiling) |
| global **langmuir**, cap 0.15 | 4.07 | 0.352 (half-sat) |

The unconstrained model predicts P_bound ≈ 0.53 — **~4× the rat-cardiac physical
ceiling** — which the linear form corrects to ≈0.14. (Langmuir uses the cap as a
half-saturation, so it clamps more gently.)

**Sweep `P_bound_max` (linear)** — bound fraction tracks the cap monotonically:

| `P_bound_max` | 0.05 | 0.10 | 0.15 | 0.20 | 0.30 |
|---|---:|---:|---:|---:|---:|
| resulting P_bound | 0.049 | 0.095 | 0.139 | 0.180 | 0.251 |
| Force | 0.57 | 1.10 | 1.61 | 2.08 | 2.90 |

**dS-invariance — the global form is grid-robust; the per-bin form is not:**

| Form | dS=0.004 (ss=82) | dS=0.002 (ss=163) |
|---|---:|---:|
| **GLOBAL** P_bound | **0.13921** | **0.13921** |
| per-bin P_bound | 0.0995 | 0.0560 |

The global form is identical to 5 digits across a 2× grid refinement (Force 1.606
vs 1.609, within integration tolerance). The per-bin form changes by ~45% — a
direct demonstration that pointwise occupancy is ill-posed in this model
(point attachment makes the injected density diverge as `dS→0`), confirming the
§3 reframe.

Static analysis (`check_matlab_code`) on all edited/created files is clean
(only pre-existing style warnings remain).

---

## 7. Full-fit A/B and FV investigation (BoundedFit basis)

Beyond the isometric mechanism check, the global term was run against the actual
multi-experiment cost using the **BoundedFit basis** (`Workbench/DriverBoundedFit.m`:
iter_17 + `params_reseeded` + SRX overlay; FV + slack; cost = `boundedOutputFn`).
Scripts: `RunBoundedFitOccupancyAB.m`, `RunFVProbe.m`. Cost = weighted
output-feature cost (sum of `plotFeatures` `cost_vec`); lower is better.

### 7.1 Occupancy alone (no force recompensation) — misleading

| config | feat cost | XTOR [3,10] | attached_ss | FV_fnorm cost |
|---|---:|---:|---:|---:|
| baseline (off) | **26.0** | 12.9 (OUT) | 0.253 | 15.98 |
| linear, cap 0.30 | 81.6 | **7.99 (IN)** | 0.158 | ~16 |
| langmuir, cap 0.30 | 47.9 | 9.41 | 0.188 | 15.7 |

Occupancy *fixes the XTOR over-bound* (12.9→8.0) but inflates total cost: it
uniformly lowers force and wrecks the *absolute*-force slack features (A, steady,
peak — each weight 50). `FV_fnorm` is normalised (scale-free) and barely moves, so
the blow-up is an **artifact of not recompensating stiffness**, not a regression.

### 7.2 Occupancy at matched force + FV shaping — substantial improvement

Recompensating with `kstiff2 ×1.6` (restores the ~0.62× force drop) and probing FV
knobs:

| config | feat cost | FV cost | XTOR | A | steady | ktr |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 26.0 | 15.98 | 12.9 | 64.0 | 79.9 | 62.4 |
| R2 reshape only | 21.8 | 11.99 | 12.9 | 64.1 | 79.9 | 63.5 |
| occ 0.30 + kstiff×1.6 | 21.8 | 12.10 | 10.1 | 64.5 | 78.6 | 75.5 |
| **occ + kstiff + R2 reshape** | **17.95** | **8.50** | **10.1** | 64.8 | 78.6 | 77.6 |

Normalised FV vs data `[1 .92 .66 .32 .20 .11]`: baseline `[1 .59 .33 .14 .08 .05]`
→ combined `[1 .64 .43 .23 .15 .10]` — **~half the FV error removed**, force level
preserved, XTOR brought to the bound edge. Total cost **26→18 (−31%)**, FV cost
**16→8.5 (−47%)**.

### 7.3 What actually shapes FV (probe results)

- **`mu` (viscous drag) is a non-lever here:** 0.015→0.008 moved FV cost
  15.98→15.94; lowering further (≤0.003) destabilises the ODE. FV steepness is
  *kinetic*, not mechanical.
- **Reducing R2 at negative strain *alone makes FV worse*** (cost 26→46): trapped
  negative-strain heads create drag and collapse Vmax. The first intuition was
  backwards.
- **The real R2 lever is detaching *earlier*** — raise the R2 multiplier at the
  onset of shortening (`PieceWiseStrainDep2Params` knot at s2=0, 1.05→5) while
  trimming the deep-negative burst, so heads release before being dragged into
  negative-force territory. ~25% FV-cost cut on its own.
- **Occupancy and the R2-onset reshape are complementary and stack** (combined FV
  cost 8.50 vs 11.99 / 12.10 singly): occupancy preferentially suppresses the
  high-force isometric states (flattening the curve at high velocity, fixing
  XTOR); the R2 reshape repairs the steep low-velocity drop and raises Vmax.

---

## 8. Verdict

- **Occupancy earns its place — at matched force.** With `kstiff2` recompensation
  it cuts total feature cost (26→22 alone, 26→18 with FV shaping), fixes the XTOR
  over-bound, improves normalised FV, and leaves absolute force essentially
  unchanged. *Without* recompensation it looks harmful — always co-vary `kstiff2`
  (or `ka`) when enabling it.
- **For the FV fit specifically**, the largest, cheapest win is the **R2
  onset-detachment reshape** (kinetic); occupancy is a complementary second lever.
  `mu` is not useful (and unsafe to lower). FV improvements carry into the slack
  features (shared kinetics), consistent with the user's note that FV is "key for
  slack as well."
- **Settings:** `UseGlobalOccupancySaturation=true`, `OccupancyForm='linear'`,
  `P_bound_max` an identified parameter in `[0.02, 0.40]` (default **relaxed to
  0.30** given uncertainty in the rat-cardiac ceiling; 0.12–0.20 is the
  literature-tight range). Re-fit `kstiff2` alongside.
- **Residual trade-off / next frontier:** `ktr` rises (62→78 vs data ~50) when
  occupancy is on — consistent with the [[mech_tradeoff]] iron law (ktr tied to
  k2eff). Closing FV *and* ktr together likely needs the 3rd low-force
  slow-detaching state, not more 2-state tuning.
- **Do not** use `UseTargetZoneSaturation` as occupancy (deprecated, §3). Local
  site competition is out-of-class — spatially-explicit models only (MUSICO,
  FiberSim, Tanner/Daniel/Regnier).

---

## 9. Files

- `Model/getParams.m` — new params (§4).
- `Model/dPUdT_CombinedTransitions.m` — global occupancy on `RD1`; per-bin form
  relabelled/made dS-invariant.
- `Model/parameterBounds.m` — bounds for `P_bound_max` (lb 0.02, ub 0.40),
  `rho_attach_max`.
- `Analysis/BindingSiteOccupancy/TestOccupancySaturation.m` — isometric mechanism
  checks (A/B, sweep, dS-invariance).
- `Analysis/BindingSiteOccupancy/RunBoundedFitOccupancyAB.m` — full-fit occupancy
  A/B against the BoundedFit basis (§7.1).
- `Analysis/BindingSiteOccupancy/RunFVProbe.m` — FV-shaping probe: occupancy +
  kstiff, R2 reshape, mu (§7.2–7.3). Saves `fv_probe_results.mat`.

## Key references

Stewart, Murthy, Dugan & Baker 2021 (JBC) · Mijailovich et al. 2016 (JGP, MUSICO)
· Smith & Mijailovich 2008 (Ann Biomed Eng) · Tanner, Daniel & Regnier 2007
(PLoS Comput Biol) · Rice, Wang, Bers & de Tombe 2008 (Biophys J) · Rice, Tu,
Poggesi & de Tombe 2008 (Pac Symp Biocomput) · Linari et al. 1998 (Biophys J) ·
Huxley 1957 (Prog Biophys) · Razumova, Bukatina & Campbell 1999 (J Appl Physiol)
· Tewari et al. 2016 (J Mol Cell Cardiol) · Land et al. 2017 (J Mol Cell Cardiol)
· Campbell et al. 2022 (FiberSim, Biophys J).
