# The Attached-Pool / ATPase / k_tr / V_max Paradox

*Generated 2026-06-04. Companion to `mechanism-evaluation.md` and `low-ATP-force-enhancement-analysis.md`. Model measurements from `ModelOptParams/ModelParamsSlackKtrOpt.m`; diagnostic tools in `Workbench/{mechProbe,isoProbe,ktrProbe}.m`. Revised after temperature/isoform scrutiny.*

---

## 1. Scope

At **8 mM ATP, saturating Ca²⁺, isometric** (the "baseline"), reconcile four targets for **mouse α-MHC** myocardium:

| Target | Originally wanted | Physiologically corrected (see §4) | Current model |
|---|---|---|---|
| Attached pool (p1+p2) | > 0.30–0.40 | 0.30–0.40 | 0.078 |
| ATP turnover (XTOR) | < 5 s⁻¹ | **≈ 8–12 s⁻¹** | 4.96 |
| k_tr | ≈ 50 s⁻¹ (mouse, 20 °C) | ≈ 50 s⁻¹ | 51.9 |
| V_max | hold (~1.5 ML/s) | hold | 1.51 |

The headline result: **the turnover target was the misset one.** For mouse α-MHC at k_tr = 50 with a duty ratio ≥ 0.3, the ATP turnover is *necessarily* ~10–12 s⁻¹, not < 5 — and this follows from a temperature-free identity, not from temperature hand-waving. Once corrected, the four are nearly consistent and the residual gap is a concrete, fixable model defect (a throttled hydrolysis step), not a deep paradox.

---

## 2. The temperature-free identity that couples the three isometric observables

Two-state cycle, detached **D** ⇌ attached **A**, attachment rate `f` (D→A), detachment rate `g` (A→D). At steady state:

```
duty        r   = f/(f+g)            (= fraction of cycling heads attached)
k_tr            = f + g              (Brenner 1988; force redevelopment)
ATPase/head J   = g·r = f·g/(f+g)    (one ATP per completed cycle; flux = g × attached)
```

Eliminating `f, g` (using `g = k_tr(1−r)`, `f = k_tr·r`):

```
┌─────────────────────────────────────────────┐
│   ATPase_per_head  =  r · (1 − r) · k_tr      │
└─────────────────────────────────────────────┘
```

**This is exact, parameter-free, and temperature-free.** It checks against measurement to 2 %:
- rat α-MHC, 20 °C (PMC2711735): r = 0.29, k_tr = 23.9 → 0.29·0.71·23.9 = **4.9** s⁻¹; measured ATPase ≈ 4.8 s⁻¹. ✓

Because `r(1−r) ≤ 0.25` (maximum at r = 0.5) and is **flat** over the plausible range — `r ∈ [0.3, 0.7] ⟹ r(1−r) ∈ [0.21, 0.25]` — the identity pins the turnover once any *two* of {duty, k_tr, ATPase} are fixed:

```
attached r ≥ 0.30  AND  k_tr = 50   ⟹   ATPase = r(1−r)·50 ≥ 10.5 s⁻¹   (NOT < 5)
attached r = 0.30  AND  ATPase = 5  ⟹   k_tr   = 5 / (0.3·0.7) ≈ 24 s⁻¹  (NOT 50)
```

There are three observables and two rate constants. **You can set any two; the third is determined.** The wanted triple (≥0.3, <5, ≈50) is over-determined by ~2× in ATPase.

---

## 3. Why this is temperature-robust (the part I got hand-wavy before)

The identity only helps across conditions if `r(1−r)` and the ATPase/k_tr ratio don't drift with temperature. **De Tombe et al. 2007** (rat myocardium, skinned, saturating Ca²⁺, *direct* 15–25 °C series) settles this:

| Quantity | 15 °C | 20 °C | 25 °C | Q₁₀ |
|---|---|---|---|---|
| k_tr (max) | 14.5 | **31.4** | 49.4 s⁻¹ | **3.3** |
| Ca²⁺-activated ATPase | — | — | — | **3.4** |
| Max tension | — | — | — | 1.4 |
| Tension cost (∝ g_app) | — | — | — | 2.5 |
| **Fraction attached** | 0.61 | 0.67 | 0.59 | **1.1 (≈ none)** |

Two facts make the identity temperature-safe:
1. **k_tr and ATPase share the same Q₁₀ (~3.3–3.4)** → their *ratio* (= r(1−r)) is temperature-invariant. So I am **not** scaling a rate by an assumed Q₁₀; I am using a measured temperature-*invariant* ratio.
2. **The attached fraction is temperature-independent** (Q₁₀ ≈ 1.1) → r is a structural property, not a thermal one.

(Strictly, tension cost has Q₁₀ 2.5 < 3.3, implying `f` is slightly more temperature-sensitive than `g`, so `r` creeps up mildly with temperature — consistent with the small 0.61→0.67 rise. This does not change the conclusion: `r(1−r)` stays in [0.21, 0.25].)

**Crucially, the user's k_tr = 50 is mouse at 20 °C, and needs no temperature scaling at all.** De Tombe's *rat* gives k_tr = 31.4 at 20 °C; mouse α-MHC is faster than rat (≈1.5–2×; mouse ≈ 1.5× rat, and α/β detach differ ~4×). 31.4 × 1.6 ≈ 50. **The mouse↔rat difference is species (isoform speed), not temperature.** My earlier "scale to warm mouse" framing was wrong in mechanism even though the number was right — the correct statement is species-based and the identity is temperature-free.

---

## 4. Mouse α-MHC turnover — the corrected number

Putting the pieces together at the **measured** operating point (mouse α, 20 °C, k_tr = 50, duty ≈ 0.29–0.3 from rat α and from the model fits):

```
ATPase_per_head  =  r(1−r)·k_tr  ≈  0.29 · 0.71 · 50  ≈  10.3 s⁻¹
```

So **mouse α-MHC isometric turnover is ~10–12 s⁻¹ per cycling head, not < 5.** The original < 5 figure corresponds to a **slower isoform / cooler / different denominator**:
- rat or human **β-MHC** (V3): duty 0.22, k_tr 11–25, ATPase ≈ **2–5 s⁻¹** (PMC2711735; this is where "< 5" lives);
- rat **α** at 20 °C: ATPase ≈ 5 s⁻¹ **but** k_tr only 24 — *not* 50;
- a **per-total-myosin** figure after SRX dilution (§7.4) — different denominator.

Independent solution data agree the heart is not that economical per cycling head when fast: actin-activated steady-state kcat ≈ 2–10 s⁻¹ for cardiac S1/HMM at low temperature, scaling up with temperature/α-content. **Bottom line: the corrected mouse-α target is (attached ≈ 0.3, k_tr ≈ 50, ATPase ≈ 10), which obeys the identity — there is no contradiction, just a number that was off by ~2×.**

---

## 5. The α-myosin hydrolysis step — and the model's throttle

The user asked specifically about the **hydrolysis step** (M·ATP → M·ADP·Pi). Literature (Dictyostelium and general myosin-II): forward rate **≈ 100–160 s⁻¹**, equilibrium constant ~1–13, i.e. a **fast, near-equilibrium** step — *not* rate-limiting. The cycle is rate-limited downstream (Pi release coupled to the power stroke, or ADP release / ATP-gated detachment), with steady-state kcat of only a few s⁻¹. The fast α isoform's hydrolysis is at least this fast.

**The model violates this.** Its `kah` ("ATP hydrolysis state change") is `kah = 18.8` × `xrate 0.435` = **8.2 s⁻¹ effective** — ~12–20× *slower* than the real hydrolysis step. Consequence: heads pile up pre-hydrolysis. The live baseline has **PT = 0.61, PD = 0.31 — 91 % of heads detached**, attached only 0.078, *because* the throttled `kah` makes hydrolysis the rate-limiting step (it should be fast). `mechanism-evaluation.md` already flagged the post-`xrate` effective rates as sub-literature; this is the most damaging instance.

**Model demonstration (this analysis).** Raising `kah` toward physiological values drains PT and lifts the attached pool exactly as predicted:

| Config | attached | XTOR | k_tr | PT | comment |
|---|---|---|---|---|---|
| baseline (`kah_eff`=8) | 0.078 | 5.0 | 52 | 0.61 | PT hoards 61 % of heads |
| `kah_eff`≈82, `k2`×0.4 | 0.30 | 8.3 | 19 | 0.10 | PT drained, attached → 0.3 |
| +`ka`×2 (fast attachment) | 0.42 | 11 | **28** | 0.13 | k_tr lifts off the ~19 floor |

So fixing the forward limb yields **(attached 0.3–0.4, XTOR ≈ 8–11, k_tr ≈ 28)** — essentially **rat α-MHC at 20 °C** (0.29, ~5–10, 24–31). The turnover lands at ~10, confirming §4. k_tr plateaus near ~28 (the attachment-limited ceiling, ≈ `ka_eff × overlap`); pushing to mouse-α's 50 needs still-faster forward kinetics, which the current ODE resists numerically (instability) — pointing to the near-equilibrium-power-stroke / cooperative routes (§7.3, §7.5).

---

## 6. Physiological critique of the "low-force third state" (the doubt — upheld)

A third attached state that **bears little force and detaches slowly**, padding the attached pool while a fast force-bearing state preserves k_tr, is **a kinetic device, not biology:**

- The post-power-stroke states (**AM·ADP**, **AM rigor**) are exactly where the lever arm has swung and **elastic strain is stored — the *high-force* states.** A strongly-bound, low-force species has no clean biochemical correlate. (`mechanism-evaluation.md` §1.1 already flags this for the model's P3.)
- The only genuinely low-force attached species is the **weakly/non-stereospecifically bound pre-stroke head (A·M·ADP·Pi)** — but that is a **fast-exchanging transient**, not a slow parking state. Repurposing it as a long-dwell buffer to inflate the attached count is fitting, not mechanism.
- It "works" mathematically only because a low-force state's slow filling doesn't drag the force transient. The moment the parking state carries real force (as a true post-stroke state must), its slow filling re-enters k_tr and k_tr collapses to the detachment rate — the 2-state wall again.

**Verdict: reject the low-force third state as a physiological explanation.** If a third state is used, it must be the **force-bearing AM·ADP with strain-dependent ADP release** (the Cooke–Pate / Beard low-ATP mechanism), not a low-stiffness park. The legitimate cousin is the **weakly-bound population** (§7.6), which is real but answers a different question ("total attached" vs "force-bearing").

---

## 7. Mechanisms that genuinely add a degree of freedom

Ranked for the warm, saturating-ATP, isometric mouse-α regime.

1. **Fix the hydrolysis throttle (§5).** Not a new mechanism — a correction. Raise effective `kah` to physiological (~100 s⁻¹), audit `xrate`. Drains PT, raises attached. *Necessary first step; demonstrated above.*
2. **ATP/ADP-gated, strain-dependent detachment.** ATP binding gates AM→M·ATP detachment; low ATP / ADP rebinding (AM·ADP trap) slows `g` → higher duty, higher force, lower k_tr, lower power. This is the Beard 2022 / Cooke–Pate axis and explains the k_tr(ATP) slope (50→25) **as motion along the frontier**. Re-enable `UseAtpK2`/`UseAtpKah` (currently off). *Reproduces the HF axis; does not move off the frontier.*
3. **Near-equilibrium (reversible) power stroke (Huxley–Simmons).** If weak⇌strong is fast in *both* directions, k_tr ≈ k₊+k₋ of that isomerization and can exceed `g`, while ATPase stays limited by slow irreversible detachment — the formal way to lift k_tr off `g`. Cardiac headroom is limited (k_tr = f+g with g dominant), but the model's power stroke is a one-way valve (k1 = 270 ≫ k_1 = 17); making it reversible is the physically-correct lever and is needed to reach k_tr 50 at high duty.
4. **SRX/DRX sequestration → whole-tissue economy.** Parking 40–60 % of heads OFF (ATPase ~1/200 s) lowers ATPase **per total myosin** without changing the cycling subset. If "< 5" were ever meant *per total head*, SRX delivers it from a per-cycling-head ~10. Caveat: it equally dilutes attached-per-total. Force-recruited (mechanosensing; Campbell 2018) — connects to `srx_reparam.md`.
5. **Cooperative / ensemble k_tr.** If k_tr is set by thin-filament activation or XB–XB cooperative recruitment rather than per-head `g`, the ensemble redevelops faster than any head's f+g — a true extra DOF, at the cost of a regulatory-unit model.
6. **Weakly-bound padding.** If "attached > 0.3" is the **total** (weak + strong) X-ray fraction (~0.6 in de Tombe), it's satisfied cheaply by low-force, ~ATP-free, fast-exchanging pre-stroke heads that don't gate k_tr. Pin the **denominator** before treating this as a target.
7. **Strain-dependent detachment for V_max (separate axis).** V_max is invariant to isometric `k2` in the model (−1.51 ML/s throughout) because it is set by the **negative-strain slip detachment** (`k2_R`, pwsd2 at s<0). V_max is **not** part of this conflict — keep it on its own strain-dependent knob.

---

## 8. Pinpointing the paradox & consequences

> **There is no deep paradox once two things are corrected.** (a) The kinetic identity `ATPase = r(1−r)·k_tr` is exact, temperature-free, and obeyed by real muscle (rat α: 0.29·0.71·24 = 4.9 ≈ measured). It says three isometric observables ride on two rate constants. (b) The **turnover target was misset**: mouse α-MHC at k_tr = 50 and duty ≥ 0.3 *must* turn over ~10–12 s⁻¹, not < 5; the < 5 belongs to β-MHC / cooler / per-total-myosin. With the corrected triple **(attached ≈ 0.3, k_tr ≈ 50, ATPase ≈ 10)** the targets are mutually consistent and sit *on* the physiological frontier.
>
> What actually ails the **model** is narrower and fixable: an **unphysiologically slow hydrolysis step** (`kah_eff` ≈ 8 vs ~100 s⁻¹) parks 91 % of heads in PT, starving the attached pool; and an **attachment-limited k_tr ceiling** (~28) that keeps it at the rat-α corner rather than the faster mouse-α one. Neither requires a low-force parking state.

**Recommendations**
1. **Re-target turnover to ~8–12 s⁻¹** for mouse α (or state the denominator if "< 5" is per-total-myosin with SRX).
2. **Fix the forward limb:** raise effective `kah` toward ~100 s⁻¹ (audit `xrate`); this alone moves attached 0.08 → 0.3 and is independently justified by hydrolysis-step literature. Then rebalance `kstiff` to hold isometric force ~85 (attached rose, so force will overshoot otherwise).
3. **Make the power stroke reversible** (raise `k_1` toward `k1`) and/or add cooperative recruitment to push the k_tr ceiling 28 → 50 at high duty.
4. **Re-enable ATP-gated detachment** (`UseAtpK2`/`UseAtpKah`) for the k_tr(ATP) slope; keep **strain-dependent slip detachment** for V_max; use **SRX** for whole-tissue economy.
5. **Drop the low-force third state**; if a third state is used, make it force-bearing AM·ADP with strain-dependent ADP release.

---

## Sources
- Rat cardiac myosin isoform kinetics, **20 °C** (f_app, g_app, k_tr, ATPase, duty 0.29): https://pmc.ncbi.nlm.nih.gov/articles/PMC2711735/
- **De Tombe et al. 2007**, temperature dependence of cross-bridge kinetics, rat myocardium (Q₁₀: k_tr 3.3, ATPase 3.4, tension 1.4, tension cost 2.5; attached-fraction Q₁₀ 1.1; k_tr 14.5/31.4/49.4 at 15/20/25 °C): https://pmc.ncbi.nlm.nih.gov/articles/PMC2277159/
- Beard et al. 2022, reduced cardiac power with low ATP (k_tr 34.9–40.7 → 25.8 s⁻¹, 8 vs 2 mM; model misses MgATP dependence): https://pmc.ncbi.nlm.nih.gov/articles/PMC9463691/
- Myosin ATP **hydrolysis step** ≈ 100–160 s⁻¹, near-equilibrium (not rate-limiting): https://pmc.ncbi.nlm.nih.gov/articles/PMC3853302/ ; https://www.med.upenn.edu/ostaplab/assets/user-content/documents/mie-reprint.pdf
- Myosin isoform ATPase / tension cost ∝ g_app (α vs β; mouse ≈ 1.5× rat): https://onlinelibrary.wiley.com/doi/10.1111/j.1440-1681.1995.tb02034.x ; https://pmc.ncbi.nlm.nih.gov/articles/PMC8835009/
- Fraction of myosin motors attached in isometric contraction at near-physiological temperature: https://pmc.ncbi.nlm.nih.gov/articles/PMC3136785/
- Mouse α/β detachment ~4× difference; ktr by slack-restretch (Stelzer/Moss): https://pmc.ncbi.nlm.nih.gov/articles/PMC4257873/
- SRX lifetime/fraction (cardiac, partial activation): https://pmc.ncbi.nlm.nih.gov/articles/PMC3077696/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC5425749/
- Internal: `mechanism-evaluation.md`, `low-ATP-force-enhancement-analysis.md`, `sources.md`
