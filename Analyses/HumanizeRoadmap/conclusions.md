# Humanizing the cross-bridge model: where we are, and which data can take it to the intact adult human myocyte

**Status:** 2026-10-27, review + literature survey, no new simulations. Target: a **healthy adult human intact
ventricular cardiomyocyte at 37 °C**, starting from a core identified on skinned rodent trabeculae (~22 °C,
saturating Ca, slack–restretch / k_tr / FV / staircase).

Most of the humanization work is **not on this branch**. It lives in worktrees and side branches (paths in §1).
Read §1 for what exists, §2 for what it adds up to, §3–4 for new data, §5 for the plan.

---

## 1. What already exists, and which thread did it

| # | analysis (date) | where it lives | session | what it established |
|---|---|---|---|---|
| A | **CalciumTwitch** (2026-09-23/24) | `ATP-wt-ca` (branch `CaActivation`, now prunable); copied into the worktrees below | before the retained session history | Ca layer design: a factorized gate on `ka_eff` (not legacy `NP`). The XB core is fast enough for a twitch. Literature thin filament gives a 5 % twitch, so amplitude is one steep, unidentified knob. Production switches `CaModel` / `CaInput`, `TestCaActivation` T2–T7 |
| B | **EHT identification report** (2026-10-02) | `Analyses/CalciumTwitch/EHT_identification_report.md` (in the worktrees) | before the retained session history | Mouse core + free thin filament + free Ca transient fits the Plate 6 control hiPSC-EHT twitch at 0.81 % rms (fit C). The fit is **not unique**: exact `k_on·Ca` symmetry, and Ca-transient vs XB-speed trade-off |
| C | **HumanizeXB** (2026-10-06) | `Analyses/HumanizeXB/` (committed in `deb4103`, also on main tree) | before the retained session history; its commit is discussed in [ATP-sensitive XB model fitting status](#local_31b9012d-1401-418d-8299-9caa4730ac53) | Humanizing is an **isoform (α→β-MHC) retune of rates, not of mechanics**. Literature β/α ratios ≈ 0.1 (k_tr, ADP release, hydrolysis, attachment). `params/human_v0_seed.m` (k2,k_2,kah,kamh,ka,kd ×0.1; k1,k_1 ×0.3). **Temperature is the dominant unknown** (rig °C never recorded). 400-line literature log |
| D | **HumanTwitch** (2026-10-08) | `C:\home\git\ATP-wt-human` (branch `HumanTwitch`, commit `08f7b1c`, plus uncommitted deck figures) | [Twitch modeling and T-world analysis](#local_14db0cff-90ee-4450-b896-893934650284) | Fixed literature hiPSC-EHT Ca transient (de Lange 2021). Tier `t1_isoT`: **net ×0.45** on k2,k_2,kah,kamh,ka,kd (isoform × Q10 inseparable, r = −1.00) gives 0.77 % rms. The answer swings ×0.21–×3.8 with the Ca-transient **shape** at the literature edges. Freeing SRX lets the fit escape into a non-human "fast cycle + SRX reservoir" regime. `params/human_twitch_t1_isoT.m` |
| E | **DrugTargets** (2026-10-12) | `ATP-wt-human/Analyses/DrugTargets/` (uncommitted) | [Twitch modeling and T-world analysis](#local_14db0cff-90ee-4450-b896-893934650284) | Inverts the ARPA force-only screen into lever magnitudes on the `t1_isoT` / `s3_free` cores |
| F | **TWorldCoupling** (2026-10-14) | commits `4e72ca7`, `93d45df`; `Docs/reference/tworld-ca-source.md` | [Twitch modeling and T-world analysis](#local_14db0cff-90ee-4450-b896-893934650284) | T-World (human myocyte, 93 states) as a Ca source: couplings `signal`, `land`, `xb`. EHT twitch with T-World's Ca: **net ×0.22–0.24** (vs 0.45 with de Lange's Ca) |
| G | **TnC linked state** (2026-10-22) | worktree `tnc-feedback-coupling-637355` (branch `claude/tnc-feedback-coupling-637355`, `65d5630`, `e39a49f`) | [Tworld feedback coupling and Ca buffering](#local_7b403c49-1f7b-471b-9aea-8c184db46f3a) → [Propose TnC state…](#local_590c1c30-86b1-4651-a7e1-e852eaaff41b) | Coupling `'xbtnc'`: separate TnC-Ca and permissive coordinates, so the T-World buffer is right for any filament set. **Length-dependent activation is already emergent** (ΔpCa50 0.08–0.11 from bound-head feedback); Land's β₁ shift would double-count. Closed-loop refit not run |
| H | **EHTRelativeTransfer** (2026-10-22) | `C:\home\git\ATP-wt-eht` (branch `EHTRelativeTransfer`, uncommitted) | [EHT data model analysis approach](#local_c9272069-eb0d-4a5b-bc7d-dadfc44da043) | A relative drug effect on the EHT **cannot be ratio-transferred** to an adult model. Activation level (EHT 84 % vs adult 25–60 %) flips the sign of RT50/RT90 effects. Mechanism-informed transfer works only with ≥ 4 shape metrics |
| I | **HumanizeLand + posts** (2026-10-22 → 27) | worktree `humanize-xb-model-231445` (branch `claude/humanize-xb-model-231445`, **uncommitted**: `Analyses/HumanizeLand/`, `params/human_land_g4valid.m`, Model changes for `CaModel 'land'` and `c_post`) | [Humanizing the XB model](#local_c5e641b1-41ad-4469-9778-33e5fd6a3454); post geometry from [HiPSC twitch measurement constraints](#local_1e2b3746-65a9-4a68-9dd6-58efc1b68198) | Land-in-T-World is calibrated **only** on 37 °C isometric twitch timing + force–Ca; its k_tr (32 s⁻¹) and FV are not data. With Land's thin filament ported, the rodent core matches Land's force–Ca and twitch fraction untuned. `g4_valid` (twitch + human k_tr 18): ka,kd ×0.57, k2,k_2 ×0.34, k1 & kah groups **unidentified**. EHT post compliance (1.92 N/m, L = 5.5 mm, ~4 % shortening) biases isometric fits ×1.2–1.5 too slow and exposes a wrong shortening response in the model |

Snapshots that matter: `params/human_v0_seed.m` (C), `params/human_twitch_t1_isoT.m` (D),
`params/human_land_g4valid.m` (I, uncommitted).

**Housekeeping risk.** Items E, H, I and part of D are uncommitted, across four worktrees that branch from the
same `93d45df`. Several of them edit the same Model files (`calciumActivation.m`, `getParams.m`,
`dPUdT_CombinedTransitions.m`). They need merging before any further human campaign; that is what
[Branch consolidation and AI code review](#local_a81c890d-f861-4b70-9e8d-acf00ec81b14) is for.

---

## 2. What it adds up to

### 2.1 One number has been estimated four ways

All are "how much slower the human cycle is than the mouse incumbent, at 37 °C", net of isoform × temperature,
on the detachment / attachment group (k2, ka, kah, kd):

| source | data | Ca / thin filament | net factor | caveat |
|---|---|---|---|---|
| Land 37 °C twitch + human k_tr 18 (`g4_valid`, I) | adult, model-derived twitch | Land filament, T-World Ca | k2 ×0.34, ka ×0.57 | only 2 directions identified; k_tr anchored on one lab |
| hiPSC-EHT, de Lange Ca (`t1_isoT`, D) | hiPSC-EHT | fixed literature Ca, free filament | ×0.45 (0.3–0.7 within de Lange's SD) | isometric; posts ignored |
| hiPSC-EHT, T-World Ca (`tw1`, F) | hiPSC-EHT | adult T-World Ca, free filament | ×0.22–0.24 | adult Ca on an immature tissue |
| hiPSC-EHT with posts (I) | hiPSC-EHT | de Lange Ca | ×0.45 × (1.19–1.38) ≈ 0.55–0.6 | the model's shortening response is wrong at L 5.5 mm |
| literature prior (C) | ratios | – | isoform 0.1 × Q10 factor 3.7–8.4 = 0.37–0.84 | rig temperature unknown |

**They agree to within a factor ~2 (0.22–0.6), and the spread is driven by the Ca transient, not by the
cross-bridges.** Every twitch-based estimate inherits the shape of the assumed transient. That is the central
identifiability fact of the whole thread: **a force twitch alone cannot separate Ca handling from cycle speed**.

### 2.2 What is identified, and what is not

| quantity | status | what pinned it |
|---|---|---|
| specific force, step size, stiffness | left at rodent values, defensible | literature: isoform-independent |
| thin filament force–Ca at 37 °C | matched (Land port) | Land/T-World steady state — itself model-derived |
| length-dependent activation | emergent (ΔpCa50 ~0.1) | G; not yet checked against human data |
| ADP release / detachment (k2) | identified (twitch + k_tr) | one human k_tr (18 s⁻¹, Milani-Nejad 2015) |
| attachment (ka·A0) | weakly identified | twitch |
| Pi release / power stroke (k1, k_1) | **not identified** | twitch is blind to it |
| hydrolysis (kah, kamh) | **not identified** | needs energetics |
| SRX / OFF↔ON | prior only; freeing it opens a non-physiological regime | needs a resting SRX or ATPase datum on the same tissue |
| FV curvature, V0, shortening deactivation, registration gate | **rodent values, unconstrained for human** | no human intact FV in the log |
| temperature (Q10) | folded into net factors | no explicit block |
| Pi / ADP in the intact cell | absent (model runs 0 Pi, 0 ADP) | – |
| passive / titin (human 30 % N2BA) | rodent | – |

### 2.3 Lessons to carry forward

1. **Do not fit Land's k_tr, FV or step responses.** They are time-compressed skinned structure, not data.
2. **Isometric twitch shape is dominated by Ca and the thin filament.** Pair every twitch with a measured
   Ca transient from the same preparation, or the XB speed floats by ×2–4.
3. **hiPSC-EHT is a different muscle** (activation 84 % vs adult 25–60 %, compliant posts, immature Ca
   handling). It is useful for the method, not as the human target.
4. **Temperature must become an explicit block** (Q10 per rate group). Without it, the 15 °C skinned human
   data and the 37 °C intact human data cannot be scored in one fit.
5. **The model's response to shortening (FV, registration gate, re-lengthening detachment) is rodent and
   untested in human.** The post analysis is the first place it visibly fails. An intact human heart muscle
   shortens; this is the biggest structural risk for the end goal.

---

## 3. New human data that can help (survey, 2026-10-27)

This extends `HumanizeXB/literature_log.md` and `HumanizeLand/literature_extension.md`; it does not repeat them.
Two delegated searches, levels as before: **V** = full text read, **U** = abstract or snippet only, **lead** =
from memory, not yet checked. DOIs from the search were inconsistent (one DOI was given for two different
Milani-Nejad papers), so **check every citation before quoting it**.

### 3.1 The honest summary first

- **No force–velocity, isotonic V_max, ATPase or X-ray data exist for intact human ventricle at 37 °C.**
  Not found by either search, nor by the earlier logs.
- **No non-failing intact human myocyte study with simultaneous Ca and sarcomere length at 37 °C was found.**
  All human living-slice data located are from end-stage failing hearts.
- **No raw human traces are publicly deposited** (Zenodo, Figshare, Dryad, PhysioNet: nothing usable).
- So the intact adult human myocyte **cannot be fitted directly**. It has to be reached in layers (§5): a human
  skinned 37 °C core, then the intact Ca and thin-filament layer on top, scored on intact twitch-level data.

### 3.2 Data that add something new

| # | dataset | prep, T | what it gives | why it matters here | level |
|---|---|---|---|---|---|
| 1 | **Lewalle, Milburn, Bell, Campbell, Niederer, JMCC 2025** (PMC13307878) | skinned human LV + LA fibres, 9 donors, 67 LV preps, **37 °C**, SL 1.94 / 2.23 µm (Kentucky) | **LV k_tr 15 ± 3 s⁻¹** (LA 40 ± 13), force–pCa at 2 SL (pCa50 5.79, nH 2.2), **quick ±1 % length steps**, energy cost 0.67 s⁻¹ kPa⁻¹; MCMC calibration of a Land-type model, which **overshoots k_tr 2×** and says so | The only human 37 °C dataset with k_tr, length dependence and step kinetics in the same preps. It scores the **core without any Ca-handling assumption**, which is the trade-off that sank every twitch fit (§2.1). Not deposited: ask Campbell for traces | V |
| 2 | **Land 2017 skinned human myocytes** (already in the extension, but not yet *used*) | skinned donor LV myocytes, 37 °C, 1.2 mM Pi | **ramp shortening 10 % in 100/200/300 ms → 0.2/0.4/0.5 of isometric force**; 0.5/1/2 % stretches in 10 ms with dip-and-recover over hundreds of ms | The only human 37 °C **shortening** data. Directly tests the model's FV / registration-gate response, which §2.3 item 5 flags as the biggest structural risk. Caveat: Land suspects glue compliance slowed the cells | V |
| 3 | **Milani-Nejad 2016** (AJP Heart) | intact RV trabeculae, 40 hearts, 37 °C, 0.5 Hz, L90 → L_opt | Frank–Starling and **length dependence of XB cycling (k_tr faster at L90)**, not altered by disease | The intact check for the emergent length-dependent activation of item G, and for the model's k_tr(SL) slope, which is already wrong-signed on the rodent ladder (`ktr_ladder_length_dependence`) | U (get full text) |
| 4 | **Chung 2019** Table 1 | intact non-failing trabeculae, 37 °C, 0.5–3 Hz | force 14.1 mN/mm², TTP 216, RT50 140 ms, dF/dt/F 7 s⁻¹ at 0.5 Hz; RT50/RT90 at higher rates | A **frequency series**: T-World produces a different Ca transient at each rate, so several rates give several Ca shapes for one core. That is the cheapest way to break the Ca/XB trade-off with existing intact data | V |
| 5 | **Chung 2018** (Physiol Rep) | intact RV trabeculae, 6 non-failing hearts, 0.5–3 Hz, **SR blocked** (CPA + ryanodine) | force–frequency stays positive without the SR; TTP and RT90 slow after block | A second, very different Ca condition on the same muscle. A held-out test of the Ca layer + core: T-World with SERCA/RyR off must reproduce it | V |
| 6 | **Holubarsch 1998** | intact human strips, isotonic vs isometric | **maximal work at ~40 % of isometric force**; work ×1.98 from 0.5 to 2 Hz | The only intact human afterloaded datum. It bounds FV curvature (max power at ~0.3–0.4 F0 is the textbook hyperbola; rat-β a/F0 0.05 in the log would put it lower) | U |
| 7 | **Milani-Nejad 2015** temperature series (in the log, but now with a role) | intact, K⁺ contracture | k_tr 6.6 / 14.1 / 18.1 s⁻¹ at 27 / 32 / 37 °C → **Q10 ≈ 2.8** | The anchor for an explicit temperature block. Note rat skinned Q10 ≈ 3.5 (de Tombe & Stienen 2007); Rice 2008's Q10 6.25 is a model constant, not a measurement | V |
| 8 | **Narolska 2005** (JMRC 26:39) | human skinned, 20 °C | ATPase and tension cost | The energetics anchor for the unidentified hydrolysis group (§2.2). Possibly already behind the log's tension cost; check | U |
| 9 | in vivo strain | echo (105 healthy adults), tagged / DENSE MRI | circumferential strain rate 1.3–1.6 s⁻¹, end-systolic fibre strain ≈ −0.14 | Bounds the velocities an intact human sarcomere actually sees: **~0.5–1.5 L/s, i.e. the low-velocity end of the FV**, where the registration gate and FV curvature live | U |
| 10 | **Bierhuizen 2025, Fischer 2019** living slices | failing LV, 37 °C, auxotonic, ~1 mN preload | TTP 166, TTR 270 ms; preload ≈ 83 % of L_max | Failing only, systematic offset vs trabeculae. Sanity check, not a target | V |

**Leads not yet checked (from memory, verify before use):**

- **Piacentino et al. 2003, Circ Res**: non-failing and failing isolated human myocytes, Ca transient and
  shortening at 37 °C.
- **Coppini 2013**: control (non-hypertrophic surgical) myocytes may include sarcomere shortening alongside the Ca
  transient we already use.
- **Pieske / Hasenfuss aequorin + force** in non-failing human trabeculae.

If any of these has non-failing Ca **and** shortening in the same cells, it is the closest thing to the end
target that exists.

---

## 4. Human-calibrated models beyond Land / T-World

| model | human data behind it | use for us | level |
|---|---|---|---|
| **Lewalle 2025** (Land-type, MCMC, 20 params) | item 1 | its posterior and its failure (k_tr 25 vs 15) are a benchmark: our core should do better on the same data | V |
| **Lewalle 2024, BJ** (Land + thick-filament OFF state with force feedback); code `github.com/AlexLewalle/OFFstateFfeedback` | Awinda 2020 donor myocardium, 37 °C, force–pCa at SL 1.9 / 2.3 | Says only **total-force feedback onto the OFF state** reproduces human length-dependent activation. Compare with our emergent LDA (item G) and our SRX force sensitivity `sigma1`; same idea, independently derived | V |
| **Regazzoni 2020 RDQ20** (Zenodo 3992553, `github.com/FrancescoRegazzoni/cardiac-activation`) | human: force–Ca and twitch kinetics only | Human k_basic 13 s⁻¹ vs rat 24–28 s⁻¹ (room temperature): its human 37 °C rates are *slower* than rat at 22 °C. That agrees with our net factors < 1 (§2.1). Twitch-tuned only; no XB identification | V |
| **T-World** | Margara ranges; no documented calibration of ν 10.16, μ 3.94, T_ref 80, koff, ca50 | Supplement section 4 ("Contractility representation") was not retrievable. Treat its contraction constants as undocumented; to learn more, diff `model_TWorld.m` against Margara or ask Tomek | V |
| **Rice 2008** Q10 table | rat | XB Q10 6.25 (g_app 2.5), Ca on 1.5, off 1.3: use the *structure* of a per-group Q10 block, not these numbers (human intact 2.8) | V |
| MUSICO (Prodanovic 2022) | known | – | V |
| Longobardi 2020, Strocchi 2022, Rodero 2022 (GP emulators) | rat / whole-heart / EP | method only: emulator-based history matching, if a human campaign gets expensive | U |

Not found: Negroni/Lascano human parameters, Marzban/Beard human, Tewari/Beard human, UT-Heart XB parameters,
and any human XB core calibrated on FV or ATPase.

---

## 5. Recommended path to the intact adult human myocyte

The data in §3 split cleanly into **core-level** data (skinned, 37 °C, saturating or clamped Ca: no Ca-handling
assumption) and **intact-level** data (twitches, where Ca, thin filament and core are entangled). Every attempt so
far has tried to identify the core from intact-level data alone. That is why the answer floated by ×2–4 with the
Ca transient. The fix is to identify in layers.

### Layer 0: decide the frame (no simulation)

1. **Get the Baker rig temperature.** It is still unrecorded, and it is a ×1.8 lever per 5 °C on every rate.
   One email.
2. **Add an explicit temperature block**: one Q10 per rate group (cycle, thin filament, SRX, passive viscosity).
   Defaults: cycle 2.8 (human intact) to 3.5 (rodent skinned), thin filament from Rice's structure. Then the
   rodent 22 °C fit, human 15 °C skinned data, Lewalle/Land 37 °C skinned data and 37 °C intact data can all be
   scored by one parameter set plus an isoform factor. Without it, isoform and temperature stay correlated −1.00.
3. **Consolidate the branches** (§1, housekeeping). Layers 1–2 need `CaModel 'land'`, `c_post`, `xbtnc` and the
   T-World Ca source in one tree.

### Layer 1: the human core at 37 °C, skinned (no Ca handling)

Score `g4_valid`, `t1_isoT` and `human_v0_seed` × Q10 against:

| target | source | protocol in our code |
|---|---|---|
| k_tr 15 ± 3 s⁻¹ (LV) | Lewalle 2025 | `ktrProbe` / `ktrP_k80` with their 20 % / 20 ms release |
| ±1 % step responses | Lewalle 2025 (traces needed) | new: quick-step protocol (the slack machinery can do it) |
| ramp shortening: 0.2 / 0.4 / 0.5 F0 after 10 % in 100 / 200 / 300 ms | Land 2017 | new: constant-velocity ramp, same as `runFVExperiment` but read at ramp end |
| 2 % stretch: dip below new steady state, recovery over hundreds of ms | Land 2017 | restretch machinery, reading `rsK80` |
| force–pCa at SL 1.9 / 2.3 | Lewalle 2025 / Awinda 2020 | `CaModel 'land'` or `twoState`, steady state |
| tension cost | Narolska 2005, Lewalle 2025 energy cost | `R2T / F0` |

**What this identifies that the twitch could not:** the detachment and attachment groups (k_tr + steps), the FV
curvature and registration-gate scaling (ramps), and with tension cost the hydrolysis group. The Pi step stays
weak unless a Pi-jump or Pi-dependence dataset is added (Stehle 2024 cites human myocyte k_tr ×1.5–2.4 per
+10 mM Pi).

**Cheapest first experiment (hours, all machinery exists):** run the Land 2017 ramp and stretch protocols and the
Lewalle k_tr protocol on the three existing human snapshots, as a forward prediction with nothing refitted.
If the ramp-end forces are far from 0.2 / 0.4 / 0.5, the rodent shortening response (§2.3 item 5) is the
blocker and needs work before anything intact is attempted.

### Layer 2: intact 37 °C (Ca handling on, core frozen from layer 1)

Run with T-World Ca and the `xbtnc` closed loop. Free only the thin filament (and the intact-vs-skinned shift
in Ca sensitivity). Score:

- Twitch timing over **frequency** (Chung 2019; T-World's own Ca at each rate): TTP, RT50, RT90, dF/dt/F.
- Intact k_tr 18 s⁻¹ and its **length dependence** (Milani-Nejad 2015/2016).
- Twitch / T_max ≈ 0.25–0.3, normalised force–frequency (Pieske, Mulieri, Holubarsch).
- **Held out:** SR-blocked twitch (Chung 2018), afterloaded maximal work at ~40 % F_iso (Holubarsch 1998), in vivo
  strain-rate range.

If the core from layer 1 cannot meet these with only the thin filament free, the gap is skinned→intact physiology
(lattice, Pi, ADP). That is then a named, testable difference, not something absorbed into "net factors".

### Layer 3: hiPSC-EHT as a test bench, not a target

Keep the EHT for the drug question and for method development, with posts in (`c_post`, L 5.5 mm). Use the
adult model from layer 2 as the healthy reference, and transfer drug effects by mechanism, never by ratio (item H).

### Data asks, in order of value

1. **Lewalle 2025 raw traces from the Campbell lab** (k_tr, steps, force–pCa at 37 °C, human LV). Highest
   leverage per effort; same type of protocols as our rodent battery.
2. **Baker lab:** rig temperature; and whether a β-MHC preparation (human donor, or pig LV as Patel 2023 did) can
   go through the slack–restretch / FV / k_tr battery at two temperatures. That would turn the isoform factor
   into a measurement.
3. **Full texts:** Milani-Nejad 2016, Holubarsch 1998, Piacentino 2003, Coppini 2013 (control myocyte shortening).
4. **Tomek:** the origin of T-World's contraction constants.
