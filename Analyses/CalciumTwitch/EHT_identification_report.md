# Identifying the Ca layer on paced control EHT twitches (Travis, Plate 6)

**Status:** final position of 2026-10-02. **Branch:** `CaActivation` (worktree `ATP-wt-ca`), uncommitted.
**Pipeline:** [`extractEHTPlate6.py`](extractEHTPlate6.py) → [`PrepareEHTControls.m`](PrepareEHTControls.m)
→ [`FitCEHTTwitch.m`](FitCEHTTwitch.m) (fit C, §3) → [`NonUniqueEHTTwitch.m`](NonUniqueEHTTwitch.m) (§4)
→ [`FitMyoSimEHTTwitch.m`](FitMyoSimEHTTwitch.m) (the MyoSim form, §8).
The path that led to fit C ([`FitEHTTwitch.m`](FitEHTTwitch.m), [`ProfileEHTTwitch.m`](ProfileEHTTwitch.m),
[`TestKoffEHTTwitch.m`](TestKoffEHTTwitch.m)) is in Appendix A. All simulate with [`simEHTTwitch.m`](simEHTTwitch.m).
**Data:** `data/travis/` (not tracked); processed CSVs in `data/travis/processed/`.

Numbers are measured unless marked *literature* or *estimate*. Literature values were recalled,
not re-checked for this report, and should be verified before they are cited.

## 0. The two preparations: what the model was built on, and what it is fitted to here

The cross-bridge core (`params/ExampleHighATP.m`) was identified on one kind of muscle, and this
report asks it to explain a very different one. Every result below should be read against this
table.

| | **mouse trabecula: what the XB core was fitted to** | **EHT: what is fitted here** |
|---|---|---|
| cells | adult mouse right-ventricular myocytes (Beard 2022: male C57BL/6J, 11–14 weeks; the 2026 Baker-lab recordings are said to be mouse, not yet confirmed in writing) | **human iPSC-derived cardiomyocytes** (lab information, 2026-10-06; not stated in the files themselves). Line, age and purity unknown; usually with fibroblasts in a hydrogel (*estimate*) |
| tissue | a dissected RV trabecula, about 140 µm wide × 100 µm thick (manuscript §3.1) | a tissue cast between two posts in a Mantarray 24-well plate; geometry not in the file |
| membrane and Ca handling | **demembranated (skinned):** no sarcolemma, no SR, no Ca transient. Ca is set by the bath | **intact:** the tissue beats spontaneously (0.13–0.24 Hz in the control wells), and paced twitches are evoked by field stimulation. Its own excitation, Ca entry, SR release and reuptake make the Ca transient, which is unmeasured |
| activation | **clamped at saturating Ca** (pCa 4.51) and held on throughout each run | **field-stimulated twitches**, 1 Hz, biphasic 2 + 2 ms pulses. The exported stimulus trace is **not synchronised** with the force record (lab information), so the stimulus-to-onset delay is unknown and is fitted as an offset (§2) |
| energy supply | 8 mM (and 2 mM) MgATP in the bath, buffered by creatine phosphate and creatine kinase | cell metabolism (mitochondria, glycolysis); [ATP] not controlled |
| myosin and thin-filament isoforms | adult mouse ventricle: mainly **α-MHC** (fast), cardiac TnI, adult titin (*literature*, recalled) | hiPSC-CM: mainly **β-MHC** (slower, 2–3× lower ATPase), much slow-skeletal TnI (Ca-sensitizing), fetal N2BA titin, immature and partly disarrayed sarcomeres (*literature*, recalled) |
| mechanics | motor-controlled length with a force transducer; protocols: slack/re-stretch train, ktr, stretch staircase, isovelocity ramps (force–velocity) | contraction against a **compliant post** ("Skeletal Muscle (12x)" stiffness factor). Force comes from magnet displacement, so the twitch is auxotonic, not isometric (platform description, recalled) |
| sarcomere length | no optical SL; nominal SL = 2 × muscle length, **SL 2.2 µm** in the core | unknown; the model here runs at **SL 2.0 µm** ([Paced1Hz_report §0](Paced1Hz_report.md)) |
| force readout | absolute stress: **F_ref ≈ 69 kPa** maximal active (n = 3, rundown-corrected) | **µN, active only**: Pulse3D subtracts the baseline, so diastolic and passive force are not in the data |
| temperature | bath temperature **not documented** (open TODO in manuscript §2.1.1) | presumably 37 °C (not in the file) |
| what the fit constrained | cross-bridge kinetics, SRX, the attachment gates and the passive/visco-elastic elements, **all at saturating Ca only**. The thin-filament (Ca) regulation was never calibrated | the Ca layer (transient + thin filament) and one force scale σ, with the core frozen |

**Species: human iPSC-derived** (lab information, 2026-10-06). The files do not state it. An earlier
version of this report called it undocumented. The slow twitch (10 %→peak 176 ms, RT50 152 ms) and
slow spontaneous beating (0.13–0.24 Hz) fit hiPSC-EHTs. The hiPSC-CM reference ranges used below
therefore apply.

**Three consequences follow for this report:**

1. **The Ca layer is new territory.** Nothing in the core's data says how force depends on Ca.
   The whole Ca→force coupling therefore comes from the EHT twitch alone. That is why it is so
   poorly identified (§4, §A.2).
2. **The core's absolute kinetics are likely too fast for these cells, and the twitch cannot tell (§4.2: a 2× slower or faster core fits equally well).** Species (mouse vs human),
   myosin isoform (α- vs β-MHC), maturity and temperature all differ; a β-MHC hiPSC core would be slower. At high activation the cross-bridge rates *do* shape the twitch
   (§6: `k2` ×0.5 slows the rise by 30 ms). So the fitted Ca layer may be partly compensating for
   a core that is too fast or too slow.
3. **Amplitude is not comparable, only shape is.** The core's 64 kPa (SL 2.0) is adult skinned
   muscle at full activation. The EHT's 429 µN depends on cross-section, myocyte content and
   maturity, none of which is known. One scalar, σ (µN per model kPa), absorbs all of it.

## Final position

| | |
|---|---|
| **Target** | Plate 6, 1 Hz field stimulation, **15 control-like wells**, beat-averaged (60 beats each), then averaged over wells (§1). The tissue is **human iPSC-derived** (§0). |
| **Best fit found: fit C** | Six parameters of the Ca layer are refitted and one force scale is solved; the cross-bridge core is unchanged (every changed value is listed in §2). The result: `tf_k_off` **3.8 s⁻¹**, `tf_k_on` 4.15 µM⁻¹s⁻¹, `K_coop` 20 (its upper bound), a Ca transient peaking 120 ms before force, σ 7.9 µN/kPa. **RMS 0.81 %** of peak (data SEM 0.35 %). Timing metrics are within 3 ms except RT90 (255 against 243 ms). The late tail beyond 550 ms is missed (§3.2). |
| **The fit is not unique (§4)** | **Exactly:** scaling both Ca levels by any λ while dividing `k_on` by λ, and swapping the rise and decay time constants, give identical force, states and ATP (difference exactly 0). Absolute Ca and `k_on` are therefore not identified, only their product. **Practically:** refits under other fixed assumptions fit as well or nearly as well: diastolic Ca 0.05 µM (RMS 0.74 %), a cross-bridge core with all cycling rates halved (0.75 %) or doubled (0.95 %), and diastolic Ca 0.20 µM (1.16 %). They do so with different parameters, internal states and energetics: the ATP per beat halves or doubles with the core speed, and the diastolic ATP rate differs 4× between the diastolic-Ca alternatives. |
| **What the data do pin down, given the model** | the attachment-gate time course; a Ca peak about 1.6× Ca₅₀ (activation 75–87 % of saturating force in every alternative); a steep effective nonlinearity (n_H 4.6–7, with `K_coop` on its bound in most); a slow thin-filament step (effective off-rate tens of s⁻¹); and a 90–145 ms delay from Ca peak to force peak (§4.4). |
| **What it does not pin down** | absolute Ca (diastolic and peak); `k_on` separately from Ca; rise versus decay of the Ca transient; the split of Ca→force kinetics between Ca handling and the thin filament (§A.3); per-head ATP (30× range across fits); diastolic activation and diastolic ATP; the cross-bridge core's own speed, and with it **the ATP cost**: tissue-level ATP per beat (AUC × σ) is 22–30 for every Ca-layer variant, but 12 for a 2× slower core and 45 for a 2× faster one, with equally good fits (§4.2). |
| **Plausibility (§5)** | no hard limit crossed. At or beyond the edge: n_H 7.0 with `K_coop` on its bound; k_off 3.8 s⁻¹, plausible only as a lumped activation step; 86 % activation per twitch; diastolic ATP about 5× basal (set by the assumed diastolic Ca). |
| **Cross-bridge core untouched (§2.2)** | A field-by-field audit of fits B, C and M against the snapshot: no cross-bridge, SRX, passive or series-element parameter differs. The changes are the Ca layer, SL0 2.2 → 2.0 (with its derived initial state), the isometric protocol settings and the output scale σ. |
| **MyoSim form (§8)** | Our `twoState` thin filament is MyoSim's kinetics; the two differ only in how partial overlap enters. The MyoSim form is now available as `CaModel = 'myosim'`. Refitted (fit M), it gives the same twitch (RMS 0.82 %), states and ATP to within 1 %, with thin-filament constants within 7 % of fit C. |
| **To make it unique (§7)** | a **relative Ca transient** in the same wells (fixes the input shape and the Ca-versus-filament split); a **[Ca]ₒ series** (fit C predicts +13 % force for a 1.5× larger transient); a Ca calibration (breaks `k_on·Ca`); and the platemap, species, post stiffness and geometry from Travis. |

## 1. The data and what to fit

### 1.1 What is in `data/travis`

| item | content | usable for fitting? |
|---|---|---|
| `Plate 6_24 hrs_ march 25-26_stimulation_….xlsx` | Mantarray / Pulse3D 3.2.2 "full" export: **continuous active twitch force (µN) at 100 Hz for 24 wells, 63.7 s**, stimulus pulse trace (biphasic 2 + 2 ms, 1 Hz, first pulse at 6.9 ms; **not synchronised** with the force record, per the lab), per-twitch and aggregate metrics | **yes: the target** |
| `Plate 6_24 hours march 25-26_spontaneous….xlsx` | same plate, same day, unpaced: 0.13–0.61 Hz, irregular | not for a first fit (no stimulus time, rate differs per well); used only to classify wells |
| `EHT dose response trials…/*.prism` (9 drugs) | per-well **% change from each tissue's own baseline**, ANOVA and Dunnett tables | not as waveforms. Useful later as drug-effect *targets* (e.g. mavacamten −40 / −78 / −98 %) |
| `drug-vs-control-report/` | summary of the Prism files | context |

Metadata of the stimulated file: plate `ML26063409-6`, recorded 2026-03-26 15:46 UTC,
post stiffness factor "Skeletal Muscle (12x)", **platemap "N/A"**. The file itself does not
say which well received what.

### 1.2 Choosing the control wells

The plate splits cleanly into two groups (figure, panel C). The well-to-drug map is inferred from
the pattern, **not from a platemap**:

| wells | paced peak→RT50 (ms) | spontaneous rate (Hz) | reading |
|---|---|---|---|
| A1–D3 (12 wells), D4, D5, D6 | 141–154 | 0.13–0.24 | one homogeneous group: **controls** |
| A4, B4, C4 | 129–141 | 0.26–0.27 | faster relaxation and beating, weakest |
| A5, B5, C5 | 122–126 | 0.34–0.38 | stronger |
| A6, B6, C6 | 99–101 | 0.52–0.61 | strongest |

Rows A–C of columns 4–6 show a **graded, dose-like lusitropic and chronotropic effect**,
which is the isoproterenol signature. The Iso Prism file ("March 19–26", 0.1 / 0.5 / 1 µM)
fits a 3 × 3 layout. Row D of columns 4–6 behaves like columns 1–3. Whether columns 1–3 are
DMSO or an inactive compound cannot be told from the file. The recording is "24 h", so the
controls are most likely 24 h DMSO wells. **Ask Travis for the Plate 6 platemap before these are
used as anything more than an initial fit.**

![control wells](results/eht_controls.png)

- **A.** Beat-averaged force of all 24 wells: blue are the controls, orange the Iso-signature wells.
- **B.** The 15 controls, each normalized to its peak, with mean ± SD. The **shape is nearly
  identical across wells**. The SD band is barely visible.
- **C.** The classification.

### 1.3 The target

Each well's force is cut into beats from −50 to +950 ms around each stimulus. The first 3 beats
are dropped, leaving 60 beats. The beats are averaged and re-zeroed on the pre-stimulus window.

| control template | value |
|---|---|
| peak active force | **429 µN** (well SD 48, range 359–503) |
| beat-to-beat CV of peak | 3.2 % (median over wells) |
| stimulus → peak | 230 ms |
| 10 % → peak | 176 ms (Pulse3D: 182) |
| peak → 50 % / 90 % relaxation | 152 / 243 ms (Pulse3D: 147 / 239) |
| width at 50 % | 275 ms (Pulse3D: 275) |

The small differences from Pulse3D come from averaging on a 10 ms grid before the metrics are
taken. Force is **active only**: Pulse3D removes the baseline, so the recording carries no
diastolic or passive force and none was fitted. Amplitude varies well to well by 11 % (CV), and
shape hardly at all. That split is what a per-tissue force scale (σ) is for.

## 2. Model, and every parameter that was changed

| layer | setting |
|---|---|
| cross-bridge core | `params/ExampleHighATP.m` as of commit **46094fb** (2026-09-04; the main branch changed it later, at dbea4ee), **frozen** (only SL0, below); isometric, fixed ends, series elastic element in the loop |
| thin filament | `CaModel = 'twoState'`: `Jon = k_on·Ca·N·(1 + K_coop(1−N))`, `Joff = k_off·Pu·(1 + K_coop·N)`, attachment gate `Pu/(1 − bound)` on `ka` |
| Ca input | `rice2008` bi-exponential (`caTransient.m`), period 1 s; the whole model beat is shifted by a fitted offset relative to the recorded stimulus, which is not synchronised with force |
| protocol per evaluation | 1 s at diastolic Ca, then 2 beats, the 2nd analysed (fit C checked over 5 beats: drift 2e-4 kPa) |
| objective | `lsqnonlin` on (σ·F_model − F_data)/F_peak,data over 0–950 ms (96 points). σ is solved in closed form at every evaluation, so the fit is on **shape** |

### 2.1 Changed parameters: original against optimized

"Original" is the value the model would use without this work: the `ExampleHighATP` snapshot, or
the code default where the snapshot does not set it (Rice 2008 for the transient, Marzban 2020 for
the thin filament). Fits A and B are earlier fits (§A.1) and are shown for comparison.

| group | parameter | meaning | original (source) | fit A | fit B | **fit C (final)** | fit C / original | status in fit C |
|---|---|---|---|---|---|---|---|---|
| switches | `CaModel` | thin-filament model | `'none'` (saturating Ca) | `'twoState'` | `'twoState'` | `'twoState'` | — | set |
| | `CaInput`, `CaTransientKind`, `CaPeriod` | Ca source | — | analytic, `rice2008`, 1 s | = | = | — | set |
| geometry | `SL0` | sarcomere length (µm) | 2.2 (snapshot) | 2.0 | 2.0 | **2.0** | ×0.91 | set (Paced1Hz_report §0) |
| Ca transient | `tau1` | rise time constant (ms) | 20 (Rice 2008) | 114.8 | 150.0 (bound) | **111.9** | ×5.6 | **fitted**; exchangeable with `tau2` (§4.1) |
| | `tau2` | decay time constant (ms) | 110 (Rice 2008) | 124.3 | 176.4 | **108.3** | ×0.98 | **fitted**; exchangeable with `tau1` |
| | `t0` / beat offset | Ca onset relative to the recorded stimulus (ms) | +5 (Rice 2008 onset delay) | +55.8 | +18.7 | **−2.7** | — | **fitted**, ±3.7 ms (bounds −150…+60) |
| | `Cadia` | diastolic Ca (µM) | 0.09 (Rice 2008) | 0.10 | 0.10 | **0.10** | ×1.11 | **assumed** (fixed); not identified (§4.2) |
| | `Camax` | peak Ca (µM) | 1.45 (Rice 2008) | 1.0 | 1.0 | **1.0** | ×0.69 | **assumed** (fixed); not identifiable (§4.1) |
| thin filament | `tf_k_on` | Ca-binding/activation rate (µM⁻¹s⁻¹) | 101.185 (Marzban 2020) | 182.8 | 657.7 | **4.146** | ×0.041 | **fitted**; identified only as `k_on·Ca` and as a ratio to `k_off` |
| | `tf_k_off` | deactivation rate (s⁻¹) | 723.852 (Marzban 2020) | 723.852 (fixed) | 723.852 (fixed) | **3.824** | ×0.0053 | **fitted**, ×/÷1.21 locally; moves with the core speed (§4.2) |
| | `k_on/k_off` | steady sensitivity (µM⁻¹) | 0.140 | 0.253 | 0.909 | **1.084** | ×7.8 | derived |
| | `tf_K_coop` | cooperativity | 9.6846 (Marzban 2020) | 13.77 | 20.0 (bound) | **19.97** | ×2.06 | **fitted, on its upper bound** (20) |
| scale | σ | µN of EHT force per model kPa | — (new) | 155.5 | 7.43 | **7.887** | — | solved in closed form |
| derived | pCa₅₀ / n_H (steady, SL 2.0) | | 5.33 / 3.33 with the Marzban constants (SL 2.2 probe, conclusions §2.2) | 5.59 / 5.2 | 6.15 / 5.9 | **6.21 / 7.0** | Ca₅₀ ×0.13 (the original is from a different snapshot at SL 2.2) | derived |

**Unchanged:** every other parameter of `ExampleHighATP.m`. That covers all cross-bridge rates
(`kah`, `kamh`, `ka`, `kd`, `k1`, `k_1`, `k2` and their strain dependences), the stiffnesses,
SRX (`ksr0`, `kmsrd`, …), the attachment gates (`A_reg`, `A_def`), passive titin (`k_pas`, `gamma`,
`Lsc0`), the Maxwell dashpot, the series element (`kSE`) and filament overlap. They were frozen
in every fit. §4.2 shows what happens if the core's cycling rates are allowed to differ.

**Bounds used for fit C:** τ_rise 10–150 ms, τ_decay 80–300 ms, offset −150…+60 ms, `tf_k_on`
0.1–3000, `tf_K_coop` 1–20 (`parameterBounds.m`), `tf_k_off` 0.5–2000 s⁻¹.

### 2.2 Audit: no cross-bridge parameter was tuned

**Confirmed.** Every fit in this report kept the cross-bridge core exactly as calibrated. The
one exception is the force scale σ, which is not a model parameter but a factor applied to the
model's output (µN of EHT force per model kPa). To check this, the parameter struct the simulator
actually used for fits B, C and M (`simEHTTwitch`'s `S.P`) was compared, field by field, with
`getParams` of the unmodified `params/ExampleHighATP.m` ([`CompareCaTransientsEHT.m`](CompareCaTransientsEHT.m);
listing in `results/eht_param_audit.txt`).

| | fields compared | fields that differ | cross-bridge, SRX, passive or series-element fields among them |
|---|---|---|---|
| fit B | 312 | 12 | **0** |
| fit C | 312 | 13 | **0** |
| fit M | 312 | 13 | **0** |

The differing fields are all of one of three kinds:

| kind | fields | why |
|---|---|---|
| **Ca layer (fitted or set)** | `CaModel` (`'none'` → `'twoState'` / `'myosim'`), `CaInput` (`'const'` → `'analytic'`), `CaFirstBeat` (0 → 1 s), `CaTransientOpt` (the transient's τ, levels and onset), `tf_k_on`, `tf_k_off` (not in fit B), `tf_K_coop`, and the derived struct `ca` | this is what was fitted (§2.1) |
| **sarcomere length (set, not fitted)** | `SL0` 2.2 → 2.0 µm, and what `getParams` derives from it: the strain-grid pivot `LXBpivot` and the initial state `PU0` (which also gains one entry for the Ca state) | the twitch runs at SL 2.0 ([Paced1Hz_report §0](Paced1Hz_report.md)). This changes overlap and passive force through the core's own, unchanged length dependence |
| **protocol** | `Vums = 0`, `v = 0` | isometric (fixed ends) |

Two qualifications:

- The thin-filament gate multiplies the attachment rate `ka_eff` (§2). So the Ca layer *acts on*
  the cross-bridge cycle, but no cross-bridge rate constant was changed.
- Cross-bridge parameters **were** varied, but only as probes, never in a reported fit. These
  are the one-at-a-time sensitivities (§6) and the "core cycling rates ×0.5 / ×2" alternatives
  of §4.2, which show that the twitch cannot constrain them.

## 3. Fit C: the best fit found

Script: [`FitCEHTTwitch.m`](FitCEHTTwitch.m). Results: `results/eht_fitC.mat`, `results/eht_fitC.png`.

### 3.1 The close-up

![fit C close-up](results/eht_fitC.png)

- **A.** Data (mean ± SD of 15 wells) and fit C.
- **B.** The residual, with ±2 SEM of the template in grey.
- **C.** The same on a log scale, which shows the late tail.
- **D.** The chain Ca → permissive sites → attachment gate → force, normalized.
- **E.** The effective site off-rate `k_off(1 + K_coop·N)` and its time constant.
- **F.** The myosin head pools.
- **G.** Net ATP turnover, with its AUC and the steady rate at saturating Ca.
- **H.** Steady force–pCa with the fitted constants; the twitch peak is marked.
- **I.** The force–Ca loop, as a relative Ca recording would see it.
- **J.** One-at-a-time sensitivities (§6).

### 3.2 Fit quality, and the one thing it misses

| | data | **fit C** | fit B | fit A |
|---|---|---|---|---|
| RMS error (% of peak) | SEM 0.35 | **0.81** | 3.8 | 8.5 |
| stimulus → peak (ms) | 230 | **228** | 208 | 198 |
| 10 % → peak | 176 | **175** | 140 | 124 |
| peak → RT50 | 152 | **153** | 173 | 166 |
| peak → RT90 | 243 | **255** | 314 | 405 |
| width at 50 % | 275 | **272** | 274 | 255 |
| peak force (µN) | 429 | 429 (σ × 54.4 kPa) | | |

**The late tail is the one visible systematic miss** (panel C). From about 550 ms the data keep
1–3 % of peak force and decay slowly to about 0.3 % by 900 ms. The model's force falls away
faster below about 3 %. This tail carries most of the structured residual (panel B, lag-1
autocorrelation 0.93). Candidates:
- a slow tail of the real Ca transient (Na/Ca exchange, *estimate*), which a single-τ
  bi-exponential cannot produce;
- visco-elastic re-lengthening of the tissue on the compliant post;
- an artefact of Pulse3D's baseline subtraction.

A relative Ca recording would separate the first from the other two.

### 3.3 Internal states

| quantity | diastole (at stimulus) | peak or minimum | at (ms, data time) |
|---|---|---|---|
| Ca (µM) | 0.100 (0.104 at 950 ms) | 1.0 | **108** |
| permissive sites `1−N` | 0.010 | 0.82 | ~200 |
| attachment gate | 0.0067 | 0.76 | 210 |
| attached heads `p1+p2` | 0.0032 | **0.23** | 205 |
| force (developed, kPa) | 0 (absolute 0.81, of which passive 0.04) | 54.4 | **228** |
| primed `PD` | 0.78 | 0.64 (min) | ~200 |
| `PT` (M·ATP) | 0.015 | 0.081 | ~200 |
| SRX (`P_SR + P_SRD`) | 0.20 | **0.04** (min) | ~200 |
| effective site off-rate | 79 s⁻¹ (τ 13 ms) | **18 s⁻¹** (τ 55 ms), min | ~200 |

- **There is a 100 ms lag from Ca peak to gate peak** (108 → 210 ms). That is the slow thin
  filament. Force follows the gate within 20 ms, because the cross-bridge core is fast. Ca peaks
  **120 ms before force**. This interval, unlike times from the stimulus, does not depend on the
  synchronisation.
- **The relaxation shape is a cooperative switch.** While the filament is on, sites switch off
  at only 18 s⁻¹. Once some have switched off, the `K_coop·N` term accelerates the rest to
  79 s⁻¹. The result is the shoulder followed by the accelerating fall of the data (panels D, E).
- **The twitch runs near saturation.** Peak force is 86 % of saturating developed force
  (panel H). Ca peaks at 1.6 × Ca₅₀.
- **Heads:** at peak 23 % are attached, 64 % are still primed in `PD`, and SRX is drawn down from
  20 % to 4 %. The force-sensitive SRX exit empties the OFF pool every beat. The model's SRX is a
  fast-exchanging pool (dwell about 30 ms, [Paced1Hz_report obs. 5](Paced1Hz_report.md)), so
  this is not a statement about mant-ATP-type SRX.
- **Diastole carries some activation:** gate 0.0067 and 0.77 kPa of cross-bridge force (1.2 % of
  saturating, or about 7 µN after scaling). This is invisible in baseline-subtracted data and
  set mainly by the fixed `Cadia` (§6).

### 3.4 ATP turnover

**Definition.** Net hydrolysis flux, in ATP per head per s:
`(PT→PD − PD→PT) + (SR→SRD − SRD→SR)`, i.e. rates `(1 − 2) + (13 − 14)` of
`dPUdT_CombinedTransitions`. The reverse steps are ATP resynthesis and are subtracted. The AUC is
taken over one full 1 s beat, so at 1 Hz **ATP per head per beat equals the mean rate per
second**.

| | **fit C** | fit B | fit A |
|---|---|---|---|
| **AUC: ATP per head per beat (= per s at 1 Hz)** | **2.90** | 3.16 | 0.19 |
| peak rate | 10.1 s⁻¹ (at ~210 ms) | 10.6 s⁻¹ | 0.56 s⁻¹ |
| steady rate at saturating Ca (pCa 4.51) | 11.6 s⁻¹ | | |
| diastolic rate | 0.139 s⁻¹ | 0.163 s⁻¹ | 0.042 s⁻¹ |
| basal share of the beat (diastolic rate × 1 s / AUC) | 5 % | 5 % | 22 % |
| economy: developed force per (ATP/head/beat) | 18.7 kPa | 17.7 kPa | 14.8 kPa |
| **tissue level: AUC × σ** | **22.9** | 23.5 | 29.6 |

- **The twitch runs the myosin at 87 % of its saturating ATPase at peak** and at a mean of 25 %
  over the beat.
- **Tissue-level turnover (AUC × σ) is robust only to the Ca layer.** It is 22–30 across every
  Ca-layer variant in this report, while per-head ATP varies 30×. But it scales almost linearly
  with the speed of the cross-bridge core, which the twitch cannot determine (12 at half speed,
  45 at double speed, with equally good fits; §4.2). Its absolute scale is inherited from the
  mouse core.
- **The DRX/SRX split is bookkeeping.** Fit C gives −1.86 (DRX step) and +4.76 (SRX step)
  ATP/head/beat. The SRX loop PT→SR→SRD→PD→PT hydrolyses at SR→SRD and resynthesizes at PD→PT,
  with no net ATP. Only the sum means anything.
- **The diastolic rate, 0.14 s⁻¹, is high:** about 5× the basal DRX turnover in mant-ATP chase
  experiments (about 0.03 s⁻¹, *literature*). It follows from the residual diastolic activation
  (§3.3), and so from the fixed `Cadia`. At `Cadia` 0.05 µM the beat's ATP falls 6.2 % (§6).

### 3.5 Force scale in plausible units

σ = 7.9 µN/kPa: the model's full-power stress acting on 0.0079 mm² would give 429 µN. For
assumed EHT cross-sections:

| assumed EHT cross-section (mm²) | 0.05 | 0.1 | 0.25 | 0.5 |
|---|---|---|---|---|
| measured peak stress (kPa) | 8.6 | 4.3 | 1.7 | 0.86 |
| XB power scale vs the mouse core (σ / (1000·A)) | 0.16 | 0.079 | 0.032 | 0.016 |

hiPSC-EHT active stresses are usually reported at about 0.5–5 mN/mm² (*estimate*), which puts
the scale at a few percent: low myofibril density and immature sarcomeres. σ also absorbs the
11 % well-to-well amplitude CV, so a **per-well σ with a shared shape** is the natural next model.

### 3.6 The fitted Ca transients against the untuned ones

![Ca transients](results/eht_ca_transients.png)

- **A.** Absolute Ca.
- **B.** Normalized to (Ca − Ca_dia)/(Ca_peak − Ca_dia), which is what the twitch fit actually
  determines (§4.1).

"Original" means the transient each model would use without tuning:
- for `twoState` (and `'myosim'` in this code), the `caTransient` default `rice2008`;
- for MyoSim itself, the pCa trace of MATMyoSim's `simple_twitch` demo
  (`code/demos/getting_started/simple_twitch/protocol.txt`, made by `generate_twitch_protocol.m`,
  a two-compartment SR pulse). Its stimulus at 0.1 s is aligned to t = 0.

| transient | diastolic (µM) | peak (µM) | time to peak (ms) | peak → 50 % decay (ms) | status |
|---|---|---|---|---|---|
| original `twoState` input: Rice 2008 defaults (rat, 24 °C) | 0.090 | 1.45 | 48 | 97 | untuned |
| original MyoSim demo (`simple_twitch`) | 0.190 | 0.747 | 38 | 71 | untuned; trace ends at 288 ms |
| fit B (`twoState`, fast filament) | 0.100 | 1.0 | 180 | 275 | fitted |
| **fit C** (`twoState`) | 0.103 | 1.0 | **108** | **185** | fitted |
| fit M (`'myosim'`) | 0.103 | 1.0 | 105 | 188 | fitted |

Times are from the recorded stimulus, which is not synchronised with force (§0). Fits C and M
place the onset 3–4 ms before it, fit B 19 ms after it.

- **Fits C and M have the same transient.** The MyoSim form changes the transient by no more
  than 3 ms in either metric.
- **Fitted transients are 2–4× slower than either untuned input:** the time to peak is about
  105 ms against 38–48 ms, and the half decay 185 ms against 71–97 ms. Both untuned transients
  are rodent-like (rat 24 °C, and a generic MyoSim demo); the EHT twitch is slow (§0).
- **Fit B's transient is slower still** (180 / 275 ms). With a fast filament, the slowness of the
  twitch has to come from Ca. With a slow filament (fits C, M) it is shared (§A.3). A measured
  Ca transient would decide between them (§7).
- **Absolute levels are not comparable.** The fitted levels 0.10 → 1.0 µM are a gauge fixed
  before fitting (§4.1). Only panel B is a result.

## 4. The fit is not unique

Script: [`NonUniqueEHTTwitch.m`](NonUniqueEHTTwitch.m). Results: `results/eht_nonunique.mat`,
`results/eht_nonunique.png`; refits cached in `results/eht_nonunique_fits.mat`.

Fit C is the best fit *under its assumptions*. It is not unique, for two different reasons:
exact symmetries of the model, and fixed assumptions that the twitch cannot test.

### 4.1 Exact symmetries: infinitely many parameter sets give the same twitch

Each of the following was simulated against fit C and compared sample by sample:

| transformation of fit C | max relative change in force | in gate | in bound heads | in ATP per beat |
|---|---|---|---|---|
| both Ca levels ×0.5 (0.05 → 0.5 µM), `tf_k_on` ×2 | **0** | **0** | **0** | **0** |
| both Ca levels ×2 (0.2 → 2.0 µM), `tf_k_on` ×0.5 | **0** | **0** | **0** | **0** |
| `tau1` ↔ `tau2` swapped | **0** | **0** | **0** | **0** |

1. **`k_on·Ca` (exact, continuous).** The ON flux is `Jon = k_on·Ca·N·(1 + K_coop(1−N))`, so Ca
   enters only multiplied by `k_on`. Any Ca amplitude λ with `k_on`/λ is the same model, down to
   the internal states and ATP. Absolute Ca and `k_on` are therefore **not identifiable from force
   at all**. Only the Ca₅₀ *relative to the transient* is meaningful: the peak is 1.6 × Ca₅₀. The
   fixed Ca levels (0.10 → 1.0 µM) are a choice of gauge, not a result.
2. **`tau1 ↔ tau2` (exact, discrete).** The peak-normalized bi-exponential is symmetric in its two
   time constants. Every optimum sits at τ₁ ≈ τ₂ (fit C: 112 / 108 ms), where it is an alpha
   function `t·e^(−t/τ)`. One Ca time constant (about 110 ms) is identified; which one is "rise"
   and which is "decay" is not.

### 4.2 Equally good fits under different fixed assumptions

Fit C fixes several things the twitch cannot test: the diastolic Ca, the speed of the
cross-bridge core (adult mouse, §0), and the bounds. Each alternative below changes **one**
assumption and refits **all six** Ca/thin-filament parameters (plus σ), starting from fit C.

![non-uniqueness](results/eht_nonunique.png)

- **A.** The fitted twitches against the data.
- **B.** The Ca transients that produce them.
- **C.** The attachment gates.
- **D, E.** Fitted parameters and derived quantities relative to fit C (log axis; colours as in A).
- **F.** RMS error, with the data SEM dashed.

| alternative | RMS (% of peak) | τ (ms) | offset (ms) | `k_on` | `K_coop` | `k_off` (s⁻¹) | pCa₅₀ / n_H | σ | activation | gate peak | bound heads | SRX min | ATP/head/beat | diastolic ATP (s⁻¹) | **AUC × σ** | Ca→force lag (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **fit C** | **0.81** | 112/108 | −2.7 | 4.15 | 20.0 | 3.82 | 6.21 / 7.0 | 7.89 | 86 % | 0.76 | 0.234 | 0.041 | 2.90 | 0.139 | 22.9 | 120 |
| `Cadia` 0.05 µM | **0.74** | 112/110 | −3.7 | 4.15 | 19.4 | 3.72 | 6.22 / 6.9 | 7.76 | 86 % | 0.77 | 0.236 | 0.040 | 2.87 | **0.071** | 22.3 | 120 |
| `Cadia` 0.20 µM | 1.16 | 107/103 | +2.5 | 4.31 | 20.0 | 4.15 | 6.20 / 7.0 | 8.21 | 84 % | 0.74 | 0.230 | 0.042 | 2.98 | **0.311** | 24.4 | 118 |
| core cycling rates ×0.5 | **0.75** | 107/106 | −3.2 | **5.34** | 19.5 | **4.81** | 6.22 / 6.6 | 7.70 | 86 % | 0.80 | 0.250 | 0.033 | **1.53** | 0.080 | **11.8** | 135 |
| core cycling rates ×2 | 0.95 | 114/116 | −3.4 | 3.77 | 20.0 | 3.58 | 6.19 / 7.3 | 8.21 | 83 % | 0.72 | 0.220 | 0.049 | **5.51** | 0.243 | **45.2** | 112 |
| `K_coop` fixed at 10 | 2.00 | 98/100 | +14.2 | 6.00 | 10 (fixed) | 5.49 | 6.22 / 4.6 | 8.85 | 78 % | 0.66 | 0.214 | 0.047 | 2.74 | 0.259 | 24.2 | 112 |
| `k_off` fixed at 0.5 × fit C | 1.96 | 75/80 | +5.4 | 2.60 | 20.0 | 1.91 (fixed) | 6.31 / 5.2 | 9.10 | 75 % | 0.62 | 0.206 | 0.050 | 2.58 | 0.175 | 23.4 | 143 |
| `k_off` fixed at 2 × fit C | 1.30 | 141/139 | −4.2 | 7.27 | 20.0 | 7.65 (fixed) | 6.17 / 6.4 | 7.72 | 87 % | 0.79 | 0.239 | 0.039 | 2.99 | 0.137 | 23.1 | 88 |

("core cycling rates" scales `kah`, `kamh`, `ka`, `kd`, `k1`, `k_1`, `k2` together, as `updateRates`
does; σ refitted in every row.)

**Three alternatives fit as well as fit C or better: `Cadia` 0.05 µM (0.74 %), a core twice as
slow (0.75 %) and a core twice as fast (0.95 %).** All three are visually indistinguishable in
panel A. The differences in RMS (≤ 0.14 percentage points) are far below the systematic residual
that every fit shares (lag-1 autocorrelation 0.93, §3.2), so they do not rank the alternatives.
What changes between these equally good fits:

- **ATP cost: per head and per tissue.** A 2× slower core halves the ATP per beat (1.53 against
  2.90 per head; AUC × σ 11.8 against 22.9). A 2× faster core nearly doubles it (5.51; 45.2). The
  thin filament absorbs the difference: `k_on` and `k_off` rise 29 % and 26 % for the slower
  core. **The twitch therefore does not determine the energetics.** Tissue-level ATP is robust
  only to the Ca-layer assumptions (22–30 across every other fit in this report); it scales
  almost linearly with the core's cycling speed.
- **Diastolic ATP and diastolic activation.** These follow the assumed diastolic Ca: 0.071 s⁻¹ at
  0.05 µM against 0.311 s⁻¹ at 0.20 µM, a 4× range.
- **The thin-filament rates.** `k_off` ranges 3.6–4.8 s⁻¹ across the equally good fits, and its
  partner `k_on` ranges 3.8–5.3.

The **other three alternatives are clearly worse** (1.3–2.0 %):
- `K_coop` held at a more physiological 10 (n_H 4.6);
- `k_off` pinned at half or double fit C's.

Within the stated assumptions, the steep cooperativity and the order of magnitude of k_off are
therefore constrained. Even these worse fits keep the same qualitative state: 75–87 % activation,
gate peak 0.62–0.79, Ca peak 88–143 ms before force.

### 4.3 Local identifiability *given* all assumptions (Jacobian at fit C)

These intervals hold **only with every assumption of §2 fixed**; §4.2 shows how far the values move when an assumption changes. At fit C the residual is 2.3× the data SEM (fit B: 11×), so local intervals mean more than they
did for fit B. They are still optimistic: the residual keeps a lag-1 autocorrelation of 0.93
(§3.2). Condition number of J: 8 × 10³.

| parameter | fit C | 95 % interval (local) | strongest correlation |
|---|---|---|---|
| τ_rise / τ_decay | 112 / 108 ms | ±354 / ±345 ms | with each other, −1.00 |
| onset offset `t0` (vs recorded stimulus) | −2.7 ms | ±3.7 ms | `K_coop`, −0.69 |
| `tf_k_on` | 4.15 µM⁻¹s⁻¹ | ×/÷ 1.18 | `k_off`, **+0.98** |
| `tf_k_off` | 3.82 s⁻¹ | ×/÷ 1.21 | `k_on`, **+0.98** |
| `tf_K_coop` | 20.0 (bound) | ±2.9 | τ, ±0.85 |

- **The ratio `k_on/k_off` is sharply identified; each rate less so.** The two move together
  (r = 0.98), keeping the steady-state sensitivity. Their common scale sets the thin-filament
  speed, and that is still determined to about ±20 %.
- **The onset is identified** to a few ms (−2.7 ± 3.7 ms), even though the stimulus record is not
  synchronised. It is set by the twitch's own rise, not by the time stamp. A Ca onset just
  before the recorded stimulus is consistent with an unsynchronised record.
- `K_coop` is on its upper bound (20). The data want **at least** this much cooperativity (§5).


### 4.4 What can be trusted, and what not

| tier | quantity | evidence |
|---|---|---|
| **pinned by the data** (stable across every alternative of §4.2 and §A.3–A.4) | the attachment-gate time course (peak ~200 ms); the Ca peak at about 1.6 × Ca₅₀; **activation 75–87 %** of saturating force; a steep effective nonlinearity (n_H 4.6–7.3); a slow Ca→force step (effective thin-filament off-rate tens of s⁻¹); Ca peaking 90–145 ms before force; the onset offset (−4 to +14 ms) | same twitch, same gate (§4.2, panel C) |
| **fixed by an assumption the twitch cannot test** | **ATP per beat, per head and per tissue** (core speed: 12–45 for AUC × σ); diastolic activation and diastolic ATP (diastolic Ca: 4×); `k_on`, `k_off` to about ±30 % (core speed); how the Ca→force delay splits between Ca handling and the thin filament (fits B and C, §A.3) | §4.2, §A.3 |
| **not identifiable at all** | absolute Ca (diastolic and peak), `k_on` separately from Ca; which Ca time constant is the rise | exact symmetries, §4.1 |

This **supersedes** the earlier claim (§A.2) that tissue-level ATP per beat is robust. It is
robust to the Ca layer, but not to the speed of the cross-bridge core.

### 4.5 Why fit C is still only a working reference

1. **The core is an adult mouse skinned core** (§0). A core 2× slower or 2× faster fits the twitch
   equally well (§4.2), with a 4× different ATP cost. If the cells are immature, or β-MHC if
   hiPSC-derived, a slower core is plausible.
2. **Force cannot separate Ca kinetics from filament kinetics** (§A.3). Fit B (slow Ca, fast
   filament) and fit C (fast Ca, slow filament) drive nearly the same gate. Fit C wins partly
   because the bi-exponential cannot make a slow-rising, fast-falling transient.
3. **The Ca levels are a gauge choice** (§4.1), and the diastolic Ca is a free assumption (§4.2).
4. **The boundary is isometric; the EHT pulls on a compliant post.** This was probed only crudely
   (§6).
5. **The late tail is missed** (§3.2), so some process is missing altogether.

**Recommendation:** use fit C as the working reference for twitch *shape* and for
drug-mechanism *directions*. Do not quote its ATP numbers, per head or per tissue, as predictions
until the core speed is constrained independently (an O₂-consumption or ktr-type measurement on
the same tissue, §7). Treat k_off, n_H and the activation level as hypotheses for the
measurements of §7.

## 5. Plausibility: does anything cross the line?

| item | fit C | physiological reference | verdict |
|---|---|---|---|
| Ca transient: time to peak, decay τ | 110 ms after its onset, about 110 ms | hiPSC-CM about 50–250 ms and 100–400 ms (*estimate*) | **ok** (fast side of the decay range) |
| Ca at next stimulus | 0.104 µM (diastolic 0.10) | returns to diastole | **ok**: the reset artefact of fit B (§A.3) is gone |
| Ca onset vs recorded stimulus | −2.7 ± 3.7 ms (free offset) | the stimulus record is not synchronised, so any offset is admissible | **ok**, and identified |
| Ca₅₀ relative to peak Ca | 0.62 µM at a 1.0 µM peak (pCa₅₀ 6.21) | skinned hiPSC-CM pCa₅₀ about 5.9–6.3 (ssTnI, *estimate*); intact steady-state force–[Ca]ᵢ in rat trabeculae Ca₅₀ ≈ 0.6 µM (*literature*, recalled) | **ok** |
| steepness | **n_H 7.0**; `K_coop` 20.0, **on its bound** | skinned myocardium n_H 2–4; intact steady-state force–[Ca]ᵢ n_H about 4–6 (*literature*, recalled) | **crosses for skinned muscle, at or beyond the edge for intact muscle.** The bound is active, so the data want even more |
| `tf_k_off` | **3.82 s⁻¹** (effective 18–79 s⁻¹, τ 13–55 ms) | isolated cTnC 700–1000 s⁻¹; thin filament + cycling S1 ≈ 13 s⁻¹ at 15 °C (Davis 2007) | **crosses as a troponin rate at 37 °C; ok as a lumped** Ca→force activation step |
| activation per twitch | **86 %** of saturating force; Ca×1.5 → +13 % | adult myocardium twitches are well below maximal Ca-activated force (*estimate*); EHT inotropy is usually available ([Ca]ₒ, Iso) | **borderline**; testable with a [Ca]ₒ series (§7) |
| heads at peak | 23 % attached; ATPase at 87 % of saturating | | plausible |
| SRX | 0.20 → 0.04 every beat | mant-ATP SRX lifetimes of order 100 s (*literature*) | **not interpretable**: the model's SRX is a fast pool |
| diastolic ATP rate | 0.139 s⁻¹ | basal DRX turnover ≈ 0.03 s⁻¹ (*literature*) | **borderline high**; set by the fixed `Cadia` |
| core speed (fixed) | adult mouse skinned core | a 2× slower or faster core fits equally well (§4.2) | **not testable from force**; sets the ATP cost |
| force scale | σ 7.9 µN/kPa, 0.0079 mm² at full power | a few % of the adult core's stress per area for any plausible EHT cross-section (§3.5) | **ok** |
| late tail (> 550 ms) | missed (1–3 % of peak) | | a missing slow process (§3.2) |

**Summary:** no hard physical limit is crossed, but two thin-filament parameters stand in for
more than troponin:
- the n_H of 7.0 with `K_coop` on its bound;
- the 3.8 s⁻¹ k_off, which is plausible only as a lumped activation rate.

The high activation level and the diastolic ATP rate are borderline, and both depend on the
fixed Ca levels. The ATP cost depends on the core speed, which the twitch cannot test (§4.2).

## 6. Sensitivities at fit C

One change at a time from fit C. "Shape RMS" is the RMS error after refitting σ alone. ΔF and ΔATP
are relative to fit C (54.4 kPa; 2.90 ATP/head/beat).

| change | ΔF_dev | ΔATP per beat | stim→peak | 10 %→peak | RT50 | RT90 (ms) | shape RMS | SRX min |
|---|---|---|---|---|---|---|---|---|
| **fit C** | — | — | 228 | 175 | 153 | 255 | 0.8 % | 0.041 |
| `Cadia` 0.05 µM | 0 % | −6 % | 228 | 174 | 143 | 239 | 2.0 % | 0.041 |
| `Cadia` 0.20 µM | −1 % | +15 % | 230 | 180 | 174 | 290 | 4.3 % | 0.040 |
| Ca amplitude ×1.5 | **+13 %** | +44 % | 200 | 158 | 243 | 341 | 10.3 % | 0.032 |
| Ca amplitude ÷1.5 | **−58 %** | −61 % | 212 | 167 | 115 | 242 | 9.0 % | 0.112 |
| `k_off` ×2 alone (sensitivity halved) | −83 % | −85 % | 165 | 136 | 112 | 282 | 19.8 % | 0.165 |
| `k_on`, `k_off` ×2 (same steady state) | +7 % | −1 % | 188 | 145 | 151 | 240 | 11.1 % | 0.036 |
| `k_on`, `k_off` ÷2 (same steady state) | −26 % | −17 % | 268 | 207 | 164 | 302 | 11.6 % | 0.064 |
| `K_coop` 10 | −16 % | +1 % | 242 | 195 | 180 | 327 | 8.1 % | 0.053 |
| `kSE` ×0.1 (compliant series element) | −37 % | +10 % | 230 | 186 | 167 | 255 | 6.9 % | 0.086 |
| `ka` ×0.5 | −46 % | −55 % | 220 | 166 | 120 | 215 | 7.5 % | 0.093 |
| `k2` ×0.5 | +33 % | −17 % | 260 | **206** | 164 | 278 | 8.0 % | 0.015 |
| `k1` ×0.5 | −10 % | −21 % | 225 | 173 | 152 | 252 | **1.2 %** | 0.046 |
| `kmsrd` ×3 | +1 % | +4 % | 225 | 175 | 158 | 263 | 1.3 % | 0.032 |
| `ksr0` ×3 | −2 % | −4 % | 230 | 177 | 148 | 248 | 1.0 % | 0.086 |

- **The thin filament's speed is now a real lever.** Co-scaling `k_on` and `k_off` by 2 either way
  (same steady state) shifts the peak by −40 / +40 ms and costs 11 % RMS. That is why its
  common scale is identified (§4.3).
- **Near saturation, Ca amplitude acts asymmetrically:** ×1.5 gives only +13 % force, ÷1.5 gives
  −58 %. This is the §7 prediction at fit C. The Iso wells' +15 % force (Prism) is consistent
  with it, though PKA also desensitizes the filament.
- **The cross-bridge core matters one rate at a time, but not jointly.** `k2` ×0.5 alone slows the
  rise by 31 ms and `ka` ×0.5 alone shortens RT90 by 40 ms, each costing about 7–8 % RMS. Scaling
  *all* cycling rates together and refitting the Ca layer (§4.2) fits as well (0.75 % at ×0.5,
  0.95 % at ×2). The thin filament absorbs the change, while the ATP cost halves or doubles.
- **Changes that leave the shape alone:** `k1` ×0.5 (Pi release) keeps the shape (1.2 % RMS) but
  saves 21 % ATP for 10 % force. The SRX levers (`kmsrd`, `ksr0`) barely matter at this
  activation.
- **Compliance:** a softer SE lowers force 37 % and *raises* ATP 10 %. Economy on compliant posts
  is not isometric economy.

## 7. What would make the fit unique: measurements and questions

### 7.1 Measurements that break the non-uniqueness

Each measurement below removes one non-uniqueness of §4. With a fast filament the internal state could flip from "few heads, foot of the curve" to "almost all heads, saturation" (§A.2); the first two observations separate these regimes, and **neither needs absolute Ca:**

1. **Relative Ca transient** in the same wells. A normalized recording fixes the input shape:
   fed in as the input (`CaInput = 'signal'`, already implemented), it removes both τ, the
   bi-exponential shape-family error and the Ca-versus-filament trade-off (§A.3) at once. It
   leaves `k_on·ΔCa`, `K_coop`, `k_off` and the diastolic-to-peak Ca ratio, and these are then
   constrained by force. Two signatures discriminate directly:
   - **the Ca-peak-to-force-peak lag:** 120 ms in fit C, 88–143 ms across the alternatives of
     §4.2, but only about 30 ms in fit B (fast filament). The lag does not need the stimulus to
     be synchronised;
   - **the force–Ca loop** (§3.1, panel I): concave and flat-topped near saturation (fit C),
     convex and power-law on the foot of the curve (fit A).
   - *Caveat:* with a single-wavelength dye (Fluo-4, Cal-520, Kd ≈ 0.3–0.35 µM, *estimate*),
     ΔF/F₀ is not linear in Ca at 1 µM, so it compresses the peak. A low-affinity dye
     (Fluo-4FF, Cal-520FF, Fura-FF) or a ratiometric one (Fura-2) keeps the relative signal
     usable. GCaMP lines have slow kinetics that would bias the time course.
2. **Ca amplitude changes without calibration:** an extracellular [Ca] series, e.g.
   0.5–3.6 mM, or a Ca-channel blocker at low dose (the diltiazem plate already exists). The last
   row of the §A.2 table and §6 give sharp predictions: **amplitude ×1.5 changes force by only
   +13 % at fit C (and ÷1.5 by −58 %)**, against +85 to +322 % in the low-activation regime. A force–[Ca]ₒ curve that is
   already flat at the standard medium [Ca] would place the tissue near saturation. The mapping
   [Ca]ₒ → transient amplitude is itself nonlinear, which is why a relative Ca recording in
   parallel helps.
   - A hint already in the data: isoproterenol raised force by only +15 % (Prism). In the fit A
     regime, a larger Ca transient would roughly double it. That points towards saturation, but
     PKA also desensitizes the myofilaments, so it is not decisive.
3. **Skinned-EHT or myofibril force–pCa from the same cells.** This gives absolute pCa₅₀, n_H
   and maximal active stress per area. Together with an absolute Ca amplitude it fixes the
   operating point outright. Maximal stress per area also turns σ into a measured
   XB-power scale.

4. **An energetic or kinetic anchor for the core speed** (new, from §4.2). Paced O₂ consumption
   of the same tissues, or a ktr / force-redevelopment measurement on skinned EHT, would fix the
   cross-bridge cycling speed. Only that turns the model's ATP-per-beat into a prediction.

### 7.2 Questions for the experimenters

Ordered roughly by value per effort.

**Data they probably already have**

1. **The platemap** for Plate 6: what is in columns 1–3, is row D of columns 4–6 DMSO, the DMSO
   %, and when dosing happened relative to the recording.
2. **Pre-dose baseline recordings** of every plate. The Prism "% change from baseline" implies
   they exist. They would give many more control tissues, paired pre/post waveforms for each drug
   (mavacamten, isoproterenol, diltiazem), and day-to-day variability.
3. **Raw Pulse3D/h5 exports** with post displacement and the absolute baseline. The "active
   twitch force" is baseline-subtracted, which removes diastolic (passive) tension.
4. **Post stiffness** in µN/µm for the "Skeletal Muscle (12x)" posts, and how force is computed
   from magnet displacement. This gives tissue shortening during the twitch, so the auxotonic
   boundary can be modelled (§6, §4.5).
5. **Recording conditions:** medium temperature (37 °C?), [Ca]ₒ and composition, time since
   medium change, and the stimulation settings (the "40" amplitude units, capture threshold).
6. **Tissue:** cell line (what "GeDi lines" means), differentiation and tissue age, CM purity, cell number per tissue, any
   cross-section or width/thickness measurement, resting length, and sarcomere length if it
   was ever imaged.

**New experiments, ranked**

1. **Relative Ca transient in the same wells at 1 Hz** (a low-affinity or ratiometric dye).
   Replaces the Ca shape assumption (§7.1) and separates Ca kinetics from thin-filament
   kinetics (§A.3). Sharp predictions, stated relative to the force peak because the stimulus
   record is not synchronised: fit C puts the Ca peak **120 ms before the force peak**, fit B only
   about 30 ms before. The force–Ca loop of fit C is the one in §3.1, panel I.
2. **Force–[Ca]ₒ series** (e.g. 0.5, 1, 1.8, 3.6 mM) and/or a **pacing-frequency series**
   (0.5–2 Hz, plus post-rest potentiation). Needs no dye. The [Ca]ₒ series discriminates the
   operating point directly (§7.1, item 2). Fit C predicts saturation: a 1.5× larger transient gives
   only +13 % force, a 1.5× smaller one −58 % (§6). The frequency series constrains Ca-handling time constants
   and the between-beat Ca carry-over that the current transient gets wrong.
3. **Ca calibration** (Fura-2 with ionomycin Rmin/Rmax, or a calibrated FF dye). Breaks the
   `k_on·Ca` degeneracy, so pCa₅₀ becomes absolute.
4. **Skinned EHT or myofibril mechanics from the same line:** force–pCa, maximal stress, kACT and
   kREL. Pins the thin filament and the XB power and kinetics (the `ka`/`k2` levers of §6).
5. **Oxygen consumption of paced EHTs** (Seahorse or a closed chamber) at 1 Hz, and with
   mavacamten. The twitch cannot fix the core speed, and the core speed sets the ATP cost (a 4×
   range for equally good fits, §4.2). A measured ATP per beat would pin it.
6. **SRX in the EHT's cardiomyocytes** (mant-ATP chase). The model's diastolic SRX fraction (0.20) and its
   exchange time scale were never identified. Depending on the operating point, the twitch
   either leaves SRX alone or empties it (0.20 → 0.03).
7. **Mavacamten, isoproterenol and diltiazem with waveforms, and ideally with Ca.** Mavacamten
   acts on the cross-bridge side, diltiazem on Ca amplitude, and isoproterenol on both. With Ca
   recorded, each drug separates one model layer.

### 7.3 Next steps (modelling)

1. **Ca input:** a transient with separate rise and decay shapes (or a measured relative
   transient via `CaInput = 'signal'`), superposed across beats.
2. **Core kinetics:** a core-rate multiplier is **not identifiable from the twitch** (§4.2:
   ×0.5 and ×2 fit as well). Fit one only together with an anchor for the core speed, such as
   O₂ consumption, ktr or myofibril kinetics (§7.1, item 4).
3. **Boundary:** the post stiffness from Travis, as a real series compliance.
4. **Hierarchy:** a shared shape with per-well σ, then the Iso wells (Ca side + `tf_k_off`), then
   mavacamten from the Prism % changes.

## 8. Is this the MyoSim thin filament? Yes in its kinetics, not in how overlap enters

Script: [`FitMyoSimEHTTwitch.m`](FitMyoSimEHTTwitch.m). Results: `results/eht_myosim.mat`,
`results/eht_myosim.png`; refits cached in `results/eht_myosim_fits.mat`. Model change:
`CaModel = 'myosim'` in [`Model/calciumActivation.m`](../../Model/calciumActivation.m), checked by
T7 of [`TestCaActivation.m`](TestCaActivation.m).

### 8.1 The two forms

The reference is MATMyoSim (`../MATMyoSim/code/@half_sarcomere/update_3state_with_SRX.m`
lines 77–96). The project's Beard port, `update_beard_atp.m` lines 173 and 186–190, uses the same
equations. MyoSim's state is `N_on`, the permissive fraction of all sites; `twoState`'s is
`N = 1 − N_on`.

| | MyoSim (`'myosim'`) | this report's `'twoState'` |
|---|---|---|
| ON flux | `k_on·Ca·(N_ov − N_on)·(1 + K_coop·N_on/N_ov)` | `k_on·Ca·N·(1 + K_coop(1 − N))` |
| OFF flux | `k_off·(N_on − N_bound)·(1 + K_coop·(N_ov − N_on)/N_ov)` | `k_off·(1 − N − bound)·(1 + K_coop·N)` |
| sites that can switch on | **only those in the overlap zone** (`N_ov` = 0.85 at SL 2.0) | all sites |
| cooperativity measured relative to | the overlap zone | all sites |
| attachment ∝ | free ON sites, `N_on − N_bound` | `(1 − N − bound)/(1 − bound)` × the RHS overlap factor `(N_ov − bound)/unattached` |
| gate implemented here | `(N_on − bound)/(N_ov − bound)`, so that with the RHS overlap factor the attachment is MyoSim's `(N_on − N_bound)` | `(1 − N − bound)/(1 − bound)` |
| bound sites | all attached heads (`p1 + p2`, one site per head) | the same |

**At full overlap (`N_ov` = 1) the two are the same model**, term for term. They differ only
where overlap is partial. With `UseOverlapFactor` off, the `'myosim'` gate is `(N_on − bound)/N_ov`,
so the attachment is again MyoSim's site term.

### 8.2 Checks

| check | result |
|---|---|
| `'myosim'` against `'twoState'` with overlap switched off (Rice twitch, fit C parameters) | max relative difference **4e-5** (solver tolerance) |
| saturating Ca (pCa 4.51, `fastgate_best`) | none 78.84, twoState 77.83, **myosim 77.84 kPa**: the calibration is kept |
| SL 2.0 with fit C's parameters, unchanged | twitch +2.9 %, shape differs by 3.25 % RMS, gate peak 0.80 against 0.76 |
| `CaModel = 'none'` path | untouched (the RHS change only keeps the geometric overlap in a new variable); T2–T6 of `TestCaActivation.m` pass unchanged (T1, the byte-identity check, needs `*_ref.m` copies of HEAD and was skipped) |
| MEX + `'myosim'` | rejected by the existing guard |

### 8.3 Fit M: fit C refitted in the MyoSim form

All six parameters were refitted with `'myosim'` from three starts (fit C's values; `k_on` and
`k_off` ×1.5; and τ ±10 ms with `K_coop` −3). The starts agree (SSR 0.00642–0.00654).

![MyoSim form](results/eht_myosim.png)

| | fit C (`twoState`) | **fit M (`'myosim'`)** | change |
|---|---|---|---|
| RMS (% of peak) | 0.81 | **0.82** | |
| stimulus→peak / 10 %→peak / RT50 / RT90 (ms) | 228 / 175 / 153 / 255 | 228 / 175 / 153 / 254 | |
| τ_rise / τ_decay (ms) | 112 / 108 | 109 / 111 | (exchangeable) |
| onset offset (ms) | −2.7 | −3.9 | −1.2 ms |
| `tf_k_on` (µM⁻¹s⁻¹) | 4.146 | 4.207 | +1.5 % |
| `tf_k_off` (s⁻¹) | 3.824 | 4.085 | +6.8 % |
| `tf_K_coop` | 19.97 (bound) | 19.99 (bound) | |
| pCa₅₀ / n_H (steady, SL 2.0) | 6.21 / 7.0 | 6.21 / 7.1 | |
| σ (µN/kPa) | 7.89 | 7.85 | −0.5 % |
| activation, gate peak, bound heads | 86 %, 0.760, 0.234 | 86 %, 0.766, 0.235 | |
| **permissive sites at peak** | 0.82 (of all sites) | **0.70** (of all sites; at most 0.85 = `N_ov`) | |
| ATP/head/beat (diastolic; saturating) | 2.90 (0.139; 11.6 s⁻¹) | 2.91 (0.133; 11.6 s⁻¹) | |
| AUC × σ | 22.9 | 22.8 | |
| Ca peak → force peak | 120 ms | 123 ms | |
| local 95 %: log₁₀ `k_on`, log₁₀ `k_off`, `K_coop`, offset | ±0.071, ±0.083, ±2.9, ±3.7 ms | ±0.068, ±0.089, ±2.1, ±3.9 ms | r(`k_on`, `k_off`) 0.98 → 0.99 |

**Conclusion.** The MyoSim form fits the twitch exactly as well, and every observable and every
energetic quantity is unchanged to within 1 %. The thin-filament constants move by less than 7 %,
well inside the non-uniqueness of §4. The visible difference is internal bookkeeping:
- In the MyoSim form only overlap sites switch on, so the permissive fraction of *all* sites
  peaks at 0.70, not 0.82.
- The attachment gate, which is what force sees, is the same (panel D).

Either form can carry fit C's results. **Use `'myosim'` when parameters or states are to be
compared with MyoSim or FiberSim models, or when overlap will change** (length steps, work loops,
Frank–Starling). The two forms differ exactly there, and only the MyoSim form keeps sites outside
the overlap from switching on.

### 8.4 Against the published model (Campbell, Janssen & Campbell 2018)

Checked against the paper's PDF: *Biophys J* 115:543–553 (2018), doi 10.1016/j.bpj.2018.07.006. The
parameter tables (S1–S3) are in its online supplement, which was not available. Values below are
from the main text only.

| part | Campbell 2018 (equation) | this model (fit M) | same? |
|---|---|---|---|
| thin filament ON | `J_on = k_on[Ca](N_overlap − N_on)(1 + k_coop·N_on/N_overlap)` (Eq. 1) | identical (`'myosim'`) | **yes** |
| thin filament OFF | `J_off = k_off(N_on − N_bound)(1 + k_coop(N_overlap − N_on)/N_overlap)` (Eq. 2) | identical | **yes** |
| initial state | all sites off, all heads OFF | all sites off (`N_on` = 0); heads from the core's own initial state | partly |
| attachment | `J3 = k3·exp(−k_cb x²/2k_BT)·M_ON·(N_on − N_bound)` (Eq. 5) | `ka·PD·(N_on − bound)/(1 − bound)` × the core's strain kernel; attaches from the primed M·ADP·Pi pool only | site term yes; extra `1/(1 − bound)` from the core's overlap factor (≈ +30 % at the twitch peak) |
| detachment | `(k4,0 + k4,1(x − x_ps)⁴)·M_FG` (Eq. 6) | the core's strain-dependent multi-state cycle | no (different core) |
| OFF → ON | `(k1 + H·k_MLCP)(1 + k_force·F_total)·M_OFF`, **linear** in force (Eq. 3) | `kmsrd·exp(F_total/24.5 kPa)` from SRX·ADP·Pi | force-dependent in both; linear vs exponential |
| ON → OFF | `k2·M_ON`, constant (Eq. 4) | `ksr0·exp(−F_total/26.6 kPa)` into SRX·ATP, suppressed by force | **no** |
| OFF pools | one | two (SRX·ATP, SRX·ADP·Pi) linked by hydrolysis | no |
| force for the thick filament | `F_total = F_active + F_passive` | `F_total` (passive included) | **yes** |
| series elasticity in the twitch fit | yes (`k_series`) | yes (`kSE`) | yes |

**Parameter values (main text):**

| | Campbell 2018 | fit M |
|---|---|---|
| `k_on` | 2.08 × 10⁷ M⁻¹s⁻¹ = **20.8 µM⁻¹s⁻¹** | 4.21 µM⁻¹s⁻¹ |
| `k_off` | **100 s⁻¹**, fixed (from Siddiqui et al.) | 4.08 s⁻¹ |
| `k_on/k_off` | 0.21 µM⁻¹ | 1.03 µM⁻¹ (≈ 5× more sensitive) |
| `k_coop` | **5.70** | 19.99 (at bound) |
| `k_force` | 95 % CI 8.8 × 10⁻⁵ – 2.1 × 10⁻⁴ Pa⁻¹ = 0.09–0.21 kPa⁻¹ | — (exponential form; ×9 at the twitch's 54 kPa, against ×5.8–12 for the linear form) |
| `k2` (ON → OFF) | 200 s⁻¹, fixed | 345 s⁻¹ × exp(−F/26.6) |

Their thin-filament constants come from **blebbistatin-treated, permeabilized rat myocardium**
(troponin-probe orientation, pCa 6.2–4.2). Their twitch fit (rabbit trabecula, 37 °C) re-adjusted
11 parameters (Table S3, not seen here) and was **driven by measured, smoothed bis-fura-2 Ca
transients**.

**What this means for fit M:**

1. **The Ca layer is the published Campbell 2018 thin filament,** term for term.
2. **Fit M's thin filament is far from their permeabilized-myocardium values:** 25× slower `k_off`,
   5× higher sensitivity relative to the transient, 3.5× stronger cooperativity. Two caveats:
   - Our Ca levels are a gauge (§4.1), so only `k_on/k_off` relative to the peak Ca compares, not
     `k_on` itself.
   - Their twitch-fit values (Table S3) may differ from the permeabilized ones; the supplement is
     needed to compare like with like.
3. **The thick filament is this project's own, not theirs.** Both are force-driven with total
   force, and at the twitch peak the recruitment boost is of similar size (×9 against ×6–12). But ours is
   exponential, has force-suppressed re-entry and two OFF pools.
4. **Methodologically, they did what §7 recommends:** they drove the twitch with a measured Ca
   transient. That removes the Ca-shape and Ca-versus-filament non-uniqueness (§4, §A.3) that fit M
   still carries.

## Appendix A. How the fit evolved (fits A → B → C)

These sections document the path to fit C. Their conclusions are superseded where §3–4 say so.

### A.1 Fits A and B, and what separated them

The two fits are **not two physiological hypotheses.** They share every assumption (below) and
differ in one thing only: the upper bound on `tf_k_on`.

- Fit A used the bound in `parameterBounds.m` (≤ 200 µM⁻¹s⁻¹, from isolated-cTnC kinetics).
- Fit B released it (≤ 3000).

Because only `k_on·Ca` matters (§4.1), fit A's bound was in effect an arbitrary cap on
sensitivity *relative to a Ca amplitude I fixed myself* (1.0 µM). It was the natural first test,
not a physiological constraint. The shared assumptions are also constraints on the answer:

| assumption | value | could it move the answer? |
|---|---|---|
| Ca shape | bi-exponential, periodized with `mod(t, T)` | **yes**: the shape family is too rigid (§4.1) and sets how saturation is used |
| Ca levels | 0.10 → 1.0 µM | only through `k_on·Ca` (exact) and the diastolic/peak ratio (weak, §6) |
| bounds | τ_rise ≤ 150 ms, `K_coop` ≤ 20, `t0` ≤ 60 ms | **yes**: fit B sits on the first two |
| cross-bridge core | adult mouse, 8 mM ATP, unchanged | yes: at high activation `ka` and `k2` shape the twitch (§6) |
| mechanics | isometric, fixed ends, model SE, SL 2.0 µm | yes: the posts are compliant |
| `tf_k_off` | 724 s⁻¹ | **yes**, dynamically: freeing it gives fit C (§A.3) |


![fits A and B](results/eht_fit.png)

Fits A and B side by side (`FitEHTTwitch.m`): A data and both fits; B Ca transients; C ATP turnover; D gate and bound heads; E force–pCa; F shape sensitivities at fit B.

### A.2 Profile over sensitivity with a fast filament

> **Read with §A.4.** This profile holds `tf_k_off` at 724 s⁻¹. With `k_off` free (§A.3) the
> low-activation regime is no longer competitive: starts from it return to fit C. The tiers below
> describe what force alone could say within the fast-filament family. Their claim that ATP × σ is
> robust holds only for the Ca layer; it scales with the core speed (§4.2, §4.4).

Script: [`ProfileEHTTwitch.m`](ProfileEHTTwitch.m). Results: `results/eht_profile.mat`,
`results/eht_profile.png`, fit cache `results/eht_profile_fits.mat`.


`tf_k_on` was fixed on 9 values from 158 to 1259 (the rest refitted). Two further fits fixed
`K_coop` at 9.68 (Marzban) and 5 (n_H ≈ 3) and freed `k_on`. Each fit used **3 starts** (between
A and B, from B's other parameters, from A's) and kept the best. The starts agree to within
0.02 SSR, so this is the profile of this model, not an optimizer artefact. Fit A's single start
had stopped slightly short: the profile at `k_on` = 200 reaches SSR 0.668, against fit A's 0.686.

![profile](results/eht_profile.png)

| `k_on` | 158 | 200 | 251 | 316 | 398 | 501 | **658 (B)** | 891 | 1259 | `K_coop` 9.68 | `K_coop` 5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| RMS error (% of peak) | 8.5 | 8.3 | 8.1 | 7.7 | 7.2 | 6.5 | **3.8** | 4.9 | 6.1 | 5.0 (`k_on` 765) | 6.3 (`k_on` 834) |
| pCa₅₀ | 5.53 | 5.62 | 5.72 | 5.82 | 5.92 | 6.02 | 6.15 | 6.26 | 6.43 | 6.21 | 6.25 |
| n_H | 5.1 | 6.8 | 5.3 | 6.8 | 5.4 | 5.8 | 5.9 | 5.3 | 5.1 | 4.4 | **2.8** |
| `K_coop` | 17.9 | 19.9 | 20.0 | 19.7 | 20.0 | 15.1 | 20.0 | 14.9 | 12.9 | 9.68 | 5 |
| `t0` (ms) | 58 | 55 | 52 | 44 | 23 | **0** | 19 | 54 | 60 | 50 | 60 |
| τ_rise / τ_decay (ms) | 117/118 | 121/121 | 121/133 | 138/138 | 150/175 | 150/239 | 150/176 | 115/115 | 93/94 | 125/126 | 108/110 |
| σ (µN/kPa) | 244 | 189 | 127 | 76 | 38 | 12.9 | 7.4 | 6.8 | 6.5 | 7.6 | 8.1 |
| activation (F_dev / F_sat) | 3 % | 4 % | 6 % | 9 % | 19 % | 57 % | 88 % | 94 % | 96 % | 87 % | 84 % |
| gate peak | 0.02 | 0.02 | 0.03 | 0.06 | 0.11 | 0.41 | 0.81 | 0.89 | 0.93 | 0.79 | 0.75 |
| bound heads, peak | 0.008 | 0.010 | 0.016 | 0.026 | 0.052 | 0.156 | 0.244 | 0.258 | 0.264 | 0.242 | 0.233 |
| SRX minimum (diastole 0.20) | 0.20 | 0.19 | 0.19 | 0.18 | 0.15 | 0.07 | 0.04 | 0.03 | 0.03 | 0.04 | 0.04 |
| ATP per head per beat | 0.12 | 0.16 | 0.22 | 0.36 | 0.69 | 2.0 | 3.2 | 3.4 | 3.7 | 3.2 | 3.3 |
| **ATP × σ** (per tissue, arbitrary units) | **29.8** | 29.2 | 28.3 | 27.5 | 26.5 | 25.7 | 23.5 | 23.3 | 24.2 | 24.4 | 26.6 |
| force change for Ca amplitude ×1.5 / ÷1.5 | +85 / −41 % | +105 / −44 % | +146 / −48 % | +276 / −54 % | +322 / −65 % | +57 / −74 % | **+9 / −59 %** | +3 / −12 % | +2 / −4 % | | |

#### How much to trust what (as concluded at the time; see §4.4)

**The data are far more precise than the model.** The control template's standard error is
0.35 % of peak (well-to-well shape SD 1.3 %, n = 15). Fit B's residual is 3.8 %, **11× larger**,
with lag-1 autocorrelation 0.94. The misfit is systematic: it is a model-structure error, not
noise. Two consequences:

- Fit B's Jacobian intervals (e.g. `k_on` ± 11 %; §4.3 now reports fit C's) are formally meaningless. They assume the
  residual is noise.
- The sharp profile minimum at `k_on` ≈ 660 (RMS 3.8 % against 4.9–6.5 % either side) says which
  **deficiency of the model is least costly**, not where the truth is. The minimum is exactly where
  τ_rise and `K_coop` both sit on their bounds. So its location is partly set by those bounds and
  by the bi-exponential shape. Moving either bound was not tested.

Sorted by how far they can be trusted:

| tier | quantity | why |
|---|---|---|
| **robust** (little model dependence) | the twitch features themselves; **tissue-level ATP per beat (ATP × σ)**, which varies only **±12 %** (23.3–29.8) while σ varies 37× and per-head ATP 30× | force and ATP both scale with the number of cycling heads; this holds as long as the core's ATP-per-force-time is right |
| robust, qualitatively | the twitch needs a **steep effective Ca→force nonlinearity**: every profile point settles at n_H 5–7, and forcing n_H to 2.8 (`K_coop` 5) costs RMS 3.8 → 6.3 % | the rise-plateau-sigmoid-fall shape |
| **conditional** (set by model structure) | the operating point: `k_on·Ca_max` relative to the Ca₅₀, i.e. how far up the force–pCa curve the peak Ca goes; all reported pCa₅₀ values are only meaningful at Ca_max = 1 µM | the best fit prefers 88 % activation, but at 8 % RMS every value from 3 % to 96 % gives a twitch that looks roughly right |
| conditional, and it follows the operating point | σ (37×), per-head ATP (30×), gate peak (0.02–0.93), bound heads (0.008–0.26), SRX recruitment (none to near-total) | each is a monotone function of the operating point (panels C, D) |
| **not identified** | `t0` (0–62 ms across the profile), τ_rise vs τ_decay, absolute Ca, `tf_k_off`, `Cadia` | trade-offs (§4) |

The robust ATP × σ result is **not** an ATP measurement. It is a model prediction whose scale
comes from the mouse core's ATP per force-time. Its value is that it does not depend on the
unidentified operating point. **Superseded in part by §4.2:** it does depend on the core speed.

### A.3 `tf_k_off` is redundant with `tf_k_on` only in steady state

> This section records how fit C was found (K3, with the onset bounded at ≥ 0). The final fit C
> frees the onset offset (§A.4) and differs from K3 by < 2 % in everything downstream.

Script: [`TestKoffEHTTwitch.m`](TestKoffEHTTwitch.m). Results: `results/eht_koff.mat`,
`results/eht_koff.png`; fit caches `results/eht_koff_K2.mat`, `results/eht_koff_bound.mat`.

The first draft of §6 (at fit B) called `k_off` "fully redundant with `k_on`". That was a steady-state statement. Scaling both by
the same λ keeps every steady state exactly (`Jon = Joff` is invariant), but it changes how fast the
thin filament follows Ca.

#### Co-scaling at fixed shape parameters

λ was scaled from 3 to 0.01, with σ refitted and nothing else. There were 3 beats, so slow
deactivation could carry over into the next beat.

| λ (k_off, s⁻¹) | 3 (2172) | 1 (724) | 0.3 (217) | 0.1 (72) | 0.03 (22) | 0.01 (7) |
|---|---|---|---|---|---|---|
| fit B: RMS | 3.8 % | 3.8 % | 3.7 % | 3.9 % | 6.4 % | 15.4 % |
| fit B: 10 %→peak / RT50 / RT90 (ms) | 140 / 172 / 314 | 140 / 173 / 314 | 140 / 175 / 315 | 143 / 174 / 312 | 151 / 178 / 308 | **188** / 172 / 295 |
| fit B: force at 700 ms (% of peak) | 2.3 | 2.3 | 2.3 | 2.3 | 2.4 | 2.8 |
| fit B: diastolic force (kPa) / gate at stimulus | 0.930 / 0.0060 | 0.930 / 0.0066 | 0.930 / 0.0063 | 0.932 / 0.0068 | 0.935 / 0.0068 | 0.950 / 0.0072 |
| A side (`k_on` 200): RMS | 8.3 % | 8.3 % | 8.3 % | 8.3 % | 8.2 % | 8.7 % |

- **The diastolic phase does not see k_off**, even at 7 s⁻¹. The tail, the diastolic force and
  gate, and beat-to-beat carry-over (drift < 1e-3 kPa) are all unchanged. The site off-flux is
  `k_off·Pu·(1 + K_coop·N)`, so the effective rate grows about 20-fold as the filament switches
  off. It is fastest exactly in diastole.
- **The rise does see it.** A slow filament delays activation (10 %→peak 140 → 188 ms). The rise
  was fit B's largest misfit, so this suggested refitting.
- On the A side (foot of the curve) nothing changes. There the gate is small and follows Ca
  proportionally.

#### Refit with k_off free: fit C

Fit B was refitted with `tf_k_off` as a 6th free parameter (log scale). The two starts kept fit B's
`k_on/k_off` at λ = 0.1 and 0.03. Both converged to SSR 0.0084 on the lower bound k_off = 5 s⁻¹.
With the bound lowered to 0.5 s⁻¹, the refit moved inside the range: **fit C**.

| | fit B | **fit C** | data |
|---|---|---|---|
| SSR / RMS | 0.138 / 3.8 % | **0.0073 / 0.87 %** | SEM 0.35 % |
| `tf_k_off` / `tf_k_on` | 724 (fixed) / 658 | **3.8 s⁻¹ / 4.2** | |
| `k_on/k_off` (steady sensitivity) | 0.91 | 1.10 | |
| `K_coop`, `t0` | 20 (bound), 19 ms | 18.9, **0 (bound)** | |
| τ_rise / τ_decay; Ca time to peak | 150 / 176 ms; 180 ms | 106 / 108 ms; **108 ms** | |
| Ca at 950 ms (diastolic 0.10 µM) | 0.147 µM | 0.103 µM (the reset artefact is gone) | |
| stimulus→peak / 10 %→peak / RT50 / RT90 (ms) | 208 / 140 / 173 / 314 | **228 / 174 / 153 / 256** | 230 / 176 / 152 / 243 |
| pCa₅₀ / n_H (steady) | 6.15 / 5.9 | 6.22 / 6.8 | |
| activation, gate peak, bound heads peak | 88 %, 0.81, 0.24 | 85 %, 0.76, 0.23 | |
| SRX diastole → minimum | 0.20 → 0.04 | 0.20 → 0.04 | |
| σ; ATP per head per beat; ATP × σ | 7.4; 3.16; 23.5 | 7.9; 2.90; **23.0** | |
| effective site off-rate over the beat | 2900–15000 s⁻¹ | **17–76 s⁻¹** | |
| periodicity (5 beats) | — | drift 5e-4 kPa | |

![k_off](results/eht_koff.png)

- **Panels A–F:** the co-scaling scan (§A.3).
- **Panel G:** fit C against fit B.
- **Panel H:** Ca transients (solid/dashed) and attachment gates (dotted).
- **Panel I:** the effective site off-rate `k_off(1 + K_coop·N)` over the beat.

**What fit C shows:**

1. **A slow, cooperative thin filament reproduces the twitch almost to noise level.** The residual
   is 2.5× the data SEM, down from 11×. The shoulder-then-accelerating relaxation comes from the
   cooperative off-rate: about 17 s⁻¹ while the filament is fully on, about 76 s⁻¹ once it is off
   (panel I).
2. **The gate is identified; how it arises is not.** Fits B and C drive almost the same attachment
   gate (panel H, dotted) from very different inputs: B uses a slow Ca transient (peak 180 ms) and
   an instantaneous filament; C uses a fast Ca transient (peak 108 ms) and a filament with a
   13–60 ms time constant. **Ca kinetics and thin-filament kinetics trade off.** The force record
   constrains only their combination, and neither the steady-state constraint nor diastole
   separates them.
3. **Everything downstream of the gate is unchanged:** activation (85 % against 88 %), bound heads,
   SRX recruitment, σ, and tissue ATP × σ (23.0 against 23.5). This is the §A.2 pattern again:
   quantities set by the gate are as (un)certain as the operating point, while tissue-level ATP is
   robust to the Ca layer (but not to the core speed, §4.2).
4. **Plausibility.** Fit C's Ca transient (peak at 108 ms, decay τ about 108 ms) is closer to
   typical hiPSC-CM values than fit B's (*estimate*). A 3.8 s⁻¹ cTnC-level k_off is far below
   isolated-cTnC values (700–1000 s⁻¹). But thin-filament Ca dissociation *with* cycling myosin has
   been measured at about 13 s⁻¹ at 15 °C (Davis 2007, in `parameterBounds.m`), and the effective
   rate here is 17–76 s⁻¹. So a slow filament is not implausible. In a lumped two-state model,
   however, `k_off` absorbs **any** slow step between Ca and force (tropomyosin shift, cooperative
   spreading, thick-filament activation). Read it as "the Ca→force coupling has a 15–60 ms time
   constant", not as a troponin rate.
5. Limits: `t0` is on its lower bound (0) and `K_coop` is near 20. The basin search of §A.4 later
   showed that starts from the A side return to fit C. A full k_off profile was not run.

**Consequence for §7:** a relative Ca transient is now even more valuable. With the measured Ca
shape as input, the Ca-versus-filament trade-off disappears, and `tf_k_off` (a rate, in s⁻¹,
independent of the `k_on·Ca` degeneracy) becomes identifiable from force. Fits B and C make
different, testable predictions for it: Ca peaking about 30 against 120 ms before the force peak.

### A.4 Basin search, and the free onset offset

Two sets of refits, all with `k_off` free:

- **Low-activation and generic starts** (with the onset still bounded at ≥ 0). Two start from
  the low-activation side of the §A.2 profile: its `k_on` 200 and 398 points, with `k_on` and
  `k_off` both divided by 50 (the same steady state with a slow filament). One is generic.
- **Free onset offset** (allowed −150 to +60 ms, because the stimulus record is not
  synchronised). Starts from the previous fit C at offsets 0, −30 and −80 ms. The best of these
  is **fit C** throughout this report.

| start | SSR | τ (ms) | onset offset | `k_on` | `K_coop` | `k_off` (s⁻¹) | σ | gate peak | ATP/head/beat | AUC × σ |
|---|---|---|---|---|---|---|---|---|---|---|
| previous fit C (K3, §A.3; offset ≥ 0) | 0.00725 | 106/108 | 0 (bound) | 4.24 | 18.9 | 3.84 | 7.9 | 0.76 | 2.90 | 23.0 |
| from `k_on` 200 (offset ≥ 0) | 0.00741 | 108/107 | 0 (bound) | 4.31 | 18.6 | 3.94 | 7.9 | 0.76 | 2.90 | 23.0 |
| from `k_on` 398 (offset ≥ 0) | 0.00762 | 108/110 | 0 (bound) | 4.44 | 18.4 | 4.08 | 7.9 | 0.76 | 2.90 | 23.0 |
| generic (offset ≥ 0) | 0.00673 | 104/105 | 0 (bound) | 3.98 | 20.0 (bound) | 3.54 | 7.9 | 0.77 | 2.91 | 22.9 |
| free offset, started at 0 | 0.00668 | 106/108 | −1.4 ms | 4.03 | 19.7 | 3.65 | | | | |
| free offset, started at −30 ms | 0.00630 | 104/117 | −2.8 ms | 4.15 | 20.0 (bound) | 3.84 | | | | |
| **free offset, started at −80 ms = fit C** | **0.00628** | 112/108 | **−2.7 ms** | **4.15** | **20.0 (bound)** | **3.82** | 7.9 | 0.76 | 2.90 | 22.9 |

**All starts land in one basin.**
- With a slow filament available, the low-activation regime is no longer competitive: it
  returns to 76 % gate opening.
- Freeing the onset moves it only 2.7 ms before the recorded stimulus, well inside the new
  bounds. The 13 % SSR gain comes mostly from `K_coop` reaching its bound.
- Everything downstream changes by less than 2 % (σ, activation, ATP, AUC × σ). The previous
  fit C's numbers (§A.3) remain valid to that precision.

## Files

| file | role |
|---|---|
| [`extractEHTPlate6.py`](extractEHTPlate6.py) | xlsx → `data/travis/processed/plate6_stim_{waveforms,times,wellmetrics}.csv` |
| [`PrepareEHTControls.m`](PrepareEHTControls.m) | beat cutting, control selection, template → `results/eht_controls.mat`, `.png` |
| [`simEHTTwitch.m`](simEHTTwitch.m) | one paced twitch for given Ca and thin-filament parameters; returns force, Ca, gate, N, head pools, force parts and net ATP turnover on the data grid |
| [`FitEHTTwitch.m`](FitEHTTwitch.m) | F0 degeneracy check; fits A and B (cached in `results/eht_fit_F1.mat`; `REFIT = true` redoes them, ~1 h); diagnostics and figure → `results/eht_fit.mat`, `.png` |
| [`ProfileEHTTwitch.m`](ProfileEHTTwitch.m) | §A.2: `k_on` profile and `K_coop`-fixed fits, 3 starts each (cached `results/eht_profile_fits.mat`); Ca-amplitude predictions (P3) → `results/eht_profile.mat`, `.png`, `eht_profile_caamp.mat` |
| [`TestKoffEHTTwitch.m`](TestKoffEHTTwitch.m) | §A.3: `k_on`/`k_off` co-scaling scan; refits with `k_off` free, K2 (cached `eht_koff_K2.mat`) and K3 = **fit C** (cached `eht_koff_bound.mat`) → `results/eht_koff.mat`, `.png` |
| [`CompareCaTransientsEHT.m`](CompareCaTransientsEHT.m) | §2.2 parameter audit (→ `results/eht_param_audit.txt`) and §3.6 Ca-transient comparison (→ `results/eht_ca_transients.png`, `.mat`) |
| [`FitMyoSimEHTTwitch.m`](FitMyoSimEHTTwitch.m) | §8: equivalence check, refit of fit C with `CaModel = 'myosim'` (cached `results/eht_myosim_fits.mat`) → `results/eht_myosim.mat`, `.png` |
| [`NonUniqueEHTTwitch.m`](NonUniqueEHTTwitch.m) | §4: exact symmetries and refits under alternative assumptions (cached `results/eht_nonunique_fits.mat`) → `results/eht_nonunique.mat`, `.png` |
| [`FitCEHTTwitch.m`](FitCEHTTwitch.m) | §3, §4.3, §5–6, §A.4: fit C close-up and the free-offset refit; Jacobian; force–pCa and saturating ATP rate; sensitivities; basin search (cached `eht_fitC_basins.mat`) → `results/eht_fitC.mat`, `results/eht_fitC.png` |

`eht_fit_F1.mat` was seeded from the optima printed in `results/eht_fit_log_F1.txt` (4
significant digits), because the first diagnostics pass crashed after the fits had finished.
Residuals and Jacobians were recomputed at those points; the SSR agrees with the optimizer's to
3 digits. `eht_koff_K2.mat` was seeded the same way from the K2 run's saved results.
