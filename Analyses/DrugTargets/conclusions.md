# Drug targets tested on the human-EHT twitch model: which mechanism each ARPA compound would have to act through, and whether the twitch-force screen can tell

**Status:** 2026-10-12, first pass (single-lever dose sweeps on two cores, 12 multi-lever combinations). Worktree
`C:\home\git\ATP-wt-human` (branch `HumanTwitch`). **Premise (decided 2026-10-12):** no further data will come: no Ca
transient, no waveforms under drug, no kinetics. The only drug observable is the Prism **% change of twitch force
from each tissue's own baseline, minus the matched DMSO control** (`data/travis/drug-vs-control-report/`). The
workflow is therefore inverted: a drug's hypothesised target is written as a model lever, the lever's dose curve is
simulated on the fitted human-EHT twitch, and the observed force change is read back as the lever magnitude the
hypothesis needs. Plausibility is then judged against pharmacology (recalled literature, marked *lit.*; re-check
before quoting) and against what the model says the twitch *shape* would do (shape is in the Mantarray export of
the control plate, not in the drug files, so it is a prediction, not a test).

**Scripts:** `drugBase.m` (base = a HumanTwitch tier optimum + the fixed hiPSC Ca transient + both ATP gates on),
`drugLevers.m` (13 levers with dose grids), `RunDrugLevers.m` (resumable sweep → `results/levers_<core>.{mat,csv}`),
`RunDrugCombos.m` (→ `results/combos.csv`), `PlotDrugLevers.m` (→ `results/levers_doseresponse.png`).
One twitch = 2 s diastole + 2 beats ≈ 6 s; a full sweep of one core ≈ 2 min on the 5-worker pool.

## 1. The base: which Ca transient and which humanised cross-bridge core

| layer | choice | why | alternative kept as a check |
|---|---|---|---|
| Ca transient | **fixed hiPSC-EHT set** (`caTransient` `'risedecay'`: 0.15 → 0.60 µM, TTP 93 ms, peak→50 % 186 ms, peak→90 % 448 ms; de Lange 2021 Fura-2, 3D hiPSC tissue, 37 °C, 1 Hz) | the only paced-tissue 37 °C measurement whose *force* twitch matches ours to 10 ms (TTP 229/RT50 160/RT90 262 vs our 230/152/243); the EHT report's freely fitted transient already had this decay | levels are irrelevant (k_on·Ca symmetry); shape edges (TTP 70–120 ms) move the fitted core 0.2–4× but, as §3 shows, hardly move the *drug* responses |
| cross-bridge core | **`t1_isoT`** (`params/human_twitch_t1_isoT.m`): the mouse incumbent `robust1_seed` with ADP release, hydrolysis, attachment and weak detachment at **×0.45** (isoform 0.1 × Q10 factor 4.5) and the Pi step at ×2.1; thin filament k_on 9.2 µM⁻¹s⁻¹, k_off 5.1 s⁻¹, K_coop 20 | the only fit that stays inside the humanising corridor (0.77 % rms); predicts ktr 23 s⁻¹ at 37 °C, F0 unchanged (54 kPa developed), 1.8 ATP/head/beat, 84 % activation, 33 % resting SRX | **`s3_free`** (0.30 % rms, noise floor): attachment ×0.35, series element ×0.5, F0 halved, 0.4 ATP/head/beat, 41 % resting SRX. Same twitch shape, different muscle; used here to mark which drug conclusions depend on the core |

Both bases reproduce the Plate 6 control twitch (429 µN with their own σ; TTP 240/230, RT50 140/151, RT90 247/246 ms
against 230/152/243 data). Both ATP gates (`UseAtpDetachR2`, K_T2 2.5 mM on the ATP-binding detachment;
`UseAtpKmsrd`, K_srx 9.3 mM on the SRX-ADP return) are switched on with the snapshot values; they are normalised to
1 at 8 mM, so the fitted twitch is untouched and only the MgATP lever sees them. No fit in this project has
constrained these two K values on human tissue; the MgATP rows below inherit that.

## 2. The compounds, their targets, and the model lever each becomes

The panel (Prism files, March–April 2026; Δ = drug − matched DMSO control, % points of baseline twitch force):

| compound | class / clinical signal | myocyte-level target (*lit.*) | model lever(s) | observed Δ force |
|---|---|---|---|---|
| **Mavacamten** | cardiac myosin inhibitor (HCM) | binds the myosin motor, slows Pi release (actin-activated ATPase IC50 ~0.3–0.5 µM), stabilises the autoinhibited IHM/SRX state; no effect on ADP release | `k1` ↓ (Pi release); `SRXstab` (ksr0, ksrd ↑, kmsrd ↓); `ka` ↓ (fewer heads available to attach: the SRX effect seen from the actin side) | −33 / −71 / −91 at 0.2 / 0.3 / 1 µM (single time) |
| **Isoproterenol** | β-agonist | PKA: phospholamban (faster SERCA, larger SR load → larger, shorter Ca transient), L-type Ca current ↑, MyBP-C (heads released from SRX), **cTnI Ser23/24 (faster Ca off-rate: absent in hiPSC-CM, which express ssTnI without the PKA site)**, titin | `CaAmp` ↑, `CaDecay` ↓, `SRXrel` ↑, (`tf_koff` ↑ only if cTnI) | +17 / +16 / +15 at 0.1 / 0.5 / 1 µM (flat: saturated) |
| **Diltiazem** | L-type Ca channel blocker | ICaL block (IC50 low µM in hiPSC-CM) → smaller Ca entry and SR load → smaller transient | `CaAmp` ↓ | +8 / +9 (n.s.) at 0.3 / 0.5 µM; −26 / −31 at 3 / 5 µM |
| **Clozapine** | antipsychotic; myocarditis, cardiomyopathy | acute in vitro: ICaL and hERG block (low–mid µM), muscarinic/adrenergic antagonism, mitochondrial effects; the clinical myocarditis is immune-mediated and not in a dish | `CaAmp` ↓ (ICaL) is the only acute lever; chronic: head loss (σ) | −21 / −52 at 7 / 10 µM (24 h); −31 / −50 (48 h) |
| **Sunitinib** | multi-TKI; LV dysfunction | off-target AMPK and RSK inhibition, mitochondrial dysfunction and ATP loss; reported Ca-handling and ICaL effects in hiPSC-CM | `MgATP` ↓ (both gates); `CaAmp` ↓; `tf_kon` ↓ (desensitisation by Pi/H⁺) | −25 / −52 at 6 / 10 µM, **the same at 24 and 48 h** (fast, non-progressive) |
| **Doxorubicin** | anthracycline; dose-cumulative cardiomyopathy | topoisomerase IIβ, mitochondrial ROS and ATP loss, RyR2/SERCA dysfunction, calpain/titin and sarcomere degradation, apoptosis | `CaAmp` ↓, `CaDecay` ↑ (SERCA), `MgATP` ↓, head loss (σ), `kstiff` ↓ (disarray) | 0 / −30 / −86 at 1 / 3 / 5 µM (48 h; progressive, steep) |
| **Carfilzomib** | proteasome inhibitor; HF signal | proteotoxic/ER stress, **PP2A activation** (dephosphorylates PLB → slower SERCA, MyBP-C → more SRX, AMPK off), mitochondrial dysfunction; eventually death | early: `CaDecay` ↑, `SRXstab`, `MgATP` ↓; late: σ loss | **+26 at 0.1 µM (24 h)**; −65 / −77 / −79 at 0.2 / 0.3 / 0.4 µM (48 h; all-or-none) |
| **5-FU** | antimetabolite; vasospasm/ischaemia | endothelial and vascular (NO, vasospasm), not the myocyte contractile machinery at 24–48 h | none (negative control) | +7…+11 vs a low control; read as 0 |
| **AZT** | NRTI; mitochondrial myopathy | mtDNA polymerase-γ inhibition: mtDNA depletion needs weeks | none at 24–48 h (negative control) | 0 |

The two negative controls are consistent with the model's premise: without a myocyte-level target there is no lever
and the predicted change is 0, which is what was measured.

## 3. Single-lever dose curves on the two cores (`results/levers_*.csv`)

![dose response](results/levers_doseresponse.png)

ΔF = change of peak developed force (%), shape changes in ms against the untreated twitch of the same core.

| lever | t1_isoT ΔF | s3_free ΔF | shape signature (t1_isoT; TTP, RT50) | robust across cores? |
|---|---|---|---|---|
| **CaAmp** ×0.9 / 0.8 / 0.7 / 0.6 / 0.5 | −9 / −24 / −46 / −65 / −79 | −10 / −27 / −48 / −66 / −79 | TTP −10…−40, **RT50 −13…−27** (relaxation faster), RT90 rises again below ×0.6 | **yes** (±3 pp) |
| CaAmp ×1.15 / 1.3 / 1.5 | +6 / +9 / **+11** | +7 / +12 / **+15** | TTP −10, **RT50 +18…+66** (relaxation slower) | yes; **inotropic ceiling** set by the 84 % activation |
| CaDecay ×0.7 / 0.85 / 1.2 / 1.5 / 2 | −10 / −4 / +3 / +6 / +7 | −11 / −4 / +4 / +6 / +7 | TTP and RT50 follow the decay (×0.7: RT50 −34; ×1.5: +71) | yes; **force-weak, shape-strong** |
| CaRise ×0.7 / 1.5 / 2 | 0 / +1 / +1 | 0 / +1 / +1 | TTP ±10 only | yes (inert) |
| **ka** ×0.7 / 0.5 / 0.35 / 0.2 / 0.1 | −25 / −46 / −65 / −82 / −92 | −30 / −52 / −67 / −82 / −91 | **TTP −10…−40, RT50 −13…−42**: a weaker *and faster* twitch | **yes** (±6 pp) |
| **k1** ×0.7 … 0.1 (Pi release) | 0 … −6 | 0 … −2 | none | yes: **the twitch does not see the Pi step** |
| **SRXstab** ×1.5 / 2 / 3 / 5 / 10 | −5 / −14 / −54 / −87 / −97 | −19 / −42 / −72 / −90 / −97 | **TTP +10…+20** at ×1.5–3 (reservoir slows the rise), RT50 −12…−44; resting SRX 0.5 → 0.98 | **no** below ×3 (the s3 core has the larger SRX reserve) |
| SRXrel ×1.5 / 2 / 3 | +1 / +2 / +2 | +7 / +10 / +11 | RT50 +11…+19 | **no** |
| k2 ×0.35 / 0.5 / 0.7 / 1.5 | +37 / +29 / +17 / −23 | +12 / +8 / +5 / −8 | TTP +20…+80, RT50 +13…+36 | **no** (3× weaker on s3) |
| **MgATP** 4 / 2 / 1 / 0.5 / 0.25 mM | **+9 / +19 / +26 / +10 / −53** | **−7 / −24 / −51 / −74 / −87** | t1: TTP +20…+160 (slower everything); s3: RT50 −11…−33 | **opposite sign**: the ATP gate on detachment wins on t1, the SRX-parking gate wins on s3 |
| tf_koff ×0.7 / 1.5 / 2 / 3 | +8 / −69 / −89 / −96 | +13 / −68 / −88 / −95 | TTP −50…−100 ms, gate peak 0.79 → 0.2–0.03 | yes |
| tf_kon ×0.7 / 0.85 / 1.2 / 1.5 | −66 / −25 / +8 / +8 | −66 / −28 / +10 / +14 | ×0.85: RT50 −23; ×1.2: RT50 +56 | yes |
| kd ×2 / 4 | −13 / −32 | −46 / −73 | RT50 −12 | **no** |
| kstiff ×0.7 / 0.5 | −32 / −54 | −34 / −56 | RT50 +7…+14, TTP −10 | yes (≈ per-head force scaling; a pure σ loss is exactly shape-neutral) |

What this says about the base choice:

1. **Every lever on the Ca side and on the attachment limb is core-independent** (Ca amplitude and decay, ka, k1, the
   thin filament, cross-bridge stiffness). For those drug hypotheses the "most plausible combination" question is
   moot: t1_isoT and s3_free give the same answer, and so, by the case-B analysis (HumanTwitch §5.2), does every
   other core that fits the twitch. The thin filament and the activation level, which the twitch identifies, carry
   these responses.
2. **The levers that depend on the core are exactly the energetic ones: SRX stabilisation/release, ADP release and
   MgATP.** On the realistic core, lowering ATP *raises* twitch force (+26 % at 1 mM) until the SRX gate takes over
   below 0.5 mM; on the possible core, which has a larger SRX reserve and a 0.4 ATP/beat economy, force falls from
   the first mM (−24 % at 2 mM). This is the same k2-gate-versus-SRX-gate question the skinned-muscle 2 mM campaign
   is about, now sitting under the three "ATP-depletion" compounds (sunitinib, doxorubicin, carfilzomib). The
   twitch-force screen cannot settle it, because the sign of the force response is the hypothesis itself.
3. **Mavacamten's canonical target, Pi release, is invisible in this model's twitch** (≤ 6 % at k1 ×0.1, on both cores,
   consistent with HumanTwitch §2.1 where the twitch could not identify k1,k_1). A mavacamten model here must act
   through the SRX pool or through attachment availability. That is not necessarily wrong for an intact twitch
   (the SRX/IHM part of mavacamten's action is the one that scales with occupancy in vivo), but it means the
   dose curve is read as "fraction of heads parked", not as "Pi release slowed".

## 4. Reading each compound's force change back as a lever magnitude

Interpolated on the t1_isoT curves (s3_free in brackets where it differs). "Plausible" judges the implied magnitude
against the pharmacology; "signature" is the shape the model attaches to that reading, which is not in the Prism
files but is in a Mantarray export.

| compound (dose) | observed Δ | lever reading | plausible? | model signature |
|---|---|---|---|---|
| Mavacamten 0.2 / 0.3 / 1 µM | −33 / −71 / −91 | `ka` ×0.62 / ×0.30 / ×0.11 ‖ `SRXstab` ×2.5 / ×3.5 / ×6 [×1.8 / ×3 / ×5] | both routes reproduce the whole curve; the 0.2→0.3 µM jump (−33 → −71) is steeper than a Hill-1 binding to the motor (ATPase IC50 0.3–0.5 µM *lit.*) and fits an SRX mechanism with cooperativity, or the tissue's Ca sensitivity sitting on the steep part of its curve (the model's K_coop 20) | **ka**: TTP −20…−40, RT50 −20…−40 ms, resting SRX unchanged. **SRXstab**: TTP *+10…+20* at the two lower doses, RT50 −26…−42, resting SRX 0.33 → 0.8–0.9. The TTP sign separates the two on t1 (on s3 both are negative). ATP/beat falls with force on either route (economy preserved) |
| Isoproterenol 0.1–1 µM | +15…+17 (saturated) | `CaAmp` ×1.5 alone +11 [+15]; + `CaDecay` ×0.8 +9.7 [+12.8]; + `SRXrel` ×2 +9.6 [**+20.8**]; + `tf_koff` ×1.5 (cTnI) **−5.5** [+3.2] | **the observed +15 % sits at the inotropic ceiling of the realistic core**: 84 % activation leaves +11 % headroom for any Ca increase, and faster Ca decay (the PLB effect) takes some back. It is reached only when the core has an SRX reserve to release (s3 +21 %) or activation is lower than fitted. The cTnI line confirms what the hiPSC biology says: with a PKA-sensitive TnI the response would invert, so the positive inotropy is itself evidence of ssTnI | Ca ×1.5 + decay ×0.8: TTP −30 ms, RT50 +25 ms (slower relaxation!); the lusitropic effect Saleem 2020 reports (TTP −17 %, decay −14 %) needs Ca decay ×0.7 or faster, which costs force. The model cannot give +15 % force *and* faster relaxation with Ca alone; SRX release does both on s3 |
| Diltiazem 3 / 5 µM | −26 / −31 | `CaAmp` ×0.79 / ×0.76 [×0.80 / ×0.77] | yes: 60–70 % ICaL block at these concentrations (*lit.*) with SR-load autoregulation damping the transient to ~−25 % | RT50 −15 ms, TTP −10 ms: a smaller, slightly faster twitch |
| Diltiazem 0.3 / 0.5 µM | +8 / +9 (n.s.) | no single Ca lever gives a gain; `CaAmp` ×0.9 + `CaDecay` ×1.3 → −2 | not reproducible by Ca entry alone; most likely noise (n.s.) or a rate/spontaneous-activity effect the paced twitch model has no handle on | — |
| Clozapine 7 / 10 µM | −21 / −52 (24 h) | `CaAmp` ×0.82 / ×0.68 (ICaL) ‖ σ ×0.79 / ×0.48 (head loss) | ICaL block in the 5–10 µM range is documented (*lit.*); the 24 h = 48 h stability favours a functional (Ca) over a cumulative (death) mechanism | Ca route: RT50 −15…−27 ms; σ route: no shape change. These two are separable only by kinetics |
| Sunitinib 6 / 10 µM | −25 / −52 (24 h = 48 h) | `CaAmp` ×0.80 / ×0.66 ‖ `tf_kon` ×0.85 / ×0.76 ‖ `MgATP`: **impossible on t1** (force rises until 0.5 mM), 2 / 1 mM on s3 | the ATP hypothesis is viable only on the SRX-rich core, and 1 mM bulk ATP in a beating, non-dying tissue is a strong claim; Ca-handling or desensitisation readings are modest and the non-progressive time course fits them | Ca route: faster relaxation; ATP route on s3: RT50 −26…−29 ms, resting SRX 0.64–0.77 (a mant-ATP-visible change) |
| Doxorubicin 3 / 5 µM (48 h) | −30 / −86 | `CaAmp` ×0.77 / ×0.45 ‖ σ ×0.70 / ×0.14 ‖ `kstiff` ×0.72 / — ; `CaDecay` ×1.5 (SERCA) gives +6, wrong sign; `MgATP` wrong sign on t1 | the steep, progressive 48 h collapse is what cell death and sarcomere loss look like; the model can only say that neither slower SERCA nor ATP depletion (on t1) can produce a loss, so the −86 % is Ca amplitude and/or heads | σ: shape-neutral; Ca: RT50 −27 then RT90 +47 (a small, long-tailed twitch at ×0.5) |
| Carfilzomib 0.1 µM (24 h) | **+26** | `k2` ×0.6 on t1 (+17…+29) ‖ `MgATP` 1–2 mM on t1 (+19…+26) ‖ `CaDecay` ×2 +7, `CaAmp` ×1.5 +11: too small | **the one positive signal only an energetic lever on the realistic core can produce**: slower ATP-dependent detachment (early mitochondrial/ATP decline, or PP2A-driven dephosphorylation slowing the cycle) raises force transiently before the 48 h collapse. On s3 no lever gives +26 %. Weak evidence (p 0.016, one dose, one time) but it is the only drug observation that prefers t1 over s3 | TTP +50…+100 ms and RT50 +20 ms: a markedly slower twitch at 24 h, then collapse |
| Carfilzomib 0.2–0.4 µM (48 h) | −65…−79, flat | σ ×0.2–0.35 (death) ‖ `SRXstab` ×3–4 [×2.5–3] | all-or-none across a 2× dose range is a death threshold, not a graded mechanism | σ: shape-neutral |

## 5. Where this leaves "the most plausible combination"

- **Base:** the de Lange hiPSC-EHT Ca transient with the corridor core `t1_isoT` is the defensible primary
  (literature-anchored transient, every humanising factor inside its corridor, a 37 °C ktr of 23 s⁻¹). `s3_free`
  stays as the mandatory companion for any lever that touches SRX, ADP release or ATP, because those answers flip
  with the core and the twitch cannot rank the cores.
- **Drug hypotheses the screen can test with either core:** mavacamten (as attachment loss or SRX parking, not as
  Pi release), diltiazem and clozapine as ICaL block, doxorubicin as amplitude/head loss, and the two negatives.
  The force readings are self-consistent with the pharmacology at the tested doses.
- **Drug hypotheses the screen cannot test:** anything energetic. ATP depletion raises force on t1 and lowers it on
  s3; sunitinib's loss and carfilzomib's early gain are each "explained" by ATP on one core and impossible on the
  other. The carfilzomib +26 % is the only observation that points one way (t1), and it is weak.
- **Isoproterenol is the sharpest constraint in the set** despite being a reference drug: +15 % is at the
  realistic core's inotropic ceiling, is reached comfortably only with an SRX reserve to release, and would turn
  negative with a PKA-sensitive TnI. It argues for ssTnI (as expected) and mildly for a core with an SRX reserve
  (s3-like) or for a lower activation than the fits give.

## 6. Next steps (no new data)

1. **Dose–response mapping** for the two drugs with clean single targets: fit a Hill occupancy (IC50, n) between
   concentration and the lever magnitude for mavacamten (ka and SRXstab routes) and diltiazem (CaAmp) and compare
   IC50/n with the literature binding numbers. If the mavacamten IC50 lands at 0.3–0.5 µM with n ≈ 1 on the SRX
   route and needs n > 2 on the ka route (or vice versa), that is a mechanism statement from force alone.
2. **Isoproterenol as a constraint on activation:** profile the fitted activation level (Ca₅₀ relative to peak Ca,
   via tf_k_on at fixed twitch fit) against the +15 % requirement, with and without SRX release; report the
   activation window compatible with both the control twitch and the iso response.
3. **Energetic drugs:** run sunitinib/carfilzomib readings as 2-D maps (MgATP × CaAmp, MgATP × SRXstab) on both cores
   and tabulate which (ATP, Ca) pairs reproduce the observed Δ at 24 h and 48 h, with the implied ATP/beat and SRX
   fraction, so that the eventual choice of core can be stated as a prediction of tissue ATP.
4. **Shape predictions table** for the lab: for each compound and dose, the TTP/RT50/RT90 change each hypothesis
   implies, so that if the Mantarray kinetics of the drug plates are ever exported (they exist in the
   per-well metrics of the control plate), each hypothesis is pre-registered.
5. Add `k_on·Ca` desensitisation (Pi, H⁺) as an explicit ischaemia lever (tf_kon) in the energetic maps; it is the
   one metabolic effect that is core-independent and force-lowering.

## Files

| file | role |
|---|---|
| `drugBase.m` | base model: tier optimum (`t1_isoT` default, `s3_free`, …) + fixed hiPSC Ca transient + ATP gates on |
| `drugLevers.m` | the 13 single-mechanism levers, their dose grids and the compounds they stand in for |
| `RunDrugLevers.m` → `results/levers_<core>.{mat,csv}` | resumable single-lever sweeps (force, shape, ATP/beat, SRX, gate) |
| `RunDrugCombos.m` → `results/combos.{mat,csv}` | multi-lever hypotheses (isoproterenol, low-dose diltiazem, ATP-depletion drugs, mavacamten) |
| `PlotDrugLevers.m` → `results/levers_doseresponse.png` | dose curves of every lever on both cores + the relaxation signature |
