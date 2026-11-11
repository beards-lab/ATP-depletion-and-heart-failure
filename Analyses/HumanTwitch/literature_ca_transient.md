# Cytosolic Ca2+ transient for a 1 Hz hiPSC-EHT twitch — literature survey

Purpose: choose a defensible Ca2+-transient SHAPE and absolute LEVELS to drive the cross-bridge model for a
1 Hz, 37 °C, field-stimulated hiPSC-derived EHT (Mantarray / Curi Bio, 24-well, compliant posts, ~1.8 mM Ca).
Measured twitch to match: stim→peak ~230 ms, 10 %→peak 176 ms, peak→50 % 152 ms, peak→90 % 243 ms.

Provenance: compiled 2026-10-08 from primary sources fetched via Europe PMC / PMC / publisher pages. Every number is
tagged **[M]** measured, **[S]** simulated (model), or **[U]** recalled-but-unverified (full text not retrievable
in this session). SEM/SD as given by the source. Where a source was not reachable this is stated explicitly.

---

## 1. hiPSC-CM Ca2+ transients (single cell, monolayer, 3D EHT)

| Prep / source | Indicator, T, pacing | Maturation | Diastolic | Amplitude / peak | TTP | Decay | Notes |
|---|---|---|---|---|---|---|---|
| **3D hiPSC-ECT**, fibrin, hiPSC-CM+fibroblasts — de Lange et al. 2021, AJP Heart 320:H1670 (PMC8260387), doi:10.1152/ajpheart.00941.2020 | **Fura-2 AM 5 µM**, Krebs-Henseleit, **37 °C**, **1 Hz** (0.5–2.5 Hz tested) | "early" 14–26 d vs "late" 42–57 d in 3D | ratio only; not calibrated | amplitude ~unchanged with culture age, marginal ↑ with frequency [M] | **83 ± 5 ms** [M] | **peak→50 % = 185 ± 15 ms** [M]; stim→75 % decay shortens with age | Force at 1 Hz: **TTP 229 ± 6, RT50 160 ± 3, RT90 262 ± 4 ms** [M] — essentially identical to our Mantarray twitch (230/152/243). **Best anchor.** t-tubules reported. |
| **hiPSC-EHT**, Hamburg (Saleem, Mannhardt, Braren, Denning, Eschenhagen, Hansen 2020, Stem Cell Reports 14:312, PMC7013237, doi:10.1016/j.stemcr.2019.12.011) | **GCaMP6f** (lentiviral), modified Tyrode 1.8 mM Ca, 37 °C, paced (0.75–2.5 Hz; baseline 1.5–2 Hz) | EHT day 20–70 | n/a (GECI) | AU only; CaT +97 % from 0.75→2 Hz (force +96 %) [M] | **TTP−80 % = 116 ± 4 ms** vs force 145 ± 5 ms [M] | **DT80 (peak→20 %) = 173 ± 7 ms** vs force RT80 149 ± 3 ms [M]; TTP/RT fall 10–30 % from 0.75→2.5 Hz | GCaMP6f τoff ~ 100s of ms → decay overestimated; rise OK-ish. "CaT precedes twitch". Isoprenaline 10 nM: TTP −17 %, DT −14 %. |
| **Cardiopatch** (Shadrin et al. 2017, Nat Commun 8:1825, PMC5705709, doi:10.1038/s41467-017-01946-x) | gCaMP6, 37 °C, 1–2.5 Hz | 3–5 weeks | n/a | dF/F stable [M] | not tabulated in retrievable text | — | force rise time 98 ± 2 ms (0.5MM patch, 3 wk) [M]; 22 mN/mm². CaT kinetics not extractable here. |
| **Intensity-trained tissue** (Ronaldson-Bouchard et al. 2018, Nature 556:239, PMC6039110, doi:10.1038/s41586-018-0016-3) | Fluo-4 NW + 5 µM blebbistatin, 37 °C, 150 fps (from patent US11299714B2) | 4 wk; 2→6 Hz ramp | n/a | n/a | — | "decay time (90 %)" and FWHM reported in figures; numbers not retrievable (PMC/Nature blocked) [U] | Qualitatively: faster CaT, positive FFR, t-tubules. Use only qualitatively. |
| **Monolayer hiPSC-CM, 3 labs** (Hwang et al. 2015, JMCC 85:79, PMC4530041, doi:10.1016/j.yjmcc.2015.05.003) | **Fura-2 AM, room temperature, 0.5 Hz** | d15/21/30 post-induction | "not different from rabbit" [M, qualitative] | amplitude & SR content comparable to adult rabbit [M] | **slower than rabbit/mouse** [M] | **slower than rabbit/mouse** [M]; SERCA 61 % / NCX 31 % of removal (rabbit 60/32, mouse 89/9) | Table with nM values not retrievable (PMC captcha, Europe PMC 500). Abstract: cytosolic buffering = rabbit; cell volume 3960 µm³. |
| **Single hiPSC-CM, micropatterned** (Pioner et al. 2022, Front Physiol 13:1030920, PMC9676373) | Cal-630, 37 °C, 1 Hz & 2 Hz, 1.8 mM Ca | d60 → d90 | — | amplitude ↑ d60→d90 [M] | TTP ↓ with age [M] | RT50 ↓ with age [M] | ms values only in figures (not extracted). |
| **Spontaneous single hiPSC-CM, Fluo-4** (experimental dataset in Paci et al. 2018, Biophys J 114:2503, PMC6028769) | Fluo-4 AM 4 µM, 35–37 °C, **spontaneous 0.70 ± 0.38 Hz** | immature | n/a | n/a (ΔF/F) | **TPEAK 270 ± 108 ms**; RT10-50 83 ± 51; RT10-90 167 ± 70 ms [M] | **DT90-10 410 ± 100 ms**; CaT duration 805 ± 188 ms [M] | Spontaneous, unpaced, single/clustered cells → much slower rise than paced 3D tissue. |
| **Optical syncytia, paced** (Paci et al. 2020, Biophys J 118:2596, PMC7231889) | OptoDyCE, **21 °C**, spont. & 0.5 Hz | immature | n/a | n/a | tRise0→peak 97–474 ms (paced) [M] | CTD50 511–1168; CTD90 863–1803 ms (paced) [M] | Room temperature: ~2–3× slower than 37 °C; do not use for levels/kinetics. |
| Rat/mouse/human EHT, Fura-2 (Stoehr et al. 2014, AJP Heart 306:H1353, PMC4116534) | Fura-2, perfused, paced | — | ratio | — | — | — | Human EHT slowest kinetics, lowest Ca sensitivity; **Fura-2 itself buffers Ca and shifts force–[Ca]o** [M]. Numbers not retrievable. |
| Mantarray / Pulse3D / Nautilus | — | — | — | — | — | — | **No peer-reviewed Ca-transient data on Mantarray tissues found** (only press releases on Nautilus optical mapping). |

Indicator caveats: Fura-2 ratio is near-linear up to ~1 µM (Kd ≈ 225 nM at 37 °C) but buffers Ca and slows decay at high loading;
Fluo-4 (Kd ≈ 345 nM) saturates/compresses peaks and ΔF/F is not convertible to µM without calibration; GCaMP6f adds
~50–100 ms rise and ~200–400 ms decay filtering (so Saleem's DT80 173 ms is an upper bound). Room-temperature data
(Hwang, Paci 2020) are ~2–3× slower than 37 °C. Pacing 0.5→2 Hz shortens TTP/decay by ~10–30 % (Saleem) and raises amplitude (positive CaT–frequency relation in EHT).

## 2. Adult human ventricular myocytes / trabeculae

| Source | Prep, indicator, T, rate | Diastolic | Peak / amplitude | TTP | Decay |
|---|---|---|---|---|---|
| Beuckelmann, Näbauer, Erdmann 1992, Circulation 85:1046, doi:10.1161/01.cir.85.3.1046 | isolated LV myocytes, **Fura-2** (pipette or AM), voltage clamp (5 mM Na pipette); temperature & rate not in abstract [U: likely RT] | **95 ± 47 nM** nonfailing; 165 ± 61 nM failing [M] | **peak 746 ± 249 nM** nonfailing; 367 ± 109 nM failing [M] | — | decline slower in failing [M] |
| Coppini et al. 2013, Circulation 127:575, doi:10.1161/CIRCULATIONAHA.112.134932 (26 HCM vs non-failing controls) | septal myocytes, Fluo-4, 0.2/0.5/1 Hz, 35–37 °C | HCM > control [M] | **control amplitude 350 (330–370) nM** as digitised by Tomek 2019 [M via S] | HCM longer [M]; control values only in figures | HCM prolonged [M] |
| Piacentino et al. 2003, Circ Res 92:651 (Fluo-3, 37 °C, 0.5 Hz) | nonfailing vs failing myocytes | — | failing smaller [M] | — | failing slower; reduced SR uptake [M] |
| Vahl et al. 1994, Circ Res 74:952 (trabeculae, 60 bpm) | working human myocardium | — | — | TTP delayed with lower afterload [M] | — |
| Mashali et al. 2021, JMCC 156:7, doi:10.1016/j.yjmcc.2021.03.007 | human trabeculae, near-physiological | force/kinetics only, **no Ca** | — | — | relaxation slower in failing [M] |
| ToR-ORd (Tomek et al. 2019, eLife 8:e48890) calibrated to Coppini | 1 Hz, endo | [U: ~0.08–0.1 µM, recalled] | **312 nM amplitude** baseline; 360 nM when APD → 351 ms [S] | **54 ms** [S], "within SD of Coppini" | **CaTD90 (TTP + peak→90 %) = 459 ms** [S], within SD of Coppini; ORd overestimates duration |

## 3. hiPSC-CM vs adult — mechanisms and magnitudes

| Feature | hiPSC-CM / EHT | Adult human | Source |
|---|---|---|---|
| Rise (TTP) | 80–120 ms in paced 3D tissue at 37 °C (de Lange; Saleem GCaMP) ; 270 ms spontaneous single cells (Paci 2018) | ~50–70 ms (ToR-ORd 54 ms within Coppini SD) | rise slowed by absent/sparse t-tubules → centripetal propagation (~50–100 ms centre–periphery delay, Koivumäki 2018) and larger reliance on ICaL/NCX entry vs RyR release |
| Decay | peak→50 % 185 ms (Fura-2 EHT, 37 °C); DT90-10 410 ms single cell | CaTD90 ~460 ms incl. TTP → peak→90 % ≈ 400 ms (ToR-ORd/Coppini) | SERCA share of removal equal to rabbit (61 %) but SERCA flux slower (Hwang 2015); lower SERCA2a/RyR2/triadin/junctin expression (Lee 2011) |
| Diastolic Ca | "not different from rabbit" (Hwang, RT); tends higher in immature/spontaneous cells | ~95 nM (Beuckelmann) | hiPSC lines with spontaneous activity have incomplete diastole |
| Amplitude | comparable to rabbit (Hwang); Paci2020 model 160 nM (low) | 350 nM (Coppini), 650 nM (Beuckelmann, clamp) | wide, indicator-dependent |
| FFR | positive CaT & force FFR 0.75–2 Hz in EHT (+~100 %) | positive | Saleem 2020 |
| Myofilament | ssTnI ≈ 100 %, β-MHC predominant, MLC-2a 81 % (see §5) | cTnI, β-MHC, MLC-2v | Iorga 2018, Pioner 2016, Bedada 2014 |

## 4. Modelled hiPSC / adult Ca transients ("model-consensus" numbers)

| Model | Rate | Diastolic | Peak / amplitude | TTP (10 %→peak) | Decay | Source |
|---|---|---|---|---|---|---|
| Paci 2018 optimised hiPSC-CM | spont. 0.65 Hz | not stated | not stated | 166 ms [S] (exp 270 ± 108) | DT90-10 368 ms [S] (exp 410 ± 100); duration 618 ms (exp 805 ± 188) | PMC6028769 |
| Paci 2020 hiPSC-CM | spont. 1712 ms CL | not stated | **amplitude 160 nM** [S] | 184 ms [S] | DT90-10 341 ms; duration 691 ms [S] | PMC7231889 |
| Forouzandehmehr 2021 hiPSC-CM-CE (electro-mechanics) | spont. 1644 ms CL; 1 Hz paced for mechanics | not stated | — | 163 ms; RT10-50 46 ms; RT10-90 102 ms [S] | DT90-10 343 ms; duration 693 ms [S]; shortening RT50 161 ms (exp 146–170) | PMC8617339 |
| Kernik 2019 hiPSC-CM | spont. 62 bpm (exp 30 ± 13) | fitted to Wu-lab Fluo-4 data at 37 °C; absolute nM not retrievable [U: ~0.1 µM diastolic, ~1 µM peak per fetch summary — treat as unverified] | — | τdecay/TTP ratio 0.77 [S] (exp 0.83) | — | PMC6767694 |
| Koivumäki 2018 hiPSC-CM (spatial) | spont. | qualitative | qualitative | "substantially slower" than adult; 50–100 ms centre–periphery delay [S] | slower | PMC5808345 |
| Akwaboah 2021 GA hiPSC-CM | spont. 2.26 Hz | — | — | 66 ms; RT10-90 31 ms [S] | DT90-10 127 ms [S] | fast, atypical; doi:10.3389/fphys.2021.675867 |
| **ToR-ORd adult endo** | 1 Hz | [U ~0.08–0.1 µM] | **312 nM amp** (360 nM at APD 351 ms) [S] | **54 ms** [S] | **CaTD90 459 ms** [S] | eLife 8:e48890 |
| ORd 2011 adult | 1 Hz | [U] | [U] | — | "slightly overestimates CaT duration" vs Coppini [S] | per Tomek 2019 |

Take-away: hiPSC single-cell models are tuned to spontaneous, Fluo-4 single-cell data (TTP 160–270 ms, DT90-10 340–410 ms) and
under-predict the speed of paced 3D tissue at 37 °C (TTP ~85 ms). Adult human models agree on amplitude 0.3–0.4 µM, TTP ~50–60 ms, CaTD90 ~450 ms.

## 5. Thin-filament Ca2+ sensitivity in hiPSC-CMs

| Quantity | hiPSC/hESC-CM | Adult human | Source |
|---|---|---|---|
| pCa50 (myofibrils, 15 °C) | **6.04 ± 0.08** (hESC-CM, d35–56) | **5.80 ± 0.05** (adult LV) | Iorga et al. 2018, Front Physiol 8:1111, doi:10.3389/fphys.2017.01111 [M] |
| Hill nH | "similar" to adult (value not extracted) | — | Iorga 2018 |
| Fmax | 42 ± 10 kPa | 94 ± 25 kPa | Iorga 2018 [M]; Pioner 2016: 18.6 ± 2.5 mN/mm² hiPSC-CM myofibrils d80–100 [M] |
| kACT / kTR (pCa 4.0, 15 °C) | 0.66 ± 0.14 / 0.74 ± 0.10 s⁻¹ (Iorga); kACT 0.51 ± 0.04 s⁻¹ (Pioner) | 0.68 ± 0.14 / 0.68 ± 0.10 s⁻¹ | both labs [M] |
| Relaxation (15 °C) | slow phase 0.30 ± 0.04 s⁻¹, fast 2.66 ± 0.3 s⁻¹ (Pioner); kREL ~4–5 s⁻¹ (Iorga) | kREL ~4–5 s⁻¹ | [M] |
| TnI isoform | **ssTnI ~100 %, cTnI below detection** (hESC-CM); hiPSC-CM keep fetal TNNI1 signature even in long-term culture; cTnI:ssTnI ratio proposed as maturation index | cTnI only | Iorga 2018; Bedada et al. 2014, Stem Cell Reports 3:594 [M]; Pioner 2016 "both isoforms" |
| MHC | β-MHC predominant (hiPSC & hESC) | β-MHC | Pioner 2016; Iorga 2018 [M] |
| Other | MLC-2a 81.5 %/MLC-2v 18.2 %; MLC-1a 86 %; cTnT3 + ssTnT; fetal titin | MLC-2v, MLC-1v, cTnT3 | Iorga 2018 [M] |
| Tn Ca2+ off-rate, cTnI vs ssTnI | ssTnI → **higher Ca affinity and slower Ca dissociation** from cTnC "in all systems" (Tn, thin filament, TF+S1); fold values in paper tables (not retrievable) | koff cTnC: Tn complex ~35 s⁻¹; thin filament ~105 s⁻¹; TF + S1 ~13 s⁻¹ (15 °C); myofibrils ~39 s⁻¹ (20 °C) | Davis et al. 2007, Biophys J 92:3195, doi:10.1529/biophysj.106.095406 [M]; Davis & Tikunova 2008, Cardiovasc Res 77:619, doi:10.1093/cvr/cvm098 [M] |
| ssTnI functional signature | blunted pH-desensitisation; ↑pCa50 ~0.2–0.3 units; slowed relaxation (transgenic mice) | — | Davis & Tikunova 2008 review; Westfall/Metzger 1997 PNAS |
| Reference tissue values | slow skeletal (rabbit soleus) pCa50 6.00 ± 0.03, nH 1.61 ± 0.07; porcine cardiac 5.66 ± 0.01, nH 2.74 ± 0.19 | — | Veltri et al. 2017, Front Physiol 8:221 [M] |

Implication: hiPSC-EHT thin filament is ~0.2–0.25 pCa units more Ca-sensitive than adult (ssTnI), with slower Tn Ca off-rate;
if the cross-bridge model's activation uses an adult pCa50 (~5.8), either shift it (~6.0) or scale the prescribed Ca down accordingly.

---

## (b) RECOMMENDED hiPSC-EHT Ca2+ transient, 1 Hz, 37 °C, 1.8 mM Ca (field-stimulated, compliant posts)

| Parameter | Recommended | Plausible range | Confidence / basis |
|---|---|---|---|
| Diastolic [Ca]i | **0.15 µM** | 0.10–0.25 µM | LOW–MED. No calibrated EHT value found. hiPSC monolayer "= rabbit" (Hwang, RT); adult human 95 nM (Beuckelmann); immature cells tend higher; Kernik ~0.1 µM [U]. |
| Peak [Ca]i | **0.6 µM** (amplitude 0.45 µM) | 0.4–1.0 µM | LOW–MED. Adult amplitude 0.31–0.35 µM (Coppini/ToR-ORd); Beuckelmann peak 0.75 µM; Hwang amplitude "= rabbit"; Paci2020 160 nM is an outlier (spontaneous). Fura-2 EHT amplitudes not calibrated anywhere. |
| Stimulus→CaT onset | ~5–10 ms | 2–20 ms | MED (AP upstroke + release; field stimulation). |
| Time to peak | **85 ms** | 70–120 ms | **HIGH-MED**: Fura-2 EHT 83 ± 5 ms (de Lange, same force kinetics as ours); GCaMP 80 % at 116 ms (Saleem) consistent. |
| 10–90 % rise | ~45 ms | 35–70 ms | MED (derived from waveform fits). |
| Peak→50 % decay | **185 ms** | 140–230 ms | **HIGH-MED**: de Lange 185 ± 15 ms (Fura-2, 37 °C, 1 Hz). |
| Peak→80 % decay | ~340 ms | 250–400 ms | MED: Saleem GCaMP 173 ms is for 80 % recovery at 1.5–2 Hz with GECI lag — treat as lower bound. |
| Peak→90 % decay | **~430 ms** | 350–550 ms | MED-LOW: single-cell DT90-10 410 ± 100 ms (Paci exp); adult peak→90 % ~400 ms. |
| Mono-exp τ (if forced) | ~250 ms | 200–330 ms | LOW: decay is sigmoidal (DT50/DT90 ≈ 0.43 vs 0.30 for exponential); a single τ cannot fit both. |
| CaTD90 total | ~520 ms | 430–650 ms | MED. |

Suggested analytic form (checked numerically, normalised shape then scaled):
`Ca(t) = Ca_dia + A · f(t)/max f`, `f(t) = (1 − exp(−t/τr)) / (1 + (t/t50)^h)`, t from stimulus.
hiPSC-EHT: **τr = 25 ms, t50 = 265 ms, h = 3.0** → TTP 88 ms, 10–90 % rise 44 ms, peak→50 % 188 ms, peak→80 % 344 ms, peak→90 % 477 ms, residual 2 % at 1000 ms.
(Alternative τr = 22, t50 = 270, h = 3.2 → TTP 86 / 192 / 338 / 459 ms.) A bi-exponential `exp(−t/τd) − exp(−t/τr)` with τr 35, τd 220 ms gives TTP 76, DT50 190 but DT90 545 ms — too long a tail; gamma-variate forms decay far too fast (DT50 < 100 ms). Avoid both.
Cross-check with the twitch: CaT peak (85 ms) precedes force peak (230 ms) by ~145 ms and CaT 50 %-decay (~270 ms after stim) lands just after force peak — consistent with de Lange's "CaT precedes force" and Saleem's loops.

## (c) Adult human LV, 1 Hz, 37 °C (for contrast)

| Parameter | Recommended | Range | Confidence / basis |
|---|---|---|---|
| Diastolic | **0.10 µM** | 0.08–0.15 | MED: Beuckelmann 95 ± 47 nM (Fura-2); ToR-ORd/ORd ~0.08–0.1 [U]. |
| Peak | **0.45 µM** (amp 0.35) | 0.35–0.75 | MED: Coppini amp 350 (330–370) nM; ToR-ORd 312–360 nM; Beuckelmann peak 746 nM under clamp. |
| TTP | **55–65 ms** | 45–90 | MED: ToR-ORd 54 ms within Coppini SD. |
| Peak→50 % | ~150 ms | 110–200 | MED-LOW (derived; Coppini ms values not extractable). |
| Peak→90 % | ~400 ms | 320–480 | MED: ToR-ORd CaTD90 459 ms incl. TTP, within Coppini SD. |
| Form | τr = 18 ms, t50 = 210 ms, h = 2.8 → TTP 65, DT50 155, DT80 290, DT90 407 ms | | |

Contrast: hiPSC-EHT at 37 °C/1 Hz is ~1.4× slower to peak and ~1.2× slower in early decay than adult, but late decay (peak→90 %) is similar;
diastolic and amplitude are probably ~1.3–1.5× adult but are the least constrained quantities. ssTnI raises thin-filament pCa50 by ~0.2 units.

## (d) Key primary sources

1. de Lange WJ, Farrell ET, Kreitzer CR, Jacobs DR, Lang D, Glukhov AV, Ralphe JC. Human iPSC-engineered cardiac tissue platform faithfully models important cardiac physiology. AJP Heart 2021;320:H1670. doi:10.1152/ajpheart.00941.2020 (PMC8260387).
2. Saleem U, Mannhardt I, Braren I, Denning C, Eschenhagen T, Hansen A. Force and calcium transients analysis in human EHT reveals positive force-frequency relation at physiological frequency. Stem Cell Reports 2020;14:312. doi:10.1016/j.stemcr.2019.12.011 (PMC7013237).
3. Hwang HS et al. (Knollmann lab). Comparable calcium handling of human iPSC-derived cardiomyocytes generated by multiple laboratories. JMCC 2015;85:79. doi:10.1016/j.yjmcc.2015.05.003 (PMC4530041).
4. Beuckelmann DJ, Näbauer M, Erdmann E. Intracellular calcium handling in isolated ventricular myocytes from patients with terminal heart failure. Circulation 1992;85:1046. doi:10.1161/01.cir.85.3.1046.
5. Coppini R et al. Late sodium current inhibition reverses electromechanical dysfunction in human HCM. Circulation 2013;127:575. doi:10.1161/CIRCULATIONAHA.112.134932.
6. Tomek J et al. Development, calibration, and validation of a novel human ventricular myocyte model (ToR-ORd). eLife 2019;8:e48890. doi:10.7554/eLife.48890.
7. Paci M et al. Automatic optimization of an in silico model of hiPSC-CMs recapitulating calcium handling abnormalities. Biophys J 2018;114:2503 (PMC6028769); Paci M et al. All-optical electrophysiology refines populations of in silico hiPSC-CMs. Biophys J 2020;118:2596. doi:10.1016/j.bpj.2020.03.018.
8. Iorga B et al. Differences in contractile function of myofibrils within hESC-CMs vs adult ventricular myofibrils are related to distinct sarcomeric protein isoforms. Front Physiol 2018;8:1111. doi:10.3389/fphys.2017.01111; Pioner JM et al. Stem Cell Reports 2016;6:885. doi:10.1016/j.stemcr.2016.04.006; Bedada FB et al. Stem Cell Reports 2014;3:594.
9. Davis JP et al. Effects of thin and thick filament proteins on calcium binding and exchange with cardiac troponin C. Biophys J 2007;92:3195. doi:10.1529/biophysj.106.095406; Davis JP, Tikunova SB. Cardiovasc Res 2008;77:619. doi:10.1093/cvr/cvm098.

Not retrievable in this session (flagged [U] above): Hwang 2015 Table 1 nM values; Kernik 2019 CaT biomarker table; Ronaldson-Bouchard 2018 CaT kinetics; Coppini 2013 control TTP/decay in ms; Stoehr 2014 human-EHT Fura-2 numbers; ToR-ORd diastolic Ca; Beuckelmann recording temperature. These should be checked against the PDFs before any number is quoted in a manuscript.
