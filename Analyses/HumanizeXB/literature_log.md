# Humanizing the cross-bridge model — literature log

2026-10-06. Companion to [`conclusions.md`](conclusions.md) (the proposed human parametrisation and
what the model does with it). This file is the evidence: what is known about how human ventricular
cross-bridge kinetics differ from the rodent kinetics the model is fitted to, which numbers are
trustworthy, which are superseded or assay artefacts, and which human measurements the model can
be run against.

Verification key used throughout: **V** = read from the paper's own text/tables; **R** = taken from
a review or a citing paper (Walklate, Ujfalusi & Geeves 2021 *CMLS* 78:7309 unless stated); **U** =
abstract-only or search-snippet, unconfirmed. Values in parentheses are ± SEM/SD as reported.

---

## 0. Where the model stands (the thing being humanised)

| item | current state | source |
|---|---|---|
| preparation fitted | skinned RV trabeculae, Baker lab (UCSF); **mouse per the abstract text, rat per its method list** — the species is not recorded consistently in this repo | `Docs/presentations/2026 Myofilament-meeting-abstract.md` lines 15 and 37 |
| isoform | α-MHC (mouse ~100 %, adult rat 80–97 %) | `Model/parameterBounds.m` header |
| temperature | **not recorded** anywhere in the repo; bounds and notes assume 20–25 °C ("room temperature"); some 2024 sessions are tagged `_10C` | `Docs/notes/labdiary_boundedfit.md`, `Analyses/CalciumTwitch/conclusions.md` §3.4 |
| activation | saturating Ca (`Ca = 1000`, `UseCa = 0`); 8 mM MgATP, 0 ADP, 0 Pi | `params/robust1_seed.m` |
| incumbent | `params/robust1_seed.m` = `refit9_opt` + robust reader, cost 0.424 on the pooled 8 mM target | `Analyses/LowATP_Refit9K2Sweep/conclusions.md` |

Incumbent read-out with the probe written for this analysis (`humanProbe.m`, 8 mM, SL 2.2 µm, see
`results/explore1.txt`):

| observable | incumbent | note |
|---|---|---|
| isometric force F0 (active) | 66 kPa (71 with passive) | pooled data 65.5 |
| k_tr (Brenner-style, instant detachment) | 45 s⁻¹ | data 44–46 |
| attached fraction p1+p2 | 0.28 | |
| SRX fraction at full activation / at rest | 0.03 / 0.27 | `SRXDuringContraction` |
| gross ATP turnover per head (XTOR) / net | 19.3 / 14.2 s⁻¹ | objective bound 3–15; both are high for mouse-α (`AttachedPool_vs_ATPase`) |
| effective post-stroke detachment k2eff / mean attached dwell | 80 s⁻¹ / 15 ms | |
| F/F0 at 0.5, 1, 2 ML/s | 0.69, 0.47, 0.21 | (active force, SL 2.2 start) |
| V_max (Hill extrapolation; model force never reaches 0 below 3 ML/s) | ≈ 4 ML/s | not kinetic — see §7 |
| base rates (s⁻¹) | ka 292, kd 93, k1 161, k_1 31, k2 117, kah 176, kamh 8.8; ksr0 367, kmsrd 27, ksr2srd 130, ksrd2sr 9.2 | |

Two things follow immediately. (i) Humanisation is first an **isoform** problem (α → β), second a
**temperature** problem (nearly all human data are at 15 °C, the model sits at ~20–25 °C), and only
third a species problem. (ii) Because the model's own temperature is unknown, every "human at the
model's temperature" number below carries a ×1.5–2 uncertainty from that alone. **Ask the Baker lab
for the rig temperature and the species before any human fit is scored.**

---

## 1. Myosin heavy-chain isoform — the first-order fact

| species / chamber | % α-MHC (protein) | source |
|---|---|---|
| human LV, non-failing | 7 ± 3 (mRNA ~35 %); failing ≈ 0 | Miyata 2000 *Circ Res* 86:386 **V**; Reiser 2001 *AJP* 280:H1814 **V** (1–10 %) |
| human LA, non-failing | ~90 (Reiser); 75 (Narolska 2005); 67 ± 4 (Noguchi 2003) | **V/R** |
| mouse ventricle | > 94–100 | Milani-Nejad & Janssen 2014 *Pharmacol Ther* table **V** |
| rat ventricle, adult | 80–97 (strain, age) | same; Fitzsimons 1998 (~80 %) **V** |
| rabbit / pig / dog | 5–15 / 3–13 / 0–20 | same |
| PTU (hypothyroid) rat | 100 % β | Fitzsimons 1998 **V** |

Consequences that matter for the parametrisation:

- **The α/β difference is kinetic, not mechanical.** Unitary step (7–10 nm) and force per head do not
  differ between rat V1/V3 or rabbit V1/V3 (Sugiura 1998, Palmiter 1999, **V** via `parameterBounds.m`);
  in large mammals β produces ~2× the force per head of α in motility (Malmqvist 2004, VanBuren 1995 **R**),
  but mouse α = mouse β, and at the fibre level specific tension is the same across species (§4).
- **Isoform mixtures are non-linear.** In motility a 1:9 V1:V3 mix moves at the V3 speed (Noguchi 2003 **V**),
  while 12 % α raises rat myocyte power by +52 % (Herron & McDonald 2002 **V**). So the 5–10 % α in
  human ventricle can be ignored for velocity but not necessarily for power. Model as 100 % β.
- **β-rat is not a human surrogate.** Rat β-S1/HMM releases ADP 1.8× (Johnson 2021, 20 °C) to 3.8×
  (Pereira 2001 vs pig β, 15 °C) faster than human/pig β; PTU-rat k_Ca at 15 °C is 3.5 s⁻¹ vs human
  0.7–1.3 (§2). Use PTU-rat data only for the *direction* of α→β effects.

---

## 2. Force-redevelopment rate k_tr (and kACT, kREL)

### 2.1 Measurements

| paper | prep | species (MHC) | T °C | kACT | **k_tr** (s⁻¹) | slow kREL | fast kREL | ver. |
|---|---|---|---|---|---|---|---|---|
| Piroddi 2007 *Pflügers* 454:63 | myofibrils | human V (β) / A (75 % α) | 15 | 0.73 / 3.73 | **0.68 ± 0.03 / 3.55 ± 0.19** | 0.15 / 0.52 | 2.9 / 16 | R |
| Piroddi 2019 *JGP* 151:18 | myofibrils, 5 donors | human V | 15 | 0.85 ± 0.02 | **0.55 ± 0.03** | 0.29 (186 ms) | 4.18 | V |
| Vitale 2019/21 *JMRCM* | myofibrils | human V donor / **mouse** V | 15 | 0.84 / 7.35 | **0.78 ± 0.04 / 7.12 ± 0.61** | 0.30 / 1.45 | 4.65 / 28.3 | V |
| Vikhorev 2017 *Sci Rep* | myofibrils | human V donor | 17 | 0.78 | **0.64 ± 0.03** | 100 ms phase | 17.1 (different definition) | V |
| Stehle 2002 *Basic Res Cardiol* | myofibril bundles | human / guinea pig / mouse | **10** | 0.3 / 1.7 / 5.2 | – | 0.15 / 0.6 / 1.8 | 4.6 / 11 / 18 | V (abstract) |
| Edes 2007 *AJP Regul* 293:R20 | skinned myocytes | human donor / **mouse** / pig | 15 | – | **1.33 ± 0.11 / 7.44 ± 0.15 / 1.02 ± 0.05** | – | – | V |
| Narolska 2005 *Cardiovasc Res* 65:221 | skinned myocytes | human V / A | 15 | – | **0.87 ± 0.06 / 7.76 ± 1.54** | – | – | R |
| van Dijk 2014 *Pflügers* | skinned strips | human donor / HCM | 15 | – | **0.75 ± 0.09 / 0.63 ± 0.03** (stretch-activation k2 0.74 / 0.77) | – | – | V |
| Milani-Nejad 2020 *Front Physiol* 11:468 | myocyte fragments | human donor (7 hearts) | 15 | – | **0.60 ± 0.05** | – | – | V |
| Lewalle 2025 *JMCC* 211:64 | skinned fibre bundles | human LV / LA, 9 donors | **37** | – | **15 ± 3 / 40 ± 13** (from t½) | – | – | V |
| Milani-Nejad 2015 *AJP* | **intact** RV trabeculae, K⁺ contracture | human donor | 27 / 32 / 37 | – | **6.6 / 14.1 / 18.1** (Q10 ≈ 2.8, n = 4) | – | – | V |
| Milani-Nejad 2013 *JGP* | intact trabeculae | **rat** | 37 | – | 27.7 (Lopt), 45–47 (L90) | – | – | V |
| Palmer & Kentish 1998 *Circ Res* | skinned trabeculae | rat / guinea pig | 22 | 14.4 / 2.6 | **13.0 ± 0.6 / 2.7 ± 0.3** | – | – | V |
| Morotti 2025 *IJMS* 26:469 | skinned trabeculae | rat | 10/15/20/25/30 | – | 5.9 / 8.3 / 12.9 / 18.8 / 29.3 (Q10 2.2) | – | – | V |
| de Tombe & Stienen 2007 *J Physiol* 584:591 | skinned trabeculae | rat | 15–25 | – | 14.5 / 31.4 / 49.4 (**Q10 3.5**) | – | – | V (values via `AttachedPool_vs_ATPase`) |
| Fitzsimons, Patel & Moss 1998 | skinned, flash photolysis | rat α / rat β (PTU) | 15 | k_Ca 8.9 / 3.5 | – | t½ 65 / 139 ms | – | V |
| Korte & McDonald 2007 | skinned myocytes | rat α / β | 12 | – | 2.4 / 1.2 | – | – | V |

### 2.2 Reading it critically

- **Human ventricular k_tr at 15 °C, max Ca = 0.55–0.85 s⁻¹ in myofibrils, 0.6–1.3 s⁻¹ in skinned
  cells/strips** (seven labs). Between-study scatter is ±30 %; Piroddi's own kTR drifted 0.68 → 0.55
  across 12 years while kACT went 0.73 → 0.85. **Working value 0.8 s⁻¹ at 15 °C.**
- **Same-lab rodent at 15 °C: mouse 7.1–7.4** (Vitale, Edes). **Human/mouse-α ratio ≈ 0.10–0.12** at
  matched temperature and prep. Rat trabeculae 8.3 (Morotti) or 14.5 (de Tombe) at 15 °C give
  0.06–0.1. Rat myofibrils (3.7–5.6, Rao 2014 **U**) and rat 12 °C myocytes (2.4) are prep-depressed.
- kACT ≈ kTR in human myofibrils → the Ca switch is not rate-limiting in β ventricle; the model's
  Ca-free, saturating-activation framing is fine for human too.
- **Intact 37 °C numbers compress the species ratio and should not set it.** Human intact 18 s⁻¹ vs rat
  intact 28 s⁻¹ (×0.65) — but rat skinned 20 °C is 13–31 s⁻¹ and Q10 3.5 would put 37 °C rat at
  110–260 s⁻¹, far above the 28 measured in intact muscle. K⁺-contracture k_tr in a 0.3–0.5 mm intact
  trabecula at 37 °C is limited by activation/diffusion and by the fit (t½ based), not by cross-bridge
  kinetics. Likewise Lewalle's 37 °C fibre LV/LA ratio (2.7×) is smaller than the 5× myofibril ratio.
  Use these only as *ceilings* for a 37 °C human model.
- **Q10.** No human temperature series below 27 °C. Projection human 15 °C (0.75) → 37 °C with Q10 3.5
  gives 12 s⁻¹ and with 2.8 gives ~8, vs 15–18 measured in intact/fibre at 37 °C: constant-Q10
  extrapolation from 15 °C under-predicts ~2×. Arrhenius curvature plus the cold-destabilised OFF state
  (Park-Holohan 2021, midpoint 16–18 °C) are the likely reasons. **For the model's 20–25 °C: human
  k_tr ≈ 1.5–3 s⁻¹** (Q10 3–3.5 from 15 °C), vs the incumbent's 45. Ratio 0.04–0.07 if the rig is at
  25 °C, 0.1 if at 20 °C. The isoform ratio (0.1) is the number to trust; the extra factor is the
  unknown rig temperature.
- Human HF: no donor-vs-HFrEF skinned k_tr number exists; HCM is 15 % slower (van Dijk). Chung 2019
  found slower intact k_tr in failing trabeculae only at 3 Hz.

---

## 3. Unloaded shortening velocity, force–velocity, power

| paper | prep | species | T °C | value | ver. |
|---|---|---|---|---|---|
| Fitzsimons 1998 | rat myocytes, slack test | rat α / rat β | 15 | V0 **4.72 / 0.92 ML s⁻¹** (unit header ambiguous, context says ML/s) | V |
| Hofmann & Lange 1994 *Circ Res* | rat myocytes, slack test | rat α | (15) | 2.88 ML s⁻¹ | V |
| Pereira 2001 *JBC* 276:4409 | myocytes, slack test | **rat β / pig β** | 15 | **1.1 / 0.35 ML s⁻¹**; ADP release 135 / 35 s⁻¹ | V |
| Korte & McDonald 2007 | rat myocytes, force clamp | rat α / β | 12 | V_opt 0.66 / 0.43 µm s⁻¹ sarc⁻¹; **a/P0 0.10 / 0.05**; peak norm. power 0.14 / 0.075 | V |
| Milani-Nejad 2020 | human myocyte fragments, force clamp | human β | 15 | peak power 54 ± 12 pW, 1.23 µW mg⁻¹ (rat α 12 °C: 1.49); **V_max and a/P0 not tabulated** | V |
| de Tombe & ter Keurs 1990/92 | intact / skinned rat trabeculae | rat α | 20–25 | V0 ~10 µm s⁻¹ hs⁻¹ at 20 °C, ~25 at 25 °C; **Q10 4.6** | V |
| Noguchi 2003 *JMCC* | motility, native myosin | human NF V / failing / A | – | NF = failing; atrium 2× speed, ½ force; 1:9 V1:V3 ≡ V3 | V (qualitative) |
| Johnson 2019/2021 *PLoS Biol* | motility, S1 | human β / rat-chimera | 20 | 0.49 / 0.90 µm s⁻¹; literature human β 0.38 vs rat β 0.62 at 25 °C | V |
| Walklate 2021 (review) | motility, full-length | human α / β | 30 | 2.0 / 1.2–1.5 µm s⁻¹; α:β 2.1–2.4 in every species | R |
| Ishii 2023 *JGP* | motility, porcine β | – | 23 / 32 / 37 / 40 | 0.4 / 1.0 / 1.5 / 2.0 µm s⁻¹ (Q10 3.1–3.5) | V |

Reading it:

- **There is no human ventricular skinned-fibre V0 in the accessible literature** (flag). Best estimate
  at 15 °C: pig β 0.35 ML/s scaled by the human/pig k_tr ratio (1.33/1.02) → **0.4–0.5 ML/s**;
  motility cross-check (human β ≈ 0.6× rat β, rat β V0 ≈ 1 ML/s) → 0.6. Rat α at 15 °C is 2.9–4.7
  ML/s. **Human/rat-α V0 ratio ≈ 0.10–0.15 at the fibre level**, vs 0.2–0.3 from motility
  (α/β 2–2.5 within species × rodent-β/human-β 1.7–3). The fibre number is the one that maps onto
  the model's V_max; the motility number is contaminated by surface/assay effects.
- The chain human β < pig β < rat β < rat α is robust and tracks ADP release (35 / 135 / > 1000 s⁻¹).
- FV curvature: rat α a/P0 0.10 at 12 °C (Korte); rat β 0.05; **no human a/P0 anywhere** — the
  McDonald-lab human myocytes would have it in raw form. Expect human a/P0 ≤ 0.1 (more curved
  than the incumbent's fitted 0.39 — which itself is a 20–25 °C α number from the Baker FV set).
- Temperature: V0 Q10 4.6 in rat trabeculae (de Tombe 1990) — larger than k_tr's 3.5. A human
  model at 37 °C needs V0 ×~30 relative to 15 °C; at the rig's 20–25 °C ×2–5: **human V0 at
  20–25 °C ≈ 1–2 ML/s** vs the incumbent's ~4 (Hill) / > 3 (zero crossing not reached).

---

## 4. Maximal isometric force — not a species discriminator

| prep | human | rodent | source |
|---|---|---|---|
| single myofibrils, 15 °C | 101–111 kPa (Piroddi 2019, Vikhorev); 108 (2007) | mouse/guinea-pig bundles 140–160 at 10 °C (Stehle) | V/R |
| skinned myocytes/strips, 15–20 °C | 35–47 kN m⁻² (Papp, van Dijk, Milani-Nejad, Piroddi strips 38.7) | rat trabeculae 37–53 at 20 °C (Ebus & Stienen) | V |
| multicellular strips, 15 °C | 21 (Narolska) | – | R |
| skinned fibres 37 °C | 16.8 ± 1.1 (Lewalle, ≤ 30 % rundown accepted) | – | V |
| intact, 37 °C | 25–44 mN mm⁻² (Alpert lab); 14–18 twitch / 36 K⁺ contracture (Janssen) | rat 28, mouse 40 (Janssen table) | V |

Verdict: force per cross-section is set by preparation quality, compliance and rundown, not by
isoform. **Keep the model's stiffness, step and F0 as they are**; put the whole human difference
into kinetics. Q10 of Tmax is only ~1.3 (de Tombe & Stienen 2007).

---

## 5. ATPase and tension cost

| paper | prep | T °C | ATPase | force | tension cost (µmol l⁻¹ s⁻¹ per kN m⁻²) | ver. |
|---|---|---|---|---|---|---|
| Narolska 2005 *CVR* | human donor V / A strips | 20 | 0.051 / 0.260 mM s⁻¹ | 21 / 33 | **2.4 ± 0.3 / 11.4 ± 1.4** | R |
| Piroddi 2019 | human donor strips (2 hearts) | 20 | 91 ± 5 µM s⁻¹ (basal 10.9) | 38.7 | **2.3 ± 0.1** | V |
| Witjas-Paalberends 2014 | HCM-smn / R403Q strips | 20 | – | – | 1.78 / 2.93 | V |
| Ebus & Stienen 1994/96/2001 | rat trabeculae | 20 | 0.43–0.49 mM s⁻¹ | 45–53 | **≈ 9–10** | V |
| Rossmanith 1995 | rat V1 / V3 | 24 | – | – | ratio 2.3 | V |
| Nguyen 1996 *Circ Res* | purified human myosin | ? | **0.3 s⁻¹ head⁻¹**, unchanged in HF | – | – | V (abstract) |
| de Tombe & Stienen 2007 | rat | 15–25 | Q10 ATPase 3.5, tension cost 2.5 | | | V |
| Narolska 2005 *JMRCM* | donor vs failing human V | 20 | economy **preserved** in end-stage HF | | | V (abstract) |

Reading it:

- **Human tension cost 2.3–2.4 vs rat 9–10, same lab and method, 20 °C: ×0.25.** Absolute human
  ATPase disagrees 1.8× between the two Amsterdam datasets because force differs (21 vs 39 kN m⁻²);
  tension cost is the robust number.
- Per head (assuming 0.15–0.2 mM heads, not measured in these papers): human **0.45–0.6 s⁻¹** at
  20 °C, max Ca; rat α 2.4–3.2. The model's incumbent is 19 (gross) / 14 (net) s⁻¹ — already
  3–5× too high for mouse-α (`AttachedPool_vs_ATPase`, `project_xtor_is_gross`). A human version
  at 20 °C should sit near **0.5–1 s⁻¹ net**, i.e. ×0.05 of the incumbent; at 25 °C ×0.1.
- The atrial/ventricular 5× (Narolska) is the cleanest in-species α/β energetics number and matches
  the 5× kACT ratio in human myofibrils (Piroddi 2007). Rat V1/V3 is only 2.3× — β-rat again.
- HF: economy and myosin ATPase per head are *not* reduced in human HF (Narolska, Nguyen); Hasenfuss
  1992 instead finds the cross-bridge force–time integral **+33 %** with tension-dependent heat
  −61 %: *longer attachment per ATP*, not slower myosin. That is exactly the k2-slowing picture
  this project already uses for low ATP — worth keeping in mind when a human HF scenario is set.

---

## 6. Solution and single-molecule kinetics: human β vs rodent α

Conditions: **G** = Geeves lab, 100 mM KCl, 20 °C; **G25** = Geeves 25 mM KCl 20 °C; **S** = Spudich
lab, 25 mM KCl, 23 °C, short S1 (sS1); **Y** = Yengo lab, 20 mM KCl, 1 mM Mg, 25 °C. All "human" myosin
is C2C12-expressed; rodent/porcine/bovine are tissue-purified.

| step | human β | rodent α (or large-mammal β where α is absent) | ratio β/α | confidence |
|---|---|---|---|---|
| ADP release from acto·S1 (model **k2**, saturating ATP) | **93 ± 5** (Deacon 2012, G); 64 (Bloemink 2014, G); 59 (Johnson 2021, G25); 72 (Sung 2015, S); 87–102 single-molecule k0 (Sung, Liu 2018); **174–332 (Tang 2021, Swenson 2017, Y — outlier)**; bovine β HMM 35 at 25 °C (Rohde 2017); pig β 35, rat β 135 at 15 °C (Pereira) | **> 1000–1250 s⁻¹** (Deacon 2012, mouse and human α, unresolved); Lowey 2013 mouse α-S1 **~350** (via `parameterBounds.m`) | **≤ 0.1 (0.07–0.25 depending on which α number)** | direction high; magnitude a bound |
| ADP affinity K_AD | 6–21 µM | 150–250 µM (α) | 0.1 | medium |
| ATP binding to acto·S1 K1k+2 | 1.1–1.6 µM⁻¹ s⁻¹ (G); ~3 (HMM trap, 25 mM KCl) | 2.0–2.5 | 0.5–0.7 | medium |
| max ATP-induced dissociation | 1000–1450 s⁻¹ | > 1000–1670 | ~1 | medium |
| hydrolysis k+3+k−3 (model **kah**) | **17 ± 2** (Deacon), 14 (Bloemink); Yengo recovery-stroke 182 (discordant) | **150–168** (mouse α tissue, human α) | **~0.1** | high within Geeves lab |
| Pi release / power stroke (model **k1**) | 12–20 s⁻¹ (Mijailovich 2017 kPi 16; Rohde bovine 11.6; Tang 17 at 25 °C — note Tang's Pi release 3 → 19 s⁻¹ from 25 → 35 °C) | **not measured for α**; modelled 32 (Mijailovich 2017), "~5× faster" assumed (Mijailovich 2022) | ~0.2–0.5 (assumed) | low |
| kcat (actin-activated) | 4.5–8.8 (6 typical at 20–23 °C) | human α 18 (G25); same-species α/β 2–3 | 0.2–0.33 | medium |
| duty ratio, unloaded | 0.08–0.11 | 0.11–0.14 | 0.8–1 | medium |
| unloaded velocity | 0.5–0.8 µm s⁻¹ (23 °C, S1) | α/β 2–3 × rodent/human 1.7–3 | 0.2–0.3 | medium |
| step size / stiffness / force | 5.4–7 nm; ≥ 2 pN nm⁻¹ (Woody 2019); 1.3–1.4 pN (Spudich) | 7 nm, same force | ~1 | medium |
| load sensitivity δ of detachment | 0.8–1.3 nm (Sung, Liu, Woody, Greenberg porcine 0.97) | not measured for α | – | – |
| temperature | k−ADP Ea 90–106 kJ mol⁻¹ (**Q10 3.4–4**); human β k−ADP 84 → 575 s⁻¹ from 20 → 37 °C (Deacon); motility Q10 3.1–3.5 | rat V0 Q10 4.6 | | medium |

Reading it:

- **ADP release is the one step where human β and rodent α differ by an order of magnitude**, and it
  is load-sensitive (δ ≈ 1 nm) — this is the model's k2 and its negative-strain PCHIP. Mijailovich 2017
  makes the point that for β the rate-limiting ADP step is the preceding isomerisation kD\* = 59 s⁻¹
  (k−D\* 354), not ADP dissociation itself (1000 s⁻¹ in both isoforms): map the model's k2 to kD\*.
- **Hydrolysis/recharge is 10× slower in β** — but the model's kah is already its slowest-sensitivity
  lever on k_tr (elasticity −0.09) and only moves F0 (+0.55). It mostly sets how much of the pool is
  parked in PT.
- **α-cardiac Pi release / power-stroke rate has never been measured**; every "α is 3–5× faster at k1"
  is a modelling assumption. ADP release for α is only a lower bound (> 1000 s⁻¹) — there is no
  resolved modern value, so the β/α ratio for k2 is likewise a bound.
- Inter-lab spread is large: human β k−ADP 59–332 s⁻¹, kcat 4.5–8.8, velocity 0.47–1.8 µm/s, driven
  by KCl (kcat drops ~60 % 20 → 50 mM, 80–90 % at 75–100 mM; Ge/Yengo 2026), Mg, temperature and
  construct (sS1 lacks the lever; tissue bovine HMM is 2–3× slower than recombinant human S1).
  Muscle at 180 mM ionic strength sits closest to the Geeves 100 mM KCl numbers.
- Attachment rate ka is not solution-measurable. Mijailovich 2022's fitted k+A 330 (mouse) vs 60 (human)
  is a model value. The constraint on ka comes from k_tr and the duty ratio: with duty ratio ~equal
  across isoforms and all exit rates ×0.1, **ka must scale ×~0.1 too**.
- Ranatunga/Linari-type tension transients: no human phase-2 data; skip.

---

## 7. The super-relaxed state: mant-ATP numbers vs what a mechanical model needs

### 7.1 mant-ATP chase measurements (cardiac)

| study | species (MHC) | prep, T | T1 (s) | T2 (s) | slow-phase amplitude "SRX" |
|---|---|---|---|---|---|
| Hooijman, Stewart & Cooke 2011 | rabbit V (β) | skinned bundles, 26 °C | 14 ± 1 | 144 ± 10 | 0.27 ± 0.01 (reverse chase 0.29; mant-ADP 0.09; pCa 5.7 0.26) |
| McNamara 2017 *PLoS One* | **human** LV donor / HCM-smn / MYBPC3 | skinned strips, RT | 14.3 / 13 / 14 | 224 / 211 / **160** | **0.276 / 0.281 / 0.223** |
| Anderson 2018 *PNAS* | human fibres (± mavacamten); human β 25-hep, 2-hep HMM, **sS1** | – | 22 s (0.046 s⁻¹) | 172 s (0.0058 s⁻¹) | sS1 (no IHM possible) still biphasic |
| Walklate 2022 *JBC* | porcine V (β) / bovine masseter / rabbit psoas / human iPSC-CM | myofibrils, 20 °C | 15 / 37 / 11 | 455 / 278 / 143 | 0.28 / 0.26 / 0.42 / 0.42 |
| Ma 2021 *Circ Res* | porcine LV | skinned, 22 °C, SL 2.0 → 2.3 | 23 → 17 | 156 | 0.154 → 0.120 (−3.4 pp) |
| *JGP* 2023 (PMC9884578) | **mouse** (α) | myofibrils, single-molecule BODIPY-ATP | 25 | 161 | **0.57 ± 0.05** (C-zone 0.71, D-zone 0.45; cMyBP-C-null 0.35) |
| Ma, Chaami, Ochala 2026 *JACC BTS* | human NF / DCM / ICM septum (792 strips) | 25 °C | shorter in DCM | shorter in DCM | lower in DCM/ICM (figures only) |

Unverified beyond secondary sources: Stewart 2010, McNamara 2016 (mouse), Toepfer 2020 numbers.

### 7.2 Why these numbers are not model rates (the prior holds, with one exception)

| finding | source |
|---|---|
| "ghosted" fibres with 95 % of myosin removed still show a fast phase comparable to DRX → the fast phase is diffusion/non-specific binding | Montesel 2025 *PNAS* |
| correcting for bleaching and FRET removes the slow phase; a single exponential suffices | Cail 2025 (via *Biophys Rev* 2026, PMC13486409) |
| purified porcine cardiac HMM shows a single exponential; "not a reliable assay to quantify the SRX population" | Mohran 2024 *JBC* |
| SRX forms within 0.2–1 s of mant-ATP binding to rigor (entry ~10 s⁻¹) yet turns over at 0.002–0.005 s⁻¹: two pools interconverting in < 1 s cannot be kinetically isolated over 100 s | Mohran 2024, Walklate 2022 |
| fitted slow-phase amplitude under-estimates true occupancy (27–28 % fitted ↔ 36 % true; 22–23 % ↔ 42 %; 58 % with bleaching) because incubation (0.2–120 s) ≪ 5–10 min equilibration; fitted k_DRX→SRX = 4.3–10 s⁻¹, k_SRX→DRX = 2–6 × 10⁻³ s⁻¹ | **Ježek, Han, Vander Roest & Beard 2025**, bioRxiv 2025.08.06.667936 — this group's own model |
| omecamtiv and piperine move the X-ray OFF/ON markers dramatically with **no change in mant SRX/DRX** → biochemical SRX ≠ structural OFF | Jani 2024 *PNAS Nexus* |
| sS1 (single-headed) is biphasic → the slow phase does not require the IHM; may be slow Pi release of a blocked head | Anderson 2018, Chu/Muretta/Thomas 2021, Walklate 2022 |
| rat trabeculae at pCa 6.6, 26 °C: stretch recruits folded motors at **~5 s⁻¹** (RLC probes); at diastolic Ca stretch alone barely moves them; folded state half-destabilised at 16–18 °C | Park-Holohan 2021 *PNAS* |
| a human-LV LDA model uses k_OFF→ON = 0.88 s⁻¹ × (1 + 1.44 kPa⁻¹ · F), k_ON→OFF = 12.6 s⁻¹ | PMC11428202 |

Verdict for this model:

1. **Exit rate (SRX→DRX).** 1/T2 (0.002–0.007 s⁻¹) is the ATPase of an unloaded, cool, resting ensemble.
   The mechanically relevant OFF→ON transition is ~1 s⁻¹ unloaded rising to 5–90 s⁻¹ under
   stress. The incumbent's `kmsrd` 27 s⁻¹ × e^(F/22.8 kPa) is in that regime. **Do not import T2.**
2. **Entry rate (DRX→SRX).** The chase-implied 4–13 s⁻¹ is compatible with mechanical ON→OFF rates
   (12.6 s⁻¹). The incumbent's `ksr0` 367 × e^(−F/27 kPa) is ~30× faster at zero force; that is a
   modelling choice (`SRXDuringContraction` shows the data need a *fast* exchange), not a human/rodent
   issue.
3. **Fraction.** Human resting SRX ≈ 0.28 fitted, **0.35–0.5 true** after correction. Rodent α-MHC
   fibre chase data essentially do not exist (the one mouse value, 0.57, is a different single-molecule
   assay with C-zone weighting). **No defensible human-vs-rodent SRX-fraction difference**: β preps
   from four species all sit at 0.15–0.30 fitted. The incumbent's rest fraction 0.27 is consistent
   with human as it stands.
4. **Temperature** is unresolved and contradictory (Walklate 2022: flat 5–30 °C; *Biophys Rev* 2026:
   rising 12–35 °C; RLC probes: folded state collapses below ~16 °C). All cardiac chase data are at
   20–26 °C. A 37 °C human model probably needs *more* OFF state than the 15–26 °C data suggest.
5. cMyBP-C is the strong modulator in both species (mouse null 0.57 → 0.35; human MYBPC3 0.28 → 0.22,
   T2 224 → 160 s); its phosphorylation falls to 23 % of normal in human HF (Copeland 2010).

So for humanisation: **leave the SRX/OFF rates and force sensitivities alone** (no species evidence
either way); consider scaling them with the cycle only as a sensitivity arm.

---

## 8. Human myocardium at 37 °C, metabolites, heart failure

### 8.1 Intact human kinetics (what a 37 °C human model must reproduce)

| quantity | human non-failing | failing | rat / mouse | source |
|---|---|---|---|---|
| twitch TTP / RT50 at 0.5 Hz, 37 °C | 220–235 / 140–153 ms | 211 / 153 | rat 54 / 29; mouse 36 / 20 ms | Milani-Nejad & Janssen 2014; 2015; Chung 2019 |
| RT50 / RT90 at 3 Hz | 90 / 157 ms | 101 / 172 | – | Chung 2019 |
| dF/dt/F | 7–11 s⁻¹ | – | rat 30, mouse 48 | same |
| k_tr (K⁺ contracture) | **18 s⁻¹** (frequency-independent) | slower only at 3 Hz | rat 28 | Milani-Nejad 2015, Chung 2019 |
| developed force | 14–18 (twitch), 36 (K⁺), 25–44 (Alpert lab) mN mm⁻² | −46 to −80 % at high rate | rat 28, mouse 40 | |
| force–frequency | +76–110 % 60 → 150 min⁻¹ (Alpert lab) vs +27 % at 2 Hz (Chung) | flat/negative | – | discordant between labs |
| XB force–time integral (myothermal) | reference | **+33 %**; tension-dependent heat −61 % | FTI ∝ %V3 across species (r 0.97) | Hasenfuss 1991/92 |

Twitch rates are 4–5× slower than rat but dominated by Ca handling; the cross-bridge-level number
is k_tr 18 vs 28 (×0.65 — compressed, §2.2).

### 8.2 Metabolites (31P-MRS, human in vivo)

| quantity | normal | HF / DCM | source |
|---|---|---|---|
| ATP | 5.7 ± 1.0 mmol kg⁻¹ ww | **−35 %** in DCM (3.7) — single 10-patient SLOOP study, never replicated absolutely | Beer 2002 *JACC* |
| PCr/ATP | 1.6–2.1 (method-dependent) | 1.3–1.6 | Neubauer 1997; Valkovič 2019 |
| free ADP | **45–92 µM** | **lower**, 45–60 µM (creatine loss) | Weiss 2005 *PNAS*; Hirsch 2012 |
| Pi | ≈ 1.2–1.4 mM (Pi/PCr 0.10–0.11) | ~2× in HCM; no HFrEF 7 T value | Valkovič 2019, Apps 2021 |
| pH | 7.09–7.12 | – | same |

For a "human HF" scenario: ATP 5.7 → ~3.7 mM, Pi 1 → 2–3 mM, ADP 50–100 µM, pH unchanged. The
2 mM ATP of the Baker protocol is below any resting human HF value; it is an ischaemic/stress regime.
Human β's 10× tighter ADP affinity (K_AD 6–21 µM vs 150–250 for α) means 50–100 µM ADP already sits
at 3–10 × K_AD for human but ≪ K_AD for rodent — **ADP is a human lever that the rodent fit cannot
see** (and the model's g1 ADP path is dead, `atp-sensitivity-mechanism`).

### 8.3 Human HF myofilament changes relevant to the model

| change | NF → HF | source | model relevance |
|---|---|---|---|
| MHC isoform | ~95 % β → ~100 % β (no shift) | Miyata 2000, Narolska 2005 | none |
| myosin ATPase / economy | unchanged (0.3 s⁻¹ head⁻¹; tension cost preserved) | Nguyen 1996; Narolska 2005 *JMRCM* | do not slow myosin for HF |
| XB force–time integral | +33 % (longer attachment) | Hasenfuss 1992 | k2-slowing-like |
| cMyBP-C phosphorylation | 2.7 → 0.6 mol Pi/mol (23 %) | Copeland 2010 | SRX/OFF stabilisation; k_tr ↓ |
| TnI phosphorylation | −80 % (donor baseline inflated by catecholamines / ex-vivo incubation) | Messer 2007/09; Cai 2018 | Ca sensitivity, not XB |
| RLC phosphorylation | contradictory (↑ Toepfer 2013; ↓ van der Velden 2001/03) | | unresolved |
| titin N2BA:N2B | 0.43 → 0.7–0.9; passive tension −25–30 % | Makarenko 2004; Neagoe 2002 | passive; rat ≈ 0 N2BA |
| SRX (mant) | ↓ in DCM/ICM, T1/T2 shorter | Ma 2026 | opposite sign to MyBP-C dephosphorylation expectation; figures only |
| Pi / pH sensitivity | 30 mM Pi → 0.33 F; pH 6.2 → 0.64 (15 °C) | van der Velden 2001 | Pi lever (f1/f2 not wired in active path) |

---

## 9. Conflicts, superseded values, and things I could not verify

- **Species of the Baker data** is inconsistent in the repo (mouse vs rat); the rig temperature is
  unrecorded. Both change the human/rodent ratio by ~1.5–2×.
- **Little 2012** ("human 40 s⁻¹, rat 150 s⁻¹", cited in `sources.md` and `parameterBounds.m` for `kd`)
  is a cTnC-reported *cross-bridge detachment* in myofibrils at 15 °C; it is closer to the model's
  post-stroke k2 than to its pre-stroke kd. Ratio 0.26, consistent with ADP-release-limited detachment.
- **Human β ADP release**: Yengo-lab 174–332 s⁻¹ (25 °C, low salt, mant-ADP) vs Geeves/Spudich 59–102
  (20–23 °C). Temperature explains ×1.8 of it. Use 60–100 at 20 °C.
- **α ADP release** has only a lower bound (> 1000); Lowey 2013's 350 is the one resolved rodent-α
  value and is what `parameterBounds.m` uses. The β/α ratio is therefore 0.07–0.25; I use 0.1.
- **Human k_tr Q10**: constant-Q10 extrapolation from 15 °C under-predicts 37 °C by ~2×. Use
  measured anchors at each temperature instead of extrapolating across > 10 °C.
- **Intact 37 °C k_tr** (human 18, rat 28) compresses the ratio to 0.65 — activation-limited; not a
  cross-bridge number.
- **SRX temperature dependence** is contradicted between Walklate 2022 and the 2026 review.
- **Mant-ATP slow phase** is partly a bleaching/diffusion artefact (Montesel, Cail, Mohran); amplitudes
  are biased low (Ježek/Beard 2025). Numbers in §7.1 are "as fitted".
- Not obtained (paywalled/captcha): Hasenfuss 1991 per-species cycling rates; de Tombe & Stienen 2007
  absolute values (taken from an earlier analysis in this repo); Awinda 2020 / Donaldson 2012 human
  37 °C 2πb/2πc; Narolska 2005 absolute ATPase (via review); Belus 2008 control values; Piroddi 2007
  ventricular tissue source; Siemankowski 1985 species table; Stewart 2010 / Reconditi 2017 /
  Brunello 2020 numbers.
- **No Kawai-lab human ventricular sinusoidal analysis exists** (search returned mouse/bovine only).
- Donor hearts are not healthy hearts: brain death (catecholamines), age (Hasenfuss FTI r = 0.86 with
  age), heart weights 271–896 g (Chung 2019), and TnI phosphorylation that rises during ex-vivo storage.

---

## 10. Human measurements the model can be run against

Ordered by how directly they map onto an existing protocol/feature of this code.

| # | measurement | human value | T °C | maps to | status |
|---|---|---|---|---|---|
| 1 | k_tr, max Ca, skinned (myofibril / cell / strip) | 0.55–0.85 / 0.6–1.3 s⁻¹ (use 0.8) | 15 | `ktrProbe`, `ktrP_k80` | best-constrained human number |
| 2 | k_tr intact, K⁺ contracture | 6.6 / 14.1 / 18.1 s⁻¹ | 27 / 32 / 37 | same | ceiling for a 37 °C model |
| 3 | stretch-activation k2 (delayed rise) | 0.74 s⁻¹ | 15 | restretch recovery `rsK80` | one lab (van Dijk 2014) |
| 4 | tension cost | 2.3–2.4 µmol l⁻¹ s⁻¹ per kN m⁻² (≈ 0.5 s⁻¹ per head at 40 kPa) | 20 | `XTOR` net / F0 | robust (two datasets) |
| 5 | slow kREL (≈ isometric detachment) | 0.29–0.30 s⁻¹ | 15 | k2eff at isometric | three datasets |
| 6 | fast kREL (relaxation after sarcomere give) | 4.2–4.7 s⁻¹ | 15 | lengthening-side detachment | definition-dependent (Vikhorev 17) |
| 7 | V0 | 0.4–0.6 ML/s (inferred) | 15 | `Vmax_hill` / zero crossing | **no direct human datum** |
| 8 | peak power, myocyte fragments | 54 pW; 1.23 µW mg⁻¹ (rat α 12 °C 1.49) | 15 | `Pmax` | needs raw FV from McDonald lab |
| 9 | FV curvature a/P0 | ≤ 0.1 (rat β 0.05, rat α 0.10 at 12 °C) | 12–15 | `a_over_F0` | no human value |
| 10 | F0 skinned | 35–47 kN m⁻² (cells/strips), 100–110 (myofibrils) | 15–20 | `steady`, `F0` | prep-dependent; do not fit |
| 11 | resting SRX (mant, fitted / corrected) | 0.28 / 0.35–0.5 | RT | rest `SRX` | assay-limited (§7) |
| 12 | Pi and pH sensitivity | 30 mM Pi → 0.33 F0; pH 6.2 → 0.64 | 15 | f1/f2 (not wired) | for the HF scenario |
| 13 | twitch TTP / RT50 at 37 °C | 220 / 140 ms (0.5 Hz), RT50 90 ms (3 Hz) | 37 | needs Ca activation (`CalciumTwitch`) | out of reach today |
| 14 | HF: XB force–time integral +33 %, economy unchanged | | 37 | k2 ↓ with kah/ka fixed | scenario |

Data that would make the human fit a real fit rather than a rescaling: a human skinned-trabecula
slack–restretch and FV battery at a known temperature (e.g. 15 and 25 °C), ideally from the same
rig as the rodent data. Lewalle 2025 (37 °C human skinned fibres, k_tr by t½) and the McDonald-lab
human myocyte FV are the two existing datasets closest to this model's protocols.

---

## 11. What the model does not represent and the human literature says matters

- **Thin-filament cooperativity / Ca activation** — absent (saturating Ca). Human kACT ≈ kTR says this is
  safe at maximal Ca, but every human twitch/relaxation number (§8.1) is unreachable without it.
- **cMyBP-C** — no C-zone, no phosphorylation state; the largest human HF myofilament change (§8.3).
- **Titin isoform** — human 30 % N2BA vs rat ~0; the model's passive is identified on rodent.
- **Temperature** — no Q10 anywhere in the code; every rate is a number at an unrecorded temperature.
- **ADP** — the ADP-release inhibition path is dead; human β's K_AD (6–21 µM) makes physiological ADP a
  real lever that the rodent fit never exercised.
