# Literature background for issue #82 (skinned-trabecula XB model into a P21 rat TriSeg host, 37 C, HR ~325 bpm)

Status: web-search check, 2026-11-09. Tags: [V] = number read from the source page/abstract this session;
[S] = from a search snippet only; [M] = from my memory, NOT verified (DOIs especially); [calc] = my own arithmetic.
PMC full texts were blocked (bot check) so most Granzier/Janssen numbers could not be read. Nothing below is a
substitute for opening the papers.

## 1. In-situ sarcomere length (rat/mouse LV)
- Rat LV in situ, open chest: SL rises from 87% to 96% of optimum (optimum 2.07 um in papillary muscle, i.e. ~1.80 -> ~1.99 um) as intraventricular pressure goes 0 -> 24 cmH2O, then plateaus; K-arrested, thorax intact: 92-95% of optimum (~1.9-1.97 um). [S] Authors/journal not captured (1968-era study, PMC2098739 area; a companion study shows SL rising almost linearly endo->epi at 6/12/24 cmH2O). Find and cite before use.
- Pearson JT, Shirai M, Tsuchimochi H, Schwenke DO, Ishida T, Kangawa K, Suga H, Yagi N. Biophys J 2007;93:4319-4329. doi:10.1529/biophysj.107.111740 [V]. Rat (SD, 9-10 wk) in situ X-ray: SL NOT measured; inferred ED epicardial SL 1.75-1.90 um. EDP 9.2 +/- 2.7 mmHg, ESP 113.5 +/- 20 mmHg (baseline). No ES SL.
- Mouse in-vivo SL: nothing found. End-systolic SL in rat/mouse: nothing found (do not quote a number).
- Implication: ED SL in rat is ~1.8-2.0 um, well below the 2.2-2.3 um of skinned-prep calibration, and below the rodent optimum ~2.07-2.2 um. Marzban 2020 (below) uses SL_rest 1.8 um, consistent.

## 2. LV end-diastolic wall stress (rat/mouse)
- No primary value found. EDP rat ~9-10 mmHg (1.2-1.4 kPa) [V, Pearson 2007].
- [calc] With Marzban 2020 rat sham geometry (LV wall volume 0.56-0.85 mL [V]) and an assumed EDV ~0.3 mL (assumption, not sourced), thin-sphere Laplace sigma = P*(Vlv/Vw)/3 ~ 0.2 kPa fibre stress; i.e. order 0.1-0.5 kPa, much lower than the ~2-3 kPa typical of canine/human. Treat as order-of-magnitude; check against the TriSeg model's own sigma_f at ED.

## 3. Passive stress, titin vs ECM, Ca effect
- Granzier H, Irving TC. Biophys J 1995;68:1027-1044 (doi 10.1016/S0006-3495(95)80278-6 [M]). Rat cardiac. Over the working range ~1.9-2.2 um titin dominates at short SL, collagen at long SL; [S] no kPa numbers read.
- Granzier H, Helmes M, Trombitas K. Biophys J 1996;70:430-442. doi:10.1016/S0006-3495(96)79586-3 [V]. Titin I-band elastic segment straightens 1.85-2.0 um, domains unfold above ~2.0 um: passive force rises steeply above 2.0 um.
- Granzier-lab mouse LV, titin vs ECM (J Mol Cell Cardiol 2011) [S]: titin dominates passive pressure within physiological volumes (~80% quoted), ECM dominates at larger volumes; midwall SL range 1.8-2.2 um. Exact citation/DOI and kPa not verified.
- Skinned passive stress as a bound: skinned/permeabilized strips retain titin and ECM but lose intracellular-pressure/cell-swelling and the in-vivo prestress; in situ ECM load is reported to dominate at high volume [S above]. Whether skinned is a true lower bound is NOT established by anything I read. Note preps are normally measured at 12-25 C (Fukuda 2005: 12/22/25 C) and passive force is temperature-sensitive.
- Ca effect: Fujita, Labeit, Gerull, Labeit, Granzier, "Titin isoform-dependent effect of calcium on passive myocardial tension", Am J Physiol Heart Circ Physiol (pp. H2528-H2534; year/DOI unverified). [S] Ca raises passive tension in N2BA-rich bovine LA significantly but little/not in N2B-rich rat RV; size not read. Rodent (N2B-dominant, ~90:10 N2B:N2BA in rat RV [V, Fukuda 2005]) => small Ca-titin effect expected.
- Fukuda N, Wu Y, Nair P, Granzier HL. J Gen Physiol 2005;125:257-271. doi:10.1085/jgp.200409177 [V]: rat/bovine, PKA phosphorylation lowers passive force most in N2B-rich RV.
- Cazorla O, Wu Y, Irving TC, Granzier H. Circ Res 2001;88:1028-1035 [S]: titin-based modulation of Ca sensitivity (active, mouse skinned).

## 4. Vmax, Q10, in-vivo sarcomere velocity
- de Tombe PP, ter Keurs HEDJ. Circ Res 1990;66:1239-1254 (doi 10.1161/01.RES.66.5.1239 [M]). Rat trabeculae, 15-35 C: Vo (unloaded sarcomere shortening velocity) linear 18-33 C, Q10 = 4.6 [S, abstract].
- de Tombe PP, Stienen GJM. J Physiol 2007;584:591-600. doi:10.1113/jphysiol.2007.138693 [V]. Permeabilized rat trabeculae 15-25 C: Q10 of ktr and of Ca-activated ATPase at saturating Ca ~3.5; Tmax ~1.3; tension cost ~2.5. ktr at ~5% activation has lower Q10 (= tension cost).
- Practical: 22 C -> 37 C with Q10 3.5 gives ktr x ~5.1 (4.6 gives Vo x ~8); Tmax x ~1.5. Extrapolating a 15-25 C Q10 to 37 C overestimates if Q10 falls at higher T; no mouse-specific Q10 found.
- Peak sarcomere shortening velocity in intact rat/mouse myocytes at 37 C or in vivo: NOT found; do not quote. (Permeabilized single rat myocytes, Vmax ~2.8 ML/s, temperature unknown [S]. Dataset note: our sawtooth Vmax 7.5 ML/s, from project memory.)

## 5. Host model parameter sets
- Marzban B, Lopez R, Beard DA. "Computational modeling of coupled energetics and mechanics in the rat ventricular myocardium." Physiome 2020, doi:10.36903/physiome.12964970 [V]; primary paper Function 2020, zqaa018 [V per Physiome page]. Code: github.com/beards-lab/Rat-Cardiac-Energetic. Rat, 37 C, TriSeg (Lumens) geometry; sham HR 318.6-367 bpm (TAC 294-349) [V]; LV Vw 0.564-0.853 mL, septal 1/3 of LV, RV 0.282-0.373 mL; LV midwall ref area 1.69-2.47 cm2 [V]. Myofilament: Campbell-type Ca activation (kon=101.2 /uM/s, koff=101.2 /s, kcoop 1.857; Ca50 5.89 / 6.00 uM-type values at SL 1.9/2.3, n 4.6/4.5) [V]. SL_rest 1.8 um, collagen threshold 2.25 um, k_passive 25 mmHg/um, K_SE 5e4 mmHg/um [V]. Ca transient re-estimated to Janssen 2002 rat trabeculae at 37 C, SL 2.2 [V]; its amplitude/time-to-peak/tau were not on the page, read the code.
- Lumens J, Delhaas T, Kirn B, Arts T. Ann Biomed Eng 2009;37:2234-2255 (doi 10.1007/s10439-009-9774-2 [M]): TriSeg.
- Tewari 2016 mouse multiscale model: not found/verified. A 2023 bioRxiv preprint (PMC9900958) reports a mouse multiscale model coupled to TriSeg with baseline and acute-ischaemia data [S]; authors/citation not read: check as a mouse 37 C host.
- Colebank/Olufsen P21 rat hypoxia TriSeg: not found by search. Not verified that it exists; ask the authors/check Olufsen-lab publications.
- Beard lab 'multiscale-ventricular-interactions-model': not searched successfully; open the GitHub repo.
- CircAdapt: no mouse/rat parameter set found; nothing linking it to the above.

## 6. Mouse LV Ca transient at 37 C
- Stull LB, Leppo M, Marban E, Janssen PML. J Mol Cell Cardiol 2002;34:1367-1376 (doi not verified) [S]: ultra-thin mouse RV trabeculae at 37 C, iontophoretic fura-2; Ca transient amplitude tracks force; kinetics numbers not read.
- Gao WD, Perez NG, Marban E. [S]: 37 C mouse trabecula experiments failed for indicator loss, most data at 20-22 C (full citation not captured).
- No verified mouse time-to-peak / tau / amplitude at 37 C, 8-10 Hz. Rat intact-trabecula decay tau 52-63 ms in WKY/SHR at 1 Hz, temperature unknown [S]. Do not use as mouse values.
- Read the Marzban repo's Ca transient for the rat 37 C shape, and the Stull 2002 full text for mouse.

## Gaps to close before relying on this
No in-vivo ES SL; no wall-stress primary source; no peak in-vivo sarcomere velocity; no mouse Ca transient numbers; no Colebank P21 rat paper; Q10 for mouse unavailable.
