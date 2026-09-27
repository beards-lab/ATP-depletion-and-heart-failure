# Can free [MgATP] change the hydrolysis rate (kah)?

## Question
H1 gates PT→PD (`kah`, 176 s⁻¹ at 8 mM) and its reverse (`kamh`) by (MgATP/8)^p, with p_kah = 0.30, so kah is ×0.66 at 2 mM. Is that defensible? `AtpEffectHypotheses/conclusions.md` calls kah "unambiguous ATP chemistry"; the `RunPlausibleSets.m` header says "NOT PLAUSIBLE".

## Chemistry
Hydrolysis (recovery stroke plus cleavage) is a **unimolecular step on nucleotide that is already bound**. Free ATP enters only upstream, in M + ATP ⇌ M·ATP ⇌ M*·ATP (K1, k+2) [1,2]. The transient k_obs rises hyperbolically with [ATP] only because binding comes first, then plateaus at k+3 + k−3.

Human S1, 20 °C [3]:

| | K1k+2 (µM⁻¹s⁻¹) | k+2 (s⁻¹) | k+3+k−3 (s⁻¹) |
|---|---|---|---|
| β | 1.5 | 158 | 17 |
| α | 2.2 | >200 | 168 (mouse α: 150) |

The β slow phase is roughly ATP-independent above about 50 µM [4]. Beard 2022 used mouse trabeculae (α-MHC) [5], so kah ≈ 176 s⁻¹ is consistent with α-S1.

- **Detached heads:** binding at 2 mM is capped by k+2, not [ATP]. It changes by less than 5% between 2 and 8 mM.
- **Actin·S1 (rigor):** 1/K′1 ≈ 0.2–1.1 mM, k′+2 ≈ 1000–1500 s⁻¹ [3]. That makes binding 0.73–0.93× as fast at 2 mM, still far faster than ADP release (β 93 s⁻¹).
- **Fibres (ferret, sinusoidal):** K1a ≈ 1 mM⁻¹ on *attached* heads [6], so occupancy is 0.66 vs 0.89 (ratio 0.75). That scheme gives hydrolysis no MgATP dependence.

## Evidence
- **Rabbit skinned fibres:** caged-ATP plus rapid freezing gave a cleavage burst of about one per head within 50 ms, at any Ca²⁺ or overlap [7].
- **Low ATP (2.5 µM–5 mM):** Cooke & Bialek explained every effect through ATP-induced detachment [8].
- **Beard 2022** (ionic strength 200 mM with KCl, 12 mM CrP): a model with ATP only at detachment failed, and they proposed an SRX/thick-filament nucleotide mechanism [5].
- **I found no report of a millimolar-[MgATP] dependence of k+3/k−3.**

## Indirect routes
- **(a) Second site:** one catalytic site per head. I found no accepted allosteric ATP site.
- **(b) Solution composition:** ionic strength is balanced [5]. Free Mg²⁺ inhibits β-cardiac myosin by slowing ADP release [9], not hydrolysis.
- **(c) Products:** K3 is an on-enzyme equilibrium (≈3.5 for cardiac S1 [10]), independent of free ATP, ADP or Pi. The cycle's ΔG shift (RT ln4 ≈ 3.4 kJ/mol) belongs in the binding step.
- **(d) SRX:** mant-ATP chase measures turnovers of bound nucleotide (about 0.005 vs 0.05 s⁻¹) [11,12]. SRX occupancy depends on nucleotide *type* (GTP, dATP) [12,13]. No 2 vs 8 mM data found.
- **(e) Thick filament:** dATP moves relaxed heads toward actin [13]. A low-ATP analogue is unverified.
- **(f) Hidden apo pool in PT:** the effective priming rate would be 1/(1/k_bind + 1/kah).
  - With solution binding rates, the loss at 2 mM is ≈2–3%.
  - With fibre-like binding (0.66 × 270 s⁻¹) it is ≈13%, and only if all of PT is apo.
  - Reproducing ×0.66 needs all of PT apo, binding still unsaturated at 8 mM, and about 0.1 µM⁻¹s⁻¹. That is 15× below the measured rate and contradicts [3,6].

## Verdict
**Gating kah by free [MgATP] is not physically defensible.** Hydrolysis acts on bound ATP, and the [ATP] dependence sits in binding, which saturates below about 1 mM. "Unambiguous ATP chemistry" is wrong; the RunPlausibleSets header is right. Gating **kamh** is even less defensible: reverse cleavage depends on bound products, and gating the two unequally would shift K3 without cause. A hidden apo pool could explain at most a 5–13% apparent effect.

## What the model is compensating for
Slowing kah slows re-priming, so it slows the attachment ("f") limb. That is the forward-limb partner k2 lacks, which explains why ktr drops 37% while F0 rises 18%. Plausible real carriers:
1. Missing strain dependence of ADP release and ATP binding: restretched heads detach more slowly at low ATP, lowering ktr.
2. Missing nucleotide- or strain-dependent thick-filament OFF/SRX recruitment [5].
3. A missing apo pool after forced detachment (small; see (f)).
4. The k2 lumping of ADP release with ATP binding.

Test (1) or (2) next, with kah fixed.

## References (DOIs verified via Europe PMC)
1. Bagshaw et al. 1974 Biochem J 141:351. doi:10.1042/bj1410351
2. Bagshaw & Trentham 1974 Biochem J 141:331. doi:10.1042/bj1410331; Johnson & Taylor 1978 Biochemistry. doi:10.1021/bi00610a002
3. Deacon et al. 2012 Cell Mol Life Sci 69:2261 (Table 1). doi:10.1007/s00018-012-0927-3
4. Bloemink et al. 2014 J Biol Chem. doi:10.1074/jbc.M113.511204
5. Beard et al. 2022 Biophys J. doi:10.1016/j.bpj.2022.07.029
6. Kawai, Saeki & Zhao 1993 Circ Res. doi:10.1161/01.res.73.1.35
7. Ferenczi, Homsher & Trentham 1984 J Physiol. doi:10.1113/jphysiol.1984.sp015311
8. Cooke & Bialek 1979 Biophys J. doi:10.1016/S0006-3495(79)85174-7
9. Swenson et al. 2014 J Biol Chem. doi:10.1074/jbc.M114.562231
10. Hazzard & Cusanovich 1986 Biochemistry. doi:10.1021/bi00373a003
11. Hooijman, Stewart & Cooke 2011 Biophys J. doi:10.1016/j.bpj.2011.02.061
12. Stewart et al. 2010 PNAS. doi:10.1073/pnas.0909468107
13. Ma et al. 2020 J Physiol. doi:10.1113/JP280402
- Siemankowski & White 1984 J Biol Chem. doi:10.1016/S0021-9258(17)42953-X (ADP release limits cardiac detachment)

## Could not verify
- β hydrolysis assignment ("provisional" per [3]).
- Temperature and free [Mg²⁺] of the Baker-lab and Beard solutions.
- 2 vs 8 mM SRX or X-ray data.
- Absence of a second nucleotide site (not searched exhaustively).
- McNamara 2015, Malnasi-Csizmadia and Nag/Spudich values (not retrieved).
