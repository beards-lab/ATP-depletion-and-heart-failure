# Abstract Draft — Myofilament Meeting

---

**A Mechanistic Cross-Bridge Model of Cardiac Muscle Mechanics Under Energetic Stress**

Filip Jezek¹, Anthony Baker², Dan Beard¹
¹ University of Michigan Medical School, Ann Arbor, MI, USA
² San Francisco VA Health Care System and Department of Medicine, San Francisco, CA, USA

---

Most forms of heart failure are accompanied by energetic deficit — impaired mitochondrial ATP production, elevated free ADP, and decreased PCr/ATP ratio — that compromises myosin cross-bridge cycling. Computational models linking mitochondrial metabolism to sarcomere mechanics [1] offer a path to understanding this coupling. A strain-discretized cross-bridge model by Beard et al. [2] applied on distinct ATP force-velocity data does not reproduce ATP-dependent contractile differences. To address this, we extend the model by including additional physical mechanisms and furthermore test the model with new data sets that include slack-restretch protocols across ATP levels.

Experiments on mouse right-ventricular cardiac trabeculae provide force-velocity and slack-restretch data at 8 mM (saturating), 2 mM (compromised), and 0.2 mM (depleted) ATP. Model calibration targets discrete features, including unloaded shortening velocity (Vmax), isometric force (F0), force-velocity profile, force redevelopment rate (ktr) at multiple sarcomere lengths, peak restretch force, viscoelastic force peak, strain-related detachment signatures and others. Feature analysis across ATP levels demonstrates systematic differences in contractile behavior, while the multi-condition dataset helps discriminate between mechanistic hypotheses.

We present a six-state strain-discretized probabilistic model. States comprise two disordered-relaxed (DRX) pools by nucleotide occupancy (ATP and ADP·Pi), two strain-sensitive attached states (pre- and post-power-stroke), and two super-relaxed (SRX) pools by nucleotide state. This explicit representation of ATP catabolism intermediates lets [ATP] and [ADP] modulate transition rates biochemically throughout the cycle. Asymmetric strain-dependent rate modifiers on the attached states capture velocity dependence of the power stroke and detachment. Force-dependent SRX mobilization encodes thick filament mechanosensing and titin passive mechanics incorporates Ca²⁺-dependent viscoelasticity [3].

The extended model qualitatively reproduces all key features of both force-velocity and slack-restretch protocols, establishing a framework capable of simulating contractile behavior across conditions. This provides the necessary foundation for mechanistic identification of ATP-dependent contractile differences — the immediate next step — and for future coupling to mitochondrial ATP production models [1] to evaluate metabolic interventions in heart failure.

[1] Collins NL, Dasika S, Van den Bergh F, Bazil J, Beard DA. Systems Analysis of Carboxylate Transport and Oxidation Pathways in Cardiac Mitochondria. bioRxiv 2026.
[2] Beard DA et al. Reduced cardiac muscle power with low ATP simulating heart failure. Biophys. J. 2022.
[3] Jezek F et al. Theoretical analysis of power-law stress relaxation and calcium-dependent passive mechanics in cardiac muscle. J. Physiol. 2026.

---

## Key Points (for reference)

### ¶1 Introduction
- ATP depletion in heart failure: energetic deficit is upstream of contractile dysfunction
- Beard et al. 2022 as parent model; established the strain-discretized cross-bridge framework
- Synergy with mitochondrial models [1]: sarcomere model as mechanical load side of ATP supply/demand axis
- Beard 2022 can reproduce FV, but fails to distinguish contractile phenotype across ATP concentrations
- Motivates the present work: extend the model to close that gap

### ¶2 Data
- Rat RV cardiac trabeculae, Baker laboratory, UCSF
- Protocols: force-velocity (FV) and slack-restretch
- ATP levels: 8 mM (saturating), 2 mM (compromised), 0.2 mM (severe, select datasets)
- Feature targets: Vmax, F0, FV profile, ktr at multiple SLs, peak restretch force, viscoelastic force peak, strain-related detachment signatures
- Feature-based fitting (not waveform fitting) — more robust and mechanistically interpretable
- SL-dependence of ktr constrains length-dependent activation mechanisms

### ¶3 Model
- Six-state strain-discretized probabilistic ODE framework
- DRX_ATP, DRX_ADP — disordered relaxed pools by nucleotide occupancy
- P1 (pre-power-stroke), P2 (post-power-stroke) — strain-sensitive attached states
- SRX_T, SRX_D — super-relaxed pools by nucleotide state
- ATP catabolism explicitly encoded: [ATP] and [ADP] modulate rates throughout cycle
- Asymmetric strain-dependent rate modifiers on attached states
- Force-dependent SRX mobilization (thick filament mechanosensing)
- Ca²⁺-dependent titin viscoelasticity [3]: systolic Ca²⁺ stiffens titin, shapes restretch transient

### ¶4 Conclusions
- Qualitative reproduction of all features across ATP conditions achieved
- No claimed quantitative advance over Beard 2022 yet — optimization in progress
- Future: couple to mitochondrial ATP supply model [1] for energetics-mechanics axis
- Clinical outlook: platform for metabolic interventions in heart failure

---

## Notes / To-Do

- [x] Six-state count: DRX_ATP, DRX_ADP, P1, P2, SRX_T, SRX_D
- [x] ATP concentrations: 8 mM, 2 mM, 0.2 mM (select)
- [x] Reference [1]: Collins NL et al. bioRxiv 2026
- [x] Filip Jezek affiliation: omit
- [ ] Confirm full author list for Beard 2022 [2]
- [ ] Final word count check before submission
