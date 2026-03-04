# Mechanism Evaluation: Literature-Guided Analysis

*Generated 2026-03-03. Based on `Model/dPUdT_CombinedTransitions.m` and `ModelOptParams/ModelParamsSlackKtrOpt.m`.*

---

## Preamble

This document provides a systematic, literature-guided evaluation of every physical mechanism implemented in the current cardiac sarcomere model. The model implements a strain-discretized Fokker–Planck-like ODE for the probability distribution of myosin cross-bridge states (Beard et al. 2022). The primary experimental targets are **force-velocity (FV)** curves and **slack-restretch (ktr)** protocols recorded from rat cardiac myocytes.

The central thesis motivating this document is that the model's failure to reproduce experimental data simultaneously for both protocols is likely caused by **incorrectly implemented or missing physics**, not merely parameter mismatch. A poorly parameterized but physically correct model can be rescued by optimization; a model missing a key mechanism cannot.

The evaluation addresses three questions for each mechanism:
1. Is the mechanism physically real and correctly captured by the ODE formulation?
2. Are the parameter values consistent with direct experimental measurements?
3. Does the mechanism materially influence the experimental observables targeted (FV curve shape, Vmax, ktr rate, force redevelopment transient, slack-phase force profile)?

All kinetic rates cited below are referenced to physiological temperature (37°C for human cardiac muscle; 25–35°C for most rat cardiac data). Temperature corrections via the Q10 rule (Q10 ~ 2–3 for actomyosin ATPase) are noted where relevant.

The **effective values** column in Part 2 accounts for the global rate multiplier `xrate = 0.435`, which is applied before each ODE integration to kah, kamh, ka, kd, k1, k_1, and k2.

---

## Part 1 — Mechanism-by-Mechanism Evaluation

### 1.1 Three-State Cross-Bridge Kinetic Scheme

**What it does.** The model implements a modified Lymn–Taylor four-step cycle collapsed into three attached states plus one free-head pool. The free (disordered relaxation) pool PD represents an ATP-bound, non-SRX head that is ready to attach. The P1 state is a weakly bound pre-power-stroke state, P2 is the strongly bound post-power-stroke (ADP-Pi released) state, and P3 is an ADP-bound "over-stroke" state. Transition rates are: PD → P1 at rate `ka × N_overlap`, P1 → P2 at `k1 × f(s)`, P2 → P3 at `k2 × f(s)`, P3 → PD at `k3`, and P1 → PD at `kd × f(s)`. The ODE propagates the full strain-resolved probability distributions dp1(s)/dt, dp2(s)/dt, dp3(s)/dt.

**Biophysical plausibility.** The Lymn–Taylor scheme (Lymn and Taylor 1971) is the mechanistic backbone of virtually all quantitative cross-bridge models. The existence of weakly and strongly bound states is established by X-ray crystallography, cryo-EM, and single-molecule studies (Rayment et al. 1993; Holmes 1997). The 3-state formulation used here (P1 weak, P2 strong, P3 ADP-bound) maps reasonably well onto the biochemically defined states: M·ATP (detached), M·ADP·Pi (pre-stroke weakly bound), M·ADP (post-stroke strongly bound), and M (rigor-like). The strain-dependent rate functions applied in the model (see Section 1.2) account for the elastic energy stored in the cross-bridge and its influence on transition rates.

**Prior model use.** The 3-state scheme (weak → strong → ADP-bound) has been used widely: Land et al. (2017) use exactly this topology for a human cardiac model at 37°C. Powers et al. (2018) use a 3-state XB in their spatially explicit titin model. Chin et al. (2006) present a 7-state extension that adds explicit ADP and ATP binding steps. Ranatunga (2010) uses a 5-state linear scheme.

**Evidence for/against.** The existence of a metastable ADP-bound state (P3 in this notation) distinct from P2 is supported by measurements showing a slow second phase of force redevelopment after rapid force steps (Piazzesi et al. 2002; Decostre et al. 2005). In cardiac myosin, ADP rebinding and release have been measured as important determinants of velocity: the ADP release rate in human cardiac myosin-II is approximately 300–400 s⁻¹ at 37°C (Malmqvist et al. 2004; Siemankowski et al. 1985), faster than the effective k2 (133 s⁻¹ in ModelParamsSlackKtrOpt, 58 s⁻¹ effective). The current k3 = 44 s⁻¹ is plausible as a rate-limiting ATP-rebinding step under saturating ATP (kcat ~ 5–30 s⁻¹ for cardiac myosin ATPase, see Section 2). However, k3 likely also captures Pi release, which is much faster (Dantzig et al. 1992: Pi release ~ 200–400 s⁻¹ in rabbit psoas). The assignment of mechanochemical steps to kinetic states in this simplified scheme should be carefully documented.

One potential problem with the current 3-state formulation is that the P3 (ADP-bound) state carries force through `kstiff3`, but by definition an ADP-bound head after the power stroke has already delivered its work. If P3 represents a post-stroke state where the head is still attached and under compressive strain from sarcomere shortening, its mechanical contribution is correctly modeled. If it represents a genuinely separate slow detachment pathway, the force assignment remains physically consistent. The ambiguity should be resolved by choosing explicit biochemical correlates for each state.

**Verdict.** ⚠️ **Questionable** — the 3-state topology is broadly correct, but the biochemical identity of P3 is ambiguous in the current implementation (it mixes ADP binding with a "slow ripping" detachment), and effective rate values are below literature ranges after the xrate correction (see Part 2).

---

### 1.2 Strain-Dependent Kinetics (Piecewise Cubic Hermite)

**What it does.** Each major rate transition (P1→P2 via `R12`, P2→P3 via `R2`, P1→PD via `R1D`, P2→P1 via `R21`) is multiplied by a strain-dependent factor computed from a piecewise cubic Hermite (PCHIP) spline evaluated at the current strain `s`. The spline is defined by five control points (`PieceWiseStrainDepX`, `PieceWiseStrainDepParams`) with extreme values of 50 at the boundaries (s = ±0.1 μm) serving as computational clamps. At zero strain, the multiplier is approximately 1.0. The multiplier rises steeply for both large positive strain (cross-bridge in tension) and large negative strain (compressed cross-bridge), with a slightly asymmetric profile.

**Biophysical plausibility.** Strain-dependent kinetics is one of the most well-established principles in cross-bridge biochemistry. The energy stored in a strained cross-bridge elastic element modifies the activation energy for each transition (Huxley 1957; Eisenberg and Hill 1978). In the thermal ratchet picture, the rate of a transition with rate constant k and transition strain s* is k × exp(−κs²/2kBT) for harmonic springs (where κ is cross-bridge stiffness). For the power stroke (P1→P2), positive strain (tension) inhibits attachment and the forward stroke; negative strain promotes it. This is consistent with the asymmetric PCHIP profiles used here, and is qualitatively supported by measurements of tension-dependent kinetics in both skeletal (Ford et al. 1977; Piazzesi et al. 1997) and cardiac muscle (Tewari, cited in sources.md).

The use of PCHIP splines rather than explicit exponential or Boltzmann functions provides more fitting flexibility but also more degrees of freedom. Each spline adds 4–5 free parameters that must be constrained from limited experimental data. The boundary multipliers of 50 (clamping detachment at extreme strain) are numerical artifacts — they prevent probability from accumulating at boundaries but have no biophysical correlate. If probabilities genuinely reach the boundary during simulation (which the `UseSpaceExtension` flag attempts to prevent), this is a sign of numerical instability rather than a correctly resolved physical process.

**Prior model use.** Strain-dependent kinetics appear in essentially all continuum cross-bridge models since Huxley (1957). Beard et al. (2022) (the parent model) use exponential strain dependence. Land et al. (2017) use a Boltzmann function for the weak→strong transition. The PCHIP parameterization in this model is a generalization that was likely introduced to obtain more flexibility in fitting asymmetric FV curves.

**Evidence for/against.** The qualitative shape of the strain-dependent multiplier is consistent with the data of Tewari (as cited): positive strain inhibits A1→A2 and A2→A3, while negative strain accelerates them. The specific functional form (PCHIP with 5 control points) is not directly measured but is used as an empirical fitting function, which is an acceptable modelling convention. The main concern is that the effective multiplier reaches 50 at the boundaries — there is no experimental measurement supporting such a large enhancement. In single-molecule studies, the ratio of forward to backward rates at extreme strains rarely exceeds 10 (Guilford et al. 1997). Boundary values of 50 are likely numerical conveniences rather than physical measurements.

**Verdict.** ⚠️ **Questionable** — the qualitative form is physically justified, but the specific PCHIP parameterization has more free parameters than can be independently constrained, and boundary multipliers of 50 are numerical artifacts rather than biophysical observations. The asymmetric profile direction is consistent with Tewari, but quantitative magnitudes need better anchoring.

---

### 1.3 Super-Relaxed State (SRX)

**What it does.** A fraction of myosin heads occupies the super-relaxed (SRX) state P_SR, which represents heads folded into the thick filament backbone in the Interacting Heads Motif (IHM) conformation (Woodhead et al. 2005; Alamo et al. 2008). These heads are biochemically inhibited (low ATPase activity) and do not contribute to force. The P_SR → PT (free ATP-bound pool) transition rate is `kmsr × exp(F_total / sigma1)`, implementing force-dependent activation. The reverse rate PT → P_SR is `ksr0 × exp(−F_total / sigma2)`. With current parameters: ksr0 = 2.4 s⁻¹ and kmsr = 1.0 s⁻¹ (base), sigma1 = 22 kPa (if F in kPa units), sigma2 = 40 kPa.

**Biophysical plausibility.** The SRX state is well-established in both skeletal and cardiac muscle (Stewart et al. 2010; McNamara et al. 2015; Hooijman et al. 2011). In the SRX, the ATP hydrolysis lifetime of a head is approximately 230 s (McNamara 2015), compared to approximately 30–100 ms in the active disordered relaxation (DRX) state. The fraction of heads in SRX at rest in cardiac muscle is debated: roughly 30–60% by mant-ATP fluorescence exchange kinetics (Hooijman et al. 2011; McNamara 2015). The force-dependent mobilization from SRX (used in Campbell 2018 to explain length-dependent activation) has experimental support from thick filament X-ray reflections that respond to sarcomere length changes and to passive stretch (Brunello et al. 2020; Irving 2017).

**Prior model use.** SRX states appear in Campbell (2018), who implements force-dependent SRX→DRX transition to reproduce Frank-Starling behavior. The present model follows that framework. Land et al. (2017) do not include SRX. Beard et al. (2022) include SRX as a major feature.

**Evidence for/against.** The mant-ATP fluorescence exchange rate for the SRX lifetime in cardiac muscle implies an effective hydrolysis cycle rate of approximately 1/230 s⁻¹ ≈ 0.004 s⁻¹ (McNamara 2015). The model's ksr0 = 2.4 s⁻¹ and kmsr = 1.0 s⁻¹ are both approximately 100-fold faster than this biochemical SRX lifetime, which means the model's SRX is not truly biochemically inhibited — it is simply a "slow pool" relative to the active cycle rates. This discrepancy may be intentional (treating SRX as a regulatory reservoir rather than a persistent off-state), but it means the modeled SRX cannot correctly reproduce the slow mant-ATP exchange kinetics. Notably, McNamara (2015) reports that cooperative interactions between adjacent myosin heads are present in skeletal but NOT cardiac thick filaments, which has implications for how SRX is recruited.

The force-dependence parameterization (sigma1 = 22 kPa, sigma2 = 40 kPa) means that at forces typical during contraction (30–100 kPa), the SRX dissolution rate increases modestly. Campbell (2018) showed this force dependence is important for explaining length-dependent activation, consistent with the implementation. However, the sigma values are not directly measured parameters — they are fitting parameters.

**Verdict.** ⚠️ **Questionable** — the SRX mechanism is biophysically real and force-dependent activation is supported by Campbell (2018) and Brunello (2020). However, the SRX exchange rates in the model (ksr0, kmsr ~1–2 s⁻¹) are 100× faster than the biochemically measured SRX lifetime (~0.004 s⁻¹), meaning the modeled SRX is a fast regulatory buffer, not a true SRX state.

---

### 1.4 Super-Relaxed ADP State (SRD)

**What it does.** The P_SRD state represents a second class of super-relaxed heads that are specifically ADP-bound while in the IHM conformation. Exchange between P_SR and P_SRD is governed by `ksr2srd = 32.1 s⁻¹` (SR→SRD) and `ksrd2sr = 2.9 s⁻¹` (SRD→SR). Exit from P_SRD to the free pool P_D occurs at `kmsrd × exp(F_SR / sigma_srd1)`. Entry from PD to P_SRD is `ksrd × exp(−F_SR / sigma_srd2)`.

**Biophysical plausibility.** The existence of an ADP-bound super-relaxed state is a theoretical extrapolation from the known biochemistry of the SRX. In the IHM structure, the blocked head (S1 inhibited by its partner) can carry various nucleotide states. Cryo-EM data show that the IHM state is stabilized by ADP binding to the blocked head: the ADP·Pi form promotes IHM release and transition toward an active conformation (Woodhead et al. 2005; Trivedi et al. 2018). The SRD formulation attempts to capture this. However, direct kinetic measurements of an ADP-specific SRX population — distinct from the general SRX — have not been reported in cardiac muscle. The rate ksr2srd = 32.1 s⁻¹ is much faster than any published SRX exchange rate and lacks direct experimental support.

The original motivation for SRD in this model appears to be reproducing the second (slow) phase of force redevelopment after slack-restretch. A slow ADP-bound pool that is mobilized during contraction could contribute to the rising phase of force. However, this is speculative, and there are alternative physical explanations (viscoelasticity, titin-based contributions).

**Prior model use.** The SRD state does not appear in canonical published cardiac XB models (Land 2017, Campbell 2018, Beard 2022). It is an extension specific to this model.

**Evidence for/against.** The concept of an ADP-stabilized IHM state is consistent with structural data (Trivedi et al. 2018; Nag et al. 2017), but the kinetic rates (ksr2srd = 32.1 s⁻¹, ksrd2sr = 2.9 s⁻¹) are orders of magnitude faster than measured SRX exchange rates. The ADP dissociation rate from thick filament heads in SRX has not been directly measured. Without direct kinetic evidence, the SRD state functions as a second adjustable time constant for the slow phase of contraction dynamics, which is a modelling convenience rather than a physically anchored mechanism.

**Verdict.** ❌ **Likely missing physical basis** — the concept is structurally motivated, but the kinetic rates are unconstrained by experiment, the state does not appear in peer-reviewed cardiac models, and its primary role in the current model appears to be providing an additional slow time constant. This mechanism should be critically scrutinized and either properly anchored to measured ADP-SRX kinetics or removed.

---

### 1.5 Filament Overlap

**What it does.** When `UseOverlap = 1`, the attachment rate to P1 is scaled by a dimensionless overlap factor `N_overlap = (L_ov × 2) / (L_thick − L_hbare)`, where `L_ov = min(L_thick/2, SL/2) − max((SL−L_thin), L_hbare/2)` represents the length of the single-overlap zone per half-sarcomere. Parameters: `L_thick = 1.14 μm` (thick filament half-length), `L_hbare = 0.10 μm` (bare zone half-length), `L_thin = 1.377 μm` (thin filament length). The overlap is linear in the computed overlap length.

**Biophysical plausibility.** Filament overlap is one of the most firmly established determinants of sarcomere force (Gordon et al. 1966). The active force-length relation is proportional to the fraction of myosin heads in the overlap zone, and this geometric relationship is a cornerstone of sarcomere mechanics. The linear overlap → force relationship is correct for constant cross-bridge density along both filaments.

**Prior model use.** Overlap is included in virtually all spatially explicit sarcomere models (Powers 2018; Prostsenko 2020; Gordon et al. 1966). Land (2017) uses a sigmoidal approximation rather than explicit geometry.

**Evidence for/against — filament lengths.** This is where a significant problem arises.

- **Thin filament length L_thin = 1.377 μm** — the current value exceeds SL0/2 = 1.1 μm, meaning that at the operating sarcomere length (SL0 = 2.2 μm), the thin filaments extend past the center of the sarcomere and would geometrically overlap with each other from both half-sarcomeres. Published thin filament lengths in rat cardiac muscle are 1.10–1.16 μm (Tonino et al. 2017; Witt et al. 2006) or 1.0–1.2 μm by electron microscopy (Robinson and Winegrad 1979; Akella et al. 1997). Human cardiac thin filament length is approximately 1.16 μm (Tonino et al. 2017). A value of 1.377 μm is incompatible with reported electron microscopy data and would imply double-overlap at SL = 2.2 μm, which is geometrically incorrect.

- **Thick filament half-length L_thick = 1.14 μm** — the total thick filament length would thus be 2.28 μm, which is within published range (1.6–1.7 μm total, so 0.8–0.85 μm half-length by conventional measurement; Tregear and Squire 1973; Luther et al. 2008). Wait — the code uses L_thick as a half-length (L_T_HS1 = min(L_thick × 0.5, SL × 0.5)), so the total thick filament is 2 × 1.14 = 2.28 μm. Published values for cardiac thick filament total length are approximately 1.65 μm (0.825 μm half), not 2.28 μm. At L_thick = 1.14 μm (half), the total = 2.28 μm, which is far above any reported measurement. This strongly suggests the overlap geometry is parameterized incorrectly.

To reconcile: if the code interprets L_thick as the full thick filament length (and divides by 2 internally), then L_thick = 1.14 μm total is far too short (half of a normal thick filament). The code line `L_T_HS1 = min(L_thick*0.5, SL*0.5)` suggests L_thick is used as total length and halved for the half-sarcomere. In that case, 1.14 μm total is much shorter than literature values (~1.65 μm), and L_thin = 1.377 μm is similarly interpreted as full thin filament length (halved to 0.689 μm per half-sarcomere). Published thin filament length of ~1.15 μm total would then halve to 0.575 μm per half-sarcomere, which is close to 0.689 μm within uncertainty. The thick filament discrepancy remains severe.

In either interpretation, at least one of the two filament lengths is outside the published experimental range, and the resulting overlap function will be incorrect at physiological SL values.

**Verdict.** ⚠️ **Questionable** — filament overlap is a necessary and well-supported mechanism, but the specific geometric parameters (L_thick and L_thin) appear to be inconsistent with published electron microscopy values, and the geometric interpretation of the overlap formula needs clarification. The overlap function likely operates in a geometrically incorrect regime, potentially explaining why the model needs large parameter compensations elsewhere.

---

### 1.6 Passive Force / Titin

**What it does.** When `UsePassive = 1` and `UseTitinIdentifiedPassive = 0` (current active configuration), passive force is: `F_passive = k_pas × max(SL − LSE − Lsc0, 0)^gamma`, where `Lsc0 = 1.51 μm`, `k_pas = 77.3`, `gamma = 2.67`. When `UseTitinIdentifiedPassive = 1` (currently OFF), a separately identified function is used: `F_passive = 0.4 × (SL/2 − x0)^7.9` with a high-power term for very short lengths.

**Biophysical plausibility.** Titin is the primary source of passive sarcomere stiffness in the physiological SL range (Granzier and Labeit 2004; Linke 2018). The titin spring behavior in the range 2.0–2.4 μm SL is dominated by the extensible I-band regions (N2B and PEVK in cardiac isoforms). The passive force-length relationship in this range is commonly described as an exponential or power law, with published exponent values for cardiac muscle of approximately 7–10 when fit as `F = k (SL − SL0)^n` (Granzier and Irving 1995; Wu et al. 2000; Terui et al. 2008). The value of gamma = 2.67 is significantly lower than this range.

The `UseTitinIdentifiedPassive` pathway uses gamma = 7.9, which is within the published range. The fact that this flag is currently OFF means the model uses an underpowered passive stiffness (gamma = 2.67 vs. literature ~7.9), which will produce too-gradual an increase in passive force with extension and too-little passive force at SL > ~2.3 μm. This may affect the slack-restretch protocol, where the sarcomere is shortened and then stretched back, as the passive restoring force contributes to the force transient.

Additional passive force contributions mentioned in Freundt (2019) and Herzog (2022) include titin oxidation, phosphorylation (Hamdani 2017), and Ca²⁺-dependent stiffness increase (~20% with Ca²⁺). The current model captures none of these secondary effects.

**Prior model use.** Powers (2018) use a nonlinear titin spring in a spatially explicit model. Prostsenko (2020) report passive force-length matching the current model's general range. Land (2017) do not include passive titin explicitly (they model cell-level data where passive stiffness is built into boundary conditions).

**Evidence for/against.** Wu et al. (2000) measured cardiac titin passive stiffness with gamma ≈ 8–10 in rat ventricle. Terui et al. (2008) found similar values for rabbit cardiac muscle. The value of gamma = 2.67 would predict much more gradual passive stiffening and lower passive force at physiological SL (2.2–2.4 μm). The identified passive law (gamma = 7.9) stored in the code as the `UseTitinIdentifiedPassive` pathway is more consistent with the literature, yet it is disabled. The Lsc0 = 1.51 μm slack length is approximately correct for rat cardiac (literature: 1.5–1.85 μm, Granzier and Irving 1995). The high-power term in the identified passive law (`0.5e9 × (x − 0.95)^13`) appears to be a numerical stiffening wall to prevent sarcomere collapse, which is a legitimate numerical convention.

**Verdict.** ❌ **Incorrectly parameterized** — passive force exponent gamma = 2.67 is far below the measured range (~7.9) for cardiac titin. The `UseTitinIdentifiedPassive` pathway with gamma = 7.9 is more physically correct but is currently OFF. This should be a high-priority fix. Ca²⁺-dependent titin stiffening (Herzog 2022) is also not accounted for.

---

### 1.7 Serial Stiffness

**What it does.** The experimental measurement apparatus introduces a serial elastic element (SE) between the muscle and the force transducer or length controller. This is modeled by adding a spring of stiffness kSE = 3203 (units: force/length, kPa/μm or equivalent) with a power-law extension `Force = sign(LSE) × |kSE × LSE|^ekSE`, where `ekSE = 0.925 ≈ 1`. The SE length LSE is an additional state variable evolving according to `dLSEdt = vel − velHS`, where velHS is the half-sarcomere velocity computed from force balance.

**Biophysical plausibility.** Serial stiffness in cardiac muscle experiments has two sources: (a) the experimental rig (glass fibers, carbon fibers, transducer compliance) and (b) myofilament lattice compliance (Huxley et al. 1994; Irving et al. 1992). Tombe (2016) emphasizes that the rig compliance causes SL to change by up to 20% even when the motor position is held fixed, which fundamentally alters the force-velocity relationship measured at the muscle level vs. the sarcomere level. Tombe (1991) modeled an internal viscous element that limits unloaded shortening velocity in isolated cardiac myocytes, which this serial element partially captures.

**Prior model use.** Land (2017) use a dashpot (viscous element) in series with the contractile element, combined with a linear spring, and report that this serial compliance substantially affects apparent Vmax. The present model follows the same rationale. Tombe (2016) provides experimental evidence that this correction is essential for interpreting skinned fiber data.

**Evidence for/against.** The nonlinearity (ekSE = 0.925) is a minor correction and biophysically plausible for a nonlinear fiber stiffness. The magnitude kSE = 3203 is difficult to compare directly to literature without knowing the force units (kPa) and the spatial units (μm per half-sarcomere). If force is in kPa and SL in μm, then kSE × LSE gives kPa when LSE is in μm. The rig stiffness of cardiac fiber rigs is typically 0.5–5 kPa/μm (Irving et al. 1992), suggesting kSE = 3203 might be specified in units inconsistent with other parameters, or the SE is much stiffer than a typical experimental rig. This needs clarification.

The `dLSEdt` term implements the idea that the SE acts as a velocity-dependent partitioner: at slow shortening, most length change occurs in the XB population; at fast shortening (slack), the SE can lengthen rapidly. This is the correct physics. The force-balance equation `velHS = (Force − F_total)/mu` is essentially Kelvin–Voigt for the viscous element. The sign convention and force balance should be checked against the convention that velHS is the sarcomere shortening velocity.

**Verdict.** ✅ **Well-supported** — serial compliance is real and well-documented (Tombe 2016; Land 2017). The implementation is physically sound. Parameter values should be confirmed against the specific experimental rig used by the Baker lab. The nonlinearity (ekSE) is minor.

---

### 1.8 Negative Stiffness Asymmetry

**What it does.** When `UseNegativeKstiff = 1`, the active force is split by strain sign. For P2 heads: `F_active = kstiff2 × p2_1_pos + kstiff2_n × p2_1_neg`, where `p2_1_pos = ∫(s > −dr) s f(s) p2(s) ds` and `p2_1_neg = ∫(s < −dr) s f(s) p2(s) ds`. Similarly for P1. With `kstiff1 = 16000`, `kstiff1_n = 160` (ratio 100×), and `kstiff2 = 97600`, `kstiff2_n = 704` (ratio ~139×). This means heads at negative strain (compressed) contribute to active force at approximately 1/100 to 1/140 the stiffness of heads at positive strain.

**Biophysical plausibility.** The elastic energy stored in a cross-bridge under positive strain (the head pulled in the direction of force generation) is different from the energy stored under negative strain (the head being pushed backward). In the standard Huxley (1957) model, the cross-bridge spring is treated as purely tensile: cross-bridges under compressive strain (negative s) are assumed to rapidly detach or contribute negligibly to force. This motivates a lower stiffness for negative-strain heads.

However, the experimental basis for a factor-of-100 asymmetry is unclear. Single-molecule measurements of myosin stiffness under negative vs. positive load show some asymmetry: myosin-II heads under compressive load do exhibit lower effective stiffness and faster detachment rates (Veigel et al. 2003; Mehta et al. 1999; Kaya and Higuchi 2010). The detachment rate increases strongly with negative strain (Guilford et al. 1997), which means the negative-strain population is dynamically depleted rather than carrying negative force. The conventional treatment is to make kd(s) large for s < 0 (strain-dependent detachment) rather than to reduce kstiff. The two approaches have different predictions: the current implementation allows P1 and P2 heads at negative strain to exist (probability is not zero) but contributes little force, whereas a fast detachment at negative strain would actively remove them from the distribution.

The 100× ratio has no direct experimental measurement support. It appears to be an empirical correction used to limit the "negative force" contribution of compressed heads and thus improve the FV curve shape. This is a modelling convention rather than a biophysical measurement.

**Prior model use.** Huxley (1957) set f(s) = 0 for s < 0 and g(s) = g₂ >> g₁ for s < 0, effectively implementing fast detachment under negative strain. Land (2017) use the same strain-independent detachment but with a velocity-dependent strain reset. Beard et al. (2022) do not appear to use asymmetric stiffness. The asymmetric stiffness approach is unusual in the literature.

**Evidence against.** The 100× stiffness ratio is not supported by any cited biophysical measurement. Single-molecule data suggest myosin stiffness under negative load (compression) is comparable to positive load (Mehta et al. 1999; Kaya and Higuchi 2010). The correct physical mechanism for limiting negative-strain force contributions is fast detachment (high kd at s < 0), not a reduced stiffness constant. Implementing negative stiffness asymmetry as a stiffness rescaling is an unphysical shortcut that will predict incorrect force transients during rapid length changes.

**Verdict.** ❌ **Likely missing or should be revised** — the physical motivation (limiting negative-strain force) is correct, but the implementation (reducing kstiff by 100×) is not supported by experimental measurements and is unconventional. The standard approach is to increase the detachment rate kd at negative strain, which is physically equivalent to rapid removal of negative-strain heads from the population. The high negative-strain detachment is already partially encoded in the PCHIP multiplier for R1D, suggesting redundancy.

---

### 1.9 A2 Attachment Shift (Actin Site Hopping)

**What it does.** When `UseA2AttachmentShift = 1`, strongly bound P2 heads at high positive or negative strain hop to the next actin binding site, a distance `d_actin = 5.5 nm` (= 0.0055 μm) closer to zero strain. The hopping rate is `RA_K = slope/dr × max(0, s − s_threshold_R) + slope/dr × max(0, −s − s_threshold_L)`, with `slope = 16496 s⁻¹/μm`, `d_actin = 5.5 nm`, `s_threshold_R = 4.6 nm`, `s_threshold_L = 5500 nm` (effectively infinity, so only positive-strain hopping is active). The mass leaving bin s is redistributed to bin (s − d_actin) for positive strain.

**Biophysical plausibility.** The actin periodicity (half-repeat) in cardiac thin filaments is approximately 5.46 nm (Tregear and Squire 1973; Egelman and DeRosier 1992), corresponding to the half-helical repeat of F-actin. The concept of a cross-bridge hopping from one actin monomer to the next has been discussed in theoretical contexts as a way to explain force-velocity data. Specifically, during rapid shortening, a strongly bound head at high positive strain may detach and reattach to the next actin site, releasing stored elastic energy and resetting its strain. This is sometimes called "strain-dependent rebinding" or "actin hopping."

Direct experimental evidence for in-cycle actin site reassignment is limited. Some structural studies (Dominguez et al. 1998) suggest myosin can contact different actin monomers depending on the power stroke state, and the cross-bridge cycle is thought to involve a sequential actin binding. However, a kinetic model where the rate of site-to-site transfer is proportional to strain above a threshold is an approximation with no direct single-molecule measurement to anchor the slope parameter (16496 s⁻¹/μm is a large number — at s = 10 nm above threshold, the rate is 0.01 μm × 16496/0.01 = 16496 s⁻¹, comparable to ATP hydrolysis rates).

**Prior model use.** Actin site hopping does not appear in standard published cardiac XB models (Land 2017, Campbell 2018, Beard 2022). It is a non-standard extension. Similar concepts have appeared in skeletal muscle models that allow multiple binding site access (Duke 1999; Lan and Sun 2005), but not as a implemented PCHIP-like redistribution at the rate constants used here.

**Evidence for/against.** The actin repeat distance (5.5 nm) is physically correct (Egelman and DeRosier 1992; Holmes et al. 1990). The functional motivation is plausible: strain-induced redistribution could act as a mechanism to enhance unloaded shortening velocity by converting high-strain P2 heads to near-zero-strain states, from which detachment is faster. However, the rate constant slope = 16496 s⁻¹/μm has no experimental measurement to support it, and the `s_threshold_L = 5500 nm` effectively disables negative-strain hopping. The current implementation has the feel of a tunable free parameter rather than a constrained biophysical mechanism.

**Verdict.** ⚠️ **Questionable** — the geometric rationale (5.5 nm actin repeat) is correct and the concept is physically motivated, but the kinetic rate (slope = 16496 s⁻¹/μm) lacks experimental support. The mechanism is not present in comparable published models. It should be treated as a hypothesis that requires additional justification from structural or mechanical data, and its role in the model should be carefully tested by disabling it and observing the effect on FV curves.

---

### 1.10 Viscous Drag

**What it does.** Internal viscous drag is modeled as a force opposing sarcomere shortening velocity: the half-sarcomere velocity is computed from `velHS = (Force − F_total) / mu` when `Force > F_total`, and `velHS = (Force − F_total) / mu_neg` when `Force < F_total` (i.e., when the sarcomere is being extended or pulled). With `mu = mu_neg = 0.633` (symmetric). A velocity-squared term `mu2 ≈ 0` is effectively disabled. The viscous element acts in series with the serial stiffness element.

**Biophysical plausibility.** Internal viscosity in cardiac myocytes is experimentally supported. Tombe (1991) directly measured an internal viscous element by analyzing the force-velocity relationship at very low loads and high velocities in isolated rat myocardium at 25°C. They found that an internal viscous term (linear in velocity) is required to explain the deviation of the FV curve from a simple Hill hyperbola at high shortening velocities. The physical origin of this viscosity may include (a) cross-bridge head movement against the aqueous medium, (b) myosin rod-rod interactions in the thick filament, (c) actin filament bending stiffness and hydrodynamic drag, and (d) viscoelastic contributions from titin.

The symmetric viscosity (mu = mu_neg) is a simplification. In reality, the viscosity during shortening and lengthening may differ due to different cross-bridge population dynamics (e.g., catch bonds, yield force asymmetry). The mu2 (velocity-squared) term is disabled, limiting the model to purely linear (Newtonian) drag.

**Prior model use.** Land (2017) use a dashpot (viscous element) explicitly and report it substantially affects apparent Vmax. Tombe (1991) provide the experimental data motivating this. The Kelvin–Voigt type coupling of viscous and elastic elements in series is standard.

**Evidence for/against.** Tombe (1991) measured viscous drag in rat myocardium at 25°C and found a linear viscous coefficient. Qualitatively, mu = 0.633 in the model units is consistent with a viscous element that limits Vmax to ~10 μm/s (as set by the `vmax` parameter). The symmetric mu = mu_neg is likely an approximation — during the slack phase (rapid shortening to zero load), the effective viscous force would be high, which is physically correct. The absence of a velocity-squared term (mu2 ≈ 0) means nonlinear drag at very high velocities is not modeled, which may contribute to inaccurate slack protocol force transients.

**Verdict.** ✅ **Well-supported** — internal viscous drag is experimentally measured (Tombe 1991) and correctly motivates the inclusion of a viscous element. The symmetric implementation is an approximation. The absence of velocity-squared drag is unlikely to cause major discrepancies at the velocities studied.

---

### 1.11 Calcium / Thin Filament Activation

**What it does.** Calcium-regulated thin filament activation is implemented through cooperative troponin-C binding. With `Ca = 1000` (effectively saturating), the thin filament is always fully ON in the current FV and ktr protocols. The cooperative binding is parameterized by `K_coop = 5.7`, `k_on = 59.7 μM⁻¹s⁻¹`, `k_off = 307 s⁻¹`. The `UseCa = 0` flag means that the troponin activation model is not actively solved as a time-varying ODE — Ca is held at a constant value, locking troponin C in a nearly saturated state.

**Biophysical plausibility.** Cooperative troponin activation is well-established (Brenner 1988; Rice et al. 2003; Dobesh et al. 2002). At saturating Ca²⁺, the thin filament should be in the "ON" state with high regulatory protein availability. The parameters K_coop = 5.7, k_on = 59.7 μM⁻¹s⁻¹, k_off = 307 s⁻¹ are within the range reported for cardiac troponin C in skinned fiber studies (Robertson et al. 1981; Morimoto 1991). The cooperativity coefficient K_coop = 5.7 is consistent with literature values of 3–8 (Dobesh et al. 2002).

Since Ca is saturating in the FV and ktr protocols, the detailed Ca-troponin dynamics do not substantially affect results. The thin filament effect is simply that all heads are permissible for attachment, which is the correct boundary condition for maximal activation experiments.

**Prior model use.** Rice et al. (2003) developed a detailed cooperative thin filament activation model. Land (2017) use a Hill-type activation model. The current model's Ca-troponin scheme (when active) follows the general cooperative binding framework.

**Evidence for/against.** With Ca = 1000 and UseCa = 0, the model effectively bypasses thin filament dynamics. This is appropriate for force-velocity and Ktr measurements at saturating Ca²⁺ (pCa 4.5–5.0, typical for maximal activation in rat cardiac skinned fibers). The k_on and k_off values yield an apparent dissociation constant of k_off/k_on = 307/59.7 ≈ 5.1 μM, which is within the physiological range for cardiac TnC (Robertson et al. 1981).

**Verdict.** ✅ **Well-supported** — for the current saturating-Ca protocols, Ca = 1000 and bypassing troponin dynamics is appropriate. The troponin parameters are within published ranges. Ca dynamics would become relevant for submaximal Ca protocols, stairs experiments, or physiological Ca transients.

---

### 1.12 Proposed New Mechanisms

The following mechanisms are absent from the current model but are supported by literature and are likely to contribute to the discrepancies between model predictions and FV/ktr data.

#### 1.12.1 Cross-Bridge Distortion Reset at State Transitions

**Physical motivation.** During shortening, strongly bound P2 heads accumulate at negative strain (compressed state) as the sarcomere slides. These heads resist shortening (negative force generation). In the Land (2017) mean-field ODE model, transitions between cross-bridge states reset the mean distortion variable — newly entering cross-bridges are assumed to bind with zero distortion, and the mean distortion of the destination state is updated via a weighted average. This prevents accumulation of negative-strain strongly bound heads.

**Important clarification.** Land et al. (2017) present the mean-distortion-variable formulation — and the associated distortion reset — as a **computational simplification**, not as a physiologically motivated mechanism. The approach replaces the full strain distribution (PDE or strain-discretized ODE system) with a small set of scalar ODEs tracking the mean strain per state, making the model tractable for embedding in whole-organ finite-element simulations. Land explicitly notes that this simplification has **minor impact on model results** within the context of his mean-field framework. The distortion "reset" is a bookkeeping consequence of the mean-field approximation (incoming cross-bridges enter with zero distortion, diluting the mean), not a separate biophysical process.

**Applicability to the current strain-discretized model.** The current model tracks the full probability distribution p(s) over a strain grid, not mean distortion variables. Implementing a distortion "reset" in this framework would require redistributing probability mass from one strain bin to another during state transitions (e.g., depositing P2 probability at strain `s + dr` rather than at the same s when the P1→P2 transition occurs). This redistribution would introduce **numerical smearing** (artificial diffusion) in the strain distribution: the probability deposited at `s + dr` would not, in general, align exactly with a grid point, requiring interpolation and spreading probability across adjacent bins. Over many time steps, this numerical diffusion would broaden the strain distribution beyond its physical width, corrupting the force calculation and strain-dependent kinetics. The mean-field approach of Land avoids this problem because it tracks only a scalar (the mean), which has no grid resolution to smear.

A physically motivated alternative that avoids numerical smearing would be to increase the strain-dependent detachment rate at negative strain (i.e., the classical Huxley approach of high g(s) for s < 0), which actively removes negative-strain heads from the distribution rather than repositioning them.

**Relevant protocols.** FV curve (especially curvature and Vmax), ktr rise rate.

**Literature support.** Land et al. (2017) — computational simplification for organ-scale models, minor impact on results. The concept of zero-distortion binding is consistent with single-molecule observations that the power stroke occurs at the moment of strong binding (Capitanio et al. 2006; Guilford et al. 1997). Rice et al. (2008) use a similar mean-field ODE approach for cooperative activation in cardiac muscle.

**Suggested implementation.** Rather than implementing a strain redistribution (which would cause numerical smearing), the preferred approach is to ensure that the strain-dependent detachment rates at negative strain are sufficiently large to prevent accumulation of compressed P2 heads. This is already partially captured by the PCHIP multipliers for R1D, but the magnitude may need to be increased. See also Section 1.8 (negative stiffness asymmetry) for the related discussion.

#### 1.12.2 Nonlinear (Power-Law) Passive Stiffness with Correct Exponent

**Physical motivation.** As discussed in Section 1.6, the current gamma = 2.67 underpredicts cardiac passive stiffness, particularly at the longer SLs reached during diastole (SL > 2.3 μm). The correct exponent (gamma ≈ 7.9, from the `UseTitinIdentifiedPassive` pathway) should be activated. Passive force contributes to total force during ktr measurement (the muscle is at a finite SL, and titin contributes passive tension), and to the restretch phase of the slack protocol.

**Relevant protocols.** Slack-restretch force transient (restretch to SL0 from shortened state), ktr measured at 2.2 μm vs. 2.0 μm SL.

**Literature support.** Granzier and Irving (1995), Wu et al. (2000), Terui et al. (2008): cardiac titin passive stiffness exponent ~ 7–10 in the physiological SL range.

**Suggested implementation.** Activate `UseTitinIdentifiedPassive = 1`. Re-run optimization to allow k_pas to re-fit.

#### 1.12.3 Ca²⁺-Dependent Titin Stiffening

**Physical motivation.** The stiffness of cardiac titin increases by approximately 20% in the presence of physiological Ca²⁺ (Herzog 2022). This is attributed to binding of PEVK titin to actin, which effectively shortens the free spring length. In maximal Ca²⁺ protocols (Ca = 1000), this would increase passive stiffness by ~20% compared to the passive (resting) measurement.

**Relevant protocols.** Force-length, ktr, slack force transient at high SL.

**Literature support.** Herzog (2022), Freundt (2019), Hamdani (2017).

**Suggested implementation.** Add a Ca-dependent multiplicative factor to k_pas: `k_pas_eff = k_pas × (1 + 0.2 × Ca / (Ca + K_Ca_titin))`. For saturating Ca, k_pas_eff ≈ 1.2 × k_pas. This is a one-parameter addition.

#### 1.12.4 Catch Bond / Cooperative XB Detachment During Lengthening

**Physical motivation.** During sarcomere stretch (as in the restretch phase of the slack-restretch protocol), cross-bridges are under positive strain and detach more slowly. This "catch bond" behavior (Guo and Bhattacharya 2012; Nishizaka et al. 1995) means that force rises faster during restretch than predicted by a simple model where all rate constants are symmetric.

**Cooperative XB detachment mechanism.** Cross-bridge cooperativity refers to the phenomenon where binding or unbinding of one cross-bridge influences the kinetics of its neighbors. Four categories of cooperativity have been identified (Razumova et al. 1999, 2000):

1. **RU → XB:** A regulatory unit (troponin–tropomyosin) in the "on" state promotes cross-bridge attachment by exposing actin binding sites.
2. **XB → RU:** A bound cross-bridge holds the regulatory unit in the "on" conformation, increasing Ca²⁺ affinity of troponin C and keeping the thin filament locally activated. This creates a positive feedback loop: once one head binds, nearby heads can bind more easily.
3. **RU → RU:** Adjacent tropomyosin units on the thin filament are structurally coupled — one shifting to the "on" position drags its neighbor toward "on." Each tropomyosin covers ~7 actin monomers (~38 nm), and the end-to-end coupling means activation spreads cooperatively along the thin filament.
4. **XB → XB:** A bound cross-bridge directly enhances the attachment rate (or reduces the detachment rate) of its nearest-neighbor cross-bridge. The mechanism is mechanical: bound heads stiffen the local myofilament lattice, reducing compliance and improving the geometric alignment of neighboring myosin heads with actin sites. Conversely, detachment of one cross-bridge increases local compliance, straining neighboring heads and accelerating their detachment — a cooperative detachment cascade.

The cooperative detachment aspect (XB → XB) is particularly relevant during the restretch phase. When the sarcomere is rapidly restretched, the population of bound cross-bridges is simultaneously strained. If the strain on any head exceeds a critical value, it detaches. This detachment increases the strain on neighboring heads (because the load is now shared among fewer heads), potentially triggering a cascade of detachments. Conversely, during isometric contraction or slow stretch, groups of bound cross-bridges mutually stabilize each other against detachment (cooperative recruitment). This cooperativity contributes to the steep force–calcium relationship observed experimentally (Hill coefficient ~5–7 in cardiac muscle), which cannot be reproduced by simple mass-action kinetics (predicted Hill coefficient ~1–2). It also affects ktr and relaxation kinetics.

The current model does not include explicit XB–XB cooperativity. All cross-bridges cycle independently, interacting only through the global force balance. Adding nearest-neighbor interactions would require either (a) an Ising-type model on the thin filament regulatory units (Campbell 2003; Campbell et al. 2010), or (b) a spatially explicit model where each myosin head has a defined position along the filament (Daniel et al. 1998).

**Relevant protocols.** Slack-restretch force rise slope, early overshoot in force transient, force–calcium relationship (Hill coefficient), ktr.

**Literature support.** Razumova et al. (1999, 2000) — stiffness-distortion sarcomere model with cooperative XB binding and systematic analysis of four cooperativity types. Campbell (2003) — Ising model of cardiac thin filament activation with nearest-neighbor cooperative interactions. Campbell et al. (2010) — coupling of adjacent tropomyosins enhances cross-bridge-mediated cooperative activation. Brenner (1988), Tewari (as cited).

**Suggested implementation.** Two approaches: (a) **Minimal:** Increase the positive-strain asymmetry of kd(s) — reduce detachment rate exponentially for s > 0 (catch bond term in R1D). The current PCHIP multiplier for R1D already does this to some extent; the issue is whether the magnitude is sufficient. (b) **Full cooperativity:** Implement a nearest-neighbor coupling on the thin filament regulatory units following the Ising model of Campbell (2003). This adds significant computational complexity (Markov chain over regulatory unit states) but would correctly capture the steep force–calcium relationship and cooperative detachment cascades.

#### 1.12.5 Maxwell Viscoelastic Dashpot (already in code, currently OFF)

**Physical motivation.** Titin and other structural elements contribute viscoelastic (not purely elastic) behavior: upon rapid shortening, titin force drops (creep), and upon restretch, it transiently overshoots. A Maxwell element (spring and dashpot in series) correctly captures this time-dependent behavior. The code already implements this (`UseMaxwellDashpot` flag) with `kSE_M = 100` and `eta_M = 0.1`, but it is currently OFF.

**Relevant protocols.** Slack-restretch force transient shape, particularly the slow rising and falling phases.

**Literature support.** Hamdani (2017), Freundt (2019), Linke (2018): titin behaves as a viscoelastic spring, not a purely elastic spring.

**Suggested implementation.** Activate `UseMaxwellDashpot = 1`. Tune kSE_M and eta_M by fitting the passive force transient during slack.

#### 1.12.6 Lattice Spacing Effects on Cross-Bridge Kinetics

**Physical motivation.** The myofilament lattice is a three-dimensional hexagonal array where thick and thin filaments are separated by a radial distance (the "lattice spacing," typically ~20–25 nm surface-to-surface at rest). This spacing changes with sarcomere length due to the constant-volume constraint of muscle fibers: at longer sarcomere lengths, the lattice compresses radially, bringing thick and thin filaments closer together. The radial distance that a myosin head must bridge to reach actin directly modulates:

- **Attachment rate:** As lattice spacing increases, the myosin head must extend further radially, reducing the probability of reaching an actin binding site and slowing attachment. The axial offset window in which the powerstroke can occur narrows by more than 5 nm over the physiological range of lattice spacings.
- **Detachment rate / attachment duration:** Decreased thick-to-thin filament surface distance prolongs cross-bridge attachment time by reducing the strain needed to reach the detachment threshold.
- **Force direction:** Cross-bridges exert both axial (along the filament) and radial force components. At typical lattice spacings, the radial component can be 30–50% of the axial component, depending on cross-bridge geometry and lever arm angle.
- **Length-dependent activation:** Changes in lattice spacing with sarcomere length contribute to the Frank–Starling mechanism: at longer SL, reduced lattice spacing brings filaments closer, enhancing cross-bridge formation — independent of, and complementary to, the overlap mechanism.

The current model is strictly one-dimensional (axial strain only) and ignores lattice spacing entirely. All cross-bridges are assumed to be at zero radial distance from actin. This simplification may contribute to the model's difficulty in simultaneously fitting FV and ktr data at different sarcomere lengths.

**Relevant protocols.** Length-tension relationship, ktr dependence on SL, FV curves at different SL.

**Literature support.** Tanner et al. (2007) — spatially explicit 3D model showing axial and radial forces depend on lattice spacing (*PLoS Comput Biol* 3: e115). Williams et al. (2012) — elastic energy storage and radial forces depend on sarcomere length (*Biophys J* 103: 677–686). Williams et al. (2013) — the length-tension curve in muscle depends on lattice spacing (*Proc R Soc B* 280: 20130697). Irving (2017) — thick filament strain and inter-filament spacing changes with SL.

**Suggested implementation.** Add a lattice-spacing-dependent scaling factor to the attachment rate ka: `ka_eff = ka × f(d10(SL))`, where `d10` is the (1,0) lattice spacing computed from the constant-volume relationship `d10(SL) = d10_ref × sqrt(SL_ref / SL)`, and `f` is a decreasing function (e.g., Gaussian or sigmoid) of the radial distance. This is a low-cost addition (one multiplicative correction) that captures the first-order lattice spacing effect without requiring a full 3D model.

#### 1.12.7 Rotational / Angular Myosin Spring Model

**Physical motivation.** Traditional cross-bridge models (including the current one) treat the myosin head as a single linear (Hookean) spring oriented parallel to the filament axis. This is a substantial simplification of the actual myosin structure, which consists of a globular motor domain, a converter domain (hinge), and a lever arm (light-chain domain). Multi-spring models replace the single linear spring with:

- A **torsional (angular) spring** at the converter domain, representing the hinge that rotates during the powerstroke. The powerstroke is modeled as torque applied at the converter, rotating the lever arm through ~70° (Rayment et al. 1993).
- An **extensional spring** along the S2 segment connecting the head to the thick filament backbone.
- In the most detailed formulation (the 4sXB model of Tanner et al. 2012), four springs correspond to: (1) the S2 link, (2) the light-chain domain/lever arm, (3) the converter domain (torsional), and (4) the catalytic domain connection to actin.

The force-extension relationship of a rotating lever arm is inherently **nonlinear and asymmetric**: the same angular displacement produces different axial forces depending on the current lever arm angle. This naturally produces both axial and radial force components and predicts different effective stiffnesses in pre- vs. post-powerstroke states — without needing an ad hoc stiffness asymmetry parameter (cf. Section 1.8).

**Relevance to current model.** The current model uses separate stiffness values for P1 (kstiff1 = 16000) and P2 (kstiff2 = 97600), with an empirical ratio of ~6×. A rotational spring model would derive this ratio from the lever arm geometry and the change in angle between states, providing a physically constrained relationship between pre- and post-stroke stiffness rather than fitting them independently. It would also naturally explain the asymmetric stiffness at negative vs. positive strain (Section 1.8), potentially eliminating the need for the kstiff_n parameters.

**Relevant protocols.** FV curvature, quick force recovery after length steps, stiffness measurements.

**Literature support.** Huxley and Simmons (1971) — original rotating cross-bridge / tilting head concept (*Nature* 233: 533–538). Daniel et al. (1998) — spatially explicit model with compliant filaments and multi-spring cross-bridges (*Biophys J* 74: 1611–1621). Tanner et al. (2012) — defines the 2sXB and 4sXB multi-spring models, showing that axial and radial forces depend on lattice spacing (*PLoS Comput Biol* 8: e1001018).

**Suggested implementation.** Replace the linear spring force `F = kstiff × s` with a geometry-derived nonlinear force `F(s) = (k_angular / L_lever) × sin(theta_0 + s/L_lever)`, where `L_lever` is the lever arm length (~8.5 nm), `theta_0` is the resting angle, and `k_angular` is the torsional stiffness. This replaces kstiff1, kstiff2, kstiff1_n, kstiff2_n with two physical parameters (k_angular, L_lever) and derives the state-dependent and strain-sign-dependent stiffness from geometry.

#### 1.12.8 Spatial Distribution of Myosin Crowns

**Physical motivation.** In vertebrate striated muscle, myosin heads are not uniformly distributed along the thick filament. They are arranged in **crowns** — groups of 3 pairs of heads at each axial level — spaced **14.3 nm** apart along the thick filament backbone. Three successive crowns form a **~43 nm repeat**, with each crown rotated ~40° relative to the previous one, forming a 3-start helix (Squire 2017). Key consequences:

- **Vernier mismatch.** The actin binding sites repeat every ~5.5 nm (half-helical repeat of F-actin, ~36 nm for the full helical repeat). The mismatch between the 14.3 nm myosin crown spacing and the 5.5 nm actin site spacing means that at any given moment, only a fraction of heads are in a favorable axial position to bind. The probability of a head being within binding range depends on the strain window of the attachment rate function — with a typical binding window of ±5 nm, roughly 50–70% of heads are geometrically competent to bind.
- **Rotational mismatch.** The 40° rotation per crown means that within one 43 nm repeat, the three crowns project heads toward different azimuthal positions around the thick filament. Not all heads face the same target thin filament. In the hexagonal lattice, each thick filament is surrounded by six thin filaments, and each crown's heads preferentially interact with two of these six.
- **Cooperative realignment.** Filament compliance provides a coupling mechanism between discrete binding sites. When one head binds and generates force, it shifts the thick filament axially relative to the thin filament. This shift can bring neighboring (previously misaligned) heads into register with actin sites, creating a cooperative recruitment cascade (Daniel et al. 1998).

The current model assumes a **continuous distribution** of cross-bridges (the Huxley 1957 convention), which overestimates the fraction of heads that can simultaneously bind and does not capture the vernier mismatch effects. A spatially explicit model with discrete crown positions would naturally produce lower duty ratios and cooperative binding effects without additional fitting parameters.

**Relevant protocols.** X-ray diffraction patterns during contraction, duty ratio, force per head, ktr.

**Literature support.** Daniel et al. (1998) — spatially explicit model with discrete binding sites and compliant filaments (*Biophys J* 74: 1611–1621). Piazzesi et al. (2007) — X-ray interference revealing discrete crown behavior in intact skeletal muscle; concluded that force modulation occurs via number of attached motors, not motor force or stroke size (*Cell* 131: 784–795). Squire (2017) — comprehensive review of myosin filament structure including crown spacing and IHM packing (*Biophys Rev* 9: 197–207). Huxley et al. (1994) — X-ray diffraction measurements of extensibility of actin and myosin filaments in contracting muscle (*Biophys J* 67: 2411–2421).

**Suggested implementation.** Full spatially explicit implementation would require tracking individual myosin heads at defined positions along the thick filament — a major computational increase incompatible with the current strain-distribution framework. A pragmatic intermediate approach is to apply a **geometric availability factor** to the attachment rate: `ka_eff = ka × p_aligned(dS)`, where `p_aligned` is the fraction of heads whose axial offset from the nearest actin site falls within the binding window, computed from the vernier geometry of 14.3 nm crown spacing vs. 5.5 nm actin repeat.

#### 1.12.9 Thick Filament Mechanosensing (OFF/ON Transition)

**Physical motivation.** In relaxed muscle, the majority of myosin heads are folded back against the thick filament backbone in a helically ordered "OFF" state — the structural basis of the super-relaxed state (SRX, Section 1.3). The structural motif is the **interacting heads motif (IHM)**, where the two heads of each myosin molecule interact with each other and with the S2 tail, parking themselves against the backbone with very low ATPase activity (~10× lower than disordered relaxed heads).

Activation involves **two parallel regulatory systems**:
1. **Thin filament regulation** (classical): Ca²⁺ binds troponin C → tropomyosin shifts → actin binding sites exposed.
2. **Thick filament regulation** (mechanosensing): Mechanical stress on the thick filament backbone — from initial cross-bridge binding or from titin-transmitted passive tension — triggers a structural transition from the OFF (IHM) to ON (disordered) state, releasing additional myosin heads and making them available for binding.

Linari et al. (2015) demonstrated this directly: at very low tensions (~0.1 T₀), the OFF-state X-ray signatures were partially lost, but imposing rapid shortening to maintain zero tension preserved the OFF state — demonstrating that the transition is **force-dependent**, not simply calcium-dependent. In cardiac muscle, Brunello et al. (2020) showed that thick filament structural changes during contraction are a key regulatory mechanism.

The current model implements SRX (Section 1.3) with force-dependent mobilization rates, which is a simplified version of this mechanism. However, the model's SRX exchange rates (ksr0 = 2.4 s⁻¹) are 600× faster than the biochemical SRX lifetime, meaning the current implementation captures the regulatory buffering role but not the true mechanosensing dynamics.

**Relevant protocols.** Length-dependent activation (Frank–Starling), ktr at different SL, force-calcium relationship.

**Literature support.** Linari et al. (2015) — force generation controlled by mechanosensing in myosin filaments (*Nature* 528: 276–279). Brunello et al. (2020) — myosin filament-based regulation of contraction dynamics in heart muscle (*Proc Natl Acad Sci USA* 117: 8177–8186). Irving (2017) — regulation of contraction by the thick filaments in skeletal muscle (*Biophys J* 113: 2579–2594). Fusi et al. (2016) — thick filament mechanosensing is a calcium-independent regulatory mechanism (*Nat Commun* 7: 13281).

**Suggested implementation.** The existing SRX framework (Section 1.3) already captures this mechanism qualitatively via force-dependent ksr/kmsr rates. To improve physical fidelity: (a) reduce ksr0 to a biochemically realistic baseline (~0.01 s⁻¹) and make the force-dependent acceleration steeper, so that at zero force the SRX pool is genuinely stable but under contraction forces it mobilizes rapidly; (b) couple titin passive tension to the SRX mobilization rate, providing a length-dependent activation pathway independent of calcium.

#### 1.12.10 Titin-Based Thick Filament Activation

**Physical motivation.** Titin spans from the Z-line to the M-line. Its I-band segment acts as a molecular spring providing passive tension (Section 1.6). Critically, titin transmits passive force **directly to the thick filament backbone**, and this tension can contribute to the OFF-to-ON transition of myosin heads (Section 1.12.9). At longer sarcomere lengths, higher titin-based passive tension promotes thick filament activation even **before calcium arrival**, providing a mechanism for length-dependent activation (the Frank–Starling mechanism) that is independent of lattice spacing effects (Section 1.12.6).

Marcucci et al. (2017) showed that incorporating titin-mediated thick filament activation into a mathematical model of rat trabeculae introduced sarcomere-length dependencies that matched experimental observations — including the steep length–tension relationship and the SL-dependence of ktr.

**Relevant protocols.** Length-tension relationship, ktr dependence on SL, diastolic function.

**Literature support.** Marcucci et al. (2017) — titin-mediated thick filament activation introduces sarcomere-length dependencies in mathematical models of rat trabecula and whole ventricle (*Sci Rep* 7: 5546). Ait-Mou et al. (2016) — titin strain contributes to the Frank–Starling law by structural rearrangements of both thin- and thick-filament proteins (*Proc Natl Acad Sci USA* 113: 2306–2311). Fusi et al. (2016) — thick filament mechanosensing is calcium-independent (*Nat Commun* 7: 13281).

**Suggested implementation.** Couple the passive titin force (Section 1.6) to the SRX mobilization rate: `kmsr_eff = kmsr × exp((F_active + alpha × F_passive) / sigma1)`, where the factor `alpha` controls the relative contribution of passive (titin-transmitted) tension vs. active (cross-bridge-transmitted) tension to thick filament activation. Currently, only F_active (or F_total) drives the SRX mobilization. Adding the passive contribution with a single parameter alpha would introduce titin-mediated length-dependent activation.

#### 1.12.11 Myosin Binding Protein C (MyBP-C)

**Physical motivation.** Cardiac myosin binding protein C (cMyBP-C) is a thick filament-associated protein localized to the **C-zone** — the central third of each half-thick filament (approximately 9 of the ~50 myosin crowns per half-thick filament). It has dual regulatory roles:

- **When dephosphorylated:** cMyBP-C binds to the myosin S2 segment, stabilizing heads in the OFF/SRX conformation and acting as a **brake** on cross-bridge cycling. It also bridges to actin, slowing thin filament sliding. The net effect is reduced attachment rate and prolonged attachment duration in the C-zone.
- **When phosphorylated** (by PKA during β-adrenergic stimulation): cMyBP-C releases myosin heads from the thick filament surface, accelerating cross-bridge recruitment and the rate of force development. PKA phosphorylation of cMyBP-C is a major inotropic mechanism in the heart — it accelerates both contraction and relaxation.

This provides a mechanism for modulating contractility that is **independent of calcium levels** — directly relevant to inotropic responses and potentially to the different rates of force development observed at different activation levels.

The current model does not distinguish between C-zone and non-C-zone regions of the thick filament. All myosin heads are treated identically, and there is no cMyBP-C regulatory layer. For skinned fiber experiments at saturating calcium (the current target protocols), cMyBP-C phosphorylation state may still influence ktr and the kinetics of force redevelopment.

**Relevant protocols.** ktr (rate of force development), FV at submaximal activation, relaxation kinetics.

**Literature support.** Previs et al. (2016) — phosphorylation of cMyBP-C releases myosin heads from the surface of cardiac thick filaments (*Proc Natl Acad Sci USA* 113: 3235–3240). Kampourakis et al. (2014) — cMyBP-C activates thin filaments and inhibits thick filaments in heart muscle cells (*Proc Natl Acad Sci USA* 111: 18763–18768). Razumova et al. (2006) — cardiac MyBP-C regulates the rate and force of contraction in mammalian myocardium (*Circ Res*).

**Suggested implementation.** Split the myosin head pool into C-zone (~30%) and non-C-zone (~70%) subpopulations. Apply a reduced effective attachment rate and increased SRX stability to C-zone heads: `ka_Czone = ka × (1 − f_brake)`, where `f_brake` depends on cMyBP-C phosphorylation state (0 = fully dephosphorylated/braked, 1 = fully phosphorylated/released). For skinned fiber experiments, f_brake would be set to the experimental phosphorylation condition.

---

## Part 2 — Parameter Range Comparison to Literature

The following table presents the key kinetic and mechanical parameters. "Effective value" accounts for the global multiplier `xrate = 0.435` applied to kah, kamh, ka, kd, k1, k_1, k2. Force stiffness values (kstiff1, kstiff2) are reported in internal model units (kPa·μm⁻¹·head⁻¹ or equivalent). Temperatures for literature values are noted; most cardiac mechanical data are at 15–25°C or 37°C (human).

| Parameter | Current Value (opt) | Effective Value (×xrate) | Literature Range | Primary Citation | Status |
|---|---|---|---|---|---|
| **kah** (ATP hydrolysis, s⁻¹) | 18.79 | 8.18 | 10–100 s⁻¹ (cardiac, 37°C) | Siemankowski et al. 1985; Malmqvist et al. 2004 | ⚠️ Below range |
| **kamh** (reverse hydrolysis, s⁻¹) | 2.44 | 1.06 | ~1–10 s⁻¹ | Eisenberg and Hill 1978 | ⚠️ Below range |
| **ka** (attachment rate, s⁻¹) | 67.0 | 29.2 | 50–500 s⁻¹ | Brenner 1988; Land 2017 | ⚠️ Below range |
| **kd** (P1→PD detachment, s⁻¹) | 9.56 | 4.16 | 5–200 s⁻¹ | Little 2012 (human 40 s⁻¹, rat 150 s⁻¹) | ⚠️ Below range |
| **k1** (P1→P2 power stroke, s⁻¹) | 270.8 | 117.9 | 100–1000 s⁻¹ (estimated from tension transients) | Ford et al. 1977; Piazzesi et al. 2002 | ✅ In range |
| **k_1** (P2→P1 reverse, s⁻¹) | 17.1 | 7.4 | < k1; 10–100 s⁻¹ | Eisenberg and Hill 1978 | ⚠️ Low end |
| **k2** (P2→P3 ADP release, s⁻¹) | 133.5 | 58.1 | 100–400 s⁻¹ (cardiac ADP release) | Siemankowski et al. 1985; Malmqvist et al. 2004 | ⚠️ Low end |
| **k3** (P3→PD ATP binding / final detach, s⁻¹) | 44.3 | 44.3 | 5–50 s⁻¹ (ATPase kcat at saturating ATP) | Siemankowski et al. 1985; Lymn and Taylor 1971 | ✅ In range |
| **gamma** (passive force exponent) | 2.67 | — | 7–10 (cardiac N2B PEVK) | Granzier and Irving 1995; Wu et al. 2000; Terui et al. 2008 | ❌ Too low |
| **L_thin** (thin filament, μm, as total length) | 1.377 | — | 1.0–1.2 μm (rat/human cardiac EM) | Tonino et al. 2017; Robinson and Winegrad 1979 | ❌ Too long |
| **L_thick** (thick filament, μm, as total length) | 1.14 | — | 1.60–1.65 μm (cardiac EM) | Tregear and Squire 1973; Luther et al. 2008 | ❌ Too short |
| **kstiff1** (P1 stiffness, internal units) | 16000 | — | ~1–3 pN/nm per XB = 1000–3000 in model units | Veigel et al. 2003; Mehta et al. 1999 | ✅ Order-of-magnitude plausible |
| **kstiff2** (P2 stiffness, internal units) | 97600 | — | ~2–4 pN/nm per XB = 2000–4000; ratio kstiff2/kstiff1 ~ 6× in model | Kaya and Higuchi 2010; Linari et al. 2009 | ⚠️ Ratio kstiff2/kstiff1 = 6; single-molecule ratio ~1.5–3 |
| **kstiff1_n** (P1 negative-strain stiffness) | 160 | — | No measurement; convention | — | ❌ Not experimentally supported |
| **kstiff2_n** (P2 negative-strain stiffness) | 704 | — | No measurement; convention | — | ❌ Not experimentally supported |
| **ksr0** (SRX off rate, s⁻¹) | 2.4 | — | 0.004 s⁻¹ (1/230 s SRX lifetime) | McNamara 2015 | ❌ 600× too fast |
| **kmsr** (SRX on rate under force, s⁻¹) | 1.0 | — | Not directly measured | — | ⚠️ No constraint |
| **ksr2srd** (SR→SRD, s⁻¹) | 32.1 | — | Not measured | — | ❌ No experimental basis |
| **mu** (viscosity coefficient) | 0.633 | — | Internal viscosity from Tombe 1991 (25°C rat cardiac) | Tombe 1991 | ✅ Qualitatively consistent |
| **xrate** (global rate multiplier) | 0.435 | — | Should be 1.0 if base rates are correct | — | ❌ Indicates base rates are systematically too high by ~2× |
| **dr** (power stroke size, μm) | 0.010 | — | 5–11 nm (single molecule) | Guilford et al. 1997; Capitanio et al. 2006; Mehta et al. 1999 | ✅ In range (10 nm) |
| **vmax** (set Vmax, μm/s per HS) | 10 | — | 10 μm/s at 25°C (Tombe 1991); 1.5 ML/s = ~3 μm/s at 12–15°C (McDonald 1998; Difee 2003) | Tombe 1991; McDonald 1998 | ⚠️ Temperature-dependent; correct for 25°C |

**Notes on xrate.** The global multiplier xrate = 0.435 reduces all major cycling rates by approximately 2.3-fold before each ODE integration. If the base rates in the parameter file were correctly calibrated to physiological values, xrate would equal 1. An xrate < 1 indicates that the unmodified base rates give a cycle time that is too fast (force rises too rapidly, ktr too high, Vmax too high). The effective kah (8.18 s⁻¹) is below the lower end of the published ATPase range (10 s⁻¹) for cardiac myosin at 37°C. This suggests the model is running at effectively below-physiological ATPase rates even at 37°C, implying either (a) the protocol is designed to match lower-temperature data (consistent with McDonald 1998 at 12°C and Difee 2003 at 15°C), or (b) the rate constants need revision, or (c) xrate serves as a temperature correction factor (Q10 correction for ~25°C data with base rates parameterized for 37°C).

---

## Part 3 — Protocol Sensitivity Analysis

### 3.1 FV Curvature and Vmax

**Physical interpretation.** The force-velocity curve is generated by running the ODE at constant sarcomere shortening velocity `v` while force settles to a new steady state, then plotting F(v). The Hill hyperbola `(F + a)(v + b) = (F0 + a)b` is the classical fit, where `a/F0` is the curvature parameter (lower a/F0 = more curved, less mechanical efficiency).

**Mechanisms controlling Vmax (velocity at zero load).**
- The primary determinant is the rate-limiting step in the cross-bridge cycle when force = 0. For most cross-bridge models, this is the ADP release rate (k2) or the net cycling rate (kcat). Vmax ≈ kcat × dr × N_XB, where N_XB is the number of attached heads per half-sarcomere length. The current effective k2 = 58 s⁻¹ combined with dr = 10 nm predicts Vmax ≈ 58 × 0.01 × N_XB μm/s. For N_XB ~ 30–50%, Vmax ~ 0.3–0.6 μm/s, far below vmax = 10 μm/s. This discrepancy is resolved by the internal viscous drag (mu) and serial stiffness acting as a governor, but the mechanical inconsistency indicates that either k2 or mu is functioning as a free parameter rather than a physically constrained value.
- The serial stiffness (kSE, Section 1.7) affects the apparent Vmax measured at the muscle ends vs. the sarcomere level. At high shortening velocities, the SE extends, causing apparent velocity at the fiber ends to exceed sarcomere shortening velocity.
- The viscous drag (mu) directly limits Vmax: `velHS = (Force − F_total)/mu`, so at zero load and zero F_total, velHS = Force/mu. This is the dominant governor of Vmax in the current model.
- Strain-dependent detachment (k1, k2 PCHIP) affects the width of the attached population during shortening, which governs the drag force at high velocities.

**Mechanisms controlling curvature (a/F0).**
- The ratio of attached P2 heads at near-zero strain to P2 heads at high negative strain determines how rapidly force falls with velocity. If negative-strain heads accumulate (kstiff2_n × p2_1_neg is large), the curve becomes less curved (more linear). The negative stiffness asymmetry (kstiff2_n = 704 << kstiff2 = 97600) specifically acts to reduce the contribution of negative-strain P2 heads, making the curve more curved (lower a/F0).
- The A2 attachment shift (Section 1.9) also affects FV curvature by redistributing high-strain P2 heads toward lower strain, reducing the negative-force contribution and potentially changing the curvature.
- ADP release rate (k2) and its strain dependence (PieceWiseStrainDep2) directly modulate the force-velocity relationship.

### 3.2 ktr and Force Redevelopment

**Physical interpretation.** The ktr protocol involves rapidly shortening the muscle to near-zero force (slack), then rapidly restretching to the original length, and measuring the rate of force redevelopment. The rate constant ktr is measured by fitting F(t) = F∞ × (1 − exp(−ktr × t)) to the rising phase. ktr reflects the rate of cross-bridge cycling in a muscle that has had all attached heads rapidly detached (by the slack shortening).

**Mechanisms controlling ktr.**
- The dominant determinants of ktr are the attachment rate (ka, R12 = k1) and the global cycling rate (kcat = kah if ATP hydrolysis is rate-limiting). Brenner (1988) showed that ktr ≈ fapp + gapp, where fapp and gapp are the apparent attachment and detachment rate constants. In the current model, ktr will be sensitive primarily to ka (29 s⁻¹ effective), k1 (117.9 s⁻¹ effective), and kd (4.16 s⁻¹ effective).
- The SRX state modulates ktr: heads in P_SR or P_SRD cannot be rapidly recruited during restretch (their mobilization rates ksr0 = 2.4 s⁻¹, kmsr = 1.0 s⁻¹ are slow relative to ktr ~ 10–30 s⁻¹). Thus the SRX fraction effectively reduces the pool of heads available for rapid force redevelopment, lowering ktr.
- Filament overlap determines how many heads are in the overlap zone and available for attachment.
- The thin filament activation state (Ca saturation) is assumed to be maximal, so Ca does not limit ktr in these protocols.

**Expected discrepancy.** If ktr measured experimentally is ~10–30 s⁻¹ at 25–37°C (Brenner 1988; Wolff et al. 1995), and the effective attachment rate in the model is ka ~ 29 s⁻¹, ktr should be in the right range. However, if the force plateau is too high (because passive stiffness is too low and the overlap factor is incorrect), the apparent ktr will be distorted. Additionally, if the SRD pool (ksr2srd = 32.1 s⁻¹) rapidly exchanges during the ktr protocol, it will add a fast component to force redevelopment that may not match experimental data.

### 3.3 Slack-Restretch Transients

**Physical interpretation.** In the slack-restretch protocol, the muscle is rapidly shortened by a preset length (slack), held until force falls to zero or near zero (allowing attached XBs to detach), then rapidly restretched. The force transient has several phases: (a) an initial rapid force drop during shortening (governed by serial stiffness, viscous drag, cross-bridge detachment), (b) a slow fall to zero if any residual force remains after slack, (c) a rapid force rise during restretch (governed by attachment rate), and (d) a slow plateau to the original force level.

**Mechanisms controlling slack transients.**
- The shortening phase is dominated by the negative-strain detachment rate (kd at s < 0), the viscous drag (mu), and the serial stiffness (kSE).
- The zero-force phase is controlled by the residual attached fraction and the spontaneous detachment rate at near-zero strain (kd at s ~ 0).
- The restretch force rise is governed by ka, k1 (attachment and power stroke), and the passive force contribution (titin). Passive force during restretch provides an elastic restoring force that contributes to the early rising phase — this is where gamma = 2.67 (too low) will cause the model to underestimate passive restoring force and predict a slower rise.
- The SRX/SRD pool dynamics affect the slow plateau phase: if SRX mobilization is slow, the plateau will be below the pre-slack force for several seconds, which is experimentally observed.

### 3.4 Mechanism-to-Observable Mapping Table

| Mechanism | FV Curvature | FV Vmax | ktr Rate | Slack Duration | Restretch Rise Rate | Restretch Plateau |
|---|---|---|---|---|---|---|
| ATPase rates (kah, xrate) | Moderate | High | High | Low | Low | Moderate |
| Attachment rate (ka) | Low | Moderate | High | Low | High | Moderate |
| Detachment rate (kd) | Moderate | Moderate | Moderate | High | Low | Low |
| Power stroke rate (k1) | Low | Moderate | High | Low | Moderate | Moderate |
| ADP release (k2) | Low | High | Moderate | Low | Moderate | Moderate |
| Strain-dependent PCHIP | High | Moderate | Moderate | Moderate | Low | Low |
| SRX (ksr0, kmsr) | Low | Low | Moderate | Low | Low | High |
| SRD (ksr2srd, kmsrd) | Low | Low | Low | Low | Low | High |
| Filament overlap (L_thick, L_thin) | Moderate | Low | High | Low | Moderate | High |
| Passive force (k_pas, gamma) | Low | Low | Low | Moderate | High | Moderate |
| Serial stiffness (kSE) | High | High | Low | High | Moderate | Low |
| Negative stiffness asym. (kstiff_n) | High | Moderate | Low | Low | Low | Low |
| A2 attachment shift (slope) | Moderate | Moderate | Low | Low | Low | Low |
| Viscous drag (mu) | Low | High | Low | High | Low | Low |
| Maxwell dashpot (currently OFF) | Low | Low | Low | High | Moderate | Moderate |
| Lattice spacing (not implemented) | Moderate | Low | Moderate | Low | Low | Moderate |
| Rotational myosin spring (not impl.) | High | Moderate | Low | Low | Low | Low |
| Crown spacing / vernier (not impl.) | Low | Low | Moderate | Low | Moderate | Moderate |
| Thick filament mechanosensing | Low | Low | Moderate | Low | Moderate | High |
| Titin-based TF activation (not impl.) | Low | Low | High | Low | Moderate | Moderate |
| MyBP-C regulation (not impl.) | Low | Low | Moderate | Low | Moderate | Low |
| Cooperative XB detachment (not impl.) | Moderate | Low | Moderate | Low | High | Moderate |

*Impact levels: High = primary determinant; Moderate = significant secondary effect; Low = minor effect.*

---

## Part 4 — Categorized Summary and Prioritized Action List

### Table A: Well-Supported Mechanisms

| Mechanism | Switch / Parameter | Evidence Summary | Citation |
|---|---|---|---|
| Three-state XB kinetic scheme | `NumberOfStates = 3` (topology) | Lymn–Taylor cycle established by X-ray, cryo-EM, single-molecule; 3-state simplification widely used | Lymn and Taylor 1971; Rayment et al. 1993; Land 2017 |
| Serial stiffness | `UseSerialStiffness`, kSE, ekSE | Rig compliance causes up to 20% SL change during contraction; essential for correct FV | Tombe 2016; Irving et al. 1992; Land 2017 |
| Viscous drag | mu, mu_neg | Internal viscous element directly measured in rat myocardium at 25°C; required to fit Vmax | Tombe 1991 |
| Strain-dependent kinetics (qualitative) | `UsePieceWiseStrainDep`, PCHIP | Positive strain inhibits A1→A2 and A2→A3; negative strain promotes — experimentally confirmed | Tewari (cited); Ford et al. 1977; Piazzesi et al. 1997 |
| Ca²⁺ / thin filament activation | K_coop, k_on, k_off, Ca = 1000 | Troponin C cooperativity well-established; saturating Ca appropriate for maximal activation protocols | Robertson et al. 1981; Dobesh et al. 2002; Brenner 1988 |
| Filament overlap (concept) | `UseOverlap` | Frank–Starling mechanism established; active force proportional to overlap zone | Gordon et al. 1966; Prostsenko 2020 |
| Power stroke size dr = 10 nm | dr | Single-molecule measurements: 5–11 nm for myosin II | Guilford et al. 1997; Mehta et al. 1999; Capitanio et al. 2006 |

### Table B: Questionable Mechanisms / Parameters

| Mechanism | Switch / Parameter | Evidence Summary | Citation |
|---|---|---|---|
| SRX exchange rates | ksr0 = 2.4 s⁻¹, kmsr = 1.0 s⁻¹ | Biochemical SRX lifetime ~230 s (rate ~0.004 s⁻¹); model SRX is 600× too fast to be a true SRX | McNamara 2015; Hooijman et al. 2011 |
| P2/P1 stiffness ratio | kstiff2/kstiff1 = 6 | Single-molecule data suggest pre- and post-stroke stiffness differ by ~1.5–3×, not 6× | Kaya and Higuchi 2010; Veigel et al. 2003 |
| Strain-dependent kinetics (quantitative) | PCHIP boundary values = 50 | Boundary multipliers lack experimental measurement; single-molecule strain-force coupling suggests factor ~5–10, not 50 | Guilford et al. 1997 |
| Filament lengths (L_thick, L_thin) | L_thick = 1.14 μm, L_thin = 1.377 μm | L_thin as total filament length (1.377 μm) exceeds published values (1.0–1.2 μm); L_thick too short | Tonino et al. 2017; Robinson and Winegrad 1979; Luther et al. 2008 |
| A2 attachment shift rate | slope = 16496 s⁻¹/μm | No experimental measurement; mechanism not in published cardiac models | — |
| Three-state scheme (P3 identity) | k3 = 44 s⁻¹, k3m = 0 | Biochemical identity of P3 state is ambiguous; mixes ADP release and slow detachment | Siemankowski et al. 1985 |
| Global xrate = 0.435 | xrate | A 2.3× scaling factor applied to all main rates indicates base rates are systematically incorrect; should not be a tuning parameter | — |

### Table C: Likely Missing Mechanisms

| Mechanism | Switch / Parameter | Evidence Summary | Citation |
|---|---|---|---|
| Correct passive exponent (gamma ~ 7.9) | `UseTitinIdentifiedPassive` = 0 (OFF) | Cardiac titin passive force exponent ~ 7–10; current gamma = 2.67 severely underestimates passive stiffness at SL > 2.3 μm | Granzier and Irving 1995; Wu et al. 2000; Terui et al. 2008 |
| Negative stiffness asymmetry (replace with fast detachment) | `UseNegativeKstiff`, kstiff1_n, kstiff2_n | No measurement supports 100× stiffness reduction at negative strain; correct physics is fast kd at s < 0 | Huxley 1957; Guilford et al. 1997 |
| Distortion reset at P1→P2 transition (Land 2017) | Not implemented | Computational simplification in Land 2017 (minor impact on results); would cause numerical smearing in strain-discretized models; prefer fast negative-strain detachment instead | Land 2017; Rice 2008 |
| Lattice spacing effects | Not implemented | Radial distance modulates attachment rate and contributes to length-dependent activation; 1D model ignores entirely | Tanner 2007/2012; Williams 2013 |
| Rotational myosin spring | Not implemented | Multi-spring geometry derives state-dependent and strain-sign-dependent stiffness from lever arm angle; would replace ad hoc kstiff_n | Huxley and Simmons 1971; Tanner 2012 |
| Discrete myosin crown spacing | Not implemented | 14.3 nm crown spacing creates vernier mismatch with 5.5 nm actin sites; continuous distribution overestimates bound fraction | Piazzesi 2007; Daniel 1998; Squire 2017 |
| Thick filament mechanosensing | Partially (SRX) | Force-dependent OFF→ON transition of IHM; current SRX rates 600× too fast for true mechanosensing | Linari 2015; Brunello 2020 |
| Titin-based thick filament activation | Not implemented | Passive tension activates thick filament independently of calcium; key for Frank–Starling | Marcucci 2017; Ait-Mou 2016 |
| MyBP-C regulation | Not implemented | C-zone brake/release of myosin heads; modulates ktr and contractility independent of calcium | Previs 2016; Kampourakis 2014 |
| Maxwell viscoelastic dashpot (titin viscoelasticity) | `UseMaxwellDashpot = 0` (OFF) | Titin is viscoelastic; slow dashpot force affects slack transient shape and restretch rising phase | Hamdani 2017; Freundt 2019; Linke 2018 |
| Ca²⁺-dependent titin stiffening | Not implemented | Titin stiffness increases ~20% in Ca²⁺; relevant for maximal activation protocols | Herzog 2022; Freundt 2019 |
| SRD mechanism with realistic rates | ksr2srd = 32.1 s⁻¹, ksrd2sr = 2.9 s⁻¹ | No experimental ADP-SRX kinetics in cardiac; current rates are orders of magnitude too fast | — |

---

### 4.4 Prioritized Action List

The following eight actions are listed in order of expected impact on reducing model-data discrepancy, based on the analysis in Parts 1–3.

**Action 1 — Activate correct passive force exponent (gamma = 7.9)**
- **Action**: Set `params0.UseTitinIdentifiedPassive = 1`. Re-optimize k_pas.
- **Expected effect**: Passive force at SL > 2.2 μm will increase significantly, improving the restretch-phase force transient and the passive component of ktr at elevated SL. This is a direct cause of underpredicted passive force in the slack-restretch protocol.
- **Affected protocol**: Slack-restretch (passive restoring force), ktr at long SL.
- **Relevant parameter**: `UseTitinIdentifiedPassive`, `k_pas`, `gamma`.

**Action 2 — Fix filament geometry (L_thick, L_thin)**
- **Action**: Set `L_thick` and `L_thin` to published values for rat cardiac: L_thin (total) ≈ 1.15 μm, L_thick (total) ≈ 1.65 μm (or L_thick half = 0.825 μm if the code uses half-lengths consistently). Verify the geometric formula in the overlap calculation to confirm whether the code interprets these as half-lengths or full lengths.
- **Expected effect**: Correcting filament geometry will fix the overlap vs. SL relationship, potentially requiring re-optimization of ka and k_pas. A correctly parameterized overlap will improve the length-dependence of force and ktr.
- **Affected protocol**: FV (force level at SL0), ktr (force plateau), length-tension.
- **Relevant parameter**: `L_thick`, `L_thin`, `L_hbare`.

**Action 3 — Replace negative stiffness asymmetry with fast negative-strain detachment**
- **Action**: Disable `UseNegativeKstiff`. Instead, increase the PCHIP multiplier for R1D and R12 at s < 0 to implement fast detachment at negative strain. Alternatively, increase the kd value and make the PCHIP for R1D asymmetric with a steep rise at s < −dr.
- **Expected effect**: Removes an unphysical stiffness asymmetry and replaces it with the standard Huxley mechanism (fast g at s < 0). Should improve FV curvature behavior and reduce artifacts in force transients during rapid shortening.
- **Affected protocol**: FV curvature, slack shortening force transient.
- **Relevant parameter**: `UseNegativeKstiff`, `kstiff1_n`, `kstiff2_n`, PCHIP (R1D).

**Action 4 — Implement distortion reset at P1→P2 (Land 2017 mechanism)**
- **Action**: In the dp2 / dp1 flux computation, when R12 removes probability from P1 at strain s, deposit it into P2 at strain `s + dr` (shifted by the power stroke size) rather than at the same strain bin. This is a one-line modification to the indexing in the dp2 accumulation.
- **Expected effect**: Prevents accumulation of post-stroke P2 heads at negative strain during shortening, directly improving FV curve shape at high velocities and reducing the artificial need for the kstiff2_n correction.
- **Affected protocol**: FV curve (Vmax, curvature), ktr (early force rise).
- **Relevant parameter**: dr = 0.01 μm (already exists). Modify dp2 accumulation at the R12 transition.

**Action 5 — Correct base rate values and eliminate xrate as a tuning parameter**
- **Action**: Using published literature values as anchors, revise the base rates (kah, ka, kd, k1, k2) so that xrate is close to 1.0. Start from ka ≈ 50–200 s⁻¹, kd ≈ 40–150 s⁻¹ (matching Little 2012 for human/rat), kah ≈ 10–30 s⁻¹, k2 ≈ 100–400 s⁻¹. Then set xrate = 1.0 and re-optimize.
- **Expected effect**: Establishes physically realistic rate constants, removing the systematic bias introduced by xrate. Will likely require full re-optimization of other parameters.
- **Affected protocol**: ktr rate, Vmax, FV plateau force.
- **Relevant parameter**: kah, ka, kd, k1, k2, xrate.

**Action 6 — Fix SRX rates to biochemically realistic values or remove SRX**
- **Action**: If SRX is retained as a model element, set ksr0 to be consistent with the biochemical SRX lifetime (~0.004–0.01 s⁻¹) and treat SRX mobilization as force-dependent with a faster effective rate that is still >> 0.004 s⁻¹ but justified by the force-dependence (Campbell 2018). Alternatively, if the SRX is meant to be a fast regulatory buffer (not true SRX), rename it and document clearly that it is not the biochemical SRX.
- **Expected effect**: Will change the ktr protocol behavior (SRX mobilization contributes to force redevelopment plateau). Correctly parameterized SRX may improve the slow plateau phase of force redevelopment.
- **Affected protocol**: ktr force plateau, slack force fall to zero.
- **Relevant parameter**: ksr0, kmsr, sigma1, sigma2.

**Action 7 — Activate Maxwell dashpot for titin viscoelasticity**
- **Action**: Set `UseMaxwellDashpot = 1`. Tune kSE_M and eta_M by fitting the passive force response during the slack phase (the slow exponential recovery of passive force after rapid shortening). Anchor eta_M to published titin viscoelastic measurements (Linke 2018; Hamdani 2017).
- **Expected effect**: Will add a viscoelastic time constant to the passive force, improving the shape of the slack-restretch force transient (slow fall and slow rise of passive force), and potentially reducing the mismatch in the restretch force peak.
- **Affected protocol**: Slack-restretch transient shape (both slack and restretch phases).
- **Relevant parameter**: `UseMaxwellDashpot`, kSE_M, eta_M.

**Action 8 — Remove or constrain the SRD state**
- **Action**: Remove `UseSuperRelaxedADP = 1` (set to 0) and test whether the remaining mechanisms can fit the ktr and slack data. If a slow secondary pool is needed, introduce it as a simple two-state relaxed/active model with rates constrained by the biochemical ADP-SRX lifetime, if that measurement becomes available.
- **Expected effect**: Reduces model complexity and removes an unconstrained mechanism. If the fits worsen substantially, it indicates that a second slow pool is genuinely needed, providing motivation for future experimental measurement of ADP-bound SRX kinetics in cardiac.
- **Affected protocol**: ktr (slow plateau), slack (slow zero-force phase).
- **Relevant parameter**: `UseSuperRelaxedADP`, ksr2srd, ksrd2sr, kmsrd, ksrd.

---

## References

- Ait-Mou, Y. et al. (2016). "Titin strain contributes to the Frank-Starling law of the heart by structural rearrangements of both thin- and thick-filament proteins." *Proc Natl Acad Sci USA* 113: 2306–2311.
- Akella, A.B. et al. (1997). "Diminished Ca²⁺ sensitivity of skinned cardiac muscle contractility coincident with troponin T-band shifts in the diabetic rat." *Biochim Biophys Acta* 1316: 109–116.
- Alamo, L. et al. (2008). "Three-dimensional reconstruction of tarantula myosin filaments suggests how phosphorylation may regulate myosin activity." *J Mol Biol* 384: 780–797.
- Beard, D.A. et al. (2022). "Reduced cardiac muscle power with low ATP simulating heart failure." *Biophys J* 121: 3213–3223. https://www.sciencedirect.com/science/article/pii/S0006349522006026
- Brenner, B. (1988). "Effect of Ca²⁺ on cross-bridge turnover kinetics in skinned single rabbit psoas fibers: implications for regulation of muscle contraction." *Proc Natl Acad Sci USA* 85: 3265–3269.
- Brunello, E. et al. (2020). "Myosin filament-based regulation of the dynamics of contraction in heart muscle." *Proc Natl Acad Sci USA* 117: 8177–8186.
- Campbell, K.S. (2018). "Force-dependent recruitment from the myosin off state contributes to length-dependent activation." *Biophys J* 115: 543–553. https://www.sciencedirect.com/science/article/pii/S0006349518307707
- Campbell, K.S. (2003). "Ising model of cardiac thin filament activation with nearest-neighbor cooperative interactions." *Biophys J* 84: 3897–3908.
- Campbell, S.G. et al. (2010). "Coupling of adjacent tropomyosins enhances cross-bridge-mediated cooperative activation in a Markov model of the cardiac thin filament." *Biophys J* 98: 2254–2264.
- Capitanio, M. et al. (2006). "Two independent mechanical events in the interaction cycle of skeletal muscle myosin with actin." *Proc Natl Acad Sci USA* 103: 87–92.
- Chin, L. et al. (2006). "Mathematical simulation of muscle cross-bridge cycle and force-velocity relationship." *Biophys J* 91: 3653–3663. https://www.sciencedirect.com/science/article/pii/S000634950672077X
- Daniel, T.L. et al. (1998). "Compliant realignment of binding sites in muscle: transient behavior and mechanical tuning." *Biophys J* 74: 1611–1621.
- Dantzig, J.A. et al. (1992). "Reversal of the cross-bridge force-generating transition by photogeneration of phosphate in rabbit psoas muscle fibres." *J Physiol* 451: 247–278.
- Decostre, V. et al. (2005). "Effect of 2,3-butanedione monoxime on the elementary steps of the cross-bridge cycle in frog skeletal muscle fibres." *J Physiol* 567: 767–785.
- Difee, G.M. et al. (2003). "Altered single cell force-velocity and power properties in exercise-trained rat myocardium." *J Appl Physiol* 94: 1941–1948. https://journals.physiology.org/doi/full/10.1152/japplphysiol.00889.2002
- Dobesh, D.P. et al. (2002). "Cooperative activation in cardiac muscle: impact of sarcomere length." *Am J Physiol Heart Circ Physiol* 282: H1271–H1282.
- Dominguez, R. et al. (1998). "Crystal structure of a vertebrate smooth muscle myosin motor domain and its complex with the essential light chain: visualization of the pre-power stroke state." *Cell* 94: 559–571.
- Duke, T.A.J. (1999). "Molecular model of muscle contraction." *Proc Natl Acad Sci USA* 96: 2770–2775.
- Egelman, E.H. and DeRosier, D.J. (1992). "Image analysis of the actin filament." *Ultramicroscopy* 40: 197–208.
- Eisenberg, E. and Hill, T.L. (1978). "A cross-bridge model of muscle contraction." *Prog Biophys Mol Biol* 33: 55–82.
- Ford, L.E. et al. (1977). "Tension responses to sudden length change in stimulated frog muscle fibres near slack length." *J Physiol* 269: 441–515.
- Fusi, L. et al. (2016). "Thick filament mechano-sensing is a calcium-independent regulatory mechanism in skeletal muscle." *Nat Commun* 7: 13281.
- Freundt, J.K. and Linke, W.A. (2019). "Titin as a force-generating muscle protein under regulatory control." *J Appl Physiol* 126: 1474–1482. https://journals.physiology.org/doi/full/10.1152/japplphysiol.00865.2018
- Gordon, A.M. et al. (1966). "The variation in isometric tension with sarcomere length in vertebrate muscle fibres." *J Physiol* 184: 170–192.
- Granzier, H.L. and Irving, T.C. (1995). "Passive tension in cardiac muscle: contribution of collagen, titin, microtubules, and intermediate filaments." *Biophys J* 68: 1027–1044.
- Granzier, H.L. and Labeit, S. (2004). "The giant protein titin: a major player in myocardial mechanics, signaling, and disease." *Circ Res* 94: 284–295.
- Guilford, W.H. et al. (1997). "Smooth muscle and skeletal muscle myosins produce similar unitary forces and displacements in the laser trap." *Biophys J* 72: 1006–1021.
- Guo, B. and Bhattacharya, S. (2012). "Catch bond mechanism in actomyosin." *J Chem Phys* 136: 015103.
- Hamdani, N. et al. (2017). "Tampering with springs: phosphorylation of titin affecting the mechanical function of cardiomyocytes." *Biophys Rev* 9: 225–237. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5498327/
- Herzog, W. (2022). "What can we learn from single sarcomere and myofibril preparations?" *Front Physiol* 13: 837611. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9092595/
- Holmes, K.C. et al. (1990). "Atomic model of the actin filament." *Nature* 347: 44–49.
- Hooijman, P. et al. (2011). "A new state of cardiac myosin with very slow ATP turnover: a potential cardioprotective mechanism in the heart." *Biophys J* 100: 1969–1976.
- Huxley, A.F. and Simmons, R.M. (1971). "Proposed mechanism of force generation in striated muscle." *Nature* 233: 533–538.
- Huxley, A.F. (1957). "Muscle structure and theories of contraction." *Prog Biophys Biophys Chem* 7: 255–318.
- Huxley, H.E. et al. (1994). "X-ray diffraction measurements of the extensibility of actin and myosin filaments in contracting muscle." *Biophys J* 67: 2411–2421.
- Irving, M. et al. (1992). "Myosin head movements are synchronous with the elementary force-generating process in muscle." *Nature* 357: 156–158.
- Irving, T.C. (2017). "Thick-filament strain and inter-filament spacing in cardiac muscle: changes with sarcomere length and calcium activation." *Biophys J* 112: 1698–1708.
- Kampourakis, T. et al. (2014). "Myosin binding protein-C activates thin filaments and inhibits thick filaments in heart muscle cells." *Proc Natl Acad Sci USA* 111: 18763–18768.
- Kaya, M. and Higuchi, H. (2010). "Nonlinear elasticity and an 8-nm working stroke of single myosin molecules in myofilaments." *Science* 329: 686–689.
- Lan, G. and Sun, S.X. (2005). "Dynamics of myosin-driven skeletal muscle contraction: I. Steady-state force generation." *Biophys J* 88: 4107–4117.
- Land, S. et al. (2017). "A model of cardiac contraction based on novel measurements of tension development in human cardiomyocytes." *J Mol Cell Cardiol* 106: 68–83. https://www.sciencedirect.com/science/article/abs/pii/S0022282817300639
- Linari, M. et al. (2015). "Force generation by skeletal muscle is controlled by mechanosensing in myosin filaments." *Nature* 528: 276–279.
- Linke, W.A. (2018). "Titin gene and protein functions in passive and active muscle." *Annu Rev Physiol* 80: 389–411.
- Little, S.C. et al. (2012). "The rates of Ca²⁺ dissociation and cross-bridge detachment from ventricular myofibrils as reported by a fluorescent cardiac troponin C." *J Biol Chem* 287: 27580–27589. https://www.jbc.org/article/S0021-9258(20)47816-0/fulltext
- Luther, P.K. et al. (2008). "Direct visualization of myosin-binding protein C bridging myosin and actin filaments in intact muscle." *Proc Natl Acad Sci USA* 108: 11423–11428.
- Lymn, R.W. and Taylor, E.W. (1971). "Mechanism of adenosine triphosphate hydrolysis by actomyosin." *Biochemistry* 10: 4617–4624.
- Marcucci, L. et al. (2017). "Titin-mediated thick filament activation, through a mechanosensing mechanism, introduces sarcomere-length dependencies in mathematical models of rat trabecula and whole ventricle." *Sci Rep* 7: 5546.
- Malmqvist, U.P. et al. (2004). "Cardiac beta-myosin heavy chain isoform changes are coupled with altered coenzyme A metabolites in diabetic rats." *J Biol Chem* 279: 16507–16513.
- McDonald, K.S. et al. (1998). "Force-velocity and power-load curves in rat skinned cardiac myocytes." *J Physiol* 511: 519–531. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2231141/
- McNamara, J.W. et al. (2015). "The role of super-relaxed myosin in skeletal and cardiac muscle." *Biophys Rev* 7: 5–14. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5425749/
- Mehta, A.D. et al. (1999). "Myosin-V is a processive actin-based motor." *Nature* 400: 590–593.
- Morimoto, S. (1991). "Interaction of troponin T with troponin C and the effect of tropomyosin on this interaction." *J Biochem* 109: 120–126.
- Nag, S. et al. (2017). "Contractility parameters of human beta-cardiac myosin with the hypertrophic cardiomyopathy mutation R403Q show loss of motor function." *Sci Adv* 1: e1500511.
- Nishizaka, T. et al. (1995). "Unbinding force of a single motor molecule of muscle measured using optical tweezers." *Nature* 377: 251–254.
- Piazzesi, G. et al. (2007). "Skeletal muscle performance determined by modulation of number of myosin motors rather than motor force or stroke size." *Cell* 131: 784–795.
- Piazzesi, G. et al. (1997). "Tension transients during steady lengthening of tetanized muscle fibres of the frog." *J Physiol* 499: 595–612.
- Piazzesi, G. et al. (2002). "Mechanism of force generation by myosin heads in skeletal muscle." *Nature* 415: 659–662.
- Powers, J.D. et al. (2018). "A spatially explicit model shows how titin stiffness modulates muscle mechanics and energetics." *Integr Comp Biol* 58: 186–193. https://academic.oup.com/icb/article/58/2/186/5036269
- Prostsenko, O.M. et al. (2020). "Changes in rat myocardium contractility under subchronic intoxication with lead and cadmium salts administered alone or in combination." *Toxicol Rep* 7: 626–633. https://www.sciencedirect.com/science/article/pii/S2214750020300445
- Previs, M.J. et al. (2016). "Phosphorylation of cardiac myosin binding protein C releases myosin heads from the surface of cardiac thick filaments." *Proc Natl Acad Sci USA* 113: 3235–3240.
- Ranatunga, K.W. (2010). "Cross-bridge mechanism(s) examined by temperature perturbation studies on muscle." In: *Advances in Experimental Medicine and Biology* 682: 247–266. https://link.springer.com/chapter/10.1007/978-1-4419-6366-6_14
- Razumova, M.V. et al. (1999). "Stiffness-distortion sarcomere model for muscle simulation." *J Appl Physiol* 87: 1861–1876.
- Razumova, M.V. et al. (2000). "Different myofilament nearest-neighbor interactions have distinctive effects on contractile behavior." *Biophys J* 78: 3120–3137.
- Razumova, M.V. et al. (2006). "Cardiac myosin-binding protein C regulates the rate and force of contraction in mammalian myocardium." *Circ Res*.
- Rayment, I. et al. (1993). "Three-dimensional structure of myosin subfragment-1: a molecular motor." *Science* 261: 50–58.
- Rice, J.J. et al. (2008). "Approximate model of cooperative activation and crossbridge cycling in cardiac muscle using ordinary differential equations." *Biophys J* 95: 2368–2390.
- Rice, J.J. et al. (2003). "Approximate model of cooperative activation and crossbridge cycling in cardiac muscle using ordinary differential equations." *Biophys J* 84: 2149–2166.
- Robertson, S.P. et al. (1981). "The effect of calcium and S-1 on the number of sites available to react with antibody on actin in regulated actomyosin." *J Biol Chem* 256: 12866–12871.
- Robinson, T.F. and Winegrad, S. (1979). "The measurement and dynamic implications of thin filament lengths in heart muscle." *J Physiol* 286: 607–619.
- Squire, J.M. (2017). "Muscle myosin filaments: cores, crowns and couplings." *Biophys Rev* 9: 197–207.
- Siemankowski, R.F. et al. (1985). "ADP dissociation from actomyosin subfragment 1 is sufficiently slow to limit the unloaded shortening velocity in vertebrate muscle." *Proc Natl Acad Sci USA* 82: 658–662.
- Stewart, M.A. et al. (2010). "Myosin ATP turnover rate is a mechanism involved in thermogenesis in resting skeletal muscle fibers." *Proc Natl Acad Sci USA* 107: 430–435.
- Tanner, B.C.W. et al. (2007). "Axial and radial forces of cross-bridges depend on lattice spacing." *PLoS Comput Biol* 3: e115.
- Tanner, B.C.W. et al. (2012). "Sarcomere lattice geometry influences cooperative myosin binding in muscle." *PLoS Comput Biol* 8: e1001018.
- Terui, T. et al. (2008). "Regulatory mechanism of length-dependent activation in skinned porcine ventricular muscle: role of thin filament cooperativity." *J Gen Physiol* 131: 241–250.
- Tombe, P.P. de and Ter Keurs, H.E.D.J. (1991). "An internal viscous element limits unloaded velocity of sarcomere shortening in rat myocardium." *J Physiol* 454: 619–642. https://physoc.onlinelibrary.wiley.com/doi/pdfdirect/10.1113/jphysiol.1992.sp019283
- Tombe, P.P. de et al. (2016). "Cardiac muscle mechanics: Sarcomere length matters." *J Mol Cell Cardiol* 91: 148–150. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5457809/
- Tonino, P. et al. (2017). "The giant protein titin regulates the length of the striated muscle thick filament." *J Mol Biol* 429: 3584–3595.
- Tregear, R.T. and Squire, J.M. (1973). "Myosin content and filament structure in smooth and striated muscle." *J Mol Biol* 77: 279–290.
- Trivedi, D.V. et al. (2018). "Direct observations of beta-cardiac myosin moving along thin filaments of cardiac muscle and the regulatory role of myosin-binding protein C." *Proc Natl Acad Sci USA* 115: E7188–E7197.
- Veigel, C. et al. (2003). "Load-dependent kinetics of force production by smooth muscle myosin measured with optical tweezers." *Nat Cell Biol* 5: 980–986.
- Witt, C.C. et al. (2006). "Nebulin regulates thin filament length, contractility, and Z-disk structure in vivo." *EMBO J* 25: 3843–3855.
- Wolff, M.R. et al. (1995). "Rate of tension development in cardiac muscle varies with level of activator calcium." *Circ Res* 76: 154–160.
- Woodhead, J.L. et al. (2005). "Atomic model of a myosin filament in the relaxed state." *Nature* 436: 1195–1199.
- Williams, C.D. et al. (2012). "Elastic energy storage and radial forces in the myofilament lattice depend on sarcomere length." *Biophys J* 103: 677–686.
- Williams, C.D. et al. (2013). "The length-tension curve in muscle depends on lattice spacing." *Proc R Soc B* 280: 20130697.
- Wu, Y. et al. (2000). "Titin isoform-dependent modulation of the amount of myosin binding to thick filaments in myocardium." *Arch Biochem Biophys* 381: 279–290.
