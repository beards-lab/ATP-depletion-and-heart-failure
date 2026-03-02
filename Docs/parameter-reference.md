# Parameter & Mechanism Reference

Current values from `ModelOptParams/ModelParamsSlackKtrOpt.m` (Feb 26 2026, most recent parametrization).
ODE function: `Model/dPUdT_CombinedTransitions.m`.

---

## Table 1 — Kinetic & Mechanical Parameters

### 1.1 Cross-Bridge Cycling Rate Constants

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `kah` | ATP hydrolysis rate: M·ATP → M·ADP·Pi (s⁻¹) | 18.79 | 10–100 s⁻¹ | Controls how fast the head cocks; cardiac β-MyHC is slow |
| `kamh` | Reverse hydrolysis: M·ADP·Pi → M·ATP (s⁻¹) | 2.44 | 1–20 s⁻¹ | Usually ≈ 0.05–0.15 × kah |
| `ka` | Attachment: PD → P1 (s⁻¹) | 67.0 | 50–500 s⁻¹ | Rate of weak-to-strong binding; Ca-regulated in reality |
| `kd` | Detachment from P1 (strain-independent base, s⁻¹) | 9.56 | 5–200 s⁻¹ | Unloaded detachment of loosely attached heads |
| `k1` | Power stroke forward: P1 → P2 (s⁻¹) | 270.8 | 100–2000 s⁻¹ | Rate-limiting for force development at short times |
| `k_1` | Power stroke reverse: P2 → P1 (s⁻¹) | 17.1 | 0–200 s⁻¹ | Back-stepping; low in cardiac muscle |
| `k2` | Post-stroke forward: P2 → P3 (s⁻¹) | 133.5 | 50–1000 s⁻¹ | Transition to the final bound state; strain-sensitive |
| `k_2` | P3 → P2 reverse (s⁻¹) | 2.79 | 0–50 s⁻¹ | Negligible in most conditions |
| `k3` | ATP binding + detachment: P3 → PT (s⁻¹) | 44.3 | 20–200 s⁻¹ | Limits duty ratio; cardiac ATPase ~1–3 s⁻¹ per head |
| `k3m` | P3 → P2 reverse re-entry (s⁻¹) | 0 | 0–20 s⁻¹ | Usually neglected |
| `k2rip` | Force-induced P2 detachment (s⁻¹) | 0 | 0–500 s⁻¹ | Off in current config |
| `xrate` | Global rate multiplier (dimensionless) | 0.435 | 0.1–2 | Single scalar scaling all turnover rates; compensates for temperature or isoform |
| `vmax` | Maximum unloaded shortening velocity (μm/s per half-sarcomere) | 10 | 4–25 μm/s | β-MyHC ~0.4 FL/s; ~6–10 μm/s at SL=2 μm |

### 1.2 Super-Relaxed (SRX) State Transitions

Myosin heads in the SRX state are folded back onto the thick filament (IHM conformation), with negligible ATPase activity. Force can recruit them back to the available (disordered-relaxed / PT) pool.

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `ksr0` | Base rate PT → SR (s⁻¹) | 2.40 | 0.1–20 s⁻¹ | Spontaneous SRX entry; cardiac SRX t½ ≈ 30–100 s at rest |
| `kmsr` | Base rate SR → PT (s⁻¹) | 1.0 | 0.1–10 s⁻¹ | Force-independent SRX exit |
| `sigma1` | Force sensitivity of SR → PT (kPa⁻¹) | 22.0 | 5–100 kPa⁻¹ | Higher force recruits more heads from SRX |
| `sigma2` | Force sensitivity of PT → SR (kPa⁻¹) | 40.4 | 5–100 kPa⁻¹ | Force suppresses SRX entry |

### 1.3 Super-Relaxed ADP (SRD) State Transitions

An ADP-bound variant of the SRX state (heads chelating ADP in the folded position). Governed by biochemical ADP/ATP competition and force.

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `ksrd` | PD → SRD (s⁻¹) | 1.0 | 0.1–20 s⁻¹ | ADP-bound head entering SRX-like state |
| `kmsrd` | SRD → PD (s⁻¹) | 1.07 | 0.1–20 s⁻¹ | Exit from ADP-SRX to available pool |
| `ksr2srd` | SR → SRD (s⁻¹) | 32.1 | 1–100 s⁻¹ | ADP rebinding to SRX heads |
| `ksrd2sr` | SRD → SR (s⁻¹) | 2.92 | 0.1–50 s⁻¹ | ADP release from SRD |
| `sigma_srd1` | Force sensitivity of SRD → PD (kPa⁻¹) | 25.3 | 5–100 kPa⁻¹ | Analogous to sigma1 for SRD branch |
| `sigma_srd2` | Force sensitivity of PD → SRD (kPa⁻¹) | 35.4 | 5–100 kPa⁻¹ | Analogous to sigma2 for SRD branch |

### 1.4 Stiffness Coefficients

Cross-bridge stiffness determines the force generated at a given strain. Values are in model force units per μm.

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `kstiff1` | Stiffness of P1 (weakly-attached) state | 16000 | 2000–30000 | Loosely attached heads; single XB ~0.2–1 pN/nm |
| `kstiff2` | Stiffness of P2 (strongly-attached) state | 97600 | 20000–150000 | Strongly attached; single XB ~1–3 pN/nm |
| `kstiff3` | Stiffness of P3 state (post-power-stroke) | 14000 | 5000–100000 | Off when UseKstiff3=0; defaults to kstiff2 |
| `kstiff1_n` | P1 stiffness at negative strain | 160 | 50–5000 | Separate negative-strain stiffness; heads resist reverse motion |
| `kstiff2_n` | P2 stiffness at negative strain | 704 | 100–5000 | Enables asymmetric force-strain relationship |
| `kSE` | Serial elastic element stiffness | 3203 | 500–20000 | Compliance of myofilaments + apparatus |
| `ekSE` | Exponent for nonlinear serial element | 0.925 | 0.5–2 | 1 = linear spring; <1 = softening |
| `k_pas` | Passive (titin) stiffness coefficient | 77.3 | 10–500 | Scales passive restoring force; titin isoform-dependent |
| `gamma` | Passive force exponent | 2.67 | 1–10 | Nonlinearity of titin force-extension; identified ≈ 7.9 |

### 1.5 Power Stroke & Geometry

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `dr` | Power stroke size (μm) | 0.010 | 0.005–0.015 μm | 5–15 nm; classic estimate ~10 nm for β-MyHC |
| `drp3` | Moment arm offset for P3 state (μm) | 0.010 | 0–0.02 μm | Strain reference shift for post-stroke state |
| `L_thick` | Thick filament half-length (μm) | 1.140 | 0.8–1.7 μm | From M-line to tip; cardiac ~0.8 μm (half-length from M-line) |
| `L_thin` | Thin filament length (μm) | 1.377 | 1.0–1.3 μm | From Z-disk; cardiac actin ~1.0–1.2 μm |
| `L_hbare` | Bare zone half-length (μm) | 0.10 | 0.08–0.12 μm | Region near M-line with no myosin heads |
| `SL0` | Initial sarcomere length (μm) | 2.20 | 1.8–2.4 μm | Optimal length for cardiac sarcomere ≈ 2.2–2.3 μm |
| `d_actin` | Actin site repeat distance (μm) | 0.0055 | 0.0055–0.0074 μm | Actin helical repeat ~5.5 nm; pairs at 36 nm |

### 1.6 Strain Grid

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `Slim_l` | Left (negative) strain boundary (μm) | 1.70 | 1.5–2.0 μm | — |
| `Slim_r` | Right (positive) strain boundary (μm) | 2.21 | 2.0–2.5 μm | — |
| `dS` | Strain grid step size (μm) | 0.002 | 0.001–0.01 μm | Finer = more accurate, slower; 2 nm resolution |

### 1.7 Strain Sensitivity of Transition Rates

Transition rates are modulated by cross-bridge strain. These parameters control the shape of the strain-rate function (exponential, piecewise, or Gaussian).

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `alpha1` | Exponential strain sensitivity of P1→P2 (μm⁻¹) | 100 | 10–500 μm⁻¹ | Higher = steeper rate increase with strain |
| `alpha_1` | Exponential strain sensitivity of P2→P1 (μm⁻¹) | 0 | 0–200 μm⁻¹ | Off here; reverse power stroke less strain-sensitive |
| `alpha2` | Symmetric strain sensitivity of P2→P3 (μm⁻¹) | 31.5 | 5–200 μm⁻¹ | Modulates detachment as function of strain |
| `alpha2_L` | P2→P3 rate sensitivity at negative strain (μm⁻¹) | 26.3 | 0–200 μm⁻¹ | Leftward (re-stretched) heads detach faster |
| `alpha2_R` | P2→P3 rate sensitivity at positive strain (μm⁻¹) | 1.0 | 0–200 μm⁻¹ | Positive-strain heads detach slowly (force-generating) |
| `alpha0` | Strain sensitivity of P1→PU detachment (μm⁻¹) | 37.5 | 0–200 μm⁻¹ | Currently overridden by piecewise function |
| `alphaRip` | Strain sensitivity of force-induced ripping (μm⁻¹) | 100 | 0–500 μm⁻¹ | Only active when k2rip > 0 |
| `StrainExp` | Exponent for generic strain-dep. rate formula | 2 | 1–4 | 2 = quadratic; 1 = linear strain sensitivity |
| `estiff` | Exponent for P2 force calculation (s+dr)^estiff | 1.001 | 0.5–2 | ~1 = nearly linear force-strain |
| `e2L` | Exponent for P2→P3 rate at negative strain | 2 | 1–4 | — |
| `e2R` | Exponent for P2→P3 rate at positive strain | 2 | 1–4 | — |
| `k2_L` | Base P2→P3 rate at negative strain (s⁻¹) | 200 | 0–2000 s⁻¹ | Fast detachment of negatively-strained heads |
| `k2_R` | Base P2→P3 rate at positive strain (s⁻¹) | 8000 | 0–20000 s⁻¹ | Rapid detachment at large positive strain (ripping) |
| `dr2` | Strain offset for P2→P3 (μm) | 0.027 | 0–0.05 μm | Shifts zero-crossing of strain-dep. detachment |
| `dr2_L` | P2→P3 strain offset at negative strain (μm) | 0 | −0.02–0.02 μm | — |
| `dr2_R` | P2→P3 strain offset at positive strain (μm) | 0.002 | 0–0.02 μm | — |
| `dr1` | Strain offset for P1→P2 rate (μm) | 0 | −0.02–0.02 μm | — |
| `dr_1` | Strain offset for P2→P1 rate (μm) | 0 | −0.02–0.02 μm | — |
| `dr0` | Strain offset for P1→PU detachment (μm) | 0 | −0.02–0.02 μm | — |
| `dr3` | Strain offset for force-induced ripping (μm) | −0.010 | −0.05–0 μm | — |

### 1.8 Piecewise Strain-Dependent Functions

When `UsePieceWiseStrainDep = 1`, the strain sensitivity of transitions is given by piecewise cubic Hermite interpolants rather than exponentials. The `*X` arrays are strain breakpoints (μm), the `*Params` arrays are rate multiplier values at those breakpoints.

| Parameter | Meaning | Current Value |
|-----------|---------|---------------|
| `PieceWiseStrainDepX` | Strain breakpoints for P1→P2 multiplier (μm) | [0.1, 0.0079, 0, −0.005, −1] |
| `PieceWiseStrainDepParams` | Rate multipliers at breakpoints (dimensionless) | [0.5, 1.155, 2.882, 58.8, 50] |
| `PieceWiseStrainDep2X` | Strain breakpoints for P2→P3 multiplier (μm) | [0.1, 0.026, 0.013, 0, −0.008, −0.008, −0.01] |
| `PieceWiseStrainDep2Params` | P2→P3 multipliers | [50, 2.91, 0.94, 1.04, 50.7, 50, 50] |
| `PieceWiseStrainDepR1DX` | Strain breakpoints for P1→PU detachment (μm) | [−0.05, −0.005, 0.006, 0.05] |
| `PieceWiseStrainDepR1DParams` | P1→PU detachment multipliers | [50, 1.15, 1.60, 50] |
| `PieceWiseStrainDepR21X` | Strain breakpoints for P2→P1 (μm) | [−0.058, −0.005, 0.006, 0.05] |
| `PieceWiseStrainDepR21Params` | P2→P1 multipliers | [50, 0.994, 0.930, 50] |

Physiological context: exponential strain sensitivity is a simplification; measured single-molecule and ensemble data suggest more complex, asymmetric shapes, hence the piecewise representation.

### 1.9 Viscosity & Damping

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `mu` | Viscous drag coefficient (positive velocity) | 0.633 | 0.01–5 | Internal sarcomere viscosity; phenomenological |
| `mu_neg` | Viscous drag at negative force/velocity | 0.633 | 0.01–5 | Asymmetric damping during re-stretch |
| `mu2` | Additional velocity-squared drag term | 1×10⁻⁹ | 0–0.1 | Essentially off |
| `kSE_M` | Maxwell dashpot spring stiffness | 100 | 10–10000 | Only active when UseMaxwellDashpot=1 |
| `eta_M` | Maxwell dashpot viscosity | 0.1 | 0.01–10 | Only active when UseMaxwellDashpot=1 |

### 1.10 A2 Attachment Shift (Actin Site Hopping)

After the power stroke, myosin can shift its strain by hopping to an adjacent actin site (spaced ~5.5 nm apart), redistributing the probability distribution in strain space.

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `slope` | Hopping rate slope (s⁻¹ μm⁻¹) | 16496 | 1000–50000 | Rate per unit distance of strain from threshold |
| `d_actin` | Actin site spacing (μm) | 0.0055 | 0.0055–0.0074 μm | 5.5 nm = first actin repeat |
| `s_threshold_R` | Positive strain threshold for hopping (μm) | 0.0046 | 0.001–0.015 μm | Heads strained beyond this hop to next site |

### 1.11 Calcium & Thin Filament (Troponin)

These parameters govern Ca-dependent thin filament activation. Note: the current model runs at saturating Ca (`Ca = 1000`), so these are not rate-limiting in the current configuration.

| Parameter | Meaning | Current Value | Plausible Range | Notes |
|-----------|---------|---------------|-----------------|-------|
| `Ca` | Free Ca²⁺ concentration (arbitrary units) | 1000 | 0–1000 | 1000 ≡ saturating; physiological systolic ~10 μM, diastolic ~0.1 μM |
| `K_coop` | Cooperative activation parameter | 5.7 | 1–20 | Nearest-neighbour cooperativity on thin filament |
| `k_on` | Troponin C Ca binding rate (μM⁻¹ s⁻¹) | 59.7 | 10–200 | — |
| `k_off` | Troponin C Ca off rate (s⁻¹) | 307 | 50–1000 | k_on/k_off = affinity |

---

## Table 2 — Mechanism Switches (`Use...` flags)

Current configuration: `ModelParamsSlackKtrOpt.m`.

### 2.1 Myosin Head State Mechanisms

| Switch | Current | Physiological Rationale | Notes |
|--------|---------|------------------------|-------|
| `UseSuperRelaxed` | **ON** | In relaxed muscle, ~50–70% of myosin heads adopt a folded (SRX/IHM) conformation with ~10× lower ATPase. This is the primary energy reserve in diastole; disrupted in hypertrophic cardiomyopathy (HCM). Force recruits heads back to the disordered pool. | Adds scalar state `P_SR` to the ODE |
| `UseSuperRelaxedADP` | **ON** | An ADP-bound SRX variant (SRD). When ADP levels rise (as in ischaemia/heart failure), ADP can stabilise the folded conformation. Provides an additional ATP-economy lever. | Adds state `P_SRD`; requires UseSuperRelaxed |
| `UseDirectSRXTransition` | **OFF** | Couples working-stroke detachment flux (P2→P3) directly to SRX recruitment. Physically motivated if heads detach directly into the folded state rather than passing through the unattached pool. | Adds a flux term; slightly changes energy balance |
| `UsePassiveForSR` | **OFF** | Gating of SRX ↔ PT transitions uses passive force only (not total force). Relevant if SRX recruitment responds to diastolic filling pressure rather than systolic active force. | Alternative to sigma1/sigma2 on total force |
| `useHalfActiveForSR` | **OFF** | Compromise: SRX gating driven by 0.5×active + passive. | — |
| `UseLimitedSRDTransition` | **OFF** | Caps the PD→SRD rate to prevent unphysical saturation at high ADP. | — |

### 2.2 Filament Geometry & Overlap

| Switch | Current | Physiological Rationale | Notes |
|--------|---------|------------------------|-------|
| `UseOverlap` | **ON** | At short sarcomere lengths (<~1.9 μm in cardiac), the zone of thick-thin filament overlap shrinks, reducing available actin sites and thus force. Frank-Starling mechanism partly operates through this. | Scales attachment rate by N_overlap(SL) |
| `UseOverlapFactor` | **OFF** | Alternative: normalise all states by available overlap fraction rather than using it as an absolute count. Changes whether force scales with available sites or with total probability. | — |

### 2.3 Passive Force (Titin)

| Switch | Current | Physiological Rationale | Notes |
|--------|---------|------------------------|-------|
| `UsePassive` | **ON** | Titin is the third filament of the sarcomere, providing passive restoring force. Dominates at SL > 2.2 μm; contributes to diastolic stiffness. Titin stiffness is increased in heart failure. | Adds k_pas·(SL−SL0)^gamma to total force |
| `UseTitinIdentifiedPassive` | **OFF** | Uses experimentally fitted titin parameters (nonlinear power-law, γ≈7.9). OFF = simpler linear model with current gamma=2.67. | — |
| `UseTitinInterpolation` | **OFF** | Reads titin force from a pre-computed lookup table (from a separate titin model). More accurate but requires external simulation. | — |

### 2.4 Elastic & Mechanical Elements

| Switch | Current | Physiological Rationale | Notes |
|--------|---------|------------------------|-------|
| `UseSerialStiffness` | **ON** | In vitro, the experimental apparatus (glass fibre, hook) introduces compliance in series with the preparation. In vivo, myofilament compliance and the Z-disk play this role. Without it, force-velocity and Ktr simulations are qualitatively wrong. | Adds ODE state `LSE`; stiffness = kSE |
| `UseMaxwellDashpot` | **OFF** | A viscoelastic (Maxwell) element in parallel with the sarcomere. Physiologically motivated by the complex time-dependent stiffness seen in rapid stretches (stiffness-frequency relationship). | Adds state `x_dash`; requires kSE_M, eta_M |

### 2.5 Cross-Bridge Strain Mechanics

| Switch | Current | Physiological Rationale | Notes |
|--------|---------|------------------------|-------|
| `UseNegativeKstiff` | **ON** | Pre-power-stroke and backward-strained heads have different elastic properties than post-power-stroke heads. Enables asymmetric force-extension: stiff in shortening, compliant in re-stretch. | Separate kstiff1_n, kstiff2_n for s < 0 |
| `UseWDetachment` | **OFF** | W-shaped detachment rate: fast at both very positive AND very negative strains. Physically: heads that are over-extended or compressed both detach rapidly. | Alternative to directional k2_L / k2_R |
| `UseNegativeForceRip` | **OFF** | Augments the P2→P3 (detachment) rate when the cross-bridge develops a force opposing shortening (negative force contribution). Related to XB behaviour during eccentric contractions. | — |
| `UseA2AttachmentShift` | **ON** | Post-power-stroke (P2) heads can hop to an adjacent actin monomer (~5.5 nm step), redistributing their strain. Reflects the spatial periodicity of actin binding sites and mechanical compliance of the working stroke. | Adds strain-redistribution flux to P2 |
| `UseA2Popping` | **OFF** | P2 heads that are being pulled to large negative strain (re-stretch) detach abruptly ("pop"). Relevant for eccentric contractions. | — |
| `UseA2Reattaching` | **OFF** | P2 heads can re-attach to a new actin site at a preferred strain. Incomplete implementation. | Not fully implemented |
| `UseA2MechanicalRecocking` | **OFF** | Mechanical force can rcock a P2 head backward to a pre-power-stroke conformation (P2D). Relevant for lengthening-induced sarcomere dysfunction. | — |
| `UseStrictDetachmentAt` | **OFF** | Forces detachment at extreme strains (|s| > threshold). Prevents unrealistic probability accumulation at grid edges. | Computational safeguard |

### 2.6 Attachment Distribution

| Switch | Current | Physiological Rationale | Notes |
|--------|---------|------------------------|-------|
| `UseA1AttachmentKernel` | **OFF** | Actin target zones are not point-like: they span several nanometres due to the helical arrangement of actin monomers. A Gaussian/triangular kernel distributes newly attached heads over a range of strains around zero. | Width set by A1AttachmentWidth = 0.006 μm |
| `UseSpaceInterpolation` | **ON** | When the strain distribution is advected (shifted by velocity), sub-grid velocities require interpolation rather than integer bin shifts. Increases accuracy at low velocities. | Computational accuracy switch |
| `UseSpaceExtension` | **ON** | If the probability distribution reaches a grid boundary (e.g. during rapid shortening), the strain grid is extended dynamically. Prevents probability loss or reflection artefacts. | Computational robustness switch |

### 2.7 Experimental Protocol Switches

| Switch | Current | Physiological Rationale | Notes |
|--------|---------|------------------------|-------|
| `UseSlack` | **ON** | Simulates the slack test: rapid release to zero force, then re-stretch. Used to extract Vmax (unloaded shortening velocity) and sarcomere dynamics during unloading. | Requires slack velocity table |
| `UseKtrProtocol` | **ON** | Simulates the Ktr (rate of force redevelopment) protocol: rapid shortening → ramp re-extension → isometric force recovery. Ktr reflects the rate of cross-bridge cycling. | — |
| `UsePieceWiseStrainDep` | **ON** | Use piecewise cubic Hermite interpolation for strain-dependent rates rather than simple exponentials. More biologically realistic: allows arbitrary strain-rate shapes. | See PieceWiseStrainDep* parameters |
| `UseStrainDep4R1D` | **OFF** | Use a Gaussian function (rather than piecewise) for the P1→PU (loosely-attached detachment) strain dependence. | Alternative to piecewise for R1D only |

---

*Generated 2026-03-02 from `ModelOptParams/ModelParamsSlackKtrOpt.m`.*
