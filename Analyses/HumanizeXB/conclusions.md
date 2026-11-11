# Humanizing the cross-bridge model — targets, what the model does, and a v0 parametrisation

2026-10-06. Evidence and per-paper critique in [`literature_log.md`](literature_log.md).
Scripts: `humanProbe.m` (species-agnostic read-out with windows long enough for β-MHC),
`applyScaling.m`, `RunHumanize.m` (resumable sweep driver), `RunExplore1.m` / `RunExplore2.m`
(the two sweeps), `WriteHumanSeed.m` (writes `params/human_v0_seed.m`). Results in
`results/explore1.*` (mechanistic sweep) and `results/explore2.*` (literature-ratio candidates).

## 1. Short answer

1. **Humanising this model is an isoform retune (α → β-MHC), not a species retune.** Human
   ventricle is ≥ 90–95 % β; mouse ~100 % α, rat 80–97 % α. The isoform changes *kinetics*
   (ADP release and hydrolysis ~10× slower, Pi release/power stroke ~3× slower, attachment rate
   scaling with them), **not mechanics** (step 7–10 nm, stiffness, force per head, specific tension
   are the same within the resolution of the data). Put the whole human difference into rates and
   leave `dr`, `kstiff*`, passive and series elements alone.
2. **The firm human numbers at matched temperature are ratios of ~0.1**: k_tr human/mouse 0.10–0.12
   (same lab, 15 °C), V0 0.10–0.15, tension cost 0.25, ATP turnover per head 0.1–0.2, ADP release
   ≤ 0.1. The model's base rates must come down by about an order of magnitude, with attachment and
   detachment moving together so that F0 and the attached fraction stay where they are.
3. **Temperature is the dominant uncertainty, not biology.** Human data are at 15 °C (mechanics) or
   20 °C (ATPase); the rig temperature of the rodent data is not recorded in this repo (assumed
   20–25 °C) and even the species is written both as mouse and rat in the abstract. Q10 for rates is
   3–3.5 (V0 4.6), so a 5 °C error is a ×1.8 error on every rate target. **Get the rig temperature and
   the species from the Baker lab before scoring any human fit.**
4. **The mant-ATP SRX numbers do not transfer and do not need to.** The slow-phase lifetime
   (140–450 s) is a resting, unloaded ATPase; the mechanically relevant OFF↔ON exchange is 1–90 s⁻¹
   (RLC probes, X-ray, LDA models), which is where the incumbent already sits. No defensible
   human-vs-rodent difference in SRX fraction exists (all β preparations 0.15–0.30 fitted,
   0.35–0.5 after this group's own correction; rodent α fibre data essentially absent). The
   model's SRX kinetics are left untouched in v0; scaling them with the cycle is a sensitivity arm.
5. **Proposed v0** (§4): `k2, k_2, kah, kamh, ka, kd ×0.1; k1, k_1 ×0.3`, everything else as
   `robust1_seed` (`params/human_v0_seed.m`). It gives F0 +8 % over the incumbent, k_tr 5.5 s⁻¹
   (×0.12), V0 ≈ 1.0 ML/s with a/F0 ≈ 0.05, ATP use 1.9 s⁻¹ per head (×0.11) — every kinetic
   observable at the literature's ~0.1 isoform ratio with force preserved, all driven by literature
   ratios rather than a fit. It is a *seed*, not a fit: there are no human slack–restretch or FV data
   in this repo to fit to, and the absolute human numbers still hinge on the unknown rig temperature.

## 2. Human targets, expressed at the model's (assumed) temperature

Human values from `literature_log.md` §§2–7, projected from 15 °C with Q10 3.5 (rates), 4.6 (V0),
1.3 (force), 2.5 (tension cost). Two columns because the rig temperature is unknown. "Ratio" is the
same-preparation human/rodent-α ratio at matched temperature — this is the number the parametrisation
should actually reproduce, independent of the Q10 guesswork.

| observable | human, 15 °C (measured) | human @ 20 °C | human @ 25 °C | ratio human/rodent-α | incumbent (α) | confidence |
|---|---|---|---|---|---|---|
| k_tr, max Ca | 0.55–1.3 s⁻¹ (use 0.8) | 1.5 | 2.8 | **0.10–0.12** | 45 | high |
| V0 | 0.4–0.6 ML/s (inferred, no direct datum) | 1.0–1.3 | 2–3 | 0.10–0.15 | > 3 (Hill ≈ 4) | medium-low |
| a/F0 | ≤ 0.1 (rat β 0.05, rat α 0.10 at 12 °C) | | | ≤ 0.5 | 0.63 (fitted to Baker FV) | low |
| peak normalised power | rat-α-like per mg at 12–15 °C | | | ~0.3–0.5 (power ∝ V0 × curvature) | 0.56 | low |
| ATP turnover per head, max Ca, isometric | 0.45–0.6 s⁻¹ (20 °C, assumes 0.15–0.2 mM heads) | 0.5 | 1.0 | 0.1–0.2 | 14 net / 19 gross | medium |
| tension cost | 2.3–2.4 µmol l⁻¹ s⁻¹ (kN m⁻²)⁻¹ (20 °C) | 2.3 | 3.6 | **0.25** | – | high |
| isometric detachment (slow kREL) | 0.29–0.30 s⁻¹ | 0.55 | 1.0 | ~0.1 | k2eff 80 | medium |
| F0, attached fraction | prep-dependent; **no isoform effect** | | | ~1 | 66 kPa, 0.28 | high |
| resting SRX fraction | 0.28 fitted / 0.35–0.5 corrected | | | ~1 | 0.27 | medium |
| ADP release (solution, 20 °C) | 60–93 s⁻¹ | | | ≤ 0.1 | – | high (direction) |
| hydrolysis (solution, 20 °C) | 14–17 s⁻¹ | | | 0.1 | – | high |
| Pi release / stroke (solution) | 12–20 s⁻¹ | | | 0.2–0.5 (α unmeasured) | – | low |

## 3. What the model does under rate scaling (`results/explore1.txt`)

Probe: isometric at SL 2.2 (3 s), Brenner-style k_tr (instant detachment, 2 s window), FV read
after 0.1 ML of shortening from the isometric state (same rule as `runFVExperiment`), active force.

| set | F0 (kPa) | attached | k_tr | V0 (zero crossing) | a/F0 | Pmax | XTOR gross / net | SRX |
|---|---|---|---|---|---|---|---|---|
| incumbent | 66 | 0.28 | 45 | > 3 (Hill 3.8) | 0.63 | 0.56 | 19.3 / 14.2 | 0.030 |
| all cycle rates ×0.25 | 68 | 0.29 | 13.1 | 2.9 | 0.10 | 0.15 | 3.6 / 2.3 | 0.021 |
| …+ SRX rates ×0.25 | 66 | 0.28 | 11.2 | 2.9 | 0.12 | 0.11 | 4.8 / 3.6 | 0.030 |
| …+ v_ref_reg ×0.25 | 68 | 0.29 | 13.1 | 2.9 | 0.23 | 0.19 | 3.6 / 2.3 | 0.021 |
| v_ref_reg ×0.25 only | 66 | 0.28 | 45 | > 3 | 1.6 | 0.65 | 19.3 / 14.2 | 0.030 |
| k2 ×0.1, rest ×0.25 | 93 | 0.41 | 4.6 | 1.0 | 0.055 | 0.07 | 3.0 / 1.8 | 0.006 |
| ka ×0.1, rest ×0.25 | 39 | 0.17 | 17.7 | 2.9 | 0.12 | 0.15 | 1.7 / 0.02 | 0.045 |
| kah ×0.1, rest ×0.25 | 66 | 0.28 | 14.9 | 2.9 | 0.11 | 0.15 | 1.9 / 1.4 | 0.029 |

What this says:

- **F0 is a ratio, as Brenner says.** Uniform slowing leaves F0 and the attached fraction untouched;
  slowing only detachment (k2) raises F0 +40 %, slowing only attachment (ka) cuts it −40 %. Since
  specific tension is not an isoform property, **ka and k2 must scale together** in the human set.
- **k_tr scales sub-linearly** (×0.25 on all rates → ×0.29 on k_tr) because the Brenner rise is
  limited by the serial element and the registration gate as well as by f + g.
- **The zero-force velocity is not kinetic in this model.** Any uniform slowdown leaves the zero
  crossing at 2.9 ML/s; only a specific, large k2 drop moves it (×0.1 → 1.0 ML/s). What *does*
  respond to uniform slowing is the FV curvature (a/F0 0.63 → 0.10) and the power (×0.26).
  Mechanistically the high-velocity tail is carried by a small positive-strain population whose
  lifetime is set by the strain-dependent (not the base) part of the detachment curve and by the
  registration gate; that is the "V_max is geometric" finding of `atp-sensitivity-mechanism` seen
  from the other side. For a human set the right V0 therefore comes out of the k2 ratio (which the
  literature puts at ≤ 0.1 anyway), not out of a global multiplier.
- **The registration gate has a fixed velocity scale** (`v_ref_reg` = 1.02 ML/s). Scaling it with
  the kinetics (×0.25) makes the FV *less* curved (a/F0 0.10 → 0.23) by letting attachment recover
  at low velocity; leaving it fixed makes the gate nearly inactive at human velocities. Neither is
  derived from human data; v0 leaves it fixed (gives the more human-like a/F0 ≤ 0.1) and flags it.
- **The SRX rates do little at full activation** (3 % pool), ±15 % on k_tr; they matter for the
  rest state and the restretch transients, where there are no human data.

## 4. Literature-ratio candidates (`results/explore2.txt`)

ATP use per head is the P2 → PT exit flux `R2T` (every post-stroke detachment binds a fresh ATP);
`RTD − RDT` is *not* it, because the SRX path also hydrolyses (`project_xtor_is_gross`). `R2T` was
added to the probe after the first four sets ran; for those it equals the v0 value to within 5 %
(same k2, same p2). Incumbent R2T ≈ 17 s⁻¹ (from `LowATP_Refit9K2Sweep`).

| set | scaling | F0 | attached | k_tr | t50 (ms) | V0 cross / Hill | a/F0 | Pmax | ATP use R2T | SRX | dwell (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| incumbent | – | 66 | 0.28 | 45 | 15 | > 3 / 3.8 | 0.63 | 0.56 | ≈ 17 | 0.030 | 15 |
| lit_A | k2,k_2,kah,kamh ×0.1; k1,k_1 ×0.3; ka ×0.1; kd ×0.25 | 60 | 0.26 | 7.9 | 65 | 1.38 / 1.6 | 0.072 | 0.069 | ≈ 1.6 | 0.017 | 336 |
| lit_A, ka ×0.15 | as above, ka ×0.15 | 74 | 0.32 | 6.4 | 62 | 1.0 / 1.8 | 0.068 | 0.071 | ≈ 1.9 | 0.013 | 274 |
| lit_B | lit_A + v_ref_reg ×0.1 | 60 | 0.26 | 7.9 | 65 | 1.42 / 1.0 | 0.30 | 0.108 | ≈ 1.6 | 0.017 | 336 |
| lit_C | lit_B + SRX rates ×0.1 | 56 | 0.25 | 5.6 | 124 | 1.0 / 0.8 | 0.24 | 0.071 | – | 0.042 | 143 |
| lit_A, k1,k_1 ×0.2 | | 60 | 0.26 | 7.8 | 66 | 1.0 / 1.5 | 0.077 | 0.067 | 1.62 | 0.017 | 344 |
| **v0 = lit_A, kd ×0.1** | k2,k_2,kah,kamh,ka,kd ×0.1; k1,k_1 ×0.3 | **71** | **0.30** | **5.5** | 74 | **1.0 / 1.8** | **0.054** | 0.061 | **1.93** | 0.014 | 279 |
| v0, k1,k_1 ×0.2 | | 71 | 0.30 | 5.5 | 75 | 1.0 / 1.6 | 0.059 | 0.059 | 1.90 | 0.014 | 284 |
| uniform ×0.1 (+ v_ref_reg) | all cycle rates ×0.1 | 70 | 0.30 | 5.4 | 78 | 1.0 / 0.9 | 0.23 | 0.083 | 1.82 | 0.014 | 299 |
| uniform ×0.05 (+ v_ref_reg) | all cycle rates ×0.05 | 71 | 0.31 | 2.7 | 130 | 0.5 / 0.4 | 0.25 | 0.043 | 0.93 | 0.009 | 910 |

Against the targets of §2:

- **F0 / attached fraction**: the ka = k2 pairing holds F0 within ±10 % in every set, as Brenner
  predicts. ka ×0.15 overshoots F0 for no human reason; keep ka = 0.1. Taking kd from 0.25 to 0.1
  raises F0 +8 % (fewer weak heads lost before the stroke) — acceptable, and physically motivated:
  β binds actin ~5× tighter than α (Deacon 2012), so the weak-binding off-rate should fall at least
  as much as the on-rate.
- **k_tr**: lit_A's 7.9 s⁻¹ is 0.18× the incumbent; the isoform ratio is 0.10–0.12. **k1 is not the
  lever** (×0.3 → ×0.2 leaves k_tr at 7.8): the Brenner rise is limited by how many weakly bound
  heads survive to stroke (k1/(k1+kd)) and by the serial/registration elements. **kd is**: kd ×0.1
  gives 5.5 s⁻¹ = ×0.12, on the isoform ratio. Uniform ×0.1 gives the same k_tr (5.4) but a far less
  curved FV, and uniform ×0.05 (the "rig at 25 °C" case) gives 2.7 s⁻¹.
- **V0 1.0 ML/s (zero crossing; Hill 1.6–1.8), a/F0 0.05, Pmax 0.06 (×0.11)**: all consistent with
  human β projected to 20–25 °C (V0 1–3 ML/s, a/F0 ≤ 0.1, power ×0.1–0.3). The zero crossing is set
  by the k2 ×0.1 step and is the same in every k2 ×0.1 set.
- **ATP use**: v0 1.9 s⁻¹ per head = ×0.11 of the incumbent's ≈ 17. In absolute terms the human
  target is 0.5–1 s⁻¹ per head at 20–25 °C (tension cost 2.3 µmol l⁻¹ s⁻¹ per kN m⁻² at 40 kPa
  and 0.17 mM heads), so v0 is ~2× high — exactly the incumbent's own 2–3× excess over mouse-α
  (`AttachedPool_vs_ATPase`), inherited through the ratio. The ratio is right; the offset is the
  rodent fit's.
- **Dwell time** 280 ms (v0): human isometric detachment from slow kREL 0.3 s⁻¹ at 15 °C projects to
  0.55–1 s⁻¹ at 20–25 °C (dwell 1–1.8 s), but kREL is a relaxation-phase number (deactivating thin
  filament), so 0.3 s is the right order.
- **lit_C** (SRX rates scaled) also reaches k_tr 5.6 but through a slow SRX refill (t50 124 ms) and a
  larger active SRX pool; no human evidence supports it. Keep as a sensitivity arm.
- **v_ref_reg scaling** (lit_B, uniform sets) makes the FV markedly less curved (a/F0 0.23–0.30 vs
  0.05–0.08). Rat β a/P0 = 0.05 argues against scaling it; v0 leaves it fixed.

### 4.1 The v0 parametrisation

Applied to `params/robust1_seed.m` (all other fields unchanged); written to
`params/human_v0_seed.m` by `WriteHumanSeed.m`.

| parameter | incumbent | ×factor | v0 | basis (literature_log §) |
|---|---|---|---|---|
| k2 (post-stroke detachment / ADP-release isomerisation) | 117 | 0.1 | 11.7 | β 60–93 vs α > 1000 (≥ 350) s⁻¹ at 20 °C; Little 2012 fibre 41 vs 159 (§6, §9) |
| k_2 | 4.36 | 0.1 | 0.44 | tied to k2 (ADP rebinding scales with K_AD too, but ADP = 0 here) |
| kah (hydrolysis) | 176 | 0.1 | 17.6 | 17 vs 150–168 s⁻¹ (Deacon 2012) — the one number with a direct human/mouse pair (§6) |
| kamh | 8.8 | 0.1 | 0.88 | keeps K_hyd; no isoform datum |
| k1, k_1 (power stroke / Pi release) | 161, 31 | 0.3 | 48, 9.4 | β Pi release 12–20 s⁻¹; α unmeasured, "3–5× faster" assumed (§6); k_tr is insensitive to it here (§4) |
| ka (attachment) | 292 | 0.1 | 29 | not measurable; fitted human/mouse 60/330 (Mijailovich 2022); must track k2 to hold F0 and duty ratio (§3) |
| kd (pre-stroke detachment) | 93 | 0.1 | 9.3 | β actin affinity ~5× tighter than α (Deacon 2012) → off-rate falls at least as much as ka; Little 2012's fibre 0.26 is a post-stroke number. kd is the k_tr lever that puts the set on the 0.1 isoform ratio (§4) |
| SRX / OFF rates, σ's | – | 1 | – | no species evidence; mechanosensing is a thick-filament property (§7) |
| registration gate v_ref_reg, tau_reg; load-realign k_cr, tau_cr, F_cr | – | 1 | – | no human data; flagged (§3) |
| dr, kstiff1/2, kstiff*_n, estiff | – | 1 | – | step/stiffness isoform-independent (§1, §4 of log) |
| passive (k_pas, γ, Lsc0, Maxwell), kSE, ekSE | – | 1 | – | human titin is 30 % N2BA (softer) — a later, separate identification |
| MgATP, MgADP, Pi | 8, 0, 0 | – | – | human in-vivo: ATP 5.7 (DCM 3.7) mM, ADP 50–90 µM, Pi 1–1.4 mM (§8.2) |

What v0 is *not*: it is not temperature-corrected (the rig temperature is unknown), not fitted to
any human trace, and does not touch the five things the model does not represent (thin-filament
activation, cMyBP-C, titin isoform, temperature, ADP-release inhibition). It is the parameter set to
start a human campaign from once human protocol data exist.

## 5. Limits found on the way

- **The incumbent's own turnover is 2–3× too high for mouse-α** (R2T ≈ 17 s⁻¹ vs 4–8; see
  `AttachedPool_vs_ATPase`, `project_xtor_is_gross`). The ratio carries that offset into v0
  (1.9 vs a human 0.5–1 s⁻¹ per head). A human campaign should carry an explicit R2T bound
  (0.3–1.5 s⁻¹ per head at 20–25 °C) rather than inherit the rodent excess.
- **V0 is not a kinetic read-out of this model.** Only the k2 step moves the zero crossing; a
  human FV fit will therefore pin k2, which is also the ATP lever. That is convenient for the ATP
  question but means V0 and the 2 mM behaviour are not independent constraints.
- **The registration gate's velocity scale (1 ML/s) is a rodent number** baked into the FV shape.
  Whether it should scale with kinetics (lit_B, a/F0 0.30) or not (lit_A, 0.07) can only be
  decided with a human FV curve; the available evidence (rat β a/F0 0.05) favours not scaling it.
- **No human slack–restretch data exist.** The restretch-transient features that dominate the
  rodent objective (peak1, valley, peak2, rsK80) have no human target; the one human
  stretch-activation number (van Dijk 2014 delayed-rise k2 = 0.74 s⁻¹ at 15 °C ≈ k_tr) is a
  consistency check only.
- **Human β's ADP affinity (K_AD 6–21 µM)** puts physiological ADP (50–90 µM) at 3–10 × K_AD; the
  model's ADP-release inhibition path is dead (`atp-sensitivity-mechanism`), so a human ATP/ADP
  scenario needs that path wired before it means anything.

## 6. Next steps (in order of value per hour)

1. **Ask the Baker lab**: rig temperature, species (mouse or rat), and whether a human (or at least
   β-MHC: pig, or PTU-rat as a stepping stone) preparation can be run on the same slack–restretch
   and FV protocols. Even one human prep at two temperatures would turn v0 from a rescaling into a fit.
2. Run v0 through the *existing* 8 mM battery at the rodent protocol timing to see which feature
   readers break (most windows are sized for k_tr ≈ 45; a human k_tr of 5 needs 5–10× longer holds),
   then build a `protocol_human_*` velocity table with holds ≥ 1 s.
3. Add a net-turnover output bound and a V0 bound to the human objective; drop the restretch
   features until human restretch data exist.
4. Wire the ADP-release inhibition (`k2_eff = k2·[ATP]/([ATP] + K_T·(1+[ADP]/K_AD))` with
   K_AD ≈ 10–20 µM for human) before any human HF scenario; then run ATP 5.7 → 3.7 mM,
   ADP 50–90 µM, Pi 1 → 3 mM — not 2 mM ATP.
5. Temperature: either record the rig temperature and keep the model there, or add a Q10 block
   (3–3.5 on cycle rates, 1.3 on force) so human 15 °C myofibril and 37 °C intact data can both
   be scored; the human intact k_tr of 18 s⁻¹ at 37 °C is the ceiling such a block must respect.
