# The Force-Velocity Shoulder / Relative Isometric Force Reduction

**Why the model's FV is too steep near v=0, what mechanism the data implies, and how to add it without breaking ktr or the restretch.**

Repo: `ATP-depletion-and-heart-failure` · Model: `dPUdT_CombinedTransitions` · Basis: BoundedFit (`params/params_reseeded.m`, p1-force config) · Date: 2026-05

---

## 1. The phenomenon

The 8 mM force-velocity data has a **flat shoulder** near isometric:

| v (ML/s) | 0 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| **data** F/F(0) | 1 | **0.919** | 0.664 | 0.316 | 0.111 |
| **model** F/F(0) | 1 | **0.51–0.65** | 0.32 | 0.15 | 0.06 |

Force falls only ~8% at v=0.5 ML/s in the data, but ~40–50% in the model. Equivalently:
the **model's isometric force is *not* depressed** relative to the cross-bridge force capacity,
whereas the **data's isometric is depressed** — the same depression that shows up as
`FV-isometric (≈56 kPa) < slack-isometric (≈66 kPa)` at the *same* SL≈2.0 in the data.

This is **not a Hill hyperbola** (Hill is steepest at v=0). A flat-shoulder/near-isometric
plateau is the classic **high-force deviation / "double-hyperbolic" FV** (Edman 1988; Edman,
Mulieri & Scubon-Mulieri 1976): above ~0.8 P0 the shortening velocity is lower than Hill
predicts, i.e. in F-vs-v space the curve sits *above* Hill near v=0. So the shoulder is a
documented, real feature — the model is missing the mechanism that produces it.

## 2. Why the model can't make it (established this campaign)

`FV_fnorm @v=0.5` is stuck **0.5–0.65 across every configuration tried** — invariant to
`ka`, `kstiff`, `kSE`, `mu`, occupancy `P_bound_max`, the R2 strain knots, and even the new
`dr1` (p1-force). Raising `ka` 50→250 (hypothesis: faster reattachment maintains shortening
force) left `FV_fnorm` at 0.51, only doubling force. Force-scale knobs are scale-invariant
under normalization; the kinetic knobs (R2 detachment shape) move Vmax and the high-v tail
but not the low-v shoulder. **The shoulder is structurally absent**: in the current model the
isometric attached fraction is *maximal* (all available heads that can bind, bind), so there
is no headroom for force to be *maintained* (rather than fall) when shortening begins.

## 3. The mechanism the data implies — velocity/registration-dependent availability

The user's physiological hypothesis (correct, and standard in the spatially-explicit
literature): **the myosin-head axial period (~14.3/43 nm) is mismatched with the actin
target-zone period (~36 nm helical pseudo-repeat).** At a *fixed* (isometric) register only a
fraction of heads face an available target zone → reduced attachment → **lower isometric
force**. During **sliding (shortening)** the register continuously refreshes, so more heads
get a binding opportunity per unit time → **higher availability** → force is *maintained* as
shortening begins (the shoulder), until the velocity-induced detachment eventually wins and
the curve falls. This is:
- **Tanner, Daniel & Regnier 2007** (PLoS CB) — sarcomere lattice geometry governs cooperative
  myosin binding; binding is target-zone-limited.
- **Daniel, Trimble & Chase 1998** (Biophys J) — *compliant realignment of binding sites*:
  filament compliance + sliding realign heads with sites, tuning binding.
- **Mijailovich et al. 2016** (MUSICO, JGP) — explicit 3-D target zones; mass-action models
  "neglect discrete positions and site-number constraints."
- **Caremani et al. 2015** (J Physiol) — the number of attached motors is *not* maximal at
  isometric and changes through shortening; consistent with sub-maximal isometric occupancy.

So the "relative isometric reduction" = **sub-maximal isometric attachment because of limited
target-zone registration, relieved by sliding.** It is the same idea as the model's existing
`UseVernierVelocity` and `UseVelGaussAttachment` flags — but those are mis-implemented (next).

## 4. Why the existing velocity-gated availability breaks ktr and the restretch

Both flags gate attachment on the **instantaneous half-sarcomere velocity** `velHS`:
- `UseVelGaussAttachment` (dPUdT:308): `ka_eff *= 1 − amp·exp(−(|velHS|−c)²/2σ²)` (dip at v=0).
- `UseVernierVelocity` (dPUdT:564): `f_sat = α + (1−α)/(1 + v_ref/|velHS|)` (sigmoid; α at v=0 → 1 at high v).

**`velHS = (Force − F_total)/mu`** (dPUdT:236,258) — it is the *internal* shortening velocity
set by the series-elastic imbalance, **not** the imposed velocity. Two consequences:
1. **During the ktr redevelopment** (end clamped, "isometric"), force is rising, so
   `Force ≠ F_total` and `velHS ≠ 0` *transiently*. The gate fires erratically as the SE
   equilibrates → the redevelopment becomes multi-phase/oscillatory (observed: ktr→150,
   rmse 5.7). The single-order redevelopment is destroyed by the SE↔gate interaction.
2. **During the restretch** (fast stretch, large +velHS) attachment is transiently *boosted*
   → a huge peak1 overshoot.

The defect is **using an instantaneous kinematic signal for what is physically a slow,
history-dependent registration state.** Velocity tracks registration only at *steady state*;
during transients (the SE-coupled redevelopment, the brief restretch) they diverge.

## 5. Proposed mechanism — a slow "registration availability" STATE variable

Replace instantaneous `f(velHS)` with a state `A_reg ∈ [A0, 1]` that *integrates* sliding:

```
dA_reg/dt = ( A_inf(|v_slide|) − A_reg ) / tau_reg
ka_eff   = ka · A_reg
A_inf(u) = A0 + (1 − A0) · u/(u + v_ref)      % A0 (~0.6–0.8) at u=0, → 1 at high sliding
```
where `v_slide` is the *filament* sliding speed (use the imposed `vel`/`Velocity`, or a
low-pass of `velHS`, **not** raw `velHS`), `A0` is the residual isometric availability, and
`tau_reg` is the registration relaxation time. Behaviour:
- **Sustained isometric** → `A_reg → A0 < 1` → **depressed isometric force** (the shoulder's
  reference) and recruits from the unexploited DRX/PT reservoir only partially.
- **Right after a slack** (sliding charged `A_reg→1`) → redevelopment starts with *high*
  availability → **fast, single-order ktr**; `A_reg` then relaxes to `A0` over `tau_reg`.
- **During shortening** → `A_reg → A_inf(v) > A0` → force *maintained* → **flat shoulder**.
- **During the brief restretch** → `A_reg` (slow) barely moves → **no attachment spike**.

**Tuning `tau_reg` is the crux** (and answers "single-order redevelopment"): if `A_reg`
decays much *slower* than the ktr rise (~20 ms), the redevelopment stays single-order and the
isometric depression appears only as a small, slow secondary settling — keep the depression
modest (`A0 ≈ 0.7–0.85`, i.e. 15–30%, matching data FV-iso/slack-iso ≈ 56/66 = 0.85).
"Maybe the lower availability is just narrower" (user) maps to **small `v_ref`** (availability
recovers within v<0.5 → full shoulder) combined with **`tau_reg` slow enough to protect the
transients**. The two knobs (`A0` depth, `v_ref` width, `tau_reg` time) are separable.

### 5b. Implementation spec (design answers)
- **Size: +1 SCALAR ODE state** (`A_reg`), strain-INsensitive. Registration availability is a
  global property of the *detached* head pool seeking target zones; it multiplies the total
  attachment rate `ka_eff`, it is NOT a function of post-attachment strain. A strain-resolved
  version (+30 states) is wasteful and physically wrong.
- **Driver: imposed filament `vel`, NOT `velHS`** — the slack itself charges `A_reg` (large
  `vel`), which then decays during the hold; using imposed `vel` sidesteps the SE coupling that
  corrupted the instantaneous gates.
- **`tau_reg ≈ 1/ktr` (~15–25 ms) — the constraint AND the physiology coincide.** Too slow
  (`≫1/ktr`) → two-phase redevelopment (fast rise then slow droop), non-single-order; too fast
  (`≪1/ktr`) → collapses to the broken instantaneous gate. At ~1/ktr the redevelopment stays
  monotonic/single-order and the depression shows only as a lower plateau. This is also the
  physiological rate: availability turns over *via* cross-bridge cycling (heads attach/detach
  and re-sample target zones at the cycling rate).

```matlab
% +1 ODE state A_reg in [A0,1]; UseRegistrationAvailability (default off via A0=1)
A_inf  = A0 + (1-A0) * abs(vel)/(abs(vel) + v_ref_reg);   % A0 at isometric -> 1 at high slide
dA_reg = (A_inf - A_reg) / tau_reg;
ka_eff = ka_eff * A_reg;                                  % multiplies attachment only
% defaults: A0 ~ 0.7-0.85 (depth; data FV-iso/slack-iso~0.85), v_ref_reg ~ 0.3-0.5 ML/s
%           (narrow; recovers by v~0.5), tau_reg ~ 0.015-0.025 s (~1/ktr)
% backward-compatible: A0=1 => A_reg≡1 => identical to current model.
```

## 6. The slack-feature constraint (must not regress, or be retunable)

The state-variable design is specifically chosen so the slack features stay tunable:
- **ktr**: protected because `A_reg` is *high* (freshly charged by the slack) during the
  redevelopment — opposite of the velocity-gate which *suppressed* it. The slow `tau_reg`
  also decouples from the SE.
- **peak1/valley (restretch transient)**: protected because `A_reg` doesn't spike on the fast
  stretch; any residual mismatch is retunable with the restretch knobs that don't touch FV
  (catch-bond `k_catch_bond`, viscoelastic SE `c_SE_visc`, `s_threshold_R`/`slope`).
- **A-vs-segment (force-length)**: unaffected — that is the *radial* lattice spacing
  (`UseLatticeSpacing`, d10 vs SL), a different axis from the *axial* registration here.

## 7. Alternative / interacting mechanisms (ranked)

1. **Series compliance "internal shortening" — RULED OUT (quantitative).** The SE relaxation
   rate is `kSE/mu = 3203/0.015 ≈ 2e5 s^-1` → time constant ~5 µs. The FV ramp lasts 0.2 s, so
   the SE is quasi-steady throughout ⇒ `velHS = vel` at the readout ⇒ no compliance-induced
   velocity lag ⇒ kSE does NOT shape the shoulder (it only matters for the sub-ms restretch).
   Also kSE is passive-identified (structural). Path cancelled.
2. **Edman high-force "give"** (Edman 1988): a strain-dependent forcible detachment near P0
   that truncates the high-force end. In this model that is the R2 strain curve at the
   *isometric* working strain — but campaign sweeps show it trades into ktr (it changes duty).
   The registration state is cleaner because it acts on *attachment*, not detachment.
3. **Distributed attachment kernel — RULED OUT.** Attachment is ALREADY distributed over a
   triangular ~0.006 µm window (`A1AttachmentWidth≠0` takes the `else` branch, dPUdT:587+;
   `UseA1AttachmentKernel` only switches Gaussian↔triangular). This is a numerics/anti-oscillation
   device — it spreads *where* heads attach (strain), not *how many* or any velocity dependence,
   so it cannot produce the shoulder. (The "compliant realignment" concept is real but is the
   registration STATE of §5, not this kernel.) Path cancelled.
4. **Two-headed / second-head binding during shortening** (Brunello 2007; Caremani) — extra
   attached heads recruited transiently. Out of class for a 1-head-per-site mean field.

## 8. The low-ATP fading

User observation: the shoulder fades at low ATP. Interpretation: at low [MgATP] the
**detachment step becomes rate-limiting** (MgATP-dependent cross-bridge release; the model's
ATP-dependent detachment is currently *commented out* — see `low-ATP-force-enhancement-
analysis`). When detachment is the bottleneck, the attached fraction is high and
**attachment-availability is no longer rate-limiting**, so the registration/availability
modulation has little leverage → the shoulder flattens out / disappears. This is consistent
and *predictive*: the registration mechanism should be placed on **attachment** (so it
naturally fades when detachment dominates at low ATP), and it should be tested across 8 vs
2 mM once ATP-dependent detachment is re-enabled. If instead the shoulder were a detachment
or SE effect, it would *not* fade with ATP in this way — so the ATP behaviour is itself a
discriminating test between the candidate mechanisms.

## 9. Recommendations (do AFTER the running optimisation; no core-model edits now)

1. **Protocol (confirmed, read-only) — both FV points are at SL 2.0**: `runFVExperiment`
   integrates each shortening point for `t = 0.1/|v|` s from SL0=2.2; with `Vums = v·ML`,
   `ML=2` (evaluateModel:158), dSL/dt = 2·v, so SL displaces `2×0.1 = 0.2 µm` ⇒ ramp ends at
   **SL 2.0** at every velocity — the SAME SL as the isometric point (cold steady run at SL 2.0).
   So there is **NO SL/overlap confound**: the shoulder is a pure velocity/availability effect.
   The SE is quasi-steady on this timescale (§7.1) so it is not a quick-release-compliance
   artifact either. (The data FV-iso(56)/slack-iso(66) gap may still be a prep/history effect
   between the two experimental protocols — worth confirming against the bench FV protocol.)
2. **Lowest-risk first**: a controlled `kSE` sweep (compliance) and `UseA1AttachmentKernel`
   (compliant realignment) — both reuse existing machinery, no new state.
3. **If those are insufficient**: implement the **registration-availability state** (§5) as a
   new optional flag `UseRegistrationAvailability` with params `A0`, `v_ref_reg`, `tau_reg`,
   gating `ka_eff` (and an extra ODE state `A_reg`). Default off → fully backward-compatible.
   Drive it by imposed `vel` (or low-passed `velHS`), *never* raw `velHS`.
4. **Validate** against: FV shoulder (8 mM), clean single-order ktr, restretch peak1/valley
   retunable, A-vs-segment unchanged, and the **8-vs-2 mM shoulder fade** (mechanism test).

## 10. Verdict

The flat-shoulder FV is the **target-zone registration / availability** effect — sub-maximal
isometric attachment relieved by sliding (Edman 1988; Tanner/Daniel/Regnier 2007; Daniel/Chase
1998; MUSICO; Caremani 2015). The model already encodes the *idea* (`UseVernierVelocity`,
`UseVelGaussAttachment`) but mis-implements it as an **instantaneous `velHS` gate**, which the
series elasticity corrupts into multi-phase ktr and a restretch spike. The fix is a **slow
registration-availability state** (modest depth `A0`, narrow `v_ref`, protective `tau_reg`)
on the *attachment* path — which produces the shoulder, preserves single-order ktr (charged by
the preceding slack), leaves the restretch retunable, and *predicts* the low-ATP fade.

## References
Edman 1988 (J Physiol 404:301) · Edman, Mulieri & Scubon-Mulieri 1976 (Acta Physiol Scand 98:143)
· Hill 1938 (Proc R Soc B 126:136) · Huxley 1957 (Prog Biophys 7:255) · Huxley & Simmons 1971
(Nature 233:533) · Piazzesi & Lombardi 1995 (Biophys J 68:1966) · Daniel, Trimble & Chase 1998
(Biophys J 74:1611) · Tanner, Daniel & Regnier 2007 (PLoS Comput Biol 3:e115) · Williams, Regnier
& Daniel 2010/2012 (PLoS Comput Biol) · Mijailovich et al. 2016 (MUSICO, J Gen Physiol 148:459) ·
Caremani et al. 2015 (J Physiol 593:3313) · Brunello et al. 2007 (PNAS 104:20114).
