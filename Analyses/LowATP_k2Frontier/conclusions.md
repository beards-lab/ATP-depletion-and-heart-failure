# Low-ATP first frontier: can reducing k2 reproduce the 2 mM signature?

**Question.** ATP binding drives cross-bridge detachment, so ATP limitation should act
first on the detachment step `k2` (ADP-release / g_app). Can a *single* reduction of
`k2` reproduce the measured 2 mM-vs-8 mM cross-bridge signature?

**Method.** Current best (`params/params_reseeded_regavail_opt2.m`, feature cost ~8.86)
run on the **new 8 mM slack protocol** (`data/protocol_03_27_2026_8mM_slack.mat`), with
`k2` scaled down from ×1.0 to ×0.3. Scoring is **relative-ratio**: the model's ATP
response `model(reduced k2)/model(base k2)` is compared to the measured data ratio
`features(2mM)/features(8mM)` — so the imperfect 8 mM baseline cancels. Rundown ignored
(files as-is). Broad feature target. Script: `RunK2Frontier.m`; data: `results/k2_frontier.mat`;
figure: `results/k2_frontier.png`.

## Result — the response is real but **bimodal**: no single k2 satisfies both clusters

| feature | data ratio (2/8 mM) | k2 scale that matches |
|---|---|---|
| steady | 1.18 | **0.67** |
| A | 1.18 | **0.69** |
| Am | 1.18 | **0.71** |
| peak1_y | 1.25 | **0.65** |
| t0 (onset) | 1.31 | **0.65** |
| peak2 | 1.18 | **0.40** |
| ktr | 0.54 | **0.39** |
| vall_y | 1.20 | — (model flat ≈1.01) |
| restretchSlopeStart | 1.26 | — (max ≈1.16 at ×0.3) |
| ktr2_overshoot | 5.56 | ≈0.30 (noisy) |

Two camps, ~1.7× apart in k2:

- **Isometric / onset cluster** (steady, A, Am, peak1, t0) → matched at **k2 ≈ ×0.68**.
- **Kinetic cluster** (ktr, peak2) → matched at **k2 ≈ ×0.40**.

The best single compromise (core-6 log-ratio cost) is **k2 ≈ ×0.50 (k2 150→75)**, which
leaves **+11–15 % residuals everywhere** (force features overshoot; ktr under-slows to
×0.62 vs target ×0.54).

## Interpretation

1. **k2 IS the right lever for the kinetics.** At ×0.40, `ktr` drops 65.7→35.8 (×0.55,
   nailing the data's ×0.54) and `peak2` rises +18 % — both spot-on. Direction of every
   feature is correct (force↑, ktr↓, onset↓, redevelopment more oscillatory). ATP-limited
   detachment cleanly explains the *dynamic* signature.

2. **The model's isometric force is ~2× too k2-sensitive.** For the ktr slowdown that the
   data shows (×0.54), `steady`/`A` jump **+41 %** where the data rises only **+18 %**.
   Equivalently: the same k2 that fixes `ktr` blows `steady` past target. This is the same
   defect as the baseline **peak2-overshoot** residual — at baseline the model's `steady`
   (75) sits *below* `peak2` (92), whereas the data has `steady ≈ peak2` in both conditions;
   lowering k2 lets `steady` race up to and past `peak2`. The duty→isometric-force coupling
   is too stiff (no force/occupancy ceiling biting here).

3. **Two features are k2-inert and need a different mechanism.** `vall_y` (post-restretch
   valley) is pinned ≈63 kPa at all scales (data +20 %), and `restretchSlopeStart` under-
   responds. These are set by stretch/restretch detachment dynamics, not steady-state k2.

## Conclusion & next levers

**A single k2 reduction is necessary but not sufficient.** It reproduces the kinetic ATP
signature (ktr, peak2, overshoot, onset) at ×0.40 but over-produces isometric force gain
by ~2×. The gap is specifically the **isometric force-per-duty sensitivity**. Candidate
second levers, in priority order:

1. **Cap isometric force rise** — engage occupancy/force saturation (`UseGlobalOccupancy
   Saturation`, `P_bound_max`) so `steady` plateaus as the attached pool grows, flattening
   the force ratio toward +18 % while keeping the ktr slowdown.
2. **Co-recruit SRX** — let low ATP also shift heads *into* SRX (fewer available), offsetting
   part of the attachment gain. Physiologically motivated (heads can't refold IHM without
   ATP rebinding; see `../LowATP_ForceEnhancement/conclusions.md`).
3. **Add an ADP-trap / 3rd slow state** for `vall_y` and the restretch stiffness, which no
   k2 value reaches.

The clean takeaway for the ATP campaign: **k2 owns the kinetics; a force-saturation (or
SRX) lever owns the amplitude.** Splitting the ATP effect across `k2` + one amplitude lever
is the minimal two-knob model.

## Second lever tested: k2 × occupancy cap (`RunK2OccupancyGrid.m`)

Paired the k2 cut with the global occupancy-saturation ceiling `P_bound_max` (lowered from
the ~off value 1.024; linear form attenuates attachment by `1 − P_bound/P_bound_max`).
**The cap helps but does not close the gap** (figure `results/k2_occupancy_landing.png`,
the force-ratio × ktr-ratio landing plane):

- It **rotates the model's force–ktr trade-off toward the data**: at a k2 that over-produces
  force, tightening the cap pulls the steady-force ratio down (at k2×0.40: steady ratio
  1.41→1.20 as P_bound_max 1.02→0.25), cutting the force-overshoot-at-correct-ktr from **+41%
  to ~+24%**, and lowering the cost at deep cuts (k2×0.30: 0.54→0.41).
- **But it erodes the ktr slowdown** — tightening the cap raises *baseline* ktr (66→97), so the
  ktr ratio drifts the wrong way (0.55→0.64). The two knobs are not independent: both act on
  attachment.
- **No (k2, P_bound_max) pair reaches the target** (force 1.18, ktr 0.54). Every iso-cap curve
  passes above-right of the data star. The model's force–ktr curve is **too shallow** — it
  cannot deliver the data's "large ktr drop for small force gain." Overall best stays the
  gentle **k2×0.50, cap off (core cost 0.088)**; adding the cap does not beat it.

**Why occupancy is the wrong second knob — and what's right.** The cap lowers force by lowering
*attachment*, which also speeds redevelopment (raises ktr), fighting what k2 buys. The data
needs a lever that **slows cycling without adding force** — i.e. moves *down-left* in the
force–ktr plane:
1. **Co-reduce attachment `ka` (global cycle slowdown).** Lowering `ka` *and* `k2` together
   keeps duty (force) ~fixed while dropping ktr = f_app+g_app. `ka` and `k2` are independent
   directions and span the plane, so this pair can in principle reach the star.
2. **SRX recruitment at low ATP.** Park heads OFF (heads can't refold the IHM without ATP),
   removing force *and* cycling capacity. Physiologically the strongest candidate
   (see `../LowATP_ForceEnhancement/`).

**Updated recommendation:** the minimal low-ATP model is **k2 + a cycling lever (`ka`), not
k2 + occupancy cap and not k2 + SRX** (SRX ruled out below). `vall_y` / `restretchSlopeStart`
remain k2- and cap-inert (need a stretch-detachment / 3rd-state lever).

## Third lever tested: SRX stabilization / destabilization

The model has a built-in ATP→SRX gate (`UseAtpOnUNR` + `K_T3`, gating `kmsr`): low ATP → less
SRX mobilization → larger SRX pool → lower force/ktr (*stabilization*). We tested both the
pool-size knob (`ksr0`) and the physical gate (`g4`):

- **SRX pool size (`ksr0` 190→1000) at k2×0.40: force/ktr flat** (steady 106.1→106.6, ktr
  35.8→35.1, attached 0.319→0.320) even though `SRX_ss` rose 0.028→0.087.
- **g4 gate (`UseAtpOnUNR=1`, `kmsr` 200–600), MgATP 8 vs 2: steadyR=1.00, ktrR=0.98** — no
  effect; combined with k2×0.40 it is indistinguishable from k2 alone.
- **Destabilization (`ksr0`↓) also does nothing** to the activated features.

**Why SRX has no leverage here:** every slack feature is at **max-Ca**, where the
force-dependent parking term `exp(−F_SR/σ2)` mechanosensitively suppresses SRX to ~0.02–0.09
(near-empty). That small pool also exchanges with the large idle PD/PT "ready" reservoir, not
the attached/force-bearing pool — growing it just drains idle heads without starving
attachment. **Verdict: neither SRX direction is an effective ATP lever for these activated
features.** SRX would only register at sub-maximal Ca (force–pCa), which this dataset lacks.
(The mavacamten trace 06 shows SRX *can* abolish force — but only via a far stronger,
activation-blocking stabilization, not the modest ATP-driven shift.)

## Fourth/fifth levers + the full map: Pi (p1→p2 block) and ka (`RunLeverMap.m`)

`results/lever_map.png` places every lever in the (steady-ratio, ktr-ratio) plane vs the data
star (1.18, 0.54; peak2 1.18 annotated):

- **Pi / power-stroke block.** Pi is a product of the power stroke, so high Pi drives the
  reverse (P2→P1). The model has `f1=(Pi/K_Pi)/(1+Pi/K_Pi)` for exactly this, but it is **only
  wired into the non-default ODE branches** — the active pchip path (lines 389–398) ignores it,
  so `params.Pi` is currently inert (and `f1`→0 at Pi=0 would zero the baseline reverse stroke).
  Tested via the faithful proxy (raise reverse `k_1` / lower forward `k1`): it **trims force but
  moves LEFT-UP** (k_1×3: steady 1.41→1.21, ktr 0.55→0.70) and **inflates peak2** (→1.28). Same
  wall as occupancy. Physiologically apt — Pi *lowers* force (Cooke & Pate), which helps explain
  why the data's force gain is only +18% not +41% — but it cannot slow ktr.
- **ka (attachment).** Reducing ka **trims steady/A and barely moves ktr** — `k2×0.40 + ka×0.7`
  lands (steady 1.20, A 1.20, ktr 0.57), the **closest any 2-lever gets** to (1.18, 0.54). BUT it
  **crashes peak2 to 0.96** (data +18%): fewer heads re-attach during the restretch, so the
  second peak collapses.

**The binding residual is peak2.** The data raises steady, A AND peak2 uniformly (+18%) while
dropping ktr 46%. k2+ka delivers steady/A/ktr; no lever keeps peak2 up while doing so (k2/Pi push
peak2 too high, ka pushes it too low). This is the *same* defect as the baseline peak2-overshoot —
the steady↔peak2 (isometric↔restretch-reattachment) coupling is structurally fixed in the 2-state
model. Closing it needs the 3rd-state / restretch-attachment mechanism (see SensitivityAnalysis &
the mech-tradeoff notes), not another rate knob.

| lever (from k2×0.40) | steady | ktr | peak2 | toward target? |
|---|---|---|---|---|
| k2 deeper | ↑↑ | ↓ | ↑ | overshoots force |
| occupancy cap | ↓ | ↑ | ↓ | left-up |
| Pi (block p1→p2) | ↓ | ↑↑ | ↑↑ | left-up |
| SRX | – | – | – | null (max-Ca) |
| **ka** | **↓** | **~flat** | **↓↓** | **closest (steady/A/ktr); peak2 crashes** |

**Bottom line:** the minimal achievable 2-lever is **k2 + ka** (force amplitude + ktr); **peak2 is
the irreducible residual** and is the concrete evidence for needing a 3rd cross-bridge state.

## Optimization target (how to draw it)

The 5 slack segments are the **same protocol in both files** (SLslack/SLdiff/t_seg identical;
v_restretch within ~2%), so they pair 1:1 and per-segment ratios are valid (`CompareATPData.m`;
figure `results/atp_per_segment.png`; object `results/atp_target.mat`). Scoring is relative:
`model(2mM)/model(8mM)` vs the data ratio. Target shape, set by per-segment constancy:

| group | features | profile (CV across segments) | target |
|---|---|---|---|
| force amplitude | steady, A, Am, peak2, peak1, vall | flat ≈1.18–1.25 (1–3%) | **scalar** per feature |
| kinetics | ktr (0.64→0.46), t0 (1.01→1.41) | **slopes with slack** (12–15%) | **per-segment vector** |
| noisy | restretchSlopeStart, ktr2_overshoot | erratic (30–52%) | de-weight |

Anchor weights: steady/A high (the clean +18% scalar); the **ktr 5-vector** high (its
slack-slope is a real constraint a correct mechanism must reproduce — not just the mean); t0
medium; noisy features ~0. The emitted `target` struct carries `scalar`, `vector`, `deweight`,
and `weights`.

## Reproduce
```matlab
cd(root); addpath(genpath('.'));
run Analyses/LowATP_k2Frontier/CompareATPData       % data compare + emit target object
run Analyses/LowATP_k2Frontier/RunK2Frontier        % ~2 min  — single-knob k2 frontier
run Analyses/LowATP_k2Frontier/RunK2OccupancyGrid   % ~2 min  — k2 × occupancy-cap 2D grid
run Analyses/LowATP_k2Frontier/RunLeverMap          % ~3 min  — all levers in force–ktr plane
```
