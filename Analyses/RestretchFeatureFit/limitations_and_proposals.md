# 2-state fit at cost 4.55: limitations & proposals

*2026-07-07. Base = `params/opt2state_opt.m` (the cost-4.48 optimum from the last
commit; re-evaluates to 4.55 here at MaxRunTime=120). Evidence:
`diagResidual_report.txt`, `probeCoupling_report.txt`.*

## 1. Where the 4.55 lives (per-feature decomposition)

| term | cost | model vs data |
|---|---|---|
| **ktr** | **1.70** | 68–72 /s vs **49** /s (40 % too fast) |
| **peak1_dSL** | **0.80** | 0.034–0.039 vs **0.022–0.028** µm (peak too late / too much excursion) |
| **peak1_y** | **0.53** | 100–104 vs **90–96** kPa (spike ~8–10 % tall) |
| ovrsht_dy | 0.38 | 0.9–1.9 vs 1.2–1.5 (scattered) |
| vall_y | 0.26 | 66–70 vs **70–77** kPa (valley ~7 kPa too deep) |
| t0_crossing | 0.20 | buckling-crossing timing |
| A | 0.18 | redevelopment amplitude — essentially matched |
| FV_fnorm | 0.12 | `[1 .91 .60 .30 .11]` vs `[1 .92 .66 .32 .11]` — **solved** |
| steady | 0.03 | `[80 80 81 81 64]` — **solved** |
| physiology group | **0.000** | XTOR 7.5, XTOR_vmax 11.7, SRX 0.083, attached 0.26, PT 0.13 — all in-bounds |

**81 % of the cost is two clusters: ktr (1.70) + restretch-peak shape (~1.97).**
FV, steady force, and every physiology bound are already satisfied — there is no
headroom being spent to hold them, so both residual clusters are "clean" targets.

## 2. ktr is NOT a structural iron wall — it is a frozen lever

Coupling probe, one lever at a time off the base (cost 4.548):

| lever | ktr | steady | peak1y | XTORvmx | FV(.5,1,2,4) | cost |
|---|---|---|---|---|---|---|
| BASE | 69.4 | 77.1 | 97.6 | 11.6 | .91 .60 .30 .11 | 4.55 |
| k1 ×0.75 | 70.8 | 74.1 | 96.6 | 11.9 | .94 .61 .31 .11 | 12.7 |
| k2 ×0.75 | 58.0 | 85.4 | 112.9 | 7.4 | **.82 .53** .26 .09 | 15.5 |
| cycle ×0.80 | 58.0 | 78.8 | 103.1 | 9.4 | **.84 .53** .26 .09 | 9.1 |
| **kmsrd 30** | **57.6** | 76.8 | 96.8 | 11.5 | .90 .56 .25 .09 | **3.67** |
| kmsrd 15 | 44.3 | 76.1 | 95.2 | 11.4 | .88 .46 .17 .08 | 5.06 |

- **Core-cycle rates (k1/k2/cycle) confirm the ktr↔FV weld**: `k2↓` lowers ktr but
  inflates force and collapses the FV shoulder/tail; a global slowdown lowers ktr
  and holds *force* but still breaks FV. Within the cross-bridge cycle, ktr and FV
  are genuinely welded (Brenner f+g). `k1↓` barely moves ktr at all.
- **The wall breaks off-cycle**: `kmsrd` (SR→SRD→PD return rate) lowers ktr 69→58
  while **holding steady force, peaks, FV shoulder, and every physiology bound**, and
  the *total cost drops 4.55→3.67*. `kmsrd 15` overshoots (ktr 44 < 49) and starts
  eroding the FV tail, so the optimum is ≈ kmsrd 18–25 (ktr ≈ 49–55).

**Mechanism.** `kmsrd` sets the slow force-redevelopment timescale via the ADP-super-
relaxed return path (mechanism-evaluation.md §1.4/§3.2): a *fast* return injects a fast
component into redevelopment that inflates the single-exponential-fitted ktr. Model
`ktr_rmse` (0.28–0.64) is far below data (1.5–1.7) — the model recovery is *too purely
single-exponential*; it lacks the slow secondary phase that pulls the data's fitted ktr
down. `kmsrd` is exactly that missing slow timescale.

**Why it was stuck.** `kmsrd = 60` — 3× over its own `parameterBounds` range [0.1, 20],
and **absent from every optimization pool to date** (the 2-state pool tunes `ksr0`, the
*into*-SRX rate, but never `kmsrd`, the *out-of*-SRX rate). It was pinned at a legacy
hand value the optimizer could never pull back.

**`dr1` (Route B) ruled out.** The code's own documented ktr↔FV decoupler (make pre-stroke
p1 force-bearing) does not help here: `dr1 0.010/0.015` inflate force massively
(steady 77→86→95, peak1y 98→114→131; cost 16→48) while barely moving ktr (69→66). At this
operating point p1 is too populated, so force-duty just adds bulk force.

## 3. The restretch-peak cluster is genuinely coupled (2-state)

| lever | ktr | steady | peak1y | note |
|---|---|---|---|---|
| kstiff2 ×0.85 | 70.1 | **68.4** | 90.2 | peak1y↓ but steady↓ too — welded via stiffness |
| kSE ×0.70 | 59.0 | 76.4 | 95.8 | peak1y↓ & ktr↓, FV held — but enlarges peak1_dSL (opposite pull) |
| P_bound ×0.85 | 70.9 | 74.7 | 95.1 | mild peak1y↓, coupled; net cost↑ |
| catchbond 0.3 | — | 75.6 | NaN | **unstable** — force runaway (peak2=440), not strain-limited |

No *static* lever separates `peak1_y` (height, ~kstiff), `peak1_dSL` (excursion, ~kSE the
other way), `vall_y` (valley depth, ~R1D detachment), and `steady` (force, ~kstiff): they
share the same elastic + detachment machinery. This is the same **timescale wall** already
documented for `vall_y`/`vall2` — a single strain-dependent rate cannot make a fast-shallow
valley *and* a slow-deep undershoot.

The escape is a **restretch-only, strain-rate-gated** force term (active only for vel>0, so
it never perturbs the isometric steady state):

- **Catch bond** — literature's top lever for "restretch rise rate" (mechanism-eval §1.12.4).
  Already coded but OFF; blew up in the probe only because `k_catch_bond=0.3` fully blocks
  detachment at high restretch velocity. Needs `k_catch_bond ≈ 0.005–0.05` **and** the
  strain-limited path (`CatchBondStrainMax`) so slip-region bridges still detach.
- **Kelvin–Voigt series element** (`UseViscoelasticSE`, `c_SE_visc`) — currently 0. Adds a
  velocity-dependent SE resistance during restretch that reshapes the spike/excursion
  without touching steady force.
- **3rd cross-bridge state** (`NumberOfStates=3`) — a slow-detaching, low-force `p3` gives an
  independent slow timescale that breaks the fast-shallow/slow-deep wall *and* supplies the
  slow ktr component. Durable but higher-effort (needs the R2 negative-strain knots re-tuned
  as a p2→p3 *feed*, not 2-state detachment — the known FV-collapse blocker).

## 4. Proposals (ranked by effort/return)

1. **[running] Free the SRX-return timescale in the 2-state pool.** `RunOptim2State_v2.m`
   adds `kmsrd, ksr2srd, ksrd2sr`; seed `params_2state_v2_seed.m` (kmsrd=18). Expected:
   ktr→~49, cost 4.55→~3.3–3.6. Low risk, mechanism- and probe-backed.
2. **Add a strain-limited catch bond (or Kelvin–Voigt SE) to attack the peak cluster.**
   First sweep `k_catch_bond ∈ {0.005,0.01,0.02,0.05}` with `CatchBondStrainMax` set, confirm
   stability + direction, then add to the pool. Targets peak1_dSL/peak1_y/vall_y without the
   steady-state tradeoff.
3. **3-state** for the durable break of the timescale wall (peaks + ktr together). Requires a
   non-collapsing seed (re-tune R2 negative knots as the p2→p3 feed).

## 5. Optimizer-procedure notes

- The residual is partly a **pool-coverage** problem, not model structure. Audit hand-set,
  out-of-pool parameters against `parameterBounds.m` (`kmsrd=60` vs ub 20 was the tell) and
  free the ones with a mechanistic tie to the live residuals.
- Keep `ktr|2` and the restretch-shape terms at full weight (no weight-gaming); the cost is
  now dominated by physically meaningful residuals, and the physiology group is at 0.

## 6. Outcome (2026-07-07) — proposal #1 executed

`RunOptim2State_v2.m` (free `kmsrd, ksr2srd, ksrd2sr`) then `ContinueOptim2State_v2.m`
(+ R2 shortening knots `PieceWiseStrainDep2Params__5/6/7`, RESUME) → new best
**`params/opt2state_v2_opt.m`, cost 3.22 (from 4.55, −29 %)**.

| feature | base (4.55) | v2 free-kmsrd (3.52) | v2 + R2 knots (3.22) | data |
|---|---|---|---|---|
| ktr | 68–72 | ~50 | **56 53 52 50 51** | ~49 |
| FV_fnorm | .91 .60 .30 .11 | .91 .53 .21 .09 | **.925 .55 .22 .09** | .92 .66 .32 .11 |
| steady | matched | matched | **80 80 81 81 64** | 80 80 81 81 64 |
| kmsrd | 60 | 19.1 | **20.6** | — (bound ≤20) |

- **ktr is solved and stays solved** (68→~50 = data) with steady force, FV shoulder, and every
  physiology bound held. Confirms the ktr wall was a frozen lever, not structure.
- **ktr and FV do NOT fully decouple.** Lowering `kmsrd` costs FV mid-tail (v=1–2). The R2
  shortening knots recover it only *partially* (FV residual 0.80→0.61; shoulder perfect, mid-tail
  .55/.22 vs .66/.32). This residual `ktr`↔`FV-mid-tail` tension is modest but real — the one
  place the fix has a genuine cost, and a candidate structural motivation for the 3rd state.
- **Top residuals now:** `peak1_dSL` 0.74, `FV_fnorm` 0.61, `ovrsht_dy` 0.40, `peak1_y` 0.39
  (ktr no longer in the top 5). The restretch-peak cluster (proposal #2/#3) is untouched by this run.

Infra fix: `optimizeFeatures.m` final verification `assert` was `< 1e-6` on a non-reproducible
objective (MaxRunTime watchdog + parallel pool ⇒ ~1e-3 noise), so successful runs exited 1
(this also killed `bgcqyp7ha`). Now `max(1e-2, 2 %)` — still catches gross divergence.
