# Justifying Direct Subtraction of the High-Ca Passive Trace

## Context

For each experimental protocol we have three co-recorded force traces under an identical prescribed length `u(t)`:

- `F_a(t)` — active trace (high Ca, cycling XBs)
- `F_p_hiCa(t)` — passive trace at high Ca (XB cycling blocked, titin in Ca-stiffened state)
- `F_p_loCa(t)` — passive trace at low Ca (XB off, titin in relaxed state)

**Independent identifications:**

- Serial stiffness `k_SE ≈ 1000 kPa/µm` (approximately linear over the operating range).
- Ca stiffens titin only — not SE, not XB stiffness per se.
- Low-Ca passive steady-state ≈ 3 kPa at operating SL; high-Ca passive steady-state ≈ 10 kPa.

We want `F_XB(t)`, the cross-bridge-only contribution, from the active recording.

## Mechanical topology

```
external load ── [ SE: F_SE(Δ) ] ── node A ── [ titin ∥ XB ] ── fixed end
                    Δ = u − xc                    at xc
```

- `u(t)` — prescribed external sarcomere length
- `xc(t)` — internal contractile-element length (unknown a priori)
- `Δ(t) = u(t) − xc(t)` — SE stretch
- Force balance at node A: `F_SE = F_titin + F_XB`
- Measured quantity: `F(t) = F_SE(t)` (external)

## Why naïve subtraction is structurally biased in general

Because the SE is compliant, xc in the *active* run differs from xc in the *passive* run even when `u(t)` is identical:

```
xc_p(t) = u(t) − F_p(t)/k_SE
xc_a(t) = u(t) − F_a(t)/k_SE
xc_a − xc_p = −(F_a − F_p)/k_SE ≈ −F_XB/k_SE
```

Titin sits at different lengths in the two runs, so

```
F_a − F_p = F_XB + [ F_T(xc_a, ẋc_a) − F_T(xc_p, ẋc_p) ]
                    \_______________________________/
                                 bias
```

The bias vanishes only if the SE is rigid (Δ → 0) or titin is sufficiently weak.

## Why direct subtraction is nearly unbiased here

With `k_SE = 1000 kPa/µm`, the SE stretch is bounded in magnitude by the force:

| State | Force (kPa) | Δ = F/k_SE (µm) |
|---|---|---|
| low-Ca passive | 3 | 0.003 |
| high-Ca passive | 10 | 0.010 |
| active peak | 50 | 0.050 |

The xc shift between the high-Ca passive and the active run is at most

```
|xc_a − xc_p_hiCa| ≤ (F_a,peak − F_p_hiCa)/k_SE ≈ 0.04 µm
```

— roughly 2 % of a 2.2 µm sarcomere. Titin, evaluated that close to its passive-run operating point, is nearly unchanged.

**Why the high-Ca passive and not the low-Ca passive:** titin in the active run is in the Ca-stiffened state (the prep is at high Ca to activate the XBs). The low-Ca passive reflects a *different physical state* of titin and would systematically overestimate the XB contribution by roughly the Ca-induced titin stiffening (~7 kPa at operating SL).

## Quantitative bound on the residual bias

Let `k_T = dF_T/dxc` be the local dynamic titin stiffness at the operating point. A first-order expansion gives

```
F_XB(t) = (F_a − F_p_hiCa) / (1 − k_T/k_SE)
       ≈ (F_a − F_p_hiCa) × (1 + k_T/k_SE)        for k_T ≪ k_SE
```

For cardiac titin in the operating range (SL ≈ 2.0–2.3 µm), `k_T` is of order 20–100 kPa/µm, giving a correction of 2–10 %. This is within other uncertainties in the cross-bridge model (attachment rates, duty ratio, overlap normalization) and smaller than typical experimental scatter.

`k_T(t)` can be estimated directly from the passive trace by local finite difference,

```
k_T(t) ≈ (dF_p_hiCa/dt) / (dxc_p_hiCa/dt)
```

and applied per-sample if needed.

## Recipe

```
# 0. Inputs: u(t), F_a(t), F_p_hiCa(t);  k_SE known.

# 1. Zeroth-order (recommended default):
F_XB(t) = F_a(t) − F_p_hiCa(t)

# 2. Internal kinematics, for any downstream XB model that needs xc:
xc_a(t)  = u(t) − F_a(t) / k_SE

# 3. Optional first-order correction:
k_T(t)   = finite_difference(F_p_hiCa(t), xc_p_hiCa(t))
F_XB(t)  = (F_a − F_p_hiCa) × (1 + k_T/k_SE)
```

## Assumptions and caveats

1. **Same `u(t)` across all three protocols.** Titin is viscoelastic, so subtraction assumes the titin strain history is identical between the passive and active runs. Small timing/rate mismatches in the prescribed length will produce transient artefacts after fast length steps.
2. **Linear elastic SE.** `k_SE = 1000 kPa/µm` is treated as constant. If the SE has its own viscoelastic time constant, Δ = F/k_SE holds only quasi-statically; expect transient errors immediately after fast length steps.
3. **Ca affects titin only.** If Ca turns out to modulate SE or some other parallel element, the high-Ca passive subtraction will absorb that effect into the apparent F_XB.
4. **First-order correction assumes `k_T < k_SE`.** True for cardiac muscle at physiological SL; may fail at SL > 2.5 µm where titin stiffens sharply.

## What the low-Ca passive trace is still good for

Even though it is not used in the F_XB extraction, the low-Ca trace remains useful:

- It gives the "naked" titin viscoelastic response, useful for identifying titin constitutive models independent of Ca effects.
- The difference `F_p_hiCa(t) − F_p_loCa(t)` isolates the Ca-induced titin stiffening — relevant for length-dependent activation mechanisms and for validating that the Ca-on-titin assumption is consistent across SL.
