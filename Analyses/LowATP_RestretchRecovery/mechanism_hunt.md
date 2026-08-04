# Is the rsK level↔ratio weld an artifact, or physiology we are missing? (2026-08-04)

Three mechanism classes were built and tested. **All fail — and they fail in a way that
identifies the constraint.** Scripts: `../LowATP_AbsoluteFit/RunMechanismRatioScan.m`,
`../LowATP_AbsoluteFit/RunStrainDepGateTest.m`.

## What the data actually asks for

At 2 mM vs 8 mM the post-restretch window shows **both**:

| | 8 mM | 2 mM | |
|---|---|---|---|
| valley drop below steady | 16.8 | **18.4** | the valley is *at least as deep* at low ATP |
| recovery τ | 28.0 ms | **71.8 ms** | and the refill is **2.6× slower** |

So: heads really are lost during the stretch, **and** the pool refills slowly. Also
`rsK` (×0.36) is suppressed *more* than `ktr` (×0.55) — the post-stretch recovery is more
ATP-sensitive than the isometric one.

## Class 1 — detachment-side (slow the clearing of stretched heads): WRONG DIRECTION

- **R2 ceiling** (`k2max`): falsified at 8 mM by the ktr overstretch peak (`../ApoPoolDetachment`),
  and here it also worsens the ratio (0.56 → 0.63–0.70).
- **Strain-dependent ATP gate** — *new, built and tested here*:
  `Km_eff(s) = K_T2·(1 + kAtpStrain·max(0, s − sAtpStrain0))`, still `MgATP_ref`-normalised so
  **8 mM is bit-identical (verified `max|Δsteady| = 0`)** — it dodges the ceiling's falsifier
  entirely, which was the whole point of building it. Physiological motivation is sound
  (positive strain inhibits ADP release → the ATP step becomes rate-controlling for strained
  heads, Nyitrai & Geeves 2004).
  **Result: wrong direction.** `rsK2` 63.8 → 70.8, separation 1.11 → 1.24 (data 0.68).

| `kAtpStrain` | rsK8 | inv8 | rsK2 | rsK ratio | separation | peak1_2 |
|---|---|---|---|---|---|---|
| 0 (uniform) | 114.2 | 0 | 63.8 | 0.56 | 1.11 | 115.7 |
| 50 | 114.2 | 0 | 64.4 | 0.56 | 1.13 | 116.4 |
| 150 | 114.2 | 0 | 67.0 | 0.59 | 1.18 | 118.5 |
| 400 | 114.2 | 0 | 70.8 | 0.62 | 1.24 | 122.1 |
| *data* | *43.7* | | *15.9* | *0.36* | *0.68* | *111.9* |

**Why every detachment-side lever fails the same way:** slowing the detachment of strained heads
keeps them **attached**. Force stays high (`peak1_2` 115.7 → 122.1, the valley fills), so the
recovery has less distance to travel and the fitted first-order rate comes out **faster**. You
cannot simultaneously hold force during the stretch and recover slowly — but the data does exactly
that. **The mechanism cannot live on the detachment axis.**

## Class 2 — re-attachment side (slow the refill): also fails

| lever | rsK2 | steady2 | vall_y 2 mM | verdict |
|---|---|---|---|---|
| baseline | 63.8 | 90.0 | 74.9 | — |
| `RegAvailShorteningOnly = 0` | *unmeasurable* | 90.0 | **115.6** | fills the valley; wrong sign |
| `tau_reg` ×3 | 60.8 | 89.9 | 75.2 | negligible |
| `tau_reg` ×10 | 64.5 | 89.3 | 75.3 | non-monotonic, negligible |
| `A0` ×0.5 | 48.5 | **65.5** ✗ | 47.1 ✗ | lowers rsK but destroys force |
| *data* | *15.9* | *91.3* | *85.2* | |

The lengthening-availability idea recommended by `../ApoPoolDetachment` §3 **was tested via the
existing flag and does not work as formulated**: enabling lengthening modulation *raises*
availability, filling the valley (74.9 → 115.6 against data 85.2). `tau_reg` is inert; `A0` is
force-welded.

## Class 3 — shared cycle rate: the weld

`k2` and `eta_M` fix the 8 mM level but compress the ATP ratio (see
`../LowATP_AbsoluteFit/conclusions.md`). `eta_M` is the clean one (peak/force untouched).

---

# Verdict: substantially an ARTIFACT of the slack protocol

The repo's own [`../ATPEffectReconciliation/conclusions.md`](../ATPEffectReconciliation/conclusions.md)
already established the key fact:

> the slack-derived `ktr` is **~5× more compliance-sensitive** than the dedicated ktr protocol; the
> two disagree by 24 % (**ktr protocol ×0.72, CV 4 %** vs **slack ×0.58, CV 7 %**); they **agree at
> 8 mM and diverge only at 2 mM**. "For cross-bridge kinetics prefer ×0.72; the slack value is
> contaminated by whatever series compliance differs between the conditions."

**`rsK` is measured in that same slack protocol, in the window immediately after the restretch —
the most compliance-loaded part of the entire protocol.** If the slack `ktr` is already 24 %
contaminated at 2 mM, `rsK` is very likely contaminated *more*. Supporting points:

1. The rundown lesion itself was identified as largely a **compliance** lesion
   (`kstiff` ×0.84 + `kSE` ×0.65 — [`../RundownCorrection`](../RundownCorrection/conclusions.md)).
2. At 2 mM the fibre is **×1.36 stronger**, so the series elastic element is loaded ×1.36 more and
   post-stretch redistribution — between fibre and series element, and between non-uniform
   sarcomeres — is correspondingly slower. That is an ATP-*dependent* effect without being an
   ATP-dependent *rate constant*.
3. The model is a **single half-sarcomere**: it structurally cannot produce inter-sarcomere
   redistribution, the classic origin of slow post-stretch creep in skinned fibres.
4. Independent measurement here: post-restretch τ is **2.6×** ATP-dependent vs post-slack **1.8×**.
   The excess sits exactly where compliance loading is highest.

This reframes the failure: the model may be "failing" to reproduce something that is **not a
cross-bridge rate at all**, and forcing it to fit `rsK` would push cross-bridge parameters to
absorb a mechanical artifact — which is precisely what the level↔ratio weld and the `absFull`
force degradation (steady 87.5 vs data 91.3) already show happening.

## The decisive experiment — cheap, and not yet run

`data/protocol_03_27_2026_velocitytable_ktr_2mM.mat` exists (also 04/03 and 04/10). **The dedicated
ktr protocol contains its own overstretch → hold, i.e. its own post-stretch recovery**, but with a
far smaller compliance excursion than the slack restretch.

> **Measure the post-stretch recovery rate in the ktr protocol at 2 mM and 8 mM.**
> - ratio ≈ **0.7** (like the ktr-protocol `ktr`) ⇒ `rsK`'s ×0.36 is largely a slack-protocol
>   compliance artifact. Do not chase it: fit the recovery from the ktr protocol, or keep `rsK`
>   down-weighted.
> - ratio ≈ **0.36** ⇒ it is genuine cross-bridge physiology, and a new mechanism is required —
>   one that produces *deep valley + slow refill*, which nothing currently in the model can.

This single measurement decides whether there is anything to model. **Run it before building any
new state.**

## If it turns out to be real, the mechanism must satisfy

1. **Not** act on detachment (the valley must stay deep) — Class 1 is excluded.
2. Slow the **refill of the attachable pool specifically after lengthening** — and not via
   availability as currently formulated (Class 2 excluded as written).
3. Be ATP-gated *more strongly than the isometric cycle* (×0.36 vs ×0.55).
4. Be `MgATP_ref`-normalised (8 mM-inert) or 8 mM must be refitted.

A head ruptured under strain carries **ADP**, so it is not apo: it must release ADP, bind ATP,
hydrolyse, and re-register before rebinding. The ADP-release and hydrolysis steps (~40–350 s⁻¹) are
the right order for a ~25 ms refractory delay; ATP binding is **not** (>10 000 s⁻¹ at 8 mM,
`../ApoPoolDetachment` §1). So a post-rupture **refractory** pool is defensible, but its ATP
dependence must come from somewhere other than mass-action ATP binding — the credible candidate
being **local [ATP] depletion in the filament lattice** at high turnover, a real skinned-fibre
phenomenon and one that would naturally be worse in the most heavily loaded window.

## Practical recommendation

Keep `rsK` **down-weighted (0.1)** — which is what the current joint fit does. It keeps the
observable visible without letting a possibly-mechanical artifact drive cross-bridge parameters.

## Code status

`UseStrainDepAtpGate` / `kAtpStrain` / `sAtpStrain0` are wired in
`Model/dPUdT_CombinedTransitions.m` + `Model/getParams.m`, **default OFF and 8 mM-inert by
construction**. Kept (tested, falsified, documented) rather than deleted, per repo practice.
