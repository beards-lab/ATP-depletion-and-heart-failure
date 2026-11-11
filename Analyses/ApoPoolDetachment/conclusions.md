# ATP-limited detachment + tearing into an apo pool (2026-07-28)

Candidate mechanism for **defect 3** identified in
[`../RestretchVsKtrRecovery/conclusions.md`](../RestretchVsKtrRecovery/conclusions.md): the model's
post-restretch recovery runs at 73–100 s⁻¹ where the data says 44–59 s⁻¹, and neither the
Maxwell-dashpot fix nor the `tau_reg` fix touches it.

Baseline for everything below: TNF + `UseMaxwellTensionOnly` + `tau_reg = 2 ms`.
Reproduce: `RunOptimApoPool.m` (refit driver), `apoPoolCost.m` (objective).

## What was implemented

Two independent, **default-OFF** mechanisms in `Model/dPUdT_CombinedTransitions.m`:

**`UseR2Ceiling`** — series cap on strong-bridge detachment:

    1/R2_eff = 1/R2_strain(s) + 1/R2max,    R2max = k2max * MgATP/(K_T2c + MgATP)

**`UseD0State`** — a detached, nucleotide-free pool `D0`, fed by `XB_Ripped` (mechanical tearing),
drained by `RD0T = k_D0T * MgATP/(K_D0T + MgATP) * D0`. One scalar ODE state appended after `A_reg`.
`D0FromR2` optionally routes a share of the normal p2→PT flux through it too.

## Does it work? Yes — it is the only thing tested that reaches the data

| variant | post-slack #2 | post-restretch | peak1 | valley | E_FV | E_slack | max D0 |
|---|---|---|---|---|---|---|---|
| **DATA** | 50.0 | **49.2** | **1.154** | **0.796** | | | |
| baseline | 45.4 | 92.7 | 1.266 | 0.862 | 0.020 | **180** | — |
| ceiling 2000 | 44.5 | 106.9 | 1.341 | 0.920 | 0.045 | 366 | — |
| ceiling 800 | 43.3 | 100.0 | 1.460 | 0.965 | 0.156 | 964 | — |
| ceiling 300 | 40.3 | 184.5 | **1.620** | 1.007 | 0.315 | 3323 | — |
| ceiling 300 + rip 3000 | 41.1 | 68.9 | 1.480 | 0.919 | 0.313 | 2473 | 0.107 |
| ceiling 300 + rip 10000 | 41.7 | 60.1 | 1.379 | 0.899 | 0.310 | 2079 | 0.744 |
| ceiling 150 + rip 10000 | 38.4 | **44.8** | 1.425 | 0.887 | 0.461 | 5057 | 0.592 |

1. **The ceiling alone produces the predicted failure**: peak1 climbs 1.266 → 1.620 (data 1.154),
   the valley fills to 1.007, and post-restretch gets *worse* (92.7 → 184.5). Nothing detaches.
2. **Tearing is the necessary release valve** — it pulls the peak back down *and* delivers the slow
   recovery: 184.5 → 68.9 → 60.1 → **44.8** against a measured 49.2.
3. **Tearing without the ceiling does nothing** (post-restretch 92.7 → 86–92 across k2rip 200–8000):
   R2 has already annihilated those bridges. The two only work as a pair.

Cost at these settings: `E_slack` 180 → 2000–5000, `E_FV` 0.020 → 0.31–0.46. Nothing else was
re-optimised — hence `RunOptimApoPool.m`.

---

# The large-strain detachment is SHARED, and the ktr protocol already measures it

The unbounded strain acceleration of `R2(s)` is not an oversight — it is the mechanism that kills the
restretch peak overshoot, and it was tuned for that. But it is one function serving several
protocols, and it was identified on only one of them.

**Strain reached, and R2 there:**

| window | max ⟨s⟩ of p2 | max R2_eff | p2 before → after |
|---|---|---|---|
| slack restretch #1 (→1.10 ML) | 0.022 µm | 693 s⁻¹ | 0.191 → 0.048 |
| slack restretch #2 | 0.023 µm | 805 s⁻¹ | 0.184 → 0.044 |
| **ktr overstretch (→1.05 ML)** | **0.050 µm** | **2989 s⁻¹** | 0.012 → 0.019 |

The ktr overstretch drives the model **2.2× further in strain and 4× higher in rate** than the slack
restretch that identified `R2(s)`. The ktr protocol is extrapolating the far tail of a piecewise
function that no other scored protocol reaches. `RunForceVelocityTime = 0` in TNF, so the FV
timecourse — the other protocol that would exercise it — is not in the cost at all.

**And the model's ktr rate is entirely a product of that tail.** Removing the 1.05 ML overstretch
from the protocol (0.80 → 1.00 direct, same duration):

| | ktr k63 | p2 at hold start |
|---|---|---|
| with overstretch (real protocol) | 53.0 | 0.052 |
| without overstretch | **119.1** | 0.035 |
| *data* | *49.4* | |

So the post-fix agreement on ktr (53.0 vs 49.4) is right for a reason that has nothing to do with
cross-bridge cycling: the overstretch happens to wipe p2 to 0.0008 and redevelopment starts from a
clean slate. Change the far tail of `R2(s)` and the ktr rate moves for purely mechanical reasons.

## …which means the ktr overstretch peak is a direct measurement of the ceiling — and it falsifies it

Force during the overstretch is in the ktr trace, and already in `E_ktr`. Nobody was reading it as a
constraint on `R2(s)`, but that is exactly what it is:

| | overstretch peak (F_iso) | hold at 1.05 ML |
|---|---|---|
| **DATA 03/27 / 04/03 / 04/10** | **0.72 / 0.80 / 0.83** | 0.37 / 0.35 / 0.43 |
| model, no ceiling | **0.746** ✓ | 0.618 |
| model, ceiling 300 | 2.249 | 1.670 |
| model, ceiling 300 + rip 3000 | 1.303 | 0.814 |
| model, ceiling 300 + rip 10000 | 1.183 | 0.655 |
| model, ceiling 150 + rip 10000 | 1.586 | 0.723 |

**The uncapped model is exact (0.746 vs 0.72–0.83). Every capped variant is 1.4–3× too high**, and
tearing only partly rescues it (2.25 → 1.18 at rip 10000, still ~50% over). The same over-prediction
shows in the slack restretch peak (1.379 vs 1.154).

This is a direct experimental falsification of the ceiling, on data already in the cost function —
stronger than the literature-bounds argument below, and it does not depend on isoform or temperature.
**The muscle really does shed force that fast at 25 nm strain.** Whatever slows the post-restretch
recovery, it is not a cap on detachment.

---

# CRITICAL REVIEW OF THE PHYSIOLOGY

Checked against `Model/parameterBounds.m`, which is anchored on **mouse/rat α-MHC** (the *fast*
cardiac isoform) with citations. Three of the four justifications do not survive.

## ❌ 1. "ATP-binding ceiling" is the wrong name — by ~30×

`parameterBounds.m:266-268` (k3 entry), citing Lowey 2013 (PMID:23564459) and Deacon 2012
(PMID:22349210):

> At saturating [MgATP] (5–8 mM) the ATP-binding sub-step is effectively instantaneous for α
> (K1k+2 ≈ 1.6–2.0 µM⁻¹s⁻¹ → **>10,000 s⁻¹ apparent**).

The fitted `k2max` is 150–300 s⁻¹. **ATP binding is 30–70× faster than the cap being attributed to
it**, and at 8 mM it cannot be rate-limiting for anything. The `MgATP/(K_T2c+MgATP)` factor with
`K_T2c = 0.15 mM` is likewise inert (0.982 at 8 mM vs 0.930 at 2 mM — a 5% swing where the low-ATP
data wants 2.5×). To get the observed [ATP] effect from this term you would need `K_T2c ≈ 5–10 mM`,
against a measured Km for ATP-induced dissociation of ~0.1–0.2 mM. **Not defensible.**

## ⚠️ 2. A ceiling IS defensible — but as ADP release, and only in β territory

`parameterBounds.m:238-255` (k2 entry) is explicit that the post-stroke step is **ADP release**:
- Lowey 2013: mouse α-S1 ADP release **~350 s⁻¹** at 20 °C (β-S1 = 130)
- Deacon 2012: α-S1 ADP release **>1000 s⁻¹** — "so fast it does **not** limit detachment (unlike β)"
- Wang 2013: fibre-level post-stroke detachment **65–160 s⁻¹** at 17 °C
- bounds `[50, 1500]`

So a cap **is** warranted in principle — the model's `R2(s)` reaching **4600 s⁻¹** is 3–70× above
every α-MHC anchor, and that much is a real defect. But the settings that actually *work* are
`k2max` = 150–300 s⁻¹, i.e. at or below the fibre-level figure and **below the α-MHC isolated-S1
range entirely** — squarely in β-myosin territory, the wrong isoform for this preparation. The
ceilings that stay inside the α range (800, 2000) make post-restretch *worse*, not better
(100.0 and 106.9 s⁻¹).

**This is the sharpest objection to the mechanism: it only delivers when pushed outside the
species-appropriate kinetics.** Either the cap is not the explanation, or the model's `k2` is
absorbing something that is not ADP release.

## ❌ 3. `k_D0T = 40 s⁻¹` cannot be ATP binding — by ~250×

Same anchor: ATP binding to a rigor head is >10,000 s⁻¹ at 8 mM. A 25 ms refractory pool is not an
ATP-binding phenomenon at physiological [ATP]. If D0 is to be slow, the slow step must be named
honestly. The candidates, with the repo's own numbers:
- **re-hydrolysis** (`kah`, bounds [40,200], mouse α 150 s⁻¹ Deacon 2012) → τ ≈ 5–25 ms. **Right
  order.** But the model *already has* this step (PT → PD).
- **re-registration** after losing the target zone → this is what `A_reg` already models.

## ⚠️ 4. D0 works in the model for a bookkeeping reason, not a chemical one

`PT` (0.11) and `PD` (0.55) are large *pre-filled* buffers, so routing heads through them creates no
bottleneck. `D0` starts **empty**, so it acts as a sink with no reservoir. That — not nucleotide
chemistry — is why it slows recovery. The same effect would follow from any initially-empty
intermediate. This should be stated plainly rather than dressed as ATP dependence.

## ⚠️ 5. The fitted tearing is outside its own bounds and predicts the wrong energetics

- `bounds.k2rip.ub = 500`; the working fits use **3000–10000**, i.e. 6–20× over.
- The strain knee at `dr3 = -0.025` puts the rupture threshold at **25 nm** of half-sarcomere
  strain. The working stroke is ~10 nm and elastic rupture is expected by ~15–20 nm. `bounds.dr3`
  allows [-0.05, 0], so this is legal, but it is at the far edge.
- `max D0 = 0.59–0.74` means **59–74% of all heads are simultaneously in the torn state**, five
  cycles in a row.
- **Energetic contradiction:** those heads must re-prime (hydrolyse) before rebinding, so this
  predicts a large stretch-induced ATPase burst. Lengthening contraction is classically *cheap*
  (Abbott & Aubert) — the Fenn effect runs the other way. A mechanism that strips three quarters of
  the myosin population per stretch and pays an ATP for each is hard to reconcile with that.

  The escape hatch is real but changes the model: a bridge ruptured from `p2` (A·M·ADP) still
  carries **ADP**, so it is *not* apo. It should return to a nucleotide-bound detached state and
  re-enter near `PD`, not require ATP binding at all. Under that reading `D0` is mislabelled.

## What survives

**Sound:** the structural diagnosis. `R2(s)` conflates two physically distinct pathways —
(a) chemical cycle completion (ADP release → ATP binding → detachment, bounded, ends in M·ATP) and
(b) mechanical rupture (unbounded in strain, ATP-**in**dependent, ends nucleotide-**bound**). The
current code gives the rupture branch the chemical branch's destination. Splitting them is right,
and the unbounded 4600 s⁻¹ is a genuine defect.

**Not sound:** every place `[ATP]` was invoked. Neither the ceiling nor the D0 exit is
ATP-binding-limited at 5–8 mM, and the Km values needed to produce the low-ATP signature are 30–70×
off the measured ones.

## ⛔ VERDICT (2026-07-28): do not run the ceiling refit

The ktr overstretch peak settles it. `R2(s)` is pinned at **both** ends of the strain range by data
already in the cost — the slack restretch peak at 0.022 µm and the ktr overstretch peak at 0.050 µm —
and the uncapped function matches both. A ceiling that fixes the post-restretch recovery breaks the
overstretch peak by 1.4–3×. Add the 30–250× mismatches in the ATP justifications below, and the
mechanism is not worth a 9-parameter, ~6 h-per-run re-identification against an 11-DOF
identifiability ceiling.

`RunOptimApoPool.m` is kept as a *falsification harness*, not a fitting plan: it is seeded in bounds
and it scores the physiology bounds, so if anyone wants to re-open this, the run will show plainly
what it costs. It should be extended with the overstretch peak as an explicit target before use.

**What the exercise bought, which is not nothing:**
1. The ktr overstretch peak is a **free, unused constraint on the large-strain detachment**. It was
   sitting inside `E_ktr` as trace MSE, invisible. Promote it to an explicit feature.
2. The post-restretch defect now has a much narrower search space: the mechanism **cannot live on the
   `R2(s)` strain axis**, because that axis is now measured at both ends.
3. The obvious remaining candidate costs nothing to test: `RegAvailShorteningOnly = 1` means
   lengthening does *not* modulate registration availability — and the post-restretch window is
   precisely a post-lengthening window. A lengthening-driven availability penalty would slow
   reattachment without touching detachment, i.e. exactly in the direction the data points and
   without disturbing either measured peak. Try that before any new state.

**If the mechanism is ever revisited, reframe it first:**
1. Rename the ceiling to an **ADP-release / cycle-completion** cap; drop `K_T2c` or set it to a
   value that makes the term inert, and let `k2max` be bounded by `bounds.k2` **[50, 1500]** with
   its weight active so the optimiser cannot walk into β territory unpunished.
2. Either rename `D0` to a **post-rupture refractory** state and justify its rate by re-hydrolysis
   + re-registration (~40–150 s⁻¹, which `kah` already brackets), or route rupture to `PD` and get
   the delay from `A_reg` instead — in which case D0 may be unnecessary.
3. Restore `k2rip ≤ 500` (its own bound) and see whether the mechanism still delivers. If it needs
   6–20× over bound and 74% of heads, **that is the falsification**, not a fitting detail.
4. Score the stretch ATPase explicitly, so the energetic contradiction in §5 shows up as a cost
   rather than as an unnoticed consequence.

The refit driver applies (1) and (3) by keeping `evalPhysiologyCost` in the objective with the
α-MHC bounds active. That is deliberate: **if the mechanism can only be fitted by violating the
bounds, the run should fail loudly.**

## Implementation notes

Wiring touched three files and is inert when off — see
[`../RestretchVsKtrRecovery/conclusions.md`](../RestretchVsKtrRecovery/conclusions.md) §"Breaking-change
assessment" and `Workbench/Tests/TestD0State.m`. Regression with both flags off is exact:
`E = [0.016 0.984 0.010 183.37]`. Known gap: `Auxiliary/diagnostics/ktrProbe.m` hand-builds
`PU0 = zeros(1, 2*ss+7)` and supports neither `A_reg` nor `D0`.
