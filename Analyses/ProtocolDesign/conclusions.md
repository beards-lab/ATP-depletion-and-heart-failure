# ProtocolDesign — how many sessions, in what order, with which repeats

> ### ⚠ §1–6 are SUPERSEDED by [§8](#8-revision-the-bracketed-alternating-ladder).
> They were costed under "brackets are expensive, so put them on a subset of sessions".
> Under the opposite constraint — **a bracket in every session** — the arithmetic changes
> sign and the recommendation with it. §1–6 are kept because the rundown-budget tables
> (§2) and the power tables (§3–5) are unchanged and still do the work; only the
> **recommended plan** (§6) is withdrawn.

**Headline (revised — see §8).**

> **Run `8-4-2-8` alternating with `2-8-2`, bracket in every session, 6 preparations,
> 21 activations, sex balanced at every rung.**
>
> Stacking 4 mM *ahead* of 2 mM deepens the descending-order bias to −12.5 %, which very
> nearly cancels `2-8-2`'s +16.0 %: **a naive 1:1 average of the two designs leaves
> 1.0–2.2 % across the whole φ range**, against 4.1–9.1 % for the plain `8-2`/`2-8`
> balance that strategy.md §2.4 rejected. That makes the corrected and uncorrected
> readouts agree to ~1 %, so the rundown correction becomes a *check* rather than a
> load-bearing assumption.
>
> **But the order inside the session must be `8-4-2-8`, not `8-2-4-8`.** With 4 mM in the
> third slot its correction (18.6 %) **exceeds its own signal (15.8 %)** and the
> uncorrected ratio reads **0.943** — 4 mM would appear *weaker* than 8 mM when the truth
> is ~16 % stronger. The sign of the measurement would be an artefact.
>
> **2/8 power:** SEM 4.6 % now → 3.3 % after 3 new sessions → 2.7 % after 6.
> **4 mM power:** 89 % that it differs from 8 mM, 77 % to catch a 15 % model error, at
> n = 3 (S5). Curvature still needs ~42 preps and is still not worth attempting.

`RunDesignPower.m` → `results/design_power.png` (§2–5).
`RunDesignLadder.m` → `results/design_ladder.png` (§8).
Extends [RundownCorrection/strategy.md](../RundownCorrection/strategy.md), which settled
order and pooling but did not cost sex, sample size, or a third dose.

---

## 1. What was already decided, and what this adds

`strategy.md` settled three things and they stand: **majority 8→2** (§2.4), **keep ≥1
reversed order to detect a broken correction rather than to cancel it** (§2.4), and
**P6 — solve the ATP effect `a` and φ jointly from both orders** (§3). It left three
things open, which is what is computed here:

| open question | status before | answered here |
|---|---|---|
| **Sex** | Appears only as a day label (`03/27 M`, `04/03 F`, `04/10 M`). No analysis in the repo treats it as a variable. | §5 — balance it, don't power for it |
| **N** | "3 preps per order is the design target" — asserted, no precision requirement behind it | §4 |
| **4 mM** | Not considered anywhere | §3 |

---

## 2. The rundown budget — what each extra activation costs

Permanent damage per run is `D = φ·|slope|·T_act`, in **kPa, not per cent** — that
absolute-subtraction property is what makes the orders asymmetric (`strategy.md` §2.2).
Measured inputs: slope 0.560 (8 mM) / 1.092 (2 mM) kPa/s, `T_act` 10.4 / 11.6 s,
reference `F₈` = 50 kPa. 4 mM is interpolated (§3.1).

Damage already suffered when each run *starts*, at φ = 0.581:

| sequence | #act | 8 mM | 4 mM | 2 mM | final 8 mM |
|---|---|---|---|---|---|
| **8-2** | 2 | 0 % | — | **5.1 %** | — |
| 2-8 | 2 | 14.7 % | — | 0 % | — |
| **8-2-8** | 3 | 0 % | — | 5.1 % | 21.5 % |
| **8-4** | 2 | 0 % | **5.8 %** | — | — |
| 4-8 | 2 | 10.2 % | 0 % | — | — |
| **8-4-8** | 3 | 0 % | 5.8 % | — | 16.9 % |
| 8-4-2 | 3 | 0 % | 5.8 % | **12.6 %** | — |
| 8-4-2-8 | 4 | 0 % | 5.8 % | 12.6 % | 31.7 % |
| 8-2-4 | 3 | 0 % | **18.6 %** | 5.1 % | — |

Uncorrected bias on the ratio against the session's first 8 mM run:

| sequence | 4 mM/8 mM | 2 mM/8 mM |
|---|---|---|
| 8-2 | — | **−5.1 %** |
| 2-8 | — | **+17.3 %** |
| **8-4** | **−5.8 %** | — |
| 4-8 | +11.3 % | — |
| 8-4-2 | −5.8 % | −12.6 % |
| 8-2-4 | −18.6 % | −5.1 % |

Three things fall out:

1. **Descending order wins at every dose.** 8→4 costs −5.8 % where 4→8 costs +11.3 %,
   the same ×2 asymmetry `strategy.md` §2.3 derived for 2 mM, for the same reason: the
   lower-ATP run decays faster *and* degrades the weaker run *and* lands in the
   denominator.
2. **Whichever dose runs last pays.** `8-2-4` is the worst design in the set — it puts
   an 18.6 % correction on the new, noisiest point.
3. **A third dose is not free.** `8-4-2` more than doubles the 2 mM correction (5.1 % →
   12.6 %) — it degrades the measurement you already have, to add one you don't.

### φ is the ×2 unknown, and it multiplies all of this

The repo quotes three values of φ from three routes: **0.45** (consistent with the 03/27
bracket, used in ATPEffectReconciliation), **0.581** (strategy §2.3), **0.889** (the P6
joint fit). The correction scales linearly with it:

| sequence | φ = 0.45 | φ = 0.581 | φ = 0.889 |
|---|---|---|---|
| 8-2 → correction on 2 mM | 3.9 % | 5.1 % | 7.7 % |
| 8-4-2 → correction on 2 mM | 9.8 % | 12.6 % | **19.3 %** |

At `8-2` the φ ambiguity is worth ±1.9 % on the ratio — tolerable. At `8-4-2` it is
±4.8 %, which is the same size as the entire 4 mM signal you are trying to measure.
**This is the quantitative reason to keep sessions at two activations.**

---

## 3. The 4 mM arm

### 3.1 4 mM is the geometric mean of 2 and 8 — which is the whole problem

`√(2×8) = 4` exactly. So a log-linear dose dependence predicts `R₄ = √R₂ = √1.34 =
1.158`, and the three candidate laws, all anchored on the measured `R₂ = 1.34`, land
within **9.7 %** of each other:

| dose law | R₄ | slope₄ (kPa/s) | separation from log-linear |
|---|---|---|---|
| hyperbolic (linear in 1/[ATP], MM-like saturation) | **1.113** | 0.737 | 3.9 % |
| **log-linear** (linear in ln[ATP]) | **1.158** | 0.782 | — |
| linear in [ATP] | **1.227** | 0.915 | 5.8 % |

That is a 3.9–5.8 % signal sitting on a 4–11 % between-prep scatter. It sets the whole
answer.

*(The 4 mM decay slope is an interpolation — the only untested input here. It is
self-resolving: the within-run slope is measured directly in every record, so the first
4 mM session replaces it with data. It affects the predicted correction, not the design.)*

### 3.2 How many sessions — it depends entirely on what the point is for

n required, at the three planning CVs (80 % power, α = 0.05):

| purpose | CV 4 % | CV 8 % | CV 11 % |
|---|---|---|---|
| Place `R₄` to ±10 % | 1 | 1 | 2 |
| **Place `R₄` to ±5 %** | 1 | **3** | 5 |
| Place `R₄` to ±3.5 % | 2 | 6 | 10 |
| Catch a 25 % model prediction error | 1 | 2 | 3 |
| **Catch a 15 % model prediction error** | 1 | **4** | 7 |
| Catch a 10 % model prediction error | 2 | 7 | 14 |
| Resolve curvature vs linear-in-ATP | 5 | 19 | 36 |
| **Resolve curvature vs hyperbolic** | 11 | **42** | 79 |

**Verdict: 3 sessions.** At the CV 8 % planning value that places the 4 mM ratio to
±4.6 % (SEM) and catches any model prediction error ≥15 %, which is the size a genuine
mechanistic failure would be. It is the point where the curve in `results/design_power.png`
(bottom left) crosses from cheap to expensive.

**What 3 sessions will *not* buy, and no affordable number will:** the shape of the
dose-response. Discriminating hyperbolic saturation from log-linearity needs **42 preps
per dose** at the realistic CV, and 11 even at the optimistic one. If the reason for
wanting 4 mM is "what is the functional form of the ATP dependence", the honest answer
is that this preparation's prep-to-prep scatter forecloses it — the money buys more by
going elsewhere. If the reason is "give the model a third point it was not fitted to",
3 is right and it is a genuinely strong test.

### 3.3 Run them as `8→4`, not as `8→4→2`

| design | sessions | 4 mM correction | 2 mM correction | notes |
|---|---|---|---|---|
| 3 × `8-4-2` | 3 | 5.8 % | **12.6 %** | contaminates the 2 mM data you already have |
| **3 × `8-4`** | **3** | **5.8 %** | untouched | 4 mM and 2 mM both referenced to 8 mM |

The 2 mM/4 mM contrast is not lost by separating them — it is formed across sessions as
`(2/8) ÷ (4/8)`, at a variance cost of √2 on an already-adequate SEM. 8 mM is the common
reference in every session, which is what makes this work. Pay √2 on a secondary
contrast rather than ×2.5 on a primary one.

---

## 4. N on the 2/8 axis

Between-prep SD of the corrected ratio is **0.053** (`strategy.md` §5.5), i.e. CV 4.0 %
at the mean of 1.34 — but that is estimated from **3 preparations (2 df)**, so its own
95 % interval is enormous. Everything above is therefore planned at **CV 8 %** and
bracketed by 4 % and 11 % (the raw pooled force CV from ATPEffectReconciliation).

Current holdings: **2 × 8→2** (03/27 M, 04/03 F) and **1 × 2→8** (04/10 M); one bracket
(03/27), and it is the sole calibration of φ.

Recommended additions: **+2 on 8→2, +2 on 2→8** → 4 and 3 per order. That satisfies P6's
stated requirement of ≥3 per order while keeping 8→2 the majority as §2.4 requires, and
delivers SEM 3.0 % on the pooled ratio at CV 8 % (1.5 % if the 4 % CV holds).

**Two of the four carry a bracket, one per order.** This is the single highest-value item
in the plan and it is not mine — ATPEffectReconciliation already flags it: *"A second
end-of-session repeat on another prep is the single most valuable next measurement."*
The `2-8-2` bracket additionally tests whether φ is ATP-dependent, which §2.3 gives
direct reason to doubt (the 2 mM run decays ×1.95 faster) and which one bracket cannot
resolve.

---

## 5. Sex

Nothing in the repo treats sex as anything but a filename. n per group to detect a sex
difference in the ATP ratio:

| effect size | CV 4 % | CV 8 % | CV 11 % |
|---|---|---|---|
| 30 % | 1 | 2 | 3 |
| 20 % | 1 | 4 | 6 |
| **10 %** | 3 | **12** | 21 |

**Recommendation: balance sex, do not power for it.** Balancing is free — it only
constrains which animal goes into an already-planned session — and it buys incidental
power against a gross (≥20 %) sex effect. Powering for a 10 % effect costs 24 sessions
at the planning CV, which is more than the entire programme.

The constraint that actually matters is **sex must not become confounded with ATP
order**, which is a live risk: the reversed-order prep (04/10) is male, so naively adding
male reversed-order sessions would alias the two.

---

## 6. The plan

Existing 3 sessions plus 7 new. Every session is `relax → [runs] → PNB+Mava`; the runs
column is the active sequence only.

| # | session | runs | sex | #act | what it buys |
|---|---|---|---|---|---|
| — | 03/27 | 8-2-8 | M | 3 | *existing* — the only measured φ |
| — | 04/03 | 8-2 | F | 2 | *existing* |
| — | 04/10 | 2-8 | M | 2 | *existing* — the only reversed order |
| A1 | new | **8-2-8** | M | 3 | bracket #2 — φ on a second prep, at the majority order |
| A2 | new | **2-8-2** | F | 3 | bracket at the reversed order — tests whether φ is ATP-dependent |
| B1 | new | 8-2 | F | 2 | replication at the low-bias order |
| B2 | new | 2-8 | F | 2 | reaches 3 per order for P6; breaks the sex↔order aliasing |
| C1 | new | **8-4-8** | F | 3 | 4 mM + its own bracket; measures the 4 mM decay slope directly |
| C2 | new | 8-4 | M | 2 | 4 mM replicate |
| C3 | new | 8-4 | M | 2 | 4 mM replicate |

**Totals — 7 new sessions, 17 activations.**

| axis | allocation | check |
|---|---|---|
| 8→2 | 03/27 M, 04/03 F, A1 M, B1 F | n = 4, **2 M / 2 F** |
| 2→8 | 04/10 M, A2 F, B2 F | n = 3, 1 M / 2 F |
| 8→4 | C1 F, C2 M, C3 M | n = 3, 2 M / 1 F |
| **sex overall** | 5 M / 5 F | **balanced, and crossed with order** |
| brackets | 03/27, A1 (8→2); A2 (2→8); C1 (8→4) | 4 total, one per axis + the existing |

**Order within every session is descending in [ATP]**, except the deliberate reversed-order
sessions, which exist to validate the correction and are kept to the minimum that P6
needs.

### Two carry-overs from `strategy.md` §4 that are still pending

* **Request a hi-res `*_stiff_ktr.txt` + `*_slack.txt` on at least one bracket.** This is
  the decisive test of whether the ktr protocol is really ~5× less compliance-sensitive
  than the slack ktr — currently a strong hypothesis with a sharp prediction (protocol
  `ktr` should fall ~3 % where slack `ktr` falls ~12 %) and **no measurement**. It decides
  whether the ATP effect on `ktr` is ×0.72 or ×0.58, which is a 24 % disagreement sitting
  in the fit target.
* **`mergeLogsAndBursts` currently excludes filenames containing `repeat` from burst
  matching.** One-line change, needed before any bracket's hi-res data will merge. With
  three new brackets planned this stops being cosmetic.

---

## 7. Analyse it staged, not at the end

The CV that drives every number above is estimated from 3 preparations. After the first
three new sessions (A1, A2, B1) re-estimate it and re-run `RunDesignPower.m`:

* If CV ≤ 8 %, the plan is right as written.
* If CV → 11 %, the 4 mM arm needs 5 not 3 to hold ±5 %, and the sex balance becomes
  purely a bias-control measure with no detection power at all.
* If CV → 4 % holds up with more data, drop B2 and take the 4 mM arm to 4.

The systematic order bias does **not** shrink with n (`strategy.md` §2.4) — only the
correction does. Extra sessions buy precision, never accuracy; the brackets buy accuracy.

![Design power](results/design_power.png)

---

## 8. Revision: the bracketed, alternating ladder

`RunDesignLadder.m` → `results/design_ladder.png`.

§6 was costed under an unstated assumption — that brackets are scarce, so they go on a
subset of sessions and φ is transferred to the rest. Under the opposite constraint,
**a bracket in every session**, three things change and the recommendation inverts.

### 8.1 A bracket everywhere removes the argument against stacking

§2's case against `8-4-2` was that it raises the 2 mM correction from 5.1 % to 12.6 %,
and that a bigger correction carries a bigger error *because φ is assumed*. If every
session brackets its own first condition, φ is **measured per preparation** and never
transferred — the correction's dominant error source is gone. A 12.6 % measured
correction beats a 5.1 % assumed one. §2's tables are still right; the conclusion drawn
from them is not, once φ stops being borrowed.

### 8.2 …but the slot order matters more than it did, not less

| | 8 mM | 4 mM | 2 mM | 8 mM |
|---|---|---|---|---|
| **`8-2-4-8`** | 0 % | **18.6 %** | 5.1 % | 31.5 % |
| **`8-4-2-8`** | 0 % | **5.8 %** | 12.5 % | 31.5 % |

Correction expressed against the effect each dose is trying to measure — the metric that
matters, because 4 mM's signal (`R₄−1` = 0.158) is less than half of 2 mM's (0.34):

| sequence | 4 mM correction / signal | 2 mM correction / signal | uncorrected `R₄` reads |
|---|---|---|---|
| `8-2-4-8` | **1.18** ⚠ | 0.15 | **0.943 — sign flipped** |
| `8-4-2-8` | 0.37 | 0.37 | 1.090 (truth 1.158) |

In `8-2-4-8` the correction on 4 mM is **larger than the entire 4 mM effect**. The
uncorrected ratio reads 0.943 — 4 mM would look *weaker* than 8 mM when it is ~16 %
stronger. Every conclusion at that dose would then rest entirely on the correction being
right, which is exactly the dependence the per-session brackets were meant to remove.

**Rule: the dose with the smallest signal takes the earliest slot.** 4 mM has roughly
half of 2 mM's effect, so it goes second, always.

### 8.3 Balanced orders *do* let you skip the correction — with `8-4-2-8`, not otherwise

Residual bias on the 2 mM/8 mM ratio from a **naive 1:1 average with no correction at
all**:

| design pair | φ = 0.45 | φ = 0.581 | φ = 0.889 |
|---|---|---|---|
| `8-2` vs `2-8` (strategy.md §2.4's case) | 4.1 % | 5.5 % | 9.1 % |
| `8-2-8` vs `2-8-2` | 4.1 % | 5.5 % | 9.1 % |
| **`8-4-2-8` vs `2-8-2`** | **1.0 %** | **1.3 %** | **2.2 %** |
| `8-2-4-8` vs `2-8-2` | 4.1 % | 5.5 % | 9.1 % |

The mechanism: strategy.md §2.4 rejected count-balancing because the two orders inflict
*unequal* damage (f ≈ 0.05 descending vs ≈ 0.15 reversed), so averaging −5 % against
+15 % leaves ~+5 %. Inserting 4 mM ahead of 2 mM deepens the descending bias to −12.5 %,
which lands close to `2-8-2`'s +16.0 %, and the average very nearly cancels.

**This is balance in *damage*, not balance in *count*** — and it is a property of this
specific pairing, not of balanced designs generally. Note the third row: `8-2-4-8` does
*not* have it, because the 2 mM run sits in slot 2 and its bias is unchanged at −5.1 %.

Consequence: the corrected and uncorrected estimates should agree to ~1 %. **Their
disagreement becomes a diagnostic** — if the two routes diverge by more than ~2 %, the
correction is broken, and you find out without an external standard.

### 8.4 The ladder

Alternating, bracketed, sex balanced at every rung so the programme can stop anywhere.

| # | session | sex | #act | n(2/8) | SEM | n(4/8) | SEM | pwr `R₄`≠1 | pwr 15 % err | M/F | desc/rev |
|---|---|---|---|---|---|---|---|---|---|---|---|
| — | *existing* | — | — | 3 | 4.6 % | 0 | — | — | — | 2/1 | 2/1 |
| **S1** | `8-4-2-8` | **F** | 4 | 4 | 4.0 % | 1 | 8.0 % | 45 % | 35 % | 2/2 | 3/1 |
| **S2** | `2-8-2` | **F** | 3 | 5 | 3.6 % | 1 | 8.0 % | 45 % | 35 % | 2/3 | 3/2 |
| **S3** | `8-4-2-8` | **M** | 4 | 6 | 3.3 % | 2 | 5.7 % | 73 % | 60 % | 3/3 | 4/2 |
| **S4** | `2-8-2` | **M** | 3 | 7 | 3.0 % | 2 | 5.7 % | 73 % | 60 % | 4/3 | 4/3 |
| **S5** | `8-4-2-8` | **M** | 4 | 8 | 2.8 % | 3 | 4.6 % | **89 %** | **77 %** | 5/3 | 5/3 |
| **S6** | `2-8-2` | **F** | 3 | 9 | 2.7 % | 3 | 4.6 % | 89 % | 77 % | 5/4 | 5/4 |

*(CV 8 % planning value; power at α = 0.05 two-sided. "pwr `R₄`≠1" = power to establish
that 4 mM differs from 8 mM at all; "pwr 15 % err" = power to catch a 15 % error in the
model's predicted `R₄`, propagating the anchor's uncertainty.)*

**21 activations across 6 preparations.** Sensible stopping points:

* **After S2** (2 preps, 7 activations) — the minimum that is self-validating: one
  bracket per order, so φ is measured in both directions and the naive-average check
  is live. 4 mM is placed but not yet significant.
* **After S4** (4 preps, 14 activations) — 2/8 at SEM 3.0 %, sex balanced 4/3 and
  crossed with order, 4 mM at 73 %. **The recommended stopping point if preparations
  are scarce.**
* **After S6** (6 preps, 21 activations) — 4 mM clears 80 % power; 2/8 at SEM 2.7 %.

### 8.5 What this still does not buy

Dose-response **curvature** — hyperbolic saturation vs log-linear — needs **n = 42** at
CV 8 % and is unchanged by any of the above. 4 mM is the geometric mean of 2 and 8, so
the candidate laws separate by only 3.9–5.8 % against a 4–11 % prep scatter (§3.1). The
4 mM arm is a **model-validation point**, not a dose-response study, and should be
reported as one.

**And no other dose rescues it.** Scanning the separation between the two laws across
the whole interval (`RunDesignLadder.m` §A extension):

| test dose | separation (log units) | n required |
|---|---|---|
| 2.5 mM | 0.0227 | 122 |
| 3.0 mM | 0.0340 | 55 |
| 3.5 mM | 0.0385 | 43 |
| **3.82 mM** *(optimum)* | **0.0391** | **41** |
| **4.0 mM** | **0.0390** | **42** |
| 5.0 mM | 0.0335 | 56 |
| 6.0 mM | 0.0237 | 113 |

4 mM is within 2 % of the theoretical optimum. The curvature question is out of reach
because of the *scatter*, not because of a poor dose choice — signal-to-noise is **0.49
per preparation** at the best possible dose. Choose 4 mM (round number, geometric mean,
`R₄ = √R₂` is a clean prediction to state) and do not attempt the curvature claim.

### 8.6 What each rung licenses — the publication ladder

The binding constraint on the *current* dataset is **not sample size** — it is that φ
rests on a **single bracket** (03/27), which is the weakest link in every corrected
number and the first thing a referee will find.

| claim | needs | reached at |
|---|---|---|
| "2 mM raises force ×1.36 and slows `ktr` ×0.55" | n ≥ 3 preps, both orders | **already held** (n = 3–4) |
| "…and our rundown correction is validated, not assumed" | ≥1 bracket **per order**, + the corrected/uncorrected agreement check (§8.3) | **S2** (2 preps, 7 act) |
| "…at SEM ≤ 3 % with sex crossed against order" | n = 7 on the 2/8 axis, 4 M / 3 F | **S4** (4 preps, 14 act) |
| "The model predicts the 4 mM response it was never fitted to" | n = 3 at 4 mM → 89 % power | **S5** (5 preps, 18 act) |
| "The ATP dose-response is hyperbolic / log-linear" | n = 42 | **never — do not claim** |

**Recommended cutoff: S4 without a 4 mM claim, S5 with one.** S2 is the point below
which the methods are attackable; S6 is comfort, not necessity.
