# 0002 — What "better fit on the new data" means

- **Status:** proposed. The metric is the PI's (2026-11-09, on #16); the comparison rule and the two open points below wait for his confirmation.
- **Date:** 2026-11-11
- **Issue:** #16 (milestone A)
- **Decided by:** Filip Jezek (metric); drafted by Claude (AI agent)

## Context

Milestone A refits the model on the new 8 → 4 → 8 → 2 → 8 mM data set. The done-criterion has to be
fixed before those data arrive, or any number can be argued afterwards. The candidates in #16 were
the pooled-target χ² per feature, held-out scores (ktr-protocol k80, the staircase) and the 2/8
ratios.

## Decision

**A fit is better when its weighted cost function is lower. The cost function is the one used in
the latest low-ATP fits** (`Analyses/LowATP_Refit9K2Sweep`, 2026-09-26).

### The metric

For one parameter set, on one data set:

- **J₈**, the 8 mM cost: the feature objective against the pooled 8 mM target, each term
  `((model − data)/data)²` weighted and floored at s²/n (`evalFeatureCost` with the `fn` list of
  `RunRefit9Optim`; `evaluateBakersExp` with `OptimizeOn = 'Feats'`).
- **J_low**, the low-ATP cost, one per low ATP level (2 mM today; 4 mM when the new set has it):
  the same terms against that level's pooled target, with the floors moved onto that level's spread
  (`rescaleFnToCondition`), and each term multiplied by its **ATP relevance**
  `a = min(z, 3)/3`. z is the z-score of the measured ATP effect on that feature,
  `|r − 1| / √(s² + s_r²)`. Bounds and guards keep weight 1. This is `score_rel` of
  `RunRefit9K2Sweep.m`. It includes the in-preparation 2 ML/s F–V point (FV2).
- **J = J₈ + Σ J_low.**

The procedure is fixed here; the numbers inside it are not. Targets, spreads, floors and relevance
weights are recomputed from the new data by the same code, **once, before any refit**, and then
frozen for every parameter set that is compared.

### Baseline of the current best on the old set

`robust1_seed` (8 mM) with fit A at 2 mM (R2 × 0.54, `pAtp_k2 = log(0.54)/log(1/4)`), pooled
03/27 + 04/03 + 04/10 targets, n = 3:

| Term | Value | Source |
|---|---|---|
| J₈ | 0.4240 | `LowATP_Refit9K2Sweep/conclusions.md`; reproduced 2026-11-10 (#88) |
| J_low at 2 mM, relevance-weighted (`score_rel`) | 4.52 | `LowATP_Refit9K2Sweep/conclusions.md`, not re-run for this record |
| J_low at 2 mM, unweighted (`score_plain`) | 7.86 | same |
| **J** | **4.94** | |

### The comparison rule

1. Build the targets and weights from the new data and freeze them.
2. Score the incumbent on them **without refitting**: J_incumbent. This is the number to beat.
3. A candidate is **better** when J_candidate < J_incumbent by more than the evaluation noise of
   the objective (0.3 % of J, `Analyses/StrainWindowNoise`), on the same frozen targets.
4. Report J₈ and each J_low separately next to J, and the unweighted `score_plain`, so a gain
   bought at one ATP level by a loss at the other is visible.
5. Report, without gating on them: the ktr-protocol k80 at each ATP level (held out on purpose),
   the staircase, and the Beard2022 F–V and P–V (#75). A candidate that lowers J and worsens a
   held-out score by more than its standard error is flagged, not rejected.

## Open points for the PI

1. **How the ATP levels combine.** Proposed: a plain sum, J₈ + J₂ + J₄. On the old set the 2 mM
   term is ten times the 8 mM term (4.52 against 0.42), so the sum is driven by low ATP.
   Alternative: compare J₈ and J_low separately and require neither to get worse.
2. **A minimum improvement.** Proposed: any improvement above evaluation noise counts. Alternative:
   require the gain to exceed the spread of J across re-draws of the target (leave-one-preparation-out),
   which protects against fitting one preparation.

## Alternatives considered

- **Unweighted χ² per feature (`score_plain`).** Simpler, but features that do not respond to ATP
  dominate the low-ATP term and reward a fit that gets the ATP-insensitive part right.
- **Held-out scores as the criterion** (ktr-protocol k80, staircase). They stay held out; using
  them to decide would turn them into training data.
- **2/8 ratios only.** Scale-free, but blind to the 8 mM fit, which carries the mechanics.

## Consequences

- The scoring code path (`RunRefit9K2Sweep.m`: targets, `rescaleFnToCondition`, relevance weights,
  FV2) becomes the reference implementation and should be lifted into one function before the new
  data arrive, so the incumbent and every candidate are scored by the same call.
- A 4 mM pooled target and its relevance weights have to be added to that function.
- Changing a feature, a floor or a weight after the new targets are frozen needs a new decision
  record.
- The relevance weights depend on the measured ATP effect, which on the old set is confounded with
  running order (#88, #92). The 8 mM brackets of the new design remove that.
