# Analyses

One self-contained folder per analysis: the **script(s)** that run it, a **`results/`**
folder for its `.mat`/figure outputs, and a **`conclusions.md`** that records what it found
and what to do next. This page is the index and the cross-analysis synthesis.

> Convention: an *analysis* lives here (script + results + conclusion). The **playground**
> (new drivers, hand-tuning, optimization entry points, tests) lives in `Workbench/`.
> Reusable, non-critical helpers live in `Auxiliary/` (diagnostic probes in
> `Auxiliary/diagnostics/`). Core model code lives in `Model/`.

> ### 📋 [Rundown & pooling strategy](RundownCorrection/strategy.md) — read before designing or analysing the next dataset
> Decision document: what to correct, what not to, which running order to request,
> where to put the repeats, and how to pool. **Key calls:** correct the *second record
> of each prep* (force **and** slack `ktr`, same lesion) and fit one prep at a time;
> and **2→8 carries ≈3.2× the bias of 8→2** — with the derivation in §2 — so
> balancing the design does *not* cancel it.

## Start here — the numbers that are settled

| question | answer | where |
|---|---|---|
| What does 2 mM ATP do? | **force ×1.34–1.35**; `ktr` **×0.72 from the ktr protocol** (CV 4 %) — *not* the ×0.58 the slack gives, which is compliance-contaminated | [ATPEffectReconciliation](ATPEffectReconciliation/conclusions.md) |
| Which `ktr` measurement should I trust? | The **ktr protocol** for kinetics; the **slack** `ktr` is ~5× more compliance-sensitive and is the right probe only for rundown | [RundownCorrection](RundownCorrection/conclusions.md) §5f |
| Does the reversed-order prep reconcile? | **Yes for force** (1.65 → 1.41 vs 1.31/1.32). Not for `ktr` — its scatter is not dose-driven | [ATPEffectReconciliation](ATPEffectReconciliation/conclusions.md) |
| How do I correct for rundown? | Reference each run to its **own** activation peak, then `loss = 0.45·\|slope\|·T_dec` on **effective activated time**; correct force, **not** `ktr`. Recipe in [methods.md](RundownCorrection/methods.md) | [RundownCorrection](RundownCorrection/conclusions.md) |
| Does the fibre decay between activations? | **No — it recovers slightly (+3 %/rest).** Decay acts only while activated | [RundownCorrection](RundownCorrection/conclusions.md) §5b |
| What *is* rundown? | **~16 % fewer heads in parallel + ~35 % added series compliance** (cell tearing). Tested against the 5-segment profile, not just the mean | [RundownCorrection](RundownCorrection/conclusions.md) §4b |
| Can I correct `ktr` too? | **No.** Free scan of the coupling: every defensible value makes cross-prep CV worse (9 %→16 %). Force and `ktr` reach the lesion by different channels | [RundownCorrection](RundownCorrection/conclusions.md) §5e |
| Which slack features can I trust? | `ktr` ratio (CV 8 %), normalised `A`/`peak1_y`/`peak2`/`vall_y` (CV 4–7 %). Not `ovrsht_*`, not `ktr2_*` | [SlackDataAnalysis](SlackDataAnalysis/conclusions.md) |
| Can I pool the Baker data? | `ktr` yes; **amplitudes no** — different protocol, truncated recovery | [SlackDataAnalysis](SlackDataAnalysis/conclusions.md) |
| How many parameters are identifiable? | **11**, from the Slack/FV/Ktr feature set | [SensitivityAnalysis](SensitivityAnalysis/conclusions.md) |
| Is `ktr` a separate process from restretch recovery? | **No** — one rate, ~41–47 s⁻¹, across all three rises | [RestretchVsKtrRecovery](RestretchVsKtrRecovery/conclusions.md) |

## Index

| Analysis | Status | One line | Conclusion |
|---|---|---|---|
| **AttachedPool_vs_ATPase** | findings | Reconciles attached-pool / ATP-turnover / k_tr / V_max at baseline; the turnover *target* was misset, leaving a fixable throttled-hydrolysis gap. | [conclusions](AttachedPool_vs_ATPase/conclusions.md) |
| **BindingSiteOccupancy** | findings | Audits attachment-saturation: the strain axis is *not* a binding-site axis, so existing per-bin "occupancy" is mis-grounded; evaluates a faithful global form. | [conclusions](BindingSiteOccupancy/conclusions.md) |
| **RestretchMechanisms** | findings | Explains the last-slack restretch double-peak; Phase-1 fit overshoots peak1 (+20%); points to kSE / attachment-shift tuning. | [conclusions](RestretchMechanisms/conclusions.md) |
| **SensitivityAnalysis** | findings | 108-parameter sensitivity + SVD: an **11-DOF identifiability ceiling**; ekSE/estiff/dr2 dominate. | [conclusions](SensitivityAnalysis/conclusions.md) |
| **LastSlackIdentification** | stub | Piecewise strain-dependent detachment fit to the last-slack transient. | [conclusions](LastSlackIdentification/conclusions.md) |
| **PassiveForceID** | findings | Justifies direct subtraction of the high-Ca passive trace to isolate F_XB, from the mechanical topology. | [conclusions](PassiveForceID/conclusions.md) |
| **LowATP_ForceEnhancement** | findings | 2 mM ATP gave ~32–44% more isometric force than 8 mM in the 03/27/2026 data; cross-referenced with literature. | [conclusions](LowATP_ForceEnhancement/conclusions.md) |
| **SlackDataAnalysis** | findings | Cross-dataset audit of all 9 slack recordings, one extraction path. **Baker is a different protocol** (2–6× shorter recovery) so its amplitudes are truncated and its "low ATP loses force" is an artifact; its `ktr` is fine. **Absolute force is untransferable (CV 20 %) but shape collapses to CV 4–7 %** once each prep's strength is divided out ⇒ fit normalised shape + a per-prep scale. Ranks every feature by reliability; retires `ovrsht_*` and `ktr2_*`. | [conclusions](SlackDataAnalysis/conclusions.md) |
| **RundownCorrection** | findings | Rundown is **linear in effective ACTIVATED time** (the within-run slope is force-independent: −0.459 fresh vs −0.462 at −17.5 % force), and only **φ ≈ 0.45** of each run's within-run decline is permanent. **Decay acts only while activated — between-activation change is a small +3 % *recovery*, not decay** — so a continuity chain gives the ATP effect ×1.30–1.32 from one prep with no φ assumed. **2 mM runs down ×1.95 faster and peaks ~2 s earlier**, so a fixed protocol window silently under-reads it; referencing each run to its own peak fixes that for free. Mechanism test against the model: fewer heads, lost SL, uniform slowdown and reduced attachment each fail on ≥1 observable; **series-elastic creep (kSE ×0.65, SL −0.098 µm) reproduces force ×0.83, `ktr` ×0.88 AND the length-tension bend**. Correcting `ktr` makes consistency worse — don't. Has a paper-ready [methods.md](RundownCorrection/methods.md). | [conclusions](RundownCorrection/conclusions.md) |
| **ATPEffectReconciliation** | findings | Applies the rundown correction to all four preps (three ran 8→2, 04/10 ran 2→8). Force CV **22 % → 11 %**, giving **×1.36**; `ktr` needs no correction and is already **×0.55 at CV 8 %**. **Low ATP = stronger and slower**, force×ktr ×0.74. Supersedes the earlier ×1.18/×1.23 figures, which were uncorrected or wall-clock-corrected. | [conclusions](ATPEffectReconciliation/conclusions.md) |
| **RestretchVsKtrRecovery** | findings | All three force rises (ktr, post-slack, post-restretch) share **one rate, ~41–47 s⁻¹**, across three protocol days — ktr is *not* a different process. Isolates three model defects: a titin-dashpot force floor unique to the ktr protocol, `tau_reg` sitting inside the ktr timescale, and a post-restretch recovery ~2× too fast. | [conclusions](RestretchVsKtrRecovery/conclusions.md) |
| **ApoPoolDetachment** | **falsified** | Detachment ceiling + tearing into a nucleotide-free pool (`UseR2Ceiling`, `UseD0State`, both default-OFF, wired and regression-clean). Reaches the measured post-restretch rate, but the **ktr overstretch peak** (data 0.72–0.83 F_iso; uncapped model 0.746, every capped variant 1.2–2.2) falsifies the ceiling on data already in the cost. ATP justifications also fail by 30–250×. Verdict: do not refit; the driver is kept as a falsification harness. | [conclusions](ApoPoolDetachment/conclusions.md) |

## Summaries & recommendations (across analyses)

**The model's flexibility is bounded, not infinite.** The sensitivity analysis shows an
**11-dimensional identifiability ceiling** from the Slack/FV/Ktr feature set — adding
parameters does not add constraint. The strongest global levers are serial-element
nonlinearity (`ekSE`), structural stiffness (`estiff`), transition geometry (`dr2`/`dr`),
and ADP-release/hydrolysis kinetics (`k2`/`kah`). Fitting effort should concentrate on those
directions rather than on widening the parameter set.

**Attachment is the recurring theme.** AttachedPool_vs_ATPase and BindingSiteOccupancy both
converge on the same lesson: the apparent "low attached pool / paradoxical turnover" tension
is partly a *mis-specified target* and partly a *mis-grounded saturation mechanism* (the
strain axis is not a site axis). The concrete, fixable defect is a throttled hydrolysis step —
not a deep contradiction. The labdiary's "set p1+p2 to ~30%" goal should be pursued through
the kinetics that the sensitivity analysis flags, not through the existing per-bin occupancy flag.

**Restretch shape is a kSE / attachment-timing problem.** RestretchMechanisms and
LastSlackIdentification target the same protocol. The Phase-1 fit reproduces peak2/valley but
overshoots peak1; the recommended path is to raise `kSE` (3–5×), let the A2 attachment shift
fire at threshold strain, then tune `slope`/`s_threshold_R` for valley/second-peak sharpness.

**Recovery rate is one number, and the model splits it into three.** RestretchVsKtrRecovery shows
the data has a *single* force-redevelopment rate (~41–47 s⁻¹) shared by the ktr protocol, the
post-slack redevelopment and the post-restretch recovery, on all three protocol days. The model gets
ratios of 3.16 / 1.40 / 0.35 instead of 1.0 — so protocol-independence of the redevelopment rate is
a sharp, cheap test that the current cost function does not score at all. Two of the three causes are
fixable by flags/parameters already present; the third (post-restretch ~2× too fast) needs new
mechanism, and ApoPoolDetachment records the leading candidate together with the reasons to distrust
its stated physiology.

**The ATP effect is now measured, not estimated.** ATPEffectReconciliation settles the target:
**force ×1.36, `ktr` ×0.55**, consistent across four preps and both ATP orders once rundown is
removed. Low ATP makes the muscle *stronger and slower* (force×`ktr` = ×0.74). This supersedes
the ×1.18 in LowATP_ForceEnhancement (uncorrected 03/27 only) and the ×1.23 briefly quoted from
a wall-clock rundown model. What remains open is *mechanistic* identification — reproducing
both numbers from cross-bridge kinetics rather than fitting them (see
[`../Docs/ROADMAP.md`](../Docs/ROADMAP.md)); LowATP_k2Frontier shows a single `k2` cannot do it.

**Rundown is a mechanical lesion, and it is identifiable.** RundownCorrection shows the
preparation degrades by **series-elastic creep** (~35 % softer, ~0.1 µm longer), not by losing
myosin heads — the length–tension curve *bends*, and a pure force scale cannot bend it. Two
practical consequences: correct force on **effective activated time** (φ ≈ 0.45 of each run's
within-run decline is permanent), and **never correct `ktr`** — doing so degrades cross-prep
consistency. When representing rundown in the model, prefer a per-run *length offset* to a
per-run *force scale*.

**Methodology that's settled.** PassiveForceID justifies direct subtraction of the high-Ca
passive trace to recover the cross-bridge-only force — treat this as the standard passive
correction in downstream fits.
