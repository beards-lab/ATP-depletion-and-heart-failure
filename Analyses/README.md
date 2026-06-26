# Analyses

One self-contained folder per analysis: the **script(s)** that run it, a **`results/`**
folder for its `.mat`/figure outputs, and a **`conclusions.md`** that records what it found
and what to do next. This page is the index and the cross-analysis synthesis.

> Convention: an *analysis* lives here (script + results + conclusion). The **playground**
> (new drivers, hand-tuning, optimization entry points, tests) lives in `Workbench/`.
> Reusable, non-critical helpers live in `Auxiliary/` (diagnostic probes in
> `Auxiliary/diagnostics/`). Core model code lives in `Model/`.

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
| **LowATP_k2Frontier** | findings | Tests whether reducing `k2` alone reproduces the 2-vs-8 mM signature (relative-ratio scoring). It owns the kinetics (ktr/peak2 at k2×0.40) but over-produces isometric force ~2×; needs a second amplitude lever. | [conclusions](LowATP_k2Frontier/conclusions.md) |

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

**ATP-dependence is the open scientific question.** LowATP_ForceEnhancement documents a real,
rundown-corrected force *increase* at low ATP — the phenomenon the whole project exists to
explain. Mechanistic identification of the ATP-dependent differences is the headline next step
(see [`../Docs/ROADMAP.md`](../Docs/ROADMAP.md)).

**Methodology that's settled.** PassiveForceID justifies direct subtraction of the high-Ca
passive trace to recover the cross-bridge-only force — treat this as the standard passive
correction in downstream fits.
