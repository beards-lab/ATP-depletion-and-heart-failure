# Docs — start here

The single map of the project's knowledge. **What do you want?**

- **Run or understand the model** → repo-root [`README.md`](../README.md) and [`CLAUDE.md`](../CLAUDE.md).
- **What has been figured out (analyses)** → [`../Analyses/README.md`](../Analyses/README.md) —
  one folder per analysis (script + results + conclusion), with a cross-analysis synthesis.
- **What to do next** → [`ROADMAP.md`](ROADMAP.md) — the consolidated backlog (all the scattered
  TODOs pulled into one prioritized list).
- **Durable reference knowledge** → [`reference/`](reference/) (below).
- **Working notes / lab diary** → [`notes/`](notes/).
- **Dated experiment write-ups** → [`experiments/`](experiments/).
- **Slides, posters, abstracts** → [`presentations/`](presentations/) (old versions in `presentations/archive/`).

## `reference/` — durable knowledge you build on

| Doc | What it is |
|---|---|
| [`parameter-reference.md`](reference/parameter-reference.md) | Catalogue of model parameters and their meaning. |
| [`mechanism-evaluation.md`](reference/mechanism-evaluation.md) | The big (108 KB) literature/mechanism analysis. |
| [`fitting-strategy.md`](reference/fitting-strategy.md) | How fits are set up and prioritized. |
| [`RundownCorrection.md`](reference/RundownCorrection.md) | Correcting for preparation rundown (+ rendered PDF). |
| [`sources.md`](reference/sources.md) | Literature links and provenance. |
| [`update-mex.md`](reference/update-mex.md) | How-to: rebuild the MEX ODE core. |
| `Ensemble vs distributed elasticity model.xlsx` | Model-variant comparison table. |

## `notes/` — transient working notes

`labdiary.md`, `labdiary_boundedfit.md` — chronological; may go stale. TODOs from these have
been lifted into [`ROADMAP.md`](ROADMAP.md).

## `experiments/` — dated result write-ups

`results-0327.md` and its `figures/`, pairing with the `data/` recordings.

## Process docs

[`reorganization-plan.md`](reorganization-plan.md) and [`refactoring-plan.md`](refactoring-plan.md)
stay here as the record of how/why the repo is laid out the way it is.

> Note: analysis *conclusions* live with their analysis in `../Analyses/<Topic>/conclusions.md`,
> not here — `reference/` holds only cross-cutting, durable knowledge, so there is exactly one
> home for each finding.
