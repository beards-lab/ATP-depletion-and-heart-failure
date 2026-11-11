# Last-Slack Identification — Conclusions

**Status:** STUB — to be filled in.

## Scope

Identify the model parameters that reproduce the *last-slack* restretch protocol
(the slack release followed by restretch at the end of the isovelocity series),
using a piecewise strain-dependent detachment description.

## Scripts

- `OptimLastSlackPieceWise.m` — workhorse optimizer for the piecewise strain-dependent
  detachment parameters against the last-slack transient.
- `IdentifyLastSlack.m` — extraction/identification of the slack-onset and restretch
  features used as fit targets.

## Results

Saved in `results/`:

- `OptimLastSlackPieceWise_state.mat`, `OptimLastSlackPieceWise_stateTimeCourse.mat`,
  `OptimLastSlackPieceWise_stateTimeCourseAll.mat` — optimizer state / time-course snapshots.
- `IdentifyLastSlack_state.mat` — identification state.

## Key findings

_To be written._ See the [RestretchMechanisms](../RestretchMechanisms/conclusions.md)
analysis, which builds on the same protocol, and the consolidated next-steps in
[`Docs/ROADMAP.md`](../../Docs/ROADMAP.md).
