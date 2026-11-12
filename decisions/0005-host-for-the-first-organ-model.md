# 0005 — The host for the first organ model

- **Status:** proposed, for the PI
- **Date:** 2026-11-12
- **Issue:** #82 (milestone 5)
- **Drafted by:** Claude (AI agent). **Decided by:** —

## Context

The cross-bridge model is calibrated on skinned mouse trabeculae at rig temperature (not
recorded, #43). The host on disk is the P21 rat TriSeg model at 37 °C
(`../Multiscale-Hx-P21-Model`; Kim et al., bioRxiv 2025, doi:10.1101/2025.10.01.679825), with 19
fitted animals. `Analyses/TriSegCaPlan/report.md` (Q5, Q7, Q9) measured the mismatch.

## Decision (proposed)

1. **The P21 rat host is used for the harness and the numerics** (WP1, WP2, WP6).
2. **Nothing the coupled model produces on it is a prediction** until three packages land:
   kinetics at the host's temperature and rate (WP3), the native passive (#85, #83), and the
   activation amplitude from force–pCa (#41, #42).
3. **The host is revisited** when the rig temperature and species are confirmed (#43) and the Q10
   block exists (#44).

## Evidence

- The host's LV sarcomere shortens at up to 27.6 µm/s (12 lengths/s) at the start of ejection;
  6.2–16.9 lengths/s across the 19 fits. The model's unloaded speed is 4.6 ML/s. On the host's
  kinematics the model's force collapses within 3 ms and it is slack for 10.5 ms.
- Scaling the cycle rates by 5 (a Q10 of 3.5 from 22 to 37 °C; de Tombe and Stienen 2007)
  removes the slack but not the collapse.
- Heart rate 292–388 bpm; sarcomere lengths from 2.33 µm down to 1.41 µm, below the calibrated
  range.
- The host has the same sarcomere interface, a steady-state routine, and data for validation:
  it reproduces its own EDV, ESV and peak pressure to a median 3–5 %.

## Alternatives considered

- **`../CVS_TriSeg`:** appears to be a human clinical model; it fits the humanised model
  (milestone 2), not the mouse data.
- **A mouse 37 °C TriSeg host:** none was found in a verified source.
- **Waiting for a matched host:** blocks the numerics, which do not depend on the match.

## Consequences

- The first coupled runs are a software result. The report and any figure from them say so.
- The validation tolerance comes from this host: EDV, ESV and peak pressure within twice the
  host's own median error (8 %, 6 %, 10 %); end-diastolic pressure is reported, not gated.
