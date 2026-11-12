# 0004 — The passive law of the walls in the organ model

- **Status:** proposed, for the PI
- **Date:** 2026-11-12
- **Issue:** #82 (milestone 5); work in #85 and #83
- **Drafted by:** Claude (AI agent). **Decided by:** —

## Context

The incumbent passive law, `k_pas·(Lsc − 1.80)^7.69`, was fitted on the slack battery at high Ca.
The host has its own law, `6.67·(Ls − 1.51)^7`, fitted to pressure and volume with an assumed
sarcomere length scale. `Analyses/TriSegCaPlan/report.md` (Q3) compared both with the passive
data of seven skinned mouse trabeculae.

## Decision (proposed)

1. **The first closed loop uses the host's passive law.** It lets the harness (WP1) be finished
   and tested without the passive question.
2. **The incumbent law is not used in the organ model as it stands.**
3. **The native passive is a Ca-dependent viscoelastic law** (#85), and the host is re-scaled to
   it (#83). A Ca-dependent static law is not pursued.

## Evidence

- At the host's end-diastolic length (2.25 µm) the incumbent law gives 23.6 kPa and the host
  0.81 kPa. The measured relaxed static passive there is about 5 kPa.
- The static passive does not depend on Ca: 3.16 kPa relaxed and 3.41 kPa at pCa 4.4 at 2.2 µm on
  100 s ramps, 5.4 and 5.5 kPa 60 s into a hold. The incumbent law gives 9.67 kPa at 2.2 µm and
  54 kPa at 2.3 µm.
- The fast response does: peak stress on a 0.1 s ramp is 10.5 kPa relaxed and 41.0 kPa at pCa 4.4
  (×4.0, range 2.9–5.9 over preparations); on a 1 s ramp ×2.9; on a 100 s ramp ×1.07.
- The host's filling stretches the LV at up to 2.1 lengths/s, the rate of the 0.1 s ramp.

## Alternatives considered

- **The measured relaxed static law plus host re-scaling:** right in slow diastole, but no
  viscosity and no Ca dependence, so wrong during filling and in systole.
- **The incumbent structure with Ca-dependent coefficients:** an option inside #85 (its form c),
  to be judged there on parity against the titin model.

## Consequences

- The first closed-loop results carry the host's diastole, not the model's.
- The 8 mM battery has to be re-checked when the passive changes (#85 acceptance), because the
  incumbent's passive part is read by the slack battery.
- End-diastolic sarcomere length becomes a result to check against the in-vivo range (rat,
  inferred 1.75–1.90 µm; Pearson et al. 2007) once the host is re-scaled.
