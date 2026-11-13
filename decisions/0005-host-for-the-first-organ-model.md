# 0005 — The host for the first organ model

- **Status:** proposed, for the PI (revised 2026-11-13 after the mouse host was examined)
- **Date:** 2026-11-12, revised 2026-11-13
- **Issue:** #82 (milestone 5)
- **Drafted by:** Claude (AI agent). **Decided by:** —

## Context

The cross-bridge model is calibrated on skinned mouse trabeculae at rig temperature (not
recorded, #43). Two TriSeg hosts are available:

- the **P21 rat** model at 37 °C (`../Multiscale-Hx-P21-Model`; Kim et al., bioRxiv 2025,
  doi:10.1101/2025.10.01.679825), public, 19 animals fitted to pressure–volume loops, with a
  cross-bridge sarcomere and a Ca transient of its own;
- the lab's **mouse** model (`beards-lab/Rat_Triseg`, private), the human TriSeg re-scaled to
  mouse, with a Lumens-type sarcomere driven by an activation function of time, fitted to a
  literature-target mouse and to one mouse's echo data.

`Analyses/TriSegCaPlan/report.md` (Q5, Q7, Q7b, Q9) measured both.

## Decision (proposed)

1. **The mouse host is the target of the organ model**: validation (#99) and the native passive
   (#83) are done on it.
2. **The P21 rat host stays the numerical test bed** for the harness (#97), the velocity signal
   (#84) and the Jacobian (#54), because it already has the thin filament, a Ca transient and
   stored beats to test against.
3. **Nothing the coupled model produces on either host is a prediction** until the kinetics
   package (#98) lands, and on the mouse host until it has a Ca transient and the activation
   amplitude is set (#41, #42).
4. **Publishing on the mouse host needs its repository released**, or its equations and
   parameters restated. To settle before any figure from it leaves the lab.

## Evidence

| | P21 rat (Nx51) | Mouse, literature targets | Mouse, one echo fit |
|---|---|---|---|
| Heart rate (bpm) | 325 | 409 | 579 |
| LV sarcomere length, max / min (µm) | 2.25 / 1.97 | 2.04 / 1.74 | 2.03 / 1.72 |
| LV peak shortening speed (lengths/s) | 12.2, a spike | 4.9, a 40 ms plateau | 4.9, a plateau |
| Incumbent static passive at the host's end-diastolic length (kPa) | 23.6 | 0.2 | 0.2 |
| Host passive there (kPa) | 0.83 | 0.44 | 1.71 |
| Model slack on the host's beat (ms) | 10.5 | 30.0 | 45.5 |
| Ca transient in the host | yes | no | no |
| Data behind the fit | PV loops, 19 animals | literature | echo, pressure assumed |

- The static passive mismatch that blocks the P21 host is absent on the mouse host, because its
  end-diastolic sarcomere length is 2.04 µm, not 2.25 µm.
- The speed mismatch is on both. The model's unloaded speed is 4.6 ML/s; the mouse host ejects
  at that speed for 40 ms, so the model carries no force during ejection, also with its cycle
  rates scaled by 5 (a Q10 of 3.5 from 22 to 37 °C; de Tombe and Stienen 2007).
- On the mouse host the model's velocity-dependent passive gives 4–9 kPa in fast filling, where
  the host is near zero (#85).

## Alternatives considered

- **P21 rat for everything** (the first version of this draft): wrong species and age, sarcomere
  lengths out of the calibrated range, and a host re-scaling needed before the native passive
  can be used.
- **Mouse host for everything:** the harness would have to wait for a mouse Ca transient, and the
  parity test against a stored cross-bridge beat would be lost.
- **`../CVS_TriSeg`:** a human clinical model; it fits the humanised model (milestone 2).

## Consequences

- A mouse Ca transient at 400–600 bpm is needed (#76, #49, or a literature transient).
- The validation tolerance has to be re-derived for the mouse host: it has one fitted animal and
  no measured pressures, so the P21 tolerances (EDV 8 %, ESV 6 %, peak pressure 10 %) do not
  carry over, and #99 has to say what data it validates against.
- #83 shrinks: little or no host re-scaling at an end-diastolic length of 2.04 µm.
- #98 (kinetics) is the main blocker on either host.
