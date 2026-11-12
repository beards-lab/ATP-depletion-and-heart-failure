# 0003 — How the full cross-bridge model is coupled to TriSeg

- **Status:** proposed, for the PI
- **Date:** 2026-11-12
- **Issue:** #82 (milestone 5)
- **Drafted by:** Claude (AI agent). **Decided by:** —

## Context

The full strain-distributed cross-bridge model is to be the sarcomere of all three TriSeg walls,
closed loop with the circulation. `Analyses/TriSegCaPlan/report.md` (Q1, Q2) measured what the
options cost on a host beat (P21 rat fit Nx51, open loop).

## Decision (proposed)

1. **One monolithic DAE.** The host's circulation and geometry states plus three copies of the
   wall state, integrated by one `ode15s` call. A wrapper hands each wall its sarcomere length
   from the geometry and gets the contractile-element velocity back (`dLsc/dt = velHS`).
2. **The model files stay unchanged, and there is no `UseTriSeg` flag.** Coupling is a way of
   using the model, like `evaluateModel` is. The wrapper lives outside `Model/`.
3. **The velocity-driven mechanisms get dLs/dt**, not the contractile-element velocity (#84).
4. **The grid stays at 60 bins with the sliding window** until the coupled system shows a reason
   to widen it. The first cost optimisation is the structured Jacobian (#54).
5. **Parity tests guard the wrapper.** (a) Wrapper against `evaluateModel` on the same length
   protocol with the Ca gate open: within 0.2 kPa (measured 0.08 kPa). (b) The closed-loop
   wrapper driven by a recorded Ls(t) against the open-loop probe: within 0.2 kPa.

## Evidence

- Wrapper against `evaluateModel`: maximum difference 0.083 kPa on 33.5–72.2 kPa.
- A lag of 0.1 ms in the host's length changes the wall force by up to 2.5 kPa, 1 ms by 24.7 kPa.
  An explicit exchange between separate solvers would need macro-steps of 0.1 ms or less.
- The solver takes about 3,000 steps per beat whatever the grid, and a smoother host signal does
  not change that.
- Force is the same to 0.01 kPa at 60, 120 and 250 bins; slides lose 2·10⁻¹² of attached mass.
- 63 % of the right-hand-side calls are finite-difference Jacobian columns at 60 bins (86 % at
  250), because 25 rows are dense. Without them the Jacobian needs 4 column groups, not 129.

## Alternatives considered

- **Explicit co-simulation** reusing `evaluateModel` per wall: no saving at the macro-step the
  stiff series element demands, and a new error source.
- **A split with the cross-bridge stiffness at the interface:** removes the stiff loop from the
  exchange but is first-order in the macro-step, against rates up to 1,500 s⁻¹.
- **A wide grid from the start:** not needed for correctness on this beat, and 1.5 to 2.8 times
  the right-hand-side calls until #54 lands.

## Consequences

- A slide in one wall restarts the whole system, about 7 times per beat at 60 bins.
- The cost per closed-loop run is about 60 wall-beats (20 beats to a periodic state, 3 walls).
- If #52 (head deletion on a slide) appears on another host beat, point 4 is revisited.
