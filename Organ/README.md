# Organ: the cross-bridge model in a TriSeg heart

The full strain-distributed cross-bridge model as the sarcomere of the three walls of a
TriSeg heart with a lumped circulation, solved as one system (decision
[0003](../decisions/0003-triseg-coupling-architecture.md)). Issue #97. Nothing in `Model/`
is changed or copied.

**Status:** a harness. It runs on the P21 rat host with the host's passive law and a
calibrated activation. Its pressures and volumes are not predictions (decision
[0005](../decisions/0005-host-for-the-first-organ-model.md)).

## Which part comes from where

| File | Provenance | Source |
|---|---|---|
| `hostP21Geometry.m` | **HOST**, copied | `beards-lab/Multiscale-Hx-P21-Model`, `model/model.m` @ 80f86d3, l.219–246 |
| `hostP21Balance.m` | **HOST**, copied; the wall stress is an input instead of `-Kse*(Lsc - Ls)` | same file, l.296–327, 356–362, 431–434, 487–492 |
| `hostP21Ca.m` | **HOST**, copied; amplitude and diastolic level read from the stored fit | same file, l.195–211 |
| `hostP21Load.m` | **HOST data**: parameters, starting state and reference traces of a saved fit | `outputs/outputs_<fit>.mat` of the host |
| `xbWall.m` | **XB**: calls `Model/dPUdT_CombinedTransitions.m`, unchanged. **GLUE**: state packing, Ca gate, velocity signal | Ca gate after `Analyses/CalciumTwitch/dPUdT_CaGate.m`; thin-filament constants from the host |
| `setupXbWall.m` | **GLUE**: a `params/` snapshot through the model's `getParams`; the host's passive law placed in the model's own passive slot (decision [0004](../decisions/0004-passive-law-for-the-organ-model.md)) | host `model.m` l.262–264 for the passive law |
| `shiftXbWall.m` | **GLUE**, reproducing the strain-window slide of `Model/evaluateModel.m` | |
| `trisegXB_rhs.m` | **GLUE**: assembles host and walls into one right-hand side | |
| `setupTrisegXB.m`, `runTrisegXB.m`, `runXbWallOpenLoop.m` | **GLUE**: initial state, mass matrix, Jacobian pattern, integration, per-beat summary | steady-state idea after the host's `model_sol.m` |

Every file states its provenance in its header, and the copied blocks are fenced by
`% ---- model.m l.xxx ----` comments.

## When something upstream changes

- **The cross-bridge model changes** (`Model/`, `params/`): nothing to do here. `xbWall`
  calls the model, and `setupXbWall` reads the state layout from `getParams`. Re-run the
  tests. If the state vector's `SL` and `LSE` positions change (today `Ns*ss + 3`,
  `Ns*ss + 4`), update `setupXbWall.m`.
- **The host changes**: re-copy the fenced blocks from the new `model.m`, update the commit
  hash in the headers and in `hostP21Load.m`, and re-run part 0 of the parity test, which
  compares the host block with the host's own stored outputs.
- **Another host**: write its own `host<Name>*.m` block with the same three roles
  (geometry → sarcomere lengths; Ca(t); balance and circulation given the wall stresses).

## Run

```matlab
addpath(genpath('Model')); addpath(genpath('Auxiliary')); addpath('params'); addpath('Organ');
H = hostP21Load('Nx51');             % host repository at ../Multiscale-Hx-P21-Model
C = setupTrisegXB(H);                % defaults: robust1_seed, 60 bins, k_on x5.3, host passive
R = runTrisegXB(C, 10);              % 10 beats; R.beat is the per-beat table
```

`Analyses/TriSegHarness/RunTrisegHarness.m` runs to a periodic state and prints the
comparison with the host fit.

## Tests

| Test | What it checks | Tolerance | Measured |
|---|---|---|---|
| `Tests/TestOrganParity.m`, part 0 | host block against the host's stored outputs | strain 1e-8, pressures 0.2 % | strain exact, pressures 9e-4 kPa |
| part A | `xbWall` against `Model/evaluateModel` (gate open, native passive) | 0.2 kPa | 0.03 kPa |
| part B | `xbWall` against the open-loop probe of `Analyses/TriSegCaPlan` on the host beat | 0.2 kPa | 0.04 kPa |
| `Tests/TestTrisegClosedLoop.m` | the coupled system starts and runs; blood volume conserved | | |

## State layout

`y = [g (4); V (6); wall LV; wall SEP; wall RV]`

- `g = [xm_LV; xm_SEP; xm_RV; ym]` (cm): algebraic, the host's rows 1–4.
- `V = [V_LV; V_RV; V_SV; V_PV; V_SA; V_PA]` (cm³): the host's states 8–13.
- a wall: `[core; Lsc; Lf; N]`: the model's own state vector (its `SL` and `LSE` slots are
  unused), the contractile length `Lsc`, the low-passed sarcomere length `Lf` behind the
  velocity signal, and the thin-filament state `N`.

The host's states 5–7 (`Lsc`) are the walls' `Lsc`; its states 14–46 (cross-bridge
moments, `N`, `U_NR`) are replaced by the walls.
