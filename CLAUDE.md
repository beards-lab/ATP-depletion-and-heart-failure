# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

MATLAB codebase implementing a cardiac cross-bridge (sarcomere) model to study ATP depletion and heart failure. Based on [Beard et al. 2022](https://www.sciencedirect.com/science/article/pii/S0006349522006026). The model simulates strain-discretized probability distributions of myosin cross-bridge states and their time evolution via ODEs.

## Running MATLAB

**Use the MATLAB MCP server — it is the primary, reliable execution path.** The MCP tools run against
`C:\Program Files\MATLAB\R2023a` in `nodesktop` mode and keep the MATLAB session alive between calls:

- `mcp__matlab__run_matlab_file` — run a `.m` script/function file (preferred for anything non-trivial).
- `mcp__matlab__evaluate_matlab_code` — run inline commands / snippets.
- `mcp__matlab__check_matlab_code` — static analysis of a `.m` file.
- `mcp__matlab__run_matlab_test_file` — run a MATLAB unit-test file.

The MCP does **not** guarantee a persistent path/workspace across calls, and its parpool workers cache
compiled code — so edits to model files are silently ignored until the pool is refreshed. At the start of a
fresh session, and after editing any model/param file, run:

```matlab
cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.'));
refreshPool(5);   % clears stale compiled code + (re)starts a fresh pool of 5 workers
```

Load parameter snapshots with `loadParams` — **never** `run`/`eval` a snapshot in the base workspace, or a
stale `params0` leaks fields and default values silently override the snapshot:

```matlab
params0 = getParams(loadParams('params_2state_seed'), [], true, false);
```

Both helpers live in `Auxiliary/`.

**Shell fallback only (the Bash tool is unreliable in this environment — prefer the MCP).** MATLAB can also
be launched non-interactively from a shell; startup takes ~10–20 s and `-batch` exits when done:

```bash
"/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe" -batch "cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.')); your_command"
```

R2025b is also installed; R2023a is the configured/confirmed version.

## Running the Model

All code runs in MATLAB. There is no build step.

**Demo run** (simplest entry point):
```matlab
% Open Workbench/RunModel.m and run it, or manually:
params = getParams([], g);
[force, out] = evaluateModel(@dPUdTCa, [0 1], params);
animateStateProbabilities(out, params);
```

**Full experimental evaluation** (compare model vs. Baker lab data):
```matlab
% Set up params0, then run:
RunBakersExp  % script in Model/, uses params0 from workspace
```

**Optimization** (Lakes optimizer):
```matlab
% Workbench/RunOptimLakes.m
```

**Handtune** (interactive parameter adjustment with plots):
```matlab
% Workbench/HandtuneAlternativeModel.m
```

**Evaluate fit quality** (returns scalar error) — current path:
```matlab
[Et, E, featCost, featNames] = evaluateBakersExp(g, params0, limitNegative, plotDebug)
% Sets params0.g = g, calls LoadData + RunBakersExp, and returns the scalar
% objective. When params0.OptimizeOn == 'Feats', Et is the FEATURE cost alone
% (sum of evalFeatureCost with costExp=2); otherwise Et = sum(E).
% featCost/featNames are computed only when nargout >= 3.
```

Legacy path (does **not** call `RunBakersExp`; builds its own `params0` and calls
`evaluateModel` directly). Only `Workbench/RunOptimLakes.m` still uses it:
```matlab
[Etot, E1] = evaluateProblem(fcn, g, true, [1 1 1 1])
% 4th arg is a binary vector selecting which experiment types to include
```

## Naming Conventions

- **`FirstCapital.m`** — script; runs and leaves variables in the workspace
- **`firstLowercase.m`** — function; has explicit return value(s)

## Directory Structure

> Reorganized 2026-06-16. Two front doors to the project's knowledge:
> **`Docs/README.md`** (knowledge map) and **`Analyses/README.md`** (analyses index + synthesis).
> See `Docs/reorganization-plan.md` for the rationale.

- **`Example/`** — the onboarding package; the canonical "how to run the model" answer.
  Start at `Example/README.md`. Heavily commented, verified end-to-end:
  - `RunExampleHighATP.m` — high-ATP (8 mM) battery: slack + ktr + stairs + FV,
    on `params/ExampleHighATP.m`. Reference feature cost **7.50** on the 03/27 8 mM set.
  - `BuildExampleProtocol.m` / `RunExampleProtocol.m` — build and run a custom
    length protocol (hold–stretch–hold), with a lossless round-trip self-test.
  - `params/ExampleHighATP.m` is FROZEN (verbatim copy of `rskR2_w025_opt.m`);
    the name is chosen so optimizer output (`<tag>_opt.m`, `<tag>_iter/`) can never collide.
- **`Model/`** — core model functions ONLY (model code, nothing else):
  - `Driver.m` — 700-line historical entry point; prefer `Example/` for onboarding
  - `getParams.m` — central parameter management; call this to build/update any params struct
  - `evaluateModel.m` — runs the ODE integrator for one experimental condition
  - `RunBakersExp.m` — script that evaluates all experimental challenges (FV, Ktr, slack, stairs) and computes total error `E`
  - `dPUdT_CombinedTransitions.m` — **primary ODE function** (current default); set via `params0.modelFcn`
  - `dPUdTCaSimpleAlternative2State.m`, `dPUdTCa.m` — alternative ODE variants (still present)
  - `resolveParams.m` — resolves param fields that are expressions (prefixed with `'='`)
  - `LoadBakersExp.m` — loads and plots experimental data
  - `extractSlackAttributes.m`, `extractForceVelocityAttributes.m`, `extractPerturbAttributes.m` — feature extraction from simulation output (called by experiment runners)
  - `evalFeatureCost.m` — feature-based cost function
  - `updateRates.m` — scales all turnover rates by a single `xrate` multiplier
  - _Archived:_ dead ODE variants `dPUdTCaSimple`, `dPUdTCaSimpleAlternative`, `dPUdT_TransitionRates` moved to `_archive/dead_ode/` (no callers).
- **`DataCuration/`** — raw data → model-ready data; all scripts use `../data/` paths:
  - See README for pipeline order (Baker legacy vs. 2026 driving-signal experiments)
- **`Analyses/`** — one self-contained folder per analysis (script(s) + `results/` + `conclusions.md`).
  Start at `Analyses/README.md`. Replaces the former singular `Analysis/` folder. Topics:
  `AttachedPool_vs_ATPase`, `BindingSiteOccupancy`, `RestretchMechanisms`, `SensitivityAnalysis`,
  `LastSlackIdentification`, `PassiveForceID`, `LowATP_ForceEnhancement`.
- **`Workbench/`** — the playground: driver scripts, hand-tuning, optimization, ad-hoc tests:
  - `RunOptim.m`, `RunOptimLakes.m` — optimization entry points
  - `HandtuneAlternativeModel.m` — manual parameter tuning workflow
  - `tunableParams.m` — lists which params are eligible for optimization
  - `DriverSimple.m`, `DriverSimple_Beard2022.m`, etc. — variant entry-point scripts (playground)
  - `Tests/` — `Test*.m` validation scripts
- **`params/`** — active/current parameter snapshots (the working source of truth); also CSV/XLSX/TXT snapshots.
- **`ModelOptParams/`** — historical named param gallery + QC plots (archive of past optima). Kept separate from `params/` by design.
- **`Auxiliary/`** — generic reusable, non-critical utilities:
  - `animateStateProbabilities.m`, `plotStateFluxes.m`, `StatesInTime.m` — visualization
  - `calcSensitivities.m`, `ResidualAndJacobian.m` — sensitivity analysis tools
  - `fitRecovery.m`, `fitSlackForceOnset.m` — fitting utilities
  - `writeParamsToMFile.m` — serializes a params struct to a `.m` file
  - `diagnostics/` — reusable diagnostic probes (`srxProbe`, `ktrProbe`, `mechProbe`, `isoProbe`)
- **`data/`** — experimental `.txt`/`.mat` data files; not tracked by git
- **`Docs/`** — knowledge base (start at `Docs/README.md`): `reference/` (durable), `notes/` (working),
  `experiments/` (dated write-ups), `presentations/` (slides/posters), `ROADMAP.md`, and process docs.
- **`scripts/`** — infrastructure (`setup-gh.sh`, `rol.sbatch`).
- **`_archive/`** — disk-only graveyard for stale `.mat`/`.png` dumps, old optimizer envs, dead code (gitignored; "archive, don't delete").

## Architecture: Parameter System

Everything flows through the `params` struct. `getParams(params, g, updateInit, updateModifiers)` is the single source of truth:

1. Fills missing fields with defaults
2. Applies **modifiers** `g` (multiplicative scalars) to fields listed in `params.mods`
3. Calls `resolveParams` — resolves fields whose value is a string starting with `'='` as MATLAB expressions (e.g., `params.kamh = '=0.1*kah'`)
4. Reconstructs array elements from `fieldname__index` notation (e.g., `PieceWiseStrainDepParams__3`)
5. Computes the strain grid `params.s` from `Slim_l`, `Slim_r`, `dS`
6. Initializes the state vector `params.PU0`

**Always call `getParams` again after changing `SL0`** — it affects the initial state vector.

## Architecture: ODE State Vector

The state vector `PU` inside `dPUdT_CombinedTransitions` is structured as:
```
[p1(1..ss), p2(1..ss), [p3(1..ss),] P_SR, NP, SL, LSE, PD, [P_SRD,] [x_dash]]
```
where `ss = length(params.s)` is the number of strain bins.

## Key Behavioral Switches in params0

Active switches (see README for current active set):
- `UseSuperRelaxed`, `UseSuperRelaxedADP` — super-relaxed myosin head states
- `UseOverlap`, `UseOverlapFactor` — filament overlap corrections
- `UsePassive`, `UseTitinIdentifiedPassive` — passive titin-based force
- `UseSerialStiffness` — serial elastic element
- `UseSlack` — enable XB slacking protocol
- `UseDirectSRXTransition`, `UsePassiveForSR` — SRX transition variants
- `UseSpaceExtension` — dynamic expansion of strain grid if boundary is hit
- `UseA1AttachmentKernel` — distributed attachment probability
- `NumberOfStates` — 2 or 3 cross-bridge states

## Optimization Workflow

Parameters are optimized as multiplicative modifiers `g` applied to base values in `getParams`. To optimize specific params:
```matlab
params0.mods = {'kstiff2', 'ka', 'k1'};  % params to optimize
% g(i) multiplies params0.mods{i}
optimfun = @(g) evaluateProblem(fcn, g, false, [1 1 0 1]);
x = fminsearch(optimfun, g0, options);
writematrix(x, 'gopt.csv');  % save result
```

`evaluateProblem` internally sets `params0`, calls `RunBakersExp`, and returns total weighted error.

The alternative **feature-based** cost uses `params0.EvalFeatures = true` and `params0.fn` cell array of feature selectors (e.g., `'FV_f|FV_v'`, `'ktr|SLslack'`).
