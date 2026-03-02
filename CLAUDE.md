# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

MATLAB codebase implementing a cardiac cross-bridge (sarcomere) model to study ATP depletion and heart failure. Based on [Beard et al. 2022](https://www.sciencedirect.com/science/article/pii/S0006349522006026). The model simulates strain-discretized probability distributions of myosin cross-bridge states and their time evolution via ODEs.

## Running MATLAB

MATLAB is installed on Windows and accessible from WSL via interop. Use this command:

```bash
"/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe" -batch "your_matlab_code_here"
```

To run a script with the repo on the MATLAB path:
```bash
"/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe" -batch "cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.')); your_script_or_command"
```

Note: MATLAB startup takes ~10–20 seconds. `-batch` mode runs non-interactively and exits when done. R2025b is also installed but R2023a is confirmed working.

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

**Evaluate fit quality** (returns scalar error):
```matlab
[Etot, E1] = evaluateProblem(fcn, g, true, [1 1 1 1])
% 4th arg is binary vector selecting which experiment types to include
```

## Naming Conventions

- **`FirstCapital.m`** — script; runs and leaves variables in the workspace
- **`firstLowercase.m`** — function; has explicit return value(s)

## Directory Structure

- **`Model/`** — core model functions:
  - `getParams.m` — central parameter management; call this to build/update any params struct
  - `evaluateModel.m` — runs the ODE integrator for one experimental condition
  - `RunBakersExp.m` — script that evaluates all experimental challenges (FV, Ktr, slack, stairs) and computes total error `E`
  - `dPUdT_CombinedTransitions.m` — **primary ODE function** (current default); set via `params0.modelFcn`
  - `dPUdTCaSimpleAlternative2State.m`, `dPUdTCa.m`, etc. — alternative ODE variants
  - `resolveParams.m` — resolves param fields that are expressions (prefixed with `'='`)
  - `LoadBakersExp.m` — loads and plots experimental data
- **`Workbench/`** — driver scripts and optimization:
  - `RunOptim.m`, `RunOptimLakes.m` — optimization entry points
  - `HandtuneAlternativeModel.m` — manual parameter tuning workflow
  - `tunableParams.m` — lists which params are eligible for optimization
  - `ModelParamsInitNiceSlack_prescribedSR.m` (and variants) — complete param snapshots saved as scripts
- **`Auxiliary/`** — utility functions:
  - `animateStateProbabilities.m`, `plotStateFluxes.m`, `StatesInTime.m` — visualization
  - `calcSensitivities.m` — one-at-a-time parameter sensitivity analysis
  - `extractSlackAttributes.m`, `extractForceVelocityAttributes.m` — feature extraction from simulation output
  - `evalFeatureCost.m` — feature-based cost function
  - `updateRates.m` — scales all turnover rates by a single `xrate` multiplier
  - `writeParamsToMFile.m` — serializes a params struct to a `.m` file
- **`data/`** — experimental `.txt` data files (ATP concentrations, slack, Ktr, force-velocity)
- **`Thrash/`** — obsolete/experimental code; do not modify
- **`Docs/`** — presentations and model description documents

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
