# ATP Depletion and Heart Failure — Cross-Bridge Model

MATLAB codebase implementing a strain-discretized cardiac cross-bridge (sarcomere) model to study the effects of ATP depletion in heart failure. Based on [Beard et al. 2022](https://www.sciencedirect.com/science/article/pii/S0006349522006026).

---

## Naming Conventions

| Pattern | Meaning |
|---------|---------|
| `FirstCapital.m` | Script — runs and leaves variables in the workspace |
| `firstLowercase.m` | Function — has explicit return values |

---

## Directory Structure

```
ATP-depletion-and-heart-failure/
│
├── Model/              Core ODE engine, parameter system, experiment coordinators
│   └── experiments/    Specialized protocol runners (FV, ktr, slack, stairs)
│
├── Auxiliary/          Reusable utility functions: visualization, feature extraction,
│                       fitting, analysis helpers
│
├── Workbench/          Interactive driver scripts, optimization entry points
│   ├── DataCuration/   Velocity-table building, protocol construction, raw-data extraction
│   └── Tests/          Test and validation scripts for individual subsystems
│
├── Analysis/           High-level analysis and tuning campaigns (sensitivity, mechanism
│                       evaluation, restretch tuning, piecewise optimization)
│
├── ModelOptParams/     Auto-generated parameter snapshots from optimization runs (73+ files)
│
├── params/             Hand-curated and validated parameter sets
│
├── PassiveTitin/       Passive titin characterization scripts and Modelica model
│
├── Docs/               Presentations, meeting notes, model descriptions
├── Figures/            Generated and reference figures
└── data/               Experimental data (ATP levels, slack, ktr, force-velocity)
```

---

## Key Workflows

### 1. Demo — run the model once

```matlab
cd Workbench
RunModel          % loads data, builds params, runs ODE, plots output
```

### 2. Evaluate fit against Baker lab data

```matlab
% In any script after addpath(genpath('.'))
LoadData;                          % loads data into workspace
params0 = getParams();             % build parameters (all defaults)
ModelParamsInitNiceSlack_prescribedSR  % apply a saved param set
RunBakersExp                       % runs all active protocols, sets E (cost vector)
```

### 3. Interactive parameter tuning

```matlab
Workbench/HandtuneAlternativeModel   % loads data, lets you tweak params0 interactively
```

### 4. Run an optimisation

```matlab
% Set params0 and g0 first, then:
Workbench/RunOptimLakes    % Lakes optimizer (recommended)
Workbench/RunOptim         % fminsearch optimizer
```

### 5. Multi-dataset protocol comparison (slack, ktr, step, staircase)

```matlab
Workbench/CompareProtocols   % loads all datasets, fits all feature types, plots comparisons
```

---

## Core Execution Stack

```
RunBakersExp (coordinator script)
  └─ evaluateProblem()            wrapper: set params0, call RunBakersExp, return error
       └─ evaluateModel()         integrate ODE for one condition (ode15s)
            └─ dPUdT_CombinedTransitions()   PRIMARY ODE right-hand side
                 └─ getParams()   central parameter management (call after any SL0 change)
```

### ODE State Vector

```
PU = [p1(1..ss), p2(1..ss), [p3(1..ss),]  P_SR, NP, SL, LSE, PD, [P_SRD,] [x_dash]]
```

where `ss = numel(params.s)` = number of strain bins.

---

## Model/ — Core Files

| File | Role |
|------|------|
| `getParams.m` | **CENTRAL** — builds/updates the params struct; call after changing any field |
| `dPUdT_CombinedTransitions.m` | **PRIMARY ODE** — set via `params0.modelFcn`; all active mechanisms |
| `evaluateModel.m` | Integrates ODE for one experimental condition |
| `RunBakersExp.m` | Coordinator: runs FV, ktr, slack, staircase protocols and computes E |
| `evaluateProblem.m` | Thin wrapper for use by optimisers |
| `resolveParams.m` | Resolves `'='`-prefixed param fields as MATLAB expressions |
| `dPUdTCa*.m` | Legacy / alternative ODE variants (kept for reference) |
| `experiments/` | `runFVExperiment`, `runKtrExperiment`, `runSlackExperiment`, `runStairsExperiment` |

---

## Auxiliary/ — Utilities

| File | Role |
|------|------|
| `extractSlackAttributes.m` | Feature extraction from slack-release protocol |
| `extractForceVelocityAttributes.m` | Feature extraction from FV protocol |
| `extractPerturbAttributes.m` | Feature extraction from perturbation (ktr) protocol |
| `evalFeatureCost.m` | Feature-based cost function (data vs simulation) |
| `fitRecovery.m` | Fit force-recovery kinetics (used for ktr) |
| `plotFeatures.m` | Plot data-vs-sim features (wrapper around plotMultipleFeatures) |
| `plotMultipleFeatures.m` | Core multi-dataset feature comparison plots |
| `animateStateProbabilities.m` | Animate cross-bridge state probability distributions |
| `calcSensitivities.m` | One-at-a-time parameter sensitivity analysis |
| `updateRates.m` | Scale all turnover rates by `xrate` multiplier |
| `writeParamsToMFile.m` | Serialize params struct to a `.m` file |

---

## Workbench/ — Drivers

| File | Role |
|------|------|
| `RunModel.m` | Simplest demo: one ODE run + plots |
| `HandtuneAlternativeModel.m` | Interactive tuning loop |
| `RunOptim.m` / `RunOptimLakes.m` | Optimisation entry points |
| `Force_pCa.m` | Force-pCa curve at varying ATP concentrations |
| `CompareProtocols.m` | **Multi-dataset comparison**: slack, ktr, step, staircase across ATP levels |
| `LoadAndPlotLogs.m` | Load and visualise optimisation log files |
| `DataCuration/` | Velocity table construction, protocol building, raw log extraction |
| `Tests/` | Validation scripts for individual subsystems (`TestSR`, `TestFitRates`, etc.) |

---

## Analysis/ — Analysis & Tuning Campaigns

| File | Role |
|------|------|
| `TuneRestretch.m` | Investigate restretch double-peak discrepancy |
| `AnalyseRestretchMechanisms.m` | Systematic sweep of restretch mechanisms |
| `RunMechanismEvaluation.m` / `OptimizeMechanismEvaluation.m` | Evaluate top-K mechanism candidates |
| `RunSensitivityAnalysis.m` / `SensitivityAllSlack.m` | Parameter sensitivity analysis |
| `IdentifyLastSlack.m` / `OptimLastSlackPieceWise.m` | Identify model parameters from last-slack segment |
| `RunOptimPiecewise.m` | Piecewise block-coordinate-descent optimisation |
| `driverSA.m` / `runSA.m` | Sensitivity analysis workflow |

---

## params/ — Parameter Sets

Hand-curated and validated parameter configurations. Load a `.m` file to populate `params0` before running:

```matlab
params0 = getParams();
ModelParamsInitNiceSlack_prescribedSR    % apply saved param set
```

See `ModelOptParams/` for auto-generated snapshots from optimisation campaigns.

---

## Key Behavioural Switches in params0

| Switch | Description |
|--------|-------------|
| `UseSuperRelaxed` | Enable super-relaxed (SRX) myosin head state |
| `UseSuperRelaxedADP` | SRX with ADP-bound variant |
| `UseOverlap` / `UseOverlapFactor` | Filament overlap corrections |
| `UsePassive` / `UseTitinIdentifiedPassive` | Passive titin-based force |
| `UseSerialStiffness` | Serial elastic element |
| `UseSlack` | Enable slack-restretch protocol in simulation |
| `UseDirectSRXTransition` | Direct SRX↔IHM transition variant |
| `NumberOfStates` | 2 or 3 cross-bridge states |
| `EvalFeatures` | Use feature-based (vs time-domain) cost function |

---

## Running from WSL

Use the MCP MATLAB server (preferred — keeps session alive, no startup penalty).  
Or via CLI:

```bash
"/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe" \
  -batch "cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.')); RunModel"
```
