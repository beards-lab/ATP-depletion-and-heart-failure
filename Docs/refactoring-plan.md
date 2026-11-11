# Codebase Reorganization and Refactoring Plan

## Context

The codebase is scientifically mature but has accumulated technical debt: files are in wrong directories, documentation is sparse, magic numbers are unexplained, and two core scripts (`RunBakersExp.m` at 1043 lines, `evaluateProblem.m` at 989 lines) are too large to maintain easily. The goal is to improve readability and maintainability without changing behavior, with execution speed as a secondary objective.

---

## Phase 1: Dead Code Removal (safest, do first)

**Files to delete:**
- `Thrash/` — 26 files, zero references from active code
- `Auxiliary/fitslackOnset.m` — 19-line stub, superseded by `fitSlackForceOnset.m`, not a proper function

**Verification:** Run `grep -r "fitslackOnset\|Thrash" Model/ Workbench/ Auxiliary/ --include="*.m"` to confirm no references before deleting.

---

## Phase 2: File Relocations

Move files to directories matching their purpose. Use `git mv` to preserve history.

| File | From | To | Reason |
|------|------|----|--------|
| `Driver.m` | `Model/` | stays in `Model/` | Minimal entry point — canonical "how to run the model" example |
| `DriverSimple.m`, `DriverSimple_DAN.m`, `driverSA.m` | `Model/` | `Workbench/` | Variant/experiment drivers belong in the playground |
| `ModelParamsInitNiceSlack_prescribedSrxD.m` | `Model/` | `params/` | Parameter snapshot, not model code |
| `ModelParamsInitNiceSlack_prescribedSR_var2.m` | `Model/` | `params/` | Parameter snapshot, not model code |
| `ModelOptParamsFeaturesOvernight.m` | root | `params/` | Generated parameter snapshot |

Update `CLAUDE.md` and `README.md` after relocation to reflect new paths.

**Verification:** Run `RunBakersExp` after each move to confirm nothing breaks.

---

## Phase 3: Add Function Documentation (Help Text)

Add MATLAB-style help text to every function that currently lacks it. Format:

```matlab
function out = myFunc(a, b)
% MYFUNC  One-line summary.
%
%   OUT = MYFUNC(A, B) detailed description.
%
%   Inputs:
%     A - description (units if applicable)
%     B - description
%
%   Outputs:
%     OUT - description
%
%   See also: RELATEDFUNC
```

**Priority order (highest-impact first):**

### 3a. `Model/getParams.m`
- Add top-level doc: purpose, all 4 arguments, return value, the 4-step pipeline
- Add section headers: `%% 1. Fill defaults`, `%% 2. Apply modifiers`, `%% 3. Resolve expressions`, `%% 4. Compute strain grid`, `%% 5. Initialize state vector`
- Document the `'='` prefix convention for linked params (e.g., `params.kamh = '=0.1*kah'`)

### 3b. `Model/dPUdT_CombinedTransitions.m`
- Add top-level doc: purpose, inputs `(t, PU, params)`, outputs `(f, outputs, rates)`
- Document PU state vector layout at the top:
  ```
  % PU = [p1(1..ss), p2(1..ss), [p3(1..ss),] P_SR, NP, SL, LSE, PD, [P_SRD,] [x_dash]]
  % where ss = number of strain bins (params.ss)
  ```
- Add section headers for major blocks: `%% Unpack state`, `%% Overlap & lattice correction`, `%% Transition rates`, `%% Super-relaxed dynamics`, `%% Force & ODE derivatives`

### 3c. `Model/evaluateModel.m`
- Add top-level doc: purpose, inputs, outputs, description of multi-segment integration strategy
- Document the `out` struct fields

### 3d. `Auxiliary/` functions missing docs
- `myoutput.m` — optimizer iteration callback (explain what it prints/logs)
- `plotRates.m` — plot transition rates from `out` struct
- `updateRates.m` — scale all turnover rates by `xrate`
- `handleAndRethrowCostException.m` — explain the error accumulation pattern and message format
- `runStatesInTime.m` — explain why the wrapper exists (MATLAB workspace scope issue)
- `evalFeatureCost.m` — document the NaN cost convention and weighting

### 3e. `Model/resolveParams.m`
- Already has reasonable docs; add a usage example in the header comment.

---

## Phase 4: Name Magic Numbers as Constants

Replace unexplained literal values with named constants at the top of each function/script. Do **not** change the values.

### 4a. `Model/dPUdT_CombinedTransitions.m`
```matlab
MAX_RATE = 1e4;    % clamp transition rates to prevent numerical blow-up during stiff phases
MIN_PROB = 0;      % probability floor — clamp negatives caused by integrator overshoot
```

### 4b. `Auxiliary/fitSlackForceOnset.m`
```matlab
% Baker lab slack protocol: 5 release velocities recorded at fixed sample indices
VELOCITY_INDICES = [3, 7, 11, 15, 19];   % column indices of each velocity in slack data table
ZONE_BOUNDS = [1162, 1209; ...];          % [start, end] sample windows per velocity zone
```

### 4c. `Auxiliary/calcSensitivities.m`
```matlab
DELTA_FRACTION = 0.1;   % one-at-a-time perturbation size (10% of nominal)
```

### 4d. `Auxiliary/extractSlackAttributes.m`
```matlab
PEAK_MIN_PROMINENCE = 0.5;  % minimum peak prominence for detection (kPa)
FORCE_THRESHOLD_KPA  = 10;  % threshold to declare force onset
SL_SEARCH_DELTA      = 0.01; % sarcomere length increment for peak search (um)
```

### 4e. `Auxiliary/evalFeatureCost.m`
- Replace `NaNCost = 10` with a named constant and comment explaining the penalty magnitude choice.
- Replace `str2num` with `str2double` (`str2num` uses `eval` internally — MATLAB best practice).

---

## Phase 5: Clean Up Dead Code (Conservative)

**Rule:** Only remove blocks that are clearly replaced by active code in the same file, or are explicitly guarded by `if false` with no future value. Leave anything whose intent is uncertain.

**Priority targets:**

- `Auxiliary/calcSensitivities.m`: The `if false` block (lines 16–20) is the **right-side** (+10%) perturbation of a two-sided finite difference — intentionally disabled to save compute. **Keep it**, but add a comment:
  ```matlab
  % NOTE: right-side perturbation disabled to halve compute cost.
  % Enable if two-sided (central) differences are needed.
  if false
  ```
  Remove only the legacy `evaluateProblem` call on lines 12–13 (superseded by `evaluateBakersExp`).

- `Auxiliary/fitRecovery.m` lines 117–164: hardcoded Baker plot data stored as comments, no active path references it — remove.

- `Auxiliary/extractSlackAttributes.m` lines 101–211: 110-line commented-out original implementation, replaced by active code above — remove.

- `Auxiliary/updateRates.m` lines 23–24: commented-out alternative rate scaling, superseded by active line above — remove.

---

## Phase 6: Refactor RunBakersExp.m

`RunBakersExp.m` is a 1043-line script implementing 4+ experimental protocols. **Keep it as a script** (so it continues to set workspace variables without needing return values), but extract each protocol's implementation into a standalone function that the script calls.

**New structure:**
```
Model/
  RunBakersExp.m               -- script coordinator (~150 lines): calls protocol functions, collects E into workspace
  makeMonotonous.m             -- extract existing helper from bottom of RunBakersExp
  mergeOutStructs.m            -- extract existing helper from bottom of RunBakersExp
  experiments/
    runFVExperiment.m          -- function: returns [E_fv, out_fv] given params0 + data
    runKtrExperiment.m         -- function: returns [E_ktr, out_ktr]
    runSlackExperiment.m       -- function: returns [E_slack, out_slack]
    runStairsExperiment.m      -- function: returns [E_stairs, out_stairs]
```

**Approach:** Each protocol function takes `(params0, data, evalFlags)` and returns its local error contribution and output struct. `RunBakersExp.m` remains a script that:
1. Calls each protocol function
2. Sums contributions into `E`
3. Sets the same workspace variables (`E`, `Force_pred`, etc.) that `evaluateProblem.m` reads

This avoids any API changes to `evaluateProblem.m` while making each protocol independently testable.

---

## Phase 7: Input Validation in getParams.m

Add lightweight warnings at function entry (all `warning()`, never `error()`, to preserve exploratory flexibility):

```matlab
% g may be empty (no modifiers applied); validate only when non-empty
if ~isempty(g) && numel(g) ~= numel(params.mods)
    warning('getParams: g has %d elements but params.mods has %d entries — must match when non-empty', ...
        numel(g), numel(params.mods));
end

% Warn if SL0 is outside physiological sarcomere length range
if isfield(params, 'SL0') && (params.SL0 < 1.4 || params.SL0 > 3.0)
    warning('getParams: SL0=%.2f um is outside physiological range [1.4, 3.0]', params.SL0);
end
```

---

## Phase 8: Performance (Deferred — TODO)

Skip for now. Add a TODO comment at the top of `dPUdT_CombinedTransitions.m`:

```matlab
% TODO: profile with `profile on; RunBakersExp; profile viewer` to identify bottlenecks.
%       Candidates: overlap/lattice recomputation on every ODE RHS call.
```

Current ODE tolerances (`AbsTol=1e-4`, `RelTol=1e-2`) are loose and computations are vectorized. No changes warranted without profiling data.

---

## Verification Plan

After each phase, run in MATLAB:

```matlab
% Smoke test
params = getParams([], []);
[force, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 1], params);
disp(force(end))  % should be same as baseline

% Full regression (save E_baseline before starting Phase 6)
RunBakersExp
assert(abs(E - E_baseline) < 1e-6, 'E changed after refactor')
```

---

## Files to Modify — Summary

| File | Phase | Change |
|------|-------|--------|
| `Thrash/` (26 files) | 1 | Delete |
| `Auxiliary/fitslackOnset.m` | 1 | Delete |
| `Model/DriverSimple*.m`, `Model/driverSA.m` | 2 | Move to `Workbench/` |
| `Model/ModelParamsInit*.m` (2 files) | 2 | Move to `params/` |
| `ModelOptParamsFeaturesOvernight.m` | 2 | Move to `params/` |
| `CLAUDE.md`, `README.md` | 2 | Update paths |
| `Model/getParams.m` | 3a, 7 | Doc + validation warnings |
| `Model/dPUdT_CombinedTransitions.m` | 3b, 4a, 8 | Doc + named constants + TODO comment |
| `Model/evaluateModel.m` | 3c | Doc |
| `Auxiliary/myoutput.m` | 3d | Doc |
| `Auxiliary/plotRates.m` | 3d | Doc |
| `Auxiliary/updateRates.m` | 3d, 5 | Doc + remove dead code |
| `Auxiliary/handleAndRethrowCostException.m` | 3d | Doc |
| `Auxiliary/evalFeatureCost.m` | 3d, 4e | Doc + `str2double` |
| `Auxiliary/fitSlackForceOnset.m` | 4b | Named constants |
| `Auxiliary/calcSensitivities.m` | 4c, 5 | Named constant + remove legacy lines 12–13, annotate `if false` block |
| `Auxiliary/extractSlackAttributes.m` | 4d, 5 | Named constants + remove dead code |
| `Auxiliary/fitRecovery.m` | 5 | Remove dead code |
| `Model/RunBakersExp.m` | 6 | Split into script coordinator + per-protocol functions |

---

## Out of Scope

- `evaluateProblem.m` refactor (989 lines): High risk during active optimization campaigns; defer
- `StatesInTime.m` GUI refactor: Fragile MATLAB UI callbacks; low ROI
- Changing ODE solver or tolerances: Risk of behavior change
- Adding a test framework: Separate decision required
- Any changes inside `Thrash/` beyond deletion
