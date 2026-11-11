# Skill: update-mex — Regenerate MEX after switch/parametrization change

## When to use
Run this skill whenever the active feature-flag set changes (e.g. enabling
`UseA1AttachmentKernel`, changing `NumberOfStates`, adding a new pchip spline,
or modifying `PieceWiseStrainDepX` knot count). The MEX encodes one fixed code
path; any switch change requires regenerating all three layers below.

---

## File inventory

| File | Role |
|------|------|
| `Model/packMexVals.m` | Packs `params` → flat `double` array `mex_vals` (187 elements). Config-specific. |
| `Model/dPUdT_CombinedTransitions_mex_flat.m` | Pure-MATLAB reference with identical fixed code path. Used to validate index layout before compiling C. |
| `Model/dPUdT_core_mex.c` | C MEX. Must match `packMexVals` index layout exactly. |
| `Model/dPUdT_CombinedTransitions_mex.m` | Thin wrapper: `dPUdT_core_mex(t, PU, params.mex_vals, params.Vums, params.LXBpivot)`. Rarely changes. |
| `Model/compileMex.m` | Compiles `dPUdT_core_mex.c` via OpenModelica GCC batch file. |
| `Model/getParams.m` | Step 6c calls `packMexVals` when `UseCompiledMex && UseFastPPval`. |

Active flag configuration is in `ModelOptParams/ModelParamsSlackKtrOpt.m`.

---

## Current mex_vals layout (187 elements, 0-based in C / 1-based in MATLAB)

```
[0-42]   43 scalars  — kd k1 k_1 k2 kstiff1 kstiff2 kstiff1_n kstiff2_n
                        ka kah kamh dS ss dr dr2
                        L_thick L_hbare L_thin k_pas Lsc0 gamma
                        kSE ekSE MaxSlackNegativeForce
                        mu mu_neg mu2
                        ksr0 sigma1 sigma2 kmsr kmsrd sigma_srd1
                        ksrd sigma_srd2 ksrd2sr ksr2srd
                        slope d_actin s_threshold_R s_threshold_L
                        A1AttachmentWidth estiff
[43-102] 60 entries  — s_grid (params.s padded to N_SS_MAX=60; actual count in mv[12]=ss)
[103-106] 4 values   — pchip R1D breaks   (N_R1D=3 segments)
[107-118] 12 values  — pchip R1D coefs    (3×4 col-major)
[119-123] 5 values   — pchip pwsd breaks  (N_PWSD=4 segments, R12)
[124-139] 16 values  — pchip pwsd coefs   (4×4 col-major)
[140-143] 4 values   — pchip R21 breaks   (N_R21=3 segments)
[144-155] 12 values  — pchip R21 coefs    (3×4 col-major)
[156-162] 7 values   — pchip pwsd2 breaks (N_PWSD2=6 segments, R2)
[163-186] 24 values  — pchip pwsd2 coefs  (6×4 col-major)
```

Pchip coef matrices are col-major (MATLAB default from `unmkpp`).
`s_grid` is always padded to 60; the C code loops only `ss` times (read from `mv[12]`).

---

## Step-by-step update procedure

### 1. Identify what changed
Read the new active params file (e.g. `ModelParamsSlackKtrOpt.m`) and note:
- Which `UseXxx` flags changed
- Whether `PieceWiseStrainDepX / R1DX / R21X / Dep2X` knot counts changed
  (determines `N_R1D`, `N_PWSD`, `N_R21`, `N_PWSD2`)
- Whether `NumberOfStates` changed (affects PU indexing)
- Whether `MaxStrainArraySize` changed (N_SS_MAX)

### 2. Update `packMexVals.m`
- Add/remove scalar fields in the `scalars` block (keep count = 43, or update `N_SCALARS` everywhere)
- Update `N_SS_MAX` if `MaxStrainArraySize` changed
- Update segment count comments for pchip blocks if knot counts changed
- The sanity check at the bottom asserts `numel(mv) == expected_len` — update `expected_len`

### 3. Update `dPUdT_CombinedTransitions_mex_flat.m`
- Mirror every change from step 2 (same index constants at top of file)
- Hard-code any newly active/inactive switches in the computation body
- Test it matches `dPUdT_CombinedTransitions`:

```matlab
ModelParamsSlackKtrOpt;
params0.UseFastPPval = true; params0.UseCompiledMex = true;
params0 = getParams(params0, [], false, false);
params0.Vums = 0;
PU0 = params0.PU0(:); t = 0;
f_ref  = dPUdT_CombinedTransitions(t, PU0, params0);
f_flat = dPUdT_CombinedTransitions_mex_flat(t, PU0, params0);
fprintf('max|err| = %.2e\n', max(abs(f_ref - f_flat)));
% Must be < 1e-12
```

### 4. Update `dPUdT_core_mex.c`
- Update `#define` constants at top to match new layout:
  - `N_SCALARS`, `N_SS`, `N_R1D`, `N_PWSD`, `N_R21`, `N_PWSD2`
  - All `OFF_*` offsets (derived from above)
  - `TOTAL_MV_LEN`
- Add/remove scalar unpacking lines to match `packMexVals` scalars block
- Update PU indexing if `NumberOfStates` changed (`p3` array, extra state scalars)
- Add/remove computation blocks for newly active switches

### 5. Recompile

```matlab
cd('C:\home\git\ATP-depletion-and-heart-failure');
addpath(genpath('.'));
clear dPUdT_core_mex;   % release lock before overwriting
compileMex;
```

### 6. Verify correctness

```matlab
% At normal working point (ss=60)
ModelParamsSlackKtrOpt;
params0.UseFastPPval = true; params0.UseCompiledMex = true;
params0 = getParams(params0, [], false, false);
params0.Vums = 0;
PU0 = params0.PU0(:); t = 0;
f_ref = dPUdT_CombinedTransitions(t, PU0, params0);
f_mex = dPUdT_CombinedTransitions_mex(t, PU0, params0);
fprintf('ss=%d  max|err|=%.2e\n', params0.ss, max(abs(f_ref - f_mex)));
% Must be < 1e-12

% At FV narrow condition (ss < 60)
params0.SL0 = 2.0;
dsl = 80/params0.kSE;
params0.Slim_r = 2.0 + params0.A1AttachmentWidth + params0.dS;
params0.Slim_l = 2.0 - dsl - params0.A1AttachmentWidth - params0.dS;
p2 = getParams(params0, [], true, false); p2.Vums = 0;
f2r = dPUdT_CombinedTransitions(t, p2.PU0(:), p2);
f2m = dPUdT_CombinedTransitions_mex(t, p2.PU0(:), p2);
fprintf('ss=%d  max|err|=%.2e\n', p2.ss, max(abs(f2r - f2m)));
% Must be 0.00e+00 or < 1e-12
```

### 7. Verify E matches

```matlab
% Sequential comparison (no parpool needed)
ModelParamsSlackKtrOpt;
params0.recalculateDataFeats = false; params0.ghostSave = '';
params0.RunSlackSegments = 'All';
LoadData; RunBakersExp; E_base = sum(E);

params0.UseFastPPval = true; params0.UseCompiledMex = true;
params0.modelFcn = '@dPUdT_CombinedTransitions_mex';
RunBakersExp; E_mex = sum(E);
fprintf('dE = %.2e  (expect < 1%%)\n', (E_mex-E_base)/E_base);
```

---

## Parpool note
MEX cannot run in `parpool('threads')` workers (MATLAB restriction).
Use `parpool('local', N)` when `UseCompiledMex = true`.
`RunModel.m` already sets this.

---

## Compile environment
- GCC from OpenModelica: `C:\Program Files\OpenModelica1.26.1-64bit\tools\msys\ucrt64\bin`
- MATLAB libs: `<matlabroot>\extern\lib\win64\mingw64\` (`libmx`, `libmex`, `libmat`)
- `compileMex.m` writes a temporary batch file to handle spaces in paths

If GCC path changes, update `gcc_bin` in `compileMex.m`.
