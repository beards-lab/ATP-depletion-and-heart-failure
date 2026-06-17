# Reorg Smoke Test — task list for Claude Code

Verifies the 2026-06-16 folder reorganization didn't break the model. Run from a WSL shell with
the repo as cwd MATLAB pattern (from `CLAUDE.md`)(as porposed below) or just run the matlab MCP if available:

```bash
"/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe" -batch "cd('C:\home\git\ATP-depletion-and-heart-failure'); addpath(genpath('.')); <CODE>"
```

> Why this matters: the reorg only *moved* files. MATLAB resolves functions by name on the path,
> so a clean run is expected — these checks confirm it and catch path-shadowing / bare-name
> `load()` surprises.

In this pass just report the failures and if the proposed fix does not work, try ONE workaround / fix. If that simple one did not help, report and stop.

## Phase 1 — Path & shadowing

- [ ] **Path builds with no error/warning.** Run `addpath(genpath('.')); disp('PATH OK')`. Expect `PATH OK`, no warnings.
- [ ] **Core functions resolve to the right place.** Run
  `which getParams; which dPUdT_CombinedTransitions; which evaluateModel; which RunBakersExp`.
  Expect all under `Model/`.
- [ ] **Relocated probes resolve once** (now in `Auxiliary/diagnostics/`). Run
  `which -all srxProbe; which -all ktrProbe; which -all mechProbe; which -all isoProbe`.
  Expect exactly one hit each, under `Auxiliary/diagnostics/`.
- [ ] **Relocated param scripts resolve once.** Run
  `which -all SimplestFVOptim; which -all SimplestFVOptim2; which -all SimplestFVOptim3`.
  Expect one hit each, under `params/`.
- [ ] **Known pre-existing duplicate** (not caused by reorg). Run `which -all ModelOptParamsFeaturesOvernight`.
  Expect **2** hits (`params/` and `ModelOptParams/`). Note which wins; decide later if it matters.
- [ ] **Archived dead code is still resolvable but unused.** Run `which dPUdTCaSimple`. Expect a hit
  under `_archive/dead_ode/` (harmless; confirms nothing breaks if a stray reference exists).

## Phase 2 — Core model run

- [ ] **Build params + integrate once.** Run
  `params = getParams([], []); [force, out] = evaluateModel(@dPUdT_CombinedTransitions,[0 1],params); fprintf('force=%g, numel(out.t)=%d\n', force(end), numel(out.t));`
  Expect finite `force`, non-empty `out`, no error.
- [ ] **Default ODE variant + a legacy one both still run.** Repeat the line above with `@dPUdTCa`
  and `@dPUdTCaSimpleAlternative2State`. Expect no "undefined function" errors.

## Phase 3 — Full evaluation (the real baseline check)

- [ ] **`RunBakersExp` completes and `E` is sane.** Set up `params0` as usual, then
  `RunBakersExp; disp(E); fprintf('Etot=%.6g\n', sum(E));` Record `Etot` here: `__________`.
- [ ] **`Etot` matches the pre-reorg value.** It should be identical (only file locations changed).
  Optional rigorous compare in Phase 5.

## Phase 4 — Bare-name `load()` resolution (trouble #2)

Scripts that `load('X.mat')` by bare name rely on path search; these `.mat` now live in `_archive/`
or an `Analyses/*/results/` folder — both on `genpath('.')`. Confirm they're still found:

- [ ] `exist('IdentifyLastSlack_state.mat','file')` → **2** (moved to `Analyses/LastSlackIdentification/results/`).
- [ ] `exist('bakers_slack8mM_all.mat','file')` and `exist('bakers_ktr_8.mat','file')` → **2** (in `data/`, untouched).
- [ ] `exist('Ghost_NoReattach_slack.mat','file')` → **2** (used by `Workbench/DriverSimple.m`; now in `_archive/`).
- [ ] `exist('filenames.mat','file')` → **2** (used by `Workbench/RunAllParametrizatinos.m`; now in `_archive/`).

> If any returns 0, the file isn't on the path — repoint that script's `load` to the new location,
> or keep `_archive/` on the MATLAB path.

## Phase 5 — Optional: rigorous baseline diff

- [ ] Create a pristine pre-reorg worktree at HEAD and compute its `Etot`, then compare:
  `git worktree add ../atp-baseline HEAD` → run Phase 3 there → diff the two `Etot` values.
  They should match to numerical tolerance. Remove with `git worktree remove ../atp-baseline`.

## Phase 6 — Finish up

- [ ] **Remove git-repair leftovers** (created during the move; safe to delete on Windows):
  `rm -f .git/index.lock.* .git/index.corrupt.*`.
- [ ] **Remove the locked leftover dir** if Explorer freed it: `Analysis/batch_output_extracted_v1`
  (gitignored), and the duplicated `resources/` at root.
