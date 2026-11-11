# Project Reorganization Plan

**Status:** EXECUTED 2026-06-16 (moves done with `mv`, not yet committed — review `git status` then `git add -A && git commit`). Originally proposed 2026-06-12.
**Date:** 2026-06-12
**Scope:** Folder/file *reorganization* (where things live). This is complementary to the existing
`Docs/refactoring-plan.md`, which covers *in-file* code quality (help text, magic numbers, splitting big
scripts). The two do not overlap and can be executed independently.

---

## 1. Guiding vision (from you)

- **`Analyses/`** — one self-contained folder per analysis: the **script(s)** + its **results** (figures/.mat)
  + its **conclusions** (a written `.md`). A top-level `Analyses/README.md` rolls these up into general
  **summaries and recommendations**.
- **`Model/`** — all necessary modeling files, nothing else.
- **`Auxiliary/`** — supporting, non-critical functions.
- **`Workbench/`** — the playground: new drivers, hand-tuning, optimization entry points, and tests.
  *Not* a place where finished analyses with results live.

## 2. Key facts that shape the plan

- **The clutter is not in git.** `.gitignore` already excludes `*.mat`, `*.png`, `*.asv`, `*.pptx`,
  `*.docx`, `*.fig`, `*.exe`, `data/`, `~*`. Only **384 files are tracked** (the `.m` code, `.md` docs,
  `.c`, configs). The ~400 MB of loose `.mat`/`.png` at the root is on disk only → moving or deleting it
  carries **zero git-history risk**.
- **Decisions you confirmed:**
  1. Clutter → **archive, don't delete** (sweep into `_archive/`, keep on disk, out of the way).
  2. Cross-analysis summaries → **`Analyses/README.md`** index.
  3. Param folders → **keep `params/` and `ModelOptParams/` separate**, document the distinction.

---

## 3. Target top-level structure

```
ATP-depletion-and-heart-failure/
├── Model/                 # core model ONLY (ODE RHS, getParams, evaluateModel, RunBakersExp,
│   └── experiments/       #   cost functions, feature extraction, data loaders, MEX, experiments/)
├── Auxiliary/             # non-critical utilities (plotting, parsing, sensitivity/fit helpers,
│   └── diagnostics/       #   serialization) + reusable diagnostic probes (SRXprobe)
├── Workbench/             # drivers, hand-tuning, optimization entry points — the playground
│   └── Tests/             # Test*.m scripts
├── Analyses/              # NEW: one subfolder per analysis (script + results + conclusions.md)
│   ├── README.md          #   index + rolled-up summaries & recommendations
│   ├── AttachedPool_vs_ATPase/
│   ├── BindingSiteOccupancy/
│   ├── RestretchMechanisms/
│   ├── SensitivityAnalysis/
│   ├── LastSlackIdentification/
│   └── PassiveForceID/
├── params/                # active parameter snapshots (current source of truth)
├── ModelOptParams/        # historical named param gallery + QC plots (archive of past optima)
├── DataCuration/          # raw experimental data -> model-ready pipeline
├── data/                  # raw experimental data (gitignored)
├── Docs/                  # cross-cutting docs, presentations, notes
│   ├── presentations/     #   .pptx / posters / abstracts
│   └── notes/             #   loose working notes
├── papers/                # reference literature PDFs
├── Modelica/              # Modelica source models
├── PassiveTitin/          # self-contained passive-titin sub-project
├── Figures/               # publication-quality figures
├── scripts/               # NEW: infra (setup-gh.sh, rol.sbatch)
└── _archive/              # NEW: disk-only graveyard for stale dumps & generated output (gitignored)
```

New top-level folders: **`Analyses/`**, **`scripts/`**, **`_archive/`**. Everything else already exists.

---

## 4. The `Analyses/` design

Each analysis folder follows one convention:

```
Analyses/<Topic>/
├── <DriverScript>.m        # the analysis script(s)
├── conclusions.md          # what it found + recommendation (rename existing report .md to this,
│                           #   or create a stub to be filled in)
└── results/                # .mat / .png outputs produced by the script
```

`Analyses/README.md` is the index: a one-line description + link per analysis, then a "Summaries &
Recommendations" section synthesizing across them (the layer you asked for).

### Proposed analysis bundles

| `Analyses/` folder | Scripts (from) | Conclusions `.md` | Results to co-locate |
|---|---|---|---|
| **BindingSiteOccupancy** | already a bundle in `Analysis/BindingSiteOccupancy/` | `OccupancySaturation_Report.md` ✅ (rename → `conclusions.md`) | `fv_probe_results.mat`, `occupancy_AB_results.mat` ✅ already there |
| **AttachedPool_vs_ATPase** | uses `SRXprobe/` probes (see §6) | `attached-pool-atpase-ktr-paradox.md` ✅ | root `Ghost_*.mat` fit snapshots referenced by it |
| **RestretchMechanisms** | `AnalyseRestretchMechanisms.m`, `TuneRestretch.m`, `RunMechanismEvaluation.m`, `OptimizeMechanismEvaluation.m`, `RunPhase2.m` | merge `Docs/restretch-analysis.md` + `restretch-discrepancy.md` + `mechanism-evaluation.md` → `conclusions.md` | `AnalyseRestretchMechanisms_results.mat` |
| **SensitivityAnalysis** | `SensitivityAllSlack.m` (current), `RunSensitivityAnalysis.m`, `FeatureCorrelation.m` | `sensitivity_analysis_explanation.md` (root) → `conclusions.md` | `SensitivityAllSlack_results.mat`, `SensitivityAnalysisResults.mat`, `ResidualsAndJacobian.mat`, `sensitivites.xlsx` |
| **LastSlackIdentification** | `OptimLastSlackPieceWise.m` (workhorse), `IdentifyLastSlack.m` | create `conclusions.md` stub | `OptimLastSlackPieceWise_state*.mat`, `IdentifyLastSlack_state.mat` |
| **PassiveForceID** | `IdentifyPassive.m` | `Docs/passive-force-subtraction.md` → `conclusions.md` | `batch_output_passive/` (or archive) |

> Note: `Analysis/` (singular, current) is **renamed/absorbed into `Analyses/` (plural)**. The remaining
> loose scripts in `Analysis/` that are *not* analyses are reclassified in §5.

---

## 5. Reclassifying the current `Analysis/` loose files

| File | Verdict | Destination |
|---|---|---|
| `OptimLastSlackPieceWise.m`, `IdentifyLastSlack.m` | analysis | `Analyses/LastSlackIdentification/` |
| `SensitivityAllSlack.m`, `RunSensitivityAnalysis.m`, `FeatureCorrelation.m` | analysis | `Analyses/SensitivityAnalysis/` |
| `AnalyseRestretchMechanisms.m`, `TuneRestretch.m`, `RunMechanismEvaluation.m`, `OptimizeMechanismEvaluation.m`, `RunPhase2.m` | analysis | `Analyses/RestretchMechanisms/` |
| `IdentifyPassive.m` | analysis | `Analyses/PassiveForceID/` |
| `BatchRunAllParams.m` | utility/infra (batch-renders param sets) | `Workbench/` |
| `SumATPSlackFitPlots.m` | one figure for an abstract | `Workbench/` (or `Figures/` generators) |
| `ExperimentWAttachments.m`, `TestSinglePerturbation.m`, `TestNewOutputFeatures.m` | scratch/tests | `Workbench/Tests/` |
| `ResidualOpt.m` | helper function (firstLowercase-style wrapper) | `Auxiliary/` |
| `RunOptimPiecewise.m` | superseded (missing state file) | `_archive/` |
| `analyzeSystemRHS.m` | dead end ("no difference") | `_archive/` |
| `driverSA.m`, `runSA.m` | obsolete SA pipeline | `_archive/` |
| `batch_output/`, `batch_output_extracted_v1/`, `batch_output_passive/`, `fit_plots/` | disposable generated PNGs (reproducible) | `_archive/` |
| `*.asv` (5 files) | MATLAB autosaves | delete (gitignored junk) |

---

## 6. `Model/`, `Auxiliary/`, `Workbench/` — relocations

These three dirs are mostly correct; the issue is **stray data files mixed in with code** and a few
genuinely misplaced scripts.

### Model/ — keep core, evict the rest

| File | Action |
|---|---|
| Core: `getParams`, `resolveParams`, `evaluateModel`, `RunBakersExp`, `evaluateProblem`, `Driver`, cost fns (`evalFeatureCost`, `evalPhysiologyCost`, `parameterBounds`, `boundedOutputFn`…), feature extractors, data loaders, `experiments/` | **stay** |
| `dPUdT_CombinedTransitions.m` (+ `_mex`, `_mex_flat`), `dPUdTCa.m`, `dPUdTCaSimpleAlternative2State.m` | **stay** (active or still-referenced) |
| `dPUdTCaSimple.m`, `dPUdTCaSimpleAlternative.m`, `dPUdT_TransitionRates.m` | → `_archive/` (dead ODE variants, no callers) |
| `ModelParams_tuesdayLunch.m` | → `params/` (param snapshot, not model code) |
| `Ghost_oujaDouble*.mat`, `envOpt*.mat` | → `_archive/` (stale workspace dumps) |
| `dPUdT_core_mex.mexw64` | leave on disk; add `*.mexw64` to `.gitignore` |
| `*.asv` (5 files) | delete |
| `experiments/runFVTimecourseExperiment.m` | → `Auxiliary/` (plot-only, no cost return) — minor, optional |

### Auxiliary/ — keep utilities, evict data

| File | Action |
|---|---|
| All 28 utility `.m` files (plotting, parse/split helpers, `calcSensitivities`, `fitRecovery`, `writeParamsToMFile`, etc.) | **stay** |
| `SimplestFVOptim`, `SimplestFVOptim2` (no extension) | → `params/` (add `.m`; confirmed param scripts) |
| `SA_all.mat`, `SA_pwsd.mat`, `env_SA.mat`, `env_fmin*.mat`, `last_env.mat`, `optim_hist_feats.mat`, `Ghost_SimplestOptim_FV.mat` | → `_archive/` |
| `preprocessDrivingSignal.asv` | delete (and note: `preprocessDrivingSignal.m` itself lives correctly in `DataCuration/`) |

### Workbench/ — keep drivers/tests, evict data

| File | Action |
|---|---|
| All driver/optim/tuning/test `.m` (`RunModel`, `RunOptim*`, `Handtune*`, `DriverSimple*`, `Compare*`, SRX validation scripts, `Tests/`) | **stay** |
| `ModelParamsInitNiceSlack_prescribedSR`, `SimplestFVOptim3` (no ext) | → `params/` (add `.m`) |
| `Ghost_*.mat`, `env*.mat`, `sens*.mat`, `_params.mat`, `testenv.mat` | → `_archive/` |
| `XBBakersDataFit_*.png` | → `Figures/` (or `_archive/`) |
| `sensitivites.xlsx` | → `Analyses/SensitivityAnalysis/results/` |
| `*.asv` (incl. `LoadAndPlotLogs.asv`) | delete |

**Diagnostic probes:** `Analysis/SRXprobe/` (`srxProbe`, `ktrProbe`, `mechProbe`, `isoProbe`) is a reusable
library used by both the AttachedPool analysis *and* Workbench scripts. Proposed home: **`Auxiliary/diagnostics/`**
(reusable, non-critical functions). Alternative: keep inside `Analyses/AttachedPool_vs_ATPase/` if you
consider them single-purpose. Flagging as a judgment call.

---

## 7. Root-level cleanup

| Group | Examples | Action |
|---|---|---|
| **Scratch/junk** | `tmp.mat`, `tmp.txt`, `done.mat`, `Unnamed.mat`, `modeltesting*.mat`, `*_smoke.png`, `*.asv`, empty `.mat` | → `_archive/` (per your "archive, don't delete") |
| **Optimizer dumps** | `env*.mat` (~248 MB, incl. 118 MB `env_fmin9.mat`), `tmpOpt*.mat` (~52 MB), `forcevelocity_out.mat` (63 MB), `out.mat`, `x2/x3.mat` | → `_archive/` |
| **Result `.mat` tied to an analysis** | `*_results.mat`, `*_state*.mat`, `ResidualsAndJacobian.mat`, `Ghost_*.mat` referenced by a report | → that analysis's `results/` where identifiable, else `_archive/` |
| **Result PNGs** | `XBBakersDataFit_*.png`, `param_parallel.png`, `param_ranges_strip.png`, `Slack protocol.jpg` | → `Figures/` |
| **Param snapshots** | `params.csv`, `modifierstbl.csv`, `gaOutparams.csv`, `all params.txt`, `PARAMETERS.xlsx`, `*_params.mat`, `ModelParamsInitNiceSlack_dr01`, `SimplestFVOptim` | → `params/` |
| **Loose notes** | `problem_summary.md`, `sensitivity_analysis_explanation.md`, `literature_mechanism_evaluation_request.md` | → `Docs/notes/` (or the matching analysis) |
| **Dymola build artifacts** | `dsmodel.c`, `dsfinal/dsin/dslog.txt`, `buildlog.txt`, `dymosim.exe`, `Xbmodel.prj` | → `_archive/` (or `Modelica/`); ensure gitignored |
| **Infra** | `setup-gh.sh`, `rol.sbatch` | → `scripts/` |

---

## 8. Secondary-directory tidy (light touch)

- **`Docs/`** — sound; add `presentations/` (the `.pptx`/poster/abstract files) and `notes/` subfolders;
  delete all `~$*` Office lock files; consider archiving superseded poster versions (`v1`–`v4`, keep `v5`).
- **`params/` vs `ModelOptParams/`** — keep both. Add a one-paragraph header to each (and to `README.md`)
  stating the roles: `params/` = active/current snapshots; `ModelOptParams/` = historical named optima + QC
  plots. No file moves between them.
- **`resources/`** — empty → delete.
- **`tools/`** — single 18 MB `.exe`; add to `.gitignore` (leave on disk).
- **`PassiveTitin/`** — keep as sub-project; internal pruning (its `*.asv`, `tmp.mat`, `SoHot.mat`, old
  `DataStruct*` snapshots) is out of scope for this pass — flag for a later dedicated cleanup.
- **`Modelica/`** — keep as-is.

---

## 9. `.gitignore` additions

Add: `*.mexw64`, `_archive/`, and (if not implicitly covered) `tools/*.exe`. Everything else is already
ignored.

---

## 10. Proposed execution order (each step independently verifiable)

1. **Create empty scaffold** — `Analyses/`, `Analyses/README.md`, `_archive/`, `scripts/`,
   `Docs/presentations/`, `Docs/notes/`, `Auxiliary/diagnostics/`. (No risk.)
2. **Delete `.asv` autosaves** everywhere (gitignored junk, ~36 files). (No risk.)
3. **Sweep clutter → `_archive/`** (root dumps, generated image folders, stray `.mat` in code dirs).
   Disk-only, not in git.
4. **Relocate tracked `.m` files** with `git mv` to preserve history: the Analyses bundles (§4–5), param
   snapshots → `params/`, probes → `Auxiliary/diagnostics/`, dead ODE variants → `_archive/`.
5. **Move/rename conclusion `.md`s** into their `Analyses/<Topic>/conclusions.md`; write `Analyses/README.md`.
6. **Tidy `Docs/`** (lock files, presentations/, notes/) and update `.gitignore`.
7. **Update `CLAUDE.md` + `README.md`** to reflect new paths (notably the `Analyses/` directory and the
   `Analysis/` → `Analyses/` rename).
8. **Verify**: run the smoke test and confirm nothing broke —
   ```matlab
   params = getParams([], []);
   [force, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 1], params);
   RunBakersExp   % E should match pre-move baseline
   ```
   MATLAB resolves functions by name on the path, so moving files within `addpath(genpath('.'))` scope is
   safe — but the smoke test confirms it.

---

## 11. Out of scope (flagged for later)

- In-file refactors in `Docs/refactoring-plan.md` (help text, magic numbers, splitting `RunBakersExp.m`).
- Internal pruning of `PassiveTitin/`.
- Merging `params/` and `ModelOptParams/` (you chose to keep them separate).
- Deleting anything — this plan **archives**; a later pass can purge `_archive/` once you're confident.

---

## 12. Open judgment calls for your call

1. **SRXprobe home** — `Auxiliary/diagnostics/` (reusable) vs inside `Analyses/AttachedPool_vs_ATPase/`.
2. **`SumATPSlackFitPlots.m` / `BatchRunAllParams.m`** — `Workbench/` vs a `Figures/`-generation analysis.
3. Whether `_archive/` should eventually be a **separate git-LFS or out-of-repo store** given it holds
   the ~400 MB of dumps.

---

# Part B — Documentation Reorganization

**The pain:** the project's *knowledge* is scattered — durable reference docs sit next to one-off prompts,
analysis conclusions are buried in `Analysis/` subfolders, TODOs live in four different files, and `Docs/`
mixes 24 markdown files with 38 presentations and 11 stale lock files. There's no single entry point, so
the project feels overwhelming and finished findings are hard to build on.

**The fix, in one line:** separate **durable knowledge you build on** from **transient working notes**,
pull every analysis conclusion into its `Analyses/` bundle, consolidate all scattered TODOs into one
roadmap, and add two `README.md` "start here" maps so nothing is hidden.

## 13. Target `Docs/` structure

```
Docs/
├── README.md            # NEW — START HERE: one-page map of all knowledge (links to everything below
│                        #   + to Analyses/README.md). The single entry point that kills the overwhelm.
├── ROADMAP.md           # NEW — consolidated TODOs / next-steps / open questions (see §15)
├── reference/           # DURABLE knowledge — the docs you build on, rarely go stale
│   ├── parameter-reference.md
│   ├── mechanism-evaluation.md          (108 KB literature analysis — the big one)
│   ├── fitting-strategy.md
│   ├── RundownCorrection.md
│   ├── passive-force-subtraction.md
│   ├── sources.md
│   └── update-mex.md                    (how-to guide)
├── notes/               # TRANSIENT working notes & lab diaries (chronological, may go stale)
│   ├── labdiary.md
│   └── labdiary_boundedfit.md
├── experiments/         # dated experimental result write-ups (pair with data/)
│   └── results-0327.md
├── presentations/       # all .pptx / .docx / posters / abstracts
│   └── archive/         #   superseded versions (poster v1–v4, old meeting decks)
└── (figures handled under top-level Figures/ — see §16)
```

## 14. Per-document disposition

Every current doc, with where it goes and why:

| Doc (current path) | Class | Destination |
|---|---|---|
| `Docs/parameter-reference.md` | durable reference | `Docs/reference/` |
| `Docs/mechanism-evaluation.md` | durable reference | `Docs/reference/` |
| `Docs/fitting-strategy.md` | durable reference | `Docs/reference/` |
| `Docs/RundownCorrection.md` | durable reference | `Docs/reference/` |
| `Docs/passive-force-subtraction.md` | durable reference | `Docs/reference/` |
| `Docs/sources.md` | durable reference (lit links) | `Docs/reference/` |
| `Docs/update-mex.md` | durable how-to | `Docs/reference/` |
| `Docs/labdiary.md` | working note | `Docs/notes/` (TODOs extracted → ROADMAP) |
| `Docs/labdiary_boundedfit.md` | working note | `Docs/notes/` |
| `Docs/results-0327.md` | dated experiment result | `Docs/experiments/` |
| `Docs/attached-pool-...` is in `Analysis/` | analysis conclusion | `Analyses/AttachedPool_vs_ATPase/conclusions.md` (Part A §4) |
| `OccupancySaturation_Report.md` (in `Analysis/`) | analysis conclusion | `Analyses/BindingSiteOccupancy/conclusions.md` |
| `Docs/restretch-analysis.md` + `restretch-discrepancy.md` | analysis conclusion | merge → `Analyses/RestretchMechanisms/conclusions.md` |
| `sensitivity_analysis_explanation.md` (root) | analysis conclusion | `Analyses/SensitivityAnalysis/conclusions.md` |
| `Docs/low-ATP-force-enhancement-analysis.md` | analysis findings (script TBD) | `Analyses/LowATP_ForceEnhancement/conclusions.md` (new bundle, seed) |
| `problem_summary.md` (root) | temporary framing prompt | `_archive/` |
| `literature_mechanism_evaluation_request.md` (root) | superseded by mechanism-evaluation.md | `_archive/` |
| `Docs/2026 Myofilament-meeting-abstract.md` | publication draft | `Docs/presentations/` |
| `Docs/*.pptx / *.docx / *.pdf` (38 files) | presentations | `Docs/presentations/` (+ `archive/` for old versions) |
| `Docs/~$*` (11 lock files) | junk | delete (gitignored) |
| `README.md`, `CLAUDE.md` | project entry / agent guide | **stay at root** |
| `Docs/refactoring-plan.md`, `Docs/reorganization-plan.md` | meta/process | stay in `Docs/` (or `Docs/process/`) |

This means **every analysis conclusion ends up in exactly one place** — its `Analyses/<Topic>/conclusions.md`
— and `Docs/reference/` holds only the cross-cutting, durable knowledge. No more hunting.

## 15. Consolidate scattered TODOs → `Docs/ROADMAP.md`

TODOs and "next steps" are currently spread across at least four docs. Pull them into one prioritized
backlog (and leave the source docs as historical record, with their TODO blocks replaced by a pointer to
ROADMAP):

| Source | Items to lift |
|---|---|
| `labdiary.md` → "Running TODOs" | 5 items (vernier effect test, 30% attached-pool target, single-slack timecourse fit, right-side rate-curve sensitivity, vall2_dy second-peak rise) |
| `restretch-analysis.md` → "Next Steps" | 5 items (fill sweep table, fminsearch 10-param run, verify FV unaffected, write g to params file, update diary) |
| `restretch-discrepancy.md` → "Recommended Next Step" | increase kSE 3–5×, then tune slope / s_threshold_R |
| `2026 Myofilament-meeting-abstract.md` | "immediate next step": mechanistic ID of ATP-dependent differences |

`ROADMAP.md` groups these by theme (e.g. *Attachment kinetics*, *Restretch shape*, *ATP-dependence*,
*Housekeeping*) so you see the real shape of the work instead of fragments. This is the most direct
antidote to the overwhelm.

## 16. Figures (light touch)

Figures are already semi-organized; the goal is only to make them findable, not to re-sort all 30+.

- **`Figures/`** — keep as the figure store. Optionally split into `Figures/raw/` (`.fig` source) and
  `Figures/exported/` (`.png`/`.svg`); keep existing `TransitionPanels/` and `proposal/` subfolders as-is.
- **`Docs/figures/`** (the `fig403_*`, `fig420_*`, `fig46_*` report figures) — these pair with the
  experiment/result write-ups, so move them alongside → `Docs/experiments/figures/` and keep the relative
  links in `results-0327.md` working.
- Loose result PNGs at the repo root (`XBBakersDataFit_*.png`, `param_*.png`, `Slack protocol.jpg`) →
  `Figures/` (already noted in Part A §7).

## 17. Two `README.md` maps (the anti-overwhelm layer)

1. **`Docs/README.md`** — one page: "What do you want?" → links to reference docs, the roadmap, the
   analyses index, presentations, and the experiment log. First thing anyone (you or an agent) opens.
2. **`Analyses/README.md`** — the analyses index from Part A §4: one line + status + link per analysis,
   then the rolled-up *Summaries & Recommendations* synthesizing across them.

Together these mean the project has exactly **two front doors**, and every other document is one click away
and in a predictable place.

## 18. Docs execution steps (fold into Part A §10)

After the Part A scaffold step, additionally:

- a. Create `Docs/reference/`, `Docs/notes/`, `Docs/experiments/`, `Docs/presentations/archive/`.
- b. Delete the 11 `Docs/~$*` lock files.
- c. `git mv` the durable docs into `reference/`, working notes into `notes/`, `results-0327.md` into
  `experiments/`, presentations into `presentations/` (old versions → `archive/`).
- d. Move analysis conclusions into their `Analyses/<Topic>/conclusions.md` (renaming as needed); archive
  the two superseded prompt docs.
- e. Write `Docs/ROADMAP.md` (consolidated TODOs) and replace the source TODO blocks with a pointer.
- f. Write `Docs/README.md` and `Analyses/README.md`.
- g. Update `README.md` + `CLAUDE.md` to point at the new doc locations.

## 19. Updated open judgment calls (Part B)

4. **Analysis conclusions** — co-locate in `Analyses/<Topic>/conclusions.md` (recommended, chosen) vs keep
   a mirror copy in `Docs/reference/`. Mirroring risks drift; a link from `Docs/README.md` is cleaner.
5. **`low-ATP-force-enhancement-analysis.md`** — promote to its own `Analyses/LowATP_ForceEnhancement/`
   bundle (it reads like a real analysis) vs file under `Docs/reference/` if you consider it settled.
6. **Meta/process docs** (`refactoring-plan.md`, this file) — leave in `Docs/` vs a `Docs/process/` subfolder.
```
