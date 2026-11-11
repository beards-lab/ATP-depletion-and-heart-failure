# Project-management triage draft (2026-10-30)

> **SUPERSEDED 2026-11-02: filed on GitHub as issues #2–#70, milestones 1–8, board https://github.com/orgs/beards-lab/projects/3. Edit the issues, not this file.** Each issue footer names its draft id.

~~Nothing here is on GitHub yet.~~ Edit freely: delete issues, change labels/milestones, move
items between sections, write `KEEP` / `DROP` / `MERGE INTO <id>` next to anything. I file only
what survives, after your go-ahead. Draft ids are prefixed by milestone (`M0-3`, `A-5`, `M1-2`, `Mα-1`, … `U-` = unassigned).
They are triage handles only: on GitHub each issue gets the next `#n` and its milestone field;
"blocked by A-5" becomes a real `#n` link at filing time.

Repo: `beards-lab/ATP-depletion-and-heart-failure`.

---

## 0. Inventory: what exists today

| Where | State |
|---|---|
| GitHub issues | **0** (none ever filed) |
| GitHub PRs | 1, merged 2024 (`Creating-ATP-exp-fit-plots`) |
| Labels | GitHub defaults only (bug, documentation, duplicate, enhancement, good first issue, help wanted, invalid, question, wontfix) |
| Milestones | none |
| Projects | none for this repo. Org `beards-lab` has *Cardiovascular simulators* and an untitled one, both unrelated |
| **`main`** | **488 commits behind** the working branch `PooledTarget-8mM` and has 1 commit of its own (2025-11, "missed branch"). **All real work since 2023 lives off `main`.** |
| Working branch | `PooledTarget-8mM` (local is 1 commit ahead of origin, unpushed). Untracked: `Analyses/HumanTwitch`, `HumanizeRoadmap`, `MomentsReduction`, `PooledFitFigures`, two 2026-10 decks |
| Live unmerged branches | `lowatp-2state` (+48, low-ATP hypotheses H0–H4, worktree `ATP-wt-lowatp`); `CaActivation` (+2, Ca thin-filament layer, worktree `ATP-wt-ca`); `HumanTwitch` = `EHTRelativeTransfer` = `claude/humanize-xb-model-231445` (one commit under three names, +3); `claude/tnc-feedback-coupling-637355` (+5, contains the previous one plus the TnC closed loop) |
| Dead/superseded branches | `BackToBasics`, `FittingCaDecay`, `PlotAttached`, `PlotTitinOutputs`, `ProximalStiffeningOnly`, `rOLLbACK` (all abandoned); `LowATP`, `thrash-branch` (superseded by `lowatp-2state`); `Refactoring-backup-pre-author-rewrite` (content already in base); `fulldata` (1 small commit) |
| Already merged (safe to delete) | `Comparing-Dan's-version`, `Creating-ATP-exp-fit-plots`, `Further-optim`, `GenFIgus`, `Generating-figures`, `Refactoring`, `SRD-state`, `Testing`, `claude/pedantic-chaum-b44eea`, `floored-campaign` (its worktree `ATP-wt-floored` too) |
| Plan/status docs | `Docs/ROADMAP.md` (**stale**, 2026-06-16); `Docs/experiments/status-2026-09-16.md` §8 (the most recent task list); `Analyses/README.md` (settled-numbers table + index); ~45 `Analyses/*/conclusions.md` with "next steps"; `Analyses/HumanizeRoadmap/conclusions.md` (2026-10-27, untracked); `Analyses/PiAccumulation/plan.md`; `Docs/refactoring-plan.md` (mostly done); `Docs/manuscript/` (TODOs and REVIEW-corrections); `Docs/notes/labdiary*.md` (stale, duplicated in ROADMAP) |

**Where open work lives today, in short:** in ~45 `conclusions.md` "next steps" sections, the
09-16 status doc §8, and in my memory files. The June ROADMAP is mostly superseded. Below every
thread is either a draft issue (§3) or marked done/stale (§4).

---

## 1. Milestones

Yours, lightly reworded so each has a checkable "done". Proposals are marked **NEW**.

| id | Title | Done when | Notes |
|---|---|---|---|
| **M0** **NEW** | **Branch & release model** | (1) The branch model in §1b is agreed and written down as a decision record; (2) `main` is up to date and is the integration trunk; (3) the paper line `paper1` exists and holds only the rodent 8 + 2 mM work; (4) a regression test proves "all non-paper flags off on `main` = `paper1` numbers"; (5) dead branches are deleted (archive tags kept); (6) every open thread is an issue | Prerequisite for "one issue → one branch → one PR". PRs wait until this is settled. ~1–2 weeks |
| **A** | **New-data fit (8 + 2 mM)** | The fit to the new dataset scores better on a **pre-registered** metric than the current fit does on the old set (see A-2) | Parallel; data not yet arrived. Contains the re-stretch mechanism work, because it limits both conditions |
| **1** | **Mouse Ca transient (wrap-up)** | `CaActivation` merged, `UseCa=0` bit-identical, hand-off doc written | Mostly done per you; missing wrap-up/hand-off |
| **α** | **Ca transient ATP sensitivity** | Pre-registered test either shows an ATP effect on the Ca/thin-filament side (e.g. pCa50 shift 8→2 mM) or rejects it | Blocked on force–pCa data at 8 and 2 mM |
| **2** | **Humanized model** | Human-parameterised model reproduces a held-out set of human/literature measurements within stated tolerances (Land 2017 ramps, Lewalle k_tr 15±3 s⁻¹, …) | From `HumanizeRoadmap` |
| **3** | **Reduced cell-mechanics model** | A reduced (moments) model reproduces the full model on the battery within stated tolerance and runs inside a host (TriSeg or similar) | From `MomentsReduction` + `TriSegCoupling` |
| **4** | **Pi and pH sensitivity** | Pi (and pH) dependence calibrated against literature and attributed | Only Pi has a plan (`PiAccumulation/plan.md`); pH is unscoped |
| **P1** **NEW** | **Paper 1: rodent 8 + 2 mM fit** | Manuscript submitted, and the published model is released: tag `paper1-v1.0`, maintenance branch `release/paper1`, archived with a DOI, figures reproducible from the tag | `Docs/manuscript/` exists with 5 reviews done. Real external goal. **Reports the new dataset** (decided 2026-11-02), so A blocks P1 |

_Decided 2026-11-02:_ the re-stretch transient mechanism stays inside A, not its own milestone (yet).

### 1b. Branch & release model (proposal, to become decision record M0-12)

The model serves four goals. Each has different stability needs:

| Goal | Milestones | Who consumes it | Stability need |
|---|---|---|---|
| **Published rodent model**: high- and low-ATP fits to mouse data | A, P1 | reviewers, readers, anyone reproducing the paper | **Frozen.** Bit-reproducible from a tag; bugfixes only |
| **Ca-activated, humanized myocyte** | 1, α, 2 | you and collaborators; the next paper(s) | Under development; flags change often |
| **Reduced cell-mechanics component for organ models** (TriSeg etc.) | 3 | other groups' heart models | **A stable interface** (inputs/outputs/units), even while internals change. *This is the goal you didn't list.* It's a third consumer with its own versioning needs |
| **Metabolite extensions** (Pi, pH) | 4 | feeds both the humanized model and future papers | Development |

Proposed branches:

```
main ───●───────●──────●──────●────────●──────▶   integration trunk: everything, new features behind flags (default OFF)
         \             ▲      ▲
paper1    ●──●──●──●───┘──●───┘──(tag paper1-v1.0)──▶ release/paper1   (bugfixes only, forever)
          rodent 8+2 mM work only; merges flow paper1 → main, never main → paper1
```

- **`paper1`** (long-lived until submission): today's `PooledTarget-8mM` plus `lowatp-2state`, **and milestone A's new-data fit** (Paper 1 reports the new dataset, decided 2026-11-02). No Ca, no humanization, no reduced model. Every change merges **forward** into `main`, never backward.
- At submission: tag **`paper1-v1.0`**. The branch becomes **`release/paper1`**, bugfixes only. Each fix is a new tag (`paper1-v1.0.1`), merged forward into `main`.
- **`main`**: the integration trunk. Ca, humanize, reduced and Pi work land here behind flags. "Eventually merge using settings to discriminate" happens here.
- **The guard that makes this safe:** a regression test on `main`. With every non-paper flag at its default (off), the paper battery reproduces the `paper1` tag's numbers. If this test is broken, `main` has silently changed the published model.
- Feature branches (one per issue) fork from `paper1` for paper work and from `main` for everything else.
- **Later**, when the reduced model has an external user, it gets its own version tags (`xb-reduced-v0.x`) on `main`.

Things that are easy to miss in a release, now issues M0-13 and P1-6:
- **Data are not in git.** A tag alone can't reproduce the paper. The release needs a data manifest (file list + checksums, see `scripts/fulldata/`, branch `fulldata`) and a pinned pooled target.
- **One entry script** that regenerates every paper figure/table from the tag (`Analyses/PooledFitFigures` is close).
- **Frozen parameter snapshot** for the published fit, under a name the optimiser can never overwrite (the same rule as `Example/params/ExampleHighATP.m`).
- **A DOI**: a GitHub release archived on Zenodo (also a manuscript TODO).

---

## 2. Labels (small set)

Replace the defaults with about ten labels. Type labels are one per issue.

| label | color | meaning |
|---|---|---|
| `bug` | red (keep) | code does the wrong thing |
| `science` | purple | open scientific question or hypothesis test; the answer is a finding, not code |
| `analysis` | blue | run something and write a `conclusions.md` (well-posed; the method is known) |
| `objective` | orange | changes to what or how the fit scores (features, weights, floors); these change every cost number |
| `refactor` | grey | code structure, deletions, no behaviour change (parity test required) |
| `docs` | light blue (rename `documentation`) | write-ups, hand-offs, README, manuscript |
| `data-request` | yellow | needs bench work or data/info from the Baker lab or another external source |
| `chore` | grey | repo hygiene: branches, worktrees, untracked files |
| `needs-judgment` | dark red | needs *your* review/decision. Counts toward the WIP limit of 2 |
| `agent-ready` | green | fully specified; an agent can do it end-to-end and the acceptance check is mechanical |
| `blocked` | black | waiting on another issue or external input (say which in the body) |

Delete: `enhancement`, `good first issue`, `help wanted`, `invalid`, `question`, `wontfix` (use
"close as not planned" instead). Keep `duplicate`.

## 3. Board

An **org-level** Project v2 under `beards-lab` (so collaborators see it), named
*ATP depletion & heart failure*. `Status` field: Backlog / Next / In progress / Review / Done.
A **WIP limit cannot be set through the API**: it's a board-column setting in the web UI. I'll give
you the two-click instruction, or you can treat it as a convention. All filed issues start in
Backlog. My suggestion for **Next** is marked ★ below.

## 4. Repo additions (step 3 PR)

- `ROADMAP.md` at the repo root (~15 lines: goal, current milestone, Now / Next / Later, linking issues).
  Replace `Docs/ROADMAP.md` with a one-line pointer to it, keeping the old content in git history.
- `.github/ISSUE_TEMPLATE/task.md` (Summary / Motivation / Acceptance criteria / Out of scope / Pointers)
  and `config.yml` (blank issues off).
- `decisions/` with a `README.md` and a template. **Optional:** seed it with ~5 records of decisions
  already made: pooled target over 03/27, deadband floor over pin, ktr protocol held out unscored,
  features not traces, `estiff` = 1.25. Want these?
- _Decided 2026-11-02:_ **no PRs until the branch model (M0) is resolved.** Labels, milestones,
  issues and the board live on GitHub only and need no branch, so they can be filed now.
  `ROADMAP.md`, the issue template and `decisions/` are committed once M0-12 says where they go
  (most likely onto `paper1`, which merges forward into `main`).

---

## 5. Draft issues

Format: **<id> · title** · `labels` · milestone · relations. ★ = suggested for "Next".

### M0 — One trunk

**M0-1 ★ · Bring `main` up to date with the working line** · `chore` `needs-judgment` · M0
- **Summary:** `main` is 488 commits behind `PooledTarget-8mM` and has 1 unique commit (2025-11 "missed branch").
- **Motivation:** One-issue-one-PR needs a base branch that is the truth. Today every branch forks off `PooledTarget-8mM`.
- **Acceptance:** `git rev-list --count main..PooledTarget-8mM` = 0 after the merge. The 1 main-only commit is either merged or explicitly dropped. Local `PooledTarget-8mM` is pushed first (it's 1 ahead).
- **Out of scope:** merging the live feature branches (M0-3, M0-4, M0-5).
- **Pointers:** `origin/main` @ 2023-02 ↔ local `main` @ 2025-11.

**M0-2 · Commit the untracked analyses and decks** · `chore` `agent-ready` · M0
- **Summary:** Untracked: `Analyses/HumanTwitch`, `HumanizeRoadmap`, `MomentsReduction`, `PooledFitFigures`, `Docs/presentations/2026-10-*`, `.claude/settings.json`, plus edits to `.gitignore` and `Analyses/README.md`.
- **Acceptance:** `git status` is clean, and large binaries (pptx/png) are checked against `.gitignore` policy before committing.
- **Out of scope:** content changes.
- **Pointers:** branch `PooledTarget-8mM`.

**M0-3 ★ · Merge `lowatp-2state` (H0–H4 ATP hypotheses) into the paper line** · `chore` `needs-judgment` · M0 · blocked by M0-1
- **Summary:** 48 commits / 142 files. `AtpEffectHypotheses`, low-ATP staircase and poster figures. `AtpEffectHypotheses/` exists on both sides, so expect conflicts. Per memory, the 2 mM gate code lives **only** there and the main tree silently runs 8 mM.
- **Acceptance:** merged. `atphyp_H4_opt` reproduces its recorded 2 mM cost on trunk. Worktree `ATP-wt-lowatp` removed.
- **Out of scope:** reconciling H4 with the later `pAtp_k2` sweep (that is A-9).
- **Pointers:** `C:/home/git/ATP-wt-lowatp`, `params/atphyp_H4_opt.m`.

**M0-4 · Consolidate the humanize branches** · `chore` · M0
- **Summary:** `HumanTwitch` = `EHTRelativeTransfer` = `claude/humanize-xb-model-231445` (same commit). `claude/tnc-feedback-coupling-637355` contains them plus the TnC loop. Keep one branch (`HumanTwitch`, fast-forwarded to tnc), delete the other three names, and prune the matching worktrees.
- **Acceptance:** one humanize branch, merged to trunk or with an open PR. `git worktree list` has no prunable entries.
- **Pointers:** worktrees `ATP-wt-human`, `ATP-wt-eht`, `.claude/worktrees/*`.

**M0-5 · Delete dead branches and worktrees** · `chore` `needs-judgment` · M0
- **Summary:** Delete local+remote: the abandoned set (`BackToBasics`, `FittingCaDecay`, `PlotAttached`, `PlotTitinOutputs`, `ProximalStiffeningOnly`, `rOLLbACK`), the superseded set (`LowATP`, `thrash-branch`, `Refactoring-backup-pre-author-rewrite`), and the already-merged set (§0). Decide on `fulldata`. Tag before deleting (`archive/<name>`) so nothing is lost.
- **Acceptance:** `git branch -a` shows only trunk + live branches, and the archive tags are pushed.
- **Out of scope:** `_archive/` (M0-10).

**M0-6 · `genpath('.')` picks up model copies in `.claude/worktrees`** · `bug` `agent-ready` · M0
- **Summary:** The documented boilerplate `addpath(genpath('.'))` shadows `Model/` with worktree copies, so runs silently use the wrong code.
- **Acceptance:** after the boilerplate, `which dPUdT_CombinedTransitions` resolves to `Model/` even with a worktree present. CLAUDE.md is updated.
- **Pointers:** `Analyses/MomentsReduction/conclusions.md`, `Auxiliary/refreshPool.m`.

**M0-7 · Delete dead switches and retire knobs** · `refactor` · M0
- **Summary:** MechanismAblation shows none of these beat the 0.055 noise floor: `UseSerialStiffness` (read nowhere, yet listed in CLAUDE.md as key), `UseA2Reattaching` (+ `kA2re`/`kA2hop`/`sA2hop`), `UseMaxwellTensionOnly`, `UseP31Shift`, `UseKstiff3`, `UseTORNegShift`, `UseAtpOnUNR`, `UseMutualPairingAttachment`, `UseTitinModel`. Also the optimiser pool clean-up (−12 params) and the retired-knob registry (§9 proposal).
- **Acceptance:** `refit9_opt` (or the current best) cost is bit-identical before/after; MEX regenerated with parity; CLAUDE.md switch list updated.
- **Out of scope:** re-running ablation arms (A-12).
- **Pointers:** `Analyses/MechanismAblation/REPORT.md`, `Model/dPUdT_CombinedTransitions.m`, `dPUdT_core_mex.c`, `getParams.m`.

**M0-8 · `RunBakersExp` silently overwrites workspace `E`** · `bug` `agent-ready` · M0
- **Acceptance:** a caller's `E` survives a `RunBakersExp` call (test in `Workbench/Tests`).
- **Pointers:** `Analyses/OverstretchedHeads/conclusions.md`.

**M0-9 · Optimiser minimises L2 while reports show L1** · `bug` `objective` · M0
- **Summary:** `costExp` differs between `evaluateBakersExp` (2) and `RunBakersExp` reporting.
- **Acceptance:** one documented convention; report numbers equal optimiser numbers for one seed.
- **Pointers:** `Analyses/RestretchRecoveryFit/conclusions.md`.

**M0-10 · Decide the fate of `_archive/` (~400 MB) and `PassiveTitin/` clutter** · `chore` `needs-judgment` · M0
- From the June ROADMAP housekeeping. **Acceptance:** a decision record, then executed.

**M0-11 · Refactoring-plan leftovers** · `refactor` `agent-ready` · M0
- **Summary:** `calcSensitivities.m` still calls legacy `evaluateProblem`. Phases 3–4 (help text, named constants) are unverified. `ktrProbe.m` hand-builds `PU0` and supports neither `A_reg` nor `D0`.
- **Acceptance:** no callers of `evaluateProblem` outside `RunOptimLakes`; `ktrProbe` runs on `fastgate_best`.
- **Pointers:** `Docs/refactoring-plan.md`, `Auxiliary/diagnostics/ktrProbe.m`.

**M0-12 ★ · Decision record: branch & release model** · `docs` `needs-judgment` · M0 · blocks M0-1, M0-3, M1-1, M2-6
- **Summary:** Write §1b as `decisions/0001-branch-and-release-model.md` once you've edited it: which branch is the paper line, forward-only merges, flags default off on `main`, release tags, where repo docs live.
- **Motivation:** Every merge issue (M0-1, M0-3, M1-1, M2-6) depends on the answer, and collaborators and agents need one place that says which branch to fork from.
- **Acceptance:** the record is merged, and CLAUDE.md names the branch to fork from for paper work and for everything else.
- **Out of scope:** executing the merges.

**M0-13 · Regression guard: non-paper flags off on `main` reproduce `paper1`** · `refactor` `agent-ready` · M0 · blocked by M0-12
- **Summary:** A test in `Workbench/Tests` runs the paper battery on the frozen paper snapshot with all non-paper flags at default, and compares the cost and per-feature values against numbers stored from `paper1`.
- **Motivation:** This is what makes "merge everything, discriminate by settings" safe. Without it, Ca/humanize work can silently change the published model.
- **Acceptance:** the test passes on `paper1` and on `main`, and fails if a non-paper flag's default is flipped (checked once by hand). Tolerance: bit-identical where the code path is unchanged, else ≤ 0.1 % on cost.
- **Out of scope:** CI (run by hand or by the agent before each merge to `main`).

### A — New-data fit (8 + 2 mM)

*Data acquisition*

**A-1 · Data-curation: `mergeLogsAndBursts` skips filenames containing `repeat`** · `bug` `agent-ready` · A
- One-line fix; needed before the new data arrives. **Acceptance:** a `*repeat*` file in a test folder is merged.
- **Pointers:** `Analyses/ProtocolDesign`, `Analyses/RundownCorrection`, `DataCuration/`.

**A-2 ★ · Pre-register what "better fit on the new data" means** · `science` `needs-judgment` `docs` · A
- **Summary:** The milestone's done-criterion must be fixed **before** the data arrive, or any number can be argued. Candidates: pooled-target χ² per feature (the s²/n floor), held-out ktr-protocol k80 and staircase scores, 2/8 ratios (force ×1.36, ktr ×0.55 → ×0.72 by protocol).
- **Acceptance:** a decision record in `decisions/` with the metric, the baseline value of the current best on the old set, and the comparison rule.
- **Out of scope:** running the fit.
- **Pointers:** `Docs/experiments/status-2026-09-16.md`, `Analyses/DataPooling`, `Analyses/FeatureWeighting`.

**A-3 · Data round 2: the bracketed ATP ladder (6 sessions)** · `data-request` · A
- **Summary:** Request `8-4-2-8` / `2-8-2` sessions with a bracket in every one, sex-balanced; stop at S4 without a 4 mM claim, or at S5 with one. Record L₀/SL.
- **Acceptance:** the data are in `data/`, curated, and pooled with `WritePooledTarget`.
- **Pointers:** `Analyses/ProtocolDesign/conclusions.md` §8, `Analyses/RundownCorrection/strategy.md`, memory `protocol_design_round2`.

**A-4 · Add-on measurements to request with round 2** · `data-request` · A
- **Summary:** (a) a second end-of-session repeat on another prep, to calibrate the rundown φ=0.45 (called "the single most valuable next measurement"); (b) a hi-res `*_stiff_ktr` + `*_slack` on a bracket (prediction: protocol ktr −3 % vs slack ktr −12 %); (c) an FV point near 0.1–0.25 ML/s (model hump +21 %); (d) a third [ATP], ideally 20 mM (Hill vs power gates are unidentifiable otherwise); (e) measured SL if the rig allows.
- **Acceptance:** a list sent to the Baker lab, with what came back recorded here.
- **Pointers:** ATPEffectReconciliation, ProtocolDesign §9, KtrProtocolVsSlackKtr §7, status 09-16 §8 #6.

*The structural defect (could be milestone B)*

**A-5 ★ · Re-stretch transient: yield gate + lengthening recruitment refit (`yieldD`)** · `science` `needs-judgment` · A
- **Summary:** Weak (A1) heads are held too long on lengthening. The yield gate solves `vall_t` and the staircase nadir. The residual names lengthening recruitment.
- **Acceptance:** cost ≤ 0.529 (the control), §5 physiology kept, staircase step amplitude no worse than refit9's, and evaluated at 2 mM too.
- **Out of scope:** the slip campaign (A-6).
- **Pointers:** `Analyses/OverstretchedHeads`, seed `yieldB_opt`, `RegAvailLengtheningGain`≈0.01.

**A-6 · Re-stretch: three-arm actin-site slip campaign** · `analysis` · A · related A-5
- **Summary:** Arm a scores 0.501 before any refit. The campaign is pre-flighted and needs 16+ cores.
- **Acceptance:** `EvaluateYieldCampaigns.m` report; the winning arm is adopted or rejected in a decision record.
- **Pointers:** `RunSlipCampaign.m`, `params/slip{A,B,C}_seed.m`, `slip.sbatch`.

**A-7 · Commit the yield-gate model code with an off-switch regression test** · `refactor` `agent-ready` · A · blocks A-5, A-6
- **Summary:** `UseYieldGate`, `s_yield`, `YieldHop1`, `YieldHopShift`, `RegAvailLengtheningGain` and the `*q` vertex features. Check whether they are committed (`cb0adc0` suggests partly).
- **Acceptance:** `refit9_opt` with everything off = 0.2079.

*Objective debts (each changes the cost numbers: do them before the new-data fit, together)*

**A-8 · Re-weight the 2 mM objective from `s2` throughout** · `objective` · A
- **Summary:** `BuildWeightedObjective` still derives every 2 mM weight from the 8 mM spread. Also: no joint 2/8 ratio objective exists yet.
- **Acceptance:** the weights table shows `s2`-derived weights, and the decision on a ratio objective is recorded.
- **Pointers:** status 09-16 §8 #1, `Analyses/FeatureWeighting`.

**A-9 · Fit the low-ATP gate on the ATP-carrying features only** · `science` · A · related A-5
- **Summary:** `pAtp_k2=0.75` alone reproduces force 2/8 (1.39 vs 1.33) and ktr (0.538 vs 0.544), but is rejected by the 8 mM re-stretch block. The 8 mM defect outvotes the ATP physiology. Get an honest `k2` from `rsK80`, `ktr`, `steady`, `rsT0` only. Also in LowATP_Refit9K2Sweep: restretch timing at 2 mM (`rsTfitE` ×1.0 vs 1.68), and `steady` weight ~17× below its precision.
- **Acceptance:** a `k2` gate value with its CI, plus a decision on whether nucleotide-dependent SRX (`kmsrd`) is needed.
- **Pointers:** status 09-16 §9, `Analyses/LowATP_Refit9K2Sweep`, memory `atp_carrying_metrics`.

**A-10 · Tighten `XTOR` to 8–12 and score net (not gross) turnover** · `objective` `science` · A
- **Summary:** The bound was relaxed by hand to 15; the model sits at 14.1 with attached fraction 0.255 vs a physiological 0.12–0.15.
- **Acceptance:** a refit with `XTOR[1]|8-12` reported, with whether the 8 mM basin moves.
- **Pointers:** status 09-16 §9 option A, `Analyses/AttachedPool_vs_ATPase`, `Analyses/SRXDuringContraction`.

**A-11 · Put the staircase in the objective (features, not trace)** · `objective` · A
- **Summary:** It's a free out-of-sample test that is known wrong in three ways. Promote `extractStairsAttributes.m` to `Model/`; score the 1 % peak, nadir and recovery; let `estiff` float with `kstiff2` compulsory.
- **Acceptance:** staircase features are in `fn`, and the current best is re-scored.
- **Out of scope:** the physical exponent question (A-16).
- **Pointers:** `Analyses/StairsPrediction`.

**A-12 · Floor the two unfloored RMSE terms (`PS_rampupRMSE`, `PS_holdDecayRMSE`)** · `objective` · A
- 18–27 % of the cost, with no noise model. **Acceptance:** each has a noise model and a floor, documented.
- **Pointers:** `Analyses/RefitStrategy`, status 09-16 §8 #3.

**A-13 · Force block: floors built from between-heart spread are ~100× too wide; level undefended (+17 %)** · `objective` `science` · A
- **Summary:** Every force feature's floor uses the 23 % between-heart force spread, not the 2 % shape agreement. Separately, the absolute level sits +17 % over pooled `steady` with nothing defending it.
- **Acceptance:** a decision on the floor basis plus a level guard, and a re-scored best.
- **Pointers:** memories `force_block_deadband_defect`, `pooled8_campaign`; `Analyses/DataPooling`.

**A-14 · FV: narrow the floor and add a scored velocity near 0.25 ML/s** · `objective` · A · related A-4
- **Summary:** The FV floor drifted to 2 sd; the model FV is non-monotonic (+21 % at 0.10 ML/s) because `FV_velocities` steps from 0 to 0.5.
- **Acceptance:** the hump is ≤ +2 % at 0.1 ML/s and the FV term is inside the narrowed floor.
- **Pointers:** `Analyses/KtrProtocolVsSlackKtr` §7, memory `objective_deadband_design` (never narrow uniformly below 0.25).

**A-15 · Re-read every protocol with the robust reader and re-score all seeds** · `analysis` `agent-ready` · A
- **Summary:** FV/stairs/passive/ktr were never re-read with `FeatureReadGrid='data'`. Seeds were scored under old solver settings, and the `s_num` floors need rebuilding. Also decide `vall_t`/`peak1_t` drop, `tFall50`, `rsK80`→`rsKls`, and retire `ovrsht_*`/`ktr2_*` if still in `fn`.
- **Acceptance:** a table of every seed under the current settings, committed.
- **Pointers:** `Analyses/FeatureRobustness`, `Analyses/StrainWindowNoise`, `Analyses/SlackDataAnalysis`.

*Open science, lower priority*

**A-16 · Which `estiff` exponent is physical (≈0.90 vs 1.25)?** · `science` · A
- Needs a step-size protocol. **Acceptance:** a protocol proposal plus a model prediction that discriminates. **Pointers:** `Analyses/StairsPrediction`.

**A-17 · ktr ladder: slope has the wrong sign (model +0.21 vs data −1.85)** · `science` · A
- **Summary:** Slack segments are an SL ladder; the data fall 15 % in all 3 preps, the model rises. `P_bound_max=1.208` disables the candidate mechanism; `L_thick`=1.17 µm acts as a fitting device. The refit6 12 h log is empty.
- **Acceptance:** model slope within the data's spread, or a recorded rejection with the reason.
- **Pointers:** `Analyses/KtrSegmentSlope`, memory `ktr_ladder_length_dependence`.

**A-18 · Slack ktr feature is biased 9–11 % low on model traces (10 kPa threshold)** · `bug` · A
- **Acceptance:** bias < 2 % on synthetic traces of known rate. **Pointers:** `extractSlackAttributes.m`, KtrProtocolVsSlackKtr.

**A-19 · Apply the ×1.2 passive correction to the active runs** · `science` · A
- **Summary:** The passive protocol ran after 4 activations. Test it on 04/10, where the reversed order should flip the sign.
- **Acceptance:** an applied/rejected decision. **Pointers:** `Analyses/PassiveRundownATP`.

**A-20 · How much SRX is there during contraction? Design a protocol that loads the pool** · `science` `data-request` · A
- **Summary:** Size is not identifiable (×0.5–×2 refit equally). Candidate: a 5–80 ms release–hold series. Also SRX at low ATP and 2 vs 8 mM mant-ATP data.
- **Pointers:** `Analyses/SRXDuringContraction`, `Analyses/LowATP_Refit9K2Sweep`.

**A-21 · Sawtooth SRS predictions: request the bench test** · `data-request` · A
- **Summary:** Untested model predictions (net work out at 0.1–0.2 ML/s, post-release overshoot 1.18–1.22 F₀, 2nd-ramp stiffness ×1.99 at 5 ms). All are on `RegAvailShorteningOnly`. Also at 2 mM.
- **Pointers:** `Analyses/SawtoothSRS`.

### 1 — Mouse Ca transient (wrap-up)

**M1-1 ★ · Rebase and merge `CaActivation`** · `chore` `needs-judgment` · 1 · blocked by M0-1
- **Summary:** 2 commits, 18 files (+2.1k), 26 behind trunk. The manuscript records it as unmerged.
- **Acceptance:** merged; `UseCa=0` is **bit-identical** to the trunk on the full battery; at pCa 4.51 ≥ 98.5 % of force is kept.
- **Pointers:** worktree `ATP-wt-ca` (prunable), `Analyses/CalciumTwitch`.

**M1-2 · `sigLookup` Hermite path is dead code** · `bug` `agent-ready` · 1
- Lines 16–18 do nearest-sample then `return`. **Acceptance:** the Hermite path is used and a unit test interpolates a known signal.

**M1-3 ★ · Write the Ca-layer hand-off doc** · `docs` · 1 · blocked by M1-1
- **Summary:** One page: what the layer does, flags (`UseCa`, `CaInput`, `CaModel`), how to drive it (analytic / measured / T-World), what's validated, and known limits (twitch 5 % with literature constants; `k_on` is the steep unidentified knob).
- **Acceptance:** a newcomer runs a mouse twitch from the doc alone (agent-checkable).
- **Pointers:** `Analyses/CalciumTwitch/conclusions.md`, memory `calcium_twitch`.

**M1-4 · Pin rig temperature and the diastolic SRX fraction (0.139)** · `data-request` · 1 · related M2-1
- **Pointers:** `Analyses/CalciumTwitch`.

### α — Ca transient ATP sensitivity

**Mα-1 · Pre-register the α test and its rejection rule** · `science` `needs-judgment` · α
- **Summary:** The prediction is a leftward pCa50 shift at 2 mM, unrun. State the effect size that would count, and what result rejects ATP sensitivity of the Ca side.
- **Acceptance:** a decision record.

**Mα-2 · Model prediction: force–pCa at 8 vs 2 mM** · `analysis` · α · blocked by M0-3, M1-1
- Needs the 2 mM gate (in `lowatp-2state`) and the Ca layer on one trunk. **Acceptance:** pCa50 and Hill n at both [ATP], with uncertainty from the gate's CI.

**Mα-3 · Request force–pCa data at 8 and 2 mM** · `data-request` · α
- **Summary:** Existing sub-maximal data are PNB-inhibited, so unusable. **Acceptance:** the request is sent and its status tracked here.

### 2 — Humanized model

**M2-1 ★ · Ask the Baker lab: rig temperature, species, β-MHC/human prep feasibility** · `data-request` · 2
- The named blocker in HumanizeXB. **Acceptance:** answers recorded in `Analyses/HumanizeXB`.

**M2-2 · Explicit Q10 block per rate group** · `science` · 2 · blocked by M2-1
- Layer 0 of HumanizeRoadmap. **Acceptance:** temperature is a parameter, and Q10 values are cited per group.

**M2-3 · Forward prediction: Land 2017 ramps and Lewalle k_tr on the human snapshots** · `analysis` · 2 · blocked by M2-2
- **Acceptance:** ramp-end forces 0.2 / 0.4 / 0.5 F₀; k_tr 15 ± 3 s⁻¹. These are the cheapest first test.
- **Pointers:** `Analyses/HumanizeRoadmap`, `params/human_v0_seed.m`, `human_land_g4valid.m`.

**M2-4 · Human v0 seed through the 8 mM battery; human protocols with holds ≥ 1 s** · `analysis` · 2
- Add R2T (0.3–1.5 s⁻¹/head) and V0 bounds; drop the restretch features. **Pointers:** `Analyses/HumanizeXB`.

**M2-5 · Wire ADP-release inhibition (K_AD 10–20 µM)** · `science` · 2 · related milestone 4

**M2-6 · Merge the humanize/T-World/TnC work (`xbtnc`, `CaInput='tworld'`, `CaModel 'land'`, `c_post`)** · `chore` `needs-judgment` · 2 · blocked by M0-4, M1-1
- **Acceptance:** all off by default and bit-identical when off.

**M2-7 · Layer 2: intact twitch on T-World Ca (Chung 2019 frequency series; Chung 2018 SR block held out)** · `analysis` · 2 · blocked by M2-3, M2-6

**M2-8 · Literature data asks: Lewalle 2025 traces, DOI and lead verification** · `data-request` · 2

### 3 — Reduced cell-mechanics model

**M3-1 · Moments twin: split RHS into kernels + periphery + backends with parity tests** · `refactor` · 3
- **Acceptance:** full and moments backends share kernels; parity tests are in `Workbench/Tests`.
- **Pointers:** `Analyses/MomentsReduction`, `dMdT_Moments.m`, `momentsSetup.m`.

**M3-2 · Production PDE deletes force-bearing heads on window shift (ktr restretch: −20 % peak, −10 % k80)** · `bug` · 3 (also affects A)
- **Summary:** Also, the 2 nm `dS` biases k80 −4 %. Fix with `ClampStrainDepExtrapolation` plus a wider or redesigned `movingWindow`, then re-score the ktr validation.
- **Acceptance:** ktr k80 within 2 % of the continuum limit. **Pointers:** `PDEGridConvergence.m`.

**M3-3 · Moments twin fails on fast stretch (ktr k80 −8 %)** · `science` · 3 · related M3-2

**M3-4 · TriSeg Phase 0: regenerate MEX with parity, then an analytic Jacobian** · `refactor` · 3
- 119 s/beat today; ~16× from the Jacobian, ~20–50× from MEX. **Pointers:** `Analyses/TriSegCoupling`.

**M3-5 · TriSeg Phase 1: moment-truncation go/no-go (~1 week)** · `analysis` `needs-judgment` · 3 · blocked by M3-1, M1-1
- Also resolve the duplicated passive force.

### 4 — Pi and pH

**M4-1 · Pi Stage 0: transient diffusion in the prep (`RunPiTransient.m`)** · `analysis` · 4
- **Pointers:** `Analyses/PiAccumulation/plan.md`.

**M4-2 · Pi Stage 1: calibrate the Pi slope (re-enable the Pi term; match Ebus & Stienen)** · `science` · 4 · blocked by M4-1

**M4-3 · Pi Stage 2: attribution (Pi relief share <20 % / 20–50 % / >50 %)** · `analysis` · 4 · blocked by M4-2

**M4-4 · Pi data asks: prep dimensions, bath Pi, Pi scavenger, added Pi; identify the Pi-diffusion paper** · `data-request` · 4

**M4-5 · Scope pH sensitivity** · `science` · 4
- **Summary:** Nothing exists yet. A literature pass and which rates to gate. **Acceptance:** a scoping `conclusions.md` in the style of CalciumTwitch.

### P1 — Paper 1

**P1-1 · Manuscript TODOs in `01-introduction-and-methods.md`** · `docs` · P1 · related M2-1
- Bench details (solutions, temperature, pCa, n, species, sex), SL channel, PNB identity, the final joint 8+2 mM fit, the K_m form, the ×1.2 passive factor, and a repository DOI.

**P1-2 · Resolve REVIEW-corrections §3–§4** · `docs` `needs-judgment` · P1
- §3: non-independent features, deadband, `estiff` knife-edge, ablation power, passive circularity. §4: missing literature, uncited `Schmidt2024`, cut H4.

**P1-3 · Parameter uncertainty: profile likelihoods + per-feature unfloored χ² table** · `analysis` · P1

**P1-4 · Decide: 0.2 mM arm in or out; ktr probe (slack ×0.581 vs protocol ×0.72); credit for the Beard22 SRX conjecture** · `science` `needs-judgment` · P1
- **Pointers:** memory `manuscript_intro_methods`, `Analyses/ATPEffectReconciliation`.

**P1-5 · Commit `PooledFitFigures` and re-run after each refit** · `docs` `agent-ready` · P1 · related M0-2

**P1-6 · Release the published model: tag, maintenance branch, data manifest, DOI** · `chore` `needs-judgment` · P1 · blocked by M0-12, M0-13
- **Summary:** At submission: freeze the parameter snapshot under a protected name; write the data manifest (file list + SHA-256, pinned pooled target); add one entry script that regenerates every paper figure/table; tag `paper1-v1.0`; create `release/paper1`; make a GitHub release archived to Zenodo.
- **Motivation:** The published model must stay reproducible while `main` moves on. Data aren't in git, so a tag alone isn't enough.
- **Acceptance:** on a clean clone of the tag plus the manifest's data, the entry script reproduces the paper's numbers; the DOI is in the manuscript.
- **Out of scope:** later model versions.
- **Pointers:** `scripts/fulldata/`, branch `fulldata`, `Analyses/PooledFitFigures`, `Example/params/ExampleHighATP.m` (frozen-name precedent).

### Unassigned (no milestone yet; your call)

**U-1 · `Docs/notes/claude-setup-improvement-candidates.md`: 9 tooling items** · `chore` · —
- Either file them individually or close as a batch.

**U-2 · 04/03 `t0` is 1.5–2× the other preps; check against that day's FV; detrend recordings before extraction; per-day passive reference** · `analysis` · — (A?)
- **Pointers:** `Analyses/SlackDataAnalysis`.

**U-3 · `D0_CAP = 5` in `buildPooledTarget`: argue it, don't inherit it** · `objective` · — (A?)
- **Pointers:** `Analyses/RestretchPeakAudit`.

**U-4 · ATPEffectReconciliation: 04/10 still an outlier after correction (×1.54 vs ×1.24/×1.29)** · `science` · — (A?)

---

## 6. Done or stale: proposed NOT to file (or file and close immediately for the record)

Mark any you want filed anyway.

| Thread | Source | Verdict |
|---|---|---|
| Phase-2 fminsearch, sweep table, raise `kSE` 3–5× | June ROADMAP / RestretchMechanisms | **Stale**: superseded by OverstretchedHeads (A-5) |
| Fix `vall2_dy` with catch k=0.05 + kSE=3000 | June ROADMAP / labdiary | **Stale**: same |
| Fit the time course of a single slack from `ModelParamsSlackKtrOpt` | June ROADMAP | **Stale**: pooled target replaces it |
| Velocity-dependent vernier on attachment | June ROADMAP | **Stale?** never run; `A_reg`/`RegAvail*` covers velocity-dependent availability. Keep? |
| Reach ~30 % attached pool | June ROADMAP | **Superseded**: physiological target is 0.12–0.15 (AttachedPool identity); covered by A-10 |
| Throttled-hydrolysis gap / re-target turnover | AttachedPool_vs_ATPase | **Merged into A-10** |
| Steady-state force sensitivity to the right side of the rate-transition curve | June ROADMAP | **Stale** |
| 11-DOF ceiling; `OptimizeMechanismEvaluation.m` | SensitivityAnalysis | **Done/obsolete** |
| `LastSlackIdentification` stub | June ROADMAP | **Archive** the folder (no findings to write). Close |
| `UseA2Reattaching` lever | A2Reattaching | **Merged into M0-7** |
| D0 / apo-pool detachment | ApoPoolDetachment | **Done: falsified**. Small leftovers (ktr overstretch peak as a feature; test `RegAvailShorteningOnly=0`) → fold into A-5? |
| Load-engaged SRX recoil, fixed dwell | RestretchSRXRecoil | **Done: falsified**. Leftovers: p2→SRXT never jointly fit, bound `k_cr`. **Stale?** |
| 3-state model, catch-bond sweep | RestretchFeatureFit | **Stale**: never run; no working 3-state seed. Keep as an idea? |
| Shoulder fade at 2 mM; FV-iso 56 vs slack-iso 66 gap | FV_Shoulder | **Stale?** low priority |
| Floor test over `dr`/`dr2` + PCHIP knots | FeatureSensitivity | **Stale** |
| `rigaware_b` degeneracy check | RecoveryOscillation | **Stale**: superseded by RefitStrategy |
| Settle the post-restretch window definition | RestretchVsKtrRecovery | **Stale**: rsK80 replaced it |
| 03/27 + 04/03 features from an older extractor | RestretchRecoveryFit | **Probably superseded** by the pooled target; check with A-15 |
| `FV_SHARE` 0.70→0.40 | OptimResultReport | **Done and superseded**: the floor replaced the pin |
| Score `rsT0` | OptimResultReport | **Covered by A-15** (robust `rsTfitE`) |
| `tau_reg` fast gate; ktr protocol held out | KtrProtocolVsSlackKtr | **Done** (`fastgate_best`). Leftover "refit from `fastgate_best` + re-validate" → has this happened in refit6–9? Close if yes |
| Peak-detector bug | RestretchPeakAudit | **Done** |
| rsK80 estimator | RestretchMetricAudit | **Done** |
| Objective numerical noise | StrainWindowNoise | **Done**. Leftovers in A-15; 0.5 nm bin check + "heads attached at +66 nm / SE exponent 0.845 < 1" → stale or a `science` issue? |
| ADP-trap / LowATP force enhancement mechanism ranking | LowATP_ForceEnhancement | **Superseded** by ATPEffectReconciliation + A-9 |
| PassiveForceID: Ca-induced titin stiffening check | PassiveForceID | Minor. **Stale?** |
| Lab-diary running TODOs | `Docs/notes/labdiary.md` | **Duplicated** in the June ROADMAP; covered above |
| Delete-the-dead-switches, `_archive`, refactoring plan | status 09-16 §8 #8, ROADMAP | **Filed** as M0-7, M0-10, M0-11 |

---

## 7. Questions for you (answer inline)

1. ~~Milestone B~~: decided, stays in A. ~~PR base~~: decided, PRs wait for M0.
2. ~~Paper 1 data~~: decided, **new dataset**; A blocks P1 and A's work goes on `paper1`. Board: `beards-lab` org. Default labels: delete. Filing: now.
3. Accept **M0 Branch & release model** and **P1 Paper 1** as milestones?
4. Org project under `beards-lab` (visible to collaborators) or under your user account?
5. Seed `decisions/` with the ~5 existing decisions, or start empty?
6. Delete the default labels listed in §2?
7. Which ★ items go to **Next** (I'd keep Next at about 5)?
