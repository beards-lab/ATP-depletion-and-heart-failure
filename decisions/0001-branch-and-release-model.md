# 0001 — Branch and release model

- **Status:** accepted
- **Date:** 2026-11-02
- **Issue:** #13 (milestone M0)
- **Decided by:** Filip Jezek

## Context

Until now all work since 2023 lived on long-running working branches (most recently
`PooledTarget-8mM`); `main` was 488 commits behind and not the truth. The model now serves
goals with different stability needs:

| Goal | Milestones | Consumer | Stability need |
|---|---|---|---|
| Published rodent model: high- and low-ATP fits at saturating Ca | A, P1 | reviewers, readers, anyone reproducing the paper | **Frozen**: bit-reproducible from a tag, bugfixes only |
| Ca-activated (beating), humanized myocyte | 1, α, 2 | the group and collaborators; next paper(s) | Under development; flags change often |
| Reduced cell-mechanics component for organ models (TriSeg etc.) | 3 | other groups' heart models | A stable interface (inputs/outputs/units) while internals change |
| Metabolite extensions (Pi, pH) | 4 | the humanized model and future papers | Development |

One-issue-one-branch-one-PR needs to know which branch is the base.

## Decision

Two long-lived lines, plus short-lived issue branches.

```
main ─────●───────●──────●──────●────────●──────▶   integration trunk: everything, new features behind flags (default OFF)
           \             ▲      ▲
ca-sat-1d-mouse ●──●──●──┘──●───┘──(tag ca-sat-1d-mouse-v1.0)──▶ release/ca-sat-1d-mouse  (bugfixes only)
           saturating-Ca rodent 8 + 2 mM work only; merges flow ca-sat-1d-mouse → main, never back
```

1. **`ca-sat-1d-mouse` is the paper line (Paper 1).** Saturating Ca (maximally activated, skinned
   preparation), full 1-D strain-distribution model, mouse, 8 + 2 mM ATP. It starts from
   `PooledTarget-8mM` plus `lowatp-2state` and takes milestone A's new-data fit. No Ca dynamics,
   no humanization, no reduced model.
   - Naming pattern for future lines: `<activation>-<model>-<species>`, e.g. `ca-twitch-1d-mouse`
     (beating), `ca-twitch-1d-human`.
2. **Merges flow forward only:** `ca-sat-1d-mouse` → `main`, never `main` → `ca-sat-1d-mouse`.
3. **`main` is the integration trunk.** Ca activation, humanization, the reduced model and Pi/pH
   work land here behind flags that **default OFF**.
4. **Regression guard (#14).** On `main`, with every non-paper flag at its default, the paper
   battery reproduces the `ca-sat-1d-mouse` numbers. If that test breaks, `main` has silently
   changed the published model.
5. **Release.** At submission, tag `ca-sat-1d-mouse-v1.0`; the branch continues as
   `release/ca-sat-1d-mouse`, bugfixes only. Each fix is a new tag (`ca-sat-1d-mouse-v1.0.1`)
   and merges forward into `main`.
6. **Issue branches** are named `<n>-<slug>` and fork from `ca-sat-1d-mouse` for paper work,
   from `main` for everything else. One issue, one branch, one PR (`Closes #<n>`).
7. **Repo-process docs** (`decisions/`, `ROADMAP.md`, the issue template, CLAUDE.md's workflow
   section) are committed on the paper line while it is the most advanced branch and reach
   `main` through the normal forward merge. Docs that only concern later goals live on `main`.
8. **Later:** once the reduced model has an external user it gets its own version tags on `main`
   (`<model>-v0.x`, name decided then).

A release must also carry what a tag alone cannot (tracked in #66): a data manifest with
checksums (data are not in git), one entry script that regenerates every paper figure and table,
a frozen parameter snapshot under a name the optimiser can never overwrite, and a Zenodo DOI.

## Alternatives considered

- **Keep working on one branch** (status quo): no stable base for PRs, and the published model is
  not protected from feature work.
- **Everything on `main`, paper as a tag only:** cheapest, but bugfixes to the published model
  would drag along unreleased features, and paper work would be gated on the flag-off guard
  being complete.
- **Name `paper1`:** rejected for a name that states what the model is (saturating Ca, 1-D,
  mouse), so sibling lines (beating, human) read as a family.

## Consequences

- PRs are unblocked once `ca-sat-1d-mouse` and an up-to-date `main` exist (#2, #4, #36, #48).
- Every flag added on `main` must default OFF and must be covered by the regression guard.
- Feature work that the paper needs goes on the paper line first, then forward to `main`.
