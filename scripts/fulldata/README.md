# Sharing the model *with* the data — the private full-data repository

The public repository excludes `data/` (`.gitignore` line 28), because the
recordings are not published yet. Colleagues who need to actually run the model
get a **separate private repository** that carries the same code plus the curated
protocol files.

- **public**  `github.com/beards-lab/ATP-depletion-and-heart-failure` — code only
- **private** `github.com/beards-lab/ATP-depletion-and-heart-failure-full-data` — code + data

Everything in this folder is the machinery that keeps the second in sync with the
first.

| File | Role |
|---|---|
| `Seed-FullDataRepo.ps1` | One-time setup. Creates the local share clone and pushes both branches. |
| `sync-code.yml` | GitHub Actions workflow. **Lives on the private repo's `main`**, copied there by the seed script. |
| `../Update-FullDataPayload.ps1` | Refreshes the *data* payload when a `.mat` is regenerated. |

---

## Why not just a fork, or a cross-repo pull request?

Neither works on GitHub:

- **A fork of a public repository is always public.** A fork's visibility cannot
  be changed to private, so the data could not live there.
- **Pull requests require both repositories to be in the same fork network.** An
  independently created private repository is not a PR target for the public one.

So the private repository carries *two* branches and the pull request happens
**inside it**, which gives the same one-click review-and-merge experience.

---

## Branch layout in the private repository

| Branch | Content | Written by |
|---|---|---|
| `code` | A pure mirror of the public repo's `Further-optim`. **No data.** | The sync workflow, only. |
| `main` | `code` + the data commit + this workflow. **The default branch — what colleagues clone.** | You (data) and PR merges (code). |

The data commit touches only files under `data/`, which never exist on `code`, so
merges from `code` into `main` cannot conflict on data. Code files merge cleanly.

> **Rule: never commit code changes directly to `main`.** Code reaches `main` only
> by merging `code`. That is what keeps `main`'s code content always identical to
> `Further-optim` rather than a divergent fork.

---

## Routine update — browser only

Whenever `Further-optim` moves in the public repository:

1. Private repo → **Actions** → *Sync code from public repo* → **Run workflow**
2. Private repo → **Pull requests** → open the PR it created → **Merge**

The workflow fetches the public repo (no credentials needed — it is public),
force-updates `code` to match `Further-optim`, and opens the PR. You review the
code diff before colleagues see it.

Because `main` only ever gains merge commits, **collaborators just `git pull`** and
it fast-forwards. No force-pushes, and it works in GitHub Desktop.

### Updating the data instead

Run `scripts/Update-FullDataPayload.ps1`, or do it by hand in GitHub Desktop:
copy the changed `.mat` files into the share clone's `data/`, review, commit to
`main`, push.

---

## Safety

The share clone has the **public remote removed** — its only remote is the private
repository. It therefore has no way to push data to the public repo, which is why
no hooks or push guards are needed, and why it is safe to open in GitHub Desktop.

Do **not** add the public working copy's `origin` back to the share clone.

The payload is an explicit allowlist (`data/protocol_*.mat`, `data/bakers_*.mat`,
`data/isovelocity.mat`, minus the example's generated file — **60 files, 2.4 MB**;
see `payload-manifest.txt`). The scripts stage the resolved per-file list, not the
globs, so what you see printed is exactly what gets committed, and they refuse any
payload over 50 MB. **Never `git add -f data/`**: that
would sweep in the ~6.6 GB raw layer under `data/PassiveCaSrc2/`, including four
zips of 217–841 MB that GitHub rejects outright (100 MB per-file limit). The model
does not read the raw layer — it is input to `DataCuration/` only.
