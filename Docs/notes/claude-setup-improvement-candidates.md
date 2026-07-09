# Claude Code setup — improvement candidates

_Generated 2026-07-07 by reflecting on 15 past Claude Code session transcripts_
_(`C:\Users\fjezek\.claude\projects\C--home-git-ATP-depletion-and-heart-failure\*.jsonl`)._
_Four sub-agents mined raw signals; signals were clustered and each cluster given a decision:
**FIX** (correctness/config), **DOC** (CLAUDE.md/memory), **SKILL**, **AUTOMATION**, or **NOTHING**._

Sessions are cited by their 8-char id prefix. A signal seen in ≥3 of the 4 mining batches is marked
**(pervasive)**.

---

## Ranked candidates

### 1. CLAUDE.md "Running MATLAB" section is stale — should be MCP-first — **FIX / DOC** ⭐ highest leverage
**(pervasive — all 4 batches)**

CLAUDE.md documents the WSL interop recipe
(`"/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe" -batch ...`) as *the* way to run MATLAB.
In practice that path is dead: across sessions the MATLAB **MCP** (`mcp__matlab__evaluate_matlab_code`
/ `run_matlab_file`) was called **~49 / 60 / 64 / 279 times** and succeeded almost every time, while the
documented WSL `matlab.exe -batch` command was **0/8 successful** (all Exit 1) and `/mnt/c/` paths were
never used successfully. Agents repeatedly burned turns trying 4+ path/cd variants before pivoting to the MCP.
_Sessions:_ 96d7a07b, cda98489, 3120450a, 4beb8e99, 526b39e2, ab7721c1.

**Action:** Rewrite the "Running MATLAB" section to mandate the MATLAB MCP tools as the primary
execution path (`mcp__matlab__run_matlab_file` for scripts, `evaluate_matlab_code` for snippets).
Keep the WSL `matlab.exe -batch` recipe only as a clearly-labeled fallback, or delete it. Note the MCP
runs `nodesktop` against `R2023a` (per `~/.claude/settings.json`).

---

### 2. `git` not on PATH; Bash tool is unreliable in this environment — **FIX / DOC**
**(pervasive — all 4 batches, plus reproduced live this session)**

`git` is not resolvable from PowerShell or the Bash tool → repeated
`The term 'git' is not recognized`. Agents guessed wrong paths for ~5 turns until the user intervened
("...you have an env var with git path. We are using github desktop's git.exe"). Separately, the **Bash
tool itself fails on nearly every call here** — even `ls`, `echo`, `cat`, `grep`, `cp` return Exit 1/2
(reconfirmed live in this very session: `ls`/`cd &&` Bash calls failed; PowerShell + Grep/Glob/Read worked).

**Verified facts (this session):**
- `git` is NOT on PATH. `C:\Users\fjezek\AppData\Local\GitHubDesktop\bin` is on PATH but does not expose `git.exe`.
- A working git.exe is already pointed to by env var `CLAUDE_CODE_GIT_BASH_PATH` →
  `C:\Users\fjezek\AppData\Local\GitHubDesktop\app-3.5.12\resources\app\git\cmd\git.exe`
  (version-pinned folder — fragile; app-3.5.11 also present).

**Action:**
- Add a CLAUDE.md "Shell & git" note: **prefer PowerShell + the Grep/Glob/Read tools; avoid the Bash tool**
  in this environment. For git, use the full path via `$env:CLAUDE_CODE_GIT_BASH_PATH` (version-proof) or
  install Git for Windows and put it on PATH permanently.
- Optional: set a stable `git` alias/PATH entry in a machine profile so the version-pinned GitHubDesktop
  folder isn't referenced directly.

---

### 3. MATLAB MCP has no persistent workspace/path/compiled-code state — **FIX + small AUTOMATION**
**(pervasive)**

The MCP does not guarantee a persistent workspace between calls, so:
- ~7 calls/session re-run `cd(root); addpath(genpath('.'))` boilerplate; still hit
  `Unrecognized function or variable 'getParams'` and stale-`params0` leaks
  (`Unrecognized field name "justPlot"`, state-vector length mismatches 187 vs 203). _cda98489, 7f1a81f3._
- **Parpool workers + `evaluate_matlab_code` cache stale compiled code after edits** — the user called this
  "a common problem" and had a sonnet agent build `loadParams.m` + `refreshPool.m`
  (restart parpool + `clear functions; rehash`) to fix it. _ab7721c1, 526b39e2, 173576af._

**⚠ Verified gap:** `Auxiliary/loadParams.m` and `Auxiliary/refreshPool.m` **do not currently exist in the
repo** (memory note `lowatp_mechanisms.md` refers to them as an "infra fix"). The fix was built mid-session
and appears to have been lost / never committed.

**Action:**
- Re-create and **commit** `loadParams.m` (safe param-snapshot loader that doesn't get overridden by
  defaults / the `base=` NaN bug) and `refreshPool.m` (restart parpool + `clear functions; rehash`).
- Add a one-line `startup`-style prelude helper (`ensurePath.m`) that idempotently sets `cd`+`addpath`, and
  document in CLAUDE.md that every fresh MCP session must call it first.
- Record the "edit-then-refreshPool-or-workers-run-stale-code" gotcha as a durable memory + CLAUDE.md note.

---

### 4. Core inner loop "load params → run experiment → dump features → analyze residuals" run by hand dozens of times — **SKILL**
**(pervasive)**

The canonical loop recurs constantly, hand-typed each time:
`root=…; cd; addpath(genpath); refreshPool(5); p0=getParams(loadParams('params_…'),[],true,false);
RunBakersExp / runSlackExperiment; extract features; compare vs data`. Users re-issue
"Run the current best-fit, analyze the restretch/FV residuals, propose next mechanism" repeatedly
(U20/U10/U23 in 4beb8e99; U1 in 5d03d446; many near-identical MCP payloads in ab7721c1/526b39e2).

**Action:** Create a **skill** (or a single driver function `runFitReport.m` behind it) that takes a param
snapshot name + which experiments to run, and returns/writes: all `cost_vect` members, feature-vs-data
table, and residual summary. Bakes in `refreshPool` + path setup. This collapses the single most-repeated
manual workflow. (Fold in the `reportFeatureCost.m` printout the user always asks for — 96d7a07b U319.)

---

### 5. Figures from headless MATLAB MCP frequently don't render — **FIX / AUTOMATION**
**(3 batches)**

Users repeatedly had to re-ask for plots: "show me the figure. It has not been generated",
"First off, give me the figures - those did not render." PNG panel-merge / montage was re-implemented
bespoke each time (Pillow path dead-ended → rewritten in MATLAB). _526b39e2, ab7721c1, 96d7a07b, 4beb8e99._

**Action:** Add a reusable `exportFeaturesFigure.m` that always writes a PNG to a known results dir
(`exportgraphics`/`print`, explicit `-r`), plus a `montagePanels.m` grid helper. Make the skill in #4 always
export PNGs rather than relying on interactive rendering. Prefer `SendUserFile` on the resulting PNG.

---

### 6. Long-running optimizers can't be hosted in-session; "smoke locally, submit to grid" should be the sanctioned path — **DOC / SKILL**
**(3 batches)**

Hosting multi-hour optims from a session repeatedly failed: `run_in_background` tasks killed ~3 min in;
detached `matlab.exe -batch` launch died in 34 s (Windows launcher-handoff, empty logs); `sleep`-then-`cat`
polling is **blocked by the harness** ("use Monitor"); Monitor `until grep` loops then timed out (because
Bash is broken). User's fallback: "Run for ~30 min at most just to debug it runs fine", then run on the grid.
There is already `scripts/rol.sbatch`. _4beb8e99, 526b39e2, 96d7a07b._

**Action:** Document the sanctioned workflow in CLAUDE.md: **(a)** debug-smoke the optimizer in-session
(short eval budget) via the MCP, **(b)** submit the real run to the grid via `scripts/rol.sbatch`, **(c)**
for any in-session polling use `Monitor` with an until-loop, never `sleep`+`cat`. Consider a small skill
wrapping the smoke-test + sbatch submission.

---

### 7. Optimizer gotchas + best-params save bug — **FIX (correctness) + DOC**

- **Best-params save mismatch (real bug):** "RunRestretchOptim reported decrease to ~6, but the SNAP_BEST
  file really ends in costs of ~160" — reported cost ≠ saved best-params. _74638a09 U50._ Needs a code fix
  in the optimizer's snapshot/save logic.
- **`surrogateopt` gotchas to codify:** `UseParallel=true` silently no-ops on a **thread** pool (falls back
  to g=1 → stalled optimizer masquerades as running) → use `UseParallel=false`; `MinSurrogatePoints`
  default 20 errors when eval budget < 20. _4beb8e99._

**Action:** Fix the snapshot-save logic so the reported best and the persisted `SNAP_BEST` agree
(add an assertion). Record the two `surrogateopt` gotchas in CLAUDE.md's optimization section.

---

### 8. Encode standing domain + working-style rules — **DOC / MEMORY**
**(pervasive preferences)**

Recurring user corrections that should be standing rules rather than re-litigated each session:
- **Physiology-first:** reject non-physical fixes; α-myosin ⇒ `kah ≥ 80`, no reverse hydrolysis,
  "everything must have physiological meaning"; don't naively temperature-scale ktr; ML = SL/2. _96d7a07b._
- **Delegation preference:** "spawn cheaper/haiku/sonnet agents for grunt work; orchestrate & plan yourself"
  (stated 5+ times across batches). Follow an explicit **OODA** loop.
- **Anti-over-engineering:** merge duplicate logic, keep code lean; bounds/eval logic must live *inside* the
  features/RunBakersExp path, not as a disconnected cost term. _f46e2fc4._
- **Distrust cheap wins:** down-tuning weights ≠ progress; prove "iron wall" claims or find a real mechanism.
- **Minimal, daisy-chainable param overlays**, not full param dumps. **Reason before simulating.**
- **Cost/quota-aware:** heavy manual `/model` switching (~16 across sessions) to manage cost per phase.

**Action:** Add a "Working agreement" section to CLAUDE.md capturing the domain constraints + the
delegation/OODA/anti-over-engineering preferences; store the durable ones as `feedback`-type memories
(several already exist — extend, don't duplicate).

---

### 9. SRX steady-state fraction tuning loop — **SKILL** (can fold into #4)
**(2 batches)**

Tuning SRX to "~0.2 at max force, equilibrate ~10 ms" was hand-rebuilt and **mis-iterated ~4×** against a
simple target (agent mis-estimated the SRXT+SRXD total three times: "0.822+0.116 is almost 1, not 0.2").
_96d7a07b, cda98489, 7f1a81f3._

**Action:** A dedicated `checkSRXfractions.m` harness: `getParams` → apply SRX overlay only →
`RunSlack=0, FV_velocities=[0]` → `RunBakersExp` → report steady-state SRXT/SRXD. Removes the arithmetic
foot-gun. Natural sub-command of the #4 skill.

---

### 10. Context-window exhaustion on long multi-file campaigns — **NOTHING new**

Frequent `/compact`, "continue from where you left off", multi-summary resumes. _f46e2fc4, 96d7a07b._
Already mitigated by the user's `labdiary_*.md`, `memory/*`, and per-analysis `conclusions.md` discipline.

**Action:** None new — reinforce keeping those load-bearing notes current. (No tooling gap justifies a build.)

---

### 11. Cross-session provenance search & minor tool gotchas — **NOTHING**

- `search_session_transcripts` returned nothing useful for "which convo built parameterBounds.m" — a tool
  limitation, not something to build around. _f46e2fc4._
- Scattered singletons: Monitor schema not pre-loaded via ToolSearch; `AskUserQuestion` rejects >4 options;
  `Read` on a directory → EISDIR; Edit "String to replace not found" on unicode arrows. General
  Claude-usage patterns, low leverage individually.

**Action:** None. Noted for awareness only.

---

## Summary table

| # | Cluster | Decision | Priority |
|---|---------|----------|----------|
| 1 | CLAUDE.md MATLAB run instructions stale (WSL vs MCP) | FIX/DOC | ⭐ highest |
| 2 | git not on PATH + Bash tool unreliable | FIX/DOC | high |
| 3 | MATLAB MCP state persistence + loadParams/refreshPool (missing) | FIX + automation | high |
| 4 | "run → feature dump → residual analysis" driver | SKILL | high |
| 5 | Reliable figure/PNG export | FIX/automation | medium |
| 6 | Smoke-locally-then-grid for long optims | DOC/SKILL | medium |
| 7 | Best-params save bug + surrogateopt gotchas | FIX/DOC | medium (bug = high) |
| 8 | Standing domain + working-style rules | DOC/MEMORY | medium |
| 9 | SRX steady-state fraction harness | SKILL | medium |
| 10 | Long-session context management | NOTHING | — |
| 11 | Provenance search / minor tool gotchas | NOTHING | — |
