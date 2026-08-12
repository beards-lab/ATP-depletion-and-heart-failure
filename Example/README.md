# Example — start here

A working, self-contained entry point to the cardiac cross-bridge model: load a
frozen parametrization, run four mechanical protocols against one day's 8 mM ATP
recordings, and see how well the model fits.

Everything here is meant to be read as much as run. The scripts are commented at
tutorial density and each one calls out the traps that are not visible from the
code alone.

| File | What it does |
|---|---|
| [`RunExampleHighATP.m`](RunExampleHighATP.m) | **The main driver.** The high-ATP protocol battery: slack, ktr, staircase, force-velocity. |
| [`BuildExampleProtocol.m`](BuildExampleProtocol.m) | Builds a protocol of your own design (hold → stretch → hold) from scratch. |
| [`RunExampleProtocol.m`](RunExampleProtocol.m) | Runs that protocol two ways and checks the round trip is lossless. |

---

## 1. Quick start

MATLAB R2023a or newer. No build step, no toolbox required for a first run (the
Parallel Computing Toolbox only speeds things up later).

```matlab
cd('<repo root>'); addpath(genpath('.'));
Example\RunExampleHighATP
```

Takes **1–2.5 minutes** (the first run is slower — MATLAB is compiling). You should see:

- **figure 1** — one tile per protocol, data in black, model in blue.
- **figure 80085** — one tile per scored feature (data ○ vs model ×), plus a
  sorted *cost culprits* bar chart. Read that bar chart first: it tells you where
  the fit is actually losing.
- in the command window, a cost breakdown ending in:

```
==> OPTIMIZER OBJECTIVE (feature cost, costExp = 2): 7.502
```

**That 7.5 is the reference number.** If you get it, your environment reproduces
the run this package was built from. (`sum(E)` printed just above it is ~210 —
that is a diagnostic total, not the fitted objective. See §5.)

Then:

```matlab
Example\RunExampleProtocol
```

which should end in `PASS — build -> save -> load -> run is lossless.`

---

## 2. What the model is

A strain-discretized cross-bridge model of a cardiac sarcomere, after
[Beard et al. 2022](https://www.sciencedirect.com/science/article/pii/S0006349522006026).
Myosin heads occupy a small number of chemo-mechanical states; within each state
the population is resolved as a probability *distribution over strain*, and the
whole thing evolves as a large system of ODEs. Force is the strain-weighted sum
over the attached states, plus passive (titin) and serial-elastic contributions.

Muscle length is the input. You command a length trajectory; the model returns
force.

---

## 3. The four protocols

All four are mechanical perturbations of a maximally activated fibre. Each
isolates a different aspect of cross-bridge kinetics.

| Protocol | The manoeuvre | What it constrains | Driven by |
|---|---|---|---|
| **Slack** | Repeated cycles: fast release to slack → hold short → re-stretch → hold | The richest protocol — force levels, redevelopment rates, re-stretch peaks. Most fitted features come from here. | `velocitytableonfile` |
| **Ktr** | Fast release then re-stretch, then watch force redevelop | How fast cross-bridges re-attach from an emptied bound pool | `ktr_velocitytableonfile` |
| **Staircase** | A series of small stretches ("ramp up") | Passive plus short-range elastic response | `stairs_velocitytableonfile` |
| **Force-velocity** | Steady force at each of several constant shortening velocities | The classic Hill relation; sets the strain-dependence of detachment | `params0.FV_velocities` + reference data inlined in `Model/LoadData.m` |

A fifth block also runs by default:

| **Passive** | The slack protocol with active attachment switched off (`ka = 0`) — the experimental PNB + mavacamten condition | Titin and serial-elastic mechanics alone | `passive_velocitytableonfile` |

> **Do not switch the passive block off casually.** It contributes no cost of its
> own — it contributes five `PS_*` *features*, and the shipped feature list asks
> for all five. Turning it off while the feature list is unchanged makes those
> features **missing**, at a flat penalty of 100 each. The cost jumps by ~500 and
> looks like a physics change when it is really bookkeeping.

---

## 4. Data files

Protocol data lives in `data/` (which is **not** in the public repository — see
the note at the end). Each protocol is one `.mat` holding two tables:

**`velocitytable`** — the driving signal.

```
[ t(s), v(ML/s), v_um(bookkeeping), L(ML) ]
```

Row *i*'s velocity is held constant over `[t_i, t_i+1]`. Row 1 is a warm-start
hold that lets the cross-bridge population settle before the protocol begins.

> **The model is driven by column 2 only.** Columns 3 and 4 are bookkeeping —
> nothing reads them during a simulation, so they can silently disagree with the
> integral of column 2 if a row was edited. When you need the commanded length,
> integrate column 2 (see `iCommandedL` at the bottom of `BuildExampleProtocol.m`).

**`datatable`** — the measurement.

```
[ t(s), SL(um), F ]
```

`F` is **absolute kPa** in the slack and passive files, and **normalised by the
pre-event isometric force** in the ktr, staircase and FV2 files. The runners know
which convention theirs uses; you only need to care when building your own file.

### Available datasets

Set them in cell `[2]` of the driver. The alternatives are all present and
commented out — uncomment one block to switch. All five have been verified to run.

| Day | 8 mM ATP | 2 mM ATP | Passive (PNB+Mava) | Notes |
|---|---|---|---|---|
| 03/27/2026 Male | ✔ (**default**) | ✔ | ✔ | The set this parametrization was fitted to |
| 04/03/2026 Female | ✔ | ✔ | ✔ | Second fibre, same protocol design |
| 04/10/2026 Male 2-8 | ✔ | ✔ | ✔ | **Reversed ATP order** (2 mM run first), which separates the real ATP effect on ktr from the rundown that contaminates the force ratio. Only day with FV2. |
| Legacy Baker | ✔ | ✔ | — | Coarser; kept for continuity with older analyses |

**FV2** is an isovelocity ramp at 2 ML/s embedded in the 04/10 recordings — a
genuine sub-Vmax force-velocity point. It is carved into its own file so the slack
files keep their strict row layout, and because its datatable uses the
normalised-force convention it is driven **through `runKtrExperiment`**.

---

## 5. Reading the cost

Two numbers come out, and confusing them is the most common early mistake.

**`sum(E)` ≈ 210** — a diagnostic total. `E` is a positional row vector: one
trace-MSE term per active protocol, then one term per scored feature. It is
dominated by the slack trace MSE (~178), which is a per-sample squared error over
a long trace, scaled by 20, that **nothing ever tried to minimise**.

**Feature cost ≈ 7.50** — the real objective. The shipped snapshot sets
`params0.OptimizeOn = 'Feats'`, so the optimizer minimises the feature cost alone.
A large `sum(E)` beside a small feature cost is the expected picture.

> `E` is **positional and protocol-dependent**. `RunBakersExp` appends terms in a
> fixed order — FV, ktr, staircase, slack, force-length, then the whole feature
> vector — so switching a protocol off silently shifts the index of everything
> after it. Never hard-code an index into `E`; rebuild the labels from the flags,
> as cell `[5]` does.

Two more things that surprise everyone once:

- **ktr and staircase contribute a scalar MSE only.** They produce no features, so
  they never appear in the feature plot.
- There *is* a feature called `ktr` in that plot — but it is measured from the
  **slack** protocol's re-stretch recovery, not from the ktr protocol.

---

## 6. The parameter system

Everything flows through one struct, and `getParams` is its single source of truth.

```matlab
params = getParams(params, g, updateInit, updateModifiers)
```

| Argument | Meaning |
|---|---|
| `params` | Partially filled struct; missing fields get defaults |
| `g` | Vector of **multiplicative modifiers** applied to the fields named in `params.mods`. This is what the optimizer tunes. Pass `[]` when not optimizing. |
| `updateInit` | Recompute the strain grid `params.s` and the initial state `params.PU0`. **Always `true` after changing `SL0`** — forgetting it is the classic silent failure: the model runs, but from the wrong initial condition. |
| `updateModifiers` | Apply `g`. `false` when loading a saved snapshot, which already has its modifiers baked in. |

`getParams` also resolves **expression fields**: any value that is a string
starting with `'='` is evaluated as MATLAB, e.g. `params.kamh = '=0.1*kah'`. That
keeps derived parameters derived instead of drifting out of sync.

### Loading a snapshot

```matlab
params0 = getParams(loadParams('params/ExampleHighATP.m'), [], true, false);
```

**Always go through `loadParams`, never `run` or `eval`.** `loadParams` executes
the snapshot in its own sandboxed workspace; running it in the base workspace lets
a stale `params0` leak fields into the new one, and defaults then silently
override the snapshot.

Give the **path**, not a bare name. `loadParams` resolves bare names via `which()`,
and this repository has snapshot names that exist in more than one folder — most
notably `Workbench/ThursdayNightFever.m` shadows `params/ThursdayNightFever.m`,
and the two are **not** identical.

### The shipped parametrization

`params/ExampleHighATP.m` — a 2-state model at 8 mM MgATP. Content copied verbatim
from `params/rskR2_w025_opt.m` (2026-08-10). The name is deliberate: optimizer
drivers write `params/<tag>_opt.m` and `params/<tag>_iter/`, so a machine that
pulls this repository and runs an optimization can never collide with this file.

---

## 7. Features and the `params0.fn` list

A **feature** is a scalar (or one scalar per protocol cycle) distilled from a force
trace — a peak height, a recovery rate, a steady level. Fitting features rather
than raw traces is what keeps the optimizer honest about *which part* of the trace
it is getting wrong.

`params0.fn` is the list of what gets scored and how heavily. The grammar (full
version in `Model/evalFeatureCost.m`):

| Entry | Meaning |
|---|---|
| `'field'` | Score this feature, weight 1 |
| `'field\|w'` | ... with weight `w` |
| `'field_y\|field_x'` | Scatter *y* against *x*, e.g. `'FV_fnorm\|FV_v'` |
| `'field\|lb-ub'` | **Boundary**: penalise only outside `[lb, ub]`. Scored against the **simulation alone** — no data is read. This is how physiological plausibility is enforced on quantities nobody measured. |
| `'f1\|w1,f2\|lb-ub'` | A comma at top level makes one **grouped** joint cost, rendered as a single tile |
| `'field[sel]'` | Reduce a per-cycle vector: `[1] [end] [mean] [min] [max] [median] [first] [last]` |

Penalties: each `NaN` costs 1; a feature the model never produced costs **100**. A
single 100 in the breakdown almost always means a missing protocol, not bad physics.

Useful helpers:

- `plotFeatures(features_data, features_model, [], fn)` — figure 80085, per-feature
  tiles plus the cost-culprit bar chart.
- `reportFeatureCost(features_data, features_model, fn)` — the same as text, plus a
  raw data-vs-model table of *every* feature, including unscored ones.

---

## 8. Known residuals

The shipped parametrization is the current best, not a solved problem. On the
8 mM 03/27 dataset:

| Feature | Data | Model | Status |
|---|---|---|---|
| `ktr` | ~50 /s | ~57 /s | Force redevelopment is **too fast**. A long-standing 2-state limitation: turnover, attached fraction and ktr are welded together by the core cycle. |
| `rsK` | ~42 /s | ~110 /s | Post-re-stretch recovery is **~2.6× too fast**, and its amplitude dependence has the **wrong sign** (model +714 vs data −126). This is the open structural defect. |
| `FV_fnorm` | 0.66 at 1 ML/s | 0.57 | The force-velocity curve falls off too steeply at intermediate velocities |
| `coolDownLS` | 0 | ~1.2 | A residual the current parametrization does not clear |

`PS_rampupRMSE` and `PS_holdDecayRMSE` are model-vs-data *residuals*, so their data
value is 0 by construction — they are not features you can "match".

Background: `Analyses/README.md` indexes the investigations behind each of these.

---

## 9. Building your own protocol

`BuildExampleProtocol.m` builds a hold → stretch → hold from four numbers and runs
it, in about 15 lines of actual code:

```matlab
vt = [ -T_WARM, 0, 0, L0        % warm-start hold
        t_str,  V, 0, L0        % stretch ramp
        t_hold, 0, 0, L0 + dL   % hold at length
        t_end,  0, 0, L0 + dL ];

p = params0;
p.SL0 = 2.0;  p.Slim_l = 1.4;  p.Slim_r = 2.3;   % grid must span the excursion
p = getParams(p, p.g, true);                     % MUST re-run after SL0
p.Velocity = vt(:, 2);
[~, out] = evaluateModel(str2func(p.modelFcn), vt(:, 1), p);
```

`RunExampleProtocol.m` then saves that as a protocol file and re-runs it through
`runStairsExperiment`, which scores it against the trace the direct run produced.
The cost must come out at **exactly 0** — that is the self-test proving the
build → save → load → run round trip is lossless. Swap in real measurements and the
same cost term becomes a real fit.

**Why `runStairsExperiment` as the host?** Of the runners it is the one with no
structural expectations: plain trace MSE, no feature extraction, no assumptions
about row count. `runSlackExperiment` *cannot* host an arbitrary protocol — it
segments the table on velocities below −1 ML/s (looking for the fast release of a
slack cycle) and its internal `switch` indexes rows literally, so a table of the
wrong shape is silently mis-read.

There is a small discovery built into the example: the peak (yield) force during
the ramp is set by the stretch **velocity**, not by how far you stretch. Change
`dL` from 0.05 to 0.01 and the peak stays at ~72.9 kPa; change `V` and it moves.

---

## 10. Troubleshooting

**`Unable to read file 'data/protocol_...mat'`**
The current directory is not the repository root. Every runner addresses data as
`data/<file>`. The Example scripts `cd` to the root themselves; if you are calling
runners by hand, do the same.

**`Reference to non-existent field 'FV_velocities'` (or `fn`, or `OptimizeOn`)**
These three have **no default in `getParams`** — they exist only if a snapshot sets
them. A from-scratch `getParams()` with `RunForceVelocity = true` will error. Load
a snapshot, or set them by hand.

**A feature costs exactly 100**
It is missing, not badly fitted. The protocol that produces it did not run — most
often `RunSlackPassive` was switched off while the `PS_*` entries stayed in `fn`.

**Edits to model code appear to have no effect**
MATLAB caches compiled code, and parpool workers keep their *own* copy. Run
`refreshPool(5)` (or `refreshPool(0)` for no pool) after editing anything under
`Model/` or `params/`. This is the single most common foot-gun in this project.

**`RunSlackSegments = 'AllPar'` hangs or errors**
It needs a live parallel pool. Use `'All'` for the serial equivalent — that is what
the Example uses, and it also gives a clean contiguous trace for plotting.

**The run writes `Ghost_*.mat` files into my working directory**
"Ghosts" are saved reference traces the plots overlay for before/after comparison.
Set `params0.ghostSave = ''` and `params0.ghostLoad = ''` to disable, as the Example
does.

**Everything is slow**
`params0.MaxRunTime` caps wall clock per simulation (320 s here). Drop the number of
protocols, or switch the slack protocol to `'AllPar'` with a pool running.

---

## 11. Where to go next

| Directory | Contents |
|---|---|
| `Model/` | Core model. `getParams.m`, `evaluateModel.m`, `RunBakersExp.m`, the ODE function `dPUdT_CombinedTransitions.m`, and `experiments/run*Experiment.m` |
| `Analyses/` | One self-contained folder per investigation, each with its own `conclusions.md`. Start at `Analyses/README.md`. |
| `Docs/` | Knowledge base — start at `Docs/README.md` |
| `Workbench/` | The playground: optimizer entry points, hand-tuning, ad-hoc tests |
| `DataCuration/` | Raw recordings → the model-ready `.mat` layer |
| `params/` | Parameter snapshots |

---

## A note on the data

The `data/` directory is excluded from the public repository — these recordings are
not published yet. This package is distributed through a **private** repository,
`beards-lab/ATP-depletion-and-heart-failure-full-data`, which carries the curated
protocol files alongside the code.

If you cloned from the public repository and `data/` is empty, that is why — ask
Filip (fjezek@umich.edu) for access to the full-data repository.

If you are maintaining that repository, see
[`scripts/fulldata/README.md`](../scripts/fulldata/README.md).
