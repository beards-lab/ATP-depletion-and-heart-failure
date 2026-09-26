# Noise-fragile features: the data share the fragility — fix the reader and one definition, not the tolerance

**Status:** findings + implemented (opt-in), 2026-09-24. Model side measured on `params/refit9_opt.m` (5 solver
settings × 9 parameter points 1e-6 apart, `AllPar` semantics); data side on the three 2026 preparations at 8 and
2 mM (40 surrogate recordings each). The fix is in `Model/` behind two switches that are **off by default and
bit-identical when off** (verified on data and model traces, [TestStagedExtractor.m](TestStagedExtractor.m)); the
seed that uses them is `params/robust1_seed.m` — see [Implemented](#implemented).

The question (fjezek): a feature that jumps when the parameters move by a part in a million is probably a bad
feature, and may be just as fragile against the recordings' own noise. Which features are they, is that true, and
what removes the noise without paying more than ~2× in run time?

## Verdict

| Question | Answer |
|---|---|
| Which features carry the objective's numerical noise? | **rsT20, rsK80, ktrfit_shape** — the noisy terms that pay at refit9 — plus **restretchSlopeStart**, which a 12 % phase sensitivity flips in and out of its deadband (0.078 or 0; see below), and **peak1_t, peak1_dSL, vall_t**, fragile but parked inside their guard deadbands. `t0_crossing` is quiet but, **in every campaign evaluation, its cycle 1 is a constant 2.2 ms artefact** (§1b). |
| Are they bad features in the data too? | **Mostly yes.** The between-heart spread of rsT20, rsTfit, vall_t, rsK80 is **equal to their own measurement noise** (noise/between 0.9–1.1); peak1_* and restretchSlopeStart are 66–77 % noise. Noise also *biases* the data value: rsT20 by **−11 %** — a target a noise-free model cannot reach. The well-measured ones are the values (steady, peak1_y, vall_y: noise ≤ 2 %) and the hold rates (ktrfit 1.6 %, noise/between 0.12). |
| Where does the numerical noise come from? | **Mostly the reader, then the definitions; once both are fixed the tolerance barely matters.** `extractSlackAttributes` reads the model at 1 kHz on a grid anchored on the trace's first stored time (or on raw solver steps, depending on the median step), with windows sized in samples and times measured from each window's first sample. The data are read at 4 kHz. |
| Proof that it is the reader? | "Dense everywhere" and "dense only in the cycles" hand the cycles a **bit-identical** state (max \|ΔPU0\| = 0), yet rsT20's cost moved 0.052 → 0.035: `storeOutputs` drops the first stored point, so the 1 kHz grid's *phase* moved. Reading the same saved traces on a grid anchored on absolute time takes ktr80's noise **6.25 % → 0.00 %**, rsK80's **0.97 % → 0.20 %**. |
| Does the reader also bias the fit? | **Yes.** The same model trace read at 1 kHz instead of the data's 4 kHz shifts restretchSlopeStart **27 %**, peak1_t **19 %**, peak1_dSL **8 %**, vall_t **5.7 %**, rsT20 **4.4 %**. Under the dense-in-cycles setting `restretchSlopeStart` pays 0.078 (27 % of the reference slack cost) **entirely** from this: read like the data it is inside its deadband. |
| Can the noise go without paying > 2×? | **Yes — at 1.16× the old run time.** Reader fix + rsT20 → rsTfitE: slack-objective noise **0.37 %** at `RelTol 1e-3` with the dense grid kept *inside the cycles* (measured 1.16× wall, 1.02× critical-path RHS calls, −0.5 % from the converged cost), and **0.36 %** with no dense grid at all (today's objective at that setting: **22 %**; at the old numerics 6.8 %). The tolerance itself is cheap: 1e-3 needs **1.02–1.05×** the old critical-path RHS calls. The previous fix's 2.4–4.5× came from putting the 0.5 ms output grid on the long isometric holds — the 94 s serial init hold alone costs 17 s of post-processing on the critical path — not from the tolerance. |

## 1. The reader is the main noise source

`extractSlackAttributes` ([Model/extractSlackAttributes.m:108-126](../../Model/extractSlackAttributes.m#L108-L126)) takes
the model trace in whatever form `evaluateModel` stored it and decides from `1/median(diff(t))` of the *whole* trace
(which the 94 s pre-slack hold dominates) whether to resample; when it does, the 1 kHz grid starts at the trace's
first stored time. Then every feature is read with sample-count machinery — 5 ms movmedians of `round(5 ms/dt)`
samples, a slope from "≥ 3 samples in the first 1 ms", times measured from the first sample of each window, the
location of an argmin/argmax sample. So the value depends on where the solver happened to step, and on an arbitrary
grid phase.

Same saved traces, three reads ([AnalyzeModelNoise.m](AnalyzeModelNoise.m)); numerical noise = relative SD over the
9 points 1e-6 apart, RMS over scored cycles:

| feature | A old numerics: pipeline → anchored 4 kHz | B `1e-3`, step output: pipeline → anchored 4 kHz | C `1e-3`, dense in cycles: pipeline → anchored 4 kHz |
|---|---|---|---|
| ktr80 | 0.88 % → 0.02 % | **6.25 % → 0.00 %** | 0.00 % → 0.00 % |
| ktrfit (→ ktrfit_shape) | 0.69 % → 0.00 % | **6.49 % → 0.00 %** | 0.00 % → 0.00 % |
| rsK80 | 0.97 % → 0.20 % | 0.68 % → 0.05 % | 0.05 % → 0.05 % |
| rsT20 | 3.38 % → 2.38 % | 2.06 % → 1.62 % | 0.22 % → **1.22 %** |

The rsT20 row is the exception that proves §2: a finer read makes its argmin anchor jump *more* often.
Setting B is *worse* than the old numerics under today's reader (slack cost SD 22 %) purely because at `1e-3` the
median step drops below 1 ms, the resample is skipped, and the 5 ms windows become "however many solver steps".

### 1b. A seam artefact: `t0_crossing` cycle 1 has been a constant

`runSlackExperiment`'s parallel path (a pool is open — every campaign) merges the init hold's trace into the
slack trace (`outs{1} = out_init`); the serial path does not. The init runs isometric to `velocitytable(3,1)` —
the END of cycle 1's release — while cycle 1's chunk starts at the release START, so the merged trace holds the
init's last, still-isometric sample (LSE > 0) at exactly t = vs(2), inside cycle 1's `t0_crossing` window. The
rule "first stored sample with LSE ≥ 0" takes it: cycle 1 reads **vs(2) − vs(1) = 2.2 ms whatever the
parameters** (serial path / true crossing: 8.3 ms; data t0: 4.9 ms). Every other feature is untouched — the
trace passed to the extractor goes through `makeMonotonous`, where the later chunk wins — but `t0_crossing`
reads the raw merged `out`. Under `FeatureReadGrid = 'data'` the window is de-overlapped the same way and the
crossing interpolated; verified: parallel and serial paths then give identical features.

## 2. Which constructions are fragile — the rule

- A **level crossing on a steep flank** is well conditioned: error = noise / slope, continuous in the trace.
- An **extremum value** is well conditioned (first-order blind to where the extremum is) — but on noisy data the
  max/min of *raw samples* is biased (vall_y **−2.6 %**; a parabola vertex removes it, vall_yq **0.0 %**).
- An **extremum location** is not: on a flat top/bottom its error is √(noise/curvature), and read off samples it is
  quantized to one sample — **1 ms on the model's 1 kHz read, 0.25 ms on the data**, against a 3 ms peak1_t.
- A time measured from **the first sample of a window** inherits the read grid's quantization; measure from the
  protocol time (velocity-table row) instead.

rsT20 violates both of the last two: it is timed from the **argmin of the post-restretch valley**. peak1_t,
peak1_dSL, vall_t read argmax/argmin samples; restretchSlopeStart fits "≥ 3 samples", i.e. 1 ms on the data and
2 ms on the model.

## 3. The data side ([RunDataNoise.m](RunDataNoise.m), [RunDataNoiseRobust.m](RunDataNoiseRobust.m))

Surrogate recordings = the trace's smooth part (Savitzky–Golay 3/9, which passes 25 % of the noise variance —
compensated) + fresh noise drawn in 10 ms blocks from the residuals of the flat end of each post-restretch hold.
The noise is white (lag 1–4 autocorrelation |r| ≤ 0.05), σ 1.28–1.56 kPa. `s_noise` = SD over 40 surrogates,
relative to the pooled 8 mM mean; `s_betw` = buildPooledTarget's between-preparation relSD; `b_noise` = the shift
noise puts into the data value; `z` = 2 mM/8 mM discriminability as in buildPooledTarget.

| current | s_betw | s_noise | noise/betw | b_noise | z_ATP | → replacement | s_betw | s_noise | b_noise | z_ATP |
|---|---|---|---|---|---|---|---|---|---|---|
| rsT20 | 23.4 % | 21.3 % | 0.91 | −11.0 % | 2.27 | **rsTfitE** | 9.9 % | 7.0 % | −2.5 % | 2.76 |
| | | | | | | rsT20e | 11.3 % | 9.8 % | −4.4 % | 3.20 |
| rsK80 | 7.7 % | 8.2 % | 1.06 | −1.3 % | 6.13 | rsKls | 9.9 % | 4.8 % | −0.6 % | 4.31 |
| t0 (fitted onset, t0_crossing's target) | 37.7 % | 2.2 % | 0.06 | | 1.05 | tOn20 | 19.5 % | 2.8 % | −0.1 % | 2.67 |
| peak1_dSL | 17.9 % | 13.6 % | 0.76 | −2.2 % | 0.85 | peak1_dSLq | 20.5 % | 10.0 % | −2.3 % | 0.69 |
| peak1_t | 18.2 % | 13.3 % | 0.73 | −2.5 % | 0.76 | peak1_tq | 20.8 % | 9.6 % | −2.3 % | 0.55 |
| vall_t | 19.1 % | 18.8 % | 0.99 | +2.3 % | 0.75 | tFall50 | 14.4 % | 7.7 % | +1.4 % | 0.19 |
| restretchSlopeStart | 30.0 % | 19.7 % | 0.66 | −0.5 % | 0.67 | rsSlope1ms | 31.2 % | 13.5 % | −0.5 % | 0.64 |
| vall_y | 24.5 % | 1.1 % | 0.05 | −2.6 % | 1.22 | vall_yq | 24.1 % | 1.0 % | 0.0 % | 1.15 |
| ktrfit | 12.8 % | 1.6 % | 0.12 | −0.4 % | 2.60 | (keep) | | | | |
| steady / peak1_y | 23.6 / 23.9 % | 0.9 / 1.6 % | ≤ 0.07 | | 1.15 / 1.32 | (keep) | | | | |

Read-outs:
- **rsT20 is noise in the data as well**: its between-heart spread *is* its measurement noise, and noise shortens
  it by 11 %. Anchoring the same crossing on the **ramp end** (rsT20e) or reading the delay off the band fit
  (rsTfitE) halves the noise, quarters the bias and *raises* the ATP discriminability.
- **rsK80's spread is also pure noise, but its ATP effect (−62 %) is so large that z stays 6.** A least-squares
  rate on every raw sample of the rise (rsKls) has **40 % less** measurement noise; numerically rsK80 is already
  fine once the reader is fixed (0.05 %), so this is a data-quality option, not a noise fix.
- **Valley and peak timings carry no ATP information under any definition** (z 0.19–0.85). vall_t is pure noise.
- **The level-crossing onset tOn20 beats the fitted t0** on the only thing that matters here (z 2.67 vs 1.05);
  on the model it is noise-free under every setting (≤ 0.01 %) and needs no model-internal LSE read.
- The ramp-read rsK80 phase test on the data (decimated to 1 kHz at its 4 phases: rsK80 moves **17 %**, rsT20
  **40 %**) says the same thing from the other side: none of these would survive being *recorded* at 1 kHz.

## 4. Strategies, scored on the objective itself ([ScoreStrategies.m](ScoreStrategies.m))

The slack terms of refit9's `fn`, evaluated on every saved trace. Noise = relative SD of the slack cost over the
9 points; the cost level is comparable **within a column only** (the objective differs between columns).

| setting | S0 today | S1 reader fix | S1 + rsT20→rsT20e | S1 + rsT20→rsTfitE | S1 + full replacement set¹ |
|---|---|---|---|---|---|
| A old: `1e-2`, step output, hard rectifier | 0.177 ± **6.84 %** | 0.183 ± 2.89 % | 0.170 ± 1.69 % | 0.195 ± 1.51 % | ± 0.78 % |
| B SiLU `1e-3`, step output | 0.204 ± **21.97 %** | 0.194 ± 1.47 % | 0.176 ± 0.47 % | 0.201 ± **0.36 %** | ± **0.11 %** |
| C SiLU `1e-3`, dense in cycles | 0.285 ± 0.36 % | 0.194 ± 0.93 % | 0.180 ± 0.47 % | 0.205 ± **0.37 %** | ± **0.11 %** |
| D SiLU `1e-2`, dense in cycles | 0.283 ± 1.41 % | 0.191 ± 0.80 % | 0.178 ± 0.41 % | 0.203 ± 0.32 % | ± 0.20 % |
| E reference SiLU `1e-4`, dense in cycles | 0.287 ± 0.04 % | 0.195 ± 0.06 % | 0.181 ± 0.05 % | 0.206 ± 0.05 % | ± 0.02 % |

¹ rsK80→rsKls, rsT20→rsT20e, peak1_dSL→peak1_dSLq, vall_t→tFall50, restretchSlopeStart→rsSlope1ms; targets = the
3-preparation 8 mM means, floors `nEl·s²/3` from the data spread (no numerical-noise term), weights of the replaced
term. **It changes what the campaign fits** (level 0.92, see §5), so it is shown for its noise only.

- The **reader fix alone** removes the ktrfit_shape noise (SD 0.0415 → 0.0000 at B) and most of rsK80's; what is left
  is rsT20's argmin. Replacing that one definition brings every cheap setting to **0.3–0.5 %**, i.e. the level the
  previous fix bought with a 2.4–4.5× slower evaluation.
- **Accuracy** (the reference's own cost, rsTfitE column): B −2.3 %, D −1.2 %, C −0.5 %. The old numerics A read
  **−38 %** under today's objective but only **−6 %** once read properly (S1: 0.183 vs 0.195) — so even most of
  *that* error was the reader; the rest is the hard rectifier and unclamped rates.

## 5. What the robust features reveal (model E vs pooled data, cycles 1–4)

vall_t and peak1_dSL are refit9 **guards** by design ("force shape as guards"), their deadbands set at 1.2× the
seed's own cost (0.72 and 0.52) — so the misfit below was parked, not scored. Measured robustly it is real:
- **The post-peak force drop is ~2× too slow** in the model: tFall50 8.6 → 6.7 ms vs data 5.0 → 3.6 ms
  (model/data 1.72–2.00). vall_t showed the same (model 20–65 % late) but paid nothing.
- **The recovery delay is ~25 % short**: rsTfitE model/data 0.75 (rsT20e 0.78–0.83).
- **The re-stretch rate is ~5–20 % fast** by the least-squares read (rsKls 1.03–1.22); the two-crossing rsK80
  straddles the data (0.94–1.15).
- **The onset ladder already fits**: tOn20 model/data 0.98–1.15.
- peak1_dSLq model/data 1.17–1.35 (the overshoot peaks late), but z_ATP 0.69 — shape, not ATP.

## 6. Cost of each setting ([RunCostCount.m](RunCostCount.m), [RunInterleavedTiming.m](RunInterleavedTiming.m))

Wall time on this machine today was not a measurement (up to three other sessions' 5-worker pools on 10 cores;
the same config measured 33 s and 62 s in `StrainWindowNoise/results/costtiming*.log`). So cost is counted in RHS
calls (contention-proof), and wall-time **ratios** come from configurations run in rotation.

| setting | critical-path RHS calls (init + slowest cycle) | full objective RHS calls (solver + output) | serial slack run, wall, interleaved ×3 |
|---|---|---|---|
| A old: `1e-2`, step output, hard rectifier | 28.0 k (1.00) | 175.8 k (1.00) | 64.5 s (1.00) |
| B SiLU `1e-3`, step output | 29.4 k (1.05) | 213.4 k (1.21) | 86.7 s (**1.33**) |
| C SiLU `1e-3`, dense in cycles only | 28.6 k (1.02) | 218.7 k (1.24) | 75.9 s (**1.16**) |
| D SiLU `1e-2`, dense in cycles only | 23.4 k (0.83) | 186.7 k (1.06) | — |
| E SiLU `1e-4`, dense in cycles only | 36.1 k (1.29) | 280.3 k (1.59) | — |
| F SiLU `1e-3`, dense everywhere (the campaign's setting) | 38.1 k (1.36) | 228.2 k (1.30) | 157.8 s (**2.45**) |

- **The tolerance is cheap.** `1e-3` costs 1.02–1.05× the old critical path; the smooth rectifier even makes `1e-2`
  cheaper than the old numerics (0.83×), because the hard corner at isometric rest was costing steps.
- **The dense grid on the long holds is what is expensive**, and mostly *not* in RHS calls: F is 1.30–1.36× in
  calls but 2.45× in wall. The 94 s init hold at 0.5 ms is 189 k interpolated points, a 189 k × 128 matrix, an
  `interp1` down to 1e4 points and a row-by-row `out.PU` growth in `storeOutputs`
  ([evaluateModel.m:453-468](../../Model/evaluateModel.m#L453-L468)) — 17 s serial on the campaign's critical path,
  where a 5-worker pool cannot hide it. `odeOutputMaxSegment 5` removes it: C is **1.16×**.
- B (no dense grid at all) is *slower* in wall than C (1.33× vs 1.16×) at equal RHS work — most likely today's
  reader running the fits on raw, uneven solver steps; with the reader fix both read the same 4 kHz grid (not
  measured).

## Recommendation

1. **Fix the reader (free, removes most of the noise).** Read the model trace on the *data's own* 4 kHz grid,
   anchored on absolute time — interpolate, never "resample if the median step says so" — and measure every
   feature time from the protocol row, not from the window's first sample. Interpolate `t0_crossing`'s LSE zero
   crossing (or switch it to tOn20 on both sides). Same for `extractKtrAttributes`, which reads raw output with
   sample-count windows (not measured here; 0.16 % at the old numerics per KtrSegmentSlope).
2. **Replace rsT20 by rsTfitE** (or rsT20e): delay from the **ramp end**, read off the 20–80 % band fit.
   Numerical noise 0.06 %, data noise 21 % → 7 %, noise bias −11 % → −2.5 %, ATP z 2.27 → 2.76.
3. **Numerics: SiLU + clamp + `RelTol 1e-3`, dense grid only inside the cycles (`odeOutputMaxSegment 5`)** —
   1.16× the old run time, −0.5 % from converged. No dense grid at all is as quiet (0.36 %) but −2.3 % off.
   Never the 0.5 ms grid on the 94 s init hold (2.45×).
4. **Floors: nothing to rebuild.** refit9 already sets every ATP term to `nEl·s²/3` from the data spread alone
   (`BuildRefit9Seed.m`); the wide floors on vall_t / peak1_dSL are its guards (1.2× the seed's cost), a design
   choice. (An earlier version of this note said otherwise — corrected.) The new rsTfitE term gets the same
   construction.
5. **Decide separately (these change the fit, not its noise):** drop vall_t and peak1_t (pure noise, no ATP);
   score tFall50 if the restretch transient's shape matters (the model is 2× slow there); consider rsKls for
   rsK80 (40 % less data noise). restretchSlopeStart needs no replacement once read at the data's rate — it is
   inside its deadband.
6. **Skip the unscored oscillation fit** (`ktr2_*`, [extractSlackAttributes.m:203-248](../../Model/extractSlackAttributes.m#L203-L248)):
   2.6 s of the 9.9 s extraction per evaluation on a model trace (10.9 s of 18.9 s on a data trace), serial on the
   client after the parallel cycles. Only Workbench/analysis scripts read it.

## Implemented

All opt-in; with the switches off every existing feature is bit-identical (55/55 fields on the three 8 mM
preparations, on a serial and on a parallel-path model trace; the ktr features 3/3). Installed by staging outside
the repo, testing, then an atomic rename, because two campaigns (`OverstretchedHeads/RunYieldRefit`, yieldA/B)
were running on this worktree's `Model/`.

| where | what |
|---|---|
| `params.FeatureReadGrid = 'data'` → [extractSlackAttributes.m](../../Model/extractSlackAttributes.m) | the model trace is read on an absolute 0.25 ms grid (the data's rate), later chunk winning where chunks overlap, instead of the median-step rule; `t0_crossing` de-overlapped and interpolated (§1b) |
| same switch → [runKtrExperiment.m](../../Model/experiments/runKtrExperiment.m) | the ktr trace is read at the ktr data's own time stamps before `extractKtrAttributes` (ktrP_k80 46.249 → 46.247 at refit9) |
| `params.SkipOscillationFit = 1` → extractSlackAttributes | skips the unscored `ktr2_*` fit (fields kept, NaN). Extraction 3.9 s → 0.7 s with both switches |
| always computed, scored by nothing unless `fn` says so | `rsT20e`, `rsTfitE` (delay from the ramp end), `tOn20` (onset as a 20 % level crossing) |
| [Model/pooledTarget8mM.m](../../Model/pooledTarget8mM.m) | regenerated with the three new features ([RegeneratePooledTarget.m](RegeneratePooledTarget.m)): re-extraction bit-identical for all 9 runs; every existing value row and `Fref` bit-identical; the moderated `s` of existing features moved −1.4 … −0.6 % (the shrinkage prior now sees 27 features) — `s` is not read at run time. rsTfitE target 15.5–17.7 ms, raw spread 9.6 %, ATP z 2.40; `pool_raw.mat` backed up as `pool_raw_20260924_pre_robust.mat.bak` |
| [BuildRobustSeed.m](BuildRobustSeed.m) → `params/robust1_seed.m` | refit9_opt with numerics C + both switches, `rsT20[1:4]` → `rsTfitE[1:4]` (fl = 4 s²/3, w = λ/s², s from the 3 preparations), guards re-derived at 1.2× this point's raw cost under the new reader. **Cost 0.4240** at refit9_opt's parameters; **4 neighbours 1e-6 away: SD 0.17 % of the full objective** (the previous fix: 0.3 % at 2.45× the run time); snapshot reload bit-identical. rsTfitE is **64 %** of it (0.270): the model's delay is ~25 % short on a feature with a 9.6 % spread, so the rule gives it w = 1.15 (rsT20 had 0.22) — the next campaign will chase the re-stretch recovery delay first. Then ktrfit_shape 0.064, rsK80 0.044, steady 0.031, t0_crossing 0.014 (now the real crossing); guards 0. |

`getParams` has no defaults for the two switches (as for `odeOutputDt` / `odeOutputMaxSegment`): absent = off.

## Caveats

- n = 3 preparations: every `s_betw` has 2 dof. The *noise* numbers are well determined (40 surrogates × 3).
- One parameter point (refit9_opt). The reader/definition findings are structural; the cost levels are not.
- The surrogate noise is white by measurement, but drawn from plateau residuals only; ramp-edge artefacts are
  not in it.
- Only the slack block and `ktrfit_shape` were scored; FV, stairs, passive and the ktr protocol were not re-read.

## Files

[robustSlackFeatures.m](robustSlackFeatures.m) — current definitions (bit-identical to extractSlackAttributes) +
replacements · [ProbeDenseSegments.m](ProbeDenseSegments.m) — identical state, different phase ·
[RunDataNoise.m](RunDataNoise.m) / [RunDataNoiseRobust.m](RunDataNoiseRobust.m) — data noise ·
[RunModelNoise.m](RunModelNoise.m) → [AnalyzeModelNoise.m](AnalyzeModelNoise.m) — numerical noise per definition
and read · [ScoreStrategies.m](ScoreStrategies.m) — the objective under each strategy ·
[RunCostCount.m](RunCostCount.m) + [dPUdT_CountedTyped.m](dPUdT_CountedTyped.m) — RHS calls ·
[RunInterleavedTiming.m](RunInterleavedTiming.m) — wall-time ratios · [TestStagedExtractor.m](TestStagedExtractor.m) —
bit-identity with the switches off, parallel vs serial path, ktr · [RegeneratePooledTarget.m](RegeneratePooledTarget.m) ·
[BuildRobustSeed.m](BuildRobustSeed.m) · `results/` — every `.mat`, `.csv`, `.log`.
