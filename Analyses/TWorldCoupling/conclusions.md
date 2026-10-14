# T-World as the Ca source of the cross-bridge model

**Status:** implemented and measured, 2026-10-14. Branch `TWorldCoupling` (worktree
`C:\home\git\ATP-wt-tworld`), built on `HumanTwitch` = PooledTarget-8mM + the Ca layer
(`CaModel`/`CaInput`, `calciumActivation.m`). Off by default: `CaInput = 'tworld'` is one more
Ca source of the existing layer, `CaModel = 'none'` stays byte-identical (checked on
`ExampleHighATP` and `human_twitch_t1_isoT`: same `PU0`, max |ΔRHS| = 0 on four states;
`TestCaActivation` T2–T7 pass).

Reproduce: `run('Analyses/TWorldCoupling/RunTWorldCoupling.m')` (~2.5 min; cases cached in
`results/case_*.mat`). Table: [`results/tworld_coupling.csv`](results/tworld_coupling.csv).
Figures: [`results/tworld_coupling_human.png`](results/tworld_coupling_human.png),
[`results/tworld_coupling_marzban.png`](results/tworld_coupling_marzban.png).

---

## 1. What T-World is, and what "force feedback" it has

T-World (Tomek, Oxford; `ThirdParty/TWorld/README.md` for provenance) is a human ventricular
myocyte: 93 states (membrane voltage, Na/K/Ca/Cl handling, buffers, SR, RyR, CaMKII), with the
**Land 2017** contraction model inside (states 69–74: `XS, XW, Ca_TRPN, TmBlocked, ZETAS, ZETAW`).
Time in ms, concentrations in mM. At 1 Hz from its stored endocardial steady state:
cytosolic Ca **0.097 → 0.42 µM**, TTP **54 ms**, Land tension 19.6 kPa peaking at 168 ms;
RHS 78 µs, ~340 solver steps per beat, 0.8 s per beat standalone.

The only feedback from contraction to Ca is **troponin Ca buffering**:
`d(buffers_TnClow)/dt = Bmax_TnClow · dCa_TRPN/dt` with `Bmax_TnClow = 70 µM`, subtracted from
the cytosolic Ca balance. Land's tension `Ta` does not feed back on anything, and the length
`lambda = 1` is hard-coded (no length feedback). So "force feedback" = the Ca that the thin
filament holds.

## 2. Coupling options and what was implemented

All three are the same switch, `params.CaTworld.coupling`, resolved once by
`Model/tworldCaLayer.m` and read by `calciumActivation.m` (`params.ca.input = 4`):

| coupling | what runs where | extra states | Ca buffering by troponin | loop |
|---|---|---|---|---|
| `'signal'` | T-World runs standalone for one beat in `getParams` (`modelRunner`, 0.8 s); its Ca is a periodic 0.25 ms table (`sigLookup`, Hermite) | 0 | Land's (T-World's own) | open |
| `'land'` | T-World's 93 states are appended to `PU` after the thin-filament state and integrated by `ode15s` together with the cross-bridge ODE; `[Ca]_i` (state 8) drives the thin filament | 93 | Land's | open (validation of the embedding) |
| `'xb'` | as `'land'`, but the flux into T-World's TnC buffer is **this model's thin filament**: `Bmax_TnClow · dA/dt`, `A` = permissive fraction (`1 − N` for `twoState`, `N_on` for `myosim`). Land's `Ca_TRPN` keeps evolving (it no longer buffers; `Ta` stays available as a comparator) | 93 | the XB thin filament, including its bound-head feedback (`Pu = 1 − N − bound`) | closed |

Mechanics of the embedding: beats start at `CaFirstBeat + k·CaPeriod`; T-World gets `t` in
ms since the beat onset, so its 1 ms stimulus fires every beat; `evaluateModel` restarts the
integrator at stimulus on/off (`ca.beatKnots`). Before `CaFirstBeat` the block is **frozen** at
its pre-paced state, so beat 1 starts exactly from the standalone steady state. The vendored
files are untouched: the `'xb'` correction is applied to the returned derivative vector
(`ydot(8) += ydot(13) − J_xb; ydot(13) = J_xb`, verified `ydot(13) ≡ Bmax·ydot(71)`).
Per-state `AbsTol` (1e-6 for the T-World block, T-World's own default) is passed by
`evaluateModel`; the cross-bridge states keep theirs. The C MEX path is refused with any Ca
model, as before. Pre-pacing at another period: `CaTworld.prebeats`; T-World parameter
overrides (cell type, current multipliers, sex, stimulus): `CaTworld.param`.

## 3. Compatibility verdict

**Compatible, with two numerical caveats that are handled.**

- The Ca layer's design (a Ca source feeding a thin-filament gate on `ka_eff`, states appended
  last, beat knots) took the 93-state block without changes to the RHS.
- Units: ms/mM vs s/µM, converted in the wrapper. Scales: `[Ca]_i ≈ 1e-4 mM` is below the
  cross-bridge `AbsTol` 1e-4, hence the per-state `AbsTol`. `RelTol`: the harness uses 1e-3
  (T-World's default); `evaluateModel`'s 1e-2 default is loose for an action potential.
- Stiffness: the AP upstroke sets the step size for the whole system. Cost on this box
  (one worker, `human_twitch_t1_isoT`, ss = 61, 0.3 s rest + 3 beats, isometric):

| case | wall (s) | points |
|---|---|---|
| `signal`, human / marzban thin filament | 13.1 / 7.8 | 6616 |
| `land` | 33.5 / 17.3 | 6729 |
| `xb` | 35.1 / 22.9 | 6730 |

  i.e. ~2.5× the cost of the prescribed-transient run for the embedded block. Fine for twitch
  studies; a fit loop should use `'signal'` unless the feedback is the question.
- Correctness of the embedding: beat 1 of `'land'` reproduces the standalone T-World beat to
  **max |ΔCa| = 5e-4 µM** (different integrator/tolerances); `'signal'` matches `'land'` to
  2.4e-3 µM (table interpolation).

## 4. Both feedback cases, measured

Two thin-filament sets: **human** = the EHT-twitch fit of `human_twitch_t1_isoT`
(`tf_k_on` 9.22 /(µM s), `tf_k_off` 5.08 /s, `K_coop` 20) and **marzban** = the Marzban 2020
rat defaults (101, 724, 9.68). Last of 3 beats:

| thin filament | coupling | Ca peak (µM) | Ca TTP (ms) | F_dev peak (kPa) | F TTP (ms) | RT50 (ms) | gate max |
|---|---|---|---|---|---|---|---|
| human | signal / land | 0.425 | 53.5 | 10.1 / 10.0 | 145.5 | 105.5 | 0.145 |
| human | **xb** | **0.581** | **34** | **12.3** | 114.5 | 118 | 0.185 |
| marzban | signal / land | 0.425 | 53.5 | 0.32 | 91.5 | 144 | 0.006 |
| marzban | **xb** | **1.123** | 49.5 | 0.98 | 73 | 63.5 | 0.018 |

**Open loop (`signal`/`land`):** T-World's 0.42 µM transient opens only 15 % of the EHT-fit's
thin filament (gate max 0.145) and gives a **10 kPa** twitch, against 46 kPa with the de Lange
0.15 → 0.60 µM `risedecay` transient the fit was made with. The fitted thin filament is tuned to
that transient's level and 93 ms TTP; T-World is lower and faster (TTP 54 ms). Force TTP 146 ms
and RT50 106 ms are in the EHT range regardless.

**Closed loop (`xb`), human set:** the XB thin filament holds about the same Ca as Land's
troponin (25 vs 27 µM peak bound) but binds it **more slowly** (k_on 9.2 /(µM s) against Land's
effective k_off 78.5 /s with n_H 1.65 at 0.76 µM), so less Ca is buffered during the upstroke:
the transient peaks **earlier and 37 % higher** (0.58 µM at 34 ms), decays faster afterwards,
and the twitch grows 23 %. Permissive fraction 0.22 vs 0.18 open loop.

**Closed loop, marzban set:** the fast/weak thin filament stays <2 % permissive, so the 70 µM
TnC buffer effectively disappears and the Ca transient **more than doubles** (1.12 µM).

The last row is the thing to decide, not a bug: **in Land's model `Ca_TRPN` is at once the Ca
buffer and the activator; in this model `N` is a regulatory switch whose permissive fraction
need not equal the Ca-bound TnC fraction.** The `'xb'` mapping takes them as equal. If that is
not wanted, the alternatives are (a) keep Land's buffer and let only the gate come from the XB
thin filament (that is exactly `'land'`), or (b) add a separate stoichiometric TnC-Ca state to
the Ca layer and drive `N` from it. Which one fits the question is the user's call; both
existing couplings stay available.

## 5. What a user does

```matlab
P = loadParams('human_twitch_t1_isoT');
P.CaModel = 'twoState';  P.CaInput = 'tworld';
P.CaPeriod = 1;  P.CaFirstBeat = 0.3;
P.CaTworld = struct('coupling', 'xb');        % 'signal' | 'land' | 'xb'; prebeats, param, state, X0, absTol
P.odeRelTol = 1e-3;  P.Velocity = 0;
P = getParams(P, [], true, false);  P.Vums = 0;  P.v = 0;
[~, out] = evaluateModel(@dPUdT_CombinedTransitions, [0 3.3], P);
% out.Ca (uM), out.CaGate, out.CaState(:, 1) = N, out.CaState(:, 2:end) = the 93 T-World states
```

Nothing of `CaTworld` is written into snapshots (it is a struct); `CaInput = 'tworld'` in a
snapshot without a `CaTworld` field resolves to `'land'` from the stored 1 Hz endocardial state.

## 6. Fitting the EHT twitch with T-World's Ca (2026-10-14)

Tiers `tw1_isoT` (adult beat as granted, coupling `'signal'`) and `tw2_caH` (SERCA, ICaL, NCX
multipliers freed, pre-paced) in `Analyses/HumanTwitch/`; comparison in `CompareTwitchFits.m`,
[`results/twitch_fit_tw1_tw2.png`](results/twitch_fit_tw1_tw2.png),
[`results/twitch_fit_tw1_tw2.csv`](results/twitch_fit_tw1_tw2.csv). rms 0.97 % as granted,
0.75 % with the Ca handling freed (de Lange transient: 0.77 %); the identified net cycle factor
is 0.22–0.24 against 0.45 with de Lange's Ca. The closed-loop `tw1_isoT_xb` stalled at 6.7 %
for the §4 reason. Full account: `Docs/reference/tworld-ca-source.md` §7.
