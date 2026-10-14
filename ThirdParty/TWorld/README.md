# T-World (vendored copy)

Human ventricular cardiomyocyte model (membrane voltage, Ca handling, ionic currents,
Land 2017 contraction, beta-adrenergic signalling). Used here as an optional **Ca source**
for the cross-bridge model: `params.CaInput = 'tworld'` (see `Model/tworldCaLayer.m`,
`Model/calciumActivation.m`, and the harness in `Analyses/TWorldCoupling/`).

## Source

| | |
|---|---|
| Author | Jakub Tomek (jakub.tomek.mff@gmail.com), University of Oxford |
| Repository | https://github.com/jtmff/TWorld |
| Commit | `581e711a300357e867319e7e67151ff7c834db12` (branch `main`) |
| Commit date | 2026-08-30 17:37:17 +0100 |
| Copied | 2026-10-14, from the local clone `C:\home\git\TWorld` (the clone's working tree differed from the commit only by CRLF line endings; the files here are `git show <commit>:<file>`, byte-identical to the commit's blobs) |
| License | GPL-3.0 (`LICENSE`, verbatim). Note: this repository is MIT-licensed; the T-World files keep their own license. |
| Online simulator | https://t-world-simulator-multipage-production.up.railway.app/ |

Publications describing the model were not yet listed upstream at the time of copying
(`README_upstream.md`); the model's structure follows the authors' earlier ToR-ORd model.

## Files

Verbatim, unmodified:

- `model_TWorld.m` — the ODE right-hand side, 93 states (plus 57 optional PKA-signalling states, not used here). `output = model_TWorld(t, y, flag_ode, cellType, <18 current multipliers>, nao, cao, ko, ICaL_fraction_junc, INaCa_fraction_junc, PKA_P, Whole_cell_PP1, runSignallingPathway, constantsSig, sex, stimAmp, stimDur, vcParameters, apClamp, extraParams)`; time in **ms**, concentrations in **mM**, voltage in mV.
- `getDefaultParameters.m` — parameter struct defaults (cell type 0 = endo, multipliers 1, `stimAmp -53`, `stimDur 1` ms, `bcl`).
- `getStartingState.m` — 1 Hz steady states: `'TW_endo'`, `'TW_epi'`, `'TW_mid'` (and `_sig` variants with the signalling block).
- `modelRunner.m` — beat-by-beat `ode15s` driver (`[time, X] = modelRunner(X0, options, param, beats, ignoreFirst)`).
- `getCurrentsStructure.m` — post-processing into `currents.time / V / Ca_i / Land_Ta / ...`.
- `LICENSE`, `README_upstream.md` — upstream license and README.

Not copied: `getConstantsPKASignalling.m`, `DataReporter.m`, the sample scripts, the
population `.mat`, the CellML and Modelica exports. Pull them from the repository at the
commit above if the signalling pathway is ever needed.

## State indices the coupling relies on (from `model_TWorld.m`)

| index | state | use |
|---|---|---|
| 8 | `ions_ca_i` (mM) | cytosolic Ca → the thin filament of the cross-bridge model |
| 13 | `buffers_TnClow` (mM) | Ca bound to troponin C; `d/dt = Bmax_TnClow * dCa_TRPN/dt`, `Bmax_TnClow = 70e-3` mM |
| 69–74 | Land 2017 contraction: `XS, XW, Ca_TRPN, TmBlocked, ZETAS, ZETAW` | `Ca_TRPN` is the only contraction state that feeds back on Ca (through index 13); `Ta` (active tension) does not feed back, and `lambda = 1` is hard-coded (no length feedback) |

The `'xb'` coupling replaces the index-13 flux (and its contribution to index 8) by the flux
into the permissive fraction of this project's thin-filament state; see `tworldCaLayer.m`.
No line of the vendored files is changed for that: the correction is applied to the
returned derivative vector.
