# ATP-dependent SRX return `kmsrd` — the ktr lever for BOTH magnitude and per-segment slope (2026-07-14)

## Question
`RunSrxDecouple.m` validated that ATP-dependent `kmsrd` (20.87 @8 mM → 8 @2 mM) gives cost 0.135 and the
right ktr **mean** (0.77→0.50). Open: does it also give the DATA's per-segment ktr **SLOPE**
(0.64 0.59 0.53 0.47 0.46), or only the mean? And does its force-sensitivity `sigma_srd1` tune it?

## Result — `kmsrd` ALREADY produces the slope (`RunKmsrdSlope.m`)
Base `params_2state_lowATP` (ADP-trap); `kmsrd`(8 mM)=20.87 (baseline), `kmsrd`(2 mM) lowered = the decoupling.

| config | cost | ktr ratio (per segment) |
|---|---|---|
| ref (no decouple) | 0.473 | 0.76 0.77 0.77 0.80 0.77 (FLAT) |
| kmsrd2=12 | 0.188 | 0.64 0.63 0.61 0.63 0.58 |
| **kmsrd2=8** | **0.135** | **0.55 0.52 0.50 0.51 0.44 (SLOPES DOWN)** |
| kmsrd2=4 | 5.39 | 0.37 0.33 0.31 0.29 0.24 (over-decoupled) |
| kmsrd2=8, sigma_srd1=10 | 0.207 | 0.68 0.65 0.64 0.63 0.61 (flatter) |
| kmsrd2=8, sigma_srd1=80 | 4.55 | breaks other features |
| kmsrd2=8, sigma_srd1=200 | 10.1 | breaks |

**The validated kmsrd2=8 (cost 0.135) already gives a downward per-segment ktr slope (0.55→0.44),
matching the data trend (0.64→0.46).** Earlier reports missed this by quoting only the MEAN.

## Why (corrects the earlier "strain-blind" claim)
`kmsrd` is NOT strain-blind: the SRX-ADP return `RSRD2PD = kmsrd·exp(F_SR/sigma_srd1)·P_SRD` depends on the
restretch force `F_SR`, which differs per slack segment (bigger slack → higher force). So the return rate —
and hence ktr — is naturally **segment-dependent** → the slope. `sigma_srd1` at baseline (24.57) is
near-optimal: lowering it flattens the slope, raising it destabilizes other features.

## Consequence
The per-segment ktr slope is NOT an open problem — the ATP-dependent `kmsrd` handles BOTH the ktr magnitude
(force↔ktr weld-break) AND the slope, at cost 0.135. This is exactly why the forced-detachment parking
state (`LowATP_ParkedState/`, a *head-count* lever) was unnecessary and refuted: `kmsrd` (a *rate* lever
with F_SR force-sensitivity) already does the job. Residual: the model slope (Δ≈0.11) is ~60% of the data
slope (Δ≈0.18) and sits slightly low on seg1–2 — a fine-tune target (kmsrd between 8–12, or a mild F_SR
shape), not a missing mechanism.

## TODO — make it a real one-file mechanism
Currently the ATP-dependence is SIMULATED (set `kmsrd` per condition). Wire `kmsrd_eff = kmsrd·g([ATP])`
into `dPUdT` (flag-gated, `g=1` at 8 mM), bake into `params/params_2state_kmsrd.m`, and confirm 0.135 +
slope emerge from concentrations alone. Reproduce: `Analyses/LowATP_kmsrd/RunKmsrdSlope.m` →
`results/kmsrd_slope.mat` + `results/kmsrd_slope.png`.
