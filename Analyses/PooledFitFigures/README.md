# Pooled-fit figures (EMC-poster layout, current fit)

The EMC poster (`Docs/presentations/2026-EMC.pptx`) showed the `spreadw` fit against **one**
preparation (03/27), with the others only as an SD band. The model is now fitted to the
**3-preparation pooled target**. This folder redraws the same panels with:

- **data line = pooled mean, band = between-preparation SD (n = 3)**;
- **8 mM model = `robust1_seed`** (= `refit9_opt` + robust feature reader, the incumbent);
- **2 mM model = fit A**: the same core with the k2 power gate, R2 ×0.54
  (`Analyses/LowATP_Refit9K2Sweep`).

| file (in `figures/`) | panel |
|---|---|
| `slack_all.png` (2 cells) | slack–restretch, all five cycles, with ML trace and ladder labels |
| `slack_cycle4.png` (2 cells) | cycle 4 detail |
| `ktr_protocol.png` (2 cells) | ktr protocol (new since the poster), F/F_iso |
| `fv_force.png`, `fv_power.png` | F–V and power; fine model sweep; pooled in-prep FV2 point (n = 2) |
| `feat_<name>.png` ×12 | slack features over the five cycles (poster panels a–f + 6 more) |
| `stairs.png` (3 cells) | staircase **prediction**, never fitted. One tile per preparation (03/27 has 8×0.1 s steps, 04/03 and 04/10 have 4×0.2 s, so they are not poolable). 8 mM data only; the 2 mM curve has no data |
| `passive.png` (1.5 cells) | Ca-sensitive passive force, 03/27 (the only PNB+Mava set) |

## Caveats on the bands

- **Force features are normalised per preparation** (Fss22 anchor + rundown lift) before
  pooling. Their band is the **shape** spread. The between-heart **level** spread (CV 15 %) is not in it,
  which is why the 8 mM steady-force band is narrow.
- **F–V** has no pooled 3-prep curve. The band is the poster's: Baker2022 vs AB2021, two
  preparations. The pooled in-prep data contribute only the FV2 diamond at 2 ML/s.
- **Passive** is one preparation and is shown as such.

## Re-run

```matlab
run('Analyses/PooledFitFigures/RunPooledFigureSims.m')   % ~1 min: fine F-V + passive -> results/sims.mat
run('Analyses/PooledFitFigures/MakePooledFitFigures.m')  % seconds: all PNGs -> figures/
```

Slack, ktr and feature model values are read from the saved fit-A runs in
`Analyses/LowATP_Refit9K2Sweep/results`. Re-run `RunBestFitFigures.m` there after a refit.
Gotcha: `runFVExperiment` takes [ATP] from the **first element of its 2nd argument**, not from
`params.MgATP`.
