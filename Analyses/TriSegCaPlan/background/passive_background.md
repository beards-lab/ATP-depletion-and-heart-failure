# Passive stress, relaxed vs Ca-stiffened (PNB+Mava) skinned mouse trabeculae

Script: `../passiveBackground.py` (figure `../figures/passive_ca.png`, raw tables `tables.txt`).

## Files
`data/PassiveCaSrc2/<date>/NN_*.txt`, 4 header lines, columns **Time (s), L (Lo), F (kPa), SL (um)**, 2 kHz (0.1 s files) or 20 Hz (100 s files). F is already stress (kPa) but carries an arbitrary offset (about 0-20 kPa; ~20 kPa on 20241121 PNB files), so a **baseline = median F in the slack hold at 0.8 Lo** (after the ramp, passive about 0) is subtracted per file.
PNB+Mava preps (7, mice, male/female mixed): 20240705, 20241010, 20241121, 20241212, 20241217, 20241219, 20241220 (older 2023 dates are PNB only, not used; `04_*` files are 10 C, not used; 20240705 pCa 5.5-6.25 are "NoZDR").
Protocol per file: hold 0.95 Lo, then ramp to 1.1757 Lo (100 s, 10 s, 1 s, or 0.1 s step), hold ~60 s, release to 0.8 Lo (baseline), return. `01_*` relaxed (no drugs), `03_02-05` relaxed PNB+Mava, `03_07-11` pCa 4.4, `03_12-15` pCa 5.5, 5.75, 6, 6.25 (0.1 s step only, 60 s hold).
**Range covered: L = 0.95 to 1.1757 Lo, SL about 1.92-1.97 um to 2.25-2.31 um** (per prep min/max in tables.txt; SL max 2.25-2.31). **SL 2.3 is reached only by one prep (20241219); otherwise extrapolated, shown in brackets.**
The laser SL channel is unusable in places (all zeros on 20240705, 1.7 um in pCa 4.4 on 20241121, saturating noise on fast ramps), so SL was **mapped from L**: per prep a quadratic SL(L) fitted to the PNB-relaxed 100 s ramp (20240705 uses the pooled fit). So the SL axis is a per-prep mapping, not a per-file measurement.

## Passive stress (kPa), 100 s ramp (quasi-static-ish), mean of 7 preps; per-prep values in tables.txt
Measured = inside the SL range reached; (x) = power-law extrapolation beyond SLmax.

| condition | SL 2.0 | 2.1 | 2.2 | 2.3 |
|---|---|---|---|---|
| relaxed, no drug | 0.78 | 1.71 | 3.58 | (7.4) |
| relaxed, PNB+Mava | 0.58 | 1.43 | 3.16 | (7.0) |
| pCa 4.4, PNB+Mava | 0.67 | 1.56 | 3.41 | (7.5) |
| ratio pCa 4.4 / relaxed PNB (mean of per-prep ratios) | 1.7 (noisy, values ~0.5 kPa) | 1.20 | 1.15 | 1.09 |
| **model law** k=11297, Lsc0=1.80093, gamma=7.689 | 0.05 | 1.05 | 9.67 | 54 |

Between preps the 2.2 um value spans 1.2-5.8 kPa (relaxed) and 1.2-5.1 (pCa 4.4); the between-prep spread (factor ~4) is far larger than the Ca effect.

## Intermediate pCa (0.1 s step to 1.1757 Lo, stress 60 s into the hold, top SL about 2.25-2.31, mean 2.27)
Only this single SL is available for pCa 5.5-6.25 (no ramp data). Values per prep (kPa), order 20240705, 1010, 1121, 1212, 1217, 1219, 1220:
- relaxed PNB+Mava 1.8, 6.0, 3.5, 7.6, 4.2, 8.8, 5.7: mean 5.4
- pCa 6.25 mean 5.3; pCa 6 mean 5.3; pCa 5.75 mean 5.3; pCa 5.5 mean 5.2; pCa 4.4 mean 5.5

Within every prep the change from relaxed to pCa 4.4 is +0 to +0.5 kPa (+0-10 %), with no systematic pCa trend.

## Findings
1. **With PNB+Mava there is essentially no Ca-dependent passive stiffening in these files**: high-Ca / relaxed = 1.1 (SL 2.1-2.3), and nothing between pCa 6.25 and 4.4. The "3 kPa low-Ca vs 10 kPa high-Ca" in `Analyses/PassiveForceID/conclusions.md` is not supported by this dataset; it presumably comes from the 2026 slack data where the drug state is different (verify there; that high Ca value may include residual cross-bridge force or a different protocol). A sigma_pas(SL,Ca) law from these data would be sigma_relaxed(SL) x (1.0 to 1.1).
2. The relaxed (drug-free) and relaxed PNB+Mava curves agree within 10 %, so Mava/PNB do not change passive stress; "Relax" is truly relaxed (no active force at pCa 9 anyway).
3. Absolute level at true SL 2.2 um is 3-3.6 kPa (mean), 1-6 per prep; ~5 kPa at SL 2.27 for the 60 s hold. The model law gives 9.7 kPa at 2.2 and 54 kPa at 2.3, i.e. 3x too stiff at 2.2 and far too steep above 2.25 against these data. Caveat: the model law was fitted to the 2026 file whose SL is the nominal 2*ML = 2.2 um at ML 1.10 (that file's ML range is 0.94-1.10, SL nominal 1.88-2.20, forces up to ~5 kPa steady at 2.2). The 2024 files' true SL (2.2 at about L = 1.08 Lo) is not on the same scale, so the SL axes are not directly comparable; the comparison needs a common length definition (Lo, not SL).
4. Viscous part: the 0.1 s step peaks far above the hold (e.g. 42 kPa transient at pCa 4.4, 20241121, vs 4 kPa steady); only the end of the hold and the 100 s ramp are used here. Even 100 s ramps contain some rate-dependent stress, so the table is an upper bound of the elastic part. Ramps are the first stretch of the protocol after S2F ordering; rundown not assessed.
