# Passive stress transient: relaxed against pCa 4.4 (PNB + Mava), by ramp duration

Script: `../passiveDynamic.py`. Ramp from 0.95 to 1.1757 Lo (top SL about 2.27 um), 7 preparations.
Stress is the rise above the 0.95 Lo hold. Peak = highest stress around the end of the ramp (5 ms smoothing);
hold = stress 50 to 60 s into the hold.

| ramp duration (s) | n | peak relaxed (kPa) | peak pCa 4.4 (kPa) | pCa 4.4 / relaxed, mean of per-prep ratios (range) | peak / hold, relaxed | peak / hold, pCa 4.4 |
|---|---|---|---|---|---|---|
| 100 | 7 | 6.6 ± 2.8 | 7.0 ± 3.0 | 1.07 (1.00 to 1.14) | 1.21 | 1.25 |
| 10 | 7 | 7.5 ± 3.2 | 11.0 ± 4.6 | 1.48 (1.32 to 1.63) | 1.46 | 2.06 |
| 1 | 7 | 8.7 ± 3.7 | 25.1 ± 10.4 | 2.94 (2.29 to 4.12) | 1.71 | 4.83 |
| 0.1 | 7 | 10.5 ± 4.6 | 41.0 ± 16.8 | 4.03 (2.91 to 5.92) | 2.08 | 8.14 |

## Per preparation, 0.1 s ramp (kPa): peak relaxed, peak pCa 4.4, hold relaxed, hold pCa 4.4

- 20240705: 3.5, 12.2, 1.5, 1.6
- 20241010: 13.7, 46.0, 5.7, 5.4
- 20241121: 7.9, 37.0, 3.6, 3.7
- 20241212: 16.0, 69.2, 7.4, 7.5
- 20241217: 7.2, 42.6, 4.2, 4.1
- 20241219: 15.0, 43.8, 8.4, 8.8
- 20241220: 10.5, 36.4, 5.2, 5.6
