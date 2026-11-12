# The P21 host fits (19 animals): data, model, sarcomere kinematics, stresses

Script: `../hostFits.py`; per-animal table: `../figures/host_fits.csv`. Last beat of each saved run.
Pressures kPa (1 kPa = 7.5 mmHg), volumes mL, lengths um, stresses kPa.

| fit | HR | LV EDV data / model | LV ESV data / model | LV EDP data / model | LV ESP data, model Pmax | Ls LV max / min | Ls RV max / min | peak XB / passive / total LV |
|---|---|---|---|---|---|---|---|---|
| Hx1 | 299 | 0.048 / 0.049 | 0.011 / 0.011 | 0.67 / 0.83 | 9.7, 9.7 | 2.25 / 1.87 | 1.76 / 1.44 | 20.8 / 0.82 / 16.5 |
| Hx10 | 299 | 0.044 / 0.045 | 0.013 / 0.013 | 0.68 / 0.62 | 9.6, 10.1 | 2.26 / 1.88 | 1.77 / 1.48 | 23.5 / 0.88 / 19.4 |
| Hx4 | 321 | 0.039 / 0.038 | 0.007 / 0.007 | 0.67 / 0.87 | 11.3, 11.5 | 2.17 / 1.77 | 1.77 / 1.45 | 22.3 / 0.36 / 18.1 |
| Hx5 | 323 | 0.049 / 0.049 | 0.012 / 0.012 | 0.59 / 0.69 | 9.6, 9.3 | 2.28 / 1.88 | 1.76 / 1.48 | 24.4 / 1.07 / 20.3 |
| Hx57 | 356 | 0.062 / 0.064 | 0.030 / 0.029 | 0.43 / 0.70 | 10.7, 10.4 | 2.08 / 1.77 | 1.66 / 1.43 | 25.9 / 0.12 / 23.3 |
| Hx58 | 324 | 0.063 / 0.065 | 0.028 / 0.028 | 0.47 / 0.54 | 8.5, 8.4 | 2.28 / 1.92 | 1.75 / 1.48 | 25.1 / 1.07 / 21.1 |
| Hx59 | 337 | 0.041 / 0.042 | 0.012 / 0.013 | 0.69 / 0.69 | 10.8, 9.8 | 2.24 / 1.94 | 1.77 / 1.52 | 20.8 / 0.71 / 17.4 |
| Hx61 | 388 | 0.032 / 0.031 | 0.006 / 0.008 | 0.53 / 0.75 | 9.5, 9.1 | 2.10 / 1.81 | 1.64 / 1.44 | 20.7 / 0.17 / 17.9 |
| Hx62 | 368 | 0.047 / 0.051 | 0.025 / 0.033 | 0.81 / 0.93 | 10.7, 11.2 | 2.27 / 2.12 | 1.73 / 1.59 | 19.9 / 0.94 / 18.7 |
| Hx7 | 382 | 0.056 / 0.058 | 0.025 / 0.024 | 1.07 / 1.37 | 8.0, 8.9 | 2.33 / 2.03 | 1.73 / 1.48 | 24.1 / 1.66 / 20.7 |
| Hx8 | 353 | 0.059 / 0.063 | 0.031 / 0.028 | 0.63 / 0.71 | 7.5, 7.6 | 2.30 / 1.99 | 1.75 / 1.48 | 21.6 / 1.27 / 17.9 |
| Hx9 | 314 | 0.051 / 0.049 | 0.007 / 0.007 | 0.44 / 0.82 | 6.9, 6.7 | 2.09 / 1.69 | 1.70 / 1.41 | 15.0 / 0.14 / 11.5 |
| Nx11 | 322 | 0.033 / 0.034 | 0.006 / 0.006 | 1.07 / 1.41 | 11.2, 12.6 | 2.28 / 1.99 | 2.19 / 1.95 | 27.2 / 1.11 / 24.1 |
| Nx12 | 292 | 0.064 / 0.060 | 0.013 / 0.013 | 0.27 / 0.42 | 8.0, 7.6 | 2.03 / 1.66 | 2.00 / 1.67 | 19.9 / 0.07 / 17.3 |
| Nx51 | 325 | 0.036 / 0.038 | 0.010 / 0.010 | 0.80 / 0.76 | 9.3, 10.5 | 2.25 / 1.97 | 2.16 / 1.93 | 24.3 / 0.83 / 21.2 |
| Nx52 | 347 | 0.055 / 0.060 | 0.029 / 0.027 | 0.67 / 0.56 | 8.0, 9.2 | 1.91 / 1.71 | 2.28 / 2.01 | 12.5 / 0.01 / 10.4 |
| Nx54 | 384 | 0.054 / 0.055 | 0.029 / 0.031 | 0.44 / 0.45 | 7.1, 7.5 | 2.14 / 1.95 | 2.10 / 1.93 | 20.4 / 0.26 / 18.0 |
| Nx55 | 357 | 0.057 / 0.054 | 0.010 / 0.010 | 0.47 / 1.28 | 8.4, 8.0 | 2.03 / 1.72 | 2.18 / 1.86 | 15.6 / 0.07 / 12.8 |
| Nx56 | 367 | 0.043 / 0.043 | 0.009 / 0.009 | 0.43 / 0.50 | 8.4, 8.0 | 2.05 / 1.78 | 2.19 / 1.96 | 18.0 / 0.09 / 15.2 |

## Ranges over the 19 fits

- HR 292 to 388 bpm; LV EDV (data) 0.032 to 0.064 mL; LV ESP (data) 6.9 to 11.3 kPa; LV EDP (data) 0.27 to 1.07 kPa.
- Sarcomere length: LV max 1.91 to 2.33, LV min 1.66 to 2.12; SEP max 1.65 to 2.28, min 1.51 to 2.13; RV max 1.64 to 2.28, min 1.41 to 2.01 um.
- Shortening per beat (max - min): LV 0.14 to 0.40 um.
- Peak shortening speed: LV 6.20 to 16.91, RV 7.48 to 18.90 lengths/s (of the maximum length), at 0.13 to 0.21 of the beat (LV); median speed while shortening: LV 0.17 to 1.37 lengths/s.
- Peak cross-bridge stress: LV 12.5 to 27.2, SEP 6.1 to 23.8, RV 18.4 to 39.5 kPa.
- Peak passive stress: LV 0.01 to 1.66, SEP 0.00 to 1.11, RV 0.00 to 1.05 kPa.
- Ca transient: diastolic 0.161 to 0.161, peak 2.16 to 2.16 uM.
- Fixed parameters (same in every fit unless a range is shown): Lsref 2.00 to 2.00, Lsc0 1.51 to 1.51 um, k_pas 6.67 to 6.67 kPa, gamma 7.0 to 7.0, Kse 6667 to 6667 kPa/um, eta 0.133 to 0.133 kPa s/um.

## How well the host itself reproduces its data (the tolerance a new sarcomere inherits)

| quantity | median abs error | max abs error |
|---|---|---|
| LV EDV | 4 % | 8 % |
| LV ESV | 3 % | 32 % |
| LV EDP | 18 % | 175 % |
| LV peak P vs ESP | 5 % | 15 % |
| RV EDV | 3 % | 14 % |
| RV ESV | 3 % | 15 % |
| RV EDP | 44 % | 184 % |
| RV peak P vs ESP | 4 % | 13 % |
