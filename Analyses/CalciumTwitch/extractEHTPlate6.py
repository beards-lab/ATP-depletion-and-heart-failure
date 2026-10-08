"""extractEHTPlate6.py -- Mantarray Plate 6 (Travis, 2026-03-26) -> model-ready CSVs.

Reads the Pulse3D 'full' export of the paced (1 Hz) recording and writes, into
data/travis/processed/:
    plate6_stim_waveforms.csv   time (s) + active twitch force (uN) for 24 wells
    plate6_stim_times.csv       stimulus pulse onsets (s)
    plate6_stim_wellmetrics.csv Pulse3D per-well means (force, timing, spont. rate)
The spontaneous export is read only for the per-well spontaneous beating rate,
which is used to separate treated from control-like wells (see EHT_identification_report.md).

Run from the repo root:  python3 Analyses/CalciumTwitch/extractEHTPlate6.py
"""
import csv
import os

import openpyxl

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
SRC = os.path.join(ROOT, 'data', 'travis')
OUT = os.path.join(SRC, 'processed')
STIM = 'Plate 6_24 hrs_ march 25-26_stimulation_ML26063409-6__2026_03_26_154634__full__2026-03-26_11-47-56.xlsx'
SPON = 'Plate 6_24 hours march 25-26_spontaneousML26063409-6__2026_03_26_154435__full__2026-03-26_11-45-51.xlsx'
WELLS = [r + str(c) for c in range(1, 7) for r in 'ABCD']        # Mantarray column-major order
METRICS = {'Active Twitch Force (μN)': 'F_uN',
           'Twitch Frequency (Hz)': 'freq_Hz',
           'Time From Contraction 10 to Peak (seconds)': 'C10toPk_s',
           'Time From Contraction 50 to Peak (seconds)': 'C50toPk_s',
           'Time From Peak to Relaxation 50 (seconds)': 'PktoR50_s',
           'Time From Peak to Relaxation 90 (seconds)': 'PktoR90_s',
           'Twitch Width 50 (seconds)': 'W50_s'}


def aggregate_means(wb):
    rows = list(wb['aggregate-metrics'].iter_rows(values_only=True))
    assert list(rows[0][2:26]) == WELLS, 'unexpected well order'
    out, cur = {}, None
    for r in rows:
        cur = r[0] or cur
        if r[1] == 'Mean' and cur in METRICS:
            out[METRICS[cur]] = list(r[2:26])
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    wb = openpyxl.load_workbook(os.path.join(SRC, STIM), read_only=True, data_only=True)
    ws = wb['continuous-waveforms']
    it = ws.iter_rows(values_only=True)
    head = next(it)
    assert [h.split(' ')[0] for h in head[1:25]] == WELLS
    wave, stim = [], []
    for r in it:
        if r[0] is not None:
            wave.append(r[:25])
        # stimulus trace is (time, amplitude) pairs; a pulse starts where amplitude steps 0 -> +40
        if r[173] is not None:
            stim.append((r[173], r[174]))
    onsets = sorted({t for (t, a), (_, a0) in zip(stim[1:], stim[:-1]) if a > 0 and a0 == 0})
    with open(os.path.join(OUT, 'plate6_stim_waveforms.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t_s'] + WELLS)
        w.writerows(wave)
    with open(os.path.join(OUT, 'plate6_stim_times.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['stim_onset_s'])
        w.writerows([[t] for t in onsets])

    m = aggregate_means(wb)
    sp = aggregate_means(openpyxl.load_workbook(os.path.join(SRC, SPON), read_only=True, data_only=True))
    with open(os.path.join(OUT, 'plate6_stim_wellmetrics.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        keys = list(METRICS.values())
        w.writerow(['well'] + keys + ['spont_freq_Hz'])
        for i, well in enumerate(WELLS):
            w.writerow([well] + [m[k][i] for k in keys] + [sp['freq_Hz'][i]])
    print(f'{len(wave)} samples, {len(onsets)} stimuli, first {onsets[0]:.4f} s -> {OUT}')


if __name__ == '__main__':
    main()
