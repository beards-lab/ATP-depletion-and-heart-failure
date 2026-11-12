"""passiveDynamic.py -- does Ca stiffen the FAST passive response? (#82 Q3, background for #85)

passiveBackground.py showed that on 100 s ramps and 60 s into a hold the passive stress is
the same relaxed and at pCa 4.4 (PNB + Mava). This script looks at the transient: the peak
stress at the end of a ramp to 1.1757 Lo (top SL about 2.27 um) for ramp durations 100, 10,
1 and 0.1 s, relaxed PNB+Mava against pCa 4.4 PNB+Mava, same 7 preparations.

    python -I passiveDynamic.py     -> background/passive_dynamic.md
Raw columns: Time (s), L (Lo), F (kPa), SL (um); 4 header lines; F has an arbitrary offset,
removed with the hold at 0.95 Lo before the ramp, so every stress is the rise above 0.95 Lo
(the passive stress there is about 0.3 kPa).
"""
import glob, os
import numpy as np
from scipy.ndimage import uniform_filter1d

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', '..', 'data', 'PassiveCaSrc2')
DATES = ['20240705', '20241010', '20241121', '20241212', '20241217', '20241219', '20241220']
DUR = ['100', '10', '1', '0.1']
COND = {'relaxed': '*_%ss_Relax_PNB_Mava.txt', 'pCa 4.4': '*_%ss_pCa_4.4_PNB_Mava.txt'}


def load(dd, pat):
    g = [x for x in glob.glob(os.path.join(DATA, dd, pat)) if '10C' not in x]
    return np.loadtxt(g[0], skiprows=4) if g else None


def measure(d, dur):
    t, L, F = d[:, 0], d[:, 1], d[:, 2]
    Ltop = np.median(L[L > 1.1])
    r0 = np.where(L > L[0] + 0.003)[0][0]                               # start of the ramp
    r1 = np.where(L > Ltop - 0.003)[0][0]                               # end of the ramp
    b = np.median(F[t < t[r0] - 1])                                     # baseline: the hold at 0.95 Lo before the ramp
    dt = np.median(np.diff(t)); f = uniform_filter1d(F - b, max(1, int(round(0.005 / dt))))
    w = (t >= t[r1] - max(0.02, 0.01 * dur)) & (t <= t[r1] + max(0.05, 0.02 * dur))
    h = (t > t[r1] + 50) & (t < t[r1] + 60) & (L > Ltop - 0.003)        # 50-60 s into the hold
    hold = np.median(F[h]) - b if h.sum() > 10 else np.nan
    return f[w].max(), hold


rows = {}
for dd in DATES:
    for du in DUR:
        for c, pat in COND.items():
            d = load(dd, pat % du)
            if d is not None: rows[(dd, du, c)] = measure(d, float(du))

out = ['# Passive stress transient: relaxed against pCa 4.4 (PNB + Mava), by ramp duration', '',
       'Script: `../passiveDynamic.py`. Ramp from 0.95 to 1.1757 Lo (top SL about 2.27 um), 7 preparations.',
       'Stress is the rise above the 0.95 Lo hold. Peak = highest stress around the end of the ramp (5 ms smoothing);', 'hold = stress 50 to 60 s into the hold.', '',
       '| ramp duration (s) | n | peak relaxed (kPa) | peak pCa 4.4 (kPa) | pCa 4.4 / relaxed, mean of per-prep ratios (range) | peak / hold, relaxed | peak / hold, pCa 4.4 |',
       '|---|---|---|---|---|---|---|']
for du in DUR:
    pr = [(rows[(dd, du, 'relaxed')], rows[(dd, du, 'pCa 4.4')]) for dd in DATES if (dd, du, 'relaxed') in rows and (dd, du, 'pCa 4.4') in rows]
    if not pr: continue
    a = np.array([[p[0][0], p[0][1], p[1][0], p[1][1]] for p in pr]); ratio = a[:, 2] / a[:, 0]
    out.append('| %s | %d | %.1f ± %.1f | %.1f ± %.1f | %.2f (%.2f to %.2f) | %.2f | %.2f |' % (
        du, len(pr), a[:, 0].mean(), a[:, 0].std(ddof=1), a[:, 2].mean(), a[:, 2].std(ddof=1), ratio.mean(), ratio.min(), ratio.max(),
        np.nanmean(a[:, 0] / a[:, 1]), np.nanmean(a[:, 2] / a[:, 3])))
out += ['', '## Per preparation, 0.1 s ramp (kPa): peak relaxed, peak pCa 4.4, hold relaxed, hold pCa 4.4', '']
for dd in DATES:
    if (dd, '0.1', 'relaxed') in rows and (dd, '0.1', 'pCa 4.4') in rows:
        r, c = rows[(dd, '0.1', 'relaxed')], rows[(dd, '0.1', 'pCa 4.4')]
        out.append('- %s: %.1f, %.1f, %.1f, %.1f' % (dd, r[0], c[0], r[1], c[1]))
txt = '\n'.join(out) + '\n'
os.makedirs(os.path.join(HERE, 'background'), exist_ok=True)
open(os.path.join(HERE, 'background', 'passive_dynamic.md'), 'w', encoding='utf-8').write(txt); print(txt)
