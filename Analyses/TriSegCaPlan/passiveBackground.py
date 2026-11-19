"""Passive stress vs SL for skinned mouse trabeculae: relaxed vs Ca-stiffened (PNB+Mava), room-temp series (03_* files).
Run: python -I passiveBackground.py   (reads ../../data/PassiveCaSrc2, writes results/ and figures/)
Columns of the raw .txt: Time (s), L (Lo), F (kPa), SL (um); 4 header lines. F carries an arbitrary offset
(~20 kPa in the PNB files) -> baseline = median F in the slack hold at 0.8 Lo (passive ~0)."""
import numpy as np, glob, os
from scipy.ndimage import median_filter, uniform_filter1d
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', '..', 'data', 'PassiveCaSrc2')
DATES = ['20240705', '20241010', '20241121', '20241212', '20241217', '20241219', '20241220']
SLT = [2.0, 2.1, 2.2, 2.3]

def load(dd, pat):
    g = [x for x in glob.glob(f'{DATA}/{dd}/{pat}') if '10C' not in x]
    return np.loadtxt(g[0], skiprows=4) if g else None

def segs(d):
    L = d[:, 1]; hi = np.where(L > 1.1)[0]; i1 = hi[0]; i2 = hi[-1]
    lo = np.where(L[i2:] < 0.85)[0] + i2; j1 = lo[0]
    ris = np.where(L[j1:] > 0.85)[0]; j2 = j1 + (ris[0] if len(ris) else len(L) - j1 - 1)
    return i1, i2, j1, j2

def baseline(d, j1, j2):
    t = d[:, 0]; m = (t > t[j1] + 5) & (t < t[j2] - 1); return np.median(d[m, 2])

def slow(d):
    i1, i2, j1, j2 = segs(d); b = baseline(d, j1, j2); t = d[:, 0]; L = d[:, 1]
    r0 = np.where(L > L[0] + 0.003)[0][0]; r1 = np.where(L > L[i2] - 0.003)[0][0]
    fs = uniform_filter1d(d[:, 2] - b, 100); sl = median_filter(d[:, 3], size=401)
    return t[r0:r1], sl[r0:r1], fs[r0:r1], b, L[r0:r1]

def curve(sl, f, bin=0.02):
    ok = (sl < 2.45) & (sl > 1.9); sl, f = sl[ok], f[ok]
    e = np.arange(sl.min(), sl.max() + bin, bin); c = []; v = []
    for a in range(len(e) - 1):
        m = (sl >= e[a]) & (sl < e[a + 1])
        if m.sum() > 5: c.append(np.median(sl[m])); v.append(np.median(f[m]))
    return np.maximum.accumulate(np.array(c)), np.array(v)

def hold_end(d):
    i1, i2, j1, j2 = segs(d); b = baseline(d, j1, j2); t = d[:, 0]
    m = (t > t[j1] - 11) & (t < t[j1] - 0.5); return np.median(d[m, 2]) - b

def interp(c, v, x):
    if x < c[0] - 0.01: return np.nan, 'ext'
    if x <= c[-1] + 0.005: return float(np.interp(x, c, v)), 'meas'
    k = c > c[-1] - 0.15
    p = np.polyfit(np.log(np.maximum(c[k] - 1.8, 1e-3)), np.log(np.maximum(v[k], 1e-3)), 1)  # local power law
    return float(np.exp(np.polyval(p, np.log(x - 1.8)))), 'ext'

kp, L0, gam = 11297.1, 1.80093, 7.68907
model = lambda s: kp * np.maximum(np.asarray(s) - L0, 0) ** gam
conds = {'relaxed (no PNB/Mava)': '*100s_RelaxS2F*', 'relaxed PNB+Mava': '*100s_Relax_PNB_Mava.txt',
         'pCa 4.4 PNB+Mava': '*100s_pCa_4.4_PNB_Mava.txt'}
stepc = {'relaxed PNB+Mava': '*0.1s_Relax_PNB_Mava.txt', 'pCa 4.4': '*0.1s_pCa_4.4_PNB_Mava.txt',
         'pCa 6.25': '*0.1s_pCa_6.25_PNB*', 'pCa 6': '*0.1s_pCa_6_PNB*', 'pCa 5.75': '*0.1s_pCa_5.75_PNB*',
         'pCa 5.5': '*0.1s_pCa_5.5_PNB*'}
# The laser SL channel is unreliable per file (all zeros on 20240705, 1.7 um at pCa 4.4 on 20241121), so SL is MAPPED
# from the commanded length: per prep, quadratic SL(L) fit to the relaxed-PNB 100 s ramp; 20240705 uses the pooled mean fit.
def fitSL(dd):
    d = load(dd, conds['relaxed PNB+Mava'])
    if d is None: return None
    t, sl, f, b, L = slow(d)
    if np.median(sl) < 1.5: return None
    return np.polyfit(L, sl, 2)
fits = {dd: fitSL(dd) for dd in DATES}
Lg = np.linspace(0.95, 1.1757, 50)
pool = np.polyfit(Lg, np.mean([np.polyval(p, Lg) for p in fits.values() if p is not None], 0), 2)
for dd in DATES:
    if fits[dd] is None: fits[dd] = pool
res = {}; steps = {}
for dd in DATES:
    res[dd] = {}
    for k, p in conds.items():
        d = load(dd, p)
        if d is None: continue
        t, sl, f, b, L = slow(d); sl = np.polyval(fits[dd], L); c, v = curve(sl, f)
        res[dd][k] = dict(c=c, v=v, sl0=float(c[0]), slmax=float(c[-1]), base=float(b))
    steps[dd] = {'SL_top': float(np.polyval(fits[dd], 1.1757))}
    for k, p in stepc.items():
        d = load(dd, p)
        if d is not None: steps[dd][k] = hold_end(d)

fig, ax = plt.subplots(1, 3, figsize=(17, 5)); cols = {'relaxed (no PNB/Mava)': 'C0', 'relaxed PNB+Mava': 'C2', 'pCa 4.4 PNB+Mava': 'C3'}
sx = np.linspace(1.9, 2.4, 100)
for dd in DATES:
    for k, r in res[dd].items():
        ax[0].plot(r['c'], r['v'], color=cols[k], alpha=.5, lw=1, label=k if dd == DATES[0] else None)
ax[0].plot(sx, model(sx), 'k--', lw=2, label='model law (fastgate_best)')
ax[0].set(xlabel='SL (um)', ylabel='passive stress (kPa)', title='100 s ramps, 7 preps (thin) vs model', ylim=(-1, 25)); ax[0].legend(fontsize=7)
tab = {}
for k in conds:
    cs = [res[dd][k] for dd in DATES if k in res[dd]]
    arr = np.array([[interp(r['c'], r['v'], x)[0] for x in SLT] for r in cs]); fl = np.array([[interp(r['c'], r['v'], x)[1] for x in SLT] for r in cs])
    tab[k] = (arr, fl)
    grid = np.linspace(2.0, min(r['slmax'] for r in cs), 30)
    ax[1].plot(grid, np.mean([np.interp(grid, r['c'], r['v']) for r in cs], 0), color=cols[k], lw=2, label=k + ' (mean)')
ax[1].plot(sx, model(sx), 'k--', lw=2, label='model')
ax[1].set(xlabel='SL (um)', title='mean over preps (to min of SLmax)', ylim=(-1, 25)); ax[1].legend(fontsize=7)
pc = [('relaxed PNB+Mava', 9.0), ('pCa 6.25', 6.25), ('pCa 6', 6.0), ('pCa 5.75', 5.75), ('pCa 5.5', 5.5), ('pCa 4.4', 4.4)]
for dd in DATES:
    ax[2].plot([x for _, x in pc], [steps[dd].get(k, np.nan) for k, _ in pc], 'o-', alpha=.5, lw=1, label=dd)
mu = [np.nanmean([steps[dd].get(k, np.nan) for dd in DATES]) for k, _ in pc]
ax[2].plot([x for _, x in pc], mu, 'ks-', lw=2, label='mean'); ax[2].invert_xaxis()
ax[2].set(xlabel='pCa (9 = relaxed)', ylabel='stress 60 s after 0.1 s step (kPa)', title='stress at top SL (~2.25 um) vs pCa'); ax[2].legend(fontsize=6)
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'figures', 'passive_ca.png'), dpi=110)

L = ['| prep | cond | SL range of ramp (um) | ' + ' | '.join(str(x) for x in SLT) + ' |', '|' + '---|' * (3 + len(SLT))]
for dd in DATES:
    for k in conds:
        if k not in res[dd]: continue
        r = res[dd][k]; vals = []
        for x in SLT:
            v, fl = interp(r['c'], r['v'], x); vals.append('-' if np.isnan(v) else (f'{v:.1f}' if fl == 'meas' else f'({v:.1f})'))
        L.append(f"| {dd} | {k} | {r['sl0']:.2f}-{r['slmax']:.2f} | " + ' | '.join(vals) + ' |')
means = []
for k in conds:
    arr, fl = tab[k]
    means.append(f'| {k} | ' + ' | '.join(f'{np.nanmean(arr[:, i]):.2f} (n={np.sum(~np.isnan(arr[:, i]))}, {int((fl[:, i] == "meas").sum())} meas)' for i in range(len(SLT))) + ' |')
a2 = tab['relaxed PNB+Mava'][0]; a3 = tab['pCa 4.4 PNB+Mava'][0]
rr = [f'SL {x}: ratio pCa4.4/relaxed(PNB) per prep: ' + ', '.join(f'{u:.2f}' for u in a3[:, i] / a2[:, i]) + f' | mean-of-ratios {np.nanmean(a3[:, i] / a2[:, i]):.2f}; model {model(x):.2f} kPa' for i, x in enumerate(SLT)]
st = ['SL_top per prep: ' + ', '.join(f"{dd}:{steps[dd]['SL_top']:.3f}" for dd in DATES)]
for k, _ in pc:
    st.append(f'{k}: ' + ', '.join(f"{steps[dd].get(k, np.nan):.1f}" for dd in DATES) + f'  mean {np.nanmean([steps[dd].get(k, np.nan) for dd in DATES]):.1f}')
txt = 'PER PREP\n' + '\n'.join(L) + '\n\nMEANS\n' + '\n'.join(means) + '\n\nRATIOS\n' + '\n'.join(rr) + '\n\nSTEP HOLD-END STRESS\n' + '\n'.join(st) + '\n'
open(os.path.join(HERE, 'results', 'tables.txt'), 'w').write(txt); print(txt)
