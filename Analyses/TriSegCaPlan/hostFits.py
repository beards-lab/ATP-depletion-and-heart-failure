"""hostFits.py -- what the P21 host's saved fits look like, per animal (#82, Q7 and Q9).

Reads ../../../Multiscale-Hx-P21-Model/outputs/outputs_*.mat (the host's own fitted
runs, two beats each; the last beat is used) and tabulates, per fit:
  the data the host was fitted to (HR, EDV, ESV, EDP, ESP of both ventricles),
  what the host model gives for the same quantities,
  the sarcomere kinematics it implies (Ls at end diastole and its minimum, per wall),
  and its stresses (peak cross-bridge, peak passive, peak total).
The host's data are in kPa and its stored pressure traces in mmHg; everything here is kPa.
Volumes mL, lengths um.

    python -I hostFits.py [host root]      -> background/host_fits.md, figures/host_fits.csv
"""
import glob, os, sys
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
HOST = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', '..', '..', 'Multiscale-Hx-P21-Model')
rows = []
for f in sorted(glob.glob(os.path.join(HOST, 'outputs', 'outputs_*.mat'))):
    S = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
    o, d = S['outputs'], S['data']
    t = np.asarray(o.time, float); T = 60.0 / float(d.HR)
    k = t >= t[-1] - T - 1e-9
    Lsref = float(d.fixpars[12])
    r = dict(fit=os.path.basename(f)[8:-4], HR=float(d.HR))
    for v in ['LV', 'RV']:
        V = np.asarray(getattr(o.volumes, 'V_' + v), float)[k]; P = np.asarray(getattr(o.pressures, 'P_' + v), float)[k] / 7.50062   # stored in mmHg
        r.update({'EDV_%s_data' % v: float(getattr(d, 'EDV_' + v)), 'ESV_%s_data' % v: float(getattr(d, 'ESV_' + v)),
                  'EDP_%s_data' % v: float(getattr(d, 'EDP_' + v)), 'ESP_%s_data' % v: float(getattr(d, 'ESP_' + v)),
                  'EDV_%s' % v: V.max(), 'ESV_%s' % v: V.min(), 'EDP_%s' % v: P[np.argmax(V)], 'Pmax_%s' % v: P.max()})
    for w in ['LV', 'SEP', 'RV']:
        Ls = Lsref * np.exp(np.asarray(getattr(o.strains, 'eps_' + w), float)[k])
        v = np.gradient(Ls, t[k]); iv = np.argmin(v)
        r.update({'LsMax_' + w: Ls.max(), 'LsMin_' + w: Ls.min(), 'vShort_' + w: -v.min() / Ls.max(), 'tShort_' + w: (t[k][iv] - t[k][0]) / T,
                  'vShortMed_' + w: -np.median(v[v < 0]) / Ls.max(),
                  'sXB_' + w: np.asarray(getattr(o.stresses.active, 'sigma_XB_' + w), float)[k].max(),
                  'sPas_' + w: np.asarray(getattr(o.stresses.passive, 'sigma_pas_' + w), float)[k].max(),
                  'sTot_' + w: np.asarray(getattr(o.stresses.total, 'sigma_' + w), float)[k].max()})
    Ca = np.asarray(o.XB.Ca_i, float)[k]; r.update(CaMin=Ca.min(), CaMax=Ca.max(), Lsref=Lsref,
        Lsc0=float(d.fixpars[13]), k_pas=float(d.fixpars[17]), gamma=float(d.fixpars[19]), Kse=float(d.fixpars[26]), eta=float(d.fixpars[25]))
    rows.append(r)

keys = list(rows[0].keys())
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True); os.makedirs(os.path.join(HERE, 'background'), exist_ok=True)
with open(os.path.join(HERE, 'figures', 'host_fits.csv'), 'w') as fh:
    fh.write(','.join(keys) + '\n')
    for r in rows: fh.write(','.join('%s' % r[k] if isinstance(r[k], str) else '%.5g' % r[k] for k in keys) + '\n')

A = lambda k: np.array([r[k] for r in rows])
rng = lambda k, f='%.2f': (f + ' to ' + f) % (A(k).min(), A(k).max())
out = ['# The P21 host fits (%d animals): data, model, sarcomere kinematics, stresses' % len(rows), '',
       'Script: `../hostFits.py`; per-animal table: `../figures/host_fits.csv`. Last beat of each saved run.',
       'Pressures kPa (1 kPa = 7.5 mmHg), volumes mL, lengths um, stresses kPa.', '',
       '| fit | HR | LV EDV data / model | LV ESV data / model | LV EDP data / model | LV ESP data, model Pmax | Ls LV max / min | Ls RV max / min | peak XB / passive / total LV |',
       '|---|---|---|---|---|---|---|---|---|']
for r in rows:
    out.append('| %s | %.0f | %.3f / %.3f | %.3f / %.3f | %.2f / %.2f | %.1f, %.1f | %.2f / %.2f | %.2f / %.2f | %.1f / %.2f / %.1f |' % (
        r['fit'], r['HR'], r['EDV_LV_data'], r['EDV_LV'], r['ESV_LV_data'], r['ESV_LV'], r['EDP_LV_data'], r['EDP_LV'], r['ESP_LV_data'], r['Pmax_LV'],
        r['LsMax_LV'], r['LsMin_LV'], r['LsMax_RV'], r['LsMin_RV'], r['sXB_LV'], r['sPas_LV'], r['sTot_LV']))
rel = lambda a, b: 100 * np.abs(A(a) / A(b) - 1)
out += ['', '## Ranges over the %d fits' % len(rows), '',
        '- HR %s bpm; LV EDV (data) %s mL; LV ESP (data) %s kPa; LV EDP (data) %s kPa.' % (rng('HR', '%.0f'), rng('EDV_LV_data', '%.3f'), rng('ESP_LV_data', '%.1f'), rng('EDP_LV_data', '%.2f')),
        '- Sarcomere length: LV max %s, LV min %s; SEP max %s, min %s; RV max %s, min %s um.' % (rng('LsMax_LV'), rng('LsMin_LV'), rng('LsMax_SEP'), rng('LsMin_SEP'), rng('LsMax_RV'), rng('LsMin_RV')),
        '- Shortening per beat (max - min): LV %.2f to %.2f um.' % ((A('LsMax_LV') - A('LsMin_LV')).min(), (A('LsMax_LV') - A('LsMin_LV')).max()),
        '- Peak shortening speed: LV %s, RV %s lengths/s (of the maximum length), at %s of the beat (LV); median speed while shortening: LV %s lengths/s.' % (rng('vShort_LV'), rng('vShort_RV'), rng('tShort_LV'), rng('vShortMed_LV')),
        '- Peak cross-bridge stress: LV %s, SEP %s, RV %s kPa.' % (rng('sXB_LV', '%.1f'), rng('sXB_SEP', '%.1f'), rng('sXB_RV', '%.1f')),
        '- Peak passive stress: LV %s, SEP %s, RV %s kPa.' % (rng('sPas_LV'), rng('sPas_SEP'), rng('sPas_RV')),
        '- Ca transient: diastolic %s, peak %s uM.' % (rng('CaMin', '%.3f'), rng('CaMax', '%.2f')),
        '- Fixed parameters (same in every fit unless a range is shown): Lsref %s, Lsc0 %s um, k_pas %s kPa, gamma %s, Kse %s kPa/um, eta %s kPa s/um.' % (
            rng('Lsref'), rng('Lsc0'), rng('k_pas', '%.3g'), rng('gamma', '%.1f'), rng('Kse', '%.4g'), rng('eta', '%.3g')),
        '', '## How well the host itself reproduces its data (the tolerance a new sarcomere inherits)', '',
        '| quantity | median abs error | max abs error |', '|---|---|---|']
for nm, a, b in [('LV EDV', 'EDV_LV', 'EDV_LV_data'), ('LV ESV', 'ESV_LV', 'ESV_LV_data'), ('LV EDP', 'EDP_LV', 'EDP_LV_data'), ('LV peak P vs ESP', 'Pmax_LV', 'ESP_LV_data'),
                  ('RV EDV', 'EDV_RV', 'EDV_RV_data'), ('RV ESV', 'ESV_RV', 'ESV_RV_data'), ('RV EDP', 'EDP_RV', 'EDP_RV_data'), ('RV peak P vs ESP', 'Pmax_RV', 'ESP_RV_data')]:
    e = rel(a, b); out.append('| %s | %.0f %% | %.0f %% |' % (nm, np.median(e), e.max()))
txt = '\n'.join(out) + '\n'
open(os.path.join(HERE, 'background', 'host_fits.md'), 'w', encoding='utf-8').write(txt); print(txt)
