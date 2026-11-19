"""mouseHostFigures.py -- figures for RunMouseHost.m (#82, Q7): the mouse host against the P21 host,
and the cross-bridge wall on the mouse host's beats.   python -I mouseHostFigures.py [P21 host root]"""
import os, sys, numpy as np, scipy.io as sio
from scipy.interpolate import PchipInterpolator
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); SP = os.path.join(HERE, 'results'); FIG = os.path.join(HERE, 'figures')
P21 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', '..', '..', 'Multiscale-Hx-P21-Model')
ld = lambda f: sio.loadmat(f, squeeze_me=True, struct_as_record=False)

def beat(t, T, cols, dt=0.25e-3):
    """last full beat on a uniform grid"""
    k = t >= t[-1] - T - 1e-9; tb = t[k] - t[k][0]
    tb, iu = np.unique(tb, return_index=True); tu = np.arange(0, T, dt)
    return tu, [PchipInterpolator(tb, np.asarray(c)[k][iu])(tu) for c in cols]

def smooth(x, n):  # periodic moving mean
    return np.convolve(np.r_[x[-n:], x, x[:n]], np.ones(n) / n, 'same')[n:-n]

H = {}
S = ld(os.path.join(P21, 'outputs', 'outputs_Nx51.mat')); o, d = S['outputs'], S['data']
t = np.asarray(o.time, float); T = 60 / float(d.HR)
tu, c = beat(t, T, [2.0 * np.exp(np.asarray(getattr(o.strains, 'eps_' + w), float)) for w in ['LV', 'SEP', 'RV']] +
             [np.asarray(o.stresses.active.sigma_XB_LV, float), np.asarray(o.stresses.passive.sigma_pas_LV, float), np.asarray(o.stresses.total.sigma_LV, float),
              np.asarray(o.volumes.V_LV, float) * 1e3, np.asarray(o.pressures.P_LV, float) / 7.50062])
H['P21 rat Nx51'] = dict(T=T, HR=float(d.HR), t=tu, Ls=c[:3], act=c[3], pas=c[4], tot=c[5], V=c[6], P=c[7])
for nm, f in [('mouse, literature targets', 'mouseHost_mouse.mat'), ('mouse, one echo fit', 'mouseHost_m22146.mat')]:
    M = ld(os.path.join(SP, f))['S']; t = np.asarray(M.t, float); T = float(M.T)
    tu, c = beat(t, T, [M.Ls_LV, M.Ls_SEP, M.Ls_RV, np.asarray(M.sact_LV) / 7.50062, np.asarray(M.spas_LV) / 7.50062, np.asarray(M.s_LV) / 7.50062, np.asarray(M.V_LV) * 1e3, np.asarray(M.P_LV) / 7.50062])
    H[nm] = dict(T=T, HR=float(M.HR), t=tu, Ls=c[:3], act=c[3], pas=c[4], tot=c[5], V=c[6], P=c[7])

print('%-26s %5s %6s | %-13s %-13s %-13s | %s' % ('host', 'HR', 'T ms', 'LV Ls', 'SEP Ls', 'RV Ls', 'LV: peak shortening um/s (lengths/s) raw | 1 ms mean | 5 ms mean; peak lengthening 1 ms mean; time of shortening 10-90 % ms'))
fig, ax = plt.subplots(2, 3, figsize=(15, 7.5))
for i, (nm, h) in enumerate(H.items()):
    dt = h['t'][1] - h['t'][0]; Ls = h['Ls'][0]; v = np.gradient(Ls, dt)
    v1 = smooth(v, int(round(1e-3 / dt))); v5 = smooth(v, int(round(5e-3 / dt)))
    a = Ls.max() - Ls.min(); i0 = np.argmax(Ls); seq = np.r_[Ls[i0:], Ls[:i0]]
    t10 = np.argmax(seq <= Ls.max() - 0.1 * a) * dt; t90 = np.argmax(seq <= Ls.max() - 0.9 * a) * dt
    rng = lambda x: '%.3f-%.3f' % (x.min(), x.max())
    print('%-26s %5.0f %6.1f | %-13s %-13s %-13s | %.1f (%.2f) | %.1f (%.2f) | %.1f (%.2f); +%.1f; %.1f' % (nm, h['HR'], 1e3 * h['T'], rng(h['Ls'][0]), rng(h['Ls'][1]), rng(h['Ls'][2]),
          -v.min(), -v.min() / Ls.max(), -v1.min(), -v1.min() / Ls.max(), -v5.min(), -v5.min() / Ls.max(), v1.max(), 1e3 * (t90 - t10)))
    print('      LV stress kPa: active max %.1f, passive max %.2f (at end diastole %.2f), total max %.1f; LV volume %.1f-%.1f uL (EF %.0f %%); LV pressure %.2f-%.1f kPa' % (
        h['act'].max(), h['pas'].max(), h['pas'][np.argmax(h['V'])], h['tot'].max(), h['V'].min(), h['V'].max(), 100 * (1 - h['V'].min() / h['V'].max()), h['P'].min(), h['P'].max()))
    ph = h['t'] / h['T']
    ax[0, 0].plot(1e3 * h['t'], Ls, label=nm); ax[0, 1].plot(1e3 * h['t'], v1, label=nm); ax[0, 2].plot(1e3 * h['t'], h['tot'], label=nm)
    ax[1, 0].plot(ph, Ls); ax[1, 1].plot(1e3 * h['t'], h['pas'], label=nm); ax[1, 2].plot(h['V'] / h['V'].max(), h['P'], label=nm)
ax[0, 0].set(xlabel='t (ms)', ylabel='LV sarcomere length (um)', title='LV Ls'); ax[0, 0].legend(fontsize=8)
ax[0, 1].axhline(-4.6 * 2.2, color='k', ls=':'); ax[0, 1].text(100, -4.6 * 2.2 + 0.4, 'XB model Vmax 4.6 ML/s', fontsize=8)
ax[0, 1].set(xlabel='t (ms)', ylabel='dLs/dt (um/s), 1 ms mean', title='LV sarcomere velocity')
ax[0, 2].set(xlabel='t (ms)', ylabel='LV total stress (kPa)', title='Host wall stress')
ax[1, 0].set(xlabel='fraction of the beat', ylabel='LV Ls (um)', title='LV Ls against phase')
ax[1, 1].set(xlabel='t (ms)', ylabel='LV passive stress (kPa)', title='Host passive stress')
ax[1, 2].set(xlabel='LV volume / EDV', ylabel='LV pressure (kPa)', title='PV loop (volume normalised)')
for a in ax.ravel(): a.grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'mouse_host_kinematics.png'), dpi=120); plt.close(fig)

# ---- the cross-bridge wall on the mouse host beats
r = sio.loadmat(os.path.join(SP, 'mouseHostProbe.mat'), squeeze_me=True, struct_as_record=False)['res']
fig, ax = plt.subplots(2, 3, figsize=(16, 8))
for i, (cs, ttl) in enumerate([('mouse', 'mouse, literature targets (409 bpm)'), ('m22146', 'mouse, one echo fit (579 bpm)')]):
    L = getattr(r, 'Ls_' + cs); T = float(getattr(r, 'T_' + cs)); Hr = getattr(r, 'H_' + cs)
    a = ax[i, 0]; a.plot(1e3 * L[:, 0], L[:, 1], 'k'); a.set(xlabel='t (ms)', ylabel='LV Ls (um)', title=ttl); a.grid(alpha=.3)
    b = a.twinx(); b.plot(1e3 * L[:, 0], L[:, 3], 'C1', lw=1); b.set_ylabel('Ca (uM), borrowed from the P21 host', color='C1')
    a = ax[i, 1]; a.plot(1e3 * L[:, 0], L[:, 2], 'r', label='host dLs/dt')
    R = getattr(r, cs + '_kon3'); k = R.t >= 3 * T - 1e-9
    a.plot(1e3 * R.tk, R.velHS[k], 'b', lw=1, label='XB contractile element')
    a.axhline(-4.6 * 2.0, color='k', ls=':'); a.text(0.55 * 1e3 * T, -4.6 * 2.0 + 0.6, 'XB model Vmax 4.6 ML/s', fontsize=8)
    a.set(xlabel='t (ms)', ylabel='um/s', title='LV velocity', ylim=(-25, 40)); a.legend(fontsize=8); a.grid(alpha=.3)
    a = ax[i, 2]; a.plot(1e3 * np.asarray(Hr.t), np.asarray(Hr.LV.sigma_XB) + np.asarray(Hr.LV.sigma_pas), 'k--', lw=1.5, label='host total')
    for nm, lab in [('kon3', 'XB model, k_on x3'), ('kon10', 'k_on x10'), ('kon3_x5', 'k_on x3, cycle rates x5')]:
        R = getattr(r, cs + '_' + nm); k = R.t >= 3 * T - 1e-9; a.plot(1e3 * R.tk, R.F[k], label=lab)
    a.set(xlabel='t (ms)', ylabel='LV total stress (kPa)', title='Stress on the host beat'); a.legend(fontsize=8); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'mouse_host_probe.png'), dpi=120); print('figures written')
