"""analyzeFvid.py -- Python twin of AnalyzeFvid.m (same method, same numbers).

Written because MATLAB could not start when the Jacobian runs finished (online
licence hang). It reads the checkpoints of RunFvidJacobian.m and the data
features AnalyzeFvid.m had already computed (results/fvid.mat, or fvid_partial.mat), and
re-implements readIsoV + beardVector for the model traces. The port is checked
against the MATLAB feature vectors in that file before anything is used.

    python analyzeFvid.py        (Python 3.10, numpy, scipy, matplotlib)

Output -> figures/fvid_*.png, figures/verdicts.csv, results/fvid_py.npz, tables on stdout.
"""
import os, sys, warnings
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.simplefilter('ignore')
np.set_printoptions(linewidth=200)
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results'); FIG = os.path.join(HERE, 'figures')
LN = np.log(1.10); SIG_PRIOR = np.log(10); THR = np.log(2)
RAMP = np.array([0.5, 1, 2, 3, 5, 8, 12, 20, 35, 60, 100, 150])
REDEV = np.array([2, 5, 10, 20, 35, 60, 100, 200, 350])
ld = lambda f: sio.loadmat(f, squeeze_me=True, struct_as_record=False)
mround = lambda x: np.floor(np.asarray(x) + 0.5)          # MATLAB round for positive numbers


def movmean(x, k):
    """MATLAB movmean(x, k) with shrinking endpoints."""
    n = len(x); b, a = ((k - 1) // 2, (k - 1) // 2) if k % 2 else (k // 2, k // 2 - 1)
    c = np.concatenate(([0.0], np.cumsum(x))); i = np.arange(n)
    lo = np.maximum(i - b, 0); hi = np.minimum(i + a, n - 1)
    return (c[hi + 1] - c[lo]) / (hi - lo + 1)


def iso_plan(t, L, half_win=0.5, Lx=1.0):
    """The part of readIsoV that depends on t and L only: per cycle the sample
    masks for Fx and F0, plus tStart, tSlack, vMeas. Last entry = the v = 0 row."""
    Lhold = np.median(L[t < 400]); below = L < Lhold - 0.01
    on = np.flatnonzero(below & ~np.concatenate(([False], below[:-1])))
    keep = []
    for i0 in on:
        if not (~below[:i0]).any(): keep.append(True); continue
        j0 = np.flatnonzero(t >= t[i0] - 50)[0]
        keep.append(not below[j0:i0].any())
    on = on[np.array(keep)]; n = len(on); dt = np.median(np.diff(t))
    vS = -movmean(np.concatenate(([0.0], np.diff(L))) / dt * 1000, max(1, int(mround(1 / dt))))
    P = dict(tStart=[], tSlack=[], vMeas=[], wFx=[], wF0=[])
    for i in range(n):
        i0 = on[i]; iEnd = len(t) - 1 if i == n - 1 else on[i + 1] - 1
        seg = np.arange(i0, iEnd + 1)
        ic = seg[np.flatnonzero(L[seg] <= Lx)[0]]
        w = np.abs(t - t[ic]) <= half_win
        s = np.flatnonzero((t[seg] > t[ic]) & (vS[seg] > 30))
        tSl = t[seg[s[0]]] if s.size else t[iEnd]
        wr = seg[(t[seg] >= t[i0] + 2) & (t[seg] <= tSl - 1) & (L[seg] > Lx - 0.03)]
        wr = wr[t[wr] >= t[ic] - max(10, 0.5 * (t[ic] - t[i0]))]
        v = -np.polyfit(t[wr] / 1000, L[wr], 1)[0]
        sl = seg[(t[seg] > tSl + 3) & (L[seg] <= L[seg].min() + 0.01)]
        P['tStart'].append(t[i0]); P['tSlack'].append(tSl); P['vMeas'].append(v); P['wFx'].append(np.flatnonzero(w)); P['wF0'].append(sl)
    w = t > t[-1] - 200
    assert abs(np.median(L[w]) - Lx) < 0.01 and n == 7
    P['vMeas'].append(0.0); P['wFx'].append(np.flatnonzero(w)); P['wF0'].append(P['wF0'][n - 1])
    P['vMeas'] = np.array(P['vMeas']); P['iv'] = np.argsort(P['vMeas'], kind='stable')
    return P


def trace_plan(t, L, P, atp):
    """The sample plan of beardVector (build mode): 0-based indices and labels."""
    dt = np.median(np.diff(t)); Lh = np.median(L[t < 400]); vNom = mround(P['vMeas'] * 2) / 2
    add = lambda tt: int(mround(tt / dt))
    idx = [add(300)]; lab = [(atp, np.nan, 'hold110', 300.0)]
    for i in range(7):
        for d in RAMP[RAMP < P['tSlack'][i] - P['tStart'][i] - 0.5]:
            idx.append(add(P['tStart'][i] + d)); lab.append((atp, vNom[i], 'ramp', d))
        Lt, tNext = Lh, t[-1]
        if i < 6: tNext = P['tStart'][i + 1]
        else: Lt = np.median(L[t > 4100])
        jr = np.flatnonzero((t > P['tSlack'][i] + 2) & (L >= Lt - 0.002))[0]
        for d in REDEV[t[jr] + REDEV < tNext - 2]:
            idx.append(add(t[jr] + d)); lab.append((atp, vNom[i], 'redev', d))
    return np.array(idx), lab


def beard_vector(F3, plans, idxs, avail):
    """pts (23) and trc (361) of one evaluation; NaN where a record is missing."""
    pts, trc, Fref = [], [], np.nan
    for a in range(3):
        P = plans[a]; npt = 7 if a == 0 else 8
        if not avail[a] or not avail[0]:
            pts.append(np.full(npt, np.nan)); trc.append(np.full(len(idxs[a]), np.nan)); continue
        f = F3[a].astype(float)
        Fx = np.array([f[w].mean() for w in P['wFx']]); F0 = np.array([f[w].mean() for w in P['wF0']])
        if a == 0: Fref = Fx[-1] - F0[-1]
        p = (Fx[P['iv']] - F0[P['iv']]) / Fref
        pts.append(p[1:] if a == 0 else p)
        trc.append(((movmean(f, 21) - np.median(F0[:7])) / Fref)[idxs[a]])
    return np.concatenate(pts), np.concatenate(trc)


def param_group(nm):
    g = [('cycle rates', ['ka', 'kd', 'k1', 'k_1', 'k2', 'k_2', 'kah', 'kamh']),
         ('cross-bridge mechanics', ['kstiff1', 'kstiff2', 'kstiff1_n', 'kstiff2_n', 'estiff', 'dr', 'dr1', 'dr2', 'mu', 'mu_neg']),
         ('series elasticity', ['kSE', 'ekSE', 'kSE_M', 'eta_M', 'x_M_slack']),
         ('SRX', ['ksr0', 'kmsrd', 'ksr2srd', 'ksrd2sr', 'sigma2', 'sigma_srd1', 'ksrd', 'sigma_srd2', 'SRXFromR2HighStrain', 'sSRXrip']),
         ('regulation', ['A0', 'v_ref_reg', 'tau_reg', 'P_bound_max']), ('recoil', ['k_cr', 'tau_cr', 'F_cr']),
         ('geometry', ['L_thick', 'L_thin', 'd_optimal', 'sigma_lattice', 'Lsc0']), ('passive', ['gamma', 'k_pas']), ('ATP gate', ['pAtp_k2'])]
    if nm.startswith('PieceWise'): return 'strain-dependence shape'
    for k, v in g:
        if nm in v: return k
    return 'other'


# ---------------------------------------------------------------- 1. load
M = ld(os.path.join(RES, 'fvid.mat' if os.path.isfile(os.path.join(RES, 'fvid.mat')) else 'fvid_partial.mat'))   # data features + MATLAB vectors for the check
dPts, semPts, dTrc, semTrc = M['dPts'], M['semPts'], M['dTrc'], M['semTrc']
Vc, dBat, semBat, labBat = M['Vc'], M['dBat'], M['semBat'], np.array([str(x) for x in M['labBat']])
okB = M['okB'].astype(bool); Vc = Vc[:, okB]; dBat = dBat[okB]; semBat = semBat[okB]; labBat = labBat[okB]
recs = ld(os.path.join(RES, 'beard_records.mat'))['recs']
Eb = ld(os.path.join(RES, 'jac_beard', 'evals.mat')); Ec = ld(os.path.join(RES, 'jac_battery', 'evals.mat'))
EB, EC = Eb['E'], Ec['E']; poolC = [str(x) for x in Ec['pool']]
kindB = np.array([str(e.kind) for e in EB]); kindC = np.array([str(e.kind) for e in EC])
isname = lambda e, nm: isinstance(e.name, str) and e.name == nm
names = [str(x) for x in Eb['pool']]
names += [str(e.name) for e in EB if str(e.kind) == 'jac' and isinstance(e.name, str) and e.name not in names]   # pAtp_k2
P_ = len(names)

plans, idxs, labs = [], [], []
for a in range(3):
    t = recs[a].t.astype(float); L = recs[a].L.astype(float)
    P = iso_plan(t, L); ix, lab = trace_plan(t, L, P, int(recs[a].atp))
    assert np.array_equal(ix, np.asarray(M['S'].idx[a]).astype(int) - 1), 'sample plan differs from MATLAB (record %d)' % a
    plans.append(P); idxs.append(ix); labs += lab
trATP = np.array([l[0] for l in labs]); trV = np.array([l[1] for l in labs]); trKind = np.array([l[2] for l in labs]); trDt = np.array([l[3] for l in labs])
ptATP = np.array([8] * 7 + [4] * 8 + [2] * 8); ptV = np.array([.5, 1, 2, 3, 4, 5, 6] + [0, .5, 1, 2, 3, 4, 5, 6] * 2)
nPts = np.where(ptATP == 8, 12, 6); nTrc = np.where(trATP == 8, 12, 6)

Vp = np.full((len(EB), 23), np.nan); Vt = np.full((len(EB), len(labs)), np.nan); nfail = 0
for e in range(len(EB)):
    F3, av = [None] * 3, [False] * 3
    for a in range(3):
        f = os.path.join(RES, 'jac_beard', 'e%03d_r%d.mat' % (e + 1, a + 1))
        if not os.path.isfile(f): continue
        r = ld(f)['r']
        if r.ok: F3[a] = np.asarray(r.F); av[a] = True
        else: nfail += 1
    Vp[e], Vt[e] = beard_vector(F3, plans, idxs, av)
print('Beard evaluations complete: %d of %d (failed simulations: %d). Battery: %d of %d files, %d with every feature present.' % (
    np.isfinite(np.c_[Vp, Vt]).all(1).sum(), len(EB), nfail, np.isfinite(Vc).any(1).sum(), len(EC), np.isfinite(Vc).all(1).sum()))

# ---- port check against the MATLAB vectors (wherever MATLAB had a value)
mp, mt = M['Vp'], M['Vt']; okp, okt = np.isfinite(mp), np.isfinite(mt)
print('PORT CHECK vs MATLAB: %d point and %d trace values compared; max |difference| %.2e and %.2e (values are O(1))' % (
    okp.sum(), okt.sum(), np.abs(Vp[okp] - mp[okp]).max(), np.abs(Vt[okt] - mt[okt]).max()))
assert np.abs(Vp[okp] - mp[okp]).max() < 1e-6 and np.abs(Vt[okt] - mt[okt]).max() < 1e-6

# ---------------------------------------------------------------- 2. Jacobians
def col(V, E, kind, nm, sem, k, h):
    i = [j for j in range(len(E)) if kind[j] == k and isname(E[j], nm)]
    return (V[i[0]] - V[0]) / sem / h if i else np.full(V.shape[1], np.nan)

Jp = np.full((23, P_), np.nan); Jt = np.full((len(labs), P_), np.nan); Jc = np.zeros((Vc.shape[1], P_)); vanish = {}
for j, nm in enumerate(names):
    Jp[:, j] = col(Vp, EB, kindB, nm, semPts, 'jac', LN); Jt[:, j] = col(Vt, EB, kindB, nm, semTrc, 'jac', LN)
    if nm in poolC:
        c = col(Vc, EC, kindC, nm, semBat, 'jac', LN); lost = ~np.isfinite(c)
        if lost.any() and not lost.all():
            vanish[nm] = list(labBat[lost]); c[lost] = 0
        Jc[:, j] = c
for nm, l in vanish.items(): print('  battery, %s x1.10: %d features vanish (%s); sensitivity set to 0' % (nm, len(l), ', '.join(l)))
m8p, m8t = ptATP == 8, trATP == 8
n3 = lambda J, n: J / np.sqrt(n[:, None] / 3)
tiers = [('A', 'Beard F-V points, 8/4/2 mM', Jp), ('B8', 'Beard records, 8 mM only', np.r_[Jp[m8p], Jt[m8t]]),
         ('B8n3', 'Beard records, 8 mM, as if n = 3', np.r_[n3(Jp[m8p], nPts[m8p]), n3(Jt[m8t], nTrc[m8t])]),
         ('B', 'Beard records, 8/4/2 mM', np.r_[Jp, Jt]), ('Bn3', 'Beard records, 8/4/2 mM, as if n = 3', np.r_[n3(Jp, nPts), n3(Jt, nTrc)]),
         ('C', 'slack battery, 8 mM', Jc), ('D', 'Beard records + slack battery', np.r_[Jp, Jt, Jc])]
# Where does the information sit? Bnr drops the redevelopment samples; Bok and Cok
# keep only the elements the incumbent reproduces within 3 standard errors.
zP = (Vp[0] - dPts) / semPts; zTr = (Vt[0] - dTrc) / semTrc; zBat = (Vc[0] - dBat) / semBat
JB = np.r_[Jp, Jt]; kindEl = np.r_[np.array(['point'] * 23, dtype=object), trKind]
okzB = np.abs(np.r_[zP, zTr]) <= 3; okzC = np.abs(zBat) <= 3
tiers += [('Bnr', 'Beard records without redevelopment', JB[kindEl != 'redev']),
          ('Bok', 'Beard records, elements the model fits', JB[okzB]), ('Cok', 'slack battery, elements the model fits', Jc[okzC])]
bad = [names[j] for j in range(P_) if not all(np.isfinite(J[:, j]).all() for _, _, J in tiers)]
assert not bad, 'missing Jacobian columns: %s' % bad
keys = [k for k, _, _ in tiers]; grp = [param_group(n) for n in names]

# ---------------------------------------------------------------- 3. Fisher information, verdicts
SP = {}; V = {}
print('\nTier                                        elements  identified  sloppy  unidentified  stiff directions  log10 lambda max / median')
for k, label, J in tiers:
    F = J.T @ J; ev = np.sort(np.linalg.eigvalsh((F + F.T) / 2))[::-1]
    with np.errstate(divide='ignore'): sFix = 1 / np.sqrt(np.diag(F))
    C = np.linalg.inv(F + np.eye(P_) / SIG_PRIOR ** 2); sFree = np.sqrt(np.diag(C))
    cls = np.array(['unidentified'] * P_, dtype=object); cls[sFix <= THR] = 'sloppy'; cls[sFree <= THR] = 'identified'
    if k in ('C', 'Cok'): cls[names.index('pAtp_k2')] = 'not in tier'
    SP[k] = dict(F=F, J=J, C=C, ev=ev); V[k] = dict(sFix=sFix, sFree=sFree, cls=cls)
    print('%-5s %-37s %6d %10d %8d %12d %15d %12.1f / %.1f' % (k, label, J.shape[0], (cls == 'identified').sum(), (cls == 'sloppy').sum(),
          (cls == 'unidentified').sum(), (ev > 1 / THR ** 2).sum(), np.log10(ev[0]), np.log10(max(np.median(ev), 1e-12))))
print('\nScale-free size of the constrained space: eigen-directions within 4 / 5 / 6 decades of the stiffest one')
for k, label, _ in tiers:
    r = SP[k]['ev'] / SP[k]['ev'][0]; print('   %-5s %-40s %3d / %3d / %3d' % (k, label, (r > 1e-4).sum(), (r > 1e-5).sum(), (r > 1e-6).sum()))
tot = (JB ** 2).sum()
print('Share of the Beard-record information (trace of F) by element kind: ' + ', '.join('%s %.0f %% (%d elements, %d within 3 sem of the data)' % (
    kd, 100 * (JB[kindEl == kd] ** 2).sum() / tot, (kindEl == kd).sum(), (okzB & (kindEl == kd)).sum()) for kd in ['point', 'ramp', 'redev', 'hold110']))
print('   by [ATP]: ' + ', '.join('%d mM %.0f %%' % (a, 100 * (JB[np.r_[ptATP, trATP] == a] ** 2).sum() / tot) for a in [8, 4, 2]))
print('   median over parameters of the redevelopment share of F_jj: %.0f %%' % (100 * np.median((JB[kindEl == 'redev'] ** 2).sum(0) / np.maximum((JB ** 2).sum(0), 1e-300))))
print('   slack battery: %d of %d elements within 3 sem of the data' % (okzC.sum(), okzC.size))
order = sorted(range(P_), key=lambda j: (grp[j], j))
print('\nPer-parameter uncertainty with all others free (log units; 0.69 = factor of two, 2.30 = the prior = no information). * = identified')
print('%-34s %-24s' % ('param', 'group') + ''.join('%9s' % k for k in keys))
for j in order:
    print('%-34s %-24s' % (names[j], grp[j]) + ''.join('%8.3f%s' % (V[k]['sFree'][j], '*' if V[k]['cls'][j] == 'identified' else ' ') for k in keys))
print('\nIdentified parameters by mechanism group:')
print('%-26s %5s' % ('group', 'size') + ''.join('%6s' % k for k in keys))
for g in sorted(set(grp)):
    m = np.array([x == g for x in grp])
    print('%-26s %5d' % (g, m.sum()) + ''.join('%6d' % (V[k]['cls'][m] == 'identified').sum() for k in keys))
ident = lambda k: set(n for n, c in zip(names, V[k]['cls']) if c == 'identified')
for a, b in [('C', 'B'), ('B', 'C'), ('C', 'B8'), ('B8', 'C'), ('B', 'B8'), ('B', 'A'), ('C', 'Bok'), ('Bok', 'C'), ('C', 'Bnr'), ('B', 'Bnr')]:
    print('Identified by %s but not by %s: %s' % (a, b, ', '.join(sorted(ident(a) - ident(b), key=names.index)) or '(none)'))
print('Identified only by the two together (D, neither B nor C): %s' % (', '.join(sorted(ident('D') - ident('B') - ident('C'), key=names.index)) or '(none)'))
for c in ['sloppy', 'unidentified']:
    print('%s under both B and C: %s' % (c, ', '.join(n for j, n in enumerate(names) if V['B']['cls'][j] == c and V['C']['cls'][j] == c)))
jA = names.index('pAtp_k2')
print('pAtp_k2 (the only ATP parameter): ' + ', '.join('%s sFix %.3f sFree %.3f (%s)' % (k, V[k]['sFix'][jA], V[k]['sFree'][jA], V[k]['cls'][jA]) for k in ['A', 'B', 'Bn3', 'D']))

# ---------------------------------------------------------------- 4. step and side check
print('\nStep check: Jacobian column from x1.05 and from x1/1.10 against the one used (x1.10). ratio = norm / norm used, cos = direction agreement.')
print('%-5s %-34s %9s %10s %8s %11s %9s %5s' % ('tier', 'param', 'normJ', 'ratioHalf', 'cosHalf', 'ratioMinus', 'cosMinus', 'lost'))
SC = []
for tk, Vs, se, E, kd in [('B', np.c_[Vp, Vt], np.r_[semPts, semTrc], EB, kindB), ('C', Vc, semBat, EC, kindC)]:
    for e in np.flatnonzero(kd == 'half'):
        nm = str(E[e].name); j = col(Vs, E, kd, nm, se, 'jac', LN); jh = col(Vs, E, kd, nm, se, 'half', np.log(1.05)); jm = col(Vs, E, kd, nm, se, 'minus', -LN)
        ok3 = np.isfinite(j) & np.isfinite(jh) & np.isfinite(jm)
        if ok3.mean() < 0.5: print('%-5s %-34s  more than half of the features lost at one of the steps' % (tk, nm)); continue
        j, jh, jm = j[ok3], jh[ok3], jm[ok3]; cs = lambda a, b: a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-300)
        row = (tk, nm, np.linalg.norm(j), np.linalg.norm(jh) / np.linalg.norm(j), cs(j, jh), np.linalg.norm(jm) / np.linalg.norm(j), cs(j, jm), (~ok3).sum())
        SC.append(row); print('%-5s %-34s %9.1f %10.2f %8.3f %11.2f %9.3f %5d' % row)
for tk in ['B', 'C']:
    q = [r for r in SC if r[0] == tk and r[2] > 1]; a = lambda i: np.array([r[i] for r in q])
    print('  tier %s: %d parameters with |J| > 1. Half step: ratio median %.2f (%.2f-%.2f), cos median %.3f, cos < 0.9 in %d. Other side: ratio median %.2f (%.2f-%.2f), cos median %.3f, cos < 0.9 in %d' % (
        tk, len(q), np.median(a(3)), a(3).min(), a(3).max(), np.median(a(4)), (a(4) < 0.9).sum(), np.median(a(5)), a(5).min(), a(5).max(), np.median(a(6)), (a(6) < 0.9).sum()))
    print('     cos < 0.9 (either check): %s' % ', '.join(r[1] for r in q if min(r[4], r[6]) < 0.9))

# ---------------------------------------------------------------- 5. synthetic recovery
print('\nSynthetic recovery: 3 truths, every parameter moved at once by exp(0.10 N(0,1)); estimate = one regularised Gauss-Newton step.')
print('rmsError = noise-free error; slope = fraction of the true move recovered; rmsErrNoisy = expected error with noise at the data sem; predSd = Fisher prediction')
recov = {}
for tk, Vs, se, E, kd in [('A', Vp, semPts, EB, kindB), ('B', np.c_[Vp, Vt], np.r_[semPts, semTrc], EB, kindB), ('C', Vc, semBat, EC, kindC)]:
    J = SP[tk]['J']; dT, dH, dN = [], [], []
    for e in np.flatnonzero(kd == 'synth'):
        fe = np.isfinite(Vs[e]); print('  tier %s, truth %d: %d of %d features exist' % (tk, e + 1, fe.sum(), fe.size))
        if fe.mean() < 0.3: continue
        Jr = J[fe]; Fr = Jr.T @ Jr; Cr = np.linalg.inv(Fr + np.eye(P_) / SIG_PRIOR ** 2)
        tn = [str(x) for x in E[e].name]; tf = np.asarray(E[e].f, float)
        dtrue = np.array([np.log(tf[tn.index(n)]) if n in tn else 0.0 for n in names])
        rr = (Vs[e, fe] - Vs[0, fe]) / se[fe]; est = Cr @ (Jr.T @ rr)
        dT.append(dtrue); dH.append(est); dN.append(np.sqrt((est - dtrue) ** 2 + np.diag(Cr @ Fr @ Cr)))
        pred = Jr @ dtrue; top = np.argsort(-np.linalg.norm(Jr * dtrue, axis=0))[:3]
        print('      real feature move %.0f sem (rms), linear prediction %.0f sem, direction agreement cos %.2f; largest predicted contributions: %s' % (
            np.sqrt((rr ** 2).mean()), np.sqrt((pred ** 2).mean()), rr @ pred / (np.linalg.norm(rr) * np.linalg.norm(pred)),
            ', '.join('%s (%+.0f %%)' % (names[i], 100 * (np.exp(dtrue[i]) - 1)) for i in top)))
    if not dT: continue
    dT, dH, dN = np.array(dT).T, np.array(dH).T, np.array(dN).T; recov[tk] = (dT, dH)
    for c in ['identified', 'sloppy', 'unidentified']:
        m = V[tk]['cls'] == c
        if not m.any(): continue
        x, y = dT[m].ravel(), dH[m].ravel()
        print('  %-2s %-13s n=%2d  rmsTruth %.3f  rmsError %.3f  slope %5.2f  corr %5.2f  rmsErrNoisy %.3f  predSd %.3f' % (
            tk, c, m.sum(), np.sqrt((x ** 2).mean()), np.sqrt(((y - x) ** 2).mean()), (x @ y) / (x @ x), np.corrcoef(x, y)[0, 1],
            np.sqrt((dN[m] ** 2).mean()), np.sqrt(np.diag(SP[tk]['C'])[m].mean())))

print('\nSingle-parameter recovery: truth = ONE parameter moved (x1.05 and x1/1.10, simulations not used for the Jacobian); all %d parameters estimated.' % P_)
print('recovered = estimate / truth for the moved parameter (1 = exact); leak = rms estimate of the other parameters / |truth|')
print('%-5s %-34s %-13s %10s %8s %10s %8s' % ('tier', 'param', 'class', 'rec x1.05', 'leak', 'rec x0.91', 'leak'))
SR = []
for tk, Vs, se, E, kd in [('A', Vp, semPts, EB, kindB), ('B', np.c_[Vp, Vt], np.r_[semPts, semTrc], EB, kindB), ('C', Vc, semBat, EC, kindC)]:
    J = SP[tk]['J']
    for e in np.flatnonzero(kd == 'half'):
        nm = str(E[e].name); j = names.index(nm); out = []
        for k2 in ['half', 'minus']:
            i = [q for q in range(len(E)) if kd[q] == k2 and isname(E[q], nm)][0]
            fe = np.isfinite(Vs[i]); Jr = J[fe]; Cr = np.linalg.inv(Jr.T @ Jr + np.eye(P_) / SIG_PRIOR ** 2)
            est = Cr @ (Jr.T @ ((Vs[i, fe] - Vs[0, fe]) / se[fe])); tr = np.log(float(E[i].f))
            out += [est[j] / tr, np.sqrt((np.delete(est, j) ** 2).mean()) / abs(tr)]
        SR.append((tk, nm, V[tk]['cls'][j]) + tuple(out)); print('%-5s %-34s %-13s %10.2f %8.2f %10.2f %8.2f' % SR[-1])
for tk in ['A', 'B', 'C']:
    q = [r for r in SR if r[0] == tk and r[2] == 'identified']
    if not q: continue
    a = np.array([[r[3], r[5]] for r in q]); l = np.array([[r[4], r[6]] for r in q])
    print('  tier %s, %d identified parameters tested: recovered within 0.5-1.5 in %d of %d cases (median %.2f); leak median %.2f, max %.2f' % (
        tk, len(q), ((a > 0.5) & (a < 1.5)).sum(), a.size, np.median(a), np.median(l), l.max()))

# ---------------------------------------------------------------- 6. fit of the incumbent
zA = (Vp[0] - dPts) / semPts; zT = (Vt[0] - dTrc) / semTrc; zC = (Vc[0] - dBat) / semBat
print('\nFit of the incumbent, chi2 per element against the standard error of the mean:')
print('   A  F-V points      chi2/N = %6.1f (N = %d)' % ((zA ** 2).mean(), zA.size))
for a in [8, 4, 2]: print('        %d mM points   chi2/N = %6.1f' % (a, (zA[ptATP == a] ** 2).mean()))
print('   B  trace samples   chi2/N = %6.1f (N = %d)' % ((zT ** 2).mean(), zT.size))
for kd in ['ramp', 'redev', 'hold110']:
    for a in [8, 4, 2]:
        m = (trKind == kd) & (trATP == a); print('        %-8s %d mM  chi2/N = %6.1f (N = %d), median |model - data| = %.3f F_ref' % (kd, a, (zT[m] ** 2).mean(), m.sum(), np.median(np.abs(Vt[0, m] - dTrc[m]))))
print('   C  slack battery   chi2/N = %6.1f (N = %d)' % ((zC ** 2).mean(), zC.size))
print('   median sem: F-V points %.4f, trace samples %.4f (F_ref units); slack battery %.4f (relative)' % (np.median(semPts), np.median(semTrc), np.median(semBat)))

# ---------------------------------------------------------------- 7. profile plan (stage 2, needs MATLAB)
cand = [j for j in range(P_) if V['C']['cls'][j] == 'identified' and V['B']['cls'][j] != 'identified']
if len(cand) < 3: cand = [j for j in range(P_) if V['C']['cls'][j] == 'identified']
worst = sorted(cand, key=lambda j: -V['B']['sFree'][j])[:3]
FB = SP['B']['F'] + np.eye(P_) / SIG_PRIOR ** 2
print('\nProfile plan (stage 2): %s' % ', '.join(names[j] for j in worst))
for j in worst:
    o = [i for i in range(P_) if i != j]
    for dl in [-THR, THR]:
        d = np.zeros(P_); d[j] = dl; d[o] = -np.linalg.solve(FB[np.ix_(o, o)], FB[o, j] * dl); d = np.clip(d, -np.log(3), np.log(3))
        print('   %-30s x%.2f: predicted move (chi2 units), Beard records %.2f, slack battery %.1f; compensation clipped in %d parameters' % (
            names[j], np.exp(dl), d @ SP['B']['F'] @ d, d @ SP['C']['F'] @ d, (np.abs(d[o]) >= np.log(3) - 1e-9).sum()))

# ---------------------------------------------------------------- 8. figures, tables
cols = dict(A=(.85, .45, 0), B8=(.55, .55, .55), B8n3=(.8, .8, .8), B=(.75, 0, 0), Bn3=(.95, .6, .6), C=(.1, .35, .75), D=(0, 0, 0),
            Bnr=(.6, .3, .6), Bok=(.95, .4, .7), Cok=(.5, .7, .95))
fig, ax = plt.subplots(figsize=(7.4, 4.4))
for k, label, _ in tiers:
    ax.semilogy(np.arange(1, P_ + 1), np.maximum(SP[k]['ev'], 1e-6), '-o', color=cols[k], ms=3, lw=1.3, label='%s  %s' % (k, label))
ax.axhline(1 / THR ** 2, color='k', ls='--', lw=.8); ax.text(P_, 1 / THR ** 2 * 1.3, 'known to a factor of 2', ha='right', fontsize=8)
ax.axhline(1 / SIG_PRIOR ** 2, color='k', ls=':', lw=.8); ax.text(P_, 1 / SIG_PRIOR ** 2 * 0.35, 'the prior (one decade)', ha='right', fontsize=8)
ax.set(xlabel='parameter combination (eigen-direction), sorted', ylabel='Fisher information eigenvalue', ylim=(1e-4, 1e8), xlim=(1, P_),
       title='How many parameter combinations each data set constrains'); ax.legend(fontsize=7, loc='upper right'); ax.grid(alpha=.3)
fig.savefig(os.path.join(FIG, 'fvid_spectrum.png'), dpi=150, bbox_inches='tight'); plt.close(fig)

Mx = np.array([[V[k]['sFree'][j] for k in keys] for j in order]); Mx[[names[j] == 'pAtp_k2' for j in order], keys.index('C')] = np.nan
fig, ax = plt.subplots(figsize=(7.2, 13))
im = ax.imshow(np.log10(Mx), vmin=np.log10(0.02), vmax=np.log10(SIG_PRIOR), cmap='viridis', aspect='auto')
ax.set_xticks(range(len(keys))); ax.set_xticklabels(keys); ax.set_yticks(range(P_)); ax.set_yticklabels(['%s  [%s]' % (names[j], grp[j]) for j in order], fontsize=6)
for r, j in enumerate(order):
    for c, k in enumerate(keys):
        if V[k]['cls'][j] == 'identified': ax.plot(c, r, 'w.', ms=5)
cb = fig.colorbar(im, ax=ax, fraction=.05); cb.set_label('log10 uncertainty with all others free (log units)')
ax.set_title('Per-parameter uncertainty by data tier (dark = well known)\ndot = identified (within a factor of 2)', fontsize=9)
fig.savefig(os.path.join(FIG, 'fvid_parameters.png'), dpi=150, bbox_inches='tight'); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
cc = dict(identified=(0, .5, .2), sloppy=(.9, .6, 0), unidentified=(.7, .7, .7))
for ax, tk in zip(axs, ['B', 'C']):
    if tk in recov:
        dT, dH = recov[tk]
        for c in ['unidentified', 'sloppy', 'identified']:
            m = V[tk]['cls'] == c
            if m.any(): ax.plot(dT[m].ravel(), dH[m].ravel(), 'o', color=cc[c], ms=4, label=c)
    ax.plot([-.3, .3], [-.3, .3], 'k-', lw=.8); ax.axhline(0, color='k', ls=':', lw=.8)
    ax.set(xlabel='true change (log)', ylabel='recovered change (log)', xlim=(-.3, .3), ylim=(-.3, .3), title='Synthetic recovery, tier %s' % tk); ax.grid(alpha=.3); ax.legend(fontsize=8, loc='upper left')
fig.savefig(os.path.join(FIG, 'fvid_recovery.png'), dpi=150, bbox_inches='tight'); plt.close(fig)

with open(os.path.join(FIG, 'verdicts.csv'), 'w') as fh:
    fh.write('param,group,' + ','.join('sFix_%s,sFree_%s,class_%s' % (k, k, k) for k in keys) + '\n')
    for j in order:
        fh.write('%s,%s,' % (names[j], grp[j]) + ','.join('%.4g,%.4g,%s' % (V[k]['sFix'][j], V[k]['sFree'][j], V[k]['cls'][j]) for k in keys) + '\n')
np.savez(os.path.join(RES, 'fvid_py.npz'), names=np.array(names), Vp=Vp, Vt=Vt, Vc=Vc, Jp=Jp, Jt=Jt, Jc=Jc, ptATP=ptATP, ptV=ptV,
         trATP=trATP, trV=trV, trKind=trKind, trDt=trDt, dPts=dPts, semPts=semPts, dTrc=dTrc, semTrc=semTrc, dBat=dBat, semBat=semBat, labBat=labBat)
print('\nfigures -> %s' % FIG)
