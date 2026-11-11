"""Simple 4-fit figure: same twitch from different Ca / cross-bridge assumptions, with the well range."""
import os, sys
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
res = sys.argv[1] if len(sys.argv) > 1 else r"C:\home\git\ATP-wt-ca\Analyses\CalciumTwitch\results"
m = sio.loadmat(os.path.join(res, "eht_nonunique.mat"), squeeze_me=True, struct_as_record=False)
D = sio.loadmat(os.path.join(res, "eht_controls.mat"), squeeze_me=True, struct_as_record=False)["D"]

pick = [("fit C", "Reference fit", "#222222"),
        ("core rates x0.5", "Cross-bridges 2× slower", "#2B6CB0"),
        ("core rates x2", "Cross-bridges 2× faster", "#C05621"),
        ("k_off fixed x2", "Different Ca handling", "#2F855A")]
names = [a.name for a in m["alts"]]
R = {n: m["R"][names.index(n)] for n, _, _ in pick}
t = m["tfine"] * 1000

# well-to-well range of the 15 control wells (beat-averaged per well)
ctrl = D.isCtrl.astype(bool)
Fw = D.Fw[:, ctrl]
tb = D.tb * 1000
lo, hi = Fw.min(axis=1), Fw.max(axis=1)

plt.rcParams.update({"font.size": 12})
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), dpi=200, gridspec_kw=dict(width_ratios=[1.25, 1, 0.9]))

ax = axs[0]
ax.fill_between(tb, lo, hi, color="#BBBBBB", alpha=0.6, lw=0, label="measured: range of 15 tissues")
ax.plot(tb, D.meanF, color="#777777", lw=1.5, ls="--", label="measured: mean")
for n, lab, c in pick:
    ax.plot(t, R[n].sigma * R[n].fm, color=c, lw=2, label=lab)
ax.set(xlim=(0, 900), xlabel="time after stimulus (ms)", ylabel="twitch force (µN)",
       title="1. All four fits reproduce the twitch")
ax.legend(fontsize=9.5, frameon=False)

ax = axs[1]
for n, lab, c in pick:
    ca = np.asarray(R[n].S.Ca, float)
    ca = (ca - ca.min()) / (ca.max() - ca.min())
    ax.plot(t, ca, color=c, lw=2)
ax.set(xlim=(0, 900), ylim=(0, 1.05), xlabel="time after stimulus (ms)", ylabel="Ca²⁺ (normalized)",
       title="2. …each needs a different Ca²⁺ transient")
ax.text(880, 0.95, "not measured", ha="right", color="#C53030", fontsize=11, style="italic")

ax = axs[2]
ref = R["fit C"].ATPsig
vals = [R[n].ATPsig / ref for n, _, _ in pick]
y = np.arange(len(pick))[::-1]
ax.barh(y, vals, color=[c for _, _, c in pick], height=0.6)
for yi, v in zip(y, vals):
    ax.text(v + 0.04, yi, f"{v:.1f}×", va="center", fontsize=12)
ax.set_yticks(y); ax.set_yticklabels([lab for _, lab, _ in pick])
ax.axvline(1, color="#999", lw=1, ls=":")
ax.set(xlim=(0, 2.4), xlabel="relative to reference fit", title="3. …and implies different\nenergy use (ATP per beat)")

for a in axs:
    for s in ["top", "right"]: a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(here, "fits_simple.png"))
print([f"{lab}: rms {R[n].rms*100:.2f}%  ATP {R[n].ATPsig/ref:.2f}" for n, lab, _ in pick])
