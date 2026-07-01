# Low-ATP cross-bridge mechanism — consolidated report

*Cardiac sarcomere model (Beard 2022 lineage). Goal: explain the measured 2 mM-vs-8 mM MgATP
cross-bridge signature with physiologically-justifiable mechanisms.* Full logs:
[`labdiary.md`](labdiary.md), [`conclusions.md`](conclusions.md), and the sibling analyses
[`LowATP_k2Frontier`](../LowATP_k2Frontier/conclusions.md), [`LowATP_3rdState`](../LowATP_3rdState/labdiary.md).

---

## 1. The data — what low ATP does (03/27/2026, 2 mM vs 8 mM slack)

The 5 slack segments are **paired** between conditions (same protocol), so per-segment ratios are
valid. The signature:

- **Force amplitude flat at ×1.18** across every segment (steady, A, Am, peak2, vall) → a single
  multiplicative gain.
- **ktr slows with slack** (×0.64→0.46), **t0** rises with slack → kinetics are segment-dependent.
- restretchSlope / overshoot are noisy (de-weighted).

![ATP signature per slack segment](../LowATP_k2Frontier/results/atp_per_segment.png)

Matches Beard 2022 (+18% force, −37% ktr). Scoring throughout = **relative ratio**
`model(2mM)/model(8mM)` vs the data ratio (cancels baseline-fit imperfection).

## 2. Single-knob frontier — k2 owns kinetics but overshoots force

Reducing `k2` (detachment) reproduces the ktr slowdown but over-produces isometric force ~2×.
Mapping every lever in the force–ktr plane: **only k2 reaches the data region**, and every
force-trim either breaks ktr (occupancy / Pi) or breaks peak2 (ka). No single 2nd lever lands it.

![Lever map in the force–ktr plane](../LowATP_k2Frontier/results/lever_map.png)

## 3. A 3rd state breaks the restretch-transient wall

Adding a post-stroke `p3` state lets the model raise the **restretch transients (peak2, vall)**
that the 2-state could not — cost 0.70 → 0.39. (But that early config used an unphysical
low-stiffness p3 + ad-hoc `ka`; see §4.)

![2-state vs 3-state per feature](../LowATP_3rdState/results/thirdstate_compare.png)

## 4. Physiologically-justifiable mechanisms (the answer)

Low ATP at the myofibril = **↓ATP + ↑ADP + ↑Pi simultaneously**. Three nucleotide effects were
wired into the active model path (previously inert), behind flags, baseline-preserving at
MgADP=0/Pi=0, with physiological stiffness (rigor `kstiff3=kstiff2`, no `ka` lever):

| mechanism | biochem | lever |
|---|---|---|
| **ADP-trap** (`UseAdpTrap`) | ↑ADP rebinds, blocks ADP release → traps stiff P2 | R2·=g2 |
| **Pi reversal** (`UsePiReversal`) | ↑Pi is a power-stroke product → inhibits forward stroke (lowers force) | R12·=f2, R21·=(1+Pi/K_Pi) |
| **Rigor** (`UseAtpDetach`) | ↓ATP slows rigor (A·M) detachment → rigor accumulates | R3·=MgATP/(MgATP+K_T_detach) |

**Result — the coupled state reproduces the full force signature from concentrations alone**
(MgATP 2, **[ADP]≈1.5 mM, [Pi]≈3 mM**; cost **0.555**): steady/A/Am/**peak2/vall** all +16–24%.
Division of labour: ADP-trap → isometric force; rigor → restretch transients; **Pi is essential**,
tempering the ADP+rigor overshoot (1.41→1.16). Single mechanisms each capture only half.

![Mechanisms vs data](results/mechanisms_compare.png)

## 5. Raw fits — what matches and what doesn't

Model (blue) vs data (black) force-time traces. **Steady plateaus match well** (the amplitudes the
mechanism explains); the **restretch peaks overshoot** (sharp spikes).

![8 mM raw trace](results/raw_trace_8mM.png)
![2 mM coupled raw trace](results/raw_trace_2mM_coupled.png)

## 6. The restretch frontier (current work)

The peak overshoot is **pre-existing** (the 2-state best already overshoots peak1 ~10%) and a
**baseline trade-off**. Diagnosis: `restretchSlopeStart` (the kSE proxy) is too low (model ~1280 vs
data 1588) because **kSE is under-constrained** — it is frozen during the `ka=0` passive fit. But:

- **kSE↑** to match the slope (≈4800) **overshoots peak1** and **speeds ktr** — see below.
- **kstiff2↓** lowers peak1 but drops `steady` (it's the force scale).
- **c_SE_visc** damps peak1 but crashes the slope.
- **`UseA2AttachmentShift` AMPLIFIES** (it re-binds strained heads at a lower-strain site → traps
  force), so it is *not* a damper.

![Restretch trade-off (kSE vs peak1)](results/restretch_tradeoff.png)

→ The restretch peak is a genuine **4-way knot** (slope↑kSE / peak1↓kstiff2 / steady↑kstiff2 /
ktr). It needs a **full bounded re-optimization** (free `ka/kd/k1/k2` to compensate) — the
`Analyses/RestretchMechanisms` Variant-A recipe — not a single knob. Underneath sits the deeper,
pre-existing **ktr-too-fast** residual (model 66 vs data 49).

## 7. Side quest — classic vs slack ktr

The classic ktr (hold 2.0 µm; slack→stretch→release) runs ~5 s⁻¹ **slower** than the slack ktr
(2.2 µm → slack), both still 1st-order. **Explanation: SRX.** The two protocols start from different
SL/force histories → different SRX occupancy at redevelopment onset (higher force at 2.2 µm
mobilizes more heads *out* of SRX → larger cycling pool → faster ktr). Because SRX mobilization is
**much faster than ktr** (timescale separation), the different starting pool just **rescales the
single-exponential rate** — it equilibrates in the first few ms and never appears as a resolvable
2nd time constant, so 1st-order kinetics are preserved. No fine-balancing of SRX timing is needed.
(The 03/27 dataset contains both `8mM_stiff_ktr` and slack; the current config has `kmsr=0` so SRX
mobilization must be enabled to test this in-model.) Useful as an extra, SL-history-dependent
constraint on SRX.

## 8. Status & next steps

- ✅ Run infrastructure fixed (`loadParams` sandbox, `refreshPool`).
- ✅ Nucleotide mechanisms wired + validated; **coupled ↓ATP+↑ADP+↑Pi is the settled answer** for
  the force-amplitude signature (cost 0.555, physiological concentrations).
- 🔜 **Restretch re-fit (Part B-full):** bounded multi-parameter optimization to resolve the
  peak1/slope/steady/ktr knot, then re-validate the coupled mechanism on the clean baseline.
- 🔜 **ktr side quest (Part D):** extract classic vs slack ktr from the dataset; enable SRX
  mobilization and confirm the rescaling explanation in-model.
