# Low-ATP cross-bridge modelling — state of the investigation

*Synthesis across LowATP_k2Frontier, LowATP_3rdState, LowATP_Mechanisms, RestretchMechanisms,
PassiveForceID, and the SRX/mech-tradeoff notes. Figures: this folder's `results/` +
`../LowATP_k2Frontier/results/`.*

---

## 1. Where we are on the **8 mM (high-ATP) fit**

The current best (`params_reseeded_regavail_opt2`) reproduces the **force amplitudes** well — steady,
A, Am all match the 8 mM data; FV is reasonable. Four **residuals persist at 8 mM**, and they are the
*same* residuals that block low ATP:

| 8 mM residual | model vs data | character |
|---|---|---|
| **ktr too fast** | 66 vs 49 s⁻¹ | the Brenner "iron law" (ktr ↔ duty ↔ force tightly coupled) |
| **restretch peak overshoot** | peak1 114 vs 89; peak2 ≫ steady (the post-restretch bump) | mean-field synchrony (see §4) |
| **restretchSlope too low** | 1280 vs 1588 | kSE under-constrained (frozen in the `ka=0` passive fit) |
| **FV shoulder** | fnorm too steep at low v | sub-maximal attachment availability |

So: **amplitudes fit; kinetics (ktr) and the restretch-transient *shape* do not.**

## 2. Where we are on the **2 mM (low-ATP) fit** — the win

Low ATP at the myofibril = **↓ATP + ↑ADP + ↑Pi simultaneously**. The coupled nucleotide mechanism
(all physiologically wired, no `ka`, rigor `kstiff3=kstiff2`) **reproduces the full force-amplitude
signature from concentrations alone** (MgATP 2, [ADP]≈1.5 mM, [Pi]≈3 mM; relative-ratio cost 0.555):
steady/A/Am/**peak2/vall** all +16–24%. Division of labour: **ADP-trap → isometric force** (traps
stiff P2); **rigor → restretch transients** (stiff P3 accumulates); **Pi is essential** — it tempers
the ADP+rigor overshoot (1.41→1.16).

![mechanisms vs data](results/mechanisms_compare.png)

## 3. What **prevents the low-ATP fit** — the residuals (all kinetic/timing)

| residual | model | data | why |
|---|---|---|---|
| **ktr level** | ratio ~0.64 (flat) | 0.54 | under-slowed; the iron law |
| **ktr SLOPE** | flat | 0.64→0.46 with slack | no SL/strain-dependence in the ATP effect |
| **t0 (onset)** | ratio ~2.1 | 1.31 | rigor over-delays onset; **decoupled from kSE** (see fig) |
| **restretch peak** | 160 vs 95 | inherited 8 mM overshoot, rigor-amplified | structural (§4) |
| **restretchSlope** | ratio 1.10 | 1.26 | rigor under-stiffens the series path |

![low-ATP kSE vs t0](results/lowatp_kSE_t0.png)

**The force is solved; the kinetics are not.** Every kinetic residual traces to one of two structural
issues: the **ktr↔force iron law** and the **mean-field restretch synchrony**.

## 4. Every path we explored (the map)

| # | path | verdict |
|---|---|---|
| 1 | **k2 alone** (detachment) | ktr ✓ but force overshoots ~2× — single knob can't do both |
| 2 | occupancy cap | left-UP in force–ktr plane — breaks ktr. ✗ |
| 3 | **SRX** (de/stabilization) | null at max-Ca (mechanosensitively suppressed pool). ✗ for force; useful for the ktr side-quest |
| 4 | Pi alone | lowers force — a *tempering* term, not a driver |
| 5 | **ka↓** | down-left but **crashes peak2**, and **not ATP-justifiable** — rejected on physiology |
| 6 | 3rd state, *wobbly* p3 (kstiff3<kstiff2) + ka | broke the restretch wall (0.39) but **unphysical** — rejected |
| 7 | **coupled rigor+ADP-trap+Pi** (physiological) | **force amplitudes ✓ (0.555)** — the answer |
| 8a | restretch: **kSE↑** | matches slope but raises peak1 + speeds ktr |
| 8b | restretch: **kstiff2↓** | lowers peak1 but drops steady (it's the force scale) |
| 8c | restretch: **c_SE_visc** | damps peak1 but crashes restretchSlope |
| 8d | restretch: **A2-shift** | **amplifies** (re-binds strained heads → traps force), not a damper |
| 8e | restretch: **R2 strain-dep** (steepen high-+s) | drops peak1 but **crashes steady + explodes ktr** (cycle coupling; restretch strain overlaps isometric strain) |
| 9 | ktr side-quest (classic vs slack) | SRX timescale-separation explains classic<slack, 1st-order preserved |
| 10 | *scoped / not finished* | full bounded restretch re-optimization; catch-bond; lattice-SL ktr-slope; classic-ktr extraction; tuned viscoelastic damping |

**Two walls fell out of this map:** (a) **no single lever lands force+ktr** (resolved by the coupled
nucleotide mechanism for *force*), and (b) **no single lever fixes the restretch peak** without
breaking steady/ktr — it is structural.

## 5. What this is telling us — the (mildly dogma-breaking) hypothesis

The residuals are not random; they converge on **two assumptions of the mean-field cross-bridge
model that the low-ATP data appears to violate**:

**(a) The force ↔ ktr "iron law" is too rigid.** In a simple cycling model, force ∝ duty ratio and
ktr = f_app+g_app, so you cannot raise force +18% while slowing ktr −46% with a single rate. The
*only* way we reproduced it was to put force into states **outside the normal cycle** — a **rigor
(nucleotide-free) and an ADP-trapped state that bear force but do not cycle**. The implication is
physiologically pointed: **a substantial fraction of low-ATP force is "static" (non-cycling) rather
than actively cycling** — against the common reading that active cardiac force is carried by actively
cycling bridges. This *is* the heart-failure phenotype: extra force (↑ diastolic stiffness) with
slowed/weaker cycling (↓ power) — the force–power paradox, here emerging mechanistically from
rigor/ADP-trapped force-bearing states under energetic stress.

**(b) The restretch is over-synchronized.** The model's attached population responds to the imposed
stretch *in unison* → a sharp elastic spike + a ringing reattachment second peak. The data shows a
small spike + a smooth recovery. No rate knob removes the spike without breaking the cycle (§4, 8e),
because the mean field has no internal heterogeneity to dephase the response. The real signal — the
**ktr that slopes with slack** and the **smooth restretch** — both point to a *distributed /
asynchronous* attached population that the mean-field single-strain-grid model cannot represent.

These are testable, not hand-waving: (a) predicts the low-ATP extra force is rigor-like (stiff,
slow-detaching, ATP-rescuable) — measurable as a stiffness/force ratio and an ADP/Pi dependence;
(b) predicts the restretch shape and the ktr-slope are governed by strain heterogeneity, not by any
single rate.

## 6. Physiologically-justified next steps (priority order)

1. **Make the rigor/ADP effect SL/strain-dependent** → the ktr **slope** (0.64→0.46). At larger slack,
   more shortening → more rigor accumulation → more slowdown. Concrete: gate `k3`/ADP-trap by the
   shortening extent or strain distribution. *(Targets the #1 unfit feature.)*
2. **Decouple onset (t0) from rigor depth** → a parallel fast-attachment / mobilization path so that
   rigor adds steady-state force without delaying the *onset*. *(t0 is currently rigor-welded.)*
3. **Address restretch synchrony physiologically, not by rate tuning** — a **velocity-dependent
   slip-bond** that detaches the *synchronized* over-strained population during the ramp (removes the
   spike without touching isometric detachment), and/or a **distributed reattachment** (dephase the
   second peak). Both are more faithful than the mean-field instantaneous response.
4. **Quantify the "static force" claim** — run the rigor model across an ADP/Pi grid and report the
   non-cycling force fraction + stiffness; compare to the measured restretchSlope/force ratio (data
   1.26) and to the classic-vs-slack ktr (SRX). This is the experiment that would *substantiate the
   dogma-challenge*.
5. **Bounded multi-feature re-optimization** of the 8 mM baseline (kSE, kstiff2, ka/kd/k1/k2, the R2
   curve) to remove the pre-existing residuals so the low-ATP kinetics are judged from a clean
   reference — then re-validate the coupled mechanism on it.

**Bottom line:** the *force* side of low ATP is physiologically solved (coupled ↓ATP+↑ADP+↑Pi via
rigor/ADP-trap/Pi). The *kinetic/timing* side (ktr level+slope, t0, restretch shape) is blocked by two
structural assumptions — the rigid force↔ktr coupling and the mean-field synchrony — and the way to
break them is **force-bearing non-cycling states + strain heterogeneity**, both of which are
physiologically real and, taken together, reframe low-ATP/failing-heart contractility as **partly
static (rigor-like) rather than purely cycling.**
