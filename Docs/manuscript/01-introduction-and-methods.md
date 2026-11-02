# A mechanistic, ATP-dependent cross-bridge model of cardiac muscle mechanics under energetic stress

**Draft — Introduction and Methods.**
Markdown with inline LaTeX, intended to be pasted into Overleaf. Math uses `$…$` / `$$…$$`;
citations use `\citep{}`/`\citet{}` against the accompanying `references.bib`.

> ## Revision status (2026-11-02)
>
> §1 Introduction replaced by the Overleaf version (main.tex @ 51e26da); `references.bib`
> replaced by the Overleaf bib (a superset). The old §1.1–1.6 subsections are gone, so the
> §1.x cross-references below and in §2–3 are stale; the old §1.6 "Principal findings"
> numbers belong in Results now. Git history (deb4103) has the previous text.
>
> ## Revision status (2026-10-07)
>
> The 2026-09-18 review findings ([`REVIEW-corrections.md`](REVIEW-corrections.md)) have been
> worked into this file. Decisions taken by the authors:
>
> 1. **D1 — headline.** The paper claims that one model reproduces *every measured aspect* of
>    the 8 mM and 2 mM behaviour, approximately. Both `ktr` probes are measurements and both
>    are reported; neither carries a mechanistic claim, and the "two gates are required"
>    argument is withdrawn (§1.6).
> 2. **D2 — transport.** Pi build-up (and ATP depletion) in the trabecula core is treated as a
>    quantified limitation with its own analysis, `Analyses/PiAccumulation/plan.md` (§3.1).
> 3. **D3 — Beard 2022.** This is the follow-up to \citet{Beard2022} (D.A.B. is a PI here). The
>    contribution is to quantify, with one parameter set, the ATP dependence that the 2022
>    model could not reproduce (§1.3).
> 4. **New:** thin-filament Ca²⁺ activation after \citet{Campbell2018} is described in §2.3.12
>    (implemented on branch `CaActivation`, not used in the fits reported here).
> 5. **New:** §2.3.2–2.3.3 now state the simulated physics (an advected density in strain), not
>    only the code's moving-frame bookkeeping.
>
> Still open: bench details ([TODO]s), the per-feature unfloored χ² table, parameter
> uncertainty, and the final low-ATP fit (§1.6 numbers are from the 2026-09-26 sweep).

> **Rendering convention — please keep to it when editing.** Display-math `$$` delimiters
> sit alone on their own line (opening and closing), and an inline `$…$` never spans a line
> break. Most Markdown renderers, GitHub's included, mis-pair the delimiters otherwise and
> everything downstream of the offending line renders as garbage.
Items needing a decision or a number from the bench are marked **[TODO]**.

---

## 1. Introduction

<!-- Synced verbatim from Overleaf ("ATP Dependent Crossbridge Model", main.tex, commit 51e26da,
2026-11-02). Overleaf is the source of truth for this section: edit there, then re-sync.

Introduction outline (working notes, remove before submission)
1. ATP is reduced in the failing heart, and this may limit power.
   - Myocardial ATP falls from about 8 mM to 3-4 mM in advanced failure as the
     adenine nucleotide pool is depleted (Ingwall2002, Beard2022); also in TAC
     rats (Lopez2020).
   - Heart failure as a limitation of cardiac power output (BeardHummelJezek2021).
   - Low MgATP lowers myofilament power by 20-25 % and ktr by ~20 %, with Fmax
     and Vmax unchanged (Beard2022).
   - Existing mechano-energetic models cannot simulate this: their cross-bridge
     cycle is ATP-saturated, so energetic failure acts only through Pi
     (Lopez2020, Tewari2016).
2. No model reproduces the measured ATP sensitivity.
   - Tewari cycle + Campbell SRX recruitment fails with one parameter set, even
     for arbitrarily large K_MgATP (Beard2022, Tewari2016, Campbell2018).
   - Other models give the direction only, or are untested in the 2-8 mM range
     (Tran2010, Schmidt2024, Wang2026).
3. Two candidate mechanisms (our postulates).
   - H1, ATP-limited detachment. The post-stroke head detaches only after
     binding MgATP. Apparent Km of whole-cycle outputs is 10-150 uM (Beard2022),
     but the dissociation constant of ATP from the attached head is 0.1-1 mM by
     sinusoidal analysis (Kawai1993, Fujita2002, Kawai2024). Not a paradox; a
     discussion point.
   - H2, ATP-dependent head availability. ATP turns thick filaments ON
     (Rhodehamel2025); nucleotide-dependent availability shown for dATP
     (Teitgen2024); SRX recruitment in the parent model (Campbell2018).
   - Both proposed but not tested by Beard2022.
4. Separating H1 and H2 needs non-physiological protocols.
   - Force-velocity and a single ktr cannot separate them (Beard2022).
   - Slack/re-stretch (recruitment), re-stretch (strain-dependent detachment),
     stretch-hold staircase (thick-filament loading Linari2015, titin
     viscoelasticity Jezek2025).
   - Paired 8 and 2 mM recordings in the same trabecula.
5. Aims, hypotheses, expected result.
   - Aim 1: one parameter set reproduces all protocols at 8 mM.
   - Aim 2: changing only ATP-dependent steps reproduces 2 mM.
   - Expected: H1 dominant, H2 needed for kinetics; model can then be coupled
     to energetics models (Bazil2016, Collins2026).
-->

The pumping action of the heart is, at the molecular level, the summed output of myosin cross-bridges,
each transducing the free energy of ATP hydrolysis into force and filament sliding \citep{hill1974}.
Heart failure can be regarded as a limitation of cardiac power output
\citep{BeardHummelJezek2021}, and the failing myocardium is simultaneously an
energy-depleted one: myocardial ATP falls from approximately 8~mM in the healthy heart to
3--4~mM in advanced heart failure as the adenine nucleotide pool is depleted
\citep{Ingwall2002, Wu2009, lopez2024}, with comparable depletion observed in rats with decompensated
pressure-overload hypertrophy \citep{Lopez2020}. Whether this depletion impairs contraction
directly---by slowing cross-bridge cycling itself---or only indirectly, through the
accompanying accumulation of hydrolysis products, remains unresolved. Here we address the
ATP dependence of cross-bridge kinetics across the millimolar concentration range spanned
by the healthy and failing heart.

The cross-bridge cycle couples the hydrolysis of ATP to the production of mechanical
work: (i)~hydrolysis of ATP on the detached myosin head primes the lever arm; (ii)~the
primed head binds actin; (iii)~phosphate release drives the working stroke, generating
force and filament sliding; and (iv)~ADP dissociation followed by ATP binding dissociates
myosin from actin. Cooke and Pate~\cite{pate1989} introduced a model of cross-bridge kinetics and force
generation in which these state transitions depend explicitly on the concentrations of ATP,
ADP and inorganic phosphate (P$_{\rm i}$), showing that maximum unloaded
shortening velocity becomes ATP-limited only when ATP falls into the micromolar range. The
corollary---that the cycle is effectively saturated at the millimolar concentrations
prevailing in vivo---is supported by measurements in permeabilized mouse trabeculae, in
which maximal force and maximal shortening velocity were insensitive to changes in MgATP
over the millimolar range \citep{Beard2022}. Current models of myocardial
mechano-energetics embody this view, and in the analysis of Lopez et al.~\cite{Lopez2020} energetic
failure impaired contraction only through a concurrent rise in P$_{\rm i}$.

Maximal force and velocity, however, are not the quantities that fail in the failing heart.
In demembranated trabeculae, lowering MgATP from 8 to 4~mM reduces peak power by
20--25\,\% and the rate of tension redevelopment ($k_\mathrm{tr}$) by approximately 20\,\%,
while leaving maximal force and unloaded shortening velocity unchanged \citep{Beard2022}.
ATP depletion of the magnitude observed in heart failure may therefore contribute directly
to contractile failure. Yet existing models cannot test this proposition, because their
cross-bridge cycles are near-saturated with MgATP at millimolar concentrations. Beard et al.~\cite{Beard2022}
were unable to reproduce the observed MgATP dependence using a strain-dependent cross-bridge
model combining the cycle of Tewari et al.~\cite{Tewari2016} with force-dependent recruitment of myosin
from the super-relaxed (SRX) state \citep{Campbell2018}: a single parameter set captured the
force--velocity relation at each concentration but not its change with MgATP, and increasing
the MgATP dissociation constant without bound did not resolve the discrepancy. Other
metabolite-sensitive formulations reproduce the direction of the MgATP effect on force and
$k_\mathrm{tr}$ \citep{Tran2010} or derive ATP-dependent rates from thermodynamic principles
\citep{Schmidt2024, Wang2026}, but none have been tested against data on muscle function at ATP concentration in the rage of 2 to
8~mM. The question is therefore what feature of cross-bridge cycling confers sensitivity to
MgATP in a range where the cycle is conventionally held to be saturated.

Two mechanisms that could make the cross-bridge cycle sensitive to ATP concentration in this range
were proposed, but not tested, by \cite{Beard2022}. First, a post-stroke head detaches from
actin only after binding MgATP. The apparent $K_m$ of 10--150~$\mu$M for ATPase activity, force
and shortening velocity \citep{Beard2022} is a property of the whole cycle, whereas the
dissociation constant of MgATP from the attached head, measured by sinusoidal analysis, is
0.1--1~mM in cardiac muscle \citep{Kawai1993, Fujita2002, Kawai2024}. At the upper end of
this range, detachment would be reduced by about a quarter at and ATP concentration of 2~mM compared to 8~mM. Second, ATP shifts
myosin heads from the ordered OFF state on the thick filament to the disordered ON state
\citep{Rhodehamel2025}, so that a lower ATP may reduce the number of heads available to cycle,
as proposed to be an effect of 2$'$-deoxy-ATP \citep{Teitgen2024}.

Measurements of force-velocity relation and a standard $k_\mathrm{tr}$ estimates cannot separate these mechanisms,
because slower detachment and fewer available heads both lower power and slow redevelopment.
However, we found that assessment of force during a sequence of large rapid changes in length can potentially differentiate these potential mechanisms.
Specifically, we designed a sequential slack and re-stretch protocol to assess sarcomere length-dependent force development and maximal unloaded velocity following slack, strain dependence of ratcheting and detachment transitions and stretch--hold steps loading the thick filament \citep{Linari2015}.
The restretch part might be affected by the viscoelasticity of the titin, which has been found to be significantly augmented at high-Ca conditions, while the active contraction was completely suppressed  \citep{Jezek2025, Baker2025-rz} . The combination of force-velocity, slack-restretch, $k_\mathrm{tr}$ and high-Ca passive protocols at high (8mM) and low (2mM) ATP concentration present an information-rich dataset, suitable for identification of complex model with ATP sensitivity.

We used the data from these experiments to identify a strain-dependent cross-bridge model that matches data from all protocols at 2, 4, and 8~mM ATP with the goal of finding the simplest explanation for the dependency of muscle dynamics on ATP concentration. In practice we started by first identifying a model that reproduces all protocols at 8~mM MgATP with one parameter set. The next step was to test whether adjusting only the ATP-dependent
steps of that model reproduces the responses at 2~mM. We tested two non-exclusive hypotheses: that low ATP slows detachment of the post-stroke state, and that low ATP reduces the availability of heads through the SRX pool. We found that both hypotheses were necessary to explain the observed data on muscle dynamics as a function of MgATP concentration.

---

## 2. Methods

### 2.1 Experimental preparation and protocols

#### 2.1.1 Preparation

Demembranated right-ventricular trabeculae were prepared from adult (about 12-week-old) male
and female C57BL/6J mice (Jackson Laboratory) as described previously
\citep{Beard2022, Baker2025-rz}. The study was approved by the Animal Care and Use
Subcommittee of the San Francisco Veterans Affairs Medical Center. In brief, trabeculae were
demembranated in ice-cold relaxing solution containing 2\,\% Triton X-100 for 1 h, washed,
and stored at $-20\,^\circ$C in a 1:1 mixture of relaxing solution and glycerol. Relaxing
solution contained (in mM) 20 EGTA, 8, 4 or 2 MgATP, 12 creatine phosphate and 100 BES, pH 7.1
(KOH), ionic strength 200 mM (KCl); in the preactivating solution 19.5 mM EGTA was replaced
by HDTA, and the activating solution contained 20 mM Ca$^{2+}$-EGTA (pCa 4.51). All solutions
contained 10 IU/mL creatine kinase and 1\,\% (v/v) protease inhibitor cocktail (P-8340,
Sigma). Myosin was inhibited, where stated, with 50 $\mu$M para-nitroblebbistatin (PNB) plus
50 $\mu$M mavacamten \citep{Baker2025-rz}.

Trabeculae were mounted with aluminium T-clips between a force transducer (Model 400) and a
servo-motor in a permeabilised fibre test system (Model 1400A; Aurora Scientific) at
21\,$^\circ$C. Sarcomere length was measured by video (Model 900B, $40\times$ objective) in
relaxing solution and set to 2.0 $\mu$m; this relaxed reference length defines muscle length
$ML = 1$, and all length commands are expressed in $ML$. Force was normalised to
cross-sectional area, assuming an elliptical cross-section from the measured width and
thickness, and is reported in kPa.

#### 2.1.2 Mechanical protocols

Each activation was preceded by 60 s in preactivating solution and lasted about 15 s at
pCa 4.51. During it the motor applied the following sequence (Fig. **[TODO]**):

1. **Force redevelopment ($k_\mathrm{tr}$).** From the isometric plateau at $ML = 1.0$, a
   release to $0.8$ (held 5 ms), a stretch to $1.05$ (held 2 ms) and a return to $1.0$,
   after which force redeveloped isometrically \citep{Beard2022}.
2. **Staircase.** Four constant-velocity stretches ($0.5$–$1.5\ ML\,\mathrm{s^{-1}}$) from
   $ML = 1.0$ to $1.01$, $1.04$, $1.07$ and $1.10$, each followed by a 140–180 ms isometric
   hold.
3. **Slack/re-stretch train.** Five cycles from $ML = 1.10$, each a rapid ($\approx 2$ ms)
   release to $1.02$, $1.00$, $0.98$, $0.96$ or $0.94$, a 298 ms hold during which force
   redeveloped at the shortened length, and a 20 ms re-stretch back to $1.10$ (at
   $3.9$–$6.8\ ML\,\mathrm{s^{-1}}$); the fifth cycle returned to $1.00$ instead.
4. **Force–velocity point (FV2).** A step to $ML = 1.06$ followed by a 35 ms shortening ramp
   at $2\ ML\,\mathrm{s^{-1}}$, giving relative force during shortening in the same
   trabecula.

The full force–velocity relation was taken from the isovelocity data of \citet{Beard2022}
(8, 4 and 2 mM MgATP; $0.5$–$6\ ML\,\mathrm{s^{-1}}$). The same length sequence was repeated
after myosin inhibition with PNB plus mavacamten to measure the passive (titin-dominated)
response.

#### 2.1.3 Experimental design

Low ATP was studied as a paired comparison within each trabecula, with every low-ATP
activation bracketed by 8 mM activations so that rundown can be interpolated rather than
modelled. Each trabecula was first activated at 8 mM, after which muscle length was readjusted
to restore the relaxed sarcomere length of 2.0 $\mu$m. The sequence was then either
8 → 4 → 8 → 2 → 8 mM or 8 → 2 → 8 → 4 → 8 mM, so that each low concentration occurs in
both positions. Each order was run on one male and one female trabecula, giving four
trabeculae (two male, two female), followed by the PNB-plus-mavacamten run.

> **[NOTE for co-authors]** §2.1.3 describes the target dataset, which has not been
> acquired yet. The fits currently in the draft use three pilot trabeculae (03/27 M,
> 04/03 F, 04/10 M; 8 and 2 mM only, with a single 8 mM bracket on 03/27), and §2.2.1
> (rundown model) is written for them. Rewrite §2.2.1 once the bracketed data exist.

### 2.2 Data reduction

#### 2.2.1 Correction for preparation rundown

A permeabilised trabecula loses force progressively across a session, so a sequential
paired comparison confounds the treatment effect with decay. Two properties of the decay,
established from a same-condition bracket (an 8 mM run repeated at the end of a session),
determine the correction. First, the fibre degrades only while activated: the change
across a relaxed interval is a small *recovery* ($+2.6$ to $+4.0\,\%$ per rest), not
further decay. Second, the decay is linear in accumulated *activated* time and not in
force — two 8 mM runs $694$ s apart gave within-run slopes of $-0.4592$ and
$-0.4618\ \mathrm{kPa\,s^{-1}}$ despite the second being $17.5\,\%$ weaker, where a
force-proportional decay predicts $-0.379$.

Only a fraction $\varphi$ of each run's within-run decline is permanent:

$$
\varphi \;=\; \frac{\Delta F_\mathrm{measured}}{\sum_i |s_i|\,T_i} \;=\; \frac{11.98}{26.52} \;=\; 0.452 ,
$$

with $s_i$ the within-run force slope and $T_i$ the activation duration of run $i$.
Damage accrued between two runs is $\mathrm{dmg} = \varphi \sum_{i\in A}|s_i|\,T_{\mathrm{dec},i}$
over the intervening force-generating activations, and the *later* (degraded) run is
corrected by $f = \mathrm{dmg}/F_\mathrm{later}$, so the low/high ATP ratio becomes
$R_\mathrm{corr} = R_\mathrm{raw}(1+f)$ for $8\to2$ order and $R_\mathrm{raw}/(1+f)$ for
$2\to8$.

Two refinements are applied unconditionally. Each run is referenced to *its own* peak
rather than to a fixed protocol window, because force reaches its peak $3.6$–$5.8$ s after
activation onset and the two ATP conditions do so at different rates; and each run's own
decay phase $T_\mathrm{dec}$ replaces a common activation duration. The own-peak reference
alone, which needs no model, reduces the between-preparation CV of the ATP force ratio from
$20\,\%$ to $17\,\%$; the remaining reduction to $7\,\%$ requires $\varphi$ and the damage
model above.

Within each record, force is additionally detrended by the record's own linear slope,
fitted on the two isometric plateaus that bracket the slack train; this removes an
ATP-dependent distortion (2 mM records decay about twice as fast) of the length dependence.

Across records, the correction is applied to the *second record of each preparation* and to
force **and** slack $k_\mathrm{tr}$ together. The underlying lesion — loss of parallel
force-generating material ($\approx 16\,\%$) plus added series compliance ($\approx 35\,\%$)
— reaches both. Correcting force alone would pair a force belonging to a fresh fibre with a
$k_\mathrm{tr}$ belonging to a degraded one, a state that never existed, and the model would
distort real parameters to fit it. (An earlier recommendation to correct force only, chosen
because it minimised between-preparation CV, is withdrawn: CV conflates correctness with
biological variation, and internal consistency is the stronger requirement.) Full
derivation, validation and limitations in `Analyses/RundownCorrection/methods.md` and
`strategy.md`. **[TODO: confirm the reported $k_\mathrm{tr}$ ratios use this rule.]**

#### 2.2.2 Normalisation and pooling

Force-dimension observables were divided, per preparation, by that preparation's isometric
force measured at $SL = 2.2\ \mu\mathrm{m}$ immediately before the slack train of its 8 mM run
($F_\mathrm{ref}$) — **one divisor per preparation, shared by both ATP conditions**, so
the 2 mM/8 mM ratio is independent of it. Between-preparation force variation is one
common scale: the first principal component of the between-preparation log-variance holds
$88\,\%$ with a flat loading ($|\cos| = 0.955$), which is what justifies a single scalar
divisor; the residual is a length offset ($\pm 0.05\ \mu\mathrm{m}$ in $SL$), not a scale.

Features were extracted from each preparation separately and then averaged; traces were
never averaged before extraction, because feature extraction is non-linear (exponential
fits, peak picking, argmax). Each feature's uncertainty is the between-preparation SD of
the normalised value, moderated by an empirical-Bayes prior estimated across the feature
set \citep{Smyth2004}, since $n = 3$ leaves two degrees of freedom per estimate;
moderation compresses the inverse-variance weight range $446\times \to 15\times$.
Absolute tension, which the normalisation removes, is reported separately:
$F_\mathrm{ref} = 69.3$ kPa, CV $15.5\,\%$ ($n=3$, rundown-corrected), and re-enters the
objective as a single explicit level term (§2.5).

#### 2.2.3 Feature extraction

Records are segmented by the commanded velocity table and reduced to scalar features by
`Model/extractSlackAttributes.m`, `Model/extractKtrAttributes.m` and
`Model/extractForceVelocityAttributes.m`. Identical code paths are applied to data and to
simulation, except that experimental traces are pre-smoothed (a $5$ ms moving median to
locate the undershoot; robust LOESS, $25$ ms, for the recovery tail) whereas simulation
output is already smooth. The rapid re-stretch peak is read from the raw trace.

The principal features are:

| feature | definition |
|---|---|
| `steady` | re-developed isometric plateau after each slack cycle |
| `A` | amplitude of the post-slack force recovery |
| `ktr` | first-order rate of the post-slack force recovery |
| `t0_crossing` | delay from the end of the slack release to force onset |
| `peak1_y`, `peak1_dSL` | height and length-offset of the first re-stretch force peak |
| `vall_y`, `vall_t` | depth and timing of the post-peak valley |
| `peak2` | second (viscoelastic) force peak |
| `rsK80`, `rsT20` | post-re-stretch recovery rate and delay, read off the 20–80 % band of the rise with the amplitude taken from `steady` |
| `ktrP_k63` | redevelopment rate from the dedicated $k_\mathrm{tr}$ protocol |
| `FV_normpowerAvg` | normalised power across the force–velocity curve |
| `PS_*` | passive-run features (re-stretch peak, ramp-up and hold-decay RMSE, plateaux at $SL$ 2.0 and 2.2) |

Two extraction choices are load-bearing. `rsK80`/`rsT20` replace an earlier pair whose
amplitude was taken from the tail of their own fit window: that estimator returned a known
$\times 0.6$ time dilation as $24\,\%$ wrong and moved $77\,\%$ when the window was
truncated to 60 ms, whereas the 20–80 % band with the amplitude anchored to `steady`
recovers it to $6\,\%$ and moves $0\,\%$, dropping the between-preparation SD from
$17.2\,\%$ to $6.1\,\%$. And `rsK80` is averaged over slack cycles 1–4 only (the other
slack features use all five): cycle 5
sits off the amplitude ladder and carries $29$–$245\,\%$ per-cycle measurement noise
against $6$–$15\,\%$ for cycles 1–4, in every run at both ATP concentrations.

The measurement noise of every feature was quantified directly by injecting each run's own
measured residual noise into its own trace $200$ times and re-extracting
(Analyses/RecoveryOscillation). Applied to the force recoveries, the same procedure shows
that the residual of a single-exponential fit is $90\,\%$ white transducer noise
($\sigma = 1.19$ kPa with the cross-bridges chemically inhibited against $1.24$ kPa
active), $6\,\%$ a slow ramp-locked wander and $4\,\%$ unexplained, with no discrete
resonance. A model sitting exactly on the truth could therefore reach at most
$R^2 = 0.834$ on such a window, which is why goodness-of-fit on a recovery trace is used
only as a boundary constraint and never as a fit target.

### 2.3 Model

#### 2.3.1 State space

Parameter values and active switches quoted below are those of the reported 8 mM set,
`params/robust1_seed.m` (the `refit9` lineage), from which the 2 mM results of §1.6 are
computed. **[TODO: update to the final joint-fit set.]**

The model tracks the probability distribution of a myosin head over six pools. Two are
detached and cycling-competent (disordered-relaxed):

- $D_T$ — $\mathrm{M\cdot ATP}$, ATP bound, not yet hydrolysed;
- $D_D$ — $\mathrm{M\cdot ADP\cdot P_i}$, primed and attachment-competent.

Two are attached and resolved over cross-bridge strain $s$ (in $\mu\mathrm{m}$), as densities:

- $p_1(s,t)$ — pre-power-stroke (A1), weakly attached;
- $p_2(s,t)$ — post-power-stroke (A2), strongly attached and force-bearing.

Two are super-relaxed (folded back onto the thick filament, interacting-heads motif), with
the same nucleotide distinction:

- $S_T$ — SRX with ATP; $S_D$ — SRX with $\mathrm{ADP\cdot P_i}$.

Probability is conserved,

$$
\int p_1\,\mathrm{d}s + \int p_2\,\mathrm{d}s + D_T + D_D + S_T + S_D = 1 ,
$$

and $D_T$ is carried as the complement of the others rather than integrated:

$$
D_T \;=\; \max\!\Big(0,\; 1 - \textstyle\int p_1 - \int p_2 - D_D - S_T - S_D\Big).
$$

Conservation is therefore exact up to two numerical safeguards — negative density values
are clipped to zero inside the right-hand side, and the $\max(0,\cdot)$ above — and up to
the grid roll of §2.3.2.

Further state variables describe the mechanics and the slow availability processes: the
commanded half-sarcomere-pair length $SL$ (the motor length expressed per sarcomere), the
series-elastic length $L_\mathrm{SE}$, from which the actual contractile length
$L_c = SL - L_\mathrm{SE}$ follows, a titin viscous stress $X_v$, a registration-availability
state $A_\mathrm{reg}$ and a compliant-realignment deficit $A_\mathrm{def}$
(§2.3.8–2.3.10). A third attached state $p_3$ is implemented and may be enabled
(`NumberOfStates = 3`); the results reported here use the two-attached-state form.

#### 2.3.2 Strain, filament sliding and the transport equation

**Physics.** A head attached to actin carries a strain $s$ (in $\mu$m): the displacement of
its lever from the position at which it attached. While the head stays attached, relative
sliding of the thick and thin filaments changes $s$ one-for-one. In a half-sarcomere the
filaments slide past each other at half the rate at which the contractile length changes,
so every attached head's strain evolves as

$$
\frac{\mathrm{d}s}{\mathrm{d}t} \;=\; \frac{v_c(t)}{2},\qquad v_c \;=\; \dot L_c \;=\; v - \dot L_\mathrm{SE},
$$

where $v = \dot{SL}$ is the imposed (motor) velocity and $v_c$ the velocity of the
contractile element, which differs from it by the extension rate of the series element.
Lengthening ($v_c > 0$) strains attached heads positively; shortening carries them to
negative strain. The attached densities $p_k(s,t)$ are therefore transported along $s$, and
obey the first-order hyperbolic (Huxley-type) balance

$$
\frac{\partial p_k}{\partial t} \;+\; \frac{v_c(t)}{2}\,\frac{\partial p_k}{\partial s} \;=\; \mathcal{R}_k(s,t),\qquad k = 1,2,
$$

where $\mathcal{R}_k$ collects the attachment, transition and detachment fluxes of §2.3.3.
The advection term carries the whole effect of velocity on the distribution: it is the reason
a shortening muscle bears less force (heads are swept to negative strain before they
detach), and the reason a stretch produces a fast force transient (the whole attached
population is strained at once). The series element enters through $v_c$ only: bridge
strain follows the contractile element, not the motor.

**Numerical realisation (method of characteristics).** The model does not discretise
$\partial p/\partial s$. It integrates the equation along its characteristics: each of
$n_s = 60$ bins is a *material* bin, holding the heads attached at one strain, and the bin's
strain moves with the contractile length,

$$
s_i(t) \;=\; \bar s_i \;+\; \tfrac{1}{2}\big(L_c(t) - L_\mathrm{piv}\big),\qquad
\bar s_i = (i - 30.5)\,\Delta s,\quad \Delta s = 2\ \mathrm{nm},
$$

so that along each bin $\mathrm{d}p_{k,i}/\mathrm{d}t = \mathcal{R}_k(s_i(t),t)$ is the total
derivative of the transport equation, exact by construction and free of numerical
diffusion. $L_\mathrm{piv}$ is the contractile length at which the grid is centred
($2.2\ \mu$m initially). What moves relative to the bins is not the population but the
point $s = 0$ at which new heads attach, and every strain-dependent rate is evaluated at
the bin's current strain $s_i(t)$. The densities $p_k$ carry units of $\mu\mathrm{m}^{-1}$;
moments are $\int p_k\,\mathrm{d}s \approx \Delta s\sum_i p_{k,i}$.

**Attachment kernel.** Attachment deposits probability around $s = 0$ with a triangular
kernel of half-width $A_x = 6$ nm, representing the spread of accessible actin target
positions; it also keeps the injection smooth as the zero-strain point moves across bins:

$$
K(s_j) \;=\; \frac{w_j}{\Delta s\,\max\!\big(1,\sum_j w_j\big)},\qquad w_j = \max\!\big(0,\,1 - |s_j|/A_x\big).
$$

**Grid roll.** The 60-bin window spans $\pm 59$ nm about $L_\mathrm{piv}$. When sustained
sliding brings the zero-strain point within $A_x$ of an edge, an ODE event stops the
integration, the arrays are shifted by $20$ bins ($40$ nm) and $L_\mathrm{piv}$ is moved by
the same strain, so physical strain is preserved. The $20$ bins at the trailing edge — heads
strained by $\approx 74$–$112$ nm in the direction of travel — are discarded and the vacated
leading bins are zeroed. Because $D_T$ is the complement, the discarded mass reappears as
detached, ATP-bound heads: physically, an instantaneous forced detachment of extremely
strained bridges, whose detachment rates there exceed $5\times10^3\ \mathrm{s^{-1}}$. Logged
over the full protocol set, each roll discards at most $10^{-5}$ of the attached mass, and
the objective is converged in the bin width to within $1.8\,\%$ at $\Delta s = 2$ nm
(Analyses/StrainWindowNoise).

#### 2.3.3 Governing equations

Writing the left-hand side as the material derivative along the sliding characteristics,
$\mathrm{D}/\mathrm{D}t = \partial_t + \tfrac12 v_c\,\partial_s$, the attached states obey

$$
\begin{aligned}
\frac{\mathrm{D}p_1}{\mathrm{D}t} &= R_{D1}\,K(s) \;-\; k_d\phi_{1D}(s)\,p_1 \;-\; k_1\phi_{12}(s)\,p_1 \;+\; k_{-1}\phi_{21}(s)\,p_2,\\[2pt]
\frac{\mathrm{D}p_2}{\mathrm{D}t} &= k_1\phi_{12}(s)\,p_1 \;-\; k_{-1}\phi_{21}(s)\,p_2 \;-\; k_2\phi_2(s+\Delta_2-\Delta)\,g_\mathrm{ATP}\,p_2 \;+\; \mathbb{1}_{v>0}\Big[\big(k_\mathrm{re}D_T + H\big)\,\delta_0(s) \;-\; k_\mathrm{hop}\,(s-s_\mathrm{hop})_+\,p_2\Big],
\end{aligned}
$$

with $H = \int k_\mathrm{hop}(s-s_\mathrm{hop})_+\,p_2\,\mathrm{d}s$, and the lumped pools obey

$$
\begin{aligned}
\frac{\mathrm{d}D_D}{\mathrm{d}t} &= k_{ah}D_T - k_{amh}D_D - R_{D1} + \int k_d\phi_{1D}\,p_1\,\mathrm{d}s + R_{S_D\to D_D},\\[2pt]
\frac{\mathrm{d}S_T}{\mathrm{d}t} &= R_{D_T\to S_T} - k_{sr2srd}S_T + k_{srd2sr}S_D + \int f_\mathrm{srx}(s)\,k_2\phi_2\,g_\mathrm{ATP}\,p_2\,\mathrm{d}s,\\[2pt]
\frac{\mathrm{d}S_D}{\mathrm{d}t} &= k_{sr2srd}S_T - k_{srd2sr}S_D - R_{S_D\to D_D}.
\end{aligned}
$$

Hydrolysis and its reverse ($k_{ah}$, $k_{amh}$) connect the detached pools. Weakly bound
heads (A1) detach to the primed pool $D_D$; strongly bound heads (A2) detach after ATP
binding, to $D_T$ (implicitly, through the complement) — except for a strain-dependent
fraction
$f_\mathrm{srx}(s) = 0.219\,\mathrm{clip}\big((s - s_\mathrm{rip})/\delta_\mathrm{rip},0,1\big)$
($s_\mathrm{rip} = 13.9$ nm, $\delta_\mathrm{rip} = 5$ nm), which re-folds directly into
the super-relaxed state: heads torn off at high positive strain return to the thick-filament
backbone rather than to the cycle. The SRX fluxes are given in §2.3.6. Attachment is

$$
R_{D1} \;=\; k_a\,A_\mathrm{reg}\,(1 - A_\mathrm{def})\;N_\mathrm{ov}\; f_\mathrm{lat}\; f_\mathrm{occ}\; D_D .
$$

The bracketed term on $p_2$ is an optional mechanism that is **disabled in the reported set**
($k_\mathrm{re} = k_\mathrm{hop} = 0$); it was active in the preceding lineage and is
described because it bears on the re-stretch residual (§1.6). It acts only while the muscle
is being lengthened ($v>0$, imposed velocity). Heads strained beyond $s_\mathrm{hop} = 4.75$ nm slip to a neighbouring actin site
at a rate proportional to the excess strain ($k_\mathrm{hop} = 3.6\times10^4\
\mu\mathrm{m^{-1}s^{-1}}$), re-entering $p_2$ at zero strain (mass-conserving), and detached
ATP-bound heads re-bind directly into the strongly bound state at $k_\mathrm{re} = 17.8\
\mathrm{s^{-1}}$. Both re-entries are point injections at $s = 0$ ($\delta_0$, a linear
two-bin split), not the attachment kernel; the $k_\mathrm{re}$ flux is not subject to the
overlap, lattice or availability factors.

#### 2.3.4 Strain dependence of the attached-state transitions

The four attached-state rates carry dimensionless, strain-dependent modifiers
$\phi_\bullet(s)$ built as shape-preserving piecewise cubic (PCHIP) interpolants through a
small set of fitted knots, which lets the optimiser shape each rate profile without
committing to an exponential or Gaussian form:

$$
\begin{aligned}
R_{1D}(s) &= k_d\,\phi_{1D}(s)\,p_1(s) && \text{A1 detachment}\\
R_{12}(s) &= k_1\,\phi_{12}(s)\,p_1(s) && \text{power stroke}\\
R_{21}(s) &= k_{-1}\,\phi_{21}(s)\,p_2(s) && \text{reverse stroke}\\
R_{2}(s)  &= k_2\,\phi_2\!\big(s + \Delta_2 - \Delta\big)\,g_\mathrm{ATP}\,p_2(s) && \text{A2 detachment}
\end{aligned}
$$

Each modifier is $\phi(s) = \mathrm{PCHIP}\{x_j + \ln(o)/10,\; y_j\}(s)$: the knot abscissae
$x_j$ carry a fitted common translation $\ln(o)/10$, which is negligible for three of the
four modifiers ($-0.4$, $0$ and $+1.4$ nm) but moves the reverse-stroke profile $\phi_{21}$ by
$-20.5$ nm. The $\phi_2$ argument is shifted by $\Delta_2 - \Delta = 20.3$ nm. Beyond the
outer knots each modifier is held at its end value (`ClampStrainDepExtrapolation`): free
cubic extrapolation can turn a rate negative at large strain. The knot sets are
fitted parameters (`PieceWiseStrainDep*X` / `*Params`); at the reported solution $\phi_{12}$ rises steeply on
the negative-strain side (shortening promotes the stroke) while $\phi_2$ is V-shaped about
the post-stroke equilibrium (both compression and stretch accelerate detachment), which is
the qualitative behaviour inferred by \citet{Tewari2016}.

#### 2.3.5 Nucleotide dependence

Because the nucleotide state is explicit, $[\mathrm{MgATP}]$, $[\mathrm{MgADP}]$ and
$[\mathrm{P_i}]$ act at the steps they chemically occupy. Each gate is written as a
Michaelis factor **normalised to unity at a reference concentration**
$[\mathrm{ATP}]_\mathrm{ref} = 8$ mM, so enabling a gate leaves the 8 mM calibration
bit-for-bit unchanged for any $K_m$ and the low-ATP prediction is a genuine
out-of-sample extrapolation:

$$
g_\mathrm{ATP}^{(x)} \;=\; \frac{[\mathrm{ATP}] / \big([\mathrm{ATP}] + K_x\big)}{[\mathrm{ATP}]_\mathrm{ref} / \big([\mathrm{ATP}]_\mathrm{ref} + K_x\big)} .
$$

Gates are *available* at each ATP-dependent step; in the 8 mM reference set they are
identically 1 by construction, and they are switched on only for the low-ATP arm. The
**cycle gate** multiplies the A2 detachment rate: a post-stroke head must bind ATP before it
can release, so low ATP retains strongly-bound bridges — the force lever (H2). The
**availability gates** act on the super-relaxed pool: on the return $S_D \to D_D$
(re-priming a folded head) and on the entry $D_T \to S_T$ — the kinetic levers, in which
isometric force is nearly flat (H3). For the exploratory low-ATP fits reported in §1.6 the
gates are written as power laws, $g = ([\mathrm{ATP}]/[\mathrm{ATP}]_\mathrm{ref})^{p}$, whose
single exponent is the identified quantity; the Michaelis form above is used when a
$K_m$ is to be reported. **[TODO: fix the final form with the joint fit.]**

$[\mathrm{MgADP}]$ and $[\mathrm{P_i}]$ do **not** enter the present right-hand side. The
\citet{Beard2022} forms (MgADP competition at the nucleotide site, $K_D = 0.194$ mM,
\citealp{Yamashita1994}; $\mathrm{P_i}$ reversal of the force-generating step) remain in the
code but are disabled. This is consistent with an ADP-buffered bath; it is *not*
consistent with $\mathrm{P_i}$ accumulating in the preparation core (§3.1), and re-enabling
the $\mathrm{P_i}$ term is part of the analysis planned there.

#### 2.3.6 Thick-filament regulation

Entry to and exit from the super-relaxed pools are gated exponentially by the stress borne
by the contractile element, $F_\mathrm{tot} = F_\mathrm{XB} + F_\mathrm{pas}$ (§2.3.9) —
the model's representation of thick-filament mechanosensing, as in the parent model
\citep{Linari2015, Campbell2018, Beard2022}:

$$
\begin{aligned}
R_{D_T\to S_T} &= k_{sr0}\,e^{-F_\mathrm{tot}/\sigma_2}\,D_T, &\qquad R_{S_T\to D_T} &= k_{msr}\,e^{+F_\mathrm{tot}/\sigma_1}\,S_T,\\
R_{D_D\to S_D} &= k_{srd}\,e^{-F_\mathrm{tot}/\sigma_{srd2}}\,D_D, &\qquad R_{S_D\to D_D} &= k_{msrd}\,e^{+F_\mathrm{tot}/\sigma_{srd1}}\,g^{(\mathrm{srx})}_\mathrm{ATP}\,S_D,\\
R_{S_T\to S_D} &= k_{sr2srd}\,S_T, &\qquad R_{S_D\to S_T} &= k_{srd2sr}\,S_D .
\end{aligned}
$$

Force suppresses entry and accelerates exit, so rising stress mobilises heads — the
Frank–Starling direction. Note that the two gates point the same way, so the in/out ratio
swings steeply with load; at the reported solution it moves $\approx 60\times$ between
zero force and $0.87\,F_\mathrm{iso}$, which makes force redevelopment at high load
intrinsically faster than at low load (Analyses/RestretchSRXRecoil).

In the reported parameter set two of the six legs refine to zero ($k_{msr} = 0$,
$k_{srd} = 0$), so the SRX pools operate as a directed cycle
$D_T \to S_T \rightleftarrows S_D \to D_D$: a head folds back from the ATP-bound detached
state, hydrolyses while folded, and is released as $\mathrm{M\cdot ADP\cdot P_i}$. This is
*non-identification*, not identification — the reverse legs remain in the equations and
were re-tested with their gates artificially revived, which changed the objective by less
than the optimiser noise (Analyses/MechanismAblation), so the data cannot distinguish the
directed cycle from the reversible one. The scheme is best read as four effective pools
with two unidentified reverse rates. A third entry route, from strongly bound heads torn
off at high strain, is given in §2.3.3.

During activation the model's SRX is a fast, load-gated switch (isometric occupancy
$0.03$, mean dwell $\approx 10$ ms); it is kinetically the structural OFF state of X-ray
work, and coincides with the slow biochemical (mant-ATP) SRX only at rest
(Analyses/SRXDuringContraction).

#### 2.3.7 Filament geometry

Thin/thick filament overlap scales attachment:

$$
N_\mathrm{ov} = \frac{2 L_\mathrm{ov}}{L_\mathrm{thick} - L_\mathrm{hbare}},\qquad
L_\mathrm{ov} = \min\!\Big(\tfrac{L_\mathrm{thick}}{2},\tfrac{SL}{2}\Big) - \max\!\Big(\tfrac{L_c}{2} - (L_c - L_\mathrm{thin}),\, \tfrac{L_\mathrm{hbare}}{2}\Big).
$$

With `UseOverlapFactor`, $N_\mathrm{ov}$ is re-expressed as the fraction of *unattached*
heads lying within the overlap zone,
$N_\mathrm{ov} \leftarrow \big(N_\mathrm{ov} - \int(p_1+p_2)\big)/(S_T + S_D + D_T + D_D)$,
so that $k_a D_D N_\mathrm{ov}$ is a flux per available head rather than per head.

Myofilament lattice spacing attenuates attachment at short lengths:

$$
d_{10} = d_{10}^\mathrm{ref}\sqrt{SL_\mathrm{ref}/SL}, \qquad
f_\mathrm{lat} = \exp\!\left[-\frac{\max\!\big(0,\, d_{10}-R_\mathrm{thick}-R_\mathrm{thin}-d_\mathrm{opt}\big)^2}{2\sigma_\mathrm{lat}^2}\right].
$$

A degeneracy analysis shows that these six parameters carry only two identifiable
combinations ($d_{10}^\mathrm{ref}\sqrt{SL_\mathrm{ref}}$ as a product and
$R_\mathrm{thick}+R_\mathrm{thin}+d_\mathrm{opt}$ as a sum), reducing operationally to
*attachment is attenuated below* $SL \approx 2.09\ \mu\mathrm{m}$; they are constrained but not
individually identified.

Finally, a mean-field occupancy term attenuates the whole attachment flux with the bound
fraction,

$$
f_\mathrm{occ} = \max\!\left(0,\; 1 - \frac{\int (p_1+p_2)\,\mathrm{d}s}{P_b^\mathrm{max}}\right).
$$

It is deliberately *global*: the strain axis is not a binding-site axis, so a per-bin
occupancy term would be both a category error and grid-spacing dependent
(Analyses/BindingSiteOccupancy). The fitted $P_b^\mathrm{max} = 1.21$ lies above the
physically attainable bound fraction, so the term does not saturate; at the operating point
it acts as a soft rescale of attachment ($f_\mathrm{occ} \approx 0.76$ at 8 mM and $0.69$
at 2 mM, from bound fractions of $0.29$ and $0.37$). Note that the overlap and lattice
factors use the commanded $SL$, while the occupancy correction inside $N_\mathrm{ov}$
uses the contractile length $L_c$.

#### 2.3.8 Attachment availability

Two slow scalar states gate attachment; both are constructed to be exactly inert at steady
state, so neither can rescale $k_a$ or move the operating point — they act only on
transients.

**Registration availability $A_\mathrm{reg}$** represents the fraction of myosin heads in
register with an accessible actin target zone. It relaxes first-order toward a
velocity-set target,

$$
\dot A_\mathrm{reg} = \frac{A_\infty - A_\mathrm{reg}}{\tau_\mathrm{reg}}, \qquad
A_\infty = A_0 + (1-A_0)\,\frac{v_\mathrm{reg}}{v_\mathrm{reg} + v_\mathrm{ref}},\qquad
v_\mathrm{reg} = \max(0, -v),
$$

driven by the *imposed* velocity so that series-elastic ringing cannot corrupt it. In the
reported set $\tau_\mathrm{reg} = 0.2$ ms, so $A_\mathrm{reg}$ follows its target
essentially instantaneously and acts as a static velocity-dependent factor on attachment;
the earlier lineage identified $\tau_\mathrm{reg} = 27$ ms, and the data cannot yet
separate the two. Using
shortening speed only ($v_\mathrm{reg} = \max(0,-v)$) keeps the force–velocity "shoulder"
while leaving the post-re-stretch undershoot intact. This mechanism, identified from the
steady force–velocity relation alone, turns out to carry the whole of the model's
predicted asymmetric response to a symmetric sawtooth length input
(Analyses/SawtoothSRS) — an independent, untested prediction.

**Compliant realignment $A_\mathrm{def}$** represents the register disturbance left behind
when bound bridges are stripped *while the filament remains loaded*. Because bound bridges
set the register of neighbouring actin sites through filament compliance
\citep{Daniel1998, Tanner2007}, stripping them under load leaves sites transiently
unavailable:

$$
\dot A_\mathrm{def} = k_\mathrm{cr}\,\mathcal{L}\big(-\dot P_b\big)\,\frac{(F_\mathrm{tot})_+}{(F_\mathrm{tot})_+ + F_\mathrm{cr}} \;-\; \frac{A_\mathrm{def}}{\tau_\mathrm{cr}},\qquad
P_b = \int (p_1+p_2)\,\mathrm{d}s ,\quad A_\mathrm{def}\le 0.9,
$$

where $\dot P_b$ is the net *reaction* change of bound mass, $\int(\mathcal R_1 + \mathcal R_2)\,\mathrm{d}s$
— advection along strain moves bound heads but does not change their number —

and $\mathcal{L}(u) = u/(1+e^{-u/\epsilon})$ is a smooth (SiLU) rectifier — a hard
$\max(0,\cdot)$ on a derivative switches continuously during post-stretch ringing and
stalls the stiff solver — and $\epsilon = 0.5\ \mathrm{s^{-1}} \ll$ the physiological
stripping rate. The driver is the *net* change in bound mass, so $\dot A_\mathrm{def}=0$
at any steady state identically. The load weighting is what makes the mechanism
selective: after a re-stretch, $p_2$ is stripped while force holds at
$\approx 0.84\,F_\mathrm{iso}$ (large loss $\times$ large weight); after a slack release or a
$k_\mathrm{tr}$ manoeuvre, $p_2$ collapses but force goes to zero with it (large loss
$\times$ zero weight). It therefore slows the one recovery window that is too fast and
leaves the two that are correct.

Both gates multiply the attachment rate:
$k_a^\mathrm{eff} = k_a\,A_\mathrm{reg}\,(1 - A_\mathrm{def})$. The lengthening-only
re-binding flux $k_\mathrm{re}D_T$ of §2.3.3 bypasses both.

#### 2.3.9 Force and the mechanical boundary condition

Cross-bridge stress is the first moment of the attached distributions about each state's
own equilibrium, with an exponent $e$ on the strongly-bound state allowing a non-Hookean
bridge, and separate stiffnesses for positive and negative strain:

$$
F_\mathrm{XB} \;=\; k_{s1}^{\pm}\!\int (s+\Delta_1)\,p_1\,\mathrm{d}s \;+\; k_{s2}^{\pm}\!\int \mathrm{sgn}(s+\Delta)\,|s+\Delta|^{\,e}\,p_2\,\mathrm{d}s ,
$$

where each integral is split at its state's zero-force strain ($s = -\Delta_1$, $s = -\Delta$),
with stiffness $k^+$ on the stretched side and a softer $k^-$ on the compressed side, and
$e = 1.24$. $\Delta = 11.3$ nm is the power-stroke size and $\Delta_1 = 4.6$ nm an offset that
lets the pre-stroke state bear force, decoupling force duty from the $p_2$ detachment that
otherwise welds $k_\mathrm{tr}$ to the force–velocity relation. Because the bridge strain
$s$ is advected by $\tfrac12 v_c$, $F_\mathrm{XB}$ responds to stretch instantaneously,
through the strain of the population already attached, and to steady sliding through the
shift of the distribution toward negative strain.

The contractile element acts in series with an elastic element (end compliance plus
damaged-region compliance; the half-sarcomere does not follow the motor,
cf. \citealp{deTombe2016}). The motor imposes the total length; the stress measured at the
transducer is the stress in the series element,

$$
F \;=\; \mathrm{sgn}(L_\mathrm{SE})\,\big|\max(0,\, k_\mathrm{SE} L_\mathrm{SE})\big|^{\,e_\mathrm{SE}},
$$

which is zero whenever the series element is slack ($L_\mathrm{SE}\le 0$, during a slack
release). The contractile length moves by a viscously regularised force balance between
that stress and the stress the contractile element bears,

$$
\dot{SL} = v(t), \qquad
v_c = \dot L_c = \frac{F - F_\mathrm{tot}}{\mu^{\pm}},\qquad
\dot L_\mathrm{SE} = v(t) - v_c,\qquad
F_\mathrm{tot} = F_\mathrm{XB} + F_\mathrm{pas},
$$

where $\mu^{\pm}$ is a drag coefficient on contractile sliding: $\mu = 0.020$ when the
contractile element is pulled out ($F > F_\mathrm{tot}$) and $\mu_\mathrm{neg} = 0.016$ when
it shortens against the series element. In the limit $\mu\to0$ this is the quasi-static
balance $F = F_\mathrm{tot}$; a finite $\mu$ both regularises the differential-algebraic
structure for the stiff solver and represents the internal viscosity that limits unloaded
shortening \citep{deTombe1992}. The resulting $v_c$ is the velocity that advects the
strain densities (§2.3.2), which closes the loop between mechanics and the cross-bridge
distribution: during a slack release the series element unloads, $F = 0$, and the
contractile element shortens at the rate at which $F_\mathrm{tot}$ can be dragged through
the viscosity, sweeping the attached heads to negative strain.

The motor velocity table is in muscle lengths per second with $ML = 1$ corresponding to
$2\ \mu$m per sarcomere; $v$ above is the full-sarcomere rate in $\mu$m/s.

#### 2.3.10 Passive element

Passive stress is a non-linear elastic term in parallel with a Maxwell (viscous)
branch representing stress-mediated unfolding of serial titin domains \citep{Jezek2025}:

$$
F_\mathrm{pas} = k_\mathrm{pas}\,\max(L_c - L_{sc0},\,0)^{\gamma} \;+\; X_v ,
$$

$$
\dot X_v \;=\; k^{M}_\mathrm{SE}\,\max(0,\dot L_c) \;-\; \frac{X_v}{\eta_M} \;-\; \frac{X_v\,\max(0,-\dot L_c)}{x_{M,\mathrm{slack}}} .
$$

The first term charges the branch during lengthening, the second relaxes it with time
constant $\eta_M$, and the third — the tension-only term — discharges it over a slack
*length* constant $x_{M,\mathrm{slack}}$, encoding that titin pulls but cannot push: any
shortening pays out slack and releases the stored unfolding stress. The two rectifiers are
the smooth SiLU form of §2.3.8 ($\epsilon = 0.01$): with hard $\max(0,\cdot)$ the corner at
$\dot L_c = 0$ sits exactly at isometric rest and stalls the solver at tight tolerance. The
branch is driven by the contractile velocity $v_c$, not the motor velocity. Ca$^{2+}$-dependence of the passive limb — PEVK–actin attachment, proximal
stiffening and stabilisation of folded serial domains — is the mechanism identified in
\citet{Jezek2025}; the present model absorbs its net effect into the identified
$(k_\mathrm{pas}, \gamma, L_{sc0}, \eta_M, k^M_\mathrm{SE}, x_{M,\mathrm{slack}})$ fitted
to the PNB+mavacamten runs, because these protocols hold Ca$^{2+}$ constant.

Model passive parameters are identified on a passive protocol that was recorded *after*
four activations and therefore reads the most-rundown state; three independent lines of
evidence (rundown accounting, the 2 mM cost optimum, and the tolerance of the direct
passive fit) put the passive contribution appropriate to the active runs at
$\approx \times 1.2$ of that value. **[TODO: decide whether to apply that factor in the
reported fits or to report it as a sensitivity.]**

#### 2.3.11 Optional mechanisms

The implementation carries a number of further mechanisms behind flags, each default-off
and verified inert when off. Those active in the reported configuration are: the
lengthening-only slip and re-binding fluxes on $p_2$ (`UseA2Reattaching`, §2.3.3), which
fill the post-re-stretch valley and cap the first re-stretch peak; the strain-dependent
re-folding of strongly bound heads torn off at high strain into SRX
(`SRXFromR2HighStrain`, §2.3.3); and the tension-only Maxwell discharge above. Those tested
and rejected include a detachment ceiling with a nucleotide-free refractory pool, catch-bond
detachment suppression during stretch, *uniform* re-routing of the A2 detachment flux into
SRX (`SRXFromR2`), and fixed-dwell (occupancy-saturating) pool exits; each is documented with
its falsifying observable in `Analyses/`. A single-mechanism ablation over 34 arms (30
completed), each given a scalar calibration plus a budget-matched re-optimisation with its
own compensators made compulsory, found no addition worth its parameters at the campaign's
acceptance threshold, and four mechanisms removable at no cost. One arm that failed to
complete was tracking below the incumbent; it must be re-run before the ablation result is
reported. **[TODO]**

#### 2.3.12 Thin-filament Ca$^{2+}$ activation

All protocols above are recorded at saturating Ca$^{2+}$ (pCa $\approx 4.5$), where the thin
filament is fully permissive, and the reported fits use no Ca$^{2+}$ layer. To carry the
model into twitches — and into the beating-heart setting the energetics question is
ultimately about — a thin-filament activation layer after \citet{Campbell2018} is
implemented as an extension (branch `CaActivation`; `Model/calciumActivation.m`). Binding
sites switch between OFF and ON with cooperative rates; following the MyoSim form, only
sites within the overlap zone, $N_\mathrm{ov}$, are switchable, and sites holding a bound
head cannot switch off (cross-bridge feedback on the thin filament):

$$
\begin{aligned}
\frac{\mathrm{d}N_\mathrm{on}}{\mathrm{d}t} &= k_\mathrm{on}[\mathrm{Ca^{2+}}]\,\big(N_\mathrm{ov} - N_\mathrm{on}\big)\Big(1 + K_\mathrm{coop}\frac{N_\mathrm{on}}{N_\mathrm{ov}}\Big)
\;-\; k_\mathrm{off}\,\big(N_\mathrm{on} - P_b\big)\Big(1 + K_\mathrm{coop}\frac{N_\mathrm{ov} - N_\mathrm{on}}{N_\mathrm{ov}}\Big),
\end{aligned}
$$

with $P_b = \int(p_1+p_2)\,\mathrm{d}s$. Attachment then sees only the ON sites that are not
already occupied, through a multiplicative gate on $k_a$,

$$
g_\mathrm{Ca} \;=\; \frac{N_\mathrm{on} - P_b}{N_\mathrm{ov} - P_b},
$$

which is exactly 1 at full activation, so the saturating-Ca$^{2+}$ calibration is unchanged
by construction (verified bit-identical on the reference parameter set). The two-state
form of \citet{Tewari2016} and \citet{Campbell2018} without the overlap restriction is the
special case $N_\mathrm{ov} = 1$. The Ca$^{2+}$ input is either a constant pCa, an analytic
transient (Rice, Marzban and Tewari forms), or a measured signal.

The thin-filament constants are not identified from the present data, which contain no
submaximal-Ca$^{2+}$ records free of myosin inhibitors. With rat values from
\citet{Marzban2020} ($k_\mathrm{on} = 101\ \mu\mathrm{M^{-1}s^{-1}}$,
$k_\mathrm{off} = 724\ \mathrm{s^{-1}}$, $K_\mathrm{coop} = 9.7$) the model gives
pCa$_{50} = 5.33$ and $n_H = 3.3$ and a 1 Hz twitch of only $\sim5\,\%$ of maximal force,
so these constants do not transfer between models and must be identified for this one.
Fitted to twitches of human iPSC-derived engineered heart tissue the layer reproduces the
force time course to $0.8\,\%$ RMS, but absolute [Ca$^{2+}$] and $k_\mathrm{on}$ are exactly
degenerate in a twitch, so the identification is not unique
(Analyses/CalciumTwitch). The prediction carried into future experiments is that the
larger bound fraction at low ATP shifts pCa$_{50}$ leftward through the cross-bridge
feedback term — i.e. low ATP should raise force more at submaximal Ca$^{2+}$ than at
saturation, the direction \citet{Rhodehamel2025} report. Force–pCa relations at 8 and
2 mM are the experiment that identifies the layer.

### 2.4 Numerical implementation

The system is integrated with MATLAB's `ode15s` (variable-order BDF, appropriate to the
stiffness introduced by the fast strain-dependent rates) with relative tolerance $10^{-3}$
and absolute tolerance $10^{-6}$, and the solution is read on a fixed $0.5$ ms output grid
(features are never read off the solver's own adaptive step points). Each protocol is supplied as a velocity table:
a list of $(t_k, v_k)$ breakpoints defining constant-velocity segments, integrated
segment-by-segment with the state carried across, so commanded length changes are exact
and no interpolation of the drive signal is required inside the right-hand side. An event
function terminates a step when the zero-strain point approaches a grid boundary (grid
roll, §2.3.2) and a second event enforces a wall-clock ceiling per simulation, so a
pathological parameter draw during optimisation fails fast rather than stalling the
campaign.

The right-hand side is available both as MATLAB (`Model/dPUdT_CombinedTransitions.m`) and
as a compiled MEX core (`Model/dPUdT_core_mex.c`); the piecewise-cubic rate modifiers are
pre-unpacked to breakpoint/coefficient arrays at parameter-resolution time and evaluated
inline, avoiding per-call `ppval` overhead. Optimisation campaigns run over a MATLAB
parallel pool; because pool workers cache compiled code, the pool is refreshed after any
edit to a model file.

All parameters flow through a single constructor, `Model/getParams.m`, which fills
defaults, applies multiplicative modifiers, resolves parameters defined as expressions of
other parameters, rebuilds the strain grid and re-initialises the state vector. Parameter
sets are serialised to plain `.m` snapshots (six significant figures; the round-trip
changes the objective by a negligible amount), so every reported result is reproducible
from a file.

### 2.5 Objective function

Calibration is on features, not traces. For each feature $k$ with $n_k$ elements,
data values $\bar f_{k,j}$ and model values $f^\mathrm{mod}_{k,j}$, the unweighted cost is
a sum of squared *relative* residuals with an optional deadband:

$$
c_k \;=\; \max\!\left(0,\;\; \sum_{j=1}^{n_k} u_{k,j}\left|\frac{\bar f_{k,j} - f^\mathrm{mod}_{k,j}}{\langle \bar f_k\rangle}\right|^{2} \;-\; \phi_k\right) \;+\; N^\mathrm{NaN}_k ,
$$

where $u_{k,j}$ are optional per-element inverse-variance weights taken from the *data*
(normalised to mean 1, floored at $10\,\%$ of their own mean) and $N^\mathrm{NaN}_k$
penalises failed evaluations. The total is

$$
E \;=\; \sum_k w_k\, c_k ,\qquad w_k \;=\; \lambda\,\frac{R_k}{s_k^2},\qquad R_k = 1 + \kappa\min(z_k, 4).
$$

$s_k$ is the measured between-preparation relative SD of feature $k$ at 8 mM. Because
$c_k$ is a sum of squared relative residuals, $c_k/s_k^2$ is a $\chi^2$ in units of one
preparation-width, which is what makes terms of different physical dimension comparable.
The global constant $\lambda$ ($\approx 0.10$ per element) is chosen so the reweighted data
terms reproduce the previous objective's data subtotal, making campaigns comparable across
reweightings; the *reported* goodness of fit is therefore $c_k/s_k^2$ per feature, not
$w_k c_k$. $R_k$ is an ATP-relevance factor ($z_k$ = the standardised 2 mM/8 mM effect
size): a baseline error in an ATP-sensitive feature is later indistinguishable from an ATP
mechanism. It is **disabled ($\kappa = 0$, $R_k \equiv 1$) in every reported fit**, because
weighting the 8 mM fit by the 2 mM effect would leak the low-ATP data into the baseline
from which the low-ATP behaviour is predicted.

Two weights are not set by this rule. The force–velocity block is pinned at $40\,\%$ of the
seed objective, and no single term may exceed a $30\,\%$ share; both were imposed to stop
one smooth term from dominating the search (see the deadband below).

$\phi_k$ is a **noise floor**, set to $\phi_k = \chi^2_\mathrm{fl}\, n_k\, s_k^2$ with
$\chi^2_\mathrm{fl} = 1$, i.e. one preparation-width. Its purpose is structural. An
objective is a fixed budget, and a search will spend it wherever descent is smoothest —
which is whichever term is *already* best fitted. In a 14-hour four-profile campaign the
force–velocity block held $40\,\%$ of the objective while sitting at $\chi^2 = 0.10$
(a third of one preparation-width, well inside the disagreement between its own two source
curves) and absorbed $110\,\%$ of the total gain, while every feature the campaign existed
to improve went backwards. Re-pinning the weight does not fix this — any nonzero weight on
a smooth, nearly-solved term remains the cheapest descent — so the floor is a *deadband*,
not a deletion: a term inside its own measurement resolution stops paying, and starts
paying again the moment it drifts out, so it still defends itself against regression.
Applying it moved the objective from $2.4088$ to $0.9081$: $62\,\%$ of the cost had been
residual inside one preparation-width.

Three classes of term are not data fits and are exempt from reweighting and flooring:
(i) **physiological boundary constraints**, which cost zero inside an interval and rise
quadratically outside it, applied to ATP turnover per head, $V_\mathrm{max}$, the SRX
fraction, the attached fraction and the detached-pool fraction, so that the fit is held
inside independently measured physiology rather than being free to reach the data by any
route; (ii) **rate-constant boxes** on $k_2$, $k_{srd}$, $k_{msrd}$; and (iii)
**guardrails** penalising qualitatively wrong solutions (a spurious second re-stretch peak;
an unbounded cool-down overshoot).

The dedicated $k_\mathrm{tr}$ manoeuvre is scored in the reported set through one term,
its redevelopment rate (`ktrP_k80`); its trace is not. The staircase is not executed inside
the objective, but it is *not* a held-out test either: it was used to choose the
series-elastic exponent between budget-matched refits, so it is a model-selection set, and
only its traces on preparations never used for selection are out-of-sample. The genuinely
held-out data are the 2 mM records, which no 8 mM fit has seen.
The guardrail on the cool-down overshoot (`coolDownLS`) is in fact a data fit to one
preparation's raw trace and has since been retired.

### 2.6 Parameter estimation

The model has $\mathcal{O}(100)$ adjustable parameters but a singular-value analysis of
the sensitivity matrix over the slack/FV/$k_\mathrm{tr}$ feature set gives an
**11-dimensional identifiability ceiling**: adding parameters does not add constraint.
Optimisation is therefore organised as repeated low-dimensional searches rather than one
high-dimensional one (`Analyses/RestretchFeatureFit/optimizeFeatures.m`). Each round:

1. **draws** a subset of $\le 12$ parameters from an eligible pool, with a small
   compulsory set always included, so each sub-problem is tractable;
2. runs **`surrogateopt`** (global, bounded) on the drawn box, a fixed multiplicative
   window about the current value ($0.4$–$2.5\times$; $0.3$–$3\times$ for PCHIP knots);
   absolute physiological limits (`Model/parameterBounds.m`) enter as penalties in the
   objective, not as the search box;
3. **polishes** with `fminsearch`;
4. **accepts** the round only if it improves the incumbent, in which case the multipliers
   are baked into the running parameter set and written to disk.

After a fixed number of non-improving rounds a transient perturbation ("kick") of a *copy*
of the incumbent is retried on the same draw and kept only if it improves; it can never
leak into a later round's starting point. The incumbent is written to disk only inside
the accept branch, and the final snapshot is asserted to match the incumbent cost — a
guard introduced after an earlier run saved a stalled excursion that was worse than its own
seed.

Every round's converged $(\boldsymbol\theta, E)$ pair is retained, so the cloud of
solutions near the optimum can be examined for uniqueness: whether parameters reconverge
(identified) or different parameters reach the same cost (sloppy).

Two procedural rules were imposed after they were found to change conclusions.
**Budget-match before ranking:** the same comparison of a shape parameter reversed its
answer between a $0.5$ h and a $1.5$ h per-arm budget, because at the short budget a refit
has barely left its seed. **Never test a shape parameter by a slice, a grid, or even a
two-parameter profile:** a profile that re-optimised the obvious compensator at every value,
validated at the incumbent, and had a smooth correction was still wrong by $10\times$,
because the parameter was coupled to a second shape that the profile held fixed.

### 2.7 Objective noise floor

In earlier lineages the objective was numerically rough: it scattered by $19.5\,\%$ under
parameter changes of $10^{-6}$. The cause was traced (Analyses/StrainWindowNoise) not to
the strain grid — rolls delete at most $10^{-5}$ of attached mass — but to features read
off the solver's adaptive step points at a relative tolerance of $10^{-2}$, with tighter
tolerances stalling on the hard rectifier corner of the titin branch at isometric rest. With
the smooth rectifier, the tolerances and the fixed output grid of §2.4, the scatter is
$0.3\,\%$. Costs from before and after the fix are not comparable. The mechanism ablation
of §2.3.11 was run before the fix, with an acceptance threshold of $2\sigma = 0.110$
($\sigma = 0.055$, the larger of two phantom gains on the campaign's controls); mechanisms
worth less than about $12\,\%$ of the objective were undetectable by that design.

### 2.8 Software and reproducibility

All simulation and analysis code is MATLAB (R2023a). Parameter sets, protocol velocity
tables, feature targets and analysis outputs are versioned alongside the code; each
analysis is a self-contained folder holding its script, its results and its conclusions.
**[TODO: repository DOI / archive link on submission.]**

---

## 3. Limitations

*Drafted 2026-09-18 from five independent reviews. Most of this belongs in the Discussion;
it is parked here so it is not rediscovered late. Items marked **[OPEN]** are unresolved and
block submission. See [`REVIEW-corrections.md`](REVIEW-corrections.md) for evidence.*

### 3.1 The preparation, and what the ATP manipulation actually changes

**[OPEN, analysis planned] Inorganic phosphate accumulates in the preparation core, and
less so at low ATP.** The trabeculae are of order $140 \times 100\ \mu$m in cross-section
\citep{Beard2022} and are permeabilised, so every metabolite reaches the myofilaments only
by diffusion from the bath. Creatine phosphate and creatine kinase regenerate ATP from ADP
in place and hold ADP low, but nothing removes $\mathrm{P_i}$, which must diffuse out. An
order-of-magnitude steady radial estimate — $0.15$ mM myosin heads, turnover taken from the
model at each ATP level, $D_{\mathrm{P_i}} \approx 3\times10^{-10}\ \mathrm{m^2s^{-1}}$,
mean excess $qR^2/8D$ — gives $\approx 3$ mM mean excess $\mathrm{P_i}$ at 8 mM and
$\approx 2$ mM at 2 mM in a $150\ \mu$m preparation, rising steeply with thickness. Because
low ATP slows turnover, the low-ATP condition carries *less* $\mathrm{P_i}$, and since
$\mathrm{P_i}$ depresses isometric force with little effect on $V_\mathrm{max}$
\citep{Cooke1985, Ebus2001}, part of the measured force gain at 2 mM may be relief from
$\mathrm{P_i}$ rather than an effect of ATP. The same gradient logic applies to ATP itself
in the core at a 2 mM bath, and \citet{Ogut2003} show that creatine phosphate supply
matters for $k_\mathrm{tr}$ in trabeculae but not in myofibrils. The design of the analysis
that bounds this — a reaction–diffusion model of the cross-section coupled to the
cross-bridge model's own ATPase, with the $\mathrm{P_i}$ term re-enabled — and of the
experiments that decide it (cross-section per preparation; added $\mathrm{P_i}$ at fixed
MgATP; covariance of the effect with thickness) is in `Analyses/PiAccumulation/plan.md`.
Until it is run, the measured ATP effect is an *apparent* bath-ATP effect.

**Only one of the three metabolites was varied.** Heart failure raises $[\mathrm{P_i}]$ and
$[\mathrm{MgADP}]$ as well as lowering $[\mathrm{MgATP}]$. The model's $\mathrm{P_i}$ and ADP
terms are present but disabled (§2.3.5), and no $\mathrm{P_i}$ arm has been recorded.
The framing that survives this is positive: 2–8 mM MgATP is 10–100× above its catalytic
$K_m$, so an effect there that survives the $\mathrm{P_i}$ correction is diagnostic of a
*non-catalytic* mechanism in a way a $\mathrm{P_i}$ effect would not be.

**Sample size.** Force is $n = 3$ preparations; only $k_\mathrm{tr}$ is $n = 4$. Every
between-preparation SD carries two degrees of freedom, so its 95 % interval spans roughly
12-fold, and since every objective weight is $\lambda/s^2$, the objective's *composition* is
itself an estimate with large uncertainty. The empirical-Bayes moderation stabilises this
but assumes exchangeability and conditional independence of the per-feature variances,
neither of which strictly holds. Sex is 2 M / 1 F and is not analysed. The reversed-order
($2\to8$) arm is a single preparation, and is simultaneously the one needing the largest
rundown correction and the one that remains highest after it.

**The rundown correction rests on one bracket.** $\varphi$ is calibrated from a single
same-condition repeat in a single preparation and transferred. Two independent routes give
0.452 and 0.625 and the corrected effect is stable across that range, and the $\varphi$-free
continuity chain gives the same answer — but the *model form* (linear in activated time,
slope-proportional damage, complete inter-activation recovery) is established over one
17.5 % decrement. A second end-of-session repeat is the highest-value missing measurement.

**Sarcomere length is derived, not measured**, on the preparation the passive and slack
analyses lean on most; the later preparations do carry a recorded SL channel, which the
curation path discards. Since series compliance carries two of the model's strongest levers
and rundown adds ~35 % series compliance within a session, this is a real limitation and not
a bookkeeping one.

**[OPEN] One published study does not see our force effect at saturating Ca$^{2+}$.**
\citet{Rhodehamel2025} report maximum Ca$^{2+}$-activated tension unchanged across
2–8 mM ATP in permeabilised porcine myocardium, where we find force raised by about a third.
Three things bear on this and all three belong in the text rather than being omitted: their
tension–pCa relation shifts *rightward* with rising ATP, so at **submaximal** Ca$^{2+}$ their
data give more force at low ATP, in our direction; their cross-bridge detachment rate falls
to $\times 0.38$ at 2 mM, agreeing with our kinetic result and exceeding it; and the authors
state that ionic strength was not matched in their mechanical assays — the measurement at
issue. Two differences remain substantive: porcine ventricle is $\beta$-myosin dominant
against mouse $\alpha$-myosin, roughly fourfold apart in detachment kinetics, and their
maximal tension is ~22 kPa against our ~69 kPa. Their $V_\mathrm{max}$ also falls by about
half at 2 mM where \citet{Beard2022} found it unchanged, which is a genuine disagreement.

**The full force–velocity relation is borrowed** from \citet{Beard2022} — averages of
$6$–$12$ trabeculae per condition — and validated in-preparation at a single velocity on
two preparations. Its low-ATP limb is contradicted by our own in-preparation point
($-51\,\%$ against $-29\,\%$ at $2\ ML\,\mathrm{s^{-1}}$), so it is scored only at 8 mM;
at 2 mM the model's peak normalised power ($0.44$–$0.48$) is compared with, but not fitted
to, the borrowed curve ($0.52$). A full in-preparation force–velocity relation at both ATP
levels is the missing measurement.

**The two $k_\mathrm{tr}$ probes disagree by 24 % and diverge only at low ATP.** Both are
reported as measurements, and the model is compared with both. The slack-derived rate is
~5× more compliance-sensitive than the dedicated manoeuvre; since the low-ATP condition also
bears ~32 % more force, the contamination has the right sign to exaggerate the slowing.
No mechanistic claim in this paper rests on the difference between them.

### 3.2 The model

**Identifiability.** The model carries of order 100 adjustable parameters against an
11-dimensional identifiability ceiling from the calibration feature set. No parameter
confidence intervals or profile likelihoods are reported. Several specific degeneracies are
known and should be stated rather than discovered: the six lattice parameters carry only two
identifiable combinations, reducing operationally to *attachment is attenuated below*
$SL \approx 2.09\ \mu\mathrm{m}$; the stiffness exponent sits on a flat ridge between steep
walls and is locally unidentifiable along its own strongest direction; and the occupancy
ceiling refines to a value above the physically attainable bound fraction, so it acts as a
soft rescale of attachment rather than as a saturation.

**Two SRX legs refine to zero.** The reported configuration runs the super-relaxed pools as
a directed cycle. This is *non-identification*, not identification: reviving those legs
changes the objective by less than the search noise, so the data cannot distinguish the
directed cycle from the reversible one. The scheme should be described as four effective
pools with two unidentified reverse rates.

**SRX, IHM and the structural OFF state are treated as one variable.** The biochemical
(nucleotide-turnover) and structural (X-ray-ordered) definitions can be uncoupled, and the
structural evidence motivating the availability mechanism is not, strictly, evidence about
the biochemical pool the model carries.

**Numerical.** The strain grid is a fixed 60-bin window that is *rolled* when the
zero-strain point approaches a boundary, discarding the probability in the trailing bins;
that mass reappears as detached heads through the complement. It is at most $10^{-5}$ of
the attached mass per roll (§2.3.2).

**Two attached states** are used where the parent model carried three, and the
Ca$^{2+}$-dependence of the passive element is absorbed into fitted constants because these
protocols hold Ca$^{2+}$ constant. The thin-filament Ca$^{2+}$ layer (§2.3.12) is
implemented but not identified.

### 3.3 The objective and the identification procedure

**Correlated residuals are summed as though independent.** The scored features are extracted
from the same traces, share a single noisy per-preparation force anchor and a single rundown
factor, and in one case one feature's amplitude is taken from another feature by
construction. The effective number of independent constraints is therefore materially
smaller than the number of terms, and "$k$ of $n$ features within one SD" overstates the
number of independent successes.

**The noise-floor deadband makes the reported total uninterpretable as a goodness of fit.**
A term inside its band has exactly zero gradient, so nothing resists degrading a well-fitted
feature up to the deadband edge; at the reported operating point a majority of data terms
pay nothing, including the absolute force level — the quantity the ATP effect is largely
*about*. The unfloored per-feature $\chi^2$ table, not the floored total, is the honest fit
statistic, and the deadband should be described as a search-conditioning device in the
optimiser's loss rather than as a measure of fit.

**One weight is hand-pinned.** The force–velocity share is set by fiat rather than from
measured precision, at a value far above what its own spread would give; a share cap
polices the remainder (§2.5).

**The staircase is not a held-out prediction test** (§2.5); only the 2 mM records are
genuinely out-of-sample for the 8 mM fit.

**The optimiser is a heuristic with no convergence guarantee.** Reported optima should be
defended by multi-start reproducibility of both cost and parameters.

**Model selection is by one-at-a-time ablation** against that noisy objective, with an
acceptance threshold that makes any mechanism worth less than roughly a tenth of the total
undetectable by construction, and which is structurally blind to mechanisms that substitute
for one another. No information criterion, cross-validation or nested test is reported, and
an information criterion cannot in any case be computed from a deadbanded objective.

**The passive element is partly circular.** Its parameters are identified on a protocol
recorded after four activations, i.e. in the most rundown state; of the three arguments
offered for scaling it up in the active runs, only one is independent of the fit the passive
branch feeds.

### 3.4 Scope

The claim is reproduction, not unique identification. A single A2-detachment gate
reproduces most of the 2 mM behaviour, and an SRX-entry partner improves it; the data do not
exclude other combinations of comparable plausibility, and the 2 mM numbers in §1.6 come
from a one-dimensional sweep on the 8 mM set, not yet from a joint fit. The 0.2 mM arm is
not used. The model is calibrated at one Ca$^{2+}$ level, one temperature and one species,
and the rate constants are therefore not directly comparable to literature values obtained
under other conditions.

---

## Notes for the co-authors

1. **Species.** The poster and \citet{Beard2022} both say *mouse* right-ventricular
   trabeculae (Beard 2022: male C57BL/6J, 11–14 weeks). One internal analysis note says
   *rat*. The draft says mouse throughout — please confirm for the 2026 Baker-lab
   recordings specifically.
2. **Reference [3] year.** The passive-mechanics paper is *J Physiol* **2025**
   (doi:10.1113/JP288807, PMID 40846493), not 2026; the poster and abstract say 2026.
3. **Conclusion numbers.** The poster's "force $\approx +20\,\%$, $k_\mathrm{tr}$ $-45\,\%$"
   predates the rundown correction. §1.6 now quotes the pooled 2 mM target used to score
   the model (force $\times 1.32 \pm 0.09$, $n = 3$; slack $k_\mathrm{tr}$ $\times 0.57$;
   $k_\mathrm{tr}$ manoeuvre $\times 0.74$) and reports both $k_\mathrm{tr}$ probes as
   measurements (decision D1). The reconciliation analysis gives force $\times 1.34$–$1.36$
   and protocol $k_\mathrm{tr}$ $\times 0.72$; keep the two sets consistent in the final
   text.
4. **Beard 2022 numbers** (§1.2) were taken from the published paper (PMC9463691) by an
   automated read and should be checked against the PDF: $n = 12/6/6$ trabeculae at
   8/4/2 mM from 6 mice; peak power $38.3/29.9/23.5$ W/L; $k_\mathrm{tr}$
   $34.9/28.1/25.8\ \mathrm{s^{-1}}$. The repo's `LoadData.m` carries a different 8 mM
   isometric force (56.4 kPa) and $k_\mathrm{tr}$ (37.8 s⁻¹) — probably a different subgroup
   average; reconcile before quoting.
5. **Six states vs. what is actually running.** The reported parameter set has
   $k_{msr} = 0$ and $k_{srd} = 0$, so the SRX block runs as a directed cycle rather than
   two reversible exchanges. §2.3.6 says this explicitly.
6. **Which parameter set is described.** §2.3 now quotes `robust1_seed` (refit9 lineage),
   the core on which the 2 mM sweep was run. It differs from the earlier `floored_seed`:
   the lengthening-only slip/re-binding fluxes are off, $\tau_\mathrm{reg} = 0.2$ ms,
   $P_b^\mathrm{max} = 1.21$, and the smooth rectifier, extrapolation clamp and tight
   tolerances are on. Update when the joint 8 + 2 mM fit closes.
7. **Ca²⁺ layer** (§2.3.12) lives on branch `CaActivation` (worktree `ATP-wt-ca`), not yet
   merged; the MyoSim form and the EHT fits there are uncommitted.
8. **0.2 mM ATP.** The Myofilament abstract mentions a 0.2 mM arm; nothing in the current
   objective uses it. Left out of the draft — say if it should go back in.
9. **Missing bench detail** is marked **[TODO]** in §2.1 — solutions, temperature, pCa,
   nucleotide buffering, inhibitor concentrations, cross-sections (needed for §3.1), and $n$.
