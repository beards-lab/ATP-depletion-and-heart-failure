# A mechanistic, ATP-dependent cross-bridge model of cardiac muscle mechanics under energetic stress

**Draft — Introduction and Methods.**
Markdown with inline LaTeX, intended to be pasted into Overleaf. Math uses `$…$` / `$$…$$`;
citations use `\citep{}`/`\citet{}` against the accompanying `references.bib`.

> ## ⚠ UNDER REVISION — do not circulate
>
> Five independent reviews (2026-09-18) found errors that are **not yet fixed in this file**.
> Read [`REVIEW-corrections.md`](REVIEW-corrections.md) first. In particular:
>
> 1. **§1.6's two-gate claim is probably wrong.** It rests on the slack `ktr` (×0.55), which
>    the project's own index says is compliance-contaminated and not to be used for kinetics.
>    On the dedicated `ktr` protocol (×0.72) the product force × `ktr` ≈ 0.97 — turnover-
>    neutral — and the argument dissolves. The claim also appears nowhere in the repo.
> 2. **The phenomenon is not yet defended against diffusion/transport.** Force up, rate down,
>    V_max unchanged is also the signature of a rigor-loaded trabecula core.
> 3. **§1.4's H2 and H3 are Beard 2022's own published conjectures** and must be credited.
> 4. Nine factual errors in §2.3–2.7 about what the code does, and a retracted rundown
>    recommendation reproduced in §2.2.1. All listed in the corrections file.

> **Rendering convention — please keep to it when editing.** Display-math `$$` delimiters
> sit alone on their own line (opening and closing), and an inline `$…$` never spans a line
> break. Most Markdown renderers, GitHub's included, mis-pair the delimiters otherwise and
> everything downstream of the offending line renders as garbage.
Items needing a decision or a number from the bench are marked **[TODO]**.

---

## 1. Introduction

### 1.1 The energetic deficit of the failing heart is transduced by the cross-bridge cycle

The heart converts chemical free energy into mechanical work at a single molecular
interface: the cyclical interaction of myosin cross-bridges with actin, powered by the
hydrolysis of MgATP. Every beat, every unit of stroke work and every determinant of
relaxation passes through that interface, and the currency it spends is the free energy
of ATP hydrolysis, $\Delta G_\mathrm{ATP}$, set jointly by $[\mathrm{MgATP}]$,
$[\mathrm{MgADP}]$ and $[\mathrm{P_i}]$.

In heart failure that currency is devalued. Phosphorus magnetic resonance spectroscopy
shows a systematic fall in the myocardial phosphocreatine-to-ATP ratio — pooled
PCr/ATP $\approx 1.96$ in healthy volunteers against $1.54$ in heart failure with
reduced ejection fraction, a $21\,\%$ reduction, and $1.60$ in HFpEF — and the deficit
predicts mortality independently of ejection fraction \citep{Bray2026,
Neubauer1997}. Because the creatine kinase system buffers ATP at the expense of PCr,
falls in PCr precede falls in ATP; by the decompensated stage the adenine nucleotide
pool itself is depleted, cytosolic ADP and $\mathrm{P_i}$ rise, and $\Delta G_\mathrm{ATP}$
falls to a level at which ATP-consuming processes begin to be thermodynamically
constrained \citep{Lopez2020, Ingwall2002}.

The component of this system we address here is the last step in that causal chain: how
the concentration of MgATP at the myofilament modulates cross-bridge cycling, and
therefore sarcomere-level force, shortening velocity, power and the *kinetics* of force
development. Our target is a mechanistic, quantitative description of that modulation
that is identifiable from mechanical measurements — the mechanical half of an
ATP-supply/ATP-demand axis whose metabolic half already exists as a validated
mitochondrial model \citep{Collins2026, Bazil2016}.

### 1.2 What is currently thought to be known

Three bodies of work define the status quo.

**Nucleotide pharmacology of the cross-bridge.** In permeabilised fibres, the products
and substrate of hydrolysis act at identifiable steps. Added $\mathrm{P_i}$ reduces
isometric force approximately linearly in $\log[\mathrm{P_i}]$ while leaving $V_\mathrm{max}$
essentially unchanged, placing the main force-generating transition after $\mathrm{P_i}$
release; added MgADP raises isometric force and lowers $V_\mathrm{max}$, consistent with
competitive occupancy of the nucleotide site slowing the ATP-binding step that terminates
the working stroke \citep{Cooke1985, Gordon2000}. MgATP itself sets the rate of
cross-bridge dissociation, so the conventional reading of ATP depletion is *slower
detachment*: a larger attached fraction, higher isometric force, lower shortening
velocity, and reduced power.

**Direct measurement of the ATP effect in cardiac muscle.** \citet{Beard2022} measured
force–velocity relations in demembranated mouse right-ventricular trabeculae at
$8$, $4$ and $2$ mM MgATP and found that halving MgATP from $8$ to $4$ mM — the range
spanned by heart failure — produced a biologically substantial, statistically significant
fall in peak power and in the rate of tension development, with the attached cross-bridge
population *greater* at $2$ mM than at $8$ mM. More recently, \citet{Rhodehamel2025}
showed by small-angle X-ray diffraction that MgATP acts on the thick filament directly:
raising $[\mathrm{ATP}]$ shifts myosin heads from quasi-helically ordered OFF states on
the filament backbone to disordered ON states, accelerates both attachment and detachment
rates, and raises power and $V_\mathrm{max}$ — so $[\mathrm{ATP}]$ is not only a substrate
for the cycle but a regulator of how many heads are available to enter it.

**Thick-filament regulation and its state in disease.** Myosin heads partition between a
disordered-relaxed (DRX) pool available for actin binding and an energy-sparing
super-relaxed (SRX) pool with ATP turnover an order of magnitude slower
\citep{Hooijman2011, McNamara2015}. The partition is mechanosensitive: force borne by the
thick filament recruits heads out of the OFF/SRX population, which is now the leading
account of length-dependent activation and the Frank–Starling relation
\citep{Linari2015, Campbell2018, Irving2017}. In failing human myocardium the SRX
fraction is reduced and the DRX fraction raised, so resting myosin ATP consumption is
elevated precisely where oxidative capacity is lowest \citep{Ma2026SRX}.

**Computational models.** Strain-distribution (Huxley-type) cross-bridge models with
explicit chemistry are the established vehicle for this problem
\citep{Huxley1957, Razumova1999, Rice2008}. \citet{Tewari2016} built a metabolite-explicit
cardiac cross-bridge model identified against complex-modulus data across MgATP and
$\mathrm{P_i}$; \citet{Land2017} and \citet{Campbell2018} added series compliance and
force-dependent OFF$\to$ON recruitment; \citet{Lopez2020} coupled a cross-bridge model to
myocardial energetics to show that metabolite changes of the magnitude seen in
decompensated failure are sufficient to impair pump function.

### 1.3 The gap

The conventional account is coherent but it does not, in fact, reproduce the measurements.
The most direct test available is \citeauthor{Beard2022}'s own: their strain-discretised
model, in which $[\mathrm{MgATP}]$, $[\mathrm{MgADP}]$ and $[\mathrm{P_i}]$ enter the
transition rates through rapid-equilibrium binding factors, fits any *individual*
force–velocity curve well, but **could not capture the observed ATP-dependent reduction in
power** \citep{Beard2022}. The ATP dependence is therefore not currently *in* the
mechanism; it is absorbed by refitting parameters per condition, which is a description,
not an explanation, and cannot be transported into a multiscale energetics model.

Two things are missing.

*Mechanistically*, a model in which ATP enters only as a detachment substrate has too few
degrees of freedom. In a cycle with a single strongly-bound state, the isometric force,
the attached fraction and the rate of force redevelopment $k_\mathrm{tr}$ are welded
together by the turnover identity $\mathrm{ATPase} = r(1-r)\,k_\mathrm{tr}$ with $r$ the
duty ratio, so any change that raises force by slowing detachment necessarily changes
$k_\mathrm{tr}$ by a coupled amount. Whether the measured pair of changes is compatible
with that constraint is a testable question, and it has not been tested.

*Experimentally*, the steady force–velocity relation is a weak instrument. It is a
near-static summary of a cycle whose ATP sensitivity is kinetic, it cannot separate
attachment from detachment, and it contains no information about how quickly the muscle
re-recruits bridges after they are stripped. A dataset that constrains the mechanism has
to contain transients.

### 1.4 Hypotheses

We consider four non-exclusive accounts of the cardiac ATP effect.

- **H1 — product inhibition.** Low $[\mathrm{ATP}]$ acts through raised $[\mathrm{MgADP}]$
  occupancy, slowing the terminal detachment step (the classical reading,
  \citealp{Cooke1985}). *Scientific rationale:* it explains force up / velocity down in
  one parameter. *Weakness:* in a creatine-kinase/phosphocreatine-buffered permeabilised
  preparation, cytosolic $[\mathrm{MgADP}]$ is held near zero, so the driver at the
  myofilament should be $[\mathrm{ATP}]$ itself, not ADP.
- **H2 — substrate limitation of detachment.** The post-stroke head must bind ATP before
  it can release, so detachment carries a Michaelis factor in $[\mathrm{ATP}]$.
  *Rationale:* direct, and it is the step ATP actually occupies. *Weakness:* the reported
  $K_m$ for ATP-induced dissociation is well below $2$ mM, so on its own it predicts a
  much smaller effect than is measured.
- **H3 — thick-filament availability.** $[\mathrm{ATP}]$ sets the SRX$\leftrightarrow$DRX
  partition directly \citep{Rhodehamel2025}, so low ATP changes the *number of heads
  entering the cycle* and the gain of force-dependent recruitment, not only the rate at
  which they leave it. *Rationale:* it decouples the force leg from the kinetic leg, which
  a single detachment rate cannot do; and it is the mechanism for which there is new
  structural evidence.
- **H4 — the phenotype is kinetic and requires a transient probe.** The ATP effect is
  larger and better determined in rate observables (force redevelopment, post-re-stretch
  recovery) than in any steady-state quantity, so a model calibrated on force–velocity
  alone is under-determined in exactly the directions that matter.

### 1.5 Approach

We test these by building a cross-bridge model with enough structure to express them
separately, and by calibrating it against a dataset that can tell them apart.

The model is a six-pool, strain-discretised probabilistic framework: two detached
(disordered-relaxed) pools distinguished by nucleotide occupancy ($\mathrm{M\cdot ATP}$
and $\mathrm{M\cdot ADP\cdot P_i}$), two strain-resolved attached states (pre- and
post-power-stroke), and two super-relaxed pools with the same nucleotide distinction.
Making the nucleotide state explicit lets $[\mathrm{MgATP}]$, $[\mathrm{MgADP}]$ and
$[\mathrm{P_i}]$ modulate rates at the chemically correct steps throughout the cycle
rather than at one lumped detachment rate. Force-dependent SRX mobilisation supplies
thick-filament mechanosensing; asymmetric, piecewise strain-dependent modifiers on the
attached-state transitions supply the velocity dependence of the stroke and of
detachment; a Ca$^{2+}$-sensitive viscoelastic titin element supplies the passive
mechanics \citep{Jezek2025}.

The data are paired 8 mM and 2 mM MgATP recordings on the same demembranated
right-ventricular trabecula, each carrying four mechanical protocols — a five-cycle
slack/re-stretch train, a force-redevelopment ($k_\mathrm{tr}$) manoeuvre, a staircase of
stretch–hold steps, and in-preparation isovelocity ramps — plus a
blebbistatin/mavacamten passive run for the titin contribution. Calibration targets
*discrete features* of those records (rates, peaks, valleys, plateaux, delays) rather than
whole traces, with weights derived from the measured between-preparation spread so that
every term is a $\chi^2$ in units of one preparation-width.

Three specific questions follow. **(i)** Can one parameter set carry all four protocols
at 8 mM simultaneously — i.e. is the framework adequate before any ATP perturbation is
applied? **(ii)** Applied to that 8 mM-calibrated model, which nucleotide-dependent gate,
or which combination, reproduces the *measured* 2 mM/8 mM ratios of force and of rate?
**(iii)** Which observables actually identify that gate, and which are within the noise?

### 1.6 Principal findings

**[TODO — tighten once the low-ATP campaign closes; the statements below are the ones
currently supported.]**

A single parameter set reproduces the 8 mM behaviour of all four protocols, including
two that were never scored, with 13 of 20 features inside one between-preparation
standard deviation of a three-preparation pooled target. A systematic single-mechanism
ablation over 35 arms found nothing that improves that fit beyond the optimiser's own
noise, and four mechanisms that can be removed at no cost — so the residual at 8 mM is
not a missing mechanism.

The measured ATP effect, after correction for preparation rundown and pooled across four
protocol days, is **force $\times 1.36$ and $k_\mathrm{tr}\times 0.55$**: low ATP makes
cardiac muscle *stronger and slower*, with the product force $\times k_\mathrm{tr}$
falling $\times 0.74$. All six force *shape* features are ATP-flat, so the effect is a
scale and a slowing, not a reshaping. This is consistent with H2–H4 and, on its own,
inconsistent with the simplest form of H1, since the preparation is nucleotide-buffered.
It is also quantitatively inconsistent with a single detachment-rate change: the force and
the kinetic legs move by amounts that the two-attached-state turnover identity cannot
deliver together, so at least two gates — one on the cycle and one on head availability —
are required. **[TODO: state the final identified gate set and the residual.]**

---

## 2. Methods

### 2.1 Experimental preparation and protocols

#### 2.1.1 Preparation

Mechanical measurements were made on demembranated right-ventricular trabeculae from
adult mice, following the procedures of \citet{Beard2022}.
**[TODO — Baker lab to supply: strain, sex, age; skinning protocol and detergent;
activating and relaxing solution composition, ionic strength and pH; free
$[\mathrm{Ca^{2+}}]$/pCa of the activating solution; creatine kinase / phosphocreatine
buffering; bath temperature; how MgATP was varied at fixed ionic strength and fixed free
$[\mathrm{Mg^{2+}}]$; cross-sectional area measurement; whether sarcomere length was
measured optically or inferred from motor length.]**

Motor length is reported as *muscle length* $ML$, normalised such that $ML = 1.0$
corresponds to a nominal sarcomere length of $2.0\ \mu\mathrm{m}$; nominal sarcomere length is
$SL = 2\,ML$ throughout. Stress is reported in kPa. No optical sarcomere-length signal is
present in these recordings, so $SL$ in the data files is a derived quantity
(Analyses/PassiveRundownATP); the model therefore carries series compliance explicitly
(§2.3.9) rather than assuming $SL$ tracks $ML$.

Contractile activation was maintained throughout each run; the nucleotide condition
(8 mM or 2 mM MgATP) was changed between runs on the same preparation, so the two
conditions are paired within preparation and sequential in time.

#### 2.1.2 Mechanical protocols

Each run contained, in order:

1. **Force redevelopment ($k_\mathrm{tr}$).** From the isometric plateau at $ML=1.0$, a
   rapid release to $ML=0.8$, a brief hold to allow attached bridges to detach, a
   re-stretch to $ML=1.05$, a hold, and a return to $ML=1.0$, after which force
   redevelops isometrically. The redevelopment rate is the primary kinetic observable.
2. **Slack/re-stretch train.** Five cycles of rapid shortening (slack release) from
   $ML=1.1$, an unloaded hold, and a re-stretch at a cycle-specific velocity, with the
   re-stretch amplitude forming a ladder across cycles. Cycles 1–4 re-stretch to
   $ML=1.10$; cycle 5 re-stretches to $ML=1.00$ by the smallest amount and at the slowest
   ramp, and is therefore excluded from any amplitude-graded quantity (§2.2.3).
3. **Staircase.** A sequence of constant-velocity stretch steps separated by isometric
   holds, used here as a held-out prediction test (§2.5).
4. **Isovelocity ramps (FV2).** In-preparation shortening ramps at $2\ ML\,\mathrm{s^{-1}}$
   providing a force–velocity point measured on the same trabecula as the transients.
5. **Passive run.** After the active runs, the preparation was inhibited with
   *N*-benzyl-*p*-toluenesulphonamide (PNB) and mavacamten and the length protocol
   repeated, isolating the passive (titin-dominated) stress. **[TODO: confirm PNB and
   mavacamten concentrations and incubation.]**

The full-range force–velocity relation at $8$, $4$ and $2$ mM MgATP is taken from the
isovelocity dataset of \citet{Beard2022} (velocities $0$, $0.5$, $1$, $2$, $3$, $4$, $5$,
$6\ ML\,\mathrm{s^{-1}}$), with a second, independently fitted Baker-laboratory curve used
to bound its reproducibility. The in-preparation FV2 ramps validate the borrowed 8 mM
curve to $0.21$ SD but contradict its 2 mM limb, so the borrowed curve is used only at
8 mM (Analyses/DataPooling).

#### 2.1.3 Datasets

Four protocol days are used: three ran $8\to2$ mM and one ran $2\to8$ mM, the reversed
order being deliberate so that the rundown bias (§2.2.1) and the true ATP effect can be
separated. **[TODO: give the $n$, the animal identities and the acquisition dates in the
form intended for publication; internally these are `03/27/2026 M`, `04/03/2026 F`,
`04/10/2026 M`, plus the legacy Baker slack set.]**

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
decay phase $T_\mathrm{dec}$ replaces a common activation duration. Together these reduce
the between-preparation CV of the ATP force ratio from $20\,\%$ to $7\,\%$ at no cost in
assumptions.

Force is corrected; $k_\mathrm{tr}$ is **not**. The coupling between the two corrections
was scanned with the coefficient left free and every physically motivated value made
cross-preparation consistency worse ($\mathrm{CV}\ 9.4\,\% \to 15.8\,\%$), because the
underlying lesion — loss of parallel force-generating material ($\approx 16\,\%$) plus
added series compliance ($\approx 35\,\%$) — reaches force and $k_\mathrm{tr}$ through
different channels. Full derivation, validation and limitations in
`Analyses/RundownCorrection/methods.md`.

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
simulation, except that experimental traces are pre-smoothed
(Savitzky–Golay, $11$ ms, for the rapid re-stretch peak; a $5$ ms moving median to locate
the undershoot; robust LOESS, $25$ ms, for the recovery tail) whereas simulation output is
already smooth.

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
$17.2\,\%$ to $6.1\,\%$. And segment-averaged quantities use slack cycles 1–4 only: cycle 5
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

and $D_T$ is carried as the complement of the others rather than integrated, which removes
one redundant equation and enforces conservation exactly:

$$
D_T \;=\; \max\!\Big(0,\; 1 - \textstyle\int p_1 - \int p_2 - D_D - S_T - S_D\Big).
$$

Five further state variables describe the mechanics and the slow availability processes:
sarcomere length $SL$, series-elastic length $L_\mathrm{SE}$, a titin viscous stress $X_v$,
a registration-availability state $A_\mathrm{reg}$ and a compliant-realignment deficit
$A_\mathrm{def}$ (§2.3.8–2.3.10). A third attached state $p_3$ is implemented and may be
enabled (`NumberOfStates = 3`); the results reported here use the two-attached-state form.

#### 2.3.2 Strain coordinate and the moving frame

Strain is discretised on a uniform grid $s^{(g)}_i$ of $n_s = 60$ bins with spacing
$\mathrm{d}s = 2$ nm, spanning $\pm 59$ nm. Rather than adding an advection term
$-v\,\partial p/\partial s$, the model advects the *coordinate*: the physical strain
carried by grid bin $i$ at time $t$ is

$$
s_i(t) \;=\; s^{(g)}_i \;-\; \tfrac{1}{2}\big(L_\mathrm{piv} - L_c(t)\big), \qquad L_c = SL - L_\mathrm{SE},
$$

with $L_\mathrm{piv}$ the reference contractile-element length. Each bin's strain therefore
changes at $\dot s = \tfrac12 \dot L_c$ — half the sarcomere rate, because a half-sarcomere
takes half the length change — which is the Lagrangian form of filament sliding and is
exact by construction, with no numerical diffusion. Shortening carries attached bridges to
negative strain; lengthening carries them positive.

Attachment injects probability at $s = 0$, whose *grid index* therefore moves. Mass is
distributed over a normalised triangular kernel of half-width $A_x = 6$ nm centred on
$s=0$,

$$
K(s) \;=\; \frac{\max\!\big(0,\,1 - |s| / A_x\big)}{\displaystyle\int \max\!\big(0,\,1-|s'|/A_x\big)\,\mathrm{d}s'},
$$

which represents the spread of accessible actin target positions and keeps the injection
smooth as the zero-strain point crosses bin boundaries. When that point comes within $A_x$
of a grid edge, an ODE event terminates the step and the grid is shifted and extended
(`UseSpaceExtension`), so no probability is lost at the boundary.

#### 2.3.3 Governing equations

With $R$ denoting fluxes (strain-resolved where written as functions of $s$):

$$
\begin{aligned}
\frac{\partial p_1}{\partial t} &= R_{D1}\,K(s) \;-\; R_{1D}(s) \;-\; R_{12}(s) \;+\; R_{21}(s),\\[2pt]
\frac{\partial p_2}{\partial t} &= R_{12}(s) \;-\; R_{21}(s) \;-\; R_{2}(s) \;+\; R_\mathrm{re}(s),\\[2pt]
\frac{\mathrm{d}D_D}{\mathrm{d}t} &= R_{TD} - R_{DT} - R_{D1} + \int R_{1D}\,\mathrm{d}s + R_{S_D\to D_D} - R_{D_D\to S_D},\\[2pt]
\frac{\mathrm{d}S_T}{\mathrm{d}t} &= R_{D_T\to S_T} - R_{S_T\to D_T} + R_{S_D\to S_T} - R_{S_T\to S_D},\\[2pt]
\frac{\mathrm{d}S_D}{\mathrm{d}t} &= R_{S_T\to S_D} - R_{S_D\to S_T} - R_{S_D\to D_D} + R_{D_D\to S_D}.
\end{aligned}
$$

The detachment flux $R_2$ returns heads to $D_T$; because $D_T$ is the complement,
this requires no explicit term. Hydrolysis and its reverse connect the detached pools,

$$
R_{TD} = k_{ah}\,D_T, \qquad R_{DT} = k_{amh}\,D_D ,
$$

and attachment is

$$
R_{D1} \;=\; k_a^\mathrm{eff}\; N_\mathrm{ov}\; f_\mathrm{lat}\; f_\mathrm{occ}\; D_D .
$$

$R_\mathrm{re}$ is an optional lengthening-only re-attachment/hop flux on $p_2$ (§2.3.11).

#### 2.3.4 Strain dependence of the attached-state transitions

The four attached-state rates carry dimensionless, strain-dependent modifiers
$\phi_\bullet(s)$ built as shape-preserving piecewise cubic (PCHIP) interpolants through a
small set of fitted knots, which lets the optimiser shape each rate profile without
committing to an exponential or Gaussian form:

$$
\begin{aligned}
R_{1D}(s) &= k_d\,\phi_d(s)\,p_1(s) && \text{A1 detachment}\\
R_{12}(s) &= k_1\,\phi_1(s)\,p_1(s) && \text{power stroke}\\
R_{21}(s) &= k_{-1}\,\phi_{-1}(s)\,p_2(s) && \text{reverse stroke}\\
R_{2}(s)  &= k_2\,\phi_2\!\big(s + \Delta_2 - \Delta\big)\,g_\mathrm{ATP}(s)\,p_2(s) && \text{A2 detachment}
\end{aligned}
$$

All modifiers are clamped at $10^4\ \mathrm{s^{-1}}$. The knot sets are fitted parameters
(`PieceWiseStrainDep*X` / `*Params`); at the reported solution $\phi_1$ rises steeply on
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

Two gates are used here. The **cycle gate** multiplies the A2 detachment rate
($K_{T2}$): a post-stroke head must bind ATP before it can release, so low ATP retains
strongly-bound bridges — this is the force lever (H2). The **availability gate**
multiplies the SRX-ADP return $S_D \to D_D$ ($K_\mathrm{srx}$): re-priming a folded head
requires ATP, so low ATP slows the re-supply of attachment-competent heads — this is the
kinetic lever, and isometric force is nearly flat in it (H3). Product terms for
$[\mathrm{MgADP}]$ (competition at the nucleotide site, $K_D = 0.194$ mM;
\citealp{Yamashita1994}) and $[\mathrm{P_i}]$ (reversal of the force-generating step,
$K_{P_i} = 4.007$ mM) are implemented on the same normalised pattern and are inert in the
buffered condition used here ($[\mathrm{MgADP}] = [\mathrm{P_i}] = 0$).
**[TODO: confirm the nominal free $[\mathrm{ADP}]$ and $[\mathrm{P_i}]$ of the activating
solution; if they are non-zero, state them and enable the corresponding terms.]**

#### 2.3.6 Thick-filament regulation

Entry to and exit from the super-relaxed pools are gated exponentially by the total
sarcomere stress $F$, which is the model's representation of thick-filament
mechanosensing \citep{Linari2015, Campbell2018}:

$$
\begin{aligned}
R_{D_T\to S_T} &= k_{sr0}\,e^{-F/\sigma_2}\,D_T, &\qquad R_{S_T\to D_T} &= k_{msr}\,e^{+F/\sigma_1}\,S_T,\\
R_{D_D\to S_D} &= k_{srd}\,e^{-F/\sigma_{srd2}}\,D_D, &\qquad R_{S_D\to D_D} &= k_{msrd}\,e^{+F/\sigma_{srd1}}\,g^{(\mathrm{srx})}_\mathrm{ATP}\,S_D,\\
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
an outcome of identification, not an assumption — the reverse legs remain in the equations
and were re-tested with their gates artificially revived, which changed the objective by
less than the optimiser noise (Analyses/MechanismAblation).

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

Finally, a mean-field occupancy ceiling attenuates the whole attachment flux as the bound
fraction approaches a physical limit,

$$
f_\mathrm{occ} = \max\!\left(0,\; 1 - \frac{\int (p_1+p_2)\,\mathrm{d}s}{P_b^\mathrm{max}}\right).
$$

This is deliberately *global*: the strain axis is not a binding-site axis, so a per-bin
occupancy term would be both a category error and grid-spacing dependent
(Analyses/BindingSiteOccupancy).

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

driven by the *imposed* velocity so that series-elastic ringing cannot corrupt it. Using
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
\dot A_\mathrm{def} = k_\mathrm{cr}\,\mathcal{L}\big(-\dot P_b\big)\,\frac{F}{F+F_\mathrm{cr}} \;-\; \frac{A_\mathrm{def}}{\tau_\mathrm{cr}},\qquad
P_b = \int (p_1+p_2)\,\mathrm{d}s ,
$$

with $\mathcal{L}(u) = u/(1+e^{-u/\epsilon})$ a smooth (SiLU) rectifier — a hard
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
$k_a^\mathrm{eff} = k_a\,A_\mathrm{reg}\,(1 - A_\mathrm{def})$.

#### 2.3.9 Force and the mechanical boundary condition

Cross-bridge stress is the first moment of the attached distributions about each state's
own equilibrium, with an exponent $e$ on the strongly-bound state allowing a non-Hookean
bridge, and separate stiffnesses for positive and negative strain:

$$
F_\mathrm{XB} \;=\; k_{s1}\!\int (s+\Delta_1)\,p_1\,\mathrm{d}s \;+\; k_{s2}\!\int \mathrm{sgn}(s+\Delta)\,|s+\Delta|^{\,e}\,p_2\,\mathrm{d}s ,
$$

where, when `UseNegativeKstiff` is set, the integrals are split at the zero-force strain
and the negative branch uses $k_{s1}^-$, $k_{s2}^-$. $\Delta$ is the power-stroke size and
$\Delta_1$ an offset that lets the pre-stroke state bear force, decoupling force duty from
the $p_2$ detachment that otherwise welds $k_\mathrm{tr}$ to the force–velocity relation.

The contractile element acts in series with an elastic element (end compliance plus
damaged-region compliance; $SL$ is not the commanded length,
cf. \citealp{deTombe2016}). Measured stress is generated by that element,

$$
F \;=\; \mathrm{sgn}(L_\mathrm{SE})\,\big|\max(F_\mathrm{min},\, k_\mathrm{SE} L_\mathrm{SE})\big|^{\,e_\mathrm{SE}} \;+\; \mu_2 v,
$$

and the series length evolves by a viscously regularised force balance,

$$
\dot{SL} = v(t), \qquad
\dot L_\mathrm{SE} = v(t) - \frac{F - F_\mathrm{tot}}{\mu^{\pm}},\qquad
F_\mathrm{tot} = F_\mathrm{XB} + F_\mathrm{pas},
$$

where $\mu^{\pm} = \mu$ when $F > F_\mathrm{tot}$ and $\mu_\mathrm{neg}$ otherwise (the
preparation is stiffer in extension than in recoil). In the limit $\mu\to0$ this is the
quasi-static balance $F = F_\mathrm{tot}$; a finite $\mu$ both regularises the
differential-algebraic structure for the stiff solver and represents the internal
viscosity that limits unloaded shortening \citep{deTombe1992}. $F_\mathrm{min}$ bounds
how negative the series element may be driven during a slack release.

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
shortening pays out slack and releases the stored unfolding stress. Written as one
continuous expression, it vanishes smoothly as $\dot L_c \to 0$, so the solver sees no
switch. Ca$^{2+}$-dependence of the passive limb — PEVK–actin attachment, proximal
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
and verified inert when off. Those active in the reported configuration are: a
lengthening-only re-attachment/hop flux on $p_2$ (`UseA2Reattaching`, which fills the
post-re-stretch valley and caps the first re-stretch peak by moving high-strain heads back
to the attachment point while conserving mass); and the tension-only Maxwell discharge
above. Those tested and rejected include a detachment ceiling with a nucleotide-free
refractory pool, catch-bond detachment suppression during stretch, re-routing of the
detachment flux into SRX, and fixed-dwell (occupancy-saturating) pool exits; each is
documented with its falsifying observable in `Analyses/`. A single-mechanism ablation
over 35 arms, each given a scalar calibration plus a budget-matched re-optimisation with
its own compensators made compulsory, found no addition worth its parameters at the
$\sigma = 0.055$ noise level of the campaign, and four mechanisms removable at no cost.

### 2.4 Numerical implementation

The system is integrated with MATLAB's `ode15s` (variable-order BDF, appropriate to the
stiffness introduced by the fast strain-dependent rates) with absolute tolerance
$10^{-4}$ and relative tolerance $10^{-2}$. Each protocol is supplied as a velocity table:
a list of $(t_k, v_k)$ breakpoints defining constant-velocity segments, integrated
segment-by-segment with the state carried across, so commanded length changes are exact
and no interpolation of the drive signal is required inside the right-hand side. An event
function terminates a step when the zero-strain point approaches a grid boundary (grid
shift/extension, §2.3.2) and a second event enforces a wall-clock ceiling per simulation,
so a pathological parameter draw during optimisation fails fast rather than stalling the
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
sets are serialised to plain `.m` snapshots with full precision (`%.17g`), so every
reported result is reproducible from a file.

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
$c_k$ is a sum of squared relative residuals, dividing by $s_k^2$ makes every term a
$\chi^2$: **$w_k c_k = 1$ means the model is one preparation-width from the data**, which
is what makes terms of different physical dimension comparable. $R_k$ is an ATP-relevance
bonus ($z_k$ = the standardised 2 mM/8 mM effect size, capped at 4): the fit is performed
at 8 mM, but its purpose is to be the platform on which the ATP perturbation is read, so a
baseline error in an ATP-sensitive feature is later indistinguishable from an ATP
mechanism, whereas an error in an ATP-flat feature is not. $\lambda$ is one global constant
chosen so the reweighted data terms reproduce the previous objective's data subtotal,
making campaigns comparable across reweightings.

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

Two protocols — the staircase and the dedicated $k_\mathrm{tr}$ manoeuvre — are run but
were historically unscored, and are therefore used as **held-out prediction tests**.

### 2.6 Parameter estimation

The model has $\mathcal{O}(100)$ adjustable parameters but a singular-value analysis of
the sensitivity matrix over the slack/FV/$k_\mathrm{tr}$ feature set gives an
**11-dimensional identifiability ceiling**: adding parameters does not add constraint.
Optimisation is therefore organised as repeated low-dimensional searches rather than one
high-dimensional one (`Analyses/RestretchFeatureFit/optimizeFeatures.m`). Each round:

1. **draws** a subset of $\le 12$ parameters from an eligible pool, with a small
   compulsory set always included, so each sub-problem is tractable;
2. runs **`surrogateopt`** (global, bounded) on the drawn box, whose bounds are derived
   from absolute physiological bounds (`Model/parameterBounds.m`) mapped to multiplicative
   windows about the current value;
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

The objective is numerically rough. Repeated evaluation under identical settings gives
$\sigma \approx 0.076$ with a peak-to-peak spread of $0.25$ at 2 mM, arising from
adaptive-step and event-detection non-determinism, so improvements are accepted against a
threshold of that order and never against a difference in the fourth decimal. In the
mechanism-ablation campaign the empirical threshold was calibrated from the campaign's own
controls, both of which re-optimised to exactly $0.9081$, giving $\sigma = 0.055$.

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

**[OPEN] Diffusive supply is not excluded as the cause of the effect.** The trabeculae are
of order 140 µm wide and 100 µm thick, maximally activated, with the ATP pool buffered by
creatine phosphate and creatine kinase. Lowering bath [MgATP] fourfold lowers diffusive
supply to the core fourfold at unchanged demand, and the observed signature — force up,
rate down, $V_\mathrm{max}$ unchanged — is also what a partially rigor-loaded core with a
well-supplied shell would produce. \citet{Ogut2003} report that lowering creatine phosphate
at *constant* MgATP depresses $k_\mathrm{tr}$ by half in cardiac trabeculae, with no such
effect in myofibrils, where there is no diffusion distance. Until core [ATP] is estimated
from the measured ATPase and cross-section, or the effect is shown not to covary with
trabecula thickness, a transport explanation cannot be excluded. This is the single largest
threat to the interpretation.

**Only one of the three metabolites was varied.** Heart failure raises $[\mathrm{P_i}]$ and
$[\mathrm{MgADP}]$ as well as lowering $[\mathrm{MgATP}]$, and $\mathrm{P_i}$ has far
greater force authority than anything ATP does between 2 and 8 mM. The model's
$\mathrm{P_i}$ and ADP gates are described but **not implemented in the present
right-hand side**, and no $\mathrm{P_i}$ arm has been recorded. The defensible framing is
positive rather than apologetic: 2–8 mM MgATP is 10–100× above its catalytic $K_m$, so a
mechanical effect there is diagnostic of a *non-catalytic* mechanism in a way a
$\mathrm{P_i}$ effect would not be. A $\mathrm{P_i}$ arm is the next experiment.

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

**The full force–velocity relation is borrowed** from a separate dataset and validated
in-preparation at a single velocity. Its low-ATP limb is contradicted by our own
in-preparation ramps, so it is used only at 8 mM — which means the one protocol with full
velocity coverage is unavailable for the condition the paper is about.

**[OPEN] The two $k_\mathrm{tr}$ probes disagree by 24 % and diverge only at low ATP.** The
slack-derived rate is ~5× more compliance-sensitive than the dedicated manoeuvre; since the
low-ATP condition also bears ~36 % more force, the contamination has the right sign to
manufacture an apparent kinetic effect. Which probe is reported changes the headline
conclusion (§1.6) and must be settled before submission.

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
zero-strain point approaches a boundary, discarding the probability in the vacated bins;
that mass reappears as detached heads through the complement. The leakage is believed small
at the operating point but is not currently quantified or reported.

**Two attached states** are used where the parent model carried three, and the
Ca$^{2+}$-dependence of the passive element is absorbed into fitted constants because these
protocols hold Ca$^{2+}$ constant.

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
polices the remainder. The stated rule $w = \lambda R/s^2$ therefore does not describe the
shipped objective, and the ATP-relevance factor $R$ is set to unity in all reported fits.

**[OPEN] The staircase is not a held-out prediction test.** The series-elastic exponent
carried by the reported parameter set was selected by ranking budget-matched refits on the
staircase. A dataset used for model selection cannot also be reported as a prediction. The
staircase *traces* on preparations the model never saw remain weak out-of-sample evidence;
the dedicated $k_\mathrm{tr}$ manoeuvre is genuinely held out in this lineage, but in the
successor lineage its rate is a fenced term while only the trace stays held out.

**The optimiser is a heuristic with no convergence guarantee**, and its accept rule is
accept-if-better against an objective whose irreproducibility arises from a wall-clock
integration guard — a mechanism that biases the cost *downward* on truncation, i.e. exactly
the direction accept-if-better selects for. Reported optima should be re-evaluated after
the fact and defended by multi-start reproducibility of both cost and parameters.

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

The mechanistic identification of the ATP effect is **not closed**: the low-ATP campaign is
in progress, and the claim that the measured force and rate changes cannot arise from a
single gate is **[OPEN]** — no fit of the low-ATP data with a one-gate model and its
residual has yet been run. The 0.2 mM arm is not used. The model is calibrated at one
Ca$^{2+}$ level, one temperature and one species, and the rate constants are therefore not
directly comparable to literature values obtained under other conditions.

---

## Notes for the co-authors

1. **Species.** The poster and \citet{Beard2022} both say *mouse* right-ventricular
   trabeculae (Beard 2022: male C57BL/6J, 11–14 weeks). One internal analysis note says
   *rat*. The draft says mouse throughout — please confirm for the 2026 Baker-lab
   recordings specifically.
2. **Reference [3] year.** The passive-mechanics paper is *J Physiol* **2025**
   (doi:10.1113/JP288807, PMID 40846493), not 2026; the poster and abstract say 2026.
3. **Conclusion numbers.** The poster's "force $\approx +20\,\%$, $k_\mathrm{tr}$ $-45\,\%$"
   predates the rundown correction. The settled figures are **force $\times 1.36$**
   (CV 11 % after correction, $n=4$) and **$k_\mathrm{tr}\times 0.55$** from the slack
   train — but $\times 0.72$ from the dedicated $k_\mathrm{tr}$ protocol, which is the
   less compliance-contaminated probe. The draft uses $\times 1.36$ / $\times 0.55$ and
   flags the discrepancy; decide which to headline.
4. **Six states vs. what is actually running.** The reported parameter set has
   $k_{msr} = 0$ and $k_{srd} = 0$, so the SRX block runs as a directed cycle rather than
   two reversible exchanges. §2.3.6 says this explicitly. If you would rather present the
   model as fully reversible, the text needs changing, not the code.
5. **0.2 mM ATP.** The Myofilament abstract mentions a 0.2 mM arm; nothing in the current
   objective uses it. Left out of the draft — say if it should go back in.
6. **Missing bench detail** is marked **[TODO]** in §2.1.1, §2.1.2 and §2.3.5 — solutions,
   temperature, pCa, nucleotide buffering, inhibitor concentrations, and $n$.
