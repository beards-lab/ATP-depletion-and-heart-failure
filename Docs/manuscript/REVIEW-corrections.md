# Review findings — corrections required before `01-introduction-and-methods.md` circulates

Five independent reviews, 2026-09-18. Four completed; the citation audit failed partway
(account spend limit) and is marked incomplete below.

**Legend.** ✅ = I verified this myself against the repo or a primary source in this session.
↪ = relayed from a reviewer, not independently re-checked. ⚠ = a reviewer over-claimed;
corrected here.

---

## 0. Three decisions that block everything else

Nothing downstream is worth editing until these are settled, because each one changes what
the Introduction is allowed to say.

### D1. Which `ktr` probe headlines — and does the two-gate claim survive?

§1.6 says force ×1.36, `ktr` ×0.55, product ×0.74, and concludes the effect is
"quantitatively inconsistent with a single detachment-rate change", so "at least two gates
… are required".

- ✅ `Analyses/RundownCorrection/conclusions.md` §5f: *"the slack `ktr` is the right probe
  for **rundown** and the **wrong** probe for cross-bridge kinetics. For the ATP effect
  prefer the ktr-protocol value ×0.72 (CV 4 %) over the slack value ×0.58 (CV 7 %)."*
- ✅ Per-preparation: protocol 0.746 / 0.732 / 0.691; slack 0.591 / 0.616 / 0.537.
- ✅ The two probes **diverge only at low ATP** (slack/protocol 1.09/1.45/0.94 at 8 mM,
  0.86/1.22/0.73 at 2 mM) — the compliance-contamination signature, with the right sign to
  manufacture the effect (2 mM bears ~36 % more force → more series extension → a
  5×-more-compliance-sensitive probe reads slower).
- ✅ Algebra: with `ktr` = f+g and duty r = f/(f+g), **force × ktr ∝ f**, the *attachment*
  rate. An invariant product is what a pure detachment gate predicts.
- ✅ On the protocol probe: 1.34 × 0.72 ≈ **0.97** (per-prep 0.956 / 0.970 / 0.996). The
  departure from unity that is the whole evidence for a second gate disappears.
- ✅ The "two gates" claim appears **nowhere in the repo**. `LowATP_k2Frontier`,
  `LowATP_Mechanisms`, `LowATP_3rdState` and `ATPDependenceFrontier` each contain a
  `results/` directory and **no `conclusions.md`**. `ROADMAP.md` still lists the
  mechanistic identification as open.

**Decision:** recompute on the protocol probe. If the product is ≈0.97, withdraw the
two-gate claim; the paper is then about why the effect is turnover-neutral yet
redistributes force and rate — still interesting, different paper. Do not submit with the
claim resting on the probe the project's own index tells readers not to use.

### D2. Can the phenomenon be defended against transport?

- ✅ Beard 2022 Methods: trabeculae **137 ± 15 µm wide × 97 ± 7 µm thick**, maximally
  activated, **12 mM creatine phosphate + 10 IU/mL creatine kinase**.
- ✅ Ogut & Brozovich, *Circ Res* 2003;93:54–60: CP 25 → 2.5 mM **cut `ktr` 51 % at constant
  5 mM MgATP**; per the abstract-level summary retrieved, lowering ATP to 1 mM or raising
  ADP to 1 mM did *not* significantly decrease `ktr`. **Read the full paper before building
  on that second clause — it was not verified against the text.**
- Myofibril rebuttal (no diffusion distance, no CP effect): *Pflügers Arch* 2006, PMID 16395599.

Force up / rate down / V_max unchanged is the signature of a partially rigor-loaded core
with a well-supplied shell. §1.4 currently dismisses H1 by *assuming* the buffering works,
while §2.1.1 admits the solution composition is unknown.

**Cheapest defences, in order:** a reaction–diffusion estimate of core [ATP]/[ADP] from
measured ATPase and cross-section; a CP titration at fixed MgATP; and the decisive one —
whether the ATP effect covaries with trabecula thickness across preparations
(cross-sections may already exist).

**Reframe that turns the related Pi objection into the rationale:** Pi and ADP act where an
effect is *expected*; 2–8 mM MgATP is 10–100× above its K_m, so an effect there is
diagnostic of a *non-catalytic* mechanism in a way a Pi effect would not be. Belongs in §1.1.

### D3. How to credit Beard 2022's own conjectures

✅ Verified verbatim from the Beard 2022 Discussion:

> "Another possible mechanism by which MgATP concentration may affect cross-bridge dynamics
> is through the transitions between the super-relaxed (SRX) and non-super-relaxed (non-SRX)
> states of the myosin heads."

> "A mechanism whereby the lower MgATP concentrations are associated with a greater
> stabilization of the SRX state may also be consistent with the observed effects of
> [MgATP] on force development kinetics."

> "it may be that the detachment kinetics or strain dependence of the detachment kinetics,
> and how these vary with changes in [MgATP], are not being appropriately formulated with
> the current model."

Those are H3 and the alternative to H2 — both presented in §1.4 as ours. D.A.B. is a
co-author. ✅ Beard 2022 Appendix A also states its cross-bridge model "integrates elements
of the model of Tewari et al. with the myosin recruitment model of Campbell et al." — so
metabolite-dependent rates *and* force-dependent OFF→ON recruitment were already present in
the model being extended.

**Fix (one sentence, and it strengthens the paper):** *Beard et al. proposed, without
testing it, that [MgATP] might act through the SRX↔non-SRX transition; the present work
makes that conjecture quantitative and tests it against data it was not fitted to.*

---

## 1. Factual errors in the model description (§2.3–2.7)

All ✅ verified against `Model/dPUdT_CombinedTransitions.m` and `params/floored_seed.m`.

| # | Section | Claim as written | Reality | Fix |
|---|---|---|---|---|
| 1.1 | §2.3.5 | MgADP and Pi terms "implemented on the same normalised pattern" | Lines 344–346 commented out; line 352 hardcodes `g1=1; g2=1; f1=0; f2=1`. Those factors exist only in the legacy, unused `dPUdTCa.m` | Delete the sentence |
| 1.2 | §2.3.5 | "Two gates are used here" | `UseAtpDetachR2 = 0`, `UseAtpKmsrd = 0` in the described set | "Two gates are **available**; both are disabled in the 8 mM reference set (where they are identically 1 by construction) and enabled only for the low-ATP arm" |
| 1.3 | §2.3.11 + §2.3.3 | SRX re-routing listed under "tested and rejected" | `SRXFromR2HighStrain = 0.212905` — **live**. The uniform variant `SRXFromR2 = 0` is the one that is off | Move to the active list; add the flux to dS_T/dt |
| 1.4 | §2.3.2 | grid "shifted and extended … no probability is lost at the boundary" | `evaluateModel.m` rolls a fixed 60-bin window by 20 bins and **zeroes** the vacated bins. 40 nm of strain discarded; mass reappears as detached heads via the complement | State the roll, state the discard, quantify the leakage |
| 1.5 | §2.3.7 | occupancy ceiling "as the bound fraction approaches a physical limit" | `P_bound_max = 1.13467` — 113 % occupancy, against the code's own stated 0.12–0.15. Acts as a ~0.82 soft rescale | Report the fitted value and say it is above the physical range |
| 1.6 | §2.4 | serialised "with full precision (`%.17g`)" | `Auxiliary/writeParamsToMFile.m:33` writes `%g` — 6 s.f. `%.17g` exists only in three purpose-built writers | Correct, and note the round-trip cost was measured at ~0 on the objective |
| 1.7 | §2.6 | surrogate bounds "derived from `Model/parameterBounds.m`" | `optimizeFeatures.m` calls it then immediately does `bounds = struct();`. Search box is fixed defaults (0.4–2.5; 0.3–3 for knots) | `parameterBounds.m` supplies the physiology *penalties*, not the search box |
| 1.8 | §2.5 | `w_k c_k = 1` = one preparation-width | Arithmetically λR·n ≈ 0.10. And F–V is **hard-pinned** at 40 % of the seed total, not λ/s² | Rewrite the χ² sentence; disclose the pin and the 30 % share cap |
| 1.9 | §2.5 | ATP-relevance factor R presented as active | `KAPPA = 0` in `BuildWeightedObjective.m`; R ≡ 1 in every reported fit. The code comment says the zero is deliberate, for the leakage reason | State κ = 0 and give the reason |

↪ Not independently re-checked, from the code reviewer (30 findings total; these are the
ones with practical consequence): the 10⁴ s⁻¹ clamp does not apply to the PCHIP branch that
runs; `RunKtr = RunStairs = 0`, so those protocols are not executed inside the objective at
all; the σ = 0.055 provenance is inverted (it is the larger of two snapshot-guard phantom
gains, not a control spread, and the applied threshold is 2σ = 0.110); `F` denotes two
different forces across §2.3.6/§2.3.8/§2.3.9; the knot abscissae carry an undocumented
fitted `_logOffset` translation (−20.95 nm on φ₋₁); `UseA2Reattaching` carries two distinct
fluxes, one of which bypasses both availability gates; `UseSerialStiffness` is read nowhere
in the ODE; `coolDownLS` is a data fit, not a guardrail, and has since been retired.

**LaTeX:** `\rightleftarrows` needs `amssymb`, not `amsmath`. Everything else in scope
compiles.

---

## 2. Data and protocol errors (§1.6, §2.1, §2.2)

| # | Section | Issue | Status |
|---|---|---|---|
| 2.1 | §2.2.1 | Reproduces an **explicitly retracted** recommendation. `RundownCorrection/strategy.md` §0(b) is headed *"'Correct force, do not correct `ktr`' was wrong"*; §"For REPORTING": *"Do not correct one and not the other."* It is also in the pinned banner at the top of `Analyses/README.md` | ✅ |
| 2.2 | §2.2.3 | Describes an 11 ms Savitzky–Golay smoothing of the re-stretch peak. `SMOOTH_WIN_FAST_S` appears **once** — its own declaration. Line 297: `y_sm = y;  % raw — smoothing not applied to this fast section` | ✅ |
| 2.3 | §2.1.1 | "No optical sarcomere-length signal is present in these recordings" is false for two of three preparations. Raw headers: 03/27 `Time \| L in \| F in`; 04/03 and 04/10 `Time \| L in \| F in \| SL in`. Curation discards the column deliberately | ✅ |
| 2.4 | §2.1.2 | PNB expanded as *N*-benzyl-*p*-toluenesulphonamide. That is **BTS**, a fast-skeletal myosin inhibitor essentially inactive on cardiac myosin. The repo only ever writes "PNB" (50 µM); §1.5 of this same draft says blebbistatin | ✅ delete the expansion |
| 2.5 | §1.6 | "reproduces … all four protocols, including two that were never scored" — `MechanismAblation/REPORT.md` says of `lean_opt`: *"kʈr (never scored) \| **worse**"*; `StairsPrediction` §3 reports the staircase elastic peak 2.0–4.1× too large with an undershoot of the wrong sign | ↪ |
| 2.6 | §1.6 | "13 of 20 features inside 1 SD" is **slack-only** and computed for `pooled8_opt`, two lineage steps behind current | ↪ |
| 2.7 | §1.6 | "pooled across four protocol days" is false for force — Baker is excluded from the force number, so force is n = 3; only `ktr` is n = 4. §2.2.2 already says n = 3 | ↪ |
| 2.8 | §1.6 | "35 ablation arms" — the report header says 30 of 34 completed; and one *failed* arm was tracking 0.85, better than the 0.9081 incumbent | ↪ |
| 2.9 | §2.2.1 | "CV 20 % → 7 % at no cost in assumptions" — the assumption-free own-peak step buys 20 % → 17 %; the 7 % requires φ and the damage model | ↪ |
| 2.10 | §2.2.3 | "segment-averaged quantities use cycles 1–4" — only `rsK80[1:4]` does; every other slack feature is scored over all five | ↪ |
| 2.11 | §2.1.2 | FV2 is a four-phase load-clamp, not a ramp, and exists only on the 04-series — "each run contained" is false for a third of the dataset. The within-run order given is also wrong | ↪ |
| 2.12 | §2.1.1 | Species: the repo says mouse in some places and rat in others, *including within one file*. It is load-bearing — "rat" set `P_bound_max` ≈ 0.12–0.15, "mouse" set the `ktr` ≈ 50 s⁻¹ target | ✅ contradiction confirmed |

**⚠ Reviewer over-claim, corrected.** The methodology review stated "every current snapshot
carries `tau_reg = 0.0002`", making the ktr protocol a selection set. ✅ `floored_seed.m`
carries **`tau_reg = 0.0266779`**. The fast-gate value belongs to the `fastgate`/`ktrfence`
lineage. The leakage charge **stands for `estiff`** (`floored_seed` carries 1.24383, and
`StairsPrediction` ranks arms on the staircase before saying "adopt this: `estiff125_opt`")
but **not** for `tau_reg` in the configuration the draft documents.

Consequence for §2.5: "held-out prediction test" is wrong for the staircase (used to select
`estiff`). For the ktr protocol it is right in this lineage, but `ktrfence_seed` has
`RunKtr = 1` and a scored `ktrP_k63` fence, so the scalar is fenced while the trace stays
held out — state the distinction.

---

## 3. Methodological concerns (not errors — things a reviewer will attack)

Ranked. All ↪ from the methodology review unless marked.

1. **Feature non-independence.** The objective sums per-feature χ² as if independent, but
   `steady`, `A`, `ktr`, `peak1_y`, `vall_y`, `rsK80` come off the same traces, share one
   noisy per-preparation anchor, share one rundown factor, and `rsK80` takes its amplitude
   *from* `steady` by definition. Effective d.o.f. is nearer 5–8 than 20. `DataPooling` §1
   has already fitted the covariance model needed (PC1 = 88 %, flat loading).
2. **The deadband is not a fit statistic.** A clamped term has exactly zero gradient inside
   its band; 10 of 16 data terms currently pay nothing, including `steady`. Reporting
   "2.4088 → 0.9081" as an improvement is not defensible — nothing about the model changed.
   Report the *unfloored* per-feature χ² as the goodness-of-fit table.
3. **The level has no gradient.** `steady` sits at 0.34 preparation-widths, paying nothing,
   so the model's absolute force may drift by ~43 % of the entire ATP effect size at zero
   cost — in a paper whose effect is largely a level effect.
4. **No uncertainty anywhere.** No CIs on the ATP ratios, on F_ref, or on any parameter, for
   a paper making mechanism claims against an 11-dimensional identifiability ceiling.
   Profile likelihoods on the six parameters the claims rest on would be the single most
   persuasive addition.
5. **`estiff` sits on a knife-edge.** A 2 % move alone costs Et 0.979 → 58.75; the
   compensating `kstiff2` tube is ~1 % wide; along the ridge the objective is flat over
   Δe = ±0.05. Locally unidentifiable along its own strongest direction — disclose.
6. **The objective's "noise" is a wall-clock ODE event**, not solver non-determinism, and it
   biases *downward* while the search rule is accept-if-better. `MechanismAblation` §10
   measured the objective as bit-reproducible. One-line fix; then delete §2.7's explanation.
7. **Ablation power.** 2σ = 0.110 on a ~0.9 objective means anything worth < 12 % is
   undetectable by construction, and the design is blind to substitutable mechanisms. No
   AIC/BIC, no cross-validation, no nested test.
8. **Passive circularity.** Of the "three independent lines" for ×1.2, one is genuinely
   independent; one is the 2 mM cost that the passive branch feeds; one is an
   absence-of-exclusion. The passive terms are 27 % of the objective and unfloored.

---

## 4. Novelty and positioning

**Verdict: solid, well-executed, incremental dressed as novel.** Realistic home *Biophysical
Journal*; *JGP* only if D1 is settled and D2 is done; *PLoS Comp Biol* if the identification
methodology is admitted as the headline.

**Drop as claimed contributions:** "six-pool", and "force-dependent SRX mobilisation supplies
thick-filament mechanosensing". Beard 2022 already had recruitment; Campbell 2018 and
Lewalle 2024 own that mechanism; and §2.3.6 concedes two of six legs refine to zero, so a
reader meets "six-pool" and finds a four-pool directed cycle.

**Lead with instead:** (i) the SRX reservoir resolved **by nucleotide occupancy**, with the
ATP gate on the S_D→D_D return — the one place ATP can set head *availability* independently
of cycling rate; (ii) the reference-normalised gate construction, which leaves the 8 mM fit
bit-identical for any K_m so 2 mM is a true extrapolation; (iii) the two steady-state-inert
availability states; (iv) the identification methodology and the rundown correction, which
other permeabilised-preparation labs would use and which is currently invisible in §1.

**Missing literature (~150 words to fix):**

| Topic | Reference |
|---|---|
| The actual first metabolite-explicit cardiac XB model (predates Tewari 2016) | Tran K, Smith NP, Loiselle DS, Crampin EJ. *Biophys J* 2010;98(2):267–276. PMID 20338848 |
| The live competitor — human atrial trabeculae across ATP *and* Pi | Musgrave JH et al. *J Physiol* 2026, doi:10.1113/JP288463 |
| Force-dependent OFF state formally analysed | Lewalle/Campbell/Niederer, *Biophys J* 2024;123(18), doi:10.1016/j.bpj.2024.05.025; and PMC11427768 |
| Strain-distributed PDE + thick-filament activation, history dependence | Liu, … Walcott. *Biophys J* 2024;123(5):555–571. PMID 38291752 |
| Spatially explicit alternative being declined | FiberSim (*Biophys J* 2021, PMID 34932957); MUSICO (*JGP* 2016) |
| Myosin modulators — the clinical frame, and we *use* mavacamten without citing it | Solaro, *PNAS* 2025;122(41):e2521205122 |
| CP/diffusion (see D2) | Ogut & Brozovich, *Circ Res* 2003;93:54–60 |
| Pi force authority | Ebus & Stienen, *Pflügers Arch* 2001;443:102–111 |
| OFF/ON ≠ SRX/DRX | Ma W et al. *PNAS* 2023;120(5):e2207615120, PMID 36696446 |
| Already in our `.bib`, **uncited** | `Schmidt2024` |

**Structural:** cut H4 — it is a claim about experimental design, not about the muscle, and
it repeats §1.3. Make H1–H3 competing and falsifiable by naming, for each, the observable
that discriminates it. Rewrite §1.1 around the K_m paradox (cardiac actomyosin K_m for
MgATP is 10–150 µM; HF leaves ATP millimolar; so ATP depletion *should* be mechanically
irrelevant — yet halving it costs power) rather than the generic opening.

---

## 5. Citation audit — completed 2026-09-21 for the load-bearing entries

### 5a. Two bibliography entries were wrong; both fixed in `references.bib`

**`Rayner2024meta` → `Bray2026`.** I had invented the first author. The paper is
Bray JJH, Henry JA, Singhal A, Mureyi T, Couch LS, Ahmad M, Valkovič L, Neubauer S,
Lewis AJM, Rider OJ, *Eur Heart J Imaging Methods Pract* 2026;4(3):qyag110, PMID 42422265.
All quoted numbers check out: 176 studies, 7843 examinations 1987–2025; PCr/ATP 1.96
healthy (n = 2988), 1.54 HFrEF (n = 575, −21 %), 1.60 HFpEF (n = 146, −18 %).

**`Nag2026SRX` → `Ma2026SRX`, and this one is substantive: I cited an editorial for a
primary data claim.** The SRX-reduced-in-HFrEF finding belongs to
Ma W, Chaami C, Melhedegaard EG, Lewis CTA, Zhang C, Guo M, Flannagan H, Landim-Vieira M,
Pyle WG, Campbell KS, Ochala J, *JACC Basic Transl Sci* 2026;11(9):101651,
doi:10.1016/j.jacbts.2026.101651, PMID 42612577. What I cited (PMID 42612576) is the
accompanying **editorial**, single-authored by **Anthony J. Baker** — a co-author of this
manuscript. Both are now in the `.bib`; data claims must cite `Ma2026SRX`, and
`Baker2026ed` only as commentary, with the co-authorship visible.

### 5b. Verified correct — no change needed

`Beard2022` · `Cooke1985` · `Rhodehamel2025` · `Daniel1998` · `Tanner2007` ·
`Yamashita1994` · `Bazil2016` · `Jezek2025` · `Collins2026` · `Ogut2003`.

Claim-level checks that passed: Cooke & Pate support the draft's exact statements (Pi lowers
force without changing $V_\mathrm{max}$; force falls linearly in $\log[\mathrm{P_i}]$;
raising [MgADP]/[MgATP] lowers $V_\mathrm{max}$ and raises isometric tension). Daniel 1998
and Tanner 2007 both support the compliant-realignment rationale — Tanner explicitly
describes "a compliant network of protein filaments where there is mechanical coupling
between myosins along the thick-filament backbone and between actin monomers along the thin
filament". `Bazil2016` is bibliographically correct; whether it is the right citation for
"a validated mitochondrial model" is an editorial choice — it is an oxidative-phosphorylation
/ ROS model, and `Collins2026` may carry that sentence alone.

### 5c. Rhodehamel 2025 — REVISED 2026-09-21 after reading the Methods

**The earlier entry below over-called this as a contradiction. It is not.** Reading the
Methods rather than the abstract:

- **ATP regeneration is present in both preparations.** Theirs: 15 mM sodium phosphocreatine
  + 15 U/mL creatine phosphokinase. Beard 2022: 12 mM CP + 10 IU/mL CK. So ADP is held low
  in both, and two independent laboratories see an ATP effect where ADP cannot be the
  carrier. **This supports the H1 dismissal rather than undermining it** and should be used
  that way in §1.4. Neither study added ADP or used a Pi scavenger.
- **Their "no change" applies only at saturating Ca.** T_max at pCa ≈ 4.3, SL 2.1 µm:
  22.30 ± 1.42 / 23.25 ± 1.14 / 22.70 ± 1.85 mN/mm². But they *also* report a **rightward
  shift of the tension–pCa curve with rising ATP** and "depressed isometric force production
  at intermediate calcium levels" — i.e. **lower ATP gives more force at submaximal Ca**,
  our direction. The conflict is confined to T_max.
- **They disclose the confound on exactly that measurement:** *"While our mechanical assays
  were conducted with matched magnesium concentration through the use of MgATP, the ionic
  strength was not matched."* Ionic strength is a strong modulator of maximal force in
  skinned myocardium.
- **Their kinetics agree with us and more steeply.** Detachment 21.91 / 38.54 / 57.07 s⁻¹ at
  2 / 6 / 8 mM — low ATP at **×0.38**, a stronger slowing than either of our `ktr` probes.
- **Real differences to name, not smooth over:** porcine LV (β-MHC) vs mouse RV (α-MHC),
  ~4× apart in detachment kinetics; their T_max ~22 kPa vs our pooled ~69 kPa; and their
  $V_\mathrm{max}$ falls ×0.49 at 2 mM where \citet{Beard2022} found it unchanged — that one
  is a genuine disagreement.
- ⚠ **Internal inconsistency in the source.** Methods list 2.25 / 4.75 / 7.25 mM Na₂ATP;
  Results report "2, 6 and 8 mM". Likely nominal-added vs computed free MgATP, or two
  experiment sets. **Verify before quoting a concentration.**
- ✗ My earlier "mitigating difference: 2.25–7.25 mM not 2–8 mM" was wrong — the ranges are
  essentially the same.

**How to use it in the paper:** cite it for the structural OFF→ON result and for the
detachment-rate agreement; state the T_max difference once, with the ionic-strength caveat
and the isoform difference, rather than omitting it. Our force result is backed by
cross-bridge logic and by their own submaximal-Ca data.

---

#### Superseded version of 5c (kept for the record)

### ~~[OPEN] Rhodehamel 2025 is cited selectively, and the omission cuts against us~~

Verbatim from its abstract: *"high ATP did not alter maximum calcium-activated tension,
although increasing ATP right-shifted the tension versus calcium (Ca²⁺) curve and
accelerated both myosin attachment and detachment rates… Power output and maximum unloaded
shortening velocity also significantly increased with increased ATP concentration."*

Inverted to low ATP, that paper reports: **tension unchanged**, rates slower, power lower,
$V_\mathrm{max}$ **lower**. Against our result (force ×1.34, `ktr` down, $V_\mathrm{max}$
unchanged) it **agrees on the rate slowing and disagrees on both force and $V_\mathrm{max}$**.
§1.2 currently uses it as structural support for H3 while omitting the half that contradicts
our headline force effect. Mitigating differences — porcine not mouse, 2.25–7.25 mM not
2–8 mM, ~22 °C mechanics — are real and should be stated, but the disagreement has to appear
in the text rather than be left for a referee to find.

### 5d. Not checked at claim level (bibliographic fields only)

`Gordon2000`, `Hooijman2011`, `McNamara2015`, `Linari2015`, `Irving2017`, `Campbell2018`,
`Huxley1957`, `Razumova1999`, `Rice2008`, `Tewari2016`, `Land2017`, `deTombe1992`,
`deTombe2016`, `Smyth2004`, `Neubauer1997`, `Lopez2020`. These are standard attributions
where the risk is low, but the "does it support *this* sentence" test was not run on them.
`Schmidt2024` remains in the `.bib` and **uncited** — use it or drop it.

Previously resolved: Ježek et al. is *J Physiol* **2025**, doi:10.1113/JP288807,
PMID 40846493 — the poster and abstract say 2026.

## 6. Superseded — original "Incomplete" section

The citation audit (claim-level verification of every `\citep` against its source, plus the
`.bib` field check) **failed partway** on an account spend limit. It had reached Cooke &
Pate 1985 and the Daniel 1998 / Tanner 2007 pair. Still unverified:

- Whether each citation supports the *specific* claim it is attached to.
- Whether Rhodehamel 2025's omitted findings — maximum Ca-activated tension **unchanged**
  with ATP, and a right-shifted tension–pCa curve — make our use of it misleading, given we
  argue elsewhere that low ATP *raises* force. **This one matters and should be resolved.**
- `.bib` field accuracy for ~25 entries; the three already flagged `UNVERIFIED`/`CHECK`
  (`Bray2026` author list and its 2024-vs-2026 year, `Ma2026SRX` — and which of the
  two closely related 2026 papers is the right one — and the `Bazil2016` placeholder).
- Previously known and unresolved: Ježek et al. is *J Physiol* **2025**, doi:10.1113/JP288807,
  PMID 40846493 — the poster and abstract say 2026.
