# A2 re-attachment during restretch — analysis

**Question.** Does a restretch-gated *direct re-attachment* pathway (`UseA2Reattaching`)
improve the 2-state restretch fit, and is it physiologically plausible? Base =
`params/opt2state_v2_opt.m` (2-state, cost ~3.22). Restretch is the dominant residual
region (`peak1_dSL`, `peak1_y`, `vall_y`, `ovrsht`).

## The mechanism

During lengthening (`vel>0`), ATP-cocked heads from the unattached permissive pool `PT`
bind **directly into the strong post-stroke state `p2`** at strain `dA2re` (default 0), and
are then carried to positive strain by the ongoing stretch — so they add force through the
restretch/valley phase. Rate `RT2 = kA2re·PT`. `PT` is the complementary pool, so no
explicit outflow term is needed. Conceptually: a head **"jumps to the next available actin
site"** presented by the sliding filament during stretch.

Implementation notes (this campaign):
- The coded `UseA2Reattaching` block was **non-functional**: it wrote `dp2` and used `s_i0k`
  *before either was defined*, so switching it on only threw. Fixed: flux applied after `dp2`
  exists, correct interpolation weight (`s_i0k_a2`), and **gated to `vel>0`** so isometric /
  FV / shortening are byte-identical when off or contracting.
- Two dedicated params added (`getParams`): `kA2re` (rate, default 0 = inert) and `dA2re`
  (landing strain, default 0). The original reused `k_2`, which also sets the p2→p1 reverse
  rate — a confound now removed.

## Plausibility

**Moderate — a reduced description of a real phenomenon (stretch-induced re-attachment).**

For:
- **Eccentric / stretch recruitment.** Actively lengthening muscle bears more force than
  isometric (residual force enhancement); part is attributed to rapid cross-bridge
  re-attachment to newly presented actin sites (Edman; Herzog eccentric literature). A
  stretch-gated attachment pathway captures this directly.
- **Actin-site periodicity ("next spot").** Actin monomers repeat ~5.5 nm; during lengthening
  the filament slide continuously presents fresh sites, so a cocked head can immediately rebind
  to the next site (Squire; the same geometry motivating A2 hopping, mechanism-eval §1.9).
- **Load-favoured strong binding under stretch** (catch-bond-like strengthening; Guo &
  Bhattacharya) makes rapid strong attachment during lengthening reasonable.

Against / caveats:
- **Direct `PT→p2` bypasses the weak-bound `p1`.** Canonical attachment is weak→strong; binding
  straight to strong is a kinetic shortcut (defensible if weak→strong is ~instantaneous under
  stretch, but a simplification).
- **`kA2re` is phenomenological** — no direct measurement; `vel>0` is a hard switch, not a smooth
  strain/velocity dependence.

## Hypothesis (why it may succeed where catch bond failed)

Catch bond raised `vall_y` but dragged `peak1_y`/`peak1_dSL`/`peak2` up with it, because it
*freezes the same strained bridges* that make the spike — valley and peak are coupled. A2
re-attachment instead adds **new low-strain heads** during restretch, which bear force through
the valley/recovery *without* contributing to the initial elastic spike (that spike comes from
pre-existing high-strain bridges). So the prediction is: **`vall_y` rises toward data (~71)
with `peak1_y` held (~96)**. Being `vel>0`-gated, `steady`/FV/shortening cannot be perturbed —
so if it works, it is a "free" restretch-only lever.

## Results

`probeA2re` off vs kA2re/dA2re sweep (base cost 3.224):

| lever | vall_y | peak1y | peak2 | steady | ktr | FV | cost |
|---|---|---|---|---|---|---|---|
| BASE off | 66.9 | 96.1 | 79.0 | 77.4 | 52.4 | .93 .55 .22 .09 | 3.224 |
| **kA2re10 d0** | **70.2** | **97.5** | 81.8 | 77.4 | 51.8 | .93 .55 .22 .09 | **3.124** |
| kA2re30 d0 | 76.4 | 100.6 | 87.2 | 77.4 | 51.1 | (same) | 4.662 |
| kA2re80 d0 | 90.2 | 108.9 | 99.3 | 77.4 | 49.9 | (same) | 27.00 |
| kA2re30 d.010 | 73.2 | 98.9 | 84.4 | 77.4 | 52.8 | (same) | 4.093 |
| kA2re80 d.010 | 82.9 | 104.1 | 92.9 | 77.4 | 52.5 | (same) | 15.22 |

_(data: vall_y ~71, peak1_y ~96, peak2 ~81)_

**Hypothesis CONFIRMED.** At `kA2re≈10` the valley lifts 66.9→70.2 (≈ data 71) **with `peak1_y`
held (96.1→97.5)** and `peak2` improved (79→82 ≈ data), while `steady`/FV/ktr are unchanged
(the `vel>0` gate guarantees it) — and total cost drops 3.224→**3.124**. This is precisely the
behaviour catch bond could NOT produce (catch bond raised valley and peak together); adding
*new low-strain heads* rather than freezing existing strained ones is the mechanistic difference.

Higher rates overshoot: `kA2re=30` over-fills the valley (76.4) and inflates `peak2`; `kA2re=80`
runs away (cost 27). So there is a sharp optimum near `kA2re≈10`. The landing strain `dA2re` is a
minor second knob (raising it spreads the added force to higher strain, mildly softening the
over-fill at high `kA2re`).

## Verdict & recommendation

A **real, physically-motivated, orthogonal restretch lever** — but a **minor** one. It fixes
`vall_y` (and helps `peak2`) essentially for free (no steady/FV/ktr cost), confirming the
"new heads, not frozen heads" reasoning. It does **not** address the *dominant* restretch
residual `peak1_dSL` (0.74, peak timing/excursion) or the FV mid-tail.

## Extension: peak-relief HOP (kA2hop) — does re-attachment cap the peak?

Hypothesis (user): the "jump to next site" should also *detach* the strained peak heads
(faster detachment ramps with strain), capping the spike. Added a restretch-gated,
strain-dependent hop `R_hop = kA2hop·max(0,s−sA2hop)·p2`: high-strain p2 heads leave and
re-land at low strain (mass-conserving). `probeA2hop`:

| lever | peak1_y | peak1_dSL | vall_y | peak2 | cost |
|---|---|---|---|---|---|
| BASE off | 96.1 | 0.0346 | 66.9 | 79.0 | 3.224 |
| **hop30k** | 95.5 | 0.0345 | 70.0 | 81.0 | **2.952** |
| hop80k | 95.1 | 0.0472 | 74.3 | 84.0 | 8.21 |
| re10+hop8k | 97.4 | 0.0351 | 71.0 | 82.4 | 3.10 |

**Two findings:**
1. **The hop does NOT cap the peak.** `peak1_y` (96.1→95.5) and `peak1_dSL` (0.0346→0.0345) are
   ~inert at useful rates and *worsen* at extreme rates. `peak1_y` is the **elastic recoil of
   pre-existing strained bridges** (forms in ~2 ms, faster than any hop can relieve), and
   `peak1_dSL` is a **series-compliance/geometry** quantity (slack takeup + kSE) — neither is a
   p2-population term the hop can reach. The user's peak-relief hypothesis is *not* supported.
2. **But the hop is the best valley lever found:** `hop30k` → cost **2.952** (−0.27, better than
   the `kA2re` valley-fill's 3.12), by shifting heads from the peak region into the valley
   (`vall_y`→70, `peak2`→81, both to data). It is either/or with `kA2re` (combining over-fills).

Because `peak1_dSL` is series-compliance-set, the correct co-lever is `kSE` (`probeA2retune`):

| lever | peak1_dSL | ktr | vall_y | peak2 | cost |
|---|---|---|---|---|---|
| hop30k | 0.0345 | ~52 | 70.0 | 81.0 | 2.952 |
| hop30k kSE×1.6 | 0.0297 | 59.5 | 70.2 | 82.7 | 5.29 |
| hop30k kSE×2.0 | **0.0281** | 60.7 | 70.2 | 83.3 | 4.62 |

`kSE↑` **does** fix `peak1_dSL` (0.0346→0.0281 ≈ data 0.025) — confirming it is a series-compliance
quantity — but it **re-inflates ktr** (52→61; stiffer series transmits force faster → faster
redevelopment). The ktr penalty (+ small `peak1_y`/`peak2` rises) outweighs the `peak1_dSL` gain,
so net cost rises. `peak1_dSL` is therefore **welded to ktr through `kSE`** — a third 2-state
coupling (alongside ktr↔FV-tail and vall_y↔vall2), and `kmsrd` (the ktr lever) is already pinned
near its bound, so there is little room to compensate.

## Final verdict

- **The A2 hop is a real, keepable win: cost 3.224→2.952** (fixes `vall_y`, `peak2` to data). Its
  lever is the **valley**, not the peak. Recommended: fold `kA2hop` (~2e4–4e4, `sA2hop`~0.005,
  hop-only — `kA2re` is redundant with it) into the 2-state pool. Saved config:
  `params/params_2state_a2hop.m`.
- **The hop cannot touch `peak1_y`/`peak1_dSL`** — elastic-recoil + series-compliance, not
  detachment-limited. `peak1_dSL` is fixable via `kSE` but welded to ktr. This is the cleanest
  specific case for the **3rd state**: its independent slow ktr timescale would free `kSE` to
  stiffen for `peak1_dSL` without the ktr penalty.

Plausibility of the mechanism itself is unchanged (moderate; stretch-recruitment + actin-site
geometry). What this campaign settled is *which residual it addresses* (valley, not peak).

