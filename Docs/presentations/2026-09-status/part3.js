// PART 3: low ATP on the current fit, ATP levers, where we stand, future work, backup.
// Standalone: NODE_PATH=~/tools/pptx/node_modules node part3.js  -> part3.pptx
// Slide 15/15b read Analyses/LowATP_Refit9K2Sweep/results/* when present (placeholders otherwise).
// Override the summary JSON for testing with env K2SWEEP_JSON=/abs/path.json
const fs = require("fs");
const path = require("path");
const S = require("./style");

const SWEEP = "Analyses/LowATP_Refit9K2Sweep/results/";

// Pixel size of a PNG (null if missing/unreadable) -> keeps aspect ratio in S.image
function pngSize(rel) {
  try {
    const b = fs.readFileSync(S.asset(rel));
    if (b.toString("ascii", 1, 4) !== "PNG") return null;
    return [b.readUInt32BE(16), b.readUInt32BE(20)];
  } catch (e) { return null; }
}
function img(s, rel, box) { S.image(s, rel, box, pngSize(rel)); }

// Sweep summary (optional). Expected shape documented in the brief; every key optional.
function loadSweep() {
  const f = process.env.K2SWEEP_JSON || S.asset(SWEEP + "k2sweep_summary.json");
  try { return JSON.parse(fs.readFileSync(f, "utf8")); } catch (e) { return null; }
}
const fmt = (x, d = 2) => (typeof x === "number" && isFinite(x) ? x.toFixed(d) : "–");

function addSlides(pres, S) {
  const { C, F, M, W, TOP } = S;
  const J = loadSweep();

  // ------------------------------------------------------------------ 15
  {
    const s = pres.addSlide();
    S.header(s, "Low ATP on the current fit: only A2 detachment feels ATP");
    const LW = 6.3;
    S.card(s, { x: M, y: TOP, w: LW, h: 3.2 }, "Design: one physiological lever", [
      "8 and 2 mM MgATP; at 2 mM only R2 (A2 detachment) × g",
      "g = 1 (no effect), 0.75, 0.54 (poster), 0.44, 0.35 (earlier best), 0.30",
      "Hydrolysis not gated (bound ATP: physiological choice)",
      "Scored on pooled 2 mM data: absolute χ² on the 2 mM between-prep spread",
      "Term weight min(z,3)/3, z = |r−1|/sd: ATP-insensitive features drop out",
      "Also reported: plain 2 mM χ²; model vs data 2/8 ratios",
    ], { fontSize: 14, bar: C.navy });
    S.card(s, { x: M, y: TOP + 3.3, w: LW, h: 1.85 }, "Earlier result (older core lean5_a_opt, 09-17)", [
      "k2 alone, g ≈ 0.35: force ×1.39 (data ×1.33), ktr ×0.54 (×0.54), rsK80 term → 0",
      "But 2 mM objective 2.13 → 5.27 (vall2_dy, t0, peak1_dSL, vall_t)",
      { text: "Same re-stretch defect tops 8 mM: not a low-ATP problem", bold: true },
    ], { fontSize: 14, bar: C.warn });
    // physics strip
    s.addShape("rect", { x: M, y: TOP + 5.25, w: LW, h: 0.5, fill: { color: "FFF4CC" }, line: { color: C.maize } });
    s.addText([
      { text: "F0 ∝ k2", options: { bold: true } }, { text: "−0.78", options: { superscript: true, bold: true } },
      { text: ",  ktr ∝ k2", options: { bold: true } }, { text: "+1.09", options: { superscript: true, bold: true } },
      { text: ": force and ktr rigidly anti-coupled via k2" },
    ], { x: M + 0.1, y: TOP + 5.25, w: LW - 0.2, h: 0.5, fontFace: F, fontSize: 14, color: C.ink, valign: "middle", margin: 0 });

    // right: cost vs g + result text
    const RX = M + LW + 0.25, RW = W - M - RX;
    img(s, SWEEP + "k2sweep_cost.png", { x: RX, y: TOP, w: RW, h: 3.9 });
    const msgs = J && Array.isArray(J.messages) && J.messages.length ? J.messages.slice(0, 4) : null;
    const head = J && typeof J.best_g === "number" ? "Result: best g = " + fmt(J.best_g) : "Result";
    if (msgs) {
      S.card(s, { x: RX, y: TOP + 4.05, w: RW, h: 1.7 }, head, msgs, { fontSize: 14, bar: C.good });
    } else {
      S.placeholder(s, { x: RX, y: TOP + 4.05, w: RW, h: 1.7 }, "sweep running: result text from k2sweep_summary.json");
    }
    S.footer(s, "Sources: Analyses/LowATP_Refit9K2Sweep (conclusions.md, results/k2sweep_*); floored-campaign: Docs/experiments/status-2026-09-16.md §9; exponents: sensitivity Jacobian (earlier core)");
    s.addNotes("Only one thing changes between 8 and 2 mM: the A2 detachment rate R2 is multiplied by g. That is the step where free ATP binds. " +
      "Hydrolysis is deliberately not gated: the ATP there is already bound. We score the 2 mM model directly against the pooled 2 mM data, " +
      "weighting each feature by how strongly it responds to ATP, so waveform-shape features that do not change with ATP do not decide the answer. " +
      "On the older core, g about 0.35 gave both headline ratios almost exactly, but the 2 mM objective got worse through the re-stretch shape and timing terms, " +
      "the same defect that dominates the 8 mM fit. Physics: through k2 alone, force and ktr move in opposite directions with fixed exponents, so one knob sets both. " +
      (J ? "Sweep result on the incumbent: see the right panel." : "The sweep on the incumbent was still running when this slide was built."));
  }

  // ------------------------------------------------------------------ 15b
  {
    const s = pres.addSlide();
    S.header(s, "Low ATP on the current fit: feature by feature", "2 mM vs pooled data");
    const LW = 6.4;
    img(s, SWEEP + "k2sweep_features.png", { x: M, y: TOP, w: LW, h: 5.7 });
    const RX = M + LW + 0.25, RW = W - M - RX;

    // native ratio chart at best g
    const R = J && J.ratios && typeof J.ratios === "object" ? J.ratios : null;
    let drawn = false;
    if (R) {
      let bi = -1;
      if (Array.isArray(J.g) && typeof J.best_g === "number") {
        let dmin = Infinity;
        J.g.forEach((g, i) => { const d = Math.abs(g - J.best_g); if (d < dmin) { dmin = d; bi = i; } });
      }
      const labels = [], dat = [], mod = [];
      for (const [k, v] of Object.entries(R)) {
        if (!v || typeof v.data !== "number") continue;
        let m = null;
        if (Array.isArray(v.model)) m = bi >= 0 ? v.model[bi] : null;
        else if (typeof v.model === "number") m = v.model;
        if (typeof m !== "number" || !isFinite(m)) continue;
        labels.push(k + (typeof v.data_sd === "number" ? " (±" + fmt(v.data_sd) + ")" : ""));
        dat.push(v.data); mod.push(m);
      }
      if (labels.length) {
        const ones = labels.map(() => 1);
        s.addText("2/8 ratio at best g" + (typeof J.best_g === "number" ? " = " + fmt(J.best_g) : "") + " (label: data ± SD)", {
          x: RX, y: TOP, w: RW, h: 0.4, fontFace: F, fontSize: 13, bold: true, color: C.navy, margin: 0 });
        s.addChart([
          { type: pres.charts.BAR, data: [{ name: "data", labels, values: dat }, { name: "model", labels, values: mod }],
            options: { chartColors: [C.navy, C.maize], barGrouping: "clustered", barGapWidthPct: 60 } },
          { type: pres.charts.LINE, data: [{ name: "no effect (1.0)", labels, values: ones }],
            options: { chartColors: [C.bad], lineSize: 1.5, lineDataSymbol: "none", lineDash: "dash" } },
        ], {
          x: RX, y: TOP + 0.4, w: RW, h: 3.1, barDir: "col",
          catAxisLabelFontSize: 11, catAxisLabelFontFace: F, catAxisLabelRotate: labels.length > 5 ? -30 : 0,
          valAxisLabelFontSize: 11, valAxisLabelFontFace: F, valAxisMinVal: 0,
          showLegend: true, legendPos: "t", legendFontSize: 11, legendFontFace: F,
          valGridLine: { color: "E5E7EB", size: 0.5 },
        });
        drawn = true;
      }
    }
    if (!drawn) S.placeholder(s, { x: RX, y: TOP, w: RW, h: 3.5 }, "native ratio chart (data vs model at best g) from k2sweep_summary.json");
    img(s, SWEEP + "k2sweep_trace.png", { x: RX, y: TOP + 3.65, w: RW, h: 2.05 });
    S.footer(s, "Sources: Analyses/LowATP_Refit9K2Sweep/results (k2sweep_features.png, k2sweep_trace.png, k2sweep_summary.json). Ratios are diagnostic only; the score is the absolute 2 mM χ².");
    s.addNotes("Left: the key 2 mM features against g, with the pooled 2 mM data band; 8 mM model and data as markers. " +
      "Top right: model versus data 2/8 ratios at the best g; the dashed line is no ATP effect. Ratios are diagnostic; the fit is scored on absolute 2 mM values. " +
      "Bottom right: the 2 mM slack and FV2 traces at the best g against the pooled 2 mM data.");
  }

  // ------------------------------------------------------------------ 16
  {
    const s = pres.addSlide();
    S.header(s, "ATP-sensitive levers: tested, and still to test");
    const B = C.bad, Wn = C.warn;
    const rows = [
      ["Lever tested", "Verdict"],
      ["H1: hydrolysis gates kah / kamh (+ k2)", { text: "99.5 % of span on spreadw, but hydrolysis of bound ATP is not physiological", color: Wn }],
      ["H2: SRX gates (4)", { text: "Nothing gained; kmsrd's implied ATP Km implausible", color: B }],
      ["H3/H4: mechanical gates", { text: "Bend compliance; no nucleotide chemistry", color: B }],
      ["k_1 / Pi rebinding (P1)", { text: "Pi predicts the wrong force sign (↓ vs data ↑)", color: B }],
      ["Strain-dependent ATP gate", { text: "No leverage", color: B }],
      ["R2 ceiling + D0 apo pool", { text: "Breaks the ktr overshoot", color: B }],
      ["Parked state", { text: "Wrong direction", color: B }],
      ["2-state ADP trap", { text: "ktr ratio 0.77 vs data 0.54", color: B }],
      ["kstiff2 gate", { text: "Moves both, but ~3× too strong on force", color: B }],
      ["XTOR relaxation", { text: "Attached 0.255 vs physiological 0.12–0.15: force bought with head count", color: B }],
      ["Level-vs-ratio conflict (15 configs)", { text: "Every lever that fixes the 8 mM rsK level worsens the ATP ratio", color: B }],
    ];
    S.table(s, rows, { x: M, y: TOP, w: 7.9, h: 5.6 }, { fontSize: 11, colW: [2.9, 5.0] });
    S.card(s, { x: M + 8.15, y: TOP, w: W - 2 * M - 8.15, h: 5.6 }, "Still to test", [
      "k2 gate on the yield/slip core, after the re-stretch defect is fixed",
      "A forward-limb partner to decouple F0 and ktr",
      "ADP/Pi gates once Pi enters the RHS (Pi does nothing in dPUdT_CombinedTransitions today)",
      "Strain-dependent ATP binding site; SRX gates revisited on pooled data",
      "ktr probe: slack ktr (×0.53) vs ktr protocol (×0.74), reviewer point",
      "Joint 8 + 2 mM fit",
    ], { fontSize: 14, bar: C.maize });
    S.footer(s, "Sources: ATP-wt-lowatp/Analyses/AtpEffectHypotheses/conclusions.md; floored-campaign: Docs/experiments/status-2026-09-16.md §6, §9; project plan (levers list). Most levers tested on earlier cores.");
    s.addNotes("Left: everything we have tried as an ATP-sensitive lever, with the verdict. The only gate that fit well, hydrolysis, is not physiological. " +
      "SRX and mechanical gates either did nothing or bent compliance. Pi goes the wrong way on force. The 15-config test shows a conflict: whatever fixes the 8 mM rsK level makes the ATP ratio worse. " +
      "Caveat: most of these ran on earlier cores. Right: what is still untested. The top item needs the re-stretch defect fixed first. The ktr-probe choice matters because the slack ktr is compliance-contaminated.");
  }

  // ------------------------------------------------------------------ 17
  {
    const s = pres.addSlide();
    S.header(s, "Where we stand (honest summary)");
    const cw = (W - 2 * M - 0.5) / 3, ch = 3.9;
    S.card(s, { x: M, y: TOP, w: cw, h: ch }, "Solid", [
      "8 mM force–velocity and force level",
      "Clean numerics: objective noise 0.17 % (was 19.5 %)",
      "Measured ATP effect settled on 3 pooled preps: force ×1.33, ktr ×0.53, rsK80 ×0.38",
    ], { fontSize: 16, bar: C.good, headColor: C.good });
    S.card(s, { x: M + cw + 0.25, y: TOP, w: cw, h: ch }, "Partial", [
      "Re-stretch timing and segment ladders (ktr, rsK80)",
      "Low ATP on the current core: " + (J && typeof J.best_g === "number" ? "R2 ×" + fmt(J.best_g) + " gets force & rsK80; ktr under-slowed, re-stretch timing wrong" : "see slide 15"),
      "Mechanism identifiability: 79 pool parameters, 32 of them PCHIP knots/offsets",
    ], { fontSize: 16, bar: C.warn, headColor: C.warn });
    S.card(s, { x: M + 2 * (cw + 0.25), y: TOP, w: cw, h: ch }, "Not started / blocked", [
      "Joint 8 + 2 mM fit",
      "Ca²⁺ activation (needs force–pCa data)",
      "Coupling to mitochondrial ATP production",
    ], { fontSize: 16, bar: C.bad, headColor: C.bad });
    S.card(s, { x: M, y: TOP + ch + 0.25, w: W - 2 * M, h: 1.2 }, "Housekeeping", [
      "floored-campaign (2 mM target, ATP gates) merged — staged, not yet committed · manuscript draft under revision (Results not written) · ROADMAP stale (06-16)",
    ], { fontSize: 14, fill: "F3F4F6", headSize: 15 });
    S.footer(s, "Sources: project plan (current state of the science); Analyses/DataPooling; Analyses/LowATP_Refit9K2Sweep.");
    s.addNotes("Solid: the 8 mM force-velocity and force level, the numerics, and the measured ATP effect. " +
      "Partial: the re-stretch transient is the remaining structural defect at 8 mM and it also blocks the 2 mM fit. " +
      "The parameter count is large relative to what the features identify. Not started: a joint 8 and 2 mM fit, calcium, and mitochondrial coupling. " +
      "Housekeeping: the low-ATP branch is being merged, the manuscript Results are not written, the roadmap is stale.");
  }

  // ------------------------------------------------------------------ 18
  {
    const s = pres.addSlide();
    S.header(s, "Future work");
    const LW = 6.1, RX = M + LW + 0.25, RW = W - M - RX;
    S.card(s, { x: M, y: TOP, w: LW, h: 1.5 }, "1  Re-stretch defect", [
      "Yield-gate / A1-slip re-fit campaign (yield objective: 0.558 → 0.501 before the slip refit)",
      "Fixes the block that outvotes ATP physiology at 2 mM",
    ], { fontSize: 14, bar: C.navy });
    S.card(s, { x: M, y: TOP + 1.65, w: LW, h: 1.5 }, "2  Low ATP", [
      "Joint 8 + 2 mM fit with physiological gates only (k2 / R2)",
      "Then test a forward-limb partner to decouple F0 and ktr",
    ], { fontSize: 14, bar: C.navy });
    S.card(s, { x: M, y: TOP + 3.3, w: LW, h: 2.4 }, "4  Mitochondria and the heart wall", [
      "Couple ATP supply (mitochondrial model) to the XB core",
      "TriSeg wall: mechanical interface already matches (SE + viscous element)",
      "Cost 119 s per 0.6 s beat, mostly recoverable; blocker is scientific: no Ca activation",
    ], { fontSize: 14, bar: C.navy });
    S.card(s, { x: RX, y: TOP, w: RW, h: 2.35 }, "3  Ca²⁺ and twitch (worktree ATP-wt-ca)", [
      "Cooperative thin-filament gate on attachment; XB core is fast: t50 7.9 ms, 10 ms pulse → 56 % of F_ss",
      "Literature (Marzban) constants: twitch only 5–6 % of F_ss (1 Hz: 3.47 kPa, TTP 70 ms)",
      { text: "Blocker: no active sub-maximal force–pCa data → ask for force–pCa at 8 and 2 mM", bold: true },
    ], { fontSize: 14, bar: C.maize });
    img(s, "Analyses/CalciumTwitch/results/twitch_probe.png", { x: RX, y: TOP + 2.45, w: RW, h: 3.25 });
    S.footer(s, "Sources: Analyses/CalciumTwitch (conclusions.md, twitch_probe.png); ATP-wt-ca/.../Paced1Hz_report.md; Analyses/TriSegCoupling/conclusions.md");
    s.addNotes("One: fix the re-stretch defect first; the yield-gate and slip variants already reduce the cost without a refit. " +
      "Two: a joint 8 and 2 mM fit with only physiological gates. Three: calcium. The cross-bridge core is fast enough to follow any twitch; " +
      "but with literature thin-filament constants a physiological transient barely activates it, about 5 percent of steady force. " +
      "Amplitude is one steep knob, k_on over k_off, that no data we have constrains. The cheapest fix is force-pCa at 8 and 2 mM on the same rig. Tony, is that feasible? " +
      "Four: mitochondrial coupling and the TriSeg wall; the mechanical interface already matches, the blocker is the missing calcium activation.");
  }

  // ------------------------------------------------------------------ Backup A
  {
    const s = pres.addSlide();
    S.header(s, "Backup: ATP effect per preparation", "2 mM / 8 mM");
    const rows = [
      ["Feature", "03/27 (8→2)", "04/03 (8→2)", "04/10 (2→8)", "mean ± SD", "Pooled target"],
      ["steady, raw", "1.180", "1.296", "1.651", "1.38 ± 0.25", ""],
      [{ text: "steady, rundown-corrected", bold: true }, "1.225", "1.403", "1.412", { text: "1.35 ± 0.11", bold: true }, { text: "1.335", bold: true }],
      ["ktr (slack)", "0.537", "0.593", "0.492", "0.54 ± 0.05", "0.53"],
      ["ktr80", "0.537", "0.657", "0.542", "0.58 ± 0.07", "0.57"],
      ["rsK80", "0.369", "0.406", "0.365", "0.38 ± 0.02", "0.38"],
      ["rsK", "0.363", "0.395", "0.277", "0.35 ± 0.06", ""],
      ["t0", "1.306", "1.353", "1.657", "1.44 ± 0.19", "1.42"],
      ["rsT20", "1.479", "1.532", "2.222", "1.75 ± 0.41", "1.72"],
      ["ktr protocol", "", "", "", "×0.72, CV 4 %", "0.74"],
    ];
    S.table(s, rows, { x: M, y: TOP, w: W - 2 * M, h: 4.3 }, { fontSize: 13, colW: [3.3, 1.65, 1.65, 1.65, 2.0, 2.083] });
    S.bullets(s, [
      "Rates are not rundown-corrected (correcting ktr degrades consistency); force needs it: 04/10 ran 2 → 8 mM",
      "04/10 stays highest on force even corrected; the two 8→2 preps alone give ×1.27, so the mean may be an over-estimate",
      "Per-prep values: simple segment means from pool_raw.mat; pooled target uses the full pooling method, hence small differences",
    ], { x: M, y: TOP + 4.45, w: W - 2 * M, h: 1.3 }, { fontSize: 14, space: 3 });
    S.footer(s, "Sources: Analyses/ATPEffectReconciliation/conclusions.md (force per prep, ktr protocol); Analyses/DataPooling/results/pool_raw.mat (rate ratios, recomputed); pooled values: current pooled target");
    s.addNotes("Per-preparation 2 over 8 mM ratios. Force needs the rundown correction because 04/10 ran in reverse order; rates are used raw. " +
      "rsK80 is the tightest ATP effect across preps. The ktr-protocol value is shown only pooled; per-prep values not tabulated here.");
  }

  // ------------------------------------------------------------------ Backup B
  {
    const s = pres.addSlide();
    S.header(s, "Backup: mechanisms available but OFF", "in refit9");
    const L = [
      ["Theme / switch", "Meaning"],
      [{ text: "ATP gates", bold: true, fill: "FFF4CC" }, { text: "", fill: "FFF4CC" }],
      ["UseAtpDetachR2", "Michaelis ATP gate on R2 (K_T2), = 1 at 8 mM"],
      ["UseAtpKmsrd", "ATP gate on kmsrd (SRX-ADP return)"],
      ["UseR2Ceiling", "ATP-binding ceiling on R2 at high strain"],
      ["UseD0State", "nucleotide-free pool for ripped heads"],
      [{ text: "Calcium", bold: true, fill: "FFF4CC" }, { text: "", fill: "FFF4CC" }],
      ["UseCa", "thin-filament activation (legacy NP slot)"],
      [{ text: "A2 strong-state variants", bold: true, fill: "FFF4CC" }, { text: "", fill: "FFF4CC" }],
      ["UseA2AttachmentShift", "strained A2 hops to next actin site"],
      ["UseA2Reattaching", "PT binds directly into A2 on re-stretch"],
      ["UseA2Popping", "fast drain of A2 beyond a strain limit"],
      ["UseA2MechanicalRecocking", "strained A2 ratchets back to A1"],
      ["UseNegativeForceRip", "force-dependent R1D rip term"],
    ];
    const R = [
      ["Theme / switch", "Meaning"],
      [{ text: "Stretch response", bold: true, fill: "FFF4CC" }, { text: "", fill: "FFF4CC" }],
      ["UseCatchBond", "less R1D/R2 detachment in rapid stretch"],
      ["UseStretchActivation", "boost ka during rapid re-stretch"],
      ["Yield gate / A1 slip", "forced detachment / slip of overstretched heads (in progress)"],
      ["Lengthening recruitment", "lengthening counts in A_reg (gain = 0 now)"],
      [{ text: "Elastic / passive", bold: true, fill: "FFF4CC" }, { text: "", fill: "FFF4CC" }],
      ["UseViscoelasticSE", "dashpot parallel to SE during re-stretch"],
      ["UseDynamicPassive", "velocity-dependent titin force"],
      ["UseTitinIdentifiedPassive", "titin law identified from passive data"],
      [{ text: "Attachment / SRX variants", bold: true, fill: "FFF4CC" }, { text: "", fill: "FFF4CC" }],
      ["UseTargetZoneSaturation", "per-bin occupancy cap (deprecated)"],
      ["UseVernierVelocity / UseVelGaussAttachment", "velocity-dependent ka"],
      ["UseDirectSRXTransition", "all R2 detachment flux → SRX"],
    ];
    const tw = (W - 2 * M - 0.3) / 2;
    S.table(s, L, { x: M, y: TOP, w: tw, h: 5.5 }, { fontSize: 11, colW: [2.55, tw - 2.55] });
    S.table(s, R, { x: M + tw + 0.3, y: TOP, w: tw, h: 5.5 }, { fontSize: 11, colW: [2.55, tw - 2.55] });
    S.footer(s, "Source: Model/getParams.m defaults and comments; Model/dPUdT_CombinedTransitions.m. All ATP gates are normalised to 1 at 8 mM, so the 8 mM fit is invariant to them.");
    s.addNotes("For reference: everything implemented but switched off in the incumbent, grouped by theme. The ATP gates and calcium are all off. " +
      "The yield gate and A1 slip are the in-progress fix for the re-stretch defect.");
  }
}

module.exports = { addSlides };
if (require.main === module) S.buildStandalone(addSlides, "part3.pptx");
