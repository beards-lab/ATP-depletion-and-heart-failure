// PART 1 of the 2026-09 status deck: intro + data + model (slides 1-9).
// Standalone build: NODE_PATH=~/tools/pptx/node_modules node part1.js  -> part1.pptx
const fs = require("fs");
const S = require("./style");

// PNG pixel size from the IHDR chunk (so S.image keeps aspect ratio). null if missing.
function pngSize(rel) {
  const abs = S.asset(rel);
  if (!fs.existsSync(abs)) return null;
  const b = fs.readFileSync(abs);
  return [b.readUInt32BE(16), b.readUInt32BE(20)];
}

// Straight arrow from (x1,y1) to (x2,y2); pptxgenjs lines need flips for up/left directions.
function arrow(s, x1, y1, x2, y2, opts = {}) {
  const x = Math.min(x1, x2), y = Math.min(y1, y2);
  const w = Math.abs(x2 - x1) || 0.0001, h = Math.abs(y2 - y1) || 0.0001;
  s.addShape("line", {
    x, y, w, h, flipH: x2 < x1, flipV: y2 < y1,
    line: { color: opts.color || S.C.ink, width: opts.width || 1.75, endArrowType: "triangle", dashType: opts.dash || "solid",
            beginArrowType: opts.both ? "triangle" : undefined },
  });
}
function label(s, text, x, y, w, h, opts = {}) {
  s.addText(text, { x, y, w, h, fontFace: S.F, fontSize: opts.fontSize || 12, color: opts.color || S.C.ink, align: opts.align || "center",
                    valign: "middle", margin: 0, italic: !!opts.italic, bold: !!opts.bold });
}
function stateBox(s, x, y, w, h, name, sub, fill, ink) {
  s.addShape("roundRect", { x, y, w, h, rectRadius: 0.08, fill: { color: fill }, line: { color: S.C.navy, width: 1.25 } });
  s.addText([{ text: name, options: { bold: true, fontSize: 18, breakLine: true } }, { text: sub, options: { fontSize: 11 } }],
            { x, y, w, h, fontFace: S.F, color: ink || S.C.navy, align: "center", valign: "middle", margin: 0 });
}

function addSlides(pres, S) {
  const { C, F, W, M, TOP } = S;

  // ───────────────────────── 1. Title ─────────────────────────
  {
    const s = pres.addSlide();
    s.background = { color: C.navy };
    s.addShape("rect", { x: 0, y: 4.55, w: W, h: 0.08, fill: { color: C.maize }, line: { color: C.maize } });
    s.addText("An ATP-dependent cross-bridge model of cardiac muscle: where we stand",
      { x: 0.8, y: 1.3, w: W - 1.6, h: 2.6, fontFace: F, fontSize: 40, bold: true, color: C.white, valign: "bottom", margin: 0 });
    s.addText("Filip Jezek, Anthony Baker, Dan Beard", { x: 0.8, y: 4.85, w: W - 1.6, h: 0.5, fontFace: F, fontSize: 22, color: C.maize, margin: 0 });
    s.addText("Lab meeting · September 2026", { x: 0.8, y: 5.4, w: W - 1.6, h: 0.45, fontFace: F, fontSize: 18, color: C.white, margin: 0 });
    s.addText("Status summary: data, model, fit, and the low-ATP effect. Limitations are stated as we go.",
      { x: 0.8, y: 6.4, w: W - 1.6, h: 0.4, fontFace: F, fontSize: 14, color: "C9D3E0", italic: true, margin: 0 });
    s.addNotes("A general status summary rather than a 'latest progress' talk. The low-ATP effect is the central story. " +
      "Several numbers shown in earlier decks and posters (uncorrected force ratios, n = 4 preps, 'reproduces all features') are corrected in this deck.");
  }

  // ───────────────────────── 2. Big picture ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "Big picture: why ATP matters", "why ATP");
    // causal chain
    const chain = [
      ["Heart failure", "energetic deficit, impaired mitochondria"],
      ["↓ [ATP], ↑ [Pi]", "↓ PCr/ATP, smaller adenine pool"],
      ["Cross-bridge cycling", "directly inhibited"],
      ["Contraction", "force, velocity, power, kinetics"],
    ];
    const bw = 2.75, bh = 1.15, gap = 0.42, y0 = TOP + 0.15;
    chain.forEach(([a, b], i) => {
      const x = M + i * (bw + gap);
      s.addShape("rect", { x, y: y0, w: bw, h: bh, fill: { color: i === 2 ? C.navy : C.light }, line: { color: C.navy, width: 1 } });
      s.addText([{ text: a, options: { bold: true, fontSize: 17, breakLine: true } }, { text: b, options: { fontSize: 13 } }],
        { x: x + 0.08, y: y0, w: bw - 0.16, h: bh, fontFace: F, color: i === 2 ? C.white : C.navy, align: "center", valign: "middle", margin: 0 });
      if (i < chain.length - 1) arrow(s, x + bw + 0.04, y0 + bh / 2, x + bw + gap - 0.04, y0 + bh / 2, { color: C.navy, width: 2.5 });
    });
    // gap statement
    s.addShape("rect", { x: M, y: 2.75, w: W - 2 * M, h: 0.85, fill: { color: "FFF6D6" }, line: { color: C.maize, width: 1 } });
    s.addText([
      { text: "The gap: ", options: { bold: true, color: C.navy } },
      { text: "Beard et al. 2022 (strain-discretised XB model) fits individual force–velocity curves, but does not reproduce the ATP-dependent differences in contractile phenotype.", options: { color: C.ink } },
    ], { x: M + 0.2, y: 2.75, w: W - 2 * M - 0.4, h: 0.85, fontFace: F, fontSize: 16, valign: "middle", margin: 0 });
    // three questions
    const qs = [
      ["Q1 · Data", "How does [ATP] modulate force, velocity and cross-bridge kinetics?"],
      ["Q2 · Mechanism", "Which cross-bridge mechanisms account for the ATP sensitivity?"],
      ["Q3 · One model", "Can a single mechanistic model reproduce force–velocity and slack–restretch across ATP levels?"],
    ];
    const cw = (W - 2 * M - 2 * 0.3) / 3;
    qs.forEach(([h, t], i) => {
      const x = M + i * (cw + 0.3), y = 3.9, ch = 2.1;
      s.addShape("rect", { x, y, w: cw, h: ch, fill: { color: C.light }, line: { color: C.light } });
      s.addShape("rect", { x, y, w: 0.08, h: ch, fill: { color: C.navy }, line: { color: C.navy } });
      s.addText(h, { x: x + 0.25, y: y + 0.12, w: cw - 0.4, h: 0.45, fontFace: F, fontSize: 17, bold: true, color: C.navy, margin: 0 });
      s.addText(t, { x: x + 0.25, y: y + 0.62, w: cw - 0.4, h: ch - 0.75, fontFace: F, fontSize: 16, color: C.ink, valign: "top", margin: 0 });
    });
    s.addText("Where we stand: Q1 answered with n = 3 preparations; Q2 and Q3 partial (8 mM fit good, low-ATP fit not yet joint).",
      { x: M, y: 6.2, w: W - 2 * M, h: 0.5, fontFace: F, fontSize: 15, italic: true, color: C.navy, margin: 0 });
    S.footer(s, "Questions as on the EMC poster (Docs/presentations/2026-EMC.pdf). Beard DA et al., Biophys J 2022.");
    s.addNotes("Kept to one slide. The three questions are the ones from the EMC poster. " +
      "Honest status line at the bottom: Q1 is essentially answered for this dataset (n = 3). Q2/Q3 are partial: the 8 mM fit is good on many features, " +
      "but no joint 8 + 2 mM fit exists yet, and the poster's 'qualitatively reproduces all features' claim is withdrawn.");
  }

  // ───────────────────────── 3. Data ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "The data we fit", "Baker lab");
    const img = "Docs/presentations/2026-09-status/assets/poster/image20.png";
    S.image(s, img, { x: M + 0.6, y: TOP, w: W - 2 * M - 1.2, h: 3.55 }, pngSize(img));
    s.addText("Slack–restretch ladder, prep 03/27 (EMC poster). Grey: 2 mM, black: 8 mM data; red/blue lines: the poster-era model fit (not the current one).",
      { x: M, y: 4.72, w: W - 2 * M, h: 0.3, fontFace: F, fontSize: 11, italic: true, color: C.mute, align: "center", margin: 0 });
    S.card(s, { x: M, y: 5.08, w: 6.05, h: 1.85 }, "Mouse RV trabeculae (Baker lab): protocols", [
      "5× slack–restretch ladder (release to ML 1.02 → 0.94, restretch at 3–7 ML/s)",
      "ktr (Brenner) protocol; passive stretches",
      "FV: borrowed Baker 2022 curve at 8 mM; in-prep FV2 ramps (2 preps)",
    ], { fontSize: 14, bar: C.navy });
    S.card(s, { x: M + 6.25, y: 5.08, w: W - 2 * M - 6.25, h: 1.85 }, "3 preps (03/27, 04/03, 04/10) · 8 vs 2 mM MgATP", [
      "Pooled 8 mM target: per-prep ÷ own pre-slack Fss22, rundown-corrected, mean of features",
      "Fref = 69.35 kPa (CV 15.5 %); weights w = 1/s² from between-prep spread",
    ], { fontSize: 14, bar: C.maize });
    S.footer(s, "Sources: Analyses/DataPooling/results/pooled_target_report.txt; Docs/experiments/pooled-metrics-2026-09.md; poster image: Docs/presentations/2026-EMC.pptx");
    s.addNotes("Species: mouse right-ventricular trabeculae (confirmed; the abstract's 'rat' was an error). " +
      "n = 3 preparations (an earlier deck said 4; one prep was dropped). " +
      "Normalisation: each prep is divided by its own 8 mM pre-slack isometric force Fss22, one divisor per prep shared by both ATP runs, so the 2/8 ratio stays a measurement. " +
      "Rundown correction explains about half of the between-prep force variance (anchor CV 22 % -> 15.5 %). " +
      "Pooled = mean of per-prep features, never features of an averaged trace. Spread is variance-moderated (n = 3 is small). " +
      "Only the slack block is pooled with error bars; FV and the ktr protocol are not pooled. The FV2 2 mM data come from only 2 preps, and the FV2 ramp is a 34 ms step, not a sustained isovelocity release. " +
      "The figure is the poster's single-prep 03/27 fit and is shown for the protocol shape only.");
  }

  // ───────────────────────── 4. ATP effect in the data ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "The ATP effect in the data", "2 mM / 8 mM");
    s.addShape("rect", { x: M, y: TOP, w: W - 2 * M, h: 0.62, fill: { color: C.maize }, line: { color: C.maize } });
    s.addText("Low ATP: stronger but slower. Force scales up, rates slow, force shape unchanged.",
      { x: M + 0.2, y: TOP, w: W - 2 * M - 0.4, h: 0.62, fontFace: F, fontSize: 19, bold: true, color: C.navy, valign: "middle", margin: 0 });
    // horizontal bar chart of (ratio - 1) in %
    const feats = [
      ["force shape (÷ level): ≈ 0", 0],
      ["rsT20  (restretch 20 % time)", 72],
      ["t0  (redevelopment delay)", 42],
      ["steady force", 33],
      ["ktr protocol", -26],
      ["ktr80  (slack)", -43],
      ["ktr  (slack)", -47],
      ["rsK80  (restretch recovery)", -62],
    ];
    const labels = feats.map((f) => f[0]);
    const up = feats.map((f) => (f[1] > 0 ? f[1] : 0));
    const dn = feats.map((f) => (f[1] < 0 ? f[1] : 0));
    s.addChart(pres.charts.BAR, [
      { name: "increase at 2 mM", labels, values: up },
      { name: "decrease at 2 mM", labels, values: dn },
    ], {
      x: M, y: 1.95, w: 7.3, h: 4.95, barDir: "bar", barGrouping: "stacked", barGapWidthPct: 45,
      chartColors: [C.warn, C.navy], valAxisMinVal: -80, valAxisMaxVal: 80, valAxisMajorUnit: 20,
      catAxisLabelFontSize: 13, catAxisLabelFontFace: F, valAxisLabelFontSize: 12, valAxisLabelFontFace: F,
      catAxisLabelPos: "low", valAxisLabelFormatCode: "0\"%\"",
      valGridLine: { color: "E5E7EB", width: 0.5 }, catGridLine: { style: "none" },
      showValAxisTitle: true, valAxisTitle: "2 mM vs 8 mM: ratio − 1 (%)", valAxisTitleFontSize: 12, valAxisTitleFontFace: F,
      showLegend: false, showValue: true, dataLabelFormatCode: '+0"%";-0"%";;', dataLabelFontSize: 12, dataLabelFontFace: F, dataLabelColor: C.white, dataLabelPosition: "inEnd",
    });
    // right-hand interpretation
    S.bullets(s, [
      { text: "Steady force ×1.33 (rundown-corrected ×1.36; pooled ×1.32 ± 8.6 %, n = 3)", bold: false },
      { text: "rsK80 ×0.38: the most ATP-discriminating feature (z = 3.6)", bold: true },
      "Headline kinetics: slack ktr ×0.53 (pooled, 3 preps; ktrfit ×0.57). The ktr protocol gives ×0.74 (reviewer point: compliance)",
      "With the protocol ktr, force × ktr ≈ 0.97: roughly constant turnover",
      "Timing: t0 ×1.42, rsT20 ×1.72",
      { text: "Supersedes earlier ×1.18–1.23 / “~20 kPa” (uncorrected)", color: C.bad },
    ], { x: 8.05, y: 1.95, w: W - 8.05 - M, h: 4.3 }, { fontSize: 15, space: 7 });
    S.footer(s, "Sources: Analyses/KtrSegmentSlope/results/atp_relevance.mat (current curation); Docs/experiments/pooled-metrics-2026-09.md §2; Analyses/ATPEffectReconciliation/conclusions.md");
    s.addNotes("These are the settled measured 2/8 mM ratios on the current curation. " +
      "Force: ×1.33 (steady); ×1.36 in the rundown-reconciliation analysis; ×1.320 ± 8.6 % in the pooled target (per prep 1.22 / 1.30 / 1.44), ~4.8 SEM from 1. " +
      "Kinetics: slack ktr ×0.53 and ktr80 ×0.57. Reviewer point (09-18): the slack ktr is compliance-contaminated; the ktr (Brenner) protocol gives ×0.72-0.74. " +
      "With the protocol value force × ktr ≈ 0.97; with the slack value it would be ≈ 0.74, which is what the July deck quoted. Which ktr probe to trust is an open decision. " +
      "rsK80 (post-restretch recovery rate, 20-80 % band) ×0.38 is the single most ATP-discriminating feature (z = 3.63). " +
      "Force shape: dividing each force feature's ratio by the level gives 1.00 ± 0.05, so ATP scales force, it does not reshape it; nine force ratios are one measurement. " +
      "Earlier decks and the poster quoted uncorrected ×1.18-1.23 / '~20 kPa higher'; those are superseded. n = 3 only.");
  }

  // ───────────────────────── 5. Model framework ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "Model framework: six-state XB model", "structure");
    // grid for the scheme
    const bw = 1.55, bh = 0.85;
    const xL = 1.0, xR = 5.1;             // left / right column x
    const y1 = 1.35, y2 = 3.25, y3 = 5.15; // SRX row / detached DRX row / attached row
    stateBox(s, xL, y1, bw, bh, "SR", "SRX · M·ATP", "E6E1F5");
    stateBox(s, xR, y1, bw, bh, "SRD", "SRX · M·ADP·Pi", "E6E1F5");
    stateBox(s, xL, y2, bw, bh, "PT", "DRX · M·ATP", C.light);
    stateBox(s, xR, y2, bw, bh, "PD", "DRX · M·ADP·Pi", C.light);
    stateBox(s, xR, y3, bw, bh, "A1", "weak, pre-stroke", "FFF1C2");
    stateBox(s, xL, y3, bw, bh, "A2", "strong, post-stroke", "FFE08A");
    const cxL = xL + bw / 2, cxR = xR + bw / 2, d = 0.13;
    // PT <-> PD (hydrolysis)
    arrow(s, xL + bw, y2 + bh / 2 - d, xR, y2 + bh / 2 - d);
    arrow(s, xR, y2 + bh / 2 + d, xL + bw, y2 + bh / 2 + d);
    label(s, "kah (hydrolysis)", xL + bw, y2 + 0.02, xR - xL - bw, 0.28);
    label(s, "kamh", xL + bw, y2 + bh - 0.28, xR - xL - bw, 0.28);
    // PD -> A1 (ka) and A1 -> PD (kd·R1D)
    arrow(s, cxR - d, y2 + bh, cxR - d, y3);
    arrow(s, cxR + d, y3, cxR + d, y2 + bh);
    label(s, "ka", cxR - 0.75, y2 + bh + 0.35, 0.55, 0.3, { align: "right" });
    label(s, "kd·R1D(s)", cxR + 0.25, y2 + bh + 0.35, 1.2, 0.3, { align: "left" });
    // A1 <-> A2
    arrow(s, xR, y3 + bh / 2 - d, xL + bw, y3 + bh / 2 - d);
    arrow(s, xL + bw, y3 + bh / 2 + d, xR, y3 + bh / 2 + d);
    label(s, "k1·R12(s)  power stroke", xL + bw, y3 + 0.02, xR - xL - bw, 0.28);
    label(s, "k_1·R21(s)", xL + bw, y3 + bh - 0.28, xR - xL - bw, 0.28);
    // A2 -> PT
    arrow(s, cxL, y3, cxL, y2 + bh, { color: C.bad, width: 2.25 });
    label(s, "k2·R2(s)", cxL + 0.15, y2 + bh + 0.12, 1.2, 0.3, { align: "left", bold: true, color: C.bad });
    label(s, "ADP release + ATP binding + detach", cxL + 0.15, y2 + bh + 0.42, 1.5, 0.55, { align: "left", fontSize: 10, color: C.bad });
    // PT <-> SR
    arrow(s, cxL - d, y2, cxL - d, y1 + bh);
    arrow(s, cxL + d, y1 + bh, cxL + d, y2, { dash: "dash", color: C.mute });
    label(s, "ksr0·e^(−F/σ)", cxL + 0.25, y1 + bh + 0.2, 1.4, 0.3, { align: "left" });
    label(s, "kmsr = 0 (dashed)", cxL + 0.25, y1 + bh + 0.5, 1.5, 0.3, { align: "left", color: C.mute });
    // SR <-> SRD
    arrow(s, xL + bw, y1 + bh / 2 - d, xR, y1 + bh / 2 - d);
    arrow(s, xR, y1 + bh / 2 + d, xL + bw, y1 + bh / 2 + d);
    label(s, "ksr2srd", xL + bw, y1 + 0.02, xR - xL - bw, 0.28);
    label(s, "ksrd2sr", xL + bw, y1 + bh - 0.28, xR - xL - bw, 0.28);
    // SRD <-> PD (force dependent)
    arrow(s, cxR - d, y1 + bh, cxR - d, y2);
    arrow(s, cxR + d, y2, cxR + d, y1 + bh);
    label(s, "kmsrd·e^(+F/σ)", cxR - 1.6, y1 + bh + 0.35, 1.4, 0.3, { align: "right" });
    label(s, "ksrd·e^(−F/σ)", cxR + 0.25, y1 + bh + 0.35, 1.35, 0.3, { align: "left" });
    // strained rupture A2 -> SR (dashed, around the left)
    s.addShape("line", { x: xL - 0.35, y: y1 + bh / 2, w: 0.0001, h: y3 + bh / 2 - (y1 + bh / 2), line: { color: C.bad, width: 1.25, dashType: "dash" } });
    s.addShape("line", { x: xL - 0.35, y: y3 + bh / 2, w: 0.35, h: 0.0001, line: { color: C.bad, width: 1.25, dashType: "dash" } });
    arrow(s, xL - 0.35, y1 + bh / 2, xL, y1 + bh / 2, { color: C.bad, width: 1.25, dash: "dash" });
    s.addText("strained rupture → SR", { x: xL - 0.55 - 1.1, y: y2 + bh / 2 - 0.15, w: 2.2, h: 0.3, fontFace: F, fontSize: 10, color: C.bad, rotate: 270, align: "center", valign: "middle", margin: 0 });
    s.addText("Attached A1/A2: distributions over strain s. Detached and SRX pools: well-mixed. F = force (kPa).",
      { x: 0.4, y: 6.2, w: 6.6, h: 0.5, fontFace: F, fontSize: 12, italic: true, color: C.mute, margin: 0 });

    // right column: mechanics + numerics
    const rx = 7.35, rw = W - rx - M;
    s.addText("Mechanical arrangement", { x: rx, y: TOP + 0.05, w: rw, h: 0.4, fontFace: F, fontSize: 17, bold: true, color: C.navy, margin: 0 });
    const mimg = "Docs/presentations/2026-09-status/assets/poster/image15.png";
    S.image(s, mimg, { x: rx + 0.2, y: TOP + 0.5, w: rw - 0.4, h: 2.1 }, pngSize(mimg));
    S.bullets(s, [
      "XB element ∥ passive titin, in series with a series elastic element (SE)",
      "Strain grid: 60 bins × 2 nm (±59 nm), moving window follows the distribution",
      "Objective: features extracted from simulated traces exactly as from the data",
      "χ²-weighted by the between-prep spread (w = 1/s²) of the pooled target",
    ], { x: rx, y: 3.85, w: rw, h: 3.0 }, { fontSize: 15, space: 6 });
    S.footer(s, "Sources: Model/dPUdT_CombinedTransitions.m; params/refit9_opt.m; mechanical diagram from the EMC poster (Docs/presentations/2026-EMC.pptx)");
    s.addNotes("Six pools: PT (M·ATP) and PD (M·ADP·Pi) detached DRX; A1 weak pre-stroke but force-bearing; A2 strong post-stroke; SR (SRX with ATP) and SRD (SRX with ADP·Pi). " +
      "A2 -> PT (red) is the lumped ADP release + ATP binding + detachment step, i.e. the natural place for an ATP gate; in the incumbent no ATP gate is on. " +
      "SRX exit kmsr is 0 in the incumbent, so SRX is left only through the ADP branch (SR -> SRD -> PD). " +
      "The strained-rupture branch sends a fraction of the high-strain R2 flux into SR instead of PT, only during a restretch. " +
      "Strain grid numbers (60 × 2 nm) are from the project plan; the window is moved by an event in evaluateModel (movingWindow). " +
      "The objective is feature-based, not trace MSE; features are read from model traces by the same extractor as the data.");
  }

  // ───────────────────────── 6. Mechanisms I: states & kinetics ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "Mechanisms ON (1/2): states and kinetics", "refit9_opt");
    const rows = [
      ["Mechanism", "What it does", "Key value (refit9)"],
      ["Six pools", "PT, PD, A1, A2 (strain-resolved), SR, SRD", "—"],
      ["Hydrolysis PT ⇄ PD", "cocks the detached head", "kah 175.7, kamh 8.80 s⁻¹"],
      ["Attachment PD → A1", "triangular kernel around zero strain", "ka 292.2 s⁻¹, half-width 6 nm"],
      ["Strain-dependent R1D, R12, R21, R2", "rate = k × PCHIP through movable knots", "k: kd 92.6, k1 160.8, k_1 31.4, k2 116.8 s⁻¹"],
      ["Power stroke; pre-stroke force", "A2 shifted by dr; A1 bears force at dr1", "dr 11.26 nm, dr1 4.79 nm"],
      ["Force-dependent SRX entry", "PT → SR, suppressed by force; no direct exit", "ksr0 367.5·e^(−F/27.3 kPa), kmsr = 0"],
      ["SRX-ADP branch", "SR ⇄ SRD; SRD ⇄ PD force-dependent", "130.2 / 9.21 s⁻¹; kmsrd 26.8·e^(+F/22.8), ksrd 10·e^(−F/37.6)"],
      ["Strained rupture → SRX", "part of high-strain R2 flux lands in SR", "21.9 % above 13.9 nm"],
    ];
    S.table(s, rows, { x: M, y: TOP + 0.1, w: W - 2 * M, rowH: 0.58 }, { fontSize: 14, colW: [3.3, 4.3, 4.73] });
    S.footer(s, "Source: params/refit9_opt.m (values rounded); Model/dPUdT_CombinedTransitions.m. Rates in s⁻¹, F in kPa, strains in nm.");
    s.addNotes("All values read from params/refit9_opt.m. " +
      "The four strain-dependent attached-state rates are the most flexible part of the model: each is a base rate times a PCHIP through movable knots (next-but-one slide). " +
      "The A1 attachment is spread over a triangular kernel of half-width 6 nm (A1AttachmentWidth 0.006 µm; implemented in the inline branch of the RHS, the separate UseA1AttachmentKernel flag is off). " +
      "SRX: entry falls with force (Frank-Starling-like mechanosensing); direct exit kmsr is 0, so return goes through SRD. " +
      "Strained rupture: fraction ramps 0 -> 21.9 % over 13.9-18.9 nm; exists only while the fibre is stretched, i.e. in restretch. " +
      "Hydrolysis kah is NOT ATP-gated in the current thinking (the poster's kah gate was not physiological).");
  }

  // ───────────────────────── 7. Mechanisms II: mechanics & regulation ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "Mechanisms ON (2/2): mechanics, regulation", "refit9_opt");
    const rows = [
      ["Mechanism", "What it does", "Key value (refit9)"],
      ["Non-Hookean XB stiffness", "force ∝ strain^estiff; split +/− strain stiffness", "estiff 1.243; kstiff2 75932 / kstiff2_n 45421"],
      ["Filament overlap (+ factor)", "attachment ∝ thick–thin overlap", "L_thick 1.195, L_thin 1.174 µm"],
      ["Lattice spacing", "Gaussian penalty on filament spacing", "d_opt 23.4 nm, σ 2.69 nm"],
      ["Global occupancy ceiling", "attachment × (1 − P_bound / P_max)", "P_bound_max 1.208"],
      ["Registration availability A_reg", "sites lost during shortening, recover", "A0 0.210, v_ref 1.025 µm/s, τ 0.205 ms"],
      ["Load-realign deficit A_def", "sites unavailable after load stripping", "τ_cr 27.5 ms, F_cr 23.4 kPa"],
      ["Series element + viscosity", "softening SE, F = (kSE·LSE)^0.845, damper μ", "kSE 4723, ekSE 0.845"],
      ["Passive titin", "power law above slack length", "k_pas·(Lc − 1.80 µm)^7.67"],
      ["Titin Maxwell branch", "viscoelastic, tension-only", "η 13.2 ms"],
    ];
    S.table(s, rows, { x: M, y: TOP + 0.05, w: W - 2 * M, rowH: 0.48 }, { fontSize: 14, colW: [3.3, 4.6, 4.43] });
    s.addShape("rect", { x: M, y: 6.3, w: W - 2 * M, h: 0.6, fill: { color: "FDECEA" }, line: { color: C.bad, width: 0.75 } });
    s.addText("Available but OFF in the incumbent: all ATP/ADP/Pi gates, Ca²⁺ activation, catch bond, stretch activation, A2 reattachment, yield gate (in progress) …",
      { x: M + 0.15, y: 6.3, w: W - 2 * M - 0.3, h: 0.6, fontFace: F, fontSize: 14, color: C.bad, valign: "middle", margin: 0 });
    S.footer(s, "Source: params/refit9_opt.m; Analyses/MechanismAblation/conclusions.md; Analyses/OverstretchedHeads/conclusions.md");
    s.addNotes("Values from params/refit9_opt.m. A_reg acts on shortening only (FV shoulder). A_def is the compliant-realignment gate, ka × (1 − A_def). " +
      "Lattice: bites only on shortening (inert at SL 2.2; attenuates attachment below SL ≈ 2.09 µm); its parameters carry only two identifiable combinations (MechanismAblation §7.1). " +
      "Occupancy ceiling: linear form, f = 1 − P_bound/P_bound_max. " +
      "MechanismAblation (35 arms at 8 mM): nothing, added or removed, improved the incumbent beyond noise; the tension-only Maxwell switch and A2 reattachment were flagged as free deletions. " +
      "Everything ATP-related is off: the 8 mM incumbent has no ATP dependence at all. Low-ATP is part 3.");
  }

  // ───────────────────────── 8. Transition rates ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "Transition rates of the incumbent", "refit9_opt");
    const fig = "Analyses/LowATP_Refit9K2Sweep/results/refit9_transitions.png";
    const nw = W - 2 * M - 3.25;
    S.image(s, fig, { x: M, y: TOP, w: nw, h: 5.7 }, pngSize(fig));
    const nx = M + nw + 0.2;
    S.card(s, { x: nx, y: TOP, w: W - nx - M, h: 5.7 }, "Piecewise method", [
      "Rate = k × PCHIP through movable knots, 4–8 knots per rate",
      "Each knot set shifts via a log-offset",
      "Clamped beyond the outer knots",
      "Very flexible: 32 of 79 pool entries shape these curves",
      { text: "→ identifiability cost", bold: true, color: C.bad },
    ], { fontSize: 14, bar: C.maize });
    S.footer(s, "Figure: Analyses/LowATP_Refit9K2Sweep/results/refit9_transitions.png; knots: params/refit9_opt.m; Model/getParams.m");
    s.addNotes("Each strain-dependent attached rate (R1D, R12, R21, R2) is a base rate times a PCHIP interpolant through knots whose positions and values the optimiser can move. " +
      "Knot counts in refit9: R12 5 (two outer knots far off-grid, 3 inside), R2 8 (evaluated at s + 20.3 nm), R1D 4, R21 4. Value at s = 0: R1D 115, R12 627, R21 202, R2 88.5 s^-1; max on grid R12 9316, R2 5840 s^-1 (refit9_transitions.txt). Knot sets shift through a log-offset (applied for R2, R1D and R21; NOT for R12, see next slide). " +
      "Extrapolation beyond the outer knots is clamped at the edge values. " +
      "32 of the 79 lean-pool entries are PCHIP knot coordinates/values or offsets (28 knots + 4 offsets). The plan said about 36; my count from leanPool in refit9_opt.m is 32. " +
      "This flexibility is what lets the model fit, and it is also why individual knot values are poorly identified.");
  }

  // ───────────────────────── 9. Parameter accounting ─────────────────────────
  {
    const s = pres.addSlide();
    S.header(s, "Honest parameter accounting", "identifiability");
    S.stat(s, M, TOP + 0.05, "79", "entries in the optimiser pool (leanPool)");
    S.stat(s, M + 3.1, TOP + 0.05, "32", "of them shape the PCHIP rates (28 knots + 4 offsets)");
    S.stat(s, M + 6.2, TOP + 0.05, "12", "moved per round: 10 compulsory + 2 drawn from the pool");
    S.bullets(s, [
      "Lattice geometry: 6 parameters, only 2 identifiable combinations",
      "Knots: flexible but weakly identified",
      "With n = 3 each χ² weight is known only to a factor of ~18: structural refit now, hard identification at n ≥ 6",
      "Costs compare only within one objective: refit9_opt = 0.2079 own / 0.424 robust / 0.558 yield",
    ], { x: M, y: 2.85, w: 6.9, h: 3.2 }, { fontSize: 16, space: 10 });
    S.card(s, { x: 7.75, y: 2.85, w: W - 7.75 - M, h: 3.85 }, "One pool entry is inert", [
      "PieceWiseStrainDepX_logOffset (R12 knot offset) is never read by the code",
      "Added as a default in Apr 2026; getParams builds R12 without it",
      "The optimiser has been moving a dead parameter (refit9: 0.9647)",
      "Not fixed: wiring it in would shift R12 by −3.6 nm and change the fit",
      { text: "Checked: offset 0.9647 vs 0.5 gives identical R12 (max diff 0 s⁻¹)", bold: true },
    ], { fontSize: 14, bar: C.bad, headColor: C.bad });
    S.footer(s, "Sources: params/refit9_opt.m; Analyses/KtrSegmentSlope/RunRefit9Optim.m; MechanismAblation §7.1; pooled-metrics-2026-09.md; LowATP_Refit9K2Sweep/results/r12_offset_check.txt");
    s.addNotes("79 = length of params0.leanPool in refit9_opt.m. Zero-valued base parameters are dropped from the pool at run time. " +
      "Optimiser design (RunRefit9Optim.m): N_DRAW = 12 per round INCLUDING the 10 compulsory levers (k2, ka, ksrd, sigma_srd2, R2 log-offset, ekSE, dr2, L_thin, k_cr, tau_cr), so only 2 are drawn at random each round " +
      "(optimizeFeatures.m: k = N_DRAW − numel(compulsory)). Note: the brief said '12 sampled + 10 compulsory'; the code says 12 total. " +
      "Lattice: d10_ref·sqrt(SL_ref) is a product and R_thick + R_thin + d_optimal a sum, so 6 parameters reduce to 2 combinations (MechanismAblation §7.1). " +
      "Inert R12 offset: default 1 in getParams.m:269; R12 built as pchip(pwsdX, pwsdP) at getParams.m:474 without the offset; both RHS paths and the MEX pack read only that interpolant; R2, R1D and R21 do apply their offsets. " +
      "Added in commit 869870c (28 Apr 2026). Deliberately not fixed. Numerical check (r12_offset_check.txt): getParams with 0.9647 vs 0.5 gives bit-identical R12 caches; positive control on the R2 offset moves R2 by up to 5652 s^-1, so the check is sensitive. " +
      "Costs: 0.2079 is refit9 on its own objective; 0.424 on the robust-reader objective; 0.558 on the yield objective. They are not comparable across objectives, and earlier campaign costs (pre-noise-fix) are not comparable at all.");
  }
}

module.exports = { addSlides };
if (require.main === module) S.buildStandalone(addSlides, "part1.pptx");
