// PART 2 of the 2026-09 status deck: fitting (slides 10-14 + backup "campaign history").
// Build standalone:  NODE_PATH=~/tools/pptx/node_modules node part2.js   -> part2.pptx
const S = require("./style");

const POSTER = "Docs/presentations/2026-09-status/assets/poster/";
const KSS = "Analyses/KtrSegmentSlope/results/";
const OH = "Analyses/OverstretchedHeads/results/";

// small text helper (keeps Arial + margins consistent)
function txt(s, S, text, box, o = {}) {
  s.addText(text, { fontFace: S.F, fontSize: 14, color: S.C.ink, margin: 0.03, valign: "top", ...box, ...o });
}

// coloured status dot
function dot(s, S, x, y, color, d = 0.2) {
  s.addShape("ellipse", { x, y, w: d, h: d, fill: { color }, line: { color } });
}

function addSlides(pres, S) {
  const { C, F, M, W } = S;

  // ---------------------------------------------------------------- 10. poster-era fit
  {
    const s = pres.addSlide();
    S.header(s, "The poster-era fit (EMC poster, Jun–Sep 2026)", "Fitting · where we were");

    S.card(s, { x: M, y: 1.2, w: 4.9, h: 2.45 }, "What it was", [
      "spreadw core + ATP power gates on k2 (p = 0.44) and hydrolysis kah (p = 0.30); one prep (03/27); 8 mM invariant by construction",
      "Separated 8 vs 2 mM on F–V and power, and on the slack features qualitatively",
    ], { bar: C.navy, fontSize: 14 });

    S.card(s, { x: M, y: 3.8, w: 4.9, h: 3.1 }, "Why we moved on", [
      { text: "Single prep: between-prep force CV 15.5 % (79.8 / 70.0 / 58.3 kPa)", color: C.bad },
      { text: "Gating hydrolysis (kah) by [ATP] is not defensible chemistry", color: C.bad },
      { text: "Objective was numerically noisy (slide 14)", color: C.bad },
      { text: "\"Reproduces all features\" overstated it: post-restretch 80 ms window and staircase were not fit", color: C.bad },
    ], { bar: C.bad, fontSize: 14 });

    // poster panels: F-V + power on top, three slack features below
    const X0 = 5.65, TW = W - M - X0; // 7.18 in
    S.image(s, POSTER + "image18.png", { x: X0, y: 1.15, w: TW / 2 - 0.05, h: 2.35 }, [1364, 885]);
    S.image(s, POSTER + "image19.png", { x: X0 + TW / 2 + 0.05, y: 1.15, w: TW / 2 - 0.05, h: 2.35 }, [1368, 886]);
    const pw = (TW - 0.2) / 3;
    [["image22.png", [1364, 886]], ["image24.png", [1369, 886]], ["image26.png", [1369, 886]]].forEach(([f, px], i) => {
      S.image(s, POSTER + f, { x: X0 + i * (pw + 0.1), y: 3.75, w: pw, h: 1.65 }, px);
    });
    // legend for the feature panels (the poster panels d/f carry none)
    const ly = 5.5;
    const leg = [["000000", "data 8 mM"], ["808080", "data 2 mM"], ["1F5FBF", "model 8 mM"], ["D03A1A", "model 2 mM"]];
    leg.forEach(([col, lab], i) => {
      const x = X0 + 0.1 + i * 1.78;
      s.addShape("line", { x, y: ly + 0.14, w: 0.35, h: 0, line: { color: col, width: 3 } });
      txt(s, S, lab, { x: x + 0.42, y: ly, w: 1.35, h: 0.3 }, { fontSize: 12, color: C.mute });
    });
    txt(s, S, "Poster panels: F–V, power–velocity; steady force, slack ktr, rsK per slack cycle (bands = SD of 3 preps, not fitted). " +
      "The ATP separation is right in sign; the 2 mM rsK and ktr drops are under-reproduced (poster’s own captions).",
      { x: X0, y: 5.9, w: TW, h: 0.95 }, { fontSize: 13, color: C.ink, italic: true });

    S.footer(s, "Sources: 2026-EMC poster; ATP-wt-lowatp/Analyses/AtpEffectHypotheses/conclusions.md (H1 gates); DataPooling/results/pooled_target_report.txt (CV)");
    s.addNotes(
      "The fit on the EMC poster: the spreadw core plus two ATP power gates, k2 exponent 0.44 and kah exponent 0.30 (params atphyp_H1_opt). " +
      "Every gate is exactly 1 at 8 mM, so the 8 mM fit was unchanged by construction. On that core those two gates captured 99.5 % of the 2 mM cost span (AtpEffectHypotheses). " +
      "It was one preparation, 03/27. " +
      "Why we moved on: (1) one prep, while the three preps differ by CV 15.5 % in plateau force; " +
      "(2) kah is hydrolysis, which happens on the detached head and whose ATP dependence is not a defensible mechanism for a 2 vs 8 mM effect; " +
      "(3) the objective was noisy at the ~20 % level (see slide 14), so small cost differences at the time meant little; " +
      "(4) the poster said 'reproduces all key features' - that overstated it. The post-restretch 80 ms window (rsK) and the held-out staircase were not fit; " +
      "the poster's own captions say rsK and ktr drops at 2 mM are under-reproduced. " +
      "Species: mouse (the abstract's 'rat' was an error).");
  }

  // ---------------------------------------------------------------- 11. what changed
  {
    const s = pres.addSlide();
    S.header(s, "What changed since the poster", "Fitting · before / after");

    const x0 = M, labW = 2.2, colW = (W - 2 * M - labW - 0.2) / 2;
    const xB = x0 + labW + 0.1, xA = xB + colW + 0.1;
    // column heads
    s.addShape("rect", { x: xB, y: 1.2, w: colW, h: 0.5, fill: { color: "9CA3AF" }, line: { color: "9CA3AF" } });
    txt(s, S, "Poster era (Jun)", { x: xB, y: 1.2, w: colW, h: 0.5 }, { fontSize: 17, bold: true, color: C.white, align: "center", valign: "middle" });
    s.addShape("rect", { x: xA, y: 1.2, w: colW, h: 0.5, fill: { color: C.navy }, line: { color: C.navy } });
    txt(s, S, "Now (refit9, Sep)", { x: xA, y: 1.2, w: colW, h: 0.5 }, { fontSize: 17, bold: true, color: C.white, align: "center", valign: "middle" });

    const rows = [
      ["Target",
        "One prep (03/27), absolute kPa",
        "Pooled 3-prep 8 mM target, Fref 69.35 kPa; χ² weights 1/s² from the between-prep spread. Pooled 2 mM target: being merged in now"],
      ["ATP gating",
        "Power gates on k2 and kah (hydrolysis)",
        "Physiological only: A2 detachment (k2 via R2). Hydrolysis not gated. The 8 mM incumbent has all gates OFF"],
      ["Numerics",
        "ode15s RelTol 1e-2; features read off solver steps (noise 19.5 %)",
        "Smooth titin rectifier, RelTol 1e-3, clamped extrapolation, reader on the data grid: noise 0.17 %"],
      ["Objective",
        "Single-prep feature cost on 03/27 (weights from a moderated spread)",
        "+ ktrfit & ktr ladder, rsK80 per segment, ktr-protocol k80, rsTfitE; force shape as guards with deadbands"],
      ["Mechanisms",
        "spreadw core",
        "Registration availability (A_reg), load-realign deficit (A_def), occupancy ceiling, lattice, tension-only Maxwell titin"],
    ];
    const y0 = 1.8, rh = 0.98;
    rows.forEach(([lab, b, a], i) => {
      const y = y0 + i * (rh + 0.05);
      s.addShape("rect", { x: x0, y, w: labW, h: rh, fill: { color: C.light }, line: { color: C.light } });
      s.addShape("rect", { x: x0, y, w: 0.08, h: rh, fill: { color: C.maize }, line: { color: C.maize } });
      txt(s, S, lab, { x: x0 + 0.18, y, w: labW - 0.2, h: rh }, { fontSize: 16, bold: true, color: C.navy, valign: "middle" });
      s.addShape("rect", { x: xB, y, w: colW, h: rh, fill: { color: "F3F4F6" }, line: { color: "E5E7EB" } });
      txt(s, S, b, { x: xB + 0.1, y, w: colW - 0.2, h: rh }, { fontSize: 14, color: C.mute, valign: "middle" });
      s.addShape("rect", { x: xA, y, w: colW, h: rh, fill: { color: C.white }, line: { color: C.navy, width: 1 } });
      txt(s, S, a, { x: xA + 0.1, y, w: colW - 0.2, h: rh }, { fontSize: 14, color: C.ink, valign: "middle" });
      // arrow between columns
    });

    S.footer(s, "Sources: DataPooling/results/pooled_target_report.txt; StrainWindowNoise, FeatureRobustness conclusions.md; KtrSegmentSlope/dashboard/decisions.jsonl");
    s.addNotes(
      "Five things changed. Target: we moved from one preparation to a pooled 3-prep 8 mM target (Fref 69.35 kPa, per-prep 79.8/70.0/58.3), " +
      "with each feature weighted 1/s^2 by its between-prep spread, so every term is a chi-square. The pooled 2 mM target exists on the floored-campaign branch and is being merged now; the 2 mM fit on the current core is not done yet. " +
      "ATP gating: only A2 detachment feels ATP; hydrolysis is not gated. The 8 mM incumbent refit9 runs with every ATP gate off. " +
      "Numerics: see slide 14 - the old objective had ~20 % noise from how features were read. " +
      "Objective: new terms for the ktr ladder (band fit ktrfit + shape), rsK80 per segment, the ktr protocol k80, and a robust recovery delay rsTfitE; the ATP-flat force shape terms are guards. " +
      "Caveat on 'guards': a guard pays nothing inside a deadband set at 1.2x the seed's own cost, which hid the re-stretch timing misfit (slide 13). " +
      "Mechanisms: the list is what is ON in refit9; some of these (registration availability, occupancy ceiling, lattice) were already being added in June around the poster, so this row is 'now in the incumbent' rather than strictly 'new since'.");
  }

  // ---------------------------------------------------------------- 12. the 8 mM fit now
  {
    const s = pres.addSlide();
    S.header(s, "The 8 mM fit now: refit9", "Fitting · incumbent");
    // figure fills the left ~60 %
    S.image(s, KSS + "refit9_protocols.png", { x: 0.3, y: 1.1, w: 7.75, h: 5.85 }, [1831, 1395]);

    const xr = 8.25, wr = W - M - xr;
    S.card(s, { x: xr, y: 1.2, w: wr, h: 2.0 }, "Incumbent: params/refit9_opt.m", [
      "0.2079 own objective · 0.424 robust · 0.558 yield",
      { text: "Costs compare only within one objective", color: C.bad, bold: true },
    ], { bar: C.navy, fontSize: 14 });

    S.card(s, { x: xr, y: 3.35, w: wr, h: 3.55 }, "Panels", [
      "A slack vs pooled trace (3 preps) · B F–V (on the Baker22 curve)",
      "C–E ktr protocol, 3 days · F ktr80 ladder",
      "G rsK80 ladder (seg 5 not scored) · I cost breakdown",
      "H staircase: held out, single prep, visual check only",
    ], { bar: C.maize, fontSize: 14 });

    S.footer(s, "Figure: KtrSegmentSlope/results/refit9_protocols.png. Costs: checkin_refit9.json; FeatureRobustness (0.424), OverstretchedHeads (0.558) conclusions.md");
    s.addNotes(
      "The current 8 mM incumbent, refit9_opt, one round from refit9_seed (0.2422 -> 0.2079), and the first banked step that re-rendered to exactly the same number. " +
      "Three numbers for the same parameters: 0.2079 on its own objective; 0.424 on the robust objective (same parameters, robust reader, rsTfitE replacing rsT20, guards re-derived); " +
      "0.558 on the yield objective where the re-stretch timing trio is floored honestly instead of guarded. These are not comparable across objectives - never read 0.2079 vs 0.558 as a regression. " +
      "Panel A is the slack protocol against the pooled mean trace; B F-V against the Baker22 curve (borrowed, not in-prep); C-E the ktr protocol on each of the 3 days; F and G the two ladders that are the main misses; " +
      "H the staircase is held out and only a single prep exists, so it is a visual check.");
  }

  // ---------------------------------------------------------------- 13. scorecard + restretch diagnosis
  {
    const s = pres.addSlide();
    S.header(s, "8 mM scorecard, and what is wrong at re-stretch", "Fitting · refit9");

    const rows = [
      ["", "Feature", "Model", "Data"],
      [C.good, "F–V", "on Baker22 curve", "—"],
      [C.good, "Steady force", "71.5 kPa", "68.5–69.1"],
      [C.good, "ktr level (mean)", "45.2 s⁻¹", "45.5"],
      [C.good, "ktr protocol k80", "46.2 s⁻¹", "43.6"],
      [C.good, "rsK80 mean", "47.8 s⁻¹", "47.5"],
      [C.warn, "Force-shape guards", "pay 0 (in deadband)", "—"],
      [C.bad, "ktr ladder", "flat, −0.0 %/seg", "−4.1 %/seg"],
      [C.bad, "rsK80 across seg 1–4", "44 → 54", "flat ~47"],
      [C.bad, "rsT20 (delay)", "5.5–6.7 ms", "8.2–9.2 ms"],
      [C.bad, "Post-peak fall / peak", "2× slow; peak 1.17–1.35× late", "1×"],
      [C.bad, "Staircase (held out)", "peaks 1.15–1.28", "1.05–1.12"],
    ];
    const tx = M, ty = 1.2, colW = [0.4, 2.55, 2.6, 1.35];
    const rh = 0.44;
    const tw = colW.reduce((a, b) => a + b, 0);
    rows.forEach((r, i) => {
      const y = ty + i * rh;
      const head = i === 0;
      s.addShape("rect", { x: tx, y, w: tw, h: rh, fill: { color: head ? C.navy : (i % 2 ? C.white : C.light) }, line: { color: "D1D5DB", width: 0.5 } });
      if (!head) dot(s, S, tx + 0.11, y + 0.12, r[0]);
      let x = tx + colW[0];
      for (let k = 1; k < 4; k++) {
        txt(s, S, r[k], { x: x + 0.05, y, w: colW[k] - 0.1, h: rh }, {
          fontSize: 14, valign: "middle", bold: head || k === 1, color: head ? C.white : (k === 1 ? C.ink : C.ink),
        });
        x += colW[k];
      }
    });
    // legend
    const lgY = ty + rows.length * rh + 0.08;
    [[C.good, "fits"], [C.warn, "fits by construction"], [C.bad, "misses"]].forEach(([c, l], i) => {
      dot(s, S, tx + 0.1 + i * 2.2, lgY + 0.06, c, 0.18);
      txt(s, S, l, { x: tx + 0.35 + i * 2.2, y: lgY, w: 1.9, h: 0.3 }, { fontSize: 12, color: C.mute });
    });

    // right: diagnosis figure + text
    const xr = tx + tw + 0.25, wr = W - M - xr;
    S.image(s, OH + "restretch_dwell_refit9_vs_yieldB.png", { x: xr, y: 1.12, w: wr, h: 3.3 }, [1326, 1007]);
    txt(s, S, "Re-stretch 1 and 4, data 03/27 (grey), refit9 (blue), yield gate + rebinding (orange). Rows: force; heads beyond +15 nm; % XB force from them.",
      { x: xr, y: 4.42, w: wr, h: 0.5 }, { fontSize: 11, color: C.mute, italic: true });
    S.bullets(s, [
      "Weak A1 heads dragged to 20–30 nm for ~4 ms, only while lengthening",
      "Yield gate / A1 actin-site slip fixes valley timing and staircase nadir",
      "Slip: 0.501 vs refit9 0.558 (yield objective), before any refit",
    ], { x: xr, y: 4.95, w: wr, h: 1.5 }, { fontSize: 14, space: 3 });
    s.addShape("rect", { x: xr, y: 6.52, w: wr, h: 0.4, fill: { color: "FFF4E0" }, line: { color: C.warn, width: 1 } });
    txt(s, S, "IN PROGRESS — not the incumbent (campaign not run)", { x: xr, y: 6.52, w: wr, h: 0.4 },
      { fontSize: 13, bold: true, color: C.warn, align: "center", valign: "middle" });

    S.footer(s, "Sources: OverstretchedHeads/results/evaluate.log + conclusions.md; KtrSegmentSlope/results/render_refit9.json + conclusions.md §8; FeatureRobustness/conclusions.md §5");
    s.addNotes(
      "Scorecard on refit9_opt against the pooled 8 mM target. " +
      "Fits: F-V on the Baker22 curve (a borrowed curve, not in-prep); steady ~71.5 kPa vs 68.5-69.1 (segments 1-4, about +4 %); ktr80 mean 45.2 vs 45.5; ktr protocol k80 46.2 vs 43.6 (+6 %, within spread); rsK80 mean 47.8 vs 47.5. " +
      "The force-shape guards pay zero, but by construction: their deadbands are 1.2x the seed's own cost, which is how the re-stretch timing misfit was parked. Hence amber. " +
      "Misses: the ktr ladder is flat where the data fall -4.1 %/segment; rsK80 climbs 44 -> 54 over segments 1-4 where the data are flat ~47; rsT20 is 5.5-6.7 ms vs 8.2-9.2 ms (~25 % short); " +
      "the post-peak fall is ~2x too slow (tFall50 model/data 1.72-2.00) and the overshoot peaks late (peak1_dSLq model/data 1.17-1.35); the held-out staircase overshoots each 1 % step (1.15-1.28 vs 1.05-1.12). " +
      "Diagnosis (OverstretchedHeads): the long-held heads are weak A1 heads in a p1<->p2 shuttle, dragged to 20-30 nm for ~4 ms, only during lengthening; no existing parameter moves the timing (3 h control refit changed the timing terms to 3 digits). " +
      "A stretched-side yield gate with fast rebinding puts A1 detachment at 6-8 nm and fixes vall_t and the staircase nadir, but the gated refit (yieldB) cost 0.763 vs control 0.529: valley too low, recovery too fast, step amplitude worse. " +
      "A1 slip to the neighbouring actin site, added on top of yieldB without refitting, scores 0.501 - below the control 0.529 and refit9's 0.558 on the same objective. " +
      "A three-arm campaign is built and pre-flighted but not run. This is in progress, not the incumbent. Held-out staircase: slip/yield effect on step amplitude still open.");
  }

  // ---------------------------------------------------------------- 14. optimisation noise
  {
    const s = pres.addSlide();
    S.header(s, "Optimisation noise, and how we resolved it", "Fitting · numerics");

    const lx = M, lw = 7.35;
    const steps = [
      [C.bad, "Symptom", "refit7’s accepted 0.3355 re-scored 0.3707 from a snapshot 1.1e-6 away; 9 evaluations 1e-6 apart spanned 0.3355–0.4160. refit8 measured 19.5 % noise and stalled. Noisy terms were crossing-time features: a ~1 ms trace shift moved one rsK80 by 4 % → 2× its cost."],
      [C.mute, "Wrong first suspect", "Strain-window shifts. They delete ≤ 1e-5 of attached mass. “Wider grid breaks it” was negative cubic extrapolation + an off-centre grid."],
      [C.warn, "Real causes", "Features read off ode15s step points at hard-coded RelTol 1e-2; tighter tolerance stalled on a corner in the titin Maxwell rectifier."],
      [C.good, "Fix 1: numerics", "Smooth rectifier (ε 0.01), RelTol 1e-3, 0.5 ms output grid, clamped extrapolation: 19.5 % → 0.3 %."],
      [C.good, "Fix 2: reader", "Model was read at 1 kHz, data at 4 kHz (restretchSlopeStart biased 27 %). Read on the data grid: ktr80 noise 6.25 % → 0. rsT20 (argmin-timed, noise in the data too) → rsTfitE; vertex features. 0.17 % at 1.16× run time."],
    ];
    let y = 1.18;
    const hs = [1.25, 0.95, 0.8, 0.8, 1.25];
    steps.forEach(([col, head, body], i) => {
      const h = hs[i];
      s.addShape("rect", { x: lx, y, w: lw, h, fill: { color: C.light }, line: { color: C.light } });
      s.addShape("rect", { x: lx, y, w: 0.08, h, fill: { color: col }, line: { color: col } });
      txt(s, S, head, { x: lx + 0.18, y: y + 0.02, w: 1.75, h: h - 0.04 }, { fontSize: 15, bold: true, color: C.navy, valign: "middle" });
      txt(s, S, body, { x: lx + 1.95, y: y + 0.02, w: lw - 2.05, h: h - 0.04 }, { fontSize: 14, valign: "middle" });
      y += h + 0.07;
    });

    // native bar chart (log scale)
    const cx = lx + lw + 0.25, cw = W - M - cx;
    s.addChart(pres.charts.BAR, [{
      name: "Objective noise (rel. SD, %)",
      labels: ["Before", "Fix 1", "Fix 2"],
      values: [19.5, 0.3, 0.17],
    }], {
      x: cx, y: 1.15, w: cw, h: 3.65, barDir: "col",
      chartColors: [C.navy], barGapWidthPct: 60,
      valAxisLogScaleBase: 10, valAxisMinVal: 0.1, valAxisMaxVal: 100, valAxisCrossesAt: 0.1,
      showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: '[<1]0.0#"%";0.0"%"', dataLabelFontSize: 16, dataLabelColor: C.ink, dataLabelFontBold: true,
      showTitle: true, title: "Objective noise, 1e-6 perturbations (log)", titleFontFace: F, titleFontSize: 14, titleColor: C.ink,
      catAxisLabelFontFace: F, catAxisLabelFontSize: 14, catAxisLabelColor: C.ink,
      valAxisLabelFontFace: F, valAxisLabelFontSize: 11, valAxisLabelColor: C.mute, valAxisLabelFormatCode: '[<1]0.0"%";0"%"',
      valGridLine: { color: "E5E7EB", size: 0.5 }, catGridLine: { style: "none" }, showLegend: false,
    });
    txt(s, S, "Before / Fix 1 at refit8_opt; Fix 2 at refit9_opt params (robust1_seed), 4 neighbours.",
      { x: cx, y: 4.8, w: cw, h: 0.45 }, { fontSize: 11, color: C.mute, italic: true });

    s.addShape("rect", { x: cx, y: 5.3, w: cw, h: 1.6, fill: { color: "FDECEA" }, line: { color: C.bad, width: 1 } });
    txt(s, S, "Consequence", { x: cx + 0.15, y: 5.35, w: cw - 0.3, h: 0.35 }, { fontSize: 15, bold: true, color: C.bad });
    txt(s, S, "Pre-fix costs are not comparable: refit8_opt ~0.14 → 0.356 re-scored. refit8’s “progress” was noise.",
      { x: cx + 0.15, y: 5.72, w: cw - 0.3, h: 1.1 }, { fontSize: 14 });

    S.footer(s, "Sources: Analyses/StrainWindowNoise/conclusions.md; Analyses/FeatureRobustness/conclusions.md; KtrSegmentSlope/dashboard/decisions.jsonl (09-23 16:36, 17:18, 22:56)");
    s.addNotes(
      "Symptom: in refit7 the optimiser accepted 0.3355; the snapshot written to disk (parameters equal to 1.1e-6, %g rounding) re-scored 0.3707. Nine evaluations at 1e-6 perturbations spanned 0.3355-0.4160 - the accepted point was a winner's curse. " +
      "refit8 put measured numerical noise into the floors and acceptance, measured 19.5 % objective noise, and stalled. The noisy terms were crossing-time features: rsK80, rsT20, the ktr shape. A ~1 ms shift of the trace moved segment-4 rsK80 52.3 -> 54.3 (4 %), and the tight rsK80 term turned that into about twice the cost. " +
      "Wrong first suspect: the moving strain window. Measured: its shifts delete at most 1e-5 of attached mass. 'Widening the grid breaks the model' was two bugs: cubic extrapolation of the rate pchips going negative (R2 = -460/s at +150 nm), and a cap that built an off-centre grid. " +
      "Real causes: features were read off ode15s step points at a hard-coded RelTol 1e-2; tightening it stalled because the tension-only titin Maxwell element has a max(0,.) corner exactly at isometric rest. " +
      "Fix 1: SiLU rectifier eps 0.01, RelTol 1e-3, fixed 0.5 ms output, clamped extrapolation: 19.5 % -> 0.3 % on refit8_opt. " +
      "Fix 2 (the reader): the model was resampled at 1 kHz on a grid anchored on its own first sample while the data are read at 4 kHz; that alone biased restretchSlopeStart 27 %. Reading on the data grid took ktr80 noise 6.25 % -> 0. rsT20 is timed from an argmin and is as noisy in the data as between hearts, so it was replaced by rsTfitE; vertex (sub-sample) peak/valley timing replaced step-function sample features. Full objective 0.17 % at 1.16x the old run time. " +
      "Caveat: the three bars are not at one parameter point (refit8_opt for the first two, refit9 params for the third). Consequence: no cost from before 2026-09-24 can be compared to one after.");
  }

}

// Backup slide, exported separately so build_all can put it after the main story.
function addBackup(pres, S) {
  const { C, W, M } = S;
  // ---------------------------------------------------------------- Backup: campaign history
  {
    const s = pres.addSlide();
    S.header(s, "Backup: campaign history (refit4–refit9)", "Fitting · backup");
    const rows = [
      ["Run", "Date", "What the objective changed", "Seed → best (own objective)", "Outcome"],
      ["refit4", "≤ 09-22", "baseline", "opt 0.2909", "Sits inside force deadbands ~100× too wide; ktr ladder +0.34/seg"],
      ["refit5", "09-22", "force s = 2·s_shape; ktr_shape term", "0.7988 → 0.6885 (reload 0.6950)", "Force block improved; ktr ladder: 5 % of gap closed"],
      ["refit6", "09-23", "ktr_shape on 1× (was 2×) shape spread; SRX branches revived", "0.7843 → 0.7586 (re-score 0.7612)", "Stopped: ~4× slower evaluations; ladder flat"],
      ["refit7", "09-23", "ATP terms on raw spreads; rsT20; shape as guards", "0.4062 → 0.3355 (reload 0.3707)", "Noise discovered: accepted point was a winner’s curse"],
      ["refit8", "09-23", "ktrfit band; noise-aware floors + acceptance", "0.1665 → 0.1409", "Stalled inside 19.5 % noise; 0.356 re-scored after fix"],
      [{ text: "refit9", bold: true }, "09-24", "clean numerics; floors nEl·s²/3; guards 1.2×", { text: "0.2422 → 0.2079 (exact on reload)", bold: true }, "Incumbent; ladder flat, rsK80 climbs"],
    ];
    S.table(s, rows, { x: M, y: 1.2, w: W - 2 * M, h: 4.9 }, { fontSize: 14, colW: [1.0, 1.05, 3.6, 3.0, 3.683] });
    s.addShape("rect", { x: M, y: 6.25, w: W - 2 * M, h: 0.6, fill: { color: "FDECEA" }, line: { color: C.bad, width: 1 } });
    txt(s, S, "Each cost is on its own objective — not comparable across rows.", { x: M, y: 6.25, w: W - 2 * M, h: 0.6 },
      { fontSize: 16, bold: true, color: C.bad, align: "center", valign: "middle" });
    S.footer(s, "Sources: Analyses/KtrSegmentSlope/conclusions.md; dashboard/decisions.jsonl; results/checkin_refit6–9.json; StrainWindowNoise/conclusions.md (0.356)");
    s.addNotes(
      "refit4 to refit9 in three days (Sep 22-24). Each row changed the objective, so read each cost only against its own seed. " +
      "refit4_opt scored 0.2909 on its own objective and 0.7988 once the force deadbands were rebuilt from the shape spread (refit5 seed). " +
      "refit5: 6.3 h, force block improved, the ktr ladder closed 5 % of the gap. refit6: stopped; the revived SRX branches made each evaluation ~4x slower; the ladder stayed flat. " +
      "refit7: ATP-carrying terms scored on their raw spreads; the accepted 0.3355 re-scored 0.3707 from disk - this is where the noise investigation started. " +
      "refit8: noise-aware floors and acceptance; recorded 0.1409 but the gains were inside its own noise; under the fixed numerics refit8_opt scores 0.356 on its own objective, and 0.2419 on the refit9 objective (vs refit7's 0.2423) - no real progress. " +
      "refit9 on clean numerics: 0.2422 -> 0.2079 in one round, reproduced exactly on reload. The refit4 seed cost is not in the files I read.");
  }
}

module.exports = { addSlides, addBackup };
if (require.main === module) S.buildStandalone((pres, S) => { addSlides(pres, S); addBackup(pres, S); }, "part2.pptx");
