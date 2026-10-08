// PART 4: which second ATP lever (sensitivity on top of k2), and can [ATP] change hydrolysis?
// Inserted by part3.js right after the "ATP-sensitive levers" slide.
// Standalone: NODE_PATH=~/tools/pptx/node_modules node part4.js  -> part4.pptx
const fs = require("fs");
const S = require("./style");
const RES = "Analyses/LowATP_Refit9K2Sweep/results/";

function pngSize(rel) {
  try {
    const b = fs.readFileSync(S.asset(rel));
    if (b.toString("ascii", 1, 4) !== "PNG") return null;
    return [b.readUInt32BE(16), b.readUInt32BE(20)];
  } catch (e) { return null; }
}
function loadJson(rel) { try { return JSON.parse(fs.readFileSync(S.asset(rel), "utf8")); } catch (e) { return null; } }

function addSlides(pres, S) {
  const { C, M, W, TOP } = S;

  // ------------------------------------------------------------ second-lever sensitivity
  {
    const s = pres.addSlide();
    S.header(s, "A second ATP lever on top of k2: one at a time", "±10 % at 2 mM");
    const L = loadJson(RES + "levers_summary.json");
    S.image(s, RES + "levers_tornado.png", { x: M, y: TOP + 0.05, w: 8.4, h: 4.7 }, pngSize(RES + "levers_tornado.png"));
    const msgs = (L && Array.isArray(L.messages) && L.messages.length) ? L.messages
      : ["Sensitivity running: results land in results/levers_summary.json"];
    S.card(s, { x: M + 8.6, y: TOP + 0.05, w: W - 2 * M - 8.6, h: 4.7 }, "Reading", msgs, { fontSize: 13, bar: C.navy });
    S.bullets(s, [
      "Base: R2 ×0.54 at 2 mM (best k2-only). Each candidate ×0.9 / ×1.1 at 2 mM only (8 mM untouched by construction).",
      "Candidates in plausibility order: A ksrd2sr, ksr0 · AB kmsrd · B A0, tau_reg · C k_1, k1 (Pi) · X kah, ka, kd (controls).",
    ], { x: M, y: 6.0, w: W - 2 * M, h: 0.9 }, { fontSize: 13, space: 2 });
    S.footer(s, "Source: Analyses/LowATP_Refit9K2Sweep (RunAtpLeverSensitivity.m, results/levers.txt); tiers: Analyses/AtpEffectHypotheses/RunPlausibleSets.m");
    s.addNotes("Replicates the earlier ATP-hypothesis ladder in its cheapest form: a local sensitivity, not a fit. " +
      "Score = the same relevance-weighted absolute 2 mM chi-2 as the k2 sweep. Right panel: per key feature, blue = the lever at its better direction moves the model 2/8 ratio towards the data.");
  }

  // ------------------------------------------------------------ follow-up: SRX entry x k2
  {
    const s = pres.addSlide();
    S.header(s, "Following up the plausible lever: SRX entry (ksr0)", "k2 × ksr0 grid");
    S.image(s, RES + "ksr0grid.png", { x: M, y: TOP + 0.05, w: 8.9, h: 4.6 }, pngSize(RES + "ksr0grid.png"));
    S.card(s, { x: M + 9.1, y: TOP + 0.05, w: W - 2 * M - 9.1, h: 4.6 }, "Result", [
      "Best: R2 ×0.54 with ksr0 ×0.4, χ² 4.52 → 4.13 (−9 %)",
      "Fixes F–V over-sensitivity: FV2 ×0.60 → 0.74 (data 0.72); ktr protocol 0.68 → 0.72 (0.74); t0 1.54 → 1.35 (1.43)",
      "Slack ktr worse (0.75 → 0.80 vs 0.57); force, rsK80 and re-stretch timing unchanged",
      "Does not allow a stronger k2 gate: at g 0.44 still worse",
    ], { fontSize: 13, bar: C.good });
    S.bullets(s, [
      "Physiology: SRX entry from the ATP-bound detached state; fewer ATP-loaded heads at low ATP → less parking. Plausible, not proven (no 2 vs 8 mM SRX data).",
      "Verdict: a modest, chemically sensible partner for k2 — not the missing forward-limb lever that would decouple force from slack ktr.",
    ], { x: M, y: TOP + 4.8, w: W - 2 * M, h: 0.95 }, { fontSize: 13, space: 2 });
    S.footer(s, "Source: Analyses/LowATP_Refit9K2Sweep (RunKsr0Grid.m, results/ksr0grid.txt)");
    s.addNotes("3 x 4 grid, 12 evaluations, no fit. The ratios at the best point are diagnostic; the score is the absolute 2 mM relevance-weighted chi-2.");
  }

  addBestFit(pres, S, "k2ksr0");
  addBestFit(pres, S, "k2ksr0kmsrd");
  addFeaturePlots(pres, S);

  // ------------------------------------------------------------ can ATP change hydrolysis?
  {
    const s = pres.addSlide();
    S.header(s, "Can [ATP] change the hydrolysis rate (kah)?", "critical check");
    s.addShape("rect", { x: M, y: TOP + 0.05, w: W - 2 * M, h: 0.7, fill: { color: "FDECEA" }, line: { color: C.bad, width: 1 } });
    s.addText("No. Hydrolysis acts on ATP already in the pocket; free [ATP] enters only through binding — so the H1 gate (kah ×0.66 at 2 mM) is not physical.",
      { x: M + 0.15, y: TOP + 0.05, w: W - 2 * M - 0.3, h: 0.7, fontFace: S.F, fontSize: 16, bold: true, color: C.bad, valign: "middle" });
    S.card(s, { x: M, y: TOP + 0.95, w: 6.05, h: 3.25 }, "Chemistry and evidence", [
      "M·ATP ⇌ M·ADP·Pi is unimolecular on bound nucleotide; transient kinetics: k_obs rises with [ATP] only via binding, then plateaus at k+3 + k−3",
      "Saturated far below 2 mM (human β-S1 slow phase ATP-independent above ~50 µM)",
      "α-cardiac S1 hydrolysis ≈ 150–170 s⁻¹ (Deacon 2012) — our kah 176 s⁻¹ matches",
      "Fibres: Kawai 1993 scheme puts MgATP only in step 1; hydrolysis burst < 50 ms (Ferenczi 1984)",
    ], { fontSize: 13, bar: C.navy });
    S.card(s, { x: M + 6.25, y: TOP + 0.95, w: W - 2 * M - 6.25, h: 3.25 }, "Indirect routes do not rescue it", [
      "Ionic strength balanced; free Mg²⁺ slows ADP release, not hydrolysis",
      "kamh: reverse cleavage depends on BOUND products — gating kah and kamh differently shifts K3 without cause",
      "Hidden apo heads inside PT: at most a 2–13 % apparent effect at 2 mM, not ×0.66",
      "mant-ATP SRX chase measures bound-nucleotide turnover, not chase [ATP]",
    ], { fontSize: 13, bar: C.warn });
    S.card(s, { x: M, y: TOP + 4.35, w: W - 2 * M, h: 1.3 }, "What the H1 gate was standing in for", [
      "Slower re-priming = slower attachment: the forward-limb partner k2 lacks (ktr down more than F0 up). Real candidates: strain-dependent ADP release / ATP binding in R2; nucleotide/strain-dependent thick-filament OFF (SRX) recruitment.",
    ], { fontSize: 13, bar: C.good });
    S.footer(s, "Deacon 2012 CMLS 10.1007/s00018-012-0927-3; Bagshaw & Trentham 1974 10.1042/bj1410331; Ferenczi 1984 10.1113/jphysiol.1984.sp015311; Kawai 1993 10.1161/01.res.73.1.35; Beard 2022 10.1016/j.bpj.2022.07.029");
    s.addNotes("Full write-up: Analyses/LowATP_Refit9K2Sweep/kah_atp_dependence.md. Our own AtpEffectHypotheses note calling kah 'unambiguous ATP chemistry' is wrong; the RunPlausibleSets header ('NOT PLAUSIBLE') is right. " +
      "Could not verify: Baker/Beard solution temperature and free Mg2+; any 2 vs 8 mM SRX-fraction or X-ray measurement; Deacon marks the beta hydrolysis assignment provisional.");
  }
}


// Steady force with its spread -- called by part3.js right after the feature-by-feature slide.
function addSteady(pres, S) {
  const { C, M, W, TOP } = S;
  const s = pres.addSlide();
  S.header(s, "Steady force: the model is inside the data spread", "8 and 2 mM");
  S.image(s, RES + "steady_state_sd.png", { x: M, y: TOP + 0.05, w: 8.9, h: 4.4 }, pngSize(RES + "steady_state_sd.png"));
  S.card(s, { x: M + 9.1, y: TOP + 0.05, w: W - 2 * M - 9.1, h: 4.4 }, "We do fit the pooled level", [
    "Late plateau, pooled data vs model: 8 mM 69 vs 71.5 kPa (+3 %), 2 mM 90 vs 88 kPa (−2.5 %)",
    "The apparent gap came from plotting prep 03/27 raw: its isometric (79.8 kPa) is 15 % above the pool (69.3 kPa); trace slides now show the pooled band",
    "2/8 level: data ×1.32 ± 0.09, model ×1.24",
  ], { fontSize: 13, bar: C.navy });
  S.card(s, { x: M, y: TOP + 4.6, w: W - 2 * M, h: 1.1 }, "Case for a higher steady weight", [
    "Raw between-prep SD of steady/own isometric is 4 % (8 mM) and 8 % (2 mM); shrinkage inflates it to 17 % / 20 %, so its weight is ~17× below its own precision. Caveat: steady and the anchor are the same-length isometric, so part of that 4 % is self-normalisation.",
  ], { fontSize: 13, bar: C.warn });
  S.footer(s, "Source: Analyses/LowATP_Refit9K2Sweep/PlotSteadyStateSD.m (buildPooledTarget; model: k2sweep.mat, 2 mM = R2 ×0.54)");
  s.addNotes("Bands: dark = raw between-prep SD, light = moderated SD used by the objective (fitFDist prior, d0 capped at 5). " +
    "Model minus data in raw SDs: 8 mM +0.9 to +1.2 (segments 1-4), 2 mM -0.2 to -0.35. Prep 04/03 at 2 mM (squares) sits under the model markers; it has no cycle-5 value. " +
    "Options: exempt the level terms from shrinkage, or pin the level like FV (a fixed share), as proposed in DataPooling conclusions.");
}


// Full-slide best-fit figure (slack + ktr protocol + F-V, 8 and 2 mM) for one parameter set.
const BESTFIT = {
  k2only: { title: "Best physiological fit (A): A2 detachment only", kicker: "8 + 2 mM",
    set: ["8 mM: incumbent (refit9 parameters, robust reader)", "2 mM: R2 ×0.54 (pAtp_k2 0.44); nothing else changes"],
    plaus: ["ATP binding is part of A2 detachment: the canonical ATP step (tier A)", "Hydrolysis not gated"],
    score: ["2 mM χ² 4.52 relevance-weighted / 7.86 plain", "8 mM cost unchanged (0.424)"] },
  k2ksr0: { title: "Best extended fit (B): A2 detachment + SRX entry", kicker: "8 + 2 mM",
    set: ["8 mM: incumbent (unchanged)", "2 mM: R2 ×0.54 and ksr0 ×0.4"],
    plaus: ["ksr0: fewer ATP-loaded heads parked in SRX; plausible, unproven (tier A)", "No 2 vs 8 mM SRX data to test it"],
    score: ["2 mM χ² 4.13 / 6.73 (set A 4.52 / 7.86)", "FV2 and ktr protocol now in range; slack ktr slightly worse"] },
  k2ksr0kmsrd: { title: "Best fit found (C): + slower SRX-ADP exit", kicker: "8 + 2 mM",
    set: ["8 mM: incumbent (unchanged)", "2 mM: R2 ×0.54, ksr0 ×0.4, kmsrd ×0.4"],
    plaus: ["kmsrd (SRX-ADP exit): weak, tier AB; ATP-binding reading rejected on its Km", "A hypothesis, not a result"],
    score: ["2 mM χ² 3.33 / 6.13 (B 4.13 / 6.73)", "Slack ktr ×0.58 (data 0.57), force unchanged", "ktr protocol ×0.59 (0.74), t0 ×1.19 (1.43)"] },
};
function addBestFit(pres, S, tag) {
  const { C, M, W, TOP } = S; const B = BESTFIT[tag];
  const s = pres.addSlide();
  S.header(s, B.title, B.kicker);
  const f = RES + "bestfit_" + tag + ".png";
  S.image(s, f, { x: M, y: TOP - 0.05, w: 9.95, h: 5.75 }, pngSize(f));
  const x = M + 10.05, w = W - M - x;
  S.card(s, { x, y: TOP, w, h: 1.6 }, "Parameter set", B.set, { fontSize: 11, headSize: 14, bar: C.navy });
  S.card(s, { x, y: TOP + 1.7, w, h: 1.8 }, "Plausibility", B.plaus, { fontSize: 11, headSize: 14, bar: C.good });
  S.card(s, { x, y: TOP + 3.6, w, h: 2.15 }, "Score", B.score, { fontSize: 11, headSize: 14, bar: C.warn });
  S.footer(s, "Pooled data ± SD (n = 3; FV2 n = 2), dashed = model; slack in kPa (pooled level), ktr and F–V normalised. Analyses/LowATP_Refit9K2Sweep/RunBestFitFigures.m");
  s.addNotes("Top: five slack-restretch cycles, pooled 8 mM (navy) and 2 mM (red) bands vs model. Bottom: ktr protocol (pooled over 3 days), force-velocity (Baker 8 mM scored; Baker 2 mM open circles, not scored; FV2 in-prep diamonds, model FV2 crosses), power. " +
    "Visible misses: the ktr-protocol re-stretch overshoot (model 1.35-1.4 vs data ~0.75), slow post-restretch recovery at 2 mM, short-length redeveloped force in cycles 1-2 at 2 mM, and power at 2-4 ML/s below Baker 8 mM.");
}

// Feature tables: 8 mM data/model and 2 mM data vs sets A/B/C (results/feature_table.json).
function addFeatureTables(pres, S) {
  const { C, M, W, TOP } = S;
  const T = loadJson(RES + "feature_table.json");
  const num = (x) => (typeof x !== "number" || !isFinite(x)) ? "–"
    : (Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : Math.abs(x) >= 1 ? x.toFixed(2) : x.toFixed(3));
  const pct = (x) => (typeof x === "number" && isFinite(x)) ? " ± " + Math.round(100 * x) + "%" : "";
  const zfill = (z) => (typeof z !== "number" || !isFinite(z)) ? null : z < 1 ? "D8F0DA" : z < 2 ? "FCEBC7" : "F8D3CF";
  const groups = [
    { title: "Backup: feature table (1/2), slack", pick: (g) => g.startsWith("Slack") },
    { title: "Backup: feature table (2/2), ktr / F–V / read-outs", pick: (g) => !g.startsWith("Slack") },
  ];
  for (const G of groups) {
    const s = pres.addSlide();
    S.header(s, G.title, "A · B · C");
    if (!T) { S.placeholder(s, { x: M, y: TOP, w: W - 2 * M, h: 5.5 }, "results/feature_table.json"); continue; }
    const rows = [["Feature", "Unit", "8 mM data", "8 mM model", "2 mM data", "A", "B", "C", "Data 2/8"]];
    let lastGroup = "";
    for (const r of T.rows.filter((r) => G.pick(r.group))) {
      if (r.group !== lastGroup && !G.title.includes("slack")) {
        rows.push([{ text: r.group, bold: true, fill: "E3E8EF" }, { text: "", fill: "E3E8EF" }, { text: "", fill: "E3E8EF" }, { text: "", fill: "E3E8EF" },
                   { text: "", fill: "E3E8EF" }, { text: "", fill: "E3E8EF" }, { text: "", fill: "E3E8EF" }, { text: "", fill: "E3E8EF" }, { text: "", fill: "E3E8EF" }]);
        lastGroup = r.group;
      }
      const z = Array.isArray(r.z) ? r.z : [];
      const isRead = r.group.startsWith("Model read-outs");
      const cell = (v, zi) => { const f = zfill(z[zi]); return f ? { text: num(v), fill: f } : num(v); };
      const mx = Array.isArray(r.max) ? r.max : null;
      const readCell = (v, k) => (mx && typeof mx[k] === "number" && r.feature.startsWith("Attached")) ? num(v) + " (max " + num(mx[k]) + ")" : num(v);
      rows.push([
        r.feature, r.unit,
        isRead ? { text: "bounds " + (r.bounds || "–"), color: C.mute } : num(r.data8) + pct(r.sd8),
        isRead ? readCell(r.model8, 0) : cell(r.model8, 0),
        isRead ? "–" : num(r.data2) + pct(r.sd2),
        isRead ? readCell(r.A, 1) : cell(r.A, 1), isRead ? readCell(r.B, 2) : cell(r.B, 2), isRead ? readCell(r.C, 3) : cell(r.C, 3),
        (typeof r.ratio === "number" && isFinite(r.ratio)) ? "×" + r.ratio.toFixed(2) : "–",
      ]);
    }
    S.table(s, rows, { x: M, y: TOP - 0.05, w: W - 2 * M, h: 5.6 }, { fontSize: 11, colW: [3.35, 0.7, 1.45, 1.25, 1.45, 1.05, 1.05, 1.05, 0.98] });
    S.footer(s, "Model cells: green < 1 SD, amber 1–2 SD, red > 2 SD from the data (SD = the objective's). A: R2 ×0.54 · B: + ksr0 ×0.4 · C: + kmsrd ×0.4 (2 mM only). Source: BuildFeatureTable.m");
    s.addNotes("Data: pooled targets (slack n = 3, ktr protocol n = 3, FV2 n = 2); F-V rows are the borrowed Baker 2022 curve (8 mM scored, 2 mM not scored; no SD). " +
      "Read-outs have no data: the objective's physiological bounds are shown instead. Attached fraction: isometric steady state per cycle, mean (max over cycles). " +
      "Largest 2 mM misses in every set: valley after the re-stretch peak (model 55 vs data 82 kPa) and second peak - the post-peak collapse is twice too deep at low ATP (vall2_dy -31 vs -14.5) - plus recovery time and valley timing.");
  }
}

// Feature plots (PlotFeaturePanels.m): slack features over the five cycles; ktr protocol, F-V, read-outs.
function addFeaturePlots(pres, S) {
  const { C, M, W, TOP } = S;
  const P = [
    { f: "features_slack.png", title: "Features over the five slack cycles: 8 mM vs 2 mM", kicker: "A · B · C",
      notes: "Bands: pooled data ± the SD the objective uses (n = 3). Dashed: models. Force features, re-stretch peak and t0 follow the data; the valley and second peak collapse at 2 mM (post-peak drop -31 vs -14.5 kPa), recovery time stays flat (data +68 %), valley time moves the wrong way. Slack ktr: only C sits in the 2 mM band; rsK80 is 20-25 % high in all sets." },
    { f: "features_ktr_fv.png", title: "ktr protocol, force–velocity and model read-outs", kicker: "A · B · C",
      notes: "Top: ktr-protocol features and the in-prep FV2 point (data ± SD; 8 mM model next to the 8 mM data, sets A/B/C next to the 2 mM data). Middle: Baker F-V (8 mM scored, 2 mM not scored) and power. Bottom: model read-outs against the objective's bounds (grey); attached fraction also against the literature 0.12-0.15 (green). XTOR sits at the objective's 15 s-1 ceiling at 8 mM." },
  ];
  for (const p of P) {
    const s = pres.addSlide();
    S.header(s, p.title, p.kicker);
    S.image(s, RES + p.f, { x: M, y: TOP - 0.05, w: W - 2 * M, h: 5.8 }, pngSize(RES + p.f));
    S.footer(s, "A: R2 ×0.54 · B: + ksr0 ×0.4 · C: + kmsrd ×0.4 (2 mM only; 8 mM = incumbent). Source: Analyses/LowATP_Refit9K2Sweep/PlotFeaturePanels.m");
    s.addNotes(p.notes);
  }
}

module.exports = { addSlides, addSteady, addBestFit, addFeatureTables, addFeaturePlots };
if (require.main === module) S.buildStandalone(addSlides, "part4.pptx");
