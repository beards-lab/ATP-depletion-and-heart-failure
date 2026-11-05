// Twitch-ladder deck: fixed rat Ca -> fitted Ca (XB frozen) -> literature Ca + humanized XB -> T-World Ca.
// Mirrors the review artifact https://claude.ai/artifact/Y5S914KMgpy2tinwoGtjvB
// Build (WSL): NODE_PATH=~/tools/pptx/node_modules node build_deck_twitch.js
// Render:      powershell.exe -NoProfile -ExecutionPolicy Bypass -File render_v3.ps1 <deck.pptx> <outdir>
const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const S = require("../2026-09-status/style");
const { C, F, W, M, TOP } = S;
const HT = "Analyses/HumanTwitch/results";
const TW = "Analyses/TWorldCoupling/results";
const ab = (rel) => S.asset(rel);
function pngSize(abs) { const b = fs.readFileSync(abs); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
function img(s, rel, box) { S.image(s, ab(rel), box, pngSize(ab(rel))); }
function take(s, text, opts = {}) {
  const y = opts.y || 6.15, h = opts.h || 0.75;
  s.addShape("rect", { x: M, y, w: W - 2 * M, h, fill: { color: opts.fill || "FFF6D6" }, line: { color: opts.line || C.maize, width: 1 } });
  s.addText(text, { x: M + 0.2, y, w: W - 2 * M - 0.4, h, fontFace: F, fontSize: 15, color: C.ink, valign: "middle", margin: 0, isTextBox: true });
}
function caption(s, text, box) { s.addText(text, { ...box, fontFace: F, fontSize: 11, italic: true, color: C.mute, valign: "top", margin: 0, isTextBox: true }); }
const B = C.bad, G = C.good;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "The human EHT twitch: four ways to supply the Ca";
pres.author = "Filip Jezek";

// 1 Title
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addShape("rect", { x: 0, y: 4.55, w: W, h: 0.08, fill: { color: C.maize }, line: { color: C.maize } });
  s.addText("The human iPSC-EHT twitch: four ways to supply the Ca", { x: 0.8, y: 1.4, w: W - 1.6, h: 2.0, fontFace: F, fontSize: 40, bold: true, color: C.white, valign: "bottom", margin: 0, isTextBox: true });
  s.addText("fixed rat Ca · fitted Ca with a frozen mouse core · literature human Ca with a humanized core · the T-World myocyte", { x: 0.8, y: 3.45, w: W - 1.6, h: 0.9, fontFace: F, fontSize: 20, color: C.white, margin: 0, isTextBox: true });
  s.addText("Filip Jezek", { x: 0.8, y: 4.85, w: W - 1.6, h: 0.5, fontFace: F, fontSize: 22, color: C.maize, margin: 0, isTextBox: true });
  s.addText("October 2026", { x: 0.8, y: 5.4, w: W - 1.6, h: 0.45, fontFace: F, fontSize: 18, color: C.white, margin: 0, isTextBox: true });
}

// 2 The ladder (native chart, log axis)
{
  const s = pres.addSlide();
  S.header(s, "Fit quality at every stage", "overview");
  const labels = [
    "1  rat Rice 2008 Ca, nothing tuned",
    "2  Ca fitted, mouse XB frozen (fit C)",
    "3a de Lange Ca, literature human XB, no Q10",
    "3b de Lange Ca, humanized XB (isoform × Q10)",
    "3c de Lange Ca, all humanizing directions",
    "3d de Lange Ca, 24 XB parameters free",
    "4a T-World Ca as granted",
    "4b T-World Ca handling freed",
    "4c T-World closed loop (XB buffers TnC)",
  ];
  const vals = [39, 0.81, 10.8, 0.77, 0.38, 0.30, 0.97, 0.75, 6.7];
  // drawn log-scale bars (a native log axis crosses at 1, so sub-1 % bars would point left)
  const x0 = M + 4.3, bw = 3.6, lo = Math.log10(0.1), hi = Math.log10(100);
  const xs = (v) => x0 + bw * (Math.log10(v) - lo) / (hi - lo);
  const bad = [0, 2, 8], rowH = 0.52, y0 = TOP + 0.15;
  [0.1, 1, 10, 100].forEach((t) => {
    s.addShape("line", { x: xs(t), y: y0 - 0.05, w: 0, h: rowH * labels.length + 0.05, line: { color: "D1D5DB", width: 0.75 } });
    s.addText(String(t), { x: xs(t) - 0.4, y: y0 + rowH * labels.length + 0.02, w: 0.8, h: 0.3, fontFace: F, fontSize: 11, color: C.mute, align: "center", margin: 0, isTextBox: true });
  });
  s.addShape("line", { x: xs(0.35), y: y0 - 0.05, w: 0, h: rowH * labels.length + 0.05, line: { color: G, width: 1.25, dashType: "dash" } });
  s.addText("noise floor 0.35 %", { x: xs(0.35) - 1.7, y: y0 + rowH * labels.length + 0.32, w: 1.8, h: 0.28, fontFace: F, fontSize: 11, color: G, align: "right", margin: 0, isTextBox: true });
  s.addText("rms residual, % of data peak (log scale)", { x: xs(1) - 0.3, y: y0 + rowH * labels.length + 0.32, w: bw, h: 0.28, fontFace: F, fontSize: 12, color: C.mute, align: "center", margin: 0, isTextBox: true });
  labels.forEach((lab, k) => {
    const y = y0 + k * rowH, col = bad.includes(k) ? B : C.navy;
    s.addText(lab, { x: M, y, w: 4.15, h: rowH - 0.1, fontFace: F, fontSize: 13, color: C.ink, align: "right", valign: "middle", margin: 0, isTextBox: true });
    s.addShape("rect", { x: xs(0.1), y: y + 0.08, w: xs(vals[k]) - xs(0.1), h: rowH - 0.26, fill: { color: col }, line: { color: col } });
    s.addText((vals[k] >= 1 ? String(vals[k]) : vals[k].toFixed(2)) + " %", { x: xs(vals[k]) + 0.08, y, w: 0.9, h: rowH - 0.1, fontFace: F, fontSize: 13, bold: true, color: col, valign: "middle", margin: 0, isTextBox: true });
  });
  S.card(s, { x: 9.2, y: TOP + 0.05, w: W - M - 9.2, h: 5.6 }, "How to read it", [
    "Target: Plate 6 control EHT twitch, 15 wells × 60 beats, 1 Hz, 37 °C",
    "TTP 230, RT50 152, RT90 243 ms",
    "rms of the residual in % of peak; force scale σ solved in closed form",
    "Noise floor = template SEM, 0.35 %",
    { text: "Three failures (1, 3a, 4c) bracket the story; every other route fits to within 1 %", bold: true },
  ], { fontSize: 14, bar: C.maize });
  S.footer(s, "Analyses/CalciumTwitch/EHT_identification_report.md; Analyses/HumanTwitch/conclusions.md; Docs/reference/tworld-ca-source.md");
  s.addNotes("Every number is the rms residual of the model twitch against the mean EHT twitch, in percent of peak. The template SEM, 0.35 %, is the floor. Stages 1, 3a and 4c are the instructive failures.");
}

// 3 Stage 1
{
  const s = pres.addSlide();
  S.header(s, "1  A rodent Ca transient gives a rodent twitch", "fixed rat Ca");
  img(s, `${HT}/deck/F2_failure.png`, { x: M, y: TOP, w: W - 2 * M, h: 4.6 });
  caption(s, "Left: Rice 2008 rat transient (peak 48 ms). Right: EHT data (black, ±SD band) and the mouse model with that Ca (red, amplitude scaled).", { x: M, y: 5.75, w: W - 2 * M, h: 0.35 });
  take(s, "Tuned: nothing. Mouse core (ExampleHighATP), Marzban rat thin filament, only σ solved. Model peaks at 70 ms vs 230 ms: rms 39 %. The core is fast enough to follow any Ca (10 ms pulse → 56 % of F_ss), so shape must come from Ca and thin filament.", { fill: "FDECEA", line: B });
  S.footer(s, "Analyses/CalciumTwitch/conclusions.md; Analyses/HumanTwitch/results/deck/F2_failure.png");
}

// 4 Stage 2: fit C
{
  const s = pres.addSlide();
  S.header(s, "2  Fitting only the Ca: 0.81 %, XB untouched", "fitted Ca, frozen XB");
  img(s, `${HT}/deck/F3_fitC.png`, { x: M, y: TOP, w: 7.4, h: 2.85 });
  caption(s, "Left: the Ca the twitch demands (red, peak 108 ms) against the rat one (blue), normalised. Right: fit C vs data.", { x: M, y: 4.1, w: 7.4, h: 0.4 });
  S.table(s, [
    ["parameter", "original", "fit C", "status"],
    ["Ca rise τ₁", "20 ms", "112 ms", "fitted"],
    ["Ca decay τ₂", "110 ms", "108 ms", "fitted"],
    ["Ca onset", "+5 ms", "−2.7 ms", "fitted"],
    ["Ca dia / peak", "0.09 / 1.45 µM", "0.10 / 1.0 µM", "assumed"],
    ["tf_k_on (µM⁻¹s⁻¹)", "101", "4.15", "fitted ×0.041"],
    ["tf_k_off (s⁻¹)", "724", "3.8", "fitted ×0.0053"],
    ["tf_K_coop", "9.7", "20", { text: "upper bound", color: C.warn, bold: true }],
    ["SL0", "2.2 µm", "2.0 µm", "set"],
    ["all XB / SRX / passive / kSE", "–", "unchanged", { text: "frozen (audited)", color: G, bold: true }],
  ], { x: 8.1, y: TOP, w: W - M - 8.1, h: 4.3 }, { fontSize: 11, colW: [1.75, 1.05, 1.05, 0.88] });
  take(s, "Timing within 3 ms except RT90 (255 vs 243 ms); the tail beyond 550 ms is missed. Not unique: core ×0.5 or ×2, or k_off ×0.5/×2, refit to 1.3–2 % with a different transient. The twitch cannot split Ca shape from core speed.", { y: 5.75, h: 0.9 });
  S.footer(s, "Analyses/CalciumTwitch/EHT_identification_report.md §2–4; FitCEHTTwitch.m");
}

// 5 Stage 3: fixed literature Ca
{
  const s = pres.addSlide();
  S.header(s, "3  Human Ca fixed: the core runs at 0.45× mouse", "de Lange Ca, humanized XB");
  img(s, `${HT}/deck/F5_tier1.png`, { x: M, y: TOP, w: W - 2 * M, h: 3.15 });
  caption(s, "Left: de Lange 2021 (Fura-2, 3D hiPSC tissue, 37 °C, 1 Hz), fixed. Right: tier-1 fit, rms 0.77 %; TTP 235 / RT50 145 / RT90 252 ms.", { x: M, y: 4.38, w: W - 2 * M, h: 0.35 });
  S.card(s, { x: M, y: 4.8, w: 6.0, h: 2.05 }, "Fixed (not tuned)", [
    "Ca 0.15 → 0.60 µM, TTP 93, peak→50 % 186, peak→90 % 448 ms",
    "k1, k_1 (Pi step) at the literature ×0.3",
    "Untouched literature β-MHC core: 10.8 %, force peaks at 360 ms (rates from a 20–25 °C rig)",
  ], { fontSize: 13, bar: C.mute });
  S.card(s, { x: 6.83, y: 4.8, w: W - M - 6.83, h: 2.05 }, "Tuned (tier 1)", [
    "r_iso × r_T on k2, k_2, kah, kamh, ka, kd: 0.066 × 6.9 = net ×0.45 (×/÷ 1.12)",
    "Only the product is identified (r = −1.00): isoform ~0.1 × Q10 factor ~4.5",
    "thin filament k_on 9.2, k_off 5.1, K_coop 20; offset +1.7 ms",
  ], { fontSize: 13, bar: G });
  S.footer(s, "Analyses/HumanTwitch/conclusions.md §1–2; literature_ca_transient.md; tierConfigs.m (t0_v0, t1_isoT)");
}

// 6 Tier table
{
  const s = pres.addSlide();
  S.header(s, "3  What was tuned, tier by tier", "tier ladder");
  S.table(s, [
    ["tier", "free, besides thin filament and offset", "result", "rms"],
    ["t0_v0", "nothing in the core: literature ratios (six rates ×0.1, k1/k_1 ×0.3)", "far too slow, RT90 442 ms", { text: "10.8 %", color: B, bold: true }],
    [{ text: "t1_isoT", bold: true }, "isoform ratio r_iso (six rates), Q10 factor r_T (all eight cycle rates)", "net ×0.45 on six rates; Pi step not seen", { text: "0.77 %", bold: true }],
    ["t2a … t2e", "one extra direction each: Pi step, SRX exchange, kSE, titin k_pas, SL0", "each ≤ 0.08 pp; kSE ×0.78, SRX ×0.49, SL 1.91 µm", "0.69–0.77 %"],
    ["t2_allc", "all of the above inside the literature corridor (r_iso 0.05–0.15, r_T 3.7–8.4)", "every factor on a bound; SRX exchange ×0.30, 48 % resting SRX", "0.38 %"],
    ["t2_all, t2b_SRX", "same, corridor removed", { text: "escape: cycle ×6–7.5 faster than mouse, 17–24 ATP/head/beat", color: C.warn }, "0.31 %"],
    ["s3_free (3 starts)", "24 parameters: 8 cycle rates, SRX rates, stiffnesses, kSE, k_pas, A0, dr, SL0", "noise floor; two families: slow attachment (ka ×0.36, kSE ×0.44) or fast cycle / slow SRX", { text: "0.30 %", color: G, bold: true }],
  ], { x: M, y: TOP + 0.1, w: W - 2 * M, h: 4.4 }, { fontSize: 13, colW: [1.9, 4.6, 4.6, 1.233] });
  take(s, "Inside the humanizing corridor one number carries the fit: ADP release, hydrolysis, attachment and weak detachment together at ×0.45 of the mouse. Everything else buys at most 0.08 pp, or leaves physiology to get to the noise floor.", { y: 5.9, h: 0.9 });
  S.footer(s, "Analyses/HumanTwitch/conclusions.md §2–3; results/<tier>.log");
}

// 7 Overview tiers figure
{
  const s = pres.addSlide();
  S.header(s, "3  Seventeen fits, one shared misfit", "all tiers");
  img(s, `${HT}/overview_tiers.png`, { x: M, y: TOP, w: 9.6, h: 5.6 });
  S.card(s, { x: 10.35, y: TOP, w: W - M - 10.35, h: 5.6 }, "Reading", [
    "A: every tier except t0_v0 lies on the data",
    "B: the shared residual sits at 400–500 ms, relaxation from 30 % to 5 %",
    "C: resting SRX ranges 0.1–0.6",
    "D: ATP use spans 100× across fits of near-equal quality",
    { text: "The twitch constrains shape, not economy", bold: true },
  ], { fontSize: 13, bar: C.maize });
  S.footer(s, "Analyses/HumanTwitch/OverviewTiers.m → results/overview_tiers.png");
}

// 8 Ca assumption sensitivity
{
  const s = pres.addSlide();
  S.header(s, "3  Core speed follows the assumed Ca shape", "Ca sensitivity");
  img(s, `${HT}/deck/F7_consequence.png`, { x: M, y: TOP, w: W - 2 * M, h: 4.7 });
  take(s, "A transient 20 ms faster to peak demands a core ×0.21; 25 ms slower demands ×3.8 (faster than the mouse). Fit quality barely moves (0.77–0.86 %). Literature human iPSC Ca TTPs span 60–270 ms: a Ca recording in these wells is the first measurement to ask for.", { y: 6.0, h: 0.9 });
  S.footer(s, "Analyses/HumanTwitch/conclusions.md §1.1, §2.2 (t1_caFast, t1_caSlow, t1_caHiLvl)");
}

// 9 Muscle-level consequences
{
  const s = pres.addSlide();
  S.header(s, "3  The tier-1 core as skinned muscle at 37 °C", "predictions");
  img(s, `${HT}/deck/F4_mouse_vs_human.png`, { x: M, y: TOP, w: W - 2 * M, h: 4.7 });
  take(s, "Saturating Ca, SL 2.2 µm: F0 kept (F0 is a ratio), ktr 23 s⁻¹ at 37 °C (intact human ceiling ~18 s⁻¹), ATP use halved (8.4 vs 17 s⁻¹ per head), FV curvature toward β-MHC. These are predictions, not fits: case-B profiles show ktr and ATPase vary 4–17× across fits of equal quality.", { y: 6.0, h: 0.9 });
  S.footer(s, "Analyses/HumanTwitch/ProbeCores.m, Analyses/HumanizeXB/humanProbe.m; conclusions.md §3.1, §5");
}

// 10 Two cores
{
  const s = pres.addSlide();
  S.header(s, "3  Same twitch, two mechanisms", "non-uniqueness");
  img(s, `${HT}/deck/F9_two_cores.png`, { x: M, y: TOP, w: 8.9, h: 5.6 });
  S.card(s, { x: 9.65, y: TOP, w: W - M - 9.65, h: 2.6 }, "Humanized core", [
    "cycle ×0.45 of mouse",
    "rms 0.77 %",
    "1.8 ATP/head/beat",
  ], { fontSize: 14, bar: B });
  S.card(s, { x: 9.65, y: TOP + 3.0, w: W - M - 9.65, h: 2.6 }, "Fast core", [
    "cycle ×10, SRX exchange ×0.22",
    "rms 0.32 %",
    "24 ATP/head/beat, a 50 % SRX reservoir shapes relaxation",
  ], { fontSize: 14, bar: "6B4FA0" });
  S.footer(s, "Shared by every fit: thin filament k_off 4–5.5 s⁻¹, K_coop 16–20, Ca→force lag ~140 ms. Analyses/HumanTwitch/DeckFigureF9.m");
}

// 11 Section divider: T-World
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addText("4  T-World as the Ca source", { x: 0.8, y: 2.3, w: W - 1.6, h: 1.2, fontFace: F, fontSize: 40, bold: true, color: C.white, margin: 0, isTextBox: true });
  s.addText("Tomek's human ventricular myocyte (93 states, Land 2017 inside) supplies the Ca; the XB model turns it into the twitch.", { x: 0.8, y: 3.6, w: W - 1.6, h: 1.0, fontFace: F, fontSize: 20, color: C.white, margin: 0, isTextBox: true });
  s.addText("Open loop ('signal'): T-World's beat as a periodic Ca table. Closed loop ('xb'): the 93 states embedded, the XB thin filament buffers troponin Ca.", { x: 0.8, y: 4.8, w: W - 1.6, h: 0.9, fontFace: F, fontSize: 16, color: C.maize, margin: 0, isTextBox: true });
}

// 12 T-World time course
{
  const s = pres.addSlide();
  S.header(s, "4  One beat: AP → Ca → force", "T-World time course");
  img(s, `${TW}/tworld_timecourse.png`, { x: M, y: TOP, w: 8.9, h: 5.65 });
  S.card(s, { x: 9.65, y: TOP, w: W - M - 9.65, h: 5.65 }, "Electrical vs Ca vs force", [
    "AP upstroke 8 ms, Ca peak 62 ms, force peak 240 ms (data 230)",
    "Force peaks during the plateau; at APD90 (376 ms) force is back to ~50 %",
    "Relaxation runs mostly after repolarisation, as Ca decays",
    "T-World's own Land tension peaks at 175–194 ms: the adult model is faster than this tissue",
    { text: "EP not fitted, but freeing I_CaL ×2.76 lengthens APD90 270 → 367 ms (from upstroke; B marks it from the stimulus)", bold: true, color: C.warn },
  ], { fontSize: 13, bar: C.maize });
  S.footer(s, "Analyses/TWorldCoupling/TWorldTimecourse.m → results/tworld_timecourse.png, tworld_timecourse.csv");
}

// 13 Event timing table
{
  const s = pres.addSlide();
  S.header(s, "4  Event timing on the data's time axis", "T-World timing");
  S.table(s, [
    ["event (ms after stimulus)", "T-World as granted (tw1)", "Ca handling freed (tw2)", "EHT data"],
    ["AP upstroke (max dV/dt)", "8", "8", "–"],
    ["Ca peak", "61", "62", "–"],
    ["Land tension peak (T-World's own, unused)", "175", "194", "–"],
    ["Ca 50 % decay", "200", "235", "–"],
    [{ text: "XB force peak", bold: true }, { text: "240", bold: true }, { text: "240", bold: true }, { text: "230", bold: true }],
    ["APD50 end", "228", "326", "–"],
    ["APD90 end", "278", "376", "–"],
    ["force RT50 end", "380", "380", "382"],
    ["Ca 90 % decay", "411", "475", "–"],
    ["force RT90 end", "490", "485", "473"],
  ], { x: M, y: TOP + 0.1, w: 8.4, h: 5.0 }, { fontSize: 13, colW: [3.6, 1.8, 1.8, 1.2] });
  S.card(s, { x: 9.2, y: TOP + 0.1, w: W - M - 9.2, h: 5.0 }, "Lags", [
    "AP → Ca peak: 53 ms (both)",
    "Ca peak → force peak: 178 ms (de Lange Ca: 140 ms)",
    "T-World's faster Ca shifts more of the twitch shaping into the cross-bridges",
    "Force timing matches data within 12 ms in both fits",
  ], { fontSize: 14, bar: C.navy });
  S.footer(s, "Analyses/TWorldCoupling/results/tworld_timecourse.csv");
}

// 14 Fit with T-World Ca
{
  const s = pres.addSlide();
  S.header(s, "4  T-World Ca fits, but halves the core speed", "fitting with T-World");
  img(s, `${TW}/twitch_fit_tw1_tw2.png`, { x: M, y: TOP, w: 9.2, h: 5.65 });
  S.card(s, { x: 9.95, y: TOP, w: W - M - 9.95, h: 5.65 }, "Reading", [
    "As granted 0.97 %; Ca handling freed 0.75 % (de Lange 0.77 %)",
    "Force, attached heads and permissive fraction are nearly identical while the Ca differs (C)",
    "The twitch pins the thin-filament output, not the Ca",
    { text: "Net cycle factor ×0.22–0.24 with T-World vs ×0.45 with de Lange", bold: true },
    "±2 % ripple on the rise and at 480 ms remains in every tier",
  ], { fontSize: 13, bar: C.maize });
  S.footer(s, "Analyses/TWorldCoupling/CompareTwitchFits.m → results/twitch_fit_tw1_tw2.png; tiers tw1_isoT, tw2_caH");
}

// 15 What was tuned with T-World
{
  const s = pres.addSlide();
  S.header(s, "4  What was tuned in the XB model and in T-World", "parameters");
  const xb = "E8F0FB", tw = "FFF4CC";
  S.table(s, [
    ["", "de Lange Ca (t1_isoT)", "T-World as granted (tw1)", "T-World Ca freed (tw2)"],
    [{ text: "rms residual, % of data peak force (429 µN)", bold: true }, "0.77 (3.3 µN)", "0.97 (4.1 µN)", { text: "0.75 (3.2 µN)", bold: true }],
    [{ text: "XB  r_iso × r_T (net on k2, k_2, kah, kamh, ka, kd)", fill: xb }, { text: "0.066 × 6.9 = 0.45", fill: xb }, { text: "0.127 × 1.9 = 0.24", fill: xb }, { text: "0.154 × 1.4 = 0.22", fill: xb }],
    [{ text: "XB  k1, k_1 (Pi step)", fill: xb }, { text: "×0.3 × r_T (fixed ratio)", fill: xb }, { text: "same", fill: xb }, { text: "same", fill: xb }],
    [{ text: "thin filament k_on / k_off / K_coop", fill: xb }, { text: "9.2 / 5.1 / 20", fill: xb }, { text: "13.3 / 4.5 / 19.9", fill: xb }, { text: "12.4 / 5.9 / 20", fill: xb }],
    [{ text: "T-World  SERCA J_up multiplier", fill: tw }, { text: "–", fill: tw }, { text: "1", fill: tw }, { text: "0.81", fill: tw, bold: true }],
    [{ text: "T-World  I_CaL multiplier", fill: tw }, { text: "–", fill: tw }, { text: "1", fill: tw }, { text: "2.76", fill: tw, bold: true }],
    [{ text: "T-World  NCX I_NaCa multiplier", fill: tw }, { text: "–", fill: tw }, { text: "1", fill: tw }, { text: "2.01", fill: tw, bold: true }],
    [{ text: "T-World  all other currents", fill: tw }, { text: "–", fill: tw }, { text: "untouched (adult endo, 1 Hz)", fill: tw }, { text: "untouched; pre-paced 20 beats", fill: tw }],
    ["Ca dia → peak, TTP", "0.15 → 0.60 µM, 95 ms", "0.097 → 0.43 µM, 60 ms", "0.093 → 0.55 µM, 63 ms"],
  ], { x: M, y: TOP + 0.1, w: W - 2 * M, h: 4.4 }, { fontSize: 13, colW: [4.3, 2.65, 2.7, 2.683] });
  take(s, "The Ca levers moved in a readable direction: more entry (I_CaL), more extrusion (NCX), slightly less uptake (SERCA): a larger, longer transient than the adult default. Blue rows: XB model; yellow: T-World. Offset +1.7 / +7.3 / +7.6 ms.", { y: 6.05, h: 0.8 });
  S.footer(s, "Docs/reference/tworld-ca-source.md §7; Analyses/HumanTwitch/tierConfigs.m (tw1_isoT, tw2_caH)");
}

// 16a Open vs closed loop: feedback diagram (native shapes)
{
  const s = pres.addSlide();
  S.header(s, "4  Open vs closed loop: who buffers TnC Ca", "feedback");
  const TWc = "FFF4CC", XBc = "E8F0FB";
  const box = (x, y, w, h, t, o = {}) => {
    s.addShape("rect", { x, y, w, h, fill: { color: o.fill || C.white }, line: { color: o.line || C.navy, width: o.lw || 1.25, dashType: o.dash || "solid" } });
    s.addText(t, { x, y, w, h, fontFace: F, fontSize: o.fs || 12, bold: !!o.bold, color: o.color || C.ink, align: "center", valign: "middle", margin: 0.03, isTextBox: true });
  };
  const arr = (x1, y1, x2, y2, o = {}) => {
    const x = Math.min(x1, x2), y = Math.min(y1, y2);
    s.addShape("line", { x, y, w: Math.abs(x2 - x1) || 0.0001, h: Math.abs(y2 - y1) || 0.0001, flipH: x2 < x1, flipV: y2 < y1,
      line: { color: o.color || C.ink, width: o.width || 1.75, endArrowType: o.noEnd ? undefined : "triangle", beginArrowType: o.both ? "triangle" : undefined, dashType: o.dash || "solid" } });
  };
  const lab = (t, x, y, w, h, o = {}) => s.addText(t, { x, y, w, h, fontFace: F, fontSize: o.fs || 11, bold: !!o.bold, italic: !!o.italic, color: o.color || C.ink, align: o.align || "center", valign: "middle", margin: 0, isTextBox: true });

  const panel = (px, closed) => {
    const tx = px, tw = 2.75, xx = px + 3.2, xw = 2.55, y0 = 1.62, y1 = 4.42;
    lab(closed ? "Closed loop ('xb'): stalls at 6.7 %" : "Open loop ('signal', 'land'): fits 0.75–0.97 %", px, 1.17, 5.75, 0.38, { fs: 15, bold: true, color: closed ? B : G, align: "left" });
    s.addShape("rect", { x: tx, y: y0, w: tw, h: y1 - y0, fill: { color: TWc }, line: { color: TWc } });
    s.addShape("rect", { x: xx, y: y0, w: xw, h: y1 - y0, fill: { color: XBc }, line: { color: XBc } });
    lab("T-World myocyte", tx, y0 + 0.02, tw, 0.25, { fs: 10, bold: true, color: C.mute });
    lab("XB model", xx, y0 + 0.02, xw, 0.25, { fs: 10, bold: true, color: C.mute });
    const bx = tx + 0.2, bw = tw - 0.4, bh = 0.45;
    box(bx, 1.92, bw, bh, "AP: I_CaL trigger");
    box(bx, 2.57, bw, bh, "SR Ca release (RyR)");
    box(bx, 3.22, bw, bh, "free Ca_i", { bold: true, fill: "E3F1E9", line: G });
    arr(bx + bw / 2, 1.92 + bh, bx + bw / 2, 2.57);
    arr(bx + bw / 2, 2.57 + bh, bx + bw / 2, 3.22);
    const fx = xx + 0.2, fw = xw - 0.4;
    box(fx, 3.22, fw, bh, "thin filament N", { bold: true });
    box(fx, 3.87, fw, bh, "cross-bridges → force");
    arr(fx + fw / 2, 3.22 + bh, fx + fw / 2, 3.87);
    arr(bx + bw, 3.445, fx, 3.445, { color: C.navy });
    lab(closed ? "drives" : "reads Ca only", bx + bw - 0.35, 3.03, fx - bx - bw + 0.7, 0.18, { fs: 10, italic: true, color: C.navy });
    if (!closed) {
      box(bx, 3.87, bw, bh, "Land TnC buffer, 70 µM", { fs: 11 });
      arr(bx + bw / 2, 3.22 + bh, bx + bw / 2, 3.87, { both: true, color: G });
      box(fx, 1.92, fw, 0.95, "N is fitted to make the right force from T-World's Ca; nothing goes back", { fs: 11, fill: XBc, line: XBc, color: C.mute });
    } else {
      box(bx, 3.87, bw, bh, "TnC flux = 70 µM · d(1−N)/dt", { fs: 11, line: B, lw: 1.75 });
      arr(bx + bw / 2, 3.22 + bh, bx + bw / 2, 3.87, { both: true, color: B });
      arr(fx, 3.445 + 0.12, bx + bw, 4.095, { color: B, width: 2.25 });
      lab("feedback", fx + 0.05, 3.68, 1.0, 0.18, { align: "left", fs: 10, bold: true, color: B });
      box(fx, 1.92, fw, 0.95, "Land's buffer switched off: N is both the Ca sensor and the Ca store", { fs: 11, fill: XBc, line: XBc, color: C.mute });
    }
  };
  panel(M, false);
  panel(7.08, true);

  // the vicious circle
  lab("Why the closed loop has no solution", M, 4.55, 6, 0.3, { fs: 13, bold: true, color: C.navy, align: "left" });
  const cw = 2.78, cg = 0.4, cy = 4.92, ch = 0.85;
  const cyc = [
    "the twitch needs a mostly-off, slow filament (1−N ≈ 0.2, k_on ≈ 9 µM⁻¹s⁻¹)",
    "a mostly-off filament holds little Ca: weak buffering on the upstroke",
    "free Ca rises higher and faster: 0.42 → 0.58 µM, TTP 54 → 34 ms",
    "to undo that, the filament must switch even slower",
  ];
  cyc.forEach((t, i) => {
    const x = M + i * (cw + cg);
    box(x, cy, cw, ch, t, { fs: 11.5, line: B, fill: "FDECEA" });
    if (i < 3) arr(x + cw, cy + ch / 2, x + cw + cg, cy + ch / 2, { color: B });
  });
  const xl = M + 3 * (cw + cg) + cw / 2, xs = M + (cw + cg) + cw / 2, yr = cy + ch + 0.22;
  arr(xl, cy + ch, xl, yr, { color: B, noEnd: true });
  arr(xl, yr, xs, yr, { color: B, noEnd: true });
  arr(xs, yr, xs, cy + ch, { color: B });
  take(s, "Fix: a separate troponin–Ca occupancy state that buffers fast (like Land's Ca_TRPN) and drives the slow, cooperative switch N. Until it exists, fit in open loop.", { y: 6.2, h: 0.7 });
  S.footer(s, "Model/tworldCaLayer.m (coupling 'signal' | 'land' | 'xb'); Analyses/TWorldCoupling/conclusions.md §4; tier tw1_isoT_xb");
  s.addNotes("Troponin C is a large Ca buffer: 70 µM of sites holding about 25 µM Ca at the peak, so contraction feeds back on the free Ca. Open loop: T-World buffers with Land's own troponin and our thin filament only reads the Ca. Closed loop: our filament's permissive fraction becomes the buffer. The fit needs a mostly-off, slow filament, which buffers little, so Ca rises higher and faster, which needs a slower filament still. N is a regulatory switch; the Ca-bound troponin fraction need not equal it.");
}

// 16 Closed loop failure
{
  const s = pres.addSlide();
  S.header(s, "4  Closed loop fails: N cannot buffer TnC Ca", "closed loop 'xb'");
  img(s, `${TW}/tworld_coupling_human.png`, { x: M, y: TOP, w: 8.2, h: 5.6 });
  S.card(s, { x: 8.95, y: TOP, w: W - M - 8.95, h: 5.6 }, "Why it stalls at 6.7 %", [
    "Open loop ('signal', 'land') reproduces the standalone beat to 5×10⁻⁴ µM",
    "In 'xb' the XB filament replaces Land's 70 µM TnC buffer",
    "The twitch needs a mostly-off, slow filament; a mostly-off filament stops buffering",
    "Ca grows and sharpens: 0.42 → 0.58 µM, TTP 54 → 34 ms",
    { text: "Fix: a separate TnC–Ca occupancy state before the loop can close", bold: true, color: B },
  ], { fontSize: 13, bar: B });
  S.footer(s, "Analyses/TWorldCoupling/RunTWorldCoupling.m, conclusions.md §4; tier tw1_isoT_xb");
}

// 17 Summary
{
  const s = pres.addSlide();
  S.header(s, "Summary: what the twitch fixes and leaves open", "summary");
  S.table(s, [
    ["stage", "Ca source", "tuned in XB", "tuned in Ca layer", "rms"],
    ["1", "Rice 2008 rat, fixed", "nothing", "nothing", { text: "39 %", color: B, bold: true }],
    ["2 (fit C)", "bi-exponential, fitted", "nothing", "τ rise, τ decay, onset; k_on, k_off, K_coop", "0.81 %"],
    ["3 (tier 1)", "de Lange 2021, fixed", "r_iso, r_T (net ×0.45)", "k_on, k_off, K_coop, offset", "0.77 %"],
    ["3 (free)", "de Lange 2021, fixed", "24 parameters", "k_on, k_off, K_coop, offset", { text: "0.30 %", color: G, bold: true }],
    ["4 (tw1)", "T-World adult, as granted", "r_iso, r_T (net ×0.24)", "k_on, k_off, K_coop, offset", "0.97 %"],
    ["4 (tw2)", "T-World, Ca handling freed", "r_iso, r_T (net ×0.22)", "+ SERCA ×0.81, I_CaL ×2.76, NCX ×2.01", "0.75 %"],
  ], { x: M, y: TOP + 0.1, w: W - 2 * M, h: 3.3 }, { fontSize: 13, colW: [1.4, 2.9, 2.6, 4.333, 1.1] });
  S.card(s, { x: M, y: 4.75, w: 3.95, h: 2.1 }, "Robust in every fit", [
    "thin filament: k_off 4–6 s⁻¹, K_coop ~20",
    "activation 73–85 %",
    "Ca ×1.5 → +11–14 % force; ka ×0.5 → −49–55 %",
  ], { fontSize: 13, bar: G });
  S.card(s, { x: 4.69, y: 4.75, w: 3.95, h: 2.1 }, "Not identified", [
    "core speed: ×0.22 (T-World) to ×0.45 (de Lange), ×0.21–3.8 across literature Ca",
    "ktr (4×), ATP per beat (16×)",
  ], { fontSize: 13, bar: B });
  S.card(s, { x: 8.88, y: 4.75, w: W - M - 8.88, h: 2.1 }, "What would decide it", [
    "relative Ca recording in these wells",
    "one muscle-level number: ATPase / O₂, or ktr",
    "an APD measurement (checks I_CaL ×2.76)",
  ], { fontSize: 13, bar: C.navy });
  S.footer(s, "Worktree C:\\home\\git\\ATP-wt-human (branch HumanTwitch)");
}

pres.writeFile({ fileName: path.resolve(__dirname, "..", "2026-10-HumanTwitch_twitch_v2.pptx") }).then((f) => console.log("wrote " + f));
