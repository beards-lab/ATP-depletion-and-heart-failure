// Storyline deck v3. Build: NODE_PATH=~/tools/pptx/node_modules node build_deck_v3.js
const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const S = require("../2026-09-status/style");
const { C, F, W, M, TOP } = S;
const DECK = path.resolve(__dirname, "..", "..", "..", "Analyses", "HumanTwitch", "results", "deck");
const fig = (n) => path.join(DECK, n);
function pngSize(abs) { if (!fs.existsSync(abs)) return null; const b = fs.readFileSync(abs); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
function img(s, name, box) { S.image(s, fig(name), box, pngSize(fig(name))); }
function msg(s, text, color) { s.addText(text, { x: M, y: TOP - 0.05, w: W - 2 * M, h: 0.55, fontFace: F, fontSize: 17, bold: true, color: color || C.navy, valign: "middle", margin: 0 }); }
function take(s, text, fill, line, y) {
  const yy = y || 6.1;
  s.addShape("rect", { x: M, y: yy, w: W - 2 * M, h: 0.8, fill: { color: fill || "FFF6D6" }, line: { color: line || C.maize, width: 1 } });
  s.addText(text, { x: M + 0.2, y: yy, w: W - 2 * M - 0.4, h: 0.8, fontFace: F, fontSize: 15, color: C.ink, valign: "middle", margin: 0 });
}
function arrow(s, x1, y1, x2, y2, opts = {}) {
  const x = Math.min(x1, x2), y = Math.min(y1, y2);
  const w = Math.abs(x2 - x1) || 0.0001, h = Math.abs(y2 - y1) || 0.0001;
  s.addShape("line", { x, y, w, h, flipH: x2 < x1, flipV: y2 < y1, line: { color: opts.color || C.ink, width: opts.width || 1.75, endArrowType: "triangle", dashType: opts.dash || "solid" } });
}
function label(s, text, x, y, w, h, opts = {}) {
  s.addText(text, { x, y, w, h, fontFace: F, fontSize: opts.fontSize || 11, color: opts.color || C.ink, align: opts.align || "center", valign: "middle", margin: 0, italic: !!opts.italic, bold: !!opts.bold });
}
function stateBox(s, x, y, w, h, name, sub, fill) {
  s.addShape("roundRect", { x, y, w, h, rectRadius: 0.08, fill: { color: fill }, line: { color: C.navy, width: 1.25 } });
  s.addText([{ text: name, options: { bold: true, fontSize: 16, breakLine: true } }, { text: sub, options: { fontSize: 10 } }], { x, y, w, h, fontFace: F, color: C.navy, align: "center", valign: "middle", margin: 0 });
}
const G = C.good, Wn = C.warn, B = C.bad;
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "From the mouse cross-bridge model to the human iPSC twitch: what we found and what we need";

// 1 Title
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addShape("rect", { x: 0, y: 4.55, w: W, h: 0.08, fill: { color: C.maize }, line: { color: C.maize } });
  s.addText("From the mouse cross-bridge model to the human iPSC twitch", { x: 0.8, y: 1.6, w: W - 1.6, h: 1.8, fontFace: F, fontSize: 40, bold: true, color: C.white, valign: "bottom", margin: 0 });
  s.addText("what we found, and what we need from the experiments", { x: 0.8, y: 3.45, w: W - 1.6, h: 0.8, fontFace: F, fontSize: 26, color: C.white, margin: 0 });
  s.addText("Filip Jezek", { x: 0.8, y: 4.85, w: W - 1.6, h: 0.5, fontFace: F, fontSize: 22, color: C.maize, margin: 0 });
  s.addText("October 2026", { x: 0.8, y: 5.4, w: W - 1.6, h: 0.45, fontFace: F, fontSize: 18, color: C.white, margin: 0 });
}

// 2 The model
{
  const s = pres.addSlide();
  S.header(s, "The model: DRX and SRX heads, gated by Ca", "model");
  const bw = 1.5, bh = 0.8, xL = 1.0, xR = 4.6, y1 = 1.7, y2 = 3.35, y3 = 5.0;
  stateBox(s, xL, y1, bw, bh, "SR", "SRX · M·ATP", "E6E1F5");
  stateBox(s, xR, y1, bw, bh, "SRD", "SRX · M·ADP·Pi", "E6E1F5");
  stateBox(s, xL, y2, bw, bh, "PT", "DRX · M·ATP", C.light);
  stateBox(s, xR, y2, bw, bh, "PD", "DRX · M·ADP·Pi, ready", C.light);
  stateBox(s, xR, y3, bw, bh, "A1", "attached, pre-stroke", "FFF1C2");
  stateBox(s, xL, y3, bw, bh, "A2", "attached, post-stroke", "FFE08A");
  const cxL = xL + bw / 2, cxR = xR + bw / 2, d = 0.12;
  arrow(s, xL + bw, y2 + bh / 2 - d, xR, y2 + bh / 2 - d);  arrow(s, xR, y2 + bh / 2 + d, xL + bw, y2 + bh / 2 + d);
  label(s, "ATP hydrolysis  kah / kamh", xL + bw, y2 - 0.3, xR - xL - bw, 0.28);
  arrow(s, cxR - d, y2 + bh, cxR - d, y3);  arrow(s, cxR + d, y3, cxR + d, y2 + bh);
  label(s, "attach  ka × gate(Ca)", cxR + 0.2, y2 + bh + 0.3, 1.8, 0.3, { align: "left", bold: true, color: B });
  label(s, "detach  kd", cxR - 1.3, y2 + bh + 0.3, 1.0, 0.3, { align: "right" });
  arrow(s, xR, y3 + bh / 2 - d, xL + bw, y3 + bh / 2 - d);  arrow(s, xL + bw, y3 + bh / 2 + d, xR, y3 + bh / 2 + d);
  label(s, "power stroke, Pi release  k1 / k_1", xL + bw, y3 + bh + 0.02, xR - xL - bw, 0.28);
  arrow(s, cxL, y3, cxL, y2 + bh, { color: B, width: 2.25 });
  label(s, "ADP release + ATP binding  k2", cxL - 1.0, y2 + bh + 0.3, 1.0, 0.5, { align: "right", bold: true, color: B, fontSize: 10 });
  arrow(s, cxL - d, y2, cxL - d, y1 + bh);  label(s, "park", cxL - 0.75, y1 + bh + 0.3, 0.5, 0.3, { align: "right", fontSize: 10 });
  arrow(s, xL + bw, y1 + bh / 2 - d, xR, y1 + bh / 2 - d);  arrow(s, xR, y1 + bh / 2 + d, xL + bw, y1 + bh / 2 + d);
  label(s, "SRX hydrolysis", xL + bw, y1 - 0.3, xR - xL - bw, 0.28);
  arrow(s, cxR - d, y1 + bh, cxR - d, y2);  arrow(s, cxR + d, y2, cxR + d, y1 + bh);
  label(s, "SRX exit (force-sensitive)", cxR + 0.2, y1 + bh + 0.3, 2.0, 0.3, { align: "left", fontSize: 10 });
  s.addShape("rect", { x: 0.6, y: y1 - 0.5, w: 5.9, h: 0.3, fill: { color: "FFFFFF" }, line: { color: "FFFFFF" } });
  label(s, "SRX: parked heads, out of the cycle (slow turnover)", 0.6, y1 - 0.5, 5.9, 0.3, { bold: true, color: "6B4FA0", fontSize: 12 });
  label(s, "DRX: heads in the cycle (detached or attached)", xL + bw + 0.1, y2 - 0.62, xR - xL - bw - 0.2, 0.3, { bold: true, color: C.navy, fontSize: 11 });
  s.addText("Attached heads carry a strain distribution; force = strain × stiffness. The whole element sits in series with an elastic element (kSE) and parallel to titin.",
    { x: 0.6, y: 6.15, w: 6.2, h: 0.6, fontFace: F, fontSize: 11, italic: true, color: C.mute, margin: 0 });
  S.card(s, { x: 7.3, y: TOP, w: W - M - 7.3, h: 2.6 }, "The three layers", [
    "cross-bridge core (above): rates, stiffness, SRX; identified on mouse trabeculae at saturating Ca",
    "thin filament: Ca binds, sites switch on cooperatively, gate × ka (Campbell/Tewari form); never calibrated on data",
    "Ca transient: prescribed input (rat literature, fitted, or human literature)",
  ], { fontSize: 13, bar: C.navy });
  S.card(s, { x: 7.3, y: TOP + 2.8, w: W - M - 7.3, h: 2.6 }, "What a twitch tests", [
    "activation: how many sites open, how fast (thin filament + Ca)",
    "cycling: how fast heads attach, stroke and let go (the core rates)",
    "the SRX pool: how many heads are parked and how fast they return",
    "one amplitude scale σ (µN of tissue force per model kPa)",
  ], { fontSize: 13, bar: C.maize });
  S.footer(s, "Model/dPUdT_CombinedTransitions.m (Beard et al. 2022 lineage); Model/calciumActivation.m");
}

// 3 Mouse model reproduces mouse data
{
  const s = pres.addSlide();
  S.header(s, "1 · The mouse model reproduces the mouse data", "where we start");
  msg(s, "Pooled 8 mM (blue) and 2 mM ATP (red) data from 3 preparations; model dashed. Slack–restretch, ktr protocol, normalised force–velocity and power.");
  img(s, "bestfit_k2ksr0.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "This is the core we take to the human tissue: the same six rates, stiffnesses and SRX kinetics, identified at saturating Ca.");
  S.footer(s, "Analyses/LowATP_Refit9K2Sweep (best high + low ATP fit: k2 ×0.54 and SRX entry ×0.4 at 2 mM); pooled target: Analyses/DataPooling.");
}

// 4 EHT data and well selection
{
  const s = pres.addSlide();
  S.header(s, "2 · The EHT twitch and the wells we used", "data");
  img(s, "eht_controls.png", { x: M, y: TOP, w: 8.6, h: 3.0 });
  S.card(s, { x: M + 8.8, y: TOP, w: W - 2 * M - 8.8, h: 3.0 }, "Plate 6, 1 Hz, 24 wells", [
    "Mantarray export: active force (baseline-subtracted), 100 Hz, 64 s",
    "no platemap in the file",
    "15 wells form one tight group (blue); 9 wells (columns 4–6, rows A–C) relax and beat faster in a graded way: the isoproterenol signature (orange)",
    "we took the 15 as controls; 60 beats each, averaged",
  ], { fontSize: 12, bar: C.navy });
  S.table(s, [
    ["control template (15 wells)", "value"],
    ["peak active force", "429 µN (well-to-well SD 48)"],
    ["stimulus → peak / 10 % → peak", "230 / 176 ms"],
    ["peak → 50 % / peak → 90 %", "152 / 243 ms"],
    ["shape noise (SEM of the template)", "0.35 % of peak"],
  ], { x: M, y: TOP + 3.15, w: 6.0, h: 1.7 }, { fontSize: 12, colW: [3.4, 2.6] });
  take(s, "Why the platemap matters: the control set is an inference from the response pattern. If any of the 15 wells had a compound, the template and every fit on it move. The cross-sectional area is the other missing number: without it the model's kPa cannot be grounded in the measured µN.", "FDECEA", B);
  S.footer(s, "Analyses/CalciumTwitch/PrepareEHTControls.m; EHT_identification_report.md §1");
}

// 5 Failure
{
  const s = pres.addSlide();
  S.header(s, "3 · Mouse model + rodent Ca vs the human twitch", "total failure");
  msg(s, "A MyoSim-type rat Ca transient + published thin-filament constants drive the mouse model. Compare with the EHT twitch.");
  img(s, "F2_failure.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "The model twitch peaks at ~70 ms where the tissue peaks at 230 ms, and relaxes three times faster. Shape wrong, amplitude 5 % of maximal.", "FDECEA", B);
  S.footer(s, "Ca: Rice et al. 2008 (rat, 24 °C); thin filament: Marzban 2020 constants; core: ExampleHighATP.");
}

// 6 Fit C
{
  const s = pres.addSlide();
  S.header(s, "4 · Keep the mouse core, fit the Ca transient", "match, at a price");
  msg(s, "The twitch is matched to within 0.8 % of peak. But look at the Ca transient the twitch demands.");
  img(s, "F3_fitC.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "The fitted Ca peaks at ~110 ms and decays over ~200 ms: 2–4× slower than the rat transient. Plausible for hiPSC cells, but it is a fitted guess, and a 2× faster or slower core fits equally well with a different Ca.", "FFF6D6", C.maize);
  S.footer(s, "EHT report fit C: Ca shape, levels (assumed 0.1 → 1 µM) and 3 thin-filament constants free; cross-bridge core frozen.");
}

// 7 Human vs mouse (literature)
{
  const s = pres.addSlide();
  S.header(s, "5 · How different is the human cross-bridge?", "literature");
  S.table(s, [
    ["observable", "human / mouse (same lab, same T)", "confidence", "what it means for the model"],
    ["ktr (rate of force redevelopment)", "0.10 – 0.12", { text: "high", color: G }, "cycle ~10× slower"],
    ["unloaded shortening velocity V0", "0.10 – 0.15", { text: "medium", color: Wn }, "ADP release ~10× slower"],
    ["ATP turnover per head; tension cost", "0.1 – 0.2; 0.25", { text: "medium / high", color: G }, "economy 4–10× better"],
    ["ADP release; ATP hydrolysis (solution)", "≤ 0.1; 0.1", { text: "high", color: G }, "the β-MHC isoform: kinetics only"],
    ["Pi release / power stroke", "0.2 – 0.5", { text: "low", color: B }, "weakly constrained"],
    ["F0, step size, stiffness, attached fraction", "≈ 1", { text: "high", color: G }, "mechanics unchanged: attachment must slow with detachment"],
    ["resting SRX fraction", "≈ 1 (no defensible difference)", { text: "low", color: B }, "leave alone, test as a sensitivity"],
    ["temperature: rig 20–25 °C → tissue 37 °C", "× 3.7 – 8.4 on every rate (Q10 3–3.5)", { text: "rig T unknown", color: Wn }, "the biggest single factor, and it is not biology"],
  ], { x: M, y: TOP, w: W - 2 * M, h: 4.6 }, { fontSize: 14, colW: [3.6, 3.0, 1.6, 4.13] });
  take(s, "Human ventricle is β-MHC: same force and stiffness, ~10× slower kinetics at the same temperature. At 37 °C that slowing is partly undone by the Q10.");
  S.footer(s, "Analyses/HumanizeXB/conclusions.md and literature_log.md (same-lab, matched-temperature ratios).");
}

// 8 Fixed Ca fit
{
  const s = pres.addSlide();
  S.header(s, "6 · Human iPSC Ca fixed, human changes only", "match");
  img(s, "F5_tier1.png", { x: M, y: TOP, w: 9.4, h: 3.7 });
  S.card(s, { x: M + 9.55, y: TOP, w: W - 2 * M - 9.55, h: 3.7 }, "Why de Lange 2021", [
    "3D hiPSC tissue, Fura-2 (ratiometric), 37 °C, 1 Hz: the same preparation class",
    "their force twitch (TTP 229, RT50 160, RT90 262 ms) is ours within 10 ms",
    "absolute Ca levels are a guess: no calibrated EHT value exists. They turn out not to matter (only k_on·Ca enters)",
  ], { fontSize: 12, bar: C.navy });
  take(s, "Free: one isoform factor on ADP release, hydrolysis, attachment and weak detachment; one temperature factor on all cycle rates; the thin filament; the beat offset. Result 0.77 % of peak, as good as fitting the Ca freely.", "E8F3EA", G, 5.2);
  S.footer(s, "Analyses/HumanTwitch/conclusions.md §2 (tier 1).");
}

// 9 What changed
{
  const s = pres.addSlide();
  S.header(s, "7 · What exactly changed in that fit", "parameters");
  S.table(s, [
    ["parameter", "role", "mouse", "fitted", "factor", "how"],
    ["ka / kd", "attachment / weak detachment (s⁻¹)", "292 / 93", "133 / 42", { text: "×0.45", bold: true }, "isoform × temperature, one shared factor"],
    ["k2 / k_2", "ADP release + ATP binding (s⁻¹)", "117 / 4.4", "53 / 2.0", { text: "×0.45", bold: true }, "same factor"],
    ["kah / kamh", "ATP hydrolysis (s⁻¹)", "176 / 8.8", "80 / 4.0", { text: "×0.45", bold: true }, "same factor"],
    ["k1 / k_1", "power stroke, Pi release (s⁻¹)", "161 / 31", "335 / 65", "×2.1", "literature ×0.3 × temperature 6.9 (not seen by the twitch)"],
    ["tf_k_on, tf_k_off, K_coop", "thin filament: Ca on-rate, site off-rate, cooperativity", "101, 724, 9.7 (rat)", "9.2 µM⁻¹s⁻¹, 5.1 s⁻¹, 20", "fitted", "never calibrated before; K_coop at its bound"],
    ["beat offset; σ", "stimulus → Ca onset; µN per model kPa", "—", "+1.7 ms; 9.4", "fitted / closed form", "stimulus record not synchronised; cross-section unknown"],
    ["SL0", "sarcomere length", "2.2 µm", "2.0 µm", "set", "EHT value unknown"],
    ["stiffness, SRX rates, titin, kSE, gates", "", "mouse", "mouse", { text: "unchanged", color: G }, "no human evidence; tested one by one: ≤ 0.08 pp gain each"],
  ], { x: M, y: TOP, w: W - 2 * M, h: 4.4 }, { fontSize: 12, colW: [2.1, 3.1, 1.5, 1.9, 1.2, 2.53] });
  take(s, "The isoform and temperature factors are perfectly anti-correlated: the twitch identifies only their product, ×0.45 (×/÷ 1.12). Read as isoform 0.1 × Q10 factor 4.5, i.e. a rig at 22–24 °C. Saturating force stays at the mouse value; ktr becomes 23 s⁻¹, ATP use per head halves.");
  S.footer(s, "robust1_seed → params/human_twitch_t1_isoT.m; conclusions.md §2 and the Jacobian in results/t1_isoT_jac.mat");
}

// 10 Paced Ca span
{
  const s = pres.addSlide();
  S.header(s, "8 · Why we still need the Ca transient measured", "the literature disagrees");
  msg(s, "Three paced, 37 °C human transients from the literature: the time to peak alone spans 60 to 190 ms.");
  img(s, "F6b_ca_paced.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "The two hiPSC tissue measurements differ 2× in time to peak (Fura-2 vs GCaMP6f, which lags). Spontaneous single-cell recordings and hiPSC models spread even wider (up to 270 ms).");
  S.footer(s, "literature_ca_transient.md; Stage1b_CaVariability.m (each curve reproduces the metrics the source reports).");
}

// 11 Consequence
{
  const s = pres.addSlide();
  S.header(s, "9 · The fitted core tracks the assumed Ca", "the catch");
  msg(s, "Same twitch, same fit quality, three Ca transients from inside the literature range: the rates asked for differ 18×.");
  img(s, "F7_consequence.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "A 20 ms change in the assumed Ca time to peak is a factor 2 on every rate of the 'human' core. Absolute Ca levels do not matter, so a RELATIVE Ca recording in these tissues is enough. This is the first thing to ask for.", "FDECEA", B);
  S.footer(s, "conclusions.md §2.2");
}

// 12 Case B
{
  const s = pres.addSlide();
  S.header(s, "10 · With Ca known, cycle speed is open upwards", "and possibly some other");
  msg(s, "Ca fixed: fits degrade steeply if the core is slower than ×0.5 of the mouse, but a core ×1 … ×10 fits equally well.");
  img(s, "F8_caseB.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "Equally good fits predict ktr 23 … 97 s⁻¹ and 1.5 … 24 ATP per head per beat. What they share: the thin filament, and the force response to Ca amplitude (×1.5: +11…+14 %) and to an attachment inhibitor (ka ×0.5: −49…−55 %).");
  S.footer(s, "conclusions.md §5: profile likelihood with the Ca transient treated as measured; PredictSet.m");
}

// 13 Two cores
{
  const s = pres.addSlide();
  S.header(s, "11 · Two cores, one twitch: where the heads are", "DRX / SRX");
  img(s, "F9_two_cores.png", { x: M, y: TOP - 0.05, w: W - 2 * M, h: 4.95 });
  take(s, "The fast core keeps half its heads parked in SRX at rest and empties the pool every beat; the slow pool is what shapes its relaxation. The humanized core parks a third and uses 13× less ATP. The twitch cannot tell them apart; an ATPase, O₂ or ktr measurement can.", "FFF6D6", C.maize);
  S.footer(s, "tier-1 fit vs the Ca-known profile point at cycle ×10; both at the fixed de Lange transient. ATP per head per beat: 1.8 vs 24.");
}

// 14 Compliance and SL
{
  const s = pres.addSlide();
  S.header(s, "12 · Do compliant posts and SL matter?", "series compliance, SL");
  img(s, "F10_compliance_SL.png", { x: M, y: TOP - 0.05, w: W - 2 * M, h: 4.95 });
  take(s, "Without refitting, a 2–10× softer series element or a 0.1 µm different resting SL changes the twitch shape by 3–10 % (far more than any fit residual), so a fit absorbs them silently through other parameters. The predictions move with them: the force response to a 1.5× Ca transient spans +3 … +19 %, to an attachment inhibitor −37 … −52 %.", "FFF6D6", C.maize);
  S.footer(s, "humanized (tier-1) core, Ca fixed, no refit; kSE is the model's series element, SL0 the resting sarcomere length. Refits: kSE free → ×0.78 (rms 0.69 %); SL0 free → 1.91 µm (0.71 %).");
}

// 15 Compliance: how to assert
{
  const s = pres.addSlide();
  S.header(s, "12 · How to pin compliance and sarcomere length", "series compliance, SL");
  S.card(s, { x: M, y: TOP, w: 6.0, h: 4.7 }, "What the model needs", [
    "the post (or membrane) stiffness in µN/µm and how force is computed from deflection → the series compliance the sarcomeres work against",
    "tissue length and cross-section → tissue shortening per twitch in % of length, and kPa from µN",
    "resting sarcomere length (one image) → the overlap the model starts from",
    "ideally: force and displacement recorded together, so the model can be run auxotonically instead of isometrically",
  ], { fontSize: 13, bar: C.navy });
  S.card(s, { x: M + 6.2, y: TOP, w: W - 2 * M - 6.2, h: 4.7 }, "What the model predicts for it", [
    "contractile-element shortening per twitch: ~20 nm per half-sarcomere (1 % of SL) at the trabecula's series stiffness, ~120 nm (6 %) at a 10× softer one: measurable as post deflection or tissue shortening",
    "peak force falls 38 % and the rise slows by 80 ms (235 → 312 ms) when the series element is 10× softer; resting SL 1.9 vs 2.2 µm moves peak force 35 → 54 kPa",
    "response to a 1.5× Ca transient: +11 % (stiff) … +3 % (soft); to attachment ×0.5: −47 … −37 %. SL 1.9 → 2.2 µm: +19 … +7 % and −52 … −39 %",
    "the response of peak force to an attachment inhibitor (mavacamten-like) and to Ca amplitude is the quantity to compare across posts of different stiffness",
  ], { fontSize: 13, bar: C.maize });
  take(s, "Assertion route: measure post stiffness and tissue geometry once; image resting SL once; then the model's isometric-vs-auxotonic difference is a computed correction, not an assumption.");
  S.footer(s, "results/deck/F10_compliance.mat");
}

// 16 Asks: measurements
{
  const s = pres.addSlide();
  S.header(s, "13 · What we need from the experiments", "ask");
  S.table(s, [
    ["", "measurement", "what it settles", "the prediction it tests"],
    ["1", { text: "relative Ca transient in the same wells at 1 Hz (ratiometric or low-affinity dye; not GCaMP)", bold: true }, "the Ca-shape assumption (×0.2 … ×4 on every rate); the thin filament", "Ca peaks 130–145 ms before force"],
    ["2", { text: "cross-sectional area and length of the tissues", bold: true }, "grounds model kPa in measured µN (σ); separates 20–30 kPa cores from 50–65 kPa cores", "twitch stress ≈ 0.8 × maximal; ~35–45 kPa for the humanized core"],
    ["3", { text: "paced O₂ consumption or ATPase of the same tissues (± mavacamten)", bold: true }, "the cycle speed: equally good fits differ 16× in ATP per beat", "1.8 (humanized) … 24 (fast core) ATP per head per beat"],
    ["4", "ktr / kACT on skinned EHT at 37 °C, or on myofibrils", "cycle speed directly", "23 s⁻¹ (humanized) vs 86–99 (fast core)"],
    ["5", "[Ca]o series (extracellular Ca in the medium, e.g. 0.5–3.6 mM): no dye needed", "scales the transient amplitude → tests the near-saturation operating point", "peak force +11…+14 % for a 1.5× transient, −53…−59 % for ÷1.5"],
    ["6", "force–pCa on skinned EHT of the same line", "thin-filament sensitivity and steepness", "pCa50 ≈ 6.2 vs peak Ca; n_H ≥ 4.6"],
  ], { x: M, y: TOP, w: W - 2 * M, h: 4.4 }, { fontSize: 12, colW: [0.4, 4.6, 4.2, 3.13] });
  take(s, "Order: 1 and 2 first (nothing substitutes for them), then 3 or 4. mant-ATP would not help here: it reads the turnover of relaxed heads (the model's diastolic rate, ~0.1 s⁻¹, resolvable but not the cycle speed), and an activated tissue does not keep the relaxed SRX.");
  S.footer(s, "conclusions.md §4–5");
}

// 17 Asks: facts and optional
{
  const s = pres.addSlide();
  S.header(s, "14 · Facts to ask for, and optional extras", "ask");
  S.card(s, { x: M, y: TOP, w: 6.0, h: 4.9 }, "Facts the lab most likely already has", [
    "Plate 6 platemap: what is in columns 1–3, is row D of 4–6 DMSO, dosing time",
    "post stiffness (µN/µm), tissue length and cross-section, how force is computed",
    "cell line, tissue age, cardiomyocyte purity, fibroblast content",
    "recording temperature, medium and its [Ca]o",
    "stimulus–force synchronisation; raw exports with the absolute baseline",
    "pre-dose baseline recordings of every plate; sarcomere length if imaged",
  ], { fontSize: 14, bar: C.navy });
  S.card(s, { x: M + 6.2, y: TOP, w: W - 2 * M - 6.2, h: 4.9 }, "If available or possible", [
    "pacing-frequency series (0.5–2 Hz) and post-rest potentiation",
    "a temperature step (e.g. 32 vs 37 °C): the Q10 of the human core directly",
    "drug waveforms, with Ca where possible: mavacamten (cross-bridge side), isoproterenol, diltiazem (Ca side)",
    "Ca calibration (Fura-2 with ionomycin) if the dye work happens",
    "MHC and TnI isoform fractions of the line",
  ], { fontSize: 14, bar: C.maize });
  S.footer(s, "EHT_identification_report.md §7; conclusions.md §4");
}

// 18 Summary
{
  const s = pres.addSlide();
  S.header(s, "Summary", "take-home");
  const rows = [
    ["1", "The mouse model reproduces the pooled mouse data at 8 and 2 mM ATP."],
    ["2", "With a rodent Ca transient it fails on the human iPSC twitch: wrong timing, 5 % amplitude."],
    ["3", "Fitting the Ca transient rescues it, but the Ca must be 2–4× slower than rodent, and it is a guess."],
    ["4", "Human β-MHC is ~10× slower than mouse at the same temperature; at 37 °C part of that is undone."],
    ["5", "With a literature human iPSC Ca fixed, one isoform × temperature factor (×0.45) plus the thin filament matches the twitch (0.77 %)."],
    ["6", "But paced literature Ca transients disagree 2–3×, and the fitted core follows the assumption 18×: we need the Ca measured."],
    ["7", "Even then the twitch bounds the cycle speed only from below; cross-section plus an ATPase / O₂ or ktr measurement are the next asks."],
  ];
  rows.forEach(([n, t], i) => {
    const y = TOP + i * 0.62;
    s.addShape("ellipse", { x: M, y: y + 0.06, w: 0.45, h: 0.45, fill: { color: i >= 5 ? B : C.navy }, line: { color: i >= 5 ? B : C.navy } });
    s.addText(n, { x: M, y: y + 0.06, w: 0.45, h: 0.45, fontFace: F, fontSize: 15, bold: true, color: C.white, align: "center", valign: "middle", margin: 0 });
    s.addText(t, { x: M + 0.65, y, w: W - 2 * M - 0.65, h: 0.58, fontFace: F, fontSize: 17, color: C.ink, valign: "middle", margin: 0 });
  });
  take(s, "Bottom line: the twitch shape is reproduced to noise, but the fitted rates are not yet a statement about the human muscle. A Ca recording, the tissue geometry and one energetic or kinetic anchor would make them one.");
  S.footer(s, "Worktree C:\\home\\git\\ATP-wt-human (branch HumanTwitch); Analyses/HumanTwitch/conclusions.md");
}

pres.writeFile({ fileName: path.resolve(__dirname, "..", "2026-10-HumanTwitch_v3.pptx") }).then((f) => console.log("wrote " + f));
