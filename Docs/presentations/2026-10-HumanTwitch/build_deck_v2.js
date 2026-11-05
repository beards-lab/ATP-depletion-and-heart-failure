// Storyline deck v2 (figure-led). Build: NODE_PATH=~/tools/pptx/node_modules node build_deck_v2.js
const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const S = require("../2026-09-status/style");
const { C, F, W, M, TOP } = S;
const DECK = path.resolve(__dirname, "..", "..", "..", "Analyses", "HumanTwitch", "results", "deck");
const fig = (n) => path.join(DECK, n);
function pngSize(abs) { if (!fs.existsSync(abs)) return null; const b = fs.readFileSync(abs); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
function img(s, name, box) { S.image(s, fig(name), box, pngSize(fig(name))); }
// one-line message under the header, big and plain
function msg(s, text, color) {
  s.addText(text, { x: M, y: TOP - 0.05, w: W - 2 * M, h: 0.55, fontFace: F, fontSize: 17, bold: true, color: color || C.navy, valign: "middle", margin: 0 });
}
// bottom takeaway box
function take(s, text, fill, line) {
  s.addShape("rect", { x: M, y: 6.1, w: W - 2 * M, h: 0.8, fill: { color: fill || "FFF6D6" }, line: { color: line || C.maize, width: 1 } });
  s.addText(text, { x: M + 0.2, y: 6.1, w: W - 2 * M - 0.4, h: 0.8, fontFace: F, fontSize: 15, color: C.ink, valign: "middle", margin: 0 });
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
  s.addNotes("Storyline: 1) the mouse model reproduces the mouse data; 2) with a rodent Ca transient it fails on the human iPSC twitch; 3) fitting the Ca transient rescues it, but the Ca must be very different; 4) how different is human from mouse (literature); 5) with a literature human Ca fixed, changing only the plausibly human parameters matches; 6) therefore the Ca transient must be measured, and one more anchor; 7) the asks.");
}

// 2 Mouse model reproduces mouse data
{
  const s = pres.addSlide();
  S.header(s, "1 · The mouse model reproduces the mouse data", "where we start");
  msg(s, "One parameter set, five protocols on one preparation: force–velocity, ktr, staircase, slack–restretch, passive.");
  img(s, "F1_example_fig1.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "Black: Baker-lab data (mouse RV trabecula, skinned, saturating Ca, 8 mM ATP). Blue: model. This is the model we take to the human tissue.");
  S.footer(s, "Example/RunExampleHighATP.m with params/ExampleHighATP.m; feature cost 0.94.");
  s.addNotes("The core is identified at saturating Ca only. Nothing in it says how force depends on Ca; the thin-filament layer is new territory for the twitch.");
}

// 3 Failure
{
  const s = pres.addSlide();
  S.header(s, "2 · Mouse model + rodent Ca vs the human twitch", "total failure");
  msg(s, "A MyoSim-type rat Ca transient + published thin-filament constants drive the mouse model. Compare with the EHT twitch.");
  img(s, "F2_failure.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "The model twitch peaks at ~70 ms where the tissue peaks at 230 ms, and relaxes three times faster. Shape wrong, amplitude 5 % of maximal: a failure, as expected.", "FDECEA", B);
  S.footer(s, "Ca: Rice et al. 2008 (rat, 24 °C); thin filament: Marzban 2020 constants; core: ExampleHighATP. EHT data: Travis, Plate 6, 15 control wells, 1 Hz.");
}

// 4 Fit the Ca, keep the core
{
  const s = pres.addSlide();
  S.header(s, "3 · Keep the mouse core, fit the Ca transient", "match, at a price");
  msg(s, "Now the twitch is matched to within 0.8 % of peak. But look at the Ca transient the twitch demands.");
  img(s, "F3_fitC.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "The fitted Ca peaks at ~110 ms and decays over ~200 ms: 2–4× slower than the rat transient. Plausible for hiPSC cells, but it is a fitted guess, and the fit is not unique (a 2× faster or slower core fits equally well).", "FFF6D6", C.maize);
  S.footer(s, "EHT report fit C: Ca shape, levels (assumed 0.1 → 1 µM) and 3 thin-filament constants free; cross-bridge core frozen. Analyses/CalciumTwitch/EHT_identification_report.md");
}

// 5 How different is human from mouse
{
  const s = pres.addSlide();
  S.header(s, "4 · How different is the human cross-bridge?", "literature");
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

// 5b Mouse vs human illustrated
{
  const s = pres.addSlide();
  S.header(s, "4 · Mouse vs human: force–velocity, ktr, ATP use", "illustration");
  msg(s, "The mouse core, the literature human core at rig temperature and at 37 °C, and the core fitted to the twitch on the next slide.");
  img(s, "F4_mouse_vs_human.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "Same isometric force in all four. The human core at rig temperature is 8× slower than the mouse (blue); brought to 37 °C it is ~2× slower (orange), and that is almost exactly where the twitch fit lands (red).", "E8F3EA", G);
  S.footer(s, "humanProbe at saturating Ca, SL 2.2 µm, 8 mM ATP. Human core = mouse with k2, k_2, kah, kamh, ka, kd ×0.1 and k1, k_1 ×0.3; 37 °C = ×5 on all cycle rates (Q10 3 over 15 °C).");
}

// 6 Fix Ca from literature, vary only the human directions
{
  const s = pres.addSlide();
  S.header(s, "5 · Human iPSC Ca fixed, human parameters only", "match");
  msg(s, "Ca fixed from de Lange 2021 (Fura-2, 37 °C, 1 Hz, 3D hiPSC tissue). Free: isoform factor, temperature factor, thin filament. Rest = mouse.");
  img(s, "F5_tier1.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "Match: 0.77 % of peak, as good as fitting the Ca freely. The twitch asks for the six isoform-sensitive rates at ×0.45 of the mouse (isoform ~0.1 × temperature ~4.5). Pi release, SRX rates, series compliance, titin, SL: little gain.", "E8F3EA", G);
  S.footer(s, "Analyses/HumanTwitch/conclusions.md §2 (tier 1); the two factors are correlated −1.00: only their product is identified.");
}

// 7 Why we need the Ca transient (variability)
{
  const s = pres.addSlide();
  S.header(s, "6 · Why we still need the Ca transient measured", "because the literature disagrees");
  msg(s, "Because 'the human iPSC Ca transient' does not exist in the literature: time to peak ranges from 60 to 270 ms.");
  img(s, "F6_ca_variability.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "Sources differ in maturity, 3D vs single cell, paced vs spontaneous, and indicator (Fura-2, Fluo-4, GCaMP). Even the two paced 37 °C tissue measurements differ 2× in time to peak.");
  S.footer(s, "Analyses/HumanTwitch/literature_ca_transient.md; Stage1b_CaVariability.m");
}

// 8 Consequence
{
  const s = pres.addSlide();
  S.header(s, "7 · The fitted core tracks the assumed Ca", "the catch");
  msg(s, "Same twitch, same fit quality, three Ca transients from inside the literature range: the rates asked for differ 18×.");
  img(s, "F7_consequence.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "A 20 ms change in the assumed Ca time to peak is a factor 2 on every rate of the 'human' core. Absolute Ca levels do not matter at all, so a RELATIVE Ca recording in these tissues is enough. This is the first thing to ask for.", "FDECEA", B);
  S.footer(s, "conclusions.md §2.2: tier-1 refits with the transient at the fast / nominal / slow edge of the range.");
}

// 9 And something else
{
  const s = pres.addSlide();
  S.header(s, "8 · With Ca known, one more anchor is needed", "and possibly some other");
  msg(s, "Even with the Ca fixed, the twitch bounds the cycle speed only from below: faster cores fit as well by slowing the SRX exchange.");
  img(s, "F8_caseB.png", { x: M, y: TOP + 0.6, w: W - 2 * M, h: 4.25 });
  take(s, "Reliable anyway: the thin filament, and how peak force responds to Ca amplitude (×1.5: +11…+14 %) or an attachment inhibitor (ka ×0.5: −49…−55 %). Not predicted: ktr (23–97 s⁻¹) and ATP cost (1.5–24 per head per beat). One ATPase/O₂ or ktr measurement closes this.");
  S.footer(s, "conclusions.md §5: profile likelihood with the Ca transient treated as measured (± indicator error); PredictSet.m");
}

// 10 Asks: measurements
{
  const s = pres.addSlide();
  S.header(s, "9 · What we need from the experiments", "ask");
  S.table(s, [
    ["", "measurement", "what it settles", "the prediction it tests"],
    ["1", { text: "relative Ca transient in the same wells at 1 Hz (ratiometric or low-affinity dye; not GCaMP)", bold: true }, "the Ca-shape assumption (×0.2 … ×4 on every rate); the thin filament", "Ca peaks 130–145 ms before force"],
    ["2", { text: "paced O₂ consumption or ATPase of the same tissues (± mavacamten)", bold: true }, "the cycle speed: equally good fits differ 16× in ATP per beat", "1.5 … 24 ATP per head per beat"],
    ["3", { text: "force per cross-section, or skinned-EHT maximal force", bold: true }, "the force scale σ; separates 20–30 kPa cores from 50–65 kPa cores", "twitch ≈ 0.8 × maximal force"],
    ["4", "ktr / kACT on skinned EHT at 37 °C (or myofibrils)", "cycle speed directly", "23 s⁻¹ (humanized core) vs 86–99 (fast alternative)"],
    ["5", "force–pCa on skinned EHT of the same line", "thin-filament sensitivity and steepness", "pCa50 ≈ 6.2 vs peak Ca; n_H ≥ 4.6"],
  ], { x: M, y: TOP, w: W - 2 * M, h: 3.9 }, { fontSize: 13, colW: [0.4, 4.6, 4.2, 3.13] });
  take(s, "Order matters: 1 first (nothing else can substitute for it), then 2 or 3. Each removes one specific ambiguity of the fit.");
  S.footer(s, "conclusions.md §4–5");
}

// 11 Asks: facts and optional
{
  const s = pres.addSlide();
  S.header(s, "10 · Facts to ask for, and optional extras", "ask");
  S.card(s, { x: M, y: TOP, w: 6.0, h: 4.9 }, "Facts the lab most likely already has", [
    "Plate 6 platemap (what is in columns 1–3; is row D of 4–6 DMSO; dosing time)",
    "cell line, tissue age, cardiomyocyte purity, fibroblast content",
    "recording temperature, [Ca]o and medium",
    "post stiffness (µN/µm) and how force is computed from displacement",
    "stimulus–force synchronisation; raw exports with the absolute baseline",
    "pre-dose baseline recordings of every plate",
    "tissue cross-section and length; sarcomere length if imaged",
  ], { fontSize: 14, bar: C.navy });
  S.card(s, { x: M + 6.2, y: TOP, w: W - 2 * M - 6.2, h: 4.9 }, "If available or possible", [
    "[Ca]o series (0.5–3.6 mM): tests the near-saturation prediction (+11…+14 % for ×1.5, −53…−59 % for ÷1.5)",
    "pacing-frequency series and post-rest potentiation",
    "a temperature step (e.g. 32 vs 37 °C): the Q10 of the human core directly",
    "drug waveforms with Ca where possible: mavacamten (cross-bridge side), isoproterenol, diltiazem (Ca side)",
    "Ca calibration (Fura-2 with ionomycin) if the dye work happens",
    "mant-ATP SRX fraction; MHC and TnI isoform fractions of the line",
  ], { fontSize: 14, bar: C.maize });
  S.footer(s, "EHT_identification_report.md §7; conclusions.md §4");
}

// 12 Summary
{
  const s = pres.addSlide();
  S.header(s, "Summary", "take-home");
  const rows = [
    ["1", "The mouse model reproduces the mouse data (FV, ktr, stairs, slack)."],
    ["2", "With a rodent Ca transient it fails on the human iPSC twitch: wrong timing, 5 % amplitude."],
    ["3", "Fitting the Ca transient rescues it, but the Ca must be 2–4× slower than rodent, and it is a guess."],
    ["4", "Human β-MHC is ~10× slower than mouse at the same temperature; at 37 °C part of that is undone."],
    ["5", "With a literature human iPSC Ca fixed, moving only the plausibly human parameters matches the twitch (0.77 %)."],
    ["6", "But the literature Ca transients disagree 4×, and the fitted core follows the assumption 18×: we need the Ca measured."],
    ["7", "Even then the twitch bounds the cycle speed only from below: one ATPase / O₂ or ktr measurement is the second ask."],
  ];
  rows.forEach(([n, t], i) => {
    const y = TOP + i * 0.62;
    s.addShape("ellipse", { x: M, y: y + 0.06, w: 0.45, h: 0.45, fill: { color: i >= 5 ? B : C.navy }, line: { color: i >= 5 ? B : C.navy } });
    s.addText(n, { x: M, y: y + 0.06, w: 0.45, h: 0.45, fontFace: F, fontSize: 15, bold: true, color: C.white, align: "center", valign: "middle", margin: 0 });
    s.addText(t, { x: M + 0.65, y, w: W - 2 * M - 0.65, h: 0.58, fontFace: F, fontSize: 17, color: C.ink, valign: "middle", margin: 0 });
  });
  take(s, "Bottom line: the twitch shape is reproduced to noise, but the fitted rates are not yet a statement about the human muscle. One Ca recording and one energetic or kinetic anchor would make them one.");
  S.footer(s, "Worktree C:\\home\\git\\ATP-wt-human (branch HumanTwitch); Analyses/HumanTwitch/conclusions.md");
}

pres.writeFile({ fileName: path.resolve(__dirname, "..", "2026-10-HumanTwitch_v2.pptx") }).then((f) => console.log("wrote " + f));
