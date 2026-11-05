// Human iPSC-EHT twitch deck (Oct 2026). Build: NODE_PATH=~/tools/pptx/node_modules node build_deck.js
const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const S = require("../2026-09-status/style");
const { C, F, W, M, TOP } = S;
const RES = path.resolve(__dirname, "..", "..", "..", "Analyses", "HumanTwitch", "results");
const fig = (n) => path.join(RES, n);
function pngSize(abs) { if (!fs.existsSync(abs)) return null; const b = fs.readFileSync(abs); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
function img(s, name, box) { S.image(s, fig(name), box, pngSize(fig(name))); }
function note(s, t) { s.addNotes(t); }
function callout(s, box, text, fill, line) {
  s.addShape("rect", { ...box, fill: { color: fill || "FFF6D6" }, line: { color: line || C.maize, width: 1 } });
  s.addText(text, { x: box.x + 0.15, y: box.y, w: box.w - 0.3, h: box.h, fontFace: F, fontSize: 14, color: C.ink, valign: "middle", margin: 0 });
}
const G = C.good, Wn = C.warn, B = C.bad;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "Human iPSC-EHT twitch: what the cross-bridge model can say without a Ca measurement";

// 1 Title
{
  const s = pres.addSlide();
  s.background = { color: C.navy };
  s.addShape("rect", { x: 0, y: 4.55, w: W, h: 0.08, fill: { color: C.maize }, line: { color: C.maize } });
  s.addText("Fitting the human iPSC-EHT twitch with the mouse cross-bridge model: what we found, and what we need from the experiments",
    { x: 0.8, y: 1.1, w: W - 1.6, h: 2.9, fontFace: F, fontSize: 34, bold: true, color: C.white, valign: "bottom", margin: 0 });
  s.addText("Filip Jezek", { x: 0.8, y: 4.85, w: W - 1.6, h: 0.5, fontFace: F, fontSize: 22, color: C.maize, margin: 0 });
  s.addText("October 2026 · Plate 6 control EHTs (Travis), 1 Hz field stimulation", { x: 0.8, y: 5.4, w: W - 1.6, h: 0.45, fontFace: F, fontSize: 18, color: C.white, margin: 0 });
  s.addText("Three questions: can a literature Ca transient plus 'humanizing' changes explain the twitch? What does the twitch identify? What should we measure next?",
    { x: 0.8, y: 6.3, w: W - 1.6, h: 0.6, fontFace: F, fontSize: 14, color: "C9D3E0", italic: true, margin: 0 });
  note(s, "Deck summarises Analyses/HumanTwitch/conclusions.md (worktree ATP-wt-human, branch HumanTwitch). Audience: the people who ran or will run the EHT experiments plus the modelling group. Each 'ask' slide says what the measurement buys.");
}

// 2 Setting
{
  const s = pres.addSlide();
  S.header(s, "Mouse-calibrated model, human iPSC tissue", "context");
  S.table(s, [
    ["", "what the cross-bridge core was fitted to", "what we fit here"],
    ["cells / tissue", "adult mouse RV trabecula, skinned (no SR, no Ca transient)", "human iPSC-derived cardiomyocytes, intact 3D tissue on compliant posts (Mantarray)"],
    ["activation", "clamped saturating Ca (pCa 4.51)", "field-stimulated twitch at 1 Hz; Ca transient unmeasured; stimulus trace not synchronised with force"],
    ["isoforms, T", "α-MHC, cTnI; rig temperature not recorded (20–25 °C assumed)", "β-MHC, slow-skeletal TnI, fetal titin, immature; presumably 37 °C"],
    ["readout", "absolute stress, 69 kPa maximal", "active force in µN, baseline-subtracted (no diastolic force); peak 429 µN"],
  ], { x: M, y: TOP, w: W - 2 * M, h: 2.9 }, { fontSize: 12, colW: [1.7, 5.1, 5.53] });
  S.card(s, { x: M, y: 4.3, w: 6.0, h: 2.55 }, "The target twitch (15 control wells × 60 beats)", [
    "stimulus → peak 230 ms · 10 % → peak 176 ms",
    "peak → 50 % 152 ms · peak → 90 % 243 ms · width 275 ms",
    "shape is nearly identical across wells (template SEM 0.35 % of peak); amplitude varies 11 %",
    "one scale σ (µN per model kPa) absorbs cross-section and myocyte content",
  ], { fontSize: 13, bar: C.navy });
  S.card(s, { x: M + 6.2, y: 4.3, w: W - 2 * M - 6.2, h: 2.55 }, "The model layers", [
    "cross-bridge core: Beard 2022 strain-discretised cycle + SRX, mouse-fitted (robust1_seed)",
    "thin filament: Campbell/Tewari two-state gate on attachment (k_on, k_off, K_coop) — never calibrated on any data",
    "Ca transient: prescribed (this work: fixed from literature)",
    "boundary: isometric with the model's series element",
  ], { fontSize: 13, bar: C.maize });
  S.footer(s, "EHT report: Analyses/CalciumTwitch/EHT_identification_report.md (fit C: free Ca transient, 0.81 % rms). This work: Analyses/HumanTwitch/conclusions.md");
  note(s, "Earlier (fit C) the Ca transient was fitted freely and the core left as mouse. Now: Ca fixed from literature, core moved along human/rodent literature directions, then freed.");
}

// 3 The Ca transient
{
  const s = pres.addSlide();
  S.header(s, "Step 1 · A fixed Ca transient from the literature", "Ca transient");
  img(s, "stage1_ca_transients.png", { x: M, y: TOP, w: W - 2 * M, h: 3.4 });
  S.table(s, [
    ["transient", "dia → peak (µM)", "time to peak", "peak → 50 %", "peak → 90 %"],
    ["hiPSC-EHT, used here (de Lange 2021 anchor: Fura-2, 37 °C, 1 Hz, 3D tissue)", "0.15 → 0.60", "93 ms", "186 ms", "448 ms"],
    ["adult human LV (Coppini 2013, ToR-ORd)", "0.10 → 0.45", "70 ms", "153 ms", "391 ms"],
    ["fit C's freely fitted transient (EHT report)", "0.10 → 1.0 (assumed)", "110 ms", "185 ms", "428 ms"],
  ], { x: M, y: 4.7, w: 8.3, h: 1.6 }, { fontSize: 10, colW: [4.1, 1.3, 1.0, 1.0, 0.9] });
  S.card(s, { x: M + 8.5, y: 4.7, w: W - 2 * M - 8.5, h: 2.2 }, "Why de Lange 2021", [
    "same preparation class (3D hiPSC tissue), ratiometric dye, 37 °C, 1 Hz",
    "their force twitch (TTP 229, RT50 160, RT90 262 ms) is ours within 10 ms",
    "absolute levels are guesses (no calibrated EHT value exists) — and they turn out not to matter",
  ], { fontSize: 12, bar: C.navy });
  s.addText("The freely fitted transient of fit C already had the literature decay; only its rise was 20 ms too slow. Ca peaks ~140 ms before force.",
    { x: M, y: 6.55, w: 8.3, h: 0.4, fontFace: F, fontSize: 11, italic: true, color: C.navy, margin: 0 });
  S.footer(s, "Analyses/HumanTwitch/literature_ca_transient.md; caTransient kind 'risedecay' (exponential rise × Hill decay, beat superposition). Panel C: zero-tuning baselines.");
}

// 4 Literature variability
{
  const s = pres.addSlide();
  S.header(s, "There is no 'the' human iPSC Ca transient", "Ca variability");
  img(s, "stage1b_ca_variability.png", { x: M, y: TOP, w: W - 2 * M, h: 3.95 });
  S.table(s, [
    ["source", "preparation, indicator, T", "TTP (ms)", "peak → 50 %"],
    ["de Lange 2021", "3D hiPSC tissue, Fura-2, 37 °C, 1 Hz", "88", "185"],
    ["Saleem 2020 (Hamburg EHT)", "GCaMP6f, 37 °C, 1.5–2 Hz", "~190 (80 % at 116)", "~100"],
    ["Paci 2018 experiments", "single hiPSC-CM, Fluo-4, 37 °C, spontaneous", "270", "245"],
    ["hiPSC-CM models (Paci 2020, Forouzandehmehr 2021, Akwaboah 2021)", "spontaneous / 2.3 Hz", "66–183", "64–169"],
    ["adult human LV", "myocytes / ToR-ORd, 37 °C, 1 Hz", "60", "174"],
  ], { x: M, y: 5.25, w: 8.6, h: 1.65 }, { fontSize: 10, colW: [3.0, 3.2, 1.4, 1.0] });
  callout(s, { x: M + 8.8, y: 5.25, w: W - 2 * M - 8.8, h: 1.65 },
    "Time to peak spans 60–270 ms across sources that all call themselves human (iPSC) cardiac. Panel B: a 20–25 ms change in the assumed Ca peak moves every rate of the fitted human core 2–8× at an unchanged fit quality.");
  S.footer(s, "Stage1b_CaVariability.m: each curve reproduces one source's reported metrics with the same functional form. Indicator caveats: Fluo-4 saturates, GCaMP6f adds 100s of ms.");
  note(s, "This is the single most important slide for the 'ask'. Everything downstream inherits this assumption.");
}

// 5 Realistic fit
{
  const s = pres.addSlide();
  S.header(s, "Step 2 · Humanizing changes only: the realistic fit", "tier 1 · 0.77 %");
  img(s, "t1_isoT.png", { x: M, y: TOP, w: 8.0, h: 4.75 });
  S.card(s, { x: M + 8.2, y: TOP, w: W - 2 * M - 8.2, h: 4.75 }, "What moved, and what it means", [
    { text: "Literature β-MHC core without temperature correction: 10.8 % rms — cannot relax a 37 °C twitch", color: B },
    { text: "Free: isoform ratio on k2, k_2, kah, kamh, ka, kd + a Q10 factor on all cycle rates (k1, k_1 at literature ×0.3)", bold: true },
    "Result 0.77 % rms (fit C with a free Ca transient: 0.81 %)",
    "Identified: the PRODUCT isoform × temperature = ×0.45 (×/÷ 1.12) on the six rates; the split is not (correlation −1.00)",
    "Read as isoform 0.1 × Q10 factor 4.5 → Baker rig at 22–24 °C",
    "Thin filament: k_off 5 s⁻¹, K_coop 20 (bound); Ca peak → force peak 140 ms",
    "Core as skinned muscle: F0 unchanged (74 kPa), ktr 23 s⁻¹ at 37 °C, ATP use halved (8 vs 17 s⁻¹ per head)",
  ], { fontSize: 12, bar: C.navy });
  s.addText("Remaining misfit: the relaxation from 30 % to 5 % of peak (400–500 ms), the same residual fit C had.",
    { x: M, y: 6.05, w: 8.0, h: 0.4, fontFace: F, fontSize: 12, italic: true, color: C.navy, margin: 0 });
  S.footer(s, "robust1_seed base; on ExampleHighATP the same two directions give 1.25 % and the same net ×0.4. Residual still 2.2× the data SEM, lag-1 autocorrelation 0.95.");
}

// 6 Ladder
{
  const s = pres.addSlide();
  S.header(s, "Adding the lesser-known directions one at a time", "stage 2 ladder");
  S.table(s, [
    ["added direction (from the tier-1 optimum)", "rms (%)", "fitted value", "verdict"],
    ["— tier 1 (isoform × temperature)", "0.77", "net ×0.45", { text: "reference", color: C.navy, bold: true }],
    ["Pi release / power stroke (k1, k_1)", "0.77", "×0.25", { text: "no gain: the twitch does not see the Pi step", color: Wn }],
    ["SRX exchange rates, cycle rates held", "0.73", "×0.49", { text: "small gain", color: G }],
    ["series stiffness kSE (compliant posts)", "0.69", "×0.78", { text: "small gain", color: G }],
    ["passive titin k_pas", "0.77", "×1.03", { text: "no gain: invisible in baseline-subtracted force", color: Wn }],
    ["sarcomere length SL0", "0.71", "1.91 µm", { text: "small gain, at the short end", color: G }],
    ["SRX rates with cycle rates free (unbounded)", "0.32", "cycle ×7.5, SRX ×0.2", { text: "leaves the corridor: 24 ATP/head/beat, ktr 86 s⁻¹ — a numerical solution, not a muscle", color: B }],
    ["all directions, inside the literature corridor", "0.38", "every factor on a bound", { text: "optimistic edge: fast corner + SRX ×0.3 + 48 % resting SRX", color: Wn }],
  ], { x: M, y: TOP, w: W - 2 * M, h: 4.3 }, { fontSize: 12, colW: [4.4, 1.0, 2.2, 4.73] });
  callout(s, { x: M, y: 5.7, w: W - 2 * M, h: 1.15 },
    "Within the humanizing corridor the twitch is explained by one quantity: a net ×0.45 slowing of ADP release, hydrolysis, attachment and weak detachment. " +
    "Everything that fits better does so by making the cycle FASTER than the mouse and using a large, slowly exchanging SRX pool as the slow element — the Ca constraint alone cannot exclude it.");
  S.footer(s, "Analyses/HumanTwitch/conclusions.md §2.1; results/overview_tiers.png");
}

// 7 Possible fit
{
  const s = pres.addSlide();
  S.header(s, "Step 3 · Anything but Ca: the possible fit", "stage 3 · 0.30 %");
  img(s, "s3_free.png", { x: M, y: TOP, w: 7.6, h: 4.5 });
  S.table(s, [
    ["", "slow-attachment family (2 starts)", "fast-cycle family (1 start + stage-2 escapes)"],
    ["rms", "0.30–0.32 % (= noise floor 0.35 %)", "0.31–0.32 %"],
    ["core", "attachment ×0.35, series ×0.5, detachment limb ×0.6–0.8", "cycle ×6–7, SRX exchange ×0.25, resting SRX 50 %"],
    ["F0 (saturating)", "19–29 kPa (mouse 71)", "59–65 kPa"],
    ["ktr", "22–23 s⁻¹", "86–99 s⁻¹"],
    ["ATP / head / beat", "0.4–0.9", "17–24"],
  ], { x: M + 7.8, y: TOP, w: W - 2 * M - 7.8, h: 3.3 }, { fontSize: 11, colW: [1.2, 1.9, 1.93] });
  S.card(s, { x: M + 7.8, y: TOP + 3.45, w: W - 2 * M - 7.8, h: 2.1 }, "What every fit shares", [
    "thin filament: k_off 4–6 s⁻¹, K_coop 16–20, Ca50 ≈ 0.6 × peak Ca",
    "activation 77–85 % of saturating force; gate peak ~0.8",
    "Ca-peak → force-peak lag 128–143 ms",
  ], { fontSize: 12, bar: C.maize });
  s.addText("24 free parameters: individually unidentifiable (condition number ~10¹⁸). Shape alone cannot rank two cores that differ 60× in ATP cost.",
    { x: M, y: 5.8, w: 7.6, h: 0.6, fontFace: F, fontSize: 12, italic: true, color: C.navy, margin: 0 });
  S.footer(s, "conclusions.md §3; probes in §3.1 (humanProbe at saturating Ca, SL 2.2 µm, 8 mM MgATP)");
}

// 8 The catch
{
  const s = pres.addSlide();
  S.header(s, "The catch: the core tracks the assumed Ca", "robustness");
  S.table(s, [
    ["Ca transient assumed (fixed)", "TTP / peak→50 % / peak→90 % (ms)", "rms", "identified net factor on the six rates", "ATP / head / beat"],
    ["fast edge of the hiPSC range", "70 / 142 / 352", "0.79 %", { text: "×0.21", bold: true }, "0.9"],
    ["nominal (de Lange 2021)", "93 / 186 / 448", "0.77 %", { text: "×0.45", bold: true }, "1.8"],
    ["slow edge of the hiPSC range", "119 / 225 / 525", "0.86 %", { text: "×3.8 (faster than the mouse)", bold: true, color: B }, "12.3"],
    ["nominal shape, levels 0.25 → 1.0 µM", "93 / 186 / 448", "0.77 %", "×0.45", "1.8"],
  ], { x: M, y: TOP, w: W - 2 * M, h: 2.4 }, { fontSize: 12, colW: [3.2, 3.0, 1.0, 3.3, 1.83] });
  S.card(s, { x: M, y: 3.9, w: 6.0, h: 3.0 }, "What this says", [
    "fit quality does not change; the core that produces it changes 18× across the literature range",
    "absolute Ca levels are irrelevant (only k_on·Ca enters) — a RELATIVE transient is enough",
    "within de Lange's own error bars (TTP ± 5, decay ± 15 ms) the factor is bounded to ~0.3–0.7",
  ], { fontSize: 13, bar: C.navy });
  callout(s, { x: M + 6.2, y: 3.9, w: W - 2 * M - 6.2, h: 3.0 },
    "Consequence: without a Ca recording in the tissues we model, no statement about the human cross-bridge speed from this twitch is defensible. " +
    "A 20 ms change in Ca time to peak is worth a factor 2 on every rate of the 'human' core. This is the first measurement to ask for.", "FDECEA", B);
  S.footer(s, "conclusions.md §2.2 (tier-1 refits t1_caFast / t1_isoT / t1_caSlow / t1_caHiLvl)");
}

// 9 Case B
{
  const s = pres.addSlide();
  S.header(s, "With a measured Ca transient: identified vs. not", "case B");
  img(s, "predict_caseB.png", { x: M, y: TOP, w: 7.9, h: 4.75 });
  S.card(s, { x: M + 8.1, y: TOP, w: W - 2 * M - 8.1, h: 4.75 }, "Profile likelihood with Ca 'measured' (± indicator error)", [
    { text: "Thin filament identified in every fit: k_off 4.2–5.5 s⁻¹, K_coop 18–20", color: G, bold: true },
    { text: "Core speed: one-sided. Steep penalty below ×0.5 of the mouse (2.8 % at ×0.15), flat from ×0.8 to ×10 (0.41 → 0.32 %)", color: Wn, bold: true },
    "a fast cycle is always compensated by slower SRX exchange (×0.56 → ×0.22) and a 50 % resting SRX pool",
    "indicator error (± 5 ms TTP, ± 10 ms decay) moves the core factor 3.6–4.1 only: no longer the limiting uncertainty",
    "→ a known Ca gives the human core a LOWER BOUND, not a value",
  ], { fontSize: 12, bar: C.navy });
  s.addText("Top row: rms when one factor is fixed and everything else refitted. Bottom row: predictions along the core-speed profile (filled: within the 95 % threshold; squares: Ca-error refits).",
    { x: M, y: 6.0, w: 7.9, h: 0.6, fontFace: F, fontSize: 11, italic: true, color: C.mute, margin: 0 });
  S.footer(s, "conclusions.md §5; 12 refits (tB set: net cycle factor, Pi factor, SRX factor, series factor, thin filament, offset)");
}

// 10 Prediction validity
{
  const s = pres.addSlide();
  S.header(s, "Prediction validity with a measured Ca transient", "case B");
  S.table(s, [
    ["prediction (high ATP, 1 Hz)", "spread over fits within 0.42 % rms", "verdict"],
    ["peak force, Ca amplitude ×1.5 / ÷1.5", "+11 … +14 % / −53 … −59 %", { text: "robust", color: G, bold: true }],
    ["peak force, attachment ×0.5 (mavacamten-like)", "−49 … −55 %", { text: "robust", color: G, bold: true }],
    ["peak force at 2 Hz (same Ca shape)", "−2 … −3 %", { text: "robust, but the model has no Ca-handling frequency response", color: G }],
    ["activation (twitch / saturating force)", "77 … 84 %", { text: "robust", color: G }],
    ["F0 at saturating Ca", "50 … 64 kPa", { text: "fair (± 15 %)", color: Wn }],
    ["Brenner ktr", "23 … 97 s⁻¹", { text: "not predicted (4×)", color: B, bold: true }],
    ["ATP use at saturating Ca", "7 … 124 s⁻¹ per head", { text: "not predicted (17×)", color: B, bold: true }],
    ["ATP per head per beat", "1.5 … 24", { text: "not predicted (16×)", color: B, bold: true }],
  ], { x: M, y: TOP, w: W - 2 * M, h: 4.0 }, { fontSize: 12, colW: [4.4, 3.6, 4.33] });
  callout(s, { x: M, y: 5.45, w: W - 2 * M, h: 1.4 },
    "With the Ca transient known, twitch-level predictions (response to Ca amplitude, rate, attachment inhibitors) are trustworthy: they follow from the thin filament and the activation level, which the twitch identifies. " +
    "The muscle-level kinetics and energetics (ktr, ATPase, ATP per beat) are not: they follow the unidentified core speed. The fitted rates do NOT report how the muscle's inner properties changed until one muscle-level quantity is measured too.");
  S.footer(s, "PredictSet.m → results/predict_caseB.csv");
}

// 11 Ascertain
{
  const s = pres.addSlide();
  S.header(s, "Facts to ascertain from the experimenters", "ask · facts");
  S.table(s, [
    ["item", "why it matters to the model", "status"],
    ["Plate 6 platemap: what is in columns 1–3, is row D of 4–6 DMSO, DMSO %, dosing time", "the 15 'control' wells were inferred from an isoproterenol-like signature, not from a map", { text: "inferred", color: Wn }],
    ["cell line, differentiation / tissue age, CM purity, fibroblast content", "maturity sets β-MHC, ssTnI, t-tubules → Ca shape and core speed", { text: "unknown", color: B }],
    ["recording temperature, [Ca]o and medium, time since medium change", "a 5 °C error is ×1.8 on every rate; [Ca]o sets the transient amplitude", { text: "presumed 37 °C", color: Wn }],
    ["post stiffness (µN/µm) for 'Skeletal Muscle 12×', how force is computed from magnet displacement", "the twitch is auxotonic; every good fit wants a softer series element (×0.5–0.8)", { text: "unknown", color: B }],
    ["stimulus–force synchronisation; pulse parameters", "the offset is fitted (±4 ms); a synchronised record would remove one free parameter", { text: "not synchronised", color: Wn }],
    ["raw Pulse3D / h5 exports with the absolute baseline (not baseline-subtracted)", "diastolic force constrains diastolic activation and SRX; currently invisible", { text: "export only", color: Wn }],
    ["pre-dose baseline recordings of all plates", "many more controls, paired pre/post for each drug", { text: "implied to exist", color: Wn }],
    ["tissue geometry (cross-section, length), sarcomere length if ever imaged", "turns σ (µN per kPa) into a force per area; SL fixes overlap", { text: "unknown", color: B }],
  ], { x: M, y: TOP, w: W - 2 * M, h: 5.4 }, { fontSize: 11, colW: [4.6, 5.6, 2.13] });
  S.footer(s, "From EHT_identification_report.md §7.2, updated with this campaign. All of these are facts the lab most likely already has.");
}

// 12 Need
{
  const s = pres.addSlide();
  S.header(s, "What we need measured, ranked by what it buys", "ask · measurements");
  S.table(s, [
    ["#", "measurement", "what it removes", "model prediction to test"],
    ["1", { text: "relative Ca transient in the same wells at 1 Hz (low-affinity or ratiometric dye; not GCaMP)", bold: true }, "the Ca-shape assumption (×0.2 … ×4 on the core); fixes the thin filament; no calibration needed (only k_on·Ca enters)", "Ca peaks 128–143 ms before force; Ca at 50 % just after the force peak"],
    ["2", { text: "paced O₂ consumption or ATPase of the same tissues (also with mavacamten)", bold: true }, "the flat core-speed profile: equally good fits differ 16× in ATP per beat", "1.5 (slow core) … 24 (fast core) ATP per head per beat"],
    ["3", { text: "force per cross-section, or skinned-EHT maximal force", bold: true }, "the σ degeneracy; separates the 20–30 kPa cores from the 50–65 kPa ones", "activation 77–85 % → twitch stress ≈ 0.8 × F0"],
    ["4", "ktr or kACT on skinned EHT at 37 °C (or myofibrils)", "core speed directly", "23 s⁻¹ (realistic core) vs 86–99 (fast family)"],
    ["5", "force–pCa on skinned EHT from the same line", "thin-filament absolute sensitivity and steepness; tests the K_coop = 20 at bound", "pCa50 ≈ 6.2 relative to peak Ca, n_H ≥ 4.6"],
  ], { x: M, y: TOP, w: W - 2 * M, h: 4.6 }, { fontSize: 11, colW: [0.4, 4.0, 4.3, 3.63] });
  callout(s, { x: M, y: 6.0, w: W - 2 * M, h: 0.85 },
    "Order matters: 1 first (removes the assumption that even destroys the lower bound), then 2 or 3 (collapses the profile). Each buys a specific non-identifiability; none of the others can substitute for 1.");
  S.footer(s, "conclusions.md §4, §5.2");
}

// 13 Nice to have
{
  const s = pres.addSlide();
  S.header(s, "What else we can ask for, if available or possible", "ask · optional");
  S.card(s, { x: M, y: TOP, w: 6.0, h: 2.75 }, "Without a dye (protocol changes only)", [
    "[Ca]o series (0.5 – 3.6 mM): fit predicts +11…+14 % force for a 1.5× transient, −53…−59 % for ÷1.5 → tests the near-saturation operating point",
    "pacing-frequency series (0.5–2 Hz) and post-rest potentiation: constrains Ca-handling time constants the model lacks",
    "a temperature step (e.g. 32 vs 37 °C): Q10 of the human core directly",
  ], { fontSize: 12, bar: C.navy });
  S.card(s, { x: M + 6.2, y: TOP, w: W - 2 * M - 6.2, h: 2.75 }, "Drug waveforms (Prism % changes exist)", [
    "mavacamten: full waveforms, ideally with Ca → acts on the cross-bridge side; fit predicts −49…−55 % for ka ×0.5",
    "isoproterenol: Ca side + thin-filament desensitisation (ssTnI blunts the PKA effect)",
    "diltiazem: pure Ca-amplitude lever, a second [Ca]o-type test",
  ], { fontSize: 12, bar: C.maize });
  S.card(s, { x: M, y: TOP + 2.95, w: 6.0, h: 2.6 }, "If the dye work is possible", [
    "Ca calibration (Fura-2 with ionomycin Rmin/Rmax): breaks k_on·Ca, makes pCa50 absolute",
    "Ca and force recorded simultaneously in the same beat (removes the synchronisation offset)",
    "note indicator kinetics: Fluo-4 compresses peaks, GCaMP6f adds 100s of ms",
  ], { fontSize: 12, bar: C.navy });
  S.card(s, { x: M + 6.2, y: TOP + 2.95, w: W - 2 * M - 6.2, h: 2.6 }, "Structural / biochemical", [
    "mant-ATP SRX fraction in these cardiomyocytes: the fast family needs ~50 % resting SRX, the realistic core 33 %",
    "MHC and TnI isoform fractions of the line (β-MHC, ssTnI vs cTnI) and titin isoform",
    "sarcomere length imaging at rest",
  ], { fontSize: 12, bar: C.maize });
  S.footer(s, "EHT_identification_report.md §7; conclusions.md §4");
}

// 14 Summary
{
  const s = pres.addSlide();
  S.header(s, "Summary", "take-home");
  S.stat(s, M, TOP, "0.77 %", "realistic fit: mouse core slowed ×0.45 on the isoform rates, Ca fixed from literature", C.navy);
  S.stat(s, M + 3.3, TOP, "0.30 %", "possible fit: noise floor, reached by two cores that differ 60× in ATP cost", C.navy);
  S.stat(s, M + 6.6, TOP, "60–270 ms", "literature spread of the hiPSC Ca time to peak; 20 ms = factor 2 on every rate", B);
  S.stat(s, M + 9.9, TOP, "×0.5 … ×10", "core speed allowed even with Ca known: lower bound only", Wn);
  S.table(s, [
    ["we can already say", "we cannot say yet", "what unlocks it"],
    ["the thin filament: k_off ≈ 5 s⁻¹, K_coop ≈ 20, Ca50 ≈ 0.6 × peak Ca; activation ~80 % of saturating; Ca → force lag ~135 ms", "the human cross-bridge speed (isoform × temperature), hence ktr", "Ca transient (relative), then ktr/kACT or ATPase"],
    ["how peak force responds to Ca amplitude, rate and an attachment inhibitor (± a few %)", "ATP cost per beat (factor 16–60 open)", "O₂ consumption / ATPase of paced tissues"],
    ["that a β-MHC core at rig temperature cannot produce a 37 °C twitch", "force per area, resting SRX, SL of the EHT", "cross-section or skinned F0; mant-ATP; SL imaging"],
  ], { x: M, y: 3.15, w: W - 2 * M, h: 2.6 }, { fontSize: 12, colW: [4.6, 3.8, 3.93] });
  callout(s, { x: M, y: 6.0, w: W - 2 * M, h: 0.85 },
    "Bottom line: the twitch shape is reproduced to noise, but the fitted rates are not yet a statement about the human muscle. One Ca recording and one energetic or kinetic anchor would make them one.");
  S.footer(s, "Worktree C:\\home\\git\\ATP-wt-human (branch HumanTwitch); Analyses/HumanTwitch/conclusions.md and results/");
}

// 15 Backup: method
{
  const s = pres.addSlide();
  S.header(s, "Backup · how the fits were done", "methods");
  S.card(s, { x: M, y: TOP, w: 6.0, h: 5.6 }, "Pipeline", [
    "target: beat-averaged active force of 15 control wells (0–950 ms, 96 points); residual (σ·F_model − F_data)/F_peak, σ in closed form",
    "model: robust1_seed core + twoState thin filament; Ca kind 'risedecay' fixed; isometric, SL 2.0 µm; 2 s diastole + 2 beats per evaluation (~6 s)",
    "stage 2: group multipliers along the HumanizeXB directions (isoform, Q10, Pi, SRX, series, titin, SL) + thin filament + offset; lsqnonlin, resumable, warm-started ladder",
    "stage 3: 24 individual multipliers, 3 starts",
    "identifiability: forward-difference Jacobian; profile likelihood (parameter fixed on a grid, rest refitted) for the case-B set",
    "predictions: humanProbe at saturating Ca (F0, SRX, ATP use, Brenner ktr) + twitch responses (Ca ×1.5, ÷1.5, 2 Hz, ka ×0.5)",
  ], { fontSize: 12, bar: C.navy });
  img(s, "overview_tiers.png", { x: M + 6.2, y: TOP, w: W - 2 * M - 6.2, h: 5.6 });
  S.footer(s, "17 + 12 converged fits, logs and figures in Analyses/HumanTwitch/results/; parameter snapshots params/human_twitch_t1_isoT.m (realistic), params/human_twitch_s3_free.m (possible)");
}

pres.writeFile({ fileName: path.resolve(__dirname, "..", "2026-10-HumanTwitch.pptx") }).then((f) => console.log("wrote " + f));
