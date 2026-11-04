// Two slides: the EMC-poster figure set redrawn for the POOLED 8 + 2 mM fit.
// Figures: Analyses/PooledFitFigures (MakePooledFitFigures.m). Style: ../2026-09-status/style.js.
//
// Build (WSL):  NODE_PATH=~/tools/pptx/node_modules node build.js   -> 2026-10-PooledFit.pptx
const path = require("path");
const pptxgen = require("pptxgenjs");
const S = require("../2026-09-status/style");
const { C, F, M, TOP } = S;

const FIG = "Analyses/PooledFitFigures/results/";
const PX = {
  slack_all: [2787, 887], slack_cycle4: [2787, 887], ktr_protocol: [2789, 891], passive: [2052, 748],
  stairs: [4080, 866], fv_force: [1364, 885], fv_power: [1368, 886], feat: [1365, 886],
};
const img = (s, f, box, px) => S.image(s, FIG + f + ".png", box, px);
const MODELS = "Model: 8 mM robust1_seed (= refit9_opt); 2 mM fit A = same core, R2 ×0.54 (k2 power gate). " +
               "Data: pooled mean of 3 preparations (03/27, 04/03, 04/10); band = between-preparation SD.";

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "Pooled 8 + 2 mM fit";

// ---------------------------------------------------------------- slide 1: time courses
{
  const s = pres.addSlide();
  S.header(s, "Pooled 8 vs 2 mM fit: slack–restretch and F–V", "fitted · n = 3");
  const LW = 8.5, RX = M + LW + 0.25, RW = 13.333 - M - RX;
  img(s, "slack_all", { x: M, y: TOP, w: LW, h: LW / 3.142 }, PX.slack_all);
  img(s, "slack_cycle4", { x: M, y: TOP + 2.85, w: LW, h: LW / 3.142 }, PX.slack_cycle4);
  img(s, "fv_force", { x: RX, y: TOP, w: RW, h: 2.6 }, PX.fv_force);
  img(s, "fv_power", { x: RX, y: TOP + 2.75, w: RW, h: 2.6 }, PX.fv_power);
  s.addText("F–V band = Baker 2022 vs AB2021 (2 preps); ◆ = pooled in-prep FV2 at 2 ML/s (n = 2)",
            { x: RX, y: TOP + 5.38, w: RW, h: 0.35, fontFace: F, fontSize: 10, italic: true, color: C.mute, margin: 0, isTextBox: true });
  S.footer(s, MODELS);
  s.addNotes(
    "Same panels as the EMC poster, but the data line is now the pooled mean the model is fitted to, not the 03/27 preparation. " +
    "Bands are the between-preparation SD on the pooled kPa scale (per-prep Fss22 anchor + rundown lift, as in the objective).\n\n" +
    "Low ATP raises force about 1.3x and slows redevelopment; one ATP lever (R2 x0.54, ATP binding in A2 detachment) carries both. " +
    "At 2 mM the model's re-stretch transient is too sharp: peak then a deep valley, where the pooled data has a broad hump. " +
    "F-V: 8 mM is scored; 2 mM F-V is NOT scored, so the red F-V curve is a prediction that lands on the data.");
}

// ---------------------------------------------------------------- slide 2: features + held-out
{
  const s = pres.addSlide();
  S.header(s, "Slack features, and what was never fitted", "features · predictions");
  const LW = 8.5, RX = M + LW + 0.25, RW = 13.333 - M - RX;
  const FEATS = ["steady", "peak1_y", "t0", "A", "ktrfit", "rsK80"];
  const cw = LW / 3, ch = cw / 1.54 + 0.06;
  FEATS.forEach((f, i) => img(s, "feat_" + f, { x: M + (i % 3) * cw, y: TOP + Math.floor(i / 3) * ch, w: cw, h: ch - 0.06 }, PX.feat));
  const sy = TOP + 2 * ch + 0.12;
  img(s, "stairs", { x: M, y: sy, w: LW, h: LW / 4.711 }, PX.stairs);

  img(s, "ktr_protocol", { x: RX, y: TOP, w: RW, h: RW / 3.13 }, PX.ktr_protocol);
  img(s, "passive", { x: RX, y: TOP + 1.4, w: RW, h: RW / 2.743 }, PX.passive);
  S.card(s, { x: RX, y: TOP + 2.95, w: RW, h: 2.75 }, "Reading", [
    { text: "Fits: force level ×1.3, t₀ slopes, redeveloped force, 8 mM ktr and rsK80, ktr protocol at both [ATP]", color: C.good },
    { text: "2 mM misses: slack ktr falls over the cycles in data, flat in model; rsK80 ~20 % high", color: C.bad },
    { text: "Staircase (unscored): 8 mM plateau 1.33 vs 1.20–1.35 in data; 2 mM is a pure prediction", color: C.warn },
  ], { fontSize: 13, headSize: 16 });
  S.footer(s, MODELS + " Staircase not poolable (03/27: 8×0.1 s steps; others 4×0.2 s): one tile per prep. Passive: 03/27 only.");
  s.addNotes(
    "Six features from the poster (a-f), now pooled. Band = SD across the three preparations of the per-prep values the target is built from; " +
    "force features are normalised per preparation, so the band is the SHAPE spread, not the 15 % between-heart level spread.\n\n" +
    "Six more feature panels (valley, second peak, post-peak drop, recovery time, valley time, peak time) are in " +
    "Analyses/PooledFitFigures/results/feat_*.png: the 2 mM valley and second peak collapse in the model (post-peak drop -31 vs -14.5 kPa), " +
    "and the 2 mM recovery time stays flat where the data rise.\n\n" +
    "The k_tr protocol is held-out validation (unscored on purpose). The staircase was never fitted and has no 2 mM data.");
}

// ---------------------------------------------------------------- slide 3: the re-stretch transient misses
{
  const s = pres.addSlide();
  S.header(s, "Where 2 mM misses: the re-stretch transient", "features · misses");
  const LW = 8.7, RX = M + LW + 0.25, RW = 13.333 - M - RX;
  const FEATS = ["vall_y", "peak2", "vall2_dy", "rsTfitE", "vall_tq", "peak1_tq"];
  const cw = LW / 3, ch = cw / 1.54 + 0.12;
  FEATS.forEach((f, i) => img(s, "feat_" + f, { x: M + (i % 3) * cw, y: TOP + 0.6 + Math.floor(i / 3) * ch, w: cw, h: ch - 0.12 }, PX.feat));
  // these panels carry no legend of their own: one key under the grid
  const key = (t, col, last) => ({ text: t, options: { color: col, bold: true, breakLine: false } });
  s.addText([
    key("●— data 8 mM    ", "0A0D12"), key("━ model 8 mM    ", "0059B3"),
    key("●— data 2 mM    ", "707378"), key("━ model 2 mM    ", "CC331A"),
    { text: "band = SD, n = 3", options: { color: C.mute } },
  ], { x: M, y: TOP + 0.6 + 2 * ch + 0.05, w: LW, h: 0.4, fontFace: F, fontSize: 13, margin: 0, isTextBox: true });
  S.card(s, { x: RX, y: TOP + 0.1, w: RW, h: 5.55 }, "2 mM: re-stretch collapses", [
    { text: "Valley after the peak: data 87 → 72 kPa, model ~55 (≈35 % low)", color: C.bad },
    { text: "Second peak: data 93–97 kPa, model 64–80", color: C.bad },
    { text: "Post-peak drop: model −31 vs data −14 kPa (too deep)", color: C.bad },
    { text: "Recovery time: data 24–37 ms, model ~12 ms (no ATP effect at all)", color: C.bad },
    { text: "Valley and peak time: model late; the valley-time ATP effect has the wrong sign", color: C.bad },
    { text: "8 mM is inside the band for the force features; valley time ~3 ms late", color: C.good },
    { text: "One gate (R2 ×0.54) slows redevelopment but cannot hold force up after the re-stretch", color: C.ink },
  ], { fontSize: 13, headSize: 16 });
  S.footer(s, MODELS + " Remaining 6 slack features (slide 2 has the first 6).");
  s.addNotes(
    "These are the feature panels that were not on the poster. 8 mM sits in or near the band; the misses are the 2 mM re-stretch transient. " +
    "In the data, low ATP gives a broad, high hump after the re-stretch with a slow recovery. The model gives a sharp peak, a deep valley " +
    "and a recovery as fast as at 8 mM. See slide 1, cycle 4, for the trace. The extended sets B (+ ksr0) and C (+ kmsrd) in " +
    "Analyses/LowATP_Refit9K2Sweep do not fix this: the post-peak drop is the same in A, B and C.");
}

const out = path.join(__dirname, "2026-10-PooledFit.pptx");
pres.writeFile({ fileName: out }).then((f) => console.log("wrote " + f));
