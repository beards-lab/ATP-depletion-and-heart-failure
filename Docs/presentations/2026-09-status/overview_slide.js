// ONE-SLIDE project overview (model schematic + 8 vs 2 mM slack fit + power–velocity blow-up).
// Build: NODE_PATH=~/tools/pptx/node_modules node overview_slide.js  -> ../2026-09-ProjectOverview-1slide.pptx
const path = require("path");
const S = require("./style");
const { C, F, W, M } = S;
const POSTER = "Docs/presentations/2026-09-status/assets/poster/";

function arrow(s, x1, y1, x2, y2, o = {}) {
  s.addShape("line", {
    x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1) || 0.0001, h: Math.abs(y2 - y1) || 0.0001,
    flipH: x2 < x1, flipV: y2 < y1,
    line: { color: o.color || C.ink, width: o.width || 1.5, endArrowType: "triangle", dashType: o.dash || "solid" },
  });
}
function label(s, text, x, y, w, h, o = {}) {
  s.addText(text, { x, y, w, h, fontFace: F, fontSize: o.fontSize || 10, color: o.color || C.ink, align: o.align || "center",
                    valign: "middle", margin: 0, bold: !!o.bold, italic: !!o.italic });
}
function stateBox(s, x, y, w, h, name, sub, fill) {
  s.addShape("roundRect", { x, y, w, h, rectRadius: 0.07, fill: { color: fill }, line: { color: C.navy, width: 1.25 } });
  s.addText([{ text: name, options: { bold: true, fontSize: 15, breakLine: true } }, { text: sub, options: { fontSize: 9 } }],
            { x, y, w, h, fontFace: F, color: C.navy, align: "center", valign: "middle", margin: 0 });
}
function sectionTitle(s, text, x, y, w) {
  s.addText(text, { x, y, w, h: 0.32, fontFace: F, fontSize: 15, bold: true, color: C.navy, margin: 0 });
}

function addSlide(pres) {
  const s = pres.addSlide();
  S.header(s, "One cross-bridge model across ATP levels", "ATP depletion & heart failure");

  // ───────── left: model schematic ─────────
  sectionTitle(s, "Model: six-state cross-bridge cycle", M, 1.15, 4.3);
  const bw = 1.3, bh = 0.62, xL = 0.8, xR = 3.25;
  const y1 = 1.65, y2 = 2.85, y3 = 4.05, d = 0.1;
  const cxL = xL + bw / 2, cxR = xR + bw / 2, gap = xR - xL - bw;
  stateBox(s, xL, y1, bw, bh, "SR", "SRX · M·ATP", "E6E1F5");
  stateBox(s, xR, y1, bw, bh, "SRD", "SRX · M·ADP·Pi", "E6E1F5");
  stateBox(s, xL, y2, bw, bh, "PT", "DRX · M·ATP", C.light);
  stateBox(s, xR, y2, bw, bh, "PD", "DRX · M·ADP·Pi", C.light);
  stateBox(s, xR, y3, bw, bh, "A1", "weak, pre-stroke", "FFF1C2");
  stateBox(s, xL, y3, bw, bh, "A2", "strong, post-stroke", "FFE08A");
  // horizontal pairs
  for (const [y, top, bot] of [[y1, "", ""], [y2, "hydrolysis", ""], [y3, "power stroke", ""]]) {
    arrow(s, xL + bw, y + bh / 2 - d, xR, y + bh / 2 - d);
    arrow(s, xR, y + bh / 2 + d, xL + bw, y + bh / 2 + d);
    if (top) label(s, top, xL + bw, y - 0.02, gap, 0.22, { fontSize: 9 });
  }
  // PD <-> A1 attach / detach
  arrow(s, cxR - d, y2 + bh, cxR - d, y3);
  arrow(s, cxR + d, y3, cxR + d, y2 + bh);
  label(s, "attach", cxR + 0.2, y2 + bh + 0.14, 0.9, 0.26, { align: "left", fontSize: 9 });
  // SRD <-> PD, force-dependent
  arrow(s, cxR - d, y1 + bh, cxR - d, y2);
  arrow(s, cxR + d, y2, cxR + d, y1 + bh);
  label(s, "force-dep.", cxR + 0.2, y1 + bh + 0.14, 0.9, 0.26, { align: "left", fontSize: 9 });
  // PT -> SR, suppressed by force
  arrow(s, cxL - d, y2, cxL - d, y1 + bh);
  label(s, "SRX entry ↓ by force", cxL + 0.05, y1 + bh + 0.14, 1.2, 0.26, { align: "left", fontSize: 9 });
  // A2 -> PT : the ATP-gated step
  arrow(s, cxL, y3, cxL, y2 + bh, { color: C.bad, width: 2.5 });
  label(s, "k2 · [ATP]ᵖ", cxL + 0.12, y2 + bh + 0.1, 1.2, 0.24, { align: "left", bold: true, color: C.bad, fontSize: 11 });
  label(s, "ADP release, ATP binding, detach", cxL + 0.12, y2 + bh + 0.33, 2.2, 0.2, { align: "left", color: C.bad, fontSize: 8 });

  // mechanical arrangement + one-liners
  S.image(s, POSTER + "image15.png", { x: M + 0.2, y: 4.85, w: 2.1, h: 0.9 }, [600, 250]);
  s.addText([
    { text: "A1/A2 are distributions over XB strain (±60 nm)", options: { breakLine: true } },
    { text: "XB ∥ titin, in series with an elastic element", options: { breakLine: true } },
    { text: "Rates: strain- and force-dependent (PCHIP knots)" },
  ], { x: 2.85, y: 4.8, w: 1.95, h: 1.0, fontFace: F, fontSize: 9, color: C.ink, valign: "middle", margin: 0, paraSpaceAfter: 3 });

  // take-home box under the schematic
  s.addShape("rect", { x: M, y: 5.95, w: 4.3, h: 0.95, fill: { color: C.light }, line: { color: C.light } });
  s.addShape("rect", { x: M, y: 5.95, w: 0.07, h: 0.95, fill: { color: C.maize }, line: { color: C.maize } });
  s.addText([
    { text: "Fit at 8 mM; 2 mM enters only via ATP-gated rates", options: { bold: true, color: C.navy, breakLine: true } },
    { text: "Same parameters separate 8 vs 2 mM in force–velocity, power and the slack–restretch ladder" },
  ], { x: M + 0.2, y: 5.98, w: 4.05, h: 0.9, fontFace: F, fontSize: 11, color: C.ink, valign: "middle", margin: 0, paraSpaceAfter: 3 });

  // divider
  s.addShape("line", { x: 5.05, y: 1.2, w: 0, h: 5.7, line: { color: "D1D5DB", width: 1 } });

  // ───────── right top: power–velocity blow-up ─────────
  const RX = 5.3, RW = W - M - RX;          // 7.53 in
  const pwH = 2.95, pwW = pwH * 1368 / 886;  // ≈ 4.94 in
  S.image(s, POSTER + "image19.png", { x: RX, y: 1.15, w: pwW, h: pwH }, [1368, 886]);

  // legend + numbers next to the power plot
  const lx = RX + pwW + 0.25, lw = W - M - lx;
  sectionTitle(s, "Power at low ATP", lx, 1.2, lw);
  const leg = [["dot", "000000", "data 8 mM"], ["dot", "808080", "data 2 mM"], ["line", "1F5FBF", "model 8 mM"], ["line", "D03A1A", "model 2 mM"]];
  leg.forEach(([kind, col, lab], i) => {
    const y = 1.58 + i * 0.27;
    if (kind === "dot") s.addShape("ellipse", { x: lx + 0.1, y: y + 0.06, w: 0.15, h: 0.15, fill: { color: col }, line: { color: col } });
    else s.addShape("line", { x: lx, y: y + 0.135, w: 0.35, h: 0, line: { color: col, width: 3 } });
    label(s, lab, lx + 0.45, y, lw - 0.45, 0.27, { align: "left", fontSize: 11, color: C.ink });
  });
  s.addText([
    { text: "−28 %", options: { fontSize: 26, bold: true, color: C.bad, breakLine: true } },
    { text: "peak power at 2 mM vs 8 mM (data 0.65 → 0.47 s⁻¹)", options: { fontSize: 11, color: C.ink, breakLine: true } },
    { text: " ", options: { fontSize: 6, breakLine: true } },
    { text: "Model follows the loss across the velocity range", options: { fontSize: 11, color: C.ink } },
  ], { x: lx, y: 2.72, w: lw, h: 1.4, fontFace: F, valign: "top", margin: 0 });

  // ───────── right bottom: slack ladder ─────────
  const slH = RW * 887 / 2787;              // ≈ 2.40 in
  S.image(s, POSTER + "image20.png", { x: RX, y: 6.95 - slH, w: RW, h: slH }, [2787, 887]);
  label(s, "Slack–restretch ladder, 8 vs 2 mM: higher force but slower redevelopment at low ATP",
        RX, 4.2, RW, 0.26, { fontSize: 11, italic: true, color: C.mute, align: "left" });

  S.footer(s, "Fit: EMC-poster model (spreadw core + ATP power gates on k2 and hydrolysis), one prep (03/27), mouse myocardium; bands = SD of 3 preps. Params atphyp_H1_opt.");
  s.addNotes(
    "One-slide overview. Left: the six-state model. Detached DRX heads PT (M·ATP) and PD (M·ADP·Pi), attached A1 (weak) and A2 (strong), " +
    "and two super-relaxed pools SR and SRD. Attached states are resolved over strain. SRX entry is suppressed by force. " +
    "The red A2→PT step lumps ADP release, ATP binding and detachment, and is where ATP acts (k2 power gate, exponent 0.44). " +
    "The poster fit also gated hydrolysis kah (exponent 0.30), which we no longer consider defensible chemistry. " +
    "Right top: power–velocity. At 2 mM peak power falls by about 28 % and the model follows the data across velocities. " +
    "Right bottom: the five-cycle slack–restretch ladder. At 2 mM steady force is higher and redevelopment slower; the model captures both. " +
    "Caveats: one preparation (03/27); the post-restretch 80 ms window and the 2 mM ktr drop are under-reproduced. " +
    "The current pooled 8 mM incumbent (refit9) has no ATP gates yet; the low-ATP refit is in progress.");
}

const pptxgen = require("pptxgenjs");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "ATP-dependent cross-bridge model: overview";
addSlide(pres);
pres.writeFile({ fileName: path.resolve(__dirname, "..", "2026-09-ProjectOverview-1slide.pptx") }).then((f) => console.log("wrote " + f));
module.exports = { addSlide };
