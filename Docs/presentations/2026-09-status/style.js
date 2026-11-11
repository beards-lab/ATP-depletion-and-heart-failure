// Shared style for the 2026-09 project-status deck (EMC-poster look: UMich navy + maize, Arial).
// Every part (part1.js, part2.js, part3.js) exports addSlides(pres, S) and uses ONLY these helpers,
// so the merged deck (build_all.js) looks like one deck.
//
// Build one part standalone:  NODE_PATH=~/tools/pptx/node_modules node part1.js   -> part1.pptx
// Build the merged deck:      NODE_PATH=~/tools/pptx/node_modules node build_all.js
const path = require("path");
const fs = require("fs");

const C = {
  navy: "00274C",   // UMich blue: header bars, titles
  maize: "FFCB05",  // UMich maize: accent rule, highlights
  ink: "1F2937",    // body text
  mute: "6B7280",   // footers, captions
  light: "EEF2F7",  // card background
  good: "2E7D32",   // "fits"
  warn: "B26A00",   // "partial"
  bad: "B3261E",    // "misses"
  white: "FFFFFF",
};
const F = "Arial";
const W = 13.333, H = 7.5;            // LAYOUT_WIDE
const M = 0.5;                         // side margin
const TOP = 1.2;                       // content starts below the header bar
const BOTTOM = 6.95;                   // content ends above the footer

const ROOT = path.resolve(__dirname, "..", "..", "..");      // repo root
const asset = (rel) => path.resolve(ROOT, rel);             // repo-relative path -> absolute

// Header bar + title (+ optional small kicker on the right)
function header(s, title, kicker) {
  s.background = { color: C.white };
  s.addShape("rect", { x: 0, y: 0, w: W, h: 0.95, fill: { color: C.navy }, line: { color: C.navy } });
  s.addShape("rect", { x: 0, y: 0.95, w: W, h: 0.06, fill: { color: C.maize }, line: { color: C.maize } });
  s.addText(title, { x: M, y: 0.08, w: W - 2 * M - (kicker ? 3 : 0), h: 0.8, fontFace: F, fontSize: 28, bold: true, color: C.white, valign: "middle", margin: 0 });
  if (kicker) s.addText(kicker, { x: W - M - 3, y: 0.08, w: 3, h: 0.8, fontFace: F, fontSize: 13, color: C.maize, align: "right", valign: "middle", margin: 0 });
}

// Footer: source / caveat line (small, grey) + slide tag
function footer(s, text) {
  if (text) s.addText(text, { x: M, y: 7.0, w: W - 2 * M, h: 0.4, fontFace: F, fontSize: 10, color: C.mute, valign: "middle", margin: 0 });
}

// Bullets. items: array of strings or {text, level (0/1), bold, color}
function bullets(s, items, box, opts = {}) {
  const runs = items.map((it, i) => {
    const o = typeof it === "string" ? { text: it } : it;
    return {
      text: o.text,
      options: {
        bullet: o.level ? { indent: 18 } : true, indentLevel: o.level || 0,
        bold: !!o.bold, color: o.color || C.ink, breakLine: i < items.length - 1,
        fontSize: o.level ? (opts.fontSize || 18) - 2 : (opts.fontSize || 18),
      },
    };
  });
  s.addText(runs, { ...box, fontFace: F, valign: "top", paraSpaceAfter: opts.space ?? 6, margin: 0.05 });
}

// Image fitted inside a box, keeping aspect ratio. Needs pixel size (w,h) of the file.
function image(s, file, box, px) {
  const abs = path.isAbsolute(file) ? file : asset(file);
  if (!fs.existsSync(abs)) return placeholder(s, box, "missing figure: " + file);
  let { x, y, w, h } = box;
  if (px) {
    const r = px[0] / px[1];
    if (w / h > r) { const nw = h * r; x += (w - nw) / 2; w = nw; } else { const nh = w / r; y += (h - nh) / 2; h = nh; }
  }
  s.addImage({ path: abs, x, y, w, h });
}

// Clearly marked placeholder for a figure/result that is still being computed
function placeholder(s, box, label) {
  s.addShape("rect", { ...box, fill: { color: "FFF8DC" }, line: { color: C.warn, dashType: "dash", width: 1.5 } });
  s.addText("PLACEHOLDER — " + label, { ...box, fontFace: F, fontSize: 16, bold: true, color: C.warn, align: "center", valign: "middle" });
}

// Card: light box with a bold heading and bullet list
function card(s, box, heading, items, opts = {}) {
  s.addShape("rect", { ...box, fill: { color: opts.fill || C.light }, line: { color: opts.line || C.light } });
  if (opts.bar) s.addShape("rect", { x: box.x, y: box.y, w: 0.08, h: box.h, fill: { color: opts.bar }, line: { color: opts.bar } });
  s.addText(heading, { x: box.x + 0.2, y: box.y + 0.1, w: box.w - 0.3, h: 0.45, fontFace: F, fontSize: opts.headSize || 17, bold: true, color: opts.headColor || C.navy, margin: 0 });
  bullets(s, items, { x: box.x + 0.2, y: box.y + 0.6, w: box.w - 0.35, h: box.h - 0.7 }, { fontSize: opts.fontSize || 14, space: 4 });
}

// Table with a navy header row. rows: array of arrays of strings (first row = header)
function table(s, rows, box, opts = {}) {
  const fs_ = opts.fontSize || 12;
  const data = rows.map((r, i) => r.map((c) => {
    const cell = typeof c === "string" ? { text: c } : c;
    return {
      text: cell.text,
      options: i === 0
        ? { bold: true, color: C.white, fill: { color: C.navy }, fontSize: fs_, fontFace: F }
        : { color: cell.color || C.ink, bold: !!cell.bold, fill: { color: cell.fill || (i % 2 ? C.white : C.light) }, fontSize: fs_, fontFace: F },
    };
  }));
  s.addTable(data, { ...box, colW: opts.colW, border: { type: "solid", color: "D1D5DB", pt: 0.5 }, valign: "middle", margin: 0.05, autoPage: false });
}

// Big number callout
function stat(s, x, y, num, label, color) {
  s.addText(num, { x, y, w: 3, h: 0.8, fontFace: F, fontSize: 36, bold: true, color: color || C.navy, margin: 0 });
  s.addText(label, { x, y: y + 0.8, w: 3, h: 0.6, fontFace: F, fontSize: 13, color: C.ink, margin: 0, valign: "top" });
}

// Standalone build helper for a part
function buildStandalone(addSlides, outName) {
  const pptxgen = require("pptxgenjs");
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE";
  pres.title = "ATP-dependent cross-bridge model: where we stand";
  addSlides(pres, module.exports);
  return pres.writeFile({ fileName: path.join(__dirname, outName) }).then((f) => console.log("wrote " + f));
}

module.exports = { C, F, W, H, M, TOP, BOTTOM, ROOT, asset, header, footer, bullets, image, placeholder, card, table, stat, buildStandalone };
