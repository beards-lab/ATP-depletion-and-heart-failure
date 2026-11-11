// Merge all parts into the final deck. Run: NODE_PATH=~/tools/pptx/node_modules node build_all.js
const path = require("path");
const pptxgen = require("pptxgenjs");
const S = require("./style");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "ATP-dependent cross-bridge model: where we stand";
for (const p of ["./part1", "./part2", "./part3"]) require(p).addSlides(pres, S);
require("./part2").addBackup(pres, S);   // backups last
require("./part4").addFeatureTables(pres, S);
pres.writeFile({ fileName: path.resolve(__dirname, "..", "2026-09-ProjectStatus.pptx") }).then((f) => console.log("wrote " + f));
