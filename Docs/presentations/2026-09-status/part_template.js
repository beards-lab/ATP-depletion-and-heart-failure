// TEMPLATE for a part module. Copy to partN.js.
const S = require("./style");
function addSlides(pres, S) {
  const s = pres.addSlide();
  S.header(s, "Slide title", "kicker");
  S.bullets(s, ["point one", { text: "sub-point", level: 1 }], { x: S.M, y: S.TOP, w: 6, h: 5 });
  S.footer(s, "Source: path/to/source");
  s.addNotes("Speaker notes.");
}
module.exports = { addSlides };
if (require.main === module) S.buildStandalone(addSlides, "partN.pptx");
