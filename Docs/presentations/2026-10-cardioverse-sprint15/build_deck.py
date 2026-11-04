"""Build the 2-slide CARDIOVERSE Sprint 15 report-out deck (EHT twitch identifiability)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor

here = os.path.dirname(os.path.abspath(__file__))
NAVY, GREY, RED = RGBColor(0x1F, 0x2A, 0x44), RGBColor(0x55, 0x5B, 0x66), RGBColor(0xB0, 0x24, 0x18)

# --- crop panels A, B, E of the non-uniqueness figure (top row A,B + bottom-middle E) ---
im = Image.open(os.path.join(here, "eht_nonunique.png"))
W, H = im.size
im.crop((0, int(0.03 * H), int(0.665 * W), int(0.50 * H))).save(os.path.join(here, "crop_AB.png"))
im.crop((int(0.33 * W), int(0.50 * H), W, H)).save(os.path.join(here, "crop_EF.png"))

# --- schematic: relative drug features on a twitch ---
t = np.linspace(0, 1, 1000)
def tw(tp, a, tr):
    x = np.where(t < tp, (t / tp) ** 2 * (3 - 2 * t / tp), np.exp(-((t - tp) / tr) ** 2))
    return a * x
base, drug = tw(0.23, 1.0, 0.17), tw(0.26, 0.65, 0.21)
fig, ax = plt.subplots(figsize=(5.2, 3.3), dpi=200)
ax.plot(t * 1000, base, color="k", lw=2, label="baseline (DMSO / pre-dose, same tissue)")
ax.plot(t * 1000, drug, color="#B02418", lw=2, label="drug")
for y, tp, c in [(base, 0.23, "k"), (drug, 0.26, "#B02418")]:
    ax.axvline(tp * 1000, color=c, ls=":", lw=1)
ax.annotate("", xy=(260, 0.6), xytext=(230, 0.6), arrowprops=dict(arrowstyle="<->", color="#555"))
ax.text(245, 0.55, "ΔTTP", ha="center", fontsize=8, color="#555")
ax.annotate("Δ amplitude\n(% of own baseline)", xy=(262, 0.62), xytext=(330, 0.9), fontsize=8, color="#555",
            arrowprops=dict(arrowstyle="->", color="#555"))
ax.text(95, 0.33, "Δ max +dF/dt", fontsize=8, color="#555", rotation=55)
ax.text(420, 0.45, "Δ max −dF/dt,\nΔ RT50 / RT90", fontsize=8, color="#555")
ax.set_xlabel("time from stimulus (ms)"); ax.set_ylabel("force (rel. to baseline peak)")
ax.set_xlim(0, 800); ax.set_ylim(-0.03, 1.08)
ax.legend(fontsize=7.5, frameon=False, loc="upper right", bbox_to_anchor=(1, 0.8))
ax.set_title("schematic — not data", fontsize=8, color="#888", loc="right")
for s in ["top", "right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(here, "schematic_features.png")); plt.close(fig)

# --- deck ---
prs = Presentation(); prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
blank = prs.slide_layouts[6]

def text(slide, x, y, w, h, paras, size=14):
    tf = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)).text_frame
    tf.word_wrap = True
    for i, p in enumerate(paras):
        s, kw = (p, {}) if isinstance(p, str) else p
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        r = para.add_run(); r.text = s
        r.font.size = Pt(kw.get("size", size)); r.font.bold = kw.get("bold", False)
        r.font.color.rgb = kw.get("color", GREY); para.level = kw.get("level", 0)
        para.space_after = Pt(kw.get("after", 4))
    return tf

def header(slide, title, sub):
    text(slide, 0.5, 0.3, 12.3, 0.7, [(title, dict(size=28, bold=True, color=NAVY))])
    text(slide, 0.5, 0.95, 12.3, 0.5, [(sub, dict(size=15, color=GREY))])

# Slide 1: the negative result
s = prs.slides.add_slide(blank)
header(s, "Force alone can't tell Ca handling from cross-bridges",
       "Model fitted to the averaged twitch of 15 control tissues, under four different assumptions")
s.shapes.add_picture(os.path.join(here, "fits_simple.png"), Inches(0.3), Inches(1.6), width=Inches(12.7))
text(s, 0.6, 5.6, 6.6, 1.8, [
    "Slower or faster cross-bridges fit the twitch just as well: a different (unmeasured) Ca transient compensates",
    "The fits disagree on what we care about: energy use differs 4× between them",
], size=15)
box = s.shapes.add_shape(1, Inches(7.5), Inches(5.5), Inches(5.4), Inches(1.75))
box.fill.solid(); box.fill.fore_color.rgb = RGBColor(0xFB, 0xEC, 0xEA); box.line.fill.background()
text(s, 7.65, 5.55, 5.15, 1.7, [
    ("From the meeting (with the EHT lab)", dict(bold=True, color=NAVY, size=14)),
    ("Their setup (Curi Bio, 37 °C, DMSO control) can't record the Ca transient or absolute tension (no cross-section), now or likely in future",
     dict(size=13)),
    ("⇒ the full model can't be identified from these data", dict(bold=True, color=RED, size=14)),
], size=13)

# Slide 2: alternatives
s = prs.slides.add_slide(blank)
header(s, "Way forward: ask questions the data can answer",
       "Don't fit the full model to immature hiPSC tissue; look at how drugs change the twitch")
text(s, 0.5, 1.7, 6.2, 2.4, [
    ("A · Simpler model", dict(bold=True, color=NAVY, size=18)),
    "Keep only what the twitch pins down: activation timing, one cross-bridge speed, amplitude",
    "Anchor it with what we can vary: pacing rate and drugs with a known target (e.g. mavacamten → myosin)",
], size=15)
text(s, 0.5, 4.1, 6.2, 2.4, [
    ("B · Relative drug effects", dict(bold=True, color=NAVY, size=18)),
    "Compare each tissue to its own baseline: time to peak, rise/relaxation speed, % amplitude",
    "Each mechanism (Ca handling vs cross-bridges) leaves its own pattern of changes",
], size=15)
s.shapes.add_picture(os.path.join(here, "schematic_features.png"), Inches(7.0), Inches(1.6), width=Inches(5.9))
text(s, 7.0, 5.45, 5.9, 1.9, [
    ("A + B: simple model on controls → predict each drug's pattern → compare with measured",
     dict(bold=True, color=NAVY, size=15)),
    ("Next: Filip: model reduction · Abhijit & Kalai: more measurement types · Matt: Curi Bio post stiffness",
     dict(size=13)),
], size=13)

out = os.path.join(here, "Sprint15_EHT_twitch.pptx")
prs.save(out); print(out)
