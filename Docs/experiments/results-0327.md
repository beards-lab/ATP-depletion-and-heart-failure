# Experiment Results — 03/27/2026

**Fiber:** M  
**Conditions:** 8mM ATP (active), 2mM ATP (active), 8mM repeat (rundown), 8mM PNB+Mava (passive)  
**Protocol:** ktr, slack at SL 2.0 and 2.2

---

## Force Comparison: 8mM vs 2mM vs PNB+Mava

Merged hi-res + fine-shifted traces. The 8mM repeat (i=4) is upscaled using SL-dependent rundown compensation (see section below) and overlaid as a dashed trace to confirm alignment with i=2.

![Force comparison](../data/03%2027%202026%20M/Report_Comparison_235.png)

**Key observations:**
- 2mM (red) runs ~20–30 kPa higher than 8mM (blue) — ATP depletion slows detachment, increasing steady-state cross-bridge force
- Upscaled i=4 repeat (purple dashed) tracks i=2 closely across all protocols, validating the rundown compensation
- PNB+Mava (green) ~0–5 kPa passive baseline throughout
- **ktr** (69.5–74 s): 2mM maintains a higher post-drop redevelopment plateau; reduced ATP slows detachment and elevates isometric force
- **Slack** (74–77.5 s): 2mM maintains elevated force through each slack release; recovery rate and plateau are both higher; the upscaled repeat closely follows the first 8mM run

---

## Rundown Rescaling (i=4 → i=2)

Active force is isolated by subtracting the PNB+Mava passive baseline (i=5), scaled by a SL-dependent factor, then passive is restored:

$$F_\text{active} = F - F_\text{passive}$$

$$F_\text{corrected} = F_\text{active,rd} \cdot \bigl(r_0 + k\,(L - L_0)\bigr) + F_\text{passive}$$

with $r_0 = 68/56 = 1.214$, $L_0 = 1.0\,L_o$, $k = -0.6$. $r_0$ is the plateau force ratio i=2 : i=4 at SL=2.0; $k < 0$ reduces the correction at higher SL, where a titin/lattice contribution persists through rundown.

![Rundown rescaling](../data/03%2027%202026%20M/Report_Rundown_Rescale.png)

---

## Rundown Linear Fits

Force decline rate estimated from two SL conditions across the protocol window (70–79.5 s).

![Rundown fits](../data/03%2027%202026%20M/Report_Rundown_23.png)

Fitting method:
- **SL=2.0:** linear fit through plateau zones `[71.5–72.4 s]` and `[77.4–78.4 s]` where L ≈ 1.0 Lo
- **SL=2.2:** linear fit through 150 ms windows immediately before each L-drop from L ≥ 1.1 Lo

| Condition | SL = 2.0 slope (kPa/s) | SL = 2.2 slope (kPa/s) |
|-----------|------------------------|------------------------|
| 2mM Active (i=3) | −1.311 | −1.303 |
| 8mM Active (i=2) | −0.495 | +0.082 |

**Notes:**
- 2mM rundown is ~2.6× faster than 8mM at SL=2.0
- 8mM is essentially flat at SL=2.2, suggesting the elevated-SL component does not run down — consistent with a titin/passive contribution at high stretch
- Both SL fits agree closely in 2mM, indicating uniform rundown rate independent of SL under ATP-depleted conditions
