# Sensitivity Analysis — 108-Parameter Results & Outlook

## Section 1: Result Summary (108 Parameters)

By expanding to 108 parameters (all "tunable" fields), we have tested the model's full flexibility. Despite adding 68 new parameters, the **active subspace is still only 11 dimensions**. 

### Raw Sensitivity Top 8
1. **`ekSE`** (4895): Serial element nonlinearity is the strongest global controller.
2. **`estiff`** (3346): Passive/structural stiffness scaling.
3. **`dr2`** (1709): Transition geometry.
4. **`k2`** (1000): ADP release kinetics.
5. **`dr`** (974): Power-stroke size.
6. **`kstiff2`** (937): Post-stroke stiffness.
7. **`kSE`** (807): Serial stiffness.
8. **`kah`** (682): Hydrolysis rate.

---

## Section 2: SVD & Identifiability — The "11-DOF" Ceiling

The model has a fundamental structural constraint. No matter how many variables we add, the experimental features we are targeting (Slack, FV, Ktr) only provide enough information to constrain **11 independent combinations**.

### Top Identifiable Parameters (The "Knobs" that work)
These 12 parameters (or their linear combinations) are the ones that actually "move" the model in unique ways:
1. `ekSE`
2. `estiff`
3. `mu_neg` (viscosity)
4. `dr` / `dr2` (geometry)
5. **`L_thin`** (Filament length): **Newly Identified** — filament geometry is crucial for fitting the overlap effects in slack tests.
6. **`xrate`**: Rate scaling.
7. **`s_threshold_R`**: **Newly Identified** — the strain threshold for detachment transition.
8. `kah`, `kSE`, `k1`.

---

## Section 3: Recommendations & Outlook

### 1. The Optimization Strategy
Since we have 11 degrees of freedom, our previous "Top 11" selection is still extremely robust. However, we should now include **`L_thin`** and **`s_threshold_R`** as they provide unique information that the kinetic rates cannot replicate.

**Recommended Optimizer Set:**
`ekSE`, `estiff`, `mu_neg`, `dr`, `dr2`, `L_thin`, `s_threshold_R`, `kah`, `kSE`, `k1`, `k2`, `kstiff2`.

### 2. Outlook on "Missing Physics"
The persistent limit of 11 DOF suggests that if the current 12 parameters cannot close the gap to the experimental data, **adding more parameters won't help.** The system is essentially "saturated" in terms of what can be tuned. 

If residuals remain high after the next optimization round, we have mathematical proof that the **mechanism structure itself** (Action 4 or similar) needs correction, as we've already exhausted the parameter space's ability to respond.

### 3. Immediate Action
We are now ready to generate the final `OptimizeMechanismEvaluation.m` using this refined set of 12 parameters.
