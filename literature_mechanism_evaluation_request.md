# Comprehensive Literature-Guided Mechanism Evaluation Request

## Background

I have a deterministic mechanistic muscle model with approximately 10
major physical state variables. Two of these states are spatially
discretized, increasing the total number of numerical states to several
dozen or around one hundred.

My experimental dataset consists of:

-   **Force--velocity (FV) experiments**
-   **Slack--restretch (ktr) experiments**

These experiments are non-physiological and are specifically designed to
uncover mechanisms that may remain hidden under in vivo--like
conditions.

My primary goal is **to understand the physics**, not merely improve
predictive accuracy.

------------------------------------------------------------------------

## Part 1 --- Mechanism Evaluation

Please perform a structured, literature-anchored analysis of the
mechanisms currently implemented in the model.

For each mechanism:

1.  Assess whether it is biophysically realistic.
2.  Determine whether it has been previously used in published muscle
    models.
3.  Provide literature references supporting or contradicting its use.
4.  Identify whether the mechanism:
    -   Can likely be rejected as unrealistic or unnecessary for
        explaining FV and slack data.
    -   Is well-established and appropriate.
    -   Is partially implemented and may require refinement.

Additionally, propose **new mechanisms** that:

-   Have literature support.
-   Are specifically relevant to force--velocity and slack--restretch
    protocols.
-   Could plausibly explain discrepancies between model and experiment.
-   Include proper literature citations.

------------------------------------------------------------------------

## Part 2 --- Parameter Range Comparison

1.  Compare the model's parameter ranges to approximate ranges reported
    in literature.
2.  Identify:
    -   Parameters within typical literature bounds.
    -   Parameters that appear unrealistic, poorly justified, or
        inconsistent with experimental findings.
3.  Provide citations for reported parameter ranges wherever possible.

------------------------------------------------------------------------

## Part 3 --- Protocol Sensitivity Analysis

Because force--velocity and slack--restretch experiments probe specific
dynamic properties, please:

-   Identify which mechanisms are most sensitive to these protocols.
-   Explain how each mechanism influences:
    -   FV curvature and maximum shortening velocity (Vmax)
    -   Force redevelopment kinetics (ktr)
    -   Transient responses after slack or stretch

------------------------------------------------------------------------

## Part 4 --- Structured Output Requirements

Please provide:

-   A categorized table of mechanisms:
    -   Well-supported
    -   Questionable
    -   Likely missing
-   A literature-referenced justification for each category.
-   A prioritized list of mechanisms to test or modify.

------------------------------------------------------------------------

## Constraints

-   Focus on mechanistic interpretability.
-   Do not suggest black-box machine learning correction approaches.
-   Emphasize mechanisms supported by experimental or modeling
    literature.
-   Include citations for each literature-backed claim.

------------------------------------------------------------------------

## Assumed Model Features

You may assume the model contains:

-   Strain-dependent cross-bridge kinetics
-   Series elastic elements
-   Possibly super-relaxed states
-   Various attachment/detachment rate constants
-   Optional viscoelastic or force-dependent detachment components

------------------------------------------------------------------------

## Final Goal

Determine:

> Which mechanisms are realistic and justified, which can be rejected,
> and which are missing but necessary to explain slack--restretch and
> force--velocity data.
