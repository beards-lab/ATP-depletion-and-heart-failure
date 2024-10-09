# ATP-depletion-and-heart-failure


## Source code structure

Core functionality (Top down):

- [RunOptim](RunOptim.m) an entry point for optimizations. Calls optim on [EvaluateProblem](EvaluateProblem.m) in various ways, and calculates sensitivities using [calcSensitivities](calcSensitivities.m)(calcSensitivities.m) under various options
- [EvaluateProblem](EvaluateProblem.m) - calculates all the model evaluations against data and calculates total error E. Could be parametrized to run only certain challenges.
- [LoadData](LoadData.m) - sets all data used for comparison at [EvaluateProblem](EvaluateProblem.m)
- [EvaluateModel](EvaluateModel.m) - evaluates the model for single or multiple velocity segments, optionally stores the result, takes care of interpolating at the exact sarcomere length etc. 
- [getParams](getParams.m) - sets all the parametrization at single place. When provided with an input params struct, updates it with all other required fields.
- [dPUdTCa](dPUdTCa.m) - The main dpudt extended with Ca sensitivity and other features.


Other functions and tests:
- [RunModel](RunModel.m) - simplest way to run the model (Demo).
- [AnimateStateProbabilities](AnimateStateProbabilities.m) - animates state probabilities from the EvaluateModel output.
- [force_pCa](force_pCa.m) calculates Force-pCa curve at varying ATP levels
- [RunOptimLakes](RunOptimLakes.m) an entry point for optimizations. 
- [LoadBakersExp](LoadBakersExp.m) - Loads Experiments from Anthony Baker and plots them

