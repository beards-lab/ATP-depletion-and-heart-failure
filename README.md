# ATP-depletion-and-heart-failure


## Naming conventions:

**FirstCapital.m** a script, that leaves all the variables in workspace

**firstLowercase.m** is a function with a return value

## Source code structure
Core functionality (Top down):

- [Driver...](RunOptim.m) an entry point for runs. Loads the parameters *params0* and runs the  [RunBakersExp](RunBakersExp.m) 
- [RunBakersExp](RunBakersExp.m) evaluates model challenges, e.g. force velocity, slack etc, based on switches in *params0* and  calculates all the model evaluations against data and calculates total error E. Could be parametrized to run only certain challenges.
- [EvaluateModel](EvaluateModel.m) - evaluates the model for single or multiple velocity segments, optionally stores the result, takes care of interpolating at the exact sarcomere length etc. 
- [getParams](getParams.m) - sets all the parametrization at single place. When provided with an input params struct, updates it with all other required fields.
- dPUdT...(eg. [dPUdTCaSimpleAlternative2State](dPUdTCaSimpleAlternative2State.m)) - The core function containing the differential equations. 

The switches changing behavior include:

- UseOverlap: 1
- UsePassive: 1
- UseSuperRelaxed: 1
- UseTitinInterpolation: 0
- MaxSlackNegativeForce: 0
- UseDirectSRXTransition: 0
- UsePassiveForSR: 0
- FudgeVmax: 0
- UseSuperRelaxedADP: 1
- UseOverlapFactor: 1
- UseTitinIdentifiedPassive: 1
- UseUniformTransitionFunc: 0
  
Switches, currently inactive:

- UseSlack
- UseP31Shift
- F_act_UseP31
- UseMutualPairingAttachment
- UseSpaceDiscretization
- UseSpaceInterpolation
- UseKstiff3
- EvalAtp
- vmax 
- UseSerialStiffness


Other functions and tests:
- [RunModel](RunModel.m) - simplest way to run the model (Demo).
- [AnimateStateProbabilities](AnimateStateProbabilities.m) - animates state probabilities from the EvaluateModel output.
- [force_pCa](force_pCa.m) calculates Force-pCa curve at varying ATP levels
- [RunOptimLakes](RunOptimLakes.m) an entry point for optimizations. 
- [LoadBakersExp](LoadBakersExp.m) - Loads Experiments from Anthony Baker and plots them

