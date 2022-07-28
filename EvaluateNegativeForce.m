function E = EvaluateNegativeForce(g)

E = evaluateModel(@dPUdT, -6, 1, 2,0,0,g);