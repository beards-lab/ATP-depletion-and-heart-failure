# Online literature resources saved for later

## Calcium
- Calcium and Excitation-Contraction Coupling in the Heart (Eisner 2017)
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5497788/
  - review
- Calcium diffusion in transient and steady states in muscle (Safford and Basingwaighte, 1977)
  - https://www.cell.com/biophysj/pdf/S0006-3495(77)85539-2.pdf
  - in depth analysis, yet mathematically copmlicated

## XB detachment rates
- The Rates of Ca2+ Dissociation and Cross-bridge Detachment from Ventricular Myofibrils as Reported by a Fluorescent Cardiac Troponin C\* (Little 2012)
  - https://www.jbc.org/article/S0021-9258(20)47816-0/fulltext
  - human 40/s, rat 150/s, with decrease to about 1/3 while ADP present

## Another mechanisms
- Force-Dependent Recruitment from the Myosin Off State Contributes to Length-Dependent Activation (Campbell 2018)
  - https://www.sciencedirect.com/science/article/pii/S0006349518307707
  -  model in which the rate of the off-to-on transition increased linearly with force reproduced the length-dependent behavior of chemically permeabilized myocardium better than a model with a constant off-to-on transition rate... and ... thick-filament transitions are modulated by force... 
  
## Force velocity
- Altered single cell force-velocity and power properties in exercise-trained rat myocardium (Difee 2003)
  - https://journals.physiology.org/doi/full/10.1152/japplphysiol.00889.2002
  - unloaded shortening velocity was not significantly different... Training increased the velocity of loaded shortening and increased peak power output
  - servo-controlled force clamp (measuring at force, not at velocity)
  - significantly reduced force-velocity profile and Vmax compared to our experiments (1.43 ML/s) - HOW??? (15degC)
- Force-velocity and power-load curves in rat skinned cardiac myocytes (McDonald 1998)
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2231141/
  - again, very low Vmax (1.5 ML/s) and max tension (26kPa), low T (12degC)
  - found two-phase relation of step-down, some experiments on pCa, 
  
## Length - tension
- Changes in rat myocardium contractility under subchronic intoxication with lead and cadmium salts administered alone or in combination (Prostsenko 2020)
  - https://www.sciencedirect.com/science/article/pii/S2214750020300445
  - relative passive and active length-tension for rat papillary and trabeculae muscle at 35 degC, length -tension comparable to our values
- Cardiac muscle mechanics: Sarcomere length matters (Tombe 2016)
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5457809/
  - serial stiffness in the measurement rig affects the contraction length, i.e. the SL changes during contraction by up to 20%, although the ML is held constant
  - based on their research from 1991
  

## viscosity
- AN INTERNAL VISCOUS ELEMENT LIMITS UNLOADED VELOCITY OF
SARCOMERE SHORTENING IN RAT MYOCARDIUM (Tombe 1991)
  - https://physoc.onlinelibrary.wiley.com/doi/pdfdirect/10.1113/jphysiol.1992.sp019283
  - 25degC, similar force-velocity to ours (Vmax~10 um/s), but calculated LINEAR sarcomere force-velocity from dynamic stiffness.
  
## Titin
- Tampering with springs: phosphorylation of titin affecting the mechanical function of cardiomyocytes (Hamdani 2017)
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5498327/
  - phosphorylation at specific amino acids can reduce or increase the stretch-induced spring force of titin, depending on where the spring region is phosphorylated. Reviews only passive properties, no regard to ATP.
- Titin as a force-generating muscle protein under regulatory control (Freundt 2019)
  - https://journals.physiology.org/doi/full/10.1152/japplphysiol.00865.2018
  - reviews titin-based viscoelastic forces in skeletal muscle cells, including chaperone binding, titin oxidation, phosphorylation, Ca2+ binding, and interaction with actin filaments
- What Can We Learn from Single Sarcomere and Myofibril Preparations? (Herzog 2022)
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9092595/ review
  - uncovering of the crucial role of titin in (skeletal) muscle contraction
  - titin is a calcium dependent molecular spring; that is, the stiffness of titin changes depending on the active state of a muscle. Passive is about 20% higher with calcium
  - residual force enhancement (step over the optimum given by the overlap function)
    - dedicated paper Residual Force Enhancement Following Eccentric Contractions: A New Mechanism Involving Titin (Herzog 2016), https://journals.physiology.org/doi/full/10.1152/physiol.00049.2014
  - binding part of titin to actin during active contraction
- A Spatially Explicit Model Shows How Titin Stiffness Modulates Muscle Mechanics and Energetics (Powers 2018)
  - https://academic.oup.com/icb/article/58/2/186/5036269
  - adding titin nonlinear spring to a hexagonal symmetric matrix of sarcomeres. Three state XB model.
  

## Another computational models
- Mathematical Simulation of Muscle Cross-Bridge Cycle and Force-Velocity Relationship (Chin 2006)
  - https://www.sciencedirect.com/science/article/pii/S000634950672077X
  - seven state XB model, with velocity dependence, ADP and ATP effects included (yet seemingly far from our observations)
- A model of cardiac contraction based on novel measurements of tension development in human cardiomyocytes (Land 2017)
  - https://www.sciencedirect.com/science/article/abs/pii/S0022282817300639?casa_token=-l6ckz5kjc8AAAAA:9TVS0ZRG4TU5YSGdUeh9mazkkqFWqxOMecEfKDhJUezv4smPccz6QXpy22_E_ybRGGYhIRgm_s8#bb0315
  - human cardiomyocyte at 37deC, 3 state XB model, velocity dependent transitions, pCa - tension fitting to 1.8, 2.0 and 2.2 um SL
  - dashpot with serial stiffness
  
## super relaxed state
- The role of super-relaxed myosin in skeletal and cardiac muscle (McNamara 2015)
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5425749/
  - might play a role in hypertrophic heart disease and hypoxia / infarction
  - cooperative interaction between adjacent myosin heads is present in skeletal but not cardiac thick filaments 
  - ATP lifetime of ~230s in SRX
  