within ;
package CrossBridgeCycling
  model XB_balance
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));

  end XB_balance;

  model SR_Balance

    Real f1=75;
    Real f2=1;
    Real SR1 = 0.4;
    Real SR2 = 0.99;

    Real qSR1=ksr0*SR1*(f1/sigma1);
    Real qNSR1=kmsr*(1-SR1)*exp(-f1/sigma2);

    Real qSR2=ksr0*SR2*(f2/sigma1);
    Real qNSR2=kmsr*(1-SR2)*exp(-f2/sigma2);

  //   Real sigma1=10;
  //   Real sigma2=10;
  //   Real ksr0=10;
  //   Real kmsr=10;
    Real sigma1;
    Real sigma2;
    Real ksr0 = 10;
    Real kmsr = 1000;

  equation
  //  der(NSR) = qSR - qNSR;

    qSR1 = qNSR1;
    qSR2 = qNSR2;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end SR_Balance;

  model XBRates_SS

  Real g1 = 1, g2 = 1, f1 = 0.1, f2 = 1;

  record Params
      Real kah = 80;
      Real kadh = 0;
      Real ka, kd;

      Real k1, k_1, k2;
      Real ksr0, kmsr;

      Real k2_L = 0, k2_R = 0;
      Real sigma1 = 1e9, sigma2  = 1e9;
      Real alpha0 = 0, dr0 = 0;
      Real alpha1 = 0;
      Real alpha_1 = 0, dr_1 = 0;
  //     Real alpha2_L, alpha2_R;
  //     Real dr2_L, dr2_R;
  end Params;

  // STATES
  Real PD, p1, p2,  P_SR;

  parameter Real F_total = 60;

  Real dS = 1, s = 1;
  Params params;

  Real N_overlap = overLap.O_overlap;

  // the cycle goes: PT (ATP bound) <-> PD(ready) <-> P1 <-> P2 -> P3 -> PT
  function strainDep
   input Real alpha, dr;
   // Real output o = min(1e4, exp((alpha*(s+dr)).^2));
   output Real o = 1;
  end strainDep;

  Real RTD = g2*params.kah*PT;
  Real RD1 = params.ka*PD*N_overlap; // to loosely attachemnt state
  Real R1D = params.kd*p1; // *strainDep(params.alpha0, params.dr0)
  //// R21 = f1*params.k_1*p2.*strainDep(params.alpha_1, params.dr_1); // backward flow from p2 to p1

  Real R12 = params.k1*p1.*exp(-params.alpha1*s); // P1 to P2
  Real R21 = f1*params.k_1*p2; // *strainDep(params.alpha_1, params.dr_1); // p2 to p1
  // R2T = g2*params.k2*p2.*min(1e9, max(1, strainDep(params.alpha2, params.dr2)));

  Real kL = 0;//*min(1e4, params.k2_L*((s-params.dr2_L)).*(exp(-(s-params.dr2_L)*params.alpha2_L)));
  Real kR = 0;//*max(0, params.k2_R*(s-params.dr2_R)).^params.alpha2_R;
  Real R2T = p2.*(params.k2 + kL + kR);

  // if params.UseSuperRelaxed && params.UseDirectSRXTransition
  Real RSR2PT = params.ksr0*exp(F_total/params.sigma1)*P_SR;
  Real RPT2SR = params.kmsr*exp(-F_total/params.sigma2)*PT;

  Real dU_SR = -RSR2PT  + RPT2SR + sum(R2T)*dS;
  // governing flows
  // PT - calculated as complement of sum of all probabilities
  Real PT = 1 - P_SR - PD - p1 - p2;
  // state 0: unattached, ATP-cocked
  Real dPD = + RTD - RD1 + R1D;
  Real dp1 = RD1 - R1D -  R12 + R21; // state 1: loosely attached, just sitting&waiting
  Real dp2 = + R12 - R21  - R2T; // strongly attached, post-ratcheted: hydrolyzed ATP to ADP, producing Pi - ready to ratchet
  Real SL = 2.2;
    XB.XBMoments.Tests.OverLap overLap(SL=SL)
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  equation
    // SET the KNOWLEDGE - states at maximally contracted
    PD = 0.2;
    p1 = 0.1;
    p2 = 0.4;
    P_SR = 0.2;

    // steady state
    dp1 = 0; dp2 = 0; dPD = 0; dU_SR = 0;

    // turnover rate is 10
    //   params.k2 = 10;
    RD1 = 14;

  //   params.k1 = 30;
    // Ratio of backward flows
  //   params.ka / params.kd = 1/5;
  //   params.k_1 / params.k1 = 1/5;
    params.kmsr / params.ksr0 = 1/5;
    R21/R12 = 1/5;

  //   R2T = 10;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end XBRates_SS;

  model elasticSegment
    Modelica.Mechanics.Translational.Components.Spring spring(c=c)
      annotation (Placement(transformation(extent={{20,-10},{40,10}})));
    Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
    Modelica.Mechanics.Translational.Components.Spring spring1(c=cStiff1*att)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,56})));
    Modelica.Mechanics.Translational.Interfaces.Support support
      annotation (Placement(transformation(extent={{-10,90},{10,110}})));
    parameter Modelica.Units.SI.TranslationalSpringConstant c=1
      "Spring constant";
    parameter Modelica.Units.SI.TranslationalSpringConstant cStiff1=1
      "Spring constant";
    parameter Real att=1;
  equation
    connect(spring.flange_a, flange_a)
      annotation (Line(points={{20,0},{-100,0}}, color={0,127,0}));
    connect(spring.flange_b, flange_b)
      annotation (Line(points={{40,0},{100,0}}, color={0,127,0}));
    connect(support, spring1.flange_b)
      annotation (Line(points={{0,100},{0,66}}, color={0,127,0}));
    connect(spring1.flange_a, flange_a)
      annotation (Line(points={{0,46},{0,0},{-100,0}}, color={0,127,0}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end elasticSegment;

  model HlasSercomere
    elasticSegment elasticSegment1(
      c=ThickStiffness,
      cStiff1=HeadStiffness,
      att=4)
      annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
    Modelica.Mechanics.Translational.Components.Fixed fixed
      annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
    Modelica.Mechanics.Translational.Components.Spring spring(c=1)
      annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
    elasticSegment elasticSegment2(
      c=ThickStiffness,
      cStiff1=HeadStiffness,
      att=0)
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    elasticSegment elasticSegment3(
      c=ThickStiffness,
      cStiff1=HeadStiffness,
      att=0)
      annotation (Placement(transformation(extent={{32,-10},{52,10}})));
    Modelica.Mechanics.Translational.Sources.Position position
      annotation (Placement(transformation(extent={{104,-10},{84,10}})));
    Modelica.Blocks.Sources.Sine sine(f= 1)
      annotation (Placement(transformation(extent={{40,38},{60,58}})));
    parameter Modelica.Units.SI.TranslationalSpringConstant ThickStiffness=100
      "Spring constant";
    parameter Modelica.Units.SI.TranslationalSpringConstant HeadStiffness=10
      "Spring constant";
    elasticSegment elasticSegment4(
      c=ThickStiffness,
      cStiff1=HeadStiffness,
      att=0)
      annotation (Placement(transformation(extent={{58,-10},{78,10}})));
  equation
    connect(fixed.flange, elasticSegment1.support) annotation (Line(points={{-80,0},
            {-80,20},{-20,20},{-20,10}},
                                     color={0,127,0}));
    connect(fixed.flange, spring.flange_a) annotation (Line(points={{-80,0},{-80,20},
            {-66,20},{-66,0},{-60,0}}, color={0,127,0}));
    connect(spring.flange_b, elasticSegment1.flange_a)
      annotation (Line(points={{-40,0},{-30,0}}, color={0,127,0}));
    connect(elasticSegment2.support, elasticSegment1.support)
      annotation (Line(points={{10,10},{10,20},{-20,20},{-20,10}},
                                                               color={0,127,0}));
    connect(elasticSegment3.support, elasticSegment1.support)
      annotation (Line(points={{42,10},{42,20},{-20,20},{-20,10}},
                                                               color={0,127,0}));
    connect(elasticSegment1.flange_b, elasticSegment2.flange_a)
      annotation (Line(points={{-10,0},{0,0}}, color={0,127,0}));
    connect(elasticSegment3.flange_a, elasticSegment2.flange_b)
      annotation (Line(points={{32,0},{20,0}}, color={0,127,0}));
    connect(sine.y, position.s_ref) annotation (Line(points={{61,48},{116,48},{
            116,0},{106,0}},
                         color={0,0,127}));
    connect(position.flange, elasticSegment4.flange_b)
      annotation (Line(points={{84,0},{78,0}}, color={0,127,0}));
    connect(elasticSegment4.flange_a, elasticSegment3.flange_b)
      annotation (Line(points={{58,0},{52,0}}, color={0,127,0}));
    connect(elasticSegment4.support, elasticSegment1.support) annotation (Line(
          points={{68,10},{68,20},{-20,20},{-20,10}}, color={0,127,0}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end HlasSercomere;

  model SrxDrx
    Physiolibrary.Chemical.Components.Substance DRX(solute_start=1 - SRX_init)
      annotation (Placement(transformation(extent={{-62,-22},{-42,-2}})));
    Physiolibrary.Chemical.Components.Substance SRX(solute_start=SRX_init)
      annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
    Modelica.Blocks.Math.MultiSum multiSum(nu=2)
      annotation (Placement(transformation(extent={{24,28},{36,16}})));
    Physiolibrary.Chemical.Components.Clearance clearance(Clearance(displayUnit
          ="l/min") = 1.6666666666667e-05)
      annotation (Placement(transformation(extent={{-32,70},{-12,90}})));
    Physiolibrary.Chemical.Components.Clearance clearance1(Clearance(
          displayUnit="ml/min") = 0.000166667)
      annotation (Placement(transformation(extent={{-108,-22},{-128,-2}})));
    parameter Physiolibrary.Types.AmountOfSubstance SRX_init(displayUnit="mol")
       = 0.5 "Initial solute amount in compartment";
    parameter Physiolibrary.Types.VolumeFlowRate SolutionFlow=0
      "Volumetric flow of solution if useSolutionFlowInput=false";
  equation
    connect(SRX.solute, multiSum.u[1]) annotation (Line(points={{-44,40},{-44,
            23.05},{24,23.05}}, color={0,0,127}));
    connect(DRX.solute, multiSum.u[2]) annotation (Line(points={{-46,-22},{-46,
            -42},{-22,-42},{-22,20.95},{24,20.95}}, color={0,0,127}));
    connect(SRX.q_out, clearance.q_in) annotation (Line(
        points={{-50,50},{-64,50},{-64,80},{-32,80}},
        color={107,45,134},
        thickness=1));
    connect(clearance1.q_in, DRX.q_out) annotation (Line(
        points={{-108,-12},{-52,-12}},
        color={107,45,134},
        thickness=1));
    annotation (experiment(
        StopTime=300,
        __Dymola_NumberOfIntervals=1500,
        Tolerance=1e-05,
        __Dymola_Algorithm="Dassl"));
  end SrxDrx;

  model SrxDrx_fluxes
    extends SrxDrx(
      SolutionFlow=0,
      clearance1(Clearance(displayUnit="l/min") = DRXClearance,
          useSolutionFlowInput=true),
      clearance(Clearance=SRXClearance, useSolutionFlowInput=true));
    Physiolibrary.Chemical.Components.Stream Stream(SolutionFlow=relaxingRate)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-66,22})));
    Physiolibrary.Chemical.Components.Stream Stream1(SolutionFlow=DrxRate)
      annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=90,
          origin={-86,22})));
    Modelica.Blocks.Sources.RealExpression realExpression(y=if time > 0 then
          DRXClearance else 0)
      annotation (Placement(transformation(extent={{-156,-2},{-136,18}})));
    parameter Physiolibrary.Types.VolumeFlowRate DRXClearance(displayUnit=
          "l/min") = 0.000166667
      "Clearance of solute if useSolutionFlowInput=false";
    Modelica.Blocks.Sources.RealExpression realExpression2(y=if time > 0 then
          SRXClearance else 0)
      annotation (Placement(transformation(extent={{16,78},{-4,98}})));
    parameter Physiolibrary.Types.VolumeFlowRate SRXClearance(displayUnit=
          "l/min") = 1.6666666666667e-05
      "Clearance of solute if useSolutionFlowInput=false";
    parameter Physiolibrary.Types.VolumeFlowRate relaxingRate(displayUnit=
          "l/min") = 8.3333333333333e-07
      "Volumetric flow of solution if useSolutionFlowInput=false";
    parameter Physiolibrary.Types.VolumeFlowRate DrxRate(displayUnit="l/min")
       = 1.6666666666667e-06
      "Volumetric flow of solution if useSolutionFlowInput=false";
  equation
    connect(Stream.q_in, DRX.q_out) annotation (Line(
        points={{-66,12},{-66,-12},{-52,-12}},
        color={107,45,134},
        thickness=1));
    connect(Stream.q_out, SRX.q_out) annotation (Line(
        points={{-66,32},{-66,36},{-64,36},{-64,50},{-50,50}},
        color={107,45,134},
        thickness=1));
    connect(Stream1.q_in, clearance.q_in) annotation (Line(
        points={{-86,32},{-86,80},{-32,80}},
        color={107,45,134},
        thickness=1));
    connect(Stream1.q_out, DRX.q_out) annotation (Line(
        points={{-86,12},{-86,-12},{-52,-12}},
        color={107,45,134},
        thickness=1));
    connect(realExpression.y, clearance1.solutionFlow)
      annotation (Line(points={{-135,8},{-118,8},{-118,-5}}, color={0,0,127}));
    connect(realExpression2.y, clearance.solutionFlow) annotation (Line(points=
            {{-5,88},{-6,88},{-6,98},{-22,98},{-22,87}}, color={0,0,127}));
    annotation (experiment(
        StartTime=-3000,
        StopTime=300,
        __Dymola_NumberOfIntervals=1500,
        Tolerance=1e-05,
        __Dymola_Algorithm="Dassl"));
  end SrxDrx_fluxes;

  model UtUdSrtSrd
    Physiolibrary.Chemical.Components.Substance UT(solute_start=0.5 - SRX_init/
          2)
      annotation (Placement(transformation(extent={{-60,-80},{-40,-60}})));
    Physiolibrary.Chemical.Components.Substance UD(solute_start=0.5 - SRX_init/
          2) annotation (Placement(transformation(extent={{38,-80},{58,-60}})));
    parameter Physiolibrary.Types.AmountOfSubstance SRX_init(displayUnit="mol")
       = 0.5 "Initial solute amount in compartment";
    parameter Physiolibrary.Types.VolumeFlowRate SolutionFlow=0
      "Volumetric flow of solution if useSolutionFlowInput=false";
    Physiolibrary.Chemical.Components.Substance SD(solute_start=SRX_init/2)
      annotation (Placement(transformation(extent={{40,74},{60,94}})));
    Physiolibrary.Chemical.Components.Substance ST(solute_start=SRX_init/2)
      annotation (Placement(transformation(extent={{-62,76},{-42,96}})));
    Physiolibrary.Chemical.Components.Stream kH(SolutionFlow(displayUnit=
            "l/min") = 4.1666666666667e-06)
      annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
    Physiolibrary.Chemical.Components.Stream kS2D(SolutionFlow(displayUnit=
            "l/min") = 8.3333333333333e-05) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={50,-10})));
    Physiolibrary.Chemical.Components.Stream kL(SolutionFlow(displayUnit=
            "l/min") = 8.3333333333333e-05)
      annotation (Placement(transformation(extent={{40,40},{20,60}})));
    Physiolibrary.Chemical.Components.Stream KS2T(SolutionFlow(displayUnit=
            "l/min") = 0.00016666666666667) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-50,30})));
    Modelica.Blocks.Math.MultiSum SRX(nu=2)
      annotation (Placement(transformation(extent={{80,54},{92,66}})));
    Modelica.Blocks.Math.MultiSum DRX(nu=2)
      annotation (Placement(transformation(extent={{82,-96},{94,-84}})));
    Physiolibrary.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure
      annotation (Placement(transformation(extent={{-20,60},{0,40}})));
    Physiolibrary.Chemical.Sources.UnlimitedSolutePumpOut
      unlimitedSolutePumpOut(useSoluteFlowInput=true)
      annotation (Placement(transformation(extent={{-80,40},{-100,60}})));
    Modelica.Blocks.Sources.RealExpression realExpression(y=if time >
          decay_time then -molarFlowMeasure.molarFlowRate else 0)
      annotation (Placement(transformation(extent={{-130,56},{-110,76}})));
    parameter Modelica.Blocks.Interfaces.RealOutput decay_time=0.0
      "Value of Real output";
    Modelica.Blocks.Math.MultiSum UnmarkedATP(nu=2)
      annotation (Placement(transformation(extent={{94,6},{82,-6}})));
  equation
    connect(ST.q_out, KS2T.q_in) annotation (Line(
        points={{-52,86},{-50,86},{-50,40}},
        color={107,45,134},
        thickness=1));
    connect(KS2T.q_out, UT.q_out) annotation (Line(
        points={{-50,20},{-50,-70}},
        color={107,45,134},
        thickness=1));
    connect(UT.q_out, kH.q_in) annotation (Line(
        points={{-50,-70},{-50,-30},{-40,-30}},
        color={107,45,134},
        thickness=1));
    connect(kH.q_out, UD.q_out) annotation (Line(
        points={{-20,-30},{50,-30},{50,-70},{48,-70}},
        color={107,45,134},
        thickness=1));
    connect(UD.q_out, kS2D.q_in) annotation (Line(
        points={{48,-70},{50,-70},{50,-20}},
        color={107,45,134},
        thickness=1));
    connect(SD.q_out, kL.q_in) annotation (Line(
        points={{50,84},{50,50},{40,50}},
        color={107,45,134},
        thickness=1));
    connect(ST.solute, SRX.u[1]) annotation (Line(points={{-46,76},{-50,76},{
            -50,72},{-66,72},{-66,100},{74,100},{74,58.95},{80,58.95}}, color={
            0,0,127}));
    connect(SD.solute, SRX.u[2]) annotation (Line(points={{56,74},{56,61.05},{
            80,61.05}}, color={0,0,127}));
    connect(UT.solute, DRX.u[1]) annotation (Line(points={{-44,-80},{-44,-91.05},
            {82,-91.05}}, color={0,0,127}));
    connect(UD.solute, DRX.u[2]) annotation (Line(points={{54,-80},{54,-88.95},
            {82,-88.95}}, color={0,0,127}));
    connect(kL.q_out, molarFlowMeasure.q_out) annotation (Line(
        points={{20,50},{0,50}},
        color={107,45,134},
        thickness=1));
    connect(molarFlowMeasure.q_in, ST.q_out) annotation (Line(
        points={{-20,50},{-50,50},{-50,86},{-52,86}},
        color={107,45,134},
        thickness=1));
    connect(unlimitedSolutePumpOut.q_in, KS2T.q_in) annotation (Line(
        points={{-80,50},{-50,50},{-50,40}},
        color={107,45,134},
        thickness=1));
    connect(realExpression.y, unlimitedSolutePumpOut.soluteFlow)
      annotation (Line(points={{-109,66},{-94,66},{-94,54}}, color={0,0,127}));
    connect(SRX.y, UnmarkedATP.u[1]) annotation (Line(points={{93.02,60},{98,60},
            {98,1.05},{94,1.05}}, color={0,0,127}));
    connect(DRX.y, UnmarkedATP.u[2]) annotation (Line(points={{95.02,-90},{100,
            -90},{100,-1.05},{94,-1.05}}, color={0,0,127}));
    connect(kL.q_in, kS2D.q_out) annotation (Line(
        points={{40,50},{50,50},{50,0}},
        color={107,45,134},
        thickness=1));
    annotation (experiment(
        StartTime=-2500,
        StopTime=5000,
        __Dymola_Algorithm="Dassl"),
      Diagram(coordinateSystem(extent={{-140,-100},{100,100}})),
      Icon(coordinateSystem(extent={{-140,-100},{100,100}})));
  end UtUdSrtSrd;

  model UtUdSrtSrd_V1
    extends UtUdSrtSrd(kS2D(SolutionFlow=8.3333333333333e-07), KS2T(
          SolutionFlow=3.3333333333333e-06));
    annotation (experiment(
        StartTime=-2500,
        StopTime=2500,
        __Dymola_Algorithm="Dassl"));
  end UtUdSrtSrd_V1;

  model UtUdSrtSrd_Ud2Ut "Reverse hydrolysis"
    extends UtUdSrtSrd(kS2D(SolutionFlow=1.6666666666667e-07), KS2T(
          SolutionFlow=1.6666666666667e-07));
    Physiolibrary.Chemical.Components.Stream Stream4(SolutionFlow(displayUnit=
            "l/min") = 1.6666666666667e-06)
      annotation (Placement(transformation(extent={{24,-66},{4,-46}})));
    Physiolibrary.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure1
      annotation (Placement(transformation(extent={{-32,-46},{-12,-66}})));
    Physiolibrary.Chemical.Components.Stream Stream5(SolutionFlow(displayUnit=
            "l/min") = 1.6666666666667e-06)
      annotation (Placement(transformation(extent={{-10,64},{10,84}})));
  equation
    connect(Stream4.q_out, molarFlowMeasure1.q_out) annotation (Line(
        points={{4,-56},{-12,-56}},
        color={107,45,134},
        thickness=1));
    connect(Stream4.q_in, UD.q_out) annotation (Line(
        points={{24,-56},{50,-56},{50,-70},{48,-70}},
        color={107,45,134},
        thickness=1));
    connect(molarFlowMeasure1.q_in, UT.q_out) annotation (Line(
        points={{-32,-56},{-50,-56},{-50,-70}},
        color={107,45,134},
        thickness=1));
    connect(Stream5.q_in, KS2T.q_in) annotation (Line(
        points={{-10,74},{-50,74},{-50,40}},
        color={107,45,134},
        thickness=1));
    connect(Stream5.q_out, SD.q_out) annotation (Line(
        points={{10,74},{50,74},{50,84}},
        color={107,45,134},
        thickness=1));
  end UtUdSrtSrd_Ud2Ut;

  model UtUdSrtSrd_BackS
    extends UtUdSrtSrd(kS2D(SolutionFlow=1.6666666666667e-06));
    Physiolibrary.Chemical.Components.Stream KS1T(SolutionFlow(displayUnit=
            "l/min") = 1.6666666666667e-06) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-70,30})));
    Physiolibrary.Chemical.Components.Stream KS1D(SolutionFlow(displayUnit=
            "l/min") = 1.6666666666667e-06) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={72,-10})));
  equation
    connect(KS1T.q_out, molarFlowMeasure.q_in) annotation (Line(
        points={{-70,40},{-70,50},{-20,50}},
        color={107,45,134},
        thickness=1));
    connect(KS1T.q_in, UT.q_out) annotation (Line(
        points={{-70,20},{-70,-30},{-50,-30},{-50,-70}},
        color={107,45,134},
        thickness=1));
    connect(KS1D.q_in, kL.q_in) annotation (Line(
        points={{72,0},{72,50},{40,50}},
        color={107,45,134},
        thickness=1));
    connect(KS1D.q_out, kH.q_out) annotation (Line(
        points={{72,-20},{72,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
  end UtUdSrtSrd_BackS;

  model UtUdSrtSrd_XBCycling "Remaining XB cycling in attached"
    extends UtUdSrtSrd;
    Physiolibrary.Chemical.Components.Stream XBCycling(SolutionFlow(displayUnit
          ="l/min") = 1.6666666666667e-06)
      annotation (Placement(transformation(extent={{20,-120},{0,-100}})));
    Physiolibrary.Chemical.Sensors.MolarFlowMeasure XBCyclingMeasure
      annotation (Placement(transformation(extent={{-40,-100},{-20,-120}})));
    Physiolibrary.Chemical.Sources.UnlimitedSolutePumpOut
      unlimitedSolutePumpOut1(useSoluteFlowInput=true)
      annotation (Placement(transformation(extent={{-78,-60},{-98,-40}})));
    Modelica.Blocks.Sources.RealExpression realExpression1(y=if time >
          decay_time then -XBCyclingMeasure.molarFlowRate else 0)
      annotation (Placement(transformation(extent={{-128,-44},{-108,-24}})));
  equation
    connect(XBCycling.q_out, XBCyclingMeasure.q_out) annotation (Line(
        points={{0,-110},{-20,-110}},
        color={107,45,134},
        thickness=1));
    connect(XBCycling.q_in, UD.q_out) annotation (Line(
        points={{20,-110},{64,-110},{64,-56},{50,-56},{50,-70},{48,-70}},
        color={107,45,134},
        thickness=1));
    connect(XBCyclingMeasure.q_in, UT.q_out) annotation (Line(
        points={{-40,-110},{-66,-110},{-66,-50},{-50,-50},{-50,-70}},
        color={107,45,134},
        thickness=1));
    connect(realExpression1.y, unlimitedSolutePumpOut1.soluteFlow) annotation (
        Line(points={{-107,-34},{-92,-34},{-92,-46}}, color={0,0,127}));
    connect(unlimitedSolutePumpOut1.q_in, UT.q_out) annotation (Line(
        points={{-78,-50},{-50,-50},{-50,-70}},
        color={107,45,134},
        thickness=1));
    annotation (Diagram(coordinateSystem(extent={{-140,-120},{100,100}})), Icon(
          coordinateSystem(extent={{-140,-120},{100,100}})));
  end UtUdSrtSrd_XBCycling;

  model UtUdSrtSrd_XBCycling_BackS
    extends UtUdSrtSrd_XBCycling(
      kL(SolutionFlow=1.6666666666667e-06),
      kH(SolutionFlow=0.00016666666666667),
      XBCycling(SolutionFlow=1.6666666666667e-05),
      kS2D(SolutionFlow=1.6666666666667e-06),
      KS2T(SolutionFlow=1.6666666666667e-06));
    Physiolibrary.Chemical.Components.Stream KS1T(SolutionFlow(displayUnit=
            "l/min") = 0)                   annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-70,30})));
    Physiolibrary.Chemical.Components.Stream KS1D(SolutionFlow(displayUnit=
            "l/min") = 0)                   annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={72,-10})));
  equation
    connect(KS1T.q_in, UT.q_out) annotation (Line(
        points={{-70,20},{-70,-30},{-50,-30},{-50,-70}},
        color={107,45,134},
        thickness=1));
    connect(KS1T.q_out, molarFlowMeasure.q_in) annotation (Line(
        points={{-70,40},{-70,50},{-20,50}},
        color={107,45,134},
        thickness=1));
    connect(KS1D.q_in, kL.q_in) annotation (Line(
        points={{72,0},{72,50},{40,50}},
        color={107,45,134},
        thickness=1));
    connect(KS1D.q_out, kH.q_out) annotation (Line(
        points={{72,-20},{72,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
  end UtUdSrtSrd_XBCycling_BackS;

  model DRXOnly
    Physiolibrary.Chemical.Components.Substance UT(solute_start=0.5 - SRX_init/2)
      annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
    Physiolibrary.Chemical.Components.Substance UD(solute_start=0.5 - SRX_init/2)
             annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
    Physiolibrary.Chemical.Components.Stream kH(SolutionFlow(displayUnit=
            "l/min") = 0.0016666666666667)
      annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
    Modelica.Blocks.Math.MultiSum DRX(nu=2)
      annotation (Placement(transformation(extent={{74,-70},{86,-58}})));
    Physiolibrary.Chemical.Components.Stream XBCycling(SolutionFlow(displayUnit
          ="l/min") = RelaxedATPCycling)
      annotation (Placement(transformation(extent={{20,-90},{0,-70}})));
    Physiolibrary.Chemical.Sensors.MolarFlowMeasure XBCyclingMeasure
      annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
    Physiolibrary.Chemical.Sources.UnlimitedSolutePumpOut
      unlimitedSolutePumpOut1(useSoluteFlowInput=true)
      annotation (Placement(transformation(extent={{-80,-20},{-100,-40}})));
    Modelica.Blocks.Sources.RealExpression realExpression1(y=if time > decay_time
           then -XBCyclingMeasure.molarFlowRate else 0)
      annotation (Placement(transformation(extent={{-60,-100},{-80,-80}})));
    Modelica.Blocks.Math.MultiSum UnmarkedATP(nu=2)
      annotation (Placement(transformation(extent={{86,28},{74,16}})));
    Physiolibrary.Chemical.Components.Substance SRX(solute_start=SRX_init)
      annotation (Placement(transformation(extent={{-44,56},{-24,76}})));
    Modelica.Blocks.Math.MultiSum multiSum(nu=1)
      annotation (Placement(transformation(extent={{74,52},{86,40}})));
    Physiolibrary.Chemical.Components.Stream Stream(SolutionFlow=relaxingRate)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-50,38})));
    Physiolibrary.Chemical.Components.Stream Stream1(SolutionFlow=DrxRate)
      annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=90,
          origin={-64,38})));
    parameter Physiolibrary.Types.AmountOfSubstance SRX_init(displayUnit="mol")
      =0.5   "Initial solute amount in compartment";
    parameter Physiolibrary.Types.VolumeFlowRate SolutionFlow=0
      "Volumetric flow of solution if useSolutionFlowInput=false";
    parameter Modelica.Blocks.Interfaces.RealOutput decay_time=0.0
      "Value of Real output";
    parameter Physiolibrary.Types.VolumeFlowRate DrxRate(displayUnit="l/min")=
      1.6666666666667e-06
      "Volumetric flow of solution if useSolutionFlowInput=false";
    parameter Physiolibrary.Types.VolumeFlowRate relaxingRate(displayUnit=
          "l/min")=8.3333333333333e-07
      "Volumetric flow of solution if useSolutionFlowInput=false";

    parameter Physiolibrary.Types.VolumeFlowRate RelaxedATPCycling(displayUnit=
          "l/min") = 0.00016666666666667
      "Volumetric flow of solution if useSolutionFlowInput=false";
  equation
    connect(UT.q_out,kH. q_in) annotation (Line(
        points={{-50,-50},{-50,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(kH.q_out,UD. q_out) annotation (Line(
        points={{0,-30},{30,-30},{30,-50}},
        color={107,45,134},
        thickness=1));
    connect(UT.solute,DRX. u[1]) annotation (Line(points={{-44,-60},{-44,-65.05},
            {74,-65.05}}, color={0,0,127}));
    connect(UD.solute,DRX. u[2]) annotation (Line(points={{36,-60},{36,-62.95},
            {74,-62.95}}, color={0,0,127}));
    connect(XBCycling.q_out,XBCyclingMeasure. q_out) annotation (Line(
        points={{0,-80},{-20,-80}},
        color={107,45,134},
        thickness=1));
    connect(XBCycling.q_in,UD. q_out) annotation (Line(
        points={{20,-80},{56,-80},{56,-30},{30,-30},{30,-50}},
        color={107,45,134},
        thickness=1));
    connect(XBCyclingMeasure.q_in,UT. q_out) annotation (Line(
        points={{-40,-80},{-64,-80},{-64,-30},{-50,-30},{-50,-50}},
        color={107,45,134},
        thickness=1));
    connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
        Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
    connect(unlimitedSolutePumpOut1.q_in,UT. q_out) annotation (Line(
        points={{-80,-30},{-50,-30},{-50,-50}},
        color={107,45,134},
        thickness=1));
    connect(SRX.solute,multiSum. u[1]) annotation (Line(points={{-28,56},{-28,
            46},{74,46}},       color={0,0,127}));
    connect(Stream.q_out,SRX. q_out) annotation (Line(
        points={{-50,48},{-50,66},{-34,66}},
        color={107,45,134},
        thickness=1));
    connect(Stream.q_in, kH.q_in) annotation (Line(
        points={{-50,28},{-50,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(Stream1.q_out, kH.q_in) annotation (Line(
        points={{-64,28},{-64,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(Stream1.q_in, SRX.q_out) annotation (Line(
        points={{-64,48},{-64,66},{-34,66}},
        color={107,45,134},
        thickness=1));
    connect(multiSum.y, UnmarkedATP.u[1]) annotation (Line(points={{87.02,46},{
            92,46},{92,23.05},{86,23.05}},
                                        color={0,0,127}));
    connect(DRX.y, UnmarkedATP.u[2]) annotation (Line(points={{87.02,-64},{92,
            -64},{92,20.95},{86,20.95}},
                                    color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false), graphics={Line(
            points={{-30,-90},{-30,-90},{-56,-90}},
            color={28,108,200},
            arrow={Arrow.None,Arrow.Filled})}),
      experiment(
        StartTime=-100,
        StopTime=100,
        __Dymola_NumberOfIntervals=5000,
        __Dymola_Algorithm="Dassl"));
  end DRXOnly;

  model DRXOnly_optim
    extends DRXOnly(
      DrxRate(displayUnit="l/min") = 3.00009E-06,
      RelaxedATPCycling(displayUnit="l/min") = 4.16667E-05,
      relaxingRate(displayUnit="l/min") = 2.66677E-05);
    extends Toepfer2020;

    Optimization.Criteria.Signals.IntegratedSquaredDeviation
      integratedSquaredDeviation
      annotation (Placement(transformation(extent={{42,60},{62,80}})));
    Modelica.Blocks.Sources.RealExpression realExpression2(y=ATPFluorescence)
      annotation (Placement(transformation(extent={{2,66},{22,86}})));
  equation
    connect(UT.q_out,kH. q_in) annotation (Line(
        points={{-50,-50},{-50,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(UT.solute,DRX. u[1]) annotation (Line(points={{-44,-60},{-44,-65.05},
            {74,-65.05}}, color={0,0,127}));
    connect(UD.solute,DRX. u[2]) annotation (Line(points={{36,-60},{36,-62.95},
            {74,-62.95}}, color={0,0,127}));
    connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
        Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
    connect(unlimitedSolutePumpOut1.q_in,UT. q_out) annotation (Line(
        points={{-80,-30},{-50,-30},{-50,-50}},
        color={107,45,134},
        thickness=1));
    connect(SRX.solute,multiSum. u[1]) annotation (Line(points={{-28,56},{-28,
            46},{74,46}},       color={0,0,127}));
    connect(Stream.q_in, kH.q_in) annotation (Line(
        points={{-50,28},{-50,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(Stream1.q_out, kH.q_in) annotation (Line(
        points={{-64,28},{-64,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(multiSum.y, UnmarkedATP.u[1]) annotation (Line(points={{87.02,46},{
            96,46},{96,23.05},{86,23.05}},
                                        color={0,0,127}));
    connect(DRX.y, UnmarkedATP.u[2]) annotation (Line(points={{87.02,-64},{96,
            -64},{96,20.95},{86,20.95}},
                                    color={0,0,127}));
    connect(integratedSquaredDeviation.u1, realExpression2.y)
      annotation (Line(points={{40,76},{23,76}}, color={0,0,127}));
    connect(UnmarkedATP.y, integratedSquaredDeviation.u2) annotation (Line(points
          ={{72.98,22},{8,22},{8,64},{40,64}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StartTime=-1000,
        StopTime=200,
        __Dymola_NumberOfIntervals=5000,
        __Dymola_Algorithm="Dassl"));
  end DRXOnly_optim;

  model Toepfer2020 "Data from PMID: 31983222"

    Real ATPFluorescence = a*exp(-b*x) + c*exp(-d*x);
    Real x = max(0, time);
    parameter Real a = 0.70, b = 0.052, c=0.30, d = 0.0061;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Toepfer2020;

  model DRXOnly_optimized
    extends DRXOnly(
      DrxRate(displayUnit="l/min") = 3.00009E-06*tune_a,
      RelaxedATPCycling(displayUnit="l/min") = 4.16667E-05*tune_b,
      relaxingRate(displayUnit="l/min") = 2.66677E-05*tune_c);
    extends Toepfer2020;
    parameter Real tune_a=2.0107438322077007,
                               tune_b=1.2262857238897689,
                                           tune_c=2.5880342477940874;
    Optimization.Criteria.Signals.IntegratedSquaredDeviation
      integratedSquaredDeviation
      annotation (Placement(transformation(extent={{42,60},{62,80}})));
    Modelica.Blocks.Sources.RealExpression realExpression2(y=ATPFluorescence)
      annotation (Placement(transformation(extent={{2,66},{22,86}})));
  equation
    connect(UT.q_out,kH. q_in) annotation (Line(
        points={{-50,-50},{-50,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(kH.q_out,UD. q_out) annotation (Line(
        points={{0,-30},{42,-30},{42,-50},{30,-50}},
        color={107,45,134},
        thickness=1));
    connect(UT.solute,DRX. u[1]) annotation (Line(points={{-44,-60},{-44,-65.05},
            {74,-65.05}}, color={0,0,127}));
    connect(UD.solute,DRX. u[2]) annotation (Line(points={{36,-60},{36,-62.95},
            {74,-62.95}}, color={0,0,127}));
    connect(XBCycling.q_out,XBCyclingMeasure. q_out) annotation (Line(
        points={{0,-80},{-4,-80},{-4,-72},{-6,-72},{-6,-88},{-20,-80}},
        color={107,45,134},
        thickness=1));
    connect(XBCycling.q_in,UD. q_out) annotation (Line(
        points={{20,-80},{56,-80},{56,-34},{42,-34},{42,-50},{30,-50}},
        color={107,45,134},
        thickness=1));
    connect(XBCyclingMeasure.q_in,UT. q_out) annotation (Line(
        points={{-40,-80},{-74,-80},{-74,-28},{-50,-28},{-50,-50}},
        color={107,45,134},
        thickness=1));
    connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
        Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
    connect(unlimitedSolutePumpOut1.q_in,UT. q_out) annotation (Line(
        points={{-80,-30},{-50,-30},{-50,-50}},
        color={107,45,134},
        thickness=1));
    connect(SRX.solute,multiSum. u[1]) annotation (Line(points={{-28,56},{-28,
            46},{74,46}},       color={0,0,127}));
    connect(Stream.q_out,SRX. q_out) annotation (Line(
        points={{-50,48},{-50,52},{-70,52},{-70,66},{-34,66}},
        color={107,45,134},
        thickness=1));
    connect(Stream.q_in, kH.q_in) annotation (Line(
        points={{-50,28},{-50,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(Stream1.q_out, kH.q_in) annotation (Line(
        points={{-64,28},{-64,-30},{-20,-30}},
        color={107,45,134},
        thickness=1));
    connect(Stream1.q_in, SRX.q_out) annotation (Line(
        points={{-64,48},{-94,48},{-94,66},{-34,66}},
        color={107,45,134},
        thickness=1));
    connect(multiSum.y, UnmarkedATP.u[1]) annotation (Line(points={{87.02,46},{
            96,46},{96,23.05},{86,23.05}},
                                        color={0,0,127}));
    connect(DRX.y, UnmarkedATP.u[2]) annotation (Line(points={{87.02,-64},{96,
            -64},{96,20.95},{86,20.95}},
                                    color={0,0,127}));
    connect(integratedSquaredDeviation.u1, realExpression2.y)
      annotation (Line(points={{40,76},{23,76}}, color={0,0,127}));
    connect(UnmarkedATP.y, integratedSquaredDeviation.u2) annotation (Line(points
          ={{72.98,22},{8,22},{8,64},{40,64}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StartTime=-1000,
        StopTime=200,
        __Dymola_NumberOfIntervals=5000,
        __Dymola_Algorithm="Dassl"));
  end DRXOnly_optimized;
  annotation (uses(Modelica(version="4.0.0"), Physiolibrary(version="2.4.1"),
      Optimization(version="2.2.6")));
end CrossBridgeCycling;
