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
  annotation (uses(Modelica(version="4.0.0")));
end CrossBridgeCycling;
