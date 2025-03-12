within ;
package CrossBridgeCycling
  package XB
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
      XB.elasticSegment elasticSegment1(
        c=ThickStiffness,
        cStiff1=HeadStiffness,
        att=4)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
      Modelica.Mechanics.Translational.Components.Fixed fixed
        annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      Modelica.Mechanics.Translational.Components.Spring spring(c=1)
        annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      XB.elasticSegment elasticSegment2(
        c=ThickStiffness,
        cStiff1=HeadStiffness,
        att=0) annotation (Placement(transformation(extent={{0,-10},{20,10}})));
      XB.elasticSegment elasticSegment3(
        c=ThickStiffness,
        cStiff1=HeadStiffness,
        att=0) annotation (Placement(transformation(extent={{32,-10},{52,10}})));
      Modelica.Mechanics.Translational.Sources.Position position
        annotation (Placement(transformation(extent={{104,-10},{84,10}})));
      Modelica.Blocks.Sources.Sine sine(f= 1)
        annotation (Placement(transformation(extent={{40,38},{60,58}})));
      parameter Modelica.Units.SI.TranslationalSpringConstant ThickStiffness=100
        "Spring constant";
      parameter Modelica.Units.SI.TranslationalSpringConstant HeadStiffness=10
        "Spring constant";
      XB.elasticSegment elasticSegment4(
        c=ThickStiffness,
        cStiff1=HeadStiffness,
        att=0) annotation (Placement(transformation(extent={{58,-10},{78,10}})));
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
  end XB;

  package mantATP
    package Data
      model Toepfer2020 "Data from PMID: 31983222"

        Real ATPFluorescence = a*exp(-b*x) + c*exp(-d*x);
        Real x = max(0, time);
        parameter Real a = 0.70, b = 0.052, c=0.30, d = 0.0061;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Toepfer2020;

      model ChaseData
        "Control cardiomyocyte mant-ATP chase. Data from experiments of Alison Van Der Roest and Julia Han."
        Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
          tableOnFile=false,
          table=[18,0.704575,0.626512,0.757024,1.175008; 27.994,0.628669,0.505251,
              0.70275,0.945396; 37.988,0.532319,0.405251,0.693936,0.873019; 47.982,
              0.519323,0.329185,0.666516,0.801315; 57.976,0.463446,0.245353,
              0.641676,0.831088; 67.97,0.443831,0.246881,0.592881,0.739397; 77.964,
              0.402628,0.178039,0.563653,0.704748; 87.958,0.347676,0.151591,
              0.522787,0.684761; 97.952,0.306668,0.12387,0.47,0.646423; 107.946,
              0.254855,0.111044,0.434256,0.644883; 117.94,0.222877,0.099013,
              0.402298,0.658037; 127.935,0.186469,0.047836,0.355838,0.627751;
              137.929,0.190436,0.052196,0.320904,0.564196; 147.923,0.191677,
              0.013144,0.295857,0.594931; 157.917,0.158652,0.009866,0.272712,
              0.540969; 167.911,0.150694,0.010948,0.252316,0.525024; 177.905,
              0.135556,0.017791,0.232335,0.524222; 187.899,0.143831,0.020146,
              0.222109,0.504267; 197.893,0.121392,-0.00987,0.213691,0.515143;
              207.887,0.11779,-0.00544,0.184426,0.509432; 217.881,0.105232,0.012826,
              0.180094,0.520404; 227.875,0.111073,0.018141,0.188456,0.490504;
              237.869,0.099732,0.012317,0.183013,0.44761; 247.863,0.107471,-0.01575,
              0.182335,0.438787; 257.857,0.104673,0.016836,0.149605,0.404299;
              267.851,0.090801,-0.02244,0.152015,0.427976; 277.845,0.0734,0.011362,
              0.134407,0.42204; 287.839,0.09944,-0.00328,0.125537,0.419313; 297.833,
              0.0734,0.007829,0.136817,0.358133; 307.827,0.08019,0.012731,0.128004,
              0.347; 317.822,0.068654,0.004042,0.145424,0.319731; 327.816,0.072937,
              0.030108,0.123804,0.319795; 337.81,0.074179,0.018491,0.125838,0.31325;
              347.804,0.077659,0.003405,0.121375,0.290889; 357.798,0.093332,-0.00719,
              0.112448,0.261373; 367.792,0.055245,0.002228,0.114444,0.294514;
              377.786,0.081747,0.012826,0.094539,0.264453; 387.78,0.060088,0.012158,
              0.090603,0.20847; 397.774,0.062521,0.006652,0.107458,0.158486;
              407.768,0.057532,-0.00414,0.091525,0.167918; 417.762,0.066488,
              0.016136,0.090282,0.118993; 427.756,0.047165,0.005124,0.097721,
              0.140488; 437.75,0.076856,0.00261,0.073879,0.128264; 447.744,0.044512,
              -0.00465,0.071224,0.120212; 457.738,0.042346,-0.00185,0.090056,
              0.091274; 467.732,0.060793,-0.00815,0.090207,0.10693; 477.726,
              0.045534,0.004742,0.10322,0.056561; 487.72,0.03816,0.009612,0.10484,
              0.079724; 497.714,0.03378,0.014577,0.0658,0.041931; 507.708,0.030129,
              0.009484,0.071073,0.035194; 517.703,0.03962,0.016136,0.049228,
              0.046198; 527.697,0.010416,-0.01181,0.057608,0.034039; 537.691,
              0.034826,-0.02139,0.037721,0.000834; 547.685,0.028304,0.000764,
              0.045669,0.009721; 557.679,0.00202,0.028517,0.05194,-0.01681; 567.673,
              0.034802,0.018014,0.047514,0.001412; 577.667,0.017985,0.009007,
              0.059831,-0.01344; 587.661,-0.00431,0.004901,0.029887,0.006224;
              597.655,0.024556,0.0155,0.038079,-0.00991; 607.649,0.021514,0.014386,
              0.03823,0.018896; 617.643,-0.00358,0.014449,0.021676,-0.0232],
          tableName="tab",
          fileName="data/tables.mat")
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Modelica.Blocks.Interfaces.RealOutput
                   y[combiTimeTable.nout] "Connector of Real output signals" annotation (Placement(
              transformation(extent={{94,-10},{114,10}})));
      equation
        connect(combiTimeTable.y, y)
          annotation (Line(points={{11,0},{104,0}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                                      Rectangle(
              extent={{-100,-100},{100,100}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              textColor={0,0,255}),
          Polygon(lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid,
            points={{-80,90},{-88,68},{-72,68},{-80,90}}),
          Line(points={{-80,68},{-80,-80}},
            color={192,192,192}),
          Line(points={{-90,-70},{82,-70}},
            color={192,192,192}),
          Polygon(lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid,
            points={{90,-70},{68,-62},{68,-78},{90,-70}}),
          Rectangle(lineColor={255,255,255},
            fillColor={244,125,35},
            fillPattern=FillPattern.Solid,
            extent={{-48,-50},{2,70}}),
          Line(points={{-48,-50},{-48,70},{52,70},{52,-50},{-48,-50},{-48,-20},{52,-20},
                    {52,10},{-48,10},{-48,40},{52,40},{52,70},{2,70},{2,-51}})}),
                                                                       Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end ChaseData;

      block TimeTable_ATPChaseControl1
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 18,0.704575;
              27.994,0.628669; 37.988,0.532319; 47.982,0.519323; 57.976,
              0.463446; 67.97,0.443831; 77.964,0.402628; 87.958,0.347676;
              97.952,0.306668; 107.946,0.254855; 117.94,0.222877; 127.935,
              0.186469; 137.929,0.190436; 147.923,0.191677; 157.917,0.158652;
              167.911,0.150694; 177.905,0.135556; 187.899,0.143831; 197.893,
              0.121392; 207.887,0.11779; 217.881,0.105232; 227.875,0.111073;
              237.869,0.099732; 247.863,0.107471; 257.857,0.104673; 267.851,
              0.090801; 277.845,0.0734; 287.839,0.09944; 297.833,0.0734;
              307.827,0.08019; 317.822,0.068654; 327.816,0.072937; 337.81,
              0.074179; 347.804,0.077659; 357.798,0.093332; 367.792,0.055245;
              377.786,0.081747; 387.78,0.060088; 397.774,0.062521; 407.768,
              0.057532; 417.762,0.066488; 427.756,0.047165; 437.75,0.076856;
              447.744,0.044512; 457.738,0.042346; 467.732,0.060793; 477.726,
              0.045534; 487.72,0.03816; 497.714,0.03378; 507.708,0.030129;
              517.703,0.03962; 527.697,0.010416; 537.691,0.034826; 547.685,
              0.028304; 557.679,0.00202; 567.673,0.034802; 577.667,0.017985;
              587.661,-0.00431; 597.655,0.024556; 607.649,0.021514; 617.643,-0.00358]);
      end TimeTable_ATPChaseControl1;

      block TimeTable_ATPChaseControl2
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 18,0.626512;
              27.994,0.505251; 37.988,0.405251; 47.982,0.329185; 57.976,
              0.245353; 67.97,0.246881; 77.964,0.178039; 87.958,0.151591;
              97.952,0.12387; 107.946,0.111044; 117.94,0.099013; 127.935,
              0.047836; 137.929,0.052196; 147.923,0.013144; 157.917,0.009866;
              167.911,0.010948; 177.905,0.017791; 187.899,0.020146; 197.893,-0.00987;
              207.887,-0.00544; 217.881,0.012826; 227.875,0.018141; 237.869,
              0.012317; 247.863,-0.01575; 257.857,0.016836; 267.851,-0.02244;
              277.845,0.011362; 287.839,-0.00328; 297.833,0.007829; 307.827,
              0.012731; 317.822,0.004042; 327.816,0.030108; 337.81,0.018491;
              347.804,0.003405; 357.798,-0.00719; 367.792,0.002228; 377.786,
              0.012826; 387.78,0.012158; 397.774,0.006652; 407.768,-0.00414;
              417.762,0.016136; 427.756,0.005124; 437.75,0.00261; 447.744,-0.00465;
              457.738,-0.00185; 467.732,-0.00815; 477.726,0.004742; 487.72,
              0.009612; 497.714,0.014577; 507.708,0.009484; 517.703,0.016136;
              527.697,-0.01181; 537.691,-0.02139; 547.685,0.000764; 557.679,
              0.028517; 567.673,0.018014; 577.667,0.009007; 587.661,0.004901;
              597.655,0.0155; 607.649,0.014386; 617.643,0.014449]);
      end TimeTable_ATPChaseControl2;

      block TimeTable_ATPChaseControl3
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 18,0.757024;
              27.994,0.70275; 37.988,0.693936; 47.982,0.666516; 57.976,0.641676;
              67.97,0.592881; 77.964,0.563653; 87.958,0.522787; 97.952,0.47;
              107.946,0.434256; 117.94,0.402298; 127.935,0.355838; 137.929,
              0.320904; 147.923,0.295857; 157.917,0.272712; 167.911,0.252316;
              177.905,0.232335; 187.899,0.222109; 197.893,0.213691; 207.887,
              0.184426; 217.881,0.180094; 227.875,0.188456; 237.869,0.183013;
              247.863,0.182335; 257.857,0.149605; 267.851,0.152015; 277.845,
              0.134407; 287.839,0.125537; 297.833,0.136817; 307.827,0.128004;
              317.822,0.145424; 327.816,0.123804; 337.81,0.125838; 347.804,
              0.121375; 357.798,0.112448; 367.792,0.114444; 377.786,0.094539;
              387.78,0.090603; 397.774,0.107458; 407.768,0.091525; 417.762,
              0.090282; 427.756,0.097721; 437.75,0.073879; 447.744,0.071224;
              457.738,0.090056; 467.732,0.090207; 477.726,0.10322; 487.72,
              0.10484; 497.714,0.0658; 507.708,0.071073; 517.703,0.049228;
              527.697,0.057608; 537.691,0.037721; 547.685,0.045669; 557.679,
              0.05194; 567.673,0.047514; 577.667,0.059831; 587.661,0.029887;
              597.655,0.038079; 607.649,0.03823; 617.643,0.021676]);
      end TimeTable_ATPChaseControl3;

      block TimeTable_ATPChaseControl4
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 18,0.704575;
              27.994,0.628669; 37.988,0.532319; 47.982,0.519323; 57.976,
              0.463446; 67.97,0.443831; 77.964,0.402628; 87.958,0.347676;
              97.952,0.306668; 107.946,0.254855; 117.94,0.222877; 127.935,
              0.186469; 137.929,0.190436; 147.923,0.191677; 157.917,0.158652;
              167.911,0.150694; 177.905,0.135556; 187.899,0.143831; 197.893,
              0.121392; 207.887,0.11779; 217.881,0.105232; 227.875,0.111073;
              237.869,0.099732; 247.863,0.107471; 257.857,0.104673; 267.851,
              0.090801; 277.845,0.0734; 287.839,0.09944; 297.833,0.0734;
              307.827,0.08019; 317.822,0.068654; 327.816,0.072937; 337.81,
              0.074179; 347.804,0.077659; 357.798,0.093332; 367.792,0.055245;
              377.786,0.081747; 387.78,0.060088; 397.774,0.062521; 407.768,
              0.057532; 417.762,0.066488; 427.756,0.047165; 437.75,0.076856;
              447.744,0.044512; 457.738,0.042346; 467.732,0.060793; 477.726,
              0.045534; 487.72,0.03816; 497.714,0.03378; 507.708,0.030129;
              517.703,0.03962; 527.697,0.010416; 537.691,0.034826; 547.685,
              0.028304; 557.679,0.00202; 567.673,0.034802; 577.667,0.017985;
              587.661,-0.00431; 597.655,0.024556; 607.649,0.021514; 617.643,-0.00358]);
      end TimeTable_ATPChaseControl4;

      block TimeTable_ATPChaseMeanPos "Mean of all 10 POS data from ATP chase"
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 5.529,0.927350434;
              11.058,0.863017052; 16.587,0.805839694; 22.116,0.756884973;
              27.645,0.715639075; 33.174,0.678791754; 38.703,0.646107701;
              44.232,0.617961965; 49.761,0.592224941; 55.29,0.569170047; 60.819,
              0.548339118; 66.348,0.530670789; 71.877,0.515357682; 77.406,
              0.499134294; 82.935,0.485184333; 88.464,0.472883152; 93.993,
              0.461307102; 99.523,0.45037915; 105.052,0.438796683; 110.581,
              0.427954043; 116.11,0.418175685; 121.639,0.409314261; 127.168,
              0.400305015; 132.697,0.391518023; 138.226,0.384122444; 143.755,
              0.378596225; 149.284,0.372118237; 154.813,0.364638911; 160.342,
              0.357411109; 165.871,0.34909594; 171.4,0.341757209; 176.929,
              0.336437754; 182.458,0.328506074; 187.987,0.322778604; 193.516,
              0.316864055; 199.045,0.312173634; 204.574,0.307026678; 210.103,
              0.30187899; 215.632,0.295187427; 221.161,0.288818556; 226.69,
              0.283511288; 232.219,0.279540521; 237.748,0.273887428; 243.277,
              0.267194997; 248.806,0.262483913; 254.335,0.256754379; 259.864,
              0.2530259; 265.393,0.247998756; 270.922,0.242715944; 276.451,
              0.237986074; 281.98,0.233339736; 287.509,0.227348635; 293.038,
              0.222794651; 298.568,0.217815087; 304.097,0.212747774; 309.626,
              0.208537212; 315.155,0.204386267; 320.684,0.200575872; 326.213,
              0.195764734; 331.742,0.190628259; 337.271,0.186022245; 342.8,
              0.181696196; 348.329,0.177508906; 353.858,0.172760629; 359.387,
              0.168569088; 364.916,0.164026577; 370.445,0.159551745; 375.974,
              0.155578284; 381.503,0.151677261; 387.032,0.148467245; 392.561,
              0.143711123; 398.09,0.139688572; 403.619,0.136124753; 409.148,
              0.131243394; 414.677,0.125071105; 420.206,0.120561432; 425.735,
              0.116331157; 431.264,0.112238916; 436.793,0.108335361; 442.322,
              0.104607132; 447.851,0.102025517; 453.38,0.09779329; 458.909,
              0.09435339; 464.438,0.091729225; 469.967,0.087877381; 475.496,
              0.083563054; 481.025,0.078837463; 486.554,0.074527306; 492.083,
              0.071352958; 497.613,0.066395998; 503.142,0.062387083; 508.671,
              0.059136577; 514.2,0.05550326; 519.729,0.051741794; 525.258,
              0.048054276; 530.787,0.044747722; 536.316,0.040591287; 541.845,
              0.036577139; 547.374,0.033348047; 552.903,0.029176759; 558.432,
              0.024649809; 563.961,0.020267953; 569.49,0.017217835; 575.019,
              0.013529917; 580.548,0.010135699; 586.077,0.00682139; 591.606,
              0.003805418; 597.135,0]);
      end TimeTable_ATPChaseMeanPos;

      block TimeTable_ATPChaseMavaMeanPos
        "Mean of all four POS data from ATP chase"
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 7.818,0.913680654;
              15.636,0.841023599; 23.455,0.777326996; 31.273,0.723555152;
              39.091,0.677882253; 46.909,0.637980255; 54.727,0.603463261;
              62.545,0.572208451; 70.364,0.544891351; 78.182,0.521108872; 86,
              0.498271314; 93.818,0.479165956; 101.636,0.461753789; 109.454,
              0.447710155; 117.273,0.429692055; 125.091,0.414403914; 132.909,
              0.401354049; 140.727,0.387305477; 148.545,0.374654799; 156.363,
              0.36443113; 164.182,0.353251496; 172,0.341742903; 179.818,
              0.332276995; 187.636,0.322722698; 195.454,0.315959531; 203.272,
              0.311113656; 211.091,0.301055514; 218.909,0.293461494; 226.727,
              0.288523875; 234.545,0.276696942; 242.363,0.270475903; 250.181,
              0.263371779; 258,0.255346822; 265.818,0.251000742; 273.636,
              0.243499858; 281.454,0.238475408; 289.272,0.232521991; 297.09,
              0.224779137; 304.909,0.221336147; 312.727,0.212461788; 320.545,
              0.206149839; 328.363,0.20028084; 336.181,0.195653179; 344,
              0.189890889; 351.818,0.184592387; 359.636,0.180791237; 367.454,
              0.178934881; 375.272,0.171880773; 383.09,0.16578103; 390.909,
              0.158030529; 398.727,0.153419729; 406.545,0.148172318; 414.363,
              0.140240898; 422.181,0.134548101; 429.999,0.128192015; 437.818,
              0.119054407; 445.636,0.113645281; 453.454,0.106140443; 461.272,
              0.100521261; 469.09,0.094313316; 476.908,0.086150135; 484.727,
              0.080038679; 492.545,0.075432478; 500.363,0.06983712; 508.181,
              0.064027152; 515.999,0.057104891; 523.817,0.05264343; 531.636,
              0.045766869; 539.454,0.039557303; 547.272,0.033853193; 555.09,
              0.026413981; 562.908,0.021942752; 570.727,0.01773855; 578.545,
              0.010838229; 586.363,0.005683596; 594.181,0]);
      end TimeTable_ATPChaseMavaMeanPos;

      block TimeTable_NormATPChaseMava1
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 5.668,0.946332737;
              11.335,0.940780951; 17.003,0.8170995; 22.67,0.821664302; 28.338,
              0.782061563; 34.005,0.67873666; 39.673,0.624884338; 45.34,
              0.645981124; 51.008,0.535870705; 56.675,0.560051817; 62.343,
              0.566097095; 68.01,0.496823145; 73.678,0.440071556; 79.345,
              0.464992906; 85.013,0.427364135; 90.68,0.438652767; 96.348,
              0.391215841; 102.015,0.404355068; 107.683,0.421133798; 113.35,
              0.421874036; 119.018,0.408796496; 124.685,0.450743322; 130.353,
              0.35130467; 136.02,0.345321078; 141.688,0.373635186; 147.355,
              0.381469373; 153.023,0.344272408; 158.69,0.334402566; 164.358,
              0.335883042; 170.025,0.308494232; 175.693,0.281783974; 181.36,
              0.304669669; 187.028,0.32903584; 192.695,0.277157486; 198.363,
              0.274319906; 204.031,0.264018259; 209.698,0.269816791; 215.366,
              0.310283141; 221.033,0.256184073; 226.701,0.282400839; 232.368,
              0.254333477; 238.036,0.198383813; 243.703,0.244710382; 249.371,
              0.248041453; 255.038,0.152920856; 260.706,0.213682068; 266.373,
              0.227623219; 272.041,0.214607365; 277.708,0.246314231; 283.376,
              0.184751095; 289.043,0.220837703; 294.711,0.20868546; 300.378,
              0.228671889; 306.046,0.202085004; 311.713,0.179877861; 317.381,
              0.181234964; 323.048,0.125717106; 328.716,0.212201592; 334.383,
              0.179939547; 340.051,0.205786195; 345.718,0.112577879; 351.386,
              0.149219666; 357.053,0.152612424; 362.721,0.146505459; 368.388,
              0.175004626; 374.056,0.214237246; 379.723,0.169946333; 385.391,
              0.155264944; 391.059,0.092344704; 396.726,0.0598976; 402.394,
              0.176485103; 408.061,0.163839368; 413.729,0.190056135; 419.396,
              0.16717044; 425.064,0.060020973; 430.731,0.102893097; 436.399,
              0.058232065; 442.066,0.08728641; 447.734,0.0077725; 453.401,
              0.019677996; 459.069,0.069767442; 464.736,0.050521251; 470.404,
              0.005798532; 476.071,0.057368454; 481.739,0.066991549; 487.406,
              0.054099069; 493.074,0.059650854; 498.741,0.071988156; 504.409,
              0.04028129; 510.076,0.078588613; 515.744,0.078403553; 521.411,
              0.07050768; 527.079,0.075195855; 532.746,0.079884029; 538.414,
              0.071247918; 544.081,0.101597681; 549.749,0.03392758; 555.416,
              0.126889149; 561.084,0.063043612; 566.751,0.009376349; 572.419,
              0.043674048; 578.086,0.027635556; 583.754,0; 589.422,0.044044168;
              595.089,0.038307322]);
      end TimeTable_NormATPChaseMava1;

      block TimeTable_NormATPChaseControl1
        extends Modelica.Blocks.Sources.TimeTable(table=[0,1; 18,0.704575;
              27.994,0.628669; 37.988,0.532319; 47.982,0.519323; 57.976,
              0.463446; 67.97,0.443831; 77.964,0.402628; 87.958,0.347676;
              97.952,0.306668; 107.946,0.254855; 117.94,0.222877; 127.935,
              0.186469; 137.929,0.190436; 147.923,0.191677; 157.917,0.158652;
              167.911,0.150694; 177.905,0.135556; 187.899,0.143831; 197.893,
              0.121392; 207.887,0.11779; 217.881,0.105232; 227.875,0.111073;
              237.869,0.099732; 247.863,0.107471; 257.857,0.104673; 267.851,
              0.090801; 277.845,0.0734; 287.839,0.09944; 297.833,0.0734;
              307.827,0.08019; 317.822,0.068654; 327.816,0.072937; 337.81,
              0.074179; 347.804,0.077659; 357.798,0.093332; 367.792,0.055245;
              377.786,0.081747; 387.78,0.060088; 397.774,0.062521; 407.768,
              0.057532; 417.762,0.066488; 427.756,0.047165; 437.75,0.076856;
              447.744,0.044512; 457.738,0.042346; 467.732,0.060793; 477.726,
              0.045534; 487.72,0.03816; 497.714,0.03378; 507.708,0.030129;
              517.703,0.03962; 527.697,0.010416; 537.691,0.034826; 547.685,
              0.028304; 557.679,0.00202; 567.673,0.034802; 577.667,0.017985;
              587.661,-0.00431; 597.655,0.024556; 607.649,0.021514; 617.643,-0.00358]);
      end TimeTable_NormATPChaseControl1;
    end Data;

    package Simple
      model SrxDrx
        Bodylight.Chemical.Components.Substance DRX(solute_start=1 - SRX_init)
          annotation (Placement(transformation(extent={{-62,-22},{-42,-2}})));
        Bodylight.Chemical.Components.Substance SRX(solute_start=SRX_init)
          annotation (Placement(transformation(extent={{-60,40},{-40,60}})));
        Modelica.Blocks.Math.MultiSum multiSum(nu=2)
          annotation (Placement(transformation(extent={{24,28},{36,16}})));
        Bodylight.Chemical.Components.Clearance clearance(Clearance(displayUnit
              ="l/min") = 1.6666666666667e-05)
          annotation (Placement(transformation(extent={{-32,70},{-12,90}})));
        Bodylight.Chemical.Components.Clearance clearance1(Clearance(
              displayUnit="ml/min") = 0.000166667)
          annotation (Placement(transformation(extent={{-108,-22},{-128,-2}})));
        parameter Bodylight.Types.AmountOfSubstance SRX_init(displayUnit="mol")
           = 0.5 "Initial solute amount in compartment";
        parameter Bodylight.Types.VolumeFlowRate SolutionFlow=0
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
        extends mantATP.Simple.SrxDrx(
          SolutionFlow=0,
          clearance1(Clearance(displayUnit="l/min") = DRXClearance,
              useSolutionFlowInput=true),
          clearance(Clearance=SRXClearance, useSolutionFlowInput=true));
        Bodylight.Chemical.Components.Stream Stream(SolutionFlow=relaxingRate)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-66,22})));
        Bodylight.Chemical.Components.Stream Stream1(SolutionFlow=DrxRate)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-86,22})));
        Modelica.Blocks.Sources.RealExpression realExpression(y=if time > 0 then
              DRXClearance else 0)
          annotation (Placement(transformation(extent={{-156,-2},{-136,18}})));
        parameter Bodylight.Types.VolumeFlowRate DRXClearance(displayUnit=
              "l/min") = 0.000166667
          "Clearance of solute if useSolutionFlowInput=false";
        Modelica.Blocks.Sources.RealExpression realExpression2(y=if time > 0 then
              SRXClearance else 0)
          annotation (Placement(transformation(extent={{16,78},{-4,98}})));
        parameter Bodylight.Types.VolumeFlowRate SRXClearance(displayUnit=
              "l/min") = 1.6666666666667e-05
          "Clearance of solute if useSolutionFlowInput=false";
        parameter Bodylight.Types.VolumeFlowRate relaxingRate(displayUnit=
              "l/min") = 8.3333333333333e-07
          "Volumetric flow of solution if useSolutionFlowInput=false";
        parameter Bodylight.Types.VolumeFlowRate DrxRate(displayUnit="l/min")
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
    end Simple;

    package SRX_ADP

      model UtUdSrtSrd
        Bodylight.Chemical.Components.Substance UT(solute_start=0.5 - SRX_init/
              2)
          annotation (Placement(transformation(extent={{-60,-80},{-40,-60}})));
        Bodylight.Chemical.Components.Substance UD(solute_start=0.5 - SRX_init/
              2) annotation (Placement(transformation(extent={{38,-80},{58,-60}})));
        parameter Bodylight.Types.AmountOfSubstance SRX_init(displayUnit="mol")
           = 0.5 "Initial solute amount in compartment";
        parameter Bodylight.Types.VolumeFlowRate SolutionFlow=0
          "Volumetric flow of solution if useSolutionFlowInput=false";
        Bodylight.Chemical.Components.Substance SD(solute_start=SRX_init/2)
          annotation (Placement(transformation(extent={{40,74},{60,94}})));
        Bodylight.Chemical.Components.Substance ST(solute_start=SRX_init/2)
          annotation (Placement(transformation(extent={{-62,76},{-42,96}})));
        Bodylight.Chemical.Components.Stream kH(SolutionFlow(displayUnit=
                "l/min") = 4.1666666666667e-06)
          annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
        Bodylight.Chemical.Components.Stream kS2D(SolutionFlow(displayUnit=
                "l/min") = 8.3333333333333e-05) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={50,-10})));
        Bodylight.Chemical.Components.Stream kL(SolutionFlow(displayUnit=
                "l/min") = 8.3333333333333e-05)
          annotation (Placement(transformation(extent={{40,40},{20,60}})));
        Bodylight.Chemical.Components.Stream KS2T(SolutionFlow(displayUnit=
                "l/min") = 0.00016666666666667) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-50,30})));
        Modelica.Blocks.Math.MultiSum SRX(nu=2)
          annotation (Placement(transformation(extent={{80,54},{92,66}})));
        Modelica.Blocks.Math.MultiSum DRX(nu=2)
          annotation (Placement(transformation(extent={{82,-96},{94,-84}})));
        Bodylight.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure
          annotation (Placement(transformation(extent={{-20,60},{0,40}})));
        Bodylight.Chemical.Sources.UnlimitedSolutePumpOut
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
        extends mantATP.SRX_ADP.UtUdSrtSrd(
                           kS2D(SolutionFlow=8.3333333333333e-07), KS2T(
              SolutionFlow=3.3333333333333e-06));
        annotation (experiment(
            StartTime=-2500,
            StopTime=2500,
            __Dymola_Algorithm="Dassl"));
      end UtUdSrtSrd_V1;

      model UtUdSrtSrd_Ud2Ut "Reverse hydrolysis"
        extends mantATP.SRX_ADP.UtUdSrtSrd(
                           kS2D(SolutionFlow=1.6666666666667e-07), KS2T(
              SolutionFlow=1.6666666666667e-07));
        Bodylight.Chemical.Components.Stream Stream4(SolutionFlow(displayUnit=
                "l/min") = 1.6666666666667e-06)
          annotation (Placement(transformation(extent={{24,-66},{4,-46}})));
        Bodylight.Chemical.Sensors.MolarFlowMeasure molarFlowMeasure1
          annotation (Placement(transformation(extent={{-32,-46},{-12,-66}})));
        Bodylight.Chemical.Components.Stream Stream5(SolutionFlow(displayUnit=
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
        extends mantATP.SRX_ADP.UtUdSrtSrd(
                           kS2D(SolutionFlow=1.6666666666667e-06));
        Bodylight.Chemical.Components.Stream KS1T(SolutionFlow(displayUnit=
                "l/min") = 1.6666666666667e-06) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-70,30})));
        Bodylight.Chemical.Components.Stream KS1D(SolutionFlow(displayUnit=
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
        extends mantATP.SRX_ADP.UtUdSrtSrd;
        Bodylight.Chemical.Components.Stream XBCycling(SolutionFlow(displayUnit
              ="l/min") = 1.6666666666667e-06)
          annotation (Placement(transformation(extent={{20,-120},{0,-100}})));
        Bodylight.Chemical.Sensors.MolarFlowMeasure XBCyclingMeasure
          annotation (Placement(transformation(extent={{-40,-100},{-20,-120}})));
        Bodylight.Chemical.Sources.UnlimitedSolutePumpOut
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
        extends mantATP.SRX_ADP.UtUdSrtSrd_XBCycling(
          kL(SolutionFlow=1.6666666666667e-06),
          kH(SolutionFlow=0.00016666666666667),
          XBCycling(SolutionFlow=1.6666666666667e-05),
          kS2D(SolutionFlow=1.6666666666667e-06),
          KS2T(SolutionFlow=1.6666666666667e-06));
        Bodylight.Chemical.Components.Stream KS1T(SolutionFlow(displayUnit=
                "l/min") = 0)                   annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-70,30})));
        Bodylight.Chemical.Components.Stream KS1D(SolutionFlow(displayUnit=
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
    end SRX_ADP;

    package SRX_ATP
      model Schematics
        Bodylight.Chemical.Components.Substance DRX_ATP(solute_start=0.5 -
              SRX_init/2)
          annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
        Bodylight.Chemical.Components.Substance DRX_ADP(solute_start=0.5 -
              SRX_init/2)
          annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
        Bodylight.Chemical.Components.Stream kH(SolutionFlow(displayUnit=
                "l/min") = 0.0016666666666667)
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
        Modelica.Blocks.Math.MultiSum DRX_all(nu=2)
          annotation (Placement(transformation(extent={{44,-60},{56,-48}})));
        Bodylight.Chemical.Components.Stream XBCycling(SolutionFlow(displayUnit
              ="l/min") = RelaxedATPCycling)
          annotation (Placement(transformation(extent={{10,-14},{-10,6}})));
        Bodylight.Chemical.Components.Substance SRX(solute_start=SRX_init)
          annotation (Placement(transformation(extent={{-122,-50},{-102,-30}})));
        Bodylight.Chemical.Components.Stream ksr(SolutionFlow=relaxingRate)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-80,-4})));
        Bodylight.Chemical.Components.Stream kmsr(SolutionFlow=DrxRate)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-80,-20})));
        Bodylight.Chemical.Components.Substance A2_Post_Ratchetted(solute_start
            =SRX_init)
          annotation (Placement(transformation(extent={{-50,-90},{-30,-70}})));
      equation
        connect(DRX_ATP.q_out, kH.q_in) annotation (Line(
            points={{-40,-40},{-40,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        connect(kH.q_out, DRX_ADP.q_out) annotation (Line(
            points={{10,-20},{30,-20},{30,-40}},
            color={107,45,134},
            thickness=1));
        connect(DRX_ATP.solute, DRX_all.u[1]) annotation (Line(points={{-34,-50},
                {-34,-55.05},{44,-55.05}}, color={0,0,127}));
        connect(DRX_ADP.solute, DRX_all.u[2]) annotation (Line(points={{36,-50},
                {36,-52.95},{44,-52.95}}, color={0,0,127}));
        connect(XBCycling.q_in, DRX_ADP.q_out) annotation (Line(
            points={{10,-4},{30,-4},{30,-40}},
            color={107,45,134},
            thickness=1));
        connect(ksr.q_out, SRX.q_out) annotation (Line(
            points={{-90,-4},{-112,-4},{-112,-40}},
            color={107,45,134},
            thickness=1));
        connect(kmsr.q_out, kH.q_in) annotation (Line(
            points={{-70,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        connect(kmsr.q_in, SRX.q_out) annotation (Line(
            points={{-90,-20},{-112,-20},{-112,-40}},
            color={107,45,134},
            thickness=1));
        connect(XBCycling.q_out, DRX_ATP.q_out) annotation (Line(
            points={{-10,-4},{-40,-4},{-40,-40}},
            color={107,45,134},
            thickness=1));
        connect(ksr.q_in, kH.q_in) annotation (Line(
            points={{-70,-4},{-40,-4},{-40,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Line(
                points={{-2,4},{-26,4},{-40,20}},
                color={28,108,200},
                arrow={Arrow.None,Arrow.Filled},
                thickness=1,
                smooth=Smooth.Bezier),                               Line(
                points={{-42,-74},{-42,-64},{-42,-50}},
                color={0,140,72},
                arrow={Arrow.None,Arrow.Filled},
                thickness=1,
                smooth=Smooth.Bezier),                               Line(
                points={{-50,-68},{-78,-50},{-96,-38}},
                color={244,125,35},
                arrow={Arrow.None,Arrow.Filled},
                thickness=1,
                smooth=Smooth.Bezier),
              Text(
                extent={{-26,-68},{42,-94}},
                textColor={244,125,35},
                textString="A2 
post-ratchetted",
                horizontalAlignment=TextAlignment.Left)}));
      end Schematics;

      model SchematicsAlternative
        Bodylight.Chemical.Components.Substance DRX_ATP(solute_start=0.5 -
              SRX_init/2)
          annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
        Bodylight.Chemical.Components.Substance DRX_ADP(solute_start=0.5 -
              SRX_init/2)
          annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
        Bodylight.Chemical.Components.Stream kH(SolutionFlow(displayUnit=
                "l/min") = 0.0016666666666667)
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
        Modelica.Blocks.Math.MultiSum DRX_all(nu=2)
          annotation (Placement(transformation(extent={{44,-60},{56,-48}})));
        Bodylight.Chemical.Components.Stream XBCycling(SolutionFlow(displayUnit
              ="l/min") = RelaxedATPCycling)
          annotation (Placement(transformation(extent={{10,-14},{-10,6}})));
        Bodylight.Chemical.Components.Substance SRX_instant(solute_start=
              SRX_init) annotation (Placement(transformation(extent={{-122,-50},
                  {-102,-30}})));
        Bodylight.Chemical.Components.Stream ksr(SolutionFlow=relaxingRate)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-80,-4})));
        Bodylight.Chemical.Components.Stream kmsr(SolutionFlow=DrxRate)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-80,-20})));
        Bodylight.Chemical.Components.Substance SRX(solute_start=SRX_init)
          annotation (Placement(transformation(extent={{-122,6},{-102,26}})));
        Bodylight.Chemical.Components.Stream ksr1(SolutionFlow=relaxingRate)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-84,52})));
        Bodylight.Chemical.Components.Stream kmsr1(SolutionFlow=DrxRate)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-84,36})));
        Bodylight.Chemical.Components.Substance A2_Post_Ratchetted(solute_start
            =SRX_init)
          annotation (Placement(transformation(extent={{-48,-92},{-28,-72}})));
      equation
        connect(DRX_ATP.q_out, kH.q_in) annotation (Line(
            points={{-40,-40},{-40,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        connect(kH.q_out, DRX_ADP.q_out) annotation (Line(
            points={{10,-20},{30,-20},{30,-40}},
            color={107,45,134},
            thickness=1));
        connect(DRX_ATP.solute, DRX_all.u[1]) annotation (Line(points={{-34,-50},
                {-34,-55.05},{44,-55.05}}, color={0,0,127}));
        connect(DRX_ADP.solute, DRX_all.u[2]) annotation (Line(points={{36,-50},
                {36,-52.95},{44,-52.95}}, color={0,0,127}));
        connect(XBCycling.q_in, DRX_ADP.q_out) annotation (Line(
            points={{10,-4},{30,-4},{30,-40}},
            color={107,45,134},
            thickness=1));
        connect(ksr.q_out, SRX_instant.q_out) annotation (Line(
            points={{-90,-4},{-112,-4},{-112,-40}},
            color={107,45,134},
            thickness=1));
        connect(kmsr.q_out, kH.q_in) annotation (Line(
            points={{-70,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        connect(kmsr.q_in, SRX_instant.q_out) annotation (Line(
            points={{-90,-20},{-112,-20},{-112,-40}},
            color={107,45,134},
            thickness=1));
        connect(XBCycling.q_out, DRX_ATP.q_out) annotation (Line(
            points={{-10,-4},{-40,-4},{-40,-40}},
            color={107,45,134},
            thickness=1));
        connect(ksr.q_in, kH.q_in) annotation (Line(
            points={{-70,-4},{-40,-4},{-40,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        connect(ksr1.q_out, SRX.q_out) annotation (Line(
            points={{-94,52},{-112,52},{-112,16}},
            color={107,45,134},
            thickness=1));
        connect(kmsr1.q_out, kH.q_in) annotation (Line(
            points={{-74,36},{-40,36},{-40,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        connect(kmsr1.q_in, SRX.q_out) annotation (Line(
            points={{-94,36},{-112,36},{-112,16}},
            color={107,45,134},
            thickness=1));
        connect(ksr1.q_in, kH.q_in) annotation (Line(
            points={{-74,52},{-40,52},{-40,-20},{-10,-20}},
            color={107,45,134},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false), graphics={Line(
                points={{2,2},{-14,8},{-24,24}},
                color={28,108,200},
                arrow={Arrow.None,Arrow.Filled},
                thickness=1,
                smooth=Smooth.Bezier),
              Text(
                extent={{-24,-70},{44,-96}},
                textColor={244,125,35},
                textString="A2 
post-ratchetted",
                horizontalAlignment=TextAlignment.Left),             Line(
                points={{-40,-76},{-40,-66},{-40,-52}},
                color={0,140,72},
                arrow={Arrow.None,Arrow.Filled},
                thickness=1,
                smooth=Smooth.Bezier),                               Line(
                points={{-48,-70},{-76,-52},{-94,-40}},
                color={244,125,35},
                arrow={Arrow.None,Arrow.Filled},
                thickness=1,
                smooth=Smooth.Bezier)}));
      end SchematicsAlternative;

      model DRXOnly
        Bodylight.Chemical.Components.Substance DRX_ATP(solute_start=0.5 -
              SRX_init/2)
          annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
        Bodylight.Chemical.Components.Substance DRX_ADP(solute_start=0.5 -
              SRX_init/2)
          annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
        Bodylight.Chemical.Components.Stream Stream_kH(SolutionFlow(displayUnit
              ="l/min") = kH)
          annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
        Modelica.Blocks.Math.MultiSum DRX_all(nu=2)
          annotation (Placement(transformation(extent={{74,-70},{86,-58}})));
        Bodylight.Chemical.Components.Stream XBCycling(SolutionFlow(displayUnit
              ="l/min") = RelaxedATPCycling)
          annotation (Placement(transformation(extent={{20,-90},{0,-70}})));
        Bodylight.Chemical.Sensors.MolarFlowMeasure XBCyclingMeasure
          annotation (Placement(transformation(extent={{-40,-90},{-20,-70}})));
        Bodylight.Chemical.Sources.UnlimitedSolutePumpOut
          unlimitedSolutePumpOut1(useSoluteFlowInput=true)
          annotation (Placement(transformation(extent={{-80,-20},{-100,-40}})));
        Modelica.Blocks.Sources.RealExpression realExpression1(y=if time > decay_time
               then -XBCyclingMeasure.molarFlowRate else 0)
          annotation (Placement(transformation(extent={{-60,-100},{-80,-80}})));
        Modelica.Blocks.Math.MultiSum UnmarkedATP(nu=2)
          annotation (Placement(transformation(extent={{86,28},{74,16}})));
        Bodylight.Chemical.Components.Substance SRX(solute_start=SRX_init)
          annotation (Placement(transformation(extent={{-44,56},{-24,76}})));
        Modelica.Blocks.Math.MultiSum SRX_all(nu=1)
          annotation (Placement(transformation(extent={{74,52},{86,40}})));
        Bodylight.Chemical.Components.Stream Stream(SolutionFlow=relaxingRate)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-50,38})));
        Bodylight.Chemical.Components.Stream Stream1(SolutionFlow=DrxRate)
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-64,38})));
        parameter Bodylight.Types.AmountOfSubstance SRX_init(displayUnit="mol")
          =0.5   "Initial solute amount in compartment";
        parameter Bodylight.Types.VolumeFlowRate SolutionFlow=0
          "Volumetric flow of solution if useSolutionFlowInput=false";
        parameter Modelica.Blocks.Interfaces.RealOutput decay_time=0.0
          "Value of Real output";
        parameter Bodylight.Types.VolumeFlowRate DrxRate(displayUnit="l/min")=
          1.6666666666667e-06
          "Volumetric flow of solution if useSolutionFlowInput=false";
        parameter Bodylight.Types.VolumeFlowRate relaxingRate(displayUnit=
              "l/min")=8.3333333333333e-07
          "Volumetric flow of solution if useSolutionFlowInput=false";

        parameter Bodylight.Types.VolumeFlowRate RelaxedATPCycling(displayUnit=
              "l/min") = 0.00016666666666667
          "Volumetric flow of solution if useSolutionFlowInput=false";
        parameter Bodylight.Types.VolumeFlowRate kH(displayUnit="l/min") =
          0.0016666666666667
          "Volumetric flow of solution if useSolutionFlowInput=false";
      equation
        connect(DRX_ATP.q_out, Stream_kH.q_in) annotation (Line(
            points={{-50,-50},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream_kH.q_out, DRX_ADP.q_out) annotation (Line(
            points={{0,-30},{30,-30},{30,-50}},
            color={107,45,134},
            thickness=1));
        connect(DRX_ATP.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},
                {-44,-65.05},{74,-65.05}}, color={0,0,127}));
        connect(DRX_ADP.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},
                {36,-62.95},{74,-62.95}}, color={0,0,127}));
        connect(XBCycling.q_out,XBCyclingMeasure. q_out) annotation (Line(
            points={{0,-80},{-20,-80}},
            color={107,45,134},
            thickness=1));
        connect(XBCycling.q_in, DRX_ADP.q_out) annotation (Line(
            points={{20,-80},{56,-80},{56,-30},{30,-30},{30,-50}},
            color={107,45,134},
            thickness=1));
        connect(XBCyclingMeasure.q_in, DRX_ATP.q_out) annotation (Line(
            points={{-40,-80},{-64,-80},{-64,-30},{-50,-30},{-50,-50}},
            color={107,45,134},
            thickness=1));
        connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
            Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
        connect(unlimitedSolutePumpOut1.q_in, DRX_ATP.q_out) annotation (Line(
            points={{-80,-30},{-50,-30},{-50,-50}},
            color={107,45,134},
            thickness=1));
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
        connect(Stream.q_out,SRX. q_out) annotation (Line(
            points={{-50,48},{-50,66},{-34,66}},
            color={107,45,134},
            thickness=1));
        connect(Stream.q_in, Stream_kH.q_in) annotation (Line(
            points={{-50,28},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_out, Stream_kH.q_in) annotation (Line(
            points={{-64,28},{-64,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_in, SRX.q_out) annotation (Line(
            points={{-64,48},{-64,66},{-34,66}},
            color={107,45,134},
            thickness=1));
        connect(SRX_all.y, UnmarkedATP.u[1]) annotation (Line(points={{87.02,46},
                {92,46},{92,23.05},{86,23.05}}, color={0,0,127}));
        connect(DRX_all.y, UnmarkedATP.u[2]) annotation (Line(points={{87.02,-64},
                {92,-64},{92,20.95},{86,20.95}}, color={0,0,127}));
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

      model DRXOnly_optim "Prepared for optimization"
        extends mantATP.SRX_ATP.DRXOnly(
          DrxRate(displayUnit="l/min") = 3.00009E-06,
          RelaxedATPCycling(displayUnit="l/min") = 4.16667E-05,
          relaxingRate(displayUnit="l/min") = 2.66677E-05,
          SRX_all(nu=2));
        extends Data.Toepfer2020;

        Optimization.Criteria.Signals.IntegratedSquaredDeviation
          integratedSquaredDeviation
          annotation (Placement(transformation(extent={{46,18},{26,38}})));
        Modelica.Blocks.Sources.RealExpression realExpression2(y=ATPFluorescence)
          annotation (Placement(transformation(extent={{74,24},{54,44}})));
      equation
        connect(DRX_ATP.q_out, Stream_kH.q_in) annotation (Line(
            points={{-50,-50},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(DRX_ATP.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},
                {-44,-65.05},{74,-65.05}}, color={0,0,127}));
        connect(DRX_ADP.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},
                {36,-62.95},{74,-62.95}}, color={0,0,127}));
        connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
            Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
        connect(unlimitedSolutePumpOut1.q_in, DRX_ATP.q_out) annotation (Line(
            points={{-80,-30},{-50,-30},{-50,-50}},
            color={107,45,134},
            thickness=1));
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
        connect(Stream.q_in, Stream_kH.q_in) annotation (Line(
            points={{-50,28},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_out, Stream_kH.q_in) annotation (Line(
            points={{-64,28},{-64,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(integratedSquaredDeviation.u1, realExpression2.y)
          annotation (Line(points={{48,34},{53,34}}, color={0,0,127}));
        connect(UnmarkedATP.y, integratedSquaredDeviation.u2) annotation (Line(points={{72.98,
                22},{48,22}},                      color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StartTime=-1000,
            StopTime=200,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end DRXOnly_optim;

      model DRXOnly_optimized
        extends mantATP.SRX_ATP.DRXOnly(
          DrxRate(displayUnit="l/min") = 3.00009E-06*tune_a,
          RelaxedATPCycling(displayUnit="l/min") = 4.16667E-05*tune_b,
          relaxingRate(displayUnit="l/min") = 2.66677E-05*tune_c);
        extends Data.Toepfer2020;
        parameter Real tune_a=2.0107438322077007,
                                   tune_b=1.2262857238897689,
                                               tune_c=2.5880342477940874;
        Optimization.Criteria.Signals.IntegratedSquaredDeviation
          integratedSquaredDeviation
          annotation (Placement(transformation(extent={{42,60},{62,80}})));
        Modelica.Blocks.Sources.RealExpression realExpression2(y=ATPFluorescence)
          annotation (Placement(transformation(extent={{2,66},{22,86}})));
      equation
        connect(DRX_ATP.q_out, Stream_kH.q_in) annotation (Line(
            points={{-50,-50},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(DRX_ATP.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},
                {-44,-65.05},{74,-65.05}}, color={0,0,127}));
        connect(DRX_ADP.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},
                {36,-62.95},{74,-62.95}}, color={0,0,127}));
        connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
            Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
        connect(Stream.q_in, Stream_kH.q_in) annotation (Line(
            points={{-50,28},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_out, Stream_kH.q_in) annotation (Line(
            points={{-64,28},{-64,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(SRX_all.y, UnmarkedATP.u[1]) annotation (Line(points={{87.02,46},
                {96,46},{96,23.05},{86,23.05}}, color={0,0,127}));
        connect(DRX_all.y, UnmarkedATP.u[2]) annotation (Line(points={{87.02,-64},
                {96,-64},{96,20.95},{86,20.95}}, color={0,0,127}));
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

      model DRXOnly_optimizedWData
        extends mantATP.SRX_ATP.DRXOnly(
          DRX_ADP(solute_start=0),
          DRX_ATP(solute_start=0),
          DrxRate(displayUnit="l/min") = 3.00009E-06*tune_a,
          RelaxedATPCycling(displayUnit="l/min") = 4.16667E-05*tune_b,
          relaxingRate(displayUnit="l/min") = 2.66677E-05*tune_c,
          SRX(solute_start=0));
        extends Data.Toepfer2020;
        parameter Real tune_a=2.0107438322077007,
                                   tune_b=1.2262857238897689,
                                               tune_c=2.5880342477940874;
        Optimization.Criteria.Signals.IntegratedSquaredDeviation
          integratedSquaredDeviation
          annotation (Placement(transformation(extent={{42,60},{62,80}})));
        Modelica.Blocks.Sources.RealExpression realExpression2(y=if time >
              decay_time then timeTable_ATPChaseControl.y else 1)
          annotation (Placement(transformation(extent={{2,66},{22,86}})));
        replaceable Data.TimeTable_ATPChaseControl1 timeTable_ATPChaseControl
          constrainedby Modelica.Blocks.Sources.TimeTable
          annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
        Bodylight.Chemical.Components.Stream LoadingRate(SolutionFlow(
              displayUnit="m3/s") = 20)
          annotation (Placement(transformation(extent={{0,-18},{-20,2}})));
        Bodylight.Chemical.Components.Substance UD1(solute_start(displayUnit=
                "mol") = 1)
                 annotation (Placement(transformation(extent={{28,-18},{48,2}})));
      equation
        connect(DRX_ATP.q_out, Stream_kH.q_in) annotation (Line(
            points={{-50,-50},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(DRX_ATP.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},
                {-44,-65.05},{74,-65.05}}, color={0,0,127}));
        connect(DRX_ADP.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},
                {36,-62.95},{74,-62.95}}, color={0,0,127}));
        connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
            Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
        connect(Stream.q_in, Stream_kH.q_in) annotation (Line(
            points={{-50,28},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_out, Stream_kH.q_in) annotation (Line(
            points={{-64,28},{-64,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(SRX_all.y, UnmarkedATP.u[1]) annotation (Line(points={{87.02,46},
                {96,46},{96,23.05},{86,23.05}}, color={0,0,127}));
        connect(DRX_all.y, UnmarkedATP.u[2]) annotation (Line(points={{87.02,-64},
                {96,-64},{96,20.95},{86,20.95}}, color={0,0,127}));
        connect(integratedSquaredDeviation.u1, realExpression2.y)
          annotation (Line(points={{40,76},{23,76}}, color={0,0,127}));
        connect(UnmarkedATP.y, integratedSquaredDeviation.u2) annotation (Line(points
              ={{72.98,22},{8,22},{8,64},{40,64}}, color={0,0,127}));
        connect(LoadingRate.q_in, UD1.q_out) annotation (Line(
            points={{0,-8},{38,-8}},
            color={107,45,134},
            thickness=1));
        connect(LoadingRate.q_out, Stream_kH.q_in) annotation (Line(
            points={{-20,-8},{-24,-8},{-24,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          experiment(
            StartTime=-60,
            StopTime=600,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end DRXOnly_optimizedWData;

      package Optimized
        model DRXOnly_optimizedWData_Ctrl1
          extends DRXOnly_optimizedWData(
            tune_c=7.58,
            tune_b=0.4424,
            tune_a=1.29414);
          annotation (experiment(
              StartTime=-60,
              StopTime=600,
              __Dymola_NumberOfIntervals=5000,
              __Dymola_Algorithm="Dassl"));
        end DRXOnly_optimizedWData_Ctrl1;

        model DRXOnly_optimizedWData_Ctrl2
          extends DRXOnly_optimizedWData(
            redeclare Data.TimeTable_ATPChaseControl2 timeTable_ATPChaseControl,
            tune_a=6.49009289961257,
            tune_b=10.0,
            tune_c=10.0);

          annotation (experiment(
              StartTime=-60,
              StopTime=600,
              __Dymola_NumberOfIntervals=5000,
              __Dymola_Algorithm="Dassl"));
        end DRXOnly_optimizedWData_Ctrl2;

        model DRXOnly_optimizedWData_Ctrl3
          extends DRXOnly_optimizedWData(
            redeclare Data.TimeTable_ATPChaseControl3 timeTable_ATPChaseControl,
            tune_a=2.456799254955921,
            tune_b=10.004366201480705,
            tune_c=10.010792516281183);

          annotation (experiment(
              StartTime=-60,
              StopTime=600,
              __Dymola_NumberOfIntervals=5000,
              __Dymola_Algorithm="Dassl"));
        end DRXOnly_optimizedWData_Ctrl3;

        model DRXOnly_optimizedWData_Ctrl4
          extends DRXOnly_optimizedWData(redeclare
              Data.TimeTable_ATPChaseControl4 timeTable_ATPChaseControl);
        end DRXOnly_optimizedWData_Ctrl4;

        model DRXOnly_optimizedWData_Mava1
          extends DRXOnly_optimizedWData(
            tune_c=7.58,
            tune_b=0.4424,
            tune_a=1.29414,
            redeclare Data.TimeTable_NormATPChaseMava1
              timeTable_ATPChaseControl);
          annotation (experiment(
              StartTime=-60,
              StopTime=600,
              __Dymola_NumberOfIntervals=5000,
              __Dymola_Algorithm="Dassl"));
        end DRXOnly_optimizedWData_Mava1;
      end Optimized;
    end SRX_ATP;

    package Predictions
      model DRXOnly_optimized_PredictRunUp
        extends mantATP.SRX_ATP.DRXOnly_optimized(
          Stream_kH(SolutionFlow=hydrolysisRate),
          DRX_ADP(solute_start=0),
          DRX_ATP(solute_start=0),
          SRX(solute_start=0));
        Bodylight.Chemical.Components.Stream LoadingRate(SolutionFlow(
              displayUnit="m3/s") = LoadingRateConst)
          annotation (Placement(transformation(extent={{-10,-4},{-30,16}})));
        Bodylight.Chemical.Components.Substance UD1(solute_start(displayUnit=
                "mol") = 1)
                 annotation (Placement(transformation(extent={{18,-4},{38,16}})));
        parameter Bodylight.Types.VolumeFlowRate LoadingRateConst(displayUnit=
              "l/min") = 0.00033333333333333
          "Volumetric flow of solution if useSolutionFlowInput=false";
        parameter Bodylight.Types.VolumeFlowRate hydrolysisRate(displayUnit=
              "l/min") = 0.0016666666666667
          "Volumetric flow of solution if useSolutionFlowInput=false";
        Bodylight.Chemical.Components.Stream kH1(SolutionFlow(displayUnit=
                "l/min") = hydrolysisRateReverse)
          annotation (Placement(transformation(extent={{0,-56},{-20,-36}})));
        parameter Bodylight.Types.VolumeFlowRate hydrolysisRateReverse(
            displayUnit="l/min") = 0.00016666666666667
          "Volumetric flow of solution if useSolutionFlowInput=false";
      equation
        connect(LoadingRate.q_in, UD1.q_out) annotation (Line(
            points={{-10,6},{28,6}},
            color={107,45,134},
            thickness=1));
        connect(LoadingRate.q_out, Stream_kH.q_in) annotation (Line(
            points={{-30,6},{-50,6},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(kH1.q_out, Stream_kH.q_in) annotation (Line(
            points={{-20,-46},{-50,-46},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(kH1.q_in, DRX_ADP.q_out) annotation (Line(
            points={{0,-46},{30,-46},{30,-50}},
            color={107,45,134},
            thickness=1));
        annotation (experiment(
            StartTime=-600,
            StopTime=600,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end DRXOnly_optimized_PredictRunUp;

      model DRXOnly_optimized_KOHydrolysis
        extends SRX_ATP.DRXOnly_optimized(
          Stream_kH(SolutionFlow=0.00083333333333333),
          tune_a=2.1746486709480752,
          tune_b=1.269709432636263,
          tune_c=1.4297135376023438);
      end DRXOnly_optimized_KOHydrolysis;

      model DRXOnly_optimized_60loading
        extends SRX_ATP.DRXOnly_optimized(
          DRX_ADP(solute_start=0),
          DRX_ATP(solute_start=0),
          SRX(solute_start=0),
          tune_a=2.284922505329536,
          tune_b=1.2623659767193207,
          tune_c=6.734504210787315);
        Bodylight.Chemical.Components.Stream LoadingRate(SolutionFlow(
              displayUnit="l/min") = 0.00033333333333333)
          annotation (Placement(transformation(extent={{-10,-14},{-30,6}})));
        Bodylight.Chemical.Components.Substance UD1(solute_start(displayUnit=
                "mol") = 1)
                 annotation (Placement(transformation(extent={{18,-14},{38,6}})));
      equation
        connect(LoadingRate.q_in, UD1.q_out) annotation (Line(
            points={{-10,-4},{28,-4}},
            color={107,45,134},
            thickness=1));
        connect(LoadingRate.q_out, Stream_kH.q_in) annotation (Line(
            points={{-30,-4},{-34,-4},{-34,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        annotation (experiment(
            StartTime=-60,
            StopTime=600,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end DRXOnly_optimized_60loading;

      model DRXOnly_optimized_600loading
        extends SRX_ATP.DRXOnly_optimized(
          DRX_ADP(solute_start=0),
          DRX_ATP(solute_start=0),
          SRX(solute_start=0));
        Bodylight.Chemical.Components.Stream LoadingRate(SolutionFlow(
              displayUnit="l/min") = 0.00033333333333333)
          annotation (Placement(transformation(extent={{0,-18},{-20,2}})));
        Bodylight.Chemical.Components.Substance UD1(solute_start(displayUnit=
                "mol") = 1)
                 annotation (Placement(transformation(extent={{28,-18},{48,2}})));
      equation
        connect(LoadingRate.q_in,UD1. q_out) annotation (Line(
            points={{0,-8},{38,-8}},
            color={107,45,134},
            thickness=1));
        connect(LoadingRate.q_out, Stream_kH.q_in) annotation (Line(
            points={{-20,-8},{-50,-8},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        annotation (experiment(
            StartTime=-600,
            StopTime=600,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end DRXOnly_optimized_600loading;
    end Predictions;

    package LoadingPhase
      model LoadingTest1
        extends SRX_ATP.Optimized.DRXOnly_optimizedWData_Ctrl1(
          kH(displayUnit="m3/s") = 0.0016666666666667*tune_kH, Stream(
              SolutionFlow(displayUnit="l/min") = relaxingRate),    UD1(solute_start=
                tune_inputSolute),
          Stream_kH(SolutionFlow=kH),
          c=300.0,
          tune_c=758.0,
          tune_a=3.238208447523665,
          tune_b=0.6824655409812524);
        parameter Bodylight.Types.AmountOfSubstance tune_inputSolute=1
          "Initial solute amount in compartment";
            Real srx_rate;
            Real srx_rateErr;
            parameter Real tune_kH=20.0;
            parameter Real srx_rate_demand = 0.1;
      equation
       srx_rate = der(SRX.solute);
       when time > 0 then
         srx_rateErr = max(srx_rate - srx_rate_demand, 0)^2;
       end when;
        annotation (experiment(
            StartTime=-60,
            StopTime=600,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end LoadingTest1;

      model LoadingTest_alternativeLoading
        extends SRX_ATP.Optimized.DRXOnly_optimizedWData_Ctrl1(
          kH(displayUnit="m3/s") = 0.0016666666666667*tune_kH, Stream(
              SolutionFlow(displayUnit="l/min") = relaxingRate),    UD1(solute_start=
                tune_inputSolute),
          Stream_kH(SolutionFlow=kH),
          break connect(LoadingRate.q_out, Stream_kH.q_in));
        parameter Bodylight.Types.AmountOfSubstance tune_inputSolute=1
          "Initial solute amount in compartment";
            Real srx_rate;
            Real srx_rateErr;
            parameter Real tune_kH=1.0;
            parameter Real srx_rate_demand = 0.1;
      equation
       srx_rate = der(SRX.solute);
       when time > 0 then
         srx_rateErr = max(srx_rate - srx_rate_demand, 0)^2;
       end when;
        connect(LoadingRate.q_out, SRX.q_out) annotation (Line(
            points={{-20,-8},{-34,-8},{-34,48},{-34,66}},
            color={107,45,134},
            thickness=1,
            smooth=Smooth.Bezier));
        annotation (experiment(
            StartTime=-1,
            StopTime=600,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end LoadingTest_alternativeLoading;
    end LoadingPhase;
  end mantATP;

  annotation (uses(
      Modelica(version="4.0.0"),
      Optimization(version="2.2.6"),
      Bodylight(version="1.0")));
end CrossBridgeCycling;
