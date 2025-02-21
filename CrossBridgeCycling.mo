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
      model DRXOnly
        Bodylight.Chemical.Components.Substance UT(solute_start=0.5 - SRX_init/2)
          annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
        Bodylight.Chemical.Components.Substance UD(solute_start=0.5 - SRX_init/2)
                 annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
        Bodylight.Chemical.Components.Stream kH(SolutionFlow(displayUnit=
                "l/min") = 0.0016666666666667)
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
      equation
        connect(UT.q_out,kH. q_in) annotation (Line(
            points={{-50,-50},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(kH.q_out,UD. q_out) annotation (Line(
            points={{0,-30},{30,-30},{30,-50}},
            color={107,45,134},
            thickness=1));
        connect(UT.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},{-44,
                -65.05},{74,-65.05}}, color={0,0,127}));
        connect(UD.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},{36,
                -62.95},{74,-62.95}}, color={0,0,127}));
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
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
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
          relaxingRate(displayUnit="l/min") = 2.66677E-05);
        extends Data.Toepfer2020;

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
        connect(UT.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},{-44,
                -65.05},{74,-65.05}}, color={0,0,127}));
        connect(UD.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},{36,
                -62.95},{74,-62.95}}, color={0,0,127}));
        connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
            Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
        connect(unlimitedSolutePumpOut1.q_in,UT. q_out) annotation (Line(
            points={{-80,-30},{-50,-30},{-50,-50}},
            color={107,45,134},
            thickness=1));
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
        connect(Stream.q_in, kH.q_in) annotation (Line(
            points={{-50,28},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_out, kH.q_in) annotation (Line(
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
        connect(UT.q_out,kH. q_in) annotation (Line(
            points={{-50,-50},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(UT.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},{-44,
                -65.05},{74,-65.05}}, color={0,0,127}));
        connect(UD.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},{36,
                -62.95},{74,-62.95}}, color={0,0,127}));
        connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
            Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
        connect(Stream.q_in, kH.q_in) annotation (Line(
            points={{-50,28},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_out, kH.q_in) annotation (Line(
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
          DrxRate(displayUnit="l/min") = 3.00009E-06*tune_a,
          RelaxedATPCycling(displayUnit="l/min") = 4.16667E-05*tune_b,
          relaxingRate(displayUnit="l/min") = 2.66677E-05*tune_c,
          SRX(solute_start=0),
          UT(solute_start=0),
          UD(solute_start=0));
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
        connect(UT.q_out,kH. q_in) annotation (Line(
            points={{-50,-50},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(UT.solute, DRX_all.u[1]) annotation (Line(points={{-44,-60},{-44,
                -65.05},{74,-65.05}}, color={0,0,127}));
        connect(UD.solute, DRX_all.u[2]) annotation (Line(points={{36,-60},{36,
                -62.95},{74,-62.95}}, color={0,0,127}));
        connect(realExpression1.y,unlimitedSolutePumpOut1. soluteFlow) annotation (
            Line(points={{-81,-90},{-94,-90},{-94,-34}},  color={0,0,127}));
        connect(SRX.solute, SRX_all.u[1]) annotation (Line(points={{-28,56},{-28,
                46},{74,46}}, color={0,0,127}));
        connect(Stream.q_in, kH.q_in) annotation (Line(
            points={{-50,28},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(Stream1.q_out, kH.q_in) annotation (Line(
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
        connect(LoadingRate.q_out, kH.q_in) annotation (Line(
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
      end Optimized;
    end SRX_ATP;

    package Predictions
      model DRXOnly_optimized_PredictRunUp
        extends mantATP.SRX_ATP.DRXOnly_optimized(
          SRX(solute_start=0),
          UT(solute_start=0),
          UD(solute_start=0),
          kH(SolutionFlow=hydrolysisRate));
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
        connect(LoadingRate.q_out, kH.q_in) annotation (Line(
            points={{-30,6},{-50,6},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(kH1.q_out, kH.q_in) annotation (Line(
            points={{-20,-46},{-50,-46},{-50,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        connect(kH1.q_in, UD.q_out) annotation (Line(
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
          kH(SolutionFlow=0.00083333333333333),
          tune_a=2.1746486709480752,
          tune_b=1.269709432636263,
          tune_c=1.4297135376023438);
      end DRXOnly_optimized_KOHydrolysis;

      record OptimFun
      extends Optimization.Internal.Version.Current.ModelOptimizationSetup(
        modelName="CrossBridgeCycling.mantATP.Predictions.DRXOnly_optimized_KOHydrolysis",
        plotScript="",
        saveSetup=true,
        saveSetupFilename="OptimizationLastRunModel.mo",
        convertSetup=false,
        askForTunerReUse=true,
        tuner=Optimization.Internal.Version.Current.Tuner(
            UseTunerMatrixForDiscreteValues=false,
            tunerParameters={Optimization.Internal.Version.Current.TunerParameter(
              name="tune_a",
              active=true,
              Value=2.01074,
              scaleToBounds=false,
              min=-10.0,
              max=10.0,
              equidistant=0,
              discreteValues=fill(0.0, 0),
              unit=""),Optimization.Internal.Version.Current.TunerParameter(
              name="tune_b",
              Value=1.22629,
              discreteValues=fill(0, 0),
              unit=""),Optimization.Internal.Version.Current.TunerParameter(
              name="tune_c",
              Value=2.58803,
              discreteValues=fill(0, 0),
              unit="")},
            tunerMatrix=fill(0.0, 0, 3)),
        criteria={Optimization.Internal.Version.Current.Criterion(
            name="integratedSquaredDeviation.y1",
            active=true,
            usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
            demand=1.0,
            unit="")},
        preferences=Optimization.Internal.Version.Current.Preferences(
            optimizationOptions=
              Optimization.Internal.Version.Current.OptimizationOptions(
              method=Optimization.Internal.Version.Current.Types.OptimizationMethod.sqp,
              ObjectiveFunctionType=Optimization.Internal.Version.Current.Types.ObjectiveFunctionType.Max,
              OptTol=0.001,
              maxEval=1000,
              GridBlock=50,
              evalBestFinal=false,
              saveBest=true,
              saveHistory=true,
              listFilename="OptimizationLog.log",
              listOn=true,
              listOnline=true,
              listIncrement=100,
              numberOfShownDigits=3,
              onPlace=true,
              listTuners=true,
              GaPopSize=10,
              GaNGen=100),
            simulationOptions=
              Optimization.Internal.Version.Current.SimulationOptions(
              startTime=0.0,
              stopTime=1.0,
              outputInterval=0.0,
              numberOfIntervals=500,
              integrationMethod=Optimization.Internal.Version.Current.Types.IntegrationMethod.Dassl,
              integrationTolerance=0.001,
              fixedStepSize=0.0,
              autoLoadResults=true,
              useDsFinal=true,
              translateModel=false,
              setCriteriaSimulationFailed=true,
              CriteriaSimulationFailedValue=1000000.0,
              simulationMode=Optimization.Internal.Version.Current.Types.SimulationMode.Single,
              parallelizationMode=Optimization.Internal.Version.Current.Types.ParallelizationMode.None,
              numberOfThreads=0,
              copyFiles=fill("", 0)),
            sensitivityOptions=
              Optimization.Internal.Version.Current.SensitivityOptions(
              TypeOfSensitivityComputation=Optimization.Internal.Version.Current.Types.SensitivityMethod.ExternalDifferencesSymmetric,
              automaticSensitivityTolerance=true,
              sensitivityTolerance=1E-06)));

      end OptimFun;

      model DRXOnly_optimized_60loading
        extends SRX_ATP.DRXOnly_optimized(
          UT(solute_start=0),
          UD(solute_start=0),
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
        connect(LoadingRate.q_out, kH.q_in) annotation (Line(
            points={{-30,-4},{-34,-4},{-34,-30},{-20,-30}},
            color={107,45,134},
            thickness=1));
        annotation (experiment(
            StartTime=-60,
            StopTime=600,
            __Dymola_NumberOfIntervals=5000,
            __Dymola_Algorithm="Dassl"));
      end DRXOnly_optimized_60loading;

      record simsetup
        extends Optimization.Internal.Version.Current.ModelOptimizationSetup(
          modelName=
              "CrossBridgeCycling.mantATP.Predictions.DRXOnly_optimized_60loading",

          plotScript="",
          saveSetup=true,
          saveSetupFilename="OptimizationLastRunModel.mo",
          convertSetup=false,
          askForTunerReUse=true,
          tuner=Optimization.Internal.Version.Current.Tuner(
                    UseTunerMatrixForDiscreteValues=false,
                    tunerParameters={
                Optimization.Internal.Version.Current.TunerParameter(
                      name="tune_a",
                      active=true,
                      Value=2.01074,
                      scaleToBounds=false,
                      min=-10.0,
                      max=10.0,
                      equidistant=0,
                      discreteValues=fill(0.0, 0),
                      unit=""),
                Optimization.Internal.Version.Current.TunerParameter(
                      name="tune_b",
                      Value=1.22629,
                      discreteValues=fill(0, 0),
                      unit=""),
                Optimization.Internal.Version.Current.TunerParameter(
                      name="tune_c",
                      Value=2.58803,
                      discreteValues=fill(0, 0),
                      unit="")},
                    tunerMatrix=fill(0.0, 0, 1)),
          criteria={Optimization.Internal.Version.Current.Criterion(
                    name="integratedSquaredDeviation.y1",
                    active=true,
                    usage=Optimization.Internal.Version.Current.Types.CriterionUsage.Minimize,
                    demand=1.0,
                    unit="")},
          preferences=Optimization.Internal.Version.Current.Preferences(
                    optimizationOptions=
                Optimization.Internal.Version.Current.OptimizationOptions(
                      method=Optimization.Internal.Version.Current.Types.OptimizationMethod.sqp,
                      ObjectiveFunctionType=Optimization.Internal.Version.Current.Types.ObjectiveFunctionType.Max,
                      OptTol=0.001,
                      maxEval=1000,
                      GridBlock=50,
                      evalBestFinal=false,
                      saveBest=true,
                      saveHistory=true,
                      listFilename="OptimizationLog.log",
                      listOn=true,
                      listOnline=true,
                      listIncrement=100,
                      numberOfShownDigits=3,
                      onPlace=true,
                      listTuners=true,
                      GaPopSize=10,
                      GaNGen=100),
                    simulationOptions=
                Optimization.Internal.Version.Current.SimulationOptions(
                      startTime=0.0,
                      stopTime=1.0,
                      outputInterval=0.0,
                      numberOfIntervals=500,
                      integrationMethod=Optimization.Internal.Version.Current.Types.IntegrationMethod.Dassl,
                      integrationTolerance=0.001,
                      fixedStepSize=0.0,
                      autoLoadResults=true,
                      useDsFinal=true,
                      translateModel=false,
                      setCriteriaSimulationFailed=true,
                      CriteriaSimulationFailedValue=1000000.0,
                      simulationMode=Optimization.Internal.Version.Current.Types.SimulationMode.Single,
                      parallelizationMode=Optimization.Internal.Version.Current.Types.ParallelizationMode.None,
                      numberOfThreads=0,
                      copyFiles=fill("", 0)),
                    sensitivityOptions=
                Optimization.Internal.Version.Current.SensitivityOptions(
                      TypeOfSensitivityComputation=Optimization.Internal.Version.Current.Types.SensitivityMethod.ExternalDifferencesSymmetric,
                      automaticSensitivityTolerance=true,
                      sensitivityTolerance=1E-06)));
      end simsetup;

      record fluf
        extends Optimization.Internal.Version.Current.ReUseOfTuners(answer=true,
            saveSetup=true);
      end fluf;
    end Predictions;
  end mantATP;

  annotation (uses(Modelica(version="4.0.0"), Physiolibrary(version="2.4.1"),
      Optimization(version="2.2.6"),
      Bodylight(version="2.4.1")));
end CrossBridgeCycling;
