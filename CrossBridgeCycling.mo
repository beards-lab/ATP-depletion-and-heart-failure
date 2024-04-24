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
end CrossBridgeCycling;
