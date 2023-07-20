clear

figure(1); 
figure(2); clf; % hold on; box on;
% figure(3); clf; hold on

ks    = 400;
del_U = 0.015;

N = 24;
P = zeros(N,1);
x0 = [P; 0];

Vlist = [1 10 100 1000 5000]*0.40/100; %  half-sarcomere velocity

opts = odeset('MaxStep',0.1);

rds = [100, 10, 1, 0.1, 0.02];
colors = colormap(lines(length(rds)));
for i = 1:5
  V = Vlist(i);
  Tend_ramp = 0.4/V; % length of ramp
  [t1,x1] = ode15s(@dPdT,2+[0 Tend_ramp],x0,opts,N,V,del_U);
  [t2,x2] = ode15s(@dPdT,2+[Tend_ramp 1200],x1(end,:),opts,N,0,del_U);

  t = [t1; t2];
  x = [x1; x2];
  L = x(:,end);

  % forces
  PF    = 1 - sum(x(:,1:N)')'; % probability that whole chain is folded
  T     = ks*PF.*L.^4 + (L+0.65).^6; % force due to 100% folded fraction PLUS OFFSET

  for j = 1:length(t)
    T(j) = T(j) + ks*sum(x(j,1:N).* max( 0,L(j) - del_U*(1:N) ).^4 );
  end
  T = T + 100*L.^4;

  % viscous force to help get peak in fast ramps
  T(1:length(t1)) = T(1:length(t1)) + 0.1*V;
  
  figure(2);
  semilogx(t,T, '--', 'Color', colors(i, :), 'LineWidth',3);
  hold on;
  rd = rds(i);
    datatable = load(['../data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']).datatable;
    [~, i_peak] = max(datatable(:, 3));
    zone = i_peak:length(datatable(:, 1));
    timebase = datatable(zone, 1) - datatable(zone(1), 1);
    fbase = datatable(zone, 3);
    t_data = datatable(:, 1);
    f_data = datatable(:, 3);
    hold on;
    % curtting off the offset?
    semilogx(t_data, f_data, 'Color', colors(i, :), 'Linewidth', 2);


%   figure(3); plot(log10(t),T)

  PeakModel(i,:) = [0.4/V max(T)];

end

figure(1); plot(t,x(:,1:N),2+[Tend_ramp Tend_ramp],[0 1],'k--')
title('state probabilities versus time for last ramp')

PeakData =[
0.02	14.1969	3.926472727
0.1	    10.6611	3.862672727
1	    7.94194	3.93316
10	    5.9797	3.8093
100	    4.772521951	3.826958537];
PeakModel


figure(6); semilogx(PeakData(:,1),PeakData(:,2),'o',PeakModel(:,1),PeakModel(:,2),'r-','LineWidth',2)
title('Peak tension vs. ramp speed')
figure(7); semilogy(t,T - 100*0.4^4)
title('semilog Y plot tension vs. time')
figure(8); loglog(t,T - 100*0.4^4)
title('log log plot tension vs. time')
figure(2); axis([0 202 0 12])
title('linear plot total tension versus time')
