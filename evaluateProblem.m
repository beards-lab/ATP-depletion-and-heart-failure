function Etot = evaluateProblem(g)

LoadData;


ML = 1.1; % half sarcomere length (microns)
% Non-zero velocities
vel = (-Data_ATP(:,1)).*ML; % micron per sec
MgATP = [8 4 2];
MgADP = 0; 
Pi    = 0; 


%% force x velocity
for k = [1 3]
    for j = 1:length(vel)
        F_active(j,k) = evaluateModel(vel(j), 1, MgATP(k),Pi,MgADP,g);
    end
end



E(1) = sum(abs(F_active(:,1)-Data_ATP(:,2)).^2) + ...
     sum(abs(F_active(:,3)-Data_ATP(:,4)).^2);

%% force x iso MgATP
for k = 1:length(MgATP_iso)
  % Zero velocity:
  F_iso(k) = evaluateModel(0, 1, MgATP_iso(k),Pi,MgADP,g);
end


E(2) = sum(abs(F_iso-F_data).^2);
%% time constant of MgATP

Tspan = [0:0.001:0.12];

for k = 1:length(MgATP)

  % Tspan array returns F active array
  F_active_ktr = evaluateModel(0, Tspan, MgATP(k),Pi,MgADP,g);
  
  Frel = F_active_ktr./F_active_ktr(end);
  
  % get the time constant
  Ktr(k) = 1/interp1(Frel,Tspan,1-exp(-1)); % time constant for Frel(1/Ktr) = 1-exp(-1)
end

E(3) = sum(abs(Ktr-Ktr_mean).^2);
%% Return
Etot = sum(E);

%% plot?