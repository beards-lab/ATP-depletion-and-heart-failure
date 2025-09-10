function PU = resetSRHelper(PU, SR_in, params)
% ResetSRAt
% SRin

ss = params.ss; % space size (length of the s for each of p1-p3)
dS = params.dS;
p1 = PU(1:ss);
p2 = PU(ss+1:2*ss);
if length(PU) > 3*ss
    p3 = PU(2*ss+1:3*ss);
    Ns = 3; % number of states
else
    p3 = 0;
    Ns = 2; % number of states
end

P_SR = PU(Ns*ss+1);

PD = PU(Ns*ss + 5);

if params.UseSuperRelaxedADP
    P_SRD = PU(Ns*ss+6);
else
    P_SRD = 0;
end

% sum of all probabilities
p1_0 = dS*sum(p1); 
p2_0 = dS*sum(p2); % p2_1 = dS*sum((sign(s+params.dr).*abs(s+params.dr).^params.estiff).*p2);
p3_0 = dS*sum(p3);

% non-hydrolized ATP in non-super relaxed state
PT = max(0, 1 - (p1_0 + p2_0 + p3_0 + PD + P_SR + P_SRD)); % unattached permissive fraction


% reseting the p2
PU(ss+1:2*ss) = p2*max([0, 1/p2_0*(p2_0 - SR_in)]);

% reseting the SR
PU(Ns*ss+1) = min(p2_0+PT, SR_in);
end