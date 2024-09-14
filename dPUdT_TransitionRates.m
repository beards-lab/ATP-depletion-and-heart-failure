strainDep = @(alpha, dr) min(1e4, exp((alpha*(s+dr)).^2));
RTD = g2*params.kah*PT;
RD1 = params.ka*PD*N_overlap; % to loosely attachemnt state

R1D = params.kd*p1.*strainDep(params.alpha0, params.dr0); %(exp(-params.alpha1*s)) + params.TK*(s>params.TK0).*s.*p1; % p1 to PU - detachment rate + tearing constant
% R12 = params.k1*p1.*exp(-params.alpha1*s); % P1 to P2
% R21 = f1*params.k_1*p2.*strainDep(params.alpha_1, params.dr_1); % backward flow from p2 to p1

% DAN's XB rates
% R1D = params.kd*p1.*exp(+(params.alpha0*s).^2); %
R12 = params.k1*p1.*exp(-params.alpha1*s); % P1 to P2
R21 = f1*params.k_1*p2.*strainDep(params.alpha_1, params.dr_1); % p2 to p1
% R2T = g2*params.k2*p2.*min(1e9, max(1, strainDep(params.alpha2, params.dr2)));

% DAN's very complicated detachment rate
% lambdaR = 0.015;
% lambdaL = 0.038;
% alpha2_R = 1/0.015;

% params.alpha2_L = 1/0.038;
% params.k2_R = 8e3;
% params.k2_L = 200;
% R0 = 0.10;
% R1 = ((s+0)<=0).*(1 - exp(-(s+0)./lambdaL)).^2;
% R2 = ((s+0)>0).*(1 - exp(+(s+0)./lambdaR)).^2;
% R2T = g2*params.k2*p2.*(R0 + R1 + R2);

% lambdaL = 0.015;
% params.k2 = 1.25*20;

kL = min(1e4, params.k2_L*((s+0)<=0).*(1 - exp(-(s+0)*params.alpha2_L)).^2);

% experiMENTAL
% kL = min(1e4, params.k2_L*((s-params.dr2_L)<=0).*(exp(-(s-params.dr2_L)*params.alpha2_L)));
% kR = 200*(s>0).*(1 - exp(+(s+0)./lambdaR)).^2;
% kR = 0*(s>0).*(s./0.01);
% kR = 600*(s>0).*((s.^2)./((s.^2) + 0.01^2));
% r = 0.030;
% kR = (100e6)*( (1/6).*s.^3 + (-1*r/2).*s.^2 + 15*(r^3).*s).*(s>0.002);

% kR = max(0, params.k2_R*(s-params.dr2_R)); %.*(s>0.002);
kR = max(0, params.k2_R*(s-params.dr2_R)).^params.alpha2_R; %.*(s>0.002);
R2T = p2.*(params.k2 + kL + kR);