function [En Ftot t_sim] = evaluatePassiveTheSecond(rd, modnames, modvals)
%% test the midpoint finder

% parameter definition
pasparam.N = 100;
pasparam.L0 = 1.5/2; % half-sarcomere

pasparam.c1 = 5;
pasparam.c2 = 15;
pasparam.g1 = 3;
pasparam.g2 = 3;
pasparam.r_a = 1;
pasparam.r_d = 5;
pasparam.k_Fbr = 0.1; % breaking force coefficient

% update the params with the modifiers, e.g.
% modnames = {"c1", "c2", "g1", "g2", "r_a", "r_d","k_Fbr"};
for i = 1:length(modnames)
    pasparam.(modnames{i}) = pasparam.(modnames{i})*modvals(i);
end


% define force functions
F1 = @(L1)(pasparam.c1.*(L1-pasparam.L0)).^pasparam.g1;
F2 = @(L2)(pasparam.c2.*L2).^pasparam.g2;

balEq = @(L1, L) F1(L1) - F2(L - L1);% balance equation, define F1 = F2 for unattached


% define the space
s_ap = [];    
s = linspace(pasparam.L0, 2.4/2, pasparam.N); % halfsarcomere space
ds = (s(end)-s(1))/pasparam.N; % so that = ds*N+1.6

% test the balance finder
% tic
% for L = s
%     % Use fzero function to solve the equation
%     X = fzero(@(L1)balEq(L1, L), [0 L]);
%     s_ap = [s_ap, X];
% end
% toc
% figure(1);clf;
% plot(s, s, s, s_ap, s, F1(s_ap), '--', s, F2(s - s_ap), ':');
%     % Display the solution
%     disp(['The solution for X is: ', num2str(X)]);
% 
%     legend('Length','L1','Force1', 'force2');
% ylim([0 max(s)]);


% init at lump
L = 1.6/2; % hals sarcomere
att = zeros(1, pasparam.N); % vector of the attached

% Use fzero function to solve the force balance equation to find the midpoint
try
s_att = fzero(@(L1)balEq(L1, L), [pasparam.L0 L]);
catch
    hovno
end

% calc the midpoint index
i_att = round((s_att-pasparam.L0)/ds);

tearingForce = pasparam.k_Fbr*abs(balEq(s', L)); % force at each of our bin from the balance midpoint

[t_sim,x] = ode15s(@dattdt,[0 1000],att,[], i_att, pasparam.r_a, pasparam.r_d, tearingForce);
att = x(end, :);
p_a = sum(att);
 
% clf;hold on;
% plot(s, x(1:10:end, :), '-', 'LineWidth',0.5);
% plot(s, x(end, :), '-', 'LineWidth',1.5);

%% ramp up
% rd = 1; % ramp duration
rs = 1.6; % ramp start um
re = 2.4; % ramp end, um
dL = (re-rs)/2; % length change over the s, half-sarcomere
V = dL/rd; % ramp velocity
dt = ds/V; % numerical time step = ds/V
p_a = [];

t_sim = [];fatt = [];p_u = [];funa = [];
for i_ramp = 1:rd/dt

    L = i_ramp*V*dt + rs/2;
    s_att = fzero(@(L1)balEq(L1, L), [pasparam.L0 L]);
    % calc the index
    i_att = round((s_att - pasparam.L0)/ds);
    tearingForce = pasparam.k_Fbr *abs(balEq(s', L)); % force at each of our bin from the balance midpoint

    [~,x] = ode15s(@dattdt,[0 dt],att,[], i_att, pasparam.r_a, pasparam.r_d, tearingForce);
    att = x(end, :);
    p_a = sum(att);

    t_sim(i_ramp) = i_ramp*dt;
    fatt(i_ramp) = sum(att.*F2(L - s));
    p_u(i_ramp) = 1 - sum(att);
    funa(i_ramp) = p_u(i_ramp)*F1(s_att);
end

% figure(1);clf;
% % semilogy(t, p_u, t, funa, t, fatt, t, fatt+funa);
% plot(t_sim, p_u, t_sim, funa, t_sim, fatt, t_sim, fatt+funa);
% legend('Attached %', 'F una', 'F att', 'F total')

%% hold on after the ramp

% i_ramp, att, L and taeringForce hold from the last iteration of the ramp
[t,x] = ode15s(@dattdt,[dt min(rd*5, 200)],att,[], i_att, pasparam.r_a, pasparam.r_d, tearingForce);
% att = x(end, :);


% clf;hold on;
% plot(s, x(1:N/10:end, :), '-', 'LineWidth',0.5);
% plot(s, x(end, :), '-', 'LineWidth',1.5);

% copy into the ramp-up results
for i_rec = 1:length(t)
    t_sim(i_rec + i_ramp) = t_sim(i_ramp) + t(i_rec);
    fatt(i_rec + i_ramp) = sum(x(i_rec, :).*F2(L - s));
    p_ui = 1 - sum(x(i_rec, :));
    p_u(i_rec + i_ramp) = p_ui;
    funa(i_rec + i_ramp) = p_ui*F1(s_att);
end
%%

t_data = t_sim + 2; % time shifted to match the data times
Ftot = funa + fatt;

datastruct = load(['data/bakers_passiveStretch_' num2str(rd*1000) 'ms.mat']);
datatable = datastruct.datatable;
Ftot_int = interp1(t_data, Ftot, datatable(:,1));    
Es = (Ftot_int - datatable(:,3)).^2; % error set
Es(isnan(Es)) = 0; % zero outside bounds
En = 1e3*sum(Es)/length(Es); % normalized error

figure(round(rd*1000));clf;
% semilogy(t_sim, p_u, t_sim, funa, t_sim, fatt, t_sim, fatt+funa);
% plot(t_data, p_u, '--', t_data, funa, '-', t_data, fatt);
% plot(datatable(:,1), Ftot_int, 'x-', datatable(:,1), datatable(:,3), 'o-' ,'LineWidth',2);
semilogy(t_data, p_u, '--', t_data, funa, '-', t_data, fatt);
hold on;
semilogy(datatable(:,1), Ftot_int, 'x-', datatable(:,1), datatable(:,3), 'o-' ,'LineWidth',2);

legend('U %', 'U_F', 'A_F', 'F total', 'Data');
xlim([t_data(1) t_data(end)]);
