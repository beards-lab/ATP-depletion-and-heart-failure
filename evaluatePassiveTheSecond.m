%% test the midpoint finder
s_ap = [];    
N = 100;
s = linspace(1.6, 2.2, N);
L0 = 1.5;
ds = s(end)/N-s(1)/N; % so that = ds*N+1.6

c1 = 5;c2 = 15;g1=3;g2 = 3;

F1 = @(L1)(c1.*(L1-L0)).^g1;
F2 = @(L2)(c2.*L2).^g2;
balEq = @(L1, L) F1(L1) - F2(L - L1);% balance equation, define F1 = F2 for unattached
tic
for L = s

    % Use fzero function to solve the equation
    X = fzero(@(L1)balEq(L1, L), [0 L]);
    s_ap = [s_ap, X];
end
toc
figure(1);clf;
plot(s, s, s, s_ap, s, F1(s_ap), '--', s, F2(s - s_ap), ':');
    % Display the solution
    disp(['The solution for X is: ', num2str(X)]);

    legend('Length','L1','Force1', 'force2');
ylim([0 max(s)]);


%% init rand

att = rand(1, N);
att = att/sum(att)*0.3;

p_a = sum(att);
p_u = 1 - p_a;
r_a = 0.5;
r_d = 0.5;


f_a = att.*F2(max(s) - s);
% figure(1);clf;semilogy(s, att, 's-', s, f_a, 'v-');

dt = 1;

% init at lump
dt = 1;
L = 1.6;
i_att = 2;
% i_att = % index of attachable midpoint
tearingForce = 0.1*abs(balEq(s', L)); % force at each of our bin from the balance midpoint

[t,x] = ode15s(@dattdt,[0 dt*10],att,[], i_att, r_a, r_d, tearingForce);
att = x(end, :);
p_a = sum(att)

clf;hold on;
plot(s, x(1:10:end, :), '-', 'LineWidth',0.5);
plot(s, x(end, :), '-', 'LineWidth',1.5);

%% ramp down from max instant stretch
L = 2.2;

% Use fzero function to solve the equation
s_att = fzero(@(L1)balEq(L1, L), [0 L]);
% calc the index
i_att = round((s_att - L0)/ds);
tearingForce = 0.1*abs(balEq(s', L)); % force at each of our bin from the balance midpoint

[t,x] = ode15s(@dattdt,[0 dt],att,[], i_att, r_a, r_d, tearingForce);

clf;hold on;
plot(s, x(1:N/10:end, :), '-', 'LineWidth',0.5);
plot(s, x(end, :), '-', 'LineWidth',1.5);
%%
%
fatt = zeros(length(t), 1);
funa = zeros(length(t), 1);
p_u = zeros(length(t), 1);
for i = 1:length(t)
    fatt(i) = sum(x(i, :).*F2(L - s));
    p_u(i) = 1 - sum(x(i, :));
    funa(i) = p_u(i)*F1(s_att);
end

figure(1);clf;
% semilogy(t, p_u, t, funa, t, fatt, t, fatt+funa);
plot(t, p_u, t, funa, t, fatt, t, fatt+funa);
legend('Attached %', 'F una', 'F att', 'F total')
