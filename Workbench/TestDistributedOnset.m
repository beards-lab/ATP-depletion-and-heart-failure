%% Test distributed passive onset model
clf; hold on;
a = 1;
% V1
y = @(t0, t) max(0, a*(t-t0)) + (t-t0 > 0)*1;

% V2
% y = @(t0, t) max(1, a*(t)).*(t > t0);

t = -2:0.1:10;
ytot = zeros(size(t));
sm = -1:0.1:1;
for t0 = sm
    ytot = ytot + y(t0, t);
    plot(t, y(t0, t), '--', LineWidth=0.5);
end
ytot = ytot / length(sm);

plot(t, y(0, t),'-',  LineWidth=2);hold on;
plot(t, ytot, '-', LineWidth=2);