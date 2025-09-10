%% Fix and Fudge Velocity table

datastruct = load('data/bakers_slack8mM_all.mat');
velocitytable =    datastruct.velocitytable;

%% Fix velcity table to have equal speeds
segments = find(velocitytable(:, 2) < 0)';
dl = velocitytable(segments, 4) - velocitytable(segments+1, 4);
% dl = 0.2802;
velocities = 160:40:320;
velocitytable(segments, 3) = -velocities;
velocitytable(segments+1, 1) = velocitytable(segments, 1) + dl./velocities';

dl_r = velocitytable(segments+2, 4) - velocitytable(segments+3, 4);
% dl_r = dl;
velocitytable(segments+3, 1) = velocitytable(segments+2, 1) + abs(dl_r./velocitytable(segments+2, 3));
velocitytable(segments, 2)  = velocitytable(segments, 3)/2;

% save('../data/bakers_slack8mM_all_fix.mat', 'velocitytable', 'datatable');
% plot that and compare
velocitytable(2:end, 5) = (velocitytable(1:end-1, 3).*(velocitytable(2:end, 1) - velocitytable(1:end-1, 1)))
velocitytable(:, 6) = cumsum([velocitytable(:, 5)])' + 2.2;
plot(velocitytable(:, 1), velocitytable(:, 4), velocitytable(:, 1), velocitytable(:, 6));

%% fudge velcity table to have the same length
segments = find(velocitytable(:, 2) < 0)';
% dl = velocitytable(segments, 4) - velocitytable(segments+1, 4);
dl = 0.2802;
velocities = 160:40:320;
velocitytable(segments, 3) = -velocities;
velocitytable(segments+1, 1) = velocitytable(segments, 1) + dl./velocities';

% dl_r = velocitytable(segments+2, 4) - velocitytable(segments+3, 4);
dl_r = dl;
velocitytable(segments+3, 1) = velocitytable(segments+2, 1) + abs(dl_r./velocitytable(segments+2, 3));

velocitytable(segments, 2)  = velocitytable(segments, 3)/2;
% datastruct.velocitytable = velocitytable;

% save('../data/bakers_slack8mM_all_fudge.mat', 'velocitytable', 'datatable');
% plot that
velocitytable(2:end, 5) = (velocitytable(1:end-1, 3).*(velocitytable(2:end, 1) - velocitytable(1:end-1, 1)))
velocitytable(:, 6) = cumsum([velocitytable(:, 5)])' + 2.2;
plot(velocitytable(:, 1), velocitytable(:, 4), velocitytable(:, 1), velocitytable(:, 6));
% velocitytable(:, 4) = velocitytable(:, 6)
%% extend the slackened time to have a better ktr fit
datastruct = load('../data/bakers_slack8mM_all_fix.mat');
velocitytable =    datastruct.velocitytable;
velocitytable(1, 1) = -40;

dur = 0.12;
t_ext = zeros(size(velocitytable, 1), 1);
t_ext(segments+2) = dur - (velocitytable(segments+2, 1)  - velocitytable(segments+1, 1));
velocitytable(:, 7) = t_ext
velocitytable(:, 1) = velocitytable(:, 1) + cumsum(t_ext);

% save('../data/bakers_slack8mM_all_extendedSlack.mat', 'velocitytable', 'datatable');

%% plot that
velocitytable(2:end, 5) = (velocitytable(1:end-1, 3).*(velocitytable(2:end, 1) - velocitytable(1:end-1, 1)))
velocitytable(:, 6) = cumsum([velocitytable(:, 5)])' + 2.2;
plot(velocitytable(:, 1), velocitytable(:, 4), velocitytable(:, 1), velocitytable(:, 6));