% CreateProtocolVelocityTable.m
dataDir = '../data/03 27 2026 M/';
fpath = [dataDir '06_Merged_8mM_Active_PNB_Mava.txt'];
M = readmatrix(fpath);

% We only want t >= 0
idx0 = find(M(:,1) >= -0.1, 1);
t = M(idx0:end, 1);
L = M(idx0:end, 2);

F = M(idx0:end, 3);

% Ramer-Douglas-Peucker line simplification to find turning points
epsilon = 0.001; % 0.1% length tolerance
keep = rdp(t, L, epsilon);

t_vt = t(keep);
L_vt = L(keep);

% Remove points that are extremely close in time to prevent ODE stiffness
% If dt < 0.001, merge them
clean_td = [true; diff(t_vt) > 0.0005];
t_vt = t_vt(clean_td);
L_vt = L_vt(clean_td);

% Compute velocities (Lo/s)
v_vt = diff(L_vt) ./ diff(t_vt);
v_vt = [v_vt; 0]; % last segment velocity is 0

% velocitytable format: [Time(s), Vel(ML/s), Vel_um(um/s), L(ML)]
% (um/s isn't used actively but we will populate it as v_vt * 2 or similar)
vel_um = v_vt * 2; % placeholder scaling

velocitytable = [t_vt, v_vt, vel_um, L_vt];

% Generate datatable (downsampled full trace for cost function evaluation)
dsf = 20; % downsample factor
datatable = [downsample(t, dsf), downsample(L, dsf), downsample(F, dsf)];

% Save it
outPath = '../data/protocol_03_27_2026_velocitytable.mat';
save(outPath, 'velocitytable', 'datatable');
fprintf('Saved new velocitytable (%dx4) and datatable (%dx3) to %s\n', ...
        size(velocitytable,1), size(datatable,1), outPath);

figure(50); clf;
plot(t, L, 'Color', [0.7 0.7 0.7]); hold on;
plot(velocitytable(:,1), velocitytable(:,4), 'ro-', 'MarkerSize', 4);
title('Extracted Velocity Table over Original L trace');
ylabel('Length (Lo)'); xlabel('Time (s)');
legend('Raw L', 'Velocity Table Segments');

function keep = rdp(x, y, epsilon)
    % Iterative RDP algorithm
    dmax = 0;
    index = 0;
    end_val = length(x);
    
    % Line from start to end: a*x + b*y + c = 0
    % Line defined by p1=(x(1),y(1)) and p2=(x(end),y(end))
    dx = x(end) - x(1);
    dy = y(end) - y(1);
    
    % calculate distance of all points to the line
    % Distance formula: |dy*x0 - dx*y0 + x(end)*y(1) - y(end)*x(1)| / sqrt(dx^2 + dy^2)
    c = x(end)*y(1) - y(end)*x(1);
    norm_factor = sqrt(dx^2 + dy^2);
    
    if norm_factor == 0
        dists = sqrt((x - x(1)).^2 + (y - y(1)).^2);
    else
        dists = abs(dy.*x - dx.*y + c) / norm_factor;
    end
    
    [dmax, index] = max(dists);
    
    if dmax > epsilon
        % Recursive call
        keep1 = rdp(x(1:index), y(1:index), epsilon);
        keep2 = rdp(x(index:end), y(index:end), epsilon);
        
        % Combine results
        keep = [keep1(1:end-1); keep2 + index - 1];
    else
        keep = [1; end_val];
    end
end
