dataDir = '../data/03 27 2026 M/';
fpath = [dataDir '02_Merged_8mM_Active.txt'];
M = readmatrix(fpath);
t = M(:,1); 
L = M(:,2);

% Reduce smoothing to preserve fast micro-events
L_sm = movmean(L, 3);
v = diff(L_sm)./diff(t);

is_moving = abs(v) > 0.005;
% find rising and falling edges of is_moving
starts = find(diff(is_moving) == 1) + 1;
ends = find(diff(is_moving) == -1);

fprintf('Detected movements:\n');
for i=1:length(starts)
    if i > length(ends); break; end
    dt = t(ends(i)) - t(starts(i));
    if dt < 0.005; continue; end % skip noise shorter than 5ms

    
    v_mean = mean(v(starts(i):ends(i)));
    dL = L_sm(ends(i)) - L_sm(starts(i));
    fprintf('Move: t = [%7.3f, %7.3f], dt = %6.3f, dL = %7.3f, v = %7.3f\n', ...
        t(starts(i)), t(ends(i)), dt, dL, v_mean);
end
