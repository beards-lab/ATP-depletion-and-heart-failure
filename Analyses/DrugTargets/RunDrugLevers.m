%% RunDrugLevers.m -- dose sweeps of every single-mechanism lever on the human-EHT twitch base.
% For each core (t1_isoT realistic; s3_free possible) and each lever of drugLevers, the twitch is
% simulated at every grid value and summarised: dF (% change of peak developed force vs untreated),
% shape (TTP, 10 %->peak, RT50, RT90, W50 in ms, and their change vs untreated), ATP per head per
% beat, diastolic SRX, peak bound fraction, gate peak. Resumable: results/levers_<core>.mat is
% updated per lever and finished levers are skipped. Output also results/levers_<core>.csv.
here = fileparts(mfilename('fullpath'));  res = fullfile(here, 'results');  if ~exist(res, 'dir'), mkdir(res); end
if ~exist('CORES', 'var'), CORES = {'t1_isoT', 's3_free'}; end
L = drugLevers();
for ic = 1:numel(CORES)
    core = CORES{ic};  out = fullfile(res, ['levers_' core '.mat']);
    if exist(out, 'file'), S = load(out);  T = S.T;  else, T = table(); end
    [base, th, tb, Fd] = drugBase(core);
    [F0, S0] = simHumanTwitch(base, th, tb);  m0 = twitchShape(tb, F0);
    ref = row(core, 'none', 1, F0, S0, tb, m0, m0);
    md = twitchShape(tb, Fd);
    fprintf('%s: untreated peak %.1f kPa (sigma %.2f -> %.0f uN, data %.0f), TTP %.0f/%.0f RT50 %.0f/%.0f RT90 %.0f/%.0f ms (model/data), ATP/beat %.2f\n', ...
        core, m0.Fpk, th.sigma, th.sigma*m0.Fpk, md.Fpk, 1e3*m0.TTP, 1e3*md.TTP, 1e3*m0.RT50, 1e3*md.RT50, 1e3*m0.RT90, 1e3*md.RT90, trapz(tb, S0.R2T));
    if isempty(T) || ~any(strcmp(T.lever, 'none')), T = [T; ref]; end
    for il = 1:numel(L)
        lv = L(il);
        if ~isempty(T) && any(strcmp(T.lever, lv.name)), fprintf('  %-8s done, skipped\n', lv.name); continue; end
        g = lv.grid;  R = cell(numel(g), 1);  t0 = tic;
        parfor ig = 1:numel(g)
            [b, t] = lv.apply(base, th, g(ig));
            try
                [F, S] = simHumanTwitch(b, t, tb);  m = twitchShape(tb, F);
                R{ig} = row(core, lv.name, g(ig), F, S, tb, m, m0);
            catch e
                warning('%s f=%g failed: %s', lv.name, g(ig), e.message);
            end
        end
        T = [T; vertcat(R{:})];  save(out, 'T');
        r = vertcat(R{:});
        fprintf('  %-8s %5.0fs  f = [%s]\n           dF%% = [%s]\n', lv.name, toc(t0), sprintf('%g ', g), sprintf('%+.0f ', r.dF));
    end
    writetable(T, fullfile(res, ['levers_' core '.csv']));
end

function r = row(core, lever, f, F, S, tb, m, m0)
r = table({core}, {lever}, f, m.Fpk, 100*(m.Fpk/m0.Fpk - 1), 1e3*m.TTP, 1e3*m.C10toPk, 1e3*m.RT50, 1e3*m.RT90, 1e3*m.W50, ...
    1e3*(m.TTP - m0.TTP), 1e3*(m.RT50 - m0.RT50), 1e3*(m.RT90 - m0.RT90), trapz(tb, S.R2T), S.SRX(1), max(S.bound), max(S.gate), max(S.Ca), S.Fdia, ...
    'VariableNames', {'core', 'lever', 'f', 'Fpk', 'dF', 'TTP', 'C10toPk', 'RT50', 'RT90', 'W50', 'dTTP', 'dRT50', 'dRT90', 'ATPbeat', 'SRXdia', 'boundPk', 'gatePk', 'CaPk', 'Fdia'});
end
