function [pts, trc, S] = beardVector(recs, F, Fref, S)
% beardVector  Feature vector of one set of Beard2022-protocol force traces.
%
%   [pts, trc, S] = beardVector(recs, F, Fref)        build the sample plan S too
%   [pts, trc]    = beardVector(recs, F, Fref, S)     reuse it
%
%   recs   1x3 records (8, 4, 2 mM) carrying t (ms) and L (Lo)
%   F      1x3 cell of force traces (kPa) on those records' time grids, model or data
%   Fref   reference force (kPa): this preparation's (or this parameter set's)
%          baseline-subtracted 8 mM isometric force at 1.00 Lo. Empty = take it
%          from F{1} (the model case: one force scale per parameter set).
%
%   pts    the F-V points tier: (F - F0)/Fref at the 1.00 Lo crossing, 8 mM at
%          v = 0.5..6 (v = 0 is the normalisation), 4 and 2 mM at v = 0..6
%   trc    the full-record tier: the baseline-subtracted, 1 ms smoothed trace
%          sampled at protocol-locked times: after each cycle onset (the step and
%          the ramp), after each re-stretch (the redevelopment), and in the first
%          hold. Slack periods are not sampled.
%   S      the plan: S.pt / S.tr are tables (ATP, v, kind, dt) labelling pts / trc

RAMP = [0.5 1 2 3 5 8 12 20 35 60 100 150];             % ms after the cycle onset
REDEV = [2 5 10 20 35 60 100 200 350];                  % ms after the end of the re-stretch
build = nargin < 4 || isempty(S);
if build, S = struct('pt', table(), 'tr', table(), 'idx', {cell(1, 3)}, 'R', {cell(1, 3)}, 'iv', []); end
pts = []; trc = [];
for a = 1:3
    t = recs(a).t; L = recs(a).L; f = double(F{a}(:));
    R = readIsoV(t, L, f); z = median(R.F0(1:7));
    if a == 1 && isempty(Fref), Fref = R.Fx(end) - R.F0(end); end
    [vS, iv] = sort(R.vMeas);                            % 0, 0.5, 1, ..., 6
    p = (R.Fx(iv) - R.F0(iv))' / Fref;
    if a == 1, p = p(2:end); vS = vS(2:end); end        % drop the normalisation point
    pts = [pts, p]; %#ok<AGROW>
    if build
        S.pt = [S.pt; table(repmat(recs(a).atp, numel(vS), 1), round(vS(:) * 2) / 2, 'VariableNames', {'ATP', 'v'})];
        idx = []; lab = table();
        dt = median(diff(t)); Lh = median(L(t < 400)); vNom = round(R.vMeas * 2) / 2;
        add = @(tt) round(tt / dt) + 1;
        idx(end+1) = add(300); lab = [lab; table(recs(a).atp, NaN, "hold110", 300, 'VariableNames', {'ATP', 'v', 'kind', 'dt'})]; %#ok<AGROW>
        for i = 1:7
            for d = RAMP(RAMP < R.tSlack(i) - R.tStart(i) - 0.5)
                idx(end+1) = add(R.tStart(i) + d); %#ok<AGROW>
                lab = [lab; table(recs(a).atp, vNom(i), "ramp", d, 'VariableNames', {'ATP', 'v', 'kind', 'dt'})]; %#ok<AGROW>
            end
            Lt = Lh; tNext = t(end);
            if i < 7, tNext = R.tStart(i + 1); else, Lt = median(L(t > 4100)); end
            jr = find(t > R.tSlack(i) + 2 & L >= Lt - 0.002, 1);
            for d = REDEV(t(jr) + REDEV < tNext - 2)
                idx(end+1) = add(t(jr) + d); %#ok<AGROW>
                lab = [lab; table(recs(a).atp, vNom(i), "redev", d, 'VariableNames', {'ATP', 'v', 'kind', 'dt'})]; %#ok<AGROW>
            end
        end
        S.idx{a} = idx; S.tr = [S.tr; lab];
    end
    y = (movmean(f, 21) - z) / Fref;
    trc = [trc, y(S.idx{a})']; %#ok<AGROW>
end
end
