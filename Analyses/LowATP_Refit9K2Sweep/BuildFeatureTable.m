% BuildFeatureTable.m -- one table of everything we measure or can read out,
% 8 mM (data, incumbent) and 2 mM (data, sets A / B / C), from the saved runs.
% No simulation. Output: results/feature_table.json (+ .txt).
%
% Model sources:  8 mM incumbent      levers_ref.mat (fm8L; F-V from the 8 mM battery)
%                 2 mM A (R2 x0.54)   k2sweep.mat RK2(g=0.54)          + bestfit_k2only.mat (F-V, FV2)
%                 2 mM B (+ksr0 x0.4) ksr0grid.mat GR(g .54, f .4)     + bestfit_k2ksr0.mat
%                 2 mM C (+kmsrd x0.4) extB2.mat RL(kmsrd x0.4)        + bestfit_k2ksr0kmsrd.mat
% Data: pooled targets (slack: pooledTarget8mM/2mM; ktr protocol: pooledKtrTarget;
% FV2: fv2Target); F-V curve: Baker 2022 (8 mM scored, 2 mM not scored).
% Colour code (JSON field z): |model - data| in units of the SD the objective uses.
rootT = 'C:\home\git\ATP-depletion-and-heart-failure'; cd(rootT); addpath(genpath(rootT));
rs = fullfile(rootT,'Analyses','LowATP_Refit9K2Sweep','results');
Lr = load(fullfile(rs,'levers_ref.mat'));  fm8 = Lr.fm8L; fd8 = Lr.fd8L;
K2 = load(fullfile(rs,'k2sweep.mat'), 'RK2');  iA = find(abs([K2.RK2.g]-0.54) < 1e-9, 1);
fmA = K2.RK2(iA).fm; fd2 = K2.RK2(iA).fd;
G  = load(fullfile(rs,'ksr0grid.mat'), 'GR'); iB = find(abs([G.GR.g]-0.54) < 1e-9 & abs([G.GR.f]-0.4) < 1e-9, 1); fmB = G.GR(iB).fm;
E2 = load(fullfile(rs,'extB2.mat'), 'RL'); iC = find(abs([E2.RL.f]-0.4) < 1e-9, 1); fmC = E2.RL(iC).fm;
BF = {load(fullfile(rs,'bestfit_k2only.mat')), load(fullfile(rs,'bestfit_k2ksr0.mat')), load(fullfile(rs,'bestfit_k2ksr0kmsrd.mat'))};
PT8 = pooledTarget8mM(); PT2 = pooledTarget2mM(); KT8 = pooledKtrTarget(8); KT2 = pooledKtrTarget(2); FV2 = fv2Target();
LoadData; D = Data_ATP;   % Baker F-V: v, 8 mM, (4 mM), 2 mM

rows = {};
% ---- slack block: {field, label, unit, scale, segments}
SL = { 'steady','Steady force (post-restretch)','kPa',1,1:5
       'A','Redeveloped force (short length)','kPa',1,1:5
       'peak1_y','Re-stretch peak','kPa',1,1:5
       'vall_y','Valley after peak','kPa',1,1:5
       'peak2','Second peak','kPa',1,1:5
       'ktrfit','Slack k_tr (band fit)','1/s',1,1:5
       'ktr80','Slack k_tr80','1/s',1,1:5
       'rsK80','Post-restretch rate rsK80','1/s',1,1:4
       't0','Redevelopment delay t0','ms',1000,1:5
       'rsTfitE','Recovery time rsTfitE','ms',1000,1:4
       'vall_tq','Valley time','ms',1000,1:5
       'peak1_tq','Peak time','ms',1000,1:5 };
for i = 1:size(SL,1)
    nm = SL{i,1}; sg = SL{i,5}; sc = SL{i,4};
    j8 = find(strcmp(PT8.name, nm), 1); j2 = find(strcmp(PT2.name, nm), 1);
    d8 = mean(PT8.val(j8, sg))*sc; d2 = mean(PT2.val(j2, sg))*sc;
    s8 = PT8.s(j8); s2 = PT2.s(j2);
    m = cellfun(@(f) iM(f, nm, sg)*sc, {fm8, fmA, fmB, fmC});
    rows(end+1,:) = iRow('Slack–restretch', SL{i,2}, SL{i,3}, d8, s8, m(1), d2, s2, m(2:4), PT8.r(j8), PT8.sr(j8)); %#ok<SAGROW>
end
% ---- ktr protocol
KP = { 'ktrP_k80','k_tr protocol k80','1/s'; 'ktrP_k63','k_tr protocol k63','1/s'; 'ktrP_F20','Force at 20 ms (F/F_iso)','–' };
for i = 1:size(KP,1)
    j = find(strcmp(KT8.name, KP{i,1}), 1);
    m = cellfun(@(f) iM(f, KP{i,1}, 1), {fm8, fmA, fmB, fmC});
    r = KT2.val(j)/KT8.val(j);
    rows(end+1,:) = iRow('k_tr protocol', KP{i,2}, KP{i,3}, KT8.val(j), KT8.s(j), m(1), KT2.val(j), KT2.s(j), m(2:4), r, r*sqrt(KT8.s(j)^2+KT2.s(j)^2)); %#ok<SAGROW>
end
% ---- force-velocity (Baker curve; model from runFVExperiment)
fvd8 = iFV(D(:,1), D(:,2)/D(1,2)); fvd2 = iFV(D(:,1), D(:,4)/D(1,4));
fvm8 = iFV(fm8.FV_v(:), fm8.FV_fnorm(:));
fvm2 = cellfun(@(b) iFV(b.fv2.FV_v(:), b.fv2.FV_fnorm(:)), BF, 'uni', 0);
FVR = { 'f1','F/F_0 at 1 ML/s','–'; 'f2','F/F_0 at 2 ML/s','–'; 'pmax','Peak power (F/F_0·v)','ML/s'; 'vpk','Velocity at peak power','ML/s' };
for i = 1:size(FVR,1)
    k = FVR{i,1}; m2 = cellfun(@(x) x.(k), fvm2);
    rows(end+1,:) = iRow('Force–velocity (Baker)', FVR{i,2}, FVR{i,3}, fvd8.(k), NaN, fvm8.(k), fvd2.(k), NaN, m2, fvd2.(k)/fvd8.(k), NaN); %#ok<SAGROW>
end
fvp2 = cellfun(@(b) b.fvp2, BF);
rows(end+1,:) = iRow('Force–velocity (Baker)', 'FV2 in-prep P/P_0, 2 ML/s (n=2)', '–', FV2.m8, FV2.s8/FV2.m8, BF{1}.fvp8, FV2.m2, FV2.s2/FV2.m2, fvp2, FV2.r, FV2.sr);
% ---- physiology read-outs (no data; objective bounds shown instead)
PH = { 'XTOR','Isometric ATP turnover XTOR','1/s','3–15'
       'XTOR_vmax','Turnover at v_max','1/s','5–25'
       'attached_ss','Attached fraction (isometric)','–','0.10–0.45'
       'SRX_ss','Super-relaxed fraction','–','0.01–0.40'
       'PT_ss','Detached M·ATP fraction','–','0.01–0.70' };
for i = 1:size(PH,1)
    m = cellfun(@(f) iM(f, PH{i,1}, 1:5), {fm8, fmA, fmB, fmC});
    mx = cellfun(@(f) iMax(f, PH{i,1}), {fm8, fmA, fmB, fmC});
    rw = iRow('Model read-outs (bounds, no data)', PH{i,2}, PH{i,3}, NaN, NaN, m(1), NaN, NaN, m(2:4), NaN, NaN);
    rw{end+1} = PH{i,4}; rw{end+1} = mx; %#ok<AGROW>
    rows(end+1,1:numel(rw)) = rw; %#ok<SAGROW>
end
% ---- write
hdr = {'group','feature','unit','data8','sd8','model8','data2','sd2','A','B','C','ratio','ratio_sd','z','bounds','max'};
J = struct('columns', {hdr}, 'sets', struct('A','R2 x0.54','B','R2 x0.54 + ksr0 x0.4','C','R2 x0.54 + ksr0 x0.4 + kmsrd x0.4'), 'rows', {{}});
fid = fopen(fullfile(rs,'feature_table.txt'),'w');
fprintf(fid, '%-38s %8s %8s | %8s %8s %8s %8s | %6s\n', 'feature', 'data8', 'model8', 'data2', 'A', 'B', 'C', 'r data');
for i = 1:size(rows,1)
    rw = rows(i,:); o = struct();
    for c = 1:numel(hdr); if c <= numel(rw); o.(hdr{c}) = rw{c}; else; o.(hdr{c}) = []; end; end
    J.rows{end+1} = o;
    fprintf(fid, '%-38s %8.4g %8.4g | %8.4g %8.4g %8.4g %8.4g | %6.3f\n', rw{2}, rw{4}, rw{6}, rw{7}, rw{9}, rw{10}, rw{11}, rw{12});
end
fclose(fid);
fw = fopen(fullfile(rs,'feature_table.json'),'w'); fwrite(fw, jsonencode(J, 'PrettyPrint', true)); fclose(fw);
type(fullfile(rs,'feature_table.txt'));

function r = iRow(grp, lab, unit, d8, s8, m8, d2, s2, m2, rat, rsd)
    z = [abs(m8-d8)/(abs(d8)*s8), abs(m2-d2)./(abs(d2)*s2)];   % in units of the objective's SD
    r = {grp, lab, unit, d8, s8, m8, d2, s2, m2(1), m2(2), m2(3), rat, rsd, z};
end
function v = iM(st, nm, sg)
    v = NaN;
    if isfield(st, nm); x = double(st.(nm)); x = x(:)'; if numel(x) >= max(sg); x = x(sg); end
        x = x(isfinite(x)); if ~isempty(x); v = mean(x); end; end
end
function v = iMax(st, nm)
    v = NaN; if isfield(st, nm); x = double(st.(nm)); x = x(isfinite(x)); if ~isempty(x); v = max(x); end; end
end
function o = iFV(v, fn)
    v = v(:); fn = fn(:); vv = linspace(0, max(v), 400)'; ff = interp1(v, fn, vv, 'pchip');
    p = vv.*ff; [o.pmax, k] = max(p); o.vpk = vv(k);
    o.f1 = interp1(v, fn, 1, 'pchip'); o.f2 = interp1(v, fn, 2, 'pchip');
end
