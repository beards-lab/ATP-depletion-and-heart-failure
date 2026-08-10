% RunRecoveryAnatomy.m — WHAT carries the force recovery in each window?
%
% The lever scan found that the series-viscoelastic class (eta_M, kSE_M,
% mu_neg) moves rsK while leaving ktr and isometric force untouched. That is
% suspicious: it suggests the model's post-restretch force recovery is partly
% a MECHANICAL re-equilibration of the series element rather than a KINETIC
% re-population of cross-bridges. If so, the model reproduces the right
% recovery AMPLITUDE by the wrong mechanism, which is exactly why its RATE is
% wrong and why no cross-bridge rate constant fixes it.
%
% This decomposes each recovery window into its physical carriers:
%   Force        - what the estimator sees
%   FXB active   - cross-bridge force  (population x strain x stiffness)
%   passive      - titin + Maxwell
%   LSE          - series-element extension (the mechanical store)
%   p1+p2        - bound fraction        (the kinetic state)
%   PD, SR+SRD   - primed reservoir and the SRX trap
%
% and reports, for each, the fraction of its total excursion completed by the
% time FORCE has completed 63 % of its own. A carrier that is already done
% when force is only 63 % recovered is NOT rate-limiting; the one that tracks
% force is.

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));

if ~exist('SNAP','var'); SNAP = 'params/rskR2_w025_opt.m'; end
if ~exist(fullfile(here,'results'),'dir'); mkdir(fullfile(here,'results')); end

p = getParams(loadParams(SNAP), [], true, false);
p.MaxRunTime = 600;
fcn = str2func(p.modelFcn);

ds = load(fullfile('data', p.velocitytableonfile));
vt = ds.velocitytable;

% serial run so out.* is one continuous trajectory with all diagnostics
pr = p; pr.Velocity = vt(:,2);
fprintf('running slack protocol (serial) at %s ...\n', SNAP); tic;
[~, out] = evaluateModel(fcn, vt(:,1), pr);
fprintf('  %.1f s\n', toc);

nr = makeMonotonous(out.t);
T  = out.t(nr);
get = @(f) reshape(out.(f)(nr), [], 1);

F    = get('Force');
LSE  = get('LSE');
p1   = get('p1_0');  p2 = get('p2_0');
PD   = get('PuATP');
SR   = get('SR');
SRD  = zeros(size(SR)); if isfield(out,'SRD'); SRD = get('SRD'); end
FPAS = get('FXBPassive');
FACT = F - FPAS;

iS = find(vt(:,2) < -1);  iR = find(vt(:,2) > 1);

CAR = {'Force','FACT','FPAS','LSE','bound','PD','SRX'};
fprintf(['\nFraction of each carrier''s own excursion completed at the moment\n' ...
         'FORCE reaches 63%% of its recovery (1.00 = already finished, i.e. NOT\n' ...
         'rate-limiting; ~0.63 = tracks force, i.e. IS the rate-limiting carrier)\n\n']);

for kind = 1:2
    if kind == 1; nm = 'postRestretch'; else; nm = 'postSlack'; end
    fprintf('--- %s ---\n%-8s', nm, 'cycle');
    fprintf('%10s', CAR{:}); fprintf('\n');

    for c = 1:numel(iR)
        if kind == 1
            ta = vt(iR(c)+1,1);
            if c < numel(iS); tb = vt(iS(c+1),1); else; tb = max(T); end
            m  = T >= ta & T <= tb;
            tw = T(m); Fw = F(m);
            mv = tw <= tw(1)+0.080; [~,iv] = min(Fw(mv));   % anchor at vall2
            idx = find(m); idx = idx(iv:end);
        else
            ta = vt(iS(c)+1,1); tb = vt(iR(c),1);
            idx = find(T >= ta & T <= tb);
        end
        if numel(idx) < 20; continue; end

        sig = struct('Force',F(idx),'FACT',FACT(idx),'FPAS',FPAS(idx), ...
                     'LSE',LSE(idx),'bound',p1(idx)+p2(idx),'PD',PD(idx), ...
                     'SRX',SR(idx)+SRD(idx));
        tw = T(idx);

        % time at which FORCE completes 63 % of its excursion
        f0 = sig.Force(1); f1 = median(sig.Force(end-round(0.15*numel(idx)):end));
        target = f0 + (1-exp(-1))*(f1-f0);
        if f1 > f0; k63 = find(sig.Force >= target, 1); else; k63 = find(sig.Force <= target, 1); end
        if isempty(k63); continue; end

        fprintf('%-8d', c);
        for ic = 1:numel(CAR)
            s  = sig.(CAR{ic});
            s0 = s(1); s1 = median(s(end-round(0.15*numel(idx)):end));
            if abs(s1-s0) < 1e-9
                fprintf('%10s', '  ~flat');
            else
                fprintf('%10.2f', (s(k63)-s0)/(s1-s0));
            end
        end
        fprintf('   [t63 %5.1f ms]\n', 1e3*(tw(k63)-tw(1)));
    end
    fprintf('\n');
end

%% ---- how much of the recovered force is mechanical vs kinetic? ----------
fprintf(['Composition of the recovered force (end - start of window):\n' ...
         '%-16s %10s %10s %10s %10s\n'], 'window', 'dForce', 'dFactive', 'dFpassive', 'dBound%');
for kind = 1:2
    if kind == 1; nm = 'postRestretch'; else; nm = 'postSlack'; end
    for c = 1:numel(iR)
        if kind == 1
            ta = vt(iR(c)+1,1);
            if c < numel(iS); tb = vt(iS(c+1),1); else; tb = max(T); end
            m = T >= ta & T <= tb; tw = T(m); Fw = F(m);
            mv = tw <= tw(1)+0.080; [~,iv] = min(Fw(mv));
            idx = find(m); idx = idx(iv:end);
        else
            ta = vt(iS(c)+1,1); tb = vt(iR(c),1);
            idx = find(T >= ta & T <= tb);
        end
        if numel(idx) < 20; continue; end
        e = @(v) median(v(idx(end-round(0.15*numel(idx))):idx(end))) - v(idx(1));
        fprintf('%-16s %10.2f %10.2f %10.2f %10.1f\n', sprintf('%s c%d',nm,c), ...
            e(F), e(FACT), e(FPAS), 100*e(p1+p2));
    end
end

save(fullfile(here,'results','recovery_anatomy.mat'), 'SNAP');
fprintf('\nSaved results/recovery_anatomy.mat\n');
