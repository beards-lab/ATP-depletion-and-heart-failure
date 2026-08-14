% RunParadoxProbe.m — "parking heads in SRX must mean less cycling, so how can
% the recovery get FASTER?"
%
% The question has two halves and this probe measures both directly, on the real
% sequential protocol (RunSlackSegments = 'All', so cycle-to-cycle carryover is
% real):
%
%   (1) HOW MUCH mass is actually diverted?  int RdetSRX dt over the restretch,
%       compared against the size of the PD reservoir the recovery draws from.
%       If PD is a big pre-filled buffer, diverting a little mass out of it
%       cannot slow the force rise — it can only lower the plateau.
%
%   (2) HOW LONG is the detour, at the force where it happens?  The normal route
%       out of p2 pays 1/kah through PT. The SRX route pays 1/ksr2srd +
%       1/(kmsrd*exp(F/sigma_srd1)). The second term is force-ACCELERATED, and
%       the post-restretch window is the only recovery window that runs at high
%       force. So the comparison must be made at the measured F, not at F = 0.
%
% Also records, for the "why is post-restretch fast in the first place" question,
% the state of the machine at each window start: how much bound mass SURVIVED,
% how full PD is, and how the SRX gates are set by the prevailing force.
%
% Variants: baseline, B2 (rip -> SRXD), D0 s0=.010 k=15 (force-flat slow pool).

clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
if isempty(gcp('nocreate')); parpool('local', 5); end

S   = load(fullfile(root,'data','protocol_03_27_2026_8mM_slack.mat'));
vt  = S.velocitytable;
iS  = find(vt(:,2) < -1);
iR  = find(vt(:,2) > 1);

pBase = getParams(loadParams('params/rskR2_w025_opt.m'), [], true, false);
pBase.MaxRunTime = 1200;
pBase.RunSlackSegments = 'All';     % sequential: real carryover
pBase.EvalFeatures = 1; pBase.PlotEachSeparately = 0; pBase.PlotFeatureFitting = 0;

% NB 'rip->SRXT' is the route originally asked for: p2 -> P_SR under stretch. It is
% the SLOWER of the two SRX routes to PD (13.3 ms vs 4.2 ms at F0, because it must
% cross SR -> SRD at ksr2srd = 110/s first), so it is the one that ought to work —
% which is exactly why it needs measuring rather than arguing about.
VAR = { 'baseline',        struct()
        'rip->SRXT .015',  struct('SRXFromR2HighStrain',1,'sSRXrip',0.015,'dSRXrip',0.005,'SRXRecoilToD',false)
        'rip->SRXD .015',  struct('SRXFromR2HighStrain',1,'sSRXrip',0.015,'dSRXrip',0.005,'SRXRecoilToD',true)
        'rip->D0 .010 k15',struct('SRXFromR2HighStrain',1,'sSRXrip',0.010,'dSRXrip',0.005, ...
                                  'SRXRecoilToD0',true,'UseD0State',true,'k_D0T',15,'K_D0T',0.15) };

P = struct();
for v = 1:size(VAR,1)
    nm = VAR{v,1}; over = VAR{v,2};
    p = pBase; f = fieldnames(over);
    for i = 1:numel(f); p.(f{i}) = over.(f{i}); end
    p = getParams(p, [], true, false);

    fprintf('\n================== %s ==================\n', nm);
    tic; [~, o] = runSlackExperiment(p); fprintf('(%.0f s)\n', toc);

    key = matlab.lang.makeValidName(nm);
    P.(key).name = nm; P.(key).p = p;

    Fiso = interp1(o.t, o.Force, vt(iS(1),1)-0.02, 'linear', 'extrap');
    has0 = isfield(o,'D0');
    hasR = isfield(o,'RdetSRX');

    fprintf('%-4s %7s %7s %7s %7s %7s %7s %7s %8s %8s\n', 'cyc', 'F0/Fis', ...
            'bound0', 'PD0', 'PT0', 'SR0', 'SRD0', 'D0_max', 'divert', 'k63');
    for c = 1:4
        % window = vall2 (min of the 80 ms after the restretch ramp ends) -> next release
        t_r1 = vt(iR(c)+1,1);
        t_end = vt(iS(c+1),1);
        m = o.t >= t_r1 & o.t <= t_end;
        j = find(m); tt = o.t(j); FF = o.Force(j)/Fiso;
        mv = tt <= tt(1)+0.080; [~, iv] = min(FF(mv));
        i0 = j(iv); i1 = j(end);

        % k63 with the shared crossing rule
        seg = i0:i1; B = o.Force(i0)/Fiso;
        Pl = median(o.Force(seg(o.t(seg) >= o.t(i1)-0.15*(o.t(i1)-o.t(i0))))/Fiso);
        y  = (movmedian(o.Force(seg)/Fiso, max(3, round(0.002/median(diff(o.t(seg)))))) - B)/max(Pl-B, eps);
        i63 = find(y >= 1-exp(-1), 1);
        k63 = NaN; if ~isempty(i63); k63 = 1/(o.t(seg(i63)) - o.t(i0)); end

        % diverted mass: integral of the routed flux over the RESTRETCH RAMP
        mr = o.t >= vt(iR(c),1) & o.t <= t_r1 + 0.010;
        div = 0;
        if hasR; div = trapz(o.t(mr), o.RdetSRX(mr)); end
        d0m = NaN; if has0; d0m = max(o.D0(m | mr)); end

        fprintf('%-4d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %8.4f %8.1f\n', c, B, ...
            o.p1_0(i0)+o.p2_0(i0), o.PuR(i0), o.PuATP(i0), o.SR(i0), o.SRD(i0), d0m, div, k63);
    end

    % transit-time arithmetic at the measured window-start force
    F0 = o.Force(i0);
    kOutSRD = p.kmsrd * exp(F0/p.sigma_srd1);
    fprintf('\n  at F0 = %.1f kPa (%.2f Fiso):\n', F0, F0/Fiso);
    fprintf('    normal   p2->PT->PD          1/kah      = %6.2f ms\n', 1000/p.kah);
    fprintf('    SRX.ADP  p2->SRD->PD         1/%.0f      = %6.2f ms   (kmsrd*exp(F/sg)=%.0f/s)\n', ...
            kOutSRD, 1000/kOutSRD, kOutSRD);
    fprintf('    SRX.ATP  p2->SR->SRD->PD                = %6.2f ms\n', ...
            1000/p.ksr2srd + 1000/kOutSRD);
    fprintf('    at F = 0 the same SRX.ADP exit is %6.2f ms  (kmsrd = %.1f/s)\n', ...
            1000/p.kmsrd, p.kmsrd);
    fprintf('    SRX ENTRY from PT: ksr0*exp(-F/sg2) = %5.0f/s at F0, %5.0f/s at F=0\n', ...
            p.ksr0*exp(-F0/p.sigma2), p.ksr0);

    % Where does the total detached mass sit through the restretch? If the routed
    % heads are genuinely PARKED, the sum PT+SR+SRD must rise and PD must fall.
    % If instead the SRX route is just a different corridor to the same place, PD
    % holds up and nothing is actually withheld from attachment.
    mrw = o.t >= vt(iR(1),1) - 0.005 & o.t <= vt(iR(1)+1,1) + 0.060;
    tw  = o.t(mrw) - vt(iR(1),1);
    fprintf('\n  through restretch cycle 1 (t rel. ramp start, ms):\n');
    fprintf('    %6s %7s %7s %7s %7s %8s %8s\n', 't', 'PD', 'PT', 'SR', 'SRD', 'detach', 'bound');
    for tq = [0 5 10 20 40 60]
        [~, iq] = min(abs(tw - tq/1000)); iq = find(mrw, 1) - 1 + iq;
        det = o.PuR(iq) + o.PuATP(iq) + o.SR(iq) + o.SRD(iq);
        fprintf('    %6d %7.3f %7.3f %7.3f %7.3f %8.3f %8.3f\n', tq, o.PuR(iq), ...
            o.PuATP(iq), o.SR(iq), o.SRD(iq), det, o.p1_0(iq)+o.p2_0(iq));
    end

    P.(key).out = struct('t', o.t, 'Force', o.Force, 'PuR', o.PuR, 'PuATP', o.PuATP, ...
                         'SR', o.SR, 'SRD', o.SRD, 'p1_0', o.p1_0, 'p2_0', o.p2_0, ...
                         'Fiso', Fiso);
    if hasR; P.(key).out.RdetSRX = o.RdetSRX; end
    if has0; P.(key).out.D0 = o.D0; end
end

save(fullfile(here,'results','paradox_probe.mat'), 'P', 'vt', 'iS', 'iR');
fprintf('\nsaved results/paradox_probe.mat\n');
