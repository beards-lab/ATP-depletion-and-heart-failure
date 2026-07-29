% RunFluxLedger.m
% Tabulate windowFluxLedger over every captured window. Zero simulation —
% reads results/window_capture.mat written by captureWindows.m.
%
% The question: what fraction of each recovery is NEW ATTACHMENT vs
% REDISTRIBUTION of cross-bridges that survived the preceding event, and does
% the SRX trap explain the post-slack / post-restretch rate split?

root = 'C:\home\git\ATP-depletion-and-heart-failure';
cd(root); addpath(genpath(root));
resdir = fullfile(root,'Analyses','RestretchVsKtrRecovery','results');
C = load(fullfile(resdir,'window_capture.mat'));

L = struct([]);
for q = 1:numel(C.W)
    w = C.W(q);
    if strcmp(w.src,'ok'); o = C.ok; pp = C.pk; else; o = C.os; pp = C.p; end
    Lq = windowFluxLedger(o, w, pp);
    if isempty(L); L = Lq; else; L(end+1) = Lq; end %#ok<SAGROW>
end

fprintf('\n================= POOL STATE AT WINDOW START =================\n');
fprintf('%-16s %7s | %6s %6s %6s %6s %6s | %6s\n', ...
        'window','k63','F/Fiso','PD','PT','p1_0','p2_0','SRX');
for q = 1:numel(L)
    s = L(q).state; fi = C.W(q).Fiso;
    fprintf('%-16s %7.1f | %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f%s\n', ...
        L(q).name, L(q).k63, s.F/fi, s.PD, s.PT, s.p1, s.p2, L(q).srx.at_i0, ...
        tern(C.W(q).excluded,'  (excl)',''));
end

fprintf('\n================= FLUX LEDGER over [t0, t63] =================\n');
fprintf('%-16s %6s | %7s %7s %7s %7s %7s | %6s %6s %6s | %-10s %s\n', ...
        'window','dur_ms','N_att','N_det1','N_12','N_21','N_det2','nu','I','T','closure','hold');
for q = 1:numel(L)
    c = L(q).closure; g = L(q).guard;
    fprintf('%-16s %6.1f | %7.4f %7.4f %7.4f %7.4f %7.4f | %6.3f %6.2f %6.2f | %-10s v=[%+.1e %+.1e]%s\n', ...
        L(q).name, L(q).dur_ms, L(q).N_att, L(q).N_det1, L(q).N_12, L(q).N_21, L(q).N_det2, ...
        L(q).nu, L(q).I, L(q).T, ...
        tern(c.ok, 'ok', sprintf('*** %.0e', max(abs([c.p1 c.p2])))), ...
        g.vmin, g.vmax, tern(g.hopRisk, ' hop?', ''));
end

fprintf('\n================= SRX SUB-LEDGER (the H_SRX test) =================\n');
fprintf('%-16s | %8s %8s | %8s %8s | %7s %7s %7s\n', ...
        'window','entry/s','exit/s','N_PT2SR','N_SRD2PD','SRX_t0','SRX_pk','excurs');
for q = 1:numel(L)
    s = L(q).srx;
    fprintf('%-16s | %8.1f %8.1f | %8.4f %8.4f | %7.3f %7.3f %+7.3f\n', ...
        L(q).name, s.entryRate0, s.exitRate0, s.N_PT2SR, s.N_SRD2PD, ...
        s.at_i0, s.peak, s.excursion);
end

%% ---- headline comparison: the SL-matched pair + ktr ----
sel = @(nm) L(strcmp({L.name}, nm));
a = sel('postSlack2'); b = sel('postRestretch2'); k = sel('ktr');
fprintf('\n================= HEADLINE =================\n');
fprintf('%-22s %10s %10s %10s\n','', 'postSlack2','postRestr2','ktr');
row = @(lbl,f) fprintf('%-22s %10.3f %10.3f %10.3f\n', lbl, f(a), f(b), f(k));
fprintf('%-22s %10.1f %10.1f %10.1f\n','k63 (1/s)', a.k63, b.k63, k.k63);
row('new-head share nu',      @(x) x.nu);
row('internal traffic I',     @(x) x.I);
row('turnover T',             @(x) x.T);
row('N_att',                  @(x) x.N_att);
row('SRX entry rate @t0',     @(x) x.srx.entryRate0);
row('SRX exit  rate @t0',     @(x) x.srx.exitRate0);
row('N_PT2SR (into trap)',    @(x) x.srx.N_PT2SR);
row('SRX excursion',          @(x) x.srx.excursion);

save(fullfile(resdir,'flux_ledger.mat'),'L');
fprintf('\nsaved results/flux_ledger.mat\n');

function o = tern(c,a,b); if c; o=a; else; o=b; end; end
