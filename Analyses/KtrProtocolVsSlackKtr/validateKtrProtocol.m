function V = validateKtrProtocol(tag, verbose)
%VALIDATEKTRPROTOCOL Held-out check of a paramset against the ktr protocol.
%
%   V = validateKtrProtocol(TAG)            e.g. validateKtrProtocol('fastgate_best')
%   V = validateKtrProtocol(TAG, false)     silent
%
%   The ktr protocol is deliberately NOT in the objective (RunKtr = 0 and no
%   ktrP_* in params0.fn), so it is a genuine out-of-sample test. This runs it
%   on all three protocol days and reports the redevelopment rate, the trace
%   error, and the trace error in units of the DATA'S OWN NOISE FLOOR -- the RMS
%   residual of the measurement about its best single exponential, which is the
%   only scale on which "the model fits this" means anything.
%
%   It also reports the post-slack #2 rate, because the point of the comparison
%   is that both windows redevelop at SL 2.00 um from zero and the data gives
%   them ONE rate (ratio 0.97). See conclusions.md.
%
%   V fields: .day .k63 .k63_data .Ektr .rmse .noise .ratio .kS2 .kS2_data

    if nargin < 2; verbose = true; end
    root = 'C:\home\git\ATP-depletion-and-heart-failure';
    DAYS = {'03_27','04_03','04_10'};

    p = getParams(loadParams(tag), [], true, false);
    p.PlotEachSeparately = 0; p.ghostSave = ''; p.ghostLoad = '';
    p.EvalFeatures = 0; p.BreakOnODEUnstable = 0; p.SaveBest = 0;
    p = getParams(p, p.g, true);

    % --- the slack side, once (the protocol file the paramset is fitted to) ---
    svt = load(fullfile(root,'data',p.velocitytableonfile)).velocitytable;
    iS  = find(svt(:,2) < -1);
    [~, os] = runSlackExperiment(p);
    FisoM = median(os.Force(os.t >= svt(iS(2),1)-0.05 & os.t <= svt(iS(2),1)-0.005));
    ws = recoveryWindows(os.t, os.Force, svt, 'slack', FisoM);
    ws = ws(strcmp({ws.type},'postSlack'));
    kS2 = ws(2).k63;

    V = struct('day',{},'k63',{},'k63_data',{},'Ektr',{},'rmse',{},'noise',{},'ratio',{});
    for d = 1:numel(DAYS)
        K = load(fullfile(root,'data',sprintf('protocol_%s_2026_velocitytable_ktr.mat',DAYS{d})));
        kvt = K.velocitytable;
        pk = p; pk.SL0 = 2*kvt(1,4); pk = getParams(pk, pk.g, true); pk.Velocity = kvt(:,2);
        [~, ok] = evaluateModel(str2func(pk.modelFcn), kvt(:,1), pk);
        FisoK = interp1(ok.t, ok.Force, kvt(2,1), 'linear', 'extrap');

        fm = extractKtrAttributes(ok.t, ok.Force/FisoK, kvt);
        fd = extractKtrAttributes(K.datatable(:,1), K.datatable(:,3), kvt);

        Fi   = interp1(ok.t, ok.Force/FisoK, K.datatable(:,1), 'linear','extrap');
        Ektr = mean((K.datatable(:,3) - Fi).^2) * 10;

        m = K.datatable(:,1) >= kvt(end-1,1) & K.datatable(:,1) <= kvt(end,1);
        rmse  = sqrt(mean((Fi(m) - K.datatable(m,3)).^2));
        noise = iNoiseFloor(K.datatable(m,1), K.datatable(m,3));

        V(d) = struct('day',DAYS{d},'k63',fm.ktrP_k63,'k63_data',fd.ktrP_k63, ...
                      'Ektr',Ektr,'rmse',rmse,'noise',noise,'ratio',fm.ktrP_k63/kS2);
    end
    [V.kS2] = deal(kS2); [V.kS2_data] = deal(50.0);   % 03/27 post-slack #2, k63

    if verbose
        fprintf('\n===== HELD-OUT ktr protocol: %s  (tau_reg = %.4g ms, RunKtr = %d) =====\n', ...
                tag, p.tau_reg*1e3, p.RunKtr);
        fprintf('%-8s %9s %9s %9s %9s %11s\n','day','k63 model','k63 data','ratio','traceMSE','rmse/noise');
        for d = 1:numel(V)
            fprintf('%-8s %9.1f %9.1f %9.2f %9.4f %10.2fx\n', V(d).day, V(d).k63, V(d).k63_data, ...
                    V(d).k63/V(d).k63_data, V(d).Ektr, V(d).rmse/V(d).noise);
        end
        fprintf('%-8s %9.1f %9.1f %9.2f %9.4f\n','MEAN', mean([V.k63]), mean([V.k63_data]), ...
                mean([V.k63])/mean([V.k63_data]), mean([V.Ektr]));
        fprintf('post-slack #2 (SL 2.00 um, IN the objective): model %.1f, data %.1f\n', kS2, 50.0);
        fprintf('k_ktr / k_postSlack2 = %.2f   (data 0.99 -- one length, one rate)\n', ...
                mean([V.k63])/kS2);
    end
end

function r = iNoiseFloor(t, y)
% RMS residual of the measurement about its own best B + A(1-exp(-k(t-t0))).
    t = t(:) - t(1); y = y(:); B = min(y); A0 = max(y) - B;
    obj = @(q) sum((y - (B + abs(q(1)).*(1-exp(-abs(q(2)).*max(t-max(q(3),0),0))))).^2);
    q = fminsearch(obj, [A0, 45, 0.005], optimset('Display','off','MaxFunEvals',8e3,'MaxIter',8e3));
    pred = B + abs(q(1)).*(1-exp(-abs(q(2)).*max(t-max(q(3),0),0)));
    r = sqrt(mean((y - pred).^2));
end
