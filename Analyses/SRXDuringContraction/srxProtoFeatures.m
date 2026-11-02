% srxProtoFeatures  the candidate identification protocols for one parameter set:
%   B  release 0.15 ML, hold T, restretch (srxProtocol 'B'), T = TB -> k63(T), F20(T)
%   A  activation step from the relaxed steady state                -> t50, t90
% Runs the 6 simulations in parallel (needs a pool). Starts from q's own isometric SRX.
%   P = struct('k63', 1xnT, 'F20', 1xnT, 't50', t90, 'srx', 'F', 'TB')
% See RunSrxTwin.m.
function P = srxProtoFeatures(q, TB)
    s = srxIsoState(q); q.SRXT_0 = s.SR; q.SRXD_0 = s.SRD;
    q.MaxRunTime = 1200;              % the short-hold runs are stiff for some knob sets
    jobs = [num2cell(TB), {0}]; kinds = [repmat({'B'}, 1, numel(TB)), {'A'}];
    R = cell(1, numel(jobs));
    parfor j = 1:numel(jobs)
        try
            R{j} = srxProtocol(q, kinds{j}, jobs{j});
        catch
            R{j} = struct('k63', NaN, 'F20', NaN, 't50', NaN, 't90', NaN);   % excluded downstream
        end
    end
    nT = numel(TB);
    P = struct('k63', cellfun(@(r) r.k63, R(1:nT)), 'F20', cellfun(@(r) r.F20, R(1:nT)), ...
               't50', R{end}.t50, 't90', R{end}.t90, 'srx', s.SRX, 'F', s.F, 'TB', TB);
end
