function logWindowShift(fcn, ts, ie, nds, PU0, params)
%LOGWINDOWSHIFT  Record what one strain-window shift in evaluateModel deletes.
%
%   Called from evaluateModel just before the arrays are shifted, only when
%   params.LogWindowShifts is true. Appends one record to the global
%   WINDOW_SHIFT_LOG: the attached mass and the active force carried by the
%   nds bins about to be zeroed, and the strains those bins sit at.
%   Serial runs only -- parallel workers keep their own globals.
%
%   See also: evaluateModel
global WINDOW_SHIFT_LOG
ss = params.ss;  Ns = params.NumberOfStates;  dS = params.dS;
PU0 = PU0(:);
SL = PU0(Ns*ss + 3);  LSE = PU0(Ns*ss + 4);
if params.LegacyStrainFlipping
    s = flipud(-(params.s(:) + (-(SL - LSE) + params.LXBpivot)/2));
else
    s = params.s(:) - (-(SL - LSE) + params.LXBpivot)/2;
end
if ie == 1, del = 1:nds; else, del = ss-nds+1:ss; end
P = reshape(PU0(1:Ns*ss), ss, Ns);
bound = dS*sum(P(:));
mass  = dS*sum(sum(P(del, :)));
PUz = PU0;  for k = 1:Ns, PUz((k-1)*ss + del) = 0; end
[~, o]  = fcn(ts, PU0, params);
[~, oz] = fcn(ts, PUz, params);
pb = sum(P, 2);
occ = find(pb > 1e-3*max(pb));
rec = struct('t', ts, 'ie', ie, 'nds', nds, 'bound', bound, 'mass', mass, ...
    'massFrac', mass/max(bound, eps), 'Fact', o(2), 'dFact', o(2) - oz(2), ...
    'sDel', [min(s(del)) max(s(del))], 'sMassMean', sum(s(del).*sum(P(del,:),2))/max(sum(sum(P(del,:))), eps), ...
    'sOccupied', [min(s(occ)) max(s(occ))]);
if isempty(WINDOW_SHIFT_LOG), WINDOW_SHIFT_LOG = rec; else, WINDOW_SHIFT_LOG(end+1) = rec; end
end
