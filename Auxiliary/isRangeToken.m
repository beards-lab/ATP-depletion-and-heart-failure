function tf = isRangeToken(tok)
%ISRANGETOKEN  True if TOK is a 'lb-ub' range string (both ends parse as finite numbers).
%
%   Uses a lookbehind to find a dash that is preceded by a digit, so that
%   negative lower bounds (e.g. '-1-5') are handled correctly.
%
%   See also: parseRange, evalFeatureCost, plotMultipleFeatures

dash = regexp(tok, '(?<=\d)-', 'once');
if isempty(dash); tf = false; return; end
lb = str2double(tok(1:dash-1));
ub = str2double(tok(dash+1:end));
tf = isfinite(lb) && isfinite(ub);
end
