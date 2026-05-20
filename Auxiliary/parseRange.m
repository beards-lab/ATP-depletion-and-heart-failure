function [lb, ub] = parseRange(tok)
%PARSERANGE  Parse a 'lb-ub' range token into two scalars.
%
%   [LB, UB] = PARSERANGE(TOK) splits TOK on the first digit-preceded dash
%   and returns the lower and upper bound as doubles.
%
%   See also: isRangeToken, evalFeatureCost, plotMultipleFeatures

dash = regexp(tok, '(?<=\d)-', 'once');
lb   = str2double(tok(1:dash-1));
ub   = str2double(tok(dash+1:end));
end
