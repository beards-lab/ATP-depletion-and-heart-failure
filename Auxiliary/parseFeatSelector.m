function [name, selector, rest] = parseFeatSelector(token)
% PARSEFEATSELECTOR  Strip an optional [selector] suffix from a feature-name token.
%
%   [NAME, SELECTOR, REST] = PARSEFEATSELECTOR(TOKEN) parses TOKEN of the form
%   'feat[selector]' or 'feat'. If a bracket is present, NAME is the bare
%   feature name, SELECTOR is the contents of [...] and REST is anything after
%   the closing bracket. If no bracket is present, NAME=TOKEN, SELECTOR='',
%   REST=''.
%
%   Examples:
%     parseFeatSelector('XTOR[3]')     -> 'XTOR',  '3',  ''
%     parseFeatSelector('XTOR[mean]')  -> 'XTOR', 'mean',''
%     parseFeatSelector('XTOR')        -> 'XTOR',  '',   ''
%
%   See also: applyFeatSelector, evalFeatureCost, plotMultipleFeatures.

m = regexp(token, '^([A-Za-z_]\w*)\[([^\]]+)\](.*)$', 'tokens', 'once');
if isempty(m)
    name = token;
    selector = '';
    rest = '';
else
    name = m{1};
    selector = m{2};
    rest = m{3};
end
end
