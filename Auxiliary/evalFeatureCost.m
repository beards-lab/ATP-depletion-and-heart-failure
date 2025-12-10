
function [cost_total, weight, cost] = evalFeatureCost(feats_data, feats_sim, fn, costExp)
% Fn is a cell array of string features. It can combine multiple ones by |.
% Last one can weight or disable its cost, e.g.
% feature1|features2|0 - evaluates feature 1, but does not count becuase 0 weight

if nargin < 4
    costExp = 1;
end
for i_feat = 1:size(fn, 2)

    featnames = split(fn{i_feat}, '|');
    feat_y{i_feat} = featnames{1};
    
    % if length(featnames ) >= 2 && featnames{2} ~= "_" && ~isempty(str2num(featnames{2}))
    %     weight(i_feat) = str2num(featnames{2});
    % elseif length(featnames) >= 3 && ~isempty(str2num(featnames{3}))
    %     weight(i_feat) = str2num(featnames{3});
    % else
    if isempty(str2num(featnames{end}))
        weight(i_feat) = 1;
    else
        weight(i_feat) = str2num(featnames{end});
    end

    fd = [feats_data.(feat_y{i_feat})];
    fs = [feats_sim.(feat_y{i_feat})];

    % normalization factor
    n_F = mean(fd);

    % cost(i_feat) = sum((fd/n_F - fs/n_F).^2);
    % cost(i_feat) = sum(abs(fd/n_F - fs/n_F));
    cost(i_feat) = sum(abs(fd/n_F - fs/n_F).^costExp);

end

cost_total = cost.*weight;