
function [cost_total, weight, cost] = evalFeatureCost(feats_data, feats_sim, fn)

for i_feat = 1:size(fn, 2)

    featnames = split(fn{i_feat}, '|');
    feat_y{i_feat} = featnames{1};
    
    if length(featnames ) >= 2 && featnames{2} ~= "_" && ~isempty(str2num(featnames{2}))
        weight(i_feat) = str2num(featnames{2});
    elseif length(featnames) >= 3 && ~isempty(str2num(featnames{3}))
        weight(i_feat) = str2num(featnames{3});
    else
        weight(i_feat) = 1;
    end

    fd = [feats_data.(feat_y{i_feat})];
    fs = [feats_sim.(feat_y{i_feat})];

    % normalization factor
    n_F = mean(fd);

    cost(i_feat) = sum((fd/n_F - fs/n_F).^2);
    cost(i_feat) = sum(abs(fd/n_F - fs/n_F));

end

cost_total = cost.*weight;