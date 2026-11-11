% srxResiduals  model features -> residuals vs the data in data-SD units.
%   [z, names, labels, m, d, sd] = srxResiduals(fm, fd)
% Slack block: pooled 8 mM target (segment means, relative between-prep SD);
% ktr protocol: pooledKtrTarget(8); F-V: normalised power vs the Baker sets,
% SD = max(across-set SD, 13.6 %). See RunSrxCompensation.m.
function [z, names, labels, m, d, sd] = srxResiduals(fm, fd)
    persistent PT8 KT8
    if isempty(PT8); PT8 = pooledTarget8mM(); KT8 = pooledKtrTarget(8); end
    FL = { 'steady',  'steady force',            1:5, 'P'
           'A',       'redeveloped force A',     1:5, 'P'
           'peak1_y', 're-stretch peak',         1:5, 'P'
           'vall_y',  'valley',                  1:5, 'P'
           'peak2',   'second peak',             1:5, 'P'
           'ktrfit',  'slack k_{tr} (fit)',      1:5, 'P'
           'ktr80',   'slack k_{tr}80',          1:5, 'P'
           't0',      'slack delay t_0',         1:5, 'P'
           'rsK80',   'post-restretch rate',     1:4, 'P'
           'vall_t',  'valley time',             1:5, 'P'
           'ktrP_k80','k_{tr} protocol k80',     1,   'K'
           'ktrP_k63','k_{tr} protocol k63',     1,   'K'
           'ktrP_F20','k_{tr} protocol F(20 ms)',1,   'K'
           'FVp1',    'F-V power 0.5 ML/s',      1,   'F'
           'FVp2',    'F-V power 1 ML/s',        2,   'F'
           'FVp3',    'F-V power 2 ML/s',        3,   'F'
           'FVp4',    'F-V power 4 ML/s',        4,   'F' };
    n = size(FL, 1); m = nan(n, 1); d = nan(n, 1); sd = nan(n, 1);
    for k = 1:n
        f = FL{k, 1}; sg = FL{k, 3};
        switch FL{k, 4}
            case 'P'
                j = find(strcmp(PT8.name, f), 1); d(k) = mean(PT8.val(j, sg)); sd(k) = abs(d(k)) * PT8.s(j);
                if isfield(fm, f); x = fm.(f); m(k) = mean(x(sg(sg <= numel(x)))); end
            case 'K'
                j = find(strcmp(KT8.name, f), 1); d(k) = KT8.val(j); sd(k) = abs(d(k)) * KT8.s(j);
                if isfield(fm, f); m(k) = fm.(f)(1); end
            case 'F'
                if isfield(fd, 'FV_normpowerAvg')
                    d(k) = fd.FV_normpowerAvg(sg);
                    sd(k) = max(sqrt(fd.FV_normpowerVar(sg)), 0.136 * abs(d(k)));
                end
                if isfield(fm, 'FV_normpowerAvg'); m(k) = fm.FV_normpowerAvg(sg); end
        end
    end
    z = (m - d) ./ sd; names = FL(:, 1); labels = FL(:, 2);
end
