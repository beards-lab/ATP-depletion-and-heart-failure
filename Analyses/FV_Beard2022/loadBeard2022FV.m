function D = loadBeard2022FV(dataDir, useCache)
% loadBeard2022FV  Load the Beard/Baker 2022 (Biophys J) isovelocity workbooks.
%
%   D = loadBeard2022FV()            reads ../../data/Beard2022/*/*isovelocity*.xlsx
%   D = loadBeard2022FV(dataDir)     reads another folder of prep subfolders
%   D = loadBeard2022FV(dataDir, false) re-reads the xlsx instead of the cache
%
% One element of D per (prep, ATP sheet). Each workbook sheet is ONE continuous
% motor record (Time ms, L in Lo, F in kPa; 0.05 ms) holding the whole
% isovelocity sequence: isometric at ~1.10 Lo, ramp-shorten at v ML/s through
% 1.00 Lo, re-stretch to 1.10, repeat for each v, and finally a step to 1.00 Lo
% for the isometric (v = 0) point. The lab's own summary (columns F..J) gives,
% per v: ML/s, iso-V force (kPa, read as the ramp passes iso-V Length), zero
% (force baseline), iso-V Length.
%
% Fields: prep (folder name), date, order ('8-4' ...), first (ATP of the first
% sheet run), atp (mM), slot (1 = first run of the prep, 2 = second), t, L, F
% (column vectors), summary table (v, Fiso_v, zero, Liso), file.

if nargin < 1 || isempty(dataDir)
    here = fileparts(mfilename('fullpath'));
    dataDir = fullfile(here, '..', '..', 'data', 'Beard2022');
end
if nargin < 2, useCache = true; end
cacheFile = fullfile(fileparts(mfilename('fullpath')), 'results', 'beard2022_fv_raw.mat');
if useCache && isfile(cacheFile)
    S = load(cacheFile, 'D'); D = S.D; return
end

preps = dir(dataDir); preps = preps([preps.isdir] & ~startsWith({preps.name}, '.'));
D = struct('prep', {}, 'date', {}, 'order', {}, 'atp', {}, 'slot', {}, ...
    't', {}, 'L', {}, 'F', {}, 'summary', {}, 'file', {});
for ip = 1:numel(preps)
    pdir = fullfile(dataDir, preps(ip).name);
    f = dir(fullfile(pdir, '*isovelocity*.xlsx'));
    assert(numel(f) == 1, 'expected one isovelocity workbook in %s', pdir);
    fn = fullfile(pdir, f.name);
    % ATP order from the folder name, e.g. "2021 05 11 - 8 mM then 4 mM"
    tok = regexp(lower(preps(ip).name), '(\d+)\s*mm\s*then\s*(\d+)\s*mm', 'tokens', 'once');
    assert(~isempty(tok), 'cannot parse ATP order from %s', preps(ip).name);
    order = str2double(tok);
    dt = regexp(preps(ip).name, '^(\d{4}) (\d{2}) (\d{2})(B?)', 'tokens', 'once');
    sheets = sheetnames(fn);
    for is = 1:numel(sheets)
        a = str2double(regexp(lower(sheets(is)), '^\s*(\d+)\s*mm', 'tokens', 'once'));
        if isnan(a), continue; end
        M = readmatrix(fn, 'Sheet', sheets(is), 'Range', 'A5:C200000');
        M = M(all(~isnan(M), 2), :);
        Sm = readmatrix(fn, 'Sheet', sheets(is), 'Range', 'F2:J40');
        Sm = Sm(~isnan(Sm(:, 1)) & ~isnan(Sm(:, 2)), :);
        k = numel(D) + 1;
        D(k).prep = preps(ip).name;
        D(k).date = sprintf('%s-%s-%s%s', dt{:});
        D(k).order = sprintf('%d-%d', order);
        D(k).atp = a;
        D(k).slot = find(order == a, 1);
        D(k).t = M(:, 1); D(k).L = M(:, 2); D(k).F = M(:, 3);
        D(k).summary = table(Sm(:, 1), Sm(:, 2), Sm(:, 3), Sm(:, 5), ...
            'VariableNames', {'v', 'Fiso_v', 'zero', 'Liso'});
        D(k).file = fn;
    end
end
if ~isfolder(fileparts(cacheFile)), mkdir(fileparts(cacheFile)); end
save(cacheFile, 'D', '-v7.3');
end
