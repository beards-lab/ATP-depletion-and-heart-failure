function outFile = buildHtmlReport(templateFile, outFile, subs)
%BUILDHTMLREPORT  Fill a {{PLACEHOLDER}} HTML template and inline every figure.
%
%   OUT = BUILDHTMLREPORT(TEMPLATEFILE, OUTFILE, SUBS)
%
%   SUBS is a struct whose field names are placeholder names without the braces:
%   field FOO replaces every occurrence of {{FOO}}. Values may be
%
%     char / string   inserted verbatim (already-escaped HTML is fine)
%     numeric scalar  inserted with %g
%     a path to an existing image file  -> read and inlined as a
%                                          data:image/<ext>;base64,... URI
%
%   Inlining is what makes the output a SINGLE self-contained file that can be
%   mailed or dropped in a shared folder without its figures going missing.
%   The cost is ~1.35x the PNG bytes; a 10-figure report lands around 3-5 MB.
%
%   Placeholders left unfilled are reported (not silently blanked) so a typo in
%   a template shows up as an error rather than a hole in the page.
%
%   Example:
%       s.DATE = datestr(now,'yyyy-mm-dd');
%       s.FIG_SLACK = 'results/fit_slack_full.png';
%       s.TOTAL = 1.6615;
%       buildHtmlReport('report_template.html','results/report.html', s);
%
%   See also: Analyses/OptimResultReport/BuildReport.m

    if ~isfile(templateFile)
        error('buildHtmlReport:noTemplate','template not found: %s', templateFile);
    end
    html = fileread(templateFile);

    IMG_EXT = {'.png','.jpg','.jpeg','.gif','.svg','.webp'};
    MIME    = struct('png','image/png','jpg','image/jpeg','jpeg','image/jpeg', ...
                     'gif','image/gif','svg','image/svg+xml','webp','image/webp');

    f = fieldnames(subs);
    nImg = 0; imgBytes = 0;
    for i = 1:numel(f)
        v = subs.(f{i});
        if isnumeric(v) && isscalar(v)
            rep = num2str(v, '%g');
        elseif isstring(v) || ischar(v)
            v = char(v);
            [~,~,ext] = fileparts(v);
            if any(strcmpi(ext, IMG_EXT)) && isfile(v)
                d = dir(v); imgBytes = imgBytes + d.bytes; nImg = nImg + 1;
                rep = iDataUri(v, MIME.(lower(ext(2:end))));
            else
                rep = v;
            end
        else
            error('buildHtmlReport:badValue', ...
                  'field %s must be char, string or numeric scalar (got %s)', f{i}, class(v));
        end
        % replace() is literal -- no regexprep, so $ and \ in a value are safe
        html = replace(html, ['{{' f{i} '}}'], rep);
    end

    left = unique(regexp(html, '\{\{[A-Za-z0-9_]+\}\}', 'match'));
    if ~isempty(left)
        error('buildHtmlReport:unfilled', ...
              'template has %d unfilled placeholder(s): %s', numel(left), strjoin(left, ' '));
    end

    d = fileparts(outFile);
    if ~isempty(d) && ~isfolder(d); mkdir(d); end
    fid = fopen(outFile, 'w', 'n', 'UTF-8');
    if fid < 0; error('buildHtmlReport:cannotWrite','cannot write %s', outFile); end
    fwrite(fid, unicode2native(html, 'UTF-8')); fclose(fid);

    di = dir(outFile);
    fprintf('buildHtmlReport: %s  (%.1f MB, %d figures inlined from %.1f MB of PNG)\n', ...
            outFile, di.bytes/1e6, nImg, imgBytes/1e6);
end

function uri = iDataUri(file, mime)
    fid = fopen(file, 'r');
    raw = fread(fid, Inf, '*uint8'); fclose(fid);
    uri = ['data:' mime ';base64,' matlab.net.base64encode(raw)];
end
