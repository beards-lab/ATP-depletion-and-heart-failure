function md2html(mdfile, htmlfile, pageTitle, embedImages)
%MD2HTML  Convert a Markdown file to a standalone, self-styled HTML page.
%
%   md2html(MDFILE) writes alongside MDFILE with the extension swapped.
%   md2html(MDFILE, HTMLFILE) writes to HTMLFILE.
%   md2html(MDFILE, HTMLFILE, PAGETITLE) also sets <title>.
%
%   There is no pandoc/python/node on this machine, so this exists to keep the
%   md -> html step reproducible in-repo. It covers the subset this project's
%   reports actually use:
%
%     # .. ######   headings          |  ```fenced``` code blocks
%     **bold**  *em*  `code`          |  > blockquotes
%     - / * bullets, 1. numbered      |  --- horizontal rules
%     | pipe | tables | with headers  |  ![alt](img)  and  [text](url)
%
%   Consecutive non-blank lines are joined into ONE paragraph before inline
%   markup is applied, as Markdown requires -- otherwise **bold** or a [link]
%   that wraps across a line break is silently left unconverted.
%
%   md2html(..., EMBEDIMAGES) controls whether local images are inlined as
%   base64 data URIs. Default TRUE, which makes the page a SINGLE self-contained
%   file that survives being moved, emailed, or opened from a sandbox that
%   cannot resolve relative paths. Pass FALSE to keep <img src="..."> links and
%   a small file, at the cost of having to ship the figure folder alongside.
%
%   Not a general Markdown implementation -- deliberately small enough to read
%   in one sitting.
%
%   See also: fileread, regexprep, matlab.net.base64encode

if nargin < 2 || isempty(htmlfile)
    [p,n] = fileparts(mdfile); htmlfile = fullfile(p,[n '.html']);
end
if nargin < 3 || isempty(pageTitle)
    [~,pageTitle] = fileparts(mdfile);
end
if nargin < 4 || isempty(embedImages); embedImages = true; end

txt   = fileread(mdfile);
lines = regexp(txt, '\r\n|\n|\r', 'split');

H = {};                     % output buffer
inCode = false;             % inside ``` fence
listKind = '';              % '', 'ul' or 'ol'
tbl  = {};                  % pending table rows
para = {};                  % pending paragraph lines

for i = 1:numel(lines)
    L = lines{i};

    % ---- fenced code -----------------------------------------------------
    if ~isempty(regexp(strtrim(L), '^```', 'once'))
        if inCode
            H{end+1} = '</code></pre>'; inCode = false; %#ok<*AGROW>
        else
            [para, H] = iFlushPara(para, H);
            [listKind, H] = iCloseList(listKind, H);
            [tbl, H] = iFlushTable(tbl, H);
            H{end+1} = '<pre><code>'; inCode = true;
        end
        continue
    end
    if inCode
        H{end+1} = iEsc(L); continue
    end

    % ---- table rows accumulate; anything else flushes them ---------------
    if ~isempty(regexp(strtrim(L), '^\|.*\|$', 'once'))
        [para, H] = iFlushPara(para, H);
        tbl{end+1} = strtrim(L); continue
    end
    [tbl, H] = iFlushTable(tbl, H);

    % ---- blank line ------------------------------------------------------
    if isempty(strtrim(L))
        [para, H] = iFlushPara(para, H);
        [listKind, H] = iCloseList(listKind, H); continue
    end

    % ---- horizontal rule -------------------------------------------------
    if ~isempty(regexp(strtrim(L), '^(-{3,}|\*{3,}|_{3,})$', 'once'))
        [para, H] = iFlushPara(para, H);
        [listKind, H] = iCloseList(listKind, H);
        H{end+1} = '<hr>'; continue
    end

    % ---- heading ---------------------------------------------------------
    tk = regexp(L, '^(#{1,6})\s+(.*)$', 'tokens', 'once');
    if ~isempty(tk)
        [para, H] = iFlushPara(para, H);
        [listKind, H] = iCloseList(listKind, H);
        n = numel(tk{1});
        H{end+1} = sprintf('<h%d>%s</h%d>', n, iInline(tk{2}), n);
        continue
    end

    % ---- blockquote (consecutive > lines merge) --------------------------
    tk = regexp(L, '^>\s?(.*)$', 'tokens', 'once');
    if ~isempty(tk)
        [para, H] = iFlushPara(para, H);
        [listKind, H] = iCloseList(listKind, H);
        if ~isempty(H) && startsWith(H{end}, '<blockquote>')
            H{end} = regexprep(H{end}, '</blockquote>$', [' ' iInline(tk{1}) '</blockquote>']);
        else
            H{end+1} = sprintf('<blockquote>%s</blockquote>', iInline(tk{1}));
        end
        continue
    end

    % ---- list items ------------------------------------------------------
    tk = regexp(L, '^\s*[-*]\s+(.*)$', 'tokens', 'once');
    if ~isempty(tk)
        [para, H] = iFlushPara(para, H);
        if ~strcmp(listKind,'ul')
            [listKind, H] = iCloseList(listKind, H);
            H{end+1} = '<ul>'; listKind = 'ul';
        end
        H{end+1} = sprintf('<li>%s</li>', iInline(tk{1})); continue
    end
    tk = regexp(L, '^\s*\d+\.\s+(.*)$', 'tokens', 'once');
    if ~isempty(tk)
        [para, H] = iFlushPara(para, H);
        if ~strcmp(listKind,'ol')
            [listKind, H] = iCloseList(listKind, H);
            H{end+1} = '<ol>'; listKind = 'ol';
        end
        H{end+1} = sprintf('<li>%s</li>', iInline(tk{1})); continue
    end

    % ---- wrapped continuation of the current list item -------------------
    % Only when no paragraph is open: an indented line under an open list item
    % belongs to that item, not to a new paragraph.
    if ~isempty(listKind) && isempty(para) && ~isempty(regexp(L, '^\s\s+\S', 'once'))
        H{end} = regexprep(H{end}, '</li>$', [' ' iInline(strtrim(L)) '</li>']);
        continue
    end

    % ---- ordinary prose: accumulate, do NOT emit yet ---------------------
    [listKind, H] = iCloseList(listKind, H);
    para{end+1} = strtrim(L);
end
[~, H] = iFlushPara(para, H);
[~, H] = iCloseList(listKind, H);
[~, H] = iFlushTable(tbl, H);
if inCode; H{end+1} = '</code></pre>'; end

if embedImages
    H = iEmbedImages(H, fileparts(mdfile));
end

fid = fopen(htmlfile, 'w');
fprintf(fid, '%s\n', iHead(pageTitle));
fprintf(fid, '%s\n', H{:});
fprintf(fid, '</article>\n</body>\n</html>\n');
fclose(fid);
d = dir(htmlfile);
fprintf('md2html: %s -> %s (%d blocks, %.0f kB)\n', mdfile, htmlfile, numel(H), d.bytes/1024);
end

% =========================================================================
function H = iEmbedImages(H, baseDir)
% Replace src="relative/path.png" with a base64 data URI, so the page is one
% self-contained file. Remote (http...) and already-inlined (data:) srcs are
% left alone; a missing file is warned about and left as a link rather than
% silently producing a broken page.
MIME = struct('png','image/png','jpg','image/jpeg','jpeg','image/jpeg', ...
              'gif','image/gif','svg','image/svg+xml');
for k = 1:numel(H)
    if isempty(strfind(H{k}, '<img')); continue; end %#ok<STREMP>
    srcs = regexp(H{k}, 'src="([^"]+)"', 'tokens');
    for s = 1:numel(srcs)
        src = srcs{s}{1};
        if startsWith(src,'http') || startsWith(src,'data:'); continue; end
        f = src; if ~isfile(f); f = fullfile(baseDir, src); end
        if ~isfile(f)
            warning('md2html:missingImage','image not found, left as a link: %s', src);
            continue
        end
        [~,~,e] = fileparts(f); e = lower(erase(e,'.'));
        if ~isfield(MIME, e)
            warning('md2html:unknownImageType','not embedding %s', src); continue
        end
        fid = fopen(f,'r'); raw = fread(fid, Inf, '*uint8'); fclose(fid);
        uri = ['data:' MIME.(e) ';base64,' matlab.net.base64encode(raw)];
        H{k} = strrep(H{k}, ['src="' src '"'], ['src="' uri '"']);
    end
end
end

% =========================================================================
function [para, H] = iFlushPara(para, H)
% Join the buffered lines into one paragraph, THEN apply inline markup, so
% emphasis and links that wrap across a line break are handled correctly.
if ~isempty(para)
    H{end+1} = sprintf('<p>%s</p>', iInline(strjoin(para, ' ')));
end
para = {};
end

function [listKind, H] = iCloseList(listKind, H)
if ~isempty(listKind); H{end+1} = ['</' listKind '>']; end
listKind = '';
end

function [tbl, H] = iFlushTable(tbl, H)
if isempty(tbl); return; end
rows = {};
for k = 1:numel(tbl)
    r = regexprep(tbl{k}, '^\||\|$', '');       % strip outer pipes
    rows{end+1} = strsplit(r, '|'); %#ok<AGROW>
end
% row 2 of a Markdown table is the ---|--- separator
hasHdr = numel(rows) >= 2 && all(cellfun(@(c) ~isempty(regexp(strtrim(c), ...
                              '^:?-{2,}:?$', 'once')), rows{2}));
H{end+1} = '<table>';
for k = 1:numel(rows)
    if hasHdr && k == 2; continue; end
    cellTag = 'td';
    if hasHdr && k == 1; cellTag = 'th'; H{end+1} = '<thead>'; end
    H{end+1} = '<tr>';
    for c = 1:numel(rows{k})
        H{end+1} = sprintf('<%s>%s</%s>', cellTag, iInline(strtrim(rows{k}{c})), cellTag);
    end
    H{end+1} = '</tr>';
    if hasHdr && k == 1; H{end+1} = '</thead><tbody>'; end
end
if hasHdr; H{end+1} = '</tbody>'; end
H{end+1} = '</table>';
tbl = {};
end

function s = iEsc(s)
s = strrep(s, '&', '&amp;');
s = strrep(s, '<', '&lt;');
s = strrep(s, '>', '&gt;');
end

function s = iInline(s)
s = iEsc(s);
% `code` first, so markup inside it is not re-interpreted
s = regexprep(s, '`([^`]+)`',                    '<code>$1</code>');
s = regexprep(s, '!\[([^\]]*)\]\(([^)]+)\)',     '<img alt="$1" src="$2">');
s = regexprep(s, '\[([^\]]+)\]\(([^)]+)\)',      '<a href="$2">$1</a>');
% lazy, and NOT [^*], so that nested emphasis (**bold *italic* bold**) works
s = regexprep(s, '\*\*(.+?)\*\*',                '<strong>$1</strong>');
s = regexprep(s, '(?<![\*\w])\*([^*]+)\*(?!\*)', '<em>$1</em>');
end

function h = iHead(ttl)
h = sprintf([ ...
'<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n' ...
'<meta name="viewport" content="width=device-width, initial-scale=1">\n' ...
'<title>%s</title>\n<style>\n' ...
' :root{--fg:#1a1a1a;--mut:#5b6570;--bg:#fff;--card:#f7f8fa;--bd:#d8dde3;--acc:#0b5cad}\n' ...
' @media (prefers-color-scheme:dark){:root{--fg:#e8eaed;--mut:#9aa4af;--bg:#15181c;' ...
'--card:#1e2228;--bd:#333a42;--acc:#68b0ff}}\n' ...
' *{box-sizing:border-box}\n' ...
' body{margin:0;background:var(--bg);color:var(--fg);\n' ...
'  font:16px/1.62 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}\n' ...
' article{max-width:1180px;margin:0 auto;padding:44px 26px 90px}\n' ...
' h1{font-size:2.0rem;line-height:1.22;margin:.2em 0 .5em;letter-spacing:-.02em}\n' ...
' h2{font-size:1.42rem;margin:2.1em 0 .6em;padding-bottom:.28em;border-bottom:2px solid var(--bd)}\n' ...
' h3{font-size:1.13rem;margin:1.7em 0 .5em;color:var(--acc)}\n' ...
' h4{font-size:1rem;margin:1.4em 0 .4em;color:var(--mut);text-transform:uppercase;letter-spacing:.05em}\n' ...
' p{margin:.75em 0} li{margin:.3em 0}\n' ...
' a{color:var(--acc)}\n' ...
' code{background:var(--card);border:1px solid var(--bd);border-radius:4px;\n' ...
'  padding:.08em .34em;font:.87em/1.4 ui-monospace,SFMono-Regular,Consolas,monospace}\n' ...
' pre{background:var(--card);border:1px solid var(--bd);border-radius:8px;\n' ...
'  padding:14px 16px;overflow-x:auto}\n' ...
' pre code{background:none;border:0;padding:0;font-size:.86rem}\n' ...
' blockquote{margin:1.1em 0;padding:.7em 1.1em;border-left:4px solid var(--acc);\n' ...
'  background:var(--card);border-radius:0 6px 6px 0;color:var(--fg)}\n' ...
' hr{border:0;border-top:1px solid var(--bd);margin:2.4em 0}\n' ...
' table{border-collapse:collapse;width:100%%;margin:1.1em 0;font-size:.93rem;display:block;overflow-x:auto}\n' ...
' th,td{border:1px solid var(--bd);padding:7px 11px;text-align:left;vertical-align:top}\n' ...
' th{background:var(--card);font-weight:650;white-space:nowrap}\n' ...
' tbody tr:nth-child(even){background:color-mix(in srgb,var(--card) 55%%,transparent)}\n' ...
' img{max-width:100%%;height:auto;display:block;margin:1.3em auto;\n' ...
'  border:1px solid var(--bd);border-radius:8px;background:#fff}\n' ...
'</style>\n</head>\n<body>\n<article>'], ttl);
end
