function docsToPdf()
%DOCSTOPDF  Render the repository's Markdown documents to printable files.
%
%   docsToPdf()          run from anywhere; writes docs/pdf/<name>.html and,
%                        when a Chromium-based browser is found (Edge on
%                        Windows, Chrome, Chromium), docs/pdf/<name>.pdf.
%
%   Markdown -> HTML is done here (headings, paragraphs, bullet and
%   numbered lists, tables, fenced code, inline code, bold, links), so the
%   tables keep their layout. The PDF step calls the browser headless:
%       msedge --headless --print-to-pdf=<pdf> <html>
%   If no browser is found the HTML files are still written: open one in
%   any browser and print it (Ctrl+P, "Save as PDF"), the print CSS is in
%   the file. ALL_DOCS combines every document in one landscape file.
%
%   To add a document, append a row to DOCS below: {source, name, landscape}.

repo   = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repo, 'docs', 'pdf');
if ~exist(outDir, 'dir'), mkdir(outDir); end

DOCS = { ...
    'README.md',               'README',          false; ...
    'docs/HOST_DECISIONS.md',  'HOST_DECISIONS',  true;  ...
    'docs/REVIEW_FINDINGS.md', 'REVIEW_FINDINGS', true;  ...
    'docs/CONCEPTS.md',        'CONCEPTS',        false};
COMBINED = 'ALL_DOCS';

browser = findBrowser();
if isempty(browser)
    fprintf('no Chromium-based browser found: HTML only (open and print it, or set the path in findBrowser)\n');
end

bodies = cell(size(DOCS, 1), 1);
for k = 1:size(DOCS, 1)
    src = DOCS{k, 1};
    bodies{k} = mdToHtml(fileread(fullfile(repo, src)));
    render(browser, outDir, DOCS{k, 2}, DOCS{k, 3}, bodies{k}, src);
end
allBody = '';
for k = 1:numel(bodies)
    allBody = [allBody '<div class="doc">' bodies{k} '</div>']; %#ok<AGROW>
end
render(browser, outDir, COMBINED, true, allBody, strjoin(DOCS(:, 1)', ', '));
end

% ------------------------------------------------------------------------
function render(browser, outDir, name, landscape, body, srcLabel)
if landscape, orient = 'landscape'; else, orient = 'portrait'; end
html = [ ...
    '<!DOCTYPE html><html><head><meta charset="utf-8"><title>' name '</title>' ...
    '<style>' pageCss(orient) '</style></head><body>' body ...
    '<div class="footer">' escapeHtml(srcLabel) ' (Spoof-Detection repository)</div></body></html>'];
htmlPath = fullfile(outDir, [name '.html']);
pdfPath  = fullfile(outDir, [name '.pdf']);
fid = fopen(htmlPath, 'w', 'n', 'UTF-8'); fwrite(fid, html, 'char'); fclose(fid);
fprintf('wrote %s\n', htmlPath);
if ~isempty(browser)
    cmd = sprintf('"%s" --headless --disable-gpu --no-sandbox --no-pdf-header-footer --print-to-pdf="%s" "%s"', ...
        browser, pdfPath, fileUrl(htmlPath));
    [status, out] = system(cmd);
    if status == 0 && exist(pdfPath, 'file')
        fprintf('wrote %s\n', pdfPath);
    else
        fprintf('PDF step failed for %s (%s); the HTML file can be printed from a browser\n', name, strtrim(out));
    end
end
end

function url = fileUrl(p)
p = strrep(p, '\', '/');
if ~strncmp(p, '/', 1), p = ['/' p]; end     % Windows drive letter: file:///C:/...
url = ['file://' p];
end

function css = pageCss(orient)
css = [ ...
    '@page { size: A4 ' orient '; margin: 14mm 12mm 16mm 12mm; }' ...
    'body { font-family: Helvetica, Arial, sans-serif; font-size: 9.5pt; line-height: 1.35; color: #111; }' ...
    'h1 { font-size: 17pt; margin: 0 0 6pt 0; border-bottom: 1.5px solid #333; padding-bottom: 3pt; }' ...
    '.doc { page-break-before: always; } .doc:first-child { page-break-before: auto; }' ...
    'h2 { font-size: 13pt; margin: 14pt 0 5pt 0; page-break-after: avoid; }' ...
    'h3 { font-size: 11pt; margin: 10pt 0 4pt 0; page-break-after: avoid; }' ...
    'p { margin: 4pt 0; } ul, ol { margin: 3pt 0 3pt 16pt; padding: 0; } li { margin: 1.5pt 0; }' ...
    'code { font-family: Menlo, Consolas, "DejaVu Sans Mono", monospace; font-size: 8.5pt; background: #f2f2f2; padding: 0 2px; border-radius: 2px; }' ...
    'pre { background: #f5f5f5; border: 1px solid #ddd; padding: 6pt; font-size: 8pt; white-space: pre-wrap; page-break-inside: avoid; }' ...
    'pre code { background: none; padding: 0; }' ...
    'table { border-collapse: collapse; width: 100%; margin: 6pt 0 8pt 0; font-size: 8.3pt; }' ...
    'thead { display: table-header-group; } tr { page-break-inside: avoid; }' ...
    'th, td { border: 1px solid #999; padding: 3pt 4pt; vertical-align: top; text-align: left; }' ...
    'th { background: #e8e8e8; font-weight: bold; } a { color: #114488; text-decoration: none; }' ...
    '.footer { font-size: 7.5pt; color: #666; margin-top: 10pt; border-top: 1px solid #ccc; padding-top: 3pt; }'];
end

% ------------------------------------------------------------------------
function browser = findBrowser()
browser = '';
cand = { ...
    getenv('CHROME'), ...
    'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe', ...
    'C:\Program Files\Microsoft\Edge\Application\msedge.exe', ...
    'C:\Program Files\Google\Chrome\Application\chrome.exe', ...
    'C:\Program Files (x86)\Google\Chrome\Application\chrome.exe', ...
    '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', ...
    '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge', ...
    '/usr/bin/chromium', '/usr/bin/chromium-browser', '/usr/bin/google-chrome'};
d = dir('/opt/pw-browsers/chromium-*/chrome-linux/chrome');          % Playwright install (Linux)
for k = 1:numel(d), cand{end+1} = fullfile(d(k).folder, d(k).name); end %#ok<AGROW>
for k = 1:numel(cand)
    if ~isempty(cand{k}) && exist(cand{k}, 'file') == 2
        browser = cand{k};
        break;
    end
end
end

% ------------------------------------------------------------------------
function html = mdToHtml(text)
% Line-based Markdown subset: headings, paragraphs, - / * bullets, 1. lists
% (indented continuation lines join the item), | tables |, ``` fenced code,
% and inline `code`, **bold**, [text](url).
text  = strrep(text, sprintf('\r\n'), sprintf('\n'));
lines = strsplit(text, sprintf('\n'));
out   = {};
para  = {};      % paragraph lines being collected
i = 1; n = numel(lines);
while i <= n
    ln = lines{i};
    if strncmp(strtrim(ln), '```', 3)                        % fenced code
        out = flushPara(out, para); para = {};
        code = {}; i = i + 1;
        while i <= n && ~strncmp(strtrim(lines{i}), '```', 3)
            code{end+1} = escapeHtml(lines{i}); i = i + 1; %#ok<AGROW>
        end
        out{end+1} = ['<pre><code>' strjoin(code, sprintf('\n')) '</code></pre>']; %#ok<AGROW>
        i = i + 1;
    elseif ~isempty(regexp(ln, '^#{1,6}\s', 'once'))         % heading
        out = flushPara(out, para); para = {};
        lvl = numel(regexp(ln, '^#+', 'match', 'once'));
        txt = regexprep(ln, '^#+\s*', '');
        out{end+1} = sprintf('<h%d>%s</h%d>', lvl, inline(txt), lvl); %#ok<AGROW>
        i = i + 1;
    elseif strncmp(strtrim(ln), '|', 1)                      % table
        out = flushPara(out, para); para = {};
        rows = {};
        while i <= n && strncmp(strtrim(lines{i}), '|', 1)
            rows{end+1} = lines{i}; i = i + 1; %#ok<AGROW>
        end
        out{end+1} = tableHtml(rows); %#ok<AGROW>
    elseif ~isempty(regexp(ln, '^\s*([-*]|\d+\.)\s+', 'once'))   % list
        out = flushPara(out, para); para = {};
        if isempty(regexp(ln, '^\s*\d+\.', 'once')), tag = 'ul'; else, tag = 'ol'; end
        items = {};
        while i <= n
            cur = lines{i};
            if ~isempty(regexp(cur, '^\s*([-*]|\d+\.)\s+', 'once'))
                items{end+1} = regexprep(cur, '^\s*([-*]|\d+\.)\s+', ''); %#ok<AGROW>
            elseif ~isempty(strtrim(cur)) && ~isempty(regexp(cur, '^\s+', 'once')) && ~isempty(items)
                items{end} = [items{end} ' ' strtrim(cur)];    % continuation line
            else
                break;
            end
            i = i + 1;
        end
        li = '';
        for k = 1:numel(items), li = [li '<li>' inline(items{k}) '</li>']; end %#ok<AGROW>
        out{end+1} = ['<' tag '>' li '</' tag '>']; %#ok<AGROW>
    elseif isempty(strtrim(ln))                              % blank: paragraph ends
        out = flushPara(out, para); para = {};
        i = i + 1;
    else
        para{end+1} = strtrim(ln); %#ok<AGROW>
        i = i + 1;
    end
end
out  = flushPara(out, para);
html = strjoin(out, sprintf('\n'));
end

function out = flushPara(out, para)
if ~isempty(para)
    out{end+1} = ['<p>' inline(strjoin(para, ' ')) '</p>'];
end
end

function h = tableHtml(rows)
cells = cell(numel(rows), 1);
for r = 1:numel(rows)
    s = strtrim(rows{r});
    s = regexprep(s, '^\|', ''); s = regexprep(s, '\|$', '');
    cells{r} = strtrim(strsplit(s, '|'));
end
keep = true(numel(rows), 1);
for r = 1:numel(rows)                                        % drop |---|---| separator rows
    if all(cellfun(@(c) ~isempty(regexp(c, '^:?-+:?$', 'once')), cells{r})), keep(r) = false; end
end
cells = cells(keep);
h = '<table>';
for r = 1:numel(cells)
    if r == 1, tag = 'th'; else, tag = 'td'; end
    row = '';
    for c = 1:numel(cells{r}), row = [row '<' tag '>' inline(cells{r}{c}) '</' tag '>']; end %#ok<AGROW>
    if r == 1, h = [h '<thead><tr>' row '</tr></thead><tbody>']; else, h = [h '<tr>' row '</tr>']; end %#ok<AGROW>
end
h = [h '</tbody></table>'];
end

function s = inline(s)
s = escapeHtml(s);
s = regexprep(s, '`([^`]+)`', '<code>$1</code>');
s = regexprep(s, '\*\*(.+?)\*\*', '<strong>$1</strong>');
s = regexprep(s, '\[([^\]]+)\]\(([^)]+)\)', '<a href="$2">$1</a>');
end

function s = escapeHtml(s)
s = strrep(s, '&', '&amp;');
s = strrep(s, '<', '&lt;');
s = strrep(s, '>', '&gt;');
end
