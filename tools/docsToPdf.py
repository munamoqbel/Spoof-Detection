#!/usr/bin/env python3
"""Render the repository's Markdown documents to printable PDF.

Markdown -> HTML (python-markdown, tables + fenced code) -> PDF through
headless Chromium, so the tables keep their layout. Wide documents are
printed in landscape.

    python3 tools/docsToPdf.py            # writes docs/pdf/*.pdf

Requires: pip install markdown; a Chromium binary (CHROME env var, or the
Playwright install under /opt/pw-browsers, or chromium on PATH).
"""
import glob
import os
import shutil
import subprocess
import sys

import markdown

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "pdf")

# (source, output name, landscape)
DOCS = [
    ("README.md",               "README",          False),
    ("docs/HOST_DECISIONS.md",  "HOST_DECISIONS",  True),
    ("docs/REVIEW_FINDINGS.md", "REVIEW_FINDINGS", True),
    ("docs/CONCEPTS.md",         "CONCEPTS",        False),
]
COMBINED = "ALL_DOCS"      # the three above in one landscape file, one document per page break

CSS = """
@page { size: A4 %(orient)s; margin: 14mm 12mm 16mm 12mm; }
body { font-family: Helvetica, Arial, sans-serif; font-size: 9.5pt; line-height: 1.35; color: #111; }
h1 { font-size: 17pt; margin: 0 0 6pt 0; border-bottom: 1.5px solid #333; padding-bottom: 3pt; }
.doc { page-break-before: always; }
.doc:first-child { page-break-before: auto; }
h2 { font-size: 13pt; margin: 14pt 0 5pt 0; page-break-after: avoid; }
h3 { font-size: 11pt; margin: 10pt 0 4pt 0; page-break-after: avoid; }
p  { margin: 4pt 0; }
ul, ol { margin: 3pt 0 3pt 16pt; padding: 0; }
li { margin: 1.5pt 0; }
code { font-family: Menlo, Consolas, "DejaVu Sans Mono", monospace; font-size: 8.5pt;
       background: #f2f2f2; padding: 0 2px; border-radius: 2px; }
pre { background: #f5f5f5; border: 1px solid #ddd; padding: 6pt; font-size: 8pt;
      white-space: pre-wrap; page-break-inside: avoid; }
pre code { background: none; padding: 0; }
table { border-collapse: collapse; width: 100%%; margin: 6pt 0 8pt 0; font-size: 8.3pt;
        page-break-inside: auto; }
thead { display: table-header-group; }
tr { page-break-inside: avoid; }
th, td { border: 1px solid #999; padding: 3pt 4pt; vertical-align: top; text-align: left; }
th { background: #e8e8e8; font-weight: bold; }
a { color: #114488; text-decoration: none; }
.footer { font-size: 7.5pt; color: #666; margin-top: 10pt; border-top: 1px solid #ccc; padding-top: 3pt; }
"""

HTML = """<!DOCTYPE html><html><head><meta charset="utf-8"><title>%(title)s</title>
<style>%(css)s</style></head><body>%(body)s
<div class="footer">%(source)s (Spoof-Detection repository)</div></body></html>"""


def findChrome():
    cand = [os.environ.get("CHROME", "")]
    cand += glob.glob("/opt/pw-browsers/chromium-*/chrome-linux/chrome")
    cand += glob.glob("/opt/pw-browsers/chromium_headless_shell-*/chrome-linux/headless_shell")
    for name in ("chromium", "chromium-browser", "google-chrome", "chrome"):
        p = shutil.which(name)
        if p:
            cand.append(p)
    for c in cand:
        if c and os.path.exists(c):
            return c
    sys.exit("no Chromium binary found: set CHROME=/path/to/chrome")


def toHtml(src):
    with open(os.path.join(ROOT, src), encoding="utf-8") as f:
        text = f.read()
    return markdown.markdown(text, extensions=["tables", "fenced_code", "sane_lists"])


def render(chrome, src, name, landscape):
    if isinstance(src, list):
        body = "".join('<div class="doc">%s</div>' % toHtml(s) for s in src)
        src = ", ".join(src)
    else:
        body = toHtml(src)
    html = HTML % {
        "title": name,
        "css": CSS % {"orient": "landscape" if landscape else "portrait"},
        "body": body,
        "source": src,
    }
    htmlPath = os.path.join(OUT, name + ".html")
    pdfPath = os.path.join(OUT, name + ".pdf")
    with open(htmlPath, "w", encoding="utf-8") as f:
        f.write(html)
    cmd = [chrome, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
           "--print-to-pdf=" + pdfPath, "file://" + htmlPath]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.remove(htmlPath)
    print("wrote", os.path.relpath(pdfPath, ROOT))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    chrome = findChrome()
    for src, name, landscape in DOCS:
        render(chrome, src, name, landscape)
    render(chrome, [d[0] for d in DOCS], COMBINED, True)
