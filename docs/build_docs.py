#!/usr/bin/env python3
# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Build the ADapt documentation: docs/src/<page>.html -> docs/<page>.html.

    python3 docs/build_docs.py

Each source file holds one page's content (the inside of <main>), starting
with <h1>. This script wraps it in the shared layout: navigation, previous /
next links and footer. The pages are static HTML; open docs/index.html in a
browser, or follow the "Docs" link of a running ADapt.
"""
import html
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'src')

# (file, navigation label, section)
PAGES = [
    ('index', 'Overview', 'Start here'),
    ('start', 'Getting started', 'Start here'),
    ('tree-3d', 'The tree and the 3D view', 'Start here'),
    ('case-biomolecule', 'Biological assembly, not the ASU', 'Cases'),
    ('case-assemblies', 'Choosing between assemblies', 'Cases'),
    ('case-sidechains', 'Missing side chains', 'Cases'),
    ('case-altlocs', 'Alternate locations', 'Cases'),
    ('case-ligand', 'Ligands and other hetero groups', 'Cases'),
    ('case-metals', 'Metals, ions, unsupported atoms', 'Cases'),
    ('case-missing-residues', 'Missing residues', 'Cases'),
    ('case-mutations', 'Mutating residues', 'Cases'),
    ('output', 'Writing the PDBQT files', 'Results'),
    ('differences', 'ADapt and AGFRgui', 'Reference'),
]

TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title} &middot; ADapt documentation</title>
<link rel="stylesheet" href="style.css">
</head><body><div class="layout">
<nav class="side">
<a class="brand" href="index.html">ADapt<small>AutoDock preparation &amp; perception tool</small></a>
{nav}
</nav>
<main>
{content}
<div class="pn">{prev}{next}</div>
<footer class="page">ADapt documentation. Author: Michel Sanner, PhD. Copyright &copy; 2026.
Every screenshot is of ADapt itself, made by <code>docs/make_screenshots.py</code>; the chemistry
is ADFR's own code (ADCP distribution).</footer>
</main></div></body></html>
"""


def nav(current):
    out, section = [], None
    for f, label, sec in PAGES:
        if sec != section:
            if section is not None:
                out.append('</ul>')
            out.append('<h4>%s</h4><ul>' % html.escape(sec))
            section = sec
        cls = ' class="here"' if f == current else ''
        out.append('<li><a%s href="%s.html">%s</a></li>' % (cls, f, html.escape(label)))
    out.append('</ul>')
    return '\n'.join(out)


def main():
    for i, (f, label, sec) in enumerate(PAGES):
        with open(os.path.join(SRC, f + '.html')) as fh:
            content = fh.read().strip()
        m = re.search(r'<h1>(.*?)</h1>', content, re.S)
        title = re.sub('<[^>]+>', '', m.group(1)) if m else label
        prev = ('<a href="%s.html">&larr; %s</a>' % (PAGES[i - 1][0], html.escape(PAGES[i - 1][1]))
                if i else '<span></span>')
        nxt = ('<a href="%s.html">%s &rarr;</a>' % (PAGES[i + 1][0], html.escape(PAGES[i + 1][1]))
               if i + 1 < len(PAGES) else '<span></span>')
        page = TEMPLATE.format(title=title, nav=nav(f), content=content, prev=prev, next=nxt)
        with open(os.path.join(HERE, f + '.html'), 'w') as fh:
            fh.write(page)
        # every image a page names must exist
        for img in re.findall(r'src="(img/[^"]+)"', content):
            if not os.path.isfile(os.path.join(HERE, img)):
                raise SystemExit('%s.html: missing %s' % (f, img))
    print('built %d pages in %s' % (len(PAGES), HERE))


if __name__ == '__main__':
    main()
