# Author: Michel Sanner, PhD. Copyright (c) 2026. Portions adapted from ADCP_docking (BSD 3-Clause); see THIRD_PARTY_NOTICES.md.
"""The web page: one file, no build step, no network access at run time.

Layout and conventions (token in the URL, plain JavaScript, bundled 3Dmol.js,
light/dark switch) follow ADCP_docking/adcp_dock/page.py by Aadhil Haq
(BSD 3-Clause). The panel itself mirrors AGFRgui's PDB analyser:
one tab per biomolecule, a residue tree with receptor and ligand columns and a
params column, the rec./lig./other show buttons and the two save buttons.
"""
from . import APP, APP_LONG, COPYRIGHT, __version__

PAGE = r"""<!doctype html><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__APP__ &middot; __APP_LONG__</title>
<script>try { const t = localStorage.getItem("adapt.theme") || localStorage.getItem("agfrweb.theme");
  if (t === "light" || t === "dark") document.documentElement.dataset.theme = t;
} catch (e) {}</script>
<link rel="stylesheet" href="/static/app.css?v=__V_CSS__">
<header>
  <b>__APP__</b> <span class=dim>__APP_LONG__</span>
  <label class=btn id=openbtn title="Read a .pdb, .ent, .pdb.gz or .ent.gz file">
    Open PDB&hellip;<input type=file id=file accept=".pdb,.ent,.gz" hidden></label>
  <form id=fetchform class=fetch title="Download an entry from the RCSB (files.rcsb.org) by its PDB id">
    <input id=pdbcode placeholder="PDB id" maxlength=4 size=6 spellcheck=false autocomplete=off>
    <button id=fetchbtn>Fetch from RCSB</button></form>
  <span id=pdblabel class=dim>or drop a PDB file anywhere on the page</span>
  <span class=grow></span>
  <span id=envnote class=dim></span>
  <a class=btn id=docslink href="/docs/index.html" target=_blank rel=noopener
     title="What __APP__ does, with worked examples">Docs</a>
  <button id=theme class=plain title="light / dark">&#9681;</button>
</header>
<main>
  <section id=left>
    <p class=intro>Select biological unit, fragments to include in receptor,
    alternate locations. Mutate/complete residue sidechains.
    <span class=dim>Right-click a residue for its menu.</span></p>
    <nav id=tabs></nav>
    <div id=treehead class=row><span id=treetitle>receptor</span><span>params</span><span>lig.</span></div>
    <div id=tree><p class=dim id=empty>No structure loaded.</p></div>
    <div class=bar>
      <span class=dim>show</span>
      <button id=showrec disabled title="Show atoms included in the receptor">rec.</button>
      <button id=showlig disabled title="Show atoms included in ligand(s)">lig.</button>
      <button id=showother disabled title="Show atoms not included in receptor or ligand">other</button>
      <span class=grow></span>
      <button id=saverec class=primary disabled title="Write the receptor PDBQT (reduce hydrogens, Gasteiger charges, AD types)">receptor PDBQT</button>
      <button id=savelig class=primary disabled title="Write one ligand PDBQT per connected fragment">ligand PDBQT</button>
    </div>
  </section>
  <section id=right>
    <div id=viewer></div>
    <div id=status class=dim>&nbsp;</div>
    <details id=savesbox open><summary>Files written</summary><ul id=saves></ul>
      <p class=dim id=folder></p></details>
    <details id=logbox><summary>Log</summary><pre id=log></pre></details>
  </section>
</main>
<div id=menu class=menu hidden></div>
<div id=busy hidden><span>working&hellip;</span></div>
<footer class=dim>__APP__ __VERSION__ &middot; __COPYRIGHT__ &middot; computation by ADFR (AGFRgui's own code,
ADCP distribution) &middot; viewer 3Dmol.js (Rego &amp; Koes, Bioinformatics 2015)</footer>
<script src="/static/3Dmol-min.js?v=__V_3DMOL__"></script>
<script src="/static/app.js?v=__V_JS__"></script>
"""


def _v(name):
    """Short md5 of a static file: its URL changes whenever its content does,
    so the hour-long cache of /static/ never serves an outdated script."""
    import hashlib
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', name)
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()[:10]


def render():
    import html
    return (PAGE.replace('__VERSION__', __version__)
                .replace('__APP_LONG__', html.escape(APP_LONG))
                .replace('__COPYRIGHT__', html.escape(COPYRIGHT).replace('(c)', '&copy;'))
                .replace('__APP__', html.escape(APP))
                .replace('__V_CSS__', _v('app.css'))
                .replace('__V_JS__', _v('app.js'))
                .replace('__V_3DMOL__', _v('3Dmol-min.js')))
