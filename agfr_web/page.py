"""The web page: one file, no build step, no network access at run time.

Layout and conventions (token in the URL, plain JavaScript, bundled 3Dmol.js,
light/dark switch) follow ADCP_docking/adcp_dock/page.py by Aadhil Haq
(BSD 3-Clause). The panel itself mirrors AGFRgui's PDB analyser:
one tab per biomolecule, a residue tree with receptor and ligand columns and a
params column, the rec./lig./other show buttons and the two save buttons.
"""
from . import __version__

PAGE = r"""<!doctype html><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>AGFR &middot; PDB panel</title>
<script>try { const t = localStorage.getItem("agfrweb.theme");
  if (t === "light" || t === "dark") document.documentElement.dataset.theme = t;
} catch (e) {}</script>
<link rel="stylesheet" href="/static/app.css">
<header>
  <b>AGFR</b> <span class=dim>PDB panel</span>
  <label class=btn id=openbtn title="Read a .pdb, .ent, .pdb.gz or .ent.gz file">
    Open PDB&hellip;<input type=file id=file accept=".pdb,.ent,.gz" hidden></label>
  <span id=pdblabel class=dim>drop a PDB file anywhere on the page</span>
  <span class=grow></span>
  <span id=envnote class=dim></span>
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
<footer class=dim>agfr_web __VERSION__ &middot; computation by ADFR (AGFRgui's own code,
ADCP distribution) &middot; viewer 3Dmol.js (Rego &amp; Koes, Bioinformatics 2015)</footer>
<script src="/static/3Dmol-min.js"></script>
<script src="/static/app.js"></script>
"""


def render():
    return PAGE.replace('__VERSION__', __version__)
