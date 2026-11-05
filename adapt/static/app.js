/* Author: Michel Sanner, PhD. Copyright (c) 2026. Portions adapted from ADCP_docking (BSD 3-Clause); see THIRD_PARTY_NOTICES.md. */
/* ADapt page behaviour. Plain JavaScript, no framework and no build step
   (the ADCP_docking convention). Every change is made by the server, which
   runs AGFRgui's own code; the page only draws what the server returns, so a
   checkbox never shows a state the ADFR atoms do not have. */
'use strict';
const T = new URLSearchParams(location.search).get('t') || '';
const $ = i => document.getElementById(i);
const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g,
  c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));

let S = null;            // session info from /api/upload
let TAB = null;          // current biomolecule name
let TREE = null;         // tree json of the current tab
let NODE = {};           // id -> node json
let EL = {};             // id -> <li>
let STATES = {};         // id -> [rec, lig, clickable]
let OPEN = {};           // id -> true/false once the user toggled it
let ATOMS = '';          // per-atom state letters, atom-index order
let VIEW = null, MODEL = null, ROT = null;
let STRUCT = {};         // tab -> {pdb, n}
let SEL = null;          // selected node id
let HL = null;           // {indices, scheme}: what is highlighted, kept across repaints

/* ---------------------------------------------------------------- api */
let busyTimer = null;
function busy(on){
  clearTimeout(busyTimer);
  if(on) busyTimer = setTimeout(() => { $('busy').hidden = false; }, 250);
  else $('busy').hidden = true;
}
async function api(path, body, raw, quiet){
  const url = path + (path.includes('?') ? '&' : '?') + 't=' + encodeURIComponent(T);
  const opt = raw ? {method:'POST', body: raw}
            : body ? {method:'POST', headers:{'Content-Type':'application/json'},
                      body: JSON.stringify(body)} : {};
  // quiet: no "working" overlay (it would cover the page and swallow the
  // next wheel event while rotamers are being browsed)
  if(!quiet) busy(true);
  try {
    const r = await fetch(url, opt);
    let d;
    try { d = await r.json(); }
    catch(e){ d = {error: 'the server did not answer with JSON (' + r.status + ')'}; }
    if(d.error){ status('<span class=bad>' + esc(d.error) + '</span>');
                 if(d.trace) log([d.trace]); }
    return d;
  } catch(e){
    status('<span class=bad>the server did not answer: ' + esc(e) + '</span>');
    return {error: String(e)};
  } finally { if(!quiet) busy(false); }
}
function status(html){ $('status').innerHTML = html || '&nbsp;'; }
function log(lines){
  if(!lines || !lines.length) return;
  const pre = $('log');
  pre.textContent += lines.join('\n') + '\n';
  pre.scrollTop = pre.scrollHeight;
}

/* ------------------------------------------------------------ loading */
async function openFile(file){
  if(!file) return;
  status('reading ' + esc(file.name) + '&hellip;');
  const d = await api('api/upload?name=' + encodeURIComponent(file.name), null, file);
  if(!d.error) loaded(d);
}

/* Fetch from RCSB: the server downloads <ID>.pdb (cached) and opens it */
$('fetchform').addEventListener('submit', async e => {
  e.preventDefault();
  const code = $('pdbcode').value.trim().toUpperCase();
  if(!/^[0-9][A-Z0-9]{3}$/.test(code)){
    status('<span class=bad>a PDB id is four characters, the first a digit (e.g. 6GGD)</span>');
    return;
  }
  status('fetching ' + esc(code) + ' from the RCSB&hellip;');
  const d = await api('api/fetch?code=' + encodeURIComponent(code), {});
  if(!d.error) loaded(d);
});

function loaded(d){
  S = d; STRUCT = {}; OPEN = {}; SEL = null;
  $('pdblabel').textContent = d.label;
  $('folder').textContent = 'session folder: ' + (d.folder_display || d.folder);
  $('saves').innerHTML = '';
  for(const sv of (d.saves || [])) listSave(sv);
  showTab(d);
  status('loaded ' + esc(d.pdbid) + ': ' + d.tabs.length + ' biomolecule tab'
         + (d.tabs.length === 1 ? '' : 's')
         + (d.fetched ? ' &middot; <span class=dim>' + esc(d.fetched.url)
                        + (d.fetched.from_cache ? ' (cached)' : '') + '</span>' : ''));
}
$('file').addEventListener('change', e => { openFile(e.target.files[0]); e.target.value = ''; });

/* A page cannot choose the folder a plain <input type=file> opens in. Where
   the browser has the File System Access API (Chrome, Edge), the picker is
   opened with an `id`, and the browser then reopens it in the folder the last
   file with that id came from, across reloads and restarts. The memory is kept
   per origin (host AND port), so it holds as long as the server keeps its
   port. Elsewhere the <input> inside the button is used as before. */
if('showOpenFilePicker' in window){
  $('openbtn').addEventListener('click', async e => {
    e.preventDefault();                  // do not also open the <input>
    let handle;
    try {
      [handle] = await window.showOpenFilePicker({
        id: 'agfr-pdb',                  // kept across the rename: a new id would forget the folder
        types: [{description: 'PDB and mmCIF files (.pdb .ent .cif, gzipped or not)',
                 accept: {'chemical/x-pdb': ['.pdb', '.ent', '.gz'],
                          'chemical/x-mmcif': ['.cif', '.mmcif']}}]});
    } catch(err){
      if(err.name !== 'AbortError')      // AbortError = the dialog was cancelled
        status('<span class=bad>could not open the file picker: ' + esc(err.message) + '</span>');
      return;
    }
    openFile(await handle.getFile());
  });
}
document.addEventListener('dragover', e => { e.preventDefault(); document.body.classList.add('drag'); });
document.addEventListener('dragleave', e => { if(!e.relatedTarget) document.body.classList.remove('drag'); });
document.addEventListener('drop', e => {
  e.preventDefault(); document.body.classList.remove('drag');
  if(e.dataTransfer.files.length) openFile(e.dataTransfer.files[0]);
});

function renderTabs(){
  const nav = $('tabs');
  nav.innerHTML = '';
  for(const name of S.tabs){
    const b = document.createElement('button');
    b.textContent = name;
    if(name === TAB) b.className = 'on';
    b.onclick = async () => {
      if(name === TAB) return;
      const d = await api('api/tab', {s: S.session, tab: name});
      if(!d.error) showTab(d);
    };
    nav.appendChild(b);
  }
}

/* the address bar names the session and tab shown, so a reload (or a link)
   reopens them while the server runs */
function remember(){
  try {
    const P = new URLSearchParams(location.search);
    P.set('session', S.session); P.set('tab', TAB); P.delete('pdb');
    history.replaceState(null, '', '?' + P.toString());
  } catch(e) {}
}
function showTab(d){
  TAB = d.tab;
  remember();
  renderTabs();
  TREE = d.tree;
  renderTree();
  applyView(d);
  loadStructure(TAB);
}

/* ---------------------------------------------------------------- tree */
function renderTree(){
  const box = $('tree');
  const top = box.scrollTop;
  NODE = {}; EL = {};
  $('treetitle').textContent = TREE.header;
  const ul = document.createElement('ul');
  for(const r of TREE.roots) ul.appendChild(makeNode(r, 0));
  box.innerHTML = '';
  box.appendChild(ul);
  box.scrollTop = top;
  paintStates();
}

function isOpen(n){ return (n.id in OPEN) ? OPEN[n.id] : !!n.open; }

function makeNode(n, depth){
  NODE[n.id] = n;
  const li = document.createElement('li');
  li.className = 'kind-' + n.kind + (n.menu ? ' hasmenu' : '');
  li.dataset.id = n.id;
  EL[n.id] = li;
  const kids = n.children || [];
  li.innerHTML =
    '<div class=row><div class=name style="padding-left:' + (depth * 14) + 'px">'
    + '<span class=tw>' + (kids.length ? (isOpen(n) ? '&#9662;' : '&#9656;') : '') + '</span>'
    + '<input type=checkbox class=cb data-col=rec>'
    + '<span class=lab title="' + esc(n.label) + '">' + esc(n.label) + '</span></div>'
    + '<div class=params></div>'
    + '<div class=lig><input type=checkbox class=cb data-col=lig></div></div>';
  if(n.id === SEL) li.firstChild.classList.add('sel');
  params(li.querySelector('.params'), n);
  if(kids.length){
    if(isOpen(n)) li.appendChild(childList(n, depth));
    else li.classList.add('closed');
  }
  return li;
}

function childList(n, depth){
  const ul = document.createElement('ul');
  for(const c of n.children) ul.appendChild(makeNode(c, depth + 1));
  return ul;
}

function toggle(id){
  const n = NODE[id], li = EL[id];
  const open = !isOpen(n);
  OPEN[id] = open;
  li.querySelector('.tw').innerHTML = open ? '&#9662;' : '&#9656;';
  if(open){
    if(!li.querySelector(':scope>ul')){
      const depth = parseInt(li.querySelector('.name').style.paddingLeft) / 14;
      li.appendChild(childList(n, depth));
      paintStates();
    }
    li.classList.remove('closed');
  } else li.classList.add('closed');
}

function paintStates(){
  for(const id in EL){
    const st = STATES[id];
    if(!st) continue;
    const boxes = EL[id].firstChild.querySelectorAll('input.cb');
    for(let c = 0; c < 2; c++){
      const cb = boxes[c], s = st[c];
      cb.style.visibility = s ? 'visible' : 'hidden';
      cb.checked = s === 'checked';
      cb.indeterminate = s === 'partial';
      cb.disabled = !st[2];
    }
  }
}

/* the params column: altloc radios, manual charge, rotamer preview */
function params(cell, n){
  if(n.altloc){
    n.altloc.alts.forEach((a, i) => {
      const lab = document.createElement('label');
      lab.innerHTML = '<input type=radio name=alt' + n.id + (i === n.altloc.current ? ' checked' : '')
        + '> ' + esc(a) + ' ' + n.altloc.occ[i].toFixed(2);
      lab.firstChild.onchange = async () => {
        const d = await api('api/altloc', {s: S.session, node: n.id, index: i});
        if(!d.error){ n.altloc.current = i; applyView(d); focus(n.id); }
        else renderTree();
      };
      cell.appendChild(lab);
    });
  }
  if(n.charge != null){
    const inp = document.createElement('input');
    inp.type = 'number'; inp.min = -10; inp.max = 10; inp.step = 0.05;
    inp.value = n.charge.toFixed(2);
    inp.title = 'charge written for this ion (AGFR default shown)';
    inp.onchange = async () => {
      const d = await api('api/charge', {s: S.session, node: n.id, value: parseFloat(inp.value)});
      if(!d.error){ n.charge = d.value; inp.value = d.value.toFixed(2); applyView(d);
                    status('charge of ' + esc(n.label) + ' set to ' + d.value.toFixed(2)); }
      else inp.value = n.charge.toFixed(2);
    };
    cell.appendChild(inp);
  }
  if(n.rotamer){
    const sel = document.createElement('select');
    sel.title = 'The rotamer to build when the receptor is written (auto: the '
      + 'best-scoring one). Scroll the mouse wheel over this box to step through them.';
    let h = '<option value=-1>auto ' + (n.rotamer.auto + 1) + '</option>';
    for(let i = 0; i < n.rotamer.n; i++) h += '<option value=' + i + '>' + (i + 1) + '</option>';
    sel.innerHTML = h;
    sel.value = n.rotamer.current;
    sel.onchange = () => { n.rotamer.current = parseInt(sel.value); rotamer(n.id); };
    // The mouse wheel over the box steps through the rotamers, one per wheel
    // notch (a notch is ~100 px of deltaY); a trackpad's small deltas are summed
    // into notch-sized steps. The tree does not scroll while the pointer is here.
    let acc = 0;
    sel.addEventListener('wheel', e => {
      e.preventDefault();
      acc += e.deltaMode === 0 ? e.deltaY : e.deltaY * 40;
      if(Math.abs(acc) < 50) return;
      const k = sel.selectedIndex + Math.sign(acc);
      acc = 0;
      if(k < 0 || k >= sel.options.length) return;
      sel.selectedIndex = k;
      sel.onchange();
      if(SEL !== n.id) select(n.id);     // show the residue the rotamers belong to
    }, {passive: false});
    cell.appendChild(sel);
  }
}

/* events on the tree, delegated */
$('tree').addEventListener('click', async e => {
  const li = e.target.closest('li');
  if(!li) return;
  const id = +li.dataset.id;
  if(e.target.classList.contains('tw')){ toggle(id); return; }
  if(e.target.matches('input.cb')){
    e.preventDefault();          // the server decides the new state
    const col = e.target.dataset.col;
    const d = await api('api/click', {s: S.session, node: id, col: col});
    if(!d.error){ applyView(d); }
    return;
  }
  if(e.target.closest('.params')) return;
  if(e.target.closest('.name')) select(id);
});
$('tree').addEventListener('dblclick', e => {
  const li = e.target.closest('li');
  if(li && e.target.closest('.name') && !e.target.matches('input')) toggle(+li.dataset.id);
});
$('tree').addEventListener('contextmenu', e => {
  const li = e.target.closest('li');
  if(!li) return;
  const n = NODE[+li.dataset.id];
  if(!n.menu) return;
  e.preventDefault();
  select(n.id);
  openMenu(n, e.clientX, e.clientY);
});

function select(id){
  if(SEL != null && EL[SEL]) EL[SEL].firstChild.classList.remove('sel');
  SEL = id;
  if(EL[id]) EL[id].firstChild.classList.add('sel');
  focus(id);
  if(NODE[id].rotamer) rotamer(id);
}

/* ---------------------------------------------------------------- menu */
function closeMenu(){ $('menu').hidden = true; }
document.addEventListener('click', e => { if(!e.target.closest('#menu')) closeMenu(); });
document.addEventListener('keydown', e => { if(e.key === 'Escape') closeMenu(); });

async function openMenu(n, x, y){
  const m = $('menu');
  let items = [];
  if(n.menu === 'sc'){
    const complete = n.kind === 'category' || n.label.indexOf(' build:') >= 0;
    items = [['Complete sidechain', complete, () => sidechain(n, 'complete')],
             ['Include as is', !complete, () => sidechain(n, 'asis')]];
    m.innerHTML = '<div class=head>' + esc(n.label.split(' ')[0]) + '</div>';
  } else if(n.menu === 'substitute'){
    const done = n.label.indexOf(' -> ') >= 0;
    m.innerHTML = '<div class=head>' + (done ? esc(n.label.split(' ')[0])
      : 'replace ' + esc(n.label.split(' ')[0]) + ' by') + '</div>';
    items = done ? [['Undo the replacement', false, () => substitute(n, null)]]
      : ['ZN', 'MG', 'CA', 'FE', 'MN'].map(el => [el, false, () => substitute(n, el)]);
  } else if(n.menu === 'mutate'){
    const d = await api('api/mutateopts?s=' + S.session + '&node=' + n.id);
    if(d.error) return;
    m.innerHTML = '<div class=head>mutate ' + esc(n.label.split(' ')[0]) + ' to</div>';
    items = d.names.map(name => {
      const target = name.replace('Base Res:', '').trim();
      return [name, d.current && d.current === target, () => mutate(n, name)];
    });
  }
  for(const [label, cur, fn] of items){
    const div = document.createElement('div');
    div.textContent = label;
    if(cur) div.className = 'cur';
    div.onclick = () => { closeMenu(); fn(); };
    m.appendChild(div);
  }
  m.hidden = false;
  const r = m.getBoundingClientRect();
  m.style.left = Math.min(x, innerWidth - r.width - 8) + 'px';
  m.style.top = Math.min(y, innerHeight - r.height - 8) + 'px';
}

async function sidechain(n, act){
  const d = await api('api/sidechain', {s: S.session, node: n.id, act: act});
  if(!d.error){ TREE = d.tree; renderTree(); applyView(d); }
}
async function substitute(n, element){
  const d = await api('api/substitute', {s: S.session, node: n.id, element: element});
  if(!d.error){ TREE = d.tree; renderTree(); applyView(d); }
}
async function mutate(n, choice){
  const d = await api('api/mutate', {s: S.session, node: n.id, choice: choice});
  if(!d.error){ TREE = d.tree; renderTree(); applyView(d);
                status(esc(NODE[n.id] ? NODE[n.id].label : n.label)); }
}

/* ------------------------------------------------------------- buttons */
function applyView(d){
  if(d.log) log(d.log);
  if(d.states){ STATES = d.states; paintStates(); }
  if(d.buttons){
    $('showrec').disabled = $('saverec').disabled = !d.buttons.receptor;
    $('showlig').disabled = $('savelig').disabled = !d.buttons.ligand;
    $('showother').disabled = !d.buttons.other;
  }
  if(d.atoms != null){ ATOMS = d.atoms; restyle(); }
}

async function show(what, scheme){
  const d = await api('api/select', {s: S.session, what: what});
  if(d.error) return;
  HL = null; restyle();
  highlight(d.indices, scheme, true);
  status(d.indices.length + ' atoms: ' + what);
}
$('showrec').onclick = () => show('receptor', 'greenCarbon');
$('showlig').onclick = () => show('ligand', 'greenCarbon');
$('showother').onclick = () => show('other', 'orangeCarbon');

/* one line per save under "Files written", a download link per file. URLs are
   relative to the page, so ADapt works under a path prefix (e.g. /adapt/). */
function downloadUrl(save, file){
  return 'api/download?s=' + encodeURIComponent(S.session) + '&save=' + save
    + (file ? '&file=' + encodeURIComponent(file) : '') + '&t=' + encodeURIComponent(T);
}
function listSave(sv){
  const li = document.createElement('li');
  li.innerHTML = '#' + sv.save + ' ' + esc(sv.kind) + ' (' + esc(sv.tab) + '): '
    + sv.files.map(f => esc(f.name) + ' <a href="' + downloadUrl(sv.save, f.name) + '">download</a>')
              .join(' &middot; ')
    + (sv.files.length > 1 && sv.kind === 'ligand'
       ? ' &middot; <a href="' + downloadUrl(sv.save) + '">all (zip)</a>' : '');
  $('saves').appendChild(li);
}
async function save(kind){
  const d = await api('api/save', {s: S.session, kind: kind});
  if(d.error) return;
  applyView(d);
  listSave({save: d.save, kind: kind, tab: TAB, files: d.files});
  // the receptor downloads its PDBQT (the protonated PDB stays a link);
  // ligands download their one file, or all fragments as a zip
  const url = kind === 'receptor' ? downloadUrl(d.save, d.files[0].name) : downloadUrl(d.save);
  const a = document.createElement('a');
  a.href = url; a.download = '';
  document.body.appendChild(a); a.click(); a.remove();
  status('wrote ' + d.files.map(f => esc(f.name)).join(', '));
}
$('saverec').onclick = () => save('receptor');
$('savelig').onclick = () => save('ligand');

/* -------------------------------------------------------------- viewer */
function dark(){
  const t = document.documentElement.dataset.theme;
  return t ? t === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
}
function viewer(){
  if(!VIEW) VIEW = $3Dmol.createViewer($('viewer'), {backgroundColor: dark() ? '#12151c' : 'white'});
  return VIEW;
}

async function loadStructure(tab){
  if(!STRUCT[tab]){
    const d = await api('api/structure?s=' + S.session + '&tab=' + encodeURIComponent(tab));
    if(d.error) return;
    STRUCT[tab] = d;
  }
  if(tab !== TAB) return;
  const v = viewer();
  v.clear(); ROT = null;
  const d = STRUCT[tab];
  HL = null;                     // its indices belonged to the previous molecule
  MODEL = v.addModel(d.pdb, 'pdb', {altLoc: '*', keepH: true});
  const n = MODEL.selectedAtoms({}).length;
  if(n !== d.n) status('<span class=bad>viewer read ' + n + ' atoms of ' + d.n
                       + '; highlighting may be off</span>');
  drawGaps(d.gaps || []);
  restyle();
  v.zoomTo();
  v.render();
}

/* missing residues, drawn as AGFRgui draws them (gridGUI.addGeomsForGaps): a
   dotted orange line across each gap with a sphere at its middle labelled with
   the missing residues, and a sphere beyond each missing terminus. The server
   computes the coordinates (adapt/gaps.py:geometry). */
const GAP_COLOR = '#ffa600';
function drawGaps(gaps){
  const v = viewer();
  for(const g of gaps){
    const [p, q] = g.line.map(c => ({x: c[0], y: c[1], z: c[2]}));
    v.addCylinder({start: p, end: q, radius: 0.07, color: GAP_COLOR, dashed: true,
                   dashLength: 0.35, gapLength: 0.25, fromCap: 1, toCap: 1});
    const c = {x: g.sphere[0], y: g.sphere[1], z: g.sphere[2]};
    v.addSphere({center: c, radius: 0.3, color: GAP_COLOR});
    v.addLabel((g.kind === 'gap' ? '' : g.kind + ' term: ') + g.label,
               {position: c, fontSize: 11, fontColor: '#8a5300', backgroundColor: '#fff3d6',
                backgroundOpacity: 0.85, borderThickness: 0, inFront: true,
                alignment: 'bottomLeft'});
  }
}

/* receptor grey carbons, ligand orange sticks, the rest purple; symmetry
   copies (lower-case letters) at half opacity, as gridGUI does */
function restyle(){
  if(!MODEL || !ATOMS) return;
  const v = viewer();
  const st = a => ATOMS[a.index] || 'o';
  v.setStyle({model: MODEL}, {});
  const groups = [
    ['R', {line: {colorscheme: 'grayCarbon'}, cartoon: {color: '#9fb8e0'}}],
    ['r', {line: {colorscheme: 'grayCarbon', opacity: 0.5}, cartoon: {color: '#9fb8e0', opacity: 0.5}}],
    ['L', {stick: {colorscheme: 'orangeCarbon', radius: 0.18}, cartoon: {color: '#f0b060'}}],
    ['l', {stick: {colorscheme: 'orangeCarbon', radius: 0.18, opacity: 0.5}, cartoon: {color: '#f0b060', opacity: 0.5}}],
    ['O', {line: {colorscheme: 'purpleCarbon'}, cartoon: {color: '#d4d4d4'}}],
    ['o', {line: {colorscheme: 'purpleCarbon', opacity: 0.5}, cartoon: {color: '#d4d4d4', opacity: 0.5}}]];
  for(const [c, style] of groups)
    v.setStyle({model: MODEL, predicate: a => st(a) === c}, style);
  // atoms without bonds (ions, waters) are invisible as lines
  v.addStyle({model: MODEL, predicate: a => !a.bonds.length},
             {sphere: {radius: 0.22, colorscheme: 'default'}});
  // a repaint follows every flag change (a checkbox, an altloc, a charge, a
  // mutation); the atoms the user is looking at stay highlighted, as in AGFRgui
  if(HL) paintHighlight(HL);
  v.render();
}

// a lone atom (an ion, a metal) has no carbon for a "greenCarbon" scheme to
// colour, so it is drawn as a larger sphere in the highlight colour itself
const HL_COLOR = {greenCarbon: '#2fd13a', orangeCarbon: '#ff9f1a'};
function paintHighlight(h){
  const set = new Set(h.indices);
  const sel = {model: MODEL, predicate: a => set.has(a.index)};
  const color = HL_COLOR[h.scheme] || '#2fd13a';
  const lone = h.indices.length === 1;
  // atoms without bonds (ions, waters) show nothing as sticks: they get a
  // sphere of the highlight colour instead; "set", not "add", because their own
  // sphere style would keep the element colour
  viewer().setStyle({model: MODEL, predicate: a => set.has(a.index) && !a.bonds.length},
                    {sphere: {radius: lone ? 0.75 : 0.45, color: color}});
  // residues scattered over a zoomed-out molecule need thicker sticks to be seen
  if(h.wide === undefined) h.wide = extent(h.indices) >= 8;
  viewer().addStyle({model: MODEL, predicate: a => set.has(a.index) && a.bonds.length > 0},
                    {stick: {colorscheme: h.scheme, radius: h.wide ? 0.4 : 0.22}});
  return sel;
}

/* How far the selected atoms spread from their centre, in Angstrom. */
function extent(indices){
  const set = new Set(indices);
  const at = MODEL.selectedAtoms({predicate: a => set.has(a.index)});
  if(!at.length) return 0;
  let cx = 0, cy = 0, cz = 0;
  for(const a of at){ cx += a.x; cy += a.y; cz += a.z; }
  cx /= at.length; cy /= at.length; cz /= at.length;
  let r = 0;
  for(const a of at) r = Math.max(r, Math.hypot(a.x - cx, a.y - cy, a.z - cz));
  return r;
}

function highlight(indices, scheme, zoom){
  if(!MODEL) return;
  const v = viewer();
  HL = {indices: indices, scheme: scheme};
  const sel = paintHighlight(HL);
  if(zoom && indices.length){
    v.zoomTo(sel);
    // a residue or an ion fills the view once zoomed to: step back so its
    // surroundings stay visible. A selection spread over the molecule (a whole
    // category, all the waters) is framed as it is.
    if(!HL.wide) v.zoom(0.35);
  }
  v.render();
}

async function focus(id){
  const d = await api('api/focus?s=' + S.session + '&node=' + id, null, null, true);   // view only: quiet
  if(d.error) return;
  HL = null; restyle();
  highlight(d.indices, 'greenCarbon', true);
}

/* One preview request at a time per residue: while one is out, further
   changes only mark it stale, and the newest choice is fetched when it returns.
   A fast wheel thus costs at most two requests, and ends on the right rotamer. */
const ROTQ = {};
async function rotamer(id){
  const q = ROTQ[id] || (ROTQ[id] = {busy: false, again: false});
  if(q.busy){ q.again = true; return; }
  q.busy = true;
  try {
    do { q.again = false; await rotamerOnce(id); } while(q.again);
  } finally { q.busy = false; }
}

async function rotamerOnce(id){
  const n = NODE[id];
  if(!n) return;
  const d = await api('api/rotamer', {s: S.session, node: id, choice: n.rotamer.current},
                      null, true);
  if(d.error) return;
  const v = viewer();
  if(ROT) v.removeModel(ROT);
  ROT = v.addModel(d.pdb, 'pdb', {keepH: true});
  v.setStyle({model: ROT}, {stick: {colorscheme: 'cyanCarbon', radius: 0.2}});
  v.render();
  status(esc(n.label.split(' ')[0]) + ' rotamer ' + (d.index + 1)
         + (n.rotamer.current < 0 ? ' (auto)' : '') + ' &middot; ' + esc(d.message)
         + ' <span class=dim>(built when the receptor is written)</span>');
}

/* --------------------------------------------------------------- theme */
$('theme').onclick = () => {
  const next = dark() ? 'light' : 'dark';
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem('adapt.theme', next); } catch(e) {}   // 'agfrweb.theme' before the rename is still read
  if(VIEW){ VIEW.setBackgroundColor(dark() ? '#12151c' : 'white'); VIEW.render(); }
};

/* Opened with ?session=<id> (a session made elsewhere, e.g. by a tool calling
   /api/upload), or ?pdb=<ID> (fetched from the RCSB), and optionally
   &tab=<name> (ASU, Biomol 1, ...): show it straight away. */
(async () => {
  const P = new URLSearchParams(location.search);
  const sid = P.get('session'), code = (P.get('pdb') || '').trim(), tab = P.get('tab');
  let d = null;
  if(sid){
    status('opening session ' + esc(sid) + '&hellip;');
    d = await api('api/session?s=' + encodeURIComponent(sid));
  } else if(code){
    $('pdbcode').value = code.toUpperCase();
    status('fetching ' + esc(code.toUpperCase()) + ' from the RCSB&hellip;');
    d = await api('api/fetch?code=' + encodeURIComponent(code.toUpperCase()), {});
  }
  if(!d || d.error) return;
  loaded(d);
  if(tab && tab !== d.tab && (d.tabs || []).includes(tab)){
    const t = await api('api/tab', {s: S.session, tab: tab});
    if(!t.error) showTab(t);
  }
})();

api('api/status').then(d => {
  if(d.error) return;
  const env = d.environment || {};
  $('envnote').innerHTML = (d.missing && d.missing.length)
    ? '<span class=bad>missing: ' + d.missing.map(esc).join('; ') + '</span>'
    : esc((env.adfr_dist || 'ADFR') + ' · python ' + (env.python || ''));
});
