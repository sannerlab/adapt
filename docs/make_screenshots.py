#!/usr/bin/env python3
# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Regenerate the documentation screenshots (docs/img/*.png) from ADapt itself.

    python3 docs/make_screenshots.py [scene ...]

Starts an ADapt server (bin/adapt, port 0) and a headless Chrome, and for
each scene fetches an entry from the RCSB, clicks through the page the way a
user would, and captures the window or a part of it. Every image in the docs
is an ADapt screenshot made here, so they can be remade after the page changes.
Sessions go to ~/adapt_runs/docs_screenshots/, the Chrome profile and the
server log to tests/_out/docs_run/ (git-ignored); never /tmp.

Needs: google-chrome, network access to files.rcsb.org (entries are cached
under ~/adapt_runs/docs_screenshots/_rcsb/ after the first run).
"""
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, 'tools'))
import cdp  # noqa: E402

IMG = os.path.join(HERE, 'img')
WORK = os.path.join(REPO, 'tests', '_out', 'docs_run')          # Chrome profile, server log
# the sessions the scenes create: under the default output folder, so the
# folder the screenshots show reads as a user's would (~/adapt_runs/...)
SESSIONS = os.path.join(os.path.expanduser('~'), 'adapt_runs', 'docs_screenshots')
W, H = 1440, 960

# Counts in-flight requests, so a scene waits for the page to finish talking
# to the server instead of sleeping for a guessed time.
INFLIGHT = """(() => { const f = window.fetch; window.__inflight = 0;
  window.fetch = function(...a){ window.__inflight++;
    return f.apply(this, a).finally(() => { window.__inflight--; }); }; })();"""

# Tree navigation helpers, evaluated in the page.
HELPERS = r"""
window.DOC = {
  lab(li){ return li.querySelector(':scope>.row .lab').textContent; },
  chain(c){ return [...document.querySelectorAll('#tree > ul > li.kind-chain')].find(li => this.lab(li) === c); },
  open(li){ if(li.classList.contains('closed')) li.querySelector(':scope>.row .tw').click(); return li; },
  close(li){ if(!li.classList.contains('closed') && li.querySelector(':scope>ul')) li.querySelector(':scope>.row .tw').click(); return li; },
  child(li, prefix){ this.open(li); return [...li.querySelectorAll(':scope>ul>li')].find(x => this.lab(x).startsWith(prefix)) || null; },
  path(chain, ...prefixes){ let li = this.chain(chain); for(const p of prefixes){ if(!li) return null; li = this.child(li, p); } return li; },
  collapseAll(){ for(const li of [...document.querySelectorAll('#tree li')].reverse()) this.close(li); },
  select(li){ li.querySelector(':scope>.row .lab').click(); },
  context(li){ const l = li.querySelector(':scope>.row .lab'); const r = l.getBoundingClientRect();
    l.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, clientX: r.x + Math.min(r.width, 90), clientY: r.y + r.height - 2})); },
  menu(label){ const d = [...document.querySelectorAll('#menu div')].find(d => d.textContent === label); d.click(); },
  check(li, col){ li.querySelector(':scope>.row input[data-col=' + col + ']').click(); },
  show(li){ li.scrollIntoView({block: 'center'}); },
  tab(name){ [...document.querySelectorAll('#tabs button')].find(b => b.textContent === name).click(); },
  rect(el){ const r = el.getBoundingClientRect(); return {x: r.x, y: r.y, width: r.width, height: r.height}; },
};
"""


class Page(object):
    def __init__(self, url):
        os.makedirs(WORK, exist_ok=True)
        self.chrome, self.ws = cdp.start_chrome(os.path.join(WORK, 'chrome'), width=W, height=H)
        ws = self.ws
        ws.call('Page.enable')
        ws.call('Runtime.enable')
        ws.call('Emulation.setDeviceMetricsOverride', width=W, height=H, deviceScaleFactor=1, mobile=False)
        ws.call('Page.addScriptToEvaluateOnNewDocument', source=INFLIGHT)
        ws.call('Page.navigate', url=url)
        cdp.wait(ws, "document.readyState === 'complete' && !!window.$3Dmol", what='page load')
        self.js(HELPERS)
        self.idle()

    def js(self, expr):
        return cdp.js(self.ws, expr)

    def idle(self, settle=0.6):
        """No request in flight for a short while, then let WebGL draw."""
        quiet = 0
        t = time.time()
        while quiet < 3:
            if time.time() - t > 240:
                raise RuntimeError('the page never went idle')
            quiet = quiet + 1 if self.js('window.__inflight === 0') else 0
            time.sleep(0.15)
        time.sleep(settle)

    def fetch(self, code, tab=None):
        self.js("document.getElementById('pdblabel').textContent = ''")
        self.js("document.getElementById('pdbcode').value = %r; document.getElementById('fetchbtn').click()" % code)
        cdp.wait(self.ws, "document.getElementById('pdblabel').textContent.length > 0", 300, 'fetch ' + code)
        self.idle()
        if tab:
            self.tab(tab)

    def tab(self, name):
        self.js('DOC.tab(%r)' % name)
        self.idle()
        cdp.wait(self.ws, 'typeof STRUCT !== "undefined" && !!STRUCT[TAB] && TAB === %r' % name, 120, 'tab ' + name)
        self.idle()

    def node(self, chain, *prefixes):
        """The tree row chain > prefix > ... (opening the way there)."""
        expr = 'DOC.path(%r%s)' % (chain, ''.join(', %r' % p for p in prefixes))
        if not self.js('!!' + expr):
            raise RuntimeError('no tree row %s' % ' > '.join((chain,) + prefixes))
        self.idle(0.2)
        return expr

    def shot(self, name, part=None, pad=0):
        """Capture the window, or the left panel / 3D view / an element."""
        clip = None
        if part == 'left':
            clip = self.js("DOC.rect(document.getElementById('left'))")
        elif part == 'right':
            clip = self.js("DOC.rect(document.getElementById('right'))")
        elif part == 'center':          # the middle of the 3D view, where a focused residue sits
            v = self.js("DOC.rect(document.getElementById('viewer'))")
            cw, ch = 520, 400
            clip = {'x': v['x'] + (v['width'] - cw) / 2, 'y': v['y'] + (v['height'] - ch) / 2,
                    'width': cw, 'height': ch}
        elif part:
            clip = self.js('DOC.rect(%s)' % part)
        if clip:
            clip = {'x': max(0, clip['x'] - pad), 'y': max(0, clip['y'] - pad),
                    'width': clip['width'] + 2 * pad, 'height': clip['height'] + 2 * pad}
        path = os.path.join(IMG, name + '.png')
        cdp.screenshot(self.ws, path, clip)
        print('  %-28s %4d kB' % (name + '.png', os.path.getsize(path) // 1024), flush=True)

    def badges(self, marks):
        """Numbered callouts: [(n, css selector, fx, fy, dx, dy)] puts badge n at
        fraction (fx, fy) of the element's box, shifted by (dx, dy) pixels."""
        self.js(r"""(() => { const box = document.createElement('div'); box.id = 'docbadges';
          box.style.cssText = 'position:fixed;inset:0;pointer-events:none;z-index:99';
          for(const [n, sel, fx, fy, dx, dy] of %s){
            const el = document.querySelector(sel); if(!el) throw new Error('no element ' + sel);
            const r = el.getBoundingClientRect();
            const b = document.createElement('div'); b.textContent = n;
            b.style.cssText = 'position:absolute;width:26px;height:26px;border-radius:13px;background:#d6336c;'
              + 'color:#fff;font:bold 14px system-ui;display:flex;align-items:center;justify-content:center;'
              + 'box-shadow:0 1px 4px rgba(0,0,0,.4);border:2px solid #fff';
            b.style.left = (r.left + fx * r.width + dx - 13) + 'px';
            b.style.top = (r.top + fy * r.height + dy - 13) + 'px';
            box.appendChild(b); }
          document.body.appendChild(box); })()""" % json.dumps(marks))

    def unbadge(self):
        self.js("const b = document.getElementById('docbadges'); if(b) b.remove()")

    def close(self):
        errs = cdp.js_errors(self.ws)
        self.chrome.terminate()
        return errs


# ---------------------------------------------------------------- scenes
def overview(p):
    p.fetch('1FV1', 'ASU')
    p.js("DOC.collapseAll(); DOC.open(DOC.chain('A')); DOC.open(DOC.chain('B'))")
    p.idle()
    p.badges([(1, '#openbtn', 0.5, 1, 0, 16), (2, '#fetchform', 0.5, 1, 0, 16),
              (3, '#pdblabel', 0.5, 1, 0, 16), (4, '#tabs button:last-child', 1, 0.5, 18, 0),
              (5, '#tree li.kind-chain .cb', 0, 0.5, -60, 22), (6, '#treehead span:nth-child(2)', 0, 0.5, -18, 0),
              (7, '#treehead span:nth-child(3)', 1, 0.5, 20, 0), (8, '#showother', 1, 0.5, 18, 0),
              (9, '#savelig', 0, 0, -8, -16), (10, '#viewer', 0, 0, 24, 24),
              (11, '#status', 1, 0.5, -24, 0), (12, '#savesbox summary', 0, 0.5, 130, 0)])
    p.shot('overview')
    p.unbadge()


def link(p):
    """The tree drives the 3D view: a residue, a whole category, the waters."""
    p.fetch('1FV1', 'Biomol 1')
    p.js("DOC.collapseAll()")
    res = p.node('A', 'protein', "TRP'")
    p.js('DOC.show(%s); DOC.select(%s)' % (res, res))
    p.idle(1.0)
    p.shot('link_residue')
    p.js("DOC.collapseAll()")
    cat = p.node('B', 'sidechains to complete')
    p.js('DOC.select(%s)' % cat)
    p.idle(1.0)
    p.shot('link_category')
    p.js("DOC.collapseAll()")
    wat = p.node('A', 'water')
    p.js('DOC.select(%s)' % wat)
    p.idle(1.0)
    p.shot('link_water')


def biomolecule(p):
    p.fetch('1HHP', 'ASU')
    p.shot('1hhp_asu')
    p.tab('Biomol 1')
    p.shot('1hhp_biomol1')


def assemblies(p):
    p.fetch('1FV1', 'Biomol 1')
    p.js("DOC.collapseAll(); DOC.child(DOC.chain('A'), 'sidechains to complete'); "
         "DOC.child(DOC.chain('B'), 'sidechains to complete'); DOC.open(DOC.chain('C'))")
    p.idle()
    p.shot('1fv1_biomol1')
    p.tab('Biomol 2')
    p.js("DOC.collapseAll(); DOC.child(DOC.chain('D'), 'sidechains to complete'); "
         "DOC.child(DOC.chain('E'), 'sidechains to complete'); DOC.open(DOC.chain('F'))")
    p.idle()
    p.shot('1fv1_biomol2')


def sidechains(p):
    p.fetch('1FV1', 'Biomol 2')
    p.js("DOC.collapseAll()")
    arg = p.node('D', 'sidechains to complete', "ARG'100")
    p.js('DOC.select(%s)' % arg)
    p.idle(1.0)
    p.shot('sc_rotamer_auto')
    p.shot('sc_rotamer_auto_zoom', 'center')
    p.js("(() => { const s = %s.querySelector('.params select'); s.value = '2'; "
         "s.dispatchEvent(new Event('change')); })()" % arg)
    p.idle(1.0)
    p.shot('sc_rotamer_3_zoom', 'center')
    glu = p.node('D', 'sidechains to complete', "GLU'158")
    p.js('DOC.select(%s)' % glu)
    p.idle(1.0)
    p.js('DOC.context(%s)' % glu)
    p.idle(0.5)
    p.shot('sc_menu')
    p.js("DOC.menu('Include as is')")
    p.idle()
    p.shot('sc_asis', "DOC.chain('D')", pad=4)


def altlocs(p):
    p.fetch('1FV1', 'Biomol 2')
    p.js("DOC.collapseAll()")
    grp = p.node('D', 'altloc')
    p.js('DOC.select(%s)' % grp)
    p.idle(1.0)
    p.shot('altloc_group')
    ser = p.node('D', 'altloc', 'protein', "SER'19")
    # SER19's locations tie at 0.50, so ADFR's default varies between runs:
    # choose them by letter, B first, so the images always match the captions
    pick = ("(() => { const r = [...%s.querySelectorAll('.params label')]"
            ".find(l => l.textContent.trim().startsWith('%%s')).querySelector('input');"
            " if(!r.checked) r.click(); })()" % ser)
    p.js(pick % 'B')
    p.idle()
    p.js('DOC.show(%s); DOC.select(%s)' % (ser, ser))
    p.idle(1.0)
    p.shot('altloc_residue')
    p.shot('altloc_B_zoom', 'center')
    p.js(pick % 'A')
    p.idle(1.0)
    p.shot('altloc_A_zoom', 'center')
    p.shot('altloc_radios', ser, pad=3)


def ligand(p):
    p.fetch('6GGD', 'Biomol 1')
    p.js("DOC.collapseAll()")
    eyb = p.node('A', 'ligands', "EYB'402")
    p.js("DOC.check(%s, 'lig')" % eyb)
    p.idle()
    p.js('DOC.select(%s)' % eyb)
    p.idle(1.0)
    p.shot('ligand_checked')


def charges(p):
    p.fetch('6GGD', 'Biomol 1')
    p.js("DOC.collapseAll()")
    zn = p.node('A', 'manual charges', "ZN'401")
    p.js("DOC.check(%s, 'rec')" % zn)
    p.idle()
    p.js('DOC.select(%s)' % zn)
    p.idle(1.0)
    p.shot('charge_zn')


def unknown(p):
    p.fetch('1PLC')
    p.js("DOC.collapseAll()")
    cu = p.node('A', 'unknown atom types', 'other', "CU'")
    p.js('DOC.select(%s)' % cu)
    p.idle(1.0)
    p.shot('unknown_cu')


def gaps(p):
    p.fetch('4A1X', 'ASU')
    p.js("DOC.collapseAll()")
    for c in 'ABCD':
        p.node(c, 'missing residues')
    p.idle()
    p.shot('gaps_asu')
    p.tab('Biomol 1')
    p.js("DOC.collapseAll()")
    g = p.node('A', 'missing residues', 'gap:')
    p.node('C', 'missing residues')
    p.js('DOC.select(%s)' % g)
    p.idle(1.0)
    p.shot('gaps_biomol1')


def rebuild(p):
    """4A1X: chain A's gap GKGG rebuilt (Biomol 1); chain B's gap GK, whose two
    residues cannot close it, refused with its reason (Biomol 2)."""
    p.fetch('4A1X', 'Biomol 1')
    p.js("DOC.collapseAll()")
    g = p.node('A', 'missing residues', 'gap:')
    p.node('C', 'missing residues')
    p.js('DOC.select(%s)' % g)
    p.idle(1.0)
    p.js("%s.querySelector(':scope>.row .params input').click()" % g)
    p.idle(1.0)
    p.shot('rebuild_4a1x')
    p.tab('Biomol 2')
    p.js("DOC.collapseAll()")
    g = p.node('B', 'missing residues', 'gap:')
    p.js("%s.querySelector(':scope>.row .params input').click()" % g)
    p.idle()
    p.js('DOC.show(%s)' % p.node('B', 'missing residues', 'gap:'))
    p.idle(0.3)
    p.shot('rebuild_refused', 'left')


def mutations(p):
    p.fetch('6GGD', 'Biomol 2')
    p.js("DOC.collapseAll()")
    ser = p.node('B', 'protein', "SER'96")
    p.js('DOC.show(%s); DOC.select(%s)' % (ser, ser))
    p.idle(1.0)
    p.js('DOC.context(%s)' % ser)
    p.idle(0.8)
    p.shot('mutate_menu')
    p.js("DOC.menu('ILE')")
    p.idle()
    ser = p.node('B', 'protein', "SER'96")
    p.js('DOC.show(%s)' % ser)
    p.shot('mutate_done', 'left')


def saving(p):
    """The 6GGD session of the README: Biomol 2, ZN and EYB kept, SER96 -> ILE,
    then both files written."""
    p.fetch('6GGD', 'Biomol 2')
    p.js("DOC.collapseAll()")
    p.js("DOC.check(%s, 'rec')" % p.node('B', 'manual charges', "ZN'401"))
    p.idle()
    p.js("DOC.check(%s, 'lig')" % p.node('B', 'ligands', "EYB'402"))
    p.idle()
    p.js("DOC.context(%s)" % p.node('B', 'protein', "SER'96"))
    p.idle(0.5)
    p.js("DOC.menu('ILE')")
    p.idle()
    p.js("document.getElementById('saverec').click()")
    p.idle()
    p.js("document.getElementById('savelig').click()")
    p.idle()
    p.js("document.getElementById('logbox').open = true")
    p.idle(0.3)
    p.shot('saved_files')
    folder = p.js("document.getElementById('folder').textContent").replace('session folder: ', '')
    return folder


SCENES = [overview, link, biomolecule, assemblies, sidechains, altlocs, ligand, charges,
          unknown, gaps, rebuild, mutations, saving]


def main(argv):
    want = set(argv)
    os.makedirs(IMG, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    log = open(os.path.join(WORK, 'server.log'), 'w')
    os.makedirs(SESSIONS, exist_ok=True)
    server = subprocess.Popen([os.path.join(REPO, 'bin', 'adapt'), 'serve', '--port', '0',
                               '--outdir', SESSIONS], stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, universal_newlines=True)
    url = None
    for line in server.stdout:
        log.write(line)
        m = re.search(r'(http://127\.0\.0\.1:\d+/\?t=\S+)', line)
        if m:
            url = m.group(1)
            break
    if not url:
        server.terminate()
        raise SystemExit('the ADapt server did not print its URL')
    try:
        for scene in SCENES:
            if want and scene.__name__ not in want:
                continue
            print(scene.__name__, flush=True)
            p = Page(url)
            try:
                out = scene(p)
                if out:
                    print('  session folder:', out)
            finally:
                errs = p.close()
            if errs:
                raise SystemExit('JavaScript errors in %s: %s' % (scene.__name__, errs[:3]))
    finally:
        server.terminate()
        log.close()


if __name__ == '__main__':
    main(sys.argv[1:])
