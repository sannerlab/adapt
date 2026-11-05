# ADapt — AutoDock preparation & perception tool

Author: Michel Sanner, PhD. Copyright © 2026. See [License](#license).

**ADapt** prepares the receptor and ligand PDBQT files for AutoDock from a PDB entry, in a web
browser. It does what AGFRgui's PDB panel does. Open a PDB file (or fetch one from the RCSB) and
you get what that panel shows:

- one tab per biomolecule (ASU, Biomol N);
- per chain, the residues grouped by category: protein, nucleic, side chains to complete, missing atoms, cofactors, modified, manual charges, ligands, additives, water, other, plus the altloc and unknown-atom-type subtrees;
- receptor and ligand checkboxes, and a params column with altloc choice, manual ion charge and rotamer choice;
- a right-click menu to complete a side chain or leave it as is, to mutate a residue, and to replace an ion AutoDock cannot type;
- two save buttons that write **receptor.pdbqt** and **ligand.pdbqt**.

ADapt's computation is **ADFR's own code**: the same `runAGFR` / `PDBProfiler` / `ADLigand` calls AGFRgui makes, run in the ADCP distribution's Python (the `adcpsuite` env, Python 3.7, with **ADFR 0.0.23 or later and MolKit2 0.0.6**; `bin/adapt check` refuses an older ADFR). `tests/test_parity.py` checks that the PDBQT files are byte-identical to what ADFR's headless `--toPdbqt` path writes for the same choices; every receptor choice ADapt offers has an `agfr` option. The Qt GUI itself is not part of the tests. Only the Qt tree and its flag logic were ported (`tree.py`, `flags.py`).

ADapt's web technology follows ADCP_docking (github.com/aadhilhaq/ADCP_docking): a standard-library `ThreadingHTTPServer` bound to localhost with a random token in the URL, one HTML page with plain JavaScript and no build step, and a bundled 3Dmol.js viewer (BSD-3; its licence is in `adapt/static/`).

**Name.** The display name is **ADapt** (capital A, capital D: it encodes "AutoDock" and reads as
the word "adapt"); it is what the page, `--help`, the server's messages and the session records
show. Identifiers are lower case: the Python package `adapt` and the command `bin/adapt`.

## Run

```bash
bin/adapt check                     # "ADapt is ready", or what is missing
bin/adapt serve                     # port 8766; prints http://127.0.0.1:8766/?t=<token>
```

The printed URL carries a random token; the page answers 403 without it. The server stays up
until stopped (Ctrl-C, or `kill` its PID when it was started detached, e.g.
`setsid nohup bin/adapt serve > server.log 2>&1 &`).

`bin/adapt` sets up the ADCP distribution's `adcpsuite` env the way the `agfr` / `agfrgui`
launchers do (`ADS_ROOT`, no user site-packages); ADFR finds `reduce` through `ADS_ROOT`. It finds
the environment by itself: `$ADCPSUITE`, `$ADS_ROOT`, the environment of an `adcp` or `agfr` on the
`PATH`, `$CONDA_PREFIX`, then an `adcpsuite` env in the usual micromamba / conda locations. Set
`ADCPSUITE=/path/to/env` to choose one.

**From another machine**, run the server on the host that has the env, then forward its port:

```bash
ssh -L 8766:127.0.0.1:8766 <host>   # on your machine; then open the printed URL there
```

**Keep a fixed port** (the default, or any `--port N`). The browser keeps what it remembers for
the page per host and port: in Chrome and Edge **Open PDB…** reopens the folder the previous file
came from (File System Access API picker with an `id`), and the light/dark choice is stored the
same way. `--port 0` picks a new port on each run, and that memory starts empty every time.
Firefox and Safari lack that API and use a plain file input, whose starting folder the browser
decides.

## Try it

**Fetch from RCSB**: type a PDB id in the box at the top and press Enter. The server downloads
`https://files.rcsb.org/download/<ID>.pdb`, keeps it under `<outdir>/_rcsb/` (so the same entry
opens from the cache the next time), and opens it. An entry the RCSB distributes only as mmCIF (the
largest assemblies, entries with five-character ligand codes such as 10HY) is fetched as
`<ID>.cif`. The session record keeps the URL the file came from.

**Open PDB…** reads a local PDB file (`.pdb`, `.ent`) or mmCIF file (`.cif`, `.mmcif`), gzipped or
not; dropping one on the page does the same. An mmCIF file is converted to PDB records (ADFR's
`mmcif.py`); chain ids longer than one character and residue names longer than three are renamed,
and the files written record the renaming. **Recent…** lists the sessions of the server and of its
output folder, and reopens one, also after a restart. ADFR ships the files of its own tests; `bin/adapt examples` prints where they
are in the installation in use. The same entries can be fetched by id:

| Entry  | Tabs | Altloc residues | Side chains to complete | Manual-charge ions | Ligands (ADFR) | Missing residues      | Good for                                              |
|--------|------|-----------------|-------------------------|--------------------|----------------|-----------------------|-------------------------------------------------------|
| `6GGD` | 3    | 14              | 13                      | 2 (ZN)             | 2 (EYB)        | 44, termini           | the fullest single example                            |
| `4A1X` | 3    | 20              | 1                       | 15 (2 ZN, 13 CL)   | -              | 31, 2 gaps + termini  | missing residues (see below)                          |
| `1H6W` | 2    | 1               | -                       | -                  | -              | 161, N-terminus       | Biomol 1 is a trimer: symmetry copies at half opacity |
| `1A3R` | 1    | -               | -                       | -                  | -              | 6, 1 gap + N-terminus | insertion codes (chain L 27A, 27B, 27C) to mutate     |
| `6BHG` | 1    | 19              | 20                      | -                  | -              | 33, 1 gap + termini   | chain B as the ligand (ADFR's own GUI test)           |
| `1JFF` | 1    | -               | -                       | 2 (MG, ZN)         | 1 (TA1)        | 32, termini (1)       | metals of two kinds; a ligand ADFR knows              |
| `5HHT` | 1    | 33              | -                       | 2 (CA)             | -              | 9, termini (2)        | many altlocs; a His tag                               |
| `5FQ2` | 1    | -               | -                       | -                  | -              | 8, termini (3)        | a modified residue in the chain                       |
| `1WKW` | 1    | -               | -                       | -                  | 1 (GTA)        | -                     |                                                       |
| `1QJ6` | 1    | -               | -                       | -                  | -              | -                     |                                                       |

Counts are for the ASU tab, measured by ADapt on the RCSB files; ADFR's copies give the same.
"-" means none. **Ligands (ADFR)** are the residues ADFR's table of PDB entries names as ligands
(they start in the **ligands** category); other hetero groups go to cofactors, additives or
**other**. **Missing residues** is what the
tree's *missing residues* nodes list:
(1) 1JFF A 35-60 is listed as missing in the file (REMARK 465), but the residues on either side are
2.68 A apart, so ADFR treats the chain as continuous there and shows no gap;
(2) the missing last His of 5HHT B's His tag is shown as a missing C-terminus;
(3) MLZ B505, a methylated lysine written as HETATM, is in the file, so it is not reported missing.

**Missing residues**, with 4A1X: fetch `4A1X` and expand **missing residues** under a chain.
Each entry is a stretch the file lists as missing (REMARK 465): `gap: SER-3--GKGG--PRO2` means
GKGG is missing between SER -3 and PRO 2; `N term:` / `C term:` are missing ends (here chain A's
His tag, `GSHHHHHH`). In the 3D view each gap is a dotted orange line from the C of the residue
before it to the N of the residue after it, with a sphere at its middle labelled with the missing
residues; a missing end is a sphere 3 A beyond the terminal atom. Clicking a gap highlights the
residues on either side of it, clicking a missing end the residue it is missing from. The ASU tab shows all four chains' stretches (8 markers); Biomol 1 and
Biomol 2 show their own two chains'. Missing residues are information: nothing is rebuilt, and
the PDBQT files are the same with or without them.

A first session with 6GGD, following ADFR's `testAGFR.py` test2:

1. Fetch `6GGD` and open the **Biomol 2** tab.
2. Under **manual charges** and **ligands**, check the receptor box of `ZN'401` and `EYB'402`
   (or check `EYB'402` in the **lig.** column to write it as the ligand instead).
3. Right-click the protein residue `SER'96` and choose `ILE`.
4. Pick another alternate location of an altloc residue in the params column.
5. **receptor PDBQT** / **ligand PDBQT** write the files, download them, and list them under
   *Files written*; the session folder shown there keeps them.

After steps 1-3 (keeping ZN and EYB in the receptor), the receptor PDBQT is byte-identical to
what ADFR's command line writes for the same choices (checked on the fetched 6GGD; the
`testAGFR.py` scenarios are in `tests/test_parity.py`):

```bash
# with the adcpsuite env's python, ADS_ROOT set to that env, in a folder holding 6ggd.pdb
python -m ADFR.bin.runAGFR -r 6ggd.pdb --toPdbqt --bioMol 2 \
    --heteroToKeep "chid B resnum 402 401" --mutate "B:SER:96::ILE"
```

(Not through the `agfr` launcher: it passes its arguments as an unquoted `$@`, so
`"chid B resnum 402 401"` is split into words and rejected.)

Step 4 is `--altlocs "B:LYS:101::A"` on the command line; a chosen rotamer is `--rotamers`, a side
chain kept as is `--sidechainsAsIs`, an ion replacement `--substitute`.

## Outputs

Each loaded file gets a session folder under `--outdir` (default `$ADAPT_OUTDIR`, else `~/adapt_runs/`; `/tmp` is refused). The folder holds:

- the input PDB, under its own name, because ADFR reads the PDB id from it to look up the known ligands;
- `save_NNN_rec/` and `save_NNN_lig/` for every save, so nothing is overwritten. A receptor save holds `<id>_rec.pdbqt` and `<id>_rec_H.pdb`: the receptor as `reduce` protonated it, with every hydrogen (non-polar ones included), the PDBQT's `REMARK 850` header, then `reduce`'s output, including its `USER  MOD` lines (each His/Asn/Gln decision: protonation, flip, score, category K/F/C/X). ADFR keeps that output (`runAGFR.reduceOutput`, `reduceChainRenames`) and `adapt/reduce_log.py` parses it. The one change to `reduce`'s text: its added hydrogens all have serial number 0, so ADapt numbers every atom 1..N in file order (hexadecimal from the 100,000th, as ProDy does);
- `session.json`, which records the input md5, the ADFR and MolKit2 package versions and the md5 of their sources that ran, the ADapt version and commit, every edit in order (with a reference to its tree node), and every file written with its md5 and a checksum of the tree's state. A session is reopened from this record (`Session.restore`): the edits are replayed and each save's state checked, and a changed input or a mismatch is refused. A receptor save also records `reduce`: its version, its His/Asn/Gln `decisions` (`{"residue": "L:HIS:198:_", "decision": "no HD1", "score": 0.771, "category": "K", ...}`) and the `chain_renames` made for symmetry copies (`{"2": {"segment": "C", "chain": "A"}}`).

## Embedding ADapt

For a program that starts ADapt (as ADCP_docking does); details in `docs/start.html#embed`.

- **URL parameters** next to `t`: `pdb=6GGD` fetches and opens an entry, `session=<id>` opens a session the server holds (e.g. made through `/api/upload`) or one recorded in its output folder, `tab=Biomol%202` switches to a tab. After a load the address bar carries `session=` and `tab=`, so a reload reopens the session, also after a restart.
- **Requests** (all with `t=`): `POST api/upload?name=<file>` (file as body), `GET api/sessions` (the sessions held, newest first, then those on disk with `"live": false`), `GET api/session?s=<id>`, `POST api/restore` (`{"folder": name}`), `POST api/save` (`{"s": id, "kind": "receptor"|"ligand"}`), `GET api/download?s=<id>&save=<n>[&file=<name>]` (one file; without `file`, a several-file save is a zip).
- **Path prefix.** The page and script use relative paths only, so a page served at `/adapt/?t=…` calls `/adapt/api/…` and a proxy needs no rewriting. For a proxy that forwards the full path, `--base-path /adapt` (or `$ADAPT_BASE_PATH`) strips the prefix, redirects `/adapt?…` to `/adapt/?…` and prints the URL with it; `/` still works.

## Where ADapt differs from AGFRgui

`docs/differences.html` has the full list. In short:

- **Corrected in ADFR 0.0.23** (found while building ADapt, so `agfr` and AGFRgui benefit too, each
  with a test in ADFR's `tests/` that fails on 0.0.22): deterministic default alternate location;
  ion charges from the Chemical Component Dictionary (FE +3, FE2 +2, MN +2, BR/CL/IOD −1; FE2 and
  MN3 no longer stop a file from opening; bonded Cl/Br charged, MolKit2 0.0.6); ions of symmetry
  copies keep their charge; residues with atoms AutoDock 4.1 cannot type stay out of the receptor,
  and an unsupported ion can be replaced; the chosen rotamer is built; missing residues are
  detected (ADFR's `_gaps` never reported one); assemblies from the tenth on are the right ones;
  mmCIF input, and gzipped PDB files keep their header; the header names the biomolecule and the
  ligand's residues, reports ion charges as written, and lists only residues of the file; new
  options `--altlocs`, `--rotamers`, `--sidechainsAsIs`, `--substitute`.
- **Only in ADapt:** symmetry copies follow every choice made on their original chain
  (`adapt/copies.py`; a ligand is taken from the original only; copies on a symmetry axis are
  written once); `<id>_rec_H.pdb` with `reduce`'s His/Asn/Gln choices; the session record and
  restore; and fixes of AGFRgui's GUI defects: an edited charge is used (AGFRgui dropped the spin
  box value), a save after a tab switch writes the new tab, an altloc chosen in one tab applies in
  a tab opened later, edits made after a save take effect, clicking a label only focuses the
  view, and the mutate menu of a modified residue opens.

## Known issues

All of them, with details and links, are in `docs/known-issues.html`. Unless marked ADapt, they
are ADFR's behaviour, shared by `agfr` and AGFRgui.

- **Symmetry copies get new chain identifiers** (`0`, `1`, …) in the receptor files (PDBQT, and the atoms and `USER  MOD` lines of `_rec_H.pdb`); ADFR does not rename them back after `reduce`. ADapt records the renames in `session.json` (`reduce.chain_renames`).
- **A ligand's atom records say chain `A`** whatever chain it came from; the header's `REMARK 850 AGFR ligand from B:EYB:402:_` records the true origin.
- **(ADapt) Some choices are in the session record only:** which hetero groups were kept, and ADapt's version, are in `session.json`, not in the PDBQT header.
- **A gap whose ends meet is not shown.** A stretch listed as missing whose neighbouring residues are within 3.5 Å (C to N) is treated as continuous (1JFF GLY34-HIS61, 2.68 Å).
- **(ADapt) Symmetry copies are not in the tree**; they are drawn translucent and follow their original chain.
- **A PDB file's id comes from its file name.** The id selects the known ligands of the entry; name PDB files by their id. An mmCIF file carries its id.
- **Copies defined only by MTRIX records are not built**; assemblies come from REMARK 350.
- **The `agfr` / `agfrgui` launchers split option values containing spaces** (unquoted `$@`). Call `python -m ADFR.bin.runAGFR` directly for those. The correction belongs to the ADFR packaging, for the next ADCP release.
- **(ADapt) Needs the ADCP distribution's `adcpsuite` environment** (Python 3.7) with ADFR 0.0.23 or later.

## Documentation

`docs/` holds the user documentation as static HTML: an overview of the problem ADapt addresses,
a getting-started page, one page per case (biological assembly vs ASU, choosing between
assemblies, side chains, alternate locations, ligands, metals and unsupported atoms, missing
residues, mutations), how the files are written, and how ADapt relates to AGFRgui. A running
ADapt serves it at `/docs/`, linked from the **Docs** button; the files also open directly from
`docs/index.html`.

Every screenshot is of ADapt itself. To remake them after the page changes, and rebuild the pages:

```bash
python3 docs/make_screenshots.py [scene ...]   # headless Chrome against a fresh ADapt server
python3 docs/build_docs.py                     # docs/src/*.html -> docs/*.html
```

`make_screenshots.py` needs `google-chrome` and access to files.rcsb.org. Its sessions go to
`~/adapt_runs/docs_screenshots/`, and its Chrome profile and log to `tests/_out/docs_run/`.

## Tests

```bash
bin/adapt test [-k name]            # = tests/run_tests.py in the adcpsuite env
```

There is no pytest in the shared env; `run_tests.py` is a small runner. Outputs go to `tests/_out/<timestamp>/` (git-ignored). The suite covers:

- byte parity with `agfr --toPdbqt` for all 10 ADFR test PDBs and for the `testAGFR.py` scenarios (bioMol, heteroToKeep, mutateRes, protToRemove, ionCharges, sidechainsAsIs), and ligand parity with ADFR's `processLigands`;
- tree coverage: every supported atom lies in exactly one leaf;
- the flag semantics;
- the HTTP layer on a real socket;
- Fetch from RCSB (id checks, download + cache, an mmCIF-only entry fetched and opened, the
  server route); the network tests skip when the machine has no route to the RCSB;
- missing residues against REMARK 465 for the ten test PDBs, and the 3D markers;
- the protonated receptor PDB (same header and heavy atoms as the PDBQT, every polar H of the
  PDBQT, `reduce`'s text with only the atom numbers changed, numbered 1..N, past 99999 too;
  `reduce`'s decisions and chain renames in `session.json`);
- symmetry copies (`test_copies.py`: a residue unchecked, an altloc, an ion charge in every copy;
  a ligand from the original only; 4INS's axial zincs written once; a chain of a later REMARK 350
  block in the tree);
- the rotamer chosen is built, an ion replaced and restored, residues left out recorded
  (`test_choices.py`);
- sessions restored after a restart, with every kind of edit, byte-identical saves, a changed
  input and an ambiguous old record refused (`test_restore.py`);
- opening a session from elsewhere, and serving under `/adapt/` through a proxy that rewrites
  nothing (`--base-path`, relative paths).

## License

Author: Michel Sanner, PhD. Copyright © 2026. All rights reserved: no license has been chosen yet,
and one will be added to this repository.

Parts of ADapt come from other software, under their own licences. `THIRD_PARTY_NOTICES.md` lists
them:

- 3Dmol.js, bundled in `adapt/static/` (BSD 3-Clause);
- code adapted from ADCP_docking (BSD 3-Clause);
- code derived from ADFR / AGFRgui, which ADapt drives but does not include.
