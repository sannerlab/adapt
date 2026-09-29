# ADapt — AutoDock preparation & perception tool

Author: Michel Sanner, PhD. Copyright © 2026. See [License](#license).

**ADapt** prepares the receptor and ligand PDBQT files for AutoDock from a PDB entry, in a web
browser. It does what AGFRgui's PDB panel does. Open a PDB file (or fetch one from the RCSB) and
you get what that panel shows:

- one tab per biomolecule (ASU, Biomol N);
- per chain, the residues grouped by category: protein, nucleic, side chains to complete, missing atoms, cofactors, modified, manual charges, ligands, additives, water, other, plus the altloc and unknown-atom-type subtrees;
- receptor and ligand checkboxes, and a params column with altloc choice, manual ion charge and rotamer preview;
- a right-click menu to complete a side chain or leave it as is, and to mutate a residue;
- two save buttons that write **receptor.pdbqt** and **ligand.pdbqt**.

ADapt's computation is **ADFR's own code**: the same `runAGFR` / `PDBProfiler` / `ADLigand` calls AGFRgui makes, run in the ADCP distribution's Python (the `adcpsuite` env, Python 3.7, `adfr-0.0.22`). `tests/test_parity.py` checks that the PDBQT files are byte-identical to what ADFR's headless `--toPdbqt` path writes for the same choices; the Qt GUI itself is not part of the tests. Only the Qt tree and its flag logic were ported (`tree.py`, `flags.py`).

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
opens from the cache the next time), and opens it. Entries the RCSB distributes only as mmCIF (the
largest assemblies) are reported, not opened, because ADapt reads PDB files. The session
record keeps the URL the file came from.

**Open PDB…** reads a local `.pdb`, `.ent`, `.pdb.gz` or `.ent.gz` file; dropping one on the page
does the same. ADFR ships the files of its own tests; `bin/adapt examples` prints where they
are in the installation in use. The same entries can be fetched by id:

| Entry  | Tabs | Altloc residues | Side chains to complete | Manual-charge ions | Ligands (ADFR) | Missing residues      | Good for                                              |
|--------|------|-----------------|-------------------------|--------------------|----------------|-----------------------|-------------------------------------------------------|
| `6GGD` | 3    | 14              | 13                      | 2 (ZN)             | 2 (EYB)        | 44, termini           | the fullest single example                            |
| `4A1X` | 3    | 20              | 1                       | 2 (ZN)             | -              | 31, 2 gaps + termini  | missing residues (see below)                          |
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

The command line has no altloc option, so step 4 is a choice only the panel can make.

## Outputs

Each loaded file gets a session folder under `--outdir` (default `$ADAPT_OUTDIR`, else `~/adapt_runs/`; `/tmp` is refused). The folder holds:

- the input PDB, under its own name, because ADFR reads the PDB id from it to look up the known ligands;
- `save_NNN_rec/` and `save_NNN_lig/` for every save, so nothing is overwritten;
- `session.json`, which records the input md5, the ADFR package version and the md5 of the ADFR sources that ran, the ADapt version and commit, every edit in order, and every file written with its md5.

## Where ADapt differs from AGFRgui (all fixes of GUI defects)

- **An edited manual charge is used.** AGFRgui's `_getCharges` reads a checkbox column the charge items never have, so the spin box value was silently dropped and the default charge written. With unedited charges both write the same bytes.
- **Tab switches apply correctly.**
  - AGFRgui kept the processed molecule across a tab switch, so a save right after switching wrote the previous biomolecule. Here it writes the new one.
  - An altloc chosen in one tab is applied to the atoms of a tab opened afterwards. AGFRgui showed the choice, and wrote it in the `REMARK 850 AGFR altlocs` line, but kept the file's default coordinates.
- **Mutations and side-chain choices always take effect.** In AGFRgui, one made after a save was ignored until some checkbox changed.
- **Clicking a label only focuses the viewer.** In AGFRgui, clicking the label of a partially checked category made every residue under it receptor while the checkbox still showed partial.
- **The "unknown atom types" subtree is display-only.** Its residues are not user-checkable in AGFRgui either.
- **Missing residues are shown.** ADFR's `PDBProfiler._gaps` never reported one: it tests `format_alignment(...).split('\n').count('-')`, which counts lines that are exactly `-`, so "missing residues" never appeared in AGFRgui. ADapt uses a copy of that function (`adapt/gaps.py`, the installed ADFR is not changed) with this test fixed and four defects it had been hiding: a chain with no missing C-terminus would have reported its whole sequence as one; modified residues and caps written as HETATM (MLZ, ALY, NH2) were reported missing; a missing end of a His tag could be dropped; a chain without SEQRES raised an error. Checked against the REMARK 465 records of ADFR's ten test PDBs (`tests/test_gaps.py`). The 3D markers are a port of AGFRgui's `addGeomsForGaps`.

## Faithful, but worth knowing (ADFR behaviour, not changed)

- The rotamer dropdown is a **preview**. When the receptor is written, `mutateResidues` builds the best-scoring rotamer, as in AGFRgui.
- A stretch listed as missing whose neighbouring residues are within 3.5 A (C to N) is not shown: ADFR treats the chain as continuous there, and its source names the case it was written for (1JFF GLY34-HIS61, 2.68 A).
- **The default alternate location is not reproducible when occupancies tie.** ADFR picks the location with the most atoms, then the highest occupancy, and breaks a tie by taking the last tied letter in `set` order, which Python varies from run to run (hash randomisation). 1FV1 D:SER19 (A 0.50 / B 0.50) defaults to A in some runs and B in others, and so does the written receptor. Ties are common: all 7 altloc residues of 1FV1, 2 of 4A1X's 20, 2 of 5HHT's 33, 4 of 6BHG's 19. The `REMARK 850 AGFR altlocs` line records the letter that was used; choose tied residues explicitly when it matters. This is ADFR's behaviour (agfr and AGFRgui alike) and is not changed here.
- The `agfr` / `agfrgui` launcher scripts pass arguments as an unquoted `$@`, so any option value containing spaces (`--heteroToKeep "chid B resnum 402"`) is split and rejected. Call `python -m ADFR.bin.runAGFR` directly for those.
- The `custom charge` REMARK is written from the ASU molecule's flags, so it is absent when the ion was added to the receptor in another biomolecule's tab. The atom still carries the charge.

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

- byte parity with `agfr --toPdbqt` for all 10 ADFR test PDBs and for the `testAGFR.py` scenarios (bioMol, heteroToKeep, mutateRes, protToRemove), and ligand parity with ADFR's `processLigands`;
- tree coverage: every supported atom lies in exactly one leaf;
- the flag semantics;
- the HTTP layer on a real socket;
- Fetch from RCSB (id checks, download + cache, mmCIF-only entries, the server route); the
  network tests skip when the machine has no route to the RCSB;
- missing residues against REMARK 465 for the ten test PDBs, the 3D markers, and that ADFR's
  own `_gaps` reports none.

## License

Author: Michel Sanner, PhD. Copyright © 2026. All rights reserved: no license has been chosen yet,
and one will be added to this repository.

Parts of ADapt come from other software, under their own licences. `THIRD_PARTY_NOTICES.md` lists
them:

- 3Dmol.js, bundled in `adapt/static/` (BSD 3-Clause);
- code adapted from ADCP_docking (BSD 3-Clause);
- code derived from ADFR / AGFRgui, which ADapt drives but does not include.
