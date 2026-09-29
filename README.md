# ADapt — AutoDock preparation & perception tool

**ADapt** prepares the receptor and ligand PDBQT files for AutoDock from a PDB entry, in a web
browser. It does what AGFRgui's PDB panel does. Open a PDB file (or fetch one from the RCSB) and
you get what that panel shows:

- one tab per biomolecule (ASU, Biomol N);
- per chain, the residues grouped by category: protein, nucleic, side chains to complete, missing atoms, cofactors, modified, manual charges, ligands, additives, water, other, plus the altloc and unknown-atom-type subtrees;
- receptor and ligand checkboxes, and a params column with altloc choice, manual ion charge and rotamer preview;
- a right-click menu to complete a side chain or leave it as is, and to mutate a residue;
- two save buttons that write **receptor.pdbqt** and **ligand.pdbqt**.

ADapt's computation is **ADFR's own code**: the same `runAGFR` / `PDBProfiler` / `ADLigand` calls AGFRgui makes, run in the ADCP distribution's Python (the `adcpsuite` env, Python 3.7, `adfr-0.0.22`). `tests/test_parity.py` checks that the PDBQT files are byte-identical to what ADFR's headless `--toPdbqt` path writes for the same choices; the Qt GUI itself is not part of the tests. Only the Qt tree and its flag logic were ported (`tree.py`, `flags.py`).

ADapt's web technology follows ADCP_docking (github.com/aadhilhaq/ADCP_docking): a standard-library `ThreadingHTTPServer` bound to localhost with a random token in the URL, one HTML page with plain JavaScript and no build step, and a bundled 3Dmol.js viewer (BSD-3; its licence is in `agfr_web/static/`).

**Name.** The display name is **ADapt** (capital A, capital D: it encodes "AutoDock" and reads as
the word "adapt"); it is what the page, `--help`, the server's messages and the session records
show. The Python package (`agfr_web`) and the command (`bin/agfr-web`) keep their technical names.

## Run

```bash
bin/agfr-web check                  # "ready", or what is missing
bin/agfr-web serve                  # port 8766; prints http://127.0.0.1:8766/?t=<token>
```

The printed URL carries a random token; the page answers 403 without it. The server stays up
until stopped (Ctrl-C, or `kill` its PID when it was started detached, e.g.
`setsid nohup bin/agfr-web serve > server.log 2>&1 &`).

`bin/agfr-web` sets up the ADCP distribution's `adcpsuite` env the way the `agfr` / `agfrgui`
launchers do (`ADS_ROOT`, no user site-packages); ADFR finds `reduce` through `ADS_ROOT`. Set
`ADCPSUITE=/path/to/env` to use another installation.

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
does the same. ADFR ships the files of its own tests; `bin/agfr-web examples` prints where they
are in the installation in use. The same entries can be fetched by id:

| entry | what it exercises |
|---|---|
| `6GGD` | 3 tabs (ASU, Biomol 1, Biomol 2); altloc residues; side chains to complete (rotamer preview); a ZN manual charge; the ligand EYB 402. The fullest single example |
| `1H6W` | Biomol 1 is a trimer: the symmetry copies are drawn at half opacity |
| `4A1X` | 3 tabs, altlocs, one side chain to complete |
| `1A3R` | one tab; chain L has insertion codes (SER 27A, LEU 27B, 27C) to mutate |
| `6BHG` | make chain B the ligand (ADFR's own GUI test does this) |
| `1JFF`, `1QJ6`, `1WKW`, `5FQ2`, `5HHT` | the rest of ADFR's test set |

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

Each loaded file gets a session folder under `--outdir` (default `/home/sanner/claude/agfr_web_runs/`, on NFS; `/tmp` is refused). The folder holds:

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

## Faithful, but worth knowing (ADFR behaviour, not changed)

- The rotamer dropdown is a **preview**. When the receptor is written, `mutateResidues` builds the best-scoring rotamer, as in AGFRgui.
- `recFromPDB._gaps` never reports a gap in adfr-0.0.22: it tests `format_alignment(...).split('\n').count('-')`, which counts lines that are exactly `-`. So "missing residues" never appears in ADapt or in AGFRgui. The node is ported and tested with a synthetic record.
- The order of the letters in `PDBProfiler._altlocsRes` varies between runs, because it comes from set iteration. The chosen letter, and hence the output, does not.
- The `agfr` / `agfrgui` launcher scripts pass arguments as an unquoted `$@`, so any option value containing spaces (`--heteroToKeep "chid B resnum 402"`) is split and rejected. Call `python -m ADFR.bin.runAGFR` directly for those.
- The `custom charge` REMARK is written from the ASU molecule's flags, so it is absent when the ion was added to the receptor in another biomolecule's tab. The atom still carries the charge.

## Tests

```bash
bin/agfr-web test [-k name]         # = tests/run_tests.py in the adcpsuite env
```

There is no pytest in the shared env; `run_tests.py` is a small runner. Outputs go to `tests/_out/<timestamp>/` (git-ignored). The suite covers:

- byte parity with `agfr --toPdbqt` for all 10 ADFR test PDBs and for the `testAGFR.py` scenarios (bioMol, heteroToKeep, mutateRes, protToRemove), and ligand parity with ADFR's `processLigands`;
- tree coverage: every supported atom lies in exactly one leaf;
- the flag semantics;
- the HTTP layer on a real socket;
- Fetch from RCSB (id checks, download + cache, mmCIF-only entries, the server route); the
  network tests skip when the machine has no route to the RCSB.
