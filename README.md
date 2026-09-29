# agfr_web — AGFRgui's PDB panel in a web browser

Open a PDB file and you get what AGFRgui's PDB analyser shows:

- one tab per biomolecule (ASU, Biomol N);
- per chain, the residues grouped by category: protein, nucleic, side chains to complete, missing atoms, cofactors, modified, manual charges, ligands, additives, water, other, plus the altloc and unknown-atom-type subtrees;
- receptor and ligand checkboxes, and a params column with altloc choice, manual ion charge and rotamer preview;
- a right-click menu to complete a side chain or leave it as is, and to mutate a residue;
- two save buttons that write **receptor.pdbqt** and **ligand.pdbqt**.

The computation is **ADFR's own code**: the same `runAGFR` / `PDBProfiler` / `ADLigand` calls AGFRgui makes, run in the ADCP distribution's Python (the `adcpsuite` env, Python 3.7, `adfr-0.0.22`). The PDBQT files are therefore byte-identical to AGFRgui's for the same choices; `tests/test_parity.py` checks this against ADFR's headless `agfr --toPdbqt` path. Only the Qt tree and its flag logic were ported (`tree.py`, `flags.py`).

The web technology follows ADCP_docking (github.com/aadhilhaq/ADCP_docking): a standard-library `ThreadingHTTPServer` bound to localhost with a random token in the URL, one HTML page with plain JavaScript and no build step, and a bundled 3Dmol.js viewer (BSD-3; its licence is in `agfr_web/static/`).

## Run

```bash
bin/agfr-web check                  # "ready", or what is missing
bin/agfr-web serve --port 0         # prints http://127.0.0.1:<port>/?t=<token>
```

`bin/agfr-web` sets up the adcpsuite env the way the `agfr` / `agfrgui` launchers do (`ADS_ROOT`, no user site-packages); ADFR finds `reduce` through `ADS_ROOT`. Set `ADCPSUITE=/path/to/env` to use another installation.

From another machine, forward the port: `ssh -L 8766:127.0.0.1:8766 bali`, then run `bin/agfr-web serve --port 8766` on bali and open the printed URL locally.

Each loaded file gets a session folder under `--outdir` (default `/home/sanner/claude/agfr_web_runs/`, on NFS; `/tmp` is refused). The folder holds:

- the input PDB, under its own name, because ADFR reads the PDB id from it to look up the known ligands;
- `save_NNN_rec/` and `save_NNN_lig/` for every save, so nothing is overwritten;
- `session.json`, which records the input md5, the ADFR package version and the md5 of the ADFR sources that ran, the agfr_web commit, every edit in order, and every file written with its md5.

## Differences from AGFRgui (all fixes of GUI defects)

- **An edited manual charge is used.** AGFRgui's `_getCharges` reads a checkbox column the charge items never have, so the spin box value was silently dropped and the default charge written. With unedited charges both write the same bytes.
- **Tab switches apply correctly.**
  - AGFRgui kept the processed molecule across a tab switch, so a save right after switching wrote the previous biomolecule. Here it writes the new one.
  - An altloc chosen in one tab is applied to the atoms of a tab opened afterwards. AGFRgui showed the choice, and wrote it in the `REMARK 850 AGFR altlocs` line, but kept the file's default coordinates.
- **Mutations and side-chain choices always take effect.** In AGFRgui, one made after a save was ignored until some checkbox changed.
- **Clicking a label only focuses the viewer.** In AGFRgui, clicking the label of a partially checked category made every residue under it receptor while the checkbox still showed partial.
- **The "unknown atom types" subtree is display-only.** Its residues are not user-checkable in AGFRgui either.

## Faithful, but worth knowing (ADFR behaviour, not changed)

- The rotamer dropdown is a **preview**. When the receptor is written, `mutateResidues` builds the best-scoring rotamer, as in AGFRgui.
- `recFromPDB._gaps` never reports a gap in adfr-0.0.22: it tests `format_alignment(...).split('\n').count('-')`, which counts lines that are exactly `-`. So "missing residues" never appears here or in AGFRgui. The node is ported and tested with a synthetic record.
- The order of the letters in `PDBProfiler._altlocsRes` varies between runs, because it comes from set iteration. The chosen letter, and hence the output, does not.
- The `custom charge` REMARK is written from the ASU molecule's flags, so it is absent when the ion was added to the receptor in another biomolecule's tab. The atom still carries the charge.

## Tests

```bash
PATH=/home/sanner/micromamba/envs/adcpsuite/bin:$PATH ADS_ROOT=/home/sanner/micromamba/envs/adcpsuite \
    python tests/run_tests.py [-k name]
```

There is no pytest in the shared env; `run_tests.py` is a small runner. Outputs go to `tests/_out/<timestamp>/` (git-ignored). The suite covers:

- byte parity with `agfr --toPdbqt` for all 10 ADFR test PDBs and for the `testAGFR.py` scenarios (bioMol, heteroToKeep, mutateRes, protToRemove), and ligand parity with ADFR's `processLigands`;
- tree coverage: every supported atom lies in exactly one leaf;
- the flag semantics;
- the HTTP layer on a real socket.
