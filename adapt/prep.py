# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Writing the receptor and ligand PDBQT files.

  receptor: runAGFR.processBM() (MSE->MET, ion substitutions, mutations and
            side-chain completion with the chosen or best rotamer, altlocs
            cleared), then the missing residues the user chose to rebuild, if
            any (adapt/rebuild.py: in chain order, so reduce protonates them
            with the rest), then processReceptor(mol.select('receptor'),
            charges) (reduce, Open Babel Gasteiger charges, ion charges, AD
            types, non-polar H merged), then the file: ADFR's REMARK 850 AGFR
            provenance lines (the assembly, mutations, rotamers, alternate
            locations, ...), with a line per rebuilt gap added by ADapt before
            the protonation line, and atom lines, as
            runAGFR.saveReceptorAsPDBQT writes them. Next to it,
            <pdbid>_rec_H.pdb: the same header, then what reduce wrote (its
            USER MOD lines and every atom with all of its hydrogens), with the
            atoms numbered in order (reduce numbers the hydrogens it adds 0).
  ligand:   processBM(), then processLigands(mol.select('ligand')), one
            ADLigand per connected fragment, written with getPDBQTlines (its
            header names the source file, the assembly and the residues the
            fragment came from).

File names: <pdbid>_rec.pdbqt, <pdbid>_lig.pdbqt, or <pdbid>_lig_NNN.pdbqt when
the ligand atoms form several fragments. Each save goes into its own numbered
folder of the session so no earlier output is ever overwritten.
"""
import os


from . import flags, rebuild


class Nothing(Exception):
    """No atoms are flagged for what was asked."""


def processed_bm(session):
    """runAGFR.processBM(), kept until an edit invalidates it, with the gaps
    chosen for rebuilding inserted (Rebuilder.apply; with none chosen it
    returns processBM's molecule itself). The ligand path shares it: the new
    residues are never ligand atoms."""
    agfr = session.agfr
    mol = agfr._BMProcessed
    if mol is None:
        mol = agfr.processBM()
        try:
            mol = session.rebuilder.apply(mol, session.profile._currentBM)
        except Exception:
            agfr._BMProcessed = None    # processBM cached a molecule without them
            raise
        agfr._BMProcessed = mol
    return mol


def _write_lines(lines, filename):
    with open(filename, 'w') as f:
        for l in lines:
            f.write('%s\n' % l)


def receptor(session, folder):
    agfr = session.agfr
    mol = processed_bm(session)
    atoms = mol.select('receptor')
    if atoms is None:
        raise Nothing('no atoms selected for the receptor')
    charges = flags.charges(session.tree())
    session.profile._customCharges = charges        # PDBAnalyserWidget.getCharges
    processed = agfr.processReceptor(atoms, charges=charges)
    path = os.path.join(folder, session.profile._pdbid + '_rec.pdbqt')
    # runAGFR.saveReceptorAsPDBQT, plus a line per rebuilt gap
    from ADFR import PDBQTlinesForAtoms
    header = _with_rebuilt(session.profile.getPDBQTHeaderLines(), mol)
    for r in getattr(mol, 'rebuilt', ()):
        if not r['written']:
            print('WARNING: rebuilt gap /%s/%s is not written: the residues on both '
                  'sides of it are not both in the receptor' % (r['segment'], r['gap']))
    lines = list(header)
    lines.extend(PDBQTlinesForAtoms(processed._ag.all))
    _write_lines(lines, path)
    processed.filename = path
    agfr._processedReceptor = processed
    paths = [path]
    # the receptor as reduce protonated it: all hydrogens, non-polar ones too,
    # and reduce's choices (USER MOD lines), under the same provenance header
    out = agfr.reduceOutput
    if out:
        from . import reduce_log
        pdb = os.path.join(folder, session.profile._pdbid + '_rec_H.pdb')
        with open(pdb, 'w') as fh:
            for l in header:
                fh.write('%s\n' % l)
            out = reduce_log.renumbered(out)
            fh.write(out if out.endswith('\n') else out + '\n')
        paths.append(pdb)
    return paths


def _with_rebuilt(lines, mol):
    """The header with a line per gap rebuilt in the receptor (rebuild.remarks)
    before "protonated using reduce": the residues were built after the side
    chains and before reduce, and the header lists the steps in that order."""
    more = rebuild.remarks(mol)
    if not more:
        return lines
    lines = list(lines)
    at = next((i for i, l in enumerate(lines)
               if l.startswith('REMARK 850 AGFR protonated')), len(lines) - 1)
    lines[at:at] = more
    return lines


def ligands(session, folder):
    agfr = session.agfr
    mol = processed_bm(session)
    ligAtoms = mol.select('ligand')
    if ligAtoms is None:
        raise Nothing('no atoms selected for the ligand')
    # one ADLigand per connected fragment; its header names the residues it
    # was made of (REMARK 850 AGFR ligand from ...)
    agfr._processedLigands = agfr.processLigands(ligAtoms)
    filepath = os.path.join(folder, session.profile._pdbid + '_lig.pdbqt')
    paths = []
    if len(agfr._adLigs) > 1:
        basename = os.path.splitext(filepath)[0]
        for num, adlig in enumerate(agfr._adLigs):
            filename = '%s_%03d.pdbqt' % (basename, num + 1)
            _write_lines(adlig.getPDBQTlines(), filename)
            agfr._processedLigands[num].filename = filename
            paths.append(filename)
    else:
        _write_lines(agfr._adLigs[0].getPDBQTlines(), filepath)
        agfr._processedLigands[0].filename = filepath
        paths.append(filepath)
    return paths
