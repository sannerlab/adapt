# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""Writing the receptor and ligand PDBQT files.

  receptor: runAGFR.processBM() (MSE->MET, ion substitutions, mutations and
            side-chain completion with the chosen or best rotamer, altlocs
            cleared), then processReceptor(mol.select('receptor'), charges)
            (reduce, Open Babel Gasteiger charges, ion charges, AD types,
            non-polar H merged), then the file: ADFR's REMARK 850 AGFR
            provenance lines (the assembly, mutations, rotamers, alternate
            locations, ...) and atom lines, as runAGFR.saveReceptorAsPDBQT
            writes them. Next to it, <pdbid>_rec_H.pdb: the same header, then
            what reduce wrote (its USER MOD lines and every atom with all of
            its hydrogens), with the atoms numbered in order (reduce numbers
            the hydrogens it adds 0).
  ligand:   processBM(), then processLigands(mol.select('ligand')), one
            ADLigand per connected fragment, written with getPDBQTlines (its
            header names the source file, the assembly and the residues the
            fragment came from).

File names: <pdbid>_rec.pdbqt, <pdbid>_lig.pdbqt, or <pdbid>_lig_NNN.pdbqt when
the ligand atoms form several fragments. Each save goes into its own numbered
folder of the session so no earlier output is ever overwritten.
"""
import os


from . import flags


class Nothing(Exception):
    """No atoms are flagged for what was asked."""


def _processed_bm(agfr):
    mol = agfr._BMProcessed
    if mol is None:
        mol = agfr.processBM()
    return mol


def _write_lines(lines, filename):
    with open(filename, 'w') as f:
        for l in lines:
            f.write('%s\n' % l)


def receptor(session, folder):
    agfr = session.agfr
    mol = _processed_bm(agfr)
    atoms = mol.select('receptor')
    if atoms is None:
        raise Nothing('no atoms selected for the receptor')
    charges = flags.charges(session.tree())
    session.profile._customCharges = charges        # PDBAnalyserWidget.getCharges
    processed = agfr.processReceptor(atoms, charges=charges)
    path = os.path.join(folder, session.profile._pdbid + '_rec.pdbqt')
    # runAGFR.saveReceptorAsPDBQT
    from ADFR import PDBQTlinesForAtoms
    lines = session.profile.getPDBQTHeaderLines()
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
        header = session.profile.getPDBQTHeaderLines()
        with open(pdb, 'w') as fh:
            for l in header:
                fh.write('%s\n' % l)
            out = reduce_log.renumbered(out)
            fh.write(out if out.endswith('\n') else out + '\n')
        paths.append(pdb)
    return paths


def ligands(session, folder):
    agfr = session.agfr
    mol = _processed_bm(agfr)
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
