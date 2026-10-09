# Author: Aadhil Haq. Copyright (c) 2026.
"""Building the residues of a gap: PDBFixer places them, Amber ff14SB makes them sane.

build() is the whole step. Its input is the PDB text adapt/rebuild.py writes for it
(the observed heavy atoms near the gaps, one chain per copy of a chain, residues
numbered 1, 2, ...) and the gaps to fill. Its output is, per gap, the heavy atoms
of the new residues; or NotBuilt, with the reason, when the loop it ends with is
not sane. Nothing here knows about ADFR.

1. PDBFixer (addMissingAtoms) is asked for these gaps and nothing else:
   missingResidues from ADapt's gaps, missingAtoms and missingTerminals empty.
   It lays the residues on an arc between the two ends and relaxes them in its
   soft force field while the observed atoms are held fixed (massless).

2. That force field is soft on purpose (angle constants of 10 kJ/mol/rad^2, a
   repulsion capped at about 1 kJ/mol), so the loop it leaves is a sketch, not a
   model. Measured on the four buildable gaps of ADFR's test entries: peptide
   bonds anywhere from cis to trans, consecutive CA 1.4 to 4.1 A apart, and in
   three of the four loops at least one CA of the wrong hand (D). Next to a Cys
   in a disulfide (1A3R H CYS 128) PDBFixer has no template for the flanking
   residue, the join carries no force, and C to N ends 5.6 A apart.

   refine() rebuilds that geometry in Amber ff14SB (OpenMM's amber14 files).
   The new residues and every hydrogen move; every observed heavy atom is
   fixed. The residues on either side of the gap take part with their own
   parameters; the other observed atoms within ENV of the loop act on it
   through Lennard-Jones terms with ff14SB's values for their element. Charges
   are off: the aim is covalent geometry without clashes, not a model of the
   loop's electrostatics in vacuum. Every new peptide bond is kept trans and
   every new CA L by flat-bottomed restraints, which exert no force on a loop
   inside their wells; each new side chain is first turned about its CA so CB
   sits where an L residue's would (_cb_ideal). A soft repulsion stage (seeded
   Langevin dynamics, so the loop can pass through itself and its
   surroundings) comes before the Lennard-Jones minimisation.

3. measure() checks the loop against TOLERANCE. When it is outside, the whole
   build is tried again with the next seed of SEEDS; when no seed gives a loop
   inside, the gap is refused (NotBuilt) with what was measured, never written.
   On ADFR's test gaps (1A3R H, 4A1X A, 6BHG A) and seeds 1 to 8, 23 of 24
   first tries passed; the miss had CA to CA 3.96 A. 4A1X B (GK) fails every
   try: an exhaustive search shows two residues cannot close it.

Everything random is seeded (SEEDS) and everything runs on OpenMM's Reference
platform: the CPU platform gives another loop for another number of threads
(loops up to 7.6 A apart, measured), the Reference platform the same one
everywhere. PDBFixer makes its Context without naming a platform; OpenMM then
takes the one $OPENMM_DEFAULT_PLATFORM names, which it reads for every Context,
so the variable is set around that one call (_platform) and checked.
Modeller.addHydrogens places new hydrogens with Python's random module before
relaxing them, so that call runs with the module seeded (_seeded) and its state
restored.
"""
import contextlib
import io
import math
import os
import random
import time
from collections import OrderedDict

import numpy

SEEDS = (1, 2, 3)       # tried in turn, the first loop that passes is kept (OpenMM: never 0)
PLATFORM = 'Reference'
FORCEFIELD = 'amber14/protein.ff14SB.xml'
METHOD = 'Amber ff14SB'  # the name the REMARK and the record give the refinement

ENV = 12.0              # A: observed atoms this close to the loop take part in refine()
CUTOFF = 0.9            # nm, Lennard-Jones cutoff (switched from 0.8 nm)
SOFT_STEPS = 1000       # Langevin steps of 0.5 fs in the soft stage
KT = 400.0              # kJ/mol/rad^2: trans restraint on each new peptide bond
TRANS_SLACK = 15.0      # deg: free play around 180
KC = 400.0              # kJ/mol/rad^2: L restraint on each new CA
L_IMPROPER = 122.0      # deg: improper N-C-CA-CB of an L residue (ADFR's test entries: 114 to 129)
L_SLACK = 20.0          # deg: free play around it

# What a loop must measure to be written (A, degrees). C to N: a peptide bond is
# 1.33 A. CA to CA across a trans peptide is 3.80 A in crystal structures; minimised
# in ff14SB, an observed stretch (4A1X A 20 to 25) measures 3.84 to 3.89 A, so the
# upper bound leaves room for the force field's own minimum. A peptide bond is trans
# within 30 degrees of 180; cis within 30 degrees of 0 (CA to CA about 2.9 A). A CA
# is L when its improper N-C-CA-CB is positive, and sound within 25 degrees of
# L_IMPROPER (observed residues next to ADFR's test gaps: 114 to 129 degrees).
# Closest contact: a heavy atom of a new residue and an observed heavy atom more than
# three bonds apart; closer than 2.5 A is a clash (hydrogen bonds are 2.6 A and up).
TOLERANCE = {'C-N': (1.30, 1.37), 'CA-CA': (3.70, 3.95), 'trans': 150.0, 'improper': 25.0,
             'contact': 2.5}

# The chiral CB of Thr and Ile: (X - CB) x (Y - CB) . (CA - CB) > 0, as in every
# Thr (143) and Ile (72) of ADFR's test entries 1A3R, 4A1X, 6BHG, 6GGD and 1H6W
CHIRAL_CB = {'THR': ('OG1', 'CG2'), 'ILE': ('CG1', 'CG2')}

# Lennard-Jones sigma (nm) and epsilon (kJ/mol) of ff14SB's types CT, N, O and S,
# for observed atoms outside the residues next to the gap; other elements as CT
LJ = {'C': (0.339967, 0.457730), 'N': (0.325, 0.711280), 'O': (0.295992, 0.878640),
      'S': (0.356359, 1.046000)}


class NotBuilt(Exception):
    """The loop built for a gap is not sane; str() says what was measured."""

    def __init__(self, site, why):
        Exception.__init__(self, why)
        self.site = site


def available():
    import importlib.util
    return all(importlib.util.find_spec(m) is not None for m in ('openmm', 'pdbfixer'))


_VERSIONS = []


def versions():
    """(PDBFixer, OpenMM) versions, as their installations name them. PDBFixer
    1.8.1's own pdbfixer.__version__ still says 1.7; the package says 1.8.1."""
    if _VERSIONS:
        return _VERSIONS[0]
    import pkg_resources
    out = []
    for dist, fallback in (('pdbfixer', 'pdbfixer.pdbfixer'), ('OpenMM', 'openmm.version')):
        try:
            out.append(pkg_resources.get_distribution(dist).version)
        except Exception:                                     # noqa: BLE001
            import importlib
            mod = importlib.import_module(fallback)
            out.append(getattr(mod, '__version__', None) or getattr(mod, 'short_version', '?'))
    _VERSIONS.append(tuple(out))
    return _VERSIONS[0]


class Loop(object):
    """One gap in PDBFixer's topology: the residue before it, the new residues,
    the residue after it (OpenMM Residues), and what refine() and measure() made
    of them: the new heavy atoms [(residue position, name, element)], their
    coordinates (A), their bonds, the measured geometry, the seed, and the
    seconds the build of its group took."""
    __slots__ = ('before', 'new', 'after', 'atoms', 'coords', 'bonds', 'geometry', 'seed',
                 'seconds')


def build(text, chains, sites):
    """[Loop], one per site of `sites` (in that order), for the gaps PDBFixer is
    asked to fill in `text`. chains[k] lists the residues of chain k of `text`;
    each site names its gap's residues (site.where(chains): chain and residue
    index of the residue after the gap; site.names: the residues to build;
    site.labels: those either side and the new ones, as faults name them).
    Each seed of SEEDS is tried in turn; raises NotBuilt, with the faults of
    the first try, when no loop passes."""
    start, first = time.time(), None
    for k, seed in enumerate(SEEDS):
        t = time.time()
        fixer, loops = place(text, chains, sites, seed)
        placed = time.time() - t
        env = refine(fixer, loops, seed)
        seconds = {'pdbfixer': round(placed, 2), 'refine': round(time.time() - t - placed, 2),
                   'tries': k + 1, 'total': round(time.time() - start, 2)}
        faults = []
        for site, loop, frag in zip(sites, loops, env['fragments']):
            loop.geometry = measure(frag, env, site.labels)
            loop.seed, loop.seconds = seed, seconds
            why = verdict(loop.geometry)
            if why:
                faults.append((site, why))
        if not faults:
            return loops
        first = first or faults[0]
    raise NotBuilt(first[0], 'none of %d loops built for it passed the checks; the first had %s'
                   % (len(SEEDS), first[1]))


# ---- 1. PDBFixer ----------------------------------------------------------
def place(text, chains, sites, seed):
    """PDBFixer on `text`, asked to build the gaps of `sites` and nothing else."""
    from pdbfixer.pdbfixer import PDBFixer
    quiet_import()
    fixer = PDBFixer(pdbfile=io.StringIO(text))
    topo = [list(c.residues()) for c in fixer.topology.chains()]
    if len(topo) != len(chains):
        raise RuntimeError('PDBFixer read %d chains, %d were written' % (len(topo), len(chains)))
    missing = OrderedDict()
    for s in sites:
        k, j = s.where(chains)
        if topo[k][j].id != str(j + 1):         # rebuild._pdb numbers residues 1, 2, ...
            raise RuntimeError('PDBFixer did not read residue %d of chain %d where it was written'
                               % (j + 1, k + 1))
        missing[(k, j)] = list(s.names)
    fixer.missingResidues = missing
    fixer.missingAtoms = {}             # side chains: ADFR's rotamers and the user's choices
    fixer.missingTerminals = {}
    with _platform():
        fixer.addMissingAtoms(seed=seed)
    return fixer, _inserted(fixer, missing)


@contextlib.contextmanager
def _platform():
    """$OPENMM_DEFAULT_PLATFORM = PLATFORM for what runs inside, and a check
    that this OpenMM honours it."""
    import openmm as mm
    old = os.environ.get('OPENMM_DEFAULT_PLATFORM')
    os.environ['OPENMM_DEFAULT_PLATFORM'] = PLATFORM
    try:
        probe = mm.System()
        probe.addParticle(1.0)
        used = mm.Context(probe, mm.VerletIntegrator(0.001)).getPlatform().getName()
        if used != PLATFORM:
            raise RuntimeError('OpenMM chose the %s platform, not %s' % (used, PLATFORM))
        yield
    finally:
        if old is None:
            del os.environ['OPENMM_DEFAULT_PLATFORM']
        else:
            os.environ['OPENMM_DEFAULT_PLATFORM'] = old


def _inserted(fixer, missing):
    """[Loop] in the order of `missing`: addMissingAtoms inserts each list
    before residue j of chain k (PDBFixer._addAtomsToTopology)."""
    found = {}
    for k, chain in enumerate(fixer.topology.chains()):
        residues = list(chain.residues())
        old = new = 0
        while new < len(residues):
            if (k, old) in missing:
                m = len(missing[(k, old)])
                found[(k, old)] = (residues[new - 1], residues[new:new + m], residues[new + m])
                new += m
            new += 1
            old += 1
    out = []
    for where, names in missing.items():
        loop = Loop()
        loop.before, loop.new, loop.after = found[where]
        if [r.name for r in loop.new] != names:
            raise RuntimeError('PDBFixer built %s, not %s' % ([r.name for r in loop.new], names))
        out.append(loop)
    return out


_QUIET = []


def quiet_import():
    """PDBFixer's soft.xml imports simtk.openmm, and OpenMM 7.6 prints a
    deprecation warning the first time; make that time here, silently."""
    if not _QUIET:
        import warnings
        with warnings.catch_warnings(), contextlib.redirect_stderr(io.StringIO()), \
                contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('ignore')
            import simtk.openmm      # noqa: F401
        _QUIET.append(True)


# ---- 2. Amber ff14SB --------------------------------------------------------
_FF = []


def _forcefield():
    import openmm.app as app
    if not _FF:
        _FF.append(app.ForceField(FORCEFIELD))
    return _FF[0]


@contextlib.contextmanager
def _seeded(seed):
    state = random.getstate()
    random.seed(seed)
    try:
        yield
    finally:
        random.setstate(state)


def refine(fixer, loops, seed):
    """Refine the loops PDBFixer made, all in one system (gaps built together
    may touch). Sets each loop's atoms, coords and bonds; returns what
    measure() needs: the final positions (A) of refine()'s topology, its
    fragments (one residue list per loop, the residues either side included)
    and the observed heavy atoms around them."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm.app.element import hydrogen
    top, xyz, pieces = _fragment(fixer, loops)
    model = app.Modeller(top, [mm.Vec3(*p) for p in xyz] * unit.nanometer)
    with _seeded(seed):
        model.addHydrogens(pH=7.0, platform=mm.Platform.getPlatformByName(PLATFORM))
    system = _forcefield().createSystem(
        model.topology, nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=CUTOFF * unit.nanometer, switchDistance=(CUTOFF - 0.1) * unit.nanometer,
        constraints=None, removeCMMotion=False, ignoreExternalBonds=True)
    residues = list(model.topology.residues())
    frags = []
    for n in pieces:
        frags.append(residues[:n])
        residues = residues[n:]
    # the new residues move, and every hydrogen: those of the residues either
    # side were placed from PDBFixer's provisional atoms
    mobile = sorted(a.index for rr in frags for r in rr for a in r.atoms()
                    if r not in (rr[0], rr[-1]) or a.element is hydrogen)
    for i in set(range(system.getNumParticles())) - set(mobile):
        system.setParticleMass(i, 0.0)
    positions = [numpy.array(p) for p in model.positions.value_in_unit(unit.nanometer)]
    observed = _observed(fixer, loops, numpy.array([positions[i] for i in mobile]) * 10.)
    _nonbonded(system, mobile, observed)
    positions.extend(x / 10. for e, x in observed)
    _soft(system, model.topology, mobile)
    _restraints(system, frags)
    integrator = mm.LangevinIntegrator(300 * unit.kelvin, 5 / unit.picosecond,
                                       0.5 * unit.femtosecond)
    integrator.setRandomNumberSeed(seed)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName(PLATFORM))
    context.setPositions(positions)
    # soft stage: the loop may pass through itself and its surroundings
    mm.LocalEnergyMinimizer.minimize(context)
    integrator.step(SOFT_STEPS)
    mm.LocalEnergyMinimizer.minimize(context)
    # Lennard-Jones stage. The restraints stay: their wells are flat (TRANS_SLACK,
    # L_SLACK), so a relaxed loop feels none of them; measured on ADFR's test gaps,
    # releasing them let a strained loop invert a CA or twist a peptide bond
    context.setParameter('soft', 0.0)
    context.setParameter('lj', 1.0)
    mm.LocalEnergyMinimizer.minimize(context)
    final = context.getState(getPositions=True).getPositions(asNumpy=True)
    final = final.value_in_unit(unit.angstrom)
    for loop, rr in zip(loops, frags):
        _take(loop, rr, model.topology, final)
    n = model.topology.getNumAtoms()
    return {'final': final, 'topology': model.topology, 'fragments': frags, 'observed': final[n:]}


def _fragment(fixer, loops):
    """An OpenMM Topology with one chain per loop: the residue before the gap,
    the new residues, the residue after it, heavy atoms only, the new side
    chains turned to an ideal CB (_cb_ideal) and their chiral CB righted
    (_cb_hand); the positions (nm); the residue count of each chain."""
    import openmm.app as app
    import openmm.unit as unit
    from openmm.app.element import hydrogen
    xyz = numpy.array(fixer.positions.value_in_unit(unit.angstrom))   # as _cb_ideal needs
    top = app.Topology()
    out, pieces = [], []
    for loop in loops:
        chain = top.addChain()
        rr = [loop.before] + loop.new + [loop.after]
        pieces.append(len(rr))
        for r in rr:
            res = top.addResidue(r.name, chain, r.id, r.insertionCode)
            atoms = [a for a in r.atoms() if a.element not in (None, hydrogen)]
            pos = xyz[[a.index for a in atoms]]
            if r in loop.new:
                names = [a.name for a in atoms]
                pos = _cb_hand(r.name, names, _cb_ideal(r.name, names, pos))
            for a, p in zip(atoms, pos):
                top.addAtom(a.name, a.element, res)
                out.append(p)
    top.createStandardBonds()
    return top, numpy.array(out) / 10., pieces


def _cb_ideal(resname, names, pos):
    """The positions of a residue with its side chain turned rigidly about CA so
    that CB points where an L residue's would (from N, CA and C, in A; the
    constants place a CB 1.52 A from CA with ideal angles). PDBFixer's soft
    dynamics leave some CA of the D hand and some half flattened; hydrogens
    placed on such a CA put HA where CB belongs, and the refinement then pushes
    CB over. A rotation keeps the side chain's own geometry, a chiral CB (Thr,
    Ile) included. Not for Pro, whose ring it would open at N; a Pro left on
    the wrong side fails the checks."""
    at = dict(zip(names, pos))
    if 'CB' not in at or resname == 'PRO':
        return pos
    n, ca, c = at['N'], at['CA'], at['C']
    b, d = ca - n, c - ca
    want = -0.58273431 * numpy.cross(b, d) + 0.56802827 * b - 0.54067466 * d
    have = at['CB'] - ca
    axis = numpy.cross(have, want)
    sin, cos = numpy.linalg.norm(axis), numpy.dot(have, want)
    if sin < 1e-6 * numpy.linalg.norm(have) * numpy.linalg.norm(want):
        return pos
    axis /= sin
    angle = math.atan2(sin, cos)
    k = numpy.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    turn = numpy.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * k.dot(k)    # Rodrigues
    side = numpy.array([name not in ('N', 'CA', 'C', 'O') for name in names])
    pos = pos.copy()
    pos[side] = (pos[side] - ca).dot(turn.T) + ca
    return pos


def _cb_hand(resname, names, pos):
    """The positions of a Thr or Ile with the wrong hand at CB (an allo form;
    PDBFixer's dynamics make them) righted by swapping two of its substituents:
    OG1 and CG2, or CG1 and CG2 with CD1 carried along with CG1."""
    at = dict(zip(names, pos))
    if resname not in CHIRAL_CB or not all(n in at for n in CHIRAL_CB[resname] + ('CA', 'CB')) \
            or _cb_right(resname, at):
        return pos
    x, y = CHIRAL_CB[resname]
    k = dict((n, i) for i, n in enumerate(names))
    pos = pos.copy()
    pos[k[x]], pos[k[y]] = at[y], at[x]
    if 'CD1' in k:
        pos[k['CD1']] = at['CD1'] + at[y] - at[x]
    return pos


def _observed(fixer, loops, mobile_xyz):
    """[(element, xyz in A)] of PDBFixer's heavy atoms outside refine()'s
    fragments, within ENV of a mobile atom (mobile_xyz, A)."""
    import openmm.unit as unit
    from openmm.app.element import hydrogen
    taken = set(r for loop in loops for r in [loop.before] + loop.new + [loop.after])
    xyz = numpy.array(fixer.positions.value_in_unit(unit.angstrom))
    atoms = [a for a in fixer.topology.atoms()
             if a.residue not in taken and a.element not in (None, hydrogen)]
    if not atoms:
        return []
    pos = xyz[[a.index for a in atoms]]
    near = numpy.zeros(len(atoms), bool)
    for p in mobile_xyz:
        near |= ((pos - p) ** 2).sum(1) <= ENV * ENV
    return [(a.element.symbol, pos[k]) for k, a in enumerate(atoms) if near[k]]


def _nonbonded(system, mobile, observed):
    """ff14SB's NonbondedForce with charges off, the mobile atoms' Lennard-Jones
    scaled by the global 'lj' (0 in the soft stage), and the observed atoms
    added as fixed particles."""
    import openmm as mm
    import openmm.unit as unit
    KJ = unit.kilojoule_per_mole
    nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
    moving = set(mobile)
    nb.addGlobalParameter('lj', 0.0)
    for i in range(nb.getNumParticles()):
        q, sigma, eps = nb.getParticleParameters(i)
        if i in moving:
            nb.setParticleParameters(i, 0.0, sigma, 0.0)
            nb.addParticleParameterOffset('lj', i, 0.0, 0.0, eps.value_in_unit(KJ))
        else:
            nb.setParticleParameters(i, 0.0, sigma, eps)
    for k in range(nb.getNumExceptions()):
        a, b, qq, sigma, eps = nb.getExceptionParameters(k)
        if a in moving or b in moving:
            nb.setExceptionParameters(k, a, b, 0.0, sigma, 0.0)
            nb.addExceptionParameterOffset('lj', k, 0.0, 0.0, eps.value_in_unit(KJ))
        else:
            nb.setExceptionParameters(k, a, b, 0.0, sigma, eps)
    for element, xyz in observed:
        system.addParticle(0.0)
        sigma, eps = LJ.get(element, LJ['C'])
        nb.addParticle(0.0, sigma, eps)


def _soft(system, topology, mobile):
    """PDBFixer's kind of repulsion, stronger: 50 kJ/mol at contact, 'soft' to
    switch it off. Between a mobile atom and anything more than three bonds away."""
    import openmm as mm
    soft = mm.CustomNonbondedForce('soft*50/((r/0.25)^4+1)')
    soft.addGlobalParameter('soft', 1.0)
    for i in range(system.getNumParticles()):
        soft.addParticle([])
    soft.createExclusionsFromBonds([(a.index, b.index) for a, b in topology.bonds()], 3)
    soft.addInteractionGroup(mobile, range(system.getNumParticles()))
    system.addForce(soft)


def _restraints(system, frags):
    """Every peptide bond of a new residue trans (CA-C-N-CA within TRANS_SLACK
    of 180 degrees), every new CA L (within L_SLACK of L_IMPROPER)."""
    import openmm as mm
    trans = mm.CustomTorsionForce('0.5*kt*max(0, d-wt)^2; d=min(dt, 2*pi-dt); dt=abs(theta-pi)')
    trans.addGlobalParameter('kt', KT)
    trans.addGlobalParameter('wt', math.radians(TRANS_SLACK))
    trans.addGlobalParameter('pi', math.pi)
    hand = mm.CustomTorsionForce('0.5*kc*max(0, abs(d)-slack)^2; '
                                 'd=atan2(sin(theta-t0), cos(theta-t0))')
    hand.addGlobalParameter('kc', KC)
    hand.addGlobalParameter('slack', math.radians(L_SLACK))
    hand.addPerTorsionParameter('t0')
    for rr in frags:
        for a, b in zip(rr[:-1], rr[1:]):
            trans.addTorsion(_index(a, 'CA'), _index(a, 'C'), _index(b, 'N'), _index(b, 'CA'), [])
        for r in rr[1:-1]:
            if any(a.name == 'CB' for a in r.atoms()):
                hand.addTorsion(_index(r, 'N'), _index(r, 'C'), _index(r, 'CA'), _index(r, 'CB'),
                                [math.radians(L_IMPROPER)])
    system.addForce(trans)
    system.addForce(hand)


def _index(residue, name):
    return next(a.index for a in residue.atoms() if a.name == name)


def _take(loop, rr, topology, final):
    """The new residues' heavy atoms, coordinates and bonds, from refine()'s
    fragment (rr: its residues, before and after the gap included)."""
    from openmm.app.element import hydrogen
    atoms = [(k, a) for k, r in enumerate(rr[1:-1]) for a in r.atoms() if a.element is not hydrogen]
    local = dict((a.index, i) for i, (k, a) in enumerate(atoms))
    loop.atoms = [(k, a.name, a.element.symbol) for k, a in atoms]
    loop.coords = numpy.array([final[a.index] for k, a in atoms])
    loop.bonds = []
    for a1, a2 in topology.bonds():
        if a1.index in local and a2.index in local:
            loop.bonds.append((local[a1.index], local[a2.index]))


# ---- 3. the checks ------------------------------------------------------------
def measure(frag, env, label):
    """{'bonds': [{'pair', 'C-N', 'CA-CA', 'omega'}] for each new peptide bond
    (the two joins included), 'C-N' and 'CA-CA': [min, max] over them, 'cis':
    the pairs within 30 degrees of cis, 'CA': [{'residue', 'improper'}] for each
    new residue with a CB, 'allo': the new Thr and Ile whose CB has the wrong
    hand, 'contact': the closest heavy-atom distance to an observed atom}; A and
    degrees. `frag`: one of refine()'s fragments; `label`: its residues as
    faults name them."""
    xyz = env['final']
    pos = [dict((a.name, xyz[a.index]) for a in r.atoms()) for r in frag]
    bonds = []
    for k in range(len(frag) - 1):
        bonds.append({'pair': '%s to %s' % (label[k], label[k + 1]),
                      'C-N': round(_dist(pos[k]['C'], pos[k + 1]['N']), 3),
                      'CA-CA': round(_dist(pos[k]['CA'], pos[k + 1]['CA']), 3),
                      'omega': round(_dihedral(pos[k]['CA'], pos[k]['C'], pos[k + 1]['N'],
                                               pos[k + 1]['CA']), 1)})
    return {'bonds': bonds,
            'C-N': [min(b['C-N'] for b in bonds), max(b['C-N'] for b in bonds)],
            'CA-CA': [min(b['CA-CA'] for b in bonds), max(b['CA-CA'] for b in bonds)],
            'cis': [b['pair'] for b in bonds if abs(b['omega']) <= 180 - TOLERANCE['trans']],
            'CA': [{'residue': label[k], 'improper': round(_dihedral(
                pos[k]['N'], pos[k]['C'], pos[k]['CA'], pos[k]['CB']), 1)}
                for k in range(1, len(frag) - 1) if 'CB' in pos[k]],
            'allo': [label[k] for k in range(1, len(frag) - 1)
                     if not _cb_right(frag[k].name, pos[k])],
            'contact': round(_contact(frag, env), 2)}


def _cb_right(name, p):
    if name not in CHIRAL_CB:
        return True
    x, y = (p[a] - p['CB'] for a in CHIRAL_CB[name])
    return numpy.dot(numpy.cross(x, y), p['CA'] - p['CB']) > 0


def _contact(frag, env):
    """The closest distance between a heavy atom of a new residue and an
    observed heavy atom: of the residues either side when more than three
    bonds away, or of the surroundings refine() saw."""
    from openmm.app.element import hydrogen
    xyz = env['final']
    new = [a for r in frag[1:-1] for a in r.atoms() if a.element is not hydrogen]
    near = _within_bonds(env['topology'], set(a.index for a in new), 3)
    flank = [a.index for r in (frag[0], frag[-1]) for a in r.atoms()
             if a.element is not hydrogen and a.index not in near]
    other = numpy.concatenate([xyz[flank], env['observed']]) if len(env['observed']) else xyz[flank]
    if not len(other):
        return float('inf')
    p = xyz[[a.index for a in new]]
    return float(numpy.sqrt(((p[:, None] - other[None]) ** 2).sum(2)).min())


def _within_bonds(topology, start, n):
    bonded = {}
    for a, b in topology.bonds():
        bonded.setdefault(a.index, []).append(b.index)
        bonded.setdefault(b.index, []).append(a.index)
    seen, edge = set(start), set(start)
    for _ in range(n):
        edge = set(j for i in edge for j in bonded.get(i, ())) - seen
        seen |= edge
    return seen


def verdict(g):
    """What is wrong with a loop of geometry `g` ('; ' between faults), or None."""
    out = []
    if g['cis']:
        out.append('a cis peptide bond (%s)' % g['cis'][0])
    else:
        worst = min(g['bonds'], key=lambda b: abs(b['omega']))
        if abs(worst['omega']) < TOLERANCE['trans']:
            out.append('a peptide bond %.0f degrees from trans (%s)' % (180 - abs(worst['omega']),
                                                                       worst['pair']))
    for name, label in (('C-N', 'C to N'), ('CA-CA', 'CA to CA')):
        lo, hi = TOLERANCE[name]
        bad = [b for b in g['bonds'] if not lo <= b[name] <= hi]
        if bad:
            b = max(bad, key=lambda b: max(lo - b[name], b[name] - hi))
            out.append('%s %.2f A (%s; %.2f to %.2f A allowed)'
                       % (label, b[name], b['pair'], lo, hi))
    for ca in g['CA']:
        if ca['improper'] <= 0:
            out.append('%s of the D hand' % ca['residue'])
        elif abs(ca['improper'] - L_IMPROPER) > TOLERANCE['improper']:
            out.append('%s with its CA distorted (improper N-C-CA-CB %.0f degrees; '
                       '%.0f to %.0f allowed)'
                       % (ca['residue'], ca['improper'], L_IMPROPER - TOLERANCE['improper'],
                          L_IMPROPER + TOLERANCE['improper']))
    if g['allo']:
        out.append('%s with its CB of the wrong hand' % g['allo'][0])
    if g['contact'] < TOLERANCE['contact']:
        out.append('an atom %.2f A from an observed atom (%.1f A allowed)'
                   % (g['contact'], TOLERANCE['contact']))
    return '; '.join(out) or None


def _dist(p, q):
    return float(numpy.linalg.norm(p - q))


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / numpy.linalg.norm(b1)
    v = b0 - numpy.dot(b0, b1) * b1
    w = b2 - numpy.dot(b2, b1) * b1
    return math.degrees(math.atan2(numpy.dot(numpy.cross(b1, v), w), numpy.dot(v, w)))
