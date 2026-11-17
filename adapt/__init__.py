# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""ADapt: AutoDock Analysis and Preparation Tool - preparing PDB entries for docking.

Reads a PDB file, shows one tab per biomolecule with its residues grouped the
way AGFRgui's PDB panel groups them, and writes the receptor and ligand PDBQT
files with ADFR's own code (the ADCP distribution's adcpsuite env).

NAME: the display/brand name is ADapt (capital A, capital D: it encodes
"AutoDock" and reads as the word "adapt"). Every user-facing string takes it
from APP / APP_LONG below. Identifiers are the lower-case `adapt`: this
package and the command bin/adapt (both `agfr_web` / `bin/agfr-web` before
2026-09-29).
"""
__version__ = '0.1.0'
__author__ = 'Michel Sanner, PhD'
__copyright__ = 'Copyright (c) 2026 Michel Sanner, PhD'

APP = 'ADapt'
APP_LONG = 'AutoDock Analysis and Preparation Tool'
AUTHOR = __author__
COPYRIGHT = __copyright__
