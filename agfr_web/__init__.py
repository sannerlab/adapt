"""ADapt, the AutoDock preparation & perception tool.

Reads a PDB file, shows one tab per biomolecule with its residues grouped the
way AGFRgui's PDB panel groups them, and writes the receptor and ligand PDBQT
files with ADFR's own code (the ADCP distribution's adcpsuite env).

NAME: the display/brand name is ADapt (capital A, capital D: it encodes
"AutoDock" and reads as the word "adapt"). Every user-facing string takes it
from APP / APP_LONG below. The package (agfr_web) and the command
(bin/agfr-web) keep their technical names.
"""
__version__ = '0.1.0'

APP = 'ADapt'
APP_LONG = 'AutoDock preparation & perception tool'
