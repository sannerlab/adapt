# Author: Michel Sanner, PhD. Copyright (c) 2026.
"""The commit of an ADapt snapshot made with `git archive` (for instance the
copy an ADCPsuite installer unpacks): git fills these in when it writes the
archive (.gitattributes: export-subst). In a git checkout they stay as they
are, and env.provenance asks git instead."""
COMMIT = '$Format:%H$'
DATE = '$Format:%cI$'
