# Third-party notices

ADapt is Copyright © 2026 Michel Sanner, PhD. It includes, adapts, or is derived from the
following software, which keeps its own licence and copyright.

## 3Dmol.js (bundled)

`adapt/static/3Dmol-min.js`, the 3D viewer: 3Dmol.js 2.5.5, BSD 3-Clause. Its licence texts are
bundled next to it: `adapt/static/3Dmol.LICENSE.txt` and `adapt/static/3Dmol-min.js.LICENSE.txt`.
Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. *Bioinformatics* 31(8):1322-1324,
2015.

## ADCP_docking (adapted code)

The web server pattern (`adapt/server.py`: a standard-library server on localhost with a URL
token, JSON routes, static files), the page conventions (`adapt/page.py`,
`adapt/static/app.js`), and the RCSB fetch (`adapt/fetch.py`) were adapted from ADCP_docking,
https://github.com/aadhilhaq/ADCP_docking, under this licence:

```
BSD 3-Clause License

Copyright (c) 2026, Aadhil Haq and Sandun Fernando, HOPE Lab, Texas A&M University

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## ADFR / AGFRgui (derived code; not included)

ADapt drives ADFR, from the ADCP distribution, to do all of its chemistry. ADFR is not included
here, and each user installs it. These files are derived from ADFR (`adfr-0.0.22`):

- `adapt/tree.py`, `adapt/flags.py` and `adapt/rotamer.py` port the logic of AGFRgui's PDB panel
  (`ADFR/GUI/pdbAnalyser.py`);
- `adapt/gaps.py` contains a copy of `ADFR/recFromPDB.py:PDBProfiler._gaps`, with marked fixes.

ADFR is developed at the Scripps Research Institute; see the ADFR distribution for its licence.

## PDBFixer and OpenMM (used; not included)

`adapt/loops.py` builds the residues of a gap the user chose to rebuild with PDBFixer 1.8.1
(`PDBFixer.addMissingAtoms`) and refines them with OpenMM 7.6.0 and its Amber ff14SB parameter
files (`amber14/protein.ff14SB.xml`), all as installed in the ADCP distribution's `adcpsuite`
environment. No PDBFixer or OpenMM code is copied into ADapt; both are imported, and only when a
gap is rebuilt. PDBFixer's licence (from its installed `LICENSE`):

```
This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2013 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
```

OpenMM is distributed under the MIT licence and, for parts of its C++ library, the GNU LGPL; see
the OpenMM distribution. Eastman P, et al. OpenMM 7: rapid development of high performance
algorithms for molecular dynamics. *PLoS Comput Biol* 13(7):e1005659, 2017. Maier JA, et al.
ff14SB: improving the accuracy of protein side chain and backbone parameters from ff99SB.
*J Chem Theory Comput* 11(8):3696, 2015.
