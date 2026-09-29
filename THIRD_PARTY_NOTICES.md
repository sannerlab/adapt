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
