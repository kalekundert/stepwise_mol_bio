#!/usr/bin/env python3

"""\
Stain nucleic acid gels using GelRed.

GelRed is a sensitive, stable and environmentally safe fluorescent nucleic acid 
dye designed to replace the highly toxic ethidium bromide (EtBr) for staining 
dsDNA, ssDNA or RNA in agarose gels or polyacrylamide gels. GelRed and EtBr 
have virtually the same spectra, so you can directly replace EtBr 
with GelRed without changing your existing imaging system. In addition, 
GelRed® is far more sensitive than EtBr (Figure 2).

GelRed can be used to stain dsDNA, ssDNA or RNA, however GelRed is twice as 
sensitive for dsDNA than ssDNA or RNA. Gel staining with GelRed is compatible 
with downstream applications such as sequencing and cloning.  GelRed is 
efficiently removed from DNA by phenol/chloroform extraction and ethanol 
precipitation.

Usage:
    gelred [-a]

Options:
    -a --attach-pdf
        Attach the manufacturer's protocol to this one, so that both will be 
        printed out by `stepwise go`.
"""

import stepwise
import autoprop
from _biotium import Biotium

class GelRed(Biotium):
    product = 'GelRed'
    uv_wavelength = '302'
    attachment = 'biotium_gelgreen.pdf'

if __name__ == '__main__':
    GelRed.main(__doc__)
