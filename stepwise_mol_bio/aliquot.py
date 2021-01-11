#!/usr/bin/env python3

import stepwise, appcli, autoprop
from stepwise_mol_bio import Main

@autoprop
class Aliquot(Main):
    """\
Make aliquots

Usage:
    aliquot <volume> [<conc>]

Arguments:
    <volume>
        The volume of each individual aliquot.  No unit is implied, so you 
        should specify one.

    <conc>
        The concentration of the aliquots, if this is not made clear in 
        previous steps.  No unit is implied, so you should specify one.
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]
    volume = appcli.param('<volume>')
    conc = appcli.param('<conc>', default=None, ignore=None)

    def __init__(self, volume, conc=None):
        self.volume = volume
        self.conc = conc

    def get_protocol(self):
        p = stepwise.Protocol()

        if self.conc:
            p += f"Make {self.volume}, {self.conc} aliquots."
        else:
            p += f"Make {self.volume} aliquots."

        return p

if __name__ == '__main__':
    Aliquot.main()
