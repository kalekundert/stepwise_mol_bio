#!/usr/bin/env python3

import autoprop
from stepwise import Protocol, Quantity
from stepwise_mol_bio import Cleanup
from appcli import DocoptConfig

@autoprop
class Miniprep(Cleanup):
    """
Purify plasmid by miniprep.

Usage:
    miniprep

Database:
    The miniprep protocol can be used in the "Cleanup" column of a FreezerBox 
    database:

        miniprep

    Miniprepped plasmids will be assumed to have a concentration of 200 ng/µL, 
    unless otherwise specified.
"""
    __config__ = [
            DocoptConfig(),
    ]

    def get_protocol(self):
        p = Protocol()
        p += "Miniprep."
        return p

    def get_product_conc(self):
        return Quantity(200, 'ng/µL')

if __name__ == '__main__':
    Miniprep.main()
