#!/usr/bin/env python3

import stepwise, appcli, autoprop
from stepwise import Quantity
from stepwise_mol_bio import Cleanup
from freezerbox import MakerArgsConfig, group_by_identity
from appcli import DocoptConfig, Key

@autoprop
class Aliquot(Cleanup):
    """\
Make aliquots

Usage:
    aliquot <volume> [<conc>]

Arguments:
    <volume>
        The volume of each individual aliquot.  No unit is implied, so you 
        must specify one.

    <conc>
        The concentration of the aliquots, if this is not made clear in 
        previous steps.  No unit is implied, so you must specify one.
"""
    __config__ = [
            DocoptConfig(),
            MakerArgsConfig(),
    ]
    volume = appcli.param(
            Key(DocoptConfig, '<volume>'),
            Key(MakerArgsConfig, 'volume'),
    )
    conc = appcli.param(
            Key(DocoptConfig, '<conc>'),
            Key(MakerArgsConfig, 'conc'),
            default=None,
            ignore=None,
    )

    group_by = {
        'volume': group_by_identity,
        'conc': group_by_identity,
    }

    def __init__(self, volume, conc=None, products=None):
        self.volume = volume
        self.conc = conc
        self.products = products or []

    def get_protocol(self):
        Q = Quantity.from_string

        if self.conc:
            aliquot_info = f'{Q(self.volume)}, {Q(self.conc)}'
        else:
            aliquot_info = f'{Q(self.volume)}'

        if self.products and self.show_product_names:
            product_names = f" of: {', '.join(self.products)}"
        else:
            product_names = "."

        return stepwise.Protocol(
                steps=[f"Make {aliquot_info} aliquots{product_names}"],
        )

if __name__ == '__main__':
    Aliquot.main()
