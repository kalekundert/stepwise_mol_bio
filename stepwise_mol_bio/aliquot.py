#!/usr/bin/env python3

import stepwise, byoc, autoprop
from stepwise import Quantity
from stepwise_mol_bio import Cleanup
from freezerbox import MakerConfig, group_by_identity
from byoc import DocoptConfig, Key, Method

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
            DocoptConfig,
            MakerConfig,
    ]
    volume = byoc.param(
            Key(DocoptConfig, '<volume>'),
            Key(MakerConfig, 'volume'),
    )
    conc = byoc.param(
            Key(DocoptConfig, '<conc>'),
            Key(MakerConfig, 'conc'),
            default=None,
    )

    group_by = {
        'volume': group_by_identity,
        'conc': group_by_identity,
    }

    def __init__(self, volume, conc=None, product_tags=None):
        self.volume = volume
        if conc: self.conc = conc
        if product_tags: self.product_tags = product_tags

    def get_protocol(self):
        Q = Quantity.from_string

        if self.conc:
            aliquot_info = f'{Q(self.volume)}, {Q(self.conc)}'
        else:
            aliquot_info = f'{Q(self.volume)}'

        if self.product_tags and self.show_product_tags:
            product_tags = f" of: {', '.join(self.product_tags)}"
        else:
            product_tags = "."

        return stepwise.Protocol(
                steps=[f"Make {aliquot_info} aliquots{product_tags}"],
        )

    def get_product_conc(self):
        return Quantity.from_string(self.conc)

if __name__ == '__main__':
    Aliquot.main()
