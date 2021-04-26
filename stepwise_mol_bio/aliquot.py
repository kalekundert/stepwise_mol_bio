#!/usr/bin/env python3

import stepwise, appcli, autoprop
from stepwise import Quantity
from stepwise_mol_bio import Main
from freezerbox import MakerArgsConfig, iter_combo_makers, group_by_identity
from appcli import DocoptConfig, Key

@autoprop
class Aliquot(Main):
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

    def __bareinit__(self):
        self.products = []
        self.show_product_names = False

    def __init__(self, volume, conc=None, products=None):
        self.volume = volume
        self.conc = conc
        self.products = products or []

    @classmethod
    def make(cls, db, products):
        makers_it = iter_combo_makers(
                cls.from_params,
                map(cls.from_product, products),
                group_by={
                    'volume': group_by_identity,
                    'conc': group_by_identity,
                },
                merge_by={
                },
        )
        makers = list(makers_it)
        show_product_names = (len(makers) != 1)

        for maker in makers:
            maker.show_product_names = show_product_names
            yield maker

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
