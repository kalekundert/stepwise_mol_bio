#!/usr/bin/env python3

import autoprop
import byoc

from stepwise import Protocol, Quantity, StepwiseConfig
from stepwise_mol_bio import Cleanup
from freezerbox import ProductConfig
from byoc import Key, DocoptConfig

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
            DocoptConfig,
            ProductConfig,
            StepwiseConfig.setup(('molbio', 'miniprep')),
    ]

    product_ori = byoc.param(
            Key(ProductConfig, 'origin'),
    )
    expected_yield_ng_uL = byoc.param(
            Key(StepwiseConfig, 'yield_ng_uL'),
    )

    @classmethod
    def _combo_maker_factory(cls, db):
        # `ProductConfig` is automatically loaded for the solo makers, but not 
        # for the combo makers.  I feel like the `ProductConfig` approach just 
        # isn't the right way to do this, though.  I think the approach I 
        # normally take is to:
        #
        # - Make a nested class representing the plasmid
        # - Use `ReagentConfig` within that class to get the info I need.
        # - If product conc is not ambiguous (i.e. if there is only one 
        #   plasmid, or all the plasmids have the same ORI), return it.
        #
        # Loading the `ProductConfig` here works around the problem in the 
        # short term, though.
        app = super()._combo_maker_factory(db)
        app.load(ProductConfig)
        return app

    def get_protocol(self):
        p = Protocol()
        p += "Miniprep."
        return p

    def get_product_conc(self):
        yield_ng_uL = self.expected_yield_ng_uL

        if isinstance(yield_ng_uL, dict):
            key = getattr(self, 'product_ori', 'default')
            yield_ng_uL = yield_ng_uL.get(key, 'default')

        return Quantity(yield_ng_uL, 'ng/µL')

if __name__ == '__main__':
    Miniprep.main()
