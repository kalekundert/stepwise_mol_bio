#!/usr/bin/env python3

"""\
Dilute the given reagent in such a way that serial dilutions are avoided.

While serial dilutions are easy to perform, they have some disadvantages: 
First, they can be inaccurate, because errors made at any step are propagated 
to all other steps.  Second, they can waste material, because they produce more 
of the lowest dilution than is needed.  If either of these disadvantages are a 
concern, this protocol can be used to create the same dilution with fewer 
serial steps.

Usage:
    direct_dilution <volume> <high> <low> <steps> [options]

Arguments:
    <volume>
        The volume of reagent needed for each concentration.  A unit may be 
        optionally given, in which case it will be included in the protocol.

    <high>
        The starting concentration for the dilution.  A unit may be optionally 
        given, in which case it will be included in the protocol.

    <low>
        The ending concentration for the dilution.  A unit may be optionally 
        given, in which case it will be included in the protocol.

    <steps>
        The number of dilutions to make, including <high> and <low>.

Options:
    -m --material NAME  [default: material]
        The substance being diluted.

    -d --diluent NAME   [default: water]
        The substance to dilute into.

    -x --max-dilution FOLD  [default: 10]
        Specify the biggest dilution that can be made at any step, as larger 
        dilutions are prone to be less accurate.
"""

import stepwise
import autoprop
from inform import Error
from stepwise_mol_bio.serial_dilution import SerialDilution, format_quantity

@autoprop
class DirectDilution(SerialDilution):

    def __init__(self):
        super().__init__()
        self.max_dilution = 10

    @classmethod
    def from_docopt(cls, args):
        self = super().from_docopt(args)
        self.max_dilution = float(args['--max-dilution'])
        return self

    def get_protocol(self):
        header = [
                format_quantity('Final', f'[{self.conc_unit}]', pad='\n'),
                format_quantity('Stock', f'[{self.conc_unit}]', pad='\n'),
                format_quantity(self.material, f'[{self.volume_unit}]', pad='\n'),
                format_quantity(self.diluent, f'[{self.volume_unit}]', pad='\n'),
        ]
        rows = []

        volumes = {}
        for target_conc in reversed(self.concentrations):
            target_volume = self.volume + volumes.get(target_conc, 0)
            stock_conc = self._pick_stock_conc(target_conc)
            stock_volume = target_volume * target_conc / stock_conc

            volumes[stock_conc] = volumes.get(stock_conc, 0) + stock_volume
            rows.insert(0, [
                format_quantity(target_conc),
                format_quantity(stock_conc),
                format_quantity(stock_volume),
                format_quantity(target_volume - stock_volume),
            ])

        protocol = stepwise.Protocol()
        protocol += f"""\
Prepare the following dilutions:

{stepwise.tabulate(rows, header, alignments='>>>>')}
"""
        return protocol

    def _pick_stock_conc(self, target_conc):
        stock_concs = {x for x in self.concentrations if x > target_conc}
        stock_concs.add(self.conc_high)

        for stock_conc in sorted(stock_concs, reverse=True):
            dilution = stock_conc / target_conc
            if dilution <= self.max_dilution:
                return stock_conc

        raise Error(f"{dilution:.1f}x dilution to make {format_quantity(target_conc, self.conc_unit)} exceeds maximum ({self.max_dilution:.1f}x).")

if __name__ == '__main__':
    DirectDilution.main(__doc__)

# vim: tw=53
