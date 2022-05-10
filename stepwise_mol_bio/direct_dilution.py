#!/usr/bin/env python3

import stepwise, byoc, autoprop

from stepwise import pl, table
from stepwise_mol_bio import UsageError
from stepwise_mol_bio.serial_dilution import SerialDilution, format_quantity

@autoprop
class DirectDilution(SerialDilution):
    """\
Dilute the given reagent in such a way that serial dilutions are minimized.

While serial dilutions are easy to perform, they have some disadvantages: 
First, they can be inaccurate, because errors made at any step are propagated 
to all other steps.  Second, they can waste material, because they produce more 
of the lowest dilution than is needed.  If either of these disadvantages are a 
concern, use this protocol to create the same dilution with fewer serial steps.

Usage:
    direct_dilution (<high> to <low> | <high> / <factor> | <low> x <factor>)
         (-n <int>) (-v <volume>) [-m <material>] [-d <diluent>] [-x <fold>] 
         [-0]

Arguments:
    <high>
        The starting concentration for the dilution.  A unit may be optionally 
        given, in which case it will be included in the protocol.

    <low>
        The ending concentration for the dilution.  A unit may be optionally 
        given, in which case it will be included in the protocol.

    <factor>
        How big of a dilution to make at each step of the protocol.

Options:
    -n --num-dilutions <int>
        The number of dilutions to make, including <high> and <low>.

    -v --volume <volume>
        The volume of reagent needed for each concentration.  A unit may be 
        optionally given, in which case it will be included in the protocol.

    -m --material <name>        [default: ${app.material}]
        The substance being diluted.

    -d --diluent <name>         [default: ${app.diluent}]
        The substance to dilute into.

    -x --max-dilution <fold>    [default: ${app.max_dilution}]
        Specify the biggest dilution that can be made at any step, as larger 
        dilutions are prone to be less accurate.

    -0 --include-zero
        Include a "dilution" with no material in the protocol.
"""
    max_dilution = byoc.param('--max-dilution', cast=float, default=10)

    def get_protocol(self):
        protocol = stepwise.Protocol()
        protocol += pl(
                "Prepare the following dilutions:",
                self.dilution_table,
        )
        return protocol

    def get_dilution_table(self):
        header = [
                format_quantity('Final', self.conc_unit and f'[{self.conc_unit}]', pad='\n'),
                format_quantity('Stock', self.conc_unit and f'[{self.conc_unit}]', pad='\n'),
                format_quantity(self.material, f'[{self.volume_unit}]', pad='\n'),
                format_quantity(self.diluent, f'[{self.volume_unit}]', pad='\n'),
        ]
        rows = []
        volumes = {}
        stock_concs = pick_stock_concs(
                self.concentrations,
                self.max_dilution,
                self.conc_unit,
        )

        for target_conc in reversed(self.concentrations):
            target_volume = self.volume + volumes.get(target_conc, 0)
            stock_conc = stock_concs[target_conc]
            stock_volume = target_volume * target_conc / stock_conc

            volumes[stock_conc] = volumes.get(stock_conc, 0) + stock_volume
            rows.insert(0, [
                format_quantity(target_conc),
                format_quantity(stock_conc),
                format_quantity(stock_volume),
                format_quantity(target_volume - stock_volume),
            ])

        return table(rows, header, align='>>>>')

def pick_stock_concs(target_concs, max_dilution, conc_unit=None):
    """
    Return a dictionary indicating which stock concentration should be used 
    for each dilution.  

    More specifically, the keys are the target concentrations and the 
    values are the stock concentrations.  There will be a key for each 
    concentration associated with this dilution.  If the dilution includes 
    a "0" concentration sample, it will also be included in the dictionary 
    (albeit with an arbitrary value).
    """
    stock_concs = {}
    stock_conc = max(target_concs)
    prev_target_conc = None

    for target_conc in target_concs:
        if target_conc == 0:
            stock_concs[target_conc] = stock_conc
            break

        dilution = stock_conc / target_conc
        if dilution > max_dilution:
            stock_conc = prev_target_conc

            dilution = stock_conc / target_conc
            if dilution > max_dilution:
                err = UsageError(
                        target_conc=target_conc,
                        conc_unit=conc_unit,
                        dilution=dilution,
                        max_dilution=max_dilution,
                )
                err.brief = lambda e: f"{e.dilution:.1g}x dilution to make {format_quantity(e.target_conc, e.conc_unit)} exceeds maximum ({e.max_dilution:.1g}x)"
                raise err

        stock_concs[target_conc] = stock_conc
        prev_target_conc = target_conc

    return stock_concs

if __name__ == '__main__':
    DirectDilution.main()
