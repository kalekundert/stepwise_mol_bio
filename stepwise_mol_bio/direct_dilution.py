#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import Error
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
    direct_dilution <volume> <high> to <low> <steps> [options]
    direct_dilution <volume> <high> / <factor> <steps> [options]
    direct_dilution <volume> <low> x <factor> <steps> [options]

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

    <factor>
        How big of a dilution to make at each step of the protocol.

    <steps>
        The number of dilutions to make, including <high> and <low>.

Options:
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
    max_dilution = appcli.param('--max-dilution', cast=float, default=10)

    def get_protocol(self):
        header = [
                format_quantity('Final', f'[{self.conc_unit}]', pad='\n'),
                format_quantity('Stock', f'[{self.conc_unit}]', pad='\n'),
                format_quantity(self.material, f'[{self.volume_unit}]', pad='\n'),
                format_quantity(self.diluent, f'[{self.volume_unit}]', pad='\n'),
        ]
        rows = []
        volumes = {}
        stock_concs = self._pick_stock_concs()

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

        protocol = stepwise.Protocol()
        protocol += f"""\
Prepare the following dilutions:

{stepwise.tabulate(rows, header, align='>>>>')}
"""
        return protocol

    def _pick_stock_concs(self):
        stock_concs = {}
        stock_conc = self.conc_high
        prev_target_conc = None

        for target_conc in self.concentrations:
            if target_conc == 0:
                stock_concs[target_conc] = stock_conc
                break

            dilution = stock_conc / target_conc
            if dilution > self.max_dilution:
                stock_conc = prev_target_conc

            dilution = stock_conc / target_conc
            if dilution > self.max_dilution:
                raise Error(f"{dilution:.1g}x dilution to make {format_quantity(target_conc, self.conc_unit)} exceeds maximum ({self.max_dilution:.1g}x)")

            stock_concs[target_conc] = stock_conc
            prev_target_conc = target_conc

        return stock_concs

if __name__ == '__main__':
    DirectDilution.main()

# vim: tw=53
