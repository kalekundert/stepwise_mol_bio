#!/usr/bin/env python3

"""\
Perform a serial dilution.

Usage:
    serial_dilution <volume> <high> <low> <steps> [options]

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
"""

import stepwise
import autoprop
from inform import Error, plural
from numbers import Real
from stepwise_mol_bio import Main

@autoprop
class SerialDilution(Main):

    def __init__(self):
        self.volume = None
        self.volume_unit = None
        self.conc_high = None
        self.conc_low = None
        self.conc_unit = None
        self.steps = None

        self.material = 'material'
        self.diluent = 'water'

    @classmethod
    def from_docopt(cls, args):
        self = cls()

        self.volume, self.volume_unit = parse_quantity(args['<volume>'])
        self.conc_high, self.conc_low, self.conc_unit = parse_high_low(
                args['<high>'],
                args['<low>'],
        )
        self.steps = int(args['<steps>'])

        self.material = args['--material']
        self.diluent = args['--diluent']

        return self

    def get_protocol(self):
        factor = self.dilution_factor
        transfer = self.volume * factor / (1 - factor)
        initial_volume = self.volume + transfer

        conc_high_str = format_quantity(self.conc_high, self.conc_unit)
        material_str = f'{conc_high_str} {self.material}'.lstrip()
        conc_table = [
                [i, format_quantity(conc, self.conc_unit)]
                for i, conc in enumerate(self.concentrations, 1)
        ]

        protocol = stepwise.Protocol()

        protocol += f"""\
Perform a serial dilution [1]:

- Put {initial_volume:.2f} μL {material_str} in the first tube.
- Add {self.volume:.2f} μL {self.diluent} in the {plural(self.steps - 1):# remaining tube/s}.
- Transfer {transfer:.2f} μL between each tube.
"""

        protocol.footnotes[1] = f"""\
The final concentrations will be:
{stepwise.tabulate(conc_table, alignments='>>')}
"""
        return protocol

    def get_dilution_factor(self):
        return (self.conc_low / self.conc_high)**(1 / (self.steps - 1))

    def get_concentrations(self):
        factor = self.dilution_factor
        return [
                self.conc_high * factor**i
                for i in range(self.steps)
        ]

def parse_quantity(x):
    try:
        return stepwise.Quantity.from_string(x).tuple
    except:
        return float(x), None

def parse_high_low(high_str, low_str):
    high, high_unit = parse_quantity(high_str)
    low, low_unit = parse_quantity(low_str)

    units = {high_unit, low_unit}
    units.discard(None)

    if len(units) > 1:
        raise Error(f"units don't match: {high_unit!r}, {low_unit!r}")

    return high, low, units.pop() if units else None

def format_quantity(value, unit=None, precision=2, pad=' '):
    if isinstance(value, Real):
        value = f'{value:.{precision}f}'
    return f'{value}{pad}{unit}' if unit else value


if __name__ == '__main__':
    SerialDilution.main(__doc__)

# vim: tw=53
