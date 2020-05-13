#!/usr/bin/env python3

"""\
Perform a serial dilution.

Usage:
    serial_dilution <volume> <high> <low> <steps> [-m <name>] [-d <name>]
    serial_dilution <volume> <high> <steps> -f <factor> [-m <name>] [-d <name>]

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
    -f --factor <x>
        How big of a dilution to make at each step of the protocol.  With this 
        flag, the <low> argument no longer needs to be specified.

    -m --material <name>    [default: material]
        The substance being diluted.

    -d --diluent <name>     [default: water]
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
        self._conc_low = None
        self.conc_unit = None
        self._dilution_factor = None
        self.steps = None

        self.material = 'material'
        self.diluent = 'water'

    @classmethod
    def from_docopt(cls, args):
        self = cls()

        self.volume, self.volume_unit = parse_quantity(args['<volume>'])
        self.steps = int(args['<steps>'])

        if x := args['--factor']:
            self.conc_high, self.conc_unit = parse_quantity(args['<high>'])
            self.dilution_factor = 1 / float(x)
        else:
            self.conc_high, self.conc_low, self.conc_unit = parse_high_low(
                    args['<high>'],
                    args['<low>'],
            )

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
- Transfer {transfer:.2f} μL between each tube to make
  {self.steps} {1/factor:.2g}-fold dilutions.
"""

        protocol.footnotes[1] = f"""\
The final concentrations will be:
{stepwise.tabulate(conc_table, align='>>')}
"""
        return protocol

    def get_conc_low(self):
        return self._conc_low or self.conc_high * self.factor**(self.steps)

    def set_conc_low(self, x):
        self._conc_low = x
        self._dilution_factor = None

    def get_dilution_factor(self):
        return self._dilution_factor or \
                (self.conc_low / self.conc_high)**(1 / (self.steps - 1))

    def set_dilution_factor(self, x):
        self._dilution_factor = x
        self._conc_low = None

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
