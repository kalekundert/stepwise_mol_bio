#!/usr/bin/env python3

import stepwise, byoc, autoprop
from byoc import DocoptConfig
from inform import plural
from stepwise import pl, ul, pre
from stepwise_mol_bio import Main, UsageError

@autoprop
class SerialDilution(Main):
    """\
Perform a serial dilution.

Usage:
    serial_dilution (<high> to <low> | <high> / <factor> | <low> x <factor>)
         (-n <int>) (-v <volume>) [-m <material>] [-d <diluent>] [-0]

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

    -m --material <name>    [default: ${app.material}]
        The substance being diluted.

    -d --diluent <name>     [default: ${app.diluent}]
        The substance to dilute into.

    -0 --include-zero
        Include a "dilution" with no material in the protocol.
"""
    __config__ = [
            DocoptConfig,
    ]

    _conc_high_str = byoc.param('<high>', default=None)
    _conc_low_str = byoc.param('<low>', default=None)
    _factor = byoc.param('<factor>', cast=float, default=None)
    _volume_str = byoc.param('--volume')

    num_dilutions = byoc.param('--num-dilutions', cast=int)
    material = byoc.param('--material', default='material')
    diluent = byoc.param('--diluent', default='water')
    include_zero = byoc.param('--include-zero', default=False)

    def __bareinit__(self):
        self._volume = None
        self._volume_unit = None

        self._conc_high = None
        self._conc_low = None
        self._conc_unit = None

    def __init__(self, volume, num_dilutions):
        self.volume = volume
        self.num_dilutions = num_dilutions

    def get_protocol(self):
        transfer = self.volume * (f := self.inv_factor) / (1 - f)
        initial_volume = self.volume + transfer

        conc_high_str = format_quantity(self.conc_high, self.conc_unit, 'g')
        material_str = f'{conc_high_str} {self.material}'.lstrip()
        conc_table = [
                [i, format_quantity(conc, self.conc_unit, 'e')]
                for i, conc in enumerate(self.concentrations, 1)
        ]

        num_tubes = self.num_dilutions if self.include_zero else self.num_dilutions - 1
        each_tube = 'each tube *except the last*' if self.include_zero else 'each tube'

        protocol = stepwise.Protocol()
        protocol += pl(
                "Perform a serial dilution [1]:",
                ul(
                    f"Put {initial_volume:.2f} {self.volume_unit} {material_str} in a tube.",
                    f"Put {self.volume:.2f} {self.volume_unit} {self.diluent} in {plural(num_tubes):# adjacent tube/s}.",
                    f"Transfer {transfer:.2f} {self.volume_unit} between {each_tube} to make {self.num_dilutions} {self.factor:.2g}-fold dilutions.",
                ),
        )

        protocol.footnotes[1] = pl(
                "The final concentrations will be:",
                pre(stepwise.tabulate(conc_table, align='>>')),
                br='\n',
        )
        return protocol

    def get_concentrations(self):
        concs = [
                self.conc_high * self.factor**-i
                for i in range(self.num_dilutions)
        ]

        if self.include_zero:
            concs.append(0)

        return concs

    def get_volume(self):
        return self._volume

    def get_volume_unit(self):
        return self._volume_unit or 'µL'

    def set_volume(self, volume):
        self._volume, self._volume_unit = parse_quantity(volume)

    def get_conc_high(self):
        return self._conc_high

    def get_conc_low(self):
        return self._conc_low

    def get_conc_unit(self):
        return self._conc_unit

    def set_conc_high_low(self, high, low):
        self._check_num_dilutions()
        self._conc_high, self._conc_low, self._conc_unit = parse_high_low(high, low)
        self._factor = (self._conc_high / self._conc_low)**(1 / (self.num_dilutions - 1))

    def set_conc_high_factor(self, high, factor):
        self._check_num_dilutions()
        self._conc_high, self._conc_unit = parse_quantity(high)
        self._conc_low = self._conc_high * factor**-(self.num_dilutions - 1)
        self._factor = factor

    def set_conc_low_factor(self, low, factor):
        self._check_num_dilutions()
        self._conc_low, self._conc_unit = parse_quantity(low)
        self._conc_high = self._conc_low * factor**(self.num_dilutions - 1)
        self._factor = factor

    def get_factor(self):
        return self._factor

    def get_inv_factor(self):
        return 1 / self._factor

    @byoc.on_load(DocoptConfig)
    def _load_volume(self):
        self.volume = self._volume_str

    @byoc.on_load(DocoptConfig)
    def _load_conc(self):
        if self._conc_high_str and self._conc_low_str:
            self.set_conc_high_low(
                    self._conc_high_str,
                    self._conc_low_str,
            )

        elif self._conc_high_str and self._factor:
            self.set_conc_high_factor(
                    self._conc_high_str,
                    self._factor,
            )

        elif self._conc_low_str and self._factor:
            self.set_conc_low_factor(
                    self._conc_low_str,
                    self._factor,
            )

    def _check_num_dilutions(self):
        if self.num_dilutions < 2:
            err = UsageError(n=self.num_dilutions)
            err.brief = "must request at least 2 dilutions: <low> and <high>"
            err.blame += lambda e: f"only {plural(e.n):# dilution/s} requested"
            raise err

def parse_quantity(x):
    try:
        return x.value, x.unit
    except:
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
        raise UsageError(
                "units don't match: {high_unit!r}, {low_unit!r}",
                high_unit=high_unit,
                low_unit=low_unit,
        )

    return high, low, units.pop() if units else None

def format_quantity(value, unit=None, format='.2f', pad=' '):
    from numbers import Real
    if isinstance(value, Real):
        value = f'{value:{format}}'
    return f'{value}{pad}{unit}' if unit else value


if __name__ == '__main__':
    SerialDilution.main()

# vim: tw=53
