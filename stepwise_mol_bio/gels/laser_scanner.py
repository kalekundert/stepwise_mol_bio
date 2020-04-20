#!/usr/bin/env python3

"""\
Image a gel using a laser scanner.

Usage:
    laser_scanner <optics>...

Arguments:
    <optics>
        A laser/filter combination to use.  The most convenient way to specify 
        such a combination is to give the name of a preset.  The following 
        presets are currently available: 

        {presets}

        You can define new presets by adding blocks like the following to your 
        stepwise configuration file:

            [molbio.laser.presets.name]
            laser = <int> or <list of ints>
            filter = <str> or <list of strs>

        You can also explicitly provide laser and filter parameters using the 
        following syntax:

            <laser>/<filter>
"""

import stepwise
import autoprop
from inform import indent, plural
from stepwise_mol_bio import Main, Presets, ConfigError, UsageError

PRESETS = Presets.from_config('molbio.laser.presets')
PRESETS_DOC = PRESETS.format_briefs('{laser} nm')
__doc__ = __doc__.format(
        presets=indent(PRESETS_DOC, 8*' ', first=-1),
)

@autoprop
class LaserScanner(Main):

    def __init__(self, preset=None):
        self.optics = []

    @classmethod
    def from_docopt(cls, args):
        self = cls()
        self.optics = args['<optics>']
        return self

    @classmethod
    def from_params(cls, laser, filter):
        self = cls()
        self.optics.append({
            'laser': laser,
            'filter': filter,
        })
        return self

    def get_protocol(self):
        optics = [parse_optics(x) for x in self.optics]
        lasers = [f"{plural(optics):laser/s}:"] + [f"{x['laser']} nm" for x in optics]
        filters = [f"{plural(optics):filter/s}:"] + [x['filter'] for x in optics]

        p = stepwise.Protocol()
        p += f"""\
Image with a laser scanner:

{stepwise.tabulate([lasers, filters])}
"""
        return p

def parse_optics(optics):
    if isinstance(optics, dict):
        return optics

    try:
        return PRESETS.load(optics)
    except ConfigError:
        try: 
            laser, filter = optics.split('/')
            return {'laser': laser, 'filter': filter}
        except ValueError as err:
            raise UsageError(f"expected a preset or '<laser>/<filter>', got {optics!r}") from None


if __name__ == '__main__':
    LaserScanner.main(__doc__)
