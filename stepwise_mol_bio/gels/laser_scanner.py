#!/usr/bin/env python3

"""\
Image a gel using a laser scanner.

Usage:
    laser_scanner [<preset>] [-l <laser>] [-f <filter>]

Arguments:
    <preset>
        The name of a predefined laser/filter combination.  The following 
        presets are available:

        {presets}

Options:
    -l --laser <nm>
        The wavelength of laser to use.

    -f --filter <desc>
        The emission filter to use.
"""

import stepwise
import autoprop
from inform import indent
from stepwise_mol_bio import Main, Presets, ConfigError

PRESETS = Presets.from_config('molbio.laser.presets')
PRESETS_DOC = PRESETS.format_briefs('{laser} nm')
__doc__ = __doc__.format(
        presets=indent(PRESETS_DOC, 8*' ', first=-1),
)

@autoprop
class LaserScanner(Main):

    def __init__(self, preset=None):
        self.preset = preset
        self.params = {}

    @classmethod
    def from_docopt(cls, args):
        self = cls()
        self.preset = args['<preset>']

        if x := args['--laser']:
            self.params['laser'] = x
        if x := args['--filter']:
            self.params['filter'] = x

        return self

    @classmethod
    def from_params(cls, laser, filter):
        self = cls()
        self.params['laser'] = laser
        self.params['filter'] = filter
        return self

    def get_config(self):
        preset = PRESETS.load(self.preset) if self.preset else {}
        return {**preset, **self.params}

    def get_protocol(self):
        p = stepwise.Protocol()
        c = self.config

        try:
            p += f"""\
Image with a laser scanner:

laser: {c['laser']} nm
filter: {c['filter']}
"""
        except KeyError as err:
            raise ConfigError(f"no {err} specified") from None

        return p

if __name__ == '__main__':
    LaserScanner.main(__doc__)
