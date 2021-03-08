#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import indent, plural
from appcli import DocoptConfig
from stepwise import StepwiseConfig, PresetConfig, UsageError, pl, pre
from stepwise_mol_bio import Main, ConfigError, merge_dicts

@autoprop
class LaserScanner(Main):
    """\
Image a gel using a laser scanner.

Usage:
    laser_scanner <optics>...

<%! from stepwise_mol_bio import hanging_indent %>\
Arguments:
    <optics>
        A laser/filter combination to use.  The most convenient way to specify 
        such a combination is to give the name of a preset.  The following 
        presets are currently available: 

        ${hanging_indent(app.preset_briefs, 8*' ')}

        You can define new presets by adding blocks like the following to your 
        stepwise configuration file:

            [molbio.laser.presets.name]
            laser = <int> or <list of ints>
            filter = <str> or <list of strs>

        You can also explicitly provide laser and filter parameters using the 
        following syntax:

            <laser>/<filter>
"""
    __config__ = [
            DocoptConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.laser'),
    ]

    presets = appcli.param(
            appcli.Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset_briefs = appcli.config_attr()
    preset_brief_template = '{laser} nm'

    # This attribute is a list of:
    # - name of preset (str)
    # - laser/filter (str)
    # - {laser: filter} (dict)
    optics = appcli.param(
            appcli.Key(DocoptConfig, '<optics>'),
            default_factory=list,
    )

    def __init__(self, *optics):
        self.optics = list(optics)

    @classmethod
    def from_laser_filter_pair(cls, laser, filter):
        self = cls()
        self.optics.append({
            'laser': laser,
            'filter': filter,
        })
        return self

    def get_protocol(self):
        optics = [self.parse_optics(x) for x in self.optics]
        lasers = [f"{plural(optics):laser/s}:"] + [f"{x['laser']} nm" for x in optics]
        filters = [f"{plural(optics):filter/s}:"] + [x['filter'] for x in optics]

        p = stepwise.Protocol()
        p += pl(
                "Image with a laser scanner:",
                pre(stepwise.tabulate([lasers, filters])),
        )
        return p

    def parse_optics(self, optics):
        if isinstance(optics, dict):
            return optics

        try:
            return self.presets[optics]
        except ConfigError:
            try: 
                laser, filter = optics.split('/')
                return {'laser': laser, 'filter': filter}
            except ValueError as err:
                raise UsageError(f"expected a preset or '<laser>/<filter>', got {optics!r}") from None

if __name__ == '__main__':
    LaserScanner.main()
