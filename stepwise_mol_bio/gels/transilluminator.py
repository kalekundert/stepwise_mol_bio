#!/usr/bin/env python3

import stepwise, appcli, autoprop

from stepwise import StepwiseConfig, PresetConfig
from stepwise_mol_bio import Main
from appcli import Key, DocoptConfig

@autoprop
class Transilluminator(Main):
    """\
Image a gel using a transilluminator.

Usage:
    transilluminator [<preset>] [-w <nm>]

<%! from stepwise_mol_bio import hanging_indent %>\
Arguments:
    <preset>
        The default set of parameters to use.  The following presets are 
        available:

        ${hanging_indent(app.preset_briefs, 8*' ')}

Options:
    -w --wavelength <nm>
        The specific wavelength of light to use, in nm.
"""
    __config__ = [
            DocoptConfig,
            PresetConfig,
            StepwiseConfig.setup('molbio.transilluminator'),
    ]
    preset_briefs = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=list,
    )
    preset = appcli.param(
            Key(DocoptConfig, '<preset>'),
            Key(StepwiseConfig, 'default_preset'),
    )
    colloquial_name = appcli.param(
            Key(PresetConfig, 'colloquial_name'),
    )
    wavelength_nm = appcli.param(
            Key(DocoptConfig, '--wavelength'),
            Key(PresetConfig, 'default_wavelength_nm'),
    )

    def get_protocol(self):
        phrases = ["Image with a"]
        if self.wavelength_nm:
            phrases += [f"{self.wavelength_nm} nm"]
        phrases += [f"{self.colloquial_name} transilluminator."]

        p = stepwise.Protocol()
        p += ' '.join(phrases)
        return p

if __name__ == '__main__':
    Transilluminator.main()
