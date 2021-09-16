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
        The primary wavelength of light produced by the transilluminator, in 
        nm.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.transilluminator.presets:
        Named groups of default parameters.  See below for the various settings 
        that can be specified in each preset.

    molbio.transilluminator.presets.<name>.colloquial_name:
        How to refer to the transilluminator in the protocol, e.g. "UV" or 
        "blue-light".

    molbio.transilluminator.presets.<name>.default_wavelength_nm:
        The default value for the `--wavelength` option.
"""
    __config__ = [
            DocoptConfig,
            PresetConfig,
            StepwiseConfig.setup('molbio.transilluminator'),
    ]
    config_paths = appcli.config_attr()
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
