#!/usr/bin/env python3

import stepwise, appcli, autoprop

from stepwise import StepwiseConfig, PresetConfig, MasterMix, pl, ul
from stepwise_mol_bio import Main
from appcli import Key, DocoptConfig
from inform import plural
from copy import deepcopy

@autoprop.cache
class Kld(Main):
    """\
Circularize a linear DNA molecule using T4 DNA ligase, e.g. to reform a plasmid 
after inverse PCR.

Usage:
    kld <dna> [-p <preset>] [-n <int>]

Arguments:
    <dna>
        The name of the DNA molecule to circularize.

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p -preset <name>
        What default reaction parameters to use.  The following parameters are 
        currently available:

        ${hanging_indent(app.preset_briefs, 8)}

    -n --num-reactions <int>        [default: ${app.num_reactions}]
        The number of reactions to set up.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.kld.default_preset:
        The default value for the `--preset` option.

    molbio.kld.presets.<name>.reaction:
        A table detailing all the components of the transcription reaction, in 
        the format understood by `stepwise.MasterMix.from_string()`.  The 
        reaction must have one reagent named "DNA"; this will be replaced with 
        the name(s) specified on the command line.
"""
    __config__ = [
            DocoptConfig,
            PresetConfig,
            StepwiseConfig.setup('molbio.kld'),
    ]
    preset_briefs = appcli.config_attr()
    config_paths = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=list,
    )
    preset = appcli.param(
            Key(DocoptConfig, '--preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    base_reaction = appcli.param(
            Key(PresetConfig, 'reaction', cast=MasterMix),
    )
    dna = appcli.param(
            Key(DocoptConfig, '<dna>'),
    )
    num_reactions = appcli.param(
            Key(DocoptConfig, '--num-reactions'),
            cast=lambda x: int(eval(x)),
            default=1,
    )

    def __init__(self, dna):
        self.dna = dna

    def get_reaction(self):
        kld = deepcopy(self.base_reaction)
        kld.num_reactions = self.num_reactions
        kld.extra_percent = 15
        kld['DNA'].name = self.dna

        return kld

    def get_protocol(self):
        p = stepwise.Protocol()
        p += pl(
                f"Setup {plural(self.num_reactions):# ligation reaction/s}:",
                self.reaction,
        )
        p += "Incubate at room temperature for 1h."
        return p

if __name__ == '__main__':
    Kld.main()
