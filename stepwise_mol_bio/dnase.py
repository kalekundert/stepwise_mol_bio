#!/usr/bin/env python3

# How to handle RNA stock conc/volume?
# - Set in preset, but don't really care
# - RT needs to set without holding, but normal users probably want to hold.
# - On the other hand, the total quantity almost doesn't matter, so maybe I 
#   should just never hold.
#
#
# - Add stock conc/volume options
# - Add 'hold' option: 'conc' or 'volume'
# - Really, want a way to say don't care as long as min./max satisfied.  But 
#   not obvious how to do that.
#
# Change reaction volume
# - Add volume option
# - Additives:
#   - Need to parse to get volume/stock conc
#   - Add function for this to stepwise

import stepwise
import byoc
import autoprop

from stepwise import (
        StepwiseConfig, PresetConfig, Reactions, Reaction, pl, ul
)
from stepwise_mol_bio import (
        App, Bindable, Thermocycler, UsageError, bind, group_samples
)
from stepwise_mol_bio.thermocycler import format_thermocycler_steps
from byoc import Key, DocoptConfig

def dnase_digest(samples):
    return stepwise.Protocol.merge(
            *plan_dnase_protocol.for_group_in(samples)
    )

@group_samples(
        'reaction_prototype',
        'rna_volume',
        'incubation',
        'denature_additives',
        'denature_incubation',
)
def plan_dnase_protocol(group):
    p = stepwise.Protocol()
    p += plan_dnase_reactions(group)
    p += plan_incubation_steps(group)
    return p

@group_samples(
        'reaction_prototype',
        'rna_volume',
)
def plan_dnase_reactions(group):
    rxn = group.reaction_prototype.copy()

    if 'RNA' not in rxn:
        raise UsageError("reaction missing required 'RNA' reagent")

    # The concentration of RNA in this reaction doesn't really matter, since 
    # it's the trace DNA (in unknown quantity) that's being acted on.  The 
    # protocols I've found just call for the quantity of RNA to be between 1 pg 
    # and 10 µg of RNA, which is a huge range.
    #
    # Normally I would hold concentration while changing volume or stock 
    # concentration, but in this case that wouldn't make sense.  In fact, if 
    # the samples we're given have different stock concentrations, I just don't 
    # show them.

    if group.rna_volume:
        rxn['RNA'].volume = group.rna_volume, rxn.volume.unit

    stock_concs = {x.rna_stock_conc for x in group}
    stock_concs.discard(None)

    if stock_concs:
        if len(stock_concs) == 1:
            rxn['RNA'].stock_conc = stock_concs.pop()
        else:
            del rxn['RNA'].stock_conc

    combos = [
            {'RNA': x.rna_name}
            for x in group
    ]

    rxns = Reactions(rxn, combos)
    rxns.step_kind = 'DNase'
    rxns.split_replicates = False
    return rxns

@group_samples(
        'incubation',
        'denature_additives',
        'denature_incubation',
)
def plan_incubation_steps(group):
    if not group.denature_additives:
        if group.denature_incubation:
            return Thermocycler([
                group.incubation,
                group.denature_incubation,
            ])
        else:
            return format_thermocycler_steps(
                    group.incubation,
                    incubate_prefix=True,
            )
    else:
        p = stepwise.Protocol()
        p += format_thermocycler_steps(
                group.incubation,
                incubate_prefix=True,
        )

        if len(group.denature_additives) == 1:
            steps = ul(f"Add {group.denature_additives[0]}.")
        else:
            steps = ul(
                    pl(
                        "Add the following:",
                        ul(*group.denature_additives),
                        br='\n',
                    ),
                    br='\n\n',
            )

        steps += format_thermocycler_steps(
                group.denature_incubation,
                incubate_prefix=True,
        )
        p += pl("Denature the DNase:", steps)

        return p

def samples_from_docopt(args):
    return [Dnase.Sample(x) for x in args['<samples>']]

@autoprop
class Dnase(App):
    """\
Remove trace DNA contamination from RNA samples.

Usage:
    dnase <samples>... [-p <preset>]

Arguments:
    <samples>
        The names of the RNA samples to treat with DNase.

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>
        What set of default reaction parameters to use.  The following presets 
        are currently available:

        ${hanging_indent(sample.preset_briefs, 8)}

    -V --rna-volume <µL>
        The volume of RNA to digest.  This option does not change the stock 
        concentration of the RNA, and therefore directly affects the 
        concentration of RNA in the final reaction.

    -C --rna-stock <conc>
        The stock concentration of the RNA to digest.  This option does not 
        change the volume of the RNA, and therefore directly affects the 
        concentration of RNA in the final reaction.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(sample.config_paths, 8)}

    molbio.dnase.default_preset
        The default value for the `--preset` option.

    molbio.dnase.presets:
        Named groups of default reaction parameters.  Typically each preset 
        corresponds to a particular kit or protocol.  See below for the various 
        settings that can be specified in each preset.

    molbio.dnase.presets.<name>.brief:
        A brief description of the preset.  This is displayed in the usage info 
        for the `--preset` option.

    molbio.dnase.presets.<name>.inherit:
        Copy all settings from another preset.  This can be used to make small 
        tweaks to a protocol, e.g. "DNase I with a non-standard additive".

    molbio.dnase.presets.<name>.reaction:
        A table detailing all of the components of the DNase reaction.  This 
        table must contain 1 reagent named "RNA".

    molbio.dnase.presets.<name>.incubation:
        The time and temperature to incubate the digestion reaction, in the 
        format expected by `format_thermocycler_steps`, e.g. a mapping with 
        "time_m" and "temp_C" keys.

    molbio.dnase.presets.<name>.denature_additives:
        A list of reagents to add to the reaction before performing the 
        denaturation step, e.g. "1 µL 500 mM EDTA".

    molbio.dnase.presets.<name>.denature_incubation:
        The same as the `incubation` setting, but for the denaturation step.
"""

    class Sample(Bindable, use_app_configs=True):

        def __init__(self, rna_name, **kwargs):
            super().__init__(**kwargs)
            self.rna_name = rna_name

        __config__ = [
                PresetConfig,
                StepwiseConfig.setup(('molbio', 'dnase')),
        ]
        presets = byoc.param(
                Key(StepwiseConfig, 'presets'),
                pick=list,
        )
        preset = byoc.param(
                Key(DocoptConfig, '--preset'),
                Key(StepwiseConfig, 'default_preset'),
        )
        reaction_prototype = byoc.param(
                Key(PresetConfig, 'reaction', cast=Reaction.from_text),
        )
        rna_volume = byoc.param(
                Key(DocoptConfig, '--rna-volume'),
                default=None,
        )
        rna_stock_conc = byoc.param(
                Key(DocoptConfig, '--rna-stock'),
                default=None,
        )
        incubation = byoc.param(
                Key(PresetConfig, 'incubation'),
        )
        denature_additives = byoc.param(
                Key(PresetConfig, ('denature', 'additives')),
                default_factory=list,
        )
        denature_incubation = byoc.param(
                Key(PresetConfig, ('denature', 'incubation')),
                default=None,
        )
        config_paths = byoc.config_attr()
        preset_briefs = byoc.config_attr()

    def get_protocol(self):
        return dnase_digest(self.samples)

    __config__ = [
            DocoptConfig,
    ]
    samples = byoc.param(
            Key(DocoptConfig, samples_from_docopt),
            get=bind,
    )

if __name__ == '__main__':
    Dnase.entry_point()
