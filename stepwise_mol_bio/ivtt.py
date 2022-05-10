#!/usr/bin/env python3

import stepwise, byoc, autoprop
from inform import fatal, warn, plural
from byoc import Key, DocoptConfig
from stepwise import StepwiseConfig, PresetConfig, pl, ul
from stepwise_mol_bio import (
        Main, BindableReagent, bind, require_reagent,
)
from freezerbox import ReagentConfig, unanimous
from more_itertools import flatten
from copy import deepcopy
from operator import not_

# `Ivtt` has the ability to add user-specified reagents and instructions to the 
# reaction.  I'd like for all my reaction-centric protocols to support those 
# features, but I'll need to generalize the implementation a bit first.  The 
# basic idea will be to create a `Reaction` class that inherits from `Main` and 
# provides the basic framework for describing a reaction.  (I might want to 
# choose a different name, to avoid conflicts with `stepwise.Reaction`, but 
# nothing else comes to mind at the moment.)
#
# One problem with `Ivtt` is that it has two versions of many of its 
# attributes.  This is because the volume of an additive specified in a config 
# file should be scaled if the reaction volume is changed on the command line, 
# but the volume of an additive specified on the command line should not.  So 
# the two version are basically "scale with volume" and "don't scale with 
# volume".
#
# To get the same effect without requiring duplicate attributes, I had the 
# thought that I could basically wrap the value of a single attribute in an 
# object that describes when that attribute should be evaluated.  For example, 
# `ivtt.volume_uL = Early(5)` would mean that the volume should be set before 
# reaction is scaled, while `ivtt.volume_uL = 5` would mean after.  I would 
# still need `volume_uL` and `default_volume_uL` attributes, because the volume 
# might be set in both the config file and the command line.  But for every 
# other parameter, only the final specification matters.

def parse_templates_from_docopt(args):
    return [
            Ivtt.Template(tag)
            for tag in args['<templates>']
    ]

def del_reagent_if_present(rxn, key):
    if key in rxn:
        del rxn[key]

def del_reagents_by_flag(rxn, flag):
    reagents = list(rxn.iter_reagents_by_flag(flag))
    for reagent in reagents:
        del rxn[reagent.key]

@autoprop.cache
class Ivtt(Main):
    """\
Express proteins from purified DNA templates.

Usage:
    ivtt <templates>... [-p <name>] [-v <µL>] [-V <µL>] [-n <rxns>]
        [-x <percent>] [-c <nM>] [-C <nM>] [-mrIX] [-a <name;conc;vol;mm>]...
        [-d <name>]...  [-i <instruction>]... [-t <time>] [-T <°C>]

Arguments:
    <templates>
        The templates to express.  The number of reactions will be inferred 
        from this list.

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>
        What default reaction parameters to use.  The following parameters are 
        currently available:

        ${hanging_indent(app.preset_briefs, 8*' ')}

    -v --volume <µL>
        The volume of the reaction in µL.  By default, the volume specified by 
        the reaction table in the chosen preset will be used.

    -n --num-reactions <int>
        The number of reactions to set up.  By default, this is inferred from
        the number of templates.

    -x --extra-percent <percent>
        How much extra master mix to prepare, as a percentage of the minimum 
        required master mix volume.

    -V --template-volume <µL>
        The volume of template to add to the reaction.  This overrides the 
        `--template-conc` option.

    -c --template-conc <nM>
        The desired final concentration of template in the reaction.

    -C --template-stock <nM>
        The stock concentration of the template DNA or mRNA, in units of nM.  
        If not specified, a concentration will be queried from the PO₄ 
        database.  In this case, all templates must be in the database and must 
        have identical concentrations.

    -m --master-mix
        Include the template in the master mix.

    -r --mrna
        Use mRNA as the template instead of DNA.

    -X --no-template
        Don't include the template in the reaction, e.g. as a negative control.

    -I --no-inhibitor
        Don't include RNase inhibitor in the reaction.

    -a --additive <name;conc;vol;mm>
        Add an additional reagent to the reaction.  See `sw reaction -h` for a 
        complete description of the syntax.  This option can be specified 
        multiple times.

    -d --exclude <name>
        Remove the given reagent from the reaction.  The reagent can be 
        identified either by its name or by its flag.  This can be used to 
        remove an additive added by the preset, for example.

    -i --instruction <text>
        Add the given instruction just below the reaction table, after any 
        instructions specified by the preset.  This option can be specified 
        multiple times to add multiple instructions.

    -t --incubation-time <time>
        The amount of time to incubate the reactions.  No unit is assumed, so 
        be sure to include one.  If '0', the incubation step will be removed 
        from the protocol (e.g. so it can be added back at a later point).

    -T --incubation-temperature <°C>
        The temperature to incubate the reactions at, in °C.
"""
    __config__ = [
            DocoptConfig,
            PresetConfig,
            StepwiseConfig.setup(('molbio', 'ivtt')),
    ]
    preset_briefs = byoc.config_attr()
    preset_brief_template = '{kit}'

    class Template(BindableReagent, use_app_configs=True):
        stock_nM = byoc.param(
                Key(DocoptConfig, '--template-stock', cast=float),
                Key(ReagentConfig, 'conc_nM'),
                Key(PresetConfig, 'template_stock_nM'),
                default=None,
        )
        is_mrna = byoc.param(
                Key(DocoptConfig, '--mrna'),
                Key(ReagentConfig, 'molecule', cast=lambda x: x == 'RNA'),
                default=None,  # interpreted as "unknown"
        )

    presets = byoc.param(
            Key(StepwiseConfig, 'presets'),
            pick=list,
    )
    preset = byoc.param(
            Key(DocoptConfig, '--preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    base_reaction = byoc.param(
            Key(PresetConfig, 'reaction'),
            cast=stepwise.MasterMix.from_text,
    )
    title = byoc.param(
            Key(PresetConfig, 'title'),
            Key(PresetConfig, 'kit'),
    )
    templates = byoc.param(
            Key(DocoptConfig, parse_templates_from_docopt),
            get=bind,
    )
    volume_uL = byoc.param(
            Key(DocoptConfig, '--volume', cast=eval),
            default=None,
    )
    default_volume_uL = byoc.param(
            # The difference between `default_volume_uL` and `volume_uL` is 
            # that the default additives are applied to the reaction after the 
            # default volume is set, but before the non-default volume is set.  
            # This allows the volume of the additive to be scaled 
            # proportionally to the volume of the reaction that the additive 
            # was specified for.
            Key(PresetConfig, 'volume_uL'),
            Key(StepwiseConfig, 'default_volume_uL'),
            default=None,
    )
    num_reactions = byoc.param(
            Key(DocoptConfig, '--num-reactions', cast=eval),
            default=None,
            get=lambda self, x: x or len(self.templates),
    )
    extra_percent = byoc.param(
            Key(DocoptConfig, '--extra-percent', cast=float),
            default=10,
    )
    template_volume_uL = byoc.param(
            Key(DocoptConfig, '--template-volume', cast=float),
            default=None,
    )
    default_template_volume_uL = byoc.param(
            # See `default_volume_uL`.
            Key(PresetConfig, 'template_volume_uL'),
            default=None,
    )
    template_conc_nM = byoc.param(
            Key(DocoptConfig, '--template-conc', cast=float),
            Key(PresetConfig, 'template_conc_nM'),
            default=None,
    )
    master_mix = byoc.param(
            Key(DocoptConfig, '--master-mix'),
            default=False,
    )
    use_template = byoc.param(
            Key(DocoptConfig, '--no-template', cast=not_),
            default=True,
    )
    use_rnase_inhibitor = byoc.param(
            Key(DocoptConfig, '--no-inhibitor', cast=not_),
            default=True,
    )
    additives = byoc.param(
            Key(DocoptConfig, '--additive'),
            default_factory=list,
    )
    default_additives = byoc.param(
            # See `default_volume_uL`.
            Key(PresetConfig, 'additives'),
            default_factory=list,
    )
    exclude = byoc.param(
            Key(DocoptConfig, '--exclude', cast=set),
            default_factory=set,
    )
    setup_instructions = byoc.param(
            Key(PresetConfig, 'setup_instructions'),
            Key(DocoptConfig, '--instruction'),
            default_factory=list,
            pick=lambda x: list(flatten(x)),
    )
    setup_footnote = byoc.param(
            Key(PresetConfig, 'setup_footnote'),
            default=None,
    )
    incubation_time = byoc.param(
            Key(DocoptConfig, '--incubation-time'),
            Key(PresetConfig, 'incubation_time'),
    )
    incubation_temp_C = byoc.param(
            Key(DocoptConfig, '--incubation-temp'),
            Key(PresetConfig, 'incubation_temp_C'),
            cast=float,
    )
    incubation_footnote = byoc.param(
            Key(PresetConfig, 'incubation_footnote'),
            default=None,
    )

    @classmethod
    def from_tags(cls, tags):
        return cls([cls.Template(x) for x in tags])

    def __init__(self, templates):
        self.templates = templates

    def get_protocol(self):
        p = stepwise.Protocol()
        rxn = self.reaction

        p += pl(
                f"Setup {plural(self.num_reactions):# {self.title} reaction/s}{p.add_footnotes(self.setup_footnote)}:",
                rxn,
                ul(*self.setup_instructions),
        )
        if self.incubation_time != '0':
            p += f"Incubate at {self.incubation_temp_C:g}°C for {self.incubation_time}{p.add_footnotes(self.incubation_footnote)}."

        return p

    def get_reaction(self):

        def add_reagents(additives):
            nonlocal i
            # It would be better if there was a utility in stepwise for parsing 
            # `sw reaction`-style strings.  Maybe `Reagent.from_text()`.
            for i, additive in enumerate(additives, i):
                reagent, stock_conc, volume, master_mix = (
                        x.strip() for x in additive.split(';'))
                rxn[reagent].stock_conc = stock_conc
                rxn[reagent].volume = volume
                rxn[reagent].master_mix = {'+': True, '-': False, '': False}[master_mix.strip()]
                rxn[reagent].order = i

        def remove_reagents(exclude):
            for reagent in rxn:
                if reagent.name in exclude or reagent.flags & exclude:
                    del rxn[reagent.key]

        rxn = deepcopy(self.base_reaction)
        rxn.num_reactions = self.num_reactions

        for i, reagent in enumerate(rxn):
            reagent.order = i

        if self.default_volume_uL:
            rxn.hold_ratios.volume = self.default_volume_uL, 'µL'

        if self.default_template_volume_uL:
            for template in ['DNA', 'mRNA']:
                if template in rxn:
                    rxn[template].volume = self.default_template_volume_uL, 'µL'

        add_reagents(self.default_additives)

        if self.volume_uL:
            rxn.hold_ratios.volume = self.volume_uL, 'µL'

        add_reagents(self.additives)
        remove_reagents(self.exclude)

        if self.use_mrna:
            template = 'mRNA'
            require_reagent(rxn, 'mRNA')
            del_reagent_if_present(rxn, 'DNA')
            del_reagents_by_flag(rxn, 'dna')
        else:
            template = 'DNA'
            require_reagent(rxn, 'DNA')
            del_reagent_if_present(rxn, 'mRNA')
            del_reagents_by_flag(rxn, 'mrna')

        rxn[template].name = f"{','.join(x.tag for x in self.templates)}"
        rxn[template].master_mix = self.master_mix

        if self.template_stock_nM:
            rxn[template].hold_conc.stock_conc = self.template_stock_nM, 'nM'

        if self.template_volume_uL:
            rxn[template].volume = self.template_volume_uL, 'µL'
        elif self.template_conc_nM:
            rxn[template].hold_stock_conc.conc = self.template_conc_nM, 'nM'
        elif self.use_template:
            warn("Template concentrations must be empirically optimized.\nThe default value is just a plausible starting point.")

        # Make sure the template is added last.
        rxn[template].order = i+1

        if not self.use_rnase_inhibitor:
            del_reagents_by_flag(rxn, 'rnase')

        if not self.use_template:
            del rxn[template]
        else:
            rxn.fix_volumes(template)

        rxn.extra_percent = self.extra_percent

        return rxn

    def get_template_stock_nM(self):
        return min(
                (c for x in self.templates if (c := x.stock_nM)),
                default=None,
        )

    @property
    def use_mrna(self):
        return unanimous(
                items=(
                    r for x in self.templates
                    if (r := x.is_mrna) is not None
                ),
                default=False,
        )

if __name__ == '__main__':
    Ivtt.main()
