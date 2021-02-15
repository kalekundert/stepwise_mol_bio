#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import fatal, warn, plural
from appcli import Key, DocoptConfig
from stepwise import StepwiseConfig, PresetConfig
from stepwise_mol_bio import Main, require_reagent, merge_dicts
from copy import deepcopy
from operator import not_

class Po4Config(appcli.Config):
    # This should be a general class with the following features:
    # - Don't load the database until needed
    # - Some option to specify the aggregation function for each attribute.  
    #   Can probably default to `unanimous()`.
    # - Move to PO₄ package.
    autoload = False

    def load(self, obj):
        try:
            import po4
        except ImportError:
            return

        db = po4.load_db()

        try:
            constructs = [db[x] for x in obj.templates]
        except po4.QueryError:
            return

        values = {}
        
        try:
            values['conc_nM'] = po4.unanimous(
                    x.conc_nM for x in constructs
            )
        except (ValueError, po4.ParseError):
            pass

        try:
            values['is_rna'] = po4.unanimous(
                    x.is_rna for x in constructs
            )
        except (ValueError, po4.ParseError):
            pass

        yield appcli.Layer(values=values, location=db.name)

@autoprop
class InVitroTranslation(Main):
    """\
Express proteins from linear DNA templates using NEBExpress.

Usage:
    ivtt <templates>... [-p <name>] [-v <µL>] [-n <rxns>] [-c <nM>] [-C <nM>]
        [-rIX] [-a <name;conc;vol;mm>]... [-t <time>] [-T <°C>]

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

    -v --volume <µL>                    [default: ${app.volume_uL}]
        The volume of the reaction in µL.

    -n --num-reactions <int>
        The number of reactions to set up.  By default, this is inferred from
        the number of templates.

    -c --template-conc <nM>
        The desired final concentration of template in the reaction.  

    -C --template-stock <nM>
        The stock concentration of the template DNA or mRNA, in units of nM.  
        If not specified, a concentration will be queried from the PO₄ 
        database.  In this case, all templates must be in the database and must 
        have identical concentrations.

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

    -t --incubation-time <time>         [default: ${app.incubation_time}]
        The amount of time to incubate the reactions.  No unit is assumed, so 
        be sure to include one.

    -T --incubation-temperature <°C>    [default: ${app.incubation_temp_C}]
        The temperature to incubate the reactions at, in °C.
"""
    __config__ = [
            Po4Config(),
            DocoptConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.ivtt'),
    ]
    preset_briefs = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '--preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    base_reaction = appcli.param(
            Key(PresetConfig, 'reaction'),
            cast=stepwise.MasterMix.from_text,
    )
    title = appcli.param(
            Key(PresetConfig, 'title'),
            Key(PresetConfig, 'kit'),
    )
    templates = appcli.param(
            Key(DocoptConfig, '<templates>'),
    )
    volume_uL = appcli.param(
            Key(DocoptConfig, '--volume', cast=eval),
            Key(PresetConfig, 'volume_uL'),
            Key(StepwiseConfig, 'default_volume_uL'),
            default=None,
    )
    num_reactions = appcli.param(
            Key(DocoptConfig, '--num-reactions', cast=eval),
            default=None,
            get=lambda self, x: x or len(self.templates),
    )
    template_conc_nM = appcli.param(
            Key(DocoptConfig, '--template-conc', cast=float),
            default=None,
    )
    template_stock_nM = appcli.param(
            Key(DocoptConfig, '--template-stock', cast=float),
            Key(PresetConfig, 'template_stock_nM'),
            Key(Po4Config, 'conc_nM'),
    )
    use_mrna = appcli.param(
            Key(DocoptConfig, '--mrna'),
            Key(Po4Config, 'is_rna'),
            default=False,
    )
    use_template = appcli.param(
            Key(DocoptConfig, '--no-template', cast=not_),
            default=True,
    )
    use_rnase_inhibitor = appcli.param(
            Key(DocoptConfig, '--no-inhibitor', cast=not_),
            default=True,
    )
    setup_instructions = appcli.param(
            Key(PresetConfig, 'setup_instructions'),
            default_factory=list,
    )
    setup_footnote = appcli.param(
            Key(PresetConfig, 'setup_footnote'),
            default=None,
    )
    incubation_time = appcli.param(
            Key(DocoptConfig, '--incubation-time'),
            Key(PresetConfig, 'incubation_time'),
    )
    incubation_temp_C = appcli.param(
            Key(DocoptConfig, '--incubation-temp'),
            Key(PresetConfig, 'incubation_temp_C'),
            cast=float,
    )
    incubation_footnote = appcli.param(
            Key(PresetConfig, 'incubation_footnote'),
            default=None,
    )

    def get_protocol(self):
        p = stepwise.Protocol()
        rxn = self.reaction

        p += stepwise.Step(
                f"Setup {plural(self.num_reactions):# {self.title} reaction/s}"
                f"{p.add_footnotes(self.setup_footnote)}:",
                rxn,
                substeps=self.setup_instructions,
        )
        p += stepwise.Step(
                f"Incubate at {self.incubation_temp_C:g}°C for {self.incubation_time}"
                f"{p.add_footnotes(self.incubation_footnote)}."
        )

        return p

    def get_reaction(self):
        rxn = deepcopy(self.base_reaction)
        rxn.num_reactions = self.num_reactions or len(self.templates)

        if self.volume_uL:
            rxn.hold_ratios.volume = self.volume_uL, 'µL'

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

        rxn[template].name = f"{','.join(self.templates)}"
        rxn[template].hold_conc.stock_conc = self.template_stock_nM, 'nM'

        if self.template_conc_nM:
            rxn[template].hold_stock_conc.conc = self.template_conc_nM, 'nM'
        elif self.use_template:
            warn("Template concentrations must be empirically optimized.\nThe default value is just a plausible starting point.")

        if not self.use_template:
            del rxn[template]

        if not self.use_rnase_inhibitor:
            del_reagents_by_flag(rxn, 'rnase')

        if self.use_template:
            rxn.fix_volumes(template)

        return rxn

def del_reagent_if_present(rxn, key):
    if key in rxn:
        del rxn[key]

def del_reagents_by_flag(rxn, flag):
    reagents = list(rxn.iter_reagents_by_flag(flag))
    for reagent in reagents:
        del rxn[reagent.key]


if __name__ == '__main__':
    InVitroTranslation.main()
