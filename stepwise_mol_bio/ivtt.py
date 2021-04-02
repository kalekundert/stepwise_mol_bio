#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import fatal, warn, plural
from appcli import Key, DocoptConfig
from stepwise import StepwiseConfig, PresetConfig, pl, ul
from stepwise_mol_bio import Main, require_reagent, merge_dicts
from copy import deepcopy
from operator import not_

class Po4Config(appcli.Config):
    # This should be a general class---provided by `po4`---with the following 
    # features:
    # - Don't load the database until needed.
    # - Will query any key given in the PO₄ database.
    # - Some option to specify the aggregation function for each key.  Can 
    #   probably default to `unanimous()`.
    # - Always yield layer if PO₄ is installed, even if template names aren't 
    #   valid tags.  This leads to better error messages.
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
        except (ValueError, po4.ParseError, po4.QueryError):
            pass

        try:
            values['is_rna'] = po4.unanimous(
                    x.is_rna for x in constructs
            )
        except (ValueError, po4.ParseError, po4.QueryError):
            pass

        yield appcli.Layer(values=values, location=db.name)

@autoprop
class InVitroTranslation(Main):
    """\
Express proteins from linear DNA templates using NEBExpress.

Usage:
    ivtt <templates>... [-p <name>] [-v <µL>] [-n <rxns>] [-c <nM>] [-C <nM>]
        [-mrIX] [-a <name;conc;vol;mm>]... [-t <time>] [-T <°C>]

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

    -t --incubation-time <time>         [default: ${app.incubation_time}]
        The amount of time to incubate the reactions.  No unit is assumed, so 
        be sure to include one.  If '0', the incubation step will be removed 
        from the protocol (e.g. so it can be added back at a later point).

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
    preset_brief_template = '{kit}'

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
            default=None,
    )
    default_volume_uL = appcli.param(
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
    num_reactions = appcli.param(
            Key(DocoptConfig, '--num-reactions', cast=eval),
            default=None,
            get=lambda self, x: x or len(self.templates),
    )
    template_conc_nM = appcli.param(
            Key(DocoptConfig, '--template-conc', cast=float),
            Key(PresetConfig, 'template_conc_nM'),
            default=None,
    )
    template_stock_nM = appcli.param(
            Key(DocoptConfig, '--template-stock', cast=float),
            Key(Po4Config, 'conc_nM'),
            Key(PresetConfig, 'template_stock_nM'),
    )
    master_mix = appcli.param(
            Key(DocoptConfig, '--master-mix'),
            default=False,
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
    additives = appcli.param(
            Key(DocoptConfig, '--additive'),
            default_factory=list,
    )
    default_additives = appcli.param(
            Key(PresetConfig, 'additives'),
            default_factory=list,
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
                reagent, stock_conc, volume, master_mix = additive.split(';')
                rxn[reagent].stock_conc = stock_conc
                rxn[reagent].volume = volume
                rxn[reagent].master_mix = {'+': True, '-': False, '': False}[master_mix.strip()]
                rxn[reagent].order = i

        rxn = deepcopy(self.base_reaction)
        rxn.num_reactions = self.num_reactions

        for i, reagent in enumerate(rxn):
            reagent.order = i

        if self.default_volume_uL:
            rxn.hold_ratios.volume = self.default_volume_uL, 'µL'

        add_reagents(self.default_additives)

        if self.volume_uL:
            rxn.hold_ratios.volume = self.volume_uL, 'µL'

        add_reagents(self.additives)

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
        rxn[template].master_mix = self.master_mix
        rxn[template].hold_conc.stock_conc = self.template_stock_nM, 'nM'

        if self.template_conc_nM:
            rxn[template].hold_stock_conc.conc = self.template_conc_nM, 'nM'
        elif self.use_template:
            warn("Template concentrations must be empirically optimized.\nThe default value is just a plausible starting point.")

        if not self.use_template:
            del rxn[template]

        if not self.use_rnase_inhibitor:
            del_reagents_by_flag(rxn, 'rnase')

        # Make sure the template is added last.
        rxn[template].order = i+1

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
