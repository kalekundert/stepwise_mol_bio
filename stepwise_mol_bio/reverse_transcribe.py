#!/usr/bin/env python3

import stepwise
import byoc
import autoprop

from stepwise import (
        StepwiseConfig, PresetConfig, Reaction, Reactions,
        format_reaction, pl, ul,
)
from stepwise_mol_bio import (
        App, Bindable, UsageError,
        bind, group_samples,
)
from stepwise_mol_bio.dnase import Dnase, dnase_digest
from stepwise_mol_bio.thermocycler import (
        Thermocycler, format_thermocycler_steps
)
from freezerbox import ReagentConfig
from byoc import Key, Method, DocoptConfig
from operator import attrgetter
from functools import partial
from more_itertools import first

def reverse_transcribe(samples):
    return stepwise.Protocol.merge(
            *plan_reverse_transcription_protocol.for_group_in(samples)
    )
    
@group_samples(
        'reaction_prototype',
        'primer_type',
        'primer_aliases',
        'anneal',
        'anneal_volume',
        'anneal_incubation',
        'no_rt_control',
        'no_rt_denature',
        'incubation',
)
def plan_reverse_transcription_protocol(group):
    p = stepwise.Protocol()

    p += plan_dnase_protocol(group)

    if group.no_rt_control and group.no_rt_denature:
        p += plan_denature_protocol(group)

    elif group.anneal:
        p += plan_standard_protocol(group)

    else:
        p += plan_quick_reactions(group)

    p += plan_incubation(group)

    return p

@group_samples(
        'reaction_prototype',
)
def plan_dnase_protocol(group):
    samples = [x.dnase for x in group if x.dnase]
    if not samples:
        return stepwise.Protocol()

    # This method modifies the reaction prototype of group object it's given.  
    # I'd much rather not do this, but right now it's by far the most succinct 
    # way to get the behavior I want.  I might have to think more about how to 
    # communicate between these `@group_samples` functions, though...

    group.reaction_prototype['template RNA'].name = "DNase-treated template RNA"
    group.reaction_prototype['template RNA'].volume = max(
            x.reaction_prototype.volume
            for x in samples
    )
    group.reaction_prototype['template RNA'].stock_conc = None
    return dnase_digest(samples)

@group_samples
def plan_denature_protocol(group):
    if group.anneal:
        raise UsageError("requested both `denature` and `anneal`.  `denature` should only be used when the enzyme and primer are premixed.  `anneal` should only be used when they are separate.")

    if list(group.reaction_prototype.iter_reagents_by_flag('primer')):
        raise UsageError("`denature` requested, but reaction includes primer.  `denature` should only be used when the enzyme and primer are premixed.")

    base_rxn = make_reaction(group)

    enz_rxn = base_rxn.copy()
    with enz_rxn.hold_solvent_volume():
        del enz_rxn['template RNA']

    rt_rxn = Reaction()
    rt_rxn['RT master mix'].volume = enz_rxn.volume
    rt_rxn['RT master mix'].flags.add('RT')
    rt_rxn['template RNA'] = base_rxn['template RNA']

    combos = make_combos(
            group, rt_rxn,
            rt_controls=('denatured', 'non-denatured'),
    )

    rt_rxns = Reactions(rt_rxn, combos)
    rt_rxns.step_kind = 'reverse transcription'

    # Normally it would be advised to make a `Reactions` object for the enzyme 
    # mix, in case the user specifies some reaction that somehow needs master 
    # mix steps.  This is not possible here, though, because this function 
    # controls how the combos are generated and removes the only reagent that 
    # can be varied from the reaction.  So there will never be any master 
    # mixes, and the extra complexity isn't really called for.
    #
    # We calculate how to scale the reaction using the extra object from 
    # `rt_rxns` object, because the user can configure the default behavior of 
    # that object, and we want to use the users configuration.
    enz_scale = rt_rxns.extra.increase_scale(len(combos) / 2, enz_rxn)

    p = stepwise.Protocol()
    p += pl(
            "Prepare reverse transcriptase (RT) master mix:",
            format_reaction(enz_rxn, scale=enz_scale, show_totals=False),
    )
    p += pl(
            "Denature half of the RT master mix:",
            ul(
                f"Take {enz_rxn.volume * enz_scale / 2:.2f} of the master mix.",
                *format_thermocycler_steps(group.no_rt_denature, incubate_prefix=True),
            ),
    )
    p += rt_rxns

    return p

@group_samples
def plan_quick_reactions(group):
    rxn = make_reaction(group)
    combos = make_combos(group, rxn)

    rxns = Reactions(rxn, combos)
    rxns.step_kind = 'reverse transcription'

    rename_rt_mix(rxns)

    return rxns

@group_samples
def plan_standard_protocol(group):
    p = stepwise.Protocol()
    anneal_rxns, rt_rxns = plan_standard_reactions(group)

    p += anneal_rxns
    p += Thermocycler(group.anneal_incubation)
    p += rt_rxns

    return p

@group_samples(
        'anneal_volume',
        'no_rt_control',
        'reaction_prototype',
        'primer_type',
)
def plan_standard_reactions(group):
    if not list(group.reaction_prototype.iter_reagents_by_flag('anneal')):
        raise UsageError("`anneal` requested, but not supported for the chosen reaction")

    def make_anneal_combo(d):
        keys = 'template RNA', 'primer'
        return {k: d[k] for k in keys if k in d}

    def make_rt_combo(d):
        d = d.copy()
        template = d.pop('template RNA')
        primer = d.pop('primer', None)
        d['annealed template/primer'] = (template, primer)
        return d

    base_rxn = make_reaction(group)
    anneal_rxn = base_rxn.copy()
    rt_rxn = base_rxn.copy()

    for reagent in base_rxn.iter_nonsolvent_reagents():
        del_rxn = rt_rxn if 'anneal' in reagent.flags else anneal_rxn
        del del_rxn[reagent.key]

    if group.anneal_volume:
        anneal_rxn.volume = group.anneal_volume, anneal_rxn.volume.unit
    else:
        anneal_rxn.volume /= 2

    rt_rxn['annealed template/primer'].volume = anneal_rxn.volume

    combos = make_combos(group, base_rxn)
    anneal_combos = [make_anneal_combo(x) for x in combos]
    rt_combos = [make_rt_combo(x) for x in combos]

    anneal_rxns = Reactions(anneal_rxn, anneal_combos)
    anneal_rxns.step_kind = 'primer-annealing'
    anneal_rxns.split_replicates = False

    rt_rxns = Reactions(rt_rxn, rt_combos)
    rt_rxns.step_kind = 'reverse transcription'

    if group.no_rt_control:
        rt_rxns.combo_step = "Setup each reaction with and without reverse transcriptase."
    else:
        rt_rxns.combo_step = None

    return anneal_rxns, rt_rxns

@group_samples
def plan_incubation(group):
    incubation = group.incubation

    if isinstance(group.incubation, dict):
        primer = choose_primer(group.primer_type, group.reaction_prototype)
        incubations = {
                group.primer_aliases.get(k, k): v
                for k, v in incubation.items()
        }
        incubation = incubations[primer]

    return Thermocycler(incubation)

@group_samples(
        'reaction_prototype',
        'primer_type',
)
def make_reaction(group):
    rxn = group.reaction_prototype.copy()
    setup_template(group, rxn)
    setup_primer(group, rxn)
    return rxn

@group_samples
def setup_template(group, rxn):
    template_stocks_ng_uL = [
            x.template_stock_ng_uL
            for x in group
            if x.template_stock_ng_uL
    ]
    if template_stocks_ng_uL:
        rxn['template RNA'].hold_conc.stock_conc = \
                min(template_stocks_ng_uL), 'ng/µL'

        rxn.repair_volumes('template RNA')

@group_samples('primer_type')
def setup_primer(group, rxn):
    primers = [x.key for x in rxn.iter_reagents_by_flag('primer')]
    if not primers:
        return

    chosen_primer = choose_primer(group.primer_type, rxn)

    for primer in primers:
        if primer != chosen_primer:
            del rxn[primer]

    rxn[chosen_primer].key = 'primer'
    rxn['primer'].flags.discard('primer')
    rxn['primer'].make_combos = False

    if chosen_primer == 'gene-specific primer':
        rxn['primer'].make_combos = True

        primer_stocks_uM = [
                x.primer_stock_uM
                for x in group
                if x.primer_stock_uM
        ]
        if primer_stocks_uM:
            rxn['primer'].hold_conc.stock_conc = \
                    min(primer_stocks_uM), 'µM'

def choose_primer(choice, rxn):
    if choice:
        return choice
    else:
        return first(
                (x.key for x in rxn.iter_reagents_by_flag('primer')),
                default=None,
        )

@group_samples('no_rt_control')
def make_combos(group, rxn, rt_controls=('+', '−')):

    def make_combo(sample):
        combo = {'template RNA': sample.template}
        if 'primer' in rxn and rxn['primer'].make_combos:
            combo['primer'] = sample.primer
        return combo

    combos = [make_combo(x) for x in group]

    if group.no_rt_control and rt_controls:
        rt_reagents = [
                x.key
                for x in rxn
                if 'RT' in x.flags
        ]
        combos, base_combos = [], combos
        for combo in base_combos:
            combos += [
                    {**combo, **{k: control for k in rt_reagents}}
                    for control in rt_controls
            ]

    return combos

def rename_rt_mix(rxns):
    rt_keys = {
            x.key
            for x in rxns.base_reaction
            if 'RT' in x.flags
    }
    for mix in rxns.unprocessed_mixes:
        if rt_keys == mix.reagents:
            mix.name = 'RT'
            break

def samples_from_docopt(args):
    primer = args.get('--primer')
    return [
        ReverseTranscribe.Sample(template, primer)
        for template in args['<templates>']
    ]

def dnase_factory_from_preset(preset):
    def factory(rt):
        rxn = rt.reaction_prototype
        return Dnase.Sample(
                rt.template,
                preset=preset,
                rna_volume=rxn['template RNA'].volume.value,
                rna_stock_conc=rxn['template RNA'].stock_conc,
        )
    return factory

class TemplateConfig(ReagentConfig):
    tag_getter = attrgetter('template')

class PrimerConfig(ReagentConfig):
    tag_getter = attrgetter('primer')

@autoprop
class ReverseTranscribe(App):
    """
Synthesize DNA from a RNA template.

Usage:
    reverse_transcribe <templates>... [-p <preset>] [-P primer] [-gaRD]

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>
        What set of default reaction parameters to use.  The following presets 
        are currently available:

        ${hanging_indent(sample.preset_briefs, 8)}

    -C --template-stock
        The stock concentration of the template RNA in ng/µL.  If not 
        specified, the protocol will attempt to look up stock concentrations 
        for each template in the FreezerBox database.  If it finds any, it will 
        use the lowest value found.  Otherwise, it will use whatever stock 
        concentration is specified in the preset.

    -P --primer <name>
        The name of the primer to use.  The names that are allowed depend on 
        the preset, although 'dt' and 'hex' are common names for oligo-dT and 
        random hexamers, respectively.  Some presets (e.g. those representing 
        master mixes with primers included) don't allow any primers to be 
        specified at all.  Some presets allow "gene-specific primers", in which 
        case the name specified here can be anything as long as the `-g` flag 
        is given.  The stock concentration of gene-specific primers is looked 
        up in the FreezerBox database, and falls back on the value specified in 
        the preset.

    -g --gene-specific-primer
        Indicate that a gene-specific primer is being used.  See the `--primer` 
        option for more information.

    -a --anneal
        Anneal the primers to the templates before starting the RT reaction.  
        The default is to skip this step, which simplifies the setup of the 
        reaction at the possible expense of some yield.

    -R --no-rt-control
        Include a −RT control for each template.  This is an important control 
        for many downstream analysis steps, e.g. qPCR.

    -D --no-dnase
        Skip the DNase treatment step, if the preset specifies one.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(sample.config_paths, 8)}

    molbio.reverse_transcribe.default_preset
        The default value for the `--preset` option.

    molbio.reverse_transcribe.presets:
        Named groups of default reaction parameters.  Typically each preset 
        corresponds to a particular kit or protocol.  See below for the various 
        settings that can be specified in each preset.

    molbio.reverse_transcribe.presets.<name>.brief:
        A brief description of the preset.  This is displayed in the usage info 
        for the `--preset` option.

    molbio.reverse_transcribe.presets.<name>.inherit:
        Copy all settings from another preset.  This can be used to make small 
        tweaks to a protocol, e.g. "SuperScript with a non-standard additive".

    molbio.reverse_transcribe.presets.<name>.reaction:
        A table detailing all the components of the reverse transcription 
        reaction, in the format understood by `stepwise.Reaction.from_text()`.  
        The reaction must have a reagent named "template RNA", and the 
        following flags are important to specify:

        'anneal': Use this flag to label every reagent that should be included 
        in the annealing reaction, e.g. template, primer, and maybe some 
        buffer.  If no reagents have this flag, the `--anneal` option will be 
        disabled.

        'primer': The reaction may include multiple different primer options, 
        each labeled with this flag.  The user can select which of these 
        primers to use via the `--primer` option.  By default, the first will 
        be used.  Use the special name "gene-specific primer" to specify the 
        volume/concentration for custom-ordered oligos.

        'RT': Use this flag to label the reagents that should be left out when 
        setting up the −RT control.

    molbio.reverse_transcribe.presets.<name>.primer:
        Specify short names for each primer in the reaction.  For example, it 
        is conventional the define 'dt' and 'hex' as aliases for whatever 
        longer names are used in the reaction table for oligo-dT primers and 
        random hexamers, respectively.  This allows users to be more succinct 
        when using the `--primer` option.

    molbio.reverse_transcribe.presets.<name>.anneal:
        The default value for the `--anneal` flag.

    molbio.reverse_transcribe.presets.<name>.anneal_volume:
        The volume to make each annealing reaction, if an annealing step is 
        requested.  Don't specify units; the units are assumed to be the same 
        as for the reaction itself (typically µL).  The default is half of the 
        total reaction volume.

    molbio.reverse_transcribe.presets.<name>.anneal_incubation:
        Describe the thermocycler protocol to use for the annealing step, in 
        the format expected by `format_thermocycler_steps`.  

    molbio.reverse_transcribe.presets.<name>.no_rt_denature:
        Specify the thermocycler protocol that should be used to heat-denature 
        the −RT control.  Setting this option also indicates that the RT should 
        be denatured rather than simply left out of the −RT control.  This is 
        sometimes recommended for all-in-one master mixes where the −RT control 
        would otherwise just be the template.

    molbio.reverse_transcribe.presets.<name>.incubation:
        Describe the thermocycler protocol to use for the reverse transcription 
        reaction, in the format expected by `format_thermocycler_steps`.  This 
        option can also be a dictionary mapping primer names to thermocycler 
        protocols, in the even that you want to use different incubation steps 
        for different primers.
"""

    @autoprop
    class Sample(Bindable, use_app_configs=True):

        def __init__(self, template, primer=None, **kwargs):
            super().__init__(**kwargs)
            self.template = template
            if primer: self.primer = primer

        def _lookup_incubation(self):
            reverse_aliases = {v: k for k, v in self.primer_aliases.items()}
            key = choose_primer(self.primer_type, self.reaction_prototype)
            key = reverse_aliases.get(key, key)
            return self.incubation_by_primer[key]

        __config__ = [
                PresetConfig,
                TemplateConfig,
                PrimerConfig,
                StepwiseConfig.setup(('molbio', 'reverse_transcribe')),
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
        template_stock_ng_uL = byoc.param(
                Key(DocoptConfig, '--template-stock'),
                Key(TemplateConfig, 'conc_ng_uL'),
                default=None,
        )
        primer = byoc.param(
                Key(DocoptConfig, '--primer'),
                default=None,
        )
        primer_stock_uM = byoc.param(
                Key(PrimerConfig, 'conc_uM'),
                default=None,
        )
        primer_aliases = byoc.param(
                Key(PresetConfig, 'primers'),
                default_factory=dict,
        )
        gene_specific_primer = byoc.param(
                Key(DocoptConfig, '--gene-specific-primer'),
                default=False,
        )
        dnase = byoc.param(
                Key(DocoptConfig, '--no-dnase', cast=lambda x: None),
                Method(lambda self: self._dnase_factory(self)),
                default=None,
        )
        _dnase_factory = byoc.param(
                Key(PresetConfig, 'dnase_preset', cast=dnase_factory_from_preset),
        )
        anneal = byoc.param(
                Key(DocoptConfig, '--anneal'),
                Key(PresetConfig, 'anneal'),
                default=False,
        )
        anneal_volume = byoc.param(
                Key(PresetConfig, 'anneal_volume'),
                default=None,
        )
        anneal_incubation = byoc.param(
                Key(PresetConfig, 'anneal_incubation'),
                default=None,
        )
        no_rt_control = byoc.param(
                Key(DocoptConfig, '--no-rt-control'),
                default=False,
        )
        no_rt_denature = byoc.param(
                Key(PresetConfig, 'no_rt_denature'),
                default=None,
        )
        incubation = byoc.param(
                Method(
                    _lookup_incubation,
                    skip=(KeyError, byoc.NoValueFound),
                ),
                Key(PresetConfig, 'incubation'),
        )
        incubation_by_primer = byoc.param(
                Key(PresetConfig, 'incubation_by_primer'),
        )
        config_paths = byoc.config_attr()
        preset_briefs = byoc.config_attr()

        def get_primer_type(self):
            if self.gene_specific_primer:
                return 'gene-specific primer'
            else:
                return self.primer_aliases.get(self.primer, self.primer)

    def __init__(self, samples):
        self.samples = samples

    def get_protocol(self):
        return reverse_transcribe(self.samples)

    __config__ = [
            DocoptConfig,
    ]
    samples = byoc.param(
            Key(DocoptConfig, samples_from_docopt),
            get=bind,
    )

if __name__ == '__main__':
    ReverseTranscribe.entry_point()
