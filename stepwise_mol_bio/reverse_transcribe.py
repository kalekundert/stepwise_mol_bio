#!/usr/bin/env python3

import stepwise
import byoc
import autoprop

from stepwise import (
        StepwiseConfig, PresetConfig, Reaction, Reactions,
        iter_all_mixes, format_reaction, pl, ul, before,
)
from stepwise_mol_bio import (
        App, Bindable, UsageError,
        bind, group_samples,
)
from stepwise_mol_bio.dnase import Dnase, dnase_digest
from stepwise_mol_bio.thermocycler import (
        Thermocycler, format_thermocycler_steps
)
from freezerbox import ReagentConfig, unanimous
from byoc import Key, Method, DocoptConfig
from operator import attrgetter
from functools import partial
from collections import Counter
from more_itertools import one, first, unique_everseen as unique
from inform import plural

def reverse_transcribe(samples):
    return plan_dnase_protocol.concat(samples) + \
            plan_reverse_transcription_protocol.concat(samples)
    
@group_samples(
        'reaction_prototype',
        'anneal',
        'anneal_volume',
        'anneal_incubation',
        'denature_rt',
        'incubation',
        'incubation_by_primer',
)
def plan_reverse_transcription_protocol(group):
    p = stepwise.Protocol()

    if group.denature_rt and not all(x.include_rt for x in group):
        p += plan_denature_protocol(group)

    elif group.anneal:
        p += plan_standard_protocol(group)

    else:
        p += plan_quick_reactions(group)

    p += plan_incubation(group)

    return p

@group_samples
def plan_dnase_protocol(group):
    dnase_samples = [x.dnase for x in group if x.dnase]
    if not dnase_samples:
        return stepwise.Protocol()

    # This method modifies the reaction prototype of group object it's given.  
    # I'd much rather not do this, but right now it's by far the most succinct 
    # way to get the behavior I want.  I might have to think more about how to 
    # communicate between these `@group_samples` functions, though...

    dnase_volume = max(
            x.reaction_prototype.volume
            for x in dnase_samples
    )

    for rt_sample in group:
        if rt_sample.dnase:
            rt_sample.reaction_prototype['template RNA'].name = "DNase-treated template RNA"
            rt_sample.reaction_prototype['template RNA'].volume = dnase_volume
            rt_sample.reaction_prototype['template RNA'].stock_conc = None

    return dnase_digest(dnase_samples)

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
    n = len(combos)

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
    enz_scale = rt_rxns.extra.increase_scale(n, enz_rxn)
    denature_frac = sum(not x.include_rt for x in group) / n

    p = stepwise.Protocol()
    p += pl(
            f"Prepare enough reverse transcriptase (RT) master mix for {plural(n):# reaction/s}:",
            format_reaction(enz_rxn, scale=enz_scale, show_totals=False),
    )
    p += pl(
            f"Denature {enz_rxn.volume * enz_scale * denature_frac:.2f} of the RT master mix:",
            ul(
                *format_thermocycler_steps(group.denature_rt, incubate_prefix=True),
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
        'reaction_prototype',
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

    if not all(x.include_rt for x in group):
        rt_rxns.combos_step = "Setup each reaction with and without reverse transcriptase."
    else:
        rt_rxns.combos_step = None

    return anneal_rxns, rt_rxns

@group_samples(
        'reaction_prototype',
        'incubation',
        'incubation_by_primer',
)
def plan_incubation(group):
    default = group.incubation

    if group.incubation_by_primer:
        rxn = group.reaction_prototype
        incubations = {
                get_primer_key(k, rxn): v
                for k, v in group.incubation_by_primer.items()
        }
        primer_keys = Counter(x.primer_key for x in group)

        if len(primer_keys) == 1:
            key = one(primer_keys)
            return Thermocycler(incubations.get(key, default))

        else:
            step = pl(
                    "Run the following thermocycler protocols:",
            )
            for k in primer_keys:
                step += f"For the {rxn[k].name} {plural(primer_keys[k]):reaction/s}:"
                step += format_thermocycler_steps(incubations.get(k, default))

            return stepwise.Protocol(steps=[step])

    else:
        return Thermocycler(default)

@group_samples(
        'reaction_prototype',
)
def make_reaction(group):
    rxn = group.reaction_prototype.copy()
    setup_template(group, rxn)
    setup_primer(group, rxn)
    setup_enzyme(rxn)
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

@group_samples
def setup_primer(group, rxn):
    primers = [x.key for x in rxn.iter_reagents_by_flag('primer')]
    if not primers:
        return

    primer_keys = {x.primer_key for x in group}

    rxn.insert_reagent('primer', before(*primers))

    # Name:

    try:
        rxn['primer'].name = unanimous(rxn[k].name for k in primer_keys)
    except:
        pass

    # Stock concentration:

    if 'gene-specific primer' in rxn:
        primer_stocks_uM = [
                x.primer_stock_uM
                for x in group
                if x.primer_stock_uM
        ]
        if primer_stocks_uM:
            rxn['gene-specific primer'].hold_conc.stock_conc = \
                    min(primer_stocks_uM), 'µM'

    try:
        rxn['primer'].stock_conc = unanimous(
                rxn[k].stock_conc
                for k in primer_keys
        )
    except ValueError:
        pass

    # Volume:

    rxn['primer'].volume = max(rxn[k].volume for k in primer_keys)

    # Flags:

    flags = (rxn[k].flags for k in primer_keys)
    rxn['primer'].flags = set.intersection(*flags)
    rxn['primer'].flags.discard('primer')

    # Remove user-provided primer prototypes:

    for k in primers:
        del rxn[k]

def setup_enzyme(rxn):
    for reagent in rxn.iter_reagents_by_flag('RT'):
        reagent.flags.add('careful')

@group_samples
def make_combos(group, rxn, rt_controls=('+', '−')):

    def make_combo(sample):
        combo = {
                'template RNA': sample.template,
                **make_primer_combo(sample),
        }
        for k in rt_keys:
            combo |= make_rt_combo(k, rt_controls[not sample.include_rt])

        return combo

    def make_primer_combo(sample):
        if 'primer' not in rxn:
            return {}

        key = sample.primer_key
        rxn0 = sample.reaction_prototype
        gene_specific_primer = get_primer_key(
                'gene-specific primer', rxn0,
                required=False,
        )

        if gene_specific_primer and key == gene_specific_primer:
            name = sample.primer
        else:
            name = rxn0[key].name

        if rxn['primer'].stock_conc or not rxn0[key].stock_conc:
            return {'primer': name}

        stock = rxn0[key].stock_conc
        dilution = rxn['primer'].volume / rxn0[key].volume

        if dilution == 1:
            return {'primer': f'{stock} {name}'}

        try:
            stock = stock / dilution
        except TypeError:
            return {'primer': f'{stock} {name} (diluted {dilution:g}x)'}
        else:
            return {'primer': f'{stock} {name}'}

    rt_keys = get_rt_keys(rxn)
    include_rt = {x.include_rt for x in group}

    if all(include_rt):
        make_rt_combo = lambda k, v: {}
    elif not any(include_rt):
        make_rt_combo = lambda k, v: {k: f'{k} ({v})'}
    else:
        make_rt_combo = lambda k, v: {k: v}

    return [make_combo(x) for x in group]

def rename_rt_mix(rxns):
    rt_keys = get_rt_keys(rxns.base_reaction)
    for mix in iter_all_mixes(rxns.mix):
        if rt_keys == mix.reagents:
            mix.name = 'RT'
            rxns.refresh_names()
            break

def get_rt_keys(rxn):
    return {x.key for x in rxn.iter_reagents_by_flag('RT')}

def get_primer_key(key_or_name, rxn, required=True):
    """
    Given either the key or the name of a primer, return the corresponding key.

    If no matching primer can be found in the given reaction, an error will be 
    raised.
    """
    primer_keys = {x.key for x in rxn.iter_reagents_by_flag('primer')}

    if key_or_name in primer_keys:
        return key_or_name

    primer_names = {rxn[k].name: k for k in primer_keys}

    try:
        return primer_names[key_or_name]

    except KeyError:
        if not required:
            return None

        def format_primers(rxn):
            primer_strs = []

            for reagent in rxn.iter_reagents_by_flag('primer'):
                key, name = reagent.key, reagent.name

                if key != name:
                    primer_strs.append(key)
                else:
                    primer_strs.append(f'{key}, {name}')

            return '\n'.join(primer_strs)

        err = UsageError(primer=key_or_name, rxn=rxn)
        err.brief = "RT reaction has no primer named {key_or_name!r}"
        err.info += lambda e: "known primers:\n" + format_primers(e)
        err.hint += "different RT presets can have different primers, are you using the right preset?"
        raise err


def samples_from_docopt(args):
    rxns = []
    primer = args.get('--primer')

    for template in args['<templates>']:
        rxn = ReverseTranscribe(template, primer)
        rxns.append(rxn)

        if args.get('--no-rt-control'):
            rxn = ReverseTranscribe(template, primer)
            rxn.include_rt = False
            rxns.append(rxn)

    return rxns

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
class ReverseTranscribe(Bindable, use_app_configs=True):

    @classmethod
    def make(cls, samples):
        return reverse_transcribe(samples)


    def __init__(self, template, primer=None, **kwargs):
        super().__init__(**kwargs)
        self.template = template
        if primer: self.primer = primer

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
    gene_specific_primer = byoc.param(
            Key(DocoptConfig, '--gene-specific-primer'),
            default=False,
    )
    include_rt = byoc.param(
            default=True,
    )
    denature_rt = byoc.param(
            Key(PresetConfig, 'no_rt_denature'),
            default=None,
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
    incubation = byoc.param(
            Key(PresetConfig, 'incubation'),
    )
    incubation_by_primer = byoc.param(
            Key(PresetConfig, 'incubation_by_primer'),
            default_factory=dict,
    )
    config_paths = byoc.config_attr()
    preset_briefs = byoc.config_attr()

    def get_primer_key(self):
        """
        Get the key for looking up the primer in the reaction table.

        The `primer` attribute is directly set by the user, which makes it hard 
        to use directly:

        - It may be either the name or the key of a primer reagent in the 
          prototype reaction.  (Names are displayed and are typically long, 
          keys are used internally and are typically short.)

        - It may refer to a gene-specific primer, in which case it won't match 
          any of the reagents in the prototype reaction.  In this case, the 
          special key "gene-specific primer" should be used.

        - It may not be specified, indicating that the default primer should be 
          chosen.

        This property handles all this and either returns a valid key or raises 
        an error.
        """
        rxn = self.reaction_prototype

        if self.gene_specific_primer:
            key_or_name = 'gene-specific primer'
        else:
            key_or_name = self.primer

        if key_or_name:
            return get_primer_key(key_or_name, rxn)
        else:
            return first(
                    (x.key for x in rxn.iter_reagents_by_flag('primer')),
                    default=None,
            )

@autoprop
class ReverseTranscribeCli(App):
    """
Synthesize DNA from a RNA template.

Usage:
    reverse_transcribe <templates>... [-p <preset>] [-C <stock>] [-P <primer>] 
        [-gaRD]

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
        the preset, although 'dt' and 'hex' are common aliases for oligo-dT and 
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
        reaction, in the format expected by `format_thermocycler_steps`.  See 
        also the `incubation_by_primer` setting, which allows different 
        protocols to be chosen based on the primer being used.

    molbio.reverse_transcribe.presets.<name>.incubation_by_primer:
        A dictionary mapping primer names to thermocycler protocols, in the 
        event that you want to use different incubation steps for different 
        primers.  This setting supersedes the `incubation` setting if the 
        primer being used is present in the dictionary.
"""
    Sample = ReverseTranscribe
    samples = byoc.param(
            Key(DocoptConfig, samples_from_docopt),
            get=bind,
    )

if __name__ == '__main__':
    ReverseTranscribeCli.entry_point()


