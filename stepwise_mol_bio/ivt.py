#!/usr/bin/env python3

import stepwise, appcli, autoprop
from stepwise import (
        StepwiseConfig, PresetConfig, MasterMix,
        paragraph_list, unordered_list,
)
from stepwise_mol_bio import (
        Main, UsageError,
        merge_dicts, format_min,
)
from freezerbox import (
        ReagentConfig, MakerArgsConfig,
        parse_volume_uL, parse_time_m, parse_temp_C, convert_conc_unit,
        unanimous, group_by_identity, normalize_seq, join_lists,
)
from appcli import DocoptConfig, Key, Method
from inform import plural, warn

def parse_reaction(reaction_input):
    if isinstance(reaction_input, dict):
        return {k: MasterMix(v) for k, v in reaction_input.items()}
    else:
        return MasterMix(reaction_input)

def pick_by_short(values, is_short):
    if not isinstance(values, dict):
        return values
    else:
        return values['short' if is_short else 'long']

def affected_by_short(values):
    if not isinstance(values, dict):
        return False
    else:
        return values['short'] != values['long']

class TemplateConfig(ReagentConfig):
    tag_getter = lambda app: app.templates

@autoprop.cache
class Ivt(Main):
    """\
Synthesize RNA using an in vitro transcription reaction.

Usage:
    ivt <templates>... [options]

Arguments:
    <templates>
        The names to the DNA templates to transcribe.  If these names can be 
        looked up in a FreezerBox database, some parameters (e.g. template 
        length) will automatically be determined.

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>          [default: ${app.preset}]
        What default reaction parameters to use.  The following parameters are 
        currently available:

        ${hanging_indent(app.preset_briefs, 8)}

    -v --volume <µL>
        The volume of the transcription reaction, in µL.

    -C --template-stock <ng/µL>
        The stock concentration of the template DNA (in ng/µL).  By default, 
        the concentration specified in the preset will be used.  Changing this 
        option will not change the quantity of template added to the reaction 
        (the volume will be scaled proportionally).

    -t --template-mass <ng>
        How much template DNA to use (in ng).  The template volume will be 
        scaled in order to reach the given quantity.  If the template is not 
        concentrated enough to reach the given quantity, the reaction will just 
        contain as much DNA as possible instead.  In order to use this option, 
        the stock concentration of the template must be in units of ng/µL.

    -V --template-volume <µL>
        The volume of DNA to add to the reaction (in µL).  This will be 
        adjusted if necessary to fit within the total volume of the reaction.

    -s --short
        Indicate that all of the templates are shorter than 300 bp.  This 
        allows the reactions to be setup with less material, so that more 
        reactions can be performed with a single kit.

    -i --incubation-time <minutes>
        How long to incubate the transcription reaction, in minutes.

    -T --incubation-temp <°C>
        What temperature the incubate the transcription reaction at, in °C.

    -x --extra <percent>        [default: ${app.extra_percent}]
        How much extra master mix to create.

    -R --toggle-rntp-mix
        Indicate whether you're using an rNTP mix, or whether you need to add 
        each rNTP individually to the reaction.  This option toggles the 
        default value specified in the configuration file.

    -a --toggle-dnase-treatment
        Indicate whether you'd like to include the optional DNase treatment 
        step.  This option toggles the default value specified in the 
        configuration file.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.ivt.default_preset:
        The default value for the `--preset` option.

    molbio.ivt.rntp_mix:
        The default value for the `--toggle-rntp-mix` flag.

    molbio.ivt.dnase_treatment:
        The default value for the `--toggle-dnase-treament` flag.

    molbio.ivt.presets:
        Named groups of default reaction parameters.  Typically each preset 
        corresponds to a particular kit or protocol.  See below for the various 
        settings that can be specified in each preset.

    molbio.ivt.presets.<name>.brief:
        A brief description of the preset.  This is displayed in the usage info 
        for the `--preset` option.

    molbio.ivt.presets.<name>.inherit:
        Copy all settings from another preset.  This can be used to make small 
        tweaks to a protocol, e.g. "HiScribe with a non-standard additive".

    molbio.ivt.presets.<name>.reaction:
        A table detailing all the components of the transcription reaction, in 
        the format understood by `stepwise.MasterMix.from_string()`.  
        Optionally, this setting can be a dictionary with keys "long" and 
        "short", each corresponding to a reaction table.  This allows different 
        reaction parameters to be used for long and short templates.

        The DNA template reagent must be named "template".  The rNTP reagents 
        may be marked with the "rntp" flag.  The `--rntp-mix` flag will replace 
        any reagents so marked with a single reagent named "rNTP mix".

    molbio.ivt.presets.<name>.instructions:
        A list of miscellaneous instructions pertaining to how the reaction 
        should be set up, e.g. how to thaw the reagents, what temperature to 
        handle the reagents at, etc.

    molbio.ivt.presets.<name>.extra_percent:
        How much extra master mix to make, as a percentage of the volume of a 
        single reaction.

    molbio.ivt.presets.<name>.incubation_time_min:
        See `--incubation-time`.  This setting can also be a dictionary with 
        keys "long" and "short", to specify different incubation times for long 
        and short templates.

    molbio.ivt.presets.<name>.incubation_temp_C:
        See `--incubation-temp`.

    molbio.ivt.presets.<name>.length_threshold
        The length of a template in base pairs that separates long from short 
        templates. Template lengths will be queried from the FreezerBox 
        database if possible.

    molbio.ivt.presets.<name>.dnase.default
        The default value for the `--dnase` flag.

    molbio.ivt.presets.<name>.dnase.reaction
        A table detailing all the components of the optional DNase reaction.  
        One component must be named "transcription reaction".

    molbio.ivt.presets.<name>.dnase.incubation_time_min
        The incubation time (in minutes) for the optional DNase reaction.

    molbio.ivt.presets.<name>.dnase.incubation_temp_C
        The incubation temperature (in °C) for the optional DNase reaction.

    molbio.ivt.presets.<name>.dnase.footnotes.reaction
    molbio.ivt.presets.<name>.dnase.footnotes.incubation
    molbio.ivt.presets.<name>.dnase.footnotes.dnase
        Lists of footnotes for the reaction setup, incubation, and DNase 
        treatment steps, respectively.

Database:
    In vitro transcription reactions can appear in the "Synthesis" column of a 
    FreezerBox database.  The associated database entry will automatically be 
    considered ssRNA, e.g. for the purpose of molecular weight calculations:

        ivt template=<tag> [preset=<preset>] [volume=<µL>] [time=<min>]
            [temp=<°C>]

    template=<tag>
        See `<templates>`.  Only one template can be specified.

    preset=<preset>
        See `--preset`.

    volume=<µL>
        See `--volume`.

    time=<min>
        See `--incubation-time`.

    temp=<°C>
        See `--incubation-temp`.

Template Preparation:
    The following information is taken directly from the HiScribe and 
    MEGAscript manuals:

    Plasmid Templates:

        To produce RNA transcript of a defined length, plasmid DNA must be 
        completely linearized with a restriction enzyme downstream of the 
        insert to be transcribed.  Circular plasmid templates will generate 
        long heterogeneous RNA transcripts in higher quantities because of high 
        processivity of T7 RNA polymerase.  Be aware that there has been one 
        report of low level transcription from the inappropriate template 
        strand in plasmids cut with restriction enzymes leaving 3' overhanging 
        ends [Schendorn and Mierindorf, 1985].

        DNA from some miniprep procedures may be contaminated with residual 
        RNase A.  Also, restriction enzymes occasionally introduce RNase or 
        other inhibitors of transcription.  When transcription from a template 
        is suboptimal, it is often helpful to treat the template DNA with 
        proteinase K (100–200 μg/mL) and 0.5% SDS for 30 min at 50°C, follow 
        this with phenol/chloroform extraction (using an equal volume) and 
        ethanol precipitation.

    PCR Templates:

        PCR products containing T7 RNA Polymerase promoter in the correct 
        orientation can be transcribed.  Though PCR mixture can be used 
        directly, better yields will be obtained with purified PCR products.  
        PCR products can be purified according to the protocol for plasmid 
        restriction digests above, or by using commercially available spin 
        columns (we recommend Monarch PCR & DNA Cleanup Kit, NEB #T1030).  PCR 
        products should be examined on an agarose gel to estimate concentration 
        and to confirm amplicon size prior to its use as a template.  Depending 
        on the PCR products, 0.1–0.5 μg of PCR fragments can be used in a 20 μL 
        in vitro transcription reaction.

    Synthetic DNA Oligonucleotides:

        Synthetic DNA oligonucleotides which are either entirely 
        double-stranded or mostly single-stranded with a double-stranded T7 
        promoter sequence can be used for transcription.  In general, the 
        yields are relatively low and also variable depending upon the 
        sequence, purity, and preparation of the synthetic oligonucleotides.
"""
    __config__ = [
            DocoptConfig(),
            MakerArgsConfig(),
            TemplateConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.ivt'),
    ]
    preset_briefs = appcli.config_attr()
    config_paths = appcli.config_attr()

    def _calc_short(self):
        return all(
                x <= self.template_length_threshold
                for x in self.template_lengths
        )

    def _pick_by_short(self, values):
        return pick_by_short(values, self.short)

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '--preset'),
            Key(MakerArgsConfig, 'preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    reaction_prototype = appcli.param(
            Key(PresetConfig, 'reaction', cast=parse_reaction),
            get=_pick_by_short,
    )
    templates = appcli.param(
            Key(DocoptConfig, '<templates>'),
            Key(MakerArgsConfig, 'template', cast=lambda x: [x]),
    )
    template_seqs = appcli.param(
            Key(TemplateConfig, 'seq'),
    )
    template_lengths = appcli.param(
            Key(TemplateConfig, 'length'),
    )
    template_length_threshold = appcli.param(
            Key(PresetConfig, 'length_threshold'),
    )
    template_volume_uL = appcli.param(
            Key(DocoptConfig, '--template-volume', cast=float),
            default=None,
    )
    template_mass_ng = appcli.param(
            Key(DocoptConfig, '--template-mass', cast=float),
            default=None,
    )
    template_stock_ng_uL = appcli.param(
            Key(DocoptConfig, '--template-stock', cast=float),
            Key(TemplateConfig, 'conc_ng_uL', cast=min),
            default=None,
    )
    short = appcli.param(
            Key(DocoptConfig, '--short'),
            Method(_calc_short),
            default=False,
    )
    volume_uL = appcli.param(
            Key(DocoptConfig, '--volume-uL', cast=float),
            Key(MakerArgsConfig, 'volume', cast=parse_volume_uL),
            Key(PresetConfig, 'volume_uL'),
            default=None,
    )
    rntp_mix = appcli.toggle_param(
            Key(DocoptConfig, '--no-rntp-mix', toggle=True),
            Key(StepwiseConfig, 'rntp_mix'),
            default=True,
    )
    extra_percent = appcli.param(
            Key(DocoptConfig, '--extra-percent'),
            Key(PresetConfig, 'extra_percent'),
            cast=float,
            default=10,
    )
    instructions = appcli.param(
            Key(PresetConfig, 'instructions'),
            default_factory=list,
    )
    incubation_times_min = appcli.param(
            Key(DocoptConfig, '--incubation-time'),
            Key(MakerArgsConfig, 'time', cast=parse_time_m),
            Key(PresetConfig, 'incubation_time_min'),
    )
    incubation_temp_C = appcli.param(
            Key(DocoptConfig, '--incubation-temp'),
            Key(MakerArgsConfig, 'temp', cast=parse_temp_C),
            Key(PresetConfig, 'incubation_temp_C'),
    )
    dnase = appcli.toggle_param(
            Key(DocoptConfig, '--toggle-dnase-treatment', toggle=True),
            Key(PresetConfig, 'dnase.treatment'),
            Key(StepwiseConfig, 'dnase_treatment'),
            default=False,
    )
    dnase_reaction_prototype = appcli.param(
            Key(PresetConfig, 'dnase.reaction', cast=MasterMix),
    )
    dnase_incubation_time_min = appcli.param(
            Key(PresetConfig, 'dnase.incubation_time_min'),
    )
    dnase_incubation_temp_C = appcli.param(
            Key(PresetConfig, 'dnase.incubation_temp_C'),
    )
    footnotes = appcli.param(
            Key(PresetConfig, 'footnotes'),
            default_factory=dict,
    )

    group_by = {
        'preset': group_by_identity,
        'volume_uL': group_by_identity,
        'incubation_times_min': group_by_identity,
        'incubation_temp_C': group_by_identity,
    }
    merge_by = {
        'templates': join_lists,
    }

    def __init__(self, templates):
        self.templates = templates

    def __repr__(self):
        return f'Ivt(templates={self.templates!r})'

    def get_reaction(self):
        rxn = self.reaction_prototype.copy()
        rxn.num_reactions = len(self.templates)
        rxn.extra_percent = self.extra_percent

        if self.volume_uL:
            rxn.hold_ratios.volume = self.volume_uL, 'µL'

        if self.rntp_mix:
            rntps = []
            
            for i, reagent in enumerate(rxn):
                reagent.order = i
                if 'rntp' in reagent.flags:
                    rntps.append(reagent)

            if not rntps:
                err = ConfigError("cannot make rNTP mix", preset=self.preset)
                err.blame += "no reagents flagged as 'rntp'"
                err.hints += "you may need to add this information to the [molbio.ivt.{preset}] preset"
                raise err

            rxn['rNTP mix'].volume = sum(x.volume for x in rntps)
            rxn['rNTP mix'].stock_conc = sum(x.stock_conc for x in rntps) / len(rntps)
            rxn['rNTP mix'].master_mix = all(x.master_mix for x in rntps)
            rxn['rNTP mix'].order = rntps[0].order

            for rntp in rntps:
                del rxn[rntp.name]

        rxn['template'].name = ','.join(self.templates)

        if self.template_stock_ng_uL:
            rxn['template'].hold_conc.stock_conc = self.template_stock_ng_uL, 'ng/µL'

        if self.template_mass_ng:
            ng_uL = convert_conc_unit(rxn['template'].stock_conc, None, 'ng/µL')
            rxn['template'].volume = self.template_mass_ng / ng_uL.value, 'µL'

        if self.template_volume_uL:
            rxn['template'].volume = self.template_volume_uL, 'µL'
            if self.template_mass_ng:
                warn(f"template quantity ({self.template_mass_ng} ng) specified but overridden by volume ({self.template_volume_uL} µL)")
        
        rxn.fix_volumes('template', rxn.solvent)
        return rxn

    def get_dnase_reaction(self):
        rxn = self.dnase_reaction_prototype
        rxn.num_reactions = len(self.templates)
        rxn.extra_percent = self.extra_percent

        if self.volume_uL:
            k = (self.volume_uL, 'µL') / rxn['transcription reaction'].volume
            rxn.hold_ratios.volume *= k

        return rxn

    def get_protocol(self):
        p = stepwise.Protocol()

        ## Clean your bench
        p += stepwise.load('rnasezap')

        ## In vitro transcription
        rxn = self.reaction
        n = plural(rxn.num_reactions)
        f = self.footnotes.get('reaction', [])
        p += paragraph_list(
                f"Setup {n:# in vitro transcription reaction/s}{p.add_footnotes(*f)}:",
                rxn,
                unordered_list(*self.instructions),
        )

        f = self.footnotes.get('incubation', [])
        if self.short and affected_by_short(self.incubation_times_min):
            f += [f"Reaction time is different than usual because the template is short (<{self.template_length_threshold} bp)."]
        p += f"Incubate at {self.incubation_temp_C}°C for {format_min(pick_by_short(self.incubation_times_min, self.short))}{p.add_footnotes(*f)}."

        ## DNase treatment
        if self.dnase:
            f = self.footnotes.get('dnase', [])
            p += paragraph_list(
                    f"Setup {n:# DNase reaction/s}{p.add_footnotes(*f)}:",
                    self.dnase_reaction,
            )
            p += f"Incubate at {self.dnase_incubation_temp_C}°C for {format_min(self.dnase_incubation_time_min)}."

        return p
    def get_dependencies(self):
        return self.templates

    def get_product_seqs(self):
        return [transcribe(x) for x in self.template_seqs]

    def get_product_molecule(self):
        return 'ssRNA'

def transcribe(template_seq):
    t7_promoter = 'TAATACGACTCACTATA'
    template_seq = normalize_seq(template_seq)

    i = template_seq.find(t7_promoter)
    j = i + len(t7_promoter)

    if i < 0:
        raise ValueError(f"no T7 promoter found")

    return template_seq[j:].translate(str.maketrans('Tt', 'Uu'))


if __name__ == '__main__':
    Ivt.main()

