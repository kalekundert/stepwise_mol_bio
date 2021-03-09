#!/usr/bin/env python3

import stepwise, appcli, autoprop
from math import sqrt, ceil
from numbers import Real
from inform import plural, indent
from appcli import Key, DocoptConfig
from stepwise import UsageError, StepwiseConfig, PresetConfig, pl, ul, pre
from stepwise_mol_bio import Main, merge_dicts, comma_set, require_reagent
from more_itertools import first_true
from copy import deepcopy
from operator import not_

def temp(x):
    try:
        return float(x)
    except ValueError:
        return x

def time(x):
    return float(x)

@autoprop
class Pcr(Main):
    """\
Amplify a DNA template using polymerase chain reaction (PCR).

Usage:
    pcr <template> <fwd> <rev> (-a <°C>) (-x <sec> | -l <kb>) [options]

Arguments:
    <template>
        The name of the template(s) to amplify.

    <fwd> <rev>
        The names of the forward and reverse primers.

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>              [default: ${app.preset}]
        The default reaction and thermocycler parameters to use.  The following 
        sets of parameters are currently available:

        ${hanging_indent(app.preset_briefs, 8*' ')}

    -l --amplicon-length <bp>
        The length of the amplicon in base pairs (bp).  This can be used to 
        calculate an appropriate extension time.

    -n --num-reactions <int>        [default: ${app.num_reactions}]
        The number of reactions to set up.

    -v --reaction-volume <μL>
        The volume of the PCR reaction.

    -V --template-volume <µL>
        The volume of template to use in the reaction.  This overrides the 
        value specified by the preset.

    -T --template-stock <conc>
        The stock concentration of template to use in the reaction.  This 
        overrides the value specified by the preset, without affecting the 
        volume.  Include a unit, because none is implied.

    -m --master-mix <reagents>      [default: ${','.join(app.master_mix)}]
        Indicate which reagents should be included in the master mix.  The 
        following values are understood:

        dna:      The DNA template.
        fwd:      The forward primer.
        rev:      The reverse primer.
        primers:  Both primers; alias for 'fwd,rev'.

    -M --nothing-in-master-mix
        Don't include anything but water and polymerase in the master mix.  
        This is an alias for: -m ''

    --primer-stock <µM>             [default: ${app.primer_stock_uM}]
        The stock concentration of the primers in µM.  Both primers must have 
        the same concentration.

    -y --num-cycles <n>
        The number of denature/anneal/extend cycles to perform, e.g. 35.

    --initial-denature-temp <°C>
        The temperature of the initial denaturation step in °C, e.g. 95.

    --initial-denature-time <sec>
        The duration of the initial denaturation step in seconds, e.g. 30.

    --denature-temp <°C>
        The temperature of the denaturation step in °C, e.g. 95.

    --denature-time <sec>
        The duration of the denaturation step in seconds, e.g. 10.

    -a --anneal-temp <°C>
        The temperature of the annealing step in °C, e.g. 60.  This is 
        determined by the sequence of the primers.
        
    -g --anneal-temp-gradient <range>
        The range of annealing temperatures that should be tried in a gradient 
        PCR reaction.  The range will be centered at the indicated annealing 
        temperature, and the protocol will indicate the corresponding high and 
        low temperatures.

    --anneal-time <sec>
        The duration of the annealing step in seconds, e.g. 20.

    --extend-temp <°C>
        The temperature of the extension step in °C, e.g. 72.
        
    -x --extend-time <sec>
        The duration of the extension step in seconds.  The rule of thumb is
        30 sec/kb, perhaps longer if you're amplifying a whole plasmid.

    --no-round-extend-time 
        When calculating an extension time from a given amplicon length (e.g. 
        `-l`), do *not* round the result to the nearest 30 sec increment.  Note 
        that when an explicit extension time is given (e.g. `-x`), it is never 
        rounded.

    --final-extend-temp <°C>
        The temperature of the final extension step in °C, e.g. 72.

    --final-extend-time <sec>
        The duration of the annealing step in seconds, e.g. 120.

    --hold-temp <°C>
        The temperature in °C to hold the reaction at after it completes.

    --melt-curve-low-temp <°C>
        The temperature in °C at which to begin recording a melt curve,
        e.g. 65.  This is only relevant for qPCR protocols.

    --melt-curve-high-temp <°C>
        The temperature in °C at which to stop recording a melt curve,
        e.g. 95.  This is only relevant for qPCR protocols.

    --melt-curve-temp-step <°C>
        How much to increment the temperature in °C at each step in the melt 
        curve, e.g. 0.5°C.  This is only relevant for qPCR protocols.

    --melt-curve-time-step <sec>
        The duration in seconds of each step in the melt curve, e.g. 5.  This 
        is only relevant for qPCR protocols.

    --two-step
        Specify that the annealing and extension steps should be combined, e.g. 
        if the primer melting temperatures are very close to the extension 
        temperature or if the amplicon is very short.

    --qpcr
        Specify that this is a qPCR protocol, i.e. that fluorescence should be 
        measured after each thermocyler cycle.
"""

    __config__ = [
            DocoptConfig(usage_getter=lambda self: self.format_usage()),
            PresetConfig(),
            StepwiseConfig('molbio.pcr'),
    ]
    usage = appcli.config_attr()
    preset_briefs = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '--preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    template = appcli.param(
            Key(DocoptConfig, '<template>'),
    )
    fwd_primer = appcli.param(
            Key(DocoptConfig, '<fwd>'),
    )
    rev_primer = appcli.param(
            Key(DocoptConfig, '<rev>'),
    )
    primer_stock_uM = appcli.param(
            Key(DocoptConfig, '--primer-stock'),
            Key(PresetConfig, 'primer_stock_uM'),
            Key(StepwiseConfig, 'primer_stock_uM'),
            cast=float,
    )
    amplicon_length_bp = appcli.param(
            Key(DocoptConfig, '--amplicon-length'),
            cast=float,
            default=None,
    )
    num_reactions = appcli.param(
            Key(DocoptConfig, '--num-reactions'),
            cast=lambda x: int(eval(x)),
            default=1,
    )
    reaction_volume_uL = appcli.param(
            Key(DocoptConfig, '--reaction-volume'),
            Key(PresetConfig, 'reaction_volume_uL'),
            cast=lambda x: float(eval(x)),
            default=None,
    )
    template_volume_uL = appcli.param(
            Key(DocoptConfig, '--template-volume'),
            cast=float,
            default=None,
    )
    template_stock = appcli.param(
            Key(DocoptConfig, '--template-stock'),
            default=None,
    )
    master_mix = appcli.param(
            Key(DocoptConfig, '--nothing-in-master-mix', cast=lambda x: set()),
            Key(DocoptConfig, '--master-mix', cast=comma_set),
            Key(StepwiseConfig, 'master_mix', cast=comma_set),
            default={'dna'},
    )
    base_reaction = appcli.param(
            Key(PresetConfig, 'reagents'),
            cast=stepwise.MasterMix.from_text,
    )
    num_cycles = appcli.param(
            Key(DocoptConfig, '--num-cycles'),
            Key(PresetConfig, 'num_cycles'),
            cast=int,
    )
    initial_denature_temp_C = appcli.param(
            Key(DocoptConfig, '--initial-denature-temp'),
            Key(PresetConfig, 'initial_denature_temp_C'),
            cast=temp,
    )
    initial_denature_time_s = appcli.param(
            Key(DocoptConfig, '--initial-denature-time'),
            Key(PresetConfig, 'initial_denature_time_s'),
            cast=time,
    )
    denature_temp_C = appcli.param(
            Key(DocoptConfig, '--denature-temp'),
            Key(PresetConfig, 'denature_temp_C'),
            cast=temp,
    )
    denature_time_s = appcli.param(
            Key(DocoptConfig, '--denature-time'),
            Key(PresetConfig, 'denature_time_s'),
            cast=time,
    )
    anneal_temp_C = appcli.param(
            Key(DocoptConfig, '--anneal-temp'),
            Key(PresetConfig, 'anneal_temp_C'),
            cast=temp,
    )
    anneal_temp_gradient_C = appcli.param(
            Key(DocoptConfig, '--anneal-temp-gradient'),
            Key(PresetConfig, 'anneal_temp_gradient_C'),
            cast=float,
    )
    anneal_time_s = appcli.param(
            Key(DocoptConfig, '--anneal-time'),
            Key(PresetConfig, 'anneal_time_s'),
            cast=time,
    )
    extend_temp_C = appcli.param(
            Key(DocoptConfig, '--extend-temp'),
            Key(PresetConfig, 'extend_temp_C'),
            cast=temp,
    )
    extend_time_s_per_kb = appcli.param(
            Key(PresetConfig, 'extend_time_s_per_kb'),
            cast=float,
    )
    extend_time_func = appcli.param(
            Key(PresetConfig, 'extend_time_func'),
            cast=eval,
    )
    round_extend_time = appcli.param(
            Key(DocoptConfig, '--no-round-extend-time', cast=not_),
            Key(PresetConfig, 'round_extend_time'),
            Key(StepwiseConfig, 'round_extend_time'),
            default=True,
    )
    final_extend_temp_C = appcli.param(
            Key(DocoptConfig, '--final-extend-temp'),
            Key(PresetConfig, 'final_extend_temp_C'),
            cast=temp,
    )
    final_extend_time_s = appcli.param(
            Key(DocoptConfig, '--final-extend-time'),
            Key(PresetConfig, 'final_extend_time_s'),
            cast=time,
    )
    hold_temp_C = appcli.param(
            Key(DocoptConfig, '--hold-temp'),
            Key(PresetConfig, 'hold_temp_C'),
            cast=temp,
            default=None,
    )
    melt_curve_low_temp_C = appcli.param(
            Key(DocoptConfig, '--melt-curve-low-temp'),
            Key(PresetConfig, 'melt_curve_low_temp_C'),
            cast=temp,
    )
    melt_curve_high_temp_C = appcli.param(
            Key(DocoptConfig, '--melt-curve-high-temp'),
            Key(PresetConfig, 'melt_curve_high_temp_C'),
            cast=temp,
    )
    melt_curve_temp_step_C = appcli.param(
            Key(DocoptConfig, '--melt-curve-temp-step'),
            Key(PresetConfig, 'melt_curve_temp_step_C'),
            cast=temp,
    )
    melt_curve_time_step_s = appcli.param(
            Key(DocoptConfig, '--melt-curve-time-step'),
            Key(PresetConfig, 'melt_curve_time_step_s'),
            cast=temp,
    )
    two_step = appcli.param(
            Key(DocoptConfig, '--two-step'),
            Key(PresetConfig, 'two_step'),
            pick=first_true,
            default=False,
    )
    qpcr = appcli.param(
            Key(DocoptConfig, '--qpcr'),
            Key(PresetConfig, 'qpcr'),
            pick=first_true,
            default=False,
    )
    footnote = appcli.param(
            Key(PresetConfig, 'footnote'),
            default=None,
    )

    @appcli.param(
            Key(DocoptConfig, '--extend-time'),
            Key(PresetConfig, 'extend_time_s'),
            cast=time,
            default=None,
    )
    def extend_time_s(self, x):
        if x:
            return x

        bp = self.amplicon_length_bp

        def round(x):
            if not self.round_extend_time: return x
            if x <= 10: return 10
            if x <= 15: return 15
            return 30 * ceil(x / 30)

        if factor := self.extend_time_s_per_kb:
            return round(factor * bp / 1000)

        return round(self.extend_time_func(bp))


    def __init__(self, template, primers):
        self.template = template
        self.primers = primers

    def format_usage(self):
        return self.__doc__

    def get_primers(self):
        return self.fwd_primer, self.rev_primer

    def set_primers(self, primers):
        self.fwd_primer, self.rev_primer = primers

    def get_reaction(self):
        pcr = deepcopy(self.base_reaction)
    
        require_reagent(pcr, 'water')
        require_reagent(pcr, 'template DNA')
        require_reagent(pcr, 'forward primer')
        require_reagent(pcr, 'reverse primer')

        pcr.extra_volume = '10 µL'
        
        if self.num_reactions:
            pcr.num_reactions = self.num_reactions
        if self.reaction_volume_uL:
            pcr.hold_ratios.volume = self.reaction_volume_uL, 'µL'

        pcr['water'].order = 1

        pcr['template DNA'].order = 2
        pcr['template DNA'].name = self.template
        pcr['template DNA'].master_mix = 'dna' in self.master_mix

        if x := self.template_volume_uL:
            pcr['template DNA'].volume = x, 'µL'
        if x := self.template_stock:
            pcr['template DNA'].stock_conc = x

        # Setup the primers.  This is complicated because the primers might 
        # get split into their own mix, if the volumes that would be added 
        # to the PCR reaction are too small.

        primer_mix = None
        primer_abbrev = {
                'forward primer': 'fwd',
                'reverse primer': 'rev',
        }
        use_primer_mix = []

        for p, name in zip(primer_abbrev, self.primers or (None, None)):
            pcr[p].order = 3
            pcr[p].name = name
            pcr[p].hold_conc.stock_conc = self.primer_stock_uM, 'µM'
            pcr[p].master_mix = (
                    primer_abbrev[p] in self.master_mix or
                           'primers' in self.master_mix
            )

            primer_scale = pcr.scale if pcr[p].master_mix else 1
            primer_vol = primer_scale * pcr[p].volume

            if primer_vol < '0.5 µL':
                use_primer_mix.append(p)

        if use_primer_mix:
            pcr['primer mix'].order = 4
            pcr['primer mix'].stock_conc = '10x'
            pcr['primer mix'].volume = pcr.volume / 10
            pcr['primer mix'].master_mix = all(
                    pcr[p].master_mix
                    for p in use_primer_mix
            )

            primer_mix = stepwise.MasterMix()
            primer_mix.volume = pcr.volume

            for p in use_primer_mix:
                primer_mix[p].name = pcr[p].name
                primer_mix[p].stock_conc = pcr[p].stock_conc
                primer_mix[p].volume = pcr[p].volume
                primer_mix[p].hold_stock_conc.conc *= 10
                del pcr[p]

            primer_mix.hold_ratios.volume = '10 µL'

        return pcr, primer_mix

    def get_thermocycler_protocol(self):

        def temp(x):
            if isinstance(x, Real):
                return f'{x:g}°C'
            elif x[-1].isdigit():
                return f'{x}°C'
            else:
                return x

        def time(x):
            try:
                x = int(x)
            except ValueError:
                return x

            if x < 60:
                return f'{x}s'
            elif x % 60:
                return f'{x//60}m{x%60:02}'
            else:
                return f'{x//60} min'

        def has_step(step, params=['temp_C', 'time_s']):
            return all(
                hasattr(self, f'{step}_{param}')
                for param in params
            )

        def step(step):
            if not hasattr(self, f'{step}_temp_gradient_C'):
                t = getattr(self, f'{step}_temp_C')
            else:
                t_mid = getattr(self, f'{step}_temp_C')
                t_range = getattr(self, f'{step}_temp_gradient_C')
                t_low = round(t_mid - t_range / 2)
                t_high = round(t_low + t_range)
                t = f'{t_low}-{t_high}'

            return f"{temp(t)} for {time(getattr(self, f'{step}_time_s'))}"

        three_step = not self.two_step and has_step('extend')

        thermocycler_steps = [
                f"- {step('initial_denature')}",
                f"- Repeat {getattr(self, 'num_cycles')}x:",
                f"  - {step('denature')}",
                f"  - {step('anneal')}",
        ]
        if three_step:
            thermocycler_steps += [
                f"  - {step('extend')}",
            ]

        if self.qpcr:
            thermocycler_steps += [
                f"  - Measure fluorescence",
            ]

        if has_step('final_extend'):
            thermocycler_steps += [
                f"- {step('final_extend')}",
            ]

        if has_step('melt_curve', 'low_temp_C high_temp_C temp_step_C time_step_s'.split()):
            thermocycler_steps += [
                f"- {self.melt_curve_low_temp_C}-{self.melt_curve_high_temp_C}°C in {time(self.melt_curve_time_step_s)} steps of {self.melt_curve_temp_step_C}°C:",
            ]

        if self.qpcr:
            thermocycler_steps += [
                f"  - Measure fluorescence",
            ]

        if self.hold_temp_C:
            thermocycler_steps += [
                f"- {temp(self.hold_temp_C)} hold",
            ]

        return pre('\n'.join(thermocycler_steps))

    def get_protocol(self):
        protocol = stepwise.Protocol()
        pcr, primer_mix = self.reaction
        thermocycler = self.thermocycler_protocol
        footnotes = []

        # Primer mix (if applicable):

        footnotes.append(pl(
                f"For resuspending lyophilized primers:",
                f"{self.primer_stock_uM} µM = {1e3 / self.primer_stock_uM:g} µL/nmol",
                br='\n',
        ))

        if primer_mix:
            protocol += pl(
                    f"Prepare 10x primer mix{protocol.add_footnotes(*footnotes)}:",
                    primer_mix,
            )
            footnotes = []

        # PCR reaction setup:

        if self.footnote:
            # Before the footnote on resuspending the primers, if it hasn't 
            # been added to the protocol already.
            footnotes.insert(0, self.footnote)

        if x := pcr['template DNA'].stock_conc:
            footnotes.append(pl(
                    f"For diluting template DNA to {x}:",
                    f"Dilute 1 µL twice into {sqrt(1000/x.value):.1g}*sqrt([DNA]) µL",
                    br='\n',
            ))

        title = 'qPCR' if self.qpcr else 'PCR'
        instructions = ul()
        protocol += pl(
                f"Setup {plural(pcr.num_reactions):# {title} reaction/s}{protocol.add_footnotes(*footnotes)}:",
                pcr,
                instructions,
        )

        if pcr.volume > '50 µL':
            instructions += f"Split each reaction into {ceil(pcr.volume.value / 50)} tubes."
        if pcr.num_reactions > 1:
            instructions += "Use any extra master mix as a negative control."

        # Thermocycler protocol:

        protocol += pl(
                "Run the following thermocycler protocol:",
                thermocycler,
        )

        protocol.renumber_footnotes()
        return protocol

if __name__ == '__main__':
    Pcr.main()

