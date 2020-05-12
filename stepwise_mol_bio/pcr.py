#!/usr/bin/env python3

"""\
Amplify a DNA template using polymerase chain reaction (PCR).

Usage:
    pcr <template> <fwd_primer> <rev_primer> (-a <°C>) (-x <sec>) [options]

Arguments:
    <template>
        The name of the template(s) to amplify.

    <fwd_primer> <rev_primer>
        The name of the primers.

Options:
    -n --num-reactions <int>        [default: {pcr.num_reactions}]
        The number of reactions to set up.

    -v --reaction-volume <μL>
        The volume of the PCR reaction.

    -m --master-mix <reagents>      [default: {master_mix}]
        Indicate which reagents should be included in the master mix.  The 
        following values are understood:

        dna:      The DNA template.
        fwd:      The forward primer.
        rev:      The reverse primer.
        primers:  Both primers; alias for 'fwd,rev'.

    -M --nothing-in-master-mix
        Don't include anything but water and polymerase in the master mix.  
        This is an alias for: -m ''

    -p --preset <name>              [default: {pcr.preset}]
        The default reaction and thermocycler parameters to use.  The following 
        sets of parameters are currently available:

        {presets}

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
        
    --anneal-time <sec>
        The duration of the annealing step in seconds, e.g. 20.

    --extend-temp <°C>
        The temperature of the extension step in °C, e.g. 72.
        
    -x --extend-time <sec>
        The duration of the extension step in seconds.  The rule of thumb is
        30 sec/kb, perhaps longer if you're amplifying a whole plasmid.

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

import autoprop
import stepwise
from math import sqrt
from numbers import Real
from inform import plural, indent
from stepwise import UsageError
from stepwise_mol_bio import Main, Presets

CONFIG = stepwise.load_config()['molbio']['pcr']
PRESETS = Presets.from_config('molbio.pcr.presets')

@autoprop
class Pcr(Main):
    preset = CONFIG['default_preset']
    master_mix = {'dna'}
    num_reactions = 1
    reaction_volume_uL = None

    def __init__(self, template=None, primers=None):
        self.template = template
        self.primers = primers

        self.reagents = None
        self.thermocycler_params = {}

    @classmethod
    def from_docopt(cls, args):
        pcr = cls()
        pcr.template = args['<template>']
        pcr.primers = args['<fwd_primer>'], args['<rev_primer>']
        pcr.num_reactions = int(eval(args['--num-reactions']))
        pcr.preset = args['--preset']
        pcr.master_mix = [] if args['--nothing-in-master-mix'] else [
                x.strip()
                for x in args['--master-mix'].split(',')
        ]
        if x := args['--reaction-volume']:
            pcr.reaction_volume_uL = float(eval(x))

        def temp(x):
            try:
                return float(x)
            except ValueError:
                return x

        def time(x):
            return float(x)

        thermocycler_keys = [
                ('num-cycles', int),
                ('initial-denature-temp', temp),
                ('initial-denature-time', time),
                ('denature-temp', temp),
                ('denature-time', time),
                ('anneal-temp', temp),
                ('anneal-time', time),
                ('extend-temp', temp),
                ('extend-time', time),
                ('final-extend-temp', temp),
                ('final-extend-time', time),
                ('hold-temp', temp),
                ('melt-curve-low-temp', temp),
                ('melt-curve-high-temp', temp),
                ('melt-curve-temp-step', temp),
                ('melt-curve-time-step', time),
                ('two-step', bool),
                ('qpcr', bool),
        ]
        for key, parser in thermocycler_keys:
            arg_key = f'--{key}'
            if args[arg_key] in (None, False):
                continue

            param_parts = key.split('-')
            param_unit = {'temp': '_C', 'time': '_s'}
            param_key = '_'.join(param_parts) + param_unit.get(param_parts[-1], '')

            pcr.thermocycler_params[param_key] = parser(args[arg_key])

        return pcr

    def get_config(self):
        from configurator import Config

        config = Config(CONFIG.data)
        config.merge(PRESETS.load(self.preset))
        config.merge(self.thermocycler_params)

        if self.reagents:
            config.merge({'reagents': self.reagents})

        return config

    def get_reaction(self):
        config = self.config

        def require_reagent(pcr, reagent):
            if reagent not in pcr:
                raise UsageError(f"reagent table for preset {self.preset!r} missing {reagent!r}.")

        pcr = stepwise.MasterMix.from_text(config.reagents)
    
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
            pcr[p].hold_conc.stock_conc = config['primer_stock_uM'], 'µM'
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
            primer_mix.volume = '10 µL'

            for p in use_primer_mix:
                primer_mix[p].name = pcr[p].name
                primer_mix[p].stock_conc = pcr[p].stock_conc
                primer_mix[p].volume = pcr[p].volume
                primer_mix[p].hold_stock_conc.conc *= 10
                del pcr[p]

        return pcr, primer_mix

    def get_thermocycler_protocol(self):
        p = self.config

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

        def has_step(p, step, params=['temp_C', 'time_s']):
            return all((f'{step}_{param}' in p) for param in params)

        def step(p, step):
            return f"{temp(p[f'{step}_temp_C'])} for {time(p[f'{step}_time_s'])}"

        three_step = not p.get('two_step', False) and has_step(p, 'extend')

        thermocycler_steps = [
                f"- {step(p, 'initial_denature')}",
                f"- Repeat {p['num_cycles']}x:",
                f"  - {step(p, 'denature')}",
                f"  - {step(p, 'anneal')}",
        ]
        if three_step:
            thermocycler_steps += [
                f"  - {step(p, 'extend')}",
            ]

        if p.get('qpcr'):
            thermocycler_steps += [
                f"  - Measure fluorescence",
            ]

        if has_step(p, 'final_extend'):
            thermocycler_steps += [
                f"- {step(p, 'final_extend')}",
            ]

        if has_step(p, 'melt_curve', 'low_temp_C high_temp_C temp_step_C time_step_s'.split()):
            thermocycler_steps += [
                f"- {p['melt_curve_low_temp_C']}-{p['melt_curve_high_temp_C']}°C in {time(p['melt_curve_time_step_s'])} steps of {p['melt_curve_temp_step_C']}°C",
            ]

        if p.get('qpcr'):
            thermocycler_steps += [
                f"  - Measure fluorescence",
            ]

        if 'hold' in p:
            thermocycler_steps += [
                f"- {p['hold_temp_C']}°C hold",
            ]

        return '\n'.join(thermocycler_steps)

    def get_protocol(self):
        protocol = stepwise.Protocol()
        config = self.config

        pcr, primer_mix = self.reaction
        thermocycler = self.thermocycler_protocol

        protocol.footnotes[1] = f"""\
For resuspending lyophilized primers:
{config.primer_stock_uM} µM = {1e3 / config.primer_stock_uM:g} µL/nmol
"""
        if x := pcr['template DNA'].stock_conc:
            protocol.footnotes[2] = f"""\
For diluting template DNA to {x}:
Dilute 1 µL twice into {sqrt(1000/x.value):.1g}*sqrt([DNA]) µL
"""
        if primer_mix:
            protocol += f"""\
Prepare 10x primer mix [1]:

{primer_mix}
"""

        footnotes = list(protocol.footnotes.keys())
        if primer_mix: footnotes.remove(1)
        footnotes = ','.join(str(x) for x in footnotes)
        title = 'qPCR' if self.config.get('qpcr') else 'PCR'

        protocol += f"""\
Setup {plural(pcr.num_reactions):# {title} reaction/s}{f' [{footnotes}]' if footnotes else ''}:

{pcr}

{'- Use any extra master mix as a negative control.' if pcr.num_reactions > 1 else ''}
"""

        protocol += f"""\
Run the following thermocycler protocol:

{thermocycler}
"""

        return protocol

__doc__ = __doc__.format(
        pcr=Pcr,
        master_mix=','.join(Pcr.master_mix),
        presets=indent(PRESETS.format_briefs(), 8*' ', first=-1),
)

if __name__ == '__main__':
    Pcr.main(__doc__)

