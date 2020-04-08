#!/usr/bin/env python3

"""\
Amplify a DNA template using polymerase chain reaction (PCR).

Usage:
    pcr <template> <fwd_primer> <rev_primer>
           <num_reactions> <annealing_temp> <extension_time> [options]

Arguments:
    <template>
        The name of the template.

    <fwd_primer> <rev_primer>
        The name of the primers.

    <num_reactions>
        The number of reactions to set up.

    <annealing_temp>
        The annealing temperature for the PCR reaction (in °C).  I typically 
        use NEB's online "Tm Calculator" to determine this parameter.

    <extension_time>
        The length of the extension step in seconds.  The rule of thumb is 30 
        sec/kb, perhaps longer if you're amplifying a whole plasmid.

Options:
    -v --reaction-volume <μL>       [default: 10]
        The volume of the PCR reaction.  The recommended volumes for Q5 are 25
        and 50 μL.

    -m --master-mix <reagents>      [default: dna]
        Indicate which reagents should be included in the master mix.  The 
        following values are understood:

        dna:        The DNA template.
        fwd:        The forward primer.
        rev:        The reverse primer.
        primers:    Both primers, alias for 'fwd,rev'.

    -M --nothing-in-master-mix
        Don't include anything but water and polymerase in the master mix.  
        This is an alias for: -m ''

    -p --polymerase <name>          [default: q5]
        The name of the polymerase being used.  Different polymerases also have 
        different thermocycler parameters, as recommended by the manufacturer.  
        Currently, the following polymerases are supported:

        q5:         Q5 High-Fidelity DNA Polymerase (NEB)
        ssoadv:     SsoAdvanced™ Universal SYBR® Green Supermix (Biorad)
"""

import docopt
import autoprop
import stepwise
from math import sqrt
from inform import plural
from configurator import Config
from stepwise import UsageError

@autoprop
class Pcr:

    def __init__(self, template=None, primers=None, num_reactions=1):
        self.template = template
        self.primers = primers
        self.num_reactions = num_reactions

        self.polymerase = 'q5'
        self.reagents = None
        self.thermocycler_params = {}
        self.reaction_volume_uL = 10
        self.master_mix = ['dna']

    @classmethod
    def from_docopt(cls, *docopt_args, **docopt_kwargs):
        args = docopt.docopt(*docopt_args, **docopt_kwargs)

        pcr = cls()
        pcr.template = args['<template>']
        pcr.primers = args['<fwd_primer>'], args['<rev_primer>']
        pcr.polymerase = args['--polymerase']
        pcr.num_reactions = eval(args['<num_reactions>'])
        pcr.reaction_volume_uL = eval(args['--reaction-volume'])
        pcr.thermocycler_params['anneal_temp_C'] = args['<annealing_temp>']
        pcr.thermocycler_params['extend_time_s'] = float(args['<extension_time>'])
        pcr.master_mix = [] if args['--nothing-in-master-mix'] else [
                x.strip()
                for x in args['--master-mix'].split(',')
        ]
        return pcr

    def get_config(self):
        config = Config(stepwise.load_config()['molbio']['pcr'].data)
        config.merge(config['polymerases'][self.polymerase].data)
        config.merge(self.thermocycler_params)
        if self.reagents:
            config.merge({'reagents': self.reagents})

        return config

    def get_reaction(self):
        config = self.config

        def require_reagent(pcr, reagent):
            if reagent not in pcr:
                raise UsageError(f"reagent table for polymerase {self.polymerase!r} missing {reagent!r}.")

        pcr = stepwise.MasterMix.from_text(config.reagents)
    
        require_reagent(pcr, 'water')
        require_reagent(pcr, 'template DNA')
        require_reagent(pcr, 'forward primer')
        require_reagent(pcr, 'reverse primer')
        
        pcr.num_reactions = self.num_reactions
        pcr.hold_ratios.volume = self.reaction_volume_uL, 'µL'
        pcr.extra_volume = '10 µL'

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
            return f"{p[f'{step}_temp_C']}°C for {time(p[f'{step}_time_s'])}"

        three_step = not p.get('two_step', False) and has_step(p, 'extend')

        thermocycler_steps = [
                f"- {p['initial_denature_temp_C']}°C for {time(p['initial_denature_time_s'])}",
                f"- Repeat {p['num_cycles']}x:",
                f"  - {step(p, 'denature')}",
                f"  - {step(p, 'anneal')}",
        ]
        if three_step:
            thermocycler_steps += [
                f"  - {step(p, 'extend')}",
            ]

        if has_step(p, 'final_extend'):
            thermocycler_steps += [
                f"- {step(p, 'final_extend')}",
            ]

        if has_step(p, 'melt_curve', 'low_temp high_temp temp_step time_step'.split()):
            thermocycler_steps += [
                f"- {p['melt_curve_low_temp_C']}-{p['melt_curve_high_temp_C']}°C in {time(p['melt_curve_time_step_s'])} steps of {p['melt_curve_temp_step_C']}°C",
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
Dilute 1 µL twice into {sqrt(1000/x.value):.1g}*sqrt(DNA) µL
"""
        if primer_mix:
            protocol += f"""\
Prepare 10x primer mix [1]:

{primer_mix}
"""

        footnotes = list(protocol.footnotes.keys())
        if primer_mix: footnotes.remove(1)
        footnotes = ','.join(str(x) for x in footnotes)

        protocol += f"""\
Setup {plural(pcr.num_reactions):# PCR reaction/s} and 1 negative control{f' [{footnotes}]' if footnotes else ''}:

{pcr}
"""

        protocol += f"""\
Run the following thermocycler protocol:

{thermocycler}
"""

        return protocol

if __name__ == '__main__':
    pcr = Pcr.from_docopt(__doc__)
    print(pcr.protocol)

