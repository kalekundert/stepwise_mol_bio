#!/usr/bin/env python3

"""\
Usage:
    pcr.py <template> <fwd_primer> <rev_primer>
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
import stepwise
from stepwise import Q
from inform import plural

general_config = {
        'primer_conc': '100 µM',
}
reagent_tables = {
        'q5': """\
                Reagent              Stock    Volume  MM?
                ==============  ==========  ========  ===
                water                       to 10 µL  yes
                template DNA     100 pg/µL    0.2 µL
                forward primer       10 µM    0.5 µL
                reverse primer       10 µM    0.5 µL
                Q5 master mix           2x    5.0 µL  yes
                """,

        'ssoadv': """\
                Reagent                   Stock    Volume  MM?
                ====================  =========  ========  ===
                water                            to 20 µL  yes
                template DNA          100 pg/µL      1 µL
                forward primer            10 µM      1 µL
                reverse primer            10 µM      1 µL
                SsoAdvanced supermix         2x     10 µL  yes
                """,
}
thermocycler_protocols = {
        # This needs to be in a config file somewhere...
        'q5': {
            'initial_denature_temp': 98,
            'initial_denature_time': 30,
            'denature_temp': 98,
            'denature_time': 10,
            'anneal_temp': 60,
            'anneal_time': 20,
            'extend_temp': 72,
            'extend_time': 120,
            'final_extend_temp': 72,
            'final_extend_time': 120,
            'num_cycles': 35,
            'hold': 4,
        },
        'ssoadv': {
            'initial_denature_temp': 95,
            'initial_denature_time': 30,
            'denature_temp': 95,
            'denature_time': 10,
            'anneal_temp': 60,
            'anneal_time': 15,
            'melt_curve_low_temp': 65,
            'melt_curve_high_temp': 95,
            'melt_curve_temp_step': 0.5,
            'melt_curve_time_step': 5,
            'num_cycles': 40,
            'two_step': True,
        },
}

def make_thermocycler_steps(params):
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

    def has_step(protocol, step, params=['temp', 'time']):
        return all((f'{step}_{param}' in protocol) for param in params)

    p = params
    three_step = not p.get('two_step', False) and has_step(p, 'extend')

    thermocycler_steps = [
            f"- {p['initial_denature_temp']}°C for {time(p['initial_denature_time'])}",
            f"- Repeat {p['num_cycles']}x:",
            f"  - {p['denature_temp']}°C for {time(p['denature_time'])}",
            f"  - {p['anneal_temp']}°C for {time(p['anneal_time'])}",
    ]
    if three_step:
        thermocycler_steps += [
            f"  - {p['extend_temp']}°C for {time(p['extend_time'])}",
        ]

    if has_step(p, 'final_extend'):
        thermocycler_steps += [
            f"- {p['final_extend_temp']}°C for {time(p['final_extend_time'])}",
        ]

    if has_step(p, 'melt_curve', 'low_temp high_temp temp_step time_step'.split()):
        thermocycler_steps += [
            f"- {p['melt_curve_low_temp']}-{p['melt_curve_high_temp']}°C in {time(p['melt_curve_time_step'])} steps of {p['melt_curve_temp_step']}°C",
        ]

    if 'hold' in p:
        thermocycler_steps += [
            f"- {p['hold']}°C hold",
        ]

    return '\n'.join(thermocycler_steps)


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    polymerase = args['--polymerase']
    master_mix = [] if args['--nothing-in-master-mix'] else [
            x.strip()
            for x in args['--master-mix'].split(',')
    ]

    # Setup the PCR reaction:
    pcr = stepwise.MasterMix.from_text(reagent_tables[polymerase])
    pcr.num_reactions = eval(args['<num_reactions>'])
    pcr.hold_ratios.volume = eval(args['--reaction-volume']), 'µL'
    pcr.extra_volume = '10 µL'

    pcr['water'].order = 1

    pcr['template DNA'].order = 2
    pcr['template DNA'].name = args['<template>']
    pcr['template DNA'].master_mix = 'dna' in master_mix

    use_primer_mix = []
    primer_names = 'forward primer', 'reverse primer'
    primer_abbrev = {
            'forward primer': 'fwd',
            'reverse primer': 'rev',
    }

    for p in primer_names:
        pcr[p].order = 3
        pcr[p].name = args[f'<{primer_abbrev[p]}_primer>']
        pcr[p].hold_conc.stock_conc = general_config['primer_conc']
        pcr[p].master_mix = (
                primer_abbrev[p] in master_mix or
                       'primers' in master_mix
        )

        primer_scale = pcr.scale if pcr[p].master_mix else 1
        primer_vol = primer_scale * pcr[p].volume

        if primer_vol < '0.5 µL':
            use_primer_mix.append(p)

    # Setup the primer mix (if needed):
    if use_primer_mix:
        pcr['primer mix'].order = 4
        pcr['primer mix'].stock_conc = '10x'
        pcr['primer mix'].volume = pcr.volume / 10
        pcr['primer mix'].master_mix = all(
                pcr[p].master_mix
                for p in use_primer_mix
        )

        primers = stepwise.MasterMix()
        primers.volume = '10 µL'

        for p in use_primer_mix:
            primers[p].name = pcr[p].name
            primers[p].stock_conc = pcr[p].stock_conc
            primers[p].volume = pcr[p].volume
            primers[p].hold_stock_conc.conc *= 10
            del pcr[p]

    # Setup the thermocycler protocol:
    thermocycler_params = thermocycler_protocols[polymerase]
    thermocycler_params['anneal_temp'] = args['<annealing_temp>']
    thermocycler_params['extend_time'] = args['<extension_time>']
    thermocycler = make_thermocycler_steps(thermocycler_params)

    # Make the protocol:
    protocol = stepwise.Protocol()

    if use_primer_mix:
        protocol += f"""\
Prepare 10x primer mixes [1]:

{primers}
"""

    protocol += f"""\
Setup {plural(pcr.num_reactions):# PCR reaction/s} and 1 negative control:

{pcr}
"""

    protocol += f"""\
Run the following thermocycler protocol:

{thermocycler}
"""

    protocol.footnotes[1] = f"""\
For resuspending lyophilized primers:
{general_config['primer_conc']} = {1e3 / Q(general_config['primer_conc']).value:g} µL/nmol
"""

    print(protocol)
