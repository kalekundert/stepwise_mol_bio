#!/usr/bin/env python3

"""\
Load, run and stain PAGE gels

Usage:
    gel <preset> <samples> [options]

Arguments:
    <preset>
        What kind of gel to run.  The following presets are available:

        {presets}

    <samples>
        The names of the samples to run, separated by commas.  This can also be 
        a number, which will be taken as the number of samples to run.

Options:
    -p --percent <int>
        The percentage of polyacrylamide/agarose in the gel being run.

    -a --additive <str>
        An extra component to include in the gel itself, e.g. 1x EtBr.

    -c --sample-conc <value>
        The concentration of the sample.  This will be used to scale how much 
        sample is mixed with loading buffer, with the goal of mixing the same 
        quantity of material specified in the preset.  In order to use this 
        option, the preset must specify a sample concentration.  The units of 
        that concentration will be used for this concentration.

    -v --sample-volume <µL>
        The volume of sample to mix with loading buffer, in µL.  This does not 
        scale the concentration, and may increase or decrease the amount of 
        sample loaded relative to what's specified in the preset.

    --mix-volume <µL>
        The volume of the sample/loading buffer mix to prepare for each sample.  
        For example, if you want to run two gels, but the preset only makes 
        enough mix for one, use this option to make more.

    --mix-extra <percent>
        How much extra sample/loading buffer mix to make.

    --incubate-temp <°C>
        What temperature to incubate the sample/loading buffer at before 
        loading it onto the gel.  The incubation step will be skipped if 
        neither `--incubate-temp` nor `--incubate-time` are specified (either 
        on the command-line or via the preset).

    --incubate-time <min>
        How long to incubate the sample/loading buffer at the specified 
        temperature before loading it onto the gel.  The incubation step will 
        be skipped if neither `--incubate-temp` nor `--incubate-time` are 
        specified (either on the command-line or via the preset).

    -l --load-volume <µL>
        The volume of the sample/loading buffer mix to load onto the gel.

    --run-volts <V>
        The voltage to run the gel at.

    -r --run-time <min>
        How long to run the gel, in minutes.

    -s --stain <command>
        The name (and arguments) of a protocol describing how to stain the gel.  
        For example, this could be 'gelred' or 'coomassie -f'.
        
    -S --no-stain
        Leave off the staining step of the protocol.
        
"""

import stepwise
import autoprop
from inform import indent
from stepwise_mol_bio import Main, Presets, ConfigError

# Incorporate information from the config file into the usage text.
PRESETS = Presets.from_config('molbio.gel.presets')
PRESETS_DOC  = PRESETS.format_briefs("{gel_percent}% {gel_type}")
__doc__ = __doc__.format(
        presets=indent(PRESETS_DOC, 8*' ', first=-1),
)

@autoprop
class Gel(Main):

    def __init__(self):
        self.preset = ""
        self.params = {}
        self.num_samples = None

    @classmethod
    def from_docopt(cls, args):
        gel = cls()
        gel.preset = args['<preset>']

        try:
            gel.params['num_samples'] = int(args['<samples>'])
        except ValueError:
            gel.params['sample_name'] = name = args['<samples>']
            gel.params['num_samples'] = len(name.strip(',').split(','))

        keys = [
                ('--percent', 'gel_percent', str),
                ('--additive', 'gel_additive', str),
                ('--sample-conc', 'sample_conc', float),
                ('--sample-volume', 'sample_volume_uL', float),
                ('--mix-volume', 'mix_volume_uL', float),
                ('--mix-extra', 'mix_extra_percent', float),
                ('--incubate-temp', 'incubate_temp_C', float),
                ('--incubate-time', 'incubate_time_min', float),
                ('--load-volume', 'load_volume_uL', float),
                ('--run-volts', 'run_volts', float),
                ('--run-time', 'run_time_min', float),
                ('--stain', 'stain', str),
        ]
        for arg_key, param_key, parser in keys:
            if args[arg_key] is not None:
                gel.params[param_key] = parser(args[arg_key])

        if args['--no-stain']:
            gel.params['stain'] = None

        return gel

    def get_config(self):
        preset = PRESETS.load(self.preset)
        return {**preset, **self.params}

    def get_protocol(self):
        p = stepwise.Protocol()
        c = self.config

        def both_or_neither(c, key1, key2):
            has_key1 = has_key2 = True

            try: value1 = c[key1]
            except KeyError: has_key1 = False

            try: value2 = c[key2]
            except KeyError: has_key2 = False

            if has_key1 and not has_key2:
                raise cError(f"specified {key1!r} but not {key2!r}")
            if has_key2 and not has_key1:
                raise cError(f"specified {key2!r} but not {key1!r}")

            if has_key1 and has_key2:
                return value1, value2
            else:
                return False

        if x := c['sample_mix']:
            mix = stepwise.MasterMix.from_text(x)
            mix.num_reactions = c.get('num_samples', 1)
            mix.extra_percent = c.get('mix_extra_percent', 50)
            mix['sample'].name = c.get('sample_name')

            if y := c.get('sample_conc'):
                stock_conc = mix['sample'].stock_conc
                if stock_conc is None:
                    raise ConfigError(f"can't change sample stock concentration, no initial concentration specified.")
                mix['sample'].hold_conc.stock_conc = y, stock_conc.unit

            if y := c.get('sample_volume_uL'):
                mix['sample'].volume = y, 'µL'

            if y := c.get('mix_volume_uL'):
                mix.hold_ratios.volume = y, 'µL'

            incubate_step = ""

            p += f"""\
Prepare samples for {c.get('title', 'electrophoresis')}:

{mix}
"""
            if y := both_or_neither(c, 'incubate_temp_C', 'incubate_time_min'):
                temp_C, time_min = y
                p.steps[-1] += f"""\

- Incubate at {temp_C}°C for {time_min} min.
"""
            
        additive = f" with {c['gel_additive']}" if c.get('gel_additive') else ''
        p += f"""\
Run a gel:

- Use a {c['gel_percent']}% {c['gel_type']} gel{additive}.
- Load {c['load_volume_uL']} µL of each sample.
- Run at {c['run_volts']}V for {c['run_time_min']} min.
        """

        if x := c.get('stain'):
            p += stepwise.load(x)

        return p

if __name__ == '__main__':
    Gel.main(__doc__)
