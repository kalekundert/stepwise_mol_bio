#!/usr/bin/env python3

import stepwise, appcli, autoprop
from appcli import DocoptConfig
from stepwise import UsageError, StepwiseConfig, PresetConfig, pl, ul
from stepwise_mol_bio import Main, ConfigError, merge_dicts
from inform import plural

def parse_num_samples(name):
    try:
        return int(name)
    except ValueError:
        return len(name.strip(',').split(','))

def parse_sample_name(name):
    try:
        int(name)
        return None
    except ValueError:
        return name

@autoprop
class Gel(Main):
    """\
Load, run and stain PAGE gels

Usage:
    gel <preset> <samples> [options]

<%! from stepwise_mol_bio import hanging_indent %>\
Arguments:
    <preset>
        What kind of gel to run.  The following presets are available:

        ${hanging_indent(app.preset_briefs, 8*' ')}

    <samples>
        The names of the samples to run, separated by commas.  This can also be 
        a number, which will be taken as the number of samples to run.

Options:
    -p --percent <number>
        The percentage of polyacrylamide/agarose in the gel being run.

    -a --additive <str>
        An extra component to include in the gel itself, e.g. 1x EtBr.

    -b --buffer <str>
        The buffer to run the gel in, e.g. TAE.

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

    __config__ = [
            DocoptConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.gel'),
    ]
    preset_briefs = appcli.config_attr()

    presets = appcli.param(
            appcli.Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            appcli.Key(DocoptConfig, '<preset>'),
    )
    title = appcli.param(
            appcli.Key(PresetConfig, 'title'),
            default='electrophoresis',
    )
    num_samples = appcli.param(
            appcli.Key(DocoptConfig, '<samples>', cast=parse_num_samples),
            ignore=None,
            default=1,
    )
    sample_name = appcli.param(
            appcli.Key(DocoptConfig, '<samples>', cast=parse_sample_name),
            default=None,
    )
    sample_mix = appcli.param(
            appcli.Key(PresetConfig, 'sample_mix'),
            default=None,
    )
    sample_conc = appcli.param(
            appcli.Key(DocoptConfig, '--sample-conc'),
            appcli.Key(PresetConfig, 'sample_conc'),
            cast=float,
            default=None,
    )
    sample_volume_uL = appcli.param(
            appcli.Key(DocoptConfig, '--sample-volume'),
            appcli.Key(PresetConfig, 'sample_volume_uL'),
            cast=float,
            default=None,
    )
    ladder_name = appcli.param(
            appcli.Key(PresetConfig, 'ladder'),
            default=None,
    )
    ladder_volume_uL = appcli.param(
            appcli.Key(PresetConfig, 'ladder_volume_uL'),
    )
    mix_volume_uL = appcli.param(
            appcli.Key(DocoptConfig, '--mix-volume'),
            appcli.Key(PresetConfig, 'mix_volume_uL'),
            cast=float,
            default=None,
    )
    mix_extra_percent = appcli.param(
            appcli.Key(DocoptConfig, '--mix-extra'),
            appcli.Key(PresetConfig, 'mix_extra_percent'),
            cast=float,
            default=50,
    )
    incubate_temp_C = appcli.param(
            appcli.Key(DocoptConfig, '--incubate-temp'),
            appcli.Key(PresetConfig, 'incubate_temp_C'),
            cast=float,
    )
    incubate_time_min = appcli.param(
            appcli.Key(DocoptConfig, '--incubate-time'),
            appcli.Key(PresetConfig, 'incubate_time_min'),
            cast=int,
    )
    gel_type = appcli.param(
            appcli.Key(PresetConfig, 'gel_type'),
    )
    gel_percent = appcli.param(
            appcli.Key(DocoptConfig, '--percent'),
            appcli.Key(PresetConfig, 'gel_percent'),
    )
    gel_additive = appcli.param(
            appcli.Key(DocoptConfig, '--additive'),
            appcli.Key(PresetConfig, 'gel_additive'),
            default=None,
    )
    gel_buffer = appcli.param(
            appcli.Key(DocoptConfig, '--buffer'),
            appcli.Key(PresetConfig, 'gel_buffer'),
    )
    load_volume_uL = appcli.param(
            appcli.Key(DocoptConfig, '--load-volume'),
            appcli.Key(PresetConfig, 'load_volume_uL'),
            cast=float,
    )
    run_volts = appcli.param(
            appcli.Key(DocoptConfig, '--run-volts'),
            appcli.Key(PresetConfig, 'run_volts'),
            cast=float,
    )
    run_time_min = appcli.param(
            appcli.Key(DocoptConfig, '--run-time'),
            appcli.Key(PresetConfig, 'run_time_min'),
            cast=int,
    )
    stain = appcli.param(
            appcli.Key(DocoptConfig, '--stain'),
            appcli.Key(DocoptConfig, '--no-stain', cast=lambda x: None),
            appcli.Key(PresetConfig, 'stain'),
            default=None,
    )

    def __init__(self, preset, num_samples=None):
        self.preset = preset
        self.num_samples = num_samples

    def get_protocol(self):
        p = stepwise.Protocol()

        def both_or_neither(key1, key2):
            has_key1 = has_key2 = True

            try: value1 = getattr(self, key1)
            except AttributeError: has_key1 = False

            try: value2 = getattr(self, key2)
            except AttributeError: has_key2 = False

            if has_key1 and not has_key2:
                raise ConfigError(f"specified {key1!r} but not {key2!r}")
            if has_key2 and not has_key1:
                raise ConfigError(f"specified {key2!r} but not {key1!r}")

            if has_key1 and has_key2:
                return value1, value2
            else:
                return False

        if self.sample_mix:
            mix = stepwise.MasterMix.from_text(self.sample_mix)
            mix.num_reactions = self.num_samples
            mix.extra_percent = self.mix_extra_percent
            mix['sample'].name = self.sample_name

            if x := self.sample_conc:
                stock_conc = mix['sample'].stock_conc
                if stock_conc is None:
                    raise ConfigError(f"can't change sample stock concentration, no initial concentration specified.")
                mix['sample'].hold_conc.stock_conc = x, stock_conc.unit

            if x := self.sample_volume_uL:
                mix['sample'].volume = x, 'µL'

            if x := self.mix_volume_uL:
                mix.hold_ratios.volume = x, 'µL'

            mix.fix_volumes('sample')

            p += pl(
                    f"Prepare {plural(self.num_samples):# sample/s} for {self.title}:",
                    mix,
            )
            if x := both_or_neither('incubate_temp_C', 'incubate_time_min'):
                temp_C, time_min = x
                p.steps[-1] += ul(
                        f"Incubate at {temp_C:g}°C for {time_min:g} min."
                )
            
        additive = f" with {x}" if (x := self.gel_additive) else ""
        p += pl("Run a gel:", ul(
            f"Gel: {self.gel_percent}% {self.gel_type}{additive}",
            f"Buffer: {self.gel_buffer}",
            f"Ladder: {self.ladder_volume_uL:g} µL {self.ladder_name}"
                if self.ladder_name else None,
            f"Samples: {self.load_volume_uL:g} µL",
            f"Run: {self.run_volts:g}V for {self.run_time_min:g} min",
        ))

        if x := self.stain:
            p += stepwise.load(x)

        return p

if __name__ == '__main__':
    Gel.main()
