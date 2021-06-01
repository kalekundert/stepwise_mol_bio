#!/usr/bin/env python3

import stepwise, appcli, autoprop
from appcli import Key, Method, DocoptConfig
from stepwise import UsageError, StepwiseConfig, PresetConfig, pl, ul, dl
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

@autoprop.cache
class Gel(Main):
    """\
Load, run, and stain gels.

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

    -M --no-mix
        Don't describe how to prepare the sample/loading buffer mix.

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

    --prerun-volts
        The voltage to pre-run the gel at.  By default, this will be the same 
        as the voltage the gel is run at.

    --prerun-time
        How long to prerun the gel, in minutes.

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
        Don't describe how to stain/visualize the gel.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.gel.presets:
        Named groups of default reaction parameters.  Typically each preset 
        corresponds to a particular kit or protocol.  See below for the various 
        settings that can be specified in each preset.

    molbio.gel.presets.<name>.title:
        How to briefly describe the gel in the protocol.  The default is 
        "electrophoresis".

    molbio.gel.presets.<name>.inherit:
        Copy all settings from another preset.  This can be used to make small 
        tweaks to a gel protocol, e.g. "SDS PAGE at lower-than-usual voltage".

    molbio.gel.presets.<name>.sample_mix:
        A table describing how to prepare the samples for the gel, in the 
        format understood by `stepwise.MasterMix.from_string()`.  The table 
        must contain one reagent named "sample".  This will be replaced with 
        the actual sample name specified on the command line, if possible.  If 
        no table is specified, the sample mix step will be left out of the 
        protocol.

    molbio.gel.presets.<name>.sample_conc:
        The default value for the `--sample-conc` option.

    molbio.gel.presets.<name>.sample_volume_uL:
        The default value for the `--sample-volume` option.

    molbio.gel.presets.<name>.ladder:
        The name of the ladder to use with this gel.

    molbio.gel.presets.<name>.ladder_volume_uL:
        How much ladder to load, in µL.

    molbio.gel.presets.<name>.mix_volume_uL:
        The default value for the `--mix-volume` option.
        
    molbio.gel.presets.<name>.mix_extra_percent:
        The default value for the `--mix-extra` option.

    molbio.gel.presets.<name>.incubate_temp_C:
        The default value for the `--incubate-temp` option.

    molbio.gel.presets.<name>.incubate_time_min:
        The default value for the `--incubate-time` option.

    molbio.gel.presets.<name>.gel_type:
        What kind of gel to use, e.g. "Bis-Tris/MES SDS PAGE", "TBE/urea PAGE", 
        "TAE/agarose", etc.  Don't include the gel percentage here; use the 
        `gel_percent` setting for that.

    molbio.gel.presets.<name>.gel_percent:
        The default value for the `--percent` option.

    molbio.gel.presets.<name>.gel_additive:
        The default value for the `--additive` option.

    molbio.gel.presets.<name>.gel_buffer:
        The default value for the `--buffer` option.
        
    molbio.gel.presets.<name>.load_volume_uL:
        The default value for the `--load-volume` option.

    molbio.gel.presets.<name>.prerun_volts:
        The default value for the `--prerun-volts` option.

    molbio.gel.presets.<name>.prerun_time_min:
        The default value for the `--prerun-time` option.

    molbio.gel.presets.<name>.run_volts:
        The default value for the `--run-volts` option.

    molbio.gel.presets.<name>.run_time_min:
        The default value for the `--run-time` option.

    molbio.gel.presets.<name>.stain:
        The default value for the `--stain` option.  If unspecified, there will 
        be no staining step by default.

    molbio.gel.presets.<name>.protocol_link:
        A hyperlink to an online description of the protocol, e.g. from the gel 
        manufacturer.  This link will be included as a footnote.
"""

    __config__ = [
            DocoptConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.gel'),
    ]
    preset_briefs = appcli.config_attr()
    config_paths = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '<preset>'),
    )
    title = appcli.param(
            Key(PresetConfig, 'title'),
            default='electrophoresis',
    )
    num_samples = appcli.param(
            Key(DocoptConfig, '<samples>', cast=parse_num_samples),
            ignore=None,
            default=1,
    )
    sample_name = appcli.param(
            Key(DocoptConfig, '<samples>', cast=parse_sample_name),
            default=None,
    )
    sample_mix_str = appcli.param(
            Key(DocoptConfig, '--no-mix', cast=lambda x: None),
            Key(PresetConfig, 'sample_mix'),
            default=None,
    )
    sample_conc = appcli.param(
            Key(DocoptConfig, '--sample-conc'),
            Key(PresetConfig, 'sample_conc'),
            cast=float,
            default=None,
    )
    sample_volume_uL = appcli.param(
            Key(DocoptConfig, '--sample-volume'),
            Key(PresetConfig, 'sample_volume_uL'),
            cast=float,
            default=None,
    )
    ladder_name = appcli.param(
            Key(PresetConfig, 'ladder'),
            default=None,
    )
    ladder_volume_uL = appcli.param(
            Key(PresetConfig, 'ladder_volume_uL'),
    )
    mix_volume_uL = appcli.param(
            Key(DocoptConfig, '--mix-volume'),
            Key(PresetConfig, 'mix_volume_uL'),
            cast=float,
            default=None,
    )
    mix_extra_percent = appcli.param(
            Key(DocoptConfig, '--mix-extra'),
            Key(PresetConfig, 'mix_extra_percent'),
            cast=float,
            default=50,
    )
    incubate_temp_C = appcli.param(
            Key(DocoptConfig, '--incubate-temp'),
            Key(PresetConfig, 'incubate_temp_C'),
            cast=float,
    )
    incubate_time_min = appcli.param(
            Key(DocoptConfig, '--incubate-time'),
            Key(PresetConfig, 'incubate_time_min'),
            cast=int,
    )
    gel_type = appcli.param(
            Key(PresetConfig, 'gel_type'),
    )
    gel_percent = appcli.param(
            Key(DocoptConfig, '--percent'),
            Key(PresetConfig, 'gel_percent'),
    )
    gel_additive = appcli.param(
            Key(DocoptConfig, '--additive'),
            Key(PresetConfig, 'gel_additive'),
            default=None,
    )
    gel_buffer = appcli.param(
            Key(DocoptConfig, '--buffer'),
            Key(PresetConfig, 'gel_buffer'),
    )
    load_volume_uL = appcli.param(
            Key(DocoptConfig, '--load-volume'),
            Key(PresetConfig, 'load_volume_uL'),
            cast=float,
    )
    prerun_volts = appcli.param(
            Key(DocoptConfig, '--prerun-volts'),
            Key(PresetConfig, 'prerun_volts'),
            Method(lambda self: self.run_volts),
            cast=float,
    )
    prerun_time_min = appcli.param(
            Key(DocoptConfig, '--prerun-time'),
            Key(PresetConfig, 'prerun_time_min'),
            cast=int,
            default=None,
    )
    run_volts = appcli.param(
            Key(DocoptConfig, '--run-volts'),
            Key(PresetConfig, 'run_volts'),
            cast=float,
    )
    run_time_min = appcli.param(
            Key(DocoptConfig, '--run-time'),
            Key(PresetConfig, 'run_time_min'),
            cast=int,
    )
    stain = appcli.param(
            Key(DocoptConfig, '--stain'),
            Key(DocoptConfig, '--no-stain', cast=lambda x: None),
            Key(PresetConfig, 'stain'),
            default=None,
    )
    protocol_link = appcli.param(
            Key(PresetConfig, 'protocol_link'),
            default=None,
    )

    def __init__(self, preset, num_samples=None):
        self.preset = preset
        self.num_samples = num_samples

    def get_protocol(self):
        p = stepwise.Protocol()

        if self.sample_mix:
            p += self.prep_step

        p += self.run_step
            
        if self.stain:
            p += stepwise.load(self.stain)

        return p

    def del_protocol(self):
        pass

    def get_prep_step(self):

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

        s = pl(
                f"Prepare {plural(self.num_samples):# sample/s} for {self.title}:",
                self.sample_mix,
        )
        if x := both_or_neither('incubate_temp_C', 'incubate_time_min'):
            temp_C, time_min = x
            s += ul(
                    f"Incubate at {temp_C:g}°C for {time_min:g} min."
            )

        return s

    def del_prep_step(self):
        pass

    def get_run_step(self):
        p = stepwise.Protocol()
        additive = f" with {x}" if (x := self.gel_additive) else ""
        percent = x.replace('-', '–') if isinstance(x := self.gel_percent, str) else x
        p += pl(f"Run a gel{p.add_footnotes(self.protocol_link)}:", dl(
            ("gel", f"{percent}% {self.gel_type}{additive}"),
            ("buffer", f"{self.gel_buffer}"),
            ("ladder", self.ladder_name and f"{self.ladder_volume_uL:g} µL {self.ladder_name}"),
            ("samples", self.load_volume_uL and f"{self.load_volume_uL:.3g} µL/lane"),
            ("prerun", self.prerun_time_min and f"{self.prerun_volts:g}V for {self.prerun_time_min:g} min"),
            ("run", f"{self.run_volts:g}V for {self.run_time_min:g} min"),
        ))
        return p

    def del_run_step(self):
        pass

    def get_sample_mix(self):
        if not self.sample_mix_str:
            return None

        mix = stepwise.MasterMix.from_text(self.sample_mix_str)
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

        if mix.solvent:
            mix.fix_volumes('sample')

        return mix

    def del_sample_mix(self):
        pass

if __name__ == '__main__':
    Gel.main()
