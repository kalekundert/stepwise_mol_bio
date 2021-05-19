#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import warn
from appcli import Key, DocoptConfig
from stepwise import StepwiseConfig, PresetConfig, Quantity, oxford_comma
from stepwise_mol_bio import Cleanup, format_sec, merge_dicts
from freezerbox import MakerArgsConfig, group_by_identity, parse_volume_uL
from more_itertools import always_iterable

def ng_uL(x):
    return Quantity(x, 'ng/µL')

@autoprop
class SpinCleanup(Cleanup):
    """\
Purify a PCR reaction using a silica spin column.

Usage:
    spin_cleanup [<preset>] [-s <µL>] [-d <buffer>] [-v <µL>]

<%! from stepwise_mol_bio import hanging_indent %>\
Arguments:
    <preset>                        [default: ${app.preset}]
        The default parameters to use.  Typically these correspond to 
        commercial kits:

        ${hanging_indent(app.preset_briefs, 8*' ')}

Options:
    -s --sample-volume <µL>
        The volume of the sample, in µL.

    -d --elute-buffer <name>
        The buffer to elute in.

    -v --elute-volume <µL>
        The volume of purified DNA/RNA to elute, in µL.  The default value 
        depends on the preset, but can usually be lowered to get more 
        concentrated product.  A warning will be displayed if the requested 
        volume is lower than the minimum recommended by the kit manufacturer.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.spin_cleanup.default_preset:
        The default value for the `--preset` option.

    molbio.spin_cleanup.presets:
        Named groups of default reaction parameters.  Typically each preset 
        corresponds to a particular kit or protocol.  See below for the various 
        settings that can be specified in each preset.

    molbio.spin_cleanup.presets.<name>.protocol_name
        How to refer to the whole protocol.  Commonly this is the name of the 
        spin column kit.

    molbio.spin_cleanup.presets.<name>.protocol_link
        A link (typically minified) to the complete protocol, e.g. as published 
        by the manufacturer of the columns.  This is not required, but if 
        specified, will be included in the protocol as a footnote.

    molbio.spin_cleanup.presets.<name>.column_name
        How to refer to the specific spin column used in the protocol.

    molbio.spin_cleanup.presets.<name>.spin_speed_g
        How fast to spin the column in each centrifugation step, in units of 
        g-force.

    molbio.spin_cleanup.presets.<name>.column_capacity_ug
        The maximum binding capacity of the column, in µg.  This information is 
        added to the protocol as a footnote.

    molbio.spin_cleanup.presets.<name>.sample_type
        How to generically refer to the sample in the protocol, e.g. "DNA".

    molbio.spin_cleanup.presets.<name>.sample_volume_uL
        The volume of sample to load on the column, in µL.  Alternatively, this 
        can be a dictionary with keys 'min' and/or 'max' specifying the minimum 
        and maximum allowed sample volumes, respectively.

    molbio.spin_cleanup.presets.<name>.bind_buffer
        The name(s) of the buffer(s) to use to bind the sample to column.  This 
        can be either a string or a list of strings.  Use a list to specify 
        that multiple buffers (e.g. binding buffer and ethanol) should be mixed 
        with the sample before it is loaded on the column.  If this option is a 
        list, the `bind_volume_uL` and `bind_volume_x` options must also be 
        lists of the same length (or left unspecified).
        
    molbio.spin_cleanup.presets.<name>.bind_volume_uL
        How much `bind_buffer` to use, in µL.  This can be either a number or a 
        list of numbers; see `bind_buffer` for more details.  This takes 
        precedence over the `bind_volume_x` setting.  

    molbio.spin_cleanup.presets.<name>.bind_volume_x
        How much `bind_buffer` to use, as a multiple of the sample volume.  
        This can be a number or a list of numbers; see `bind_buffer` for more 
        details.  This is superseded by the `bind_volume_uL` setting.  

    molbio.spin_cleanup.presets.<name>.bind_spin_sec
        How long to centrifuge the column during the bind step.

    molbio.spin_cleanup.presets.<name>.bind_vacuum
        Whether or not to use a vacuum manifold for the bind step.  The default 
        is False.  If True, the `bind_spin_sec` option is ignored.

    molbio.spin_cleanup.presets.<name>.pH_buffer
        The name of the buffer to use when adjusting the pH of the sample.

    molbio.spin_cleanup.presets.<name>.pH_volume_uL
        How much `pH_buffer` to use, in µL.  This takes precedence over the 
        `pH_volume_x` setting.

    molbio.spin_cleanup.presets.<name>.pH_volume_x
        How much `pH_buffer` to use, as a multiple of the sample volume.  
        This is superseded by the `pH_volume_uL` setting.

    molbio.spin_cleanup.presets.<name>.pH_color
        The color the sample/binding buffer should be after reaching the 
        correct pH.

    molbio.spin_cleanup.presets.<name>.wash_buffer
        The name of the buffer to use when washing the column.  This can either 
        be a string or a list of strings.  Use a list to specify that there 
        should be multiple wash steps.  If this option is a list, the 
        `wash_volume_uL`, `wash_spin_sec`, and `wash_vacuum` options must also 
        be lists of the same length (or left unspecified).

    molbio.spin_cleanup.presets.<name>.wash_volume_uL
        The volume of `wash_buffer` to use, in µL.  This can either be a number 
        or a list of numbers; see `wash_buffer` for more details.

    molbio.spin_cleanup.presets.<name>.wash_spin_sec
        How long to centrifuge the column during the wash step.  This can 
        either be a number or a list of numbers; see `wash_buffer` for more 
        details.

    molbio.spin_cleanup.presets.<name>.wash_vacuum
        Whether or not to use a vacuum manifold for the wash step.  This can 
        either be a boolean or a list of booleans; see `wash_buffer` for more 
        details.  The default is False.  If True, the `wash_spin_sec` option is 
        ignored.

    molbio.spin_cleanup.presets.<name>.dry_spin_sec
        How long to centrifuge the column after the wash step(s), e.g. to 
        remove any residual ethanol.  If left unspecified, this step will not 
        be included in the protocol.

    molbio.spin_cleanup.presets.<name>.elute_buffer
        The default value for the `--elute-buffer` flag.

    molbio.spin_cleanup.presets.<name>.elute_volume_uL
        The default value for the `--elute-volume` flag.

    molbio.spin_cleanup.presets.<name>.elute_min_volume_uL
        The minimum recommended volume to elute in.  Smaller volumes can still 
        be specified, but will be accompanied by a warning.

    molbio.spin_cleanup.presets.<name>.elute_wait_sec
        How long to incubate the column with elution buffer before eluting, in 
        seconds.

    molbio.spin_cleanup.presets.<name>.elute_spin_sec
        How long to centrifuge the column when eluting.

Database:
    Spin-column cleanup protocols can appear in the "Cleanups" column of a 
    FreezerBox database:

        spin-cleanup [<preset>] [volume=<µL>] [buffer=<name>]
    
    <preset>
        See `<preset>`.

    volume=<µL>
        See `--elute-volume`.  Must specify a unit.

    buffer=<µL>
        See `--elute-buffer`.
"""
    __config__ = [
            DocoptConfig(),
            MakerArgsConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.spin_cleanup'),
    ]
    preset_briefs = appcli.config_attr()
    config_paths = appcli.config_attr()
    preset_brief_template = '{protocol_name}'

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '<preset>'),
            Key(MakerArgsConfig, 1),
            Key(StepwiseConfig, 'default_preset'),
    )
    protocol_name = appcli.param(
            Key(PresetConfig, 'protocol_name'),
    )
    protocol_link = appcli.param(
            Key(PresetConfig, 'protocol_link'),
            default=None,
    )
    column_name = appcli.param(
            Key(PresetConfig, 'column_name'),
            default='silica spin column',
    )
    spin_speed_g = appcli.param(
            Key(PresetConfig, 'spin_speed_g'),
            default=None,
    )
    column_capacity_ug = appcli.param(
            Key(PresetConfig, 'column_capacity_ug'),
            default=None,
    )
    sample_type = appcli.param(
            Key(PresetConfig, 'sample_type'),
            default='DNA',
    )
    sample_volume_uL = appcli.param(
            Key(DocoptConfig, '--sample-volume', cast=float),
            default=None,
    )
    target_sample_volume_uL = appcli.param(
            Key(PresetConfig, 'sample_volume_uL'),
            default=None,
    )
    bind_buffer = appcli.param(
            Key(PresetConfig, 'bind_buffer'),
    )
    bind_volume_uL = appcli.param(
            Key(PresetConfig, 'bind_volume_uL'),
            default=None
    )
    bind_volume_x = appcli.param(
            Key(PresetConfig, 'bind_volume_x'),
            default=None
    )
    bind_spin_sec = appcli.param(
            Key(PresetConfig, 'bind_spin_sec'),
            default=None
    )
    bind_vacuum = appcli.param(
            Key(PresetConfig, 'bind_vacuum'),
            default=False,
    )
    ph_buffer = appcli.param(
            Key(PresetConfig, 'pH_buffer'),
            default=None,
    )
    ph_volume_uL = appcli.param(
            Key(PresetConfig, 'pH_volume_uL'),
            default=None
    )
    ph_volume_x = appcli.param(
            Key(PresetConfig, 'pH_volume_x'),
            default=None
    )
    ph_color = appcli.param(
            Key(PresetConfig, 'pH_color'),
    )
    wash_buffer = appcli.param(
            Key(PresetConfig, 'wash_buffer'),
    )
    wash_volume_uL = appcli.param(
            Key(PresetConfig, 'wash_volume_uL'),
    )
    wash_spin_sec = appcli.param(
            Key(PresetConfig, 'wash_spin_sec'),
            default=None,
    )
    wash_vacuum = appcli.param(
            Key(PresetConfig, 'wash_vacuum'),
            default=False,
    )
    dry_spin_sec = appcli.param(
            Key(PresetConfig, 'dry_spin_sec'),
            default=None,
    )
    elute_buffer = appcli.param(
            Key(DocoptConfig, '--elute-buffer'),
            Key(MakerArgsConfig, 'buffer'),
            Key(PresetConfig, 'elute_buffer'),
    )
    elute_volume_uL = appcli.param(
            Key(DocoptConfig, '--elute-volume', cast=float),
            Key(MakerArgsConfig, 'volume', cast=parse_volume_uL),
            Key(PresetConfig, 'elute_volume_uL'),
    )
    elute_min_volume_uL = appcli.param(
            Key(PresetConfig, 'elute_min_volume_uL'),
            default=None,
    )
    elute_wait_sec = appcli.param(
            Key(PresetConfig, 'elute_wait_sec'),
            default=None,
    )
    elute_spin_sec = appcli.param(
            Key(PresetConfig, 'elute_spin_sec'),
    )

    group_by = {
            'preset': group_by_identity,
            'elute_buffer': group_by_identity,
            'elute_volume_uL': group_by_identity,
    }

    def __init__(self, preset=None):
        if preset is not None:
            self.preset = preset

    def get_protocol(self):
        p = stepwise.Protocol()
        pl = stepwise.paragraph_list()
        ul = stepwise.unordered_list()

        def break_if_too_long(pl, ul, n=4):
            if len(ul) > n:
                ul = stepwise.unordered_list()
                pl += ul
            return ul

        footnotes = []
        if self.protocol_link:
            footnotes.append(self.protocol_link)
        if self.column_capacity_ug:
            footnotes.append(f"Column capacity: {self.column_capacity_ug} µg")

        if self.products and self.show_product_names:
            product_names = oxford_comma(self.products) + ' '
        else:
            product_names = ''

        p += pl
        pl += f"Purify {product_names}using {self.protocol_name}{p.add_footnotes(*footnotes)}:"
        pl += ul

        if self.spin_speed_g:
            ul += f"Perform all spin steps at {self.spin_speed_g}g."

        ## Dilute
        if x := self.target_sample_volume_uL:
            v = self.sample_volume_uL

            if not isinstance(x, dict):
                target = f'{x} µL'
                skip = v and v == x
                self.sample_volume_uL = x
            elif 'min' in x and 'max' in x:
                target = f"between {x['min']}–{x['max']} µL"
                skip = v and x['min'] <= v <= x['max']
            elif 'min' in x:
                target = f"at least {x['min']} µL"
                skip = v and x['min'] <= v
            elif 'max' in x:
                target = f"at most {x['max']} µL"
                skip = v and v <= x['max']

            if not skip:
                ul += f"Ensure that the sample is {target}."

        ## Bind
        bind_params = zip_params(
                self.bind_buffer,
                self.bind_volume_x,
                self.bind_volume_uL,
        )
        for bind_buffer, bind_volume_x, bind_volume_uL in bind_params:
            bind_volume = resolve_volume(bind_volume_uL, bind_volume_x, self.sample_volume_uL)
            ul += f"Add {bind_volume} {bind_buffer} to the crude {self.sample_type}."

        if self.ph_buffer:
            ph_volume = resolve_volume(self.ph_volume_uL, self.ph_volume_x, self.sample_volume_uL)
            ul += f"If not {self.ph_color}: Add {ph_volume} {self.ph_buffer}."

        ul += f"Load on a {self.column_name}."
        ul += flush_column(self.bind_spin_sec, self.bind_vacuum)
        ul = break_if_too_long(pl, ul)

        ## Wash
        wash_params = zip_params(
                self.wash_buffer,
                self.wash_volume_uL,
                self.wash_spin_sec,
                self.wash_vacuum,
        )
        for wash_buffer, wash_volume_uL, wash_spin_sec, wash_vacuum in wash_params:
            ul += f"Add {wash_volume_uL} µL {wash_buffer}."
            ul += flush_column(wash_spin_sec, wash_vacuum)

        ## Dry
        if self.dry_spin_sec:
            ul += flush_column(self.dry_spin_sec)

        ul = break_if_too_long(pl, ul)

        ## Elute
        if self.elute_volume_uL < self.elute_min_volume_uL:
            warn(f"Elution volume ({self.elute_volume_uL} µL) is below the recommended minimum ({self.elute_min_volume_uL} µL).")

        ul += f"Add {self.elute_volume_uL} µL {self.elute_buffer}."
        if self.elute_wait_sec:
            ul += f"Wait at least {format_sec(self.elute_wait_sec)}."
        ul += flush_column(self.elute_spin_sec, keep_flowthrough=True)

        return p

def zip_params(*params):
    from itertools import repeat
    from more_itertools import always_iterable

    yield from zip(*(
            always_iterable(p or repeat(p))
            for p in params
    ))

def resolve_volume(volume_uL, volume_x, sample_volume_uL):
    if volume_uL:
        return f'{volume_uL} µL'
    elif sample_volume_uL:
        return f'{volume_x * sample_volume_uL} µL'
    else:
        return f'{volume_x} volumes'

def flush_column(spin_time_sec, use_vacuum=False, keep_flowthrough=False):
    if use_vacuum:
        return "Apply vacuum."
    else:
        if not spin_time_sec:
            raise ValueError("no spin time specified")
        return f"Spin {format_sec(spin_time_sec)}; {'keep' if keep_flowthrough else 'discard'} flow-through."


if __name__ == '__main__':
    SpinCleanup.main()

