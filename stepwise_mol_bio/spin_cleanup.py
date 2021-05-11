#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import warn
from appcli import Key, DocoptConfig
from stepwise import StepwiseConfig, PresetConfig, Quantity, oxford_comma
from stepwise_mol_bio import Cleanup, format_sec, merge_dicts
from freezerbox import MakerArgsConfig, group_by_identity

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

    molbio.spin_cleanup.presets.<name>.title
        How to refer to the whole protocol.  Commonly this is the name of the 
        spin column kit.

    molbio.spin_cleanup.presets.<name>.column_name
        How to refer to the spin column in the protocol.

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
        The name of the buffer to use to bind the sample to column.
        
    molbio.spin_cleanup.presets.<name>.bind_volume_uL
        How much `bind_buffer` to use, in µL.  This takes precedence over the 
        `bind_volume_x` setting.

    molbio.spin_cleanup.presets.<name>.bind_volume_x
        How much `bind_buffer` to use, as a multiple of the sample volume.  
        This is superseded by the `bind_volume_uL` setting.

    molbio.spin_cleanup.presets.<name>.bind_spin_sec
        How long to centrifuge the column during the bind step.

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
        The name of the buffer to use when washing the column.

    molbio.spin_cleanup.presets.<name>.wash_volume_uL
        The volume of `wash_buffer` to use, in µL.

    molbio.spin_cleanup.presets.<name>.wash_spin_sec
        How long to centrifuge the column during the wash step.

    molbio.spin_cleanup.presets.<name>.dry_buffer
        The name of the buffer to use when drying the column.  This can be left 
        empty to indicate that no buffer should be used for this step.

    molbio.spin_cleanup.presets.<name>.dry_volume_uL
        The volume of `dry_buffer` to use, in µL.

    molbio.spin_cleanup.presets.<name>.dry_spin_sec
        How long to centrifuge the column during the drying step.

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

        spin-cleanup [preset=<preset>] [volume=<µL>] [buffer=<name>]
    
    preset=<preset>
        See `--preset`.

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

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '<preset>'),
            Key(MakerArgsConfig, 'preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    protocol_link = appcli.param(
            Key(PresetConfig, 'protocol_link'),
            default=None,
    )
    title = appcli.param(
            Key(PresetConfig, 'title'),
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
    )
    dry_buffer = appcli.param(
            Key(PresetConfig, 'dry_buffer'),
            default=None,
    )
    dry_volume_uL = appcli.param(
            Key(PresetConfig, 'dry_volume_uL'),
    )
    dry_spin_sec = appcli.param(
            Key(PresetConfig, 'dry_spin_sec'),
    )
    elute_buffer = appcli.param(
            Key(DocoptConfig, '--elute-buffer'),
            Key(MakerArgsConfig, 'buffer'),
            Key(PresetConfig, 'elute_buffer'),
    )
    elute_volume_uL = appcli.param(
            Key(DocoptConfig, '--elute-volume', cast=int),
            Key(MakerArgsConfig, 'volume'),
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
            'elute_volume_uL': group_by_identity,
            'elute_buffer_uL': group_by_identity,
    }

    def __init__(self, preset=None):
        if preset is not None:
            self.preset = preset

    def get_protocol(self):
        p = stepwise.Protocol()

        footnotes = []
        if self.protocol_link:
            footnotes.append(self.protocol_link)
        if self.column_capacity_ug:
            footnotes.append(f"Column capacity: {self.column_capacity_ug} µg")

        pl = stepwise.paragraph_list()
        ul = stepwise.unordered_list()

        if self.products and self.show_product_names:
            product_names = oxford_comma(self.products) + ' '
        else:
            product_names = ''

        pl += f"Purify {product_names}using {self.title}{p.add_footnotes(*footnotes)}:"
        pl += ul

        p += pl

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
        bind_volume = resolve_volume(self.bind_volume_uL, self.bind_volume_x, self.sample_volume_uL)
        ul += f"Add {bind_volume} {self.bind_buffer} to the crude {self.sample_type}."

        if self.ph_buffer:
            ph_volume = resolve_volume(self.ph_volume_uL, self.ph_volume_x, self.sample_volume_uL)
            ul += f"If not {self.ph_color}: Add {ph_volume} {self.ph_buffer}."

        ul += f"Load on a {self.column_name}."
        ul += f"Spin {format_sec(self.bind_spin_sec)}; discard flow-through."

        ## Wash
        ul += f"Add {self.wash_volume_uL} µL {self.wash_buffer}."
        ul += f"Spin {format_sec(self.wash_spin_sec)}; discard flow-through."

        ## Dry
        if self.dry_buffer:
            ul += f"Add {self.dry_volume_uL} µL {self.dry_buffer}."
        ul += f"Spin {format_sec(self.dry_spin_sec)}; discard flow-through."

        ## Elute
        if self.elute_volume_uL < self.elute_min_volume_uL:
            warn(f"Elution volume ({self.elute_volume_uL} µL) is below the recommended minimum ({self.elute_min_volume_uL} µL).")

        ul += f"Add {self.elute_volume_uL} µL {self.elute_buffer}."
        if self.elute_wait_sec:
            ul += f"Wait at least {format_sec(self.elute_wait_sec)}."
        ul += f"Spin {format_sec(self.elute_spin_sec)}; keep flow-through."

        return p

def resolve_volume(volume_uL, volume_x, sample_volume_uL):
    if volume_uL:
        return f'{volume_uL} µL'
    elif sample_volume_uL:
        return f'{volume_x * sample_volume_uL} µL'
    else:
        return f'{volume_x} volumes'

if __name__ == '__main__':
    SpinCleanup.main()

