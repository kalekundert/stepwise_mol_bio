#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import warn
from appcli import Key, DocoptConfig
from stepwise import StepwiseConfig, PresetConfig
from stepwise_mol_bio import Main, format_sec, merge_dicts

@autoprop
class SpinCleanup(Main):
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
"""
    __config__ = [
            DocoptConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.spin_cleanup'),
    ]
    usage = appcli.config_attr()
    preset_briefs = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '<preset>'),
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
            Key(PresetConfig, 'column_capacity'),
            default=None,
    )
    sample_type = appcli.param(
            Key(PresetConfig, 'sample_type'),
            default='DNA',
    )
    sample_volume_uL = appcli.param(
            Key(DocoptConfig, '--sample-volume', cast=int),
            default=None,
    )
    bind_buffer = appcli.param(
            Key(PresetConfig, 'bind_buffer'),
    )
    bind_volume_x = appcli.param(
            Key(PresetConfig, 'bind_volume_x'),
    )
    bind_spin_sec = appcli.param(
            Key(PresetConfig, 'bind_spin_sec'),
    )
    ph_buffer = appcli.param(
            Key(PresetConfig, 'pH_buffer'),
            default=None,
    )
    ph_volume_x = appcli.param(
            Key(PresetConfig, 'pH_volume_x'),
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
            Key(PresetConfig, 'elute_buffer'),
    )
    elute_volume_uL = appcli.param(
            Key(DocoptConfig, '--elute-volume', cast=int),
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

    def get_protocol(self):
        p = stepwise.Protocol()

        footnotes = []
        if self.protocol_link:
            footnotes.append(self.protocol_link)
        if self.column_capacity_ug:
            footnotes.append(f"Column capacity: {self.column_capacity_ug} µg")

        pl = stepwise.paragraph_list()
        ul = stepwise.unordered_list()

        pl += f"Purify using {self.title}{p.add_footnotes(*footnotes)}:"
        pl += ul

        p += pl

        if self.spin_speed_g:
            ul += f"Perform all spin steps at {self.spin_speed_g}g."

        ## Bind
        bind_volume = resolve_volume(self.bind_volume_x, self.sample_volume_uL)
        ul += f"Add {bind_volume} {self.bind_buffer} to the crude {self.sample_type}."

        if self.ph_buffer:
            ph_volume = resolve_volume(self.ph_volume_x, self.sample_volume_uL)
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

def resolve_volume(volume_x, sample_volume_uL):
    if sample_volume_uL:
        return f'{volume_x * sample_volume_uL} µL'
    else:
        return f'{volume_x} volumes'


if __name__ == '__main__':
    SpinCleanup.main()

