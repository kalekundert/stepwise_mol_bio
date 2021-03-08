#!/usr/bin/env python3

import stepwise, appcli, autoprop
from appcli import Key, DocoptConfig
from stepwise import StepwiseConfig, PresetConfig, pl, ul
from stepwise_mol_bio import Main, merge_dicts, format_min, format_sec
from inform import plural
from collections import defaultdict

@autoprop
class Stain(Main):
    """\
Stain a gel.

Usage:
    stain <preset> [-t <min>] [-i <cmd> | -I]

<%! from stepwise_mol_bio import hanging_indent %>\
Arguments:
    <preset>
        The default parameters to use.  The following presets are available:

        ${hanging_indent(app.preset_briefs, 8*' ')}

Options:
    -t --time <min>
        How long to incubate the gel in the stain, in minutes.

    -i --imaging-protocol <cmd>

    -I --no-imaging
        Don't include the imaging step in the protocol (e.g. so you can provide
        a customized alternative).
"""
    __config__ = [
            DocoptConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.stain'),
    ]
    preset_briefs = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '<preset>'),
    )
    title = appcli.param(
            Key(PresetConfig),
    )
    light_sensitive = appcli.param(
            Key(PresetConfig),
            default=False,
    )
    prep_buffer = appcli.param(
            Key(PresetConfig),
    )
    prep_volume_mL = appcli.param(
            Key(PresetConfig),
            default=30,
    )
    prep_time_min = appcli.param(
            Key(PresetConfig),
            Key(PresetConfig, 'prep_time_hr', cast=lambda x: x*60),
    )
    prep_microwave = appcli.param(
            Key(PresetConfig),
            default=False,
    )
    prep_steps = appcli.param(
            Key(PresetConfig),
    )
    prep_repeats = appcli.param(
            Key(PresetConfig),
            default=1,
    )
    stain_buffer = appcli.param(
            Key(PresetConfig),
    )
    stain_volume_mL = appcli.param(
            Key(PresetConfig),
            Key(StepwiseConfig, 'default_stain_volume_mL'),
            default=30,
    )
    stain_time_min = appcli.param(
            Key(DocoptConfig, '--time'),
            Key(PresetConfig),
            Key(PresetConfig, 'stain_time_hr', cast=lambda x: x*60),
    )
    stain_microwave = appcli.param(
            Key(PresetConfig),
            default=False,
    )
    stain_steps = appcli.param(
            Key(PresetConfig),
    )
    stain_repeats = appcli.param(
            Key(PresetConfig),
            default=1,
    )
    stain_rinse_repeats = appcli.param(
            Key(PresetConfig),
            default=0,
    )
    destain_buffer = appcli.param(
            Key(PresetConfig),
    )
    destain_volume_mL = appcli.param(
            Key(PresetConfig),
            default=30,
    )
    destain_time_min = appcli.param(
            Key(PresetConfig),
            Key(PresetConfig, 'destain_time_hr', cast=lambda x: x*60),
    )
    destain_microwave = appcli.param(
            Key(PresetConfig),
            default=False,
    )
    destain_steps = appcli.param(
            Key(PresetConfig),
    )
    destain_repeats = appcli.param(
            Key(PresetConfig),
            default=1,
    )
    imaging_cmd = appcli.param(
            Key(DocoptConfig, '--imaging-protocol'),
            Key(DocoptConfig, '--no-imaging', cast=None),
            Key(PresetConfig, 'imaging_protocol'),
            default=None,
    )
    protocol_link = appcli.param(
            Key(PresetConfig),
            default=None,
    )
    footnotes = appcli.param(
            Key(PresetConfig),
            default_factory=list,
    )

    def get_protocol(self):
        p = stepwise.Protocol()

        footnotes = [x] if (x := self.protocol_link) else []
        footnotes += self.footnotes

        s = pl(f"Stain gel with {self.title}{p.add_footnotes(*footnotes)}:")

        if self.light_sensitive:
            #s += "Keep gel in the dark in all following steps."
            s += ul("Keep the stain protected from light.")
            #s += ul("In all following steps, protect the stain from light.")

        s += self._format_stage('prep')
        s += self._format_stage('stain')
        s += self._format_stage('destain')

        p += s

        if self.imaging_cmd:
            p += stepwise.load(self.imaging_cmd)

        return p

    def _format_stage(self, prefix):

        class SkipStage(Exception):
            pass

        def has(attr):
            return hasattr(self, f'{prefix}_{attr}')

        def get(attr, *args):
            return getattr(self, f'{prefix}_{attr}', *args)


        def have_repeats():
            n = get('repeats')
            return n > 1 if isinstance(n, int) else bool(n)

        def format_repeats():
            n = get('repeats')
            if isinstance(n, int):
                n = f'{n}x'
            return f"Repeat {n}:"

        def format_submerge():
            if not has('buffer'): raise SkipStage
            return f"Submerge gel in ≈{get('volume_mL', 30)} mL {get('buffer')}."

        def format_microwave():
            if get('microwave', False):
                return f"Microwave until almost boiling (≈{format_sec(get('microwave_time_s', 45))})."

        def format_incubate():
            return f"Shake gently for {format_min(get('time_min'))}."

        step_templates = get('steps', ['submerge', 'microwave', 'incubate'])
        step_formatters = {
                'submerge': format_submerge,
                'microwave': format_microwave,
                'incubate': format_incubate,
        }

        out = steps = ul()
        for template in step_templates:
            formatter = step_formatters.get(template, lambda: template)
            try:
                steps += formatter()
            except SkipStage:
                return

        if have_repeats():
            out = ul(pl(format_repeats(), steps, br='\n'))

        if n := get('rinse_repeats', False):
            out += f"Rinse {plural(n)://#x }with water."

        return out



if __name__ == '__main__':
    Stain.main()


