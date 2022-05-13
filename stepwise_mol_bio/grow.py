#!/usr/bin/env python3

import stepwise
import byoc
import autoprop

from stepwise import StepwiseConfig, PresetConfig, pl, ul, dl
from stepwise_mol_bio import (
        Main, BindableReagent, UsageError, bind, comma_list
)
from freezerbox import ReagentConfig
from byoc import Key, Attr, DocoptConfig
from itertools import groupby

@autoprop.cache
class Grow(Main):
    """
Grow cultures of microbial cells.

Usage:
    grow <strains>... [-p <preset>] [-m <media>] [-a <antibiotic>] [-v <mL>] 
        [-i <µL> | -I <ratio>] [-T <°C>] [-r <rpm>] [-l <phase>] [-d <od>] 
        [-t <hr>] [-fDO] [options]

Arguments:
    <strains>
        The names of the strains to grow.  If these names can be found in the 
        FreezerBox database, the database will automatically be queried for 
        information like which antibiotics to use.

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>
        Which set of default parameters to use.  The following presets are 
        currently available:

        ${hanging_indent(app.preset_briefs, 8)}

    -m --media <name>
        What media to grow the cells in.

    -a --antibiotics <names>
        What antibiotic to add to the media.  You can specify multiple 
        antibiotics separated by commas.  Note that this option overrides any 
        information about antibiotics that may be in the FreezerBox database.

    -f --fresh
        Specify that freshly-plated colonies (e.g. either transformed or 
        restreaked) should be used.

    -v --volume <mL>
        What volume of media to use, in mL.

    -i --inoculate-volume <µL>
        What volume of overnight culture to inoculate the day culture with, in 
        µL.

    -I --inoculate-ratio <ratio>
        What ratio of sterile media to overnight culture to use when 
        inoculating the day culture.  For example, given a 3 mL culture, a 
        ratio of 100 would mean to inoculate with 30 µL of overnight culture.
        Note that this option is overridden by the `--inoculate-volume` option.

    -T --temp <°C>
        What temperature to grow the cells at, in °C.

    -r --shaking <rpm>
        What speed to shake the cells while they're growing, in rpm.

    -l --phase <desc>
        What phase to grow the cells to, e.g. "log", "early log", "late log", 
        etc.  This can be specified in conjunction with `--density` and 
        `--time`.

    -d --density <od600>
        What density to grow the cells to, in units of OD600.  This can be 
        specified in conjunction with `--phase` and `--time`.

    -t --time <hr>
        How long to grow the cells, in hours.  If you specified either 
        `--phase` or `--density`, this time will be treated as an estimate.

    --overnight-media <name>
        What media to grow the overnight culture in.  By default, this will be 
        the same as the day culture.

    --overnight-volume <mL>
        What volume of media to use for the overnight culture.

    --overnight-temp <°C>
        What temperature to grow the overnight culture at, in °C.  By default, 
        this will be the same as the day culture.

    --overnight-shaking <°C>
        What speed to shake the overnight culture at, in rpm.  By default, this 
        will be the same as the day culture.

    -D --skip-overnight
        Skip the overnight culture step.

    -O --only-overnight
        Only include to overnight culture step.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.grow.default_preset
        The default value for the `--preset` option.

    molbio.grow.presets:
        Named groups of default protocol parameters.  Typically each preset 
        corresponds to a particular experiment or type of cell.  See below for 
        the various settings that can be specified in each preset.

    molbio.grow.presets.<name>.media
        The default value for the `--media` option.

    molbio.grow.presets.<name>.antibiotic
        The default value for the `--antibiotic` option.  Note that this 
        setting, like `--antibiotic`, also overrides any information about 
        antibiotics that may be in the FreezerBox database.  Use with caution.

    molbio.grow.presets.<name>.fresh
        The default value for the `--fresh` option.

    molbio.grow.presets.<name>.volume_mL
        The default value for the `--volume` option.

    molbio.grow.presets.<name>.inoculate_volume_mL
        The default value for the `--inoculate-volume` option.

    molbio.grow.presets.<name>.inoculate_ratio
        The default value for the `--inoculate-ratio` option.

    molbio.grow.presets.<name>.temp_C
        The default value for the `--temp` option.

    molbio.grow.presets.<name>.shaking_rpm
        The default value for the `--shaking` option.

    molbio.grow.presets.<name>.target_phase
        The default value for the `--phase` option.

    molbio.grow.presets.<name>.target_density
        The default value for the `--density` option.

    molbio.grow.presets.<name>.time_hr
        The default value for the `--time` option.

    molbio.grow.presets.<name>.overnight.media
        The default value for the `--overnight-media` option.

    molbio.grow.presets.<name>.overnight.volume_mL
        The default value for the `--overnight-volume` option.

    molbio.grow.presets.<name>.overnight.temp_C
        The default value for the `--overnight-temp` option.

    molbio.grow.presets.<name>.overnight.shaking_rpm
        The default value for the `--overnight-shaking` option.
"""

    class Strain(BindableReagent, use_app_configs=True):
        antibiotics = byoc.param(
                Key(DocoptConfig, '--antibiotics', cast=comma_list),
                Key(PresetConfig, 'antibiotics'),
                Key(ReagentConfig),
                cast=lambda xs: [
                    {'Amp': 'Carb'}.get(x, x)
                    for x in xs
                ],
        )

        @classmethod
        def from_tags(cls, tags):
            return [cls(x) for x in tags]

    class Ratio:

        def __init__(self, ratio):
            self.ratio = ratio

    def __init__(self, strains):
        self.strains = strains

    def get_protocol(self):
        p = stepwise.Protocol()
        if not self.skip_overnight:
            p += self.overnight_culture_step
        if not self.only_overnight:
            p += self.day_culture_step
        return p

    def get_overnight_culture_step(self):
        return pl(f"Grow overnight cultures:",
                self.get_strains_by_antibiotics_dl(self.overnight_media),
                ul(
                    f"Inoculate {self.overnight_volume_mL:g} mL sterile media.",
                    "Use fresh colonies (i.e. either transformed or restreaked within 1 week)." if self.fresh else None,
                    f"Incubate overnight at {self.overnight_temp_C}°C with shaking at {self.overnight_shaking_rpm} rpm."
                ),
        )

    def get_day_culture_step(self):
        step = pl("Grow day cultures:")

        if self.skip_overnight or self.media != self.overnight_media:
            step += self.strains_by_antibiotics_dl

        step += ul(
                f"Inoculate {self.volume_mL:g} mL sterile media with {self.inoculate_volume_uL:g} µL saturated overnight culture.",
                f"Incubate at {self.temp_C:g}°C with shaking at {self.shaking_rpm} rpm {self.target}.",
        )

        return step

    def get_strains_by_antibiotics(self):
        by_antibiotics = lambda x: x.antibiotics
        return [
                (k, list(g))
                for k, g in groupby(
                    sorted(self.strains, key=by_antibiotics),
                    key=by_antibiotics,
                )
        ]

    def get_strains_by_antibiotics_dl(self, media=None):
        media = media or self.media
        return dl(*(
            (format_media(media, antibiotics), format_strains(strains))
            for antibiotics, strains in self.strains_by_antibiotics
        ))

    def get_antibiotics(self):
        return [x for x, _ in self.strains_by_antibiotics]

    def get_media_with_antibiotics(self):
        return '/'.join(
                format_media(self.media, antibiotics)
                for antibiotics in self.antibiotics
        )
    def get_target(self):
        if self.target_phase:
            target = f'until the cells reach {self.target_phase} phase'

            if self.target_od:
                target += f', OD={self.target_od}'

            if self.time_h:
                target += f', ≈{self.time_h:g}h'

            return target

        if self.target_od:
            target = f'until the cells reach OD={self.target_od}'

            if self.time_h:
                target += f', ≈{self.time_h:g}h'

            return target

        if self.time_h:
            return f'for {self.time_h:g}h'

        raise UsageError("must specify how long to grow the cells for")

    __config__ = [
            DocoptConfig,
            PresetConfig,
            StepwiseConfig.setup(('molbio', 'grow')),
    ]
    presets = byoc.param(
            Key(StepwiseConfig, 'presets'),
            pick=list,
    )
    preset = byoc.param(
            Key(DocoptConfig, '--preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    strains = byoc.param(
            Key(DocoptConfig, '<strains>', cast=Strain.from_tags),
            get=bind,
    )
    fresh = byoc.param(
            Key(DocoptConfig, '--fresh'),
            Key(PresetConfig, 'fresh'),
            default=False,
    )
    media = byoc.param(
            Key(DocoptConfig, '--media'),
            Key(PresetConfig, 'media'),
    )
    volume_mL = byoc.param(
            Key(DocoptConfig, '--volume', cast=float),
            Key(PresetConfig, 'volume_mL'),
    )
    inoculate_volume_uL = byoc.param(
            Key(DocoptConfig, '--inoculate-volume', cast=float),
            Key(DocoptConfig, '--inoculate-ratio', cast=[Ratio, float]),
            Key(PresetConfig, 'inoculate_volume_uL'),
            Key(PresetConfig, 'inoculate_ratio', cast=Ratio),
            get=lambda self, x: 1000 * self.volume_mL / x.ratio if isinstance(x, Grow.Ratio) else x
    )
    temp_C = byoc.param(
            Key(DocoptConfig, '--temp', cast=float),
            Key(PresetConfig, 'temp_C'),
    )
    shaking_rpm = byoc.param(
            Key(DocoptConfig, '--shaking', cast=float),
            Key(PresetConfig, 'shaking_rpm'),
    )
    target_phase = byoc.param(
            Key(DocoptConfig, '--phase'),
            Key(PresetConfig, 'target_phase'),
            default=None,
    )
    target_od = byoc.param(
            Key(DocoptConfig, '--od', cast=float),
            Key(PresetConfig, 'target_od'),
            default=None,
    )
    time_h = byoc.param(
            Key(DocoptConfig, '--time', cast=float),
            Key(PresetConfig, 'target_time_h'),
            default=None,
    )
    overnight_media = byoc.param(
            Key(DocoptConfig, '--overnight-media'),
            Key(PresetConfig, ('overnight', 'media')),
            Attr('media'),
    )
    overnight_volume_mL = byoc.param(
            Key(DocoptConfig, '--overnight-volume', cast=float),
            Key(PresetConfig, ('overnight', 'volume_mL')),
            default=1,
    )
    overnight_temp_C = byoc.param(
            Key(DocoptConfig, '--overnight-temp', cast=float),
            Key(PresetConfig, ('overnight', 'temp_C')),
            Attr('temp_C'),
    )
    overnight_shaking_rpm = byoc.param(
            Key(DocoptConfig, '--overnight-shaking', cast=float),
            Key(PresetConfig, ('overnight', 'shaking_rpm')),
            Attr('shaking_rpm'),
    )
    skip_overnight = byoc.param(
            Key(DocoptConfig, '--skip-overnight'),
            default=False,
    )
    only_overnight = byoc.param(
            Key(DocoptConfig, '--only-overnight'),
            default=False,
    )

    preset_briefs = byoc.config_attr()
    config_paths = byoc.config_attr()

def format_media(media, antibiotics):
    sep = ' + ' if ' ' in media else '+'
    return sep.join([media, *antibiotics])

def format_strains(strains):
    return ','.join(x.tag for x in strains)


if __name__ == '__main__':
    Grow.entry_point()

