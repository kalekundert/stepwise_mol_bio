#!/usr/bin/env python3

import stepwise, appcli, autoprop, freezerbox

from stepwise import (
        StepwiseConfig, PresetConfig, pl, ul, dl,
)
from stepwise_mol_bio import (
        Main, Bindable, BindableReagent, UsageError,
        bind, comma_list,
)
from freezerbox import (
        MakerConfig, ProductConfig, ReagentConfig, Plasmid, Strain,
        parse_bool, join_lists, group_by_identity, unanimous,
)
from appcli import Key, Method, DocoptConfig
from itertools import combinations
from more_itertools import (
        one, flatten, always_iterable, unique_everseen as unique,
)
from inform import plural

def parse_transformations_from_docopt(args):
    return [
            parse_transformation_from_str(p)
            for p in args['<plasmids>']
    ]

def parse_transformations_from_freezerbox(reagent):
    fields = reagent.maker_args

    if isinstance(reagent, Strain):
        plasmid_tags = [x.tag for x in reagent.plasmids]
        strain_tag = reagent.parent_strain

    if isinstance(reagent, Plasmid):
        plasmid_tags = [reagent.tag]

        # The strain doesn't need to be specified.  If not specified, the 
        # strain from the default preset will be used.
        strain_tag = fields.get('strain')

    plasmids = [Transform.Plasmid(x) for x in plasmid_tags]
    strain = Transform.Strain(strain_tag) if strain_tag else None

    return [Transform.Reaction(plasmids, strain)]

def parse_transformation_from_str(rxn_str):
    fields = rxn_str.split('>', 1)
    plasmids = [Transform.Plasmid(x) for x in fields[0].split('+')]
    strain = None

    if len(fields) == 2:
        strain = Transform.Strain(fields[1])

    return Transform.Reaction(plasmids, strain)

@autoprop.cache
class Transform(Main):
    """\
Transform one or more plasmids into chemically competent cells.

Usage:
    transform <plasmids>... [-p <preset>] [-s <strain>] [-a <antibiotics>]

Arguments:
    <plasmids>
        The names of the plasmids to transform.  If these names can be found in 
        the FreezerBox database, a number of default parameters (e.g.  
        resistance, recovery time) will be derived from them.  You can indicate 
        that multiple plasmids should be co-transformed by separating their 
        name with plus (+) signs.  You can also indicate which strain the 
        plasmid(s) should be transformed into by specifying a greater-than sign 
        (>) followed by the strain at the end of the argument.

<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>              [default: ${app.preset}]
        The default parameters to use.  The following sets of parameters are 
        currently available:

        ${hanging_indent(app.preset_briefs, 8*' ')}

    -s --strain <name>
        The name of the strain of competent cells to transform into.  If this 
        name can be found in the FreezerBox database, any plasmids harbored by 
        this strain will be taken into account when choosing the antibiotic.

    -a --antibiotics <name,...>
        The names of the antibiotics to plate with, separated by commas.  
        Usually this is determined automatically based on the plasmids being 
        transformed, but this option provides a means to override the default 
        if necessary.
"""
    __config__ = [
            DocoptConfig,
            MakerConfig,
            ProductConfig,
            PresetConfig,
            StepwiseConfig.setup('molbio.transform'),
    ]

    class Reaction(Bindable, use_app_configs=True):

        def __init__(self, plasmids, strain=None, *, antibiotics=None):
            self.plasmids = plasmids
            if strain: self.strain = strain
            if antibiotics: self.plate_antibiotics = antibiotics

        def __repr__(self):
            plasmid_tags = [x.tag for x in self.plasmids]
            strain_tag = self.strain.tag
            return f'Transform.Reaction({plasmid_tags!r}, {strain_tag!r})'

        @classmethod
        def from_tags(cls, plasmid_tags, strain_tag=None, *, antibiotics=None):
            plasmids = [Transform.Plasmid(x) for x in plasmid_tags]
            strain = Transform.Strain(strain_tag) if strain_tag else None
            return cls(plasmids, strain, antibiotics=antibiotics)

        def on_bind(self, app, force=False):
            super().on_bind(app)
            bind(app, self.plasmids, force=force)
            bind(app, self.strain, force=force)

        def _infer_antibiotics(self):
            return infer_antibiotics(
                    self.plasmids,
                    self.strain,
            )

        strain = appcli.param(
                Key(DocoptConfig, '--strain'),
                Key(MakerConfig, 'strain'),
                Key(PresetConfig, 'strain'),
                cast=lambda x: Transform.Strain(x),
        )
        plate_media = appcli.param(
                Key(PresetConfig, 'plate_media'),
                default='LB',
        )
        plate_antibiotics = appcli.param(
                Key(DocoptConfig, '--antibiotics', cast=comma_list),
                Key(PresetConfig, 'antibiotics'),
                Method(_infer_antibiotics),
        )

    class Plasmid(BindableReagent, use_app_configs=True):
        quantity = appcli.param(
                Key(PresetConfig, 'plasmid_quantity'),
        )
        antibiotic = appcli.param(
                # It's possible for plasmids to have multiple resistance genes, 
                # but in these cases the user will just have to manually 
                # specify which ones to use.
                Key(ReagentConfig, 'antibiotics', cast=one),
        )
        ori = appcli.param(
                Key(ReagentConfig, 'ori'),
                default=None,
        )

    class Strain(BindableReagent):
        antibiotics = appcli.param(
                Key(ReagentConfig, 'resistance'),
                Method(lambda self: [p.antibiotic for p in self.plasmids]),
        )
        plasmids = appcli.param(
                Key(ReagentConfig, 'plasmids', cast=lambda xs: [
                    Transform.Plasmid(x.tag)
                    for x in xs
                ]),
                default_factory=list,
        )

    def __init__(self, transformations):
        self.transformations = transformations

    @classmethod
    def from_tags(cls, plasmid_tags, strain_tag=None, *, antibiotics=None):
        rxn = cls.Reaction.from_tags(
                plasmid_tags,
                strain_tag,
                antibiotics=antibiotics,
        )
        return cls([rxn])

    def get_protocol(self):
        p = stepwise.Protocol()

        for rxn in self.transformations:
            check_oris(rxn.plasmids, self.incompatibility_groups)

        strains = make_strain_map(self.transformations)
        media_antibiotics = make_media_antibiotic_map(self.transformations)

        if len(self.transformations) == 1:
            rxn = one(self.transformations)
            intro = f"Transform {format_plasmids(rxn.plasmids)} into {format_strain(rxn.strain)}{p.add_footnotes(*self.footnotes)}:"
        else:
            plasmids = list(unique(x.plasmids for x in self.transformations))
            intro = f"Transform the following {plural(plasmids):plasmid/s}: {format_plasmid_groups(plasmids)}{p.add_footnotes(*self.footnotes)}"

        ul1 = ul()
        pl1 = pl(ul1)
        p += pl(intro, pl1)

        if len(media_antibiotics) == 1:
            k, v = one(media_antibiotics.items())
            ul1 += f"Pre-warm {plural(v):# {format_media_antibiotics(*k)} plate/s}."
        else:
            ul1 += pl(
                    "Pre-warm selective plates:",
                    dl(*(
                        (format_media_antibiotics(*k), format_plasmid_groups(v))
                        for k, v in media_antibiotics.items()
                    )),
                    br='\n',
            )

        if len(self.transformations) > 1:
            ul2 = ul()
            pl2 = pl("For each transformation:", ul2)
            ul1 += pl2
        else:
            ul2 = ul1
            pl2 = pl1

        try:
            strain = unanimous(strains)
        except ValueError:
            ul2 += pl(
                    f"Thaw {self.cell_volume_uL} µL of the appropriate competent cells on ice:",
                    dl(*(
                        (format_strain(k), format_plasmid_groups(v))
                        for k, v in strains.items()
                    )),
                    br='\n',
            )
            ul2 = ul(); pl2 += ul2
        else:
            ul2 += f"Thaw {self.cell_volume_uL} µL competent {format_strain(strain)} cells on ice."
        
        try:
            quantity = unanimous(
                    p.quantity.strip()
                    for t in self.transformations
                    for p in t.plasmids
            )
        except ValueError:
            ul2 += pl(
                    "Add the following quantities of plasmids:",
                    dl(*(
                        (
                            format_transformation(x, strains),
                            format_quantities(p.quantity for p in x.plasmids),
                        )
                        for x in self.transformations
                    )),
                    br='\n',
            )
            ul2 = ul(); pl2 += ul2
        else:
            ul2 += f"Add {quantity} plasmid."

        ul2 += "Gently flick to mix"
        if not self.skip_rest:
            ul2 += f"Incubate on ice for {self.rest_time_min} min."

        ul2 = ul(); pl2 += ul2

        ul2 += f"Incubate at {self.heat_shock_temp_C}°C for {self.heat_shock_time_sec}s."
        ul2 += f"Incubate on ice for 2 min."
        ul2 += f"Add {self.recovery_media_volume_uL} µL {self.recovery_media}."
        if not self.skip_recovery:
            ul2 += f"Incubate at {self.recovery_temp_C}°C for {self.recovery_time_min} min with end-over-end mixing."

            ul2 = ul(); pl2 += ul2

            if self.concentrate:
                ul2 += f"Spin at {self.conc_spin_speed_g}g for {self.conc_spin_time_min} min."
                ul2 += f"Remove {self.cell_volume_uL + self.recovery_media_volume_uL - self.conc_volume_uL} µL media."
                ul2 += f"Resuspend pelleted cells."

        ul2 = ul(); pl2 += ul2

        if self.plate_dilution != 1:
            ul2 += f"Dilute cells 1:{self.plate_dilution}."

        ul2 += f"Plate {self.plate_volume_uL} µL cells."
        ul2 += f"Incubate at {self.outgrowth_temp_C}°C for {self.outgrowth_time_h}h."

        return p

    def get_dependencies(self):
        deps = set()

        for rxn in self.transformations:
            deps.add(rxn.strain.tag)
            deps.update(p.tag for p in rxn.plasmids)

        deps.difference_update(p.tag for p in self.products)
        return deps

    def _skip_rest_recovery(self):
        skip = set(self.skip_recovery_antibiotics)
        return not any(
                set(x.plate_antibiotics) - skip
                for x in self.transformations
        )

    preset_briefs = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=list,
    )
    preset = appcli.param(
            Key(DocoptConfig, '--preset'),
            Key(MakerConfig, 'preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    transformations = appcli.param(
            Key(DocoptConfig, parse_transformations_from_docopt),
            Key(ProductConfig, parse_transformations_from_freezerbox),
            get=bind,
    )
    incompatibility_groups = appcli.param(
            Key(StepwiseConfig, 'incompatibility_groups'),
            pick=join_lists,
    )
    skip_rest = appcli.param(
            Key(DocoptConfig, '--fast'),
            Key(DocoptConfig, '--skip-rest'),
            Key(PresetConfig, 'skip_rest'),
            Method(_skip_rest_recovery),
    )
    skip_recovery = appcli.param(
            Key(DocoptConfig, '--fast'),
            Key(DocoptConfig, '--skip-recovery'),
            Key(PresetConfig, 'skip_recovery'),
            Method(_skip_rest_recovery),
    )
    skip_recovery_antibiotics = appcli.param(
            Key(StepwiseConfig, 'skip_recovery_antibiotics'),
    )
    cell_volume_uL = appcli.param(
            Key(PresetConfig, 'cell_volume_uL'),
    )
    rest_time_min = appcli.param(
            Key(PresetConfig, 'rest_time_min'),
    )
    heat_shock_temp_C = appcli.param(
            Key(PresetConfig, 'heat_shock_temp_C'),
    )
    heat_shock_time_sec = appcli.param(
            Key(PresetConfig, 'heat_shock_time_sec'),
    )
    recovery_media = appcli.param(
            Key(PresetConfig, 'recovery_media'),
    )
    recovery_media_volume_uL = appcli.param(
            Key(PresetConfig, 'recovery_media_volume_uL'),
    )
    recovery_temp_C = appcli.param(
            Key(PresetConfig, 'recovery_temp_C'),
    )
    recovery_time_min = appcli.param(
            Key(PresetConfig, 'recovery_time_min'),
    )
    concentrate = appcli.param(
            Key(PresetConfig, 'concentrate', cast=parse_bool),
            default=False,
    )
    conc_spin_speed_g = appcli.param(
            Key(PresetConfig, 'conc_spin_speed_g'),
    )
    conc_spin_time_min = appcli.param(
            Key(PresetConfig, 'conc_spin_time_min'),
    )
    conc_volume_uL = appcli.param(
            Key(PresetConfig, 'conc_volume_uL'),
    )
    plate_dilution = appcli.param(
            Key(PresetConfig, 'plate_dilution'),
            default=1,
    )
    plate_volume_uL = appcli.param(
            Key(PresetConfig, 'plate_volume_uL'),
    )
    outgrowth_temp_C = appcli.param(
            Key(PresetConfig, 'outgrowth_temp_C'),
            Key(PresetConfig, 'recovery_temp_C'),
    )
    outgrowth_time_h = appcli.param(
            Key(PresetConfig, 'outgrowth_time_h'),
    )
    footnotes = appcli.param(
            Key(PresetConfig, 'footnotes'),
    )

    group_by = {
            'preset': group_by_identity,
    }
    merge_by = {
            'transformations': join_lists,
    }

def check_oris(plasmids, incompatibility_groups):
    """
    Make sure all of the given plasmids have compatible ORIs.
    """
    plasmids = [p for p in plasmids if p.ori]

    for group in incompatibility_groups:
        for a, b in combinations(plasmids, 2):
            if a.ori in group and b.ori in group:
                err = UsageError(plasmids=[a,b])
                err.brief += "can't co-transform two plasmids with incompatible ORIs"
                err.info += lambda e: '\n'.join(
                        f'{p.tag}: {p.ori}'
                        for p in e.plasmids
                )
                raise err

def infer_antibiotics(plasmids, strain):
    antibiotics = {}

    def iter_construct_antibiotics():
        for plasmid in plasmids:
            yield plasmid, plasmid.antibiotic

        for antibiotic in strain.antibiotics:
            yield strain, antibiotic

    for construct, antibiotic in iter_construct_antibiotics():
        antibiotics.setdefault(antibiotic, []).append(construct)

    for antibiotic, constructs in antibiotics.items():
        if len(constructs) > 1:
            err = UsageError(
                    antibiotic=antibiotic,
                    constructs=constructs,
            )
            err.brief = "can't use the same antibiotic to select for multiple constructs"
            err.blame += lambda e: f"the following constructs provide {e.antibiotic} resistance:\n{', '.join(x.tag for x in e.constructs)}"
            raise err

    return tuple(antibiotics.keys())

def make_media_antibiotic_map(transformations):
    out = {}

    for tfn in transformations:
        key = tfn.plate_media, tuple(tfn.plate_antibiotics)
        out.setdefault(key, []).append(tfn.plasmids)

    return out

def make_strain_map(transformations):
    strains = {}

    for tfn in transformations:
        strains.setdefault(tfn.strain, []).append(tfn.plasmids)

    return strains

def format_transformation(transformation, strains):
    plasmids_str = format_plasmids(transformation.plasmids)
    strain_str = format_strain(transformation.strain)

    if len(strains) == 1:
        return plasmids_str
    else:
        return f'{plasmids_str}→{strain_str}'

def format_strain(strain):
    return strain.tag

def format_plasmid_groups(plasmids):
    return ', '.join(format_plasmids(p) for p in plasmids)

def format_plasmids(plasmid):
    return '+'.join(p.tag for p in plasmid)

def format_quantities(quantities):
    return ', '.join(quantities)

def format_media_antibiotics(media, antibiotics):
    return '+'.join([media, *antibiotics])


if __name__ == '__main__':
    Transform.main()
