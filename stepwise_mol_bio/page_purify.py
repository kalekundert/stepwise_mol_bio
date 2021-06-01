#!/usr/bin/env python3

import stepwise
import autoprop
import appcli

from stepwise import (
        PresetConfig, StepwiseConfig, Quantity, pl, ul, oxford_comma,
)
from appcli import Key, Method, DocoptConfig
from freezerbox import (
        ReagentConfig, MakerArgsConfig, ParseError, unanimous, 
        convert_conc_unit, group_by_identity, join_lists,
        parse_volume_uL as parse_strict_volume_uL,
)
from stepwise_mol_bio import (
        Cleanup, Gel, SpinCleanup, ConfigError, merge_dicts, comma_list,
)
from operator import not_
from more_itertools import one
from inform import plural
from math import ceil, inf

def parse_conc(x):
    return Quantity.from_string_or_float(x, 'µg/µL')

def parse_volume_uL(x):
    return float(x)

def parse_samples(values):
    return [Sample.from_colon_separated_string(x) for x in values]

class Sample:
    __config__ = [
            ReagentConfig(
                db_getter=lambda self: self.app.db,
                tag_getter=lambda self: self.name,
                transform=one,
            )
    ]

    def _calc_mass_ug(self):
        try:
            stock_conc_ug_uL = convert_conc_unit(self.stock_conc, self.mw, 'µg/µL')
        except ParseError:
            err = ConfigError(sample=self)
            err.brief = "can't calculate mass in µg for sample: {sample.name!r}"
            err.info += "MW: {sample.mw}"
            err.info += "stock conc: {sample.stock_conc}"
            err.info += "volume: {sample.stock_conc}"
            raise err
        else:
            return Quantity(stock_conc_ug_uL.value * self.volume_uL, 'µg')

    name = appcli.param()
    molecule = appcli.param(
            Method(lambda self: self.app.default_molecule),
            Key(ReagentConfig, 'molecule'),
            default='RNA',
            ignore=None,
    )
    stock_conc = appcli.param(
            Method(lambda self: self.app.default_stock_conc),
            Key(ReagentConfig, 'conc'),
            ignore=None,
    )
    volume_uL = appcli.param(
            Method(lambda self: self.app.default_volume_uL),
            ignore=None,
    )
    mw = appcli.param(
            Key(ReagentConfig, 'mw'),
            default=None,
            ignore=None,
    )
    mass_ug = appcli.param(
            Method(_calc_mass_ug),
            ignore=None,
    )

    @classmethod
    def from_colon_separated_string(cls, value):
        fields = {i: x for i, x in enumerate(value.split(':'))}

        name = fields[0]
        molecule = fields.get(1) or None
        stock_conc = fields.get(2) or None
        volume_uL = fields.get(3) or None

        if stock_conc:
            stock_conc = parse_conc(stock_conc)
        if volume_uL:
            volume_uL = parse_volume_uL(volume_uL)

        return cls(
                name=name,
                molecule=molecule,
                stock_conc=stock_conc,
                volume_uL=volume_uL,
        )

    def __init__(self, name, molecule=None, stock_conc=None, volume_uL=None, mass_ug=None):
        self.app = None
        self.name = name
        self.molecule = molecule
        self.stock_conc = stock_conc
        self.volume_uL = volume_uL
        self.mass_ug = mass_ug

    def __repr__(self):
        attrs = ['name', 'molecule', 'stock_conc', 'volume_uL', 'mass_ug']
        attr_reprs = [
                f'{k}={v!r}'
                for k in attrs
                if (v := getattr(self, k, None))
        ]
        return f'Sample({", ".join(attr_reprs)})'

    def bind(self, app, force=False):
        if not self.app or force:
            self.app = app

@autoprop
class PagePurify(Cleanup):
    """\
Purify nucleic acids by PAGE.

Usage:
    page_purify <samples>... [-p <preset>] [-P <gel>] [-c <conc>] [-v <µL>]
        [-b <bands>] [-R | -D] [-C]

Arguments:
    <samples>
        A description of each crude sample to purify.  Each argument can 
        contain several pieces of information, separated by colons as follows:

            name[:molecule[:conc[:volume]]]

        name:
            The name of the product.  If this is corresponds to a tag in the 
            FreezerBox database, default values for the other parameters will 
            be read from the database.

        molecule:
            What kind of molecule the product is: either "RNA" or "DNA" 
            (case-insensitive).  Use `--rna` or `--dna` if all samples are the 
            same molecule.

        conc:
            The concentration of the product.  This may include a unit.  If no 
            unit is specified, µg/µL is assumed.  Use `--conc` if all samples 
            are the same concentration.

        volume:
            The volume of the product to load on the gel, in µL.  Do not 
            include a unit.  Use `--volume` if all samples have the same 
            volume.
            
<%! from stepwise_mol_bio import hanging_indent %>\
Options:
    -p --preset <name>          [default: ${app.preset}]
        The default parameters to use.  The following presets are available:

        ${hanging_indent(app.preset_briefs, 8*' ')}

    -P --gel-preset <name>
        The default gel electrophoresis parameters to use.  See `sw gel -h` for 
        a list of available presets.

    -c --conc <float>
        The concentration of each sample.  This may include a unit.  If no unit 
        is specified, µg/µL is assumed.  This is superseded by concentrations 
        specified in the <samples> argument, but supersedes concentrations 
        found in the FreezerBox database.

    -v --volume <µL>
        The volume of each sample to load on the gel, in µL.  Do not include a 
        unit.  This is superseded by volumes specified in the <samples> 
        argument.

    -b --bands <names>
        A comma separated list of names identifying the bands to cut out of the 
        gel.  Typically these would be the lengths of the desired products, 
        e.g. "400 nt" or "1.5 kb".

    -R --rna
        Assume that each sample is RNA.  This is superseded by the <samples> 
        argument, but supersedes the FreezerBox database.

    -D --dna
        Assume that each sample is DNA.  This is superseded by the <samples> 
        argument, but supersedes the FreezerBox database.

    -C --no-cleanup
        Don't include the final spin-column purification step in the protocol.

Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.page_purify.default_preset:
        The default value for the `--preset` option.

    molbio.page_purify.presets:
        Named groups of default reaction parameters.  Typically each preset 
        corresponds to a particular kit or protocol.  See below for the various 
        settings that can be specified in each preset.

    molbio.page_purify.presets.<name>.conc
        The default value for the `--conc` option.  Note that if this option is 
        set, concentrations will never be read from the FreezerBox database.

    molbio.page_purify.presets.<name>.volume_uL
        The default value for the `--volume` option.

    molbio.page_purify.presets.<name>.molecule
        The default value for the `--rna`/`--dna` options.  This should either 
        be "RNA" or "DNA".

    molbio.page_purify.presets.<name>.gel_preset
        The default value for the `--gel-preset` option.

    molbio.page_purify.presets.<name>.bands
        The default value for the `--bands` option.

    molbio.page_purify.presets.<name>.cleanup_preset
        The default presets to use for the spin column cleanup step after 
        recovering the DNA/RNA from the gel.  See `sw spin_cleanup -h` for a 
        list of valid presets.  This option should be a dictionary with the 
        keys 'rna' and 'dna', each specifying the preset to use with samples of 
        the corresponding type.  Alternatively, this option can be 'false' to 
        indicate that the cleanup step should be skipped (see `--no-cleanup`).

Database:
    PAGE purification protocols can appear in the "Cleanups" column of a 
    FreezerBox database:

        page-purify [preset=<name>] [gel=<name>] [conc=<conc>] [volume=<µL>]
    
    preset=<name>
        See `--preset`.

    gel=<name>
        See `--gel-preset`.

    conc=<conc>
        See `--conc`.  Must include a unit.

    volume=<µL>
        See `--volume`.  Must include a unit.
"""
    __config__ = [
            DocoptConfig(),
            MakerArgsConfig(),
            PresetConfig(),
            StepwiseConfig('molbio.page_purify'),
    ]
    Sample = Sample
    preset_briefs = appcli.config_attr()
    config_paths = appcli.config_attr()

    presets = appcli.param(
            Key(StepwiseConfig, 'presets'),
            pick=merge_dicts,
    )
    preset = appcli.param(
            Key(DocoptConfig, '--preset'),
            Key(MakerArgsConfig, 'preset'),
            Key(StepwiseConfig, 'default_preset'),
    )
    samples = appcli.param(
            Key(DocoptConfig, '<samples>', cast=parse_samples),
            Method(lambda self: [Sample(name=x.tag) for x in self.products]),
    )
    default_stock_conc = appcli.param(
            Key(DocoptConfig, '--conc', cast=parse_conc),
            Key(MakerArgsConfig, 'conc', cast=Quantity.from_string),
            Key(PresetConfig, 'conc', cast=Quantity.from_string),
    )
    default_volume_uL = appcli.param(
            Key(DocoptConfig, '--volume', cast=parse_volume_uL),
            Key(MakerArgsConfig, 'volume', cast=parse_strict_volume_uL),
            Key(PresetConfig, 'volume_uL', cast=parse_volume_uL),
            default=5,
    )
    default_molecule = appcli.param(
            Key(DocoptConfig, '--dna', cast=lambda x: 'DNA'),
            Key(DocoptConfig, '--rna', cast=lambda x: 'RNA'),
            Key(PresetConfig, 'molecule'),
    )
    gel_preset = appcli.param(
            Key(DocoptConfig, '--gel-preset'),
            Key(MakerArgsConfig, 'gel'),
            Key(PresetConfig, 'gel_preset'),
    )
    gel_percent = appcli.param(
            Key(DocoptConfig, '--gel-percent'),
            default=None,
    )
    gel_run_volts = appcli.param(
            Key(DocoptConfig, '--gel-run-volts'),
            default=None,
    )
    gel_run_time_min = appcli.param(
            Key(DocoptConfig, '--gel-run-time-min'),
            default=None,
    )
    desired_bands = appcli.param(
            Key(DocoptConfig, '--bands', cast=comma_list),
            Key(PresetConfig, 'bands'),
            default_factory=list,
    )
    rna_cleanup_preset = appcli.param(
            Key(PresetConfig, 'cleanup_preset.rna'),
    )
    dna_cleanup_preset = appcli.param(
            Key(PresetConfig, 'cleanup_preset.dna'),
    )
    cleanup = appcli.param(
            Key(DocoptConfig, '--no-cleanup', cast=not_),
            Key(PresetConfig, 'cleanup_preset', cast=bool),
    )

    group_by = {
            'preset': group_by_identity,
            'gel_preset': group_by_identity,
    }
    merge_by = {
            'samples': join_lists,
    }

    def __init__(self, samples):
        self.samples = samples

    def get_protocol(self):
        self._bind_samples()

        p = stepwise.Protocol()
        p += self.gel_electrophoresis_steps
        p += self.gel_extraction_steps
        p += self.product_recovery_steps
        return p

    def get_gel_electrophoresis_steps(self):
        p = stepwise.Protocol()
        n = plural(self.samples)
        gel = Gel(self.gel_preset)

        mix = gel.sample_mix
        k = mix.volume / mix['sample'].volume

        for sample in self.samples:
            sample.stock_conc
            sample.sample_volume_per_lane_uL = min(
                    sample.volume_uL,
                    sample.volume_uL / ceil(sample.mass_ug / '20 µg'),
                    mix['sample'].volume.value or inf,
            )
            sample.load_volume_per_lane_uL = k * sample.sample_volume_per_lane_uL
            sample.num_lanes = int(ceil(sample.volume_uL / sample.sample_volume_per_lane_uL))

        conc_vol_groups = group_by_identity(
                self.samples,
                key=lambda x: (x.stock_conc, x.volume_uL),
        )
        for (conc, volume_uL), group in conc_vol_groups:
            prep_gel = Gel(self.gel_preset)
            prep_gel.num_samples = len(group)
            mix = prep_gel.sample_mix

            mix.hold_ratios.volume = k * volume_uL, 'µL'
            mix['sample'].name = ','.join(str(x.name) for x in group)
            mix['sample'].stock_conc = conc

            p += prep_gel.prep_step

        if x := self.gel_percent:
            gel.gel_percent = x
        if x := self.gel_run_volts:
            gel.run_volts = x
        if x := self.gel_run_time_min:
            gel.run_time_min = x

        try:
            gel.load_volume_uL = unanimous(
                    x.load_volume_per_lane_uL
                    for x in self.samples
            )
        except ValueError:
            gel.load_volume_uL = 0

            p += pl(f"Load the {n:sample/s} in the gel as follows:", u := ul())
            for sample in self.samples:
                num_lanes = f" in each of {m} lanes" if (m := sample.num_lanes) != 1 else ""
                u += f"{sample.load_volume_per_lane_uL:.2g} µL {sample.name}{num_lanes}"

        p += gel.run_step

        return p

    def get_gel_extraction_steps(self):
        p = stepwise.Protocol()
        names = self._cluster_names_by_molecule()

        desired = oxford_comma(x) if (x := self.desired_bands) else 'desired'
        n = plural(max(len(self.desired_bands), len(self.samples)))
        f = "Based on Fitzy's DNA PAGE purification protocol, [Nilson2013], and [Petrov2013]."

        p += pl(f"Cut the {desired} {n:band/s} out of the gel{p.add_footnotes(f)}.", ul(
            "Place the gel over a TLC plate.",
            "Use a UV light to visualize the RNA (dark spot).",
            "Consider visualizing remaining gel to ensure that all desired RNA was excised.",
        ))
        p += pl(f"Crush the gel {n:slice/s}.", ul(
            "Poke a hole in the bottom of a 0.65 mL tube with a 27 g needle.",
            "Place gel slice inside the 0.65 mL tube.",
            "Place the 0.65 mL tube inside a 1.5 mL tube.",
            "Spin at max speed for 5 min.",
        ))

        eb = elution_buffer()
        eb.hold_ratios.volume = 500 * len(self.samples), 'µL'
        p += pl(
                "Resuspend gel in 400 µL PAGE elution buffer:",
                eb,
        )

        if names['RNA']:
            p += f"Incubate {join_if(names['RNA'], names['DNA'])}at 4°C overnight with end-over-end mixing."

        if names['DNA']:
            p += f"Incubate {join_if(names['DNA'], names['RNA'])}at 55°C overnight with 800 rpm mixing."

        return p
            
    def get_product_recovery_steps(self):
        p = stepwise.Protocol()
        names = self._cluster_names_by_molecule()

        # Filter material:
        # - Corning has a good guide on which material to select:
        #   https://www.corning.com/catalog/cls/documents/selection-guides/t_filterselectionguide.pdf
        #
        # - I want 0.22 µm, because that's the standard size for filter 
        #   sterilizing biological buffers (that's not my application here, but 
        #   I can see myself wanting to do that).
        #
        # - I want cellulose acetate filters.  Nylon and cellulose nitrate have 
        #   high DNA binding, which will cause me to lose material.  The 
        #   downside to cellulose acetate is that it has a wetting agent that 
        #   will end up in the sample.  However, this will be removed by the 
        #   Zymo column in the subsequent step.
        #
        # - Product number: 8161 (non sterile)
        #
        # Centrifugation speed:
        # - Fitzy's DNA PAGE purification protocol calls for 4 min at 7000 rpm
        # - The Corning guide (see above) includes an agarose gel purification 
        #   protocol, which calls for 13,000g for 5-20 min.  But this protocol 
        #   has no incubation step, so I gather that the spin is supposed to 
        #   pull the solvent out of the gel.  I probably don't need to go so 
        #   fast, but why not go as fast as the columns can tolerate?
        f = "Nylon and cellulose nitrate have high DNA-binding: https://tinyurl.com/3pkyc8dr"
        p += pl(
                "Remove gel debris by spin filtration:",
                ul(
                    f"Load samples onto a 0.22 µm cellose-acetate Spin-X column{p.add_footnotes(f)}.",
                    "Spin at 13,000g for 5 min."
                ),
        )

        if self.cleanup:
            if names['RNA']:
                p += cleanup(self.rna_cleanup_preset, names['RNA'], names['DNA'])

            if names['DNA']:
                p += cleanup(self.dna_cleanup_preset, names['DNA'], names['RNA'])

        return p

    def _bind_samples(self):
        assert self.samples

        for sample in self.samples:
            sample.bind(self)

    def _cluster_names_by_molecule(self):
        names = {'RNA': [], 'DNA': []}
        for sample in self.samples:
            names[sample.molecule].append(sample.name)
        return names

    @classmethod
    def _solo_maker_factory(cls, product):
        app = super()._solo_maker_factory(product)
        app._bind_samples()
        return app

def elution_buffer():
    return stepwise.MasterMix("""\
            Reagent               Stock      Volume
            ===================  ======  ==========
            nuclease-free water          to 1000 µL
            Tris, pH 7.5            1 M       10 µL
            NaCl                    5 M      100 µL
            EDTA                 500 mM        2 µL
            SDS                     10%       10 µL
    """)

def join_if(items, cond):
    if not cond:
        return ''
    else:
        return stepwise.oxford_comma(items) + ' '

def cleanup(preset, products, show_products):
    sc = SpinCleanup(preset)
    sc.sample_volume_uL = 400

    if show_products:
        sc.products = products
        sc.show_product_names = True

    return sc.protocol


if __name__ == '__main__':
    PagePurify.main()
