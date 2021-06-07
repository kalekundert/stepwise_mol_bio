#!/usr/bin/env python3

import appcli
import autoprop

from stepwise import Quantity
from stepwise_mol_bio import Main, UsageError
from freezerbox import (
        ReagentConfig, MakerArgsConfig, mw_from_length, parse_conc, 
        parse_size_bp, parse_volume_uL, convert_conc_unit, iter_combos, 
        group_by_identity, group_by_cluster, unanimous, join_lists,
)
from appcli import Key, Method, DocoptConfig
from dataclasses import dataclass
from itertools import zip_longest
from more_itertools import one
from functools import partial
from collections import defaultdict
from inform import warn
from math import inf

# fragment:
#   Parameters for a single dsDNA construct, including name, concentration, 
#   molecular weight, volume, etc.
#
# assembly:
#   A list of fragments that will be assembled into a single plasmid.
#
# assemblies:
#   A list of assemblies that may be performed simultaneously, e.g. using 
#   master mixes.

def parse_assemblies_from_docopt(args):
    assemblies = args['<assemblies>']
    concs = parse_fragment_attrs(
            args.get('--conc', []),
            partial(parse_conc, default_unit='ng/µL'),
    )
    lengths = parse_fragment_attrs(
            args.get('--length', []),
            partial(parse_size_bp, default_unit='bp'),
    )
    return parse_assemblies_from_comma_strs(
            args['<assemblies>'],
            concs=concs,
            lengths=lengths,
    )

def parse_assemblies_from_freezerbox(fields):
    return [parse_assembly_from_comma_str(fields[1])]

def parse_assemblies_from_comma_strs(assembly_strs, concs=None, lengths=None):
    return [
            parse_assembly_from_comma_str(x, concs=concs, lengths=lengths)
            for x in assembly_strs
    ]

def parse_assembly_from_comma_str(assembly_str, concs=None, lengths=None):
    if concs is None: concs = defaultdict(lambda: None)
    if lengths is None: lengths = defaultdict(lambda: None)

    assembly = []
    for name in assembly_str.split(','):
        fragment = Fragment(
                name,
                conc=concs[name],
                length=lengths[name],
        )
        assembly.append(fragment)

    return assembly

def parse_fragment_attrs(args, parse_value=lambda x: x):
    named_values = {}
    default_value = None

    for arg in args:
        fields = arg.split('=', 1)
        if len(fields) == 1:
            default_value = parse_value(arg)
        else:
            key, value = fields
            named_values[key] = parse_value(value)

    all_values = defaultdict(lambda: default_value)
    all_values.update(named_values)
    return all_values

def bind_assemblies(app, assemblies):
    for assembly in assemblies:
        for frag in assembly:
            frag.bind(app)
    return assemblies

def add_fragments_to_reaction(
        rxn,
        assemblies,
        *,
        target_pmol,
        min_pmol,
        excess_insert=1,
        order=-1,
    ):

    class FragmentReagent:

        def __init__(self, all_frags):
            self.fragments = all_frags

            # name, master_mix
            names = [x.name for x in all_frags]
            try:
                self.name = unanimous(names)
                self.master_mix = True
            except ValueError:
                self.name = ','.join(names)
                self.master_mix = False

            # conc
            frags = [x for x in all_frags if x is not null_reagent]
            self.conc_nM = inf
            self.conc = None

            for frag in frags:
                if frag.conc_nM < self.conc_nM:
                    self.conc_nM = frag.conc_nM
                    self.conc = frag.conc

            def all_close(x, threshold=0.1):
                try:
                    a, b = max(x), min(x)
                    return (a - b) / a < 0.1
                except TypeError:
                    return False

            concs = [x.conc for x in frags]
            concs_nM = [x.conc_nM for x in frags]

            if not all_close(concs_nM) or not all_close(concs):
                self.conc = self.conc_nM, 'nM'

            # volume
            self.vol_uL = uL_from_pmol(target_pmol, self.conc_nM)

    null_reagent = Fragment('−')
    frag_reagents = [
            FragmentReagent(frags)
            for frags in zip_longest(*assemblies, fillvalue=null_reagent)
    ]
    for i, reagent in enumerate(frag_reagents):
        excess = 1 if i == 0 else excess_insert
        reagent.vol_uL *= excess

    # Make sure the maximum volume is not exceeded.
    total_vol_uL = sum(x.vol_uL for x in frag_reagents)
    max_vol_uL = rxn.free_volume.value 
    if max_vol_uL < total_vol_uL:
        k = max_vol_uL / total_vol_uL
        for reagent in frag_reagents:
            reagent.vol_uL *= k

    # Warn if any fragment is below the minimum.
    for reagent in frag_reagents:
        reagent.pmol = pmol_from_uL(reagent.vol_uL, reagent.conc_nM)
        if reagent.pmol < min_pmol:
            warn(f"using {reagent.pmol:.3f} pmol of {reagent.name}, {min_pmol:.3f} pmol recommended.")

    for reagent in frag_reagents:
        rxn[reagent.name].volume = reagent.vol_uL, 'µL'
        rxn[reagent.name].stock_conc = reagent.conc
        rxn[reagent.name].master_mix = reagent.master_mix
        rxn[reagent.name].order = order
        rxn[reagent.name].flags.add('fragment')
        rxn[reagent.name].pmol = reagent.pmol

    rxn.num_reactions = len(assemblies)

def uL_from_pmol(pmol, conc_nM):
    return 1e3 * pmol / conc_nM

def pmol_from_uL(uL, conc_nM):
    return uL * conc_nM / 1e3

ARGUMENT_DOC = """\
    <assemblies>...
        The names of the DNA fragments to assemble, separated by commas.  Each 
        argument defines a single assembly reaction.  You can specify any 
        number of arguments to setup multiple assemblies at once.  Any 
        fragments that are used in every assembly will be automatically 
        included in the master mix.

        The concentration of every fragment must be known.  If the fragment 
        name corresponds to an entry in the FreezerBox database, the 
        concentration specified in the database will be used by default.  
        Otherwise, you must specify concentrations using the `--conc` option.

        The first fragment in each assembly is considered the "backbone", and 
        all the remaining fragments are considered "inserts".  This is only 
        relevant for the `--excess-insert` option."""
OPTION_DOC = """\
    -v --volume <µL>                [default: ${app.volume_uL}]
        The volume of the complete assembly reaction.  This does not affect the 
        amount of DNA that will be used; the recomendations for that amount are 
        a fixed number of pmol per reaction.  But it may be necessary to 
        increase the volume if your DNA is dilute, or if you have a large 
        number of inserts.

    -c --conc <[name=]value>
        Either the concentration of a single fragment, or the default 
        concentration for all fragments.  The former applies if a name is 
        specified (followed by an equals sign); the latter applies 
        otherwise.

        The concentration value can optionally include a unit, e.g. nM, ng/µL, 
        etc.  If no unit is specified, ng/µL is assumed.  If the concentration 
        is not specified in terms of molarity, the length of the fragment(s) in 
        question must also be specified.  This information can be queried from 
        the FreezerBox database, or explicitly given via the `--length` option.

        Although concentrations can be read from the FreezerBox database, any 
        values specified on the command-line will override any specified in the 
        database.

    -l --length <[name=]value>
        Either the length of a single fragment, or the default length for all 
        fragments.  The former applies if a name is specified (followed by an 
        equals sign); the latter applies otherwise.

        The length can optionally include a unit, e.g. bp or kb.  If no unit is 
        specified, bp is assumed.

        Lengths only need to be specified for fragments that don't have 
        concentrations in molarity units.  The lengths will be used to estimate 
        to molecular weight of each fragment.  This assumes that each fragment 
        is dsDNA.  For fragments that are present in the FreezerBox database, 
        molecular wieght will be inferred by default from any sequence or 
        length information present in the database.  Lengths specified on the 
        command-line override any values present in the database.

    -x --excess-insert <ratio>      [default: ${app.excess_insert}]
        The molar-excess of each insert relative to the backbone.  Values 
        between 1-10 (e.g. 1-10x excess) are typical."""

@autoprop.cache
class Fragment:
    __config__ = [
            ReagentConfig(transform=one, tag_getter=lambda self: self.name),
    ]

    name = appcli.param()
    conc = appcli.param(
            Key(ReagentConfig),
    )
    mw = appcli.param(
            Key(ReagentConfig),
            Method(lambda self: mw_from_length(self.length)),
            default=None,
    )
    length = appcli.param(
            Key(ReagentConfig),
    )

    def __init__(self, name=None, *, conc=None, length=None, mw=None):
        if name: self.name = name
        if conc: self.conc = conc
        if length: self.length = length
        if mw: self.mw = mw

    def __repr__(self):
        attrs = 'name', 'conc', 'length', 'mw'
        attr_strs = [
                f'{attr}={value!r}'
                for attr in attrs
                if (value := getattr(self, attr, None)) is not None
        ]
        return f'Fragment({", ".join(attr_strs)})'

    def __eq__(self, other):
        undef = object()
        attrs = 'name', 'conc', 'length'
        return all(
                getattr(self, attr, undef) == getattr(other, attr, undef)
                for attr in attrs
        )

    def bind(self, app, force=False):
        if not hasattr(self, 'app') or force:
            self.app = app

    def get_db(self):
        return self.app.db

    def get_conc_nM(self):
        return convert_conc_unit(self.conc, self.mw, 'nM').value

    def del_conc_nM(self):
        pass

@autoprop
class Assembly(Main):
    __config__ = [
            DocoptConfig(),
            MakerArgsConfig(),
    ]

    Fragment = Fragment
    target_pmol_per_frag = 0.06
    min_pmol_per_frag = 0.02

    assemblies = appcli.param(
            Key(DocoptConfig, parse_assemblies_from_docopt),
            Key(MakerArgsConfig, parse_assemblies_from_freezerbox),
            get=bind_assemblies
    )
    volume_uL = appcli.param(
            Key(DocoptConfig, '--volume', cast=float),
            Key(MakerArgsConfig, 'volume', cast=parse_volume_uL),
            default=5,
    )
    excess_insert = appcli.param(
            Key(DocoptConfig, '--excess-insert', cast=float),
            default=2,
    )

    group_by = {
            'volume_uL': group_by_identity,
            'excess_insert': group_by_identity,
    }
    merge_by = {
            'assemblies': join_lists,
    }

    def __init__(self, assemblies):
        self.assemblies = assemblies

    def get_num_fragments(self):
        return max(len(x) for x in self.assemblies)

    def get_product_conc(self):
        min_pmol = min(
                x.pmol
                for x in self.reaction.iter_reagents_by_flag('fragment')
        )
        return Quantity(1e3 * min_pmol / self.volume_uL, 'nM')

    def get_dependencies(self):
        return {
                frag.name
                for assembly in self.assemblies
                for frag in assembly
        }

    def _add_fragments_to_reaction(self, rxn, order=-1):
        rxn.hold_ratios.volume = self.volume_uL, 'µL'
        rxn.extra_min_volume = '0.5 µL'
        add_fragments_to_reaction(
                rxn, self.assemblies,
                target_pmol=self.target_pmol_per_frag,
                min_pmol=self.min_pmol_per_frag,
                excess_insert=self.excess_insert,
                order=order,
        )
        return rxn

