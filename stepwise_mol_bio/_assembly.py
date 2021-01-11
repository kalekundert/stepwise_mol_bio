#!/usr/bin/env python3

import appcli
import autoprop
from inform import warn
from appcli import Key, DocoptConfig
from dataclasses import dataclass
from stepwise_mol_bio import Main, UsageError

FRAGMENT_DOC = """\
    <backbone> <inserts>
        The DNA fragments to assemble.  For each fragment, the following 
        information can be specified:
        
        - A name (optional): This is how the fragment will be referred to in 
          the protocol, but will not have any effect other than that.  By 
          default, a generic name like "Insert #1" will be chosen.

        - A concentration (required): This is used to calculate how much of 
          each fragment needs to be added to get the ideal backbone:insert 
          ratio.  You may specify a unit (e.g. nM), but ng/µL will be assumed 
          if you don't.  You may also specify the concentration as just "PCR" 
          (no other value or units).  This indicates that the fragment is the 
          product of an PCR reaction, and is assumed to be 50 ng/μL.

        - A length (required if concentration is in ng/µL): The length of the 
          fragment in bp.  This is used to convert the above concentration into
          a molarity.

        For each individual fragment, specify whichever of the above fields are 
        relevant, separated by colons.  In other words:
        
            [<name>:]<conc>[:<length>]"""

OPTION_DOC = """\
    -n --num-reactions <int>        [default: ${app.num_reactions}]
        The number of reactions to setup.

    -m, --master-mix <bb,ins>       [default: ${app.master_mix_str}]
        Indicate which fragments should be included in the master mix.  Valid 
        fragments are "bb" (for the backbone), "ins" (for all the inserts), 
        "1" (for the first insert), "2", "3", etc.

    -v, --reaction-volume <µL>      [default: ${app.volume_uL}]
        The volume of the complete assembly reaction.  You might want larger 
        reaction volumes if your DNA is dilute, or if you have a large number 
        of inserts.

    -x, --excess-insert <ratio>     [default: ${app.excess_insert}]
        The molar-excess of each insert relative to the backbone.  Values 
        between 1-10 (e.g. 1-10x excess) are typical."""

class Assembly(Main):
    __config__ = [
            DocoptConfig(
                usage_getter=lambda self: self._get_docopt_usage(),
            ),
    ]

    fragments = appcli.param(
            lambda d: parse_fragments([d['<backbone>'], *d['<inserts>']]),
    )
    num_reactions = appcli.param(
            Key(DocoptConfig, '--num-reactions'),
            cast=int,
            default=1,
    )
    volume_uL = appcli.param(
            Key(DocoptConfig, '--reaction-volume'),
            cast=float,
            default=5,
    )
    master_mix = appcli.param(
            Key(DocoptConfig, '--master-mix'),
            cast=lambda x: frozenset(x.split(',')),
            default=frozenset(),
    )
    excess_insert = appcli.param(
            Key(DocoptConfig, '--excess-insert'),
            cast=float,
            default=2,
    )

    target_pmol_per_frag = 0.06
    min_pmol_per_frag = 0.02

    @property
    def master_mix_str(self):
        return ','.join(self.master_mix)

    def _get_docopt_usage(self):
        return self.__doc__.format(self, **globals())

    def _add_fragments_to_reaction(self, rxn):
        calc_fragment_volumes(
                self.fragments,
                target_pmol=self.target_pmol_per_frag,
                min_pmol=self.min_pmol_per_frag,
                max_vol_uL=rxn.free_volume.value,
                excess_insert=self.excess_insert,
        )

        for i, frag in enumerate(self.fragments):
            rxn[frag.name].volume = frag.vol_uL, 'µL'
            rxn[frag.name].stock_conc = frag.conc.value, frag.conc.unit
            rxn[frag.name].master_mix = self._is_frag_in_master_mix(i)
            rxn[frag.name].order = 1

        rxn.hold_ratios.volume = self.volume_uL, 'µL'
        rxn.num_reactions = self.num_reactions
        rxn.extra_min_volume = '0.5 µL'

        return rxn

    def _is_frag_in_master_mix(self, i):
        if i == 0:
            return 'bb' in self.master_mix
        else:
            return 'ins' in self.master_mix or str(i) in self.master_mix

@dataclass
class Fragment:
    name: str
    conc_nM: float
    vol_uL: float = None

@dataclass
class Concentration:
    value: float
    unit: str

def parse_fragments(frag_strs):
    """
    Parse fragments from colon-separated strings, i.e. that could 
    be specified on the command-line.

    See the usage text for a description of the syntax of this string.

    Note that if only two fields are specified, they could refer to a name and 
    a concentration, or a concentration and a size.  This is resolved by 
    attempting to interpret each field as a concentration.  Note that this will 
    break if given a name that could be interpreted as a concentration, so 
    don't do that.
    """

    fragments = []

    if len(frag_strs) < 2:
        raise UsageError("must specify at least two fragments")

    for i, frag_str in enumerate(frag_strs):
        fields = frag_str.split(':')
        frag_name = default_fragment_name(i)
        frag_size = None

        if len(fields) == 3:
            frag_name = fields[0]
            frag_conc = conc_from_str(fields[1])
            frag_size = int(fields[2])

        elif len(fields) == 2:
            # Is this (name, conc) or (conc, size)?
            try:
                frag_conc = conc_from_str(fields[0])
                frag_size = int(fields[1])

            except ValueError:
                frag_name = fields[0]
                frag_conc = conc_from_str(fields[1])

        elif len(fields) == 1:
            try:
                frag_conc = conc_from_str(fields[0])

            except ValueError:
                frag_name = fields[0]
                frag_conc = conc_from_str('PCR')

        else:
            raise UsageError("cannot parse fragment '{fragment_str}'")

        if frag_conc.unit == 'ng/µL' and frag_size is None:
            raise UsageError(f"'{frag_str}' specifies a concentration in ng/µL, so the size of the fragment must also be specified (e.g. '{frag_str}:<size>')")

        frag_nM = nM_from_conc(frag_conc, frag_size)
        frag = Fragment(frag_name, frag_nM)
        frag.conc = frag_conc
        fragments.append(frag)

    return fragments

def calc_fragment_volumes(
        frags,
        target_pmol,
        min_pmol,
        max_vol_uL,
        excess_insert=1,
):
    import numpy as np

    # Calculate the ideal amount of each fragment.
    for i, frag in enumerate(frags):
        excess = 1 if i == 0 else excess_insert
        frag.vol_uL = uL_from_pmol(excess * target_pmol, frag.conc_nM)

    # Make sure the maximum volume is not exceeded.
    total_vol_uL = sum(x.vol_uL for x in frags)
    if max_vol_uL < total_vol_uL:
        k = max_vol_uL / total_vol_uL
        for frag in frags:
            frag.vol_uL *= k

    # Warn if any fragment is below the minimum.
    for frag in frags:
        best_pmol = pmol_from_uL(frag.vol_uL, frag.conc_nM)
        if best_pmol < min_pmol:
            warn(f"using {best_pmol:.3f} pmol of {frag.name}, {min_pmol:.3f} pmol recommended.")

def default_fragment_name(i):
    return "backbone" if i == 0 else f"insert #{i}"

def conc_from_str(x):
    import re

    if x.upper().strip() == 'PCR':
        return Concentration(50, 'ng/µL')

    value_pattern = '[0-9.]+'
    unit_pattern = 'ng/[uµ]L|[muµnpf]M'
    conc_pattern = rf'({value_pattern})(\s*({unit_pattern}))?$'

    match = re.match(conc_pattern, x.strip())
    if not match:
        raise ValueError(f"could not interpret '{x}' as a concentration.")

    value, _, unit = match.groups()
    unit = (unit or 'ng/µL').replace('u', 'µ')
    return Concentration(float(value), unit)

def nM_from_conc(conc, num_bp):
    multiplier = {
            'mM': 1e9 / 1e3,
            'µM': 1e9 / 1e6,
            'nM': 1e9 / 1e9,
            'pM': 1e9 / 1e12,
            'fM': 1e9 / 1e15,
    }

    # https://www.neb.com/tools-and-resources/usage-guidelines/nucleic-acid-data
    if num_bp:
        multiplier['ng/µL'] = 1e6 / (650 * num_bp)

    try:
        return conc.value * multiplier[conc.unit]
    except KeyError:
        raise ValueError(f"cannot convert '{conc.unit}' to nM.")

def uL_from_pmol(pmol, conc_nM):
    return 1e3 * pmol / conc_nM

def pmol_from_uL(uL, conc_nM):
    return uL * conc_nM / 1e3


def test_conc_from_str():
    from pytest import approx, raises
    c = Concentration

    assert conc_from_str('1') == c(1.0, 'ng/µL')
    assert conc_from_str('10') == c(10.0, 'ng/µL')
    assert conc_from_str('1.0') == c(1.0, 'ng/µL')
    assert conc_from_str('1 ng/uL') == c(1.0, 'ng/µL')
    assert conc_from_str('1 ng/µL') == c(1.0, 'ng/µL')
    assert conc_from_str('1 fM') == c(1.0, 'fM')
    assert conc_from_str('1 pM') == c(1.0, 'pM')
    assert conc_from_str('1 nM') == c(1.0, 'nM')
    assert conc_from_str('1 uM') == c(1.0, 'µM')
    assert conc_from_str('1 µM') == c(1.0, 'µM')

    assert conc_from_str(' 1 ') == c(1.0, 'ng/µL')
    assert conc_from_str(' 1 nM') == c(1.0, 'nM')
    assert conc_from_str('1 nM ') == c(1.0, 'nM')
    assert conc_from_str(' 1  nM ') == c(1.0, 'nM')

    assert conc_from_str('PCR') == c(50.0, 'ng/µL')
    assert conc_from_str('pcr') == c(50.0, 'ng/µL')

    with raises(ValueError, match="''"):
        conc_from_str('')
    with raises(ValueError, match='xxx'):
        conc_from_str('1 xxx')

def test_nM_from_conc():
    from pytest import approx, raises
    c = Concentration

    assert nM_from_conc(c(1, 'mM'), None) == approx(1e6)
    assert nM_from_conc(c(1, 'µM'), None) == approx(1e3)
    assert nM_from_conc(c(1, 'nM'), None) == approx(1e0)
    assert nM_from_conc(c(1, 'pM'), None) == approx(1e-3)
    assert nM_from_conc(c(1, 'fM'), None) == approx(1e-6)

    assert nM_from_conc(c(1, 'ng/µL'), 100) == approx(1e6/(650 * 100))
    assert nM_from_conc(c(1, 'ng/µL'), 1000) == approx(1e6/(650 * 1000))

def test_fragments_from_strs():
    from pytest import approx, raises
    f = Fragment

    ## 0 fragments
    with raises(UsageError):
        fragments_from_strs('')

    ## 1 fragment
    with raises(UsageError):
        fragments_from_strs(['30nM'])

    ## 2 fragments
     # 3 arguments
    assert fragments_from_strs(['30nM', 'Gene:60:1000']) == [
            f('Backbone', 30),
            f('Gene', approx((60 * 1e6) / (650 * 1000))),
    ]
     # 2 arguments: name, conc
    assert fragments_from_strs(['30nM', 'Gene:60nM']) == [
            f('Backbone', 30),
            f('Gene', 60),
    ]
    with raises(UsageError):
        fragments_from_strs(['30nM', 'Gene:60'])

     # 2 arguments: conc, size
    assert fragments_from_strs(['30nM', '60:1000']) == [
            f('Backbone', 30.0),
            f('Insert #1', approx((60 * 1e6) / (650 * 1000))),
    ]
     # 1 argument: conc
    assert fragments_from_strs(['30nM', '60nM']) == [
            f('Backbone', 30),
            f('Insert #1', 60),
    ]
    with raises(UsageError):
        fragments_from_strs('30nM:60')

    ## 3 fragments
    assert fragments_from_strs(['30nM', '60nM', '61nM']) == [
            f('Backbone', 30),
            f('Insert #1', 60),
            f('Insert #2', 61),
    ]

