#!/usr/bin/env python3

import stepwise, appcli, autoprop

from stepwise import StepwiseConfig, pl, ul
from stepwise_mol_bio import (
        Main, Argument, UsageError,
        bind_arguments, merge_names,
)
from freezerbox import (
        ReagentConfig, MakerConfig, QueryError,
        parse_conc_uM, parse_volume_uL, group_by_identity, join_lists,
)
from appcli import Key, Method, DocoptConfig
from inform import plural
from more_itertools import all_equal, flatten
from collections.abc import Sequence
from functools import partial

def parse_oligo_pairs_docopt(pair_strs):
    return [
            parse_oligo_pair_docopt(x)
            for x in pair_strs
    ]

def parse_oligo_pair_docopt(pair_str):
    pair = pair_str.split(',')

    if len(pair) != 2:
        raise UsageError(
                lambda e: f"expected 2 oligos, got {len(e.pair)}: {e.pair_str}",
                pair=pair,
                pair_str=pair_str,
        )

    return Anneal.Oligo(pair[0]), Anneal.Oligo(pair[1])

def parse_oligo_pair_freezerbox(fields):
    if len(fields.by_index) != 3:
        raise QueryError(
                lambda e: f"expected 2 oligos, got {len(e.fields.by_index) - 1}: '{e.fields}'",
                fields=fields,
        )

    cmd, o1, o2 = fields.by_index
    pair = Anneal.Oligo(o1), Anneal.Oligo(o2)
    return [pair]

def parse_oligo_stock_docopt(conc_str):
    concs = [float(x) for x in conc_str.split(',')]

    if len(concs) == 1:
        return concs[0]
    if len(concs) == 2:
        return tuple(concs)

    raise UsageError(
            lambda e: f"expected 1 or 2 stock concentrations, got {len(e.concs)}: {e.conc_str}",
            concs=concs,
            conc_str=conc_str,
    )

@autoprop.cache
class Anneal(Main):
    """\
Anneal two complementary oligos.

Usage:
    anneal <oligo_1,oligo_2>... [-n <num_rxns>] [-v <µL>] [options]

Arguments:
    <oligo_1,oligo_2>
        The names of the two oligos to anneal, separated by a comma.  Any 
        number of oligo pairs can be specified.

Options:
    -n --num-rxns <num_rxns>
        The number of reactions to set up.  By default, this is inferred from 
        the number of oligo pairs specified.

    -v --volume <µL>
        The volume of each annealing reaction in µL.

    -c --oligo-conc <µM>
        The final concentration of each oligo in the reaction, in µM.  This 
        will also be the concentration of the annealed duplex, if the reaction 
        goes to completion.  The default is to use as much oligo as possible.

    -C --oligo-stock <µM[,µM]>
        The stock concentrations of the oligos, in µM.  You can optionally use 
        a comma to specify different stock concentrations for the two oligos.

<%! from stepwise_mol_bio import hanging_indent %>\
Configuration:
    Default values for this protocol can be specified in any of the following 
    stepwise configuration files:

        ${hanging_indent(app.config_paths, 8)}

    molbio.anneal.volume_uL:
        The default value for the `--volume` option.

    molbio.anneal.stock_conc_uM:
        The default value for the `--oligo-stock` option.  Unlike the above 
        flag, though, this must be just a single number.
    
FreezerBox:
    Annealing reactions can appear in the "Synthesis" column of a FreezerBox 
    database.  The associated database entry will automatically be considered 
    dsDNA, e.g. for the purpose of molecular weight calculations.  The 
    following options can be specified:

        anneal <oligo_1> <oligo_2> [conc=<µM>]

    conc=<µM>
        See `--oligo-conc`.
"""
    __config__ = [
            DocoptConfig,
            MakerConfig,
            StepwiseConfig.setup('molbio.anneal')
    ]

    class Oligo(Argument, use_main_configs=True):
        __config__ = [ReagentConfig]
        stock_uM = appcli.param(
                Key(DocoptConfig, '--oligo-stock', cast=parse_oligo_stock_docopt),
                Key(ReagentConfig, 'conc_uM'),
                Key(StepwiseConfig, 'stock_conc_uM'),
        )

    oligo_pairs = appcli.param(
            Key(DocoptConfig, '<oligo_1,oligo_2>', cast=parse_oligo_pairs_docopt),
            Key(MakerConfig, parse_oligo_pair_freezerbox),
            get=partial(bind_arguments, iter=flatten),
    )
    num_reactions = appcli.param(
            Key(DocoptConfig, '--num-reactions'),
            Method(lambda self: len(self.oligo_pairs)),
    )
    volume_uL = appcli.param(
            Key(DocoptConfig, '--volume'),
            Key(MakerConfig, 'volume', cast=parse_volume_uL),
            Key(StepwiseConfig, 'volume_uL'),
            cast=float,
    )
    oligo_conc_uM = appcli.param(
            Key(DocoptConfig, '--oligo-conc', cast=float),
            Key(MakerConfig, 'conc', cast=parse_conc_uM),
            default=None,
    )
    config_paths = appcli.config_attr()

    group_by = {
            'volume_uL': group_by_identity,
            'oligo_conc_uM': group_by_identity,
    }
    merge_by = {
            'oligo_pairs': join_lists,
    }

    @classmethod
    def from_tags(cls, tag_pairs):
        oligo_pairs = [
                (cls.Oligo(a), cls.Oligo(b))
                for a, b in tag_pairs
        ]
        return cls(oligo_pairs)


    def __init__(self, oligo_pairs):
        self.oligo_pairs = oligo_pairs

    def get_reaction(self):
        rxn = stepwise.MasterMix.from_text("""\
            Reagent  Stock     Volume  MM?
            =======  =====  =========  ===
            water           to 4.0 µL  yes
            PBS      10x       0.4 µL  yes
        """)
        rxn.num_reactions = self.num_reactions
        rxn.hold_ratios.volume = self.volume_uL, 'µL'

        oligos_1, oligos_2 = zip(*self.oligo_pairs)
        rxn['oligo1'].name = merge_names(x.tag for x in oligos_1)
        rxn['oligo2'].name = merge_names(x.tag for x in oligos_2)
        rxn['oligo1'].master_mix = all_equal(oligos_1)
        rxn['oligo2'].master_mix = all_equal(oligos_2)
        rxn['oligo1'].stock_conc = min(x.stock_uM for x in oligos_1), 'µM'
        rxn['oligo2'].stock_conc = min(x.stock_uM for x in oligos_2), 'µM'

        if self.oligo_conc_uM:
            rxn['oligo1'].hold_stock_conc.conc = self.oligo_conc_uM, 'µM'
            rxn['oligo2'].hold_stock_conc.conc = self.oligo_conc_uM, 'µM'

        else:
            V = rxn.get_free_volume_excluding('oligo1', 'oligo2')
            C1 = rxn['oligo1'].stock_conc
            C2 = rxn['oligo2'].stock_conc
            C12 = C1 + C2

            rxn['oligo1'].volume = V * (C2 / C12)
            rxn['oligo2'].volume = V * (C1 / C12)

        return rxn

    def get_protocol(self):
        protocol = stepwise.Protocol()
        n = self.num_reactions

        protocol += pl(
                f"Setup {plural(n):# annealing reaction/s}:",
                self.reaction,
        )
        protocol += pl(
                f"Perform the {plural(n):annealing reaction/s}:",
                ul(
                    "Incubate at 95°C for 2 min.",
                    "Cool at room temperature.",
                ),
        )
        return protocol

    def get_product_conc(self):
        rxn = self.reaction
        return min(rxn['oligo1'].conc, rxn['oligo2'].conc)

    def get_product_volume(self):
        return self.reaction.volume

    def get_product_molecule(self):
        return 'dsDNA'

    def get_dependencies(self):
        tags = [(a.tag, b.tag) for a, b in self.oligo_pairs]
        return set(flatten(tags))

if __name__ == '__main__':
    Anneal.main()
