#!/usr/bin/env python3

import stepwise, appcli, autoprop

from stepwise import pl, table
from stepwise_mol_bio import Cleanup
from freezerbox import ProductConfig
from appcli import Key, DocoptConfig
from more_itertools import one

def parse_reactions_docopt(rxn_strs):
    return merge_reactions(
            parse_reaction_docopt(x)
            for x in rxn_strs
    )

def parse_reaction_docopt(rxn_str):
    plasmid_strs, primer_strs = rxn_str.split(':', 1)
    plasmids = plasmid_strs.split(',')
    primers = primer_strs.split(',')
    return {
            plasmid: primers
            for plasmid in plasmids
    }

def parse_reactions_freezerbox(product):
    primers = product.maker_args.by_index[1:]
    return {product.tag: primers}

def merge_reactions(rxns):
    # Maintain insertion order.
    out = {}
    for d in rxns:
        for k, v in d.items():
            out.setdefault(k, []).extend(v)
    return out


@autoprop.cache
class Sequence(Cleanup):
    """\
Send samples for Sanger sequencing.

Usage:
    sequence <plasmids:primers>...

Arguments:
    <plasmids:primers>
        Which plasmids to sequence with which primers.  You can specify this 
        argument any number of times.  You can also specify comma-separated 
        lists of plasmids or primers on either side of the colon.  For example, 
        "p1,p2:o1,o2" would be taken to mean that the plasmids p1 and p2 should 
        both be sequenced with the primers o1 and o2.
"""
    __config__ = [
            DocoptConfig,
            ProductConfig,
    ]

    reactions = appcli.param(
            Key(DocoptConfig, '<plasmids:primers>', cast=parse_reactions_docopt),
            Key(ProductConfig, parse_reactions_freezerbox),
    )

    merge_by = {
            'reactions': merge_reactions,
    }

    def __init__(self, reactions):
        self.reactions = reactions

    def get_protocol(self):
        p = stepwise.Protocol()

        if len(self.reactions) == 1:
            plasmid, primers = one(self.reactions.items())
            p += f"Sequence {plasmid} with {', '.join(primers)}."

        else:
            p += pl(
                    "Sequence the following plasmids:",
                    table(
                        header=['Plasmid', 'Primers'],
                        rows=[
                            [k, ','.join(v)]
                            for k, v in self.reactions.items()
                        ],
                    ),
            )

        return p

if __name__ == '__main__':
    Sequence.main()


