#!/usr/bin/env python3

import stepwise
import autoprop
import appcli

from stepwise import pl, ul
from stepwise_mol_bio import Assembly
from stepwise_mol_bio._assembly import ARGUMENT_DOC, OPTION_DOC
from freezerbox import MakerConfig, group_by_identity, parse_bool
from appcli import Key, DocoptConfig
from inform import plural

@autoprop.cache
class Ligate(Assembly):
    __doc__ = f"""\
Assemble restriction-digested DNA fragments using T4 DNA ligase.

Usage:
    ligate <assemblies>... [-c <conc>]... [-l <length>]... [options]

Arguments:
{ARGUMENT_DOC}

Options:
{OPTION_DOC}

    -k --kinase
        Add T4 polynucleotide kinase (PNK) to the reaction.  This is necessary 
        to ligate ends that are not already 5' phosphorylated (e.g. annealed 
        oligos, PCR products).

Database:
    Ligation reactions can appear in the "Synthesis" column of a FreezerBox 
    database:
        
        ligate <assembly> [volume=<µL>] [kinase=<bool>]

    <assembly>
        See the <assemblies>... command-line argument.

    volume=<µL>
        See --volume.  You must include a unit.

    kinase=<bool>
        See --kinase.  Specify "yes" or "no".
"""

    excess_insert = appcli.param(
            Key(DocoptConfig, '--excess-insert', cast=float),
            default=3,
    )
    use_kinase = appcli.param(
            Key(DocoptConfig, '--kinase'),
            Key(MakerConfig, 'kinase', cast=parse_bool),
            default=False,
    )

    group_by = {
            **Assembly.group_by,
            'use_kinase': group_by_identity,
    }

    def get_reaction(self):
        rxn = stepwise.MasterMix.from_text('''\
        Reagent               Stock        Volume  Master Mix
        ================  =========   ===========  ==========
        water                         to 20.00 μL         yes
        T4 ligase buffer        10x       2.00 μL         yes
        T4 DNA ligase      400 U/μL       1.00 μL         yes
        T4 PNK              10 U/μL       1.00 μL         yes
        ''')
        if not self.use_kinase:
            del rxn['T4 PNK']

        return self._add_fragments_to_reaction(rxn)

    def get_protocol(self):
        p = stepwise.Protocol()
        rxn = self.reaction

        p += pl(
                f"Setup {plural(rxn.num_reactions):# ligation reaction/s}{p.add_footnotes('https://tinyurl.com/y7gxfv5m')}:",
                rxn,
        )
        p += pl(
                "Incubate at the following temperatures:",
                ul(
                    "25°C for 15 min",
                    "65°C for 10 min",
                ),
        )
        return p


if __name__ == '__main__':
    Ligate.main()
