#!/usr/bin/env python3

import stepwise
import appcli
import autoprop

from stepwise import pl, ul
from stepwise_mol_bio import UsageError, Assembly
from stepwise_mol_bio.digest import NebRestrictionEnzymeDatabase
from stepwise_mol_bio._assembly import ARGUMENT_DOC, OPTION_DOC
from freezerbox import MakerArgsConfig, group_by_identity
from appcli import Key, DocoptConfig
from inform import plural

@autoprop
class GoldenGate(Assembly):
    __doc__ = f"""\
Perform a Golden Gate assembly reaction.

Usage:
    golden_gate <assemblies>... [-c <conc>]... [-l <length>]... [options]

Arguments:
{ARGUMENT_DOC}

Options:
{OPTION_DOC}

    -e --enzymes <type_IIS>         [default: ${{','.join(app.enzymes)}}]
        The name(s) of the Type IIS restriction enzyme(s) to use for the 
        reaction.  To use more than one enzyme, enter comma-separated names.  
        Only NEB enzymes are supported.

Database:
    Golden gate assemblies can appear in the "Synthesis" column of a FreezerBox 
    database:
        
        golden-gate <assembly> [enzymes=<type IIS>] [volume=<µL>]

    <assembly>
        See the <assemblies>... command-line argument.

    volume=<µL>
        See --volume.  You must include a unit.

    enzymes=<type IIS>
        See --enzymes.
"""

    enzymes = appcli.param(
            Key(DocoptConfig, '--enzymes'),
            Key(MakerArgsConfig, 'enzymes'),
            cast=lambda x: x.split(','),
            default=['BsaI-HFv2'],
    )

    group_by = {
            **Assembly.group_by,
            'enzymes': group_by_identity,
    }

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = '20 µL'

        rxn['T4 ligase buffer'].volume = '2.0 μL'
        rxn['T4 ligase buffer'].stock_conc = '10x'
        rxn['T4 ligase buffer'].master_mix = True
        rxn['T4 ligase buffer'].order = 2

        enz_uL = 0.5 if self.num_fragments <= 10 else 1.0

        rxn['T4 DNA ligase'].volume = enz_uL, 'µL'
        rxn['T4 DNA ligase'].stock_conc = '400 U/μL'
        rxn['T4 DNA ligase'].master_mix = True
        rxn['T4 DNA ligase'].order = 3

        enzyme_db = NebRestrictionEnzymeDatabase()

        for enzyme in self.enzymes:
            stock = enzyme_db[enzyme]['concentration'] / 1000
            rxn[enzyme].volume = enz_uL, 'µL'
            rxn[enzyme].stock_conc = stock, 'U/µL'
            rxn[enzyme].master_mix = True
            rxn[enzyme].order = 4

        return self._add_fragments_to_reaction(rxn)

    def get_protocol(self):
        # Maybe this should be the getter function for the assembly param...

        p = stepwise.Protocol()
        rxn = self.reaction
        n = rxn.num_reactions
        n_frags = self.num_fragments

        f = 'https://tinyurl.com/yaa5mqz5'
        p += pl(
                f"Setup the Golden Gate {plural(n):reaction/s}{p.add_footnotes(f)}:",
                rxn,
        )
        if n_frags <= 2:
            p += pl(
                    "Run the following thermocycler protocol:",
                    ul(
                        "37°C for 5 min",
                    ),
                    "Or, to maximize the number of transformants:",
                    ul(
                        "37°C for 60 min",
                        "60°C for 5 min",
                    ),
            )

        elif n_frags <= 4:
            p += pl(
                    "Run the following thermocycler protocol:",
                    ul(
                        "37°C for 60 min",
                        "60°C for 5 min",
                    ),
            )

        elif n_frags <= 10:
            p += pl(
                    "Run the following thermocycler protocol:",
                    ul(
                        "Repeat 30 times:",
                        ul(
                            "37°C for 1 min",
                            "16°C for 1 min",
                        ),
                        "60°C for 5 min",
                    ),
            )

        else:
            p += pl(
                    "Run the following thermocycler protocol:",
                    ul(
                        "Repeat 30 times:",
                        ul(
                            "37°C for 5 min",
                            "16°C for 5 min",
                        ),
                        "60°C for 5 min",
                    ),
            )

        return p

if __name__ == '__main__':
    GoldenGate.main()
