#!/usr/bin/env python3

import stepwise
import autoprop

from stepwise import pl, ul
from stepwise_mol_bio import Assembly
from stepwise_mol_bio._assembly import ARGUMENT_DOC, OPTION_DOC
from inform import plural

@autoprop.cache
class Gibson(Assembly):
    __doc__ = f"""\
Perform a Gibson assembly reaction, using the NEB master mix (E2611).

Usage:
    gibson <assemblies>... [-c <conc>]... [-l <length>]... [options]

Arguments:
{ARGUMENT_DOC}

Options:
{OPTION_DOC}

Database:
    Gibson assemblies can appear in the "Synthesis" column of a FreezerBox 
    database:
        
        gibson <assembly> [volume=<µL>]

    <assembly>
        See the <assemblies>... command-line argument.

    volume=<µL>
        See --volume.  You must include a unit.
"""

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = '20 µL'

        rxn['Gibson master mix'].volume = rxn.volume / 2
        rxn['Gibson master mix'].stock_conc = "2x"
        rxn['Gibson master mix'].catalog_num = 'NEB E2611'
        rxn['Gibson master mix'].master_mix = True
        rxn['Gibson master mix'].order = 2

        return self._add_fragments_to_reaction(rxn)

    def del_reaction(self):
        pass

    def get_protocol(self):
        p = stepwise.Protocol()
        rxn = self.reaction
        n = rxn.num_reactions

        f = "https://tinyurl.com/ychbvkra"
        p += pl(
                f"Setup the Gibson assembly {plural(n):reaction/s}{p.add_footnotes(f)}:",
                rxn,
        )

        incubation_time = '15 min' if self.num_fragments <= 3 else '1h'
        p += f"Incubate at 50°C for {incubation_time}."
        return p

    def del_protocol(self):
        pass

if __name__ == '__main__':
    Gibson.main()

# vim: tw=50
