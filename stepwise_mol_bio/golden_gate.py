#!/usr/bin/env python3

import stepwise
import appcli
import autoprop
from inform import plural
from stepwise_mol_bio import UsageError
from stepwise_mol_bio.digest import NebRestrictionEnzymeDatabase
from _assembly import Assembly

@autoprop
class GoldenGate(Assembly):
    """\
Perform a Golden Gate assembly reaction.

Usage:
    golden_gate <backbone> <inserts>... [options]

Arguments:
{FRAGMENT_DOC}

        The only difference between the backbone and the fragments is that the 
        backbone is typically present at half the concentration of the 
        inserts, see --excess-insert.

Options:
{OPTION_DOC}

    -e --enzymes <type_IIS>         [default: {0.enzymes_str}]
        The name(s) of the Type IIS restriction enzyme(s) to use for the 
        reaction.  To use more than one enzyme, enter comma-separated names.  
        The default is to use a single generic name.
"""

    enzymes = appcli.param(
            '--enzymes',
            cast=lambda x: x.split(','),
            default=['BsaI-HFv2'],
    )

    def get_enzymes_str(self):
        # This is why __doc__ should be a Jinja template...
        return ','.join(self.enzymes)

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = '20 µL'

        rxn['T4 ligase buffer'].volume = '2.0 μL'
        rxn['T4 ligase buffer'].stock_conc = '10x'
        rxn['T4 ligase buffer'].master_mix = True
        rxn['T4 ligase buffer'].order = 2

        enz_uL = 0.5 if len(self.fragments) <= 10 else 1.0

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
        p = stepwise.Protocol()
        n = self.num_reactions

        p += f"""\
Setup the Golden Gate {plural(n):reaction/s} [1]:

{self.reaction}
"""
        if len(self.fragments) == 2:
            p += f"""\
Run the following thermocycler protocol:

- 37°C for 5 min

Or, to maximize the number of transformants:

- 37°C for 60 min
- 60°C for 5 min
"""
        elif len(self.fragments) <= 4:
            p += f"""\
Run the following thermocycler protocol:

- 37°C for 60 min
- 60°C for 5 min
"""
        elif len(self.fragments) <= 10:
            p += f"""\
Run the following thermocycler protocol:

- Repeat 30 times:
  - 37°C for 1 min
  - 16°C for 1 min
- 60°C for 5 min
"""
        else:
            p += f"""\
Run the following thermocycler protocol:

- Repeat 30 times:
  - 37°C for 5 min
  - 16°C for 5 min
- 60°C for 5 min
"""
        p += """\
Transform 2 µL.
"""
        p.footnotes[1] = """\
https://preview.tinyurl.com/yaa5mqz5
"""
        return p

if __name__ == '__main__':
    GoldenGate.main()

# vim: tw=50

