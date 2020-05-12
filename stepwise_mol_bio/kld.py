#!/usr/bin/env python3
# vim: tw=50

"""\
Circularize a linear DNA molecule using T4 DNA ligase, e.g. to reform a plasmid 
after inverse PCR.

Usage:
    kld <substrate> [-n <int>]

Arguments:
    <substrate>
        The name of the DNA molecule to circularize.

Options:
    -n --num-reactions <int>        [default: {0.num_reactions}]
        The number of reactions to set up.

"""

import docopt
import stepwise
import autoprop
from inform import plural
from stepwise_mol_bio import Main

@autoprop
class Kld(Main):
    num_reactions = 1

    def __init__(self, substrate):
        self.substrate = substrate

    @classmethod
    def from_docopt(cls, args):
        self = cls(args['<substrate>'])
        self.num_reactions = int(eval(args['--num-reactions']))
        return self

    def get_reaction(self):
        kld = stepwise.MasterMix.from_text('''\
        Reagent               Stock        Volume  Master Mix
        ================  =========   ===========  ==========
        water                         to 10.00 μL         yes
        T4 ligase buffer        10x       1.00 μL         yes
        T4 PNK              10 U/μL       0.25 μL         yes
        T4 DNA ligase      400 U/μL       0.25 μL         yes
        DpnI                20 U/μL       0.25 μL         yes
        DNA                50 ng/μL       1.50 μL
        ''')

        kld.num_reactions = self.num_reactions
        kld.extra_percent = 15
        kld['DNA'].name = self.substrate

        return kld

    def get_protocol(self):
        protocol = stepwise.Protocol()
        protocol += f"""\
Run {plural(self.num_reactions):# ligation reaction/s}:

{self.reaction}

- Incubate at room temperature for 1h.
"""
        return protocol

__doc__ = __doc__.format(Kld)

if __name__ == '__main__':
    Kld.main(__doc__)
