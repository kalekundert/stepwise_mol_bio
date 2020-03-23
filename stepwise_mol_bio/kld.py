#!/usr/bin/env python3
# vim: tw=50

"""\
Circularize a linear DNA molecule using T4 DNA ligase, e.g. to reform a plasmid 
after inverse PCR.

Usage:
    kld <num_reactions> 

Arguments:
    <num_reactions>
        The number of reactions to set up.
"""

import docopt
import stepwise
import autoprop
from inform import plural

@autoprop
class Kld:

    def __init__(self):
        self.num_reactions = 1

    @classmethod
    def from_docopt(cls, *docopt_args, **docopt_kwargs):
        args = docopt.docopt(*docopt_args, **docopt_kwargs)

        kld = cls()
        kld.num_reactions = eval(args['<num_reactions>'])

        return kld

    def get_reaction(self):
        kld = stepwise.MasterMix.from_text('''\
        Reagent               Stock        Volume  Master Mix
        ================  =========   ===========  ==========
        water                         to 10.00 μL         yes
        T4 ligase buffer        10x       1.00 μL         yes
        T4 PNK              10 U/μL       0.25 μL         yes
        T4 DNA ligase      400 U/μL       0.25 μL         yes
        DpnI                20 U/μL       0.25 μL         yes
        PCR product        50 ng/μL       1.50 μL
        ''')

        kld.num_reactions = self.num_reactions
        kld.extra_percent = 15

        return kld

    def get_protocol(self):
        protocol = stepwise.Protocol()
        protocol += f"""\
Run {plural(self.num_reactions):# ligation reaction/s}:

{self.reaction}

- Incubate at room temperature for 1h.
"""
        return protocol

if __name__ == '__main__':
    kld = Kld.from_docopt(__doc__)
    print(kld.protocol)
