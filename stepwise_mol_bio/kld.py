#!/usr/bin/env python3
# vim: tw=50

import stepwise, appcli, autoprop
from inform import plural
from stepwise_mol_bio import Main

@autoprop
class Kld(Main):
    """\
Circularize a linear DNA molecule using T4 DNA ligase, e.g. to reform a plasmid 
after inverse PCR.

Usage:
    kld <dna> [-n <int>]

Arguments:
    <dna>
        The name of the DNA molecule to circularize.

Options:
    -n --num-reactions <int>        [default: ${app.num_reactions}]
        The number of reactions to set up.
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]

    dna = appcli.param('<dna>')
    num_reactions = appcli.param(
            '--num-reactions',
            cast=lambda x: int(eval(x)),
            default=1,
    )

    def __init__(self, dna):
        self.dna = dna

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
        kld['DNA'].name = self.dna

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
    Kld.main()
