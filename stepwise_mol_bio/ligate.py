#!/usr/bin/env python3
# vim: tw=50

"""\
Assemble restrictions digested DNA fragments using T4 DNA ligase.

Usage:
    ligate <backbone> <inserts>... [-n <int>] [-v <µL>] [-k] [options]

Arguments:
{FRAGMENT_DOC}

Options:
{OPTION_DOC}

    -k --kinase
        Add T4 polynucleotide kinase (PNK) to the reaction.  This is necessary 
        to ligate ends that are not already 5' phosphorylated (e.g. annealed 
        oligos, PCR products).
"""

import docopt, stepwise, autoprop
from inform import plural
from _assembly import Assembly, format_docstring

@autoprop
class Ligate(Assembly):
    excess_insert = 3
    use_kinase = False

    @classmethod
    def from_docopt(cls, args):
        self = super().from_docopt(args)
        self.use_kinase = args['--kinase']
        return self

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
        p= stepwise.Protocol()
        p+= f"""\
Run {plural(self.num_reactions):# ligation reaction/s} [1]:

{self.reaction}
"""
        p += """\
Incubate at the following temperatures:

- 25°C for 15 min
- 65°C for 10 min
"""
        p += """\
Transform 2 µL.
"""
        p.footnotes[1] = """\
https://preview.tinyurl.com/y7gxfv5m
"""
        return p

__doc__ = format_docstring(Ligate, __doc__)

if __name__ == '__main__':
    Ligate.main(__doc__)
