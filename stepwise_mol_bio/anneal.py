#!/usr/bin/env python3

"""\
Anneal linker-N and mRNA prior to ligation.

Usage:
    anneal <oligo_1> <oligo_2> [-n <num_rxns>] [-v <µL>] [options]

Arguments:
    <oligo_1> <oligo_2>
        The names of the two oligos to anneal.

Options:
    -n --num-rxns <num_rxns>            [default: {0.num_reactions}]
        The number of reactions to set up.

    -v --volume <µL>                    [default: {0.volume_uL}]
        The volume of each annealing reaction in µL.

    -c --oligo-conc <µM>                [default: {0.oligo_conc_uM}]
        The final concentration of each oligo in the reaction, in µM.  This 
        will also be the concentration of the annealed duplex, if the reaction 
        goes to completion.

    -C --oligo-stock <µM>                [default: {0.oligo_stock_uM}]
        The stock concentrations of the oligos, in µM.

    -m --master-mix <reagents>          [default: {1}]
        The reagents to include in the master mix.  The following reagents are 
        understood: '1' (the first oligo), '2' (the second oligo), or the name 
        of either oligo specified on the command line. To specify both 
        reagents, separate the two names with a comma.
"""

import stepwise, docopt, autoprop
from inform import plural
from stepwise_mol_bio import Main

@autoprop
class Anneal(Main):
    num_reactions = 1
    volume_uL = 4
    oligo_conc_uM = 45
    oligo_stock_uM = 100
    master_mix = set()

    def __init__(self, oligo_1, oligo_2):
        self.oligo_1 = oligo_1
        self.oligo_2 = oligo_2

    @classmethod
    def from_docopt(cls, args):
        self = cls(args['<oligo_1>'], args['<oligo_2>'])
        self.num_reactions = int(args['--num-rxns'])
        self.volume = float(args['--volume'])
        self.oligo_conc_uM = float(args['--oligo-conc'])
        self.oligo_stock_uM = float(args['--oligo-stock'])
        self.master_mix = {*args['--master-mix'].split(',')}
        return self

    def get_reaction(self):
        rxn = stepwise.MasterMix.from_text("""\
            Reagent  Stock     Volume  MM?
            =======  =====  =========  ===
            water           to 4.0 µL  yes
            PBS      10x       0.4 µL  yes
            oligo1   10 µM     0.5 µL
            oligo2   10 µM     0.5 µL
        """)
        rxn.num_reactions = self.num_reactions
        rxn.hold_ratios.volume = self.volume_uL, 'µL'

        rxn['oligo1'].name = self.oligo_1
        rxn['oligo2'].name = self.oligo_2
        rxn['oligo1'].master_mix = bool({'1', self.oligo_1} & self.master_mix)
        rxn['oligo2'].master_mix = bool({'2', self.oligo_2} & self.master_mix)
        rxn['oligo1'].stock_conc = self.oligo_stock_uM, 'µM'
        rxn['oligo2'].stock_conc = self.oligo_stock_uM, 'µM'
        rxn['oligo1'].hold_stock_conc.conc = self.oligo_conc_uM, 'µM'
        rxn['oligo2'].hold_stock_conc.conc = self.oligo_conc_uM, 'µM'

        return rxn

    def get_protocol(self):
        protocol = stepwise.Protocol()
        n = self.num_reactions

        protocol += f"""\
Setup {plural(n):# annealing reaction/s}:

{self.reaction}
"""
        protocol += f"""\
Perform the {plural(n):annealing reaction/s}:

- Incubate at 95°C for 2 min.
- Cool at room temperature.
"""
        return protocol

__doc__ = __doc__.format(Anneal, ','.join(Anneal.master_mix))

if __name__ == '__main__':
    Anneal.main(__doc__)
