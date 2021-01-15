#!/usr/bin/env python3

import stepwise, appcli, autoprop
from inform import plural
from stepwise_mol_bio import Main, comma_set
from collections.abc import Sequence

def float_pair(x):
    if ',' in x:
        return tuple(map(float, x.split(',', maxsplit=1)))
    else:
        return float(x)

def pair(x, i):
    return x[i] if isinstance(x, Sequence) else x

@autoprop
class Anneal(Main):
    """\
Anneal two complementary oligos.

Usage:
    anneal <oligo_1> <oligo_2> [-n <num_rxns>] [-v <µL>] [options]

Arguments:
    <oligo_1> <oligo_2>
        The names of the two oligos to anneal.

Options:
    -n --num-rxns <num_rxns>            [default: ${app.num_reactions}]
        The number of reactions to set up.

    -v --volume <µL>                    [default: ${app.volume_uL}]
        The volume of each annealing reaction in µL.

    -c --oligo-conc <µM>
        The final concentration of each oligo in the reaction, in µM.  This 
        will also be the concentration of the annealed duplex, if the reaction 
        goes to completion.  The default is to use as much oligo as possible.

    -C --oligo-stock <µM[,µM]>          [default: ${app.oligo_stock_uM}]
        The stock concentrations of the oligos, in µM.  You can optionally use 
        a comma to specify different stock concentrations for the two oligos.

    -m --master-mix <reagents>          [default: ${','.join(app.master_mix)}]
        The reagents to include in the master mix.  The following reagents are 
        understood: '1' (the first oligo), '2' (the second oligo), or the name 
        of either oligo specified on the command line. To specify both 
        reagents, separate the two names with a comma.
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]
    oligo_1 = appcli.param('<oligo_1>')
    oligo_2 = appcli.param('<oligo_2>')
    num_reactions = appcli.param('--num-rxns', default=1, cast=int)
    volume_uL = appcli.param('--volume', default=4, cast=float)
    oligo_conc_uM = appcli.param('--oligo-conc', default=None, cast=float)
    oligo_stock_uM = appcli.param('--oligo-stock', default=100, cast=float_pair)
    master_mix = appcli.param('--master-mix', default={'1'}, cast=comma_set)

    def __init__(self, oligo_1, oligo_2):
        self.oligo_1 = oligo_1
        self.oligo_2 = oligo_2

    def get_reaction(self):
        rxn = stepwise.MasterMix.from_text("""\
            Reagent  Stock     Volume  MM?
            =======  =====  =========  ===
            water           to 4.0 µL  yes
            PBS      10x       0.4 µL  yes
        """)
        rxn.num_reactions = self.num_reactions
        rxn.hold_ratios.volume = self.volume_uL, 'µL'

        rxn['oligo1'].name = self.oligo_1
        rxn['oligo2'].name = self.oligo_2
        rxn['oligo1'].master_mix = bool({'1', self.oligo_1} & self.master_mix)
        rxn['oligo2'].master_mix = bool({'2', self.oligo_2} & self.master_mix)
        rxn['oligo1'].stock_conc = pair(self.oligo_stock_uM, 0), 'µM'
        rxn['oligo2'].stock_conc = pair(self.oligo_stock_uM, 1), 'µM'

        if self.oligo_conc_uM:
            rxn['oligo1'].hold_stock_conc.conc = self.oligo_conc_uM, 'µM'
            rxn['oligo2'].hold_stock_conc.conc = self.oligo_conc_uM, 'µM'

        else:
            V = rxn.get_free_volume_excluding('oligo1', 'oligo2')
            C1 = rxn['oligo1'].stock_conc
            C2 = rxn['oligo2'].stock_conc
            C12 = C1 + C2

            rxn['oligo1'].volume = V * (C2 / C12)
            rxn['oligo2'].volume = V * (C1 / C12)

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

if __name__ == '__main__':
    Anneal.main()
