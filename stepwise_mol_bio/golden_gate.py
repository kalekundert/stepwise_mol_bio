#!/usr/bin/env python3

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
    -n --num-reactions <N>  [default: 1]
        The number of reactions to setup.

    -e --enzymes <type_IIS>
        The name(s) of the Type IIS restriction enzyme(s) to use for the 
        reaction.  To use more than one enzyme, enter comma-separated names.  
        The default is to use a single generic name.

    -m, --master-mix <bb,ins>  [default: ""]
        Indicate which fragments should be included in the master mix.  Valid 
        fragments are "bb" (for the backbone), "ins" (for all the inserts), 
        "1" (for the first insert), "2", "3", etc.

    -v, --reaction-volume <µL>  [default: 10]
        The volume of the complete assembly reaction.  You might want larger 
        reaction volumes if your DNA is dilute, or if you have a large number 
        of inserts.

    -d, --dna-volume <µL>
        The combined volume of backbone and insert DNA to use, in µL.  The 
        default is use the as much volume as possible, i.e. the full reaction 
        volume (see --reaction-volume) less the volumes of any enzymes and 
        buffers.  You might want to use less DNA if you're trying to conserve 
        material.

    -x, --excess-insert <ratio>  [default: 2]
        The fold-excess of each insert relative to the backbone.  The default 
        specifies that there will be twice as much of each insert as there 
        is backbone.
"""

import stepwise
import autoprop
from inform import plural
from stepwise_mol_bio import UsageError
from _assembly import Assembly, calc_fragment_volumes, FRAGMENT_DOC

__doc__ = __doc__.format(**globals())

@autoprop
class GoldenGate(Assembly):

    def __init__(self):
        super().__init__()
        self.enzymes = ['Golden Gate enzyme']
        self.excess_insert = 2

    @classmethod
    def from_docopt(cls, args):
        self = super().from_docopt(args)
        if x := args['--enzymes']: self.enzyme = x.split(',')
        self.excess_insert = float(args['--excess-insert'])
        return self

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = self.reaction_volume_uL, 'µL'
        rxn.num_reactions = self.num_reactions or 1

        rxn['T4 ligase buffer'].volume = '1.0 μL'
        rxn['T4 ligase buffer'].stock_conc = '10x'
        rxn['T4 ligase buffer'].master_mix = True
        rxn['T4 ligase buffer'].order = 2

        rxn['T4 DNA ligase'].volume = '0.5 μL'
        rxn['T4 DNA ligase'].stock_conc = '400 U/μL'
        rxn['T4 DNA ligase'].master_mix = True
        rxn['T4 DNA ligase'].order = 3

        for enzyme in self.enzymes:
            rxn[enzyme].volume = '0.5 μL'
            rxn[enzyme].master_mix = True
            rxn[enzyme].order = 4

        calc_fragment_volumes(
                self.fragments,
                vol_uL=self.pick_dna_volume_uL(rxn.free_volume.value),
                excess_insert=self.excess_insert,
        )

        for i, frag in enumerate(self.fragments):
            rxn[frag.name].volume = self.pick_frag_volume_uL(frag), 'µL'
            rxn[frag.name].stock_conc = frag.conc.value, frag.conc.unit
            rxn[frag.name].master_mix = self.is_frag_in_master_mix(i)
            rxn[frag.name].order = 1

        return rxn

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

        p.footnotes[1] = """\
https://international.neb.com/protocols/2018/10/02/golden-gate-assembly-protocol-for-using-neb-golden-gate-assembly-mix-e1601
"""
        return p

    def get_recommended_pmol_per_frag(self):
        # These recommendations are actually for Gibson assemblies [1].  I'm 
        # using them here because I couldn't find clear recommendations 
        # specifically for Golden Gate assemblies.  NEB's Golden Gate protocol 
        # [2] uses 0.05 pmol backbone, but that is neither a minimum nor a 
        # maximum.
        #
        # [1] https://www.neb.com/protocols/2012/09/25/gibson-assembly-master-mix-assembly
        # [2] https://international.neb.com/protocols/2018/10/02/golden-gate-assembly-protocol-for-using-neb-golden-gate-assembly-mix-e1601

        if len(self.fragments) <= 3:
            return 0.5, 0.02
        else:
            return 1.0, 0.2

if __name__ == '__main__':
    GoldenGate.main(__doc__)

# vim: tw=50

