#!/usr/bin/env python3

"""\
Perform a Gibson assembly reaction, using the NEB master mix (E2611).

Usage:
    gibson_assembly <backbone> <inserts>... [options]

Arguments:
{FRAGMENT_DOC}

Options:
    -n --num-reactions <N>  [default: 1]
        The number of reactions to setup.

    -m, --master-mix <bb,ins>  [default: ""]
        Indicate which fragments should be included in the master mix.  Valid 
        fragments are "bb" (for the backbone), "ins" (for all the inserts), 
        "1" (for the first insert), "2", "3", etc.

    -v, --reaction-volume <µL>  [default: 10]
        The volume of the complete assembly reaction.  You might want larger 
        reaction volumes if your DNA is dilute, or if you have a large number 
        of inserts.

    -d, --dna-volume <µL>
        The maximum combined volume of backbone and insert DNA to use, in µL.  
        The default is use the quantities of DNA recommended by NEB, i.e. 0.2 
        pmol/fragment for 3 or fewer fragments, or as much as possible if those 
        quantities can't be reached.  You might want to use less DNA if you're 
        trying to conserve material.
"""

import stepwise
import autoprop
from inform import plural
from _assembly import Assembly, calc_fragment_volumes, FRAGMENT_DOC

__doc__ = __doc__.format(**globals())

@autoprop
class Gibson(Assembly):

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = self.reaction_volume_uL, 'µL'
        rxn.num_reactions = self.num_reactions or 1

        rxn['Gibson master mix (NEB E2611)'].volume = rxn.volume / 2
        rxn['Gibson master mix (NEB E2611)'].stock_conc = "2x"
        rxn['Gibson master mix (NEB E2611)'].master_mix = True
        rxn['Gibson master mix (NEB E2611)'].order = 2

        calc_fragment_volumes(
                self.fragments,
                vol_uL=self.pick_dna_volume_uL(rxn.free_volume.value),
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
        incubation_time = '15 min' if len(self.fragments) <= 3 else '1h'

        p += f"""\
Setup the Gibson assembly {plural(n):reaction/s} [1]:

{self.reaction}
"""
        p += f"""\
Incubate at 50°C for {incubation_time}.
"""
        p.footnotes[1] = "https://www.neb.com/protocols/2012/09/25/gibson-assembly-master-mix-assembly"
        return p

    def get_recommended_pmol_per_frag(self):
        if len(self.fragments) <= 3:
            return 0.5, 0.02
        else:
            return 1.0, 0.2



if __name__ == '__main__':
    Gibson.main(__doc__)

# vim: tw=50
