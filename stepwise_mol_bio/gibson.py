#!/usr/bin/env python3

import stepwise
import appcli
import autoprop
from inform import plural
from _assembly import Assembly

@autoprop
class Gibson(Assembly):
    """\
Perform a Gibson assembly reaction, using the NEB master mix (E2611).

Usage:
    gibson_assembly <backbone> <inserts>... [options]

Arguments:
{FRAGMENT_DOC}

Options:
{OPTION_DOC}
"""

    def get_reaction(self):
        rxn = stepwise.MasterMix()
        rxn.volume = '20 µL'

        rxn['Gibson master mix (NEB E2611)'].volume = rxn.volume / 2
        rxn['Gibson master mix (NEB E2611)'].stock_conc = "2x"
        rxn['Gibson master mix (NEB E2611)'].master_mix = True
        rxn['Gibson master mix (NEB E2611)'].order = 2

        return self._add_fragments_to_reaction(rxn)

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
        p += f"""\
Transform 2 µL.
"""
        p.footnotes[1] = """\
https://preview.tinyurl.com/ychbvkra
"""
        return p

    def get_usage(self):
        return format_docstring(self, self.__doc__)

if __name__ == '__main__':
    Gibson.main()

# vim: tw=50
