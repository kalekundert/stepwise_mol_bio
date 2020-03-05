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
from inform import plural

args = docopt.docopt(__doc__)
protocol = dirty_water.Protocol()

kld = stepwise.Reaction('''\
Reagent               Stock        Volume  Master Mix
================  =========   ===========  ==========
water                         to 10.00 μL         yes
T4 ligase buffer        10x       1.00 μL         yes
T4 PNK              10 U/μL       0.25 μL         yes
T4 DNA ligase      400 U/μL       0.25 μL         yes
DpnI                20 U/μL       0.25 μL         yes
PCR product        50 ng/μL       1.50 μL
''')

kld.num_reactions = eval(args['<num_reactions>'])
kld.extra_fraction = 15

protocol += """\
Run {plural(kld.num_reactions):# ligation reaction/s}:

{kld}

- Incubate at room temperature for 1h.
"""

print(protocol)
