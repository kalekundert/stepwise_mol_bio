#!/usr/bin/env python3
# vim: tw=50

import stepwise, appcli, autoprop
from appcli import Key, DocoptConfig
from coomassie import Coomassie

@autoprop
class SimplyBlue(Coomassie):
    """\
Stain protein gels using SimplyBlue SafeStain (Invitrogen LC6060).

SimplyBlue SafeStain is a ready-to-use, proprietary Coomassie G-250 stain that 
is specially-formulated for rapid, sensitive detection, and safe, non-hazardous 
disposal. This stain does not require the use of methanol or acetic acid, thus 
eliminating the need to dispose of hazardous waste.  Proteins stained using 
SimplyBlue SafeStain are compatible with mass spectrometry (MS) analysis.

Usage:
    simplyblue [-f | -o | -q] [-r] [-a] [-g <mm>]

Options:
    -f --fast
        Use a microwave to speed up the staining and destaining steps.

    -o --optimized
        Use an optimized protocol that takes slightly longer and requires 
        concentrated NaCl, but significantly improves sensitivity.

    -q --quiet
        Don't display any details for how to stain the gel.  Use this if you've 
        memorized the protocol and don't want to waste space.

    -r --fluorescent
        Image the gel using a near-IR fluorescence laser scanner, rather than
        a colorimetric gel imager.

    -a --attach-pdf
        Attach the manufacturer's protocol to this one, so that both will be 
        printed out by `stepwise go`.

    -g --gel-thickness <mm>  [default: 1.0]
        The thickness of the gel in millimeters.  Thicker gels sometimes 
        require more stain/destain or longer incubate times.
"""

    stain_type = appcli.param(
            Key(DocoptConfig, '--fast', cast=lambda x: 'microwave'),
            Key(DocoptConfig, '--quiet', cast=lambda x: 'quiet'),
            Key(DocoptConfig, '--optimized', cast=lambda x: 'optimal'),
            default='basic',
    )
    gel_thickness_mm = appcli.param(
            '--gel-thickness',
            cast=float,
            default=1.0,
    )
    attach_pdf = appcli.param(
            '--attach-pdf',
            default=False,
    )

    def __init__(self):
        super().__init__()
        self.default_imaging_protocol = 'fluorescent'
        self.gel_thickness_mm = 1.0
        self.attach_pdf = False

    def get_staining_protocols(self):
        return {
                'basic': self.get_basic_protocol,
                'optimal': self.get_optimal_protocol,
                'microwave': self.get_microwave_protocol,
                'quiet': self.get_quiet_protocol,
        }

    def get_basic_protocol(self):
        return self.get_standard_protocol(False)

    def get_optimal_protocol(self):
        return self.get_standard_protocol(True)

    def get_standard_protocol(self, max_sensitivity):
        p = stepwise.Protocol()
        k = self.gel_thickness_mm

        p += """\
Stain gel with SimplyBlue SafeStain:

- Rinse gel 3x for 5 min with 100 mL water [1].
- Add enough stain to cover the gel (≈20 mL).
- Incubate at 25°C for 1h with gentle shaking [2].
"""
        if max_sensitivity:
            p.steps[0] += f"""\
- Wash the gel with {100*k:.0f} mL water for 1h [3].
- Add {20*k:.0f} mL 5M NaCl to the above wash.
- Continue washing for 2-16h.
"""
        else: 
            p.steps[0] += """\
- Wash the gel with 100 mL water for 1–3h [3].
- To obtain the clearest background, wash again 
  with 100 mL water for 1h [4].
"""

        p.footnotes[1] = """\
Rinsing removes SDS and buffer salts, which 
interfere with binding of the dye to the protein.
"""
        p.footnotes[2] = """\
The gel can be stained for up to 3 hours, but 
after 3 hours, sensitivity decreases.  If you need 
to leave the gel overnight in the stain, add 2 mL 
of 5M NaCl in water for every 20 mL of stain. This 
procedure will not affect sensitivity.
"""
        p.footnotes[3] = """\
The gel can be left in the water for several days 
without loss of sensitivity. There is a small 
amount of dye in the water that is in equilibrium 
with the dye bound to the protein, so proteins 
remain blue.
"""
        p.footnotes[4] = """\
Sensitivity decreases if the gel is stored in 
water for more than 1 day. The decrease in the 
amount of free dye in the water favors 
dissociation of the dye from the protein.  If you 
need to store the gel in water for a few days, add 
20 mL of 5M NaCl.
"""
        p.prune_footnotes()

        return p

    def get_microwave_protocol(self):
        p = stepwise.Protocol()
        t = 5 if self.gel_thickness_mm == 1.0 else 10

        p += f"""\
Stain gel with SimplyBlue SafeStain:

- Repeat 3 times:
    - Add 100 mL water.
    - Microwave until almost boiling (1 min).
    - Shake gently for 1 min, then discard water.
- Add enough stain to cover the gel.
- Microwave until almost boiling (45-60s).
- Shake gently for {t} min [1].
- Wash the gel with 100 mL water for 10 min [2].
- Add 20 mL 5M NaCl and wash for 5 min [3].

"""
        p.footnotes[1] = """\
Detection limit: 20 ng BSA
"""
        p.footnotes[2] = """\
Detection limit: 10 ng BSA
"""
        p.footnotes[3] = """\
Detection limit:  5 ng BSA
"""
        return p

    def get_quiet_protocol(self):
        p = stepwise.Protocol()
        p += """\
Stain gel with SimplyBlue SafeStain.
"""
        if self.attach_pdf:
            from pathlib import Path
            p.attachments = [
                    Path(__file__).parent / 'invitrogen_simplyblue_safestain.pdf'
            ]

        return p

if __name__ == '__main__':
    SimplyBlue.main()

