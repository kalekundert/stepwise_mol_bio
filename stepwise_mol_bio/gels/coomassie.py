#!/usr/bin/env python3

"""\
Stain protein gels using Coomassie R-250.

Usage:
    coomassie [-f | -q] [-r]

Options:
    -f --fast
        Use a microwave to speed up the staining and destaining steps.
    
    -q --quiet
        Don't display any details for how to stain the gel.  Use this if you've 
        memorized the protocol and don't want to waste space.

    -r --fluorescent
        Image the gel using a near-IR fluorescence laser scanner, rather than
        a colorimetric gel imager.
"""

import stepwise
import autoprop
from _stain import Stain
from laser_scanner import LaserScanner

@autoprop
class Coomassie(Stain):

    def __init__(self):
        self.stain_type = 'basic'
        self.image_type = None
        self.default_image_type = 'colorimetric'

    @classmethod
    def from_docopt(cls, args):
        self = cls()
        if args['--fast']:
            self.stain_type = 'microwave'
        if args['--quiet']:
            self.stain_type = 'quiet'
        if args['--fluorescent']:
            self.image_type = 'fluorescent'
        return self

    def get_staining_protocols(self):
        return {
                'basic': self.get_basic_staining,
                'microwave': self.get_microwave_staining,
                'quiet': self.get_quiet_staining,
        }

    def get_imaging_protocols(self):
        return {
                'colorimetric': self.get_colorimetric_imaging,
                'fluorescent': self.get_fluorescent_imaging,
        }

    def get_basic_staining(self):
        p = stepwise.Protocol()
        p += """\
Stain gel with Coomassie:

- Submerge gel in fresh stain.
- Incubate 1-16h with gentle shaking.
- Repeat until the background is clear:
  - Submerge gel in fresh destain
  - Gently shake for 30 min.
"""
        return p

    def get_microwave_staining(self):
        p = stepwise.Protocol()
        p += """\
Stain gel with Coomassie:

- Submerge the gel in fresh stain.
- Microwave on high for 30 sec [1].
- Gently shake for 5–10 min.
- Rinse twice with water.

- Repeat until the background is clear:
  - Submerge the gel in fresh destain.
  - Microwave on high for 30 sec.
  - Place a wadded-up kimwipe in the destain.
  - Gently shake for 10 min.
"""
        p.footnotes[1] = """\
Coomassie stain contains methanol, so avoid 
breathing fumes when using the microwave.
"""
        return p

    def get_quiet_staining(self):
        p = stepwise.Protocol()
        p += """\
Stain gel with Coomassie.
"""
        return p

    def get_colorimetric_imaging(self):
        return stepwise.Protocol()

    def get_fluorescent_imaging(self):
        scan = LaserScanner.from_params(658, '710BP40')
        return scan.protocol

if __name__ == '__main__':
    Coomassie.main(__doc__)

