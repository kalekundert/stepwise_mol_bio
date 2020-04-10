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
from inform import Error, did_you_mean

@autoprop
class Coomassie:

    def __init__(self):
        self.stain_type = 'basic'
        self.image_type = None
        self.default_image_type = 'colorimetric'

    @classmethod
    def main(cls, *args, **kwargs):
        try:
            from docopt import docopt
            args = docopt(*args, **kwargs)
            coom = cls.from_docopt(args)
            print(coom.protocol)
        except Error as err:
            err.report()

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

    def get_protocol(self):
        return self.staining_protocol + self.imaging_protocol

    def get_staining_protocol(self):
        try:
            return self.staining_protocols[self.stain_type]()
        except KeyError:
            raise Error("unknown staining protocol {self.stain_type!r}, did you mean {did_you_mean(self.stain_type, self.staining_protocols)}")

    def get_staining_protocols(self):
        return {
                'basic': self.get_basic_staining,
                'microwave': self.get_microwave_staining,
                'quiet': self.get_quiet_staining,
        }

    def get_imaging_protocol(self):
        key = self.image_type or self.default_image_type
        try:
            return self.imaging_protocols[key]()
        except KeyError:
            raise Error("unknown imaging protocol {key!r}, did you mean {did_you_mean(key, self.imaging_protocols.keys())}")

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
        p = stepwise.Protocol()
        p += """\
Image with a laser scanner [1]:

laser: 658 nm
filter: 710BP40
"""
        p.footnotes[1] = """\
Coomassie seems to quench fluorophores like 
FITC/GFP.  I don't know exactly why this is.
"""
        return p

if __name__ == '__main__':
    Coomassie.main(__doc__)

