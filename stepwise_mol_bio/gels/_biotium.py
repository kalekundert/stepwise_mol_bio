#!/usr/bin/env python3

"""\
Stain nucleic acid gels using GelGreen.

Usage:
    gelred
"""

import stepwise
import autoprop
from _stain import Stain

@autoprop
class Biotium(Stain):
    product = None
    uv_wavelength = None
    attachment = None
    default_image_type = 'uv'

    def __init__(self):
        super().__init__()
        self.attach_pdf = False

    @classmethod
    def from_docopt(cls, args):
        self = super().from_docopt(args)
        self.attach_pdf = args['--attach-pdf']
        return self

    def get_staining_protocol(self):
        p = stepwise.Protocol()
        p += f"""\
Stain gel with {self.product}:

- Submerge gel in ≈50 mL 3x {self.product} [1].
- Shake gently for 30 min.
"""
        p.footnotes[1] = """\
Including 100 mM NaCl in the staining solution 
enhances sensitivity, but may promote dye 
precipitation if the gel stain is reused.

Staining solution can be reused at least 2-3x.
Store staining solution at room temperature
protected from light.
"""
        if self.attach_pdf:
            from pathlib import Path
            p.attachments = [
                    Path(__file__).parent / self.attachment,
            ]

        return p

    def get_imaging_protocols(self):
        return {
                'uv': self.get_uv_imaging,
        }

    def get_uv_imaging(self):
        p = stepwise.Protocol()
        p += f"""\
Image with a {self.uv_wavelength} nm UV transilluminator.
"""
        return p
