#!/usr/bin/env python3

import stepwise, appcli, autoprop
from appcli import Key, DocoptConfig
from stepwise_mol_bio.gels._stain import Stain

@autoprop
class Biotium(Stain):
    __config__ = [
            DocoptConfig(),
    ]

    product = None
    uv_wavelength = None
    attachment = None
    image_type = appcli.param(
            Key(DocoptConfig, '--no-imaging', cast=lambda x: None),
            default='uv',
    )
    attach_pdf = appcli.param(
            Key(DocoptConfig, '--attach-pdf'),
            default=False,
    )

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
