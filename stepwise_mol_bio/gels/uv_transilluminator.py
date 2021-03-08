#!/usr/bin/env python3

import stepwise, appcli, autoprop
from stepwise_mol_bio import Main

@autoprop
class UvTransilluminator(Main):
    """\
Image a gel using a UV transilluminator.

Usage:
    uv_transilluminator [<wavelength>]
"""
    __config__ = [
            appcli.DocoptConfig(),
    ]

    wavelength_nm = appcli.param('<wavelength>', default=300)

    def get_protocol(self):
        p = stepwise.Protocol()
        p += f"Image with a {self.wavelength_nm} nm UV transilluminator."
        return p

if __name__ == '__main__':
    UvTransilluminator.main()
