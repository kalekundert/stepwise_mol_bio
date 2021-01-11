#!/usr/bin/env python3

import stepwise, appcli, autoprop
from appcli import Key, DocoptConfig
from stepwise_mol_bio import Main, pcr, kld
from pcr import Pcr

@autoprop
class InversePcr(Pcr):
    """\
Clone a plasmid by inverse PCR.

Usage:
    {pcr}

    -L --skip-kld
        Skip the KLD reaction.  This is equivalent to just using the `pcr` 
        command.

    -P --skip-pcr
        Skip the PCR reaction.  This is equivalent to using the `kld` command 
        with a transformation step afterwards.
"""
    skip_pcr = appcli.param(
            Key(DocoptConfig, '--skip-pcr'),
            default=False,
    )
    skip_kld = appcli.param(
            Key(DocoptConfig, '--skip-kld'),
            default=False,
    )

    def __bareinit__(self):
        self.kld = kld.Kld('PCR product')

    def format_usage(self):
        return self.__doc__.format(
                pcr=super().__doc__
                        .split('Usage:', 1)[1]
                        .replace('pcr', 'invpcr', 1)
                        .strip(),
        )

    def get_protocol(self):
        self.kld.num_reactions = self.num_reactions
        p = stepwise.Protocol()
        if not self.skip_pcr:
            p += super().protocol
        if not self.skip_kld:
            p += self.kld.protocol
            p += "Transform 2 µL."
        return p

if __name__ == '__main__':
    InversePcr.main()




