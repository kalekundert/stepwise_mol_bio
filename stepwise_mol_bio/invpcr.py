#!/usr/bin/env python3

"""\
Clone a plasmid by inverse PCR.

{pcr}

    -L --skip-kld
        Skip the KLD reaction.  This is equivalent to just using the `pcr` 
        command.

    -P --skip-pcr
        Skip the PCR reaction.  This is equivalent to using the `kld` command 
        with a transformation step afterwards.
"""

import docopt
import autoprop
from stepwise import Protocol
from stepwise_mol_bio import Main, pcr, kld

@autoprop
class InversePcr(Main):

    def __init__(self):
        self.pcr = pcr.Pcr()
        self.kld = kld.Kld('PCR product')
        self.skip_pcr = False
        self.skip_kld = False

    @classmethod
    def from_docopt(cls, args):
        inv = cls()
        inv.pcr = pcr.Pcr.from_docopt(args)
        inv.skip_pcr = args['--skip-pcr']
        inv.skip_kld = args['--skip-kld']
        return inv

    def get_protocol(self):
        self.kld.num_reactions = self.pcr.num_reactions
        p = Protocol()
        if not self.skip_pcr:
            p += self.pcr.protocol
        if not self.skip_kld:
            p += self.kld.protocol
            p += "Transform 2 µL."
        return p

def copy_pcr_usage():
    i = pcr.__doc__.find('Usage:')
    return __doc__.format(pcr=pcr.__doc__[i:].replace('pcr', 'invpcr', 1).strip())

if __name__ == '__main__':
    usage = copy_pcr_usage()
    InversePcr.main(usage)




