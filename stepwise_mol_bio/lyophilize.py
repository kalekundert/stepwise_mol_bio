#!/usr/bin/env python3

import stepwise, byoc, autoprop
from stepwise import Quantity
from stepwise_mol_bio import Cleanup, UsageError
from freezerbox import MakerConfig, group_by_identity, parse_volume, parse_conc
from byoc import DocoptConfig, Key
from inform import plural

@autoprop
class Lyophilize(Cleanup):
    """\
Concentrate samples by lyophilization.

Usage:
    lyophilize [-v <µL> | -c <conc>]

Options:
    -v --volume <µL>
        The volume to bring the samples to after lyophilization.  If no unit is 
        specified, µL is assumed.

    -c --conc <conc>
        The concentration to bring the samples to after lyophilization.  
        Include a unit.  It is assumed that the experimenter will measure the 
        concentration of the sample before lyophilization and calculate the 
        correct volume of buffer to add.

Database:
    The lyophilize protocol can appear in the "Cleanup" column of a FreezerBox 
    database:
        
        lyophilize [volume=<volume>] [conc=<conc>]

    volume=<volume>
        See --volume.  You must specify a unit.

    conc=<conc>
        See --conc.  You must specify a unit.
"""
    __config__ = [
            DocoptConfig,
            MakerConfig,
    ]
    volume = byoc.param(
            Key(DocoptConfig, '--volume', cast=lambda x: parse_volume(x, default_unit='µL')),
            Key(MakerConfig, 'volume', cast=parse_volume),
            default=None,
    )
    conc = byoc.param(
            Key(DocoptConfig, '--conc'),
            Key(MakerConfig, 'conc', cast=parse_conc),
            default=None,
    )

    group_by = {
        'volume': group_by_identity,
        'conc': group_by_identity,
    }

    def __init__(self, *, volume=None, conc=None):
        if volume: self.volume = volume
        if conc: self.conc = conc

    def get_protocol(self):
        phrases = []
        n = plural(self.product_tags)

        if self.show_product_tags:
            phrases.append(f"Concentrate the following {n:sample/s}")
        elif self.product_tags:
            phrases.append(f"Concentrate the {n:sample/s}")
        else:
            phrases.append("Concentrate the sample(s)")

        if self.volume and self.conc:
            err = UsageError(volume=self.volume, conc=self.conc)
            err.brief = "cannot specify volume and concentration"
            err.info += "volume: {volume}"
            err.info += "conc: {conc}"
            raise err
        elif self.volume:
            phrases.append(f"to {self.volume}")
        elif self.conc:
            phrases.append(f"to {self.conc}")
        else:
            pass

        if self.show_product_tags:
            phrases.append(f"by lyophilization: {', '.join(map(str, self.product_tags))}")
        else:
            phrases.append("by lyophilization.")

        return stepwise.Protocol(
                steps=[' '.join(phrases)],
        )

    def get_product_conc(self):
        return self.conc

    def get_product_volume(self):
        return self.volume

if __name__ == '__main__':
    Lyophilize.main()
