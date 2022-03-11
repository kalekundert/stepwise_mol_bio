#!/usr/bin/env python3

import stepwise
import autoprop
import appcli

from stepwise_mol_bio import Cleanup
from freezerbox import MakerConfig, parse_volume_mL, parse_time_m
from appcli import Key, Method, DocoptConfig
from more_itertools import pairwise
from math import ceil

@autoprop
class Autoclave(Cleanup):
    """
    Sterilize media/buffers by autoclave.

    Usage:
        autoclave <volume_mL>
        autoclave -t <min>

    Arguments:
        <volume_mL>
            The greatest volume of liquid in any of the bottles being 
            autoclaved.

    Options:
        -t --time <min>
            The specific amount of time to autoclave for.

    Database:
        The autoclave protocol can appear in the "Cleanup" column of a 
        FreezerBox database:

            autoclave [<volume>]
            autoclave [time=<min>]

        <volume>
            See <volume_mL>.  You must include a unit.

        time=<min>
            See --time.  You must specify a unit.

    References:
       https://tinyurl.com/2bsfur3b
    """

    __config__ = [
            DocoptConfig,
            MakerConfig,
    ]

    def __init__(self, volume_mL=None, *, time_min=None):
        if volume_mL:
            self.volume_mL = volume_mL
        if time_min:
            self.time_min = time_min

    def get_protocol(self):
        p = stepwise.Protocol()
        p += f"Autoclave at 121°C for {self.time_min} min."
        return p

    def _calc_sterilization_time_min(self):
        sterilization_mL_min = [
                (75, 25),
                (250, 30),
                (500, 40),
                (1000, 45),
                (1500, 50),
                (2000, 55),
        ]
        for volume_mL, time_min in sterilization_mL_min:
            if self.volume_mL <= volume_mL:
                return time_min

        return 35 + 10 * int(ceil(self.volume_mL / 1000))

    volume_mL = appcli.param(
            Key(DocoptConfig, '<volume_mL>', cast=int),
            Key(MakerConfig, 1, cast=[parse_volume_mL, int]),
    )
    time_min = appcli.param(
            Key(DocoptConfig, '--time', cast=int),
            Key(MakerConfig, 'time', cast=[parse_time_m, int]),
            Method(_calc_sterilization_time_min),
    )
    merge_by = {
            'time_min': max,
    }

if __name__ == '__main__':
    Autoclave.main()
