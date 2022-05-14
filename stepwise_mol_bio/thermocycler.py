#!/usr/bin/env python3

import stepwise
import autoprop
import byoc

from stepwise import pl, ul
from stepwise_mol_bio import Main, UsageError
from freezerbox import parse_temp_C, parse_time_s, format_time_s
from byoc import Key, DocoptConfig

def parse_thermocycler_steps(step_strs):
    # The command-line interface only supports "regular" incubation steps.
    return [
            parse_thermocycler_step(x)
            for x in step_strs
    ]

def parse_thermocycler_step(step_str):
    fields = step_str.split('/')
    if len(fields) != 2:
        raise UsageError(f"expected 'temperature/time', not: {step_str!r}")

    temp_str, time_str = fields
    return {
            'temp_C': parse_temp_C(temp_str, default_unit='°C'),
            'time_s': parse_time_s(time_str, default_unit='s'),
    }

def format_thermocycler_steps(given):
    match given:
        case list():
            return ul(*map(format_thermocycler_steps, given))

        case {'temp_C': temp_C, 'time_m': time_m}:
            return f'{temp_C:g}°C for {format_time_s(int(60 * time_m))}'

        case {'temp_C': temp_C, 'time_s': time_s}:
            return f'{temp_C:g}°C for {format_time_s(int(time_s))}'

        case {'temp_C': temp_C, 'time': time_str}:
            return f'{temp_C:g}°C for {time_str}'

        case {'repeat': n, 'steps': steps}:
            return pl(
                    f'Repeat {n}x:',
                    format_thermocycler_steps(steps),
                    br='\n',
            )
        case err:
            raise UsageError(f"unexpected step in thermocycler protocol: {err!r}")

@autoprop
class Thermocycler(Main):
    """
Format a thermocycler protocol.

Usage:
    thermocycler <steps>...

Arguments:
    <steps>
        Each step specifies an incubation temperature and time, in that order, 
        separated by a slash.  Units are optional.  By default, temperatures 
        are in °C and times are in seconds.  For example, here are several 
        equivalent ways to specify an incubation at 72°C for 90 seconds:

            72/90
            72/90s
            72/1m30

        Note that the command-line interface does not currently support steps, 
        like loops, that are more complicated than simple incubations.
"""

    def __init__(self, steps):
        self.steps = steps

    def get_protocol(self):
        p = stepwise.Protocol()
        p += pl(
                'Run the following thermocycler protocol:',
                format_thermocycler_steps(self.steps),
        )
        return p

    __config__ = [DocoptConfig]
    steps = byoc.param(
            Key(DocoptConfig, '<steps>', cast=parse_thermocycler_steps),
    )

if __name__ == '__main__':
    Thermocycler.entry_point()
