#!/usr/bin/env python3

import stepwise
import autoprop
import byoc

from stepwise import pl, ul
from stepwise_mol_bio import Main, UsageError
from freezerbox import parse_temp_C, parse_time_s, format_time_s
from byoc import Key, DocoptConfig
from more_itertools import first

# I thought of a way for this code to accommodate arbitrary user-configured 
# thermocycler protocols in the context of PCR.  This issue is that the PCR 
# protocol needs to tweak certain parts of the thermocycler protocol, e.g. 
# setting the annealing temperature and extension time.  That would be hard to 
# do now, but easy to do if a name could be associated with each step.  Here's 
# specifically what I'm thinking:
#
# - Allow each step dictionary to have a name.
# - Replace the dictionaries with actual objects, than abstract the process of 
#   setting temperatures etc.  This is important because there's no point 
#   accessing individual steps if you can't easily mutate them.  The API can 
#   stay the same: just convert the dicts into these objects on intake.
# - Save a mapping between names and these objects.

def parse_thermocycler_steps(step_strs):
    # The command-line interface only supports "regular" incubation steps.
    return [
            parse_thermocycler_step(x)
            for x in step_strs
    ]

def parse_thermocycler_step(step_str):
    # Syntax to support repeat:
    #
    #   "35x 98/10 60/20 72/2m"
    #
    # - This would allow PCR thermocycler protocols to be specified on the 
    #   command-line, which could be very useful.
    # - Parsing: Look for \d+x, then interpret all space-separated fields to 
    #   follow as incubation steps.  Recursive repeating not allowed.

    fields = step_str.split('/')
    if len(fields) != 2:
        raise UsageError(f"expected 'temperature/time', not: {step_str!r}")

    temp_str, time_str = fields
    return {
            'temp_C': parse_temp_C(temp_str, default_unit='°C'),
            'time_s': parse_time_s(time_str, default_unit='s'),
    }

def format_thermocycler_steps(given, *, incubate_prefix=False):
    normalize_temp_time(given)

    match given:
        case list():
            return ul.from_iterable(
                    format_thermocycler_steps(
                        x,
                        incubate_prefix=incubate_prefix,
                    )
                    for x in given if x is not None
            )

        case {'temp': temp_str, 'time': time_str}:
            step = f'{temp_str} for {time_str}'
            return f'Incubate at {step}.' if incubate_prefix else step

        case {'hold_C': hold_C}:
            return f'Hold at {hold_C:g}°C'

        case {'repeat': n, 'steps': steps}:
            return pl(
                    f'Repeat {n}x:',
                    format_thermocycler_steps(steps),
                    br='\n',
            )
        # Melt curves...
        case 'fluorescence':
            return "Measure fluorescence"

        case err:
            raise UsageError("unexpected step in thermocycler protocol: {err!r}", err=err)

def normalize_temp_time(given):
    if not isinstance(given, dict):
        return

    time_formatters = {
            'time': lambda x: x,
            'time_s': lambda x: format_time_s(int(x)),
            'time_m': lambda x: format_time_s(int(60 * x)),
            'time_h': lambda x: format_time_m(int(60 * x)),
    }
    temp_formatters = {
            'temp': lambda x: x,
            'temp_C': lambda x: f'{x:g}°C'
    }

    def normalize(given, formatters):
        keys = set(given) & set(formatters)

        if not keys:
            return
        if len(keys) > 1:
            raise UsageError(f"found multiple values for single thermocycler parameter: {repr_join(k for k in formatters if k in params)}")

        key = keys.pop()
        normalized_key = first(formatters)
        given[normalized_key] = formatters[key](given.pop(key))
    
    normalize(given, time_formatters)
    normalize(given, temp_formatters)

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
        if len(self.steps) == 1:
            p += format_thermocycler_steps(self.steps[0], incubate_prefix=True)
        else:
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
