#!/usr/bin/env python3

from stepwise import repr_join
from stepwise_mol_bio import UsageError

def plan_centrifuge_step(params):
    time_formatters = {
            'time': '{}',
            'time_s': '{}s',
            'time_m': '{} min',
            'time_h': '{} hr',
    }
    speed_formatters = {
            'speed': '{}',
            'speed_g': lambda x: f'{int(x):,}g',
            'speed_rpm': lambda x: f'{int(x):,} rpm',
    }
    temp_formatters = {
            'temp': '{}',
            'temp_C': '{}°C',
    }

    def format_param(params, formatters, missing_ok=False):
        keys = set(params) & set(formatters)

        if not keys:
            if missing_ok:
                return ''
            else:
                raise UsageError(f"expected one of the following centrifugation parameters: {repr_join(formatters)}")
        if len(keys) > 1:
            raise UsageError(f"found multiple values for single centrifugation parameter: {repr_join(k for k in formatters if k in params)}")

        key = keys.pop()
        formatter = formatters[key]
        if isinstance(formatter, str):
            formatter = formatter.format

        return formatter(params[key])

    time = format_param(params, time_formatters)
    speed = format_param(params, speed_formatters)
    temp = format_param(params, temp_formatters, missing_ok=True)


    if temp:
        return f"Centrifuge at {speed} for {time} at {temp}."
    else:
        return f"Centrifuge at {speed} for {time}."

