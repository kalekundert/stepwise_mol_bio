#!/usr/bin/env python3

import sys
import appcli
import tidyexc
from appdirs import AppDirs
from inform import format_range
from more_itertools import all_equal
from pathlib import Path

app_dirs = AppDirs("stepwise_mol_bio")

class Main(appcli.App):
    usage_io = sys.stderr

    @classmethod
    def main(cls):
        self = cls.from_params()
        self.load(appcli.DocoptConfig)
        
        try:
            self.protocol.print()
        except StepwiseMolBioError as err:
            print(err, file=sys.stderr)



class StepwiseMolBioError(tidyexc.Error):
    pass

class ConfigError(StepwiseMolBioError):
    pass

class UsageError(StepwiseMolBioError):
    pass

def try_except(expr, exc, failure, success=None):
    try:
        x = expr()
    except exc:
        return failure()
    else:
        return success() if success else x

def hanging_indent(text, prefix):
    from textwrap import indent
    if isinstance(prefix, int):
        prefix = ' ' * prefix
    return indent(text, prefix)[len(prefix):]

def join_lists(lists):
    return sum(lists, [])

def merge_dicts(dicts):
    result = {}
    for dict in reversed(list(dicts)):
        result.update(dict)
    return result

def comma_list(x):
    return [x.strip() for x in x.split(',')]

def comma_set(x):
    return {x.strip() for x in x.split(',')}

def int_or_expr(x):
    return type_or_expr(int, x)

def float_or_expr(x):
    return type_or_expr(float, x)

def type_or_expr(type, x):
    if isinstance(x, type):
        return x
    else:
        return type(eval(x))

def require_reagent(rxn, reagent):
    if reagent not in rxn:
        raise UsageError(f"reagent table missing {reagent!r}")

def merge_names(names):
    names = list(names)
    if all_equal(names):
        return names[0]
    else:
        return ','.join(names)




def format_sec(x):
    if x < 60:
        return f'{x}s'

    min = x // 60
    sec = x % 60

    return f'{min}m{f"{sec:02}" if sec else ""}'

def format_min(x):
    if x < 60:
        return f'{x}m'

    hr = x // 60
    min = x % 60

    return f'{hr}h{f"{sec:02}" if min else ""}'

