#!/usr/bin/env python3

import sys
import appcli
from appdirs import AppDirs
from inform import Error, format_range
from pathlib import Path

app_dirs = AppDirs("stepwise_mol_bio")

class Main(appcli.App):
    usage_io = sys.stderr

    @classmethod
    def main(cls):
        self = cls.from_params()
        self.load()
        
        try:
            self.protocol.print()
        except Error as err:
            err.report()



class StepwiseMolBioError(Error):
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

def merge_dicts(dicts):
    result = {}
    for dict in reversed(list(dicts)):
        result.update(dict)
    return result

def comma_list(x):
    return [x.strip() for x in x.split(',')]

def comma_set(x):
    return {x.strip() for x in x.split(',')}

def require_reagent(rxn, reagent):
    if reagent not in rxn:
        raise UsageError(f"reagent table missing {reagent!r}")


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

