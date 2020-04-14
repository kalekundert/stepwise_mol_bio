#!/usr/bin/env python3

from appdirs import AppDirs
from inform import Error

app = AppDirs("stepwise_mol_bio")

class Main:

    @classmethod
    def main(cls, *args, **kwargs):
        try:
            from docopt import docopt
            args = docopt(*args, **kwargs)
            self = cls.from_docopt(args)
            print(self.protocol)
        except Error as err:
            err.report()

    @classmethod
    def from_docopt(cls, args):
        return cls()



class StepwiseMolBioError(Error):
    pass

class ConfigError(StepwiseMolBioError):
    pass

class UsageError(StepwiseMolBioError):
    pass

def hanging_indent(text, prefix):
    from textwrap import indent
    if isinstance(prefix, int):
        prefix = ' ' * prefix
    return indent(text, prefix)[len(prefix):]

