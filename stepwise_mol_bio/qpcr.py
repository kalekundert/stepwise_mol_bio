#!/usr/bin/env python3

from appcli import Key, DocoptConfig
from stepwise import UsageError, StepwiseConfig, PresetConfig
from pcr import Pcr

class Qpcr(Pcr):
    """\
Amplify a DNA template using polymerase chain reaction (PCR).

Note that this protocol is simply an alias for the PCR protocol, with a 
different default preset.

Usage:
    {pcr}
"""
    __config__ = [
            DocoptConfig(usage_getter=lambda self: self.format_usage()),
            PresetConfig(),
            StepwiseConfig('molbio.qpcr'),
            StepwiseConfig('molbio.pcr'),
    ]
    def format_usage(self):
        return self.__doc__.format(
                pcr=super().__doc__
                        .split('Usage:', 1)[1]
                        .replace('pcr', 'qpcr', 1)
                        .strip(),
        )

if __name__ == '__main__':
    Qpcr.main()




