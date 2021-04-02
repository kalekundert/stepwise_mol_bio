#!/usr/bin/env python3

"""
Protocols relating to molecular biology, e.g. PCR.
"""

__version__ = '1.14.0'

from ._utils import *
from ._presets import *

from pathlib import Path
from numbers import Real
from voluptuous import Any

class Plugin:
    protocol_dir = Path(__file__).parent
    config_defaults = protocol_dir / 'conf.toml'
