#!/usr/bin/env python3

"""
Protocols relating to molecular biology, e.g. PCR.
"""

__version__ = '1.23.1'

from ._utils import *
from ._assembly import Assembly
from .aliquot import Aliquot
from .anneal import Anneal
from .digest import RestrictionDigest
from .direct_dilution import DirectDilution
from .ethanol_precipitation import EthanolPrecipitation
from .gels.gel import Gel
from .gels.laser_scanner import LaserScanner
from .gels.stain import Stain
from .gels.transilluminator import Transilluminator
from .gibson import Gibson
from .golden_gate import GoldenGate
from .ivt import Ivt
from .ivtt import Ivtt
from .kld import Kld
from .ligate import Ligate
from .lyophilize import Lyophilize
from .miniprep import Miniprep
from .pcr import Pcr
from .serial_dilution import SerialDilution
from .sequence import Sequence
from .spin_cleanup import SpinCleanup

# Avoid circular imports
from .invpcr import InversePcr
from .qpcr import Qpcr
from .page_purify import PagePurify

import stepwise
from pathlib import Path

class Plugin:
    protocol_dir = Path(__file__).parent
    config_path = protocol_dir / 'conf.toml'
    priority = stepwise.Builtins.priority + 10

del stepwise
del Path
