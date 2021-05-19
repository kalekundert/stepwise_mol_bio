#!/usr/bin/env python3

"""
Protocols relating to molecular biology, e.g. PCR.
"""

__version__ = '1.19.0'

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
from .gels.uv_transilluminator import UvTransilluminator
from .gibson import Gibson
from .golden_gate import GoldenGate
from .ivt import Ivt
from .ivtt import InVitroTranslation
from .kld import Kld
from .ligate import Ligate
from .pcr import Pcr
from .serial_dilution import SerialDilution
from .spin_cleanup import SpinCleanup

# Avoid circular imports
from .invpcr import InversePcr
from .qpcr import Qpcr

from pathlib import Path

class Plugin:
    protocol_dir = Path(__file__).parent
    config_defaults = protocol_dir / 'conf.toml'

del Path
