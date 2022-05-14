#!/usr/bin/env python3

"""
Protocols relating to molecular biology, e.g. PCR.
"""

__version__ = '1.24.0'

from ._utils import *
from ._assembly import Assembly
from .aliquot import Aliquot
from .anneal import Anneal
from .autoclave import Autoclave
from .digest import RestrictionDigest
from .direct_dilution import DirectDilution
from .ethanol_precipitation import EthanolPrecipitation
from .gels.gel import Gel, Ladder
from .gels.laser_scanner import LaserScanner
from .gels.stain import Stain
from .gels.transilluminator import Transilluminator
from .gibson import Gibson
from .golden_gate import GoldenGate
from .grow import Grow
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
from .thermocycler import Thermocycler
from .transform import Transform
from .trizol import Trizol

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
