#!/usr/bin/env python3

"""
Protocols relating to molecular biology, e.g. PCR.
"""

__version__ = '1.2.0'

from ._utils import *
from ._presets import *

from pathlib import Path
from numbers import Real
from voluptuous import Any

class Plugin:
    protocol_dir = Path(__file__).parent
    config_defaults = protocol_dir / 'conf.toml'
    config_schema = {
            'pcr': {
                'primer_stock_uM': Real,
                'default_preset': str,
                'presets': {
                    str: {
                        'brief': str,
                        'inherit': str,
                        'reagents': str,
                        'num_cycles': int,
                        'initial_denature_temp_C': Real,
                        'initial_denature_time_s': Real,
                        'denature_temp_C': Real,
                        'denature_time_s': Real,
                        'anneal_temp_C': Real,
                        'anneal_time_s': Real,
                        'extend_temp_C': Real,
                        'extend_time_s': Real,
                        'final_extend_temp_C': Real,
                        'final_extend_time_s': Real,
                        'melt_curve_low_temp_C': Real,
                        'melt_curve_high_temp_C': Real,
                        'melt_curve_temp_step_C': Real,
                        'melt_curve_time_step_s': Real,
                        'hold_temp_C': Real,
                        'two_step': bool,
                        'qpcr': bool,
                    },
                },
            },
            'gel': {
                'presets': {
                    str: {
                        'brief': str,
                        'inherit': str,
                        'title': str,
                        'gel_type': str,
                        'gel_percent': Any(str, Real),
                        'gel_additive': str,
                        'sample_mix': str,
                        'sample_name': str,
                        'sample_conc': Real,
                        'sample_volume_uL': Real,
                        'mix_volume_uL': Real,
                        'mix_extra_percent': Real,
                        'incubate_temp_C': Real,
                        'incubate_time_min': Real,
                        'load_volume_uL': Real,
                        'run_volts': Real,
                        'run_time_min': Real,
                        'stain': str,
                    },
                },
            },
            'laser': {
                'presets': {
                    str: {
                        'brief': str,
                        'inherit': str,
                        'laser': Any(str, Real),
                        'filter': str,
                    },
                },
            },
    }
