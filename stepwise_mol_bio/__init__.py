#!/usr/bin/env python3

"""
Protocols relating to molecular biology, e.g. PCR.
"""

__version__ = '0.0.0'

from pathlib import Path
from numbers import Real

class Plugin:
    protocol_dir = Path(__file__).parent
    config_defaults = protocol_dir / 'conf.toml'
    config_schema = {
            'pcr': {
                'primer_stock_uM': Real,
                'polymerases': {
                    str: {
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
                    }
                }
            }
    }
