#!/usr/bin/env python3

import parametrize_from_file
from stepwise_mol_bio.autoclave import *
from param_helpers import *

@parametrize_from_file(
        schema=Schema({
            'volume_mL': Coerce(int),
            'expected': Coerce(int),
        }),
)
def test_calc_sterilization_time(volume_mL, expected):
    app = Autoclave()
    app.volume_mL = volume_mL
    assert app.time_min == expected

