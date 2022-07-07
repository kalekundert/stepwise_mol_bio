#!/usr/bin/env python3

from stepwise_mol_bio.centrifuge import plan_centrifuge_step
from param_helpers import *

@parametrize_from_file(schema=with_swmb.error_or('expected'))
def test_plan_centrifuge_step(params, expected, error):
    with error:
        assert plan_centrifuge_step(params) == expected
