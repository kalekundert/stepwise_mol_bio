#!/usr/bin/env python3

import stepwise
import stepwise_mol_bio.dnase as dnase
from param_helpers import *

@parametrize_from_file(
        schema=with_swmb.error_or('expected'),
)
def test_plan_dnase_reactions(group, expected, error):
    group = eval_sample_group(group)

    with error:
        rxns = dnase.plan_dnase_reactions(group)
        assert match_protocol(rxns.protocol, expected)

@parametrize_from_file
def test_plan_incubation_steps(group, expected):
    group = eval_sample_group(group)

    p = stepwise.Protocol()
    p += dnase.plan_incubation_steps(group)
    assert match_protocol(p, expected)

@parametrize_from_file
def test_plan_precipitation_steps(group, expected):
    group = eval_sample_group(group)

    p = stepwise.Protocol()
    p += dnase.plan_precipitation_steps(group)
    assert match_protocol(p, expected)

test_cli()
