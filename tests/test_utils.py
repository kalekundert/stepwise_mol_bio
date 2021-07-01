#!/usr/bin/env python3

import parametrize_from_file

from pytest import approx
from stepwise_mol_bio._utils import *
from param_helpers import *

@parametrize_from_file(
        schema=Schema({
            'given': eval,
            'length': eval,
            **error_or({
                'expected': eval,
            }),
        }),
)
def test_match_len(given, length, expected, error):
    with error:
        assert match_len(given, length) == expected

@parametrize_from_file(schema=Schema({str: eval}))
def test_int_or_expr(given, expected):
    actual = int_or_expr(given)
    assert actual == expected
    assert isinstance(actual, int)

@parametrize_from_file(schema=Schema({str: eval}))
def test_float_or_expr(given, expected):
    actual = float_or_expr(given)
    assert actual == approx(expected)
    assert isinstance(actual, float)
