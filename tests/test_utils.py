import parametrize_from_file

from pytest import approx
from stepwise_mol_bio._utils import *
from param_helpers import *

@parametrize_from_file(
        schema=[
            with_swmb.error_or('expected'),
            with_py.eval,
        ]
)
def test_match_len(given, length, expected, error):
    with error:
        assert match_len(given, length) == expected

@parametrize_from_file(schema=with_py.eval)
def test_int_or_expr(given, expected):
    actual = int_or_expr(given)
    assert actual == expected
    assert isinstance(actual, int)

@parametrize_from_file(schema=with_py.eval)
def test_float_or_expr(given, expected):
    actual = float_or_expr(given)
    assert actual == approx(expected)
    assert isinstance(actual, float)

@parametrize_from_file(schema=with_py.eval)
def test_round_down_to_1_sig_fig(given, expected):
    assert round_down_to_1_sig_fig(given) == approx(expected)

@parametrize_from_file(schema=with_py.eval)
def test_round_up_to_1_sig_fig(given, expected):
    assert round_up_to_1_sig_fig(given) == approx(expected)
