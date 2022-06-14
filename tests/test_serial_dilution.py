import pytest
from stepwise_mol_bio.serial_dilution import *
from param_helpers import *
from pytest import approx

SET_CONC_METHODS = [
        lambda app, high, low, factor: app.set_conc_high_low(high, low),
        lambda app, high, low, factor: app.set_conc_high_factor(high, factor),
        lambda app, high, low, factor: app.set_conc_low_factor(low, factor),
]

@parametrize_from_file(
        schema=[
            cast(expected=with_py.eval),
            with_swmb.error_or('expected'),
        ],
)
def test_parse_high_low(high_str, low_str, expected, error):
    with error:
        high, low, unit = parse_high_low(high_str, low_str)
        assert high == expected['high']
        assert low == expected['low']
        assert unit == expected['unit']

@pytest.mark.parametrize('set_conc', SET_CONC_METHODS)
@pytest.mark.parametrize('include_zero', [True, False])
@parametrize_from_file(schema=with_math.eval)
def test_concentrations(set_conc, high, low, factor, include_zero, expected):
    sd = SerialDilution(1, len(expected))
    sd.include_zero = include_zero
    set_conc(sd, high, low, factor)

    if include_zero:
        expected = [*expected, 0]

    assert sd.conc_high == approx(high)
    assert sd.conc_low == approx(low)
    assert sd.factor == approx(factor)
    assert sd.inv_factor == approx(1 / factor)
    assert sd.concentrations == approx(expected)

@pytest.mark.parametrize('set_conc', SET_CONC_METHODS)
def test_check_num_dilutions(set_conc):
    sd = SerialDilution(1, 1)
    error = with_swmb.error({
        'type': 'UsageError',
        'message': [
            'must request at least 2 dilutions',
            'only 1 dilution requested',
        ],
    })

    with error:
        set_conc(sd, 16, 1, 4)

test_protocol()
test_cli()
