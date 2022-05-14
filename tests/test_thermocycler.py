from param_helpers import *
from stepwise_mol_bio.thermocycler import *

test_protocol()
test_cli()

@parametrize_from_file(
        schema=Schema({
            'given': [str],
            **with_swmb.error_or({
                'expected': with_py.eval,
            }),
        }),
)
def test_parse_thermocycler_steps(given, expected, error):
    with error:
        assert parse_thermocycler_steps(given) == expected

@parametrize_from_file(
        schema=Schema({
            'steps': with_py.eval,
            **with_swmb.error_or({
                'expected': with_sw.eval,
            }),
        }),
)
def test_format_thermocycler_steps(steps, expected, error):
    with error:
        assert format_thermocycler_steps(steps) == expected

