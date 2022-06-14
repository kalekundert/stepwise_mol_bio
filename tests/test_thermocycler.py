from param_helpers import *
from stepwise_mol_bio.thermocycler import *

test_protocol()
test_cli()

@parametrize_from_file(
        schema=[
            cast(expected=with_py.eval),
            with_swmb.error_or('expected'),
        ],
)
def test_parse_thermocycler_steps(given, expected, error):
    with error:
        assert parse_thermocycler_steps(given) == expected

@parametrize_from_file(
        schema=[
            cast(steps=with_py.eval, expected=with_sw.eval),
            with_swmb.error_or('expected'),
        ],
)
def test_format_thermocycler_steps(steps, expected, error):
    with error:
        assert format_thermocycler_steps(steps) == expected

