from stepwise_mol_bio import Ladder
from stepwise_mol_bio.gels.gel import *
from param_helpers import *

test_protocol()
test_cli()

@parametrize_from_file(
        schema=Schema({
            'name': str,
            'expected': Coerce(int),
        }),
)
def test_parse_num_samples(name, expected):
    assert parse_num_samples(name) == expected

@parametrize_from_file(
        schema=Schema({
            'ladder_str': str,
            'name': str,
            'volume_uL': with_py.eval,
        }),
)
def test_ladder_from_string(ladder_str, name, volume_uL):
    ladder = Ladder.from_string(ladder_str)
    assert ladder.name == name
    assert ladder.volume_uL == volume_uL

