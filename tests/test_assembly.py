from stepwise_mol_bio._assembly import *
from param_helpers import *

@parametrize_from_file(
        schema=[
            cast(parse_value=with_py.eval, expected=with_py.eval),
            defaults(parse_value=lambda x: x),
        ],
)
def test_parse_fragment_attrs(args, parse_value, expected):
    actual = parse_fragment_attrs(args, parse_value)
    for k, v in expected.items():
        assert actual[k] == v

@parametrize_from_file(
        schema=[
            cast(args=with_py.eval, expected=with_swmb.eval),
            with_swmb.error_or('expected'),
        ],
)
def test_parse_assemblies_from_docopt(args, expected, error):
    with error:
        assert parse_assemblies_from_docopt(args) == expected

@parametrize_from_file(
        schema=[
            cast(assemblies=with_swmb.eval, kwargs=with_swmb.eval),
            defaults(kwargs={}, warning=''),
        ],
)
def test_add_fragments_to_reaction(assemblies, kwargs, expected, warning, capsys):
    rxn = stepwise.MasterMix()
    rxn['water'].volume = 'to 5 µL'

    add_fragments_to_reaction(rxn, assemblies, **kwargs)
    cap = capsys.readouterr()

    assert str(rxn) == expected
    assert warning in cap.err
