#!/usr/bin/env python3

from stepwise_mol_bio._assembly import *
from schema_helpers import *

@parametrize_from_file(
        schema=Schema({
            'args': [str],
            Optional('parse_value', default='lambda x: x'): eval,
            'expected': {str: eval},
        }),
)
def test_parse_fragment_attrs(args, parse_value, expected):
    actual = parse_fragment_attrs(args, parse_value)
    for k, v in expected.items():
        assert actual[k] == v

@parametrize_from_file(
        schema=Schema({
            'args': {str: eval_with()},
            **error_or({
                'expected': eval_swmb,
            })
        }),
)
def test_parse_assemblies_from_docopt(args, expected, error):
    with error:
        assert parse_assemblies_from_docopt(args) == expected

@parametrize_from_file(
        schema=Schema({
            'assemblies': [[eval_swmb]],
            Optional('kwargs', default={}): {str: eval_swmb},
            'expected': str,
            Optional('warning', default=''): str,
        }),
)
def test_add_fragments_to_reaction(assemblies, kwargs, expected, warning, capsys):
    rxn = stepwise.MasterMix()
    rxn.volume = 5, 'µL'

    add_fragments_to_reaction(rxn, assemblies, **kwargs)
    cap = capsys.readouterr()

    assert str(rxn) == expected
    assert warning in cap.err
