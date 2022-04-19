#!/usr/bin/env python3

import parametrize_from_file
from stepwise_mol_bio.sequence import *
from param_helpers import *

@parametrize_from_file
def test_parse_reactions_docopt(given, expected):
    assert parse_reactions_docopt(given) == expected

@parametrize_from_file(
        schema=Schema({
            'given': with_py.eval,
            'expected': with_py.eval,
        }),
)
def test_merge_reactions(given, expected):
    assert merge_reactions(given) == expected

