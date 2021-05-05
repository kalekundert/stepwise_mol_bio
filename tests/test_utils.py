#!/usr/bin/env python3

import parametrize_from_file
from stepwise_mol_bio._utils import *
from schema_helpers import *

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

