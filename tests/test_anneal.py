#!/usr/bin/env python3

import parametrize_from_file
from stepwise_mol_bio.anneal import *
from param_helpers import *

@parametrize_from_file(
        schema=Schema({
            'given': [str],
            **with_swmb.error_or({
                'expected': with_swmb.eval,
            }),
        }),
)
def test_parse_oligo_pairs_docopt(given, expected, error):
    with error:
        assert parse_oligo_pairs_docopt(given) == expected

@parametrize_from_file(
        schema=Schema({
            'given': freezerbox.parse_fields,
            **with_swmb.error_or({
                'expected': with_swmb.eval,
            }),
        }),
)
def test_parse_oligo_pair_freezerbox(given, expected, error):
    with error:
        assert parse_oligo_pair_freezerbox(given) == expected

@parametrize_from_file(
        schema=Schema({
            'given': str,
            **with_swmb.error_or({
                'expected': with_py.eval,
            }),
        }),
)
def test_parse_oligo_stock_docopt(given, expected, error):
    with error:
        assert parse_oligo_stock_docopt(given) == expected

