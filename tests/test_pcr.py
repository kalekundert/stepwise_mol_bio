#!/usr/bin/env python3

import pytest
import parametrize_from_file

from stepwise_mol_bio.pcr import Pcr, find_amplicon
from pytest import approx
from schema_helpers import *

@parametrize_from_file(
        schema=Schema({
            str: And(str, str.strip),
            'is_linear': eval,
            **error_or({
                'expected': And(str, str.strip),
            }),
        }),
)
def test_find_amplicon(template, primer_1, primer_2, is_linear, expected, error):
    with error:
        amplicon = find_amplicon(template, primer_1, primer_2, is_linear)
        assert amplicon == expected.upper()

@parametrize_from_file(schema=app_expected_error)
def test_product_seqs(app, expected, error):
    with error:
        assert app.product_seqs == [x.upper() for x in expected]

@parametrize_from_file(schema=app_expected)
def test_anneal_temp_C(app, expected):
    assert app.anneal_temp_C == approx(expected)

@parametrize_from_file(schema=app_expected)
def test_extend_time_s(app, expected):
    assert app.extend_time_s == approx(expected)

@parametrize_from_file(schema=app_expected)
def test_master_mix(app, expected):
    assert app.master_mix == expected

@parametrize_from_file(schema=db_expected)
def test_make(db, expected):
    for tag in expected:
        assert db[tag].seq == expected[tag]['seq'].upper()
        assert db[tag].dependencies == expected[tag]['dependencies']
        assert db[tag].conc_ng_uL == expected[tag]['conc_ng_uL']
        assert db[tag].molecule == 'DNA'
        assert db[tag].is_double_stranded == True
        assert db[tag].is_circular == False

