#!/usr/bin/env python3

import pytest
import parametrize_from_file

from stepwise.testing import disable_capture
from stepwise_mol_bio.pcr import *
from freezerbox.stepwise import Make
from pytest import approx
from param_helpers import *

@parametrize_from_file(
        schema=Schema({
            str: And(str, str.strip),
            'is_circular': eval,
            **error_or({
                'expected': And(str, str.strip),
            }),
        }),
)
def test_find_amplicon(template, primer_1, primer_2, is_circular, expected, error):
    with error:
        amplicon = find_amplicon(template, primer_1, primer_2, is_circular)
        assert amplicon == expected.upper()

@parametrize_from_file(
        schema=Schema({
            'given': str,
            **error_or({
                'expected': {str: str},
            }),
        })
)
def test_parse_amplicon(given, expected, error):
    with error:
        amplicon = parse_amplicon(given)
        assert amplicon.template.tag == expected['template']
        assert amplicon.fwd.tag == expected['fwd']
        assert amplicon.rev.tag == expected['rev']

@parametrize_from_file(
        schema=Schema({
            'given': str,
            **error_or({
                'expected': eval_with(),
            }),
        })
)
def test_parse_primers(given, expected, error):
    with error:
        assert parse_primers(given) == expected

@parametrize_from_file(
        schema=Schema({
            'app': exec_app,
            'expected': {str: str},
        }),
)
def test_reaction(app, expected):
    pcr, primers = app.reaction

    assert str(pcr) == expected['pcr']

    if 'primers' in expected:
        assert str(primers) == expected['primers']
    else:
        assert primers is None

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

@parametrize_from_file(
        key='test_freezerbox_make',
        schema=db_tags_expected_protocol,
)
def test_protocol_products(db, tags, expected):
    if len(tags) != 1:
        pytest.skip()

    app = Pcr.from_product(tags[0])
    app.db = db

    assert match_protocol(app, expected)

@parametrize_from_file(schema=db_tags_expected_protocol)
def test_freezerbox_make(db, tags, expected, disable_capture):
    app = Make(db, tags)
    assert match_protocol(app, expected, disable_capture)

@parametrize_from_file(schema=db_expected)
def test_freezerbox_attrs(db, expected):
    for tag in expected:
        assert db[tag].seq == expected[tag]['seq'].upper()
        assert db[tag].dependencies == expected[tag]['dependencies']
        assert db[tag].conc_ng_uL == expected[tag]['conc_ng_uL']
        assert db[tag].volume_uL == expected[tag]['volume_uL']
        assert db[tag].molecule == 'DNA'
        assert db[tag].is_double_stranded == True
        assert db[tag].is_circular == False

