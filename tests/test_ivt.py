#!/usr/bin/env python3

import parametrize_from_file
from stepwise_mol_bio import ivt
from schema_helpers import *

@parametrize_from_file(schema=app_expected_reaction)
def test_reaction(app, expected):
    assert app.reaction.format_text() == expected

@parametrize_from_file(schema=app_expected_reaction)
def test_dnase_reaction(app, expected):
    assert app.dnase_reaction.format_text() == expected

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    actual = app.protocol.format_text()
    print(actual)
    for x in expected:
        assert x in actual

@parametrize_from_file(schema=app_expected)
def test_short(app, expected):
    assert app.short == expected

@parametrize_from_file(schema=Schema({str: eval}))
def test_pick_by_short(values, is_short, expected):
    assert ivt.pick_by_short(values, is_short) == expected

@parametrize_from_file(schema=Schema({str: eval}))
def test_affected_by_short(values, expected):
    assert ivt.affected_by_short(values) == expected

@parametrize_from_file(
        schema=Schema({
            'template': str,
            **error_or({
                'expected': str,
            }),
        }),
)
def test_transcribe(template, expected, error):
    with error:
        assert ivt.transcribe(template) == expected

@parametrize_from_file(schema=db_expected)
def test_make(db, expected):
    for tag in expected:
        assert db[tag].seq == expected[tag]['seq']
        assert db[tag].dependencies == expected[tag]['dependencies']
        assert db[tag].molecule == 'RNA'
        assert db[tag].is_single_stranded == True


