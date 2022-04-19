#!/usr/bin/env python3

import parametrize_from_file
from stepwise_mol_bio import ivt
from param_helpers import *

@parametrize_from_file
def test_dnase_reaction(app, expected):
    app = exec_app(app)
    assert app.dnase_reaction.format_text() == expected

@parametrize_from_file(
        schema=Schema({
            'app': str,
            'expected': eval,
        }),
)
def test_short(app, expected):
    app = exec_app(app)
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
            **with_swmb.error_or({
                'expected': str,
            }),
        }),
)
def test_transcribe(template, expected, error):
    with error:
        assert ivt.transcribe(template) == expected

