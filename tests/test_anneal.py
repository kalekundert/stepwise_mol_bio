#!/usr/bin/env python3

import pytest
import parametrize_from_file

from stepwise.testing import check_command, disable_capture
from stepwise_mol_bio.anneal import *
from freezerbox.stepwise import Make
from param_helpers import *

@parametrize_from_file(schema=app_expected_reaction)
def test_reaction(app, expected):
    assert app.reaction.format_text() == expected

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    actual = app.protocol.format_text()
    assert match_protocol(app, expected)

@pytest.mark.slow
@parametrize_from_file(schema=cmd_stdout_stderr)
def test_cli(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout)

@parametrize_from_file(schema=db_tags_expected_protocol)
def test_freezerbox_make(db, tags, expected, disable_capture):
    app = Make(db, tags)
    assert match_protocol(app, expected, disable_capture)

@parametrize_from_file(schema=db_expected)
def test_freezerbox_attrs(db, expected):
    for tag in expected:
        assert db[tag].dependencies == expected[tag]['dependencies']
        assert db[tag].conc == expected[tag]['conc']
        assert db[tag].volume == expected[tag]['volume']
        assert db[tag].molecule == 'DNA'
        assert db[tag].is_double_stranded

@parametrize_from_file(
        schema=Schema({
            'given': [str],
            **error_or({
                'expected': eval_swmb,
            }),
        }),
)
def test_parse_oligo_pairs_docopt(given, expected, error):
    with error:
        assert parse_oligo_pairs_docopt(given) == expected

@parametrize_from_file(
        schema=Schema({
            'given': freezerbox.parse_fields,
            **error_or({
                'expected': eval_swmb,
            }),
        }),
)
def test_parse_oligo_pair_freezerbox(given, expected, error):
    with error:
        assert parse_oligo_pair_freezerbox(given) == expected

@parametrize_from_file(
        schema=Schema({
            'given': str,
            **error_or({
                'expected': eval_py,
            }),
        }),
)
def test_parse_oligo_stock_docopt(given, expected, error):
    with error:
        assert parse_oligo_stock_docopt(given) == expected

