#!/usr/bin/env python3

import pytest
import parametrize_from_file

from stepwise.testing import check_command, disable_capture
from stepwise_mol_bio.sequence import *
from freezerbox.stepwise import Make
from param_helpers import *

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    actual = app.protocol.format_text()
    assert match_protocol(app, expected)

@pytest.mark.slow
@parametrize_from_file(schema=cmd_stdout_stderr)
def test_cli(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout, stderr=stderr)

@parametrize_from_file(schema=db_tags_expected_protocol)
def test_freezerbox_make(db, tags, expected, disable_capture):
    app = Make(db, tags)
    assert match_protocol(app, expected, disable_capture)

@parametrize_from_file(
        schema=Schema({
            'given': [str],
            'expected': {str: [str]},
        }),
)
def test_parse_reactions_docopt(given, expected):
    assert parse_reactions_docopt(given) == expected

@parametrize_from_file(
        schema=Schema({
            'given': eval_py,
            'expected': eval_py,
        }),
)
def test_merge_reactions(given, expected):
    assert merge_reactions(given) == expected

