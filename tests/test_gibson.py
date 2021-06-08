#!/usr/bin/env python3

from stepwise.testing import check_command, disable_capture
from schema_helpers import *
from freezerbox.stepwise import Make

@parametrize_from_file(schema=app_expected_protocol)
def test_python_protocol(app, expected):
    assert match_protocol(app, expected)

@pytest.mark.slow
@parametrize_from_file(schema=cmd_stdout_stderr)
def test_cli_protocol(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout)

@parametrize_from_file(schema=db_tags_expected_protocol)
def test_freezerbox_protocol(db, tags, expected, disable_capture):
    app = Make(db, tags)
    assert match_protocol(app, expected, disable_capture)

@parametrize_from_file(schema=db_expected)
def test_freezerbox_attrs(db, expected):
    for tag in expected:
        assert db[tag].dependencies == expected[tag]['dependencies']
        assert db[tag].conc == expected[tag]['conc']



