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

@parametrize_from_file(schema=db_expected_protocol)
def test_freezerbox_protocol(db, expected, disable_capture):
    tags = list(db.keys())
    app = Make(db, tags)
    assert match_protocol(app, expected, disable_capture)



