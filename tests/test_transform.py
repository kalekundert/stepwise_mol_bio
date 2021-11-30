#!/usr/bin/env python3

import pytest
import parametrize_from_file

from stepwise.testing import check_command, disable_capture
from stepwise_mol_bio.transform import *
from freezerbox.stepwise import Make
from param_helpers import *

@parametrize_from_file(schema=app_expected_protocol_error)
def test_protocol(app, expected, error):
    with error:
        actual = app.protocol.format_text()
        assert match_protocol(app, expected)

@pytest.mark.slow
@parametrize_from_file(schema=cmd_stdout_stderr)
def test_cli(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout)

#@parametrize_from_file(schema=db_tags_expected_protocol)
@parametrize_from_file
def test_freezerbox_make(db, tags, expected, disable_capture):
    db = eval_db(db)
    app = Make(db, tags)
    assert match_protocol(app, expected, disable_capture)


