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
    check_command(cmd, stdout=stdout, stderr=stderr)

