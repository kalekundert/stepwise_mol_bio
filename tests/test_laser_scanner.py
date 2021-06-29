#!/usr/bin/env python3

from param_helpers import *
from stepwise.testing import check_command

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    assert match_protocol(app, expected)

@pytest.mark.slow
@parametrize_from_file(schema=cmd_stdout_stderr)
def test_cli(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout)


