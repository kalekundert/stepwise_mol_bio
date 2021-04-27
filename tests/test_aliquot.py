#!/usr/bin/env python3

import parametrize_from_file
from schema_helpers import *

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    actual = app.protocol.format_text()
    print(actual)
    for x in expected:
        assert x in actual


