#!/usr/bin/env python3

import parametrize_from_file
from schema_helpers import *

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    actual = app.protocol.format_text()
    print(actual)

    i = 0
    for x in expected:
        j = actual.find(x, i)
        assert j != -1, x
        i = j + len(x)

