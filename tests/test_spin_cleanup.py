#!/usr/bin/env python3

import parametrize_from_file
from pytest import approx
from freezerbox import QueryError
from param_helpers import *

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    assert match_protocol(app, expected)

@parametrize_from_file
def test_freezerbox_attrs(db, expected):
    db = eval_db(db)
    expected = eval_swmb(expected)

    for tag in expected:
        assert db[tag].volume == expected[tag]['volume']

        if 'conc' in expected[tag]:
            assert db[tag].conc == expected[tag]['conc']
        else:
            with pytest.raises(QueryError):
                db[tag].conc

