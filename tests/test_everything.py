#!/usr/bin/env python3

import pytest
import parametrize_from_file
import nestedtext as nt

from stepwise.testing import check_command, disable_capture
from freezerbox.stepwise import Make
from functools import partial
from param_helpers import *

def add_protocol_to_id(params, context):
    import re
    prefix = re.match('test_(.*).nt', context.path.name).group(1)
    prefix = prefix.replace('_', '-')

    for param in params:
        if 'id' in param:
            param['id'] = f'{prefix}-{param["id"]}'
        else:
            param['id'] = prefix

    return params

def load_nt_ignore_missing(path):
    suite = nt.load(path)
    optional_tests = [
            'test_reaction',
            'test_protocol',
            'test_cli',
            'test_freezerbox_make',
            'test_freezerbox_attrs',
    ]
    for test in optional_tests:
        suite.setdefault(test, [])
    return suite

parametrize_from_files = partial(
        parametrize_from_file,
        path=[
            'test_aliquot.nt',
            'test_anneal.nt',
            'test_autoclave.nt',
            'test_digest.nt',
            'test_gel.nt',
            'test_gibson.nt',
            'test_golden_gate.nt',
            'test_ivt.nt',
            'test_ivtt.nt',
            'test_kld.nt',
            'test_laser_scanner.nt',
            'test_ligate.nt',
            'test_lyophilize.nt',
            'test_miniprep.nt',
            'test_page_purify.nt',
            'test_pcr.nt',
            'test_sequence.nt',
            'test_spin_cleanup.nt',
            'test_transform.nt',
        ],
        loaders={'.nt': load_nt_ignore_missing},
        preprocess=add_protocol_to_id,
)

@parametrize_from_files()
def test_reaction(app, expected):
    app = exec_app(app)
    assert app.reaction.format_text() == expected

@parametrize_from_files(
        schema=Schema({
            'app': str,
            **with_swmb.error_or({
                'expected': [str],
                'forbidden': [str],
            }),
        }),
)
def test_protocol(app, expected, forbidden, error):
    app = exec_app(app)
    with error:
        assert match_protocol(app, expected, forbidden)

@pytest.mark.slow
@parametrize_from_files(
        schema=Schema({
            'cmd': str,
            Optional('stdout', default='^$'): str,
            Optional('stderr', default='^$'): str,
        }),
)
def test_cli(cmd, stdout, stderr):
    check_command(cmd, stdout=stdout, stderr=stderr)

@parametrize_from_files(
        schema=Schema({
            'db': {str: str},
            Optional('tags', default=[]): [str],
            'expected': [str],
        }),
)
def test_freezerbox_make(db, tags, expected, disable_capture):
    db = eval_db(db)
    tags = tags or list(db.keys())
    app = Make(db, tags)
    assert match_protocol(app, expected, capture=disable_capture)

@parametrize_from_files(
        schema=Schema({
            'db': {str: str},
            'expected': {str: {str: str}},
            Optional('errors', default={}): {str: {str: str}}
        }),
)
def test_freezerbox_attrs(db, expected, errors):
    db = eval_db(db)
    expected = with_swmb.eval(expected)

    for tag in expected:
        for attr in expected[tag]:
            assert getattr(db[tag], attr) == expected[tag][attr]

    for tag in errors:
        for attr in errors[tag]:
            with with_swmb.error(errors[tag][attr]):
                getattr(db[tag], attr)

