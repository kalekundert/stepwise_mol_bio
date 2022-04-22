import pytest
import parametrize_from_file

from stepwise_mol_bio.pcr import *
from pytest import approx
from param_helpers import *

test_freezerbox_make()
test_freezerbox_attrs()

@parametrize_from_file(
        schema=Schema({
            str: And(str, str.strip),
            'is_circular': eval,
            **with_swmb.error_or({
                'expected': And(str, str.strip),
            }),
        }),
)
def test_find_amplicon(template, primer_1, primer_2, is_circular, expected, error):
    with error:
        amplicon = find_amplicon(template, primer_1, primer_2, is_circular)
        assert amplicon == expected.upper()

@parametrize_from_file(
        schema=Schema({
            'given': str,
            **with_swmb.error_or({
                'expected': {str: str},
            }),
        })
)
def test_parse_amplicon(given, expected, error):
    with error:
        amplicon = parse_amplicon(given)
        assert amplicon.template.tag == expected['template']
        assert amplicon.fwd.tag == expected['fwd']
        assert amplicon.rev.tag == expected['rev']

@parametrize_from_file(
        schema=Schema({
            'given': str,
            **with_swmb.error_or({
                'expected': with_py.eval,
            }),
        })
)
def test_parse_primers(given, expected, error):
    with error:
        assert parse_primers(given) == expected

@parametrize_from_file
def test_reactions(app, expected):
    app = exec_app(app)
    pcr, primers = app.reaction

    assert str(pcr) == expected['pcr']

    if 'primers' in expected:
        assert str(primers) == expected['primers']
    else:
        assert primers is None

@parametrize_from_file(
        schema=Schema({
            'app': str,
            **with_swmb.error_or({
                'expected': with_py.eval,
            }),
        }),
)
def test_product_seqs(app, expected, error):
    app = exec_app(app)
    with error:
        assert app.product_seqs == [x.upper() for x in expected]

@parametrize_from_file(
        schema=Schema({
            'app': str,
            'expected': with_py.eval,
        }),
)
def test_anneal_temp_C(app, expected):
    app = exec_app(app)
    assert app.anneal_temp_C == approx(expected)

@parametrize_from_file(
        schema=Schema({
            'app': str,
            'expected': with_py.eval,
        }),
)
def test_extend_time_s(app, expected):
    app = exec_app(app)
    assert app.extend_time_s == approx(expected)

@parametrize_from_file(key='test_freezerbox_make')
def test_protocol_products(db, tags, expected):
    if len(tags) != 1:
        pytest.skip()

    app = Pcr.from_product(tags[0])
    app.db = eval_db(db)

    assert match_protocol(app, expected)

