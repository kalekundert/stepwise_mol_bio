import pytest
import parametrize_from_file

from stepwise_mol_bio.digest import *
from more_itertools import one
from warnings import catch_warnings, simplefilter
from param_helpers import *

with catch_warnings():
    simplefilter('ignore', DeprecationWarning)
    import requests_testing

test_reaction()
test_protocol()
test_freezerbox_make()
test_freezerbox_attrs()

@parametrize_from_file(
        schema=Schema({
            'seq': with_py.eval,
            'enzymes': with_py.eval,
            'is_circular': with_py.eval,
            'target_size': with_py.eval,
            **with_swmb.error_or({
                'product': with_py.eval,
                'products': with_py.eval,
            }),
        }),
)
def test_calc_digest_products(seq, enzymes, is_circular, target_size, product, products, error):
    with error:
        assert products == calc_digest_products(
                seq, enzymes,
                is_circular=is_circular,
        )
    with error:
        assert product == calc_digest_product(
                seq, enzymes,
                is_circular=is_circular,
                target_size=target_size,
        )

def test_neb_restriction_enzyme_database(tmp_path):
    cache_path = tmp_path / 'cache.json'
    db = NebRestrictionEnzymeDatabase(cache_path)

    # This doesn't test the case where the internet is inaccessible.

    assert db['EcoRI'] == db['ecori'] == {
            'amt': '50000/10000/50000/10000/5000 units',
            'bcn': 'R0101',
            'blueWhite': True,
            'buf1': 25,
            'buf2': 100,
            'buf3': 50,
            'buf4': 50,
            'buf5': 0,
            'cat': 'R0101',
            'clonedAtNEB': True,
            'concentration': 20000,
            'cpg': True,
            'dam': False,
            'dcm': False,
            'displayName': 'EcoRI',
            'engineered': False,
            'epimarkValidated': False,
            'heatInactivationTemp': 65,
            'heatInactivationTime': 20,
            'hfEnzyme': False,
            'incubateTemp': 37,
            'methylationSensitivity': [
                'cpg (Blocked by Some Combinations of Overlapping)',
            ],
            'mul': False,
            'name': 'EcoRI',
            'plainname': 'EcoRI',
            'recombinant': True,
            'recommBuffer': 'NEBuffer EcoRI/SspI',
            'reducedStarActivity': False,
            'size': 'L/S/M/T/V ',
            'star1': False,
            'star2': True,
            'star3': False,
            'star4': True,
            'star5': False,
            'supplement': {'atp': 0.0, 'bsa': 0.0, 'dtt': 0.0, 'enzact': 0.0, 'sam': 0.0},
            'thermophilic': False,
            'timeSaver': True,
            'url': 'https://www.neb.com/products/r0101-ecori',
    }

    with pytest.raises(ConfigError) as err:
        db['EcoRJ']

    assert err.match(r"no such enzyme 'EcoRJ'")
    assert err.match(r"successfully downloaded the most recent restriction enzyme data from NEB \(in case 'EcoRJ' is a new enzyme\)")
    assert err.match(r"did you mean: 'EcoRI'")

@requests_testing.activate
def test_neb_restriction_enzyme_database_offline(tmp_path):
    cache_path = tmp_path / 'cache.json'

    with pytest.raises(ConfigError) as err:
        NebRestrictionEnzymeDatabase(cache_path)

    assert err.match("failed to download")
    assert err.match("URL: http://nebcloner.neb.com/data/reprop.json")

@parametrize_from_file(
        schema=Schema({
            'enzymes': with_py.eval,
            'expected': str,
        }),
)
def test_pick_compatible_buffer(enzymes, expected):
    assert pick_compatible_buffer(enzymes) == expected

def test_reaction_unknown_supplement():
    mock_db = {
            'XyzX': {
                'name': 'XyzX',
                'displayName': 'XyzX',
                'concentration': 10_000,
                'recommBuffer': 'rCutSmart Buffer',
                "supplement": {
                    "atp": 0.0,
                    "bsa": 0.0,
                    "dtt": 0.0,
                    "enzact": 20.0, 
                    "sam": 0.0,
                },
            },
    }
    app = RestrictionDigest.from_tags(['x1'], ['XyzX'], mock_db)

    with pytest.raises(ConfigError) as err:
        app.reaction

    assert err.match(r"'XyzX' requires an unknown supplement: 'enzact'")

@parametrize_from_file(key='test_freezerbox_make')
def test_protocol_product(db, tags, expected):
    if len(tags) != 1:
        pytest.skip()

    app = RestrictionDigest.from_product(one(tags))
    app.db = eval_db(db)

    assert match_protocol(app, expected)

