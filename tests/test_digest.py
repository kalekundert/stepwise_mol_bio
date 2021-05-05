#!/usr/bin/env python3

import pytest
import parametrize_from_file
from stepwise_mol_bio.digest import *
from freezerbox import Database, Plasmid
from schema_helpers import *
from warnings import catch_warnings, simplefilter

with catch_warnings():
    simplefilter('ignore', DeprecationWarning)
    import requests_testing

@parametrize_from_file(
        schema=Schema({
            'seq': eval,
            'enzymes': eval,
            'is_circular': eval,
            'target_size': eval,
            **error_or({
                'product': eval,
                'products': eval,
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
            'amt': '50000/10000/50000/10000 units',
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
            'size': 'L/S/M/T ',
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
            'enzymes': eval_with(),
            'expected': str,
        }),
)
def test_pick_compatible_buffer(enzymes, expected):
    assert pick_compatible_buffer(enzymes) == expected

@parametrize_from_file(schema=app_expected_reaction)
def test_reaction(app, expected):
    assert app.reaction.format_text() == expected

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
    app = RestrictionDigest(['x1'], ['XyzX'], mock_db)

    with pytest.raises(ConfigError) as err:
        app.reaction

    assert err.match(r"'XyzX' requires an unknown supplement: 'enzact'")

@parametrize_from_file(schema=app_expected_protocol)
def test_protocol(app, expected):
    actual = app.protocol.format_text()
    print(actual)
    for x in expected:
        assert x in actual

@parametrize_from_file(schema=db_expected)
def test_make(db, expected):
    for tag in expected:
        assert db[tag].seq == expected[tag]['seq']
        assert db[tag].dependencies == expected[tag]['dependencies']

