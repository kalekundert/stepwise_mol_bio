import parametrize_from_file
from stepwise_mol_bio.anneal import *
from param_helpers import *

test_reaction()
test_protocol()
test_cli()
test_freezerbox_make()
test_freezerbox_attrs()

@parametrize_from_file(
        schema=[
            cast(expected=with_swmb.eval),
            with_swmb.error_or('expected'),
        ],
)
def test_parse_oligo_pairs_docopt(given, expected, error):
    with error:
        assert parse_oligo_pairs_docopt(given) == expected

@parametrize_from_file(
        schema=[
            cast(given=freezerbox.parse_fields, expected=with_swmb.eval),
            with_swmb.error_or('expected'),
        ],
)
def test_parse_oligo_pair_freezerbox(given, expected, error):
    with error:
        assert parse_oligo_pair_freezerbox(given) == expected

@parametrize_from_file(
        schema=[
            cast(expected=with_py.eval),
            with_swmb.error_or('expected'),
        ],
)
def test_parse_oligo_stock_docopt(given, expected, error):
    with error:
        assert parse_oligo_stock_docopt(given) == expected

