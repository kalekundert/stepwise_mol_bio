from stepwise_mol_bio.direct_dilution import *
from param_helpers import *

@parametrize_from_file
def test_dilution_table():
    assert app.dilution_table.rows == expected

@parametrize_from_file(
        schema=[
            cast(
                target_concs=Schema([Int]),
                max_dilution=int,
                expected=Schema({Int: Int}),
            ),
            with_swmb.error_or('expected'),
        ],
)
def test_pick_stock_concs(target_concs, max_dilution, expected, error):
    with error:
        assert pick_stock_concs(target_concs, max_dilution) == expected

test_protocol()
test_cli()


