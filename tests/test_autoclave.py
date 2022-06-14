import parametrize_from_file
from stepwise_mol_bio.autoclave import *
from param_helpers import *

test_protocol()
test_cli()
test_freezerbox_make()

@parametrize_from_file(
        schema=cast(volume_mL=int, expected=int),
)
def test_calc_sterilization_time(volume_mL, expected):
    app = Autoclave()
    app.volume_mL = volume_mL
    assert app.time_min == expected

