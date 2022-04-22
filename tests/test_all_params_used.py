import pytest, re
import nestedtext as nt
from pathlib import Path

TEST_DIR = Path(__file__).parent

@pytest.mark.parametrize('test_path_nt', TEST_DIR.glob('test_*.nt'))
def test_all_params_used(test_path_nt):
    test_suite = nt.load(test_path_nt)
    test_script_py = test_path_nt.with_suffix('.py').read_text()

    for test in test_suite:
        pattern = fr'^{test}\(.*\)|^def {test}\(.*\)'
        assert re.search(pattern, test_script_py, re.MULTILINE)

