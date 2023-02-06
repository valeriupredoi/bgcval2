"""
Run end-to-end test/dummy-data analyses.

Dummy data is produced with produce_dummy_data.py (needs iris),
but since BGCVal2 needs a LOT of data, we don't store it but
produce it on the fly. Some dummy data still needs some fixing
to fully work with various components of BGCVal2, hence some tests
here catch errors rather than run full end-to-end analyses.

"""
import contextlib
import functools
import pytest
import sys
from . import  bgcval2_test_data
import bgcval2

from bgcval2.analysis_compare import main
from unittest.mock import patch


def wrapper(f):
    @functools.wraps(f)
    def empty(*args, **kwargs):
        if kwargs:
            raise ValueError(f'Parameters not supported: {kwargs}')
        return True

    return empty


@contextlib.contextmanager
def arguments(*args):
    backup = sys.argv
    sys.argv = list(args)
    yield
    sys.argv = backup


@pytest.fixture(scope='session')
def test_create_files():
    bgcval2_test_data()


@patch('bgcval2.analysis_compare.main', new=wrapper(bgcval2.analysis_compare.main))
@pytest.mark.usefixtures('test_create_files')
def test_run_analysis_compare():
    """Patch and run the whole thing."""
    bgcval2_test_data()
    with arguments('analysis_compare', '-y', 'input_yml/integration_testing.yml'):
        main()
