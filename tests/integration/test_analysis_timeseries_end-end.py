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

from . import bgcval2_test_data
import bgcval2

from bgcval2.analysis_timeseries import main
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


@patch('bgcval2.analysis_timeseries.main', new=wrapper(bgcval2.analysis_timeseries.main))
def test_run_analysis_timeseries_debug():
    """Patch and run the whole thing."""
    with arguments('analysis_timeseries', '--jobID', 'u-cp416debug',
                       '--keys', 'debug'):
        main()
