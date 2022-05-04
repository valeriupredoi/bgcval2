"""Tests for BGCVal2 CLI.

Includes a context manager to temporarily modify sys.argv
"""
import contextlib
import copy
import functools
import sys
from unittest.mock import patch

import pytest

from bgcval2 import analysis_timeseries, bgcval
from bgcval2.analysis_timeseries import main


def wrapper(f):
    @functools.wraps(f)
    def empty(*args, **kwargs):  # noqa
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


def test_setargs():
    original = copy.deepcopy(sys.argv)
    with arguments('testing', 'working', 'with', 'sys.argv'):
        assert sys.argv == ['testing', 'working', 'with', 'sys.argv']
    assert sys.argv == original


@patch('bgcval2.analysis_timeseries.main', new=wrapper(analysis_timeseries))
def test_run_analysis_timeseries():
    """Test run command."""
    with arguments('analysis_timeseries', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0


@patch('bgcval2.bgcval.run', new=wrapper(bgcval))
def test_run_bgcval():
    """Test run command."""
    with arguments('bgcval', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
