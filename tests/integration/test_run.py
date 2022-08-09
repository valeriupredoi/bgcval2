"""Tests for BGCVal2 CLI.

Includes a context manager to temporarily modify sys.argv
"""
import contextlib
import copy
import functools
import sys
from unittest.mock import patch

import pytest

from bgcval2 import (
    analysis_timeseries,
    bgcval,
    download_from_mass,
    bgcval2_make_report,
    analysis_compare
)
from bgcval2.analysis_timeseries import main
from bgcval2.download_from_mass import main
from bgcval2.bgcval2_make_report import main
from bgcval2.analysis_compare import main


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
def test_run_analysis_timeseries_command():
    """Test run command."""
    with arguments('analysis_timeseries', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "analysis_timeseries: error: argument \
        -j/--jobID: expected at least one argument"
    with arguments('analysis_timeseries', '--jobID', '--keys'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
            assert err == pytest_wrapped_e


@patch('bgcval2.bgcval.run', new=wrapper(bgcval))
def test_run_bgcval_command():
    """Test run command."""
    with arguments('bgcval', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "bgcval: error: the following arguments are required: -j/--jobID"
    with arguments('bgcval', '--job-id'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
            assert err == pytest_wrapped_e


@patch('bgcval2.download_from_mass.main', new=wrapper(download_from_mass))
def test_download_from_mass_command():
    """Test run command."""
    with arguments('download_from_mass', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "download_from_mass: error: the following \
        arguments are required: -i/--job-id"
    with arguments('bgcval', '--job-id'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
            assert err == pytest_wrapped_e


@patch('bgcval2.analysis_compare.main', new=wrapper(analysis_compare))
def test_analysis_compare_command():
    """Test run command."""
    with arguments('analysis_compare', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "download_from_mass: error: the following \
        arguments are required: -i/--job-id"
    with arguments('bgcval', '--job-id'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()
            assert err == pytest_wrapped_e
