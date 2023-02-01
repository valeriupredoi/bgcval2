"""Tests for BGCVal2 CLI.

Includes a context manager to temporarily modify sys.argv
"""
import contextlib
import copy
import functools
import sys
from unittest.mock import patch
from io import StringIO

import pytest

from bgcval2 import (
    analysis_timeseries,
    bgcval,
    download_from_mass,
    bgcval2_make_report,
    analysis_compare
)
from bgcval2.analysis_timeseries import main as analysis_timeseries_main
from bgcval2.download_from_mass import main as download_from_mass_main
#from bgcval2.bgcval import run as bgcval_main
from bgcval2.bgcval2_make_report import main as bgcval2_make_report_main
from bgcval2.analysis_compare import main as analysis_compare_main


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


@contextlib.contextmanager
def capture_sys_output():
    capture_out, capture_err = StringIO(), StringIO()
    current_out, current_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = capture_out, capture_err
        yield capture_out, capture_err
    finally:
        sys.stdout, sys.stderr = current_out, current_err


@patch('bgcval2.analysis_timeseries.main', new=wrapper(analysis_timeseries))
def test_run_analysis_timeseries_command():
    """Test run command."""
    with arguments('analysis_timeseries', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            analysis_timeseries_main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "analysis_timeseries: error: the following arguments " \
        "are required: -j/--jobID"
    with arguments('analysis_timeseries'):
        with pytest.raises(SystemExit) as cm, capture_sys_output() \
            as (stdout, stderr):
            analysis_timeseries_main()
        assert err in str(stderr.getvalue())
    err = "--jobID: expected at least one argument"
    with arguments('analysis_timeseries', '--jobID', '--keys'):
        with pytest.raises(SystemExit) as cm, capture_sys_output() \
            as (stdout, stderr):
            analysis_timeseries_main()
        assert err in str(stderr.getvalue())


#@patch('bgcval2.bgcval.run', new=wrapper(bgcval))
#def test_run_bgcval_command():
#    """Test run command."""
#    with arguments('bgcval', '--help'):
#        with pytest.raises(SystemExit) as pytest_wrapped_e:
#            bgcval_main()
#        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "the following arguments are required: -i/--job-id\n"
    with arguments('bgcval'):
        with pytest.raises(SystemExit) as cm, capture_sys_output() \
            as (stdout, stderr):
            bgcval_main()
        assert err in str(stderr.getvalue())


@patch('bgcval2.download_from_mass.main', new=wrapper(download_from_mass))
def test_download_from_mass_command():
    """Test run command."""
    with arguments('download_from_mass', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            download_from_mass_main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "the following arguments are required: -j/--jobID"
    with arguments('download_from_mass'):
        with pytest.raises(SystemExit) as cm, capture_sys_output() \
            as (stdout, stderr):
            download_from_mass_main()
        assert err in str(stderr.getvalue())


@patch('bgcval2.analysis_compare.main', new=wrapper(analysis_compare))
def test_analysis_compare_command():
    """Test run command."""
    with arguments('analysis_compare', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            analysis_compare_main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "the following arguments are required: -y/--compare_yml"
    with arguments('analysis_compare'):
        with pytest.raises(SystemExit) as cm, capture_sys_output() \
            as (stdout, stderr):
            analysis_compare_main()
        assert err in str(stderr.getvalue())


@patch('bgcval2.bgcval2_make_report.main', new=wrapper(bgcval2_make_report))
def test_bgcval2_make_report_command():
    """Test run command."""
    with arguments('analysis_compare', '--help'):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            bgcval2_make_report_main()
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
    err = "the following arguments are required: -i/--job-id"
    with arguments('bgcval2_make_report'):
        with pytest.raises(SystemExit) as cm, capture_sys_output() \
            as (stdout, stderr):
            bgcval2_make_report_main()
        assert err in  str(stderr.getvalue())
    err = "argument -r/--report: expected one argument"
    with arguments('bgcval2_make_report', '--job-id DUM', '--report'):
        with pytest.raises(SystemExit) as cm, capture_sys_output() \
            as (stdout, stderr):
            bgcval2_make_report_main()
        assert err in str(stderr.getvalue())
