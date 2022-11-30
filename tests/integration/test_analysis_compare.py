"""Integration test for compare."""
import contextlib
import functools
import pytest
import sys

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


@patch('bgcval2.analysis_compare.main', new=wrapper(bgcval2.analysis_compare.main))
def test_run_analysis_compare():
    """Patch and run the whole thing."""
    with pytest.raises(ValueError) as exc:
        with arguments('analysis_compare', '-y', 'input_yml/debug.yml'):
            main()
    expected_exc = "operands could not be broadcast together with shapes (3,0) (0,3)"
    assert expected_exc in str(exc.value)
