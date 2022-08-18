"""Test funcs in analysis_timeseries."""
import bgcval2
import os
import pytest
import yaml

from bgcval2.analysis_timeseries import (
    load_function,
    load_key_file
)

def test_load_function(tmp_path):
    """Test load_function()."""
    dummy_func_file = tmp_path / 'myfuncfile.py'
    with open(dummy_func_file, "w") as fil:
        fil.write("import os\n\ndef myfunc():\n    return 22\n")
    convert = {"path": str(dummy_func_file), "function": "myfunc"}
    expected_func, expected_kwargs = load_function(convert)
    res = expected_func()
    assert res == 22
    assert expected_kwargs == {}


def test_load_function_kwargs(tmp_path):
    """Test load_function()."""
    dummy_func_file = tmp_path / 'myfuncfile2.py'
    with open(dummy_func_file, "w") as fil:
        fil.write("import bgcval2\n\ndef myfunc():\n    return 22\n")
    convert = {
        "path": str(dummy_func_file), "function": "myfunc",
        "cow": "moo"}
    expected_func, expected_kwargs = load_function(convert)
    res = expected_func()
    assert res == 22
    assert expected_kwargs == {"cow": "moo"}


# test globalVolumeMean
def test_load_function_internal(tmp_path):
    """Test load_function()."""
    package_path = os.path.dirname(bgcval2.__file__)
    func_file = os.path.join(package_path, 'functions', 'globalVolMean.py')
    convert = {
        "path": str(func_file), "function": "globalVolumeMean",
        "cow": "moo"}
    expected_func, expected_kwargs = load_function(convert)
    with pytest.raises(KeyError) as exc:
        res = expected_func("gt40.nc", keys={})
    assert "Needs an `areafile` in kwargs" in str(exc)
    assert expected_kwargs == {"cow": "moo"}


def test_load_key_file(tmp_path):
    """Test load_key_file() func."""
    package_path = os.path.dirname(os.path.dirname(bgcval2.__file__))
    class Object(object):
        pass
    paths = Object()
    paths.bgcval2_repo = package_path
    paths.orcaGridfn = str(tmp_path / "dummy_orca")
    with open(paths.orcaGridfn, "w") as fil:
        fil.write("free Willy!")
    paths.ModelFolder_pref = "moo"
    paths.ObsFolder = "ox"
    key = "temperature"
    jobID = "cow"
    runtime_output_dict = load_key_file(key, paths, jobID)
    expected_output_dict = dict()
    expected_output_dict['modeldetails'] = dict()
    expected_output_dict['name'] = 'Temperature'
    assert expected_output_dict['name'] == runtime_output_dict['name']
    expected_output_dict['units'] = 'degrees C'
    assert expected_output_dict['units'] == runtime_output_dict['units']
    expected_output_dict['dimensions'] = 3
    assert expected_output_dict['dimensions'] == runtime_output_dict['dimensions']
    expected_output_dict['layers'] = ['Surface']
    assert expected_output_dict['layers'] == runtime_output_dict['layers']
    expected_output_dict['gridFile'] = paths.orcaGridfn
    assert expected_output_dict['gridFile'] == runtime_output_dict['gridFile']
    expected_output_dict['modeldetails']['name'] = 'Temperature'
    assert expected_output_dict['modeldetails']['name'] == \
        runtime_output_dict['modeldetails']['name']
    expected_output_dict['modeldetails']['vars'] = ['thetao_con']
    assert expected_output_dict['modeldetails']['vars'] == \
        runtime_output_dict['modeldetails']['vars']
