"""Test funcs in analysis_timeseries."""
import bgcval2
import numpy as np
import os
import pytest
import yaml

from netCDF4 import Dataset
from bgcval2.analysis_timeseries import (
    load_function,
    load_key_file
)


def create_hyb_pres_file_without_ap(dataset, short_name):
    """Create dataset without vertical auxiliary coordinate ``ap``."""
    dataset.createDimension('time', size=1)
    dataset.createDimension('lev', size=2)
    dataset.createDimension('lat', size=3)
    dataset.createDimension('lon', size=4)
    dataset.createDimension('bnds', size=2)

    # Dimensional variables
    dataset.createVariable('time', np.float64, dimensions=('time',))
    dataset.createVariable('lev', np.float64, dimensions=('lev',))
    dataset.createVariable('lev_bnds', np.float64, dimensions=('lev', 'bnds'))
    dataset.createVariable('lat', np.float64, dimensions=('lat',))
    dataset.createVariable('lon', np.float64, dimensions=('lon',))
    dataset.variables['time'][:] = [0.0]
    dataset.variables['time'].standard_name = 'time'
    dataset.variables['time'].units = 'days since 6543-2-1'
    dataset.variables['lev'][:] = [1.0, 2.0]
    dataset.variables['lev'].bounds = 'lev_bnds'
    dataset.variables['lev'].standard_name = (
        'atmosphere_hybrid_sigma_pressure_coordinate')
    dataset.variables['lev'].units = '1'
    dataset.variables['lev_bnds'][:] = [[0.5, 1.5], [1.5, 3.0]]
    dataset.variables['lev_bnds'].standard_name = (
        'atmosphere_hybrid_sigma_pressure_coordinate')
    dataset.variables['lev_bnds'].units = '1'
    dataset.variables['lat'][:] = [-30.0, 0.0, 30.0]
    dataset.variables['lat'].standard_name = 'latitude'
    dataset.variables['lat'].units = 'degrees_north'
    dataset.variables['lon'][:] = [30.0, 60.0, 90.0, 120.0]
    dataset.variables['lon'].standard_name = 'longitude'
    dataset.variables['lon'].units = 'degrees_east'

    # Coordinates for derivation of pressure coordinate
    dataset.createVariable('b', np.float64, dimensions=('lev',))
    dataset.createVariable('b_bnds', np.float64, dimensions=('lev', 'bnds'))
    dataset.createVariable('ps', np.float64,
                           dimensions=('time', 'lat', 'lon'))
    dataset.variables['b'][:] = [0.0, 1.0]
    dataset.variables['b_bnds'][:] = [[-1.0, 0.5], [0.5, 2.0]]
    dataset.variables['ps'][:] = np.arange(1 * 3 * 4).reshape(1, 3, 4)
    dataset.variables['ps'].standard_name = 'surface_air_pressure'
    dataset.variables['ps'].units = 'Pa'
    dataset.variables['ps'].additional_attribute = 'xyz'

    # Variable
    dataset.createVariable(short_name, np.float32,
                           dimensions=('time', 'lev', 'lat', 'lon'))
    dataset.variables[short_name][:] = np.full((1, 2, 3, 4), 0.0,
                                               dtype=np.float32)
    dataset.variables[short_name].standard_name = (
        'cloud_area_fraction_in_atmosphere_layer')
    dataset.variables[short_name].units = '%'


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
    with pytest.raises(FileNotFoundError) as exc:
        res = expected_func("gt40.nc", keys={})
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

    # test bogus paths
    paths.ModelFolder_pref = "moo"
    paths.ObsFolder = "ox"
    key = "temperature"
    jobID = "cow"
    with pytest.raises(OSError) as exc:
        load_key_file(key, paths, jobID)
    assert "Base moo/cow is not a valid directory" in str(exc)

    # assemble OK paths, no files
    jobID = "cow"
    data_dir = os.path.join(str(tmp_path), jobID)
    os.mkdir(data_dir)
    paths.ModelFolder_pref = str(tmp_path)
    model_file = os.path.join(data_dir, "model")
    with open(model_file, "w") as fil:
        fil.write("cowabunga")
    paths.ObsFolder = str(tmp_path)
    obs_dir = os.path.join(str(tmp_path), "WOA", "annual")
    os.makedirs(obs_dir)
    key = "temperature"
    with pytest.raises(FileNotFoundError) as exc:
        runtime_output_dict = load_key_file(key, paths, jobID)
    assert "nemo_cowo_1y_*_grid-T.nc" in str(exc)

    # OK paths, OK files but they are not netCDF4
    model_file = os.path.join(data_dir, "nemo_cowo_1y_1990_grid-T.nc")
    with open(model_file, "w") as fil:
        fil.write("cowabunga")
    with pytest.raises(OSError) as exc:
        runtime_output_dict = load_key_file(key, paths, jobID)
    assert "Unknown file format - not netCDF" in str(exc)

    # finally create valid netCDF4 files (dummy, but valid)
    nc_path = os.path.join(data_dir, 'nemo_cowo_1y_1990_grid-T.nc')
    dataset = Dataset(nc_path, mode='w')
    create_hyb_pres_file_without_ap(dataset, 'ta')
    dataset.close()
    print(f"Saved {nc_path}")
    obs_nc_path = os.path.join(obs_dir, "woa13_decav_t00_01v2.nc")
    obs_dataset = Dataset(obs_nc_path, mode='w')
    create_hyb_pres_file_without_ap(obs_dataset, 'ta')
    obs_dataset.close()
    print(f"Saved {obs_nc_path}")
    runtime_output_dict = load_key_file(key, paths, jobID)

    expected_output_dict = dict()
    expected_output_dict['modeldetails'] = dict()
    expected_output_dict['name'] = 'Temperature'
    assert expected_output_dict['name'] == runtime_output_dict['name']
    expected_output_dict['units'] = 'degrees C'
    assert expected_output_dict['units'] == runtime_output_dict['units']
    expected_output_dict['dimensions'] = 3
    assert expected_output_dict['dimensions'] == runtime_output_dict['dimensions']
    expected_output_dict['layers'] = ['Surface', '500m']
    assert expected_output_dict['layers'] == runtime_output_dict['layers']
    expected_output_dict['gridFile'] = paths.orcaGridfn
    assert expected_output_dict['gridFile'] == runtime_output_dict['gridFile']
    expected_output_dict['modeldetails']['name'] = 'Temperature'
    assert expected_output_dict['modeldetails']['name'] == \
        runtime_output_dict['modeldetails']['name']
    expected_output_dict['modeldetails']['vars'] = ['thetao_con', 'thetao']
    assert expected_output_dict['modeldetails']['vars'] == \
        runtime_output_dict['modeldetails']['vars']
