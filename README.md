[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

BGCVal
======

This is the Python3 (Python 3.9+) version of [BGCVal](https://gmd.copernicus.org/articles/11/4215/2018/).

Environment and installation
============================

**Supported Operating Systems so far**: Linux/UNIX

Grab [miniconda3](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html), install `mamba` in the `base` environment:

```
conda install -c conda-forge mamba
```

then create the `bgcval` environment and activate it:

```
mamba env create -n bgcval -f environment.yml
conda activate bgcval
```

then install the development dependencies:

```
pip install -e .[develop]
```

Running the tool
================

The tool has a number of executables one can invoke individually, e.g.:

```
analysis_timeseries u-bc179 level2
```

or to make the summary HTML page:

```
makeReport u-bc179 2010
```

### Available executables

Executable name | What it does
:--------------:|:------------:
`analysis_compare` | runs comparison
`analysis_level3_amoc` | something
`analysis_level3_dms` | something
`analysis_level3_omz` | something
`analysis_level3_sameYear` | something
`analysis_p2p` | runs p2p
`analysis_timeseries` | runs timeseries
`bgcval` | runs everything
`makeReport` | makes the HTML report

