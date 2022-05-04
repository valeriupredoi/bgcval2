[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Github Actions Test](https://github.com/valeriupredoi/bgcval2/actions/workflows/run-tests.yml/badge.svg)](https://github.com/valeriupredoi/bgcval2/actions/workflows/run-tests.yml)

BGCVal
======

This is the Python3 (Python 3.8 and 3.9) version of [BGCVal](https://gmd.copernicus.org/articles/11/4215/2018/).

**This is a fully deployable Python3 package.**

Suport for Python 3.10 is not yet enabled due to the current use of Basemap, that is obsolete, but still usable with Python 3.8-3.9.

Environment and installation
============================

**Supported Operating Systems so far**: Linux/UNIX

To install locally:

- get the `git` file of this repository:

```
git clone https://github.com/valeriupredoi/bgcval2.git
```

- then grab [miniconda3](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html), install `mamba` in the `base` environment:

```
conda install -c conda-forge mamba
```

- then create the `bgcval2` environment and activate it:

```
cd bgcval2
mamba env create -n bgcval2 -f environment.yml
conda activate bgcval2
```

- then install the development dependencies and the tool itself:

```
pip install -e .[develop]
```

Running the tool
================

The tool has a number of executables one can invoke individually, e.g.:

```
analysis_timeseries u-bc179 level1
```

or to make the summary HTML page:

```
makeReport u-bc179 2010
```

### Available executables

Executable name | What it does | Command
:--------------:|:------------:|:------------:
`analysis_timeseries` | runs timeseries analysis for single model run. | analysis_timeseries jobID key
`analysis_p2p` | runs point to point analysis of model against observational dataset. | analysis_p2p jobID YEAR
`bgcval` | runs time series and point to point. | bgcval jobID
`makeReport` | makes the single model HTML report. | makeReport jobID
`analysis_compare` | runs comparison of multiple single jobs  | analysis_compare


Appendix
========

## Python2 to Python3 Migration

Migration from the orginal BGCVal code, which was Python2, has been done with the `2to3` tool:

- Install 2to3:

```
pip install 2to3
```
- Usage: use the 3.9 extension and write to disk option:

```
2to3-3.9 script.py -w
```

Remove the backup `.py.bak` files or stash them.

## View `html` output pages straight into GitHub

Github renders html straight into the browser, with source being any given `index.html`. For that, one needs to prepend the `https://htmlpreview.github.io/?` before the http address of the `index.html`, example:

https://htmlpreview.github.io/?https://github.com/valeriupredoi/reports-bgc-val/blob/main/bgc-val-initial-report/u-bc179/index.html
