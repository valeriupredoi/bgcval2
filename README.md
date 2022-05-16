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
analysis_p2p u-bc179 level2 2010
```
Once these have completed, a summary HTML page can be generated with the command:

```
makeReport u-bc179 2010
```
This produces an html5 mobile-friendly website which can be opened using your
browser of choice.


### Available executables

Executable name | What it does | Command
:--------------:|:------------:|:------------:
`analysis_timeseries` | runs timeseries analysis for single model run. | analysis_timeseries jobID key
`analysis_p2p` | runs point to point analysis of model against observational dataset. | analysis_p2p jobID YEAR
`bgcval` | runs time series and point to point. | bgcval jobID
`makeReport` | makes the single model HTML report. | makeReport jobID
`analysis_compare` | runs comparison of multiple single jobs  | analysis_compare


Time series analysis
--------------------

This is an analysis that investigates the time development of specific marine
physics and Biogeochemistry fields in the given model, and then compares them
against historical observations.

The command to run it is `analysis_timeseries jobID key`, where jobID is a mass
job id, such a `u-ab123`, and the key is a pre-defined key, which generates a
list of variables.

Key | What it is | Description
:--------------:|:------------:|:------------:
`kmf` | Key Metrics First | A short and quick list of the most important metrics.
`physics` | Physics | A comprehensive list of physical metrics.
`bgc` | Biogeochemistry | A comprehensive list of biogeochemical metrics.
`level1` | Level 1 | A comprehensive list of physical and biogeochemical metrics.
`debug` | Debug | A very short list of a couple keys to test code changes.
`fast` | UKESM1-fast  | A list of metrics tailored to the UKESM1-Fast model.


Note that there may be some overlap between the contents of these keys.


Point to point analysis
-----------------------

This is an analysis that compares a specific year of the model
against several climatological datasets, and produces comparison
maps, histograms and scatter plots.


The command to run it is `analysisp2p jobID YEAR key`, where jobID is a mass
job id, such a `u-ab123`, `YEAR` is the single year to investigate
and the key is a pre-defined key, which generates a list of variables.

Key | What it is | Description
:--------------:|:------------:|:------------:
`physics` | Physics | A comprehensive list of physical metrics.
`level2` | Level 2 | A comprehensive list of physical and biogeochemical metrics.
`debug` | Debug | A very short list of a couple keys to test code changes.

Note that there may be some overlap between the contents of these keys.


Single Model report
-------------------

Once an analysis has run, either time series or point to point, a report
can be generated from this output, using the command:
```
makeReport jobID year
```
This  gnerated an HTML5 mobile-friendlyt report, summarising the output of a
single model run.


Multi-model comaprison report
-----------------------------

Once several models have been analysed using the time series analysis,
their time development can be compared using the analysis_compare command.
This is currently hardwired, but work in on-going.

Documentation
=============

See available Sphinx [documentation](https://htmlpreview.github.io/?https://github.com/valeriupredoi/bgcval2/blob/main/doc/build/index.html). To build locally the documentation run:

```
sphinx-build -Ea doc doc/build
```

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
