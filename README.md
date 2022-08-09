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

Test that the tool has been installed correctly with:
```
analysis_compare -h
```
which should print the module information and instructions on how to run the tool. 


### Available executables

Executable name | What it does | Command
:--------------:|:------------:|:------------:
`analysis_timeseries` | runs timeseries analysis for single model run. | analysis_timeseries jobID key
`analysis_p2p` | runs point to point analysis of model against observational dataset. | analysis_p2p jobID YEAR
`bgcval` | runs time series and point to point. | bgcval jobID
`bgcval2_make_report` | makes the single model HTML report. | bgcval2_make_report jobID
`analysis_compare` | runs comparison of multiple single jobs  | analysis_compare



Running the tool to compare multiple jobs
=========================================

The time developmenmt of several models can be compared 
and summarized in a single comparison report html. 
This report can be generated with a single command, based on a simple yml file input:

```
analysis_compare --compare_yml comparison_recipe.yml
```

Example input yaml files exist in the `input_yml` directory.
However, there are a few key values:


In this yml file, the structure is:
```
---
name: <Analysis name string>
do_analysis_timeseries: <bool>
do_mass_download: <bool>
master_suites: <str>

jobs:
   <jobID1>:
      description: <descrption of the first job>
      colour: red
      thickness: 0.7
      linestyle: '-'
      shifttime: 0.
      suite: physics
   <jobID2>:
      description: <descrption of the second job>
      ...
```

These values are:
 - `name`: 
   - The name of the analysis.
   - This will be the path of output report.
 - `do_analysis_timeseries`: 
   -  A Boolean value to run or skip the single model timeseries. 
   -  Set to False if the single model analysis has already completed.
 - `do_mass_download`: 
   - A boolean value to run the mass download. 
   - This is not currently possible as we can only download mass file from mass-cli1 on jasmin.
   - See below for details on how to download jobs data.
 - `master_suites`:
   - A list of the type of analysis report to produce. 
   - Options are: `physics`, `bio`, `debug`.
   - Default is `['physics', 'bio',]`.
 - `jobs`:
   - A list of jobIDs, and some options on how they will appear in the final report. 
   - The options are:
     - `description`:
       - A description of job, which helps people differentiate the jobs in the final report.
     - `colour`: 
       - The colour of this jobs line in the report plots.
       - Default colour is a randomly generated hex colour.
     - `thickness`:
       - The thickness of this jobs line in the report plots.
       - default is 0.7
     - `linestyle`:
       - The linestyle of this jobs line in the report plots.
       - Accepts typical matplotlib line styles: `'solid', 'dashed', 'dotted', 'dashdot', '-'. ';', ', etc`
       - default is `-'`
     - `shiftime`:
       - A number in whole years by which the jobs line is shifted.
       - this is useful if jobs start from different initial years in spin up, for instance.
       - Default is `0.` ( no time shift ).
     - `suite`:
       - An analysis suite to run the analysis_timeseries.
       - See `analysis_timeseries` for more details.


A sample yaml exists in `input_yml/comparison_analysis_template.yml`,
which can be adapted to additional analyses.



Downloading data using MASS
===========================

Data can be downloaded and prepared for analysis using the `download_from_mass` bgcval2 tool,
with the command:
```
download_from_mass -j jobID
```
where `jobID` is one or more jobIDs.

This script will only work on jasmin's `mass-cli1` machine,
which is set up to interact with the Met Office Storate System MASS.

The optional flag `--dry-run` skips the download part of the script,
and generates a script in the `bgcval2/mass_scripts` directory.
From there, users can ssh to `mass-cli1` and execute this script:

```
# from login1.jasmin.ac.uk, ssh to the mass machine:
ssh -X  mass-cli1

# run script with:
source /path/to/bgcval2/is/here/bgcval2/mass_scripts/jobID.sh
```

Alternatively, the whole `download_from_mass` tool could be executed on the `mass-cli1` machine.

Several different keys can be included in the download if monthly data is required.
However, it's not recommended to include monthly data at this stage as that code
both to download, but also to run the monthly analysis is not currently tested.

This tool downloads the data, but also includes several functions which create symbolic links 
in the data's directory in order to accomodate incompatible changes in NEMO's output formatting.


Running the tool for a single job
=================================

The multi-job analysis described above can only do timeseries analysis.
To run an in-depth analysis of a single job, the following command can be run:

```
bgcval2 -j jobID
```

This will run a time series analysis, a point to point analysis, and 
publish the reports into a single job html5 report.


Alternatively, these tasks can be invoked individually, e.g.:

```
analysis_timeseries --jobID u-bc179 --keys kmf level1
analysis_p2p u-bc179 level2 2010
bgcval2_make_report u-bc179 2010

```
This produces an html5 mobile-friendly website which can be opened using your
browser of choice.



Time series analysis
--------------------

This is an analysis that investigates the time development of specific marine
physics and Biogeochemistry fields in the given model, and then compares them
against historical observations.

The command to run it is `analysis_timeseries --jobID jobID --keys key`, 
where jobID is one or more mass jobIDs, such a `u-ab123`.
The key is one or more pre-defined key, which generates a
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

WORK IN PROGRESS.

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

WORK IN PROGRESS.

Once an analysis has run, either time series or point to point, a report
can be generated from this output, using the command:
```
bgcval2_make_report jobID year
```
This  gnerated an HTML5 mobile-friendlyt report, summarising the output of a
single model run.


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
