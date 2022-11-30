[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Github Actions Test](https://github.com/valeriupredoi/bgcval2/actions/workflows/run-tests.yml/badge.svg)](https://github.com/valeriupredoi/bgcval2/actions/workflows/run-tests.yml)

![bgcval2logo](https://github.com/valeriupredoi/bgcval2/blob/main/doc/figures/BGCVal2-logo-2-DARK-200x177.png)

bgcval2
=======


bgcval2 is the updated and modernised version of BGC-val.
There are several updates over the previously published version of [BGCVal](https://gmd.copernicus.org/articles/11/4215/2018/).
The primary improvements are to the user interface, the ease of use, ease of installation,
and a general modernisation of the core approach, including the move to python3, 
and using a conda enviroment.

This work was funded through WP1 of the Terrafirma project.


Current version notes:

- Suport for Python 3.10 is not yet enabled due to the current use of Basemap, that is obsolete, but still usable with Python 3.8-3.9.

Environment and installation
============================

**Supported Operating Systems so far**: Linux/UNIX

To install locally:

- get the `git` file of this repository:

```
git clone https://github.com/valeriupredoi/bgcval2.git
```

- then grab [miniconda3](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
- This can be installed using (bearing in mind that your file name may differ), and follow the on screen instructions:
```
  bash Miniconda3-py39_4.12.0-Linux-x86_64.sh
```

- Once you have a local version of miniconda3 installed, install `mamba` in the `base` environment:

```
conda install -c conda-forge mamba
```

- mamba is a replacement for conda that is faster, but works exactly the same way. Ie the command “conda env” -> becomes “mamba env”. However, this step is somewhat optional. If you can’t install mamba, but have conda working, just replace “mamba” with “conda” below.


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


On the jasmin computational system, members of the esmeval working group should be able to run the debug analysis script:
```
analysis_compare -y input_yml/debug.yml
```
This script performs an analysis of two small physics-only UKESM development jobIDs,
using the debug suite.


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
auto_download: <bool>

jobs:
   <jobID1>:
      description: <descrption of the first job>
      colour: red
      thickness: 0.7
      linestyle: '-'
      shifttime: 0.
      suite: physics
      auto_download: False

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
 - `auto_download`: 
   - If True, adds the jobID to a list of paths which will download automatically over night
   - The script it renewed every time that the analysis\_compare is run.
   - If the job is not run, older jobs are removed from the nightly download after some time.
   - This boolean flag can be set at the top level or for individual jobs.
   - Please set it to false if your job has completed and all relevant has been downloaded.
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
       - One or more analysis suites to run the analysis_timeseries.
       - Options include: `bgc` (biogeochemistry), `kmf` (key metrics first), `physics`, `level1`.
       - See `analysis_timeseries` for more details, below.


A sample yaml exists in `input_yml/comparison_analysis_template.yml`,
which can be adapted to additional analyses.

Once the comparison suite has been run, members of the esmeval group workspace on JASMIN
can copy the html report to a web-visible directory, using the script:

```
./rsync_to_esmeval.sh
```

then the report will appear on the [JASMIN public facing page](https://gws-access.jasmin.ac.uk/public/esmeval/CompareReports/bgcval2/),
which is public facing but password protected.


Downloading data using MASS
===========================

Data can be downloaded and prepared for analysis using the `download_from_mass` bgcval2 tool,
with the command:
```
download_from_mass -j jobIDs
```
where `jobIDs` is one or more jobIDs.

This script will only work on JASMIN's `mass-cli1` machine,
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

Please consult the command help for more details:
```
download_from_mass -h
```


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

The command to run it is `analysis_timeseries --jobID jobID --keys suites`,
where jobID is one or more mass jobIDs, such a `u-ab123`.
The suite is one or more pre-defined suites, which is a yaml file of
list of variables.

Suite  | What it is | Description
:--------------:|:------------:|:------------:
`kmf` | Key Metrics First | A short and quick list of the most important metrics.
`physics` | Physics | A comprehensive list of physical metrics.
`bgc` | Biogeochemistry | A comprehensive list of biogeochemical metrics.
`alkalinity` | Alkalininty & CO2 | A list of alkalininty, pH, DIC, and pCO2 metrics.
`level1` | Level 1 | A comprehensive list of physical and biogeochemical metrics.
`debug` | Debug | A very short list of a couple keys to test code changes.
`fast` | UKESM1-fast  | A list of metrics tailored to the UKESM1-Fast model.

Note that there may be some overlap between the contents of these suites.

These suites are defined in the bgcval2 directory, `key_lists`.
To create a new user defined suite, the file must be present there
as a yaml and must be named with all-lower-case letters.

The contents of the suite is a list of keys. Each key in the suite
much have an associated yaml file in the `key_files` directory.
These files define how bgcval2 locates and interacts with model and
observational data.
The values in these files depend on the analysis but includes:

```
---
name: Analysis_Name
units: Units
dimensions: 1,2 or  3 # The number of dimensions after the calculations are performed
layers          : Surface
regions         : Global ignoreInlandSeas SouthernOcean ArcticOcean Equator10 Remainder NorthernSubpolarAtlantic NorthernSubpolarPacific

#paths:
modelFiles      : $BASEDIR_MODEL/$JOBID/nemo_$JOBIDo_1y_*_grid-T.nc
dataFile        : $BASEDIR_OBS/WOA/annual/woa13_decav_t00_01v2.nc
gridFile        : $PATHS_GRIDFILE

# model details
model: Model_name
model_vars      : thetao_con
model_convert   : NoChange

# Observational Data coordinates names
datasource      : WOA
data_vars       : t_an
data_convert    : NoChange
data_tdict      : ZeroToZero
```

Certain key fields in paths can be replaced with another value, allowing more
flexibility of inputs. For instance, if the path includes the username,
then `$USER` will be replaced by the users name.
These fields include:
  - `$JOBID`
  - `$USER`
  - `$YEAR`
  - `$MODEL`

In addition, some paths from the `Paths.py` can also be used:
  - `basedir_model`: The path to the model base directory, defined in Paths.paths
  - `basedir_obs`: he path to the observations base directory, defined in Paths.paths
  - `PATHS_GRIDFILE`: The path to the gridfile, defined in Paths.paths
  - `PATHs_BGCVAL2`: the path to the bgcval2 repository directory.
  
A `convert`  dictionary can be given to each model or observation data in the yml:

```
model_convert:
    function: custom_function
    path: path/to/custom/function.py
    kwarg_1: 5.
    kawrg_2: yellow
```

For instance, this example applies the function `custom_function` which is in the
`path/to/custom/function.py` file, and gives it two key word arguments.
Several example functions exists in `bgcval2/functions` which may be useful 
for how to write your own.

Simiarly, the `bgcval2/functions/standard_functions.py` contains several 
basic functions such as multiply by or add to, or `noChange`, which can all be 
called without providing the `path`, and which may have their own key word
arguments.



Clearing the Cache
------------------

In order to produce a comparison report, bgcval2 first generates the comparison images, 
then populates the report using all images in the reports path. 
This means that older images may appear in new reports, even if they were removed from the 
input yamls. If you want to "clear the cache", these images need to be deleted.

The key place to clear is set by default on jasmin to be:
```
/gws/nopw/j04/ukesm/BGC_data/$USER/bgcval2/images/TimeseriesCompare/$NAME
```
where `$USER` is your jasmin user name and `$NAME` is the name given to this analysis 
in your `input_yml` file.
This is where the comparison plots are stored.

From there, these plots are copied to
```
CompareReports2/$NAME
```
This is where the report is generated.

The third place that these plots are kept is on the public facing jasmin directory:
```
/gws/nopw/j04/esmeval/public/CompareReports/bgcval2/$USER/$NAME
```
This is where the report is hosted.


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
