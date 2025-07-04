[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Github Actions Test](https://github.com/valeriupredoi/bgcval2/actions/workflows/run-tests.yml/badge.svg)](https://github.com/valeriupredoi/bgcval2/actions/workflows/run-tests.yml)

![bgcval2logo](https://github.com/valeriupredoi/bgcval2/blob/main/doc/figures/bgcval2_logo_v_small.png)

bgcval2
=======


bgcval2 is the updated and modernised version of BGC-val.
There are several updates over the previously published version of [BGCVal](https://gmd.copernicus.org/articles/11/4215/2018/).
The primary improvements are to the user interface, the ease of use, ease of installation,
and a general modernisation of the core approach, including the move to python3, 
and using a conda enviroment.

This work was funded through WP1 of the Terrafirma project.


Current version notes:

- Suported versions of Python: 3.10, 3.11, 3.12, and 3.13
- Numpy 2+ compatible

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

- then install the development dependencies and the tool itself
  (note that all our dependencies are from `conda-forge` so there
  is no risk of mixing Conda and PyPI packages):

```
pip install -e .[develop]
```

Test that the tool has been installed correctly with:
```
analysis_compare -h
```
which should print the module information and instructions on how to run the tool.


On the JASMIN computational system (HPC), members of the `esmeval` working group should
be able to run the debug analysis script:
```
analysis_compare -y input_yml/debug.yml
```
This script performs an analysis of two small physics-only UKESM development jobIDs,
using the debug suite.


Keeping the code up to date
---------------------------

To keep the code up to date with the main branch, fetch the changes from github:
```
git fetch
```

You can merge the main branch into your local code:
```
git merge origin/main
```

If the `environment.yml` file has changed in the merge, you can update your conda environment with:
```
mamba env update -n bgcval2 -f environment.yml
pip install -e .[develop]
```

or delete your old environment and follow the instructions above to create a new one:
```
conda env remove -n ENV_NAME
```


### Available executables

Executable name | What it does | Command
:--------------:|:------------:|:------------:
`analysis_timeseries` | runs timeseries analysis for single model run. | analysis_timeseries -j jobID -k key
`analysis_p2p` | runs point to point analysis of model against observational dataset. | analysis_p2p jobID YEAR
`bgcval` | runs time series and point to point. | bgcval jobID
`bgcval2_make_report` | makes the single model HTML report. | bgcval2_make_report jobID
`analysis_compare` | runs comparison of multiple single jobs  | analysis_compare
`batch_timeseries` | Submits single job time series analysis to slurm | batch_timeseries


### Checking out development branches

Should you want to work on several git development branches, you can switch between them, making sure they always
stay up-to-date with the remote branch; right after you have checked out the `bgcval2` git repository (see above)
you find yourself on the `main` branch (the main development branch); after a while, with new updates upstream (at remote),
this branch (your local `main`) becomes out of date with remote, so you should always update it:

```
git pull origin main
```

After you have updated it, you can now check out a new branch `devel_branch` from upstream, some branch that was recently pushed and you need to test or contribute to:

```
git fetch -v
git checkout devel_branch
git pull origin devel_branch
```

If you want to create your own branch, first make sure you are happy with your current code. A good tip is to start from the main branch

```
git checkout origin/main
git branch new_branch_name
git checkout new_branch_name
```

### Updating an existing environment

Should you wish to update your `bgcval2` environment (e.g. when dependencies have changed), we recommend
deleting the old environment, and creating a brand new one from scratch. To delete the old environment just remove
it from the conda tree, after you have deactivated it first:

```
conda deactivate
rm -rf $USER/miniconda3/envs/bgcval2
```

Now you can follow the steps above to (re-)create the new environment, still called `bgcval2`, and
installing the tool inside it:

```
mamba env create -n bgcval2 -f environment.yml
conda activate bgcval2
pip install -e .[develop]
```

### Creating multiple environments from different branches

Should you need to have multiple environments, each with a different installation of bgcval2 (e.g when working
with different development branches), you can simply create them like instructed above, with the only difference
that the environments must have relevant and differing names. An example: you are working on a branch called
`bgcval2_some-feature` and would like to test it in a bespoke environment, you can create it:

```
mamba env create -n bgcval2_some-feature -f environment.yml
conda activate bgcval2_some-feature
```

then pip-install the tool there, and on you go with testing. This environment exists in parallel with the main feature `bgcval2`
environment, and you can toggle between them with `conda deactivate` current env, then `conda activate` any other environment.


> **_NOTE:_** You don't need to create new environments in parallel if you just want to test a new branch;
after checking out the new branch you can test it immediately, while still in the `bgcval2` environment,
since `pip install .` is a local installation that uses the local `bgcval2/` repository (it doesn't move
installed scripts in such locations as `/opt` or `/lib`)

### Saving and committing local changes

If you made changes to the code, and you want to bring in the latest remote changes (via `git pull`), you have two options:
- either commit your changes to the remote branch and then push your branch back to github (origin):

```
 git commit your-changed-code -m "commit message"
 git push origin your_branch
```

- or, save your changes elsewhere (in a directory that is not under bgcval2's git control), and run `git stash` to bring your
  repository to the state of the latest pull from remote (at the HEAD of the latest pushed commit)

Either way you go, please make sure your changes are either committed, or saved and stashed for a later commit!

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
strictFileCheck: <bool>
dpi: <int>
savepdf: <bool>
savejson: <bool>

jobs:
   <jobID1>:
      description: <descrption of the first job>
      label: 'Job info'
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
   -  This value is overwritten by the `analysis_compare -s` which skips new timeseries analyses.
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
 - `strictFileCheck`: 
   - A boolean which when True will raise an error if input model files are missing. 
   - Default is True, set to False to skip this check.
 - `dpi`:
   - The resolution in dots per inch (dpi) of the output image.
   - dpi=100 is fine for most purposes but 300 is needed for highres posters/publications.
   - If set to `None`, then the dpi will be set to your default value.
 - `savepdf`:
   - Output the image as a pdf in addition to the standard png. 
   - This doesn't replace the image in the report, but saves a pdf version of the image in the images directory.
   - The pdfs will be web-visible in the report image directory, but will not linked in the html.
   - To view a pdf from a report, simply click the image, and replace the `png` extension with `pdf` in the browser path.
 - `safejson`:
   - Outputs the data and plotting details of each timeseries plot as a json file.
   - json is a human readable format that is similar to a python dictionary or a shelve file. 
   - This file includes all data, time values, units, and line colors, thicknesses and styles.

 - `jobs`:
   - A list of jobIDs, and some options on how they will appear in the final report.
   - The options are:
     - `description`:
       - A description of job, which helps people differentiate the jobs in the final report.
     - `label`:
       - A short description of the job which will appear in the plot legend.
       - If you make it too long, it won't fit on the plot.  
       - Optional: Default behaviour is the jobID.
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

In order to only run the report making part of the comparison analysis 
(skip the `analysis_timeseries` part),
either set the `do_analysis_timeseries` key to `False` in the `input_yml` file,
or run `analysis_compare` with the command line argument: `-s` or `--skip-timeseries`.
To skip the analysis timeseries command, use `--no-skip-timeseries`.
Though without without either of these command line arguments, 
bgcval2 will default to the value in your `input_yml` file.
Also, the command line argument overwrites the value in `input_yml`.

Once the comparison suite has been run, members of the esmeval group workspace on JASMIN
can copy the html report to a web-visible directory, using the script:

```
./rsync_to_esmeval.sh
```

then the report will appear on the [JASMIN public facing page](https://gws-access.jasmin.ac.uk/public/esmeval/CompareReports/bgcval2/),
which is public facing but password protected.


Information about the json output:
----------------------------------

The `savejson` flag allows bgcval2 to output all data from the multi-model plots to json. 
Json files are effectively python dictionaries that are saved as human readable ASCII files. 
This means that they can contain a lot more metadata than a simple csv file. 
It also means that it's easier to deal with multiple datasets with different time ranges. 

So, in these files, everything is saved as a series of dictionaries, and the main ones are:
 - `timesD`: The time value (in decimal years)
 - `arrD`: The data value. 

There's also a lot of info used to make the bgcavl2 plots, including 
 - `colours`: the colour used by bgcval2 for each jobid
 - `linestyle`: the plot line style for each job id
 - `label`: the plot lable for each jobid

Each of these dictionaries uses the jobID as the index, so in python, these can be plotted using the code:


```
import json
from matplotlib import pyplot

json_file = 'filename.json'
with open(json_file) as json_file:
    amoc_data = json.load(json_file)

for jobID in sorted(amoc_data['timesD'].keys()):
    times = amoc_data['timesD'][jobID]
    amoc = amoc_data['arrD'][jobID]
    colour = amoc_data['colours'][jobID]
    pyplot.plot(times, amoc, c=colour, label=jobID)

pyplot.show()
```




Batch times series Analysis
===========================

The `batch_timeseries` tool can take an `analysis_compare` input yaml file,
and instead of running the time series analysis for each job on
the interactive shell terminal in series, it uses slurm to submit
each job as an independent job. 

On jasmin, users can run up to five jobs simulataneously,
so this can singnificantly boost the speed of the analysis. 

The command to run it is:
```
batch_timeseries - y comparison_recipe.yml
```

This will submit a time-series analysis for each job, using a command which looks like this:
```
sbatch -J jobID --error=logs/jobID .err --output=logs/jobID .out lotus_timeseries.sh jobID  kmf physics bgc
```
The output and error messages will be in the `logs` directory with the jobID as the file prefix.
The job name on slurm will also be the jobID, so it's easy to tell which jobs are running.
The analysis suites will be appended as a list to the end of the command.
In order to reduce the chance of analysing the same jobID twice, `batch_timeseries`
checks whether a job exists, either currently running or in the queue before submitting.
If a jobID exists, it is not re-submitted. However, this means that
if two versions of the same jobID are submitted one after the other
with different suite lists (`kmf`, `physics`, `bgc`), then only the first 
set of suites will be run. 

There is also an optional flag `-d` or `--dry_run` to test `batch_timeseries`, 
which outputs the submission command to screen but does not submit the jobs.

Note that this task does not run the `analysis_compare` suite so it will 
not generate the html report. However, the html report can be generated more quickly
with the `-s` argument to skip the `analysis_timeseries` section 
described above. 

In addition, note that this will not run the `download_from_mass`
script, so jobs added here will not be included in the automated download.
However, these jobs are added for automated download when `analysis_compare` 
is used. 


Downloading data using MASS
===========================

UKESM jobs on JASMIN automatically download data from Mass. 
The `analysis_compare` tool creates a script in the `mass_scripts` directory.
This script is copied to a shared directory in the shared ukesm/terrafirma diskspace:

```
/gws/nopw/j04/esmeval/bgcval2/shared_mass_scripts
```

A crontab job is set up to execute these scripts overnight on mass-cli1,
but only executes scripts that are younger than 30 days. 
This prevents older jobs from being run after they are no longer needed.
Very old scripts are automatically deleted. 
The crontab script is:

```
/gws/nopw/j04/esmeval/bgcval2/shared_mass_cron.sh
```

Note that cron executes this script with the `ldemora` username, so permissions may need to be changed
in some cases. The output log for this script is sent to the file:

```
/gws/nopw/j04/esmeval/bgcval2/shared_mass_cron.out
```

In the case that data is needed immediately, it's posiible to download manually as well.
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
smoothings      : DataOnly both5 both30 movingav30years 5and30 30and100


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
    kwarg_2: yellow
```

For instance, this example applies the function `custom_function` which is in the
`path/to/custom/function.py` file, and gives it two key word arguments.
Several example functions exists in `bgcval2/functions` which may be useful
for how to write your own.

Simiarly, the `bgcval2/functions/standard_functions.py` contains several
basic functions such as multiply by or add to, or `noChange`, which can all be
called without providing the `path`, and which may have their own key word
arguments.

Some of the other options for each analysis include:
  - `layers`:
    - The z axis layers that you want to use.
    - Typically include things like `Surface`, `500m`, `1000m` etc.
    - Note that this is a pre-curated list, so you can't use it to select a specific depth (like `327m` or something)
  - `regions`:
    - A list of regions to investigate.
    - This is a pre-curated list, and they are defined in so you can't.
    - These regions are defined in the file `bgcval2/bgcvaltools/makeMask.py`.
  -  `smoothings`:
    - This is the smoothing function (if any) to apply to the data before plotting it.
    - The smoothing it added before plotting, but after saving the shelve file, so it doesn't impact the data.
    - No smoothing is `DataOnly`, which is also the default behaviour. 
    - If several smoothing options are added here, bgcval2 will generate a comparison plot for each one.
    - The smoothings are defined and performed in `bgcval2/timeseries/timeseriesPlots.py`
    - Other modes exist:
      - `movingav30years`: A 30 year moving average
      - `both5`: Both no smoothing and a 5 year moving average
      - `both30`: Both no smoothing and a 30 year moving average
      - `5and30`: a 5 year moving average and a 30 year moving average.

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

**WORK IN PROGRESS**

This is an analysis that compares a specific year of the model
against several climatological datasets, and produces comparison
maps, histograms and scatter plots.

The command to run it is `analysis_p2p -j JOBID -y YEAR -k key_list`,
where jobID is a unique model run identifier, such a `u-ab123`.
The  `YEAR` is the single year to investigate
and the key is a pre-defined key, which generates a list of variables.

Key | What it is | Description
:--------------:|:------------:|:------------:
`physics_p2p` | Physics | A short list of physical metrics for p2p analysis.

Single Model report
-------------------

**WORK IN PROGRESS**

Once an analysis has run, either time series or point to point, a report
can be generated from this output, using the command:

```
bgcval2_make_report jobID year
```

This generates an HTML5 mobile-friendlyt report, summarising the output of a
single model run.

Running a job on JASMIN's batch system, LOTUS
=============================================

Instead of running these jobs interactively on the command line, 
it is possible to submit them as batch jobs to JASMIN's LOTUS machine.
[See the LOTUS documentation here](https://help.jasmin.ac.uk/article/5004-lotus-overview).

To submit a job, make a copy of the `lotus_bgcval2.sh` file:

```
rsync -av lotus_bgcval2.sh lotus_bgcval2_$USER.sh
```

From there, you'll need to edit your copy of the file `lotus_bgcval2_$USER.sh`.

In this file, you'll need to edit the following bash objects:

1. `CONDA_ENV`: The bgcavl2 conda environment (default is bgcval2)
2. `BGCVAL2_PATH`: The path to your bgcval2 directory (default is ~/bgcval2)
3. `BGCVAL2_SUITE`: Your yml comparison suite file. (no default)

You may also need to edit the job time limit, at the top of the file.
The default it set to 6 hours.

Once that is done, the script is submitted to LOTUS with:

```
sbatch  lotus_bgcval2_$USER.sh
```

You can monitor your script in the queue with:

```
squeue | grep $USER
```

You will also see a couple log files appear which will allow you to follow
the output of the process.
If you wish to send to to a specific path, you may edit the lines at the top of 
`lotus_bgcval2_$USER.sh`:

```
#SBATCH -o log_bgcval2_%J.out
#SBATCH -e log_bgcval3_%J.err
```

Documentation
=============

**WORK IN PROGRESS**

See available Sphinx [documentation](https://htmlpreview.github.io/?https://github.com/valeriupredoi/bgcval2/blob/main/doc/build/index.html). To build locally the documentation run:

```
sphinx-build -Ea doc doc/build
```

Appendix
========

## Python2 to Python3 Migration

Migration from the orginal BGCVal code, which was Python2, has been done with the `2to3` tool:

- Install `2to3` package:

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
