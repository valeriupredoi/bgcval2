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
