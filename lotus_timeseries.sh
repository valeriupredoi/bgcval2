#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --time 6:00:00
#SBATCH -o logs/log_bgcval2_ts_%J.out
#SBATCH -e logs/log_bgcval3_ts_%J.err

# slurm job name should be set in the command line with -J option to the jobID

######################
# This script runs a single time series job.
#
# Runs:
#    sbatch -J jobID lotus_timeseries.sh jobID suite1 suite2 etc...
#
#########################


#########################
# Change this to your bgcval2 conda environment name
CONDA_ENV=bgcval2

# Change this to your bgcval2 directory - or wheever you wish to run.
#BGCVAL2_PATH=$PWD

# Add one or more input_yml files here
#BGCVAL2_SUITE2=input_yml/my_bgcval2_suite2.yml

##########################
# Source global definitions
if [ -f ~/.bashrc ]; then
    echo 'source ~/.bashrc'
    source ~/.bashrc
fi

##########################
# load your conda env:
#conda activate bgcval2
echo conda activate $CONDA_ENV
conda activate ${CONDA_ENV}

##########################
# Change directory to your bgcval2 directory:
#echo cd $BGCVAL2_PATH
#cd $BGCVAL2_PATH
#pwd

##########################
# Load command line arguments:
args=$@
jobID=$1
suites="${@:2}"
echo $suites

echo "analysis_timeseries -j $jobID -k $suites"
analysis_timeseries -j ${jobID} -k ${suites}
#########################
# Add one or more input_yml files here
#cho analysis_timeseries -y $BGCVAL2_SUITE
#nalysis_timeseries -y $BGCVAL2_SUITE 
