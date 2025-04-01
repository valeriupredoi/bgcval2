#!/bin/bash
#SBATCH --job-name=bgcval2
#SBATCH --partition=standard
#SBATCH --time 6:00:00
#SBATCH --qos=standard
#SBATCH --account=esmeval
#SBATCH -o logs/log_bgcval2_ts_%J.out
#SBATCH -e logs/log_bgcval3_ts_%J.err

# Note that the 
# slurm job name should be set in the command line with -J option to the jobID
# and the output and error filers are also set at the command line.
#

######################
# This script runs a single time series job.
#
# Runs:
#    sbatch -J jobID lotus_timeseries.sh jobID suite1 suite2 etc...
# Note that that batch_timeseries command also
# adds specific out and job scripts,
# and outputs them to ./logs directory..  
#########################


#########################
# Change this to your bgcval2 conda environment name
CONDA_ENV=bgcval2


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
# Load command line arguments:
args=$@
jobID=$1
suites="${@:2}"
echo $suites

##########################
# Run single jog timeseries analysis.
echo "analysis_timeseries -j $jobID -k $suites"
analysis_timeseries -j ${jobID} -k ${suites}


