#!/bin/bash
#SBATCH --job-name=$1
#SBATCH --partition=short-serial
#SBATCH --time 6:00:00
#SBATCH -o logs/log_bgcval2_ts_%J.out
#SBATCH -e logs/log_bgcval3_ts_%J.err


######################
# First make a copy of this file:
# rysnc -av lotus_bgcval2.sh lotus_bgcval2_$USER.sh
#
# Then edit your new copy (lotus_bgcval2_$USER.sh):
#
#  1. Check the name of your conda environment: CONDA_ENV
#     "bgcval2" is the default, but may have a different one.
#
#  2. Check the path to your bgcval2 directory: BGCVAL2_PATH
#
#  3. Change the analyses to run: BGCVAL2_SUITE
#     You can add several suites here if needed.
#
#  4. Submit this script with with:
#     sbatch lotus_bgcval2_$USER.sh
#
#  5. Monitor progress with:
#     squeue | grep $USER
#
#########################



#########################
# Edit these bits:

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

jobID=$1
keys=$2

echo "analysis timeseries -j $jobID -k $2"
#analysis_timeseries -j $jobID -k $2
#########################
# Add one or more input_yml files here
#cho analysis_timeseries -y $BGCVAL2_SUITE
#nalysis_timeseries -y $BGCVAL2_SUITE 
