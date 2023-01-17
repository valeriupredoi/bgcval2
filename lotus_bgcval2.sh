#!/bin/bash
#BSUB -q short-serial
#BSUB -n 1
#BSUB -W 24:00
#BSUB -o "logs/bgcval2%J.out"
#BSUB -e "logs/bgcval2%J.err"
#BSUB -J bgcval2
#BSUB -u EMAIL@email.com    
#BSUB -B
#BSUB -N


######################
# First make a copy of this file:
# rysnc -av lotus_bgcval2.sh lotus_bgcval2_$USER.sh
#
# Then edit your new copy (lotus_bgcval2_$USER.sh):
#  1. Change your email above (#BSUB -u email@email.com)
#
#  2. Check the name of your conda environment: CONDA_ENV 
#     "bgcval2" is the default, but may have a different one.
#
#  3. Change the path to your bgcval2 directory: BGCVAL2_PATH
#
#  4. Change the analyses to run: BGCVAL2_SUITE
#     You can add several suites here if needed. 
#
#  5. Submit this script with with:
#     bsub <  lotus_bgcval2_$USER.sh
#
#  6. Monitor progress with:
#     bjobs | grep $USER
#########################


# Edit these bits:

# Change this to your bgcval2 conda environment name
CONDA_ENV=bgcval2
# Change this to your bgcval2 directory - or wheever you wish to run.
BGCVAL2_PATH=/home/users/$USER/bgcval2
# Add one or more input_yml files here
BGCVAL2_SUITE=input_yml/my_bgcval2_suite.yml
#BGCVAL2_SUITE2=input_yml/my_bgcval2_suite2.yml


##########################
# script

# load the conda program:
module load jaspy

# load your conda env:
conda activate $CONDA_ENV

# Change directory to your bgcval2 directory:
cd $BGCVAL2_PATH



#########################
# Add one or more input_yml files here 
analysis_compare -y $BGCVAL2_SUITE
#analysis_compare -y $BGCVAL2_SUITE2



#########################
# Rsync report to web facing directory:
# (shouldn't need to change this!)
Outpath=/gws/nopw/j04/esmeval/public/CompareReports/bgcval2/$USER
mkdir -p $Outpath
rsync -av CompareReports2/* $Outpath/.


