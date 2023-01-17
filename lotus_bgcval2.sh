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
# What to edit:
#  1. Change your email above (#BSUB -u email@email.com)
#
#  2. Check the name of your conda environment. 
#     "bgcval2" is the default, but may have a different one.
#
#  3. Change the path to your bgcval2 directory.
#
#  4. Change the analyses to run below.
#     You can add several suites here if needed. 
#
#  5. Submit this script with with:
#     bsub <  lotus_bgcval2_$USER.sh
#
#  6. Monitor progress with:
#     bjobs | grep $USER
#########################

# load the conda program:
module load jaspy



# Change this to your bgcval2 directory - or wheever you wish to run.
cd /home/users/$USER/bgcval2

# Change this to your bgcval2 conda environment name
conda activate bgcval2


#########################
# Add one or more input_yml files here 
analysis_compare -y input_yml/my_bgcval2_suite.yml



#########################
# Rsync report to web facing directory:
# (shouldn't need to change this!)
Outpath=/gws/nopw/j04/esmeval/public/CompareReports/bgcval2/$USER
mkdir -p $Outpath
rsync -av CompareReports2/* $Outpath/.


