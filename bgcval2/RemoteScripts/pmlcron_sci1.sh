#!/bin/bash
source /users/modellers/ledm/.bashrc
. ~/.keychain/$HOSTNAME-sh

export XAUTHORITY=/users/modellers/ledm/.Xauthority
export DISPLAY=':0'

#####
# parsing job id
jobid=${1:-u-ad980}
echo jobid=$jobid
export jobid=$jobid

# ssh -A -X ldemora@jasmin-login1.ceda.ac.uk "/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/helloOnMass.sh $jobid"

# ssh -A -X ldemora@jasmin-login1.ceda.ac.uk "/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/runWholePackageOnSci1.sh $jobid"
ssh -A -X ldemora@login1.jasmin.ac.uk "/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/runWholePackageOnSci2.sh $jobid"

