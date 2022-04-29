!/bin/bash
source /users/modellers/ledm/.bashrc
. ~/.keychain/$HOSTNAME-sh

export XAUTHORITY=/users/modellers/ledm/.Xauthority
export DISPLAY=':0'
#ssh -A -X ldemora@jasmin-login1.ceda.ac.uk '/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/helloOnMass.sh'

#####
# parsing job id
#jobid=${1:-compare}
#echo jobid=$jobid
#export jobid=$jobid

#ssh -A -X ldemora@jasmin-login1.ceda.ac.uk "/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/helloOnMass.sh $jobid"

ssh -A -X ldemora@jasmin-login1.ceda.ac.uk "/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/runCompareOnSci1.sh"

