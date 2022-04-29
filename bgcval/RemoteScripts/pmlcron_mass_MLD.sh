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


#eval `ssh-agent -s` 
#ssh -A -X ldemora@jasmin-login1.ceda.ac.uk '/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/helloOnMass.sh' 
ssh -A -X ldemora@jasmin-login1.ceda.ac.uk "/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/downloadOnMassCli_MLD.sh $jobid"
#kill $SSH_AGENT_PID

