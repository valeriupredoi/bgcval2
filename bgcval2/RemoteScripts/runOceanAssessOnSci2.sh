#!/bin/bash
echo "uname -a"
uname -a 
echo "USER: " $USER


echo 'SSH_AUTH_SOCK:' $SSH_AUTH_SOCK  
echo 'SSH_CLIENT:' $SSH_CLIENT 
echo 'SSH_CONNECTION:' $SSH_CONNECTION 
echo 'SSH_TTY:'  $SSH_TTY   

#####
# parsing job id
jobid=${1:-u-ad980}
echo runOceanAsssessOnSci2 jobid=$jobid
export jobid=$jobid

python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py $jobid

ssh -X -A jasmin-sci2 "cd /home/users/ldemora/workspace/ocean_assess; ./ocean_assess_jobid.sh $jobid" 	# supposed to be jobid in script filename

echo "The end of runOceanAssessOnSci2.sh"
