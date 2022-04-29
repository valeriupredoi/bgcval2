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
#jobid=${1:-u-ad980}
#echo jobid=$jobid
#export jobid=$jobid

#python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py $jobid
#ssh -X -A jasmin-sci1 "rsync -av /home/users/ldemora/workspace/ukesm-validation/images/TimeseriesCompare/* /group_workspaces/jasmin4/esmeval/public/Intercomparison/."

ssh -X -A jasmin-sci1 "cd /home/users/ldemora/workspace/ukesm-validation; pwd; python2.7 /home/users/ldemora/workspace/ukesm-validation/analysis_compare.py; rsync -av /home/users/ldemora/workspace/ukesm-validation/CompareReports /group_workspaces/jasmin4/esmeval/public/."
echo "The end of runCompareOnSci1.sh on login lander"
