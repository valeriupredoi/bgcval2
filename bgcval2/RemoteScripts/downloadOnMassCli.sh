#!/bin/bash
echo "uname -a"
uname -a 
echo "USER: " $USER
#eval "$(ssh-agent)"

echo 'SSH_AUTH_SOCK:' $SSH_AUTH_SOCK  
echo 'SSH_CLIENT:' $SSH_CLIENT 
echo 'SSH_CONNECTION:' $SSH_CONNECTION 
echo 'SSH_TTY:'  $SSH_TTY        

#####
# parsing job id
jobid=${1:-u-ad980}
echo jobid=$jobid
export jobid=$jobid

 
python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py $jobid


#printf "\n\nssh -A -X mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'\n"
#ssh -X -A mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'

ssh -X -A mass-cli.jasmin.ac.uk "python /home/users/ldemora/workspace/ukesm-validation/bgcvaltools/downloadFromMass.py $jobid"
#ssh -X -A mass-cli.jasmin.ac.uk "mkdir -p gws/nopw/j04/ukesm/BGC_data/u/group_workspaces/jasmin2/ukesm/BGC_data/$jobid/1y; ln -s  /group_workspaces/jasmin2/ukesm/BGC_data/$jobid/*1y*grid_[UVWT]*.nc /group_workspaces/jasmin2/ukesm/BGC_data/$jobid/1y/.; ls -lhrt /group_workspaces/jasmin2/ukesm/BGC_data/$jobid"
ssh -X -A mass-cli.jasmin.ac.uk "mkdir -p /gws/nopw/j04/ukesm/BGC_data/$jobid/1y; ln -s  /gws/nopw/j04/ukesm/BGC_data/$jobid/*1y*grid_[UVWT]*.nc /gws/nopw/j04/ukesm/BGC_data/$jobid/1y/.; ls -lhrt /gws/nopw/j04/ukesm/BGC_data/$jobid"

#####
# MLD:
#ssh -X -A mass-cli1 "python /home/users/ldemora/workspace/ukesm-validation/bgcvaltools/downloadFromMass.py $jobid mld"

echo "The end of downloadOnMassCli.sh $jobid"

