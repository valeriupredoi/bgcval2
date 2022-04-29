#!/bin/bash
echo "uname -a"
uname -a 
echo "USER: " $USER
#eval "$(ssh-agent)"

echo 'SSH_AUTH_SOCK:' $SSH_AUTH_SOCK  
echo 'SSH_CLIENT:' $SSH_CLIENT 
echo 'SSH_CONNECTION:' $SSH_CONNECTION 
echo 'SSH_TTY:'  $SSH_TTY         

#rintf "\n\nssh-add -l\n"
#sh-add -l

#####
# parsing job id
jobid=${1:-u-ad980}
echo jobid=$jobid
export jobid=$jobid

python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py $jobid

#printf "\n\nssh -v -X -i /home/users/ldemora/.ssh/id_rsa ldemora@mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'\n"
#ssh -v -X -A -i /home/users/ldemora/.ssh/id_rsa ldemora@mass-cli1 'source /home/users/ldemora/.bashrc; python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'

#printf "\n\nssh -v -X ldemora@mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'\n"
#ssh -v -X -A ldemora@mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'

printf "\n\nssh -A -X mass-cli1 'python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py'\n"
ssh -X -A mass-cli1 "python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py $jobid"

ssh -X -A jasmin-sci1 "python /home/users/ldemora/workspace/ukesm-validation/RemoteScripts/hello.py $jobid"




