#!/bin/bash
source /users/modellers/ledm/.bashrc
. ~/.keychain/$HOSTNAME-sh

export XAUTHORITY=/users/modellers/ledm/.Xauthority
export DISPLAY=':0'
#ssh -A -X ldemora@jasmin-login1.ceda.ac.uk '/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/helloOnMass.sh'

ssh -A -X ldemora@jasmin-login1.ceda.ac.uk '/home/users/ldemora/workspace/ukesm-validation/RemoteScripts/runReportOnSci1.sh'
cd /users/modellers/ledm/ImagesFromJasmin
tar xfvz report-u-ad980.tar.gz
