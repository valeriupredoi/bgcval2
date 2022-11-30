# Run this script on mass-cli1.jasmin.ac.uk
# from login1.jasmin.ac.uk, ssh to the mass machine:
#     ssh -X  mass-cli
# run script with:
# source /home/valeriu/bgcval2/mass_scripts/u-cp416debug.sh
# moo passwd -r # if mass password is expired
moo get --fill-gaps moose:/crum/u-cp416debug/ony.nc.file/*.nc /home/valeriu/bgcval2/local_test/BGC_data/u-cp416debug/
