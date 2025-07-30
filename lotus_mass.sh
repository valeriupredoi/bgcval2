#!/bin/bash
#SBATCH --job-name=massbgcval
#SBATCH --partition=mass
#SBATCH --time 24:00:00
#SBATCH --qos=mass
#SBATCH --account=mass
#SBATCH -o /home/users/ldemora/mass2_logs/log_mass2.out
#SBATCH -e /home/users/ldemora/mass2_logs/log_mass2.err

# Submit this script to mass-cli2.jasmin.ac.uk
#  sbatch lotus_mass.sh

# run script with:
# moo passwd -r # if mass password is expired


##########################
# Run mass download.

# Set your target directory
TARGET_DIR="/gws/nopw/j04/esmeval/bgcval2/shared_mass_scripts/"

# Create an array of filenames without extensions
file_list=()
for file in "$TARGET_DIR"*; do
  basename=$(basename "$file")
  name="${basename%.*}"  # Remove extension
  file_list+=("$name")
done

# Iterate through the list and echo the file names
for fname in "${file_list[@]}"; do
  echo "Filename: $fname"
done


# Shuffle the list and iterate through it
for JOBID in $(printf "%s\n" "${file_list[@]}" | shuf); do
  mkdir -p /gws/nopw/j04/ukesm/BGC_data/$JOBID/
  moo get --fill-gaps moose:/crum/$JOBID/ony.nc.file/*.nc /gws/nopw/j04/ukesm/BGC_data/$JOBID/
done



