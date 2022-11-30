#!/bin/bash

# single script to automatically download bgcval2 users data.
# 
# cron tab line:
# 00 19 * * * /gws/nopw/j04/esmeval/bgcval2/shared_mass_cron.sh >  /gws/nopw/j04/esmeval/bgcval2/shared_mass_cron.out 2>&1 # bgcval2 automatic download 


#self-destruct:
{
    sleep 1380m  # Kill script after 23 hours.
    echo "Script ran out of time (23 hour limit)"
    kill $$
} &


# delete files older than 2 weeks
find /gws/nopw/j04/esmeval/bgcval2/shared_mass_scripts -mindepth 1 -mtime +15 -print -delete

# run all files:
find /gws/nopw/j04/esmeval/bgcval2/shared_mass_scripts/*.sh -exec bash {} \;

