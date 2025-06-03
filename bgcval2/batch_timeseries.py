#!/usr/bin/env python
#
# Copyright 2015, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license.

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with bgc-val.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk
#
"""
.. module:: batch_timeseries
   :platform: Unix
   :synopsis: A script to submit slurm scripts time series.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
import argparse
import subprocess
import os
import sys

from getpass import getuser

from bgcval2.analysis_compare import load_comparison_yml


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-y',
                        '--compare_yml',
                        nargs='+',
                        type=str,
                        help='One or more Comparison Analysis configuration file, for examples see bgcval2 input_yml directory.',
                        required=True,
                        )

    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'default-bgcval2-config.yml'),
                        help='User configuration file (for paths).',
                        required=False)

    parser.add_argument('--dry_run',
                        '-d',
                        default=False,
                        help='When True: Do not submit the jobs to lotus.',
                        action=argparse.BooleanOptionalAction,
                        required=False)

    args = parser.parse_args()
    return args


def submits_lotus(compare_yml, config_user, dry_run=False):
    """
     Loads the yaml file and submits individual time series to sbatch.
     """
    # Load details from yml file
    details = load_comparison_yml(compare_yml)

    # list of job IDS
    jobs = details['jobs']

    # username
    user = getuser()

    # Load current on-going list of this users slurm jobs:
    out = str(subprocess.check_output(["squeue", "--user="+user]))

    # loop over jobs:
    for job in jobs:
        # Check whether there's already a job running for this jobID
        if out.find(job) > -1:
            print("That job exists already: skipping", job)
            continue

        # Get list of suites for each job
        suites = details['suites'][job]

        # Make it a list:
        if isinstance(suites, str):
            suites = suites.split(' ')

        # prepare the command
        command_txt = ['sbatch', 
             '-J', job, 
              ''.join(['--error=logs/', job,'.err']),
              ''.join(['--output=logs/', job,'.out']),
             'lotus_timeseries.sh', job]
        for suite in suites:
            command_txt.append(suite)

        # Send it!
        if dry_run:
            print('Not submitting (dry-run):', ' '.join(command_txt))
        else:
            # Submit job:
            print('Submitting:', ' '.join(command_txt))
            #command1 = subprocess.Popen(command_txt)
            command1 = subprocess.Popen(
                  command_txt,
                  stdout=subprocess.PIPE,
                  stderr=subprocess.STDOUT,
            )
      

def main():

    """Run the main routine."""
    args = get_args()

    # This has a sensible default value.
    config_user=args.config_file

    # This shouldn't fail as it's a required argument.
    compare_ymls = args.compare_yml

    for compare_yml in compare_ymls:
        print(f"analysis_timeseries: Comparison config file {compare_yml}")

        if not os.path.isfile(compare_yml):
            print(f"analysis_timeseries: Could not find comparison config file {compare_yml}")
            sys.exit()
        dry_run = args.dry_run
        submits_lotus(compare_yml, config_user, dry_run)


if __name__ == "__main__":
    from ._version import __version__
    print(f'BGCVal2: {__version__}')
    main()

