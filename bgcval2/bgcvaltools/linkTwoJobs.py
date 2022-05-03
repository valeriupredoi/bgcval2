#!/usr/bin/python

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
"""
.. module:: linkTwoJobs
   :platform: Unix
   :synopsis: Links two jobs, used when a job fails and is restarted with a new jobID.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from sys import argv
import subprocess
from socket import gethostname
import os
import shutil
from glob import glob
from shelve import open as shopen
from datetime import datetime

from ..UKESMpython import folder
from getpass import getuser


def linkTwoJobs(jobID1, jobID2, copyShelves=False):
    """
		command:
		linkTwoJobs jobID1 jobID2.
			job 1 is the source job
			job 2 is the new job, where the links will be put.
		
	"""

    if '' in [
            jobID1,
            jobID2,
    ]: return

    ########
    # Two aspects of this:
    #	link the netcdf files
    #	copy the shelve files.

    linkNetcdfs = True
    #copyShelves = False
    linkchain = True

    #####
    # Netcdf stuff:
    if linkNetcdfs:
        machine = gethostname()
        knownmachine = False
        if machine.find('ceda') > -1:
            knownmachine = True

            netcdfFold1 = "/group_workspaces/jasmin2/ukesm/BGC_data/" + jobID1
            netcdfFold2 = netcdfFold1.replace(jobID1, jobID2)
            if not os.path.exists(netcdfFold1):
                print(
                    "This folder doesn't exist. Do you have the right jobID 1?",
                    jobID1)
                print(netcdfFold1)
                return
            if not os.path.exists(netcdfFold2):
                print("This folder doesn't exist yet.")
                print("making", folder(netcdfFold2))
                #return
            shelvefold1 = "/group_workspaces/jasmin2/ukesm/BGC_data/" + getuser(
            ) + "/shelves/*/" + jobID1
            shelvefold2 = "/group_workspaces/jasmin2/ukesm/BGC_data/" + getuser(
            ) + "/shelves/*/" + jobID2

        if machine.find('monsoon') > -1:
            knownmachine = True
            netcdfFold1 = "/projects/ukesm/ldmora/UKESM/" + jobID1
            netcdfFold2 = netcdfFold1.replace(jobID1, jobID2)
            if not os.path.exists(netcdfFold1):
                print(
                    "This folder doesn't exist. Do you have the right jobID 1?",
                    jobID1)
                return
            if not os.path.exists(netcdfFold2):
                print(
                    "This folder doesn't exist. Do you have the right jobID 2?",
                    jobID2)

                return

            shelvefold1 = "shelves/timeseries/" + jobID1
            shelvefold2 = folder(shelvefold1.replace(jobID1, jobID2))

        if not knownmachine:
            print("Are you running this on the correct machine?")
            print(
                "\tYou should be on mass-cli1.ceda.ac.uk at jasmin or on monsoon at the MO"
            )
            print("\tBut you're at", machine)
            return

        #####
        # link the files:
        physicsonly = 0  #True
        bgconly = 1  #False
        for fn1 in sorted(glob(netcdfFold1 + '/*' + jobID1 + '*')):
            print("linking", fn1)
            fn2 = fn1.replace(jobID1, jobID2)
            if os.path.exists(fn2):
                #print "Already exists:\t",fn2
                continue

        # only does physics
            if physicsonly:
                if fn1.lower().replace('-', '').replace('_',
                                                        '').find('diadt') > -1:
                    continue
                if fn1.lower().replace('-', '').replace('_',
                                                        '').find('ptrct') > -1:
                    continue

        # only does bgc.
            if bgconly:
                if fn1.lower().replace('-', '').replace('_',
                                                        '').find('grid') > -1:
                    continue

            try:
                os.symlink(fn1, fn2)
                print("linking:\t", fn1, '--->', fn2)
            except:
                print("linking:\t", fn1, '--->', fn2, 'FAILED')

    #####
    # shelve stuff: (same location on both machines!)
    if copyShelves:

        for fn1 in sorted(glob(shelvefold1 + '/*')):
            fn2 = fn1.replace(jobID1, jobID2)

            if os.path.exists(fn2):
                print("Already exists:\t", fn2)
                continue

            sh2 = folder(os.path.dirname(fn2))
            try:
                shutil.copy2(fn1, fn2)
                print("Copying:\t", fn1, '--->', fn2)
            except:
                print("Copying:\t", fn1, '--->', fn2, 'FAILED')

            if os.path.exists(fn2) and fn2.find('insitu') == -1:
                try:
                    s = shopen(fn2)
                    if 'readFiles' not in list(s.keys()):
                        s.close()
                        print("Shelve has no readFiles:", fn2)
                        continue
                except:
                    continue
                readFiles = s['readFiles']
                for i, fnshel in enumerate(readFiles):
                    newfnshe = fnshel.replace(jobID1, jobID2)
                    readFiles[i] = newfnshe
                    print("changing shelve contents from ", fnshel, "to",
                          newfnshe)
                s['readFiles'] = readFiles
                s.close()

    if linkchain:
        fn = netcdfFold1.replace(jobID1, 'linking.log')
        #		try:
        #			f = open(fn,'r')
        #			txt = f.read()
        #			f.close()
        #		except:
        #			txt = ''
        newline = str(datetime.now()) + ':\t' + jobID1 + ' -> ' + jobID2 + '\n'

        #	if txt.find(newline)>-1:
        #		txt+=newline

        f = open(fn, 'a')
        f.write(newline)
        f.close()


if __name__ == "__main__":

    try:
        jobID1 = argv[1]
        jobID2 = argv[2]
    except:
        print("Please provide two jobIDs")
        jobID1 = ''
        jobID2 = ''
        os.exit()

    copyShelves = False
    try:
        key = argv[3]
        if str(key) in ['True', 1, '1', 'copyShelves']:
            copyShelves = True
    except:
        copyShelves = False
    linkTwoJobs(jobID1, jobID2, copyShelves=copyShelves)
