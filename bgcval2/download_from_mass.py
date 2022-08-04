#!/usr/bin/python2.7

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
.. module:: downloadFromMass
   :platform: Unix
   :synopsis: A set of tools to download the UKESM model run data from MASS.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
#####
# Load Standard Python modules:
from sys import argv, stdout
import subprocess
from socket import gethostname
import os
from glob import glob
from re import findall

#####
# Load specific local code:
from .Paths import paths


"""
This module includes a series of tools to download the UKESM model run data from MASS.

When run as a script, the command is::

	./downloadFromMass.py jobID

This tool will only work on machines that have mass enabled.
 
"""


def folder(name):
    """ This snippet takes a string, makes the folder and the string.
            It also accepts lists of strings.
        """
    if type(name) == type(['a', 'b', 'c']):
        name = '/'.join(name, )
    if name[-1] != '/':
        name = name + '/'
    if os.path.exists(name) is False:
        os.makedirs(name)
        print('makedirs ', name)

####
# ensure that permissions are : drwxrwsr-x+
#os.chmod(name , 02775)
    return name


def mnStr(month):
    """ 
        :param month: An int between 1 and 100.
        Returns a 2 digit number string with a leading zero, if needed.
        """
    mn = '%02d' % month
    return mn


def getYearFromFile(fn):
    """ 
	Takes a file anem, and looks for 8 consequetive numbers, then removes those that are months, and returns the year.
	"""
    a = findall(r'\d\d\d\d\d\d\d\d', fn)
    a.reverse()  # prefer second year.
    datestrs = [
        '1130',
    ]
    datestrs.extend([mnStr(i) + '01' for i in range(1, 13)])
    for i in a:
        if i[-4:] in datestrs:
            yr = i[:4]
            return yr

    return False


def rebaseSymlinks(fn, dryrun=True, debug=False):
    """ 
	:param fn: A full path to a filename. It should be a symbolic link.
	:param dryrun: A boolean switch to do a trial run of this function.

	This function reduces a chain of symbolic links down to one. It takes a full path, 
	checks whether it is a sym link, then checks whether the real path is the  target of the link.
	If not, it replaces the target with the real path.
	
	"""

    #####
    # fn is a link file
    #if not os.path.exists(fn):
    #       print "rebaseSymlinks:\tfile does not exist.",fn
    #       return
    if not os.path.islink(fn):
        if debug: print("download_from_mass:\trebaseSymlinks:\tfile is not a symlink.", fn)
        return

#####
# Real path and first target:
    realpath = os.path.realpath(fn)  # The final end of the link chin
    linkpath = os.readlink(fn)  # The first target in the symlink chain

    if realpath == linkpath: return

    print("download_from_mass:\trebaseSymlinks:\tdeleting and re-linking ", fn, '-->', realpath)
    if dryrun: return
    os.remove(fn)
    os.symlink(realpath, fn)


def findLastFinishedYear(jobID, dividby=1, numberfiles=6):
    """
	:param jobID: The job ID, as elsewhere.
	:param 	dividby: Outputs every "dividby" years.
	:param numberfiles: The expected number of files per a year. (usually 6, but sometimes 4)

	This tool find the best year to have a close look at the model, by searching through the files
	and guessing which years are finished.
	
	"""
    if jobID == '': return

    machine = gethostname()
    if machine.find('jasmin') > -1:
        #		outputFold = "/group_workspaces/jasmin2/ukesm/BGC_data/"+jobID+'/'
        outputFold = "/gws/nopw/j04/ukesm/BGC_data/" + jobID + '/'
    if machine.find('monsoon') > -1:
        outputFold = "/projects/ukesm/ldmora/UKESM/" + jobID + '/'

    if gethostname().find('pmpc') > -1:
        outputFold = "/data/euryale7/scratch/ledm/UKESM/MEDUSA/" + jobID + '/'

    fnDict = {}
    files = sorted(glob(outputFold + jobID + 'o_1y_*_????_?.nc'))

    for fn in files:
        yr = getYearFromFile(fn)
        print("download_from_mass:\tgetYearFromFile:", fn, yr)
        try:
            fnDict[yr] += 1
        except:
            fnDict[yr] = 1

    years = sorted(fnDict.keys())
    years.reverse()

    print("download_from_mass:\t",years, fnDict)

    if len(years) == 0:
        print("download_from_mass:\tfindLastFinishedYear:\tNo files found.\t")
        return False

    if len(years) < dividby:
        print("download_from_mass:\tfindLastFinishedYear:\tLess than", dividby,
              "years of model run, returning first year:", years[-1])
        return years[0]

    for y in years:
        if int(y) % dividby != 0: continue

        print(y, ':', fnDict[y])
        if fnDict[y] >= numberfiles: return y

    print(
        "No correct year, there's probably a problem here findLastFinishedYear(",
        jobID, ")")
    print("Machine", machine)
    print("outputFold:", outputFold)
    return False
    #assert 0


def downloadField(jobID,
                  keys,
                  extension='grid[-_]T',
                  timeslice='m',
                  name='',
                  timerange='*',
                  dryrun=False,
                  starttime=0,
                  stoptime=1E20):
    """
	:param jobID: The job ID
	:param keys: a list of fields as they are saved in the Netcdf. (can also be a single string)
	:param timeslice: The time granularity (monthly or yearly)
	:param dryrun: does not download files, just prints.
	:param extension: Nemo style file extension
	:param name: Name of the analysis group, used for the folder.
	
	This tool takes the jobID, the field name, and using the known structure of universally similar MASS and the local filesystem structure
	from paths.py, downloads the monthly jobID data for the field requested to the local file structure.
	
	This tool will only work on machines that have mass enabled.
	
	"""

    if jobID == '': return

    #####
    # verify time granularity
    timeslice = str(timeslice).lower()
    if timeslice in ['monthly', 'month', '1m', 'm']: ts = 'm'
    elif timeslice in ['yearly', 'year', '1y', 'y']: ts = 'y'
    else: ts = 'y'

    #####
    # verify keys
    if type(keys) == type('string'): keys = [
            keys,
    ]
    if len(keys) == 1 and name == '':
        name = keys[0]

    #####
    # Verify output folder:
    outputFold = folder(paths.ModelFolder_pref + "/" + jobID + "/" + name)
    ####
    # ensure that permissions are : drwxrwsr-x+
    os.chmod(paths.ModelFolder_pref + "/" + jobID, 0o2775)
    os.chmod(outputFold, 0o2775)

    print("downloadField:", name, jobID, keys, timeslice, 'being saved to:',
          outputFold)

    #####
    # make query file:
    querytxt = '-v ' + ','.join(keys)
    queryfile = folder('queryfiles/') + name + '.txt'
    qf = open(queryfile, 'w')
    qf.write(querytxt)
    qf.close()
    print("downloadField:\tquery text", querytxt)

    #####
    # moose file path:
    massfp = "moose:/crum/" + jobID + "/on" + ts + ".nc.file/*_1" + ts + '_' + timerange + '_' + extension + ".nc"
    print("downloadField:\tmoose path:", massfp)

    ######
    # print files
    #bashCommand = "moo passwd -r"
    #print "running the command:",bashCommand
    #process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    #output = process.communicate()[0]
    #print "output",output
    #assert 0

    ######
    # print files
    bashCommand = "moo ls " + massfp
    print("running the command:", bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = process.communicate()[0]
    print("output", output)

    #####
    # create a bash command to download the files:
    shell = '/bin/bash'
    bashCommand = '#!' + shell + ' \n'
    filesToDL = 0
    for l, line in enumerate(output.split('\n')):
        if line in [
                '',
                ' ',
        ]: continue
        outfn = outputFold + os.path.basename(line)
        if os.path.exists(outfn): continue
        yr = int(getYearFromFile(line))
        if not starttime < yr < stoptime:
            print("File outside requireste time:", line, 'not ', starttime,
                  '<', yr, '<', stoptime)
            continue
        bashCommand += "\nmoo filter --fill-gaps " + queryfile + " " + line + " " + outfn
        filesToDL += 1
    bashCommand += "\necho \"The End of " + jobID + ' ' + name + "\"\n"
    print("running the command:\n######\n", bashCommand)
    bashfile = folder('queryfiles/') + jobID + '-' + name + '.sh'
    print("Saved script as ", bashfile)
    qf = open(bashfile, 'w')
    qf.write(bashCommand)
    qf.close()

    #####
    # Run the bash command to download the files:
    if not dryrun and filesToDL > 0:
        script = shell + " " + bashfile
        print(script)
        process = subprocess.Popen(
            script.split(),
            stdout=subprocess.PIPE,
        )  #shell=True, executable=shell)
        output = process.communicate()
        print("bash out", output)

    fixFilePaths(outputFold, jobID)


######
# Some spefici wrappers for the downloadField
def nemoMonthlyIce(jobID):
    downloadField(jobID, [
        'soicecov',
    ],
                  extension='grid[-_]T',
                  timeslice='m',
                  name='monthlyIce',
                  dryrun=False)


def nemoMonthlyMLD(jobID, starttime=0, stoptime=1E20):
    downloadField(jobID, [
        'somxl010',
    ],
                  extension='grid[-_]T',
                  timeslice='m',
                  name='monthlyMLD',
                  dryrun=False,
                  starttime=starttime,
                  stoptime=stoptime)


def monthlyChl(jobID, months=['01', '02', '06', '07', '08', '12']):
    for month in months:  #['01','02','06','07','08','12']: # They want JJA and DJF
        ts = '????' + month + '01-??????01'
        downloadField(jobID, ['CHD', 'CHN'],
                      extension='ptrc[-_]T',
                      timeslice='m',
                      timerange=ts,
                      name='monthlyCHL',
                      dryrun=False)


def medusaMonthlyexport(jobID):
    downloadField(jobID, [
        'SDT__100',
        'FDT__100',
        'SDT__200',
        'FDT__200',
        'SDT__500',
        'FDT__500',
        'SDT_1000',
        'FDT_1000',
        'PRD',
        'PRN',
    ],
                  extension='diad[-_]T',
                  timeslice='m',
                  name='monthlyExport',
                  dryrun=False)


def download_from_mass(jobID, doMoo=True):
    """
	:param jobID: The job ID
	
	This tool takes the jobID, and using the known structure of universally similar MASS and the local filesystem structure
	from paths.py, downloads the jobID data to the local file structure.
	
	This tool will only work on machines that have mass enabled.
	
	"""
    if jobID == '': return

    machine = gethostname()
    knownmachine = False
    if machine.find('jasmin') > -1:
        knownmachine = True
        # outputFold = "/group_workspaces/jasmin2/ukesm/BGC_data/"+jobID
        outputFold = "/gws/nopw/j04/ukesm/BGC_data/" + jobID + '/'

        if not os.path.exists(outputFold):
            print("Making ", outputFold)
            os.makedirs(outputFold)
    ####
    # ensure that permissions are : drwxrwsr-x+
        os.chmod(outputFold, 0o2775)

    if machine.find('monsoon') > -1:
        knownmachine = True
        outputFold = "/projects/ukesm/ldmora/UKESM/" + jobID

        if not os.path.exists(outputFold):
            print("Making ", outputFold)
            os.makedirs(outputFold)

    if not knownmachine:
        print("Are you running this on the correct machine?")
        print(
            "\tYou should be on mass-cli1.ceda.ac.uk at jasmin or on monsoon at the MO"
        )
        print("\tBut you're at", machine)
        print(
            "\tTo skip this warning, use the \"anymachine\" option at the command line"
        )
        return

    deleteBadLinksAndZeroSize(outputFold, jobID)

    fixFilePaths(outputFold, jobID)
    deleteBadLinksAndZeroSize(outputFold, jobID)

    # Set up a file to save command to a new file.
    download_script_path = ''.join([folder('mass_scripts/'), jobID,'.sh'])
    header_lines = ['# Run this script on mass-cli1.jasmin.ac.uk\n',]
    header_lines.append('# from login1.jasmin.ac.uk, ssh to the mass machine:\n#     ssh -X  mass-cli\n')
    header_lines.append(''.join(['# run script with:\n# source ', os.path.abspath(download_script_path),'\n']))
    header_lines.append('# moo passwd -r # if mass password is expired\n')
    download_script_txt = ''.join(header_lines)

    # moo ls: 
    bashCommand = "moo ls moose:/crum/" + jobID + "/ony.nc.file/*.nc"
    download_script_txt = ''.join([download_script_txt, bashCommand, '\n'])

    if not doMoo:
        print("download_from_mass:\tthe command is (dry-run): \n", bashCommand)
        output = ''
    else:
        print("download_from_mass:\trunning the command:", bashCommand)
        stdout.flush()
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = process.communicate()[0]

    # moo get:
    if len(output.split('\n')) > 6000:
        failed = 0
        process1 = {}
        output1 = {}
        for i in range(100):
            bashCommand = "moo get --fill-gaps moose:/crum/" + jobID + "/ony.nc.file/*_1y_??" + mnStr(
                i) + "*.nc " + outputFold
            print("download_from_mass:\trunning the command:", bashCommand)
            download_script_txt = ''.join([download_script_txt, bashCommand, '\n'])

            stdout.flush()
            try:
                process1[i] = subprocess.Popen(bashCommand.split(),
                                               stdout=subprocess.PIPE)
                process1[i].wait()
                output1[i] = process.communicate()[0]
                print("download_from_mass:\t",i, output1[i])
            except:
                failed += 1
                print("Failed", i, 'total fails:',failed, '\t', bashCommand)
    else:
        print("download_from_mass:\tDownloading at the following files:")
        bashCommand = "moo get --fill-gaps moose:/crum/" + jobID + "/ony.nc.file/*.nc " + outputFold
        download_script_txt = ''.join([download_script_txt, bashCommand, '\n'])
        if not doMoo:
            print("download_from_mass:\tthe command is (dry-run): \n", bashCommand)
        else:
            print("download_from_mass:\trunning the command:", bashCommand)
            process = subprocess.Popen(bashCommand.split(),
                                       stdout=subprocess.PIPE)
            output = process.communicate()[0]

    print('writing file:',download_script_path, '\nfile contents:\n', download_script_txt)
    outfile = open(download_script_path, 'w')
    outfile.write(download_script_txt)
    outfile.close()

    fixFilePaths(outputFold, jobID)
    deleteBadLinksAndZeroSize(outputFold, jobID)


def fixFilePaths(outputFold, jobID, debug=False):
    #####
    # The coupled model looses the first two characters of the name in the netcdf file.
    fns = glob(outputFold + "/*" + jobID[2:] + "*.nc")
    print("downloadFromMass:\tfixFilePaths:\tLooking for",
          outputFold + "/" + jobID[2:] + "*.nc")
    fns.extend(
        glob(outputFold +
             '/MetOffice*'))  # Because ocean assess might use the lisence?
    for fn in sorted(fns):
        #####
        correctfn = fn.replace(jobID[2:], jobID)
        correctfn = correctfn.replace('u-u-', 'u-')

        if os.path.exists(correctfn):
            if debug:
                print("downloadFromMass:\tfixFilePaths:\tcorrect path exists.",
                      correctfn)
            continue
        if correctfn == fn: continue
        print("downloadFromMass:\tfixFilePaths:\tFixing file prefix", fn,
              '-->', correctfn)
        try:
            os.symlink(fn, correctfn)
        except:
            print("Unable to make link:", correctfn)
            continue
#	        print "downloadFromMass:\tfixFilePaths:\t", correctfn

#####
# Some runs have nemo/medusa as a preface to the file name.
    for pref in ['nemo_', 'medusa_']:
        #nemo_u-ai886o_1y_26291201-26301201_grid-V.nc
        fns = glob(outputFold + "/" + pref + jobID + "*.nc")
        print("downloadFromMass:\tfixFilePaths:\tLooking for new prefix:",
              pref, outputFold + "/" + pref + jobID + "*.nc")
        for fn in sorted(fns):
            #####
            correctfn = os.path.dirname(fn) + '/' + os.path.basename(
                fn).replace(pref, '')
            if os.path.exists(correctfn):
                if debug:
                    print(
                        "downloadFromMass:\tfixFilePaths:\tcorrect path exists.",
                        correctfn)
                continue
            print("downloadFromMass:\tfixFilePaths:\tFixing file prefix",
                  pref,
                  end=' ')
            os.symlink(fn, correctfn)
            print("downloadFromMass:\tfixFilePaths:\t", correctfn)

#####
# Some runs have nemo/medusa as a preface to the file name.
    suffDict = {
        'grid-T': 'grid_T',
        'grid-U': 'grid_U',
        'grid-V': 'grid_V',
        'grid-W': 'grid_W',
        'diad-T': 'diad_T',
        'ptrc-T': 'ptrc_T',
        #'diaptr':'diad_T', # diaptr is not the same as BGC diad-T
    }
    for badsuff, suff in list(suffDict.items()):
        #nemo_u-ai886o_1y_26291201-26301201_grid-V.nc
        fns = glob(outputFold + "/" + jobID + "*" + badsuff + ".nc")
        print("downloadFromMass:\tfixFilePaths:\tLooking for new suff:",
              badsuff, outputFold + "/" + jobID + "*" + badsuff + ".nc")
        for fn in sorted(fns):
            #####
            correctfn = os.path.dirname(fn) + '/' + os.path.basename(
                fn).replace(badsuff, suff)
            if os.path.exists(correctfn):
                if debug:
                    print(
                        "downloadFromMass:\tfixFilePaths:\tcorrect path exists.",
                        correctfn)
                continue
            print("downloadFromMass:\tfixFilePaths:\tFixing file suffix",
                  badsuff,
                  '->',
                  suff,
                  end=' ')
            if correctfn == fn: continue

            try:
                os.symlink(fn, correctfn)
            except:
                continue
            print("downloadFromMass:\tfixFilePaths:\t", correctfn)

    #####
    # This code looks at symoblic links and points them at their ultimate source, removing the long link chains.
    for fn in glob(outputFold + '/*'):
        rebaseSymlinks(fn, dryrun=False)


def deleteBadLinksAndZeroSize(outputFold, jobID):

    bashCommand1 = "find " + outputFold + "/. -size 0 -print -delete"
    bashCommand2 = "find -L " + outputFold + "/. -type l -delete  -print"

    print("deleteBadLinksAndZeroSize:\t", bashCommand1)

    process1 = subprocess.Popen(bashCommand1.split(), stdout=subprocess.PIPE)
    output1 = process1.communicate()[0]

    print("deleteBadLinksAndZeroSize:\t", bashCommand2)

    process2 = subprocess.Popen(bashCommand2.split(), stdout=subprocess.PIPE)
    output2 = process2.communicate()[0]


def main():
    try:
        jobID = argv[1]
    except:
        print("Please provide a jobID")
        jobID = ''
    try:
        keys = argv[2:]
    except:
        keys = []

    #####
    # All yearly files

    if keys == []:
        download_from_mass(jobID, doMoo=True)
        return

    if 'noMoo' in keys or 'dryrun' in keys or '--dry-run' in keys:
        download_from_mass(jobID, doMoo=False)
        return  
    #####
    # Monthly Ice files
    if 'ice' in keys or 'soicecov' in keys:
        nemoMonthlyIce(jobID)
    #####
    # Monthly MLD
    if 'mld' in keys  or 'MLD' in keys:
        nemoMonthlyMLD(jobID, starttime=2570, stoptime=2610)
#####
# Monthly chl
    if 'chl' in keys:
        monthlyChl(jobID,
                   months=[
                       '01',
                       '02',
                       '03',
                       '04',
                       '05',
                       '06',
                       '07',
                       '08',
                       '09',
                       '10',
                       '11',
                       '12',
                   ])

    if 'export' in keys:
        medusaMonthlyexport(jobID)

    #####
    # Other specific monthly files.
#    else:
#        downloadField(jobID, keys, timeslice='m', dryrun=0)

if __name__ == "__main__":
    main()
