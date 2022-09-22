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
.. module:: download_from_mass
   :platform: Unix
   :synopsis: A set of tools to download the UKESM model run data from MASS.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
#####
# Load Standard Python modules:
import argparse

from sys import stdout
import subprocess
from socket import gethostname
import shutil
import os
from glob import glob
from re import findall

#####
# Load specific local code:
from bgcval2._runtime_config import get_run_configuration
from bgcval2.Paths.paths import paths_setter


"""
This module includes a series of tools to download the UKESM model run data from MASS.

When run as a script, the command is::

	./download_from_mass.py jobID

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

def get_paths(config_user):
    # get runtime configuration
    if config_user:
        paths_dict, config_user = get_run_configuration(config_user)
    else:
        paths_dict, config_user = get_run_configuration("defaults")
    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)
    return paths

def getYearFromFile(fn):
    """
	Takes a file name, and looks for 8 consecutive numbers, then removes those that are months, and returns the year.
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
        if debug: 
            print("download_from_mass:\trebaseSymlinks:\tfile is not a symlink.", fn)
        return

#####
# Real path and first target:
    realpath = os.path.realpath(fn)  # The final end of the link chin
    linkpath = os.readlink(fn)  # The first target in the symlink chain

    if realpath == linkpath: return

    if debug:
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

    outputFold = folder([paths.ModelFolder_pref,  jobID,] )

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
                  stoptime=1E20,
                  config_user=None):
    """
	:param jobID: The job ID
	:param keys: a list of fields as they are saved in the Netcdf. (can also be a single string)
	:param timeslice: The time granularity (monthly or yearly)
	:param dryrun: does not download files, just prints.
	:param extension: Nemo style file extension
	:param name: Name of the analysis group, used for the folder.

	This tool takes the jobID, the field name, and using the known structure of universally similar MASS and the local filesystem structure
	from paths.py, downloads the monthly jobID data for the field requested to the local file structure.

	This tool will only work on machines that have connection to MASS enabled.

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
    paths = get_paths(config_user)
    outputFold = folder([paths.ModelFolder_pref, jobID, name])
    ####
    # ensure that permissions are : drwxrwsr-x+
    #os.chmod(paths.ModelFolder_pref + "/" + jobID, 0o2775)
    #os.chmod(outputFold, 0o2775)

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
    do_ls = False
    if do_ls:
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
def nemoMonthlyIce(jobID, dryrun=False, config_user=None):
    downloadField(jobID, [
        'soicecov',
    ],
                  extension='grid[-_]T',
                  timeslice='m',
                  name='monthlyIce',
                  dryrun=dryrun, 
                  config_user=config_user)  


def nemoMonthlyMLD(jobID, starttime=0, stoptime=1E20, dryrun=False, config_user=None):
    downloadField(jobID, [
        'somxl010',
    ],
                  extension='grid[-_]T',
                  timeslice='m',
                  name='monthlyMLD',
                  dryrun=dryrun,
                  starttime=starttime,
                  stoptime=stoptime,
                  config_user=config_user)


def monthlyChl(jobID, months=['01', '02', '06', '07', '08', '12'], dryrun=False, config_user=None):
    for month in months:  #['01','02','06','07','08','12']: # They want JJA and DJF
        ts = '????' + month + '01-??????01'
        downloadField(jobID, ['CHD', 'CHN'],
                      extension='ptrc[-_]T',
                      timeslice='m',
                      timerange=ts,
                      name='monthlyCHL',
                      dryrun=dryrun,
                      config_user=config_user)


def medusaMonthlyexport(jobID, dryrun=False, config_user=None):
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
                  dryrun=dryrun,
                  config_user=config_user)


def download_from_mass(
        jobID, 
        doMoo=True, 
        mass_shared_path=True,
        config_user=None
    ):
    """
	:param jobID: The job ID

	This tool takes the jobID, and using the known structure of universally similar MASS and the local filesystem structure
	from paths.py, downloads the jobID data to the local file structure.

	This tool will only work on machines that have mass enabled.

	"""
    if jobID == '': return

    paths = get_paths(config_user)
    outputFold = folder([paths.ModelFolder_pref,  jobID,] )

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
    do_ls=False
    if do_ls:
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
    else:
        output=''

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

    if mass_shared_path:
        shared_file_path = os.path.join(paths.shared_mass_scripts, os.path.basename(download_script_path))
        print('writing file in shared path', shared_file_path)
        shutil.copy(download_script_path, shared_file_path)

    fixFilePaths(outputFold, jobID, debug=False,)
    deleteBadLinksAndZeroSize(outputFold, jobID, debug=False,)


def fixFilePaths(outputFold, jobID, debug=False):
    #####
    # The coupled model looses the first two characters of the name in the netcdf file.
    fns = glob(outputFold + "/*" + jobID[2:] + "*.nc")
    if debug:
        print("download_from_mass:\tfixFilePaths:\tLooking for",
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
                print("download_from_mass:\tfixFilePaths:\tcorrect path exists.",
                      correctfn)
            continue
        if correctfn == fn: continue
        if debug:
            print("download_from_mass:\tfixFilePaths:\tFixing file prefix", fn,
                  '-->', correctfn)
        try:
            os.symlink(fn, correctfn)
        except:
            if debug:
                print("Unable to make link:", correctfn)
            continue

    #####
    # Some runs have nemo/medusa as a preface to the file name.
    for pref in ['nemo_', 'medusa_']:
        fns = glob(outputFold + "/" + pref + jobID + "*.nc")
        if debug:
            print("download_from_mass:\tfixFilePaths:\tLooking for new prefix:",
                  pref, outputFold + "/" + pref + jobID + "*.nc")
        for fn in sorted(fns):
            #####
            correctfn = os.path.dirname(fn) + '/' + os.path.basename(
                fn).replace(pref, '')
            if os.path.exists(correctfn):
                if debug:
                    print(
                        "download_from_mass:\tfixFilePaths:\tcorrect path exists.",
                        correctfn)
                continue
            if debug:
                print("download_from_mass:\tfixFilePaths:\tFixing file prefix",
                      pref,
                      end=' ')
            os.symlink(fn, correctfn)
            if debug:
                print("download_from_mass:\tfixFilePaths:\t", correctfn)

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
        if debug:
            print("download_from_mass:\tfixFilePaths:\tLooking for new suff:",
                  badsuff, outputFold + "/" + jobID + "*" + badsuff + ".nc")
        for fn in sorted(fns):
            #####
            correctfn = os.path.dirname(fn) + '/' + os.path.basename(
                fn).replace(badsuff, suff)
            if os.path.exists(correctfn):
                if debug:
                    print(
                        "download_from_mass:\tfixFilePaths:\tcorrect path exists.",
                        correctfn)
                continue
            if debug:
                print("download_from_mass:\tfixFilePaths:\tFixing file suffix",
                      badsuff,
                      '->',
                      suff,
                      end=' ')
            if correctfn == fn: continue

            try:
                os.symlink(fn, correctfn)
            except:
                continue
            if debug:
                print("download_from_mass:\tfixFilePaths:\t", correctfn)

    #####
    # This code looks at symoblic links and points them at their ultimate source, removing the long link chains.
    for fn in glob(outputFold + '/*'):
        rebaseSymlinks(fn, dryrun=False, debug=False)


def deleteBadLinksAndZeroSize(outputFold, jobID, debug=True):

    bashCommand1 = "find " + outputFold + "/. -size 0 -print -delete"
    bashCommand2 = "find -L " + outputFold + "/. -type l -delete  -print"

    if debug: print("deleteBadLinksAndZeroSize:\t", bashCommand1)

    process1 = subprocess.Popen(bashCommand1.split(), stdout=subprocess.PIPE)
    output1 = process1.communicate()[0]

    if debug: print("deleteBadLinksAndZeroSize:\t", bashCommand2)

    process2 = subprocess.Popen(bashCommand2.split(), stdout=subprocess.PIPE)
    output2 = process2.communicate()[0]


def pop_keys(keys, remove_keys):
   for k in remove_keys:
       if k in keys:
           keys.remove(k)

   return keys


def perform_download(jobID, keys, doMoo, config_user=None):
    """
    Single model download.
    """
    #####
    # Default behaviour is to download annual files
    if not keys or 'annual' in keys:
        download_from_mass(jobID, doMoo=doMoo, config_user=config_user)
        keys = pop_keys(keys, ['annual', ])

    dryrun = not doMoo
    #####
    # Monthly Ice files
    if 'ice' in keys or 'soicecov' in keys:
        nemoMonthlyIce(jobID, dryrun=dryrun, config_user=config_user)
        keys = pop_keys(keys, ['ice', 'soicecov'])

    #####
    # Monthly MLD
    if 'mld' in keys  or 'MLD' in keys:
        nemoMonthlyMLD(jobID, starttime=0, stoptime=5000,dryrun=dryrun, config_user=config_user)
        keys = pop_keys(keys, ['mld', 'MLD'])

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
        keys = pop_keys(keys, ['chl', 'CHL'])

    if 'export' in keys:
        medusaMonthlyexport(jobID, dryrun=dryrun, config_user=config_user)
        keys = pop_keys(keys, ['export',])

    #####
    # Other specific monthly files.
    if keys:
        downloadField(jobID, keys, timeslice='m', dryrun=dryrun, config_user=config_user)


def get_args():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-j',
                        '--jobID',
                        nargs='+', type=str,
                        help='One or more JobIDs to download',
                        required=True)
    parser.add_argument('-k',
                        '--keys',
                        default=['annual', ],
                        nargs='+', type=str,
                        help='Download keys - default options are: annual (which downloads all the annual files), '
                             'or chl, mld, ice, export, which downkoads monthly files for these fields. '
                             'Note that monthly files download is unstable and slow.',
                        required=False)
    parser.add_argument('-d',
                        '--dry-run',
                        action='store_true', 
                        help='Dry run: Do not download any files.',
                        )

    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'default-bgcval2-config.yml'),
                        help='User configuration file (for paths).',
                        required=False)

    args = parser.parse_args()

    return args


def main():
    """Run the main routine."""
    args = get_args()

    jobIDs = args.jobID
    keys = args.keys
    dryrun = args.dry_run
    doMoo = not dryrun
    config_user = args.config_file

    if keys in [None, '', [],]:
        keys = []
    if keys:
        keys = [str(k) for k in keys]

    print(f"Running with job_ids: {jobIDs} and keys {keys}")

    for jobID in jobIDs:
        perform_download(jobID, keys, doMoo, config_user=config_user)


if __name__ == "__main__":
    main()
