
"""
In this script, we remove years of data from a shelve
"""
from shelve import open as shOpen
import glob
import os
import argparse


from bgcval2.bgcvaltools import bv2tools as bvt
from bgcval2._runtime_config import get_run_configuration
from bgcval2.Paths.paths import paths_setter


def get_paths(
        config_user=None
    ):
    # get runtime configuration
    if config_user:
        paths_dict, config_user = get_run_configuration(config_user)
    else:
        paths_dict, config_user = get_run_configuration("defaults")

    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)
    return paths.shelvedir


def load_all_datatypes(shelvedir, jobID):
    """
    Returns a 
    """ 
    wildcards = shelvedir + '_'.join([
            jobID,
            '*',
        ]) + '.shelve.dat'
    files = glob.glob(wildcards)
    datatypes = []
    for fn in files:
        if fn.find('insitu') > -1:
            continue
        basename = os.path.basename(fn)
        basename = basename.replace(jobID+'_', '')
        basename = basename.replace('.shelve.dat', '')
        datatypes.append(basename)
    return datatypes



def remove_data(jobID, 
        year, 
        month=None,
        config_user=None,
        dataTypes= ['MA_AMOC_26N', ],
        dry_run = True,
        ):
    """
    Remove all data from the year and month from shelve files.
    """

    path_shelvedir = get_paths(config_user=config_user) 

    shelvedir = bvt.folder([path_shelvedir, "timeseries", jobID])
   
    if dataTypes == ['all', ]:
        dataTypes = load_all_datatypes(shelvedir, jobID)

    for dataType in dataTypes:
        shelvefn = shelvedir + '_'.join([
            jobID,
            dataType,
        ]) + '.shelve'

        if glob.glob(shelvefn+'*'):
            sh = shOpen(shelvefn)
            print('Shelve loads okay:', shelvefn +'*')
            sh = shOpen(shelvefn)
            readFiles       = sh['readFiles']
            modeldataD      = sh['modeldata']
            sh.close()

        # generate a time key from year/month.
        if month:
            if isinstance(month, str) and len(month) == 1:
                month = int(month)
            if isinstance(month, int):
                month = bvt.mnStr(month)
            
            time_key = year + month
        else:
            time_key = year

        changes = 0
        remove_files = [] 
        for readFile in readFiles:
            if readFile.find(time_key) > -1:
                remove_files.append(readFile)
                changes +=1
        
        if not changes:
            print('Nothing to remove')
            continue

        # Remove files from list
        for remove_file in remove_files:
            print('Removing', remove_file)
            if dry_run:
                pass
            else:
                readFiles.remove(remove_file)

        # remove processed data from file.
        for (r, l, m), values in modeldataD.items(): 
            key_removes = []
            for time in values.keys():
                if int(time) == int(year): 
                    # year matches:
                    if not month: 
                        # No month (delete all entries from this year
                        key_removes.append(time)
                        continue
                    # search for month:
                    mn = int((time - int(year)) *12)
                    if mn == int(month):
                        # found same month
                        key_removes.append(time)

            for key_remove in key_removes:
                print('Removing', key_remove)
                if dry_run:
                    pass
                else:
                    del values[key_remove]

        # Save file.
        if not dry_run:
            print('Saving:', shelvefn)
            sh = shOpen(shelvefn)
            sh['readFiles'] = readFiles
            sh['modeldata'] = modeldataD
            sh.close()
        else:
           print('Not saving (dry_run):', shelvefn)


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.getcwd(),
                                             'config-user.yml'),
                        help='User configuration file',
                        required=False)

    parser.add_argument('-j',
                        '--jobids',
                        default=None,
                        nargs='+',
                        type=str,
                        help='One or more jobIDs (required)',
                        required=True)

    parser.add_argument('-y',
                        '--years',
                        default=None,
                        nargs='+',
                        type=str,
                        help='One or more years (required)',
                        required=True)
    parser.add_argument('-m',
                        '--months',
                        default=[None, ],
                        nargs='+',
                        type=str,
                        help='One or more months - default is all months',
                        required=False)

    parser.add_argument('-k',
                        '--keys',
                        default=['all',],
                        nargs='+',
                        type=str,
                        help='One or more datasets - default is everything (all available keys). ',
                        required=False)

    parser.add_argument('-d',
                        '--dry-run',
                        action='store_true',
                        help='Dry run: Do not edit any files.',
                        required=False
                        )


    args = parser.parse_args()

    return args



def main():
    args = get_args()
    jobIDs = args.jobids
    years = args.years
    months = args.months
    dataTypes = args.keys
    dry_run = args.dry_run

    for jobID in jobIDs:
        for year in years:
            for month in months:
                print('start: remove_data', jobID, year,month, dataTypes, dry_run)
                remove_data(jobID, year, month=month,dataTypes=dataTypes, dry_run=dry_run) 



if __name__ == "__main__":
    main()

