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
.. module:: makeReport
   :platform: Unix
   :synopsis: A script to produce an html5 document summarising a jobs performance.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

#####
# Load Standard Python modules:
from glob import glob
from sys import argv
import os
import shutil

#####
# Load specific local code:
from .UKESMpython import folder, shouldIMakeFile, round_sig
from .html5 import html5Tools, htmltables
from .bgcvaltools.pftnames import getLongName
from .timeseries.analysis_level0 import analysis_level0, analysis_level0_insitu

#####
# User defined set of paths pointing towards the datasets.
#from .Paths.paths import paths_setter

#from .Paths.paths import imagedir

#    # filter paths dict into an object that's usable below
#    paths = paths_setter(paths_dict)


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            try:
                shutil.copytree(s, d, symlinks, ignore)
            except:
                pass
        else:
            try:
                shutil.copy2(s, d)
            except:
                pass


def addImageToHtml(fn, imagesfold, reportdir, debug=True):
    #####
    # Note that we use three paths here.
    # fn: The original file path relative to here
    # newfn: The location of the new copy of the file relative to here
    # relfn: The location of the new copy relative to the index.html

    newfn = imagesfold + os.path.basename(fn)
    relfn = newfn.replace(reportdir, './')

    if debug:
        print('addImageToHtml, fn:', fn, os.path.isdir(fn))
        print('addImageToHtml, imagesfold:', imagesfold)
        print('addImageToHtml, reportdir:', reportdir)
        print('addImageToHtml, newfn:', newfn, os.path.isdir(newfn))
        print('addImageToHtml, relfn:', relfn)

    # it's a directory:
    if os.path.isdir(fn) and not os.path.exists(newfn):
        os.mkdir(newfn)
        return relfn

    if not os.path.exists(newfn):
        if debug: 
            print("cp", fn, newfn)
        shutil.copy2(fn, newfn)
    else:
        ####
        # Check if the newer file is the same one from images.
        if os.path.getmtime(fn) == os.path.getmtime(newfn): return relfn
        ####
        # Check if file is newer than the one in images.
        if shouldIMakeFile(
                fn,
                newfn,
        ):
            if debug: print("removing old file", fn)
            os.remove(newfn)
            shutil.copy2(fn, newfn)
            if debug: print("cp", fn, newfn)
    return relfn


def html5Maker(
    jobID='u-ab749',
    reportdir='reports/tmp',
    year='*',
    clean=False,
    doZip=False,
    physicsOnly=False,
):

    if clean:
        #####
        # Delete old files
        print("Removing old files from:", reportdir)
        try:
            shutil.rmtree(reportdir)
        except:
            pass

    reportdir = folder(reportdir)
    year = str(year)

    ####
    # Copy all necceasiry objects and templates to the report location:
    print("Copying html and js assets to", reportdir)
    basedir = os.path.dirname(__file__)
    copytree(os.path.join(basedir, 'html5/html5Assets'), reportdir)
    indexhtmlfn = reportdir + "index.html"
    try:
        os.rename(reportdir + 'index-template.html', indexhtmlfn)
    except:
        pass

    imagesfold = folder(reportdir + 'images/')

    def newImageLocation(fn):
        return imagesfold + os.path.basename(fn)

    #####
    #
    descriptionText = 'Validation of the job: ' + jobID
    if year != '*': descriptionText += ', in the year: ' + year

    html5Tools.writeDescription(
        indexhtmlfn,
        descriptionText,
    )

    #####
    # Two switches to turn on Summary section, and groups of plots of field and region.
    Level0 = True
    Level1 = True
    Level1Regional = True
    Level1Profiles = True
    level2Horizontal = True
    level2Physics = False
    summarySections = False
    level3OMZ = False
    Level3Salinity = False
    regionMap = True

    #####
    # A list of caveats linked to specific datasets or regions, or jobs.
    ListofCaveats, ListofCaveats_regions = {}, {}
    ListofCaveats[
        'ExportRatio'] = 'Note that there is no historic data set for this variable.'
    ListofCaveats[
        'Iron'] = 'Note that there is no suitable historic data set for this variable.'
    ListofCaveats[
        'MLD'] = 'Note that the Model MLD is calculated with based on a sigma_0 difference of 0.01 with the surface where as data uses as \
			sigma_0 difference of +/- 0.2 degrees from a depth on 10m.'

    ListofCaveats[
        'Nitrate'] = 'Note that MEDUSA uses Dissolved Inorganic Nitrogen (DIN) rather than nitrate. We assume that non-nitrate parts of DIN are of relatively minor importance and so assume that WOA nitrate is comparable to model DIN.'
    ListofCaveats[
        'TotalDust'] = "This field has the Mediterranean, Red, Caspian, Black and Baltic seas excluded in both the model and data."
    if jobID == 'u-ad371':
        ListofCaveats[
            'Chlorophyll_cci'] = 'Note that the Non-diatom chlorophyll failed in run:' + jobID
        ListofCaveats[
            'IntegratedPrimaryProduction_OSU'] = 'Note that the Non-diatom chlorophyll does not contribute to IntPP in this run:' + jobID

    if Level0:
        # Level 0 metrics
        # x	1. Globally integrated primary production
        # x	2. Global average export fraction (i.e. globally integrated export / globally integrated primary production)
        # x	3. Globally integrated CO2 flux
        # x	4. Globally averaged surface DIN concentration
        # x	5. Globally averaged surface silicic acid concentration
        # x	6. Globally averaged surface DIC concentration
        #	7. Globally averaged surface alkalinity concentration
        #	8-11. As 4-7, but for the surface Southern Ocean [*]
        # x	12. AMOC
        # x	13. Drake transport
        # 	14-15. Volume of the Arctic / Antarctic sea-ice
        #	16-17. Extent of the Arctic / Antarctic sea-ice

        SectionTitle = 'Level 0'
        href = 'level0-table'

        table_data = []

        Description = 'A set of metrics to describe the state of the run.'

        Caption = ''

        realData_dict = {
            'TotalIntegratedPrimaryProduction': 'Target: 40-60',
            'ExportRatio': 'Target: 0.15 - 0.25',
            'TotalAirSeaFluxCO2': 'Target: -0.1 - 0.1',
            'Nitrate': '',
            'Silicate': '',
            'DIC': '',
            'Alkalinity': '',
            'AMOC_26N': 'Target: 10-20',
            'DrakePassageTransport': '136.7 	&#177; 7.8 ',  #	&#177; is +/-
            'NorthernTotalIceExtent': '14.9 	&#177; 0.3',
            'SouthernTotalIceExtent': '19.8 	&#177; 0.6',
            'TotalOMZVolume': '',
        }
        realData_source = {
            'TotalIntegratedPrimaryProduction': '',
            'ExportRatio': '',
            'TotalAirSeaFluxCO2': '',
            'Nitrate': 'World Ocean Atlas',
            'Silicate': 'World Ocean Atlas',
            'DIC': 'GLODAP',
            'Alkalinity': '',
            'AMOC_26N': '',
            'DrakePassageTransport':
            'Cunningham, S. A., S. G. Alderson, B. A. King, and M. A. Brandon (2003), Transport and variability of the Antarctic Circumpolar Current in Drake Passage, J. Geophys. Res., 108, 8084, doi:10.1029/2001JC001147, C5.',
            'NorthernTotalIceExtent':
            'HadISST: Rayner, N. A.; et al (2003) Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century J. Geophys. Res.Vol. 108, No. D14, 4407 10.1029/2002JD002670',
            'SouthernTotalIceExtent':
            'HadISST: Rayner, N. A.; et al (2003) Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century J. Geophys. Res.Vol. 108, No. D14, 4407 10.1029/2002JD002670',
            'TotalOMZVolume': 'World Ocean Atlas',
        }
        units_dict = {
            'TotalIntegratedPrimaryProduction': 'Gt/yr',
            'ExportRatio': '',  # no units
            'TotalAirSeaFluxCO2': 'Pg C/yr',
            'Nitrate': 'mmol N/m&#179;',  # &#179; is superscript 3
            'Silicate': 'mmol Si/m&#179;',
            'DIC': 'mmol C/m&#179;',
            'Alkalinity': 'meq/m&#179;',
            'AMOC_26N': 'Sv',
            'DrakePassageTransport': 'Sv',
            'NorthernTotalIceExtent':
            'x 1E6 km&#178;',  # &#178; is superscript 2
            'SouthernTotalIceExtent': 'x 1E6 km&#178;',
            'TotalOMZVolume': 'm&#179;',
        }

        fields = [
            'TotalIntegratedPrimaryProduction',
            'ExportRatio',
            'TotalAirSeaFluxCO2',
            'Nitrate',
            'Silicate',
            'DIC',
            'Alkalinity',
            'TotalOMZVolume',
            'AMOC_26N',
            'DrakePassageTransport',
            'NorthernTotalIceExtent',
            'SouthernTotalIceExtent',
        ]
        physFields = [
            'AMOC_26N',
            'DrakePassageTransport',
            'NorthernTotalIceExtent',
            'SouthernTotalIceExtent',
        ]
        timestrs = []
        for field in fields:
            if physicsOnly and field not in physFields: continue
            if field in [
                    'Nitrate', 'Silicate', 'DIC', 'Alkalinity', 'Temperature',
                    'Salinity'
            ]:
                for (r, l, m) in [('Global', 'Surface', 'mean'),
                                  ('SouthernOcean', 'Surface', 'mean'),
                                  ('AtlanticSOcean', 'Surface', 'mean')]:

                    #####
                    # Data column:
                    try:
                        rdata = realData_dict[field]
                    except:
                        rdata = ''
                    if rdata == '':
                        rdata = analysis_level0_insitu(jobID=jobID,
                                                       field=field,
                                                       region=r,
                                                       layer=l,
                                                       metric=m)
                        if rdata == False: rdata = ''
                        else: rdata = round_sig(rdata, 4)

                    try:
                        u = ' ' + units_dict[field]
                    except:
                        u = ''
                    try:
                        source = realData_source[field]
                    except:
                        source = ''
                    datcol = str(rdata) + u

                    name, mdata, timestr = analysis_level0(jobID=jobID,
                                                           field=field,
                                                           region=r,
                                                           layer=l,
                                                           metric=m)
                    longname = getLongName(name, debug=1)
                    if False in [name, mdata, timestr]:
                        table_data.append([longname, '', datcol])
                        continue

                    if timestr not in timestrs: timestrs.append(timestr)

                    modcol = str(round_sig(mdata, 4)) + u

                    table_data.append([longname, modcol, datcol])

                if len(source):
                    Caption += '<br><b>' + getLongName(
                        field) + '</b>: The data was taken from: ' + source

            else:
                #####
                # Data column:
                try:
                    rdata = realData_dict[field]
                except:
                    rdata = ''
                if rdata == '':
                    rdata = analysis_level0_insitu(jobID=jobID,
                                                   field=field,
                                                   region='regionless',
                                                   layer='layerless',
                                                   metric='metricless')
                    if rdata == False: rdata = ''
                    else: rdata = round_sig(rdata, 4)

                try:
                    u = ' ' + units_dict[field]
                except:
                    u = ''
                try:
                    source = realData_source[field]
                except:
                    source = ''
                datcol = str(rdata) + u
                name, mdata, timestr = analysis_level0(
                    jobID=jobID,
                    field=field,
                )  #region='regionless', layer='layerless', metric='metricless')
                longname = getLongName(name, debug=1)
                if False in [name, mdata, timestr]:
                    table_data.append([longname, '', datcol])
                    continue
                if timestr not in timestrs: timestrs.append(timestr)

                longname = getLongName(name, debug=1)
                modcol = str(round_sig(mdata, 4)) + u

                #			Caption+= '<br><b>'+longname+':</b> Model is the mean of the range '+timestr +'. '+\
                #							'The data was taken from: '+source
                table_data.append([longname, modcol, datcol])

                if len(source):
                    Caption += '<br><b>' + longname + '</b>: The data was taken from: ' + source
        if len(timestrs):
            Caption += '<br><b> Model</b> is the mean of the annual means in the range ' + timestrs[
                0] + '. '

        l0htmltable = htmltables.table(
            table_data,
            header_row=['Property', 'Model', 'Data'],
            col_align=['left', 'center', 'center'],
        )

        html5Tools.AddTableSection(indexhtmlfn,
                                   href,
                                   SectionTitle,
                                   Description=Description,
                                   Caption=Caption,
                                   tablehtml=l0htmltable)

    if Level1:
        level1Fields = [
            'TotalIntegratedPrimaryProduction',
            'ExportRatio',
            'AirSeaFluxCO2',
            'Nitrate',
            'DIC',
            'pH',
            'Alkalinity',
            'Chlorophyll',
            'Chlorophyll_cci',
            'DiaFrac',
            'TotalDust',
            #'DiatomChlorophyll',
            #'NonDiatomChlorophyll',
            'TotalOMZVolume',
            'Temperature',
            'GlobalMeanTemperature',
            'IcelessMeanSST',
            'Salinity',
            'GlobalMeanSalinity',
            'TotalIceArea',
            'TotalIceExtent',
            'DrakePassageTransport',
            'AMOC_26N',
            'ADRC_26N',
            'ZonalCurrent',
            'MeridionalCurrent',
            'VerticalCurrent',
            'MLD',
            'MaxMonthlyMLD',
            'MinMonthlyMLD'
        ]
        lev1physFields = [
            'Temperature',
            'GlobalMeanTemperature',
            'IcelessMeanSST',
            'Salinity',
            'GlobalMeanSalinity',
            'TotalIceArea',
            'TotalIceExtent',
            'WeddelIceExtent',
            'DrakePassageTransport',
            'AMOC_26N',
            'ADRC_26N',
            'MIZ',
            'ZonalCurrent',
            'MeridionalCurrent',
            'VerticalCurrent',
            'MLD',
            'MaxMonthlyMLD',
            'MinMonthlyMLD',
        ]
        SectionTitle = 'Level 1'
        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}

        #region = 'Global'
        for key in level1Fields:
            #print "Make Report\tLevel 1:",key
            if physicsOnly and key not in lev1physFields: continue
            #print "Make Report\tLevel 1:",key
            #####
            # href is the name used for the html
            href = 'L1' + key + '-global'
            hrefs.append(href)
            #print "Make Report\tLevel 1:",key, href
            #####
            # Title is the main header, SidebarTitles is the side bar title.
            Titles[href] = getLongName(key)
            SidebarTitles[href] = getLongName(key)
            #print "Make Report\tLevel 1:",key, Titles[href]

            #####
            # Descriptions is a small sub-header
            desc = ''
            if key in list(ListofCaveats.keys()):
                desc += ListofCaveats[key] + '\n'
            #if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'
            Descriptions[href] = desc
            #print "Make Report\tLevel 1:",key, desc

            #####
            # A list of files to put in this group.
            FileLists[href] = {}
            #####
            # Determine the list of files:
            vfiles = glob(imagedir + '/' + jobID +
                          '/timeseries/*/percentiles*' + key + '*' +
                          'Global*10-90pc.png')
            #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
            #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))
            vfiles.extend(
                glob(imagedir + '/' + jobID + '/timeseries/*/sum*' + key +
                     '*' + 'Global*sum.png'))
            vfiles.extend(
                glob(imagedir + '/' + jobID + '/timeseries/*/mean*' + key +
                     '*' + 'Global*mean.png'))
            vfiles.extend(
                glob(imagedir + '/' + jobID + '/timeseries/*/*' + key + '*' +
                     'regionless*metricless.png'))

            #####
            # Exceptions:
            if key in [
                    'AirSeaFluxCO2',
            ]:
                vfiles.extend(
                    glob(
                        imagedir + '/' + jobID +
                        '/timeseries/*/sum*NoCaspianAirSeaFluxCO2_ignoreCaspian_layerless_sum.png'
                    ))

        #vfiles.extend(glob(imagedir+'/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
        #vfiles.extend(glob(imagedir+'/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad.png'))
        #vfiles.extend(glob(imagedir+'/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*scatter.png'))
        #vfiles.extend(glob(imagedir+'/'+jobID+'/Targets/'+year+'/*'+key+'*/BGCVal/*.png'))

        #####
        # Create plot headers for each file.
            for fn in vfiles:
                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                ####
                # WOA fields that also produce transects, etc.
                if key in [
                        'Nitrate', 'Silicate', 'Temperature', 'Salinity',
                        'Oxygen', 'DIC', 'Alkalinity'
                ] and fn.lower().find('surface') < 0:
                    continue
                if key in ['ExportRatio'] and fn.find('_' + key) < 0:
                    continue  # make sure it's the global one, not the local one.
                #####
                # Create custom title by removing extra bits.
                #title = filenameToTitle(relfn)

                FileLists[href][relfn] = html5Tools.fnToTitle(relfn)
                print("Adding ", relfn, "to script")
        #print "Make Report\tLevel 1:",key, FileLists[href]

        html5Tools.AddSubSections(
            indexhtmlfn,
            hrefs,
            SectionTitle,
            SidebarTitles=SidebarTitles,  #
            Titles=Titles,
            Descriptions=Descriptions,
            FileLists=FileLists)

    if Level1Regional:
        l1regions = [
            'Global',
            'SouthernOcean',
            'AtlanticSOcean',
            'NorthernSubpolarAtlantic',
            'NorthernSubpolarPacific',
            'Equator10',
            'ArcticOcean',
            'Remainder',
            'ignoreInlandSeas',
        ]

        regionalFields = [
            'AirSeaFluxCO2',
            'Nitrate',
            'Silicate',
            'Iron',
            'IntegratedPrimaryProduction_OSU',
            'Chlorophyll',
            'Chlorophyll_cci',
            'pH',
            'DiaFrac',
            'Dust',
            #'DiatomChlorophyll',
            #'NonDiatomChlorophyll',
            'OMZThickness',
            'OMZMeanDepth',
            'Temperature',
            'Salinity',
            'ZonalCurrent',
            'MeridionalCurrent',
            'VerticalCurrent',
            'MaxMonthlyMLD',
            'MinMonthlyMLD',
            #'TotalIceArea'
        ]
        physregionalFields = [
            'Temperature',
            'Salinity',
            'TotalIceArea',
            'TotalIceExtent',
            'ZonalCurrent',
            'MeridionalCurrent',
            'VerticalCurrent',
            'MaxMonthlyMLD',
            'MinMonthlyMLD',
        ]
        SectionTitle = 'Level 1 - regional'
        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
        FileOrder = {}
        for key in regionalFields:
            if physicsOnly and key not in physregionalFields: continue
            #if key not in ['Alkalinity','Nitrate']: continue

            href = 'L1region' + key  #+'-'+region

            desc = ''
            if key in list(ListofCaveats.keys()):
                desc += ListofCaveats[key] + '\n'
            #if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'

            hrefs.append(href)
            Titles[href] = getLongName(key)
            SidebarTitles[href] = getLongName(key)
            Descriptions[href] = desc
            FileLists[href] = {}
            FileOrder[href] = {}

            #####
            # Determine the list of files:
            # It preferentially plots

            vfiles = []
            for region in l1regions:
                regfiles = glob(imagedir + '/' + jobID +
                                '/timeseries/*/percentiles*' + key + '*' +
                                region + '*10-90pc.png')
                print(
                    "Adding",
                    imagedir + '/' + jobID + '/timeseries/*/percentiles*' +
                    key + '*' + region + '*10-90pc.png')
                vfiles.extend(regfiles)

                if len(regfiles): continue

                vfiles.extend(
                    glob(imagedir + '/' + jobID + '/timeseries/*/mean*' + key +
                         '*' + region + '*mean.png'))
                print(
                    "Adding", imagedir + '/' + jobID + '/timeseries/*/mean*' +
                    key + '*' + region + '*mean.png')
        #####
        # Create plot headers for each file.
            count = 0
            for fn in vfiles:
                #####
                # Skip transects, they'll be added below.
                if fn.find('Transect') > -1: continue
                #if fn.lower().find('surface')<0 or fn.lower().find('layerless')<0:continue
                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        html5Tools.AddSubSections(
            indexhtmlfn,
            hrefs,
            SectionTitle,
            SidebarTitles=SidebarTitles,  #
            Titles=Titles,
            Descriptions=Descriptions,
            FileLists=FileLists,
            FileOrder=FileOrder)

    if Level1Profiles:
        #for plottype in ['profile','profilehov']:	# with Hovs
        for plottype in [
                'profile',
        ]:  # without hovs.
            l1regions = [
                'Global',
                'SouthernOcean',
                'AtlanticSOcean',
                'NorthernSubpolarAtlantic',
                'NorthernSubpolarPacific',
                'Equator10',
                'ArcticOcean',
                'Remainder',
                'ignoreInlandSeas',
            ]

            regionalFields = [
                'Nitrate',
                'Silicate',
                'Iron',
                'DIC',
                'Chlorophyll',
                'pH',
                'Alkalinity',
                'Oxygen',
                'Temperature',
                'Salinity',
                'ZonalCurrent',
                'MeridionalCurrent',
                'VerticalCurrent',
            ]
            physregionalFields = [
                'Temperature',
                'Salinity',
                'ZonalCurrent',
                'MeridionalCurrent',
                'VerticalCurrent',
            ]

            if plottype == 'profile': SectionTitle = 'Level 1 - Profiles'
            if plottype == 'profilehov':
                SectionTitle = 'Level 1 - Hovmoeller plots'
            hrefs = []
            Titles = {}
            SidebarTitles = {}
            Descriptions = {}
            FileLists = {}
            FileOrder = {}
            for key in regionalFields:
                if physicsOnly and key not in physregionalFields: continue
                #if key not in ['Alkalinity','Nitrate']: continue

                href = 'L1' + plottype + '-' + key  #+'-'+region

                desc = ''
                if key in list(ListofCaveats.keys()):
                    desc += ListofCaveats[key] + '\n'
                #if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'

                hrefs.append(href)
                Titles[href] = getLongName(key)
                SidebarTitles[href] = getLongName(key)
                Descriptions[href] = desc
                FileLists[href] = {}
                FileOrder[href] = {}
                #####
                # Determine the list of files:
                vfiles = []
                for region in l1regions:
                    #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png'))
                    if plottype == 'profile':
                        vfiles.extend(
                            glob(imagedir + '/' + jobID +
                                 '/timeseries/*/profile_*' + key + '*' +
                                 region + '*mean.png'))
                    if plottype == 'profilehov':
                        vfiles.extend(
                            glob(imagedir + '/' + jobID +
                                 '/timeseries/*/profilehov_*' + key + '*' +
                                 region + '*mean.png'))
            #####
            # Create plot headers for each file.
                count = 0
                for fn in vfiles:
                    #####
                    # Skip transects, they'll be added below.
                    if fn.find('Transect') > -1: continue
                    #if fn.lower().find('surface')<0:continue

                    #####
                    # Copy image to image folder and return relative path.
                    relfn = addImageToHtml(fn, imagesfold, reportdir)

                    #####
                    # Create custom title by removing extra bits.
                    title = html5Tools.fnToTitle(relfn)

                    FileLists[href][relfn] = title
                    FileOrder[href][count] = relfn
                    count += 1
                    print("Adding ", relfn, "to script")

            html5Tools.AddSubSections(
                indexhtmlfn,
                hrefs,
                SectionTitle,
                SidebarTitles=SidebarTitles,  #
                Titles=Titles,
                Descriptions=Descriptions,
                FileLists=FileLists,
                FileOrder=FileOrder)

    if level2Horizontal:
        l2Fields = [
            'Nitrate',
            'Silicate',
            'DIC',
            'pH',
            'Alkalinity',
            'Oxygen',
            'Chlorophyll_cci',
            #'TotalIntegratedPrimaryProduction',
            'IntegratedPrimaryProduction_OSU',
            'AirSeaFluxCO2',
            'Dust',
            #'TotalOMZVolume',
            #'TotalAirSeaFluxCO2' ,
            'Temperature',
            'Salinity',
            'MLD',
            'ZonalCurrent',
            'MeridionalCurrent',
            'VerticalCurrent',
        ]
        physl2Fields = [
            'Temperature',
            'Salinity',
            'MLD',
            'ZonalCurrent',
            'MeridionalCurrent',
            'VerticalCurrent',
        ]
        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
        SectionTitle = 'Level 2'
        region = 'Global'
        slices = [
            'Surface',
            '1000m',
            'Transect',
        ]
        FileOrder = {}

        for key in l2Fields:
            if physicsOnly and key not in physl2Fields: continue
            #if key not in ['Alkalinity','Nitrate']: continue

            href = 'l2-' + key + '-' + region

            desc = ''
            if key in list(ListofCaveats.keys()):
                desc += ListofCaveats[key] + '\n'
            if region in list(ListofCaveats_regions.keys()):
                desc += ListofCaveats_regions[key] + '\n'

            hrefs.append(href)
            Titles[href] = getLongName(region) + ' ' + getLongName(key)
            SidebarTitles[href] = getLongName(key)
            Descriptions[href] = desc
            FileLists[href] = {}
            FileOrder[href] = {}
            #####
            # Determine the list of files:
            vfiles = []
            #vfiles = glob(imagedir+'/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
            #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
            #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))
            #vfiles.extend(glob(imagedir+'/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
            for s in slices:
                if s in [
                        'Surface',
                        '1000m',
                ]:
                    print(
                        "looking for", imagedir + '/' + jobID +
                        '/P2Pplots/*/*' + key + '*/*/*' + s + '*' + region +
                        '*' + key + '*' + year + '*robinquad.png')
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*/*/*' + s + '*' + region + '*' + key + '*' +
                             year + '*robinquad.png'))
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*/*/*' + s + '*' + region + '*' + key + '*' +
                             year + '*robinquad-cartopy.png'))
                if s in [
                        'Transect',
                ]:
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*Transect/*/*' + s + '*' + region + '*' + key +
                             '*' + year + '*transect.png'))
            if key in [
                    'Chlorophyll_cci',
                    'IntegratedPrimaryProduction_OSU',
                    'AirSeaFluxCO2',
                    'MLD',
            ]:
                vfiles.extend(
                    glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                         '*/*/*' + region + '*' + key + '*' + year +
                         '*robinquad.png'))
                vfiles.extend(
                    glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                         '*/*/*' + region + '*' + key + '*' + year +
                         '*robinquad-cartopy.png'))

        #####
        # Create plot headers for each file.
            count = 0
            for fn in vfiles:
                #####
                # Skip transects, they'll be added below.
                #if fn.find('Transect') >-1: continue

                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        html5Tools.AddSubSections(
            indexhtmlfn,
            hrefs,
            SectionTitle,
            SidebarTitles=SidebarTitles,  #
            Titles=Titles,
            Descriptions=Descriptions,
            FileLists=FileLists,
            FileOrder=FileOrder)

    if level2Physics:
        l2Fields = [
            'Temperature',
            'Salinity',
            'MLD',
            'ZonalCurrent',
            'MeridionalCurrent',
            'VerticalCurrent',
        ]
        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
        SectionTitle = 'Level 2 - Physics'
        region = 'Global'
        slices = [
            'Surface',
            '1000m',
            'Transect',
        ]
        FileOrder = {}

        for key in sorted(l2Fields):
            #if key not in ['Alkalinity','Nitrate']: continue

            href = 'l2p-' + key + '-' + region

            desc = ''
            if key in list(ListofCaveats.keys()):
                desc += ListofCaveats[key] + '\n'
            if region in list(ListofCaveats_regions.keys()):
                desc += ListofCaveats_regions[key] + '\n'

            hrefs.append(href)
            Titles[href] = getLongName(region) + ' ' + getLongName(key)
            SidebarTitles[href] = getLongName(key)
            Descriptions[href] = desc
            FileLists[href] = {}
            FileOrder[href] = {}
            #####
            # Determine the list of files:
            vfiles = []
            #vfiles = glob(imagedir+'/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
            #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
            #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))
            #vfiles.extend(glob(imagedir+'/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
            for s in slices:
                if s in [
                        'Surface',
                        '1000m',
                ]:
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*/*/*' + s + '*' + region + '*' + key + '*' +
                             year + '*robinquad.png'))
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*/*/*' + s + '*' + region + '*' + key + '*' +
                             year + '*robinquad-cartopy.png'))
                if s in [
                        'Transect',
                ]:
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*Transect/*/*' + s + '*' + region + '*' + key +
                             '*' + year + '*transect.png'))
            if key in [
                    'Chlorophyll_cci',
                    'IntegratedPrimaryProduction_OSU',
                    'AirSeaFluxCO2',
                    'MLD',
            ]:
                vfiles.extend(
                    glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                         '*/*/*' + region + '*' + key + '*' + year +
                         '*robinquad.png'))
                vfiles.extend(
                    glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                         '*/*/*' + region + '*' + key + '*' + year +
                         '*robinquad-cartopy.png'))

            #####
            # Create plot headers for each file.
            count = 0
            for fn in vfiles:
                #####
                # Skip transects, they'll be added below.
                #if fn.find('Transect') >-1: continue

                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        html5Tools.AddSubSections(
            indexhtmlfn,
            hrefs,
            SectionTitle,
            SidebarTitles=SidebarTitles,  #
            Titles=Titles,
            Descriptions=Descriptions,
            FileLists=FileLists,
            FileOrder=FileOrder)

    if level3OMZ:

        l3omzFields = [
            'ExtentMaps'
            #			  'O2',
            #			  'OMZ',
            #                          'ZonalCurrent','MeridionalCurrent','VerticalCurrent',
        ]
        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
        SectionTitle = 'Level 3 - Oxygen Minimum Zone'
        region = 'Global'
        slices = [
            '500m',
            '1000m',
            'Transect',
            'Surface',
        ]
        FileOrder = {}

        for key in sorted(l3omzFields):
            #if key not in ['Alkalinity','Nitrate']: continue

            href = 'l3omz-' + key + '-' + region

            desc = ''
            if key in list(ListofCaveats.keys()):
                desc += ListofCaveats[key] + '\n'
            if region in list(ListofCaveats_regions.keys()):
                desc += ListofCaveats_regions[key] + '\n'

            hrefs.append(href)
            Titles[href] = getLongName(region) + ' ' + getLongName(key)
            SidebarTitles[href] = getLongName(key)
            Descriptions[href] = desc
            FileLists[href] = {}
            FileOrder[href] = {}
            #####
            # Determine the list of files:
            vfiles = []
            #	for s in slices:
            #	    if s in ['Surface','1000m','500m']:
            vfiles.extend(
                glob(imagedir + '/' + jobID +
                     '/Level3/OMZ/ExtentMaps/*/*_Global.png'))

            #		vfiles.extend(glob(imagedir+'/'+jobID+'/Level3/OMZ/ExtendMaps/*/*.png'))
            #			    if s in ['Transect',]:
            #				vfiles.extend(glob(imagedir+'/'+jobID+'/Level3/OMZ/*'+key+'*Transect/*/*'+s+'*'+region+'*'+key+'*'+year+'*transect.png'))
            #if key in [	'Chlorophyll_cci',
            #	  	'IntegratedPrimaryProduction_OSU',
            #		'AirSeaFluxCO2',
            #		'MLD',
            #	  ]:
            #	vfiles.extend(glob(imagedir+'/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad.png'))
            #	vfiles.extend(glob(imagedir+'/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad-cartopy.png'))

            #####
            # Create plot headers for each file.
            count = 0
            for fn in vfiles:
                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        html5Tools.AddSubSections(
            indexhtmlfn,
            hrefs,
            SectionTitle,
            SidebarTitles=SidebarTitles,  #
            Titles=Titles,
            Descriptions=Descriptions,
            FileLists=FileLists,
            FileOrder=FileOrder)

    if Level3Salinity:
        l3sal_regions = [
            'NordicSea',
            'LabradorSea',
            'NorwegianSea',
        ]

        regionalFields = [
            'Salinity', 'Temperature', 'MLD', 'sowaflup', 'sohefldo',
            'sofmflup', 'sosfldow', 'soicecov', 'MaxMonthlyMLD'
        ]
        profileFields = [
            'Salinity',
            'Temperature',
        ]
        SectionTitle = 'Level 3 - Salinity time series'
        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
        FileOrder = {}
        for key in regionalFields:
            #if physicsOnly and key not in physregionalFields:continue
            #if key not in ['Alkalinity','Nitrate']: continue
            href = 'L3nassalinity' + key  #+'-'+region

            desc = ''
            if key in list(ListofCaveats.keys()):
                desc += ListofCaveats[key] + '\n'

            hrefs.append(href)
            Titles[href] = getLongName(key)
            SidebarTitles[href] = getLongName(key)
            Descriptions[href] = desc
            FileLists[href] = {}
            FileOrder[href] = {}

            #####
            # Determine the list of files:
            # It preferentially plots

            vfiles = []
            for region in l3sal_regions:
                regfiles = glob(imagedir + '/' + jobID +
                                '/timeseries/*/percentiles*' + key + '*' +
                                region + '*10-90pc.png')
                print(
                    "Adding",
                    imagedir + '/' + jobID + '/timeseries/*/percentiles*' +
                    key + '*' + region + '*10-90pc.png')
                vfiles.extend(regfiles)

                if len(regfiles): continue

                vfiles.extend(
                    glob(imagedir + '/' + jobID + '/timeseries/*/mean*' + key +
                         '*' + region + '*mean.png'))
                print(
                    "Adding", imagedir + '/' + jobID + '/timeseries/*/mean*' +
                    key + '*' + region + '*mean.png')
            #####
            # Create plot headers for each file.
            count = 0
            for fn in vfiles:
                #####
                # Skip transects, they'll be added below.
                if fn.find('Transect') > -1: continue
                #if fn.lower().find('surface')<0 or fn.lower().find('layerless')<0:continue
                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        html5Tools.AddSubSections(
            indexhtmlfn,
            hrefs,
            SectionTitle,
            SidebarTitles=SidebarTitles,  #
            Titles=Titles,
            Descriptions=Descriptions,
            FileLists=FileLists,
            FileOrder=FileOrder)

        #for plottype in ['profile','profilehov']:	# with Hovs
        for plottype in [
                'profile',
        ]:  # without hovs.
            if plottype == 'profile':
                SectionTitle = 'Level 3 - Salininty Profiles'
            hrefs = []
            Titles = {}
            SidebarTitles = {}
            Descriptions = {}
            FileLists = {}
            FileOrder = {}
            for key in profileFields:
                #if key not in ['Alkalinity','Nitrate']: continue

                href = 'L3nas' + plottype + '-' + key  #+'-'+region

                desc = ''
                if key in list(ListofCaveats.keys()):
                    desc += ListofCaveats[key] + '\n'
                #if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'

                hrefs.append(href)
                Titles[href] = getLongName(key)
                SidebarTitles[href] = getLongName(key)
                Descriptions[href] = desc
                FileLists[href] = {}
                FileOrder[href] = {}
                #####
                # Determine the list of files:
                vfiles = []
                for region in l3sal_regions:
                    #vfiles.extend(glob(imagedir+'/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png'))
                    if plottype == 'profile':
                        vfiles.extend(
                            glob(imagedir + '/' + jobID +
                                 '/timeseries/*/profile_*' + key + '*' +
                                 region + '*mean.png'))
                    if plottype == 'profilehov':
                        vfiles.extend(
                            glob(imagedir + '/' + jobID +
                                 '/timeseries/*/profilehov_*' + key + '*' +
                                 region + '*mean.png'))
                #####
                # Create plot headers for each file.
                count = 0
                for fn in vfiles:
                    #####
                    # Skip transects, they'll be added below.
                    if fn.find('Transect') > -1: continue
                    #if fn.lower().find('surface')<0:continue

                    #####
                    # Copy image to image folder and return relative path.
                    relfn = addImageToHtml(fn, imagesfold, reportdir)

                    #####
                    # Create custom title by removing extra bits.
                    title = html5Tools.fnToTitle(relfn)

                    FileLists[href][relfn] = title
                    FileOrder[href][count] = relfn
                    count += 1
                    print("Adding ", relfn, "to script")

            html5Tools.AddSubSections(
                indexhtmlfn,
                hrefs,
                SectionTitle,
                SidebarTitles=SidebarTitles,  #
                Titles=Titles,
                Descriptions=Descriptions,
                FileLists=FileLists,
                FileOrder=FileOrder)

    if regionMap:
        vfiles = []
        vfiles.extend(glob('html5/html5Assets/images/*Legend*.png'))
        relfns = [addImageToHtml(fn, imagesfold, reportdir) for fn in vfiles]
        print(relfns)
        href = 'regionMap_default'
        html5Tools.AddSubSections(
            indexhtmlfn,
            [
                href,
            ],
            'Legends',
            SidebarTitles={
                href: 'Regional and Transect Legends',
            },
            Titles={
                href: 'Regional and Transect Legends',
            },
            Descriptions={
                href:
                'Maps showing the boundaries of the regions and transects used in this report.',
            },
            FileLists={
                href: relfns,
            },
        )

    tar = "tar cfvz  report-" + jobID + ".tar.gz " + reportdir

    print("-------------\nSuccess\ntest with:\nfirefox", indexhtmlfn)
    print("To zip it up:\n", tar)
    if doZip:
        import subprocess
        subprocess.Popen(tar.split())


def comparehtml5Maker(
    jobIDs=['u-ab749', 'u-ad371'],
    reportdir='reports/tmp',
    files=[],
    clean=False,
    doZip=False,
    jobDescriptions={},
    jobColours={},
):

    if clean:
        #####
        # Delete old files
        print("Removing old files from:", reportdir)
        try:
            shutil.rmtree(reportdir)
        except:
            pass

    reportdir = folder(reportdir)

    ####
    # Copy all necceasiry objects and templates to the report location:
    print("Copying html and js assets to", reportdir)
    copytree('html5/html5Assets', reportdir)
    indexhtmlfn = reportdir + "index.html"
    try:
        os.rename(reportdir + 'index-compare-template.html', indexhtmlfn)
    except:
        pass

    imagesfold = folder(reportdir + 'images/')

    def newImageLocation(fn):
        return imagesfold + os.path.basename(fn)

    descriptionText = 'Comparison of the jobs: ' + ', '.join(jobIDs)

    html5Tools.writeDescription(
        indexhtmlfn,
        descriptionText,
    )

    #####
    if jobColours != {}:
        SectionTitle = 'Legend'
        Description = ''
        Caption = ''

        href = 'legendTable'
        table_data = []
        for jobID in jobIDs:
            colourhtml = "<span style=\"color: " + jobColours[
                jobID] + "\">" + jobColours[jobID].title() + "</span>"
            try:
                jobdesc = jobDescriptions[jobID]
            except:
                print("comparehtml5Maker:\tjobID:", jobID,
                      " not in jobDescriptions.")
                jobdesc = jobID
            jobIDlink = '<a href="http://gws-access.ceda.ac.uk/public/esmeval/bgcval/' + jobID + '">' + jobID + '</a>'
            table_data.append([jobIDlink, colourhtml, jobdesc])
            print("comparehtml5Maker:\tadded:",
                  [jobIDlink, colourhtml, jobdesc])

        htmltable = htmltables.table(
            table_data,
            header_row=[
                'Job ID',
                'Colour',
                'Description',
            ],
            col_align=[
                'center',
                'center',
                'left',
            ],
        )

        html5Tools.AddTableSection(indexhtmlfn,
                                   href,
                                   SectionTitle,
                                   Description=Description,
                                   Caption=Caption,
                                   tablehtml=htmltable)

    physicsKM = [
        'AMOC_26N',
        'ADRC_26N',
        'DrakePassage',
        'GlobalMeanTemperature',
        'GlobalMeanSalinity',
        'NorthernTotalIceExtent',
        'SouthernTotalIceExtent',
        'Temperature_Global_Surface',
        'Salinty_Global_Surface',
        'FreshwaterFlux_Global',
        'TotalHeatFlux',
    ]

    bgcKM = [
        'TotalAirSeaFluxCO2',
        'NoCaspianAirSeaFluxCO2',
        'TotalIntegratedPrimaryProduction',
        'TotalDust',
        'Chlorophyll_Global_Surface',
        'TotalOMZVolume',
        'ExportRatio',
        'Nitrate_Global_Surface',
        'DIC_Global_Surface',
        'pH_Global_Surface',
        'Alkalinity_Global_Surface',
        'Silicate_Global_Surface',
        'Iron_Global_Surface',
    ]
    categories = {
        'Physics Key Metrics': [],
        'BGC Key Metrics': [],
        'Other Plots': [],
    }
    extrafolds = []

    for fn in files:
        found = False
        for key in physicsKM:
            if found: continue
            sectionTitle = 'Physics Key Metrics'
            if fn.find(key) > -1:
                try:
                    categories[sectionTitle].append(fn)
                except:
                    categories[sectionTitle] = [
                        fn,
                    ]
                found = True
        for key in bgcKM:
            if found: continue
            sectionTitle = 'BGC Key Metrics'
            if fn.find(key) > -1:
                try:
                    categories[sectionTitle].append(fn)
                except:
                    categories[sectionTitle] = [
                        fn,
                    ]
                found = True
        if found: continue
        for extracat in extrafolds:
            if found: continue
            if fn.find('/' + extracat + '/') > -1:
                try:
                    categories[extracat].append(fn)
                except:
                    categories[extracat] = [
                        fn,
                    ]

        if found: continue
        try:
            categories['Other Plots'].append(fn)
        except:
            categories['Other Plots'] = [
                fn,
            ]

    categoryOrder = []
    if len(categories['Physics Key Metrics']):
        categoryOrder.append('Physics Key Metrics')
    if len(categories['BGC Key Metrics']):
        categoryOrder.append('BGC Key Metrics')
    for exf in extrafolds:
        if exf not in list(categories.keys()): continue
        if len(categories[exf]): categoryOrder.append(exf)

    if len(categories['Other Plots']): categoryOrder.append('Other Plots')
    categories['Other Plots'] = sorted(categories['Other Plots'])

    for cat in categoryOrder:

        #####
        # These plots get added below.
        if cat in [
                'Other Plots',
        ]:
            continue

        catfiles = categories[cat]
        if not len(catfiles): continue

        href = cat.replace(' ', '')
        Title = cat

        #####
        # Copy image to image folder and return relative path.
        relativeFiles = [
            addImageToHtml(catfn, imagesfold, reportdir) for catfn in catfiles
        ]

        ####
        # sort alphabeticaly.
        relativeFiles = sorted(relativeFiles)

        html5Tools.AddSection(indexhtmlfn,
                              href,
                              Title,
                              Description='',
                              Files=relativeFiles)

    if len(categories['Other Plots']):
        otherFilenames = files[:]  #categories['Other Plots'][:]
        SectionTitle = 'All Plots'

        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
        FileOrder = {}

        names = [
            'Chlorophyll',
            'MLD',
            'Nitrate',
            'Salinity',
            'Temperature',
            'Current',
            #'so',
            'Ice',
            'DIC',
            'pH',
            'DMS',
            'DiaFrac',
            'Dust',
            'Iron',
            'Silicate',
            'Alkalinity',
            'AMOC',
            'ADRC',
            'DrakePassage',
            'AirSeaFlux',
            'DTC',
            'Oxygen',
            'OMZ',
            'Production',
            'Export',
            'FreshwaterFlux',
            'HeatFlux',
            'soga',
            'scvoltot',
            'thetaoga',
            'scalarHeatContent',
        ]
        for key in sorted(names):
            #####
            # Determine the list of files:
            vfiles = []
            for ofn in otherFilenames:
                if ofn.find(key) > -1:
                    vfiles.append(ofn)

            for ofn in vfiles:
                otherFilenames.remove(ofn)

            if len(vfiles) == 0: continue

            href = 'OtherPlots-' + key
            desc = ''

            hrefs.append(href)
            Titles[href] = getLongName(key)
            SidebarTitles[href] = getLongName(key)
            Descriptions[href] = desc
            FileLists[href] = {}
            FileOrder[href] = {}

            #####
            # Create plot headers for each file.
            count = 0
            for fn in sorted(vfiles):
                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        if len(otherFilenames):
            href = 'OtherPlots-others'

            hrefs.append(href)
            Titles[href] = 'Everything Else'
            SidebarTitles[href] = 'Everything Else'
            Descriptions[href] = ''
            FileLists[href] = {}
            FileOrder[href] = {}

            #####
            # Create plot headers for each file.
            count = 0
            for fn in sorted(otherFilenames):
                #####
                # Copy image to image folder and return relative path.
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        print(hrefs)
        print(SectionTitle)
        print(SidebarTitles)
        print(Titles)
        print(Descriptions)
        print(FileLists)
        print(FileOrder)
        #		hrefs=[]
        if len(hrefs):
            html5Tools.AddSubSections(indexhtmlfn,
                                      hrefs,
                                      SectionTitle,
                                      SidebarTitles=SidebarTitles,
                                      Titles=Titles,
                                      Descriptions=Descriptions,
                                      FileLists=FileLists,
                                      FileOrder=FileOrder)

    legend = True
    if legend:
        vfiles = []
        vfiles.extend(glob('html5/html5Assets/images/*Legend*.png'))
        relfns = [addImageToHtml(fn, imagesfold, reportdir) for fn in vfiles]
        print(relfns)
        href = 'regionMap_default'
        html5Tools.AddSubSections(
            indexhtmlfn,
            [
                href,
            ],
            'Legends',
            SidebarTitles={
                href: 'Regional and Transect Legends',
            },
            Titles={
                href: 'Regional and Transect Legends',
            },
            Descriptions={
                href:
                'Maps showing the boundaries of the regions and transects used in this report.',
            },
            FileLists={
                href: relfns,
            },
        )

    print("-------------\nSuccess\ntest with:\nfirefox", indexhtmlfn)


def main():
    try:
        jobID = argv[1]
    except:
        print("Please provide a jobID next time")
        exit()

    #defaults:
    clean = False
    physicsOnly = False
    year = '*'
    reportdir = folder('reports/' + jobID)

    for i, arg in enumerate(argv):
        if i <= 1: continue

        try:
            y = int(arg)
            year = y
            continue
        except:
            pass

        if arg == 'clean':
            clean = True
            continue

        if arg == 'physics':
            physicsOnly = True
            continue

        reportdir = arg

    html5Maker(
        jobID=jobID,
        reportdir=reportdir,
        year=year,
        clean=clean,
        physicsOnly=physicsOnly,
    )


if __name__ == "__main__":
    main()
