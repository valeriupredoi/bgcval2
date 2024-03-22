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
.. module:: bgcval2_make_report
   :platform: Unix
   :synopsis: A script to produce an html5 document summarising a jobs performance.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

#####
# Load Standard Python modules:
import argparse
import sys

from glob import glob
import os
import shutil
import yaml

#####
# Load specific local code:
from bgcval2.bgcvaltools.bv2tools import folder, shouldIMakeFile, round_sig
from bgcval2.html5 import html5Tools, htmltables
from bgcval2.bgcvaltools.pftnames import getLongName
from bgcval2.timeseries.analysis_level0 import analysis_level0, analysis_level0_insitu
from bgcval2._runtime_config import get_run_configuration
from bgcval2.Paths.paths import paths_setter


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
        print("addImageToHtml:\tfn:", fn, 
            "\n\timagesfold:", imagesfold, 
            "\n\treportdir:", reportdir,
            "\n\tnewfn:", newfn,
            "\n\trelfn:", relfn)

    # sending a pdf as well.
    pdf_fn = fn.replace('.png', '.pdf')
    if os.path.exists(pdf_fn):
        new_pdf_fn = newfn.replace('.png', '.pdf')
        # Check if file is newer than the one in images.
        if shouldIMakeFile(
                pdf_fn,
                new_pdf_fn,
        ):
            if debug:
                print("addImageToHtml:\tAdding new pdf to report.", new_pdf_fn)
            if os.path.exists(new_pdf_fn):
                os.remove(new_pdf_fn)
            shutil.copy2(pdf_fn, new_pdf_fn)

    if not os.path.exists(newfn):
        if debug: 
            print("addImageToHtml:\tcopytree", fn, newfn)
        basedir = folder(os.path.dirname(newfn))
        if os.path.isdir(fn):
            shutil.copytree(fn, newfn, symlinks, ignore)
        else:
            shutil.copy2(fn, newfn)
    else:
        ####
        # Check if the newer file is the same one from images.
        if os.path.getmtime(fn) == os.path.getmtime(newfn): 
            return relfn
        ####
        # Check if file is newer than the one in images.
        if shouldIMakeFile(
                fn,
                newfn,
        ):
            if debug: 
                print("addImageToHtml:\tremoving old file", fn)
            os.remove(newfn)
            shutil.copy2(fn, newfn)
            if debug: 
                print("addImageToHtml:\t copy2", fn, newfn)
    return relfn


def html5Maker(
    jobID='u-ab749',
    reportdir='reports/tmp',
    year='*',
    clean=False,
    physicsOnly=False,
    paths=None,
    config_user=None,
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
    copytree(os.path.join(basedir, os.path.join(paths.bgcval2_repo,'bgcval2/html5/html5Assets')), reportdir)
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
    level2_auto = True
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
            'Phosphate',
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
                    'Nitrate', 'Phosphate', 'Silicate', 'DIC', 'Alkalinity', 'Temperature',
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
                                                       metric=m,
                                                       paths=paths)
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
                                                           metric=m,
                                                           paths=paths)
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
                                                   metric='metricless',
                                                   paths=paths)
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
                    paths=paths,
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
            'Phosphate',
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
            imagedir = paths.imagedir
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
                glob(imagedir + '/' + jobID + '/timeseries/*/min*' + key +
                     '*' + 'Global*min.png'))
            vfiles.extend(
                glob(imagedir + '/' + jobID + '/timeseries/*/max*' + key +
                     '*' + 'Global*max.png'))

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
                print('create plot headers:', fn)
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                ####
                # WOA fields that also produce transects, etc.
                if key in [
                        'Nitrate', 'Phosphate', 'Silicate', 'Temperature', 'Salinity',
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
            'Phosphate',
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
                print('adding Level1 regional plot:', jobID, fn)
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
                'Phosphate',
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
                    print('adding Level1 Profile plot:', jobID, fn)

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
            'Phosphate',
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
                print('adding Level2  plot:', jobID, fn)
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





    if level2Physics or level2_auto:
        if level2Physics: 
            l2Fields = [
                'Temperature',
                'Salinity',
                'MLD',
                'ZonalCurrent',
                'MeridionalCurrent',
                'VerticalCurrent',
            ]
            slices = [
                'Surface',
                '1000m',
                'Transect',
            ]
            SectionTitle = 'Level 2 - Physics'
            region = 'Global'

        if level2_auto:
            l2Fields = glob(imagedir + '/' + jobID + '/P2Pplots/*/*')
            l2Fields = [os.path.basename(fn) for fn in sorted(l2Fields)]
            levels = ['Surface', '4000m', '2000m', '1000m', '750m','500m','200m', '100m', '50m', 'Transect']
            outdict = {}
            outlevels = {}
            for i, fn in enumerate(l2Fields):
                for lv in levels:
                    if fn.find(lv) > -1:
                        outlevels[lv] = True
                    fn = fn.replace(lv, '')
                outdict[fn] = True
                
            l2Fields = [key for key, v in outdict.items()]
            slices = [key for key, v in outlevels.items()]
            SectionTitle = 'Level 2'
            region = '*'

        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
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
                        '*',
                ]:
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*/*/*' + s + '*' + region + '*' + key + '*' +
                             year + '*robinquad.png'))
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*/*/*' + s + '*' + region + '*' + key + '*' +
                             year + '*robinquad-cartopy.png'))
                    vfiles.extend(
                        glob(imagedir + '/' + jobID + '/P2Pplots/*/*' + key +
                             '*/*/*' + s + '*' + region + '*' + key + '*' +
                             year + '*.png'))
                        
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
                print('adding Level2 Feilds plot:', jobID, fn)
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
                print('adding Level3 plot:', jobID, fn)
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
                print('adding Level3 regional plot:', fn)
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
                    print('adding Level3 Hovmoeller:',  fn)
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
        vfiles.extend(glob(os.path.join(paths.bgcval2_repo,'bgcval2/html5/html5Assets/images/*Legend*.png')))
        print('adding Legend')

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


def comparehtml5Maker(
    jobIDs=[],
    reportdir='reports/tmp',
    files=[],
    clean=False,
    jobDescriptions={},
    jobColours={},
    paths = {},
    analysisKeys=[]
):
    """
    Generates the multi-job comparison report.
    
    jobIDs: list of job IDs
    reportdir: output job direction
    files: a list of paths to images
    clean: bool to remove old report.
    jobDescriptions: dict to describe eahc job
    jobColours: dict for each jobs colour in the legend.
    paths: bgcval.paths
    analysisKeys: list of analysis keys (generated from the input_yml).
    """
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
    #html5Assets_dir = 

    copytree(os.path.join(paths.bgcval2_repo,'bgcval2/html5/html5Assets'), reportdir)
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

    # A curated list of specific plots to include in the headlines
    # of the comparison report
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
        'MA_SST_Global_Surface',
        'MA_SSS_Global_Surface',
        'MA_Drake',
        'MA_AMOC_26N',
        'MA_AEU',
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
        'MA_Nitrate_Global_Surface',
        'MA_Phosphate_Global_Surface',
        'MA_TotalIntNPP', 
        'MA_TotalIntGPP',
        'MA_TotalPhytoC_Global_Surface',
        'MA_TotalZooC_Global_Surface',
    ]

    categories = {
        'Physics Key Metrics': [],
        'BGC Key Metrics': [],
        'Other Plots': [],
    }

    
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
        print('adding compare plots')
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
        otherFilenames = files[:]  
        SectionTitle = 'All Plots'

        hrefs = []
        Titles = {}
        SidebarTitles = {}
        Descriptions = {}
        FileLists = {}
        FileOrder = {}

        for key in sorted(analysisKeys):
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
                print('Adding ', key, 'compare plots:', fn)
                relfn = addImageToHtml(fn, imagesfold, reportdir)

                #####
                # Create custom title by removing extra bits.
                title = html5Tools.fnToTitle(relfn)

                FileLists[href][relfn] = title
                FileOrder[href][count] = relfn
                count += 1
                print("Adding ", relfn, "to script")

        if len(otherFilenames):
            # I think this is never happens anymore.
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
                print('Adding ', key, 'compare plots:', fn)
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
        vfiles.extend(glob(os.path.join(paths.bgcval2_repo,'bgcval2/html5/html5Assets/images/*Legend*.png')))
        print('Adding compare plots legend')
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
    parser.add_argument('-i',
                        '--job-id',
                        default=None,
                        help='Job ID',
                        required=True)
    parser.add_argument('-y',
                        '--year',
                        default=None,
                        help='Year',
                        required=False)
    parser.add_argument('-a',
                        '--clean',
                        action='store_true',
                        help='Clean or not',
                        required=False)
    parser.add_argument('-p',
                        '--physics',
                        action='store_true',
                        help='Physics or not',
                        required=False)
    parser.add_argument('-r',
                        '--report',
                        default=None,
                        help='Report repo',
                        required=False)

    args = parser.parse_args()

    return args


def main():
    """Run the html maker for a single job ID."""
    from ._version import __version__
    print(f'BGCVal2: {__version__}')

    args = get_args()
    jobID = args.job_id 

    if args.config_file:
        config_user = os.path.join(os.getcwd(), args.config_file)
        print(f"analysis_timeseries: Using user config file {config_user}")
    else:
        config_user = os.path.join(os.getcwd(), "bgcval2-config-user.yml")
        print(f"analysis_timeseries: Using user default file {config_user}")
    if not os.path.isfile(config_user):
        print(f"analysis_timeseries: Could not find configuration file {config_user}")
        config_user = None

    #defaults:
    clean = False
    physicsOnly = False
    year = '*'
    reportdir = folder('reports/' + jobID)

    if args.year:
        try:
            year = int(args.year)
        except ValueError:
            print("analysis_timeseries: Invalid input for year - must be an integer, got {args.year}")
    if args.clean:
        clean = True
        print("analysis_timeseries: Running with Clean option!")
    if args.physics:
        physicsOnly = True
        print("analysis_timeseries: Running with Physics option!")
    if args.report:
        reportdir = os.path.abspath(args.report)

    # get runtime configuration; not implemented yet
    if config_user:
        paths_dict, config_user = get_run_configuration(config_user)
    else:
        paths_dict, config_user = get_run_configuration("defaults")

    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)

    html5Maker(
        jobID=jobID,
        reportdir=reportdir,
        year=year,
        clean=clean,
        physicsOnly=physicsOnly,
        paths=paths,
        config_user=config_user
    )


if __name__ == "__main__":
    main()
