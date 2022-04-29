#!/usr/bin/ipython

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
.. module:: html5Tools
   :platform: Unix
   :synopsis: A swiss army knife of tools for html5 making.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

import os
from ..bgcvaltools.pftnames import getLongName


def AddtoFile(filepath, linenumber, text):
    """	Adds "text" at line "linenumber" in the file "filepath"
	"""
    f = open(filepath, "r")
    contents = f.readlines()
    f.close()

    contents.insert(linenumber, text)

    f = open(filepath, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()


def locateLineNumber(filepath, linekey):
    """	Locates the line above the first instance of "linekey" in the file "filepath."
	"""

    f = open(filepath, "r")
    contents = f.readlines()
    f.close()
    for l, line in enumerate(contents):
        if line.find(linekey) >= 0:
            return l


def writeSideBar(filepath, href, text, option='head', debug=False):
    """	Addes a sidebar to the file "filepath."
	"""
    print("Adding sidebar", href, text)

    if option == 'sub':
        newline = '\n\t\t\t\t\t\t\t\t<li><a href="#' + href + '" id="' + href + '-link" class="skel-layers-ignoreHref"><span class="icon fa-square">' + text + '</span></a></li>\n'
    if option == 'head':
        newline = '\n\t\t\t\t\t\t\t  <li><a href="#' + href + '" id="' + href + '-link" class="skel-layers-ignoreHref"><span class="icon fa-th">' + text + '</span></a></li>\n'

    if option == 'startSub': newline = '\n\t\t\t\t\t\t\t\t<ul>\n'
    if option == 'endSub': newline = '\n\t\t\t\t\t\t\t\t</ul>\n'

    linenumber = locateLineNumber(filepath, 'AddSideBar') - 1
    if debug:
        print("Adding sidebar at ", linenumber, "line:\n",
              newline.replace('\t', ''))
    AddtoFile(filepath, linenumber, newline)


def writeDescription(filepath, text, debug=False):
    """	Addes a sidebar to the file "filepath."
	"""
    newline = '\t\t\t\t\t\t\t<p>' + text + '</p>\n'
    linenumber = locateLineNumber(filepath, 'HeaderDescription')
    if debug:
        print("Adding descriptionat ", linenumber, "line:\n",
              newline.replace('\t', ''))
    AddtoFile(filepath, linenumber, newline)


def removeDuplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def fnToTitle(fn, jobID='u-a'):
    """	This script takes a filename and returned a usable description of the plot title.
	"""
    titlelist = os.path.basename(fn).replace('.png',
                                             '').replace('_', ' ').split(' ')

    titlelist = removeDuplicates(titlelist)

    IgnoreList = [
        'percentiles',
        'BGCVal',
        'MEDUSA',
        '',
        '10-90pc',
        '5-95pc',
    ]

    title = ''
    for i, t in enumerate(titlelist):
        # remove redundent versus field
        if t.find('vs') > -1: continue
        if t.find(jobID) == 0: continue
        if t in IgnoreList: continue
        if t in ['sum', 'Sum']: continue

        ln = getLongName(t)
        if ln in IgnoreList: continue
        title += ln + ', '
    return title[:-2]


def addImagesText(imagePath, title=''):
    """ Creates the appropriate text to insert an image onto the page.
	"""
    basedir = os.path.dirname(__file__)
    f = open(os.path.join(basedir, "figure-template.html"), "r")
    contents = f.readlines()
    f.close()
    if title == '': title = fnToTitle(imagePath)

    #####
    # Halves the title if it is too long
    if len(title) > 40:
        line = title.split(' ')
        line.insert(int(len(line) / 2), '<br>')
        title = ' '.join(line)

    for l, line in enumerate(contents):
        if line.find('imagefilename') >= 0:
            contents[l] = contents[l].replace('imagefilename', imagePath)
        if line.find('PlotTitle') >= 0:
            contents[l] = contents[l].replace('PlotTitle', title)

    outtxt = '\n'.join(contents)
    return outtxt


def AddSection(filepath, href, Title, Description='', Files=[]):
    """	Addes a section to the file "filepath."
	"""
    #####
    # Add a link to the side bar
    writeSideBar(filepath, href, Title)

    #####
    # Copy the template and add the images.
    f = open("html5/section-template.html", "r")
    contents = f.readlines()
    f.close()
    if type(Files) == type([
            'a',
            'list',
    ]):
        imagesTxt = '\n'.join([addImagesText(f) for f in Files])

    if type(Files) == type({
            'a': 'dict',
    }):
        imagesTxt = '\n'.join(
            [addImagesText(f, title=Files[f]) for f in sorted(Files.keys())])

    #####
    # Add Title, description and figures to the copied template
    for l, line in enumerate(contents):
        if line.find('inserthref') >= 0:
            contents[l] = contents[l].replace('inserthref', href)
        if line.find('Title') >= 0:
            contents[l] = contents[l].replace('Title', Title)
        if line.find('Description') >= 0:
            contents[l] = contents[l].replace('Description', Description)
        if line.find('Figures') >= 0: contents[l + 1] += imagesTxt

    #####
    # Convert the list into a string
    outtxt = '\n'.join(contents)

    #####
    # Add this into the template file.
    linenumber = locateLineNumber(filepath, 'AddSectionHere') - 1

    AddtoFile(filepath, linenumber, outtxt)


def AddTableSection(filepath,
                    href,
                    Title,
                    Description='',
                    Caption='',
                    tablehtml=[]):
    """	Addes a section to the file "filepath."
	"""
    #####
    # Add a link to the side bar
    writeSideBar(filepath, href, Title)

    #####
    # Copy the template and add the images.
    basedir = os.path.dirname(__file__)
    f = open(os.path.join(basedir, "section-template.html"), "r")
    contents = f.readlines()
    f.close()

    #####
    # Add Title, description and figures to the copied template
    for l, line in enumerate(contents):
        if line.find('inserthref') >= 0:
            contents[l] = contents[l].replace('inserthref', href)
        if line.find('Title') >= 0:
            contents[l] = contents[l].replace('Title', Title)
        if line.find('Description') >= 0:
            contents[l] = contents[l].replace('Description', Description)
        if line.find('Table') >= 0: contents[l + 1] += tablehtml
        if len(Caption) and line.find('TabCaption') >= 0:
            contents[l + 1] += '<p>' + Caption + '</p>'

    #####
    # Convert the list into a string
    outtxt = '\n'.join(contents)

    #####
    # Add this into the template file.
    linenumber = locateLineNumber(filepath, 'AddSectionHere') - 1

    AddtoFile(filepath, linenumber, outtxt)


def AddSubSections(filepath,
                   hrefs,
                   SectionTitle,
                   SidebarTitles={},
                   Titles={},
                   Descriptions={},
                   FileLists={},
                   FileOrder={}):
    """	Addes a section and a series of nested subsectionto the file "filepath."
	"""
    #####
    # Add a leading link to the side bar
    if len(hrefs) == 0:
        print("Empty AddSubSections:", hrefs, SectionTitle, SidebarTitles,
              Titles, Descriptions)
        assert 0
        return
    writeSideBar(filepath, hrefs[0], SectionTitle, option='head')
    writeSideBar(filepath, '', '', option='startSub')

    for href in hrefs:
        SidebarTitle = SidebarTitles[href]
        Title = Titles[href]
        Description = Descriptions[href]
        Files = FileLists[href]  #

        writeSideBar(filepath, href, SidebarTitle, option='sub')

        #####
        # Copy the template and add the images.
        basedir = os.path.dirname(__file__)
        f = open(os.path.join(basedir, "section-template.html"), "r")
        contents = f.readlines()
        f.close()
        if type(Files) == type([
                'a',
                'list',
        ]):
            try:
                Order = FileOrder[href]
            except:
                Order = {i: f for i, f in enumerate(sorted(Files))}
            imagesTxt = '\n'.join([
                addImagesText(Order[count]) for count in sorted(Order.keys())
            ])

        if type(Files) == type({
                'a': 'dict',
        }):
            try:
                Order = FileOrder[href]
            except:
                Order = {i: f for i, f in enumerate(sorted(Files.keys()))}
            imagesTxt = '\n'.join([
                addImagesText(Order[count], title=Files[Order[count]])
                for count in sorted(Order.keys())
            ])

        #####
        # Add Title, description and figures to the copied template
        for l, line in enumerate(contents):
            if line.find('inserthref') >= 0:
                contents[l] = contents[l].replace('inserthref', href)
            if line.find('Title') >= 0:
                contents[l] = contents[l].replace('Title', Title)
            if line.find('Description') >= 0:
                contents[l] = contents[l].replace('Description', Description)
            if line.find('Figures') >= 0: contents[l + 1] += imagesTxt

        #####
        # Convert the list into a string
        outtxt = '\n'.join(contents)

        #####
        # Add this into the template file.
        linenumber = locateLineNumber(filepath, 'AddSectionHere') - 1

        AddtoFile(filepath, linenumber, outtxt)
    writeSideBar(filepath, '', '', option='endSub')


#def AddImage(filepath, href, text):
#	"""	Addes a sidebar to the file "filepath."
#	"""
#	newline = '<li><a href="#'+href+'" id="'+href+'-link" class="skel-layers-ignoreHref"><span class="icon fa-th">'+text+'</span></a></'
#	linenumber = locateLineNumber(filepath, 'AddSideBar')
#	print "Adding at ", linenumber,"line:\n",newline
#	AddtoFile(filepath,linenumber,newline)
