.. BGC-Val documentation master file, created by
   sphinx-quickstart on Fri Sep 23 10:33:38 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
===================================
Biogeochemistry Validation toolkit
===================================

The Biogeochemistry Validation toolkit (BGC-Val) is a toolkit built to assist with the validation of marine ocean models. 
The tools were originally built for analysis of UKESM annual fields, to help identify 
potential issues during the spin up phase. 



Accessing the toolkit
======================

The code itself is available upon request from the 
`PML-hosted gitlab server. <https://gitlab.ecosystem-modelling.pml.ac.uk/UKESM/bgc-val>`_

Access can be requested via the `Code registration form <http://www.pml.ac.uk/Modelling_at_PML/Access_Code>`_.
The PML-contact is "Lee de Mora, ledm@pml.ac.uk", please select "Other" in the Model code/Tutorial, and please put
"Access to BGC-val toolkit" as the main purpose of use.  Note that processing may take a few days.



SSH keys on gitlab
-----------------------
Before using gilab, you need to set up your ssh keys. 
**You will not be able to access the git repositories if you skip this step.**
On your working machine (where you want to use the code, but on the login node).
make a ssh key (leave passphrase blank)

	* Instructions `here <https://gitlab.ecosystem-modelling.pml.ac.uk/help/ssh/README>`_.
	
	* Copy paste the sshkey.pub file into your gitlab settings `here <https://gitlab.ecosystem-modelling.pml.ac.uk/profile/keys>`_.


There are some key points required before you can run git.
You'll need to set your name and email address. (git will prompt you to do this the first time)::

	git config --global user.name "your name"
	git config --global user.email "your_email@example.com"
	
You may also want to change your default editor. It's used mostly for writing commit messages::

	git config --global core.editor vim
	

Cloning the git repository
---------------------------
Once access to the gitlab server has been established, and you have set up your ssh keys, the toolkit can be
downloade with the command::

	git clone git@gitlab.ecosystem-modelling.pml.ac.uk:UKESM/bgc-val.git
        git clone git@gitlab.ecosystem-modelling.pml.ac.uk:ledm/netcdf_manip.git

The first command makes a local copy of the BGC-val toolkit.
You also will also need the netcdf_manip toolkit. Netcdf_manip is a python toolkit designed for manipulating netcdf files in python. 
The ukesm validation tools occassionally use them, and you'll need them to run the BGC validation code.
Each of the tools have a different role:

	:pruneNC: allows you to copy netcdf files while reducing the number of fields
	:mergeNC: merges multiple files along the time axis.
	:changeNC: does a lot of stuff
	:convertto1D: converts all non-masked values in a netcdf to a 1D array.

Once you have downloaded the Netcdf_manip repositories, you'll need to link it to pthon by adding some lines to your ~/.bashrc file, (or ~/.cshrc, ~/.tshrc if you use that shell).
Now, add the following lines to ~/.bahrc or ~/.cshrc ::

        export PYTHONPATH=$PYTHONPATH:/home/username/etc/pathToWorkingSpace/netcdf_manip 


Required python modules
------------------------

The following standard python modules are required to run BGC-Val. 
These modules should already be available, but can be installed using pip, 
or a similar package manager.

	* numpy
	* scipy
	* shelve
	* os
	* sys
	* socket
	* calendar
	* shutil
	* glob
	* getpass
	* string
	* math
	* itertools

The following modules are also required, but may be more complicated to install. 
Please see their websites for more details on how to install these packages. 
Note that these are likely to be pre-installed on post processing machines,
such as jasmin's  sci1 nodes and monsoon post-processing machines. 

	* `cartopy <http://scitools.org.uk/cartopy/>`_.
	* `Matplotlib, mpl_toolkits <http://matplotlib.org>`_.
	* `NetCDF4 <http://unidata.github.io/netcdf4-python/>`_.

The `sphinx module <http://www.sphinx-doc.org>`_ was used to create this documentation, but is not required.



BGC-Val Structure
======================

The rest of the cover page shows how to set up the analysis package. 
The biogeochemical validation toolkit is build with flexibility in mind.

The executable top level structure of this module includes the scripts:
	
	:theWholePackage.py: This script runs a specific
	:makeReport.py: This script creates an html report.
	:analysis_timeseries.py: This script runs the time series analysis.
	:analysis_p2p.py: This script runs the point to point analysis.
	:download_from_mass.py: This script is used to download UKESM files from mass.
	
These files are set up to run the validation is run for the UK Earth System Model in the eORCA1 grid.
They can be used as templates for running other models. 


Folders:
	 
	 :timeseries/: The tools required for the time series analysis.
	 :p2p/: The point to point analysis tools.
	 :bgcvaltools/: Miscellaneous tools used around the BGC-val toolkit.
	 :data/: These are data, masks, and coastline files required for some models. They do not include any observational data.
	 :html5/: The tools and assets required for the makeReport toolkits.
	 :Paths/: A set of paths that is used by many of the analyses to locate various datasets.
	 :RemoteScripts/: A set of scripts that execute the toolkit on jasmin from PML.
	 :LISENCES/: The software lisence(s).

	
Paths
---------------------------
This is a set of paths that is used by many of the analyses to locate various datasets.
The idea is that it is set once per machine, and then used in multiple analyses, ultimately saving time and effort.
It also sets the paths for storing the post processed data. 

Before running the analysis, the paths.py must be created in the main directory::

	cp Paths/paths_template.py paths.py

paths.py contains a list of paths to various locations, saved as strings.
These paths will be unique for each machine.
In order to run the analysis, edit your version of paths.py so that the path strings reflect your machines file system.


The important paths to set are:
	:shelvedir: A location for post processed data, saved in pythons shelve format.
	:esmvalFolder: A path for the model data.
	:ObsFolder: The base folder for many observational datasets. These will vary depending on your path location. 

The Path folder contains some other examples of paths.py files that have been used for various file systems.



Standard nomencature
=====================
There are some standard definitions used accross the software:

	:jobID: The jobID is used as a unique identifyer of each run. In the case of UKESM, jobID is automatically created by the ROSE suite. A unique jobId for each run is required to create unique paths for the post processed data.
	:region: The regions of the global ocean which to perform the analysis. These are usually assigned as a list of regions. 'Global' indicates the entire earth, but many other regions are available. Regions are defined in the bv2tools.py file in the function makeMask().
	:layer:  The depth layer used in the analysis.
	:suite:  The suite is a set of analysises that 
	:debug: A boolean flag for printing extra run-time debugging statements.
	
These definitions are laregly applicable to marine models. 


Autovivification
-----------------
An Autovivification is an object that behaves like a dictionairy in a dictionairy. 
It simifies the process of creating a nested dictionairy. ie::

	from bv2tools import Autovivification
	av = Autovivification()
	av['three']['level']['dictionairy'] = 'a string'

is similar to::
	
	av = {}
	av['three'] = {}
	av['three']['level'] = {}
	av['three']['level']['dictionairy'] = 'a string'	

Except the the Autovivification object does not require the user to specific that the
contents of the dictionary is another dictionairy. 


Coordinates dictionairies
--------------------------------------

Coordinates dictionairies are tools that are built into the times series, level0 and point to point analysis toolkits.
They are dictionairies with very specific contents that allow the analysis package to open, load and save the coordinate fields.
A coordinate dictionairy is required for each model and for each dataset type to be analysed.

The coordinates dictionairy includes the names of each coordinate as they are used in the netcdf file.
The contents of the coordinates dictionairy are:
	:t: The name used in the netcdf for the time coordinate.
	:z: The name used in the netcdf for the depth coordinate.
	:lat: The name used in the netcdf for the latitude coordinate.
	:lon: The name used in the netcdf for the longitude coordinate.
	:cal: The calendar used in the netcdf, this is important for converting certain run time fields (such as a 360 day year) into a gregorian date. 
	:tdict: The time-dictiorary, a dictionairy used to map the time coordinate onto the data.

An exmaple of this dictionairy for the NEMO model is::

	coords 	= {
		't':'time_counter', 
		'z':'nav_lev', 
		'lat': 'nav_lat', 
		'lon': 'nav_lon', 
		'cal': '360_day',
		'tdict':bv2tools.tdicts['ZeroToZero']
	}	


The 'tdict' field is a dictionairy mapping the time index onto the array index, Ie some time fields will have the day of year
in the time field so that a 360 day year looks like::
 
  	nc('time')[:] = [15,45,75, ..., 345]


This dictionairy allows the analysis suite to load any model or observational data in the same way, increasing
the re-usability of the overall work package.



Details dictionairies
--------------------------------------

Similar to the coordinates dictionairies, details dictionaries are tools that are built into the times series, level0 and point to point analysis toolkits.
They are dictionairies with very specific contents that allow the analysis package to open, load, and interpret the data fields.
A details dictionairy is required for each model and for each dataset type to be analysed.

	
The contents of this dictionairy is:
	:name: The name of the field, used in plotting, saving fields and elsewhere.
	:vars: A list holding the variable name (or names) as they appear in the netcdf.
	:convert: A function that can be used to either the field in vars and apply a simple numerical function to the field. (see below)
	:units: A string holding the units **after** the function in convert.

An example would be for chlorophull in MEDUSA::

	modeldetails = {
		'name': 	'Chlorophyll', 
		'vars':		['CHN','CHD'], 
		'convert': 	bv2tools.sums,
		'units':	'mg C/m^3'
	} 

The 'name' and 'units' fields are self evident in this exmaple.
The 'vars' are the chlorophyll concentrations of two plankton functional types in MEDUSA.
The sum of these two fields is the total chlorophyll concentration.
The 'convert' function is the sums function in bv2tools, which produces the sum of all the fields  in the 'vars' item.
In this case, the sums function is::

	def sums(nc,keys):		
		a = nc.variables[keys[0]][:]
		for k in keys[1:]:a += nc.variables[k][:]
		return a 



bv2tools holds many other useful functions to fulfil this role, such as:
	:NoChange: Loads keys[0] from the netcdf.
	:mul1000: Loads keys[0] from the netcdf, then multiplies by 1000.
	:div1000: Loads keys[0] from the netcdf, but divides by 1000.
	:sums: Loads keys[0] from the netcdf,  then sums the other keys.
	
It is also possible to write a custom function for the "convert" item of this dictionairy, many exmaples exist in the
analysis_timeseries.py module.


The modules of BGC-val
=======================

Here is a 


The Whole Package
--------------------


.. automodule:: theWholePackage
    :members: theWholePackage


Timeseries Analysis 
--------------------
.. automodule:: analysis_timeseries
    :members: analysis_timeseries



Level 0 Analysis 
--------------------
.. automodule:: analysis_level0
    :members: analysis_level0


Point to point Analysis 
--------------------
.. automodule:: analysis_p2p
    :members: analysis_p2p
    
    
    
.. toctree::
   :maxdepth: 2
   
        
    
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

