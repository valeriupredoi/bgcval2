#!/usr/bin/python


from ConfigParser import ConfigParser
import os

from collections import defaultdict


package_directory = os.path.dirname(os.path.abspath(__file__))
headerline =  "\n\n####"


def checkConfig(Config,debug=False):
	"""
	If it's a string, it opens it as a ConfigParser.
	"""
	if type(Config) == type('string'):
		if debug: print "Reading", Config
		Config1 = ConfigParser()
		Config1.read(Config)
		Config = Config1
	return Config

def parseList(Config,section,option):
	"""
	This tool loads an string from config file and returns it as a list.
	"""
	Config = checkConfig(Config)
	try:	list1 = Config.get(section, option)
	except: 
		print "No option ",option," in section: ",section
		return ''
		
	list1 = list1.replace(',', ' ')
	list1 = list1.replace('  ', ' ')
	list1 = list1.replace('\'', '')
	list1 = list1.replace('\"', '')	
	while list1.count('  ')>0: 
		list1 = list1.replace('  ', ' ')
	return list1.split(' ')




def addCompare(txt=''):
	txt += headerline
	txt += "\n# analysis_compare"
	txt += "\n1 5   * * * "+package_directory+"/pmlcron_compare.sh >  /users/modellers/ledm/ImagesFromJasmin/cron-compare.log 2>&1"
#        txt += "\n1 21  * * * "+package_directory+"/pmlcron_compare.sh >  /users/modellers/ledm/ImagesFromJasmin/cron-compare.log 2>&1"
	return txt


def loadJobs(cp):
	cp = checkConfig(cp)
	jobs = cp.options('jobs')
	jobs = {j:True for j in jobs} #remove duplicates
	jobs = sorted(jobs.keys())
	jobs.reverse()
	return jobs


def hourcheck(hour):
	if 0 <= hour < 24:
		return hour
	while hour < 0: 
		hour += 24
	while hour > 23:
		hour -= 24
		if hour <0: assert 0
	return hour



def get_times(jobs, startinghour = 19, lasthour = 8, ):

	""" Figure out the times to set the runs.
	returns two dictionaries: minutes and hours
	
	"""
	hours = {}
	minutes = {}
	hours_available = (24 - startinghour ) + lasthour
	if len(jobs) < 	hours_available:
		for i, job in enumerate(jobs):
		    minutes[job] = ' 1'
		    hours[job] = hourcheck(i + startinghour)
		    
	elif len(jobs) <= 2 * hours_available:	
		for i, job in enumerate(jobs):
			hour = hourcheck(i + startinghour)
			if i <= hours_available:
			    minutes[job] = ' 1'
		    	    hours[job] = hour	    
			elif hours_available < i <= 2 * hours_available:
			    minutes[job] = '31'
		    	    hours[job] = hourcheck(hour - hours_available)
		    	else: assert 0
		    	
	elif len(jobs) <= 3 * hours_available:	
		minutes = defaultdict(lambda:' 1')
		for i, job in enumerate(jobs):
			hour = hourcheck(i + startinghour)
			if i <= hours_available:
			    minutes[job] = '1'
		    	    hours[job] = hour			    
			elif hours_available < i <= 2 * hours_available:
			    minutes[job] = '21'
		    	    hours[job] = hourcheck(hour - (1*hours_available))
			elif 2 * hours_available < i <= 3 * hours_available:
			    minutes[job] = '41'
		    	    hours[job] = hourcheck(hour - (2*hours_available))
		    	else: assert 0
		    	    	    			    
	elif len(jobs) <= 4 * hours_available:	
		minutes = defaultdict(lambda:' 1')
		for i, job in enumerate(jobs):
			hour = hourcheck(i + startinghour)
			if i <= hours_available:
			    minutes[job] = '1'
		    	    hours[job] = hour			    
			elif hours_available < i <= 2 * hours_available:
			    minutes[job] = '16'
		    	    hours[job] = hourcheck(hour - (1*hours_available))
			elif 2 * hours_available < i <= 3 * hours_available:
			    minutes[job] = '31'
		    	    hours[job] = hourcheck(hour - (2*hours_available))
			elif 3 * hours_available < i <= 4 * hours_available:
			    minutes[job] = '46'
		    	    hours[job] = hourcheck(hour - (3*hours_available))
		    	else: assert 0		    	    
	else:	
		assert 0	
	return hours, minutes



def addMassjobs(configfn='',startinghour = 17, lasthour = 8,):
	cp = checkConfig(configfn,)
	jobs = loadJobs(cp)
	
	massscripts = ''
	massscripts += headerline	
	massscripts += "\n# DownloadFromMass"	

	hours, minutes = get_times(jobs, startinghour = startinghour, lasthour = lasthour, )
		
	for i, job in enumerate(jobs):
		hour = str(int(hours[job]))
		minute = minutes[job]

		options = 	parseList(cp,'jobs',job)
		if 'standard' in options or 'physics' in options:	
			massscripts	+= '\n'+minute+'  '+hour+\
					'   * * * '+package_directory+'/pmlcron_mass.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron-mass_'+\
					job + '.log 2>&1'
		
		if 'MLD' in options:	
			massscripts	+= '\n'+minute+' '+hour+\
					'   * * * '+package_directory+'/pmlcron_mass_MLD.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron-mass_mld_'+\
					job + '.log 2>&1'
	return massscripts	



def addSci1jobs(configfn='',startinghour = 18, lasthour = 8, ):
	cp = checkConfig(configfn,)
	jobs = loadJobs(cp)
	
	scijobs = ''
	scijobs += headerline
	scijobs += "\n# Run evaluation"	
	
	hours, minutes = get_times(jobs, startinghour = startinghour, lasthour = lasthour, )
	
	for i, job in enumerate(jobs):
		hour = str(int(hours[job]))
		minute = minutes[job]
		options = parseList(cp,'jobs',job)
		if 'standard' in options:	
			scijobs	+= '\n'+ minute+'  '+hour+\
					'   * * * '+package_directory+'/pmlcron_sci1.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron-sci1_'+\
					job + '.log 2>&1'
		
		if 'physics' in options:	
			scijobs	+= '\n'+ minute+'   '+hour+\
					'   * * * '+package_directory+'/pmlcron_sci1_physOnly.sh '+\
					job + ' >  /users/modellers/ledm/ImagesFromJasmin/cron_sci1_phys_'+\
					job + '.log 2>&1'
	return scijobs	


def runNow(configfn=''):
	#####
	# this script 
	cp = checkConfig(configfn,)
	jobs = loadJobs(cp)
	scijobs = ''
	scijobs += "\n# Run evaluation\n"			
	for i, job in enumerate(jobs):
		options = 	parseList(cp,'jobs',job)
		if 'standard' in options:	
			scijobs	+= package_directory+'/pmlcron_mass.sh '+job + '; '
		
		if 'MLD' in options:	
			scijobs	+= package_directory+'/pmlcron_mass_MLD.sh '+job + '; '
			
			
	for i, job in enumerate(jobs):
		options = 	parseList(cp,'jobs',job)
		if 'standard' in options:	
			scijobs	+= package_directory+'/pmlcron_sci1.sh '+job + '; '
		
		if 'physics' in options:	
			scijobs	+= package_directory+'/pmlcron_sci1_physOnly.sh '+job + '; '
		
	return scijobs	
	
		
def main():
	fn = str(package_directory+'/jobids_config.ini')

	crontxt = ''

	crontxt+= headerline+"\nMAILTO=\"\""
	
	crontxt+= headerline+"\n# This crontable was made by: "+str(os.path.abspath(__file__))
	
	crontxt+= headerline+"\n# using the file: "+str(os.path.abspath(fn))
	

	#crontxt += addCompare()
	
	crontxt += addMassjobs(fn)
	
	#crontxt += addSci1jobs(fn)
	
	crontxt+= headerline
	print crontxt
	save =True
	if save:
		cfn = 'crons/crontab.txt'
		f = open(cfn,'w')
		f.write(crontxt)
		f.close()
		print "saved as :", cfn
		print "install with: \ncrontab "+cfn
	
	print runNow(fn)
	
main()
