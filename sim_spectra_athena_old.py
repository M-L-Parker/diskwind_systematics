#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from xspec import *
import os
import itertools
from subprocess import call

os.environ['RELLINE_TABLES']='/home/mlparker/programs/xspec_models/relxill/'
AllModels.lmod('relxill','/home/mlparker/programs/xspec_models/relxill/')
# AllModels.lmod('absmodel','/home/mlparker/programs/xspec/tbnew/')
cwd=os.getcwd()

Xset.chatter=0

def read_xcm(xcmfilename):
	'''Read data, model etc. from an xspec .xcm file
	Note that this is not complete. Check the list below.
	'''

	cwd = os.getcwd()

	print('\nReading .xcm file...')
	xcmfile = open(xcmfilename)

	model_index = -1

	pars = {}
	newpars = {}
	n_datasets=1
	for i, row in enumerate(xcmfile):
		temp = row.strip().split()

		# print row
		if len(temp) > 0:
			# Change directory
			if temp[0] == 'cd':
				os.chdir(temp[1])

			# Set abundances and cross-sections
			elif temp[0] == 'abund':
				Xset.abund = temp[1]
			elif temp[0] == 'xsect':
				Xset.xsect = temp[1]

			# Load data
			elif temp[0] == 'data':
				print('Loading data:', ' '.join(temp[1:]))
				AllData(' '.join(temp[1:]))
				n_datasets = int(temp[1][0])

			# Ignore specified channels
			elif temp[0] == 'ignore':
				AllData.ignore(' '.join(temp[1:]))

			# Find model definition
			elif temp[0] == 'model':
				print('Loading model:', ' '.join(temp[1:]))
				model_index = i
				modelstr = row.strip()[7:] #7
				model = Model(modelstr)

			elif temp[0] == "bayes":
				pass

			elif model_index > 0 and i > model_index:
				if temp[0] == 'newpar':
					parnum = int(temp[1])
					newpars[parnum] = ' '.join(temp[2:])
				else:
					pars[i - model_index] = ' '.join(temp)

	n_pars = len(pars)
	if n_pars % n_datasets != 0:
		print(n_pars, n_datasets, n_pars % n_datasets)
		print("Number of parameters should be a multiple of the number of datasets. This isn't.")

	else:
		n_pars_per_dset = n_pars / n_datasets
		for p in pars:
			# print(p, pars[p])
			mnum = int(1 + (p - 1) / n_pars_per_dset)
			pnum = int(p - (mnum - 1) * n_pars_per_dset)# -1 # added -1
			# print(mnum, pnum)
#             print '{0} {1}'.format('p', p)
#             print '{0} {1}'.format('pars', pars)
#             print '{0} {1}'.format('pnum', pnum)
			AllModels(mnum).setPars({pnum: pars[p]})
			# print(pnum, pars[p])
		for p in newpars:
			# print p, newpars[p]
			mnum = int(1 + (p - 1) / n_pars_per_dset)
			pnum = int(p - (mnum - 1) * n_pars_per_dset)

			AllModels(mnum).setPars({pnum: newpars[p]})

	return pars, n_datasets



#### Define parameter ranges ####

## Reflection parameters ##
a=[0.,0.5,0.98]
h=[3,6,10]
A_Fe=[1,3,5]

## Disk wind parameters ##
m_dot=[0.1,0.3,0.5]
fv=[0.5,1,1.5]   # Needs to be set to a grid point for interpolation reasons
lx=[0.5,1,2]

## Joint parameters ##
gamma=[2]
mu=[0.875,0.725,0.475]   # Needs to be set to a grid point for interpolation reasons
# i=[np.arccos(x)/2/np.pi*360 for x in mu]    # arccos(mu) in degrees






#### Load .xcm file to get some pn data and the model ####
read_xcm("baseline_hres.xcm")



#### Parameter combinations (for skipping nested for loops)
parsets=itertools.product(a,h,A_Fe,m_dot,fv,lx,gamma,mu)

m=AllModels(1)
# print(m.pds_hres_all.parameterNames)
# exit()
for parset in parsets:
	print("Parameter set:")
	print(" ".join([str(x) for x in parset]))
	filestem="_".join([str(x) for x in parset])
	filename="fakespec_xifu_"+filestem+".pha"
	backfilename="fakespec_xifu_"+filestem+"_bkg.pha"
	AllData.clear()
	if not os.path.exists("simulated_spectra/"+filename):

		# Set Relxill parameters
		m.relxilllp.a=parset[0]
		m.relxilllp.h=parset[1]
		m.relxilllp.Afe=parset[2]

		# Set disk wind parameters
		m.pds_hres_all.Mdot=parset[3]
		m.pds_hres_all.fv=parset[4]
		m.pds_hres_all.LxLEdd=parset[5]

		m.relxilllp.gamma=parset[6]

		m.pds_hres_all.mu=parset[7]
		i=np.arccos(parset[7])/2/np.pi*360
		m.relxilllp.Incl=i

		# Define settings for FakeIt
		fakeitsettings=FakeitSettings(
									response="XIFU_CC_BASELINECONF_2018_10_10.rmf",\
									arf="XIFU_CC_BASELINECONF_2018_10_10.arf",\
									background="Total_pointsources_XIFU_CC_BASELINECONF_2018_10_10.pha",\
									exposure=100000,\
									backExposure=100000,\
									correction=1,\
									fileName=filename
									)

		# Generate fake spectrum
		AllData.fakeit(1,fakeitsettings)

		# ##### Need to rebin spectrum here ####
		# binned_file=specgroup_call(filename)
		# Can skip this if we just copy the binning from original file, much faster

		os.rename(filename,"simulated_spectra/"+filename)
		# os.rename(binned_file,"simulated_spectra/"+binned_file)
		os.rename(backfilename,"simulated_spectra/"+backfilename)

	binnedfilename='.'.join(filename.split('.')[:-1])+'_grp25.pha'
	print(binnedfilename)
	if not os.path.exists('simulated_spectra/'+binnedfilename):
		print('./grppha_quick.sh %s %s' % ('simulated_spectra/'+filename, 'simulated_spectra/'+binnedfilename))
		call('./grppha_quick.sh %s %s' % ('simulated_spectra/'+filename, 'simulated_spectra/'+binnedfilename),shell=True)
