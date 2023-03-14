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



#### Load .xcm file to get some pn data and the model ####
read_xcm("control_dw.xcm")



def to_gridpoint(x,increment):
	scaling=1./increment
	return np.round(x*scaling)/scaling
	# return 0

#### Draw N parameter combinations
N=1000
# h=np.random.rand(N)*7+3
# h=3
mu=np.random.rand(N)*0.75+0.2
gamma=2

m_dot=np.random.rand(N)*0.4+0.1
fv=np.random.rand(N)*1.25+0.25
lx=np.random.rand(N)*1.5+0.5

mu=to_gridpoint(mu,0.025)
fv=to_gridpoint(fv,0.25)
# print(fv)
# exit()


m=AllModels(1)
# print(m.pds_hres_all.parameterNames)
# exit()

parset_file=open("simulated_spectra_2/dw_control/parameters.txt",'w')

for i in range(0,N):
	parset=[m_dot[i],fv[i],lx[i],gamma,mu[i]]

	print("\nParameter set %s:" % str(i+1))
	print(" ".join([str(x) for x in parset]))

	parset_file.write(str(i)+" "+" ".join([str(x) for x in parset])+"\n")

	filestem=str(i)
	filename="fakespec_"+filestem+"_control_dw.pha"
	backfilename="fakespec_"+filestem+"_control_dw_bkg.pha"
	if not os.path.exists("simulated_spectra_2/dw_control/"+filename):

		# Set Relxill parameters
		# m.relxilllp.a=parset[0]
		# m.relxilllp.h=parset[1]
		# m.relxilllp.Afe=parset[2]

		# Set disk wind parameters
		m.pds_hres_all.Mdot=parset[0]
		m.pds_hres_all.fv=parset[1]
		m.pds_hres_all.LxLEdd=parset[2]

		m.powerlaw.PhoIndex=parset[3]

		m.pds_hres_all.mu=parset[4]
		i=np.arccos(parset[4])/2/np.pi*360
		# m.relxilllp.Incl=i

		# Define settings for FakeIt
		fakeitsettings=FakeitSettings(
									response="pn_1.rmf",\
									arf="pn_1.arf",\
									background="pn_bkg_1.pha",\
									exposure=100000,\
									backExposure=100000,\
									correction=1,\
									fileName=filename
									)

		# Generate fake spectrum
		AllData.fakeit(1,fakeitsettings)


		print('Binning spectrum...')
		binnedfilename='.'.join(filename.split('.')[:-1])+'_binned.pha'
		grppha_str='ftgrouppha\
					infile=%s\
					backfile=%s\
					outfile=%s\
					grouptype=opt\
					respfile=pn_1.rmf' % (filename, backfilename, binnedfilename)
		call(grppha_str, shell=True)


		os.rename(filename,"simulated_spectra_2/dw_control/"+filename)
		os.rename(binnedfilename,"simulated_spectra_2/dw_control/"+binnedfilename)
		os.rename(backfilename,"simulated_spectra_2/dw_control/"+backfilename)


parset_file.flush()
parset_file.close()
