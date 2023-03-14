#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from xspec import *
import os
import itertools
from glob import glob

os.environ['RELLINE_TABLES']='/home/mlparker/programs/xspec_models/relxill/'
AllModels.lmod('relxill','/home/mlparker/programs/xspec_models/relxill/')
# AllModels.lmod('absmodel','/home/mlparker/programs/xspec/tbnew/')

Xset.chatter=0
Xset.parallel.leven=8
Xset.parallel.error=8

# Plot.device='/xw'

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




Fit.query='y'

### Relxill only fits:
print("Fitting with relxill, no disk wind:")


# parsets=itertools.product(a,h,A_Fe,gamma,mu)

read_xcm("control_ref.xcm")
m=AllModels(1)
os.chdir("simulated_spectra_2/ref_control")



if not os.path.exists("fits_relxill_control_3to10.dat"):
	outputfile1=open('fits_relxill_control_3to10.dat','w')
	titlerow="filename,a,frac,Afe,gamma,mu,fit_a,fit_frac,fit_Afe,fit_gamma,fit_i,chi2,dof\n"
	# titlerow="filename,a,h,Afe,mdot,fv,lx,gamma,mu,fit_a,fit_h,fit_Afe,fit_gamma,fit_i,chi2,dof\n"
	outputfile1.write(titlerow)
	prev_runs=[]

else:
	prev_runs=np.loadtxt('fits_relxill_control_3to10.dat',skiprows=1,delimiter=',',dtype='str')[:,0]
	outputfile1=open('fits_relxill_control_3to10.dat','a')


# outputfile1=open("fits_relxill_control_3to10.dat",'w')
# titlerow="filename,a,h,Afe,gamma,mu,fit_a,fit_h,fit_Afe,fit_gamma,fit_i,chi2,dof\n"
# outputfile1.write(titlerow)

# spectra=glob("fakespec_*_control_ref.pha")

N=1000
parsets=np.loadtxt("parameters.txt")

for i,p in enumerate(parsets):
	parset=p[1:]
	# print(parset)
	# exit()
	parsetstr=','.join([str(x) for x in parset])
	AllData.clear()
	print("Parameter set %s:" % str(i+1))
	print(" ".join([str(x) for x in parset]))
	filestem=str(i)
	filename="fakespec_"+filestem+"_control_ref_binned.pha"

	if filename not in prev_runs:
		# print(filename)
		# filename=spectrum
		print(filename)
		# backfilename="fakespec_"+filestem+"_control_ref_bkg.pha"
		# backfilename=filename.replace(".pha","_bkg.pha")
		AllData(filename)
		AllData.ignore("0.-3.,10.-**")

		# Set Relxill parameters
		m.relxill.a=parset[0]
		m.relxill.a.frozen=True
		m.relxill.refl_frac=parset[1]
		m.relxill.refl_frac.frozen=True
		m.relxill.Afe=parset[2]
		m.relxill.Afe.frozen=True

		# # Set disk wind parameters
		# m.pds_hres_all.Mdot=parset[3]
		# m.pds_hres_all.fv=parset[4]
		# m.pds_hres_all.LxLEdd=parset[5]

		m.relxill.gamma=parset[3]
		# m.relxilllp.gamma.frozen=True

		# m.pds_hres_all.mu=parset[7]
		i=  np.arccos(parset[4])/2/np.pi*360
		m.relxill.Incl=i
		m.relxill.Incl.frozen=True
		m.relxill.norm="1 -1"

		print("Running initial fit...")
		Fit.perform()

		print("Running secondary fit...\n")
		m.relxill.a.frozen=False
		m.relxill.Afe.frozen=False
		m.relxill.Incl.frozen=False
		# m.relxill.refl_frac.frozen=False
		# m.relxilllp.h.frozen=False
		Fit.perform()
		# Plot("ld rat")
		# print(Fit.statistic,Fit.dof)
		# exit()
		# print("Running errors...")
		# Fit.error("1. 1-16")
		# print("Running secondary fit...")
		# Fit.perform()

		outputstring=filename+','+parsetstr+','+\
								','.join([str(m.relxill.a.values[0]),\
								str(m.relxill.refl_frac.values[0]),\
								str(m.relxill.Afe.values[0]),\
								str(m.relxill.gamma.values[0]),\
								str(m.relxill.Incl.values[0]),\
								str(Fit.statistic),\
								str(Fit.dof)+'\n'
								])
		# print(outputstring)
		outputfile1.write(outputstring)
		outputfile1.flush()
	else:
		print('file name',filename,'already fit, skipping')

outputfile1.close()

os.chdir("..")

# ### DW only fits:
# print("Fitting with disk wind, no reflection:")
#
#
# parsets=itertools.product(m_dot,fv,lx,gamma,mu)
#
# read_xcm("control_dw.xcm")
# m=AllModels(1)
# os.chdir("simulated_spectra")
#
# outputfile2=open("fits_dw_control_3to10.dat",'w')
# titlerow="filename,mdot,fv,lx,gamma,mu,fit_mdot,fit_fv,fit_lx,fit_mu,fit_gamma,chi2,dof\n"
# outputfile2.write(titlerow)
#
# for parset in parsets:
# 	parsetstr=','.join([str(x) for x in parset])
# 	AllData.clear()
# 	print("Parameter set:")
# 	print(" ".join([str(x) for x in parset]))
# 	filestem="_".join([str(x) for x in parset])
# 	filename="fakespec_"+filestem+"_control_dw.pha"
# 	backfilename="fakespec_"+filestem+"_control_dw_bkg.pha"
# 	AllData(filename)
# 	AllData.ignore("0.-3.,10.-**")
#
# 	# Set Relxill parameters
# 	# m.relxilllp.a=parset[0]
# 	# m.relxilllp.h=parset[1]
# 	# m.relxilllp.Afe=parset[2]
#
# 	# Set disk wind parameters
# 	m.pds_hres_all.Mdot=parset[0]
# 	m.pds_hres_all.fv=parset[1]
# 	m.pds_hres_all.fv.frozen=True
# 	m.pds_hres_all.LxLEdd=parset[2]
#
# 	m.powerlaw.PhoIndex=parset[3]
#
# 	m.pds_hres_all.mu=parset[4]
# 	m.pds_hres_all.mu.frozen=True
# 	# i=  np.arccos(parset[7])/2/np.pi*360
# 	# m.relxilllp.Incl=i
# 	# m.relxilllp.Incl.frozen=False
# 	m.powerlaw.norm="1 -1"
#
# 	print("Running initial fit...")
# 	Fit.perform()
# 	# print(Fit.statistic,Fit.dof)
# 	# exit()
# 	# print("Running errors...")
# 	# Fit.error("1. 1-16")
# 	# print("Running secondary fit...")
# 	# Fit.perform()
#
# 	outputstring=filename+','+parsetstr+','+\
# 							','.join([str(m.pds_hres_all.Mdot.values[0]),\
# 							str(m.pds_hres_all.fv.values[0]),\
# 							str(m.pds_hres_all.LxLEdd.values[0]),\
# 							str(m.pds_hres_all.mu.values[0]),\
# 							str(m.powerlaw.PhoIndex.values[0]),\
# 							str(Fit.statistic),\
# 							str(Fit.dof)+'\n'
# 							])
# 	print(outputstring)
# 	outputfile2.write(outputstring)
# 	outputfile2.flush()
#
# outputfile2.close()
