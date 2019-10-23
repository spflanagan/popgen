'''
Running 2D model for LAFW vs ALST  pops
Run this from outside the dadi directory
'''

#start with ipython -pylab from ~/Research/popgen/fwsw_results/dadi_analysis

# Numpy is the numerical library dadi is built upon
import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime

#use dportik's functions
#get the optimize functions
execfile("../../programs/dadi_pipeline-master/Two_Population_Pipeline/Optimize_Functions.py")
execfile("../../programs/dadi_pipeline-master/Two_Population_Pipeline/Models_2D.py")
execfile("../../scripts/rougeux_models.py")
#%load ../../programs/dadi_pipeline-master/Two_Population_Pipeline/Optimize_Functions.py
#%load ../../programs/dadi_pipeline-master/Two_Population_Pipeline/Models_2D.py
#%load ../../scripts/250_custom_dadi_models.py

#### SET THESE ####
plot=False
homogeneous=True
heterogeneous=False

######################
# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )
#projections is sample size of alleles
#need to use MINIMUM projections

#pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
#projs = [70,      60,     72,     70,    72,    46,     61]

spect = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'LAFW','ALST' ],projections =[72,70] ,polarized = False )  #polarized = False creates folded spectrum

os.chdir("LA2D")

#=================================================================================================#
#										PLOT SPECTRA	 										  #
#=================================================================================================#
if plot is True:
	dadi.Plotting.plot_single_2d_sfs(spect,vmin=0.00001)


#=================================================================================================#
#										LOOP TO OPTIMIZE 										  #
#=================================================================================================#
if homogeneous is True:
	pts = [ 220,240 ]
	rounds=4
	#define the lists for optional arguments
	#you can change these to alter the settings of the optimization routine
	reps = [10,20,30,40]
	maxiters = [3,5,10,15]
	folds = [3,2,2,1]
	fs_folded = True
	prefix = "la2d"

	for i in range(1,3):
		prefix = "LA2D_{}".format(i)
			
		# Split with complete isolation
		Optimize_Routine(spect, pts, prefix, "SI", SI, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, Ts, O",in_upper=[15,15,10,1],in_lower=[0.1,0.1,0.1,0.00001])

		# Split into two populations with migration during divergence
		Optimize_Routine(spect, pts, prefix, "IM", IM, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, Ts, O",in_upper=[15,15,10,10,10,1],in_lower=[0.1,0.1,0.1,0.01,0.1,0.00001])
		
		# Split into two populations with ancient migration 
		Optimize_Routine(spect, pts, prefix, "AM", AM, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, Ts, Tam, O",in_upper=[15,15,10,10,10,10,1],in_lower=[0.1,0.1,0.1,0.01,0.1,0.1,0.00001])

		# Split into two populations followed by secondary contact 
		Optimize_Routine(spect, pts, prefix, "SC", SC, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, Ts, Tsc, O",in_upper=[15,15,10,10,10,10,1],in_lower=[0.1,0.1,0.1,0.01,0.1,0.1,0.00001])
		

#=================================================================================================#
#							LOOP TO OPTIMIZE HETEROGENEOUS MODELS								  #
#=================================================================================================#
if heterogeneous is True:
	pts = [ 220,240 ]
	rounds=4
	#define the lists for optional arguments
	#you can change these to alter the settings of the optimization routine
	reps = [10,20,30,40]
	maxiters = [3,5,10,15]
	folds = [3,2,2,1]
	fs_folded = True
	prefix = "la2Dhet"

	for i in range(1,2):
		prefix = "AL2Dhet_Number_{}".format(i)
		
		# Split into two populations, with continuous asymmetric migration.
		Optimize_Routine(spect, pts, prefix, "asym_mig_O", asym_mig_O, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, O",in_upper=[15,15,10,10,10,1],in_lower=[.1,.1,0.01,0.01,0.01,0])
		
		# Split into two populations, with continuous asymmetric migration.
		Optimize_Routine(spect, pts, prefix, "asym_mig_2N", asym_mig_2N, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, O, Q, hrf",in_upper=[15,15,10,10,10,1,1,1],in_lower=[.1,.1,0.01,0.01,0.01,0,0,0.0001])

		# Split into two populations, with continuous asymmetric migration and genomic islands.
		Optimize_Routine(spect, pts, prefix, "asym_mig_2m", asym_mig_2m, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, me12, me21, O, P",in_upper=[15,15,10,10,10,10,10,1,1],in_lower=[.1,.1,0.01,0.01,0.01,0.01,0.01,0,0])

		# Split into two populations, with continuous asymmetric migration.
		Optimize_Routine(spect, pts, prefix, "asym_mig_2N2m", asym_mig_2N2m, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, me12, me21, O, P, Q, hrf",in_upper=[15,15,10,10,10,10,10,1,1,1,1],in_lower=[.1,.1,0.01,0.01,0.01,0.01,0.01,0,0,0,0])
