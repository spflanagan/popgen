'''
Running 2D model for FL pops
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
execfile("../../scripts/250_custom_dadi_models.py")
execfile("../../scripts/dadi_selection_mods.py")
#%load ../../programs/dadi_pipeline-master/Two_Population_Pipeline/Optimize_Functions.py
#%load ../../programs/dadi_pipeline-master/Two_Population_Pipeline/Models_2D.py
#%load ../../scripts/250_custom_dadi_models.py

#### SET THESE ####
plot=False
homogeneous=False
heterogeneous=True

######################
# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )
#projections is sample size of alleles
#need to use MINIMUM projections

#pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
#projs = [70,      60,     72,     70,    72,    46,     61]

fl = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLLG','FLCC' ],projections =[70,60] ,polarized = False )  #polarized = False creates folded spectrum

os.chdir("FL2D")

#=================================================================================================#
#										PLOT SPECTRA	 										  #
#=================================================================================================#
if plot is True:
	dadi.Plotting.plot_single_2d_sfs(fl,vmin=0.01)


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
	prefix = "fl"

	for i in range(4,6):
		prefix = "V1_Number_{}".format(i)
			
		# Split into two populations, no migration.
		Optimize_Routine(fl, pts, prefix, "no_mig", no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T",in_params=[0.16,9.79,0.13],in_upper=[15,10,10],in_lower=[0.1,0.1,0.01])

		# Split into two populations, with continuous symmetric migration.
		Optimize_Routine(fl, pts, prefix, "sym_mig", sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T",in_params=[0.18,5.98,0.33,0.18],in_upper=[15,15,5,10],in_lower=[0.1,1,0.01,0.01])
		
		# Split into two populations, with continuous asymmetric migration.
		Optimize_Routine(fl, pts, prefix, "asym_mig", asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T",in_params=[0.15,9.39,0.44,0.06,0.15],in_upper=[15,15,10,10,10],in_lower=[.1,.1,0.01,0.01,0.01])

		# Split into two pops, isolation with migration model
		Optimize_Routine_Extrap(fl, pts, prefix, "IM", dadi.Demographics2D.IM, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "s, nu1, nu2, T, m12, m21",in_params=[1.83,0.41,5.36,1.01,0.16,5.05],in_upper=[10,15,15,10,10,10],in_lower=[0.01,0.1,.1,.1,0.01,0.01])
		
		# Split into two pops, growth in one pop and two epoch in another
		Optimize_Routine(fl, pts, prefix, "growth_twoep_no", growth_twoep_no_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2,T,Tc",in_params=[1.01,4.99,0.55,0.74],in_upper=[15,5,10,1],in_lower=[1,1,0.01,0.01])

		# Split into two pops, growth in one pop and two epoch in another
		Optimize_Routine(fl, pts, prefix, "growth_twoep_sym", growth_twoep_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2,m, T,Tc",in_params=[1.01,11.39,0.17,7.76,0.39],in_upper=[15,15,10,10,1],in_lower=[1,1,0.01,0.01,0.01])

		# Split into two pops, growth in one pop and two epoch in another
		Optimize_Routine(fl, pts, prefix, "growth_twoep_asym", growth_twoep_asym_mig, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2,m12,m21, T,Tc",in_params=[1.01,12.13,0.06,0.14,0.42,0.63],in_upper=[15,15,10,10,10,1],in_lower=[1,1,0.01,0.01,0.01,0.01])
		
		"""
		# Split with no migration, then instantaneous size change with no migration.
		Optimize_Routine(fs, pts, prefix, "no_mig_size", no_mig_size, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, T1, T2")

		# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
		Optimize_Routine(fs, pts, prefix, "sym_mig_size", sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")

		# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
		Optimize_Routine(fs, pts, prefix, "asym_mig_size", asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")

		# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.
		Optimize_Routine(fs, pts, prefix, "anc_sym_mig_size", anc_sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2")

		# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
		Optimize_Routine(fs, pts, prefix, "anc_asym_mig_size", anc_asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2")
		"""

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
	prefix = "flhet"

	for i in range(2,6):
		prefix = "V1_Number_{}".format(i)
		
		# Split into two populations, with continuous asymmetric migration.
		Optimize_Routine(fl, pts, prefix, "asym_mig_O", asym_mig_O, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, O",in_upper=[15,15,10,10,10,1],in_lower=[.1,.1,0.01,0.01,0.01,0])
		
		# Split into two populations, with continuous asymmetric migration.
		Optimize_Routine(fl, pts, prefix, "asym_mig_2N", asym_mig_2N, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, O, Q, hrf",in_upper=[15,15,10,10,10,1,1,1],in_lower=[.1,.1,0.01,0.01,0.01,0,0,0.0001])

		# Split into two populations, with continuous asymmetric migration and genomic islands.
		Optimize_Routine_Extrap(fl, pts, prefix, "asym_mig_2m", asym_mig_2m, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, me12, me21, O, P",in_upper=[15,15,10,10,10,10,10,1,1],in_lower=[.1,.1,0.01,0.01,0.01,0.01,0.01,0,0])

		# Split into two populations, with continuous asymmetric migration.
		Optimize_Routine(fl, pts, prefix, "asym_mig_2N2m", asym_mig_2N2m, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T, m12, m21, me12, me21, O, P, Q, hrf",in_upper=[15,15,10,10,10,10,10,1,1,1,1],in_lower=[.1,.1,0.01,0.01,0.01,0.01,0.01,0,0,0,0])
