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
#%load ../../programs/dadi_pipeline-master/Two_Population_Pipeline/Optimize_Functions.py
#%load ../../programs/dadi_pipeline-master/Two_Population_Pipeline/Models_2D.py
#%load ../../scripts/250_custom_dadi_models.py


# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )
#projections is sample size of alleles
#need to use MINIMUM projections

#pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
#projs = [70,      61,     72,     70,    72,    46,     61]

fl = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLLG','FLCC' ],projections =[70,61] ,polarized = False )  #polarized = False creates folded spectrum

dadi.Plotting.plot_single_2d_sfs(fl,vmin=0.01)


#=================================================================================================#
#										LOOP TO OPTIMIZE 										  #
#=================================================================================================#
pts = [ 200,220,240 ]
rounds=4
#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]
fs_folded = True
prefix = "fl"

for i in range(1,6):
	prefix = "V5_Number_{}".format(i)
	"""	
	# Split into two populations, no migration.
	Optimize_Routine(fl, pts, prefix, "no_mig", no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T",in_upper=[15,5,10],in_lower=[1,1,0])

	# Split into two populations, with continuous symmetric migration.
	Optimize_Routine(fl, pts, prefix, "sym_mig", sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T",in_upper=[15,5,5,10],in_lower=[1,1,0,0])
	
	# Split into two populations, with continuous asymmetric migration.
	Optimize_Routine(fl, pts, prefix, "asym_mig", asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T",in_upper=[15,5,10,10,10],in_lower=[1,1,0,0,0])

	# Split into two pops, isolation with migration model
	Optimize_Routine(fl, pts, prefix, "IM", dadi.Demographics2D.IM, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "s, nu1, nu2, T, m12, m21",in_upper=[10,15,5,10,10,10],in_lower=[0,0.1,1,1,0,0])
	"""
	# Split into two pops, growth in one pop and two epoch in another
	Optimize_Routine(fl, pts, prefix, "growth_twoep_no", growth_twoep_no_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2,T,Tc",in_upper=[15,5,10,1],in_lower=[1,1,0,0])

	# Split into two pops, growth in one pop and two epoch in another
	Optimize_Routine(fl, pts, prefix, "growth_twoep_sym", growth_twoep_sym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2,m, T,Tc",in_upper=[15,15,10,10,1],in_lower=[1,1,0,0,0])

	# Split into two pops, growth in one pop and two epoch in another
	Optimize_Routine(fl, pts, prefix, "growth_twoep_asym", growth_twoep_asym_mig, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2,m12,m21, T,Tc",in_upper=[15,15,10,10,10,1],in_lower=[1,1,0,0,0,0])
	
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

	