'''
Running 2D model for TX pops
Run this from outside the dadi directory
'''

#start with ipython -pylab from ~/Research/popgen/fwsw_results/dadi_analysis/TX2D

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
execfile( "../../programs/dadi_pipeline-master/Two_Population_Pipeline/Models_2D.py")
execfile("../../scripts/250_custom_dadi_models.py")


# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )
#projections is sample size of alleles
#need to use MINIMUM projections

#pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
#projs = [70,      61,     72,     70,    72,    46,     61]

tx = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'TXFW','TXCC' ],projections =[46,61] ,polarized = False )  #polarized = False creates folded spectrum

os.chdir("TX2D")
#=================================================================================================#
#										PLOT SPECTRA	 										  #
#=================================================================================================#
dadi.Plotting.plot_single_2d_sfs(tx,vmin=0.01)


#=================================================================================================#
#										LOOP TO OPTIMIZE 										  #
#=================================================================================================#
pts = [ 220,240 ]
rounds=4
#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]
fs_folded = True
prefix = "tx"

for i in range(1,2):
	prefix = "V1_Number_{}".format(i)
	# Split into two populations, no migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "no_mig", no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T",in_params=[1.01,9.73,0.25],in_upper=[15,15,10],in_lower=[1,1,0.01])

	# Split into two populations, with continuous symmetric migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "sym_mig", sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T",in_params=[1.01,9.56,0.55,0.95],in_upper=[15,15,10,10],in_lower=[1,1,0.01,0.01])

	# Split into two populations, with continuous asymmetric migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "asym_mig", asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T",in_params=[1.01,9.22,0.33,1.37,0.88],in_upper=[15,15,10,10,10],in_lower=[1,1,0.01,0.01,0.01])

	# Split with no migration, then instantaneous size change with no migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "no_mig_size", no_mig_size, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, T1, T2",in_params=[1.1,7.57,1.01,9.29,0.05,0.2],in_upper=[15,15,15,15,10,10],in_lower=[1,1,1,1,0.01,0.01])

	# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "sym_mig_size", sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2",in_params=[12.03,1.65,1.01,9.47,0.58,0.04,0.83],in_upper=[15,15,15,15,10,10,10],in_lower=[1,1,1,1,0.01,0.01,0.01])

	# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "asym_mig_size", asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2",in_params=[1.01,8.90,1.01,9.92,0.63,0.51,0.16,1.15],in_upper=[15,15,15,15,10,10,10,10],in_lower=[1,1,1,1,0.01,0.01,0.01,0.01])

	# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "anc_sym_mig_size", anc_sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2",in_params=[1.48,6.91,1.01,13.94,6.2,1,0.21],in_upper=[15,15,15,15,10,10,10],in_lower=[1,1,1,1,0.01,0.01,0.01])

	# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
	Optimize_Routine_Extrap(tx, pts, prefix, "anc_asym_mig_size", anc_asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2",in_params=[4.29,5.86,1.01,14.93,1.11,8.58,1.0,0.23],in_upper=[15,15,15,15,10,10,10,10],in_lower=[1,1,1,1,0.01,0.01,0.01,0.01])

