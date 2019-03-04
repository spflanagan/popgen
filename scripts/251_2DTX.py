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
# #=================================================================================================#
# #										PLOT SPECTRA	 										  #
# #=================================================================================================#
# dadi.Plotting.plot_single_2d_sfs(tx,vmin=0.01)


# #=================================================================================================#
# #										LOOP TO OPTIMIZE 										  #
# #=================================================================================================#
# pts = [ 200,220,240 ]
# rounds=4
# #define the lists for optional arguments
# #you can change these to alter the settings of the optimization routine
# reps = [10,20,30,40]
# maxiters = [3,5,10,15]
# folds = [3,2,2,1]
# fs_folded = True
# prefix = "tx"

# for i in range(1,6):
# 	prefix = "V5_Number_{}".format(i)
# 	# Split into two populations, no migration.
# 	Optimize_Routine(tx, pts, prefix, "no_mig", no_mig, rounds, 3, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, T",in_upper=[15,15,10],in_lower=[1,1,0])

# 	# Split into two populations, with continuous symmetric migration.
# 	Optimize_Routine(tx, pts, prefix, "sym_mig", sym_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m, T",in_upper=[15,15,10,10],in_lower=[1,1,0.01,0])

# 	# Split into two populations, with continuous asymmetric migration.
# 	Optimize_Routine(tx, pts, prefix, "asym_mig", asym_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, T",in_upper=[15,15,10,10,10],in_lower=[1,1,0.01,0.01,0])

# 	# Split with no migration, then instantaneous size change with no migration.
# 	Optimize_Routine(tx, pts, prefix, "no_mig_size", no_mig_size, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, T1, T2",in_upper=[15,15,15,15,10,10],in_lower=[1,1,1,1,0,0])

# 	# Split with symmetric migration, then instantaneous size change with continuous symmetric migration.
# 	Optimize_Routine(tx, pts, prefix, "sym_mig_size", sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2",in_upper=[15,15,15,15,10,10,10],in_lower=[1,1,1,1,0.01,0,0])

# 	# Split with different migration rates, then instantaneous size change with continuous asymmetric migration.
# 	Optimize_Routine(tx, pts, prefix, "asym_mig_size", asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2",in_upper=[15,15,15,15,10,10,10,10],in_lower=[1,1,1,1,0.01,0.01,0,0])

# 	# Split with continuous symmetrical gene flow, followed by instantaneous size change with no migration.
# 	Optimize_Routine(tx, pts, prefix, "anc_sym_mig_size", anc_sym_mig_size, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m, T1, T2",in_upper=[15,15,15,15,10,10,10],in_lower=[1,1,1,1,0.01,0,0])

# 	# Split with continuous asymmetrical gene flow, followed by instantaneous size change with no migration.
# 	Optimize_Routine(tx, pts, prefix, "anc_asym_mig_size", anc_asym_mig_size, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1a, nu2a, nu1b, nu2b, m12, m21, T1, T2",in_upper=[15,15,15,15,10,10,10,10],in_lower=[1,1,1,1,0.01,0.01,0,0])

	
#=================================================================================================#
#									FIT MODEL TO EMPIRICAL DATA 								  #
#=================================================================================================#
#Make sure to define your extrapolation grid size.
pts = [ 200,220,240 ]

#Provide best optimized parameter set for empirical data.
#These will come from previous analyses you have already completed (above)
emp_params = [1.01,3.5,1.01,14.5,0.5,0.85,0.4]

#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True

#Fit the model using these parameters and return the folded model SFS (scaled by theta).
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model function,
#but everything else can stay as it is. See above for argument explanations.
scaled_tx = Optimize_Functions_GOF.Optimize_Empirical(tx, pts, "Empirical", "sym_mig_size", sym_mig_size, emp_params, fs_folded=fs_folded)


#=================================================================================================#
#										PERFORM SIMULATIONS 								      #
#=================================================================================================#

#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True

#Set the number of simulations to perform here. This should be ~100 or more.
sims = 100

#Enter the number of parameters found in the model to test.
p_num = 7

#Set the number of rounds here.
rounds = 3

#I strongly recommend defining the lists for optional arguments to control the settings 
#of the optimization routine for all the simulated data.
reps = [20,30,50]
maxiters = [5,10,20]
folds = [3,2,1]

#Execute the optimization routine for each of the simulated SFS.
Optimize_Functions_GOF.Perform_Sims(sims, scaled_tx, pts, "growth_twoep_sym", growth_twoep_sym_mig, rounds, p_num, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds)
