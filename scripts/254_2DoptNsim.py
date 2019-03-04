'''
Optimizing and simulating 2D models
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


# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )
#projections is sample size of alleles
#need to use MINIMUM projections

#pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
#projs = [70,      61,     72,     70,    72,    46,     61]

#=================================================================================================#
#											FLORIDA				 								  #
#=================================================================================================#

fl = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLLG','FLCC' ],projections =[70,61] ,polarized = False )  #polarized = False creates folded spectrum

os.chdir("FL2D")



#===================================FIT MODEL TO EMPIRICAL DATA ==================================#

#Make sure to define your extrapolation grid size.
pts = [ 200,220,240 ]

#Provide best optimized parameter set for empirical data.
#These will come from previous analyses you have already completed (above)
emp_params = [1.01,11.6,0.1725,8.5,0.32]

#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True

#Fit the model using these parameters and return the folded model SFS (scaled by theta).
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model function,
#but everything else can stay as it is. See above for argument explanations.
scaled_fl = Optimize_Functions_GOF.Optimize_Empirical(fl, pts, "Empirical", "growth_twoep_sym", growth_twoep_sym_mig, emp_params, fs_folded=fs_folded)


#=======================================PERFORM SIMULATIONS=======================================#

#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = True

#Set the number of simulations to perform here. This should be ~100 or more.
sims = 100

#Enter the number of parameters found in the model to test.
p_num = 5

#Set the number of rounds here.
rounds = 3

#I strongly recommend defining the lists for optional arguments to control the settings 
#of the optimization routine for all the simulated data.
reps = [20,30,50]
maxiters = [5,10,20]
folds = [3,2,1]

#Execute the optimization routine for each of the simulated SFS.
Optimize_Functions_GOF.Perform_Sims(sims, scaled_fl, pts, "growth_twoep_sym", growth_twoep_sym_mig, rounds, p_num, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds)




#=================================================================================================#
#											TEXAS				 								  #
#=================================================================================================#

os.chdir("TX2D")
	
#===================================FIT MODEL TO EMPIRICAL DATA ==================================#

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


#=======================================PERFORM SIMULATIONS=======================================#

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
