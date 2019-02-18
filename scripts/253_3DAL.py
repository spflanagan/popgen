'''
Running 3D model for AL and LA pops
Run this from outside the dadi directory
'''

#start with ipython -pylab from ~/Research/popgen/fwsw_results/dadi_analysis/AL3D

# Numpy is the numerical library dadi is built upon
import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime

#use dportik's functions
#get the optimize functions
execfile("../../programs/dadi_pipeline-master/Three_Population_Pipeline/Optimize_Functions.py")
execfile( "../../programs/dadi_pipeline-master/Three_Population_Pipeline/Models_3D.py")


# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )
#projections is sample size of alleles
#need to use MINIMUM projections

#pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
#projs = [70,      61,     72,     70,    72,    46,     61]

al = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'ALFW','ALST','LAFW' ],projections =[72,70,72] ,polarized = False )  #polarized = False creates folded spectrum

os.chdir("AL3D")
#dadi.Plotting.plot_3d_sfs(al,vmin=0.01)


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
prefix = "al"

for i in range(1,6):
	prefix = "V5_Number_{}".format(i)
	# Model with split between pop 1 and (2,3), gene flow does not occur. Period of symmetric secondary contact occurs between adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete. 
	Optimize_Routine(al, pts, prefix, "refugia_adj_1", refugia_adj_1, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, m1, m2, m3, T1, T2, T3")

	# Model with split between pop 1 and (2,3), then split between 2 and 3. Symmetric migration
	Optimize_Routine(al, pts, prefix, "split_symmig_all", split_symmig_all, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2")

	# Model with split between pop 1 and (2,3), then split between 2 and 3
	Optimize_Routine(al, pts, prefix, "split_symmig_adjacent", refugia_adj_1, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2")

	# Model with split between pop 1 and (2,3), with gene flow, which then stops. Split between pops 2 and 3, gene flow does not occur at all.
	Optimize_Routine(al, pts, prefix, "ancmig_adj_3", ancmig_adj_3, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA,m1, m2,  T1a, T1b, T2")
