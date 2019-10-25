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

#### GET ARGUMENTS ####
snps_file=sys.argv[1]
popi=sys.argv[2]
prji=int(sys.argv[3])
popj=sys.argv[4]
prjj=int(sys.argv[5])
x=int(sys.argv[6])
y=int(sys.argv[7])

print prji
print prjj
######################
# Load the data
dd = dadi.Misc.make_data_dict ( snps_file )
#projections is sample size of alleles
#need to use MINIMUM projections

spect = dadi.Spectrum.from_data_dict(dd , pop_ids =[ popi, popj ],projections =[ prji, prjj ] ,polarized = False )  #polarized = False creates folded spectrum

pathname=popi+"_"+popj
if not os.path.exists(pathname):
    os.makedirs(pathname)
os.chdir(pathname)

#=================================================================================================#
#										PLOT SPECTRA	 										  #
#=================================================================================================#
#figname=pathname+".png"
#dadi.Plotting.plot_single_2d_sfs(spect,vmin=0.00001)
#pylab.savefig(name)


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
prefix = pathname

for i in range(x,y):
	prefix = prefix+"_{}".format(i)
		
	# Split with complete isolation
	Optimize_Routine(spect, pts, prefix, "SI", SI, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, Ts, O",in_upper=[100, 100, 10, 0.99],in_lower=[0.01, 0.01, 0, 0.01])

	# Split into two populations with migration during divergence
	Optimize_Routine(spect, pts, prefix, "IM", IM, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, Ts, O",in_upper=[120, 120, 40, 40, 10, 0.99],in_lower= [0.01, 0.01, 0, 0, 0, 0.01])
	
	# Split into two populations with ancient migration 
	Optimize_Routine(spect, pts, prefix, "AM", AM, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, Ts, Tam, O",in_upper=[15,15,10,10,10,10,1],in_lower=[0.1,0.1,0.1,0.01,0.1,0.1,0.00001])

	# Split into two populations followed by secondary contact 
	Optimize_Routine(spect, pts, prefix, "SC", SC, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, Ts, Tsc, O",in_upper=[15,15,10,10,10,10,1],in_lower=[0.1,0.1,0.1,0.01,0.1,0.1,0.00001])
	
	
	#================================ HETEROGENEOUS MODELS =====================================#

	# SI derivatives
	Optimize_Routine(spect, pts, prefix, "SI2N", SI2N, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, Ts, nr, bf, O",in_upper=[100, 100, 12, 1, 0.999 , 0.99],in_lower=[0.01, 0.01, 0, 0, 0.001, 0.01])
	Optimize_Routine(spect, pts, prefix, "SIG", SIG, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, Ts, O",in_upper=[80, 80, 100, 100, 10, 0.99],in_lower=[0.01, 0.01, 0.01, 0.01, 0, 0.01])
	Optimize_Routine(spect, pts, prefix, "SI2NG", SI2NG, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, hrf, Ts, Q, O",in_upper=[100, 100, 100, 100, 1, 10, 0.5, 0.99],in_lower=   [0.01, 0.01, 0.01, 00.1, 0.01, 0.01, 0.01, 0.01])
	
	# IM derivatives
	Optimize_Routine(spect, pts, prefix, "IMG", IMG, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, m12, m21, Ts, O",in_upper=[120, 120, 100, 100, 40, 40, 10, 0.99],in_lower=[0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0.01])
	Optimize_Routine(spect, pts, prefix, "IM2N", IM2N, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, hrf, m12, m21, Ts, Q, O",in_upper=[120, 120, 1, 40, 40, 10, 0.5, 0.99],in_lower=[0.01, 0.01, 0.1, 0, 0, 0, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "IM2m", IM2m, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, me12, me21, Ts, P, O",in_upper=[250, 250, 150, 150, 45, 45, 10, 0.95, 0.99],in_lower= [0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.01])
	Optimize_Routine(spect, pts, prefix, "IM2NG", IM2NG, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, hrf, m12, m21, Ts, Q, O",in_upper=[120, 120, 80, 80, 1, 40, 40, 10, 0.5, 0.99],in_lower=[0.01, 0.01, 0.01, 0.01, 0.1, 0, 0, 0, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "IM2mG", IM2mG, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, P, O",in_upper=[120, 120, 200, 200, 40, 40, 30, 30, 10, 0.95, 0.99],in_lower= [0.01, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.01])
	
	# AM derivatives
	Optimize_Routine(spect, pts, prefix, "AM2N", AM2N, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, hrf, m12, m21, Tam, Ts, Q, O",in_upper=[80, 80, 1, 30, 30, 10, 10, 0.5, 0.99],in_lower=[0.01, 0.01, 0.1, 0, 0, 0, 0, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "AMG", AMG, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, m12, m21, Tam, Ts, O",in_upper=[120, 120, 100, 100, 30, 30, 10, 2, 0.99],in_lower=[0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0.01])
	Optimize_Routine(spect, pts, prefix, "AM2m", AM2m, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, me12, me21, Ts, Tam, P, O",in_upper=[250, 250, 150, 150, 25, 25, 10, 2, 0.95, 0.99],in_lower=[0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01])
	Optimize_Routine(spect, pts, prefix, "AM2NG", AM2NG, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, hrf, m12, m21, Tam, Ts, Q, O",in_upper=[90, 90, 60, 60, 1, 50, 50, 10, 10, 0.5, 0.99],in_lower=[0.01, 0.01, 0.1, 0, 0, 0, 0, 0, 0, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "AM2N2m", AM2N2m, rounds, 12, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q, O",in_upper=[100, 100, 1, 30, 30, 20, 20, 5, 10, 0.95, 0.5, 0.99],in_lower=[0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "AM2mG", AM2mG, rounds, 12, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, m12, m21, me12, me21, Tam, Ts, P, O",in_upper=[90, 90, 40, 40, 30, 30, 20, 20, 2, 10, 0.95, 0.99],in_lower=[0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01])
	Optimize_Routine(spect, pts, prefix, "AM2N2mG", AM2N2mG, rounds, 14, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q, O",in_upper=[90, 90, 50, 50, 1, 30, 30, 20, 20, 5, 10, 0.95, 0.5, 0.99],in_lower= [0.01, 0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01, 0.01])
	
	# SC derivatives
	Optimize_Routine(spect, pts, prefix, "SCG", SCG, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, m12, m21, Ts, Tsc, O",in_upper=[120, 120, 150, 150, 30, 30, 10, 2, 0.99],in_lower=[0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.01])
	Optimize_Routine(spect, pts, prefix, "SC2N", SC2N, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, hrf, m12, m21, Ts, Tsc, Q, O",in_upper=[120, 120, 1, 80, 80, 10, 2, 0.5, 0.99],in_lower=[0.01, 0.01, 0.1, 0, 0, 0, 0, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "SC2m", SC2m, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P, O",in_upper=[250, 250, 150, 150, 30, 30, 10, 2, 0.95, 0.99],in_lower= [0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01])
	Optimize_Routine(spect, pts, prefix, "SC2NG", SC2NG, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, hrf, m12, m21, Ts, Tsc, Q, O",in_upper=[120, 120, 150, 150, 1, 80, 80, 10, 8, 0.5, 0.99],in_lower=[0.01, 0.01, 0.01, 0.01, 0.1, 0, 0, 0, 0, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "SC2N2m", SC2N2m, rounds, 12, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q, O",in_upper=[120, 120, 1, 80, 80, 10, 10, 10, 2, 0.95, 0.5, 0.99],in_lower=[0.01, 0.01, 0.1, 0, 0, 0, 0, 0, 0, 0.5, 0.01, 0.01])
	Optimize_Routine(spect, pts, prefix, "SC2mG", SC2mG, rounds, 12, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, Tsc, P, O",in_upper=[120, 120, 150, 150, 80, 80, 10, 10, 10, 2, 0.95, 0.99],in_lower=[0.01, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.05, 0.01])
	Optimize_Routine(spect, pts, prefix, "SC2N2mG", SC2N2mG, rounds, 14, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q, O",in_upper= [120, 120, 200, 200, 1, 80, 80, 10, 10, 10, 5, 0.95, 0.5, 0.99],in_lower= [0.01, 0.01, 0.1, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0.5, 0.01, 0.01])