'''
Making plots with filtered data
Run this from ~/Research/popgen/fwsw_results/dadi_results/
'''

#start with ipython -pylab in dadi output directory 

# Numpy is the numerical library dadi is built upon
import sys
import os
import numpy
import dadi
import matplotlib
matplotlib.use('Agg')
import pylab
from datetime import datetime

execfile("../../scripts/rougeux_models.py")
execfile("../../programs/dadi_pipeline-master/Plotting/Plotting_Functions.py")


# get the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )

# create a spectrum
pops = ['FLLG', 'FLCC', 'ALFW', 'ALST', 'LAFW', 'TXFW','TXCC' ]
projs = [  70, 60 , 72, 70, 72, 46, 61 ]
pts = [ 220,240 ]


#Fit the model using these parameters and return the model SFS.
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model function,
#but everything else can stay as it is. See above for argument explanations.


# FLLG-FLCC
prefix="FLLG_FLCC"
vmin_val=float(0.0001)
fs=dadi.Spectrum.from_data_dict(dd , pop_ids = [ "FLLG", "FLCC" ],projections = [ projs[0], projs[1] ],polarized = False )  #polarized = False creates folded spectrum

### best model
emp_params = [100.7513, 0.2839, 0.7258,  0.0679, 0.3924, 4.0179, 0.0985, 11.3563, 0.9802, 0.2197, 0.1026, 0.2420]
model_fit = Fit_Empirical(fs, pts, prefix, "SC2mG", SC2mG, emp_params, fs_folded=True)
dadi.Plotting.plot_2d_comp_multinom(data=fs,model=model_fit,vmin=vmin_val,show=False)
name = "../../figs/dadi/FLLG_FLCC_SC2mG_best.png"
pylab.savefig(name)
pylab.close()

dadi.Plotting.plot_2d_comp_Poisson(data=fs,model=model_fit,vmin=vmin_val,show=False)
name = "../../figs/dadi/FLLG_FLCC_SC2mG_bestPois.png"
pylab.savefig(name)
pylab.close()

### Next best
emp_params = [24.4124, 1.6570, 0.7252,  4.7761, 7.8046, 0.1623, 0.0872, 13.5789, 0.3007, 0.3189, 0.1001, 0.9411]
model_fit2 = Fit_Empirical(fs, pts, prefix, "SC2mG", SC2mG, emp_params, fs_folded=True)
dadi.Plotting.plot_2d_comp_multinom(data=fs,model=model_fit2,vmin=vmin_val,show=False)
name = "../../figs/dadi/FLLG_FLCC_SC2mG_2nd.png"
pylab.savefig(name)
pylab.close()

dadi.Plotting.plot_2d_comp_Poisson(data=fs,model=model_fit2,vmin=vmin_val,show=False)
name = "../../figs/dadi/FLLG_FLCC_SC2mG_2ndPois.png"
pylab.savefig(name)
pylab.close()








