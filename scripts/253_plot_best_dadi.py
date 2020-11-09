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

execfile("../../scripts/rougemont_models_folded.py")
execfile("../../programs/dadi_pipeline-master/Plotting/Plotting_Functions.py")

###########################################
## Generic function
###########################################

# get the info from arguments
pop1=sys.argv[1]
pop2=sys.argv[2]
model=sys.argv[3]
emp_params = [] 
for i in range(4,len(sys.argv)):    
    emp_params.append(float(sys.argv[i]))

print( "Making plot for "+pop1+"vs"+pop2+", model "+model+".")
print("Empirical parameters are: ")
print(emp_params)


# create a spectrum
pops = ['FLLG', 'FLCC', 'ALFW', 'ALST', 'LAFW', 'TXFW','TXCC' ]
projs = [  30, 30 , 30, 30, 30, 30, 30 ]
pts = [ 220,240 ]
vmin_val=float(0.0001)


# setup
prefix=pop1+"_"+pop2
pop1_index=pops.index(pop1)
pop2_index=pops.index(pop2)

# get the data
dd = dadi.Misc.make_data_dict ( prefix+".dadi.snps" )
fs=dadi.Spectrum.from_data_dict(dd , pop_ids = [ pop1, pop2 ],projections = [ projs[pop1_index], projs[pop2_index] ],polarized = False )  
model_fit = Fit_Empirical(fs, pts, prefix, model, globals()[model], emp_params, fs_folded=True)


# plot multinomial
dadi.Plotting.plot_2d_comp_multinom(data=fs,model=model_fit,vmin=vmin_val,show=False)
name = "../../figs/dadi/"+pop1+"_"+pop2+"_"+model+"multinom.png"
pylab.savefig(name)
pylab.close()

# plot poisson
dadi.Plotting.plot_2d_comp_Poisson(data=fs,model=model_fit,vmin=vmin_val,show=False)
name = "../../figs/dadi/"+pop1+"_"+pop2+"_"+model+"_pois.png"
pylab.savefig(name)
pylab.close()

###########################################


#Fit the model using these parameters and return the model SFS.
#Here, you will want to change the model arguments to match your model function,
#and use the empirical parameters from the best-fitting model,
#but everything else can stay as it is. See above for argument explanations.


# # FLLG-FLCC
# prefix="FLLG_FLCC"
# vmin_val=float(0.0001)
# fs=dadi.Spectrum.from_data_dict(dd , pop_ids = [ "FLLG", "FLCC" ],projections = [ projs[0], projs[1] ],polarized = False )  #polarized = False creates folded spectrum

# ### best model
# emp_params = [100.7513, 0.2839, 0.7258,  0.0679, 0.3924, 4.0179, 0.0985, 11.3563, 0.9802, 0.2197, 0.1026, 0.2420]
# model_fit = Fit_Empirical(fs, pts, prefix, "SC2mG", SC2mG, emp_params, fs_folded=True)
# dadi.Plotting.plot_2d_comp_multinom(data=fs,model=model_fit,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_FLCC_SC2mG_best.png"
# pylab.savefig(name)
# pylab.close()

# dadi.Plotting.plot_2d_comp_Poisson(data=fs,model=model_fit,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_FLCC_SC2mG_bestPois.png"
# pylab.savefig(name)
# pylab.close()

# ### Next best
# emp_params = [24.4124, 1.6570, 0.7252,  4.7761, 7.8046, 0.1623, 0.0872, 13.5789, 0.3007, 0.3189, 0.1001, 0.9411]
# model_fit2 = Fit_Empirical(fs, pts, prefix, "SC2mG", SC2mG, emp_params, fs_folded=True)
# dadi.Plotting.plot_2d_comp_multinom(data=fs,model=model_fit2,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_FLCC_SC2mG_2nd.png"
# pylab.savefig(name)
# pylab.close()

# dadi.Plotting.plot_2d_comp_Poisson(data=fs,model=model_fit2,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_FLCC_SC2mG_2ndPois.png"
# pylab.savefig(name)
# pylab.close()


# # FLLG-ALFW
# prefix="FLLG_ALFW"
# vmin_val=float(0.0001)
# fs=dadi.Spectrum.from_data_dict(dd , pop_ids = [ "FLLG", "ALFW" ],projections = [ projs[0], projs[1] ],polarized = False )  #polarized = False creates folded spectrum

# ### best model
# #nu1, nu2, b1, b2, hrf, m12, m21, Ts, Q, O
# emp_params = [1.5709, 3.2861, 0.2746,  23.8587, 0.2227, 0.2811, 1.2932,0.4282,0.979]
# model_fit = Fit_Empirical(fs, pts, prefix, "IM2NG", IM2NG, emp_params, fs_folded=True)
# dadi.Plotting.plot_2d_comp_multinom(data=fs,model=model_fit,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_ALFW_IM2NG.png"
# pylab.savefig(name)
# pylab.close()

# dadi.Plotting.plot_2d_comp_Poisson(data=fs,model=model_fit,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_ALFW_IM2NG_pois.png"
# pylab.savefig(name)
# pylab.close()

# # FLLG-ALST
# prefix="FLLG_ALFW"
# vmin_val=float(0.0001)
# fs=dadi.Spectrum.from_data_dict(dd , pop_ids = [ "FLLG", "ALFW" ],projections = [ projs[0], projs[1] ],polarized = False )  #polarized = False creates folded spectrum

# ### best model
# #nu1, nu2, b1, b2, hrf, m12, m21, Ts, Q, O
# emp_params = [1.5709, 3.2861, 0.2746,  23.8587, 0.2227, 0.2811, 1.2932,0.4282,0.979]
# model_fit = Fit_Empirical(fs, pts, prefix, "IM2NG", IM2NG, emp_params, fs_folded=True)
# dadi.Plotting.plot_2d_comp_multinom(data=fs,model=model_fit,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_ALFW_IM2NG.png"
# pylab.savefig(name)
# pylab.close()

# dadi.Plotting.plot_2d_comp_Poisson(data=fs,model=model_fit,vmin=vmin_val,show=False)
# name = "../../figs/dadi/FLLG_ALFW_IM2NG_pois.png"
# pylab.savefig(name)
# pylab.close()




