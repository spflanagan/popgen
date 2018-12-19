'''
Testing dadi with 1D demographic scenarios
Run this from outside the dadi directory
'''

#start with ipython -pylab from ~/Research/popgen/fwsw_results/dadi_analysis

# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi



# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw.dadi.pruned.snps" )
#full: pop_ids =[ 'TXSP','TXCC','TXFW','TXCB','LAFW','ALST','ALFW','FLSG','FLKB','FLFD','FLSI','FLAB','FLPB','FLHB','FLCC','FLLG' ]
#   projections=[][124,  82,  62,  72,  96,  94,  96,  88,  84,  80,  90,  84,  86,  82,  82,  94] #projections is sample size of alleles
#need to use MINIMUM projections
#let's start with FLFW
fs = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLLG' ],projections =[70] ,polarized = False )  #causes FutureWarning
ns = fs.sample_sizes

# set parameters 
nu = 0.5 #ratio of contemporary to ancient population size
T = 500 #2NA generations
nuB = 0.1 #ratio of pop size after instantaneous change to ancient pop size
nuF = 0.7 #ratio of pop size now to ancient pop size
TB = 100 #length of bottleneck in 2Na generations
TF = 100 #length of time since bottleneck recovery



fllg = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLLG' ],projections =[70] ,polarized = False )  
dadi.Plotting.plot_1d_fs(fllg)

flcc = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLCC' ],projections =[60] ,polarized = False )  
dadi.Plotting.plot_1d_fs(flcc)

alfw = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'ALFW' ],projections =[72] ,polarized = False )  
dadi.Plotting.plot_1d_fs(alfw)

alst = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'ALST' ],projections =[70] ,polarized = False )  
dadi.Plotting.plot_1d_fs(alst)

lafw = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'LAFW' ],projections =[72] ,polarized = False )  
dadi.Plotting.plot_1d_fs(lafw)

txfw = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'TXFW' ],projections =[46] ,polarized = False )  
dadi.Plotting.plot_1d_fs(txfw)

txcb = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'TXCB' ],projections =[54] ,polarized = False )  
dadi.Plotting.plot_1d_fs(txcb)


# Single population models

snm =  dadi.Demographics1D.snm(fllg,ns = ns,pts= 60)

two = dadi.Demographics1D.two_epoch(params = (nu,T), ns = ns, pts = 60)

growth = dadi.Demographics1D.growth(params = (nu,T), ns = ns, pts = 60)

bottle = dadi.Demographics1D.bottlegrowth(params = (nuB,nuF,T), ns = ns, pts = 60)

three = dadi.Demographics1D.three_epoch(params = (nuB,nuF,TB,TF), ns = ns, pts = 60)

# Compare models to data



import pylab
dadi.Plotting.plot_1d_comp_multinom(snm,fllg)
dadi.Plotting.plot_1d_comp_Poisson(snm,fllg)

dadi.Plotting.plot_1d_comp_multinom(two,fllg)

dadi.Plotting.plot_1d_comp_multinom(growth,fllg)

dadi.Plotting.plot_1d_comp_multinom(bottle,fllg)

dadi.Plotting.plot_1d_comp_multinom(three,fllg)


### I think I need to optimize parameters before continuing.