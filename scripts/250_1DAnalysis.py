'''
Testing dadi with 1D demographic scenarios
Run this from outside the dadi directory
'''

#start with ipython -pylab from ~/Research/popgen/fwsw_results/dadi_analysis

# Numpy is the numerical library dadi is built upon
from numpy import array
import os
import dadi

#use dportik's functions
#get the optimize functions
%load ../../programs/dadi_pipeline-master/Optimize_Functions.py



# Load the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )
#projections is sample size of alleles
#need to use MINIMUM projections
pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
projs = [70,61,72,70,72,46,61]


#=================================================================================================#
#									CHECK THE FREQUENCY SPECTRA 								  #
#=================================================================================================#

#save each figure to a file
fllg = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLLG' ],projections =[70] ,polarized = False )  
dadi.Plotting.plot_1d_fs(fllg)

flcc = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'FLCC' ],projections =[61] ,polarized = False )  
dadi.Plotting.plot_1d_fs(flcc)

alfw = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'ALFW' ],projections =[72] ,polarized = False )  
dadi.Plotting.plot_1d_fs(alfw)

alst = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'ALST' ],projections =[70] ,polarized = False )  
dadi.Plotting.plot_1d_fs(alst)

lafw = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'LAFW' ],projections =[72] ,polarized = False )  
dadi.Plotting.plot_1d_fs(lafw)

txfw = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'TXFW' ],projections =[46] ,polarized = False )  
dadi.Plotting.plot_1d_fs(txfw)

txcc = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'TXCC' ],projections =[61] ,polarized = False )  
dadi.Plotting.plot_1d_fs(txcc)


#=================================================================================================#
#										LOOP TO OPTIMIZE 										  #
#=================================================================================================#




for i in range(len(pops)):
	print(pops[i])
	if(os.path.exists(pops[i])):
		os.chdir(pops[i])
	else:	
		os.mkdir( pops[i])
		os.chdir(pops[i])
	fs = dadi.Spectrum.from_data_dict(dd , pop_ids = [ pops[i] ],projections = [ projs[i] ],polarized = False ) 
	pts = projs[i]+10

	p_labels = "nu, T"
	Optimize_Routine(fs, pts, pops[i], "two_epoch", dadi.Demographics1D.two_epoch, 5,2, fs_folded=True)
	Optimize_Routine(fs, pts, pops[i], "growth", dadi.Demographics1D.growth, 5,2, fs_folded=True)
	p_labels = "nuB, nuF, T"
	Optimize_Routine(fs, pts, pops[i], "bottlegrowth", dadi.Demographics1D.bottlegrowth, 5,3, fs_folded=True)
	p_labels = "nuB, nuF, TB, TF"
	Optimize_Routine(fs, pts, pops[i], "bottlegrowth", dadi.Demographics1D.three_epoch, 5,4, fs_folded=True)
	os.chdir("../")



#This seems to work.



#=================================================================================================#
#								ONE AT A TIME THE EXPLORATORY WAY 								  #
#=================================================================================================#
# set parameters 
nu = 0.5 #ratio of contemporary to ancient population size
T = 500 #2NA generations
nuB = 0.1 #ratio of pop size after instantaneous change to ancient pop size
nuF = 0.7 #ratio of pop size now to ancient pop size
TB = 100 #length of bottleneck in 2Na generations
TF = 100 #length of time since bottleneck recovery

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


prefix = "FLLG1_1"
pts = 60

p_labels = "ns, pts"

Optimize_Routine(fllg, 60, prefix, "snm", dadi.Demographics1D.snm, 3,2, fs_folded=True)

#or with a loop
for i in range(1,6):
    prefix = "FLLG1_1_Number_{}".format(i)
    Optimize_Routine(fllg, pts, prefix, "snm", dadi.Demographics1D.snm, 3, 2, fs_folded=True, param_labels = p_labels)

