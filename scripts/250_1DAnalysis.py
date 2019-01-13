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
	Optimize_Routine(fs, pts, pops[i], "two_epoch", dadi.Demographics1D.two_epoch, rounds=5,param_number=2, fs_folded=True, param_labels = "nu, T")
	Optimize_Routine(fs, pts, pops[i], "growth", dadi.Demographics1D.growth, rounds=5,param_number=2, fs_folded=True, param_labels = "nu, T")
	p_labels = "nuB, nuF, T"
	Optimize_Routine(fs, pts, pops[i], "bottlegrowth", dadi.Demographics1D.bottlegrowth, rounds=5,param_number=3, fs_folded=True, param_labels = "nuB, nuF, T")
	p_labels = "nuB, nuF, TB, TF"
	Optimize_Routine(fs, pts, pops[i], "three_epoch", dadi.Demographics1D.three_epoch, rounds=5,param_number=4, fs_folded=True, param_labels = "nuB, nuF, TB, TF")
	os.chdir("../")



#This seems to work. Extracted best fits in R

#=================================================================================================#
#									COMPARE MODELS TO DATA										  #
#=================================================================================================#
import pylab

#pops = ['FLLG', 'FLCC', 'ALFW','ALST','LAFW','TXFW','TXCC']
#projs = [70,      61,     72,     70,    72,    46,     61]

#----------------------------------------------ALFW-----------------------------------------------#
pts = 80
nuB=0.0135
nuF=10.5357
T=0.2692
alfw_bottle = dadi.Demographics1D.bottlegrowth(params = (nuB,nuF,T), ns=[72], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(alfw_bottle,alfw)

#----------------------------------------------ALST-----------------------------------------------#
pts = 80
nuB=10.7382
nuF=28.716
T=0.9324
alst_bottle = dadi.Demographics1D.bottlegrowth(params = (nuB,nuF,T), ns=[70], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(alst_bottle,alst)

#----------------------------------------------FLCC-----------------------------------------------#
pts = 80
nu=29.7
T=0.0881
flcc_twoep = dadi.Demographics1D.two_epoch(params = (nu,T), ns=[61], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(flcc_twoep,flcc)

pts = 80
nu=29.0771
T=0.119
flcc_growth = dadi.Demographics1D.growth(params = (nu,T), ns=[61], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(flcc_growth,flcc)

#----------------------------------------------FLLG-----------------------------------------------#
pts = 80
nuB=0.0114
nuF=5.3033
T=0.0586
fllg_bottle = dadi.Demographics1D.bottlegrowth(params = (nuB,nuF,T), ns=[70], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(fllg_bottle,fllg)

#----------------------------------------------LAFW-----------------------------------------------#
pts = 80
nuB=0.0368
nuF=20.0107
T=0.5108
lafw_bottle = dadi.Demographics1D.bottlegrowth(params = (nuB,nuF,T), ns=[72], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(lafw_bottle,lafw)

#----------------------------------------------TXCC-----------------------------------------------#
pts = 80
nu=16.3039
T=1.1496
txcc_growth = dadi.Demographics1D.growth(params = (nu,T), ns=[61], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(txcc_growth,txcc)

#----------------------------------------------TXFW-----------------------------------------------#
pts = 80
nu=0.602
T=0.0198
txfw_growth = dadi.Demographics1D.growth(params = (nu,T), ns=[46], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(txfw_growth,txfw)

pts = 80
nu=0.5634
T=0.0319
txfw_twoep = dadi.Demographics1D.two_epoch(params = (nu,T), ns=[46], pts = pts)
dadi.Plotting.plot_1d_comp_multinom(txfw_twoep,txfw)






#=================================================================================================#
#								ONE AT A TIME THE EXPLORATORY WAY 								  #
#=================================================================================================#
# # set parameters 
# nu = 0.5 #ratio of contemporary to ancient population size
# T = 500 #2NA generations
# nuB = 0.1 #ratio of pop size after instantaneous change to ancient pop size
# nuF = 0.7 #ratio of pop size now to ancient pop size
# TB = 100 #length of bottleneck in 2Na generations
# TF = 100 #length of time since bottleneck recovery

# # Single population models

# snm =  dadi.Demographics1D.snm(fllg,ns = ns,pts= 60)

# two = dadi.Demographics1D.two_epoch(params = (nu,T), ns = ns, pts = 60)

# growth = dadi.Demographics1D.growth(params = (nu,T), ns = ns, pts = 60)

# bottle = dadi.Demographics1D.bottlegrowth(params = (nuB,nuF,T), ns = ns, pts = 60)

# three = dadi.Demographics1D.three_epoch(params = (nuB,nuF,TB,TF), ns = ns, pts = 60)

# # Compare models to data



# import pylab
# dadi.Plotting.plot_1d_comp_multinom(snm,fllg)
# dadi.Plotting.plot_1d_comp_Poisson(snm,fllg)

# dadi.Plotting.plot_1d_comp_multinom(two,fllg)

# dadi.Plotting.plot_1d_comp_multinom(growth,fllg)

# dadi.Plotting.plot_1d_comp_multinom(bottle,fllg)

# dadi.Plotting.plot_1d_comp_multinom(three,fllg)


# ### I think I need to optimize parameters before continuing.


# prefix = "FLLG1_1"
# pts = 60

# p_labels = "ns, pts"

# Optimize_Routine(fllg, 60, prefix, "snm", dadi.Demographics1D.snm, 3,2, fs_folded=True)

# #or with a loop
# for i in range(1,6):
#     prefix = "FLLG1_1_Number_{}".format(i)
#     Optimize_Routine(fllg, pts, prefix, "snm", dadi.Demographics1D.snm, 3, 2, fs_folded=True, param_labels = p_labels)

