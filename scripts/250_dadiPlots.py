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


# get the data
dd = dadi.Misc.make_data_dict ( "fwsw75.dadi.snps" )

# create a spectrum
pops = ['FLLG', 'FLCC', 'ALFW', 'ALST', 'LAFW', 'TXFW','TXCC' ]
projs = [  70, 60 , 72, 70, 72, 46, 61 ]

# loop through for 1D 
for i in range(0,len(pops)):
    spect = dadi.Spectrum.from_data_dict(dd , pop_ids = [ pops[i] ],projections = [ projs[i] ],polarized = False )  #polarized = False creates folded spectrum
    name = "../../figs/dadi/{}_1D.png".format(pops[i])
    dadi.Plotting.plot_1d_fs(spect,show=False)
    pylab.savefig(name)


# loop through for 2D
for i in range(0,(len(pops)-1)):
    for j in range(i+1,len(pops)):
        spect = dadi.Spectrum.from_data_dict(dd , pop_ids = [ pops[i], pops[j] ],projections = [ projs[i], projs[j] ],polarized = False )  #polarized = False creates folded spectrum
        name = "../../figs/dadi/%s-%s_2D.png" % (pops[i],pops[j])
        dadi.Plotting.plot_single_2d_sfs(spect,vmin=0.0001)
        pylab.savefig(name)
        pylab.close()

