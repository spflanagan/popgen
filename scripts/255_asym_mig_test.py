'''
Running asym_mig to compare frequency spectra
'''

import os
import numpy
import dadi
import pylab
from datetime import datetime
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum



def asym_mig(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1, nu2, m12, m21, T = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    


# Set the parameters (from a run that resulted in lots of warnings)
nu1 = 1.01
nu2 = 15.4299
m12 = 0.9824
m21 = 1.1719
T = 3.6611
ns = [46,60]

# First set of points
pts = 200
xx = Numerics.default_grid(pts)
phi = PhiManip.phi_1D(xx)
phi = PhiManip.phi_1D_to_2D(xx, phi)
phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
fs200 = Spectrum.from_phi(phi, ns, (xx,xx))

# Second set of points
pts = 250
xx = Numerics.default_grid(pts)
phi = PhiManip.phi_1D(xx)
phi = PhiManip.phi_1D_to_2D(xx, phi)
phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
fs250 = Spectrum.from_phi(phi, ns, (xx,xx))

# Third set of points
pts = 300
xx = Numerics.default_grid(pts)
phi = PhiManip.phi_1D(xx)
phi = PhiManip.phi_1D_to_2D(xx, phi)
phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
fs300 = Spectrum.from_phi(phi, ns, (xx,xx))

