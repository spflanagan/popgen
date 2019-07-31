"""
Testing ms models compared to dadi models
"""

import os
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


# Full model - ms does not produce a sfs that looks right
nu1 = 0.1317
nu2 = 8.4225
m12 = 0.7777
m21 = 0.0561
T = 0.1441
ns = [46,60]
pts = 300

asym_mig_fs = asym_mig(params=(nu1,nu2,m12,m21,T),ns=ns,pts=pts)
asym_mig_folded = dadi.Spectrum.fold(asym_mig_fs)
dadi.Plotting.plot_single_2d_sfs(asym_mig_folded,vmin=0.000001)

core="-m 1 2 1.5554 -m 2 1 0.1122 -n 1 0.1317 -n 2 8.4225 -ej 0.07205 2 1"
command=dadi.Misc.ms_command(theta=0.1026112, ns=(60,70), core=core, iter=12103)
ms_fs=dadi.Spectrum.from_ms_file(os.popen(command))
folded_sfs=dadi.Spectrum.fold(ms_fs)
dadi.Plotting.plot_single_2d_sfs(folded_sfs,vmin=0.000001)



#### What if I try it without some of these components?
# equal migration
asym_mig_fs = asym_mig(params=(nu1,nu2,m12,m12,T),ns=ns,pts=pts)
asym_mig_folded = dadi.Spectrum.fold(asym_mig_fs)
dadi.Plotting.plot_single_2d_sfs(asym_mig_folded,vmin=0.000001)

core="-m 1 2 1.5554 -m 2 1 1.5554 -n 1 0.1317 -n 2 8.4225 -ej 0.07205 2 1"
command=dadi.Misc.ms_command(theta=0.1026112, ns=(60,70), core=core, iter=12103)
ms_fs=dadi.Spectrum.from_ms_file(os.popen(command))
folded_sfs=dadi.Spectrum.fold(ms_fs)
dadi.Plotting.plot_single_2d_sfs(folded_sfs,vmin=0.000001)

# no migration
asym_mig_fs = asym_mig(params=(nu1,nu2,0,0,T),ns=ns,pts=pts)
asym_mig_folded = dadi.Spectrum.fold(asym_mig_fs)
dadi.Plotting.plot_single_2d_sfs(asym_mig_folded,vmin=0.000001)

core="-m 1 2 0 -m 2 1 0 -n 1 0.1317 -n 2 8.4225 -ej 0.07205 2 1"
command=dadi.Misc.ms_command(theta=0.1026112, ns=(60,70), core=core, iter=12103)
ms_fs=dadi.Spectrum.from_ms_file(os.popen(command))
folded_sfs=dadi.Spectrum.fold(ms_fs)
dadi.Plotting.plot_single_2d_sfs(folded_sfs,vmin=0.000001)

# Without the ej term
asym_mig_fs = asym_mig(params=(nu1,nu2,m12,m12,0),ns=ns,pts=pts)
asym_mig_folded = dadi.Spectrum.fold(asym_mig_fs)
dadi.Plotting.plot_single_2d_sfs(asym_mig_folded,vmin=0.000001)

core="-m 1 2 1.5554 -m 2 1 0.1122 -n 1 0.1317 -n 2 8.4225"
command=dadi.Misc.ms_command(theta=0.1026112, ns=(60,70), core=core, iter=12103)
ms_fs=dadi.Spectrum.from_ms_file(os.popen(command))
folded_sfs=dadi.Spectrum.fold(ms_fs)
dadi.Plotting.plot_single_2d_sfs(folded_sfs,vmin=0.000001)

