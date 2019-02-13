"""
    Custom demographic models
"""

import numpy
import dadi
import sys



#=================================================================================================#
#                                           GROWTH                                                #
#=================================================================================================#

#####################################################
#MODELS WITH GROWTH IN BOTH POPS
def  growth_no_mig(params,ns,pts):
    """
        Split with symmetric migration, followed by a bottleneck in one population with ongoing symmmetric migration
        nu1  : Size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2f : Final size of the population 2
        T    : Total time since the split
    """

    # get the params 
    nu1,nu1f,nu2,nu2f,T = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #modified from founder_*() eqs in Models_2D.py
    #exponential growth
    nu1_func = lambda t: nu1 * (nu1f/nu1)**(t/T)
    nu2_func = lambda t: nu2 * (nu2f/nu2)**(t/T)
    #make it happen
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=0, m21=0)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs


def  growth_sym_mig(params,ns,pts):
    """
        Split with symmetric migration, followed by a bottleneck in one population with ongoing symmmetric migration
        nu1  : Size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2f : Final size of the population 2
        m    : Migration rate
        T    : Total time since the split
    """

    # get the params 
    nu1,nu1f,nu2,nu2f,m,T = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #modified from founder_*() eqs in Models_2D.py
    #exponential growth
    nu1_func = lambda t: nu1 * (nu1f/nu1)**(t/T)
    nu2_func = lambda t: nu2 * (nu2f/nu2)**(t/T)
    #make it happen
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m, m21=m)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def  growth_sym_mig(params,ns,pts):
    """
        Split with symmetric migration, followed by a bottleneck in one population with ongoing symmmetric migration
        nu1  : Size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2f : Final size of the population 2
        m12  : Migration rate from pop 1 to pop 2
        m21  : Migration rate from pop 2 to pop 1
        T    : Total time since the split
    """

    # get the params 
    nu1,nu1f,nu2,nu2f,m12,m21,T = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #modified from founder_*() eqs in Models_2D.py
    #exponential growth
    nu1_func = lambda t: nu1 * (nu1f/nu1)**(t/T)
    nu2_func = lambda t: nu2 * (nu2f/nu2)**(t/T)
    #make it happen
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m12, m21=m21)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

#####################################################
#MODELS WITH GROWTH IN ONE POP, TWO-EPOCH IN OTHER
def  growth_twoep_no_mig(params,ns,pts):
    """
        Split with no migration
        Exponential growth in one population and instantaneous size change in the other
        
        nu1: Ratio of contemporary to population 1 size after split
        nu2: Ratio of contemporary to population 2 size after split
        T    : Total time since the split
        Tc   : Time at which the size change happened as a proportion of time since split
    """

    # get the params 
    nu1,nu2,T,Tc = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    #first epoch
    phi = Integration.two_pops(phi, xx, T, 1, 1, m12=0, m21=0)
    #pop 1 has exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    #second epoch
    phi = Integration.two_pops(phi, xx, T*Tc, nu1_func, nu2, m12=0, m21=0)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs



def  growth_twoep_sym_mig(params,ns,pts):
    """
        Split with symmetric migration
        Exponential growth in one population and instantaneous size change in the other
        
        nu1: Ratio of contemporary to population 1 size after split
        nu2: Ratio of contemporary to population 2 size after split
        m  : Migration rate
        T    : Total time since the split
        Tc   : Time at which the size change happened as a proportion of time since split
    """

   # get the params 
    nu1,nu2,m,T,Tc = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    #first epoch
    phi = Integration.two_pops(phi, xx, T, 1, 1, m12=m, m21=m)
    #pop 1 has exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    #second epoch
    phi = Integration.two_pops(phi, xx, T*Tc, nu1_func, nu2, m12=m, m21=m)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs



def  growth_twoep_asym_mig(params,ns,pts):
    """
        Split with asymmetric migration
        Exponential growth in one population and instantaneous size change in the other
        
        nu1: Ratio of contemporary to population 1 size after split
        nu2: Ratio of contemporary to population 2 size after split
        m12: Migration rate from pop 1 to pop 2
        m21: Migration rate from pop 2 to pop 1
        T    : Total time since the split
        Tc   : Time at which the size change happened as a proportion of time since split
    """

   # get the params 
    nu1,nu2,m12,m21,T,Tc = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    #first epoch
    phi = Integration.two_pops(phi, xx, T, 1, 1, m12=m12, m21=m21)
    #pop 1 has exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    #second epoch
    phi = Integration.two_pops(phi, xx, T*Tc, nu1_func, nu2, m12=m12, m21=m21)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs



#=================================================================================================#
#                                   GROWTH & BOTTLENECKS                                          #
#=================================================================================================#
#MODELS WITH GROWTH IN ONE POP, BOTTLENECK IN OTHER
def  growth_bottle_no_mig(params,ns,pts):
    """
        Split with no migration
        Exponential growth in population 1
        instantaneous size change followed by exponential growth in pop 2
        
        nu1: Ratio of contemporary to population 1 size after split
        nu2b: Ratio of bottleneck population size to split population size for pop 2
        nu2f: Ratio of contemporary to split population size for pop 2
        T    : Total time since the split
        Tc   : Time at which the size change happened as a proportion of time since split
    """ 

    # get the params 
    s,nu1,nu2,T,Tc = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    #first epoch
    phi = Integration.two_pops(phi, xx, T, 1, 1, m12=0, m21=0)
    #pop 1 has exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    #pop 2 has bottleneck
    nu2_func = lambda t: nu2b*numpy.exp(numpy.log(nu2f/nu2b) * t/T)
    #second epoch
    phi = Integration.two_pops(phi, xx, T*Tc, nu1_func, nu2_func, m12=0, m21=0)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs


def  growth_bottle_sym_mig(params,ns,pts):
    """
        Split with symmetric migration
        Exponential growth in population 1
        instantaneous size change followed by exponential growth in pop 2
        
        nu1: Ratio of contemporary to population 1 size after split
        nu2b: Ratio of bottleneck population size to split population size for pop 2
        nu2f: Ratio of contemporary to split population size for pop 2
        m: migration rate
        T    : Total time since the split
        Tc   : Time at which the size change happened as a proportion of time since split
    """

    # get the params 
    s,nu1,nu2,m,T,Tc = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    #first epoch
    phi = Integration.two_pops(phi, xx, T, 1, 1, m12=m, m21=m)
    #pop 1 has exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    #pop 2 has bottleneck
    nu2_func = lambda t: nu2b*numpy.exp(numpy.log(nu2f/nu2b) * t/T)
    #second epoch
    phi = Integration.two_pops(phi, xx, T*Tc, nu1_func, nu2_func, m12=m, m21=m)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs


def  growth_bottle_asym_mig(params,ns,pts):
    """
        Split with asymmetric migration
        Exponential growth in population 1
        instantaneous size change followed by exponential growth in pop 2
        
        nu1: Ratio of contemporary to population 1 size after split
        nu2b: Ratio of bottleneck population size to split population size for pop 2
        nu2f: Ratio of contemporary to split population size for pop 2
        m12: migration rate from pop 1 to pop 2
        m21: migration rate from pop 2 to pop 1
        T    : Total time since the split
        Tc   : Time at which the size change happened as a proportion of time since split
    """

    # get the params 
    s,nu1,nu2,m12,m21,T,Tc = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    #first epoch
    phi = Integration.two_pops(phi, xx, T, 1, 1, m12=m12, m21=m21)
    #pop 1 has exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    #pop 2 has bottleneck
    nu2_func = lambda t: nu2b*numpy.exp(numpy.log(nu2f/nu2b) * t/T)
    #second epoch
    phi = Integration.two_pops(phi, xx, T*Tc, nu1_func, nu2_func, m12=m12, m21=m21)

    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs



#=================================================================================================#
#                                        BOTTLENECKS                                              #
#=================================================================================================#
"""
    These 2D models are modified from models found on the dadi google group 
    (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    They are modified to work with dportik's dadi_pipeline
"""

#####################################################
# MODELS WITH BOTTLENECKS IN ONE POP
def  bottle1_sym_mig(params,ns,pts):
    """
        Split with symmetric migration, followed by a bottleneck in one population with ongoing symmmetric migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        m    : Migration rate
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,m,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m, m21=m)
    # bottleneck in pop 1
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2, m12=m, m21=m)
    # recovery from bottleneck in pop 1
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2, m12=m, m21=m)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def  bottle1_asym_mig(params,ns,pts):
    """
        Split with asymmetric migration, followed by a bottleneck in one population with ongoing asymmmetric migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        m12  : Migration rate from pop 1 to pop 2
        m21  : Migration rate from pop 2 to pop 1
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,m12,m21,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    # bottleneck in pop 1
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2, m12=m12, m21=m21)
    # recovery from bottleneck in pop 1
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2, m12=m12, m21=m21)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def  bottle1_no_mig(params,ns,pts):
    """
        Split with no migration, followed by a bottleneck in one population with no migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=0, m21=0)
    # bottleneck in pop 1
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2, m12=0, m21=0)
    # recovery from bottleneck in pop 1
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2, m12=0, m21=0)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

#####################################################
# MODELS WITH BOTTLENECKS IN BOTH POPS

def  bottle2_sym_mig(params,ns,pts):
    """
        Split with symmetric migration, followed by a bottleneck in both populations with ongoing symmmetric migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2b : Bottleneck size of pop 2
        nu2f : Final size of pop 2
        m    : Migration rate
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,nu2b,nu2f,m,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m, m21=m)
    # bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2b, m12=m, m21=m)
    # recovery from bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2f, m12=m, m21=m)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def  bottle2_asym_mig(params,ns,pts):
    """
        Split with asymmetric migration, followed by a bottleneck in both populations with ongoing asymmmetric migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2b : Bottleneck size of pop 2
        nu2f : Final size of pop 2
        m12  : Migration rate from pop 1 to pop 2
        m21  : Migration rate from pop 2 to pop 1
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,nu2b,nu2f,m12,m21,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    # bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2b, m12=m12, m21=m21)
    # recovery from bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2f, m12=m12, m21=m21)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def  bottle2_no_mig(params,ns,pts):
    """
        Split with no migration, followed by a bottleneck in one population with no migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2b : Bottleneck size of pop 2
        nu2f : Final size of pop 2
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,nu2b,nu2f,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=0, m21=0)
    # bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2b, m12=0, m21=0)
    # recovery from bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2f, m12=0, m21=0)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs


#####################################################
# MODELS WITH HISTORICAL MIGRATION THEN BOTTLENECK IN ONE POP

def  bottle1_hist_sym_mig(params,ns,pts):
    """
        Split with symmetric migration, followed by a bottleneck in one population with no migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        m    : Migration rate
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,m,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m, m21=m)
    # bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2, m12=0, m21=0)
    # recovery from bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2, m12=0, m21=0)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def  bottle1_hist_asym_mig(params,ns,pts):
    """
        Split with asymmetric migration, followed by a bottleneck in one population with no migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        m12  : Migration rate from pop 1 to pop 2
        m21  : Migration rate from pop 2 to pop 1
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,m12,m21,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    # bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2, m12=0, m21=0)
    # recovery from bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2, m12=0, m21=0)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

#####################################################
# MODELS WITH HISTORICAL MIGRATION THEN BOTTLENECK IN ONE POP

def  bottle2_hist_sym_mig(params,ns,pts):
    """
        Split with symmetric migration, followed by a bottleneck in both populations with no migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2b : Bottleneck size of pop 2
        nu2f : Final size of pop 2
        m    : Migration rate
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover

        modified from models found on the dadi google group 
        (https://groups.google.com/forum/#!searchin/dadi-user/2d$20bottleneck%7Csort:date/dadi-user/gUYz8QzT18c/F8cmzI0ZCQAJ)
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,nu2b,nu2f,m,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m, m21=m)
    # bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2b, m12=0, m21=0)
    # recovery from bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2f, m12=0, m21=0)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs

def  bottle2_hist_asym_mig(params,ns,pts):
    """
        Split with asymmetric migration, followed by a bottleneck in both populations with no migration
        nu1  : Size of the population 1
        nu1b : Bottleneck size of the population 1
        nu1f : Final size of the population 1
        nu2  : Size of the population 2
        nu2b : Bottleneck size of pop 2
        nu2f : Final size of pop 2
        m12  : Migration rate from pop 1 to pop 2
        m21  : Migration rate from pop 2 to pop 1
        T    : Total time since the split
        TB   : Time since bottleneck
        TF   : Time since recover
    """

    # get the params 
    nu1,nu1b,nu1f,nu2,nu2b,nu2f,m12,m21,T,TB,TF = params

    # grid
    xx = dadi.Numerics.default_grid(pts)
    # phi of the ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # divergence, with gene flow
    phi = dadi.Integration.two_pops(phiX, xx, T, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    # bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TB, nu1=nu1b, nu2=nu2b, m12=0, m21=0)
    # recovery from bottleneck in both pops
    phi = dadi.Integration.two_pops(phiX, xx, TF, nu1=nu1f, nu2=nu2f, m12=0, m21=0)
    
    # computing fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs
