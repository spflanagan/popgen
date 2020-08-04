#!/usr/bin/env python
# -*- coding: utf-8 -*- 

"""
module modeledemo with demographic models of divergence
downloaded from https://github.com/QuentinRougemont/DemographicInference/blob/master/02-modifs_v2/folded/modeledemo_new_models_folded.py
on 3 Aug 2020
methods from Rougemont et al 2020 https://www.biorxiv.org/content/10.1101/732750v4.full.pdf
"""

import numpy
import dadi


def SI(params, (n1,n2), pts):
    nu1, nu2, Ts = params
    """
    Model with split and complete isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    ### Sum the two spectra in proportion O
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def SI2N(params, (n1,n2), pts):
    nu1, nu2, Ts, nr, bf = params
    """
    Model with split and complete isolation, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: The scaled time of the split
    n1,n2: Size of fs to generate.
    nr: Proportion of non/low-recombining regions
    bf : Background factor (to which extent the effective population size is reduced in the non-recombining regions)
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # Spectrum of non-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2	
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    # Finally, calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    # Spectrum of recombining regions
    # phi for the equilibrium ancestral population
    phir = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    # We set the population sizes after the split to nu1 and nu2
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    # oriented
    fsr = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    ### Sum the two spectra in proportion O
    fs= (nr*fsnr + (1-nr)*fsr) 

    return fs


def SIG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, Ts = params
    
    """ 
    Model with split and complete isolation.

    nu1: Size of population 1 at split.
    nu2: Size of population 2 at split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split independantly in each population
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def SI2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, Ts, Q = params
    """
    Model with split and complete isolation, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    Ts: The scaled time of the split
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We start the population size change after the split independantly in each population
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) and integrate the hrf for low-recombining regions
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    ### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnr + Q*fslr) 
    return fs

def PAN(params, (n1), pts): 
    
    """
    Panmictic model

    nu1: Size of population 1
    n1 Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    nu, T = params
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    # phi for the equilibrium panmictic population
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T, nu)
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, n1, (xx,))
    return fs

def PANG(params, (n1), pts): 
    nu1, T, = params
    """
    Model with exponential growth some tme ago

    nu1: Ratio of contemporary to ancient population size
    T: Time in the past at which growth began (in units of 2*Na 
       generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # phi for the equilibrium ancestral population
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    #bnu1_func = lambda t: nu1 * b1**(t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, bnu1_func)
    ###
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1), (xx,))
    return fs

def PANG2N(params, (n1), pts): 
    nu1, hrf, T, Q = params
    """
    Model with exponential growth some tme ago

    ns = (n1,)
    nu1: Ratio of contemporary to ancient population size
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    T: Time in the past at which growth began (in units of 2*Na 
       generations) 
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """

    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    
    # phi for the equilibrium ancestral population
    # We start the population size change after the split and set the migration rates to m12 and m21
    #bnu1_func = lambda t: nu1 * b1**(t/T)
    bnu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/T)
    phinr = dadi.Integration.one_pop(phinr, xx,T, bnu1_func, )
    ###
    ## Finally, calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1), (xx,))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    
    # phi for the equilibrium ancestral population
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func_hrf = lambda t: (numpy.exp(numpy.log(nu1) * t/T) ) * hrf
    #bnu1_func_hrf = lambda t: (nu1 * b1**(t/T) )* hrf
    philr = dadi.Integration.one_pop(philr , xx,T, bnu1_func_hrf, )
    ###
    ## Finally, calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1), (xx,))

    ### Sum the two spectra in proportion Q
    fs= ((1-Q)*fsnr + Q*fslr) 
    return fs

def PANGb(params, ns, pts): 
    nuB,nuF,T = params
    """
    Instantanous size change followed by exponential growth. #=bottlegrowth

    params = (nu,T)
    ns = (n1,)
    nu1: Ratio of contemporary to ancient population size
    T: Time in the past at which growth began (in units of 2*Na 
       generations) 
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #bottleneck and growth
    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, nu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

def PANGb2N(params, ns, pts): 
    nuB, nuF, hrf, T, Q = params
    """
    Instantanous size change followed by exponential growth. #=bottlegrowth

    params = (nu,T)
    ns = (n1,)
    nu1: Ratio of contemporary to ancient population size
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    T: Time in the past at which growth began (in units of 2*Na 
       generations) 
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    # phi for the equilibrium ancestral population
    #NON recombining region:
    philr = dadi.PhiManip.phi_1D(xx)
    
    #bottleneck and growth
    nu_func_hrf = lambda t: (nuB*numpy.exp(numpy.log(nuF/nuB) * t/T) )* hrf
    philr = dadi.Integration.one_pop(philr, xx, T, nu_func_hrf)
    #spectrum
    fslr = dadi.Spectrum.from_phi(philr, ns, (xx,))
    
    # phi for the equilibrium ancestral population
    #recombining region:
    phinr = dadi.PhiManip.phi_1D(xx)
    
    #bottleneck and growth
    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
    phir = dadi.Integration.one_pop(phinr, xx, T, nu_func)
    #spectrum
    fsnr = dadi.Spectrum.from_phi(phinr, ns, (xx,))

    ### Sum the two spectra in proportion O (and Q)
    fs = ((1-Q)*fsnr + Q*fslr)  
    return fs

def IM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts = params
    """
    Model with migration during the divergence.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # Finally, calculate the spectrum.
    # oriented
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    ### Sum the two spectra in proportion O
    return fs

def IMG(params, (n1,n2), pts): 
    nu1, nu2, b1, b2, m12, m21, Ts = params
    """
    Model with migration during the divergence.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # mis-oriented
    return fs

def IM2N(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnr + Q*fslr) 
    return fs

def IM2N2m(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, me12, me21, Ts, P, Q = params

    """
    Model of semi permeability with split, ongoing migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts,nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the spectra
    fs = (P*fsN+(1-P)*fsI+(1-Q)*fsnr+Q*fslr) 
    return fs

def IM2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) & integrate the hrf for low-recombining regions and set the migration rates to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b1**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the two spectra in proportion O (and Q)
    fs= ((1-Q)*fsnr + Q*fslr) 
    return fs

def IM2mG(params, (n1,n2), pts): 
    nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, P = params
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Finally, calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## Finally, calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion O (and P)
    fs = (P*fsN + (1-P)*fsI) 
    return fs

def IM2N2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, P, Q = params

    """
    Model of semi permeability with split, ongoing migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    #### Spectrum of low-recombining regions

    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the spectra
    fs = (P*fsN+(1-P)*fsI+(1-Q)*fsnr+Q*fslr) 
    return fs


def AM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts = params

    """
    Model with split, ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21 
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


def AMG(params, (n1,n2), pts): 
    nu1, nu2, b1, b2, m12, m21, Tam, Ts = params
    """
    Model with split, ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
  	b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split independantly in each population and set the migration rate to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tam)
    bnu2_func = lambda t: nu2 * b2**(t/Tam)
    phi = dadi.Integration.two_pops(phi, xx, Tam, bnu1_func, bnu2_func, m12=m12, m21=m21)
    # We continue the population size change after ancient migration (until present) independantly in each population and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def AM2N(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, Tam, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration.
    Ts: The scaled time between the end of ancient migration and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to hrf*nu1 and hrf*nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the two spectra in proportion O and 1-O
    fs= ((1-Q)*fsnr + Q*fslr)
    return fs
    
    
def AM2N2m(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rate to zero 
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to hrf*nu1 and hrf*nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the spectra
    fs = (P*fsN+(1-P)*fsI+(1-Q)*fsnr+Q*fslr) 
    return fs


def AM2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, Tam, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration.
    Ts: The scaled time between the end of ancient migration and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We start the population size change after the split independantly in each population and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We start the population size change after the split independantly (bnu{1,2}_func) & integrate the hrf for low-recombining regions in each population and set the migration rates to zero
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b1**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the two spectra in proportion O and 1-O
    fs= ((1-Q)*fsnr + Q*fslr) 
    return fs
    

def AM2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, m12, m21, me12, me21, Tam, Ts, P = params
    """
    Model of semi permeability with split, ancient migration with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    P: The porportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes nu1 & nu2 after the split and set the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We start the population reduction after the split and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes nu1 & nu2 after the split and set the migration rate to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    # We start the population reduction after the split and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion O (and P)
    fs =(P*fsN + (1-P)*fsI) 
    return fs
    

def AM2N2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tam, Ts, P, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We SET the population size nu1 & nu2 after the split and set the migration rates to m12 & m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We set the population sizes change independently in each population after ancient migration to bnu1_func and bnu2_func and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population size after the split in each population and set the migration rates to me12 & me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes change after ancient migration to bnu1_func and bnu2_func and set the migration rates to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to me12 and me21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We set the population sizes after the split to hrf*bnu1_func and hrf*bnu2_func and set the migration rates to zero
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
 
    #### Sum the spectra
    fs = (P*fsN + (1-P)*fsI + (1-Q)*fsnr + Q*fslr) 
    return fs

def SC(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc = params
    """
    Model with split, complete isolation, followed by secondary contact

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def SC2N(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, Ts, Tsc, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = ((1-Q)*fsnr + Q*fslr)
    return fs

def SCG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, m12, m21, Ts, Tsc = params
   
    """
    Model with split, complete isolation, followed by secondary contact with exponential growth

    nu1: Size of population 1 at split.
    nu2: Size of population 2 at split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population reduction after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Calculate the spectrum
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs
    
def SC2N2m(params, (n1,n2), pts):
    nu1, nu2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = (Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI)
    return fs


def SC2NG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, Ts, Tsc, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We set the population sizes changes independently between pops after the split and isolation to nu1 and nu2 and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)   
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We set the population sizes changes independently between pops after the split and isolation to bnu1 and bnu2 (bnu{1,2}_func) & integrate the hrf and set the migration rates to m12 & m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf   
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = ((1-Q)*fsnr+Q*fslr) 
    return fs


def SC2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, m12, m21, me12, me21, Ts, Tsc, P = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates and exponential growth

    nu1: Size of pop 1 after split.
    nu2: Size of pop 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The porportion of the genome evolving neutrally.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population size change after the split and set the migration rates to me12 and me21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def SC2N2mG(params, (n1,n2), pts):
    nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Ts, Tsc, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population sizes change independently in each populations after the split to bnu1 and bnu2 and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes change after the split to bnu1 and bnu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes change independently in each populations after the isolation to bnu1 and bnu2 and set the migration rates to m12 & m21
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI 
    return fs


def IM2m(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, P = params
    
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def AM2m(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tam, P = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def SC2m(params, (n1,n2), pts):
    nu1, nu2, m12, m21, me12, me21, Ts, Tsc, P  = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def PIM2m(params, (n1,n2), pts):
    nu1, nu2, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts, Tam, Tsc, P = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    mA12: Migration from pop 2 to pop 1 (2*Na*m12).
    mA21: Migration from pop 1 to pop 2.
    meA12: Effective migration from pop 2 to pop 1 in genomic islands.
    meA21: Effective migration from pop 1 to pop 2 in genomic islands.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scaled time between the ancient migration and present.
    Tsc: The scaled time between the secondary contact and present (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
       ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
  
    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    #Set population sizes after the split to nu1 and nu2 and the migration rate to meA12 and meA21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=meA12, m21=meA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def PSC2m(params, (n1,n2), pts):
    nu1, nu2, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts, Tsc1, Tam, Tsc, P = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    mA12: Migration from pop 2 to pop 1 (2*Na*m12).
    mA21: Migration from pop 1 to pop 2.
    meA12: Effective migration from pop 2 to pop 1 in genomic islands.
    meA21: Effective migration from pop 1 to pop 2 in genomic islands.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tsc1: The scaled time between the first secondary contact and Tam (in units of 2*Na generations).
    Tam: The scaled time between the ancient migration (SC) and stop of gene flow.
    Tsc: The scaled time between the secondary contact and present (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc1, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    #Set population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to meA12 and meA21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=meA12, m21=meA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P 
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def PIM2N(params, (n1,n2), pts):
    nu1, nu2, hrf, mA12, mA21, m12, m21, Ts, Tam, Tsc, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    mA12: Migration from pop 2 to pop 1 (2*Na*m12).
    mA21: Migration from pop 1 to pop 2.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scaled time between the ancient migration and present.
    Tsc: The scaled time between the secondary contact and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    ### Calculate the neutral spectrum in normally recombining regions
    # phi for the equilibrium ancestral population
    phiNr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiNr = dadi.PhiManip.phi_1D_to_2D(xx, phiNr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phiNr = dadi.Integration.two_pops(phiNr, xx, Ts, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiNr = dadi.Integration.two_pops(phiNr, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiNr = dadi.Integration.two_pops(phiNr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phiNr, (n1,n2), (xx,xx))
  
    ### Calculate the genomic island spectrum low-recombining regions
    # phi for the equilibrium ancestral population
    phiLr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiLr = dadi.PhiManip.phi_1D_to_2D(xx, phiLr)
    #Set population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rate to meA12 and meA21
    phiLr = dadi.Integration.two_pops(phiLr, xx, Ts, nu1*hrf, nu2*hrf, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiLr = dadi.Integration.two_pops(phiLr, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiLr = dadi.Integration.two_pops(phiLr, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(phiLr, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (1-Q)*fsnr+Q*fslr 
    return fs

def PSC2N(params, (n1,n2), pts):
    nu1, nu2, hrf, mA12, mA21, m12, m21, Ts, Tsc1, Tam, Tsc, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    mA12: Migration from pop 2 to pop 1 (2*Na*m12).
    mA21: Migration from pop 1 to pop 2.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tsc1: The scaled time between the first secondary contact and Tam (in units of 2*Na generations).
    Tam: The scaled time between the ancient migration (SC) and stop of gene flow.
    Tsc: The scaled time between the secondary contact and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum in normally recombining region
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc1, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc1, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)

    ## calculate the spectrum.
    # oriented
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    
    ### Sum the spectra 
    fs = (1-Q)*fsnr + Q*fslr 
    return fs


########################### PSC2N2m PIM2N2m #############################################

def PIM2N2m(params, (n1,n2), pts):
    nu1, nu2, hrf, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts, Tam, Tsc, P, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    mA12: Migration from pop 2 to pop 1 (2*Na*m12).
    mA21: Migration from pop 1 to pop 2.
    meA12: Effective migration from pop 2 to pop 1 in genomic islands.
    meA21: Effective migration from pop 1 to pop 2 in genomic islands.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scaled time between the ancient migration and present.
    Tsc: The scaled time between the secondary contact and present (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
       ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))
  
    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    #Set population sizes after the split to nu1 and nu2 and the migration rate to meA12 and meA21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=meA12, m21=meA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Calculate the neutral spectrum in normally recombining regions
    # phi for the equilibrium ancestral population
    phiNr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiNr = dadi.PhiManip.phi_1D_to_2D(xx, phiNr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phiNr = dadi.Integration.two_pops(phiNr, xx, Ts, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiNr = dadi.Integration.two_pops(phiNr, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiNr = dadi.Integration.two_pops(phiNr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phiNr, (n1,n2), (xx,xx))
  
    ### Calculate the genomic island spectrum low-recombining regions
    # phi for the equilibrium ancestral population
    phiLr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiLr = dadi.PhiManip.phi_1D_to_2D(xx, phiLr)
    #Set population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rate to meA12 and meA21
    phiLr = dadi.Integration.two_pops(phiLr, xx, Ts, nu1*hrf, nu2*hrf, m12=meA12, m21=meA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiLr = dadi.Integration.two_pops(phiLr, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiLr = dadi.Integration.two_pops(phiLr, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(phiLr, (n1,n2), (xx,xx))


    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) + (Q*fslr +(1-Q)*fsnr) 
    return fs

def PSC2N2m(params, (n1,n2), pts):
    nu1, nu2, hrf, mA12, mA21, meA12, meA21, m12, m21, me12, me21, Ts, Tsc1, Tam, Tsc, P, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    mA12: Migration from pop 2 to pop 1 (2*Na*m12).
    mA21: Migration from pop 1 to pop 2.
    meA12: Effective migration from pop 2 to pop 1 in genomic islands.
    meA21: Effective migration from pop 1 to pop 2 in genomic islands.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tsc1: The scaled time between the first secondary contact and Tam (in units of 2*Na generations).
    Tam: The scaled time between the ancient migration (SC) and stop of gene flow.
    Tsc: The scaled time between the secondary contact and present (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc1, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    #Set population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to meA12 and meA21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=meA12, m21=meA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Calculate the neutral spectrum in normally recombining region
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and the migration rate to mA12 and mA21
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc1, nu1, nu2, m12=mA12, m21=mA21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc1, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ## We now start a new migration event at time Tsc:
    ## We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)

    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    
    ### Sum the two spectra in proportion P 
    fs = (P*fsN+(1-P)*fsI) + (Q*fslr +(1-Q)*fsnr)
    return fs


####### models with ancestral pop size change: ##########################
def AMA(params, (n1,n2), pts):
    nuA, nu1, nu2, m12, m21, Tp, Ts, Tam = params

    """
    Model with split, ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nuA)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21 
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=0, m21=0)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


def AMAG(params, (n1,n2), pts): 
    nuA, nu1, nu2, b1, b2, m12, m21, Tam, Tp, Ts = params
    """
    Model with split, ancient migration

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nuA)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split independantly in each population and set the migration rate to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tam)
    bnu2_func = lambda t: nu2 * b2**(t/Tam)
    phi = dadi.Integration.two_pops(phi, xx, Tam, bnu1_func, bnu2_func, m12=m12, m21=m21)
    # We continue the population size change after ancient migration (until present) independantly in each population and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

def AMA2N(params, (n1,n2), pts):
    nuA, nu1, nu2, hrf, m12, m21, Tam, Tp, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration.
    Ts: The scaled time between the end of ancient migration and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)

    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to hrf*nu1 and hrf*nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)

    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the two spectra in proportion O and 1-O
    fs= (1-Q)*fsnr + Q*fslr
    return fs
    
    
def AMA2N2m(params, (n1,n2), pts):
    nuA, nu1, nu2, hrf, m12, m21, me12, me21, Tam, Tp, Ts, P, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rate to zero 
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    # oriented
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We keep the population sizes after the split to hrf*nu1 and hrf*nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the spectra
    fs = P*fsN+(1-P)*fsI+(1-Q)*fsnr+Q*fslr
    return fs


def AMA2NG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, hrf, m12, m21, Tam, Tp, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Tam: The scaled time between the split and the end of ancient migration.
    Ts: The scaled time between the end of ancient migration and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We start the population size change after the split independantly in each population and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We start the population size change after the split independantly (bnu{1,2}_func) & integrate the hrf for low-recombining regions in each population and set the migration rates to zero
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b1**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)

    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the two spectra 
    fs= (1-Q)*fsnr + Q*fslr
    return fs


def AMA2mG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, m12, m21, me12, me21, Tam, Tp, Ts, P = params
    """
    Model of semi permeability with split, ancient migration with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Tam: The scaled time between the split and the end of ancient migration (in units of 2*Na generations).
    Ts: The scaled time between the end of ancient migration and present.
    P: The porportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes nu1 & nu2 after the split and set the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We start the population reduction after the split and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes nu1 & nu2 after the split and set the migration rate to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    # We start the population reduction after the split and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion O (and P)
    fs = P*fsN + (1-P)*fsI
    return fs


def AMA2N2mG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tam, Tp, Ts, P, Q = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We SET the population size nu1 & nu2 after the split and set the migration rates to m12 & m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We set the population sizes change independently in each population after ancient migration to bnu1_func and bnu2_func and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population size after the split in each population and set the migration rates to me12 & me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes change after ancient migration to bnu1_func and bnu2_func and set the migration rates to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))


    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to me12 and me21 
    philr = dadi.Integration.two_pops(philr, xx, Tam, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    # We set the population sizes after the split to hrf*bnu1_func and hrf*bnu2_func and set the migration rates to zero
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
 
    #### Sum the spectra
    fs = (P*fsN + (1-P)*fsI + (1-Q)*fsnr + Q*fslr)
    return fs






def SCA(params, (n1,n2), pts):
    nuA, nu1, nu2, m12, m21, Tp, Ts, Tsc = params
    """
    Model with split, complete isolation, followed by secondary contact

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nuA)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


def SCA2N(params, (n1,n2), pts):
    nuA, nu1, nu2, hrf, m12, m21, Tp, Ts, Tsc, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = ((1-Q)*fsnr + Q*fslr)
    return fs


def SCAG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, m12, m21, Tp, Ts, Tsc = params
   
    """
    Model with split, complete isolation, followed by secondary contact with exponential growth

    nu1: Size of population 1 at split.
    nu2: Size of population 2 at split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nuA)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population reduction after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## Calculate the spectrum
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs
    
    
def SCA2N2m(params, (n1,n2), pts):
    nuA, nu1, nu2, hrf, m12, m21, me12, me21, Tp, Ts, Tsc, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    philr = dadi.Integration.two_pops(philr, xx, Tsc, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
 
    ### Sum the spectra 
    fs = (Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI)
    return fs


def SCA2NG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, hrf, m12, m21, Tp, Ts, Tsc, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We set the population sizes changes independently between pops after the split and isolation to nu1 and nu2 and set the migration rates to zero
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)   
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    
    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We set the population sizes changes independently between pops after the split and isolation to bnu1 and bnu2 (bnu{1,2}_func) & integrate the hrf and set the migration rates to m12 & m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf   
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = ((1-Q)*fsnr+Q*fslr)
    return fs


def SCA2mG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, m12, m21, me12, me21, Tp, Ts, Tsc, P = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates and exponential growth

    nu1: Size of pop 1 after split.
    nu2: Size of pop 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The porportion of the genome evolving neutrally.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population size change after the split and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population size change after the split and set the migration rates to me12 and me21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P 
    fs = (P*fsN+(1-P)*fsI) 
    return fs


def SCA2N2mG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tp, Ts, Tsc, P, Q = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We start the population sizes change independently in each populations after the split to bnu1 and bnu2 and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Tsc)
    bnu2_func = lambda t: nu2 * b2**(t/Tsc)
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes change after the split to bnu1 and bnu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes change independently in each populations after the isolation to bnu1 and bnu2 and set the migration rates to m12 & m21
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split and isolation to nu1 and nu2 and set the migration rates to zero
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # We keep the population sizes after the split and isolation to nu1 and nu2 and set the migration rate to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Tsc)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Tsc)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Tsc, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs =(Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI)
    return fs


def IMA(params, (n1,n2), pts):
    nuA, nu1, nu2, m12, m21, Tp, Ts = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nuA
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ## calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    return fs


def IMAG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, m12, m21, Tp, Ts = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nuA: anc pop after size chage
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    return fs


def IMA2N(params, (n1,n2), pts):
    nuA, nu1, nu2, hrf, m12, m21, Tp, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    #### Sum the two spectra
    fs= (1-Q)*fsnr + Q*fslr
    return fs


def IMA2NG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, hrf, m12, m21, Tp, Ts, Q = params

    """
    Model with split, ancient migration, heterogenous effective population size (with 2 classes of loci shared by the two populations = Hill-Robertson effects)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)
    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)
    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We start the population size change after the split independantly in each population and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ## Finally, calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) & integrate the hrf for low-recombining regions and set the migration rates to m12 and m21
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b1**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    #### Sum the two spectra
    fs= (1-Q)*fsnr + Q*fslr 
    return fs


def IMA2mG(params, (n1,n2), pts): 
    nuA, nu1, nu2, b1, b2, m12, m21, me12, me21, Tp, Ts, P = params
    """
    Model with migration during the divergence with two type of migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    P: The porportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We start the population size change after the split, independently in each population, and set the migration rates to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ## Finally, calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We start the population size change after the split, independently in each population, and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## Finally, calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the spectra
    fs = P*fsN + (1-P)*fsI
    return fs

def IMA2N2m(params, (n1,n2), pts):
    nuA, nu1, nu2, hrf, m12, m21, me12, me21, Tp, Ts, P, Q = params

    """
    Model of semi permeability with split, ongoing migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts,nu1, nu2, m12=m12, m21=m21)
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)

    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    ## calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    philr = dadi.Integration.two_pops(philr, xx, Ts, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI
    return fs


def IMA2N2mG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, hrf, m12, m21, me12, me21, Tp, Ts, P, Q = params

    """
    Model of semi permeability with split, ongoing migration with 2 migration rates, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    #### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, bnu1_func, bnu2_func, m12=me12, m21=me21)
    ###
    ## calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    #### Calculate the spectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2 and the migration rates to m12 and m21 
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=m12, m21=m21)
    ###
    ## calculate the spectrum.
    #fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We set the population sizes after the split to hrf*nu1 and hrf*nu2 and the migration rates to m12 and m21 
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=m12, m21=m21)

    ## calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))

    ### Sum the spectra 
    fs = Q*fslr+(1-Q)*fsnr+P*fsN+(1-P)*fsI 
    return fs



def IMA2m(params, (n1,n2), pts):
    nuA, nu1, nu2, m12, m21, me12, me21, Tp, Ts, P = params
    
    """
    Model with migration during the divergence with two type of migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """    

    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    # oriented
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def AMA2m(params, (n1,n2), pts):
    nuA, nu1, nu2, m12, m21, me12, me21, Tp, Ts, Tam, P = params

    """
    Model of semi permeability with split, ancient migration with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    Tam: The scale time between the ancient migration and present.
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))


    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=me12, m21=me21)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tam, nu1, nu2, m12=0, m21=0)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI)
    return fs

def SCA2m(params, (n1,n2), pts):
    nuA, nu1, nu2, m12, m21, me12, me21, Tp, Ts, Tsc, P = params
    """
    Model of semi permeability with split, complete isolation, followed by secondary contact with 2 migration rates

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    me12: Effective migration from pop 2 to pop 1 in genomic islands.
    me21: Effective migration from pop 1 to pop 2 in genomic islands.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    Tsc: The scale time between the secondary contact and present.
    P: The proportion of the genome evolving neutrally
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    ### Calculate the neutral spectrum
    # phi for the equilibrium ancestral population
    phiN = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiN = dadi.Integration.one_pop(phiN, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiN = dadi.PhiManip.phi_1D_to_2D(xx, phiN)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiN = dadi.Integration.two_pops(phiN, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phiN = dadi.Integration.two_pops(phiN, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    # calculate the spectrum.
    fsN = dadi.Spectrum.from_phi(phiN, (n1,n2), (xx,xx))

    ### Calculate the genomic island spectrum
    # phi for the equilibrium ancestral population
    phiI = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phiI = dadi.Integration.one_pop(phiI, xx, Tp, nu=nuA)

    # Now do the divergence event
    phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to zero
    phiI = dadi.Integration.two_pops(phiI, xx, Ts, nu1, nu2, m12=0, m21=0)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to me12 and me21
    phiI = dadi.Integration.two_pops(phiI, xx, Tsc, nu1, nu2, m12=me12, m21=me21)
    # calculate the spectrum.
    fsI = dadi.Spectrum.from_phi(phiI, (n1,n2), (xx,xx))

    ### Sum the two spectra in proportion P (and O)
    fs = (P*fsN+(1-P)*fsI) 
    return fs

def SIA(params, (n1,n2), pts):
    nuA, nu1, nu2, Tp, Ts = params
    """
    Model with split and complete isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nuA)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    ### Sum the two spectra in proportion O
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs



def SIA2N(params, (n1,n2), pts):
    nuA, nu1, nu2, hrf,Tp, Ts, Q = params
    """
    Model with split and complete isolation, heterogenous effective population size (2 classes, shared by the two populations = background selection)
    nuA: Size of ancestral pop after change
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    Ts: The scaled time of the split
    n1,n2: Size of fs to generate.
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # Spectrum of non-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We set the population sizes after the split to nu1 and nu2
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*hrf, nu2*hrf, m12=0, m21=0)
    # Finally, calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    
    # Spectrum of recombining regions
    # phi for the equilibrium ancestral population
    phir = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phir = dadi.Integration.one_pop(phir, xx, Tp, nu=nuA)

    # Now do the divergence event
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    # We set the population sizes after the split to nu1 and nu2
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    fsr = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))

    ### Sum the two spectra
    fs= (1-Q)*fsnr + Q*fsr
    return fs


def SIAG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, Tp, Ts = params
    
    """ 
    Model with split and complete isolation.

    nu1: Size of population 1 at split.
    nu2: Size of population 2 at split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    Ts: The scaled time between the split and present (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nuA)

    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We start the population size change after the split independantly in each population
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ## Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


def SIA2NG(params, (n1,n2), pts):
    nuA, nu1, nu2, b1, b2, hrf, Tp, Ts, Q = params
    """
    Model with split and complete isolation, heterogenous effective population size (2 classes, shared by the two populations = background selection)

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    b1: Population growth coefficient of population 1
    b2: Population growth coefficient of population 2
    hrf: Hill-Robertson factor, i.e. the degree to which Ne is locally reduced due to the effects of background selection and selective sweep effects
    Ts: The scaled time of the split
    Q: The proportion of the genome with a reduced effective size due to selection at linked sites
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    #### Calculate the pectrum in normally-recombining regions
    # phi for the equilibrium ancestral population
    phinr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    phinr = dadi.Integration.one_pop(phinr, xx, Tp, nu=nuA)

    # Now do the divergence event
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # We start the population size change after the split independantly in each population
    bnu1_func = lambda t: nu1 * b1**(t/Ts)
    bnu2_func = lambda t: nu2 * b2**(t/Ts)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, bnu1_func, bnu2_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    fsnr = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))

    #### Spectrum of low-recombining regions
    # phi for the equilibrium ancestral population
    philr = dadi.PhiManip.phi_1D(xx)

    #Now do the population growth event.
    philr = dadi.Integration.one_pop(philr, xx, Tp, nu=nuA)

    # Now do the divergence event
    philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
    # We start the population size change after the split independantly in each population (bnu{1,2}_func) and integrate the hrf for low-recombining regions
    bnu1hrf_func = lambda t: (nu1 * b1**(t/Ts)) * hrf
    bnu2hrf_func = lambda t: (nu2 * b2**(t/Ts)) * hrf
    philr = dadi.Integration.two_pops(philr, xx, Ts, bnu1hrf_func, bnu2hrf_func, m12=0, m21=0)
    ###
    ## Finally, calculate the spectrum.
    fslr = dadi.Spectrum.from_phi(philr, (n1,n2), (xx,xx))
    
    ### Sum the two spectra 
    fs= (1-Q)*fsnr + Q*fslr
    return fs




