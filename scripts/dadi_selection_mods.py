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


	
def asym_mig_O(params, ns, pts):
    """
    Split into two populations, with different migration rates.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	O: The proportion of accurate orientation
	"""
    nu1, nu2, m12, m21, T, O = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
	## Calculate the spectrum for correctly oriented loci
	fsO = Spectrum.from_phi(phi, ns, (xx,xx))
	## Calculate the spectrum for incorrectly oriented loci
	fsM = dadi.Numerics.reverse_array(fsO)
	
	# Sum the two spectra in proportion O
    fs = O*fsO+(1-O)*fsM
    return fs

def asym_mig_2N(params, ns, pts):
    """
    Split into two populations, with different migration rates.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	O: The proportion of accurate orientation
	Q: The proportion with reduced effective size due to selection at linked sites
	hrf: Hill-Robertson factor, the degree to which Ne is reduced due to background selection
	"""
    nu1, nu2, m12, m21, T, O, Q, hrf = params
    xx = Numerics.default_grid(pts)
   
	# Spectrum for normally recombining regions
	## phi for equilibrium ancestral population
    phinr = PhiManip.phi_1D(xx)
	## divergence event
    phinr = PhiManip.phi_1D_to_2D(xx, phinr)
    ## set population sizes after split with different migration rates
    phinr = Integration.two_pops(phinr, xx, T, nu1, nu2, m12=m12, m21=m21)
	## Calculate the spectrum for correctly oriented loci
    fsnrO = Spectrum.from_phi(phinr, ns, (xx,xx))
    ## Calculate the spectrum for incorrectly oriented loci
	fsnrM = dadi.Numerics.reverse_array(fsnrO)
	
	# Spectrum for low-recombining regions
	## phi for equilibrium ancestral population
	philr = dadi.PhiManip.phi_1D(xx)
	## divergence event
	philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
	## set population sizes after split with different migration rates
	philr = Integration.two_pops(philr, xx, T, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
	## Calculate the spectrum for correctly oriented loci
	fslrO = Spectrum.from_phi(philr, ns, (xx,xx))
	## Calculate the spectrum for incorrectly oriented loci
	fslrM = dadi.Numerics.reverse_array(fslrO)
	
	# Sum the two spectra in proportion O
    fs= O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
	
    return fs
	

def asym_mig_2m(params, ns, pts):
    """
    Split into two populations, with different migration rates.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	me12: Effective migration from pop 2 to pop 1 in genomic islands
	me21: Effective migration from pop 1 to pop 2 in genomic islands
	O: The proportion of accurate orientation
	P: The proportion of the genome evolving neutrally
	"""
    nu1, nu2, T, m12, m21, me12, me21, O, P = params
    xx = Numerics.default_grid(pts)
   
	# Spectrum for neutral spectrum
	## phi for equilibrium ancestral population
    phiN = PhiManip.phi_1D(xx)
	## divergence event
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)
    ## set population sizes after split with different migration rates
    phiN = Integration.two_pops(phiN, xx, T, nu1, nu2, m12=m12, m21=m21)
	## Calculate the spectrum for correctly oriented loci
    fsNO = Spectrum.from_phi(phiN, ns, (xx,xx))
    ## Calculate the spectrum for incorrectly oriented loci
	fsNM = dadi.Numerics.reverse_array(fsNO)
	
	# Spectrum for genomic islands
	## phi for equilibrium ancestral population
	phiI = dadi.PhiManip.phi_1D(xx)
	## divergence event
	phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
	## set population sizes after split with different migration rates
	phiI = Integration.two_pops(phiI, xx, T, nu1, nu2, m12=me12, m21=me21)
	## Calculate the spectrum for correctly oriented loci
	fsIO = Spectrum.from_phi(phiI, ns, (xx,xx))
	## Calculate the spectrum for incorrectly oriented loci
	fsIM = dadi.Numerics.reverse_array(fsIO)
	
	# Sum the two spectra in proportion O
    fs = O*(P*fsNO + (1-P)*fsIO) + (1-O)*(P*fsNM + (1-P)*fsIM)
	
    return fs
	
	
def asym_mig_2N2m(params, ns, pts):
    """
    Split into two populations, with different migration rates.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	me12: Effective migration from pop 2 to pop 1 in genomic islands
	me21: Effective migration from pop 1 to pop 2 in genomic islands
	O: The proportion of accurate orientation
	P: The proportion of the genome evolving neutrally
	Q: The proportion of the genome with reduced effective size due to linked sites
	hrf: Hill-Robertson factor
	"""
    nu1, nu2, T, m12, m21, me12, me21, O, P, Q, hrf = params
    xx = Numerics.default_grid(pts)
   
	# Spectrum for neutral loci
	## phi for equilibrium ancestral population
    phiN = PhiManip.phi_1D(xx)
	## divergence event
    phiN = PhiManip.phi_1D_to_2D(xx, phiN)
    ## set population sizes after split with different migration rates
    phiN = Integration.two_pops(phiN, xx, T, nu1, nu2, m12=m12, m21=m21)
	## Calculate the spectrum for correctly oriented loci
    fsNO = Spectrum.from_phi(phiN, ns, (xx,xx))
    ## Calculate the spectrum for incorrectly oriented loci
	fsNM = dadi.Numerics.reverse_array(fsNO)
	
	# Spectrum for genomic islands
	## phi for equilibrium ancestral population
	phiI = dadi.PhiManip.phi_1D(xx)
	## divergence event
	phiI = dadi.PhiManip.phi_1D_to_2D(xx, phiI)
	## set population sizes after split with different migration rates
	phiI = Integration.two_pops(phiI, xx, T, nu1, nu2, m12=me12, m21=me21)
	## Calculate the spectrum for correctly oriented loci
	fsIO = Spectrum.from_phi(phiI, ns, (xx,xx))
	## Calculate the spectrum for incorrectly oriented loci
	fsIM = dadi.Numerics.reverse_array(fsIO)
	
		# Spectrum for normally recombining regions
	## phi for equilibrium ancestral population
    phinr = PhiManip.phi_1D(xx)
	## divergence event
    phinr = PhiManip.phi_1D_to_2D(xx, phinr)
    ## set population sizes after split with different migration rates
    phinr = Integration.two_pops(phinr, xx, T, nu1, nu2, m12=m12, m21=m21)
	## Calculate the spectrum for correctly oriented loci
    fsnrO = Spectrum.from_phi(phinr, ns, (xx,xx))
    ## Calculate the spectrum for incorrectly oriented loci
	fsnrM = dadi.Numerics.reverse_array(fsnrO)
	
	# Spectrum for low-recombining regions
	## phi for equilibrium ancestral population
	philr = dadi.PhiManip.phi_1D(xx)
	## divergence event
	philr = dadi.PhiManip.phi_1D_to_2D(xx, philr)
	## set population sizes after split with different migration rates
	philr = Integration.two_pops(philr, xx, T, nu1*hrf, nu2*hrf, m12=m12, m21=m21)
	## Calculate the spectrum for correctly oriented loci
	fslrO = Spectrum.from_phi(philr, ns, (xx,xx))
	## Calculate the spectrum for incorrectly oriented loci
	fslrM = dadi.Numerics.reverse_array(fslrO)
	
	# Sum the spectra
    fs = O*(Q*fslrO+(1-Q)*fsnrO+P*fsNO+(1-P)*fsIO) + (1-O)*(Q*fslrM+(1-Q)*fsnrM+P*fsNM+(1-P)*fsIM)
	
    return fs


