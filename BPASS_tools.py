import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u

# Convert BPASS model to microJansky/M_solar
def flux_conv(wav,L,z):
    # Speed of light in angstrom/sec
    c       = 2.99792458e18
    # Solar luminosity in ergs/s
    L_solar = 3.846e33
    # Lumninosity distance at z in cm
    DL = cosmo.luminosity_distance(z).to(u.cm).value

    # Luminosity in ergs/s/A for 1 solar mass
    out = (L*L_solar)/(1.e6)
    # convert to flux at z in ergs/s/cm^2/A
    out /= (4.*np.pi*DL*DL)

    # convert to microJansky
    out *= (3.34e4)*(wav*wav)*(1.e6)
    
    return out

# Get BPASS fluxes in each photometric filter
def BPASS_phot(wav,spec,filts,z):

    pwavs = []
    pflxs = []
    tfit  = []
    for i,k in enumerate(filts.keys()):

        pwavs.append(np.average(filts[k]['wav'],weights=filts[k]['T']))
        tms = np.interp(filts[k]['wav'],wav*(1.+z),spec)
        pflxs.append(np.average(tms,weights=filts[k]['T']))

        if np.min(filts[k]['wav']/(1.+z)) > 1216.:
            tfit.append(i)

    return np.array(pwavs),np.array(pflxs),tfit

# Build grid of BPASS models at fixed metallicity
def make_grid(tmd,filts,z,nSED=33):
    da = np.loadtxt(tmd)
    grid_out = np.empty((nSED,len(list(filts.keys()))))
    for i in range(nSED):
        filt_wavs,grid_out[i,:],tft = BPASS_phot(da[:,0],flux_conv(da[:,0],da[:,i+1],z),filts,z)

    return filt_wavs,grid_out,tft

# Create the output spectrum after running linear regression
def out_spec(theta,tmod,z,ll,bdc,ebv,Rv=2.74,nSED=33):
    da   = np.loadtxt(tmod)
    spec = np.empty(da.shape[0])
    dc = np.interp(da[:,0],ll,np.exp((-1.)*(bdc*ebv*Rv/1.086)))
    for ii in range(nSED):
        spec+=theta[ii]*flux_conv(da[:,0],da[:,ii+1],z)*dc
    return da[:,0],spec
