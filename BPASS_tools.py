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
def BPASS_phot(wav,spec,fnames,filts,z):

    pflxs = []
    tfit  = []
    for i,k in enumerate(fnames):
        fwav = filts[k]['wav']
        fT = filts[k]['T']
        t = np.where(fT > 0.01)[0]

        pwav = (np.average(fwav,weights=fT))/(1.+z)
        tms = np.interp(fwav,wav*(1.+z),spec)
        pflxs.append(np.average(tms,weights=fT))
        minw,maxw = np.min(fwav[t])/(1.+z),np.max(fwav[t])/(1.+z)

        # Leave bands below LAF and beyond 1.5 micron out
        if  minw > 1216.:
            if minw < 5007. and maxw > 5007.:
                pass
            else:
                tfit.append(i)

    return np.array(pflxs),tfit

# Build grid of BPASS models at fixed metallicity
def make_grid(tmd,filts,fnames,z,nSED=33):
    da = np.loadtxt(tmd)
    grid_out = np.empty((nSED,len(fnames)))
    phot_out = np.empty((2,len(fnames)))
    for i in range(nSED):
        grid_out[i,:],tft = BPASS_phot(da[:,0],flux_conv(da[:,0],da[:,i+1],z),fnames,filts,z)
            
    return grid_out,tft

def make_big_grid(ffs,filts,fnames,z,nSED=33):

    grid_out = np.empty((nSED*len(ffs),len(fnames)))
    for j,fn in enumerate(ffs):
        da = np.loadtxt(fn)
        for i in range(nSED):
            grid_out[j*nSED + i,:],tft = BPASS_phot(da[:,0],flux_conv(da[:,0],da[:,i+1],z),fnames,filts,z)
            
    return grid_out,tft

# Create the output spectrum after running linear regression
def out_spec(theta,tmod,z,ll,bdc,ebv,Rv=2.74,nSED=33):
    da   = np.loadtxt(tmod)
    spec = np.empty(da.shape[0])
    dc = np.interp(da[:,0],ll,np.exp((-1.)*(bdc*ebv*Rv/1.086)))
    for ii in range(nSED):
        spec+=theta[ii]*flux_conv(da[:,0],da[:,ii+1],z)*dc
    return da[:,0],spec

def out_spec_big(tOPT,mods,z,ll,bdc,ebv,fesc,Rv=2.74,nSED=33):
    
    da   = np.loadtxt(mods[0])
    spec = np.empty(da.shape[0])
    dc = np.interp(da[:,0],ll,np.exp((-1.)*(bdc*ebv*Rv/1.086)))
    topt = np.reshape(tOPT,(len(mods),nSED))
    for i,mod in enumerate(mods):
        td = np.loadtxt(mod)
        for j in range(nSED):
            spec += topt[i][j]*flux_conv(td[:,0],td[:,j+1],z)

    return da[:,0],spec
