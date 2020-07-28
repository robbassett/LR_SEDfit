
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

# Build grid of BPASS models at fixed metallicity
def make_grid(tmd,filts,z,tm,nSED=33):
    da = np.loadtxt(tmd)
    grid_out = np.empty((nSED,len(tm)))
    for i in range(nSED):
        filt_wavs,grid_out[i,:],tft = BPASS_phot(da[:,0],flux_conv(da[:,0],da[:,i+1],z),filts,z,tm)

    return filt_wavs,grid_out,tft
