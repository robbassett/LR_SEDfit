import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as units
from astropy.cosmology import z_at_value

def SFhist_plot(dat,k):
    topt = dat[k]['theta_opt']
    z = dat[k]['z']
    lages = np.zeros(len(topt))
    for i in range(len(lages)): lages[i] = 6.0+0.1*i
    ages = 10.**lages
    ages = ages*units.year
    uage = cosmo.age(z).to_value(units.year)
    uage = uage*units.year
    zs = np.zeros(len(topt))
    for i in range(len(zs)):
        tuage = (uage-(ages[i])).to_value(units.Gyr)
        zs[i] = z_at_value(cosmo.age,tuage*units.Gyr)
        
    sfr = np.zeros(len(topt))
    for i in range(len(ages)):
        if i == 0:
            lage = 0.
        else:
            lage = ages[i-1]
        dage = ages[i]-lage
        sfr[i] = topt[i]/dage.value

    F = plt.figure()
    ax= F.add_subplot(211)
    ax.loglog(ages,topt,'ko')
    ax.plot(ages,topt,'k-')
    ax.set_xlabel('Log(age)')
    ax.set_ylabel('Stellar Mass')

    ax= F.add_subplot(212)
    ax.plot(zs,sfr,'ko')
    ax.plot(zs,sfr,'k-')
    ax.set_xlabel('redshift')
    ax.set_ylabel('SFR')
    plt.show()

def SED_plot(dat,k):

    tmo = dat[k]
    swav,out_spec = tmo['wav'],tmo['spec']
    dust_curve = tmo['attenuation_curve']
    int_spec = out_spec/dust_curve
    pht = tmo['obs_phot']
    pwav,pflx,perr = pht['wav'],pht['flx'],pht['err']
    pft = pht['fit']
    redshift = tmo['z']
        
    F = plt.figure()
    ax = F.add_subplot(111)
    ax.plot(swav,int_spec,'g-',lw=.5,label='Intrinsic SED')
    ax.plot(swav,out_spec,'k-',lw=.5,label='Attenuated SED')
    ax.plot(pwav,pflx,'ro',mfc='None',label='Photometry (not fit)')
    ax.errorbar(pwav[pft],pflx[pft],yerr=perr[pft],fmt='o',ms=5,c='r',label='Photometry (fit)')
    ax.set_xlim(600,13999)
    ax.set_ylim(0,1.1*np.max(pflx[pft]))
    ax.legend(ncol=2)
    plt.show()
