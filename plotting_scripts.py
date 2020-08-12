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
    lbt = np.zeros(len(topt))
    for i in range(len(zs)):
        tuage = (uage-(ages[i])).to_value(units.Gyr)
        
        lbt[i] = uage.value-ages[i].value
        zs[i] = z_at_value(cosmo.age,tuage*units.Gyr)

    age_ticks = [10.**_e for _e in [6.,7.,8.,9.]]
    zticks = np.zeros(len(age_ticks))
    for i in range(len(zticks)):
        tuage = (uage-(age_ticks[i]*units.year)).to_value(units.Gyr)
        zticks[i] = round(100.*z_at_value(cosmo.age,tuage*units.Gyr))/100.
    
        
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
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Stellar Mass')

    ax= F.add_subplot(212)
    ax.semilogx(ages,sfr,'ko')
    ax.plot(ages,sfr,'k-')
    xa = ax.twiny()
    xa.semilogx(ages,sfr,'ko')
    xa.set_xticks(age_ticks)
    xa.set_xticklabels(zticks)
    xa.set_xlabel('redshift')
    ax.set_xlabel('lookback time (years)')
    ax.set_ylabel('SFR')
    plt.tight_layout()
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

if __name__ == '__main__':
    dat = np.load('test_output.npy',allow_pickle=True).item()
    for k in dat.keys():
       SFhist_plot(dat,k)
        
