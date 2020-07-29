import numpy as np

def base_dust_curve():
    ll = np.arange(0.06,.56,.0001)
    tv = np.where((ll > .545)&(ll < .555))[0]
    tb = np.where((ll > .440)&(ll < .450))[0]

    raw_dc = dust_mod(ll)
    nv = np.mean(raw_dc[tb])-np.mean(raw_dc[tv])
    return ll*(1.e4),raw_dc/nv

def get_dust_mod(dmn):

    global dust_mod
    
    if dmn == 'reddy16':
        from dust_curves.reddy16 import reddy16 as dust_mod
    elif dmn == 'wild':
        from dust_curves.wild_av import wild_av as dust_mod
    elif dmn == 'calz00':
        from dust_curves.calz00b import calz00b as dust_mod
    elif dmn == 'SMC':
        from dust_curves.SMCdc import SMCdc as dust_mod
    else:
        print('\n\n')
        print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
        print('Invalid Dust Model Name, Please Choose From:')
        print('[reddy16,wild,calz00,SMC]')
        print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=')
        print('\n\n')

    
    print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
    print(f' Dust Model Set to {dmn}')
    print('-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')

class DustCurve(object):
    
    def __init__(self,mod_name):
        get_dust_mod(mod_name)
        self.wav,self.base_curve = base_dust_curve()

    def tmp(self):
        import matplotlib.pyplot as plt
        F = plt.figure()
        ax= F.add_subplot(111)
        ax.plot(self.wav,self.base_curve,'k-')
        plt.show()

def add_spec(sw,spec,dw,bdc,ebv,Rv=2.74):

    dc = np.exp((-1.)*(np.interp(sw,dw,bdc)*ebv*Rv/1.086))
    return spec*dc

def add_phot(pw,phot,dw,bdc,ebv,z,Rv=2.74):

    dc = np.exp((-1.)*(bdc*ebv*Rv/1.086))
    phot_out = np.copy(phot)
    for i in range(len(pw)):
        t = np.where(np.abs(dw-(pw[i]/(1.+z))) == np.min(np.abs(dw-(pw[i]/(1.+z)))))[0]
        phot_out[i] *= dc[t]

    return phot_out

def add_grid(pw,grd,dw,bdc,ebv,z,Rv=2.74):

    dc = np.exp((-1.)*(bdc*ebv*Rv/1.086))
    grid_out = np.copy(grd)
    for i in range(len(pw)):
        t = np.where(np.abs(dw-(pw[i]/(1.+z))) == np.min(np.abs(dw-(pw[i]/(1.+z)))))[0]
        grid_out[:,i] *= dc[t]

    return grid_out
