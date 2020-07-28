import numpy as np

def spec(sw,spec,dw,bdc,ebv,Rv=2.74):

    dc = np.exp((-1.)*(np.interp(sw,dw,bdc)*ebv*Rv/1.086))
    return spec*dc

def phot(pw,phot,dw,bdc,ebv,z,Rv=2.74):

    dc = np.exp((-1.)*(bdc*ebv*Rv/1.086))
    phot_out = np.copy(phot)
    for i in range(len(pw)):
        t = np.where(np.abs(dw-(pw[i]/(1.+z))) == np.min(np.abs(dw-(pw[i]/(1.+z)))))[0]
        phot_out[i] *= dc[t]

    return phot_out

def grid(pw,grd,dw,bdc,ebv,z,Rv=2.74):

    dc = np.exp((-1.)*(bdc*ebv*Rv/1.086))
    grid_out = np.copy(grd)
    for i in range(len(pw)):
        t = np.where(np.abs(dw-(pw[i]/(1.+z))) == np.min(np.abs(dw-(pw[i]/(1.+z)))))[0]
        grid_out[:,i] *= dc[t]

    return grid_out
