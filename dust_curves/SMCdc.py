import numpy as np

def Drude(x,x0,gamma):
    num = x**2.
    den = (((x**2.)-(x0**2.))**2.) + ((x**2.)*(gamma**2.))
    return num/den

def UVfx(x):
    return (0.5392*((x-5.9)**2.))+(0.05644*((x-5.9)**3.))


def SMCdc(lam):
    SMCpar = [-.96,1.18,2.57,4.71,1.0,0.1]
    xx     = lam**(-1.)
    klam   = np.zeros(len(lam))

    for i in range(len(lam)):
        Dt     = Drude(xx[i],SMCpar[3],SMCpar[4])
        if xx[i] < 5.9:
            Ft = UVfx(xx[i])
        else:
            Ft = 0.
        klam[i]=(SMCpar[0] + SMCpar[1]*xx[i] + SMCpar[1]*Dt + SMCpar[5]*Ft)

    return klam
