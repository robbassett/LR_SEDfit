import numpy as np

# CONVERT ANYTHING TO MICRO-JANSKY

def AB_mag(mag):
    fnu = 10.**((-0.4)*(mag-8.9))
    return fnu*1.e6

def ZFOURGE_flx(flx):
    return 0.3631*flx

def Jy(flx):
    return flx*1.e6
