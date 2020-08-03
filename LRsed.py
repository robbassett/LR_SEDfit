import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import glob
import sys

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

import BPASS_tools as bpt
import LR_mod as lrm
import dust_curves as dm
from load_photometry import Photometry as pht

def filter_wavs(filts,fnames):
    fwav = np.empty(len(fnames))
    for i,k in enumerate(fnames):
        fwav[i] = np.average(filts[k]['wav'],weights=filts[k]['T'])
    return fwav

def fit_photometry(tphot,fnames):
    flux,err = np.empty(len(fnames)),np.empty(len(fnames))
    for i,k in enumerate(fnames):
        flux[i] = tphot[k][0]
        err[i]  = tphot[k][1]

    return flux,err

def fit_one_SED(gal_ID,base_dust_curve,filters,bpass_folder,photometry,EBVres=0.0025,nages=33,Rv=2.74):

    # BPASS SPECTRA FILE NAMES/LOCATIONS
    bpass_spectra = glob.glob(f'{bpass_folder}/spectra*.gz')
    
    # CREATE EBV AND AGE ARRAYS
    EBVs = np.arange(0.,0.2,EBVres)
    ages = np.arange(0,nages,1).astype(int)+1

    # SETTING BOUNDS ON RETURNED MASSES FOR OUTPUT MODEL,
    # PREVENTS NEGATIVE MASS MODELS BEING INCLUDED
    bnd = [0.,1.e13]
    bounds = []
    for i in range(33):
        bounds.append(tuple(bnd))
    bounds = tuple(bounds)

    # DETERMINE WHICH FILTERS ARE IN USE AND SET UP PHOTOMETRY FOR FITTING
    fit_filters = photometry.get_fit_filters(gal_ID,filters)
    redshift = photometry.data[str(gal_ID)]['redshift']
    wav = filter_wavs(filters,fit_filters)
    flux,error = fit_photometry(photometry.data[str(gal_ID)]['phot'],fit_filters)
    
    bpmn=5 # JUST FIT ONE MODEL NOW, WILL INCLUDE LOOP OVER ALL MODELS LATER
    grid,tft = bpt.make_grid(bpass_spectra[bpmn],filters,fit_filters,redshift)

    tFit = list(set(tft))
    theta = np.zeros(grid.shape[0])*1.e8
    losses = np.copy(EBVs)*0.
    topts = np.empty((len(EBVs),len(theta)))
    for i,EBV in enumerate(EBVs):

        dust_grid = dm.add_grid(wav,grid,base_dust_curve.wav,base_dust_curve.base_curve,EBV,redshift,Rv=Rv)
        Res = opt.minimize(fun = lrm.cost,
                           x0 = theta,
                           args = (dust_grid[:,tFit],flux[tFit],error[tFit]),
                           method = 'TNC',
                           jac = lrm.grad,
                           bounds = bounds
              )

        theta_opt = Res.x
        SED_opt = lrm.h(theta_opt,dust_grid)
        topts[i,:] = theta_opt
        losses[i] = lrm.SED_loss(SED_opt,flux[tFit],error[tFit])

    tbest = np.where(losses == np.min(losses))[0]
    tEBV = EBVs[tbest]
    tOPT = topts[tbest][0,:]
    swav,out_spec = bpt.out_spec(tOPT,
                                 bpass_spectra[bpmn],
                                 redshift,
                                 base_dust_curve.wav,
                                 base_dust_curve.base_curve,
                                 tEBV
                    )
    bpf,tf2 = bpt.BPASS_phot(swav,out_spec,fit_filters,filters,redshift)

    F = plt.figure()
    ax = F.add_subplot(111)
    ax.plot(swav,out_spec,'k-',lw=.5)
    ax.plot(wav[tFit]/(redshift+1.),bpf[tFit],'g^')
    ax.plot(wav/(redshift+1.),bpf,'g^',mfc='None')
    ax.plot(wav/(redshift+1.),flux,'bo',mfc='None')
    ax.plot(wav[tFit]/(redshift+1.),flux[tFit],'bo')
    ax.set_xlim(600,13999)
    plt.show()
    
