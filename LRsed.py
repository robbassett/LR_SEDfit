from astropy.cosmology import WMAP9 as cosmo
import astropy.units as units
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

def fit_one_SED(gal_ID,base_dust_curve,filters,bpass_folder,photometry,EBVres=0.0025,Rv=2.74,maxage=5e8):

    # Determine which ages to use based on age of the universe at z=redshift
    redshift = photometry.data[str(gal_ID)]['redshift']
    ages = 10.**(np.array([6. + 0.1*_i for _i in range(51)]))
    cage = cosmo.age(redshift).to_value(units.year)
    tage = np.where(cage-maxage > ages)[0]
    nages = len(tage)
    
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
    wav = filter_wavs(filters,fit_filters)
    flux,error = fit_photometry(photometry.data[str(gal_ID)]['phot'],fit_filters)

    grid,tft = bpt.make_grid(bpass_spectra[0],filters,fit_filters,redshift)
    alosses = np.zeros(len(bpass_spectra))
    aopts   = np.zeros((len(bpass_spectra),grid.shape[0]))
    aEBVs   = np.zeros(len(bpass_spectra))
    for j,bps in enumerate(bpass_spectra):
        grid,tft = bpt.make_grid(bps,filters,fit_filters,redshift)

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
        alosses[j] = losses[tbest]
        aopts[j] = tOPT
        aEBVs[j] = tEBV

    abest = np.where(alosses == np.min(alosses))[0][0]
    EBV_out = aEBVs[abest]
    OPT_out = aopts[abest]
    
    swav,out_spec = bpt.out_spec(OPT_out,
                                 bpass_spectra[abest],
                                 redshift,
                                 base_dust_curve.wav,
                                 base_dust_curve.base_curve,
                                 EBV_out
                    )
    bpf,tf2 = bpt.BPASS_phot(swav,out_spec,fit_filters,filters,redshift)

    return {'wav':swav,
            'spec':out_spec,
            'BPASS_mod':bpass_spectra[abest],
            'theta_opt':OPT_out,
            'EBV':EBV_out,
            'attenuation_curve':np.exp((-1.)*(np.interp(swav,base_dust_curve.wav,base_dust_curve.base_curve)*EBV_out*Rv/1.086)),
            'z':redshift,
            'obs_phot':{
                'wav':wav/(1.+redshift),
                'flx':flux,
                'err':error,
                'fit':tFit
            }
        }

def fit_one_SED_big(gal_ID,base_dust_curve,filters,bpass_folder,photometry,EBVres=0.0025,Rv=2.74,maxage=5e8):

    # Determine which ages to use based on age of the universe at z=redshift
    redshift = photometry.data[str(gal_ID)]['redshift']
    ages = 10.**(np.array([6. + 0.1*_i for _i in range(51)]))
    cage = cosmo.age(redshift).to_value(units.year)
    tage = np.where(cage-maxage > ages)[0]
    nages = len(tage)
    
    # BPASS SPECTRA FILE NAMES/LOCATIONS
    bpass_spectra = glob.glob(f'{bpass_folder}/spectra*.gz')
    
    # CREATE EBV AND AGE ARRAYS
    EBVs = np.arange(0.,0.2,EBVres)
    ages = np.arange(0,nages,1).astype(int)+1

    # SETTING BOUNDS ON RETURNED MASSES FOR OUTPUT MODEL,
    # PREVENTS NEGATIVE MASS MODELS BEING INCLUDED
    bnd = [0.,1.e13]
    bounds = []
    for i in range(33*len(bpass_spectra)):
        bounds.append(tuple(bnd))
    bounds = tuple(bounds)

    # DETERMINE WHICH FILTERS ARE IN USE AND SET UP PHOTOMETRY FOR FITTING
    fit_filters = photometry.get_fit_filters(gal_ID,filters)
    wav = filter_wavs(filters,fit_filters)
    flux,error = fit_photometry(photometry.data[str(gal_ID)]['phot'],fit_filters)

    grid,tft = bpt.make_big_grid(bpass_spectra,filters,fit_filters,redshift)
    tFit = list(set(tft))
    theta = np.zeros(grid.shape[0])*1.e8
    losses = np.copy(EBVs)*0.
    topts = np.empty((len(EBVs),len(theta)))
    tfescs = np.zeros(len(EBVs))
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

    tbest = np.argmin(losses)
    tEBV = EBVs[tbest]
    tOPT = topts[tbest]
    
    swav,out_spec = bpt.out_spec_big(tOPT,
                                 bpass_spectra,
                                 redshift,
                                 base_dust_curve.wav,
                                 base_dust_curve.base_curve,
                                 tEBV,
                                 1.0
     )
    
    return {'wav':swav,
            'spec':out_spec,
            'BPASS_files':bpass_spectra,
            'theta_opt':tOPT,
            'EBV':tEBV,
            'attenuation_curve':np.exp((-1.)*(np.interp(swav,base_dust_curve.wav,base_dust_curve.base_curve)*tEBV*Rv/1.086)),
            'z':redshift,
            'obs_phot':{
                'wav':wav/(1.+redshift),
                'flx':flux,
                'err':error,
                'fit':tFit
            }
        }

if __name__ == '__main__':
    ID = 14759
    dc = dm.DustCurve('reddy16')
    flt = np.load('./COSMOS_filters.npy',allow_pickle=True).item()  
    phot = pht('myphot.dat','myphot_labels.dat')

    fit_one_SED(ID,dc,flt,'../../BPASS/bp2',phot)
