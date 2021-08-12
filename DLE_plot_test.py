import argparse

import matplotlib.pyplot as plt
import numpy as np
from linetools.spectra.xspectrum1d import XSpectrum1D
from pypeit import specobjs
from pypeit.core.coadd import multi_combspec
from IPython import embed

parser = argparse.ArgumentParser()
parser.add_argument('--file', help='Spectrum Fits File')
parser.add_argument('--redshift', help='Redshift of emission lines')
parser.add_argument('--blue', help='Exten value for target in Blue Detector ')
parser.add_argument('--red', help='Exten value for target in Blue Detector ')
parser.add_argument('--width', help='Width of boxcar smoothing ')
args = parser.parse_args()


def ivarsmooth(flux, ivar, window):
    '''
    Boxcar smoothign of width window with ivar weights
    Args:
        flux:
        ivar:
        window:
    Returns:
    '''
    nflux = (flux.shape)[0]
    halfwindow = int(np.floor((np.round(window) - 1)/2))
    shiftarr = np.zeros((nflux, 2*halfwindow + 1))
    shiftivar = np.zeros((nflux, 2*halfwindow + 1))
    shiftindex = np.zeros((nflux, 2*halfwindow + 1))
    indexarr = np.arange(nflux)
    indnorm = np.outer(indexarr,(np.zeros(2 *halfwindow + 1) + 1))
    for i in np.arange(-halfwindow,halfwindow + 1,dtype=int):
        shiftarr[:,i+halfwindow] = np.roll(flux,i)
        shiftivar[:, i+halfwindow] = np.roll(ivar, i)
        shiftindex[:, i+halfwindow] = np.roll(indexarr, i)
    wh = (np.abs(shiftindex - indnorm) > (halfwindow+1))
    shiftivar[wh]=0.0
    outivar = np.sum(shiftivar,axis=1)
    nzero, = np.where(outivar > 0.0)
    zeroct=len(nzero)
    smoothflux = np.sum(shiftarr * shiftivar, axis=1)
    if(zeroct > 0):
        smoothflux[nzero] = smoothflux[nzero]/outivar[nzero]
    else:
        smoothflux = np.roll(flux, 2*halfwindow + 1) # kill off NAN’s
    return (smoothflux, outivar)


blue_exten = int(args.blue)
red_exten = int(args.red)
fits_file = args.file

sobjs = specobjs.SpecObjs.from_fitsfile(fits_file, chk_version=False)
blue_spec = sobjs[blue_exten-1].to_xspec1d(extraction='OPT', fluxed=False)
red_spec = sobjs[red_exten-1].to_xspec1d(extraction='OPT', fluxed=False)

if len(blue_spec.wavelength)>len(red_spec.wavelength):
    diff = len(blue_spec.wavelength)-len(red_spec.wavelength)

    blue_wave = blue_spec.wavelength[diff:]
    blue_flux = blue_spec.flux[diff:]
    blue_ivar = blue_spec.ivar[diff:]

    red_wave = red_spec.wavelength
    red_flux = red_spec.flux
    red_ivar = red_spec.ivar
else:
    diff = len(red_spec.wavelength)-len(blue_spec.wavelength)
    length = len(red_spec.wavelength)

    blue_wave = blue_spec.wavelength
    blue_flux = blue_spec.flux
    blue_ivar = blue_spec.ivar

    red_wave = red_spec.wavelength[:length-diff]
    red_flux = red_spec.flux[:length-diff]
    red_ivar = red_spec.ivar[:length-diff]

waves = np.zeros((len(blue_wave),2))
fluxes = np.zeros((len(blue_wave),2))
ivars = np.zeros((len(blue_wave),2))

waves[:,0] = blue_wave
waves[:,1] = red_wave
fluxes[:,0] = blue_flux
fluxes[:,1] = red_flux
ivars[:,0] = blue_ivar
ivars[:,1] = red_ivar
masks = np.ones_like(waves,dtype=int)

wgmax = np.max(red_wave.value)
wgmin = np.min(blue_wave[blue_wave.value>10].value)

new_waves, new_flux, new_ivars, new_masks = multi_combspec(waves,fluxes,ivars,masks,wave_grid_max=wgmax,wave_grid_min=wgmin)

zero_skip = new_waves> 10
new_waves = new_waves[zero_skip]
new_flux = new_flux[zero_skip]
new_ivars = new_ivars[zero_skip]
new_masks = new_masks[zero_skip]

'''
new_wavelength = np.concatenate((blue_spec.wavelength[zero_skip],red_spec.wavelength[red_spec.wavelength>blue_spec.wavelength[-1]]))
new_flux = np.concatenate((blue_spec.flux[zero_skip],red_spec.flux[red_spec.wavelength>blue_spec.wavelength[-1]]))
new_sig = np.concatenate((blue_spec.sig[zero_skip],red_spec.sig[red_spec.wavelength>blue_spec.wavelength[-1]]))
new_spec = XSpectrum1D.from_tuple((new_wavelength, new_flux, new_sig), verbose=False)
'''

width=int(args.width)
new_flux, new_ivars = ivarsmooth(new_flux,new_ivars,width)

new_sig = new_ivars/np.sqrt(new_ivars)

# Line List
redshift = float(args.redshift)
Lines = {"C III":[1175.71], "Si II":[1190,1260.42], "Lya":[1215.670],
         "N V":[1240.81], "O I":[1305.53], "C II":[1335.31]}


fig, ax = plt.subplots()

trans = ax.get_xaxis_transform()
ax.plot(new_waves,new_flux)

for tup in Lines.items():

    if np.shape(tup[1])[0]==1:
        z_wavelength = (1+redshift)*tup[1][0]
        ax.vlines(z_wavelength,new_flux.min(),new_flux.max(),'k',linestyles='--')
        plt.text(z_wavelength, .85, tup[0], transform=trans,backgroundcolor='0.75')
    else:
        for l in tup[1]:
            z_wavelength = (1 + redshift) * l
            ax.vlines(z_wavelength, new_flux.min(), new_flux.max(), 'k', linestyles='--')
            plt.text(z_wavelength, .85, tup[0], transform=trans, backgroundcolor='0.75')

plt.show()