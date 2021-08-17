import matplotlib.pyplot as plt
import numpy as np
from pypeit import specobj
from pypeit import specobjs
from pypeit.core.coadd import multi_combspec
from IPython import embed
import skycalc_ipy

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


blue_exten = 10
red_exten = 25
fits_file = '../../DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Star/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits'

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
masks = np.ones_like(waves,dtype=bool)

wgmax = np.max(red_wave.value)
wgmin = np.min(blue_wave[blue_wave.value>10].value)

#Continium Continuity

red_spec.wavelength[red_spec.wavelength.value>10].value<blue_spec.wavelength[-1].value
overlap_top = blue_spec.wavelength[-1].value
overlap_mask = red_spec.wavelength[red_spec.wavelength.value>10].value < overlap_top
overlap_num = int(np.sum(overlap_mask))
blue_sum = np.sum(blue_spec.flux.value[-1*overlap_num:])
red_sum = np.sum(red_spec.flux[red_spec.wavelength.value>10].value[overlap_mask])
ratio = blue_sum/red_sum
fluxes[:,1] *= ratio

new_waves, new_flux, new_ivars, new_masks = multi_combspec(waves,fluxes,ivars,masks,wave_grid_max=wgmax,wave_grid_min=wgmin)

zero_skip = new_waves> 10
new_waves = new_waves[zero_skip]
new_flux = new_flux[zero_skip]
new_ivars = new_ivars[zero_skip]
new_masks = new_masks[zero_skip]

width=5
new_flux, new_ivars = ivarsmooth(new_flux,new_ivars,width)

new_sig = 1/np.sqrt(new_ivars)

sobj = specobj.SpecObj.from_arrays('MultiSlit',new_waves,new_flux,new_ivars)

new_sobjs = specobjs.SpecObjs()
new_sobjs.add_sobj(sobj)

subheader = {'test':'test'}

new_sobjs.write_to_fits(subheader,'star_spec.fits')