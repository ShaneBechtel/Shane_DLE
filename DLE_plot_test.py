import argparse

import matplotlib.pyplot as plt
import numpy as np
from pypeit import specobjs
from pypeit.core.coadd import multi_combspec
from astropy.io import fits
from pypeit import utils
from IPython import embed

parser = argparse.ArgumentParser()
parser.add_argument('--file', help='Spectrum Fits File')
parser.add_argument('--lines', default='yes', help='Plot spectral lines')
parser.add_argument('--redshift', default=0.0, help='Redshift of emission lines')
parser.add_argument('--blue', help='Exten value for target in Blue Detector')
parser.add_argument('--red', help='Exten value for target in Blue Detector')
parser.add_argument('--width', help='Width of boxcar smoothing')
args = parser.parse_args()
#Example Call
#python DLE_plot_test.py --file ../../DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Both/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits --redshift 4.95 --blue 9 --red 30 --width 5

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
    halfwindow = int(np.floor((np.round(window) - 1) / 2))
    shiftarr = np.zeros((nflux, 2 * halfwindow + 1))
    shiftivar = np.zeros((nflux, 2 * halfwindow + 1))
    shiftindex = np.zeros((nflux, 2 * halfwindow + 1))
    indexarr = np.arange(nflux)
    indnorm = np.outer(indexarr, (np.zeros(2 * halfwindow + 1) + 1))
    for i in np.arange(-halfwindow, halfwindow + 1, dtype=int):
        shiftarr[:, i + halfwindow] = np.roll(flux, i)
        shiftivar[:, i + halfwindow] = np.roll(ivar, i)
        shiftindex[:, i + halfwindow] = np.roll(indexarr, i)
    wh = (np.abs(shiftindex - indnorm) > (halfwindow + 1))
    shiftivar[wh] = 0.0
    outivar = np.sum(shiftivar, axis=1)
    nzero, = np.where(outivar > 0.0)
    zeroct = len(nzero)
    smoothflux = np.sum(shiftarr * shiftivar, axis=1)
    if (zeroct > 0):
        smoothflux[nzero] = smoothflux[nzero] / outivar[nzero]
    else:
        smoothflux = np.roll(flux, 2 * halfwindow + 1)  # kill off NAN’s
    return (smoothflux, outivar)


blue_exten = int(args.blue)
red_exten = int(args.red)
fits_file = args.file

sobjs = specobjs.SpecObjs.from_fitsfile(fits_file, chk_version=False)
blue_spec = sobjs[blue_exten - 1].to_xspec1d(extraction='OPT', fluxed=False)
red_spec = sobjs[red_exten - 1].to_xspec1d(extraction='OPT', fluxed=False)

if len(blue_spec.wavelength) < len(red_spec.wavelength):

    red_wave = red_spec.wavelength
    red_flux = red_spec.flux
    red_ivar = red_spec.ivar

    blue_wave = np.zeros_like(red_wave)
    blue_flux = np.zeros_like(red_wave)
    blue_ivar = np.zeros_like(red_wave)

    diff = len(red_spec.wavelength) - len(blue_spec.wavelength)

    blue_wave[diff:] = blue_spec.wavelength
    blue_flux[diff:] = blue_spec.flux
    blue_ivar[diff:] = blue_spec.ivar


else:

    blue_wave = blue_spec.wavelength
    blue_flux = blue_spec.flux
    blue_ivar = blue_spec.ivar

    red_wave = np.zeros_like(blue_wave)
    red_flux = np.zeros_like(blue_wave)
    red_ivar = np.zeros_like(blue_wave)

    diff = len(blue_spec.wavelength) - len(red_spec.wavelength)

    red_wave[diff:] = red_spec.wavelength
    red_flux[diff:] = red_spec.flux
    red_ivar[diff:] = red_spec.ivar

waves = np.zeros((len(blue_wave), 2))
fluxes = np.zeros((len(blue_wave), 2))
ivars = np.zeros((len(blue_wave), 2))

waves[:, 0] = blue_wave
waves[:, 1] = red_wave
fluxes[:, 0] = blue_flux
fluxes[:, 1] = red_flux
ivars[:, 0] = blue_ivar
ivars[:, 1] = red_ivar
masks = ivars > 0.0

wgmax = np.max(red_wave.value)
wgmin = np.min(blue_wave[blue_wave.value > 10].value)

# Continium Continuity

overlap_top = blue_spec.wavelength[-1].value
overlap_mask = red_spec.wavelength[red_spec.wavelength.value > 10].value < overlap_top
overlap_num = int(np.sum(overlap_mask))
blue_sum = np.sum(blue_spec.flux.value[-1 * overlap_num:])
red_sum = np.sum(red_spec.flux[red_spec.wavelength.value > 10].value[overlap_mask])
ratio = blue_sum / red_sum
fluxes[:, 1] *= ratio

new_waves, new_flux, new_ivars, new_masks = multi_combspec(waves, fluxes, ivars, masks, wave_grid_max=wgmax,
                                                           wave_grid_min=wgmin)

zero_skip = new_waves > 10
new_waves = new_waves[zero_skip]
new_flux = new_flux[zero_skip]
new_ivars = new_ivars[zero_skip]
new_masks = new_masks[zero_skip]

width = int(args.width)
new_flux, new_ivars = ivarsmooth(new_flux, new_ivars, width)

new_sig = 1 / np.sqrt(new_ivars)

# Quasar Redshifts

lya = 1215.67
J14_wave = lya*(1+4.709)
J16_wave = lya*(1+3.8101)

# Atmospheric Effects

tell_hdu = fits.open('star_spec_tellcorr.fits')
tell_waves = tell_hdu[1].data['wave']
tell_spec = tell_hdu[1].data['telluric']
tell_corr = np.interp(new_waves,tell_waves,tell_spec,left=1,right=1)

flux_corr = new_flux/(tell_corr + (tell_corr == 0))
ivar_corr = (tell_corr > 0.0) * new_ivars * tell_corr * tell_corr
mask_corr = (tell_corr > 0.0) * new_masks
sig_corr = np.sqrt(utils.inverse(ivar_corr))

# Composite Spectrum
comp_hdu = fits.open('Composite/MUSYC_LBGonly_stack.fits')
comp_data = comp_hdu[0].data[0, 0, :]
comp_redshift = 3381.89 / 1216.00 - 1
comp_waves = (3200.43 + 0.86 * np.arange(len(comp_hdu[0].data[0, 0, :]))) / (1 + comp_redshift)
comp_data = (np.mean(new_flux) / np.mean(comp_data)) * comp_data

fig, ax = plt.subplots()
redshift = float(args.redshift)
trans = ax.get_xaxis_transform()
ax.step(new_waves, flux_corr, 'k', where='mid')
ax.plot(new_waves, sig_corr, 'r:')
ax.plot(tell_waves, tell_spec, 'g--', transform=trans, alpha=0.5)
ax.plot((1 + redshift) * comp_waves, comp_data, 'orange', alpha=0.5)

ax.vlines(J16_wave, new_flux.min(), new_flux.max(), 'y', linestyles='--', alpha=0.5)
plt.text(J16_wave, .85, 'J1630', transform=trans, backgroundcolor='0.75')
ax.vlines(J14_wave, new_flux.min(), new_flux.max(), 'y', linestyles='--', alpha=0.5)
plt.text(J14_wave, .85, 'J1438A', transform=trans, backgroundcolor='0.75')

ln_flag = bool(args.lines[0] == 'n')

if not ln_flag:

    line_file = open('gal_vac.lst')
    ln_lst = line_file.readlines()
    line_file.close()
    Lines = {}
    for ln in ln_lst:
        lam = float(ln[:8])
        name = ln[-9:-1]
        if name in Lines.keys():
            Lines[name] = np.concatenate([Lines[name], np.array([lam])])
        else:
            Lines[name] = np.array([lam])

    for tup in Lines.items():

        if np.shape(tup[1])[0] == 1:
            z_wavelength = (1 + redshift) * tup[1][0]
            ax.vlines(z_wavelength, flux_corr.min(), flux_corr.max(), 'b', linestyles='--', alpha=0.5)
            plt.text(z_wavelength, .85, tup[0], transform=trans, backgroundcolor='0.75')
        else:
            for l in tup[1]:
                z_wavelength = (1 + redshift) * l
                ax.vlines(z_wavelength, new_flux.min(), new_flux.max(), 'b', linestyles='--', alpha=0.5)
                plt.text(z_wavelength, .85, tup[0], transform=trans, backgroundcolor='0.75')
plt.xlim(new_waves.min(), new_waves.max())

plt.show()
