import argparse

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
import numpy as np
from glob import glob
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit import sensfunc
from pypeit.core.coadd import multi_combspec
from pypeit.core import flux_calib
from astropy.io import fits
from pypeit import utils
from scipy.interpolate import interp1d
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed

parser = argparse.ArgumentParser()
parser.add_argument('--file_path', help='Coadd directory path')
parser.add_argument('--blue', help='Exten value for target in Blue Detector')
parser.add_argument('--red', help='Exten value for target in Blue Detector')
parser.add_argument('--width', help='Width of boxcar smoothing')
parser.add_argument('--channel', default=2, help='Which channel to include for 2D image in figure')
parser.add_argument('--contnorm', default='yes', help='Normalize the continuum flux?')
args = parser.parse_args()


# Example Call
# python DLE_general_plot.py --file /home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Both/Science_coadd/ --blue 9 --red 30 --width 5 --channel 1

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
file_path=args.file_path
if file_path[-1] != '/':
    file_path += '/'
files = glob(file_path+'*')
spec1d_file = files[0]
spec2d_file = files[1]

sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=True)
blue_spec = sobjs[blue_exten - 1].to_xspec1d(extraction='OPT', fluxed=True)
red_spec = sobjs[red_exten - 1].to_xspec1d(extraction='OPT', fluxed=True)

if len(blue_spec.wavelength) < len(red_spec.wavelength):

    red_wave = red_spec.wavelength.value
    red_flux = red_spec.flux.value
    red_ivar = red_spec.ivar.value

    blue_wave = np.zeros_like(red_wave)
    blue_flux = np.zeros_like(red_wave)
    blue_ivar = np.zeros_like(red_wave)

    diff = len(red_spec.wavelength) - len(blue_spec.wavelength)

    blue_wave[diff:] = blue_spec.wavelength.value
    blue_flux[diff:] = blue_spec.flux.value
    blue_ivar[diff:] = blue_spec.ivar.value


else:

    blue_wave = blue_spec.wavelength.value
    blue_flux = blue_spec.flux.value
    blue_ivar = blue_spec.ivar.value

    red_wave = np.zeros_like(blue_wave)
    red_flux = np.zeros_like(blue_wave)
    red_ivar = np.zeros_like(blue_wave)

    diff = len(blue_spec.wavelength) - len(red_spec.wavelength)

    red_wave[diff:] = red_spec.wavelength.value
    red_flux[diff:] = red_spec.flux.value
    red_ivar[diff:] = red_spec.ivar.value

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

wgmax = np.max(red_wave)
wgmin = np.min(blue_wave[blue_wave > 10])

cont_flag = bool(args.contnorm[0] == 'n')

if not cont_flag:
    # Continium Continuity

    overlap_top = blue_wave[-1]
    overlap_mask = red_wave[red_wave > 10] < overlap_top
    overlap_num = int(np.sum(overlap_mask))
    blue_sum = np.sum(blue_flux[-1 * overlap_num:])
    red_sum = np.sum(red_flux[red_wave > 10][overlap_mask])
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
J14_wave = lya * (1 + 4.709)
J16_wave = lya * (1 + 3.8101)

# Atmospheric Effects #TODO Don't hard code in this telluric file.
tell_hdu = fits.open('star_spec_tellcorr.fits')
tell_waves = tell_hdu[1].data['wave']
tell_spec = tell_hdu[1].data['telluric']
tell_corr = np.interp(new_waves, tell_waves, tell_spec, left=1, right=1)

flux_corr = new_flux / (tell_corr + (tell_corr == 0))
ivar_corr = (tell_corr > 0.0) * new_ivars * tell_corr * tell_corr
mask_corr = (tell_corr > 0.0) * new_masks
sig_corr = np.sqrt(utils.inverse(ivar_corr))

# Composite Spectrum
comp_file = open('deimos_z4composite.dat')
comp_lines = comp_file.readlines()

comp_data = []

for i in range(len(comp_lines)):
    data_hold = [float(x) for x in comp_lines[i].split()]

    comp_data.append(data_hold)

comp_data = np.array(comp_data)

comp_waves = comp_data[:, 0]
comp_flux = comp_data[:, 1]
comp_flux = (np.median(new_flux) / np.nanmedian(comp_flux)) * comp_flux


# 2D Image
det = sobjs[blue_exten - 1].DET
spec2DObj_blue = spec2dobj.Spec2DObj.from_file(spec2d_file, det, chk_version=False)
spec2DObj_red = spec2dobj.Spec2DObj.from_file(spec2d_file, det+4, chk_version=False)
channel = int(args.channel)


wave_ind_blue = int(np.round(sobjs[blue_exten - 1].SPAT_PIXPOS))
wave_ind_red = int(np.round(sobjs[red_exten - 1].SPAT_PIXPOS))

red_flag = False
blue_flag = False

if wave_ind_blue>wave_ind_red:
    wave_diff = wave_ind_blue-wave_ind_red
    red_flag = True
    spat_buffer = np.zeros((spec2DObj_red.sciimg.shape[0],wave_diff))
    spat_max = np.min([spec2DObj_red.sciimg.shape[1]+wave_diff,spec2DObj_blue.sciimg.shape[1]])

elif wave_ind_blue<wave_ind_red:
    wave_diff = wave_ind_red - wave_ind_blue
    blue_flag = True
    spat_buffer = np.zeros((spec2DObj_blue.sciimg.shape[0], wave_diff))
    spat_max = np.min([spec2DObj_red.sciimg.shape[1],spec2DObj_blue.sciimg.shape[1]+wave_diff])


if red_flag:
    if channel == 0:
        img_data_blue = spec2DObj_blue.sciimg
        img_data_red = spec2DObj_red.sciimg
        img_data_red = np.append(spat_buffer,img_data_red,axis=1)
        img_data = np.append(img_data_blue[:,:spat_max],img_data_red[:,:spat_max],axis=0)
        vmax = 15
        vmin = -3
    elif channel == 1:
        gpm_blue = spec2DObj_blue.bpmmask == 0
        img_data_blue = (spec2DObj_blue.sciimg - spec2DObj_blue.skymodel) * gpm_blue
        gpm_red = spec2DObj_red.bpmmask == 0
        img_data_red = (spec2DObj_red.sciimg - spec2DObj_red.skymodel) * gpm_red
        img_data_red = np.append(spat_buffer, img_data_red, axis=1)
        img_data = np.append(img_data_blue[:, :spat_max], img_data_red[:, :spat_max], axis=0)
        vmax = 15
        vmin = -3
    elif channel == 2:
        gpm_blue = spec2DObj_blue.bpmmask == 0
        img_data_blue = (spec2DObj_blue.sciimg - spec2DObj_blue.skymodel) * np.sqrt(spec2DObj_blue.ivarmodel) * gpm_blue
        gpm_red = spec2DObj_red.bpmmask == 0
        img_data_red = (spec2DObj_red.sciimg - spec2DObj_red.skymodel) * np.sqrt(spec2DObj_red.ivarmodel) * gpm_red
        img_data_red = np.append(spat_buffer, img_data_red, axis=1)
        img_data = np.append(img_data_blue[:, :spat_max], img_data_red[:, :spat_max], axis=0)
        vmax = 4
        vmin = -1
    else:
        raise ValueError('Expected channel value of 0, 1, or 2')
elif blue_flag:
    if channel == 0:
        img_data_blue = spec2DObj_blue.sciimg
        img_data_red = spec2DObj_red.sciimg
        img_data_blue = np.append(spat_buffer,img_data_blue,axis=1)
        img_data = np.append(img_data_blue[:,:spat_max],img_data_red[:,:spat_max],axis=0)
        vmax = 15
        vmin = -3
    elif channel == 1:
        gpm_blue = spec2DObj_blue.bpmmask == 0
        img_data_blue = (spec2DObj_blue.sciimg - spec2DObj_blue.skymodel) * gpm_blue
        gpm_red = spec2DObj_red.bpmmask == 0
        img_data_red = (spec2DObj_red.sciimg - spec2DObj_red.skymodel) * gpm_red
        img_data_blue = np.append(spat_buffer, img_data_blue, axis=1)
        img_data = np.append(img_data_blue[:, :spat_max], img_data_red[:, :spat_max], axis=0)
        vmax = 15
        vmin = -3
    elif channel == 2:
        gpm_blue = spec2DObj_blue.bpmmask == 0
        img_data_blue = (spec2DObj_blue.sciimg - spec2DObj_blue.skymodel) * np.sqrt(spec2DObj_blue.ivarmodel) * gpm_blue
        gpm_red = spec2DObj_red.bpmmask == 0
        img_data_red = (spec2DObj_red.sciimg - spec2DObj_red.skymodel) * np.sqrt(spec2DObj_red.ivarmodel) * gpm_red
        img_data_blue = np.append(spat_buffer, img_data_blue, axis=1)
        img_data = np.append(img_data_blue[:, :spat_max], img_data_red[:, :spat_max], axis=0)
        vmax = 4
        vmin = -1
    else:
        raise ValueError('Expected channel value of 0, 1, or 2')
else:
    if channel == 0:
        img_data_blue = spec2DObj_blue.sciimg
        img_data_red = spec2DObj_red.sciimg
        vmax = 15
        vmin = -3
    elif channel == 1:
        gpm_blue = spec2DObj_blue.bpmmask == 0
        img_data_blue = (spec2DObj_blue.sciimg - spec2DObj_blue.skymodel) * gpm_blue
        vmax = 15
        vmin = -3
    elif channel == 2:
        gpm_blue = spec2DObj_blue.bpmmask == 0
        img_data = (spec2DObj_blue.sciimg - spec2DObj_blue.skymodel) * np.sqrt(spec2DObj_blue.ivarmodel) * gpm_blue
        vmax = 4
        vmin = -1
    else:
        raise ValueError('Expected channel value of 0, 1, or 2')


# Figure Plotting
img_hdu = fits.open(spec2d_file)
img_wave = img_hdu[(det - 1) * 11 + 8].data
img_wave[img_wave<10.0] = np.nan
#redshift = float(args.redshift)
#wave_lya = (1 + redshift) * lya
#wave_low = wave_lya-150
#wave_high = wave_lya+300
wave_low = J14_wave-300
wave_high = J14_wave+1000
spec_low = np.where(img_wave[:, wave_ind_blue] > wave_low)[0][0]
spec_high = np.where(img_wave[:, wave_ind_blue] < wave_high)[0][-1]

if np.sum(img_wave[:,wave_ind_blue]>wave_high) == 0:
    img_wave_red = img_hdu[(det + 4 - 1) * 11 + 8].data
    img_wave_red[img_wave_red < 10.0] = np.nan
    spec_high += np.where(img_wave_red[:, wave_ind_red] < wave_high)[0][-1]

blue_slit = sobjs[blue_exten - 1].SLITID
wave_ind = np.max([wave_ind_red,wave_ind_blue])

if channel == 1:
    # 2D Sensfunc

    sens = sensfunc.SensFunc.from_file('sens_2010sep24_d0924_0010.fits')

    spectrograph = load_spectrograph('keck_deimos')
    exptime = spectrograph.get_meta_value(files[1],'exptime')
    #exptime = 1600.0 #Obj 4219

    sens_factor = flux_calib.get_sensfunc_factor(spec2DObj.waveimg[:,wave_ind],
                                                 sens.wave.reshape(15353), sens.zeropoint.reshape(15353), exptime,
                                                 extrap_sens=True)

    sens_gpm = sens_factor < 100.0*np.nanmedian(sens_factor)
    sens_factor_masked = sens_factor*sens_gpm
    sens_factor_img = np.repeat(sens_factor_masked[:, np.newaxis], spec2DObj.waveimg[0].shape[0], #pseudo_dict['nspat']
                                            axis=1)

    img_data *= sens_factor_img
    #imgminsky_gpm = sens_gpm[:, np.newaxis] & pseudo_dict['inmask']
    vmax = 0.015
    vmin = -0.005


# 2D Spatial Range
spat_low = wave_ind - 35
spat_high = wave_ind + 35

slit_mask = img_hdu[(det - 1) * 11 + 10].data.spat_id == blue_slit
slit_low = img_hdu[(det - 1) * 11 + 10].data.left_init[slit_mask][0, 0]
slit_high = img_hdu[(det - 1) * 11 + 10].data.right_init[slit_mask][0, 0]

if (slit_low>spat_low)&(slit_high<spat_high):
    spat_low = slit_low
    spat_high = slit_high
elif slit_low>spat_low:
    pix_diff = slit_low-spat_low
    spat_low = slit_low
    spat_high += pix_diff
elif slit_high<spat_low:
    pix_diff = slit_high-spat_low
    spat_high = slit_high
    spat_low += pix_diff


# 1D Flux Range
wave_low_ind = np.where(np.abs(new_waves-wave_low)==np.min(np.abs(new_waves-wave_low)))[0][0]
wave_high_ind = np.where(np.abs(new_waves-wave_high)==np.min(np.abs(new_waves-wave_high)))[0][0]
flux_range = flux_corr[wave_low_ind:wave_high_ind+1]
flux_space = (flux_range.max()-flux_range.min())/20


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
fig, ax = plt.subplots(2, 1, figsize=(20, 16))

ax[0].imshow(img_data.transpose(), vmin=vmin, vmax=vmax, cmap='gray')
ax[0].set_ylim(spat_low, spat_high)
ax[0].set_xlim(spec_low, spec_high)
ax[0].axes.get_xaxis().set_visible(False)
ax[0].axes.get_yaxis().set_visible(False)

trans = ax[1].get_xaxis_transform()
ax[1].step(new_waves, flux_corr, 'k', linewidth=1, where='mid', label=r'\textbf{Observed Spectrum}')
ax[1].plot(new_waves, sig_corr, 'r:', linewidth=3, label=r'\textbf{Observed Uncertainty}')

ax[1].axvline(J14_wave, color='y', linestyle='--', alpha=0.5)
ax[1].text(J14_wave, .85, 'J1438', transform=trans, backgroundcolor='0.75')

ax[1].set_xlabel(r'\textbf{Wavelength (\AA)}', size=30)
ax[1].set_ylabel(r'$$\bf F_{\lambda} \quad (10^{-17} erg s^{-1} cm^{-2} \AA^{-1})$$', size=30)
ax[1].legend(prop={"size": 20})
#ax[1].set_ylim(-0.02, 0.11)
ax[1].set_ylim(flux_range.min()-flux_space,flux_range.max()+flux_space)
ax[1].set_xlim(wave_low, wave_high)
ax[1].xaxis.set_minor_locator(MultipleLocator(10))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.004))
ax[1].tick_params('both', length=20, width=2, which='major', labelsize=22)
ax[1].tick_params('both', length=10, width=1, which='minor')
#ax[0].set_title(r'\textbf{Obj 8986 Spectrum}', size=24)
plt.tight_layout(h_pad=0)
plt.subplots_adjust(hspace=-.42)
plt.savefig('spec_figure.png', bbox_inches='tight')
plt.show()
plt.close()

