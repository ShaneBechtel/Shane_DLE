import argparse

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, AutoLocator, MultipleLocator, FixedLocator
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
parser.add_argument('--redshift', default=0.0, help='Redshift of emission lines')
parser.add_argument('--exten', help='Exten value for target')
parser.add_argument('--width', help='Width of boxcar smoothing')
parser.add_argument('--flux', default=False, help='Show fluxed spectra?')
parser.add_argument('--channel', default=2, help='Which channel to include for 2D image in figure')
args = parser.parse_args()


# Example Call
# python DLE_lya_figure_mosaic.py --file_path /home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1438A/final_redux/setup_FWHM/Science_coadd/ --redshift 5.15 --exten 35 --width 5 --channel 1 --flux True

def ivarsmooth(flux, ivar, window):
    '''
    Boxcar smoothign of width window with ivar weights
    Args:
        flux:
        ivar:
        window:
    Returns:
    '''
    nflux = flux.shape[0]
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


exten = int(args.exten)
file_path=args.file_path
if file_path[-1] != '/':
    file_path += '/'
files = glob(file_path+'*')
files.sort()
spec1d_file = files[0]
spec2d_file = files[1]
#tell_file = files[2]

sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=True)

spec1d = sobjs[exten - 1]


obj_name = r'\textbf{' + spec1d.MASKDEF_OBJNAME.lstrip('0') + r'}'
obj_name = r"{}".format(obj_name)
new_waves = spec1d.OPT_WAVE

zero_skip = new_waves > 10
new_waves = new_waves[zero_skip]

flux = eval(args.flux)

if flux:
    new_flux = spec1d.OPT_FLAM[zero_skip]
    new_ivars = spec1d.OPT_FLAM_IVAR[zero_skip]
else:
    new_flux = spec1d.OPT_COUNTS[zero_skip]
    new_ivars = spec1d.OPT_COUNTS_IVAR[zero_skip]
new_masks = new_ivars > 0.0


width = int(args.width)
new_flux, new_ivars = ivarsmooth(new_flux, new_ivars, width)
new_sig = 1 / np.sqrt(new_ivars)


# Quasar Redshifts
lya = 1215.67
J14_wave = lya * (1 + 4.709)
J16_wave = lya * (1 + 3.8101)

'''
# Atmospheric Effects
tell_hdu = fits.open(tell_file)
tell_waves = tell_hdu[1].data['wave']
tell_spec = tell_hdu[1].data['telluric']
tell_corr = np.interp(new_waves, tell_waves, tell_spec, left=1, right=1)

flux_corr = new_flux / (tell_corr + (tell_corr == 0))
ivar_corr = (tell_corr > 0.0) * new_ivars * tell_corr * tell_corr
mask_corr = (tell_corr > 0.0) * new_masks
sig_corr = np.sqrt(utils.inverse(ivar_corr))
'''

flux_corr = new_flux
ivar_corr = new_ivars
mask_corr = new_masks
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
#comp_flux = (np.median(new_flux[new_flux!=0.0]) / np.nanmedian(comp_flux[comp_flux!=0.0])) * comp_flux






# 2D Image
det = sobjs[exten - 1].DET
det_num = int(det[-1])
spec2DObj = spec2dobj.Spec2DObj.from_file(spec2d_file, det, chk_version=False)
channel = int(args.channel)

wave_ind = int(np.round(sobjs[exten - 1].SPAT_PIXPOS))

# 2D Sensitivity test

if channel == 0:
    img_data = spec2DObj.sciimg
    vmax = 15
    vmin = -3
elif channel == 1:
    gpm = spec2DObj.bpmmask == 0
    img_data = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm
    vmax = 15
    vmin = -3
elif channel == 2:
    gpm = spec2DObj.bpmmask == 0
    img_data = (spec2DObj.sciimg - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) * gpm
    vmax = 4
    vmin = -1
else:
    raise ValueError('Expected channel value of 0, 1, or 2')

# Figure Plotting
img_hdu = fits.open(spec2d_file)
img_wave = img_hdu[(det_num - 1) * 12 + 8].data[:,wave_ind]
img_wave[img_wave<10.0] = np.nan
redshift = float(args.redshift)
wave_lya = (1 + redshift) * lya
wave_low = wave_lya-150
wave_high = wave_lya+300
wave_ind = int(np.round(sobjs[exten - 1].SPAT_PIXPOS))
spec_low = np.where(img_wave > wave_low)[0][0]
spec_high = np.where(img_wave< wave_high)[0][-1]

wave_low = img_wave[spec_low]
wave_high = img_wave[spec_high]

comp_low = np.where((1+redshift)*comp_waves >= wave_low)[0][0]
comp_high = np.where((1+redshift)*comp_waves <= wave_high)[0][-1]
comp_flux = (np.median(new_flux[spec_low:spec_high]) / np.nanmedian(comp_flux[comp_low:comp_high])) * comp_flux

blue_slit = sobjs[exten - 1].SLITID


# 2D Spatial Range
spat_low = wave_ind - 35
spat_high = wave_ind + 35

slit_mask = img_hdu[(det_num - 1) * 12 + 10].data.spat_id == blue_slit
slit_low = img_hdu[(det_num - 1) * 12 + 10].data.left_init[slit_mask][0, 0]
slit_high = img_hdu[(det_num - 1) * 12 + 10].data.right_init[slit_mask][0, 0]

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


if flux:
    if channel == 1:
        # 2D Sensfunc

        #sens = sensfunc.SensFunc.from_file('../DLE_auxillary/sens_2010sep24_d0924_0010.fits')
        sens = sensfunc.SensFunc.from_file('../DLE_auxillary/keck_deimos_600ZD_sensfunc.fits')

        spectrograph = load_spectrograph('keck_deimos')
        exptime = spectrograph.get_meta_value(files[1],'exptime')
        #exptime = 1600.0 #Obj 4219

        sens_factor = flux_calib.get_sensfunc_factor(spec2DObj.waveimg[:,wave_ind],
                                                     sens.wave.squeeze(), sens.zeropoint.squeeze(), exptime,
                                                     extrap_sens=True)
        sens_gpm = sens_factor < 100.0 * np.nanmedian(sens_factor)
        sens_factor_masked = sens_factor*sens_gpm
        sens_factor_img = np.repeat(sens_factor_masked[:, np.newaxis], spec2DObj.waveimg[0].shape[0], #pseudo_dict['nspat']
                                                axis=1)

        img_data *= sens_factor_img
        #imgminsky_gpm = sens_gpm[:, np.newaxis] & pseudo_dict['inmask']

        #2D Flux Range

        fwhm_low = wave_ind - 10
        fwhm_high = wave_ind + 10

        mad_std_low = utils.nan_mad_std(img_data[spec_low:spec_high,int(spat_low):fwhm_low])
        mad_std_high = utils.nan_mad_std(img_data[spec_low:spec_high,fwhm_high:int(spat_high)])
        mad_std = np.mean([mad_std_low,mad_std_high])

        vmax = 5*mad_std
        vmin = -2*mad_std


# 1D Flux Range
wave_low_ind = np.where(np.abs(new_waves-wave_low)==np.min(np.abs(new_waves-wave_low)))[0][0]
wave_high_ind = np.where(np.abs(new_waves-wave_high)==np.min(np.abs(new_waves-wave_high)))[0][0]
flux_range = flux_corr[wave_low_ind:wave_high_ind+1]


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
ax[1].plot((1 + redshift) * comp_waves, comp_flux, 'g', linewidth=2, alpha=0.5,
           label=r'\textbf{Jones et al. 2012 Composite}')

ax[1].vlines(wave_lya, new_flux.min(), new_flux.max(), 'b', linewidth=2, linestyles='--', alpha=0.5)
ax[1].text(wave_lya - 30, .9, r'$\bf Ly\alpha$', transform=trans, backgroundcolor='0.75', size=24)

ax[1].set_xlabel(r'\textbf{Wavelength (\AA)}', size=30)
ax[1].set_ylabel(r'$$\bf F_{\lambda} \quad (10^{-17} erg s^{-1} cm^{-2} \AA^{-1})$$', size=30)
ax[1].legend(prop={"size": 20})
mad_std_1D = utils.nan_mad_std(flux_range)
ax[1].set_ylim(flux_range.mean()-mad_std_1D*5,flux_range.mean()+mad_std_1D*15)
ax[1].set_xlim(wave_low, wave_high)

tick_label_size = 20
major_tick_width = 2
major_tick_length = 15
minor_tick_width = 2
minor_tick_length = 7

ax[1].tick_params(axis='x', which='both', direction='in', top=True, bottom=True)
ax[1].tick_params(axis='y', which='both', direction='in', left=True, right=True)
ax[1].tick_params(axis="x", which='major', labelsize=tick_label_size, length=major_tick_length, width=major_tick_width)
ax[1].tick_params(axis="x", which='minor', labelsize=tick_label_size, length=minor_tick_length, width=minor_tick_width)
ax[1].tick_params(axis="y", which='major', labelsize=tick_label_size, length=major_tick_length, width=major_tick_width)
ax[1].tick_params(axis="y", which='minor', labelsize=tick_label_size, length=minor_tick_length, width=minor_tick_width)
ax[1].yaxis.set_minor_locator(AutoMinorLocator())
ax[1].yaxis.set_major_locator(AutoLocator())
ax[1].xaxis.set_minor_locator(AutoMinorLocator())
ax[1].xaxis.set_major_locator(AutoLocator())

ax[0].set_title(r'\textbf{Obj 18799 Spectrum}', size=24)
#ax[0].set_title(r'\textbf{Obj \# Spectrum}', size=24)
plt.tight_layout(h_pad=0)
plt.subplots_adjust(hspace=-.42)
#plt.savefig('spec_figure.png', bbox_inches='tight')
plt.show()
plt.close()


