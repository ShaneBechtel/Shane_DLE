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
from astropy.stats import sigma_clipped_stats
from pypeit import utils
from scipy.interpolate import interp1d
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed

parser = argparse.ArgumentParser()
parser.add_argument('--file_path', help='Coadd directory path')
parser.add_argument('--exten', help='Exten value for target')
parser.add_argument('--width', help='Width of boxcar smoothing')
parser.add_argument('--flux', default=False, help='Show fluxed spectra?')
parser.add_argument('--channel', default=2, help='Which channel to include for 2D image in figure')
args = parser.parse_args()


# Example Call
# python DLE_mosaic_plot.py --file_path /home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1630B/branch_mosaic/setup_FWHM/Science_coadd/ --exten 59 --width 5 --channel 1 --flux True

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


exten = int(args.exten)

file_path=args.file_path
if file_path[-1] != '/':
    file_path += '/'
files = glob(file_path+'*')
files.sort()

spec1d_file = files[0]
spec2d_file = files[1]
#tell_file = files[2] #TODO Fix J1630B Telluric Issues
sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=True) #TODO Get OBJ NAME & OTHER DETAILS

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

# Atmospheric Effects #TODO Don't hard code in this telluric file.
#tell_hdu = fits.open(tell_file)
tell_hdu = fits.open('../DLE_auxillary/star_spec_tellcorr.fits')
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
sig_corr = new_sig
'''


# 2D Image
det = sobjs[exten - 1].DET
det_num = int(det[-1])
spec2DObj = spec2dobj.Spec2DObj.from_file(spec2d_file, det, chk_version=True)
channel = int(args.channel)

spat_max = spec2DObj.sciimg.shape[1]

wave_ind = int(np.round(sobjs[exten - 1].SPAT_PIXPOS))
waves = spec2DObj.waveimg[:,wave_ind]

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
img_wave = img_hdu[(det_num - 1) * 12 + 8].data[:,wave_ind] #TODO Make sure this indexing is still correct for Mosaic
img_wave[img_wave<10.0] = np.nan
#redshift = float(args.redshift)
#wave_lya = (1 + redshift) * lya
#wave_low = wave_lya-150
#wave_high = wave_lya+300
#wave_low = J14_wave-800
#wave_high = J14_wave+3000

# J1438 Range
#wave_low = 6000
#wave_high = 10000

# J1630 Range
wave_low = 5500
wave_high = 9000
spec_low = np.where(img_wave > wave_low)[0][0]
spec_high = np.where(img_wave < wave_high)[0][-1]

#TODO Try and see what is necessary for Mosaic here
'''
if np.sum(img_wave>wave_high) == 0:
    img_wave_red = img_hdu[(det_num + 4 - 1) * 12 + 8].data
    img_wave_red[img_wave_red < 10.0] = np.nan
    spec_high += np.where(img_wave_red[:, wave_ind_red] < wave_high)[0][-1]
'''

slit_id = sobjs[exten - 1].SLITID

# 2D Spatial Range
spat_low = wave_ind - 70
spat_high = wave_ind + 70

slit_mask = img_hdu[(det_num - 1) * 12 + 10].data.spat_id == slit_id
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

        sens = sensfunc.SensFunc.from_file('../DLE_auxillary/keck_deimos_600ZD_sensfunc.fits')

        spectrograph = load_spectrograph('keck_deimos')
        exptime = spectrograph.get_meta_value(files[1],'exptime')
        #exptime = 1600.0 #Obj 4219

        sens_factor = flux_calib.get_sensfunc_factor(waves,sens.wave.squeeze(), sens.zeropoint.squeeze(), exptime,
                                                     extrap_sens=True)

        sens_gpm = sens_factor < 100.0*np.nanmedian(sens_factor)
        sens_factor_masked = sens_factor*sens_gpm
        sens_factor_img = np.repeat(sens_factor_masked[:, np.newaxis], waves.shape[0], #pseudo_dict['nspat']
                                                axis=1)

        img_data *= sens_factor_img[:, :spat_max]
        #imgminsky_gpm = sens_gpm[:, np.newaxis] & pseudo_dict['inmask']

        #2D Flux Range

        fwhm_low = wave_ind - 10
        fwhm_high = wave_ind + 10

        mad_std_low = utils.nan_mad_std(img_data[spec_low:spec_high,int(spat_low):fwhm_low])
        mad_std_high = utils.nan_mad_std(img_data[spec_low:spec_high,fwhm_high:int(spat_high)])
        mad_std = np.mean([mad_std_low,mad_std_high])

        vmax = 5*mad_std
        vmin = -2*mad_std

# TODO USE ASTROPY sigma_clip_stats
# 1D Flux Range
wave_low_ind = np.where(np.abs(new_waves-wave_low)==np.min(np.abs(new_waves-wave_low)))[0][0]
wave_high_ind = np.where(np.abs(new_waves-wave_high)==np.min(np.abs(new_waves-wave_high)))[0][0]
flux_range = flux_corr[wave_low_ind:wave_high_ind+1]


def forceAspect(ax,aspect):
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


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
forceAspect(ax[0],aspect=1.0)

trans = ax[1].get_xaxis_transform()
ax[1].step(new_waves, flux_corr, 'k', linewidth=1, where='mid', label=r'\textbf{Observed Spectrum}')
ax[1].plot(new_waves, sig_corr, 'r:', linewidth=3, label=r'\textbf{Observed Uncertainty}')

#ax[1].axvline(J14_wave, color='y', linestyle='--', alpha=0.5)
#ax[1].text(J14_wave, .85, 'J1438', transform=trans, backgroundcolor='0.75')
ax[1].axvline(J16_wave, color='y', linestyle='--', alpha=0.5)
ax[1].text(J16_wave, .85, 'J1630', transform=trans, backgroundcolor='0.75')

ax[1].axvline(J16_wave, color='y', linestyle='--', alpha=0.5)
ax[1].text(J16_wave, .85, 'J1630', transform=trans, backgroundcolor='0.75')

ax[1].set_xlabel(r'\textbf{Wavelength (\AA)}', size=30)
ax[1].set_ylabel(r'$$\bf F_{\lambda} \quad (10^{-17} erg s^{-1} cm^{-2} \AA^{-1})$$', size=30)
ax[1].legend(prop={"size": 20})
sig_clip_1D = sigma_clipped_stats(flux_range)[2]
ax[1].set_ylim(flux_range.mean()-sig_clip_1D*5,flux_range.mean()+sig_clip_1D*8)
ax[1].set_xlim(wave_low, wave_high)
ax[1].xaxis.set_minor_locator(MultipleLocator(100))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.25))
#ax[1].yaxis.set_minor_locator(MultipleLocator(0.005))
ax[1].tick_params('both', length=20, width=2, which='major', labelsize=22)
ax[1].tick_params('both', length=10, width=1, which='minor')
ax[0].set_title(r'\textbf{(OBJ }' + obj_name + r'\textbf{)}', size=24)
plt.tight_layout(h_pad=0)
plt.subplots_adjust(hspace=-.442)
#plt.savefig('test_figure.png', bbox_inches='tight')
plt.show()
plt.close()


