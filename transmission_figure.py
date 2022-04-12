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
from scipy.optimize import curve_fit
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed


parser = argparse.ArgumentParser()
parser.add_argument('--file_path', help='Coadd directory path')
parser.add_argument('--exten', help='Exten value for target')
parser.add_argument('--redshift', default=0.0, help='Redshift of emission lines')
args = parser.parse_args()


# Example Call
# python transmission_figure.py --file_path /home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1630B/branch_mosaic/setup_FWHM/Science_coadd/ --exten 59 --redshift 4.166


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

new_waves = spec1d.OPT_WAVE

zero_skip = new_waves > 10
new_waves = new_waves[zero_skip]

new_flux = spec1d.OPT_FLAM[zero_skip]
new_ivars = spec1d.OPT_FLAM_IVAR[zero_skip]
new_masks = new_ivars > 0.0
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

'''
flux_corr = new_flux / (tell_corr + (tell_corr == 0))
ivar_corr = (tell_corr > 0.0) * new_ivars * tell_corr * tell_corr
mask_corr = (tell_corr > 0.0) * new_masks
sig_corr = np.sqrt(utils.inverse(ivar_corr))
'''

flux_corr = new_flux
ivar_corr = new_ivars
mask_corr = new_masks
sig_corr = new_sig



#TODO Fit Continuum to Power Law

#TODO Encountering Issue of trying to fit detector flux gap. Need to solve this first. Try hard coding for 8986 for now?


z = float(args.redshift)

lya_wave = lya * (1+z)

cont_spec_low = lya_wave + 30.0 #Angstroms
waves_holder = new_waves[mask_corr]
cont_spec_ind = np.where(waves_holder >= cont_spec_low)[0][0]

ind_gap_8986 = int(300/0.6)

def cont_powerlaw(x, a, b, c):
    return a * (x/c)**b

x = waves_holder[cont_spec_ind:cont_spec_ind+ind_gap_8986]
y = flux_corr[mask_corr][cont_spec_ind:cont_spec_ind+ind_gap_8986]

cont_fit = curve_fit(cont_powerlaw,x,y,maxfev=100000,bounds=([0,-np.inf,500],[np.inf,0,2000]))


a,b,c = cont_fit[0]

#TODO Normalize Spectrum by Power Law Fit

good_waves = new_waves[mask_corr]
good_flux = flux_corr[mask_corr]
good_sig = sig_corr[mask_corr]

cont_flux = cont_powerlaw(good_waves,a,b,c)

trans_flux = good_flux/cont_flux
trans_flux[trans_flux>1.0] = 1.0
trans_flux[trans_flux<0.0] = 0.0

trans_sig = sig_corr[mask_corr]/cont_flux

# J1630 Range #TODO Change to be small region centered on Quasar Lya Wavelength (want about 100 angstroms of coverage)
wave_low = J16_wave - 60
wave_high = J16_wave + 30
spec_low = np.where(good_waves > wave_low)[0][0]
spec_high = np.where(good_waves < wave_high)[0][-1]


#TODO Convert Wavelengths to Velocity with v=0 being Quasar Lya

vel_range = 299792.458 * (good_waves-J16_wave)/J16_wave #km/s
vel_low = -3000
vel_high = 1500


plt.plot(vel_range,trans_flux)
plt.plot(vel_range,trans_sig,'r:')
plt.xlim(vel_low,vel_high)
plt.ylim(-0.1,1.2)
plt.close()



#TODO Take Velocities and create distance analogue measurement



dist_range = vel_range / 67.7 #Mpc


plt.plot(dist_range,trans_flux)
plt.plot(dist_range,trans_sig,'r:')
plt.xlim(-44,22)
plt.ylim(-0.1,1.2)
plt.close()

#TODO Convert Plotting Process to Single Axis


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
fig, ax = plt.subplots(figsize=(20,12))

trans = ax.get_xaxis_transform()
ax.step(vel_range, trans_flux, 'k', linewidth=1, where='mid', label=r'\textbf{Observed Spectrum}')
ax.plot(vel_range, trans_sig, 'r:', linewidth=3, label=r'\textbf{Observed Uncertainty}')

#plt.axvline(J14_wave, color='y', linestyle='--', alpha=0.5)
#plt.text(J14_wave, .85, 'J1438', transform=trans, backgroundcolor='0.75')
ax.vlines(0,-0.1,1.2, color='y', linestyle='--', alpha=0.5)
ax.text(0, .9, 'J1630', transform=trans, backgroundcolor='0.75')

ax.set_xlabel(r'\textbf{Velocity (km $$s^{-1}$$)}', size=30)
ax.set_ylabel(r'$$Ly\alpha$$ Transmission', size=30)
ax.legend(prop={"size": 20})
ax.set_ylim(-0.1,1.2)
ax.set_xlim(vel_low, vel_high)
#ax[1].xaxis.set_minor_locator(MultipleLocator(100))
#ax[1].yaxis.set_minor_locator(MultipleLocator(0.25))
ax.tick_params('both', length=20, width=2, which='major', labelsize=22)
ax.tick_params('both', length=10, width=1, which='minor')
ax.set_title(r'\textbf{(OBJ NAME)}', size=24)
#plt.savefig('test_figure.png', bbox_inches='tight')
plt.show()
plt.close()

