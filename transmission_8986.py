import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from glob import glob
from pypeit import specobjs
from astropy.io import fits
from scipy.optimize import curve_fit
from matplotlib.ticker import AutoMinorLocator, AutoLocator, MultipleLocator, FixedLocator
from transmission_significance import trans_sig_comp
from IPython import embed


parser = argparse.ArgumentParser()
parser.add_argument('--file_path', help='Coadd directory path')
parser.add_argument('--exten', help='Exten value for target')
parser.add_argument('--redshift', default=0.0, help='Redshift of emission lines')
parser.add_argument('--qso', help='Which Quasar field is object in? (J14 or J16)')
args = parser.parse_args()


# Example Call
# python transmission_8986.py --file_path /home/sbechtel/Documents/DEIMOS_Light_Echo/Targets/J1630B/final_redux/setup_FWHM/Science_coadd/ --exten 59 --redshift 4.166 --qso J1630


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
qso = args.qso
if qso.lower() == 'j1630':
    qso_wave = lya * (1 + 3.8101)
    qso_sig = 400 # km/s
elif qso.lower() == 'j1438':
    qso_wave = lya * (1 + 4.709)
else:
    raise Exception('Please enter a valid QSO option: J1438 or J1630')

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


# Currently Hard coded for Obj 8986

z = float(args.redshift)

lya_wave = lya * (1+z)

cont_spec_low = lya_wave + 30.0 #Angstroms
waves_holder = new_waves[mask_corr]
cont_spec_ind = np.where(waves_holder >= cont_spec_low)[0][0]

ind_gap_8986 = int(300/0.6) #TODO Fix detector flux discontinuity

red_side_ind_low = np.where(waves_holder>=6800.)[0][0]
red_side_ind_high = np.where(waves_holder<=7200)[0][-1]

#TODO Manually fix Flux for 8986

def cont_powerlaw(x, a, b, c):
    return a * (x/c)**b


y = flux_corr[mask_corr][cont_spec_ind:cont_spec_ind+ind_gap_8986]

flux_stitch = np.median(y)/np.median(flux_corr[mask_corr][red_side_ind_low:red_side_ind_high])


x = waves_holder[cont_spec_ind:]


det_ind = np.where(x>=6700)[0][0]

y_corr = flux_corr[mask_corr][cont_spec_ind:]
y_corr[det_ind:] *= flux_stitch

cont_fit = curve_fit(cont_powerlaw,x,y_corr,maxfev=100000,bounds=([0,-np.inf,1300*(1+z)-1],[np.inf,0,1300*(1+z)+1]))


a,b,c = cont_fit[0]

good_waves = new_waves[mask_corr]
good_flux = flux_corr[mask_corr]
good_sig = sig_corr[mask_corr]

cont_flux = cont_powerlaw(good_waves,a,b,c)

trans_flux = good_flux/cont_flux
trans_sig = sig_corr[mask_corr]/cont_flux

vel_range = 299792.458 * (good_waves-qso_wave)/qso_wave #km/s
vel_low = -3000
vel_high = 3000

H0 = 67.7 #TODO DON'T USE LOCAL H0, use H(z)!!!!!!!!!
H_z = np.sqrt(H0**2 * (0.3*(4.83101)**3 + 0.7))
dist_range = vel_range / H_z #Mpc
comov_range = dist_range * (1+z) # Comoving distance cMpc

gamma_uvb = 1.03214425e-12 # Gamma_uvb from Nyx simulations
gamma_qso_1cmpc = 2.62575329e-09 # Gamma_qso at 1 cMpc for wqso
gamma_qso = 7.46286498e-12 # Gamma_qso from Joe Code

transverse_dist = 3.8996 # pMpc

wqso = gamma_qso_1cmpc/gamma_uvb


#gamma_vals = gamma_uvb * (1.0 + wqso/(dist_range**2 + transverse_dist**2)) # Uses Proper distance
gamma_vals = gamma_uvb * (1.0 + wqso/(comov_range**2 + transverse_dist**2)) # Uses Comoving distance

#lower_trans_ind = np.where(dist_range>-4)[0][0]
#upper_trans_ind = np.where(dist_range<4)[0][-1]
#lower_trans_ind = np.where(vel_range>-qso_sig)[0][0]
#upper_trans_ind = np.where(vel_range<qso_sig)[0][-1]

# trans_mask = (dist_range>=-4) & (dist_range<=4)
trans_mask = (vel_range>=-qso_sig) & (vel_range<=qso_sig)

#trans_sig_comp(dist_range[lower_trans_ind:upper_trans_ind],trans_flux[lower_trans_ind:upper_trans_ind],trans_sig[lower_trans_ind:upper_trans_ind],gamma_vals[lower_trans_ind:upper_trans_ind])
trans_sig_comp(dist_range[trans_mask],trans_flux[trans_mask],trans_sig[trans_mask],gamma_vals[trans_mask])


qso_uncertainty = (vel_range > -qso_sig) & (vel_range < qso_sig)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'


fig= plt.figure(figsize=(20,12))
ax = fig.add_subplot(111,label='1')
ax2 = fig.add_subplot(111,label='2',frame_on=False)

trans = ax.get_xaxis_transform()
ax.step(vel_range, trans_flux, 'k', linewidth=1, where='mid', label=r'\textbf{Observed Spectrum}')
ax.fill_between(vel_range[qso_uncertainty], 0, 1, transform=trans, color='gray', alpha=0.5)
ax.axvline(0, color='y', linestyle='--', alpha=0.5)
ax.text(0, .9, qso.upper(), transform=trans, backgroundcolor='0.75')
ax.set_xlabel(r'\textbf{Velocity (km s$^{-1}$)}', size=30)
ax.set_ylabel(r'\textbf{Transmission}', size=30)
ax.set_ylim(-0.1,1.2)
ax.set_xlim(vel_low, vel_high)

tick_label_size = 20
major_tick_width = 2
major_tick_length = 15
minor_tick_width = 2
minor_tick_length = 7

ax.tick_params(axis='x', which='both', direction='in', top=True, bottom=True)
ax.tick_params(axis='y', which='both', direction='in', left=True, right=True)
ax.tick_params(axis="x", which='major', labelsize=tick_label_size, length=major_tick_length, width=major_tick_width)
ax.tick_params(axis="x", which='minor', labelsize=tick_label_size, length=minor_tick_length, width=minor_tick_width)
ax.tick_params(axis="y", which='major', labelsize=tick_label_size, length=major_tick_length, width=major_tick_width)
ax.tick_params(axis="y", which='minor', labelsize=tick_label_size, length=minor_tick_length, width=minor_tick_width)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(AutoLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_major_locator(AutoLocator())


ax2.plot(dist_range, trans_sig, 'r:', linewidth=3, label=r'\textbf{Observed Uncertainty}')
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.set_xlim(vel_low/H0, vel_high/H0)
ax2.set_xlabel(r'\textbf{Distance (pMpc)}', size=30,labelpad=15)
ax2.set_ylim(-0.1,1.2)
ax2.tick_params('both', length=20, width=2, which='major', labelsize=22)
ax2.tick_params('both', length=10, width=1, which='minor')
ax2.xaxis.set_label_position('top')
ax2.yaxis.set_label_position('right')
ax2.get_yaxis().set_visible(False)


plt.savefig('trans_test.png')
plt.show()
plt.close()

