import numpy as np
from pypeit import specobj
from pypeit import specobjs
from pypeit.core.coadd import multi_combspec
from IPython import embed

# blue_exten = 10
# red_exten = 25
# fits_file = '../../DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Star/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits'

blue_exten = 3
red_exten = 38
fits_file = '../../DEIMOS_Light_Echo/Targets/J1630A/det_all/setup_Star/Science_coadd/spec1d_DE.20190705.25097-DE.20190705.36380-J1630A.fits'

# blue_exten = 7
# red_exten = 39
# fits_file = '../../DEIMOS_Light_Echo/Targets/J1630B/det_all/setup_Star/Science_coadd/spec1d_DE.20190706.22720-DE.20190706.32153-J1630B.fits'

sobjs = specobjs.SpecObjs.from_fitsfile(fits_file, chk_version=True)

blue_spec = sobjs[blue_exten - 1]
red_spec = sobjs[red_exten - 1]

if len(blue_spec.OPT_WAVE) < len(red_spec.OPT_WAVE):

    red_wave = red_spec.OPT_WAVE
    red_flux = red_spec.OPT_COUNTS
    red_ivar = red_spec.OPT_COUNTS_IVAR

    blue_wave = np.zeros_like(red_wave)
    blue_flux = np.zeros_like(red_wave)
    blue_ivar = np.zeros_like(red_wave)

    diff = len(red_spec.OPT_WAVE) - len(blue_spec.OPT_WAVE)

    blue_wave[diff:] = blue_spec.OPT_WAVE
    blue_flux[diff:] = blue_spec.OPT_COUNTS
    blue_ivar[diff:] = blue_spec.OPT_COUNTS_IVAR


else:

    blue_wave = blue_spec.OPT_WAVE
    blue_flux = blue_spec.OPT_COUNTS
    blue_ivar = blue_spec.OPT_COUNTS_IVAR

    red_wave = np.zeros_like(blue_wave)
    red_flux = np.zeros_like(blue_wave)
    red_ivar = np.zeros_like(blue_wave)

    diff = len(blue_spec.OPT_WAVE) - len(red_spec.OPT_WAVE)

    red_wave[diff:] = red_spec.OPT_WAVE
    red_flux[diff:] = red_spec.OPT_COUNTS
    red_ivar[diff:] = red_spec.OPT_COUNTS_IVAR

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

# Continium Continuity

new_waves, new_flux, new_ivars, new_masks = multi_combspec(waves, fluxes, ivars, masks, wave_grid_max=wgmax,
                                                           wave_grid_min=wgmin)

zero_skip = new_waves > 10
new_waves = new_waves[zero_skip]
new_flux = new_flux[zero_skip]
new_ivars = new_ivars[zero_skip]
new_masks = new_masks[zero_skip]
'''
#Masking J1438A
mask_lowind = np.where(np.abs(new_waves-6700)==np.min(np.abs(new_waves-6700)))[0][0]
mask_highind = np.where(np.abs(new_waves-7120)==np.min(np.abs(new_waves-7120)))[0][0]
new_masks[mask_lowind:mask_highind] = 0
new_ivars[mask_lowind:mask_highind] = 0
'''

'''
# Masking J1630
mask_lowind = np.where(np.abs(new_waves - 6650) == np.min(np.abs(new_waves - 6650)))[0][0]
mask_highind = np.where(np.abs(new_waves - 6850) == np.min(np.abs(new_waves - 6850)))[0][0]
new_masks[mask_lowind:mask_highind] = 0
new_ivars[mask_lowind:mask_highind] = 0
mask_lowind = np.where(np.abs(new_waves - 6500) == np.min(np.abs(new_waves - 6500)))[0][0]
mask_highind = np.where(np.abs(new_waves - 6600) == np.min(np.abs(new_waves - 6600)))[0][0]
new_masks[mask_lowind:mask_highind] = 0
new_ivars[mask_lowind:mask_highind] = 0
'''

# Edge of Sensitivity Function
new_ivars[new_waves < 5500] = 0
new_flux[new_waves < 5500] = 0

sobj = specobj.SpecObj.from_arrays('MultiSlit', new_waves, new_flux, new_ivars)

new_sobjs = specobjs.SpecObjs()
new_sobjs.add_sobj(sobj)
embed()

'''
# J1438A
subheader = {'PYP_SPEC':'keck_deimos','PYPELINE':'MultiSlit','INSTRUME':'DEIMOS  ',
             'airmass':1.08991231,'exptime':1600.0,'ra':'219.610125','dec':'43.23066666666667',
             'target':'PS1_00013','DISPNAME':'600ZD','decker':'J1438A','binning':'1,1',
             'mjd':58639.34915,'FILENAME':'J1438A_star_spec.fits'}
'''

# J1630A
subheader = {'PYP_SPEC': 'keck_deimos', 'PYPELINE': 'MultiSlit', 'INSTRUME': 'DEIMOS  ',
             'airmass': 1.07627586, 'exptime': 1800.0, 'ra': '247.7417083333333', 'dec': '4.622166666666667',
             'target': 'PS1_00211', 'DISPNAME': '600ZD', 'decker': 'J1630A', 'binning': '1,1',
             'mjd': 58669.290408, 'FILENAME': 'J1630A_star_spec.fits'}

'''
# J1630B
subheader = {'PYP_SPEC': 'keck_deimos', 'PYPELINE': 'MultiSlit', 'INSTRUME': 'DEIMOS  ',
             'airmass': 1.08212524, 'exptime': 1800.0, 'ra': '247.79883333333328', 'dec': '4.614694444444444',
             'target': 'PS1_00570', 'DISPNAME': '600ZD', 'decker': 'J1630B', 'binning': '1,1',
             'mjd': 58670.262897, 'FILENAME': 'J1630B_star_spec.fits'}
'''

new_sobjs.write_to_fits(subheader, 'J1630A_star_spec.fits', overwrite=True)
