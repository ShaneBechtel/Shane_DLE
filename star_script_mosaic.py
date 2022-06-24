import numpy as np
from pypeit import specobj
from pypeit import specobjs
from pypeit.core.coadd import multi_combspec
from IPython import embed

#TODO Alter Code for Mosaic Inputs

#exten = 24
#fits_file = '../../DEIMOS_Light_Echo/Targets/J1438A/final_redux/setup_FWHM/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits'

#exten = 5
#fits_file = '../../DEIMOS_Light_Echo/Targets/J1630A/final_redux/setup_FWHM/Science_coadd/spec1d_DE.20190705.25097-DE.20190705.36380-J1630A.fits'

exten = 30
fits_file = '../../DEIMOS_Light_Echo/Targets/J1630B/final_redux/setup_FWHM/Science_coadd/spec1d_DE.20190706.22720-DE.20190706.32153-J1630B.fits'

sobjs = specobjs.SpecObjs.from_fitsfile(fits_file, chk_version=True)

star_spec = sobjs[exten - 1]

waves = star_spec.OPT_WAVE
flux = star_spec.OPT_COUNTS
ivars = star_spec.OPT_COUNTS_IVAR
masks = ivars > 0.0

zero_skip = waves > 10
new_waves = waves[zero_skip]
new_flux = flux[zero_skip]
new_ivars = ivars[zero_skip]
new_masks = masks[zero_skip]
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
'''# Masking J1630A
mask_lowind = np.where(np.abs(new_waves - 6958) == np.min(np.abs(new_waves - 6958)))[0][0]
mask_highind = np.where(np.abs(new_waves - 7058) == np.min(np.abs(new_waves - 7058)))[0][0]
new_masks[mask_lowind:mask_highind] = 0
new_ivars[mask_lowind:mask_highind] = 0
mask_lowind = np.where(np.abs(new_waves - 6515) == np.min(np.abs(new_waves - 6515)))[0][0]
mask_highind = np.where(np.abs(new_waves - 6615) == np.min(np.abs(new_waves - 6615)))[0][0]
new_masks[mask_lowind:mask_highind] = 0
new_ivars[mask_lowind:mask_highind] = 0
'''
# Edge of Sensitivity Function
# new_ivars[new_waves < 5500] = 0

sobj = specobj.SpecObj.from_arrays('MultiSlit', new_waves, new_flux, new_ivars)

new_sobjs = specobjs.SpecObjs()
new_sobjs.add_sobj(sobj)
embed()

'''
# J1438A
subheader = {'PYP_SPEC':'keck_deimos','PYPELINE':'MultiSlit','INSTRUME':'keck_deimos',
             'airmass':1.08991231,'exptime':1600.0,'ra':'219.610125','dec':'43.23066666666667',
             'target':'PS1_00013','DISPNAME':'600ZD','decker':'J1438A','binning':'1,1',
             'mjd':58639.34915,'FILENAME':'J1438A_star_spec.fits'}
'''
'''
# J1630A
subheader = {'PYP_SPEC': 'keck_deimos', 'PYPELINE': 'MultiSlit', 'INSTRUME': 'keck_deimos',
             'airmass': 1.07627586, 'exptime': 1800.0, 'ra': '247.7417083333333', 'dec': '4.622166666666667',
             'target': 'PS1_00211', 'DISPNAME': '600ZD', 'decker': 'J1630A', 'binning': '1,1',
             'mjd': 58669.290408, 'FILENAME': 'J1630A_star_spec.fits'}

'''
# J1630B
subheader = {'PYP_SPEC': 'keck_deimos', 'PYPELINE': 'MultiSlit', 'INSTRUME': 'keck_deimos',
             'airmass': 1.08212524, 'exptime': 1800.0, 'ra': '247.79883333333328', 'dec': '4.614694444444444',
             'target': 'PS1_00570', 'DISPNAME': '600ZD', 'decker': 'J1630B', 'binning': '1,1',
             'mjd': 58670.262897, 'FILENAME': 'J1630B_star_spec.fits'}


new_sobjs.write_to_fits(subheader, 'J1630B_star_spec.fits', overwrite=True)
