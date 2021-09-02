import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from IPython import embed

star_hdu = fits.open('star_spec.fits')
star_wave = star_hdu[1].data['OPT_WAVE']
star_flux = star_hdu[1].data['OPT_COUNTS']

tell_hdu = fits.open('star_spec_tellcorr.fits')
tell_wave = tell_hdu[1].data['wave']
tell_flux = tell_hdu[1].data['flux']

tell_corr = tell_flux/star_flux

np.savez('tell_corr',tell_wave=tell_wave,tell_corr=tell_corr)