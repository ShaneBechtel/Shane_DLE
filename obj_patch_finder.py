import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from glob import glob
from IPython import embed

path = '/home/sbechtel/Downloads/unzip_hold/tobias_code/J1438/deepCoAdd/'
patch_files = 'calexp-HSC-R2*.fits'

files = glob(path+patch_files)

coords = '14:38:41.170 +43:06:36.202' #OBJ 6364
#coords = '14:38:42.323 +43:08:57.939' #OBJ 7325

c = SkyCoord(coords, unit=(u.hourangle,u.deg))

ra_val = c.ra.value
dec_val = c.dec.value

min_val = np.inf
min_ind = None

good_files = []

for i,file in enumerate(files):

    patch_fits = fits.open(file)[1]

    w = WCS(patch_fits.header)

    ra_ind, dec_ind = w.world_to_pixel(c)

    ra_flag = (ra_ind>0) & (ra_ind<patch_fits.header['NAXIS1'])
    dec_flag = (dec_ind>0) & (dec_ind<patch_fits.header['NAXIS2'])

    if ra_flag & dec_flag:
        good_files.append(file)


print(good_files)

'''
for i,file in enumerate(files):

    patch_fits = fits.open(file)[1]

    ref_ra = patch_fits.header['CRVAL1']
    ref_dec = patch_fits.header['CRVAL2']

    ra_pix_step = patch_fits.header['CD1_1']
    dec_pix_step = patch_fits.header['CD2_2']

    ra_pix_off = patch_fits.header['CRPIX1'] - (size/2)
    dec_pix_off = patch_fits.header['CRPIX2'] - (size / 2)

    patch_ra = (ra_pix_off*ra_pix_step) + ref_ra
    patch_dec = (dec_pix_off*dec_pix_step) + ref_dec

    ra_off = ra_val-patch_ra
    dec_off = dec_val-patch_dec

    dist = np.sqrt(ra_off**2 + dec_off**2)

    if dist < min_val:

        min_ind = i
        min_val = dist
        
print(files[min_ind])

'''

