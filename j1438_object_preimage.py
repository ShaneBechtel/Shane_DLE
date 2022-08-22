import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Circle
from astropy.stats import sigma_clipped_stats
from matplotlib.offsetbox import AnchoredText

from IPython import embed


mosaic_path = '/home/sbechtel/Downloads/unzip_hold/tobias_code/J1438/deepCoAdd/'
bands = ['r','i','z']

#obj = '6364'
obj = '7325'

if obj == '6364':
    rband_name = 'calexp-HSC-R2-0-6.6.fits'
    iband_name = 'calexp-HSC-I2-0-6.6.fits'
    zband_name = 'calexp-HSC-Z-0-6.6.fits'
    #coords = '14:38:41.170 +43:06:36.202' #OBJ 6364
    coords = '14:38:41.2 +43:06:33.5' #OBJ 6364 Shifted manually to nearest object

if obj == '7325':
    rband_name = 'calexp-HSC-R2-0-6.5.fits'
    iband_name = 'calexp-HSC-I2-0-6.5.fits'
    zband_name = 'calexp-HSC-Z-0-6.5.fits'
    #coords = '14:38:42.323 +43:08:57.939' #OBJ 7325
    coords = '14:38:42.40 +43:08:55.8'  # OBJ 7325 Shifted manually to nearest object

band_paths = [rband_name,iband_name,zband_name]


c = SkyCoord(coords, unit=(u.hourangle,u.deg))

coords_ra = c.ra.degree
coords_dec = c.dec.degree

image_fits = fits.open(mosaic_path+band_paths[0])[1]

ra_ref_val = image_fits.header['CRVAL1']
ra_ref_pix = image_fits.header['CRPIX1']
ra_scale = image_fits.header['CD1_1']

dec_ref_val = image_fits.header['CRVAL2']
dec_ref_pix = image_fits.header['CRPIX2']
dec_scale = image_fits.header['CD2_2']

size = image_fits.header['NAXIS1']

ra_cen_val = (ra_ref_pix - (size/2))*ra_scale + ra_ref_val
dec_cen_val = (dec_ref_pix - (size/2))*dec_scale + dec_ref_val


ra_pixel_offset = (coords_ra - ra_cen_val)/ra_scale
ra_index = int(np.round((size/2) + ra_pixel_offset))

dec_pixel_offset = (coords_dec - dec_cen_val)/dec_scale
dec_index = int(np.round((size/2) + dec_pixel_offset))


image_width = 30 #pixels

mpl.rcParams['axes.linewidth'] = 5

fig, ax = plt.subplots(1, 3,sharex=True, sharey=True,figsize=(20,5))
fig.subplots_adjust(wspace=0)

for i in range(len(band_paths)):

    if i!=0:
        image_fits = fits.open(mosaic_path+band_paths[i])[1]

    image = image_fits.data[dec_index-image_width:dec_index+image_width,ra_index-image_width:ra_index+image_width]

    sig_clip_1D = sigma_clipped_stats(image)[2]

    ax[i].imshow(image, origin='lower', cmap='binary', vmin=image.mean()-5*sig_clip_1D,vmax=image.mean()+5*sig_clip_1D)
    ax[i].tick_params(axis=u'both', which=u'both', length=0)
    ax[i].xaxis.set_ticklabels([])
    ax[i].yaxis.set_ticklabels([])

    circ = Circle((image_width,image_width),0.4*image_width,color='red',linestyle='--',linewidth=2,fill=False)
    ax[i].add_patch(circ)

    at = AnchoredText(
        bands[i], prop=dict(size=30), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax[i].add_artist(at)

ax[0].set_ylabel('OBJ 8986',size=30)
plt.show()

embed()