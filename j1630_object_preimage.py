import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from matplotlib.patches import Circle
from astropy.stats import sigma_clipped_stats
from matplotlib.offsetbox import AnchoredText

from IPython import embed


mosaic_path = '/home/sbechtel/Downloads/unzip_hold/Data/J1630/final/'

uband_name = 'mosaic_QSOJ1630+0435_bUs_med.fits'
gband_name = 'mosaic_QSOJ1630+0435_bg_med.fits'
rband_name = 'mosaic_QSOJ1630+0435_rr_med.fits'
iband_name = 'mosaic_QSOJ1630+0435_ri_med.fits'

bands = ['u', 'g', 'r', 'i']
band_paths = [uband_name, gband_name, rband_name, iband_name]


#J1630
coords = '16:31:23.852 +04:30:21.493' # 8986
#coords = '16:31:18.383 +04:38:29.937' # 22571
#coords = '16:31:08.308 +04:36:21.113' # 18799
#coords = '16:31:03.981  +04:30:33.559' # 9321
#coords = '16:31:00.714  +04:33:35.532' #PS1_65, bright star below and to left of J1630

c = SkyCoord(coords, unit=(u.hourangle,u.deg))

coords_ra = c.ra.degree
coords_dec = c.dec.degree

image_fits = fits.open(mosaic_path+band_paths[0])[0]

w = WCS(image_fits.header)

dec_index, ra_index = w.world_to_pixel(c)
dec_index = int(dec_index)
ra_index = int(ra_index)

'''
ra_center_val = image_fits.header['CRVAL1']
ra_center_pix = image_fits.header['CRPIX1']
ra_scale = image_fits.header['CD1_1']

dec_center_val = image_fits.header['CRVAL2']
dec_center_pix = image_fits.header['CRPIX2']
dec_scale = image_fits.header['CD2_2']

ra_pixel_offset = (coords_ra - ra_center_val)/ra_scale
ra_index = int(np.round(ra_center_pix + ra_pixel_offset))

dec_pixel_offset = (coords_dec - dec_center_val)/dec_scale
dec_index = int(np.round(dec_center_pix + dec_pixel_offset))
'''
pix_scale = image_fits.header['CD1_1']
pix_arc = -3600*pix_scale

image_width = int(np.round(6/pix_arc)) # Create images with width 6"

mpl.rcParams['axes.linewidth'] = 5

fig, ax = plt.subplots(1, 4,sharex=True, sharey=True,figsize=(20,5))
fig.subplots_adjust(wspace=0)

for i in range(len(band_paths)):

    if i!=0:
        image_fits = fits.open(mosaic_path+band_paths[i])[0]

    image = image_fits.data[dec_index-image_width:dec_index+image_width,ra_index-image_width:ra_index+image_width]

    sig_clip_1D = sigma_clipped_stats(image)[2]

    """ # Used for OBJ 18799
    if i!= 3:
        ax[i].imshow(image,origin='lower',cmap='binary',vmin=image.mean()-5*sig_clip_1D,vmax=image.mean()+5*sig_clip_1D)
    else:
        ax[i].imshow(image, origin='lower', cmap='binary', vmin=image.mean() - 15 * sig_clip_1D,
                     vmax=image.mean() + 15 * sig_clip_1D)
    """

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
