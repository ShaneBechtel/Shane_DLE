import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from IPython import embed


#Currently just working for J1438
color_path = '/home/sbechtel/Downloads/unzip_hold/tobias_code/HSC_forced_J1438_short.cat'

color_fits = fits.open(color_path)[1]

rmag = color_fits.data['HSC-R2_MAG_APER']
imag = color_fits.data['HSC-I2_MAG_APER']
zmag = color_fits.data['HSC-Z_MAG_APER']

i_z_mag = imag - zmag
r_i_mag = rmag - imag

i_z_mag[i_z_mag<-0.6] = -0.6
r_i_mag[r_i_mag>2.55] = 2.55

i_z_cut = 0.7
r_i_cut = 1.2

red_i_z = np.linspace((r_i_cut-1.0)/1.5,i_z_cut,100)
red_r_i = np.linspace(r_i_cut,(1.5*i_z_cut+1.0),100)

plt.scatter(i_z_mag, r_i_mag, color='b', s=0.1)
plt.plot(red_i_z,red_r_i,'r')
plt.hlines(r_i_cut,-1,(r_i_cut-1.0)/1.5,'r')
plt.vlines(i_z_cut,(1.5*i_z_cut+1.0),3,'r')
plt.xlim(-0.65,2.1)
plt.ylim(-0.1,2.6)
plt.xlabel('HSC I - HSC Z in mag')
plt.ylabel('HSC R - HSC I in mag')
plt.tick_params(axis='x', which='both', direction='in', top=True, bottom=True)
plt.tick_params(axis='y', which='both', direction='in', left=True, right=True)
plt.show()