import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from IPython import embed

color_path = '/home/sbechtel/Downloads/unzip_hold/Data/QSOJ1630+0435.cat'

color_fits = fits.open(color_path)[1]

umag = color_fits.data['u'] # 'u' or 'ua'?
gmag = color_fits.data['g']
rmag = color_fits.data['r']
imag = color_fits.data['i']

u_r_mag = umag - rmag
g_r_mag = gmag - rmag
r_i_mag = rmag - imag

u_r_mag[u_r_mag>6.9] = 6.9
g_r_mag[g_r_mag>4.1] = 4.1
r_i_mag[r_i_mag<-0.9] = -0.9

u_r_cut = 4.8
g_r_cut = 1.0
r_i_cut = 1.0

red_r_i = np.linspace((g_r_cut-0.8)/1.5,r_i_cut,100)
red_g_r = np.linspace(g_r_cut,(1.5*r_i_cut+0.8),100)

plt.scatter(r_i_mag, g_r_mag, color='b', s=0.1)
plt.plot(red_r_i,red_g_r,'r')
plt.hlines(g_r_cut,-2,(g_r_cut-0.8)/1.5,'r')
plt.vlines(r_i_cut,(1.5*r_i_cut+0.8),5,'r')
plt.xlim(-0.95,2.1)
plt.ylim(-0.2,4.15)
plt.xlabel('LBC R - LBC I in mag')
plt.ylabel('LBC G - LBC R in mag')
plt.tick_params(axis='x', which='both', direction='in', top=True, bottom=True)
plt.tick_params(axis='y', which='both', direction='in', left=True, right=True)
plt.show()


plt.scatter(r_i_mag, u_r_mag, color='b', s=0.1)
plt.hlines(u_r_cut,-2,r_i_cut,'r')
plt.vlines(r_i_cut,u_r_cut,8,'r')
plt.xlim(-0.95,2.1)
plt.ylim(-0.2,6.95)
plt.xlabel('LBC R - LBC I in mag')
plt.ylabel('LBC U - LBC R in mag')
plt.tick_params(axis='x', which='both', direction='in', top=True, bottom=True)
plt.tick_params(axis='y', which='both', direction='in', left=True, right=True)
plt.show()