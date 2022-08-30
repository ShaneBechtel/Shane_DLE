import numpy as np
import astropy.units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from IPython import embed

T = Table.read('/home/sbechtel/Downloads/unzip_hold/Data/QSOJ1630+0435.cat')

# E_BV    = Targets['E_BV'][i]

U_SN = T['u'] / T['ue']
G_SN = T['g'] / T['ge']
R_SN = T['r'] / T['re']
I_SN = T['i'] / T['ie']

#T['COL_IZ'] = T['HSC-I2_MAG_APER'] - T['HSC-Z_MAG_APER'] - E_BV * (1.698 - 1.263)
#T['COL_RI'] = T['HSC-R2_MAG_APER'] - T['HSC-I2_MAG_APER'] - E_BV * (2.285 - 1.698)
GR = T['g'] - T['r']
RI = T['r'] - T['i']

Sel = np.ones(len(T), dtype='bool')
#Sel = Sel * (T['ra'] > 19.5) * (T['ra'] < 26.5)
Sel = Sel * (R_SN>5) * (I_SN>5)

Sel1 = Sel
Sel1 = Sel1 * (GR > 1.0)
Sel1 = Sel1 * (RI < 1.0)
Sel1 = Sel1 * (GR > 1.5 * RI + 0.8)
Sel1 = Sel1 * ( T['r']>20 ) * ( T['r']<25.3)

Sel2 = Sel1 * (U_SN < 3)


plt.plot(np.fmax(-0.9, RI[Sel]), np.fmin(4.2, GR[Sel]), marker='.', lw=0, markersize=3, alpha=0.5)
plt.plot([-99, 0.2 / 1.5, 1.0, 1.0], [1.0, 1.0, 2.3, 99], lw=1.5, ls='-', color='red')
plt.xlabel("$\mathrm{LBC\:R - LBC\:I \; in \; mag}$")
plt.ylabel("$\mathrm{LBC\:G - LBC\:R \; in \; mag}$")
plt.xlim(-0.95, 2.1)
plt.ylim(-0.2, 4.25)
plt.tight_layout()

plt.plot(np.fmax(-0.9, RI[Sel1]), np.fmin(4.2, GR[Sel1]), marker='.', lw=0, markersize=3,
         color='darkorange', alpha=1.0)
plt.plot(np.fmax(-0.9, RI[Sel2]), np.fmin(4.2, GR[Sel2]), marker='o', lw=0, markersize=3,
         color='darkorange', alpha=1.0)

#plt.savefig(Target['label'].split('+')[0] + '/' + 'ColorSelection_' + Targets['label'][i] + '_SEx.pdf')
plt.show()