import numpy as np
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from MakePostageStamps import MakePostageStampLBC
from PyFindingchart import Findingchart
from IPython import embed

path = "/home/sbechtel/Downloads/unzip_hold/Data/J1630/final/"

RA = 247.7347 * u.deg
Dec = 4.5998 * u.deg



T = Table.read('/home/sbechtel/Downloads/unzip_hold/Data/QSOJ1630+0435.cat')

# E_BV    = Targets['E_BV'][i]

U_SN = T['u'] / T['ue']
G_SN = T['g'] / T['ge']
R_SN = T['r'] / T['re']
I_SN = T['i'] / T['ie']

GR = T['g'] - T['r']
RI = T['r'] - T['i']

Sel = np.ones(len(T), dtype='bool')
Sel = Sel * (T['ra'] > 19.5) * (T['ra'] < 26.5)
Sel = Sel * (R_SN>5) * (I_SN>5)

Sel1 = Sel
Sel1 = Sel1 * (GR > 1.0)
Sel1 = Sel1 * (RI < 1.0)
Sel1 = Sel1 * (GR > 1.5 * RI + 0.8)
Sel1 = Sel1 * ( T['r']>20 ) * ( T['r']<25.3)

Sel2 = Sel1 * (U_SN < 3)



#MakePostageStampLBC("QSOJ1630+0435", RA, Dec, "J1630", Filters=['r'], ImgPath=path, OutPath=path)

Targets=Table.read('/home/sbechtel/Downloads/unzip_hold/tobias_code/HSC_Hennawi_2019A_Targets_Final.fits')

Figs, figs, Catalog = Findingchart(RA, Dec, z=3.8101, M1450=-28.37*u.mag,
                                   Date=Time('2019-06-05'), Observatory='LBT', SDSS=False,
                                   Path=None, ImgPath=path, TVG_lim=17.0*u.mag,
                                   ExtCat = T[Sel2],
                                   ImgSTR='mosaic_QSOJ1630+0435_rr_med.fits',
                                   Radius=4 * u.arcsec, ExtCatLabel='$\mathrm{LBC \; g-drop}$')

