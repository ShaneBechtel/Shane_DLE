
import sys
from IPython import embed
from math import *

# Example Call
# python transverse_dist_calc.py 3.8101 0.1479

if sys.argv[1] == '-h':
    print( '''Cosmology calculator obtained from Ned Wright (www.astro.ucla.edu/~wright)
input values = redshift, theta
ouput values = transverse_distance between QSO and b/g galaxy

Options:   -h for this message ''')
    sys.exit()

z = float(sys.argv[1])  # redshift, z_QSO
theta = float(sys.argv[2]) # Angular Seperation, degrees, use https://cads.iiap.res.in/tools/angularSeparation

H0 = float(70.0)  # Hubble constant, 70 These values from Worseck Paper (J1630)
WM = float(0.27)  # Omega(matter)
WV = float(0.73)  # Omega(vacuum) or lambda

# initialize constants

WR = 0.  # Omega(radiation)
WK = 0.  # Omega curvaturve = 1-Omega(total)
c = 299792.458  # velocity of light in km/sec
Tyr = 977.8  # coefficent for converting 1/H into Gyr
DTT = 0.5  # time from z to now in units of 1/H0
DTT_Gyr = 0.0  # value of DTT in Gyr
age = 0.5  # age of Universe in units of 1/H0
age_Gyr = 0.0  # value of age in Gyr
zage = 0.1  # age of Universe at redshift z in units of 1/H0
zage_Gyr = 0.0  # value of zage in Gyr
DCMR = 0.0  # comoving radial distance in units of c/H0
DCMR_Mpc = 0.0
DCMR_Gyr = 0.0
DA = 0.0  # angular size distance
DA_Mpc = 0.0
DA_Gyr = 0.0
kpc_DA = 0.0
DL = 0.0  # luminosity distance
DL_Mpc = 0.0
DL_Gyr = 0.0  # DL in units of billions of light years
V_Gpc = 0.0
a = 1.0  # 1/(1+z), the scale factor of the Universe
az = 0.5  # 1/(1+z(object))

h = H0 / 100.
WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species, T0 = 2.72528
WK = 1 - WM - WR - WV
az = 1.0 / (1 + 1.0 * z)
age = 0.
n = 1000  # number of points in integrals
for i in range(n):
    a = az * (i + 0.5) / n
    adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
    age = age + 1. / adot

zage = az * age / n
zage_Gyr = (Tyr / H0) * zage
DTT = 0.0
DCMR = 0.0

# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
for i in range(n):
    a = az + (1 - az) * (i + 0.5) / n
    adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
    DTT = DTT + 1. / adot
    DCMR = DCMR + 1. / (a * adot)

DTT = (1. - az) * DTT / n
DCMR = (1. - az) * DCMR / n
age = DTT + zage
age_Gyr = age * (Tyr / H0)
DTT_Gyr = (Tyr / H0) * DTT
DCMR_Gyr = (Tyr / H0) * DCMR
DCMR_Mpc = (c / H0) * DCMR

# tangential comoving distance

ratio = 1.00
x = sqrt(abs(WK)) * DCMR
if x > 0.1:
    if WK > 0:
        ratio = 0.5 * (exp(x) - exp(-x)) / x
    else:
        ratio = sin(x) / x
else:
    y = x * x
    if WK < 0: y = -y
    ratio = 1. + y / 6. + y * y / 120.
DCMT = ratio * DCMR
DA = az * DCMT
DA_Mpc = (c / H0) * DA
kpc_DA = DA_Mpc / 206.264806
DA_Gyr = (Tyr / H0) * DA

theta_rad = theta * pi / 180

transverse_dist = theta_rad * DA_Mpc #Mpc

print("Transverse distance of : "+str(transverse_dist)+" Mpc")