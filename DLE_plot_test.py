import numpy as np
from pypeit import specobjs
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from linetools.spectra.xspectrum1d import XSpectrum1D
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file', help='Spectrum Fits File')
parser.add_argument('--redshift', help='Redshift of emission lines')
parser.add_argument('--exten', help='Exten value for target in Blue Detector ')
args = parser.parse_args()


blue_exten = int(args.exten)
#red_exten = int(args.red)


#fits_file = "../../DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Both/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits"
fits_file = args.file
sobjs = specobjs.SpecObjs.from_fitsfile(fits_file, chk_version=False)
blue_spec = sobjs[blue_exten-1].to_xspec1d(extraction='OPT', fluxed=False)

det = sobjs[blue_exten-1].DET
spat_pix = sobjs[blue_exten-1].SPAT_PIXPOS

for i, sobj in enumerate(sobjs):
    if (sobj.DET == det+4) & (sobj.SPAT_PIXPOS == spat_pix):
        red_ind = i
        break

red_spec = sobjs[red_ind].to_xspec1d(extraction='OPT', fluxed=False)
zero_skip = blue_spec.wavelength.value > 10

new_wavelength = np.concatenate((blue_spec.wavelength[zero_skip],red_spec.wavelength[red_spec.wavelength>blue_spec.wavelength[-1]]))
new_flux = np.concatenate((blue_spec.flux[zero_skip],red_spec.flux[red_spec.wavelength>blue_spec.wavelength[-1]]))
new_sig = np.concatenate((blue_spec.sig[zero_skip],red_spec.sig[red_spec.wavelength>blue_spec.wavelength[-1]]))
new_spec = XSpectrum1D.from_tuple((new_wavelength, new_flux, new_sig), verbose=False)

#new_spec.plot()


# Line List
redshift = float(args.redshift)
Lines = {"C III":[1175.71], "Si II":[1190,1260.42], "Lya":[1215.670],
         "N V":[1240.81], "O I":[1305.53], "C II":[1335.31]}


fig, ax = plt.subplots()

# the x coords of this transformation are data, and the
# y coord are axes
trans = ax.get_xaxis_transform()

ax.plot(new_wavelength,new_flux)

for tup in Lines.items():

    if np.shape(tup[1])[0]==1:
        z_wavelength = (1+redshift)*tup[1][0]
        ax.vlines(z_wavelength,new_flux.min(),new_flux.max(),'k',linestyles='--')
        plt.text(z_wavelength, .85, tup[0], transform=trans,backgroundcolor='0.75')
    else:
        for l in tup[1]:
            z_wavelength = (1 + redshift) * l
            ax.vlines(z_wavelength, new_flux.min(), new_flux.max(), 'k', linestyles='--')
            plt.text(z_wavelength, .85, tup[0], transform=trans, backgroundcolor='0.75')

plt.show()