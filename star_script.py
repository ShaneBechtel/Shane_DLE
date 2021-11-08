import numpy as np
from pypeit import specobj
from pypeit import specobjs
from pypeit.core.coadd import multi_combspec
from IPython import embed


#blue_exten = 10
#red_exten = 25
#fits_file = '../../DEIMOS_Light_Echo/Targets/J1438A/det_all/setup_Star/Science_coadd/spec1d_DE.20190605.30172-DE.20190605.35227-J1438A.fits'

blue_exten = 41
red_exten = 117
fits_file = '../../DEIMOS_Light_Echo/Targets/J1630B/det_all/setup_Both/Science_coadd/spec1d_DE.20190706.22720-DE.20190706.32153-J1630B.fits'

#blue_exten =
#red_exten =
#fits_file = ''

sobjs = specobjs.SpecObjs.from_fitsfile(fits_file, chk_version=False)
blue_spec = sobjs[blue_exten-1].to_xspec1d(extraction='OPT', fluxed=True)
red_spec = sobjs[red_exten-1].to_xspec1d(extraction='OPT', fluxed=True)

if len(blue_spec.wavelength)>len(red_spec.wavelength):
    diff = len(blue_spec.wavelength)-len(red_spec.wavelength)

    blue_wave = blue_spec.wavelength[diff:]
    blue_flux = blue_spec.flux[diff:]
    blue_ivar = blue_spec.ivar[diff:]

    red_wave = red_spec.wavelength
    red_flux = red_spec.flux
    red_ivar = red_spec.ivar
else:
    diff = len(red_spec.wavelength)-len(blue_spec.wavelength)
    length = len(red_spec.wavelength)

    blue_wave = blue_spec.wavelength
    blue_flux = blue_spec.flux
    blue_ivar = blue_spec.ivar

    red_wave = red_spec.wavelength[:length-diff]
    red_flux = red_spec.flux[:length-diff]
    red_ivar = red_spec.ivar[:length-diff]

waves = np.zeros((len(blue_wave),2))
fluxes = np.zeros((len(blue_wave),2))
ivars = np.zeros((len(blue_wave),2))

waves[:,0] = blue_wave
waves[:,1] = red_wave
fluxes[:,0] = blue_flux
fluxes[:,1] = red_flux
ivars[:,0] = blue_ivar
ivars[:,1] = red_ivar
masks = np.ones_like(waves,dtype=bool)

wgmax = np.max(red_wave.value)
wgmin = np.min(blue_wave[blue_wave.value>10].value)

#Continium Continuity

'''overlap_top = blue_spec.wavelength[-1].value
overlap_mask = red_spec.wavelength[red_spec.wavelength.value>10].value < overlap_top
overlap_num = int(np.sum(overlap_mask))
blue_sum = np.sum(blue_spec.flux.value[-1*overlap_num:])
red_sum = np.sum(red_spec.flux[red_spec.wavelength.value>10].value[overlap_mask])
ratio = blue_sum/red_sum
fluxes[:,1] *= ratio'''

new_waves, new_flux, new_ivars, new_masks = multi_combspec(waves,fluxes,ivars,masks,wave_grid_max=wgmax,wave_grid_min=wgmin)

zero_skip = new_waves> 10
new_waves = new_waves[zero_skip]
new_flux = new_flux[zero_skip]
new_ivars = new_ivars[zero_skip]
new_masks = new_masks[zero_skip]
#Masking strange dip feature in stellar spectrum
mask_lowind = np.where(np.abs(new_waves-6700)==np.min(np.abs(new_waves-6700)))[0][0]
mask_highind = np.where(np.abs(new_waves-7120)==np.min(np.abs(new_waves-7120)))[0][0]
new_masks[mask_lowind:mask_highind] = 0
new_ivars[mask_lowind:mask_highind] = 0

# Edge of Sensitivity Function
new_ivars[new_waves<5500] = 0
new_flux[new_waves<5500] = 0

sobj = specobj.SpecObj.from_arrays('MultiSlit',new_waves,new_flux,new_ivars)

new_sobjs = specobjs.SpecObjs()
new_sobjs.add_sobj(sobj)
embed()
subheader = {'PYP_SPEC':'keck_deimos','PYPELINE':'MultiSlit','INSTRUME':'DEIMOS  ',
             'airmass':1.08991231,'exptime':1600.0,'ra':'219.610125','dec':'43.23066666666667',
             'target':'PS1_00013','DISPNAME':'600ZD','decker':'J1438A','binning':'1,1',
             'mjd':58639.34915,'FILENAME':'star_spec.fits'}

new_sobjs.write_to_fits(subheader,'J16B_star_spec.fits',overwrite=True)