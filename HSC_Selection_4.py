

import numpy as np
import astropy.units as u 
from astropy.table import Table

import matplotlib.pyplot as plt
from astropy.time import Time
from IPython import embed

from PyFindingchart import Findingchart

def MakeDsimulatorTable( Targets, GuideStars, OutFile, CatRA='HSC-I2_ALPHAWIN_J2000', CatDec='HSC-I2_DELTAWIN_J2000', CatFilter='HSC-I2', CatMag='MAG_AUTO', GuideMag='i_mean_psf' ):
	from astropy.coordinates import SkyCoord
	f=open( OutFile, 'w' )
	Targets['ICRS']	= SkyCoord( Targets[ CatRA ].quantity, Targets[ CatDec ] )
	STR	= '# OBJNAME' + '  ' + '     RA     ' + '    DEC     ' + '  ' + ' EQX  ' + '  ' + ' MAG ' + '  ' +  'band ' + '  ' + 'PCODE' + '  '
	f.write( STR + '\n' )
	for i,target in enumerate(Targets):
		STR = ''
		STR = STR + '{:09.0f}'.format( Targets['NUMBER'][i] )	+ '  '
		STR = STR + Targets['ICRS'][i].to_string('hmsdms', precision=3).split(' ')[0].replace('h',':').replace('m',':').replace('s','')		+ '  '
		STR = STR + Targets['ICRS'][i].to_string('hmsdms', precision=3).split(' ')[1].replace('d',':').replace('m',':').replace('s','')		+ '  '
		STR = STR + '2000.0'	+ '  '
		#STR = STR + '{:05.2f}'.format( Targets[ CatFilter + '_' + CatMag ][i] )	+ '  '
		STR = STR + '  {:1s}  '.format( CatFilter[4:5] )	+ '  '
		#STR = STR + '{:5.0f}'.format( Targets['PCode'][i] )	+ '  '
		f.write( STR + '\n' )
		print(STR)
	GuideStars['ICRS']	= SkyCoord( GuideStars[ 'RA' ].quantity, GuideStars[ 'Dec' ] )
	for i,target in enumerate(GuideStars):
		STR = ''
		STR = STR + '{:9s}'.format( 'PS1_Guide' )	+ '  '
		STR = STR + GuideStars['ICRS'][i].to_string('hmsdms', precision=3).split(' ')[0].replace('h',':').replace('m',':').replace('s','')		+ '  '
		STR = STR + GuideStars['ICRS'][i].to_string('hmsdms', precision=3).split(' ')[1].replace('d',':').replace('m',':').replace('s','')		+ '  '
		STR = STR + '2000.0'	+ '  '
		STR = STR + '{:05.2f}'.format( GuideStars[ GuideMag ][i] )	+ '  '
		STR = STR + '  {:1s}  '.format( GuideMag[0:1] )	+ '  '
		STR = STR + '{:5.0f}'.format( -2 )	+ '  '
		f.write( STR + '\n' )
		print(STR)
	f.close()

R = 1

Targets=Table.read('/home/sbechtel/Downloads/unzip_hold/tobias_code/HSC_Hennawi_2019A_Targets_Final.fits')

for i, Target in enumerate( Targets ):
	if Targets['label'][i]!='J1438+4314' or Targets['Drop'][i] != 'r' or i!=7: continue
	#if Targets['label'][i]!='J1347+4956': continue
	E_BV    = Targets['E_BV'][i]


	COL	=       ['id', 'HSC-I2_classification_extendedness','HSC-I2_RA','HSC-I2_Dec', 'HSC-I2_FLAGS' ]
	COL	= COL + [ 'HSC-R2_FLUX_APER', 'HSC-R2_flux_kron', 'HSC-R2_flux_psf', 'HSC-R2_cmodel_flux' ]
	COL	= COL + [ 'HSC-I2_FLUX_APER', 'HSC-I2_flux_kron', 'HSC-I2_flux_psf', 'HSC-I2_cmodel_flux' ]
	COL	= COL + [ 'HSC-Z_FLUX_APER',  'HSC-Z_flux_kron',  'HSC-Z_flux_psf',  'HSC-Z_cmodel_flux' ]
	COL	= COL + [ 'HSC-R2_SN', 'HSC-I2_SN', 'HSC-Z_SN' ]

	T	= Table.read( Targets['label'][i].split('+')[0] + '/CatalogShort_' + Targets['label'][i].split('+')[0]  + '_4.fits' )

	#T['HSC-R2_SN']	= T['HSC-R2_cmodel_flux'] / T['HSC-R2_cmodel_flux_err']
	#T['HSC-I2_SN']	= T['HSC-I2_cmodel_flux'] / T['HSC-I2_cmodel_flux_err']
	#T['HSC-Z_SN']	= T['HSC-Z_cmodel_flux']  / T['HSC-Z_cmodel_flux_err']

	T['HSC-R2_SN']	= T['HSC-R2_FLUX_APER'] / T['HSC-R2_FLUXERR_APER']
	T['HSC-I2_SN']	= T['HSC-I2_FLUX_APER'] / T['HSC-I2_FLUXERR_APER']
	T['HSC-Z_SN']	= T[ 'HSC-Z_FLUX_APER'] / T[ 'HSC-Z_FLUXERR_APER']

	#T['HSC-R2_MAG_cmodel']	= -2.5*np.log10( np.fmax(T['HSC-R2_cmodel_flux']       , 0) ) + 2.5*np.log10( 63095734448.0194 )
	#T['HSC-I2_MAG_cmodel']	= -2.5*np.log10( np.fmax(T['HSC-I2_cmodel_flux']       , 0) ) + 2.5*np.log10( 63095734448.0194 )
	#T['HSC-Z_MAG_cmodel']	= -2.5*np.log10( np.fmax(T['HSC-Z_cmodel_flux']        , 0) ) + 2.5*np.log10( 63095734448.0194 )
	T['HSC-R2_MAG_ap']	= -2.5*np.log10( np.fmax(T['HSC-R2_FLUX_APER'], 0) ) + 2.5*np.log10( 63095734448.0194 )
	T['HSC-I2_MAG_ap']	= -2.5*np.log10( np.fmax(T['HSC-I2_FLUX_APER'], 0) ) + 2.5*np.log10( 63095734448.0194 )
	T['HSC-Z_MAG_ap' ]	= -2.5*np.log10( np.fmax(T[ 'HSC-Z_FLUX_APER'], 0) ) + 2.5*np.log10( 63095734448.0194 )

	#COL_IZ		= T['HSC-I2_MAG_cmodel'] - T['HSC-Z_MAG_cmodel']		- E_BV*(1.698-1.263)
	#COL_RI		= T['HSC-R2_MAG_cmodel'] - T['HSC-I2_MAG_cmodel']		- E_BV*(2.285-1.698)
	T['COL_IZ']	= T['HSC-I2_MAG_ap'] - T[ 'HSC-Z_MAG_ap']		- E_BV*(1.698-1.263)
	T['COL_RI']	= T['HSC-R2_MAG_ap'] - T['HSC-I2_MAG_ap']		- E_BV*(2.285-1.698)

	Sel	= np.ones( len(T), dtype='bool' )
	Sel	= Sel *	( np.all(T['HSC-I2_FLAGS']==False) == 0 )
	#Sel	= Sel *	( np.all(T['HSC-Z_FLAGS' ]==False) == 0 )
	Sel	= Sel *	( T['HSC-I2_SN'] > 5 )
	Sel	= Sel *	( T['HSC-Z_SN']  > 5 )

	#Using Selection method from Other Script
	Sel = np.ones(len(T), dtype='bool')
	Sel = Sel * (T['HSC-I2_MAG_APER'] > 20) * (T['HSC-I2_MAG_APER'] < 26.0)

	plt.close()
	plt.plot( np.fmax(-0.6, T['COL_IZ'][Sel]), np.fmin(2.6, T['COL_RI'][Sel]), marker='.', lw=0, markersize=3, alpha=0.5 )
	plt.plot( [-99, 0.2, 0.7, 0.7], [1.2,1.2,2.1,99], lw=1.5, color='red' )

	Sel1	= Sel
	Sel1	= Sel1 * ( T['COL_RI'] > 1.2 )
	Sel1	= Sel1 * ( T['COL_IZ'] < 0.7 )
	Sel1	= Sel1 * ( T['COL_RI'] > ( 1.5 * T['COL_IZ'] + 1.0 ) )
	#Sel1	= Sel1 * ( T['HSC-I2_MAG_ap']>21 ) * ( T['HSC-I2_MAG_ap']<24.7 )
	Sel1 = Sel1 * (T['HSC-I2_MAG_ap'] > 20) * (T['HSC-I2_MAG_ap'] < 25.3) # From Other Script
	#Sel1	= Sel1 * ( np.all(T['HSC-I2_FLAGS'] ) == 0 )

	plt.plot( np.fmax(-0.6, T['COL_IZ'][Sel1]), np.fmin(2.6, T['COL_RI'][Sel1]), marker='o', lw=0, markersize=3, color='darkorange', alpha=1.0 )

	plt.xlabel("$\mathrm{HSC\:I - HSC\:Z \; in \; mag}$")
	plt.ylabel("$\mathrm{HSC\:R - HSC\:I \; in \; mag}$")
	plt.xlim( -0.65, 2.1 )
	plt.ylim( -0.2, 2.7 )
	plt.tight_layout()
	plt.savefig( 'ColorSelection_' + Targets['label'][i] + '_4.pdf' )
	plt.close()

	#print(T[Sel1][COL])
	#input()

	from MakePostageStamps import MakePostageStamp as MakePostageStamp
	CatRA='HSC-R2_ALPHAWIN_J2000'; CatDec='HSC-R2_DELTAWIN_J2000'
	#CatRA='HSC-I2_RA'; CatDec='HSC-I2_Dec']

	ImgPath = "/home/sbechtel/Downloads/unzip_hold/tobias_code/J1438/deepCoAdd/"

	for j, target in enumerate( T ):
		if (Sel*Sel1)[j] == False : continue
		#print('\n', T['id'][j])
		#MakePostageStamp(	T[CatRA].quantity[j], T[CatDec].quantity[j], T['Patch'][j], str(T['id'][j]), Filters=['HSC-R2', 'HSC-I2', 'HSC-Z' ], Mags=T['HSC-R2_MAG_cmodel', 'HSC-I2_MAG_cmodel', 'HSC-Z_MAG_cmodel'][j],
		#			ImgPath = Targets['label'][i].split('+')[0] + '/deepCoAdd/', OutPath = Targets['label'][i].split('+')[0] + '/Charts4/', Size=12*u.arcsec, Radius=3*u.arcsec )
		MakePostageStamp(	T[CatRA].quantity[j], T[CatDec].quantity[j], T['Patch'][j], str(T['NUMBER'][j]), Filters=['HSC-R2', 'HSC-I2', 'HSC-Z' ], Mags=T['HSC-R2_MAG_ap', 'HSC-I2_MAG_ap', 'HSC-Z_MAG_ap'][j],
					ImgPath = ImgPath, OutPath = './Charts4/', Size=12*u.arcsec, Radius=1.5*u.arcsec )


	#T['PCode'] = 1000 - T['HSC-I2_MAG_cmodel'].quantity.to(u.mag).value * 10
	
	Figs, figs, Catalog	= Findingchart(	Targets['RA'].quantity[i], Targets['Dec'].quantity[i],
					z=Targets['Ref_Z'][i], M1450=Targets['M1450'].quantity[i],
					PA= Targets['PA'].quantity[i], OffsetSky=Targets['OffsetSky'].quantity[i], OffsetFP=Targets['OffsetFP'].quantity[i],
					Date=Time('2019-06-05'), Observatory='Keck', TVG_lim=17.0*u.mag, SDSS=False, Path=None, ImgPath=ImgPath, ImgSTR = 'calexp-HSC-I2-0-*.fits',
					ExtCat=T[Sel*Sel1], ExtCatRA='HSC-I2_RA', ExtCatDec='HSC-I2_Dec', Radius=4*u.arcsec, ExtCatLabel='$\mathrm{HSC \; r-drop}$' )

	#input()
	MakeDsimulatorTable( T[Sel], Catalog, Targets['label'][i].split('+')[0] + '_Targets.dat', CatRA='HSC-R2_ALPHAWIN_J2000', CatDec='HSC-R2_DELTAWIN_J2000', CatFilter='HSC-I2', CatMag='MAG_cmodel', GuideMag='i_mean_psf' )
	for j in range(3):
		PA 	= Targets['PA'].quantity[i][j] + 90*u.deg
		Offset	= -1 * np.array([ -1/np.cos(Targets['Dec'].quantity[i]), 1 ]) * np.dot( np.array( [ [ np.cos( PA ), -np.sin( PA ) ], [ np.sin( PA ), np.cos( PA ) ] ] ), Targets['OffsetFP'].quantity[i][j] )
		print(( Targets['RA'].quantity[i] + Offset[0] ).to(u.deg).value/360*24, ( Targets['Dec'].quantity[i] + Offset[1] ).to(u.deg).value, Targets['PA'].quantity[i][j].to(u.deg).value)
	#input()
	plt.savefig( 'DropoutSelection_' + Targets['label'][i] + '_4.pdf', dpi=800 )






















