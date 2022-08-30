

import numpy as np
import astropy.units as u 
from astropy.table import Table

import matplotlib.pyplot as plt
from astropy.time import Time
from IPython import embed

from PyFindingchart import Findingchart

if 0:
	Sel	= np.ones( len(T), dtype='bool' )
	Sel 	= Sel * ( T['HSC-R2_FLUX_APER'][:,R] / T['HSC-R2_FLUXERR_APER'][:,R] > 5.0 )
	Sel 	= Sel * ( T['HSC-I2_FLUX_APER'][:,R] / T['HSC-I2_FLUXERR_APER'][:,R] > 5.0 )

	Sel	= Sel * ( np.all(T['HSC-R2_FLAGS']==False) )
	Sel	= Sel * ( -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] ) + 2.5*np.log10( 63095734448.0194 ) < 24.7 )

	COL_RI 	= -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] )
	COL_GR	= -2.5*np.log10( T['HSC-G_FLUX_APER' ][:,R] ) - -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] )

	plt.close()
	plt.plot( COL_RI[Sel], COL_GR[Sel], marker='.', lw=0, markersize=2, alpha=0.2 )
	plt.plot( [-99, 0.2/1.5, 1, 1], [1,1,2.3,99], lw=1.5, color='red' )
	plt.xlim( -0.6, 2.1 )
	plt.ylim( -0.2, 2.7 )

	Sel = Sel * ( COL_GR > 1.0 )
	Sel = Sel * ( COL_RI < 1.0 )
	Sel = Sel * ( COL_GR > 1.5*COL_RI + 0.8 )

	plt.plot( COL_RI[Sel], COL_GR[Sel], marker='.', lw=0, markersize=4, alpha=1 )

if 0:
	Sel	= np.ones( len(T), dtype='bool' )
	Sel 	= Sel * ( T['HSC-I2_FLUX_APER'][:,R] / T['HSC-I2_FLUXERR_APER'][:,R] > 5.0 )
	#Sel 	= Sel * ( T['HSC-Z_FLUX_APER'][:,R] / T['HSC-Z_FLUXERR_APER'][:,R] > 5.0 )
	#
	Sel	= Sel * ( np.all(T['HSC-I2_FLAGS']==False) )
	Sel	= Sel * ( -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] ) + 2.5*np.log10( 63095734448.0194 ) < 24.7 )
	#
	COL_IZ 	= -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-Z_FLUX_APER' ][:,R] )
	COL_RI	= -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] )
	#
	plt.close()
	plt.plot( COL_IZ[Sel], COL_RI[Sel], marker='.', lw=0, markersize=3, alpha=0.3, label='$\mathrm{Aperture\;1.5"}$' )
	plt.plot( [-99, 0.2, 0.7, 0.7], [1.2,1.2,2.1,99], lw=1.5, color='red' )
	plt.xlim( -0.6, 2.1 )
	plt.ylim( -0.2, 2.7 )

	Sel = Sel * ( COL_RI > 1.2 )
	Sel = Sel * ( COL_IZ < 0.7 )
	Sel = Sel * ( COL_RI > 1.5*COL_IZ + 1.0 )

	#plt.plot( COL_RI[Sel], COL_GR[Sel], marker='.', lw=0, markersize=4, alpha=1 )


	#plt.close()
	#plt.plot( T['HSC-R2_RA'][Sel], T['HSC-R2_Dec'][Sel], marker='.', lw=0 )

if 0:
	print("\nI-Band SN > 5:")
	print("Aperture", np.sum( T['HSC-I2_FLUX_APER'][:,R]	/ T['HSC-I2_FLUXERR_APER'][:,R] > 5.0 ))
	print("Kron    ", np.sum( T['HSC-I2_flux_kron'] 		/ T['HSC-I2_flux_kron_err'] > 5.0 ))
	print("PSF     ", np.sum( T['HSC-I2_flux_psf'] 			/ T['HSC-I2_flux_psf_err'] > 5.0 ))
	print("cModel  ", np.sum( T['HSC-I2_cmodel_flux'] 		/ T['HSC-I2_cmodel_flux_err'] > 5.0 ))
	print("\nObjects < 24.7 mag:")
	print("Aperture", np.sum( -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] )	+ 2.5*np.log10( 63095734448.0194 ) < 24.7 ))
	print("Kron    ", np.sum( -2.5*np.log10( T['HSC-I2_flux_kron'] ) 		+ 2.5*np.log10( 63095734448.0194 ) < 24.7 ))
	print("PSF     ", np.sum( -2.5*np.log10( T['HSC-I2_flux_psf'] ) 		+ 2.5*np.log10( 63095734448.0194 ) < 24.7 ))
	print("cModel  ", np.sum( -2.5*np.log10( T['HSC-I2_cmodel_flux'] ) 		+ 2.5*np.log10( 63095734448.0194 ) < 24.7 ))
	print("\nObjects < 24.7 mag, SN > 5 and FLAGS=False:")
	print("Aperture", np.sum( (-2.5*np.log10(T['HSC-I2_FLUX_APER'][:,R])+2.5*np.log10(63095734448.0194)<24.7) * (T['HSC-I2_FLUX_APER'][:,R]/T['HSC-I2_FLUXERR_APER'][:,R]>5.0) * np.all(T['HSC-I2_FLAGS']==False) ))
	print("Kron    ", np.sum( (-2.5*np.log10(T['HSC-I2_flux_kron']) 	+2.5*np.log10(63095734448.0194)<24.7) * (T['HSC-I2_flux_kron'    ]     /T['HSC-I2_flux_kron_err']         >5.0) * np.all(T['HSC-I2_FLAGS']==False) ))
	print("PSF     ", np.sum( (-2.5*np.log10(T['HSC-I2_flux_psf']) 		+2.5*np.log10(63095734448.0194)<24.7) * (T['HSC-I2_flux_psf'     ]     /T['HSC-I2_flux_psf_err']          >5.0) * np.all(T['HSC-I2_FLAGS']==False) ))
	print("cModel  ", np.sum( (-2.5*np.log10(T['HSC-I2_cmodel_flux']) 	+2.5*np.log10(63095734448.0194)<24.7) * (T['HSC-I2_cmodel_flux'  ]     /T['HSC-I2_cmodel_flux_err']       >5.0) * np.all(T['HSC-I2_FLAGS']==False) ))

if 0:
	E_BV	= Targets['E_BV'][i]
	plt.close()
	Sel	= np.ones( len(T), dtype='bool' )
	Sel 	= Sel * ( T['HSC-I2_flux_kron'][:] / T['HSC-I2_flux_kron_err'] > 5.0 )
	#Sel 	= Sel * ( T['HSC-Z_flux_kron'][:]  / T['HSC-Z_flux_kron_err'][:]  > 5.0 )
	#
	Sel	= Sel * ( np.all(T['HSC-I2_FLAGS']==False) )
	Sel	= Sel * ( -2.5*np.log10( T['HSC-I2_flux_kron'][:] ) + 2.5*np.log10( 63095734448.0194 ) < 24.7 )
	Sel	= Sel * ( T['HSC-I2_classification_extendedness']==0 )
	#
	COL_IZ 	= -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-Z_FLUX_APER' ][:,R] ) - E_BV*(1.698-1.263)
	COL_RI	= -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] ) - E_BV*(2.285-1.698)
	plt.plot( COL_IZ[Sel] , COL_RI[Sel] , marker='.', lw=0, markersize=2, alpha=0.2, color='blue', label='_nolegend_' )
	SelC 	= Sel * ( COL_RI > 1.2 ) * ( COL_IZ < 0.7 ) * ( COL_RI > 1.5*COL_IZ + 1.0 )
	plt.plot( COL_IZ[SelC], COL_RI[SelC], marker='.', lw=0, markersize=4, alpha=0.5, color='blue', label='$\mathrm{Aperture\;1.5"}$' )
	#
	COL_IZ 	= -2.5*np.log10( T['HSC-I2_flux_kron'] ) - -2.5*np.log10( T['HSC-Z_flux_kron' ] ) - E_BV*(1.698-1.263)
	COL_RI	= -2.5*np.log10( T['HSC-R2_flux_kron'] ) - -2.5*np.log10( T['HSC-I2_flux_kron'] ) - E_BV*(2.285-1.698)
	plt.plot( COL_IZ[Sel] , COL_RI[Sel] , marker='.', lw=0, markersize=2, alpha=0.2, color='green', label='_nolegend_' )
	SelC 	= Sel * ( COL_RI > 1.2 ) * ( COL_IZ < 0.7 ) * ( COL_RI > 1.5*COL_IZ + 1.0 )
	plt.plot( COL_IZ[SelC], COL_RI[SelC], marker='.', lw=0, markersize=4, alpha=0.5, color='green', label='$\mathrm{Kron}$' )
	#
	COL_IZ 	= -2.5*np.log10( T['HSC-I2_flux_psf'] ) - -2.5*np.log10( T['HSC-Z_flux_psf' ] ) - E_BV*(1.698-1.263)
	COL_RI	= -2.5*np.log10( T['HSC-R2_flux_psf'] ) - -2.5*np.log10( T['HSC-I2_flux_psf'] ) - E_BV*(2.285-1.698)
	plt.plot( COL_IZ[Sel] , COL_RI[Sel] , marker='.', lw=0, markersize=2, alpha=0.2, color='red', label='_nolegend_' )
	SelC 	= Sel * ( COL_RI > 1.2 ) * ( COL_IZ < 0.7 ) * ( COL_RI > 1.5*COL_IZ + 1.0 )
	plt.plot( COL_IZ[SelC], COL_RI[SelC], marker='.', lw=0, markersize=3, alpha=0.5, color='red', label='$\mathrm{PSF}$' )
	#
	COL_IZ 	= -2.5*np.log10( T['HSC-I2_cmodel_flux'] ) - -2.5*np.log10( T['HSC-Z_cmodel_flux' ] ) - E_BV*(1.698-1.263)
	COL_RI	= -2.5*np.log10( T['HSC-R2_cmodel_flux'] ) - -2.5*np.log10( T['HSC-I2_cmodel_flux'] ) - E_BV*(2.285-1.698)
	plt.plot( COL_IZ[Sel] , COL_RI[Sel] , marker='.', lw=0, markersize=2, alpha=0.2, color='purple', label='_nolegend_' )
	SelC 	= Sel * ( COL_RI > 1.2 ) * ( COL_IZ < 0.7 ) * ( COL_RI > 1.5*COL_IZ + 1.0 )
	plt.plot( COL_IZ[SelC], COL_RI[SelC], marker='.', lw=0, markersize=4, alpha=0.5, color='purple', label='$\mathrm{cModel}$' )
	#
	plt.plot( [-99, 0.2/1.5, 0.7, 0.7], [1.2,1.2,2.05,99], lw=1.5, color='dimgray' )
	plt.xlabel("$\mathrm{HSC\:I - HSC\:Z \; in \; mag}$")
	plt.ylabel("$\mathrm{HSC\:R - HSC\:I \; in \; mag}$")
	plt.legend( fontsize=13, numpoints=1, loc='lower right' )
	plt.xlim( -0.6, 1.6 )
	plt.ylim( -0.2, 2.7 )
	plt.tight_layout()

	#<TableColumns names=('Source','RA','Dec','E_BV','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z','PS1_g','PS1_r','PS1_i','PS1_z','PS1_y','POSS1_B','POSS1_R','M1450','SDSS_Z','SDSS_Plate','SDSS_Fiber','SDSS_MJD','SDSS_SN','HMQ_Z','HMQ_TYPE','HMQ_CITE','HMQ_ZCITE','ELQS_Z','ELQS_reference','ELQS_BAL','Check','ID','Ref_Z','Ref_i','URL1','URL2','URL3','M1450i','M1450z','label','Time','Az','Alt','Culm_Alt','Culm_Time','visible')>


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
	if Targets['Drop'][i] != 'g' or True: continue

	##T	= Table.read( Targets['Catalog'][i] )
	T	= Table.read( Targets['label'][0].split('+')[0] + '/HSC_forced_' + Targets['label'][0].split('+')[0]  + '_short.cat' )

	Sel	= np.ones( len(T), dtype='bool' )
	###Sel 	= Sel * ( T['HSC-R2_FLUX_APER'][:,R] / T['HSC-R2_FLUXERR_APER'][:,R] > 5.0 )
	###Sel 	= Sel * ( T['HSC-I2_FLUX_APER'][:,R] / T['HSC-I2_FLUXERR_APER'][:,R] > 5.0 )
	##Sel 	= Sel * ( T['HSC-R2_flux_kron'] / T['HSC-R2_flux_kron_err'] > 5.0 )
	###Sel 	= Sel * ( T['HSC-I2_flux_kron'] / T['HSC-I2_flux_kron_err'] > 5.0 )

	Sel	= Sel * ( T['HSC-R2_FLUX_AUTO'] / T['HSC-R2_FLUXERR_AUTO'] > 5.0 )

	##Sel	= Sel * ( np.all(T['HSC-R2_FLAGS']==False) )
	###Sel	= Sel * ( -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] ) + 2.5*np.log10( 63095734448.0194 ) < 24.7 )
	##Sel	= Sel * ( -2.5*np.log10( T['HSC-R2_flux_kron'][:] ) + 2.5*np.log10( 63095734448.0194 ) < 24.7 )

	Sel	= Sel * ( T['HSC-R2_MAG_AUTO'] < 24.7 )

	##COL_RI 	= -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] )
	##COL_GR	= -2.5*np.log10( T['HSC-G_FLUX_APER' ][:,R] ) - -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] )

	COL_RI		= T['HSC-R2_MAG_AUTO'] - T['HSC-I2_MAG_AUTO']
	COL_GR		= T['HSC-G_MAG_AUTO' ] - T['HSC-R2_MAG_AUTO']


	plt.close()
	plt.plot( COL_RI[Sel], COL_GR[Sel], marker='.', lw=0, markersize=3, alpha=0.3 )
	plt.plot( [-99, 0.2/1.5, 1, 1], [1,1,2.3,99], lw=1.5, color='red' )
	plt.xlabel("$\mathrm{HSC\:R - HSC\:I \; in \; mag}$")
	plt.ylabel("$\mathrm{HSC\:G - HSC\:R \; in \; mag}$")
	plt.xlim( -0.6, 2.1 )
	plt.ylim( -0.2, 2.7 )
	plt.tight_layout()
	plt.savefig( 'ColorSelection_' + Targets['label'][i] + '.pdf' )
	#raw_input()

	Sel = Sel * ( COL_GR > 1.0 )
	Sel = Sel * ( COL_RI < 1.0 )
	Sel = Sel * ( COL_GR > 1.5*COL_RI + 0.8 )

	plt.close()
	Figs	= Findingchart(	Targets['RA'].quantity[i], Targets['Dec'].quantity[i],
				z=Targets['Ref_Z'][i], M1450=Targets['M1450'].quantity[i],
				PA= Targets['PA'].quantity[i], OffsetSky=Targets['OffsetSky'].quantity[i], OffsetFP=Targets['OffsetFP'].quantity[i],
				Date=Time('2019-06-05'), Observatory='Keck', TVG_lim=17.0*u.mag, SDSS=False, Path=None,
				ExtCat=T[Sel], ExtCatRA='HSC-R2_ALPHAWIN_J2000', ExtCatDec='HSC-R2_DELTAWIN_J2000', Radius=5*u.arcsec, ExtCatLabel='$\mathrm{HSC \; g-drop}$' )

	plt.savefig( 'DropoutSelection_' + Targets['label'][i] + '.pdf' )

	#raw_input()

for i, Target in enumerate( Targets ):
	if Targets['label'][i]!='J1111+1336' or True: continue
	E_BV    = Targets['E_BV'][i]

	T	= Table.read( Targets['label'][i].split('+')[0] + '/HSC_forced_' + Targets['label'][i].split('+')[0]  + '_short.cat' )

	Sel	= np.ones( len(T), dtype='bool' )

	#Sel	= Sel * ( T['HSC-R2_FLUX_AUTO'] / T['HSC-R2_FLUXERR_AUTO'] > 5.0 )
	Sel	= Sel * ( T['HSC-R2_MAG_AUTO'] < 24.7 )

	COL_RI		= T['HSC-R2_MAG_AUTO'] - T['HSC-I2_MAG_AUTO']		- E_BV*(2.285-1.698)
	COL_GR		= T['HSC-G_MAG_AUTO' ] - T['HSC-R2_MAG_AUTO']		- E_BV*(3.303-2.285)


	plt.close()
	plt.plot( COL_RI[Sel], COL_GR[Sel], marker='.', lw=0, markersize=3, alpha=0.5 )
	plt.plot( [-99, 0.2, 0.7, 0.7], [1.2,1.2,2.1,99], lw=1.5, color='red' )
	plt.xlabel("$\mathrm{HSC\:I - HSC\:Z \; in \; mag}$")
	plt.ylabel("$\mathrm{HSC\:R - HSC\:I \; in \; mag}$")
	plt.xlim( -0.6, 2.1 )
	plt.ylim( -0.2, 2.7 )
	plt.tight_layout()
	plt.savefig( 'ColorSelection_' + Targets['label'][i] + '.pdf' )

	Sel = Sel * ( COL_GR > 1.0 )
	Sel = Sel * ( COL_RI < 1.0 )
	Sel = Sel * ( COL_GR > 1.5*COL_RI + 0.8 )

	print(T[Sel])

	Sel1	= ( T['HSC-I2_MAG_AUTO']>18 ) * ( T['HSC-I2_MAG_AUTO']<24.7 ) * ( T['HSC-I2_FLAGS']==0 )

	from MakePostageStamps import MakePostageStamp as MakePostageStamp
	CatRA='HSC-R2_ALPHAWIN_J2000'; CatDec='HSC-R2_DELTAWIN_J2000'
	for j, target in enumerate( T ):
		if (Sel*Sel1)[j] == False : continue
		print('\n', T['NUMBER'][j])
		MakePostageStamp(	T[CatRA].quantity[j], T[CatDec].quantity[j], T['Patch'][j], str(T['NUMBER'][j]), Filters=['HSC-R2', 'HSC-I2', 'HSC-Z' ], Mags=T['HSC-R2_MAG_AUTO', 'HSC-I2_MAG_AUTO', 'HSC-Z_MAG_AUTO'][j],
					ImgPath = Targets['label'][i].split('+')[0] + '/deepCoAdd/', OutPath = Targets['label'][i].split('+')[0] + '/Charts/', Size=12*u.arcsec, Radius=3*u.arcsec )

	#raw_input()

	Figs	= Findingchart(	Targets['RA'].quantity[i], Targets['Dec'].quantity[i],
				z=Targets['Ref_Z'][i], M1450=Targets['M1450'].quantity[i],
				PA= Targets['PA'].quantity[i], OffsetSky=Targets['OffsetSky'].quantity[i], OffsetFP=Targets['OffsetFP'].quantity[i],
				Date=Time('2019-06-05'), Observatory='Keck', TVG_lim=17.0*u.mag, SDSS=False, Path=None, ImgPath=Target['label'].split('+')[0]+'/deepCoAdd/',
				ExtCat=T[Sel], ExtCatRA='HSC-R2_ALPHAWIN_J2000', ExtCatDec='HSC-R2_DELTAWIN_J2000', Radius=4*u.arcsec, ExtCatLabel='$\mathrm{HSC \; r-drop}$' )

	plt.savefig( 'DropoutSelection_' + Targets['label'][i] + '.pdf', dpi=800 )

	input()

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

	for j, target in enumerate( T ):
		if (Sel*Sel1)[j] == False : continue
		#print('\n', T['id'][j])
		#MakePostageStamp(	T[CatRA].quantity[j], T[CatDec].quantity[j], T['Patch'][j], str(T['id'][j]), Filters=['HSC-R2', 'HSC-I2', 'HSC-Z' ], Mags=T['HSC-R2_MAG_cmodel', 'HSC-I2_MAG_cmodel', 'HSC-Z_MAG_cmodel'][j],
		#			ImgPath = Targets['label'][i].split('+')[0] + '/deepCoAdd/', OutPath = Targets['label'][i].split('+')[0] + '/Charts4/', Size=12*u.arcsec, Radius=3*u.arcsec )
		MakePostageStamp(	T[CatRA].quantity[j], T[CatDec].quantity[j], T['Patch'][j], str(T['NUMBER'][j]), Filters=['HSC-R2', 'HSC-I2', 'HSC-Z' ], Mags=T['HSC-R2_MAG_ap', 'HSC-I2_MAG_ap', 'HSC-Z_MAG_ap'][j],
					ImgPath = Targets['label'][i].split('+')[0] + '/deepCoAdd/', OutPath = Targets['label'][i].split('+')[0] + '/Charts4/', Size=12*u.arcsec, Radius=1.5*u.arcsec )


	#T['PCode'] = 1000 - T['HSC-I2_MAG_cmodel'].quantity.to(u.mag).value * 10
	
	Figs, figs, Catalog	= Findingchart(	Targets['RA'].quantity[i], Targets['Dec'].quantity[i],
					z=Targets['Ref_Z'][i], M1450=Targets['M1450'].quantity[i],
					PA= Targets['PA'].quantity[i], OffsetSky=Targets['OffsetSky'].quantity[i], OffsetFP=Targets['OffsetFP'].quantity[i],
					Date=Time('2019-06-05'), Observatory='Keck', TVG_lim=17.0*u.mag, SDSS=False, Path=None, ImgPath=Target['label'].split('+')[0]+'/deepCoAdd/', ImgSTR = 'calexp-HSC-I2-0-*.fits',
					ExtCat=T[Sel*Sel1], ExtCatRA='HSC-I2_RA', ExtCatDec='HSC-I2_Dec', Radius=4*u.arcsec, ExtCatLabel='$\mathrm{HSC \; r-drop}$' )

	#input()
	MakeDsimulatorTable( T[Sel], Catalog, Targets['label'][i].split('+')[0] + '_Targets.dat', CatRA='HSC-R2_ALPHAWIN_J2000', CatDec='HSC-R2_DELTAWIN_J2000', CatFilter='HSC-I2', CatMag='MAG_cmodel', GuideMag='i_mean_psf' )
	for j in range(3):
		PA 	= Targets['PA'].quantity[i][j] + 90*u.deg
		Offset	= -1 * np.array([ -1/np.cos(Targets['Dec'].quantity[i]), 1 ]) * np.dot( np.array( [ [ np.cos( PA ), -np.sin( PA ) ], [ np.sin( PA ), np.cos( PA ) ] ] ), Targets['OffsetFP'].quantity[i][j] )
		print(( Targets['RA'].quantity[i] + Offset[0] ).to(u.deg).value/360*24, ( Targets['Dec'].quantity[i] + Offset[1] ).to(u.deg).value, Targets['PA'].quantity[i][j].to(u.deg).value)
	#input()
	plt.savefig( 'DropoutSelection_' + Targets['label'][i] + '_4.pdf', dpi=800 )

for i, Target in enumerate( Targets ):
	#if Targets['label'][i]!='J1438+4314' or Targets['Drop'][i] != 'r' or i!=7: continue
	if Targets['label'][i]!='J1347+4956': continue
	if Targets['label'][i]!='J1111+1336' or True: continue
	E_BV    = Targets['E_BV'][i]

	#T	= Table.read( Targets['Catalog'][i] )
	T	= Table.read( Targets['label'][i].split('+')[0] + '/HSC_forced_' + Targets['label'][i].split('+')[0]  + '_short.cat' )

	Sel	= np.ones( len(T), dtype='bool' )
	###Sel 	= Sel * ( T['HSC-I2_FLUX_APER'][:,R] / T['HSC-I2_FLUXERR_APER'][:,R] > 4.0 )
	###Sel 	= Sel * ( T['HSC-Z_FLUX_APER' ][:,R] / T['HSC-Z_FLUX_APER_err' ][:,R] > 4.0 )
	###Sel 	= Sel * ( T['HSC-I2_flux_kron'] / T['HSC-I2_flux_kron_err'] > 5.0 )
	###Sel 	= Sel * ( T['HSC-Z_flux_kron' ] / T['HSC-Z_flux_kron_err' ] > 5.0 )
	##Sel 	= Sel * ( T['HSC-I2_flux_psf'] / T['HSC-I2_flux_psf_err'] > 5.0 )
	##Sel 	= Sel * ( T['HSC-Z_flux_psf' ] / T['HSC-Z_flux_psf_err' ] > 5.0 )

	#Sel	= Sel * ( T['HSC-I2_FLUX_AUTO'] / T['HSC-I2_FLUXERR_AUTO'] > 5.0 )

	##Sel	= Sel * ( np.all(T['HSC-I2_FLAGS']==False) )
	###Sel	= Sel * ( -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] ) + 2.5*np.log10( 63095734448.0194 ) < 24.7 ) 
	###Sel	= Sel * ( -2.5*np.log10( T['HSC-I2_flux_kron'] ) + 2.5*np.log10( 63095734448.0194 ) < 24.7 ) 
	###Sel	= Sel * ( -2.5*np.log10( T['HSC-I2_flux_psf'] ) + 2.5*np.log10( 63095734448.0194 ) < 25.5 ) 

	Sel	= Sel * ( T['HSC-I2_MAG_AUTO'] < 25.0 )

	###COL_IZ	= -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-Z_FLUX_APER' ][:,R] )
	###COL_RI	= -2.5*np.log10( T['HSC-R2_FLUX_APER'][:,R] ) - -2.5*np.log10( T['HSC-I2_FLUX_APER'][:,R] )
	##COL_IZ	= -2.5*np.log10( T['HSC-I2_flux_psf'][:] ) - -2.5*np.log10( T['HSC-Z_flux_psf' ][:] )
	##COL_RI 	= -2.5*np.log10( T['HSC-R2_flux_psf'][:] ) - -2.5*np.log10( T['HSC-I2_flux_psf'][:] )

	COL_IZ		= T['HSC-I2_MAG_AUTO'] - T['HSC-Z_MAG_AUTO']		- E_BV*(1.698-1.263)
	COL_RI		= T['HSC-R2_MAG_AUTO'] - T['HSC-I2_MAG_AUTO']		- E_BV*(2.285-1.698)


	plt.close()
	plt.plot( COL_IZ[Sel], COL_RI[Sel], marker='.', lw=0, markersize=3, alpha=0.5 )
	plt.plot( [-99, 0.2, 0.7, 0.7], [1.2,1.2,2.1,99], lw=1.5, color='red' )
	plt.xlabel("$\mathrm{HSC\:I - HSC\:Z \; in \; mag}$")
	plt.ylabel("$\mathrm{HSC\:R - HSC\:I \; in \; mag}$")
	plt.xlim( -0.6, 2.1 )
	plt.ylim( -0.2, 2.7 )
	plt.tight_layout()
	plt.savefig( 'ColorSelection_' + Targets['label'][i] + '.pdf' )
	input()

	Sel = Sel * ( COL_RI > 1.2 )
	Sel = Sel * ( COL_IZ < 0.7 )
	Sel = Sel * ( COL_RI > 1.5*COL_IZ + 1.0 )

	print(T[Sel])

	MakeDimulatorTable( T[Sel] )

	Figs	= Findingchart(	Targets['RA'].quantity[i], Targets['Dec'].quantity[i],
				z=Targets['Ref_Z'][i], M1450=Targets['M1450'].quantity[i],
				PA= Targets['PA'].quantity[i], OffsetSky=Targets['OffsetSky'].quantity[i], OffsetFP=Targets['OffsetFP'].quantity[i],
				Date=Time('2019-06-05'), Observatory='Keck', TVG_lim=17.0*u.mag, SDSS=False, Path=None, #ImgPath=Target['label'].split('+')[0]+'/deepCoAdd/',
				ExtCat=T[Sel], ExtCatRA='HSC-I2_ALPHAWIN_J2000', ExtCatDec='HSC-I2_DELTAWIN_J2000', Radius=4*u.arcsec, ExtCatLabel='$\mathrm{HSC \; r-drop}$' )

	plt.savefig( 'DropoutSelection_' + Targets['label'][i] + '.pdf', dpi=800 )

	#raw_input()






















