###############################################################################
# Py_Findingchart - Optimized for Deimos Imaging
# Version 1.0
# 2018-04-08
# Authors: Tobias Schmidt
#
#
# Py_Findingchart creates findingcharts based on DSS and SDSS imaging (if available)
# In addition, offset stars are selected from the Gaia astrometric catalog
#
###############################################################################

from __future__ import division

import os
import glob
import time
import sys
import argparse

import numpy as np

import astropy.units as u
from astropy.units import Quantity
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, SkyOffsetFrame, ICRS
# from SkyCoordinates import SkyCoordinates_from_name
from astropy.table import Table

import matplotlib as mpl
from matplotlib import pylab
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt

import aplpy

from PyStarAltFC import PyStarAlt
from IPython import embed

np.seterr(all='ignore')

#################################################################################################################################

### Important: a pure grayscale will be reduced and alpha fucked up!
cdict2 = {'red': [[0.0, 1.0, 1.0],
                  [0.9, .09, .09],
                  [1.0, 0.0, 0.0]],

          'green': [[0.0, 1.0, 1.0],
                    [1.0, 0.0, 0.0]],

          'blue': [[0.0, 1.0, 1.0],
                   [1.0, 0.0, 0.0]]
          }
cdict2['alpha'] = [[0.0, 0.0, 1.0],
                   [1.0, 1.0, 0.0]]

# plt.register_cmap('GrayAlpha', cdict2)
newcmp = LinearSegmentedColormap('GrayAlpha', segmentdata=cdict2, N=256)


#################################################################################################################################

def Get_SkyView_Image_pl(RA, Dec, Size=5 * u.arcmin, Scale=1 * u.arcsec, Survey='DSS2R'):
    # POSS image query
    # download DSS images via some PERL script (sic!) from http://skyview.gsfc.nasa.gov/current/docs/batchpage.html
    # https://skyview.gsfc.nasa.gov/current/cgi/titlepage.pl

    Pos = SkyCoord(RA, Dec, frame='icrs')
    File = "SkyView_Query_{:11.0f}.fits".format(time.time() * 100)
    Command = "./skybatch_wget.pl file=" + File + \
              " Position='" + Pos.ra.to_string(unit=u.h, sep=' ', precision=0, pad=True) + \
              ", " + Pos.dec.to_string(sep=' ', precision=0, pad=True, alwayssign=True) + "'" \
                                                                                          ", Size={:4.2f}, Pixels={:3.0f}".format(
        Size.to(u.deg).value, (Size / Scale).to('').value) + \
              ", Survey='" + Survey + "'" + ", Scaling='Linear'"
    print(Command)
    os.system(Command)

    DSS_Image = fits.open(File)[0]
    os.remove(File)

    return DSS_Image


#################################################################################################################################

def Get_SkyView_Image(RA, Dec, Size=5 * u.arcmin, Scale=1 * u.arcsec, Survey='DSS2R'):
    # query various surveys from https://skyview.gsfc.nasa.gov/current/cgi/titlepage.pl
    import time, os
    from astropy.io import fits

    Pos = SkyCoord(RA, Dec, frame='icrs')
    File = "SkyView_Query_{:11.0f}.fits".format(time.time() * 100)

    URL = 'http://skyview.gsfc.nasa.gov//current/cgi/pskcall?'
    URL += 'Position={:},'.format(
        Pos.to_string('hmsdms', precision=1).split(' ')[0].replace('h', '+').replace('m', '+').replace('s', ''))
    URL += '+{:},'.format(
        Pos.to_string('hmsdms', precision=1).split(' ')[1].replace('d', '+').replace('m', '+').replace('s', ''))
    URL += '&Size={:1.3f},'.format(Size.to(u.deg).value)
    URL += '&Pixels={:1.0f},'.format((Size / Scale).to('').value)
    URL += '&Survey={:},'.format(Survey)
    URL += '&Scaling=Linear'

    os.system('wget -q -O ' + File + ' "' + URL + '"')
    SkyView_Image = fits.open(File)[0]
    os.remove(File)

    return SkyView_Image


#################################################################################################################################

def Get_PanSTARRS(RA, Dec, Size=5 * u.arcmin, nDetections=5):
    File = "PanSTARRS_Query_{:11.0f}.fits".format(time.time() * 100)
    URL = 'https://archive.stsci.edu/panstarrs/search.php?'
    URL += 'RA={:7.5f}&'.format(RA.to(u.deg).value)
    URL += 'DEC={:7.5f}&'.format(Dec.to(u.deg).value)
    URL += 'SR={:5.3f}&'.format(Size.to(u.deg).value)
    URL += 'coordformat=dec&max_records=100000&outputformat=VOTable&action=Search'
    os.system('wget -q -O ' + File + ' "' + URL + '"')
    Catalog = Table.read(File, format='votable')[
        'ramean', 'decmean', 'ndetections', 'gmeanpsfmag', 'rmeanpsfmag', 'imeanpsfmag', 'zmeanpsfmag', 'ymeanpsfmag']
    os.remove(File)
    Catalog = Catalog[Catalog['ndetections'] >= nDetections]
    Catalog.rename_column('ramean', 'RA')
    Catalog.rename_column('decmean', 'Dec')
    Catalog['RA'].unit = u.deg
    Catalog['Dec'].unit = u.deg
    for Band in ['g', 'r', 'i', 'z', 'y']:
        Catalog.rename_column(Band + 'meanpsfmag', Band + '_mean_psf')
        Catalog[Band + '_mean_psf'].unit = u.mag
        Catalog[Band + '_mean_psf'][Catalog[Band + '_mean_psf'] == -999.0] = np.nan
    return Catalog


#################################################################################################################################

def Findingchart( 	RA, Dec, Name=None, mag=0*u.mag, Filter=None, z=None, M1450=None, TVG_lim=18*u.mag,
            PA=np.array([0])*u.deg, OffsetSky=np.array([[0.0,0.0]])*u.arcmin, OffsetFP=np.array([[+0.0,+0.0]])*u.arcmin,
            Path='./Charts', ReDo=True, Size=18*u.arcmin,
            Date=None, Observatory=None, SDSS=False, ImgPath=None, ImgSTR = 'calexp-HSC-R2-0-*.fits',
            ExtCat=None, ExtCatRA='HSC-R2_RA', ExtCatDec='HSC-R2_Dec', Radius=5*u.arcsec, ExtCatLabel='$\mathrm{HSC g-drop}$' ):

    Pos 	= SkyCoord( RA, Dec, frame='icrs' )
    Name    = 'J' + (Pos.ra -             360/24/3600/2*u.deg).to_string(unit=u.h, sep='', precision=1, pad=True)[:-2] +  \
            (Pos.dec-np.sign(Pos.dec.deg)*.5*u.arcsec).to_string(          sep='', precision=1, pad=True, alwayssign=True)[:-2] if Name==None else Name

    if Path != None and os.path.exists(Path + Name + '.pdf') and ReDo == False:
        print('Already exists: ' + Name)
    else:
        print('Processing ' + Name)

    plt.close()
    # SubPlot	= [0.14, 0.28, 0.79, 0.63]
    # Fig 	= plt.figure(figsize=(8.5,10.3))
    SubPlot = [0.12, 0.08, 0.85, 0.83]
    Fig = plt.figure(figsize=(11.0, 11.5))
    # figs	= [ aplpy.FITSFigure( 'Dummy.fits', figure=Fig, subplot=SubPlot) ]
    figs = [aplpy.FITSFigure(
        Get_SkyView_Image(RA=RA, Dec=Dec, Size=Size + .1 * u.arcmin, Scale=1.0 * u.arcsec, Survey='DSS2R'), north=True,
        figure=Fig, subplot=SubPlot)]  ## Dummy to initialize

    if ImgPath!=None:
        print(ImgPath)
        print(glob.glob(ImgPath +ImgSTR))
        # vmin, vmax	= np.nanmean( [ np.nanpercentile( fits.open(File)[1].data, [30, 99] ) for File in glob.glob(ImgPath+ImgSTR) ], axis=0 )
        HDUList		= fits.open( "/home/sbechtel/Downloads/unzip_hold/tobias_code/J1438/deepCoAdd/calexp-HSC-R2-0-5.6.fits")
        vmin, vmax      = np.nanpercentile( HDUList[1].data[ (-2000 < HDUList[1].data) * (HDUList[1].data < 99999) ] , [10, 90] ) * np.array([ 5.0, 20.18 ])
        if Name == 'J163055+043558':
            vmin = 1000
            vmax = 10000
        print(vmin, vmax)
        for i ,File in enumerate(glob.glob(ImgPath +ImgSTR)):
            try:
                Patch	= fits.open(File)[1]
            except:
                Patch	= fits.open(File)[0]
            Header	= Patch.header
            SepRA	= 	( ( Header['CRVAL1'] + Header['CD1_1'] * ( Header['NAXIS1' ] /2 - Header['CRPIX1'] ) )* u.deg - Pos.ra) * np.cos(Pos.dec)
            SepDec = ((Header['CRVAL2'] + Header['CD2_2'] * (
                        Header['NAXIS2'] / 2 - Header['CRPIX2'])) * u.deg - Pos.dec)
            print(File, SepRA.to(u.arcmin), SepDec.to(u.arcmin), np.nanpercentile(Patch.data, [30.0, 99]),
                  bool((np.abs(SepRA) > Size / 2 * 1.8) + (np.abs(SepDec) > Size / 2 * 1.8)))
            # figs[-1].show_rectangles( ( Header['CRVAL1'] + Header['CD1_1'] * ( Header['NAXIS1']/2 - Header['CRPIX1'] ) ), ( Header['CRVAL2'] + Header['CD2_2'] * ( Header['NAXIS2']/2 - Header['CRPIX2'] ) ),
            #				(24*u.arcsec).to(u.deg).value, (24*u.arcsec).to(u.deg).value, color='teal', alpha=1.0 )
            if np.abs(SepRA) > Size / 2 * 1.3 or np.abs(SepDec) > Size / 2 * 1.3:
                continue
            figs += [aplpy.FITSFigure(Patch, north=False, figure=Fig, subplot=SubPlot, downsample=2)]
            figs[-1].show_colorscale(cmap=newcmp, stretch='arcsinh', vmin=vmin, vmax=vmax)
    embed()
    ### GET SDSS Image if possible
    try:
        if SDSS == False: a = b
        print("Getting SDSS r")
        SDSS = Get_SkyView_Image(RA=RA, Dec=Dec, Size=Size + 0.5 * u.arcmin, Scale=0.396 * 4 * u.arcsec, Survey='SDSSr')
        SDSS.data[SDSS.data == 0.0] = np.nan
        figs += [aplpy.FITSFigure(SDSS, north=True, figure=Fig, subplot=SubPlot)]
        figs[-1].show_colorscale(cmap=newcmp, stretch='arcsinh', vmin=np.nanpercentile(SDSS.data, 24.0),
                                 vmax=np.nanpercentile(SDSS.data, 99.925))
        plt.gca().plot([np.nan, np.nan], ls='', lw=1.5, label='$\mathrm{SDSS\;r\;band}$')
    except:
        plt.gca().plot([np.nan, np.nan], ls='', lw=1.5, label='$\mathrm{DSS\;r\;band}$')

    if 1:
        print("Getting PanSTARRS")
        Catalog = Get_PanSTARRS(RA, Dec, Size / np.sqrt(2))
        Catalog = Catalog[Catalog['r_mean_psf'].quantity < TVG_lim]
        plt.gca().plot([np.nan], [np.nan], marker='s', ls='', markersize=9, markerfacecolor='none',
                       markeredgecolor='limegreen',
                       label='$\mathrm{{PS1 \; r < {:4.1f}\,mag}}$'.format(TVG_lim.to(u.mag).value))
    figs[-1].show_rectangles(Catalog['RA'].to(u.deg).value, Catalog['Dec'].to(u.deg).value,
                             np.ones(len(Catalog)) * (12 * u.arcsec).to(u.deg).value,
                             np.ones(len(Catalog)) * (12 * u.arcsec).to(u.deg).value, color='limegreen', alpha=0.7)

    ### Plot Catalog
    try:
        if ExtCat is not None:
            ExtCatRA = 'HSC-R2_ALPHAWIN_J2000'
            ExtCatDec = 'HSC-R2_DELTAWIN_J2000'
            figs[-1].show_circles(ExtCat[ExtCatRA].quantity.to(u.deg).value, ExtCat[ExtCatDec].quantity.to(u.deg).value,
                                  u.Quantity(np.ones(len(ExtCat)) * [Radius.to(u.deg).value]), color='DarkOrange',
                                  alpha=1.0)
            plt.gca().plot([np.nan], [np.nan], marker='o', ls='', markersize=8, lw=1.8, markerfacecolor='none',
                           markeredgecolor='DarkOrange', label=ExtCatLabel)
    except:
        if ExtCat is not None:
            ExtCatRA = 'RA'
            ExtCatDec = 'Dec'
            figs[-1].show_circles(ExtCat[ExtCatRA], ExtCat[ExtCatDec],
                                  u.Quantity(np.ones(len(ExtCat)) * [Radius.to(u.deg).value]), color='DarkOrange',
                                  alpha=1.0)
            plt.gca().plot([np.nan], [np.nan], marker='o', ls='', markersize=8, lw=1.8, markerfacecolor='none',
                           markeredgecolor='DarkOrange', label=ExtCatLabel)

    ### Mark Target
    figs[-1].show_circles(Pos.ra.deg, Pos.dec.deg, (9 * u.arcsec).to(u.deg).value, lw=1.4, color='teal')
    figs[-1].show_circles(RA.to(u.deg).value, Dec.to(u.deg).value, (8 * u.arcsec).to(u.deg).value, lw=1.4, color='red')
    figs[-1].add_label((RA - 0 * u.arcsec).to(u.deg).value, (Dec - 25 * u.arcsec).to(u.deg).value, Name, ha='center',
                       va='top', fontsize=11)

    ### Make Legend
    plt.gca().legend(loc='upper right', fontsize=10, frameon=True, numpoints=1, scatterpoints=1)

    ### Plot FoV
    Target = Quantity([RA, Dec])
    PathFoV = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
    ## https://www2.keck.hawaii.edu/inst/KSDs/40/html/ksd40-55.2.html#pgfId=771977
    # PA0		= 179*u.deg								## Statement from Carlos Alvarez, instrument angle i=179, e.g at PA=90 XIM lines up with RA and YIM line sup with Dec
    PA0 = 270 * u.deg  ## I define the PA in the standard way: PA=0 corresponds to slits in N-S, Guider in the East, Vignetted Edges in the West!
    ## This is the same as in DSIMULATOR. FoV definition has to be rotated by 270 to reach this position
    FoV = np.loadtxt(PathFoV + 'Deimos_FoV.dat') * [-1,
                                                    1] * u.arcsec  ## The downloaded files seem to be flipped in XIM but correct in YIM
    TVG = np.loadtxt(PathFoV + 'Deimos_TVG.dat') * [-1,
                                                    1] * u.arcsec  ## The downloaded files seem to be flipped in XIM but correct in YIM
    AimpointIM = Quantity([-95.20 * u.mm, 202.10 * u.mm]) / (0.728 * u.mm / u.arcsec)
    AimpointTV = Quantity([-73.25 * u.mm, 100.00 * u.mm]) / (0.728 * u.mm / u.arcsec)
    PlateScale = 0.728 * u.mm / u.arcsec
    for pa, offsetSky, offsetFP in zip(PA, OffsetSky, OffsetFP):
        AimIM = Target + offsetSky + np.array([-1 / np.cos(Dec), 1]) * np.dot(
            [[np.cos(pa + PA0), -np.sin(pa + PA0)], [np.sin(pa + PA0), np.cos(pa + PA0)]], (offsetFP).T).T
        AimTV = Target + offsetSky + np.array([-1 / np.cos(Dec), 1]) * np.dot(
            [[np.cos(pa + PA0), -np.sin(pa + PA0)], [np.sin(pa + PA0), np.cos(pa + PA0)]],
            (offsetFP - AimpointIM + AimpointTV).T).T
        PolygonFoV = Target + offsetSky + np.array([-1 / np.cos(Dec), 1]) * np.dot(
            [[np.cos(pa + PA0), -np.sin(pa + PA0)], [np.sin(pa + PA0), np.cos(pa + PA0)]],
            (offsetFP - AimpointIM + FoV).T).T
        PolygonTVG = Target + offsetSky + np.array([-1 / np.cos(Dec), 1]) * np.dot(
            [[np.cos(pa + PA0), -np.sin(pa + PA0)], [np.sin(pa + PA0), np.cos(pa + PA0)]],
            (offsetFP - AimpointIM + TVG).T).T
        figs[-1].show_polygons([PolygonFoV.to(u.deg).value], color='navy', lw=0.8, alpha=0.8)
        figs[-1].show_polygons([PolygonTVG.to(u.deg).value], color='crimson', lw=0.8, alpha=0.8)
        figs[-1].show_markers(AimIM[0].to(u.deg).value, AimIM[1].to(u.deg).value, marker='+', lw=0.8, facecolor='navy',
                              edgecolor='navy', s=12 ** 2)
        figs[-1].show_markers(AimTV[0].to(u.deg).value, AimTV[1].to(u.deg).value, marker='+', lw=0.8,
                              facecolor='crimson', edgecolor='crimson', s=12 ** 2)
        print(AimIM, AimTV)

    ### add title
    Title = '\LARGE{' + Name + '}'
    Title += '\nTarget:  ' + Pos.to_string('hmsdms', precision=0) + ' -- ' + '{:8.4f}d {:+7.4f}d'.format(Pos.ra.deg,
                                                                                                         Pos.dec.deg)
    Title += '\nMagnitudes: $\mathrm{{ m_{{ {:} }} = {:4.1f}\,mag}}$'.format(Filter.replace(' ', '\;'),
                                                                             mag.to(u.mag).value) if (
                Filter != None and mag != None) else '\nMagnitude unknown'
    Title += ', $\mathrm{{ M_{{1450}} = {:4.1f}\, mag}}$'.format(M1450.to(u.mag).value) if (M1450 != None) else ''
    Title += ' -- Redshift: $\mathrm{{ z = {:4.2f} }}$'.format(z) if (z != None) else ' -- Redshift unknown'
    plt.title(Title)

    ### adjust figures

    for i, axis in enumerate(Fig.axes):
        Fig.axes[i].patch.set_visible(False)

    for i, fig in enumerate(figs):
        fig.set_system_latex(True)
        fig.ticks.set_color('black')
        fig.axis_labels.hide()
        fig.tick_labels.hide()
        fig.tick_labels.set_xformat('hh:mm:ss')
        fig.tick_labels.set_yformat('dd:mm')
        [fig.recenter(Pos.ra.deg, Pos.dec.deg, (Size / 2 * f).to(u.deg).value) for f in
         [2, 1.8, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1]]

    figs[-1].axis_labels.show()
    figs[-1].tick_labels.show()

    ### Make StarAlt Plot
    '''
    if Date != None and Observatory != None and False:
        ax_Alt = Fig.add_axes((0.14, 0.05, 0.79, 0.14))
        T = Table();
        T['ICRS'] = [Pos]
        PyStarAlt(T, Date, Observatory, Sampling=10 * u.min, ax1=ax_Alt)
    '''
    ### safe figure
    if Path != None:
        Fig.savefig(Path + Name + '_4.pdf')
    if True:
        Fig.savefig('current_chart.pdf')
    return Fig, figs, Catalog


#################################################################################################################################

if __name__ == '__main__':

    ## Possible Call Syntax:
    ## ./PyFindingchart.py -RA=12h56m11s -Dec=-05d47m22s -Observatory=Keck -PA=90 -SDSS=True
    ## ./PyFindingchart.py -RA 00:42:44 -Dec +41:14:09 -Observatory Keck -Date 2018-04-21
    ## ./PyFindingchart.py -RA 0.0d -Dec +1.23d -Filter r -mag 17.5 -z 0.15 -Observatory Keck -Date 2018-04-21
    ## ./PyFindingchart.py -FromName 3C273 -z 0.15 -Observatory Keck -OffsetFP=-4.5,-0.7

    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('-RA', '--RA', required=False)
    parser.add_argument('-Dec', '--Dec', required=False)
    parser.add_argument('-Name', '--Name', required=False)
    parser.add_argument('-mag', '--mag', required=False)
    parser.add_argument('-Filter', '--Filter', required=False)
    parser.add_argument('-z', '--z', required=False)
    parser.add_argument('-FromName', '--FromName', required=False)
    parser.add_argument('-Observatory', '--Observatory', required=False)
    parser.add_argument('-Date', '--Date', required=False)
    parser.add_argument('-SDSS', '--SDSS', required=False)
    parser.add_argument('-PA', '--PA', required=False)
    parser.add_argument('-OffsetSky', '--OffsetSky', required=False)
    parser.add_argument('-OffsetFP', '--OffsetFP', required=False)

    args = vars(parser.parse_args())

    print('\nArguments:')
    print(args, '\n')
    '''
    ## Try resolvng Name first and set coordinates
    if args['FromName']!=None:
        ICRS	= SkyCoordinates_from_name( args['FromName'] )
        print('Using Object Name', args['FromName'])
        #ICRS.ra, ICRS.dec
        args['RA']	= '{:7.5f}d'.format(ICRS.ra.deg)
        args['Dec']	= '{:7.5f}d'.format(ICRS.dec.deg)
        args['Name']	= args['FromName']
    '''

    ## If Coordinates given or could be resolved, go for it
    if args['RA'] != None and args['Dec'] != None:
        try:
            ICRS = SkyCoord(args['RA'], args['Dec'], frame='icrs')
        except:
            ICRS = SkyCoord(args['RA'], args['Dec'], frame='icrs', unit=(u.hourangle, u.deg))
        print('Using Coordinates\n', ICRS, '\n')
    '''
    ## If no magnitude or Filter is given, try to get at least Gaia mags
    if args['mag'] == None or args['Filter'] == None:
        Catalog = Gaia_Query(ICRS.ra, ICRS.dec, 10 * u.arcsec)['ra', 'dec', 'dist', 'phot_g_mean_mag']
        Catalog.sort('dist')
        if len(Catalog):
            args['mag'] = Catalog['phot_g_mean_mag'][0]
            args['Filter'] = 'Gaia g'

        ## Make Chart
        Figs = Findingchart(ICRS.ra, ICRS.dec,
                            Name=args['Name'], Filter=args['Filter'],
                            mag=float(args['mag']) * u.mag if args['mag'] != None else None,
                            z=float(args['z']) if args['z'] != None else None,
                            PA=np.array([float(args['PA'])]) * u.deg if args['PA'] != None else np.array([0.0]) * u.deg,
                            OffsetSky=np.array([[float(f) for f in args['OffsetSky'].split(',')]]) * u.arcmin if args[
                                                                                                                     'OffsetSky'] != None else np.array(
                                [[0.0, 0.0]]) * u.arcmin,
                            OffsetFP=np.array([[float(f) for f in args['OffsetFP'].split(',')]]) * u.arcmin if args[
                                                                                                                   'OffsetFP'] != None else np.array(
                                [[0.0, 0.0]]) * u.arcmin,
                            Date=Time(args['Date']) if args['Date'] != None else Time.now(),
                            Observatory=args['Observatory'],
                            Path=None, SDSS=bool(args['SDSS']), Size=19 * u.arcmin)

        raw_input()'''
#################################################################################################################################



