import sys
import numpy as np
import ephem
from astropy import units as u
from astropy.time import Time
from astropy.table import Table, Column, vstack
import matplotlib.dates

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, search_around_sky
# from SkyCoordinates import SkyCoordinates_from_name

import matplotlib.pyplot as plt

#################################################################################

Observatories = {'LSW': EarthLocation(lat=' 49d23m55s', lon='   0h34m53s', height=590 * u.m),
                 'CAHA': EarthLocation(lat=' 37d13m25s', lon='  -2d32m46s', height=2168 * u.m),
                 'Paranal': EarthLocation(lat='-24d37m33s', lon=' -70d24m11s', height=2635 * u.m),
                 'VLT': EarthLocation(lat='-24d37m33s', lon=' -70d24m11s', height=2635 * u.m),
                 'LBT': EarthLocation(lat=' 32d42m05s', lon='-109d53m22s', height=3267 * u.m),
                 'Keck': EarthLocation(lat=' 19d49m36s', lon='-155d28m24s', height=4123 * u.m),
                 'SB': EarthLocation(lat=' 34.44074d', lon='-119.71569d', height=100 * u.m),
                 'Lick': EarthLocation(lat=' 37d20m24s', lon='-121.63666d', height=1200 * u.m) }


# SpecPhotStandards	= Table.read('SpecPhotStandards.fits')

#################################################################################



def PyStarAlt( Targets, Date, Observatory, Sampling=10*u.min, ax1=None ):

    LOC = Observatories[Observatory]
    if 'ICRS'  not in Targets.columns:	Targets['ICRS']  = SkyCoord( Targets['RA'], Targets['Dec'], unit=u.deg, frame='icrs' )

    ### Compute Sun

    solTimes	= Date + 0* u.day + np.linspace(-12, 24, 36 * 60) * u.hour
    solPTimes = matplotlib.dates.date2num(solTimes.datetime)
    sunaltazs = get_sun(solTimes).transform_to(AltAz(obstime=solTimes, location=LOC))

    Night = solTimes[sunaltazs.alt < -15 * u.deg]
    Midnight = solPTimes[np.isclose(sunaltazs.alt.deg, np.min(sunaltazs.alt.deg))][-1]
    Darktime = np.min(
        solTimes[(solPTimes > Midnight - 12 / 24) * (solPTimes < Midnight) * (sunaltazs.alt < -15 * u.deg)]), np.max(
        solTimes[(solPTimes < Midnight + 12 / 24) * (solPTimes > Midnight) * (sunaltazs.alt < -15 * u.deg)])
    Sunset = np.min(
        solTimes[(solPTimes > Midnight - 12 / 24) * (solPTimes < Midnight) * (sunaltazs.alt < 0 * u.deg)]), np.max(
        solTimes[(solPTimes < Midnight + 12 / 24) * (solPTimes > Midnight) * (sunaltazs.alt < 0 * u.deg)])
    MidNight = Sunset[0] + 0.5 * (Sunset[1] - Sunset[0]).value * u.day

    ### Compute Objects
    print('Calculating Altitudes...   ', Time.now())

    Times = Sunset[0] + np.arange(-1 / 24, (Sunset[1] - Sunset[0]).value + 1 / 24, Sampling.to(u.day).value) * u.day
    PTimes = matplotlib.dates.date2num(Times.datetime)

    sunaltaz = get_sun(Times).transform_to(AltAz(obstime=Times, location=LOC))
    AltAzm = [radec.transform_to(AltAz(obstime=Times, location=LOC)) for radec in Targets['ICRS']]

    Targets['Time'] = PTimes[np.newaxis, :]
    Targets['Az'] = np.array([altaz.az for altaz in AltAzm]) * u.deg
    Targets['Alt'] = np.array([altaz.alt for altaz in AltAzm]) * u.deg
    Targets['Culm_Alt'] = np.max(Targets['Alt'].quantity[:, sunaltaz.alt < -10 * u.deg], axis=1)
    Targets['Culm_Time'] = [np.mean(PTimes[np.isclose(Targets['Alt'][i], Targets['Culm_Alt'][i])]) for i in
                            range(len(Targets))]

    ### plot visibilities ###########################################################

    if ax1 == None: fig, ax1 = plt.subplots(figsize=(16, 9))

    ### Plot Sun

    ax1.plot(solPTimes, sunaltazs.alt, color='orange', label='Sun')
    ax1.fill_between(solPTimes, 0, 90, sunaltazs.alt < 90 * u.deg, color='skyblue', zorder=0, alpha=1.0)
    ax1.fill_between(solPTimes, 0, 90, sunaltazs.alt < 0 * u.deg, color='darkblue', zorder=0, alpha=0.8)
    ax1.fill_between(solPTimes, 0, 90, sunaltazs.alt < -5 * u.deg, color='black', zorder=0, alpha=0.4)
    ax1.fill_between(solPTimes, 0, 90, sunaltazs.alt < -10 * u.deg, color='black', zorder=0, alpha=0.4)
    ax1.fill_between(solPTimes, 0, 90, sunaltazs.alt < -15 * u.deg, color='black', zorder=0, alpha=0.4)

    ### Plot Moon

    MoonPos = SkyCoord(ephem.Moon(MidNight.datetime).ra * u.rad, ephem.Moon(MidNight.datetime).dec * u.rad,
                       frame='icrs')
    ax1.plot(PTimes, MoonPos.transform_to(AltAz(obstime=Times, location=LOC)).alt, ls='--', lw=0.7, color='yellow',
             label='Moon')

    ### Plot Targets

    for i, obj in enumerate(Targets):
        ax1.plot_date(Targets['Time'][i], Targets['Alt'].quantity[i], '-', color='w')

    ax1.set_xlim(matplotlib.dates.date2num((Sunset[0] - 0.5 * u.hour).datetime),
                 matplotlib.dates.date2num((Sunset[1] + 0.5 * u.hour).datetime))
    ax1.set_ylim(0, 90)
    ax1.set_yticks(np.linspace(0, 90, 4))
    ax1.grid(color='white')
    ax1.text(0.5, 1.02,
             'Sky on {Date:} -- Observatory:  {Obs:}  {Lat:8.4f} / {Lon:8.4f}'.format(Date=Date.iso.split()[0],
                                                                                      Obs=Observatory,
                                                                                      Lat=LOC.latitude.value,
                                                                                      Lon=LOC.longitude.value),
             ha='center', va='bottom', transform=ax1.transAxes, fontsize=12)
    ax1.set_xlabel('UTC', fontsize=12)
    ax1.set_ylabel('Altitude [deg]', fontsize=12)

    ax1.tick_params(labelsize=12)
    ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

    # plt.subplots_adjust( top=0.96, bottom=0.06, left=0.04, right=0.95 )

    return Targets

#################################################################################


