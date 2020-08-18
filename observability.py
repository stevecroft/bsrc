#!/usr/bin/env python

import numpy as np
import scipy as sp
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse

parser = argparse.ArgumentParser(description="Check observability at GBO and (default) Parkes")
parser.add_argument('source',help="Source name (e.g. LHS1140)")
parser.add_argument('-date',help="Date (format yyyy-mm-dd)")
parser.add_argument('-mwa',help="Plot MWA instead of Parkes",action="store_true", default=False)
parser.add_argument('-utc',help="Plot UTC instead of PST",action="store_true", default=False)
parser.add_argument('-pdt',help="Plot PDT instead of PST",action="store_true", default=False)
args = parser.parse_args()

#target = SkyCoord.from_name('LHS1140')
target = SkyCoord.from_name(args.source)

print("Target RA=",target.ra," Dec=",target.dec)
# Parkes and GBO lat/long

if args.mwa:
	tel2name = 'MWA'
	tel2 = EarthLocation(lat=-26.7*u.deg, lon=116.7*u.deg, height=377*u.m)
	tel2utcoffset = 8*u.hour
else:
	tel2name = 'Parkes'
	tel2 = EarthLocation(lat=-33.1*u.deg, lon=148.2*u.deg, height=324*u.m)
	tel2utcoffset = 10*u.hour # tel2 UTC offset

gbo = EarthLocation(lat=38.4*u.deg, lon=-79.8*u.deg, height=808*u.m)

pstutoffset = -8 * u.hour # Pacific time UTC offset
pdtutoffset = -7 * u.hour # Pacific daylight time UTC offset
gboutcoffset = -5*u.hour # GBO UTC offset

# midnight on the specified day
# default to tomorrow if not specified
if not args.date:
	tstart = Time.now().iso[0:10]
else:
	tstart = Time(args.date).iso[0:10]
print("Date "+str(tstart))
utcmidnight = Time(tstart+' 00:00:00')
pstmidnight = utcmidnight - pstutoffset
pdtmidnight = utcmidnight - pdtutoffset
tel2midnight = utcmidnight - tel2utcoffset
gbomidnight = utcmidnight - gboutcoffset

# divide the day into 1000 timesteps for plotting
delta_midnight = np.linspace(0, 24, 1000)*u.hour
# defaults to plotting times in PST, but replace with tel2midnight, gbomidnight, or
# utcmidnight in the following line as desired
if args.utc:
	times = utcmidnight+delta_midnight
	tzone = 'UTC'
elif args.pdt:
	times = pdtmidnight+delta_midnight
	tzone = 'PDT'
else:
	times = pstmidnight+delta_midnight
	tzone = 'PST'

# compute altitude and azimuth of the target for tel2
altazframe_tel2 = AltAz(obstime=times, location=tel2)
targetaltazs_tel2 = target.transform_to(altazframe_tel2)

# can compute Sun altitude too if desired
#sunaltazs = get_sun(times).transform_to(altazframe)
alts = targetaltazs_tel2.alt
# find times where the altitude crosses 30 degrees (elevation limit at tel2)
altz = alts.deg - 30.
altzc = delta_midnight[np.where(altz[:-1] * altz[1:] < 0)]
try:
	altz1t = min(altzc)
	altz2t = max(altzc)
	altz1 = str(TimeDelta(altz1t, format='sec').to_datetime())
	altz2 = str(TimeDelta(altz2t, format='sec').to_datetime())
	print("Target crosses 30 degrees at %s at %s %s and %s" % (tel2name, tzone, altz1,altz2))
except:
	altz1 = 0 * u.h
	altz2 = 24 * u.h
	print("Target never crosses 30 degrees at",tel2name)

fig1 = plt.figure(args.source, figsize=(10,4))
ax1 = fig1.add_subplot(121)
plt.plot(delta_midnight, alts)
#plt.xticks(range(0, 24, 1))
# can plot Sun altitude too if desired
#plt.plot(delta_midnight, sunaltazs.alt, color='y', label='Sun')
plt.xlim(0, 24) # time 0 - 24 hours
plt.ylim(0, 90) # alt 0 - 90 deg
plt.xlabel(tzone)
plt.ylabel('Altitude [deg]')
plt.title(tel2name)
ax1.add_patch(patches.Rectangle(
        (0., 0.), 24, 30,
        hatch='/', fill=False))
#ax1.add_patch(patches.Rectangle(
#        (0., 0.), altz1 / u.h, 30,
#        hatch='/', color='red', fill=True))
#ax1.add_patch(patches.Rectangle(
#        (altz2 / u.h, 0.), 24, 30,
#        hatch='/', color='red', fill=True))


altazframe_gbo = AltAz(obstime=times, location=gbo)
targetaltazs_gbo = target.transform_to(altazframe_gbo)
alts = targetaltazs_gbo.alt
altz = alts.deg - 5.
altzc = delta_midnight[np.where(altz[:-1] * altz[1:] < 0)]
try:
	altz1t = min(altzc)
	altz2t = max(altzc)
	altz1 = str(TimeDelta(altz1t, format='sec').to_datetime())
	altz2 = str(TimeDelta(altz2t, format='sec').to_datetime())
	print("Target crosses 5 degrees at GBO at %s %s and %s" % (tzone, altz1,altz2))
except:
	altz1 = 0 * u.h
	altz2 = 24 * u.h
	print("Target never crosses 5 degrees at GBO")

#fig2 = plt.figure(2)
ax2 = fig1.add_subplot(122)
plt.plot(delta_midnight, alts)
plt.xlim(0, 24)
plt.ylim(0, 90)
plt.xlabel(tzone)
plt.ylabel('Altitude [deg]')
plt.title('GBO')
ax2.add_patch(patches.Rectangle(
        (0., 0.), 24, 5,
        hatch='/', fill=False))
#ax2.add_patch(patches.Rectangle(
#        (0., 0.), altz1 / u.h, 5,
#        hatch='/', color='red', fill=True))
#ax2.add_patch(patches.Rectangle(
#        (altz2 / u.h, 0.), 24, 5,
#        hatch='/', color='red', fill=True))
plt.show()

# find times where the altitude crosses 5 degrees (elevation limit at GBO)
