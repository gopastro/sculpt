import types
from astropy.time import Time
from sculpt.utils import SculptArgumentError
from astropy.coordinates import Angle, SkyCoord, EarthLocation, FK5
import astropy.units as u

#FIX THIS: we can probably get rid of this function, if we use try-except clauses with appropriate exception handling from astropy
def _check_angle(angle):
    if isinstance(angle, Angle) 
    or type(angle) == types.StringType 
    or type(angle) == types.FloatType:
        return Angle(angle).degree
    else:
        raise IdealPyArgumentError('angle', "angle should be of type ephem.Angle or string or float")
    

def azel2radec(az, el, latitude, longitude, date=None, height=2000.0):
    """
    ==============================================================
    WARNING: This will probably break what you thought it would do,
    if you've used it before. It is untested and unverified. The
    documentation below does not reflect what will happen!!

    THIS IS HERE AS LEGACY, BUT IS CONVERTED TO USE ASTROPY 
    UTILITIES (astropy.coordinates.SkyCoords) INSTEAD

    NEW HEIGHT PARAMETER: assumes 2000m elevation.
    PUT IT AT THE END OF THE PARAMETER LIST, SINCE WE DON'T WANT
    TO BREAK THE WAY SOMEONE HAS IT CALLED IN OLDER CODE.

    WILL ACCEPT AN astropy.time.Time OBJECT, OR ACCEPTABLE STRING
    FOR INITIALIZING ONE. NO LONGER ACCEPTS datetime.date OBJECTS
    (THAT DIDN'T MAKE SENSE ANYWAY, IN CONTEXT, IMHO.)
    
    ===== BJS 12 AUG 2015 ========================================

    Given the azimuth and elevation of an object in radians
    (its horizontal coordinates), and the latitude and longitude
    (in degrees) of the observer's location, returns the
    Right Ascension and Declination in ephem.Angle units.
    Simply rendering them as __str__ will give the right stringified
    representation of RA and Dec.

    Example Usage:

       >>> ra, dec = azel2radec('135:23', '75.0', '42.38028', '-72.52361')
       >>> print ra, dec
       >>> (17:48:05.94, 31:00:06.1)
       
    @param az: Azimuth of the source given in radians typically. If
       azimuth is given as a float it is interpreted as radians. It can
       also be specified as a string like: '170:23:20'. Or it can be
       specified as ephem.Angle.
    @type az: ephem.Angle units or string or floating point (radians)
    @param el: Elevation of the source given in radians typically. If
       elevation is given as a float it is interpreted as radians. It can
       also be specified as a string like: '60:23:20'. Or it can be
       specified as ephem.Angle.
    @type el: ephem.Angle units or string or floating point (radians)
    @param latitude: Latitude of the source given in radians typically. If
       latitude is given as a float it is interpreted as radians. It can
       also be specified as a string like: '60:23:20'. Or it can be
       specified as ephem.Angle.
    @type latitude: ephem.Angle units or string or floating point (radians)    
    @param longitude: Longitude of the source given in radians typically. If
       longitude is given as a float it is interpreted as radians. It can
       also be specified as a string like: '60:23:20'. Or it can be
       specified as ephem.Angle. Longitudes west of GMT should be given as
       negative numbers.
    @type longitude: ephem.Angle units or string or floating point (radians)
    @param date: a datetime.date or datetime.datetime instance of date. This is
       an optional entry. If not given, it will assume current time
    @type date: datetime.date or datetime.datetime instance
    @return: a tuple of RA and Dec of the source for the given observer
       location. Units of RA and Dec are ephem.Angle units.
    """

    #FIX THIS: WE SHOULD REALLY BE USING BUILT-IN ANGLE OBJECTS FROM ASTROPY
    latitude = _check_angle(latitude)
    longitude = _check_angle(longitude)
    az = _check_angle(az)
    el = _check_angle(el)
    
    #FIX THIS: THROW THE ABOVE IN A TRY-EXCEPT CLAUSE TO HANDLE AND/OR 
    #          RE-RAISE EXCEPTIONS APPROPRIATELY!
    if date is not None and (isinstance(date, datetime.datetime)): 
        time = Time(date) 
    else:
        time = Time.now()

    c = SkyCoord(alt=el*u.deg, az=az*u.deg, 
                 obstime=time, frame='altaz',
                 location=EarthLocation(lat=latitude*u.deg, 
                                        lon=longitude*u.deg, height=height*u.m))
    return (c.fk5.ra.degree, c.fk5.dec.degree)


def jprecess(ra, dec, epoch="B1950.0"):
    """
    WARNING: THIS WILL PROBABLY NOT WORK ASK YOU EXPECT IT.
             UNTESTED, UNVERIFIED. PLEASE REFER TO THE ORIGINAL DOCSTRING
             BELOW FOR THE EXPECTED BEHAVIOR, AS WELL AS ASTROPY.COORDINATES
             DOCUMENTATION FOR WHERE WE SHOULD GO FROM HERE.

    Given RA and Dec in degrees for a given epoch (default 1950.0), returns
    the precessed RA and Dec in degrees in J2000 coordinates

    Example Usage:

       >>> ra = ten('4:23:43')*360/24.
       >>> dec = ten('65:42:55')
       >>> ra2000, dec2000 = jprecess(ra, dec)
       >>> print ra2000, dec2000

    @param ra: Right Ascension in degrees. Or you can give a vector or list
       of RA coordinates.
    @type ra: float in degrees
    @param dec: Declination in degrees. Or you can give a vector or list of
       Dec coordinates
    @type dec: float in degrees
    @param epoch: Scalar giving epoch of original observations, default 1950.0
    @type epoch: string or float
    @return: a tuple of RA and Dec or list of tuples (if input is a vector)
       precessed to FK5 system of J2000. Angles are in degrees.
    """
    c = SkyCoord(ra=ra, dec=dec, unit='deg', frame="fk5", equinox=epoch)
    return (c.fk5.ra.degree, c.fk5.dec.degree)

def bprecess(ra, dec, epoch="J2000.0"):
    """
    WARNING: THIS WILL PROBABLY NOT WORK ASK YOU EXPECT IT.
             UNTESTED, UNVERIFIED. PLEASE REFER TO THE ORIGINAL DOCSTRING
             BELOW FOR THE EXPECTED BEHAVIOR, AS WELL AS ASTROPY.COORDINATES
             DOCUMENTATION FOR WHERE WE SHOULD GO FROM HERE.

    Given RA and Dec in degrees for a given epoch (default 2000.0), returns
    the precessed RA and Dec in degrees in B1950 coordinates (FK4)

    Example Usage:

       >>> ra = ten('4:23:43')*360/24.
       >>> dec = ten('65:42:55')
       >>> ra1950, dec1950 = jprecess(ra, dec)
       >>> print ra1950, dec1950

    @param ra: Right Ascension in degrees. Or you can give a vector or list
       of RA coordinates.
    @type ra: float in degrees
    @param dec: Declination in degrees. Or you can give a vector or list of
       Dec coordinates
    @type dec: float in degrees
    @param epoch: Scalar giving epoch of original observations, default 2000.0
    @type epoch: string or float
    @return: a tuple of RA and Dec or list of tuples (if input is a vector)
       precessed to FK4 system of B1950. Angles are in degrees.
    """
    c = SkyCoord(ra=ra, dec=dec, unit='deg', frame="fk5", equinox=epoch)
    c1950 = c.transform_to(FK5(equinox="B1950.0"))
    return (c1950.ra.degree, c1950.dec.degree)

