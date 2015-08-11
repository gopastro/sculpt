"""
A few functions dealing with formatted RA and Dec coordinates.
"""
import types
import numpy
import decimal
import math

def adstring (*args, **kwargs):
    """
    Return RA and Dec as character string(s) in sexigesimal format.
    
    RA and Dec may be entered as either a 2 element list or as
    two separate lists (or scalars).  One can also specify the precision 
    of the declination in digits after the decimal point.
    @param ra_dec: the first parameter if it is the only parameter
        given is a 2 element list giving the Right Ascension and
        Declination in decimal degrees.
    @param ra: if args consist of two elements, ra is first element. ra can
        be numeric scalar or list of RAs.
    @param dec: if args consist of two elements, dec is second element.
        dec can be numeric scalar or list of decs.
    @param precision: optional input. Either the third element or keyword
        argument. If absent, defaults to 1. Integer scalar (0-4) giving the
        number of digits after the decimal of DEClination. The RA is
        automatically 1 digit more. It is not available for just DEC.
        If no PRECISION parameter is passed, a precision of 1 is used.
        Values of precision larger than 4 will be truncated to 4.
    @return: Character string(s) containing HR,MIN,SEC,DEC,MIN,SEC formatted
        as ( 2I3,F5.(p+1),2I3,F4.p ) where p is the PRECISION parameter.
        If only a single scalar is supplied it is converted to a
        sexigesimal string (2I3,F5.1).
       """
    precision = kwargs.get('precision', 1)
    nargs = len(args)
    if nargs == 1:
        outlist = False
        ra, dec = args[0]
    elif nargs in (2, 3):
        ra, dec = args[0], args[1]
        if nargs == 3:
            precision = args[2]
            if precision > 4:
                precision = 4
            elif precision < 0:
                precision = 0
        if type(ra) in (types.ListType, types.TupleType):
            outlist = True
            if len(ra) != len(dec):
                print "Ra and Dec lists should be of same length!"
                return None
        else:
            outlist = False
    else:
        print "Input argument list can only be of length 1, 2 or 3"
        return None

    #rprec = decimal.Decimal(("%f" % math.pow(10.0, -(precision+1))).rstrip('0'))
    #dprec = decimal.Decimal(("%f" % math.pow(10.0, -precision)).rstrip('0'))
    rprec = precision+1
    dprec = precision
    if outlist:
        #ra and dec are lists
        for i, r in enumerate(ra):
            ihr, imin, xsec, ideg, imn, xsc = radec(r, dec[i])
            outstr = "%02d:%02d:%02d.%s  %02d:%02d:%02d.%s" % (ihr, imin, int(math.floor(xsec)), ("%f" % round((xsec-math.floor(xsec)), rprec))[2:2+rprec], ideg, imn, int(math.floor(xsc)), ("%f" % round((xsc-math.floor(xsc)), dprec))[2:2+dprec])
            print outstr
            return outstr
    else:
        ihr, imin, xsec, ideg, imn, xsc = radec(ra, dec)
        outstr = "%02d:%02d:%02d.%s  %02d:%02d:%02d.%s" % (ihr, imin, int(math.floor(xsec)), ("%f" % round((xsec-math.floor(xsec)), rprec))[2:2+rprec], ideg, imn, int(math.floor(xsc)), ("%f" % round((xsc-math.floor(xsc)), dprec))[2:2+dprec])
        print outstr
        return outstr

def sixty(scalar):
    """
    Converts a decimal number to sexigesimal.

    Reverse of L{ten()} function.

    Example Usage:

       >>> sixty(23.445)
       [23, 26, 42.0]
       >>> sixty(-23.445)
       [-23, 26, 42.0]

    @param scalar: Decimal quantity
    @return: List of 3 floats of degrees, minutes, seconds in
        sexagesimal format. By default, a negative number is
        signified by making the first non-zero element of the
        output vector negative.

        """
    if scalar < 0.0:
        neg = True
    else:
        neg = False
    second = abs(3600.0*scalar)
    minute = abs(60.*scalar)
    degree = abs(scalar)
    result = []
    result.append(int(degree))
    result.append(int(minute-60.0*result[0]))
    result.append(second - 3600.*result[0] - 60.0*result[1])
    if neg:
        applied = False
        for i in range(3):
            if result[i] == 0.0:
                continue
            if not applied:
                result[i] = -result[i]
                applied = True
    return result


def ten(*args):
    """
    Converts a sexigesimal number to decimal.
    Inverse of sixty function.

    Example Usage:

       >>> ten('4:23:43')
       4.3952777777777783
       >>> ten(4, 23, 43)
       4.3952777777777783

    @param Inputs: Input can be a string of form '04:23:30' or
        can be 1 to 3 scalar floats which will be interpreted
        as hours/degrees, minutes, seconds.

    @return: Decimal (float) equivalent value of the input sexigesimal.
        A minus sign on any nonzero element of the input vector
        causes all the elements to be taken as < 0.
       """
    np = len(args)    
    vec = numpy.zeros((3,), dtype='float')
    if np == 1:
        if type(args[0]) == types.StringType:
            sp_args = map(float,args[0].split(':'))
            if len(sp_args)>=1 and len(sp_args)<=3:
                for i,a in enumerate(sp_args):
                    vec[i] = a
            else:
                print 'Incorrect format for input string'
                return
        else:
            vec[0] = args[0]
    elif np<1 or np>3:
        print "Number of arguments should be between 1 to 3"
        return
    else:
        for i,a in enumerate(args):
            vec[i] = a
    if len(numpy.where(vec<0)[0]):
        sign = -1.0
    else:
        sign = 1.0
    vec = abs(vec)
    facs  = [1., 60.0, 3600.0]
    decim = float(vec[0])
    if len(vec)>1:
        for i in range(1,len(vec)):
            decim += vec[i]/facs[i]
    return decim*sign

def radec(ra, dec, hours=False):
    """
    To convert RA and Dec  from decimal to sexigesimal units.

    The conversion is to sexigesimal hours for RA,  and sexigesimal 
    degrees for declination. If ra and dec are lists rather than scalars,
    all outputs are lists as well.

    @param ra: right ascension, scalar or vector, in DEGREES unless the
        HOURS keyword is set to True
    @param dec: declination in decimal DEGREES, scalar or vector, same
        number of elements as RA
    @param hours: Optional input. if set to True, then the input righ ascension
        should be specified in hours instead of degrees.
    @return: 6 outputs, either lists or scalars based on inputs
        ihr  - right ascension hours   (int)
        imin - right ascension minutes (int)
        xsec - right ascension seconds  (float)
        ideg - declination degrees (int)
        imn  - declination minutes (int)
        xsc  - declination seconds (float)
      """
    if type(ra) in (types.ListType, types.TupleType):
        vector = True
        if type(dec) not in (types.ListType, types.TupleType) and len(dec) != len(ra):
            print "ra and dec should have similar lengths"
            return None, None, None, None, None, None
    else:
        if type(dec) in (types.ListType, types.TupleType):
            print "ra and dec should have similar lengths"
            return None, None, None, None, None, None
        vector = False
    if vector:
        ihr = []
        imin = []
        xsec = []
        ideg = []
        imn = []
        xsc = []
        for i, r in enumerate(ra):
            if hours:
                r = r % 24
            else:
                r = (r % 360)/15
            h,m,s = sixty(r)
            ihr.append(h)
            imin.append(m)
            xsec.append(s)

            h,m,s = sixty(dec[i])
            ideg.append(h)
            imn.append(m)
            xsc.append(s)
    else:
        if hours:
            ra = ra % 24
        else:
            ra = (ra % 360)/15
        ihr, imin, xsec = sixty(ra)
        ideg, imn, xsc = sixty(dec)
    return ihr, imin, xsec, ideg, imn, xsc
