"""
Convert 4 dimensional FITS cube files 
produced by GILDAS to 3 dimensional
FITS cubes
"""
from astropy.io import fits as pyfits
from sculpt.idealpy.radio import transpose_cube

def convert_gildas_and_transpose(filename,
                                 output_filename=None):
    """
    Takes a 4d GILDAS file and converts it first 
    into a 3d xyv format cube, and then transposes
    to 3d vxy format and saves output
    """
    if output_filename is None:
        output_filename = filename.lower()
    hdulist = pyfits.open(filename)
    hdu = hdulist[0]
    header = hdu.header
    data = hdu.data
    data.shape = data.shape[1:]
    header['NAXIS'] = 3
    for attr in ('NAXIS4', 'CTYPE4', 'CRVAL4',
                 'CRPIX4', 'CROTA4', 'CDELT4'):
        del header[attr]
    hduout = transpose_cube(hdu, origin='xyv')
    hdulist = pyfits.HDUList([hduout])
    print "Writing converted FITS file to %s" % output_filename
    hdulist.writeto(output_filename)
    return
