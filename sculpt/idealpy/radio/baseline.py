from astropy.io import fits as pyfits
from sculpt.idealpy.fits import sxpar, sxdelpar, sxaddpar, getaxes, sxaddhist
from sculpt.utils import SculptArgumentError
import numpy

def baseline(hdu, chan, windows, order = 0, subtract = True, returnrms = True, kms = True):
    #this function needs to be passed a FITS HDU
    #with data in the (V,X,Y) format

    if isinstance(hdu, pyfits.hdu.image.PrimaryHDU):
        #get data and header from the hdu
	header = hdu.header
        data = hdu.data

    elif isinstance(hdu, numpy.ndarray):
        if header is None or not isinstance(header, pyfits.header.Header):
            raise SculptArgumentError('header', "Since you passed in data that is a numpy array, set header to a pyfits header type")
        data = hdu

    shape = data.shape

    #windows over which to do the test for sigma values
    #basically excluse major lines you observe
    
    lenv, lenx, leny = shape
    #defaultswindow = ((100,200),(650,750))
    sigma = numpy.zeros((lenx, leny))

    x = numpy.arange(lenv)
    c_loop = 0

# this commented section is the makermsimage way of selecting them
# if someone wants to make a change regarding one being better, do it

#    velax = getaxes(header, 1)
#    indices = numpy.zeros(velax.shape, dtype=bool)
#    for v1, v2 in window:
#        if chan:
#            indices[v1:v2] = True
#        else:
#            ind = numpy.logical_and(velax>= v1, velax<= v2)
#            indices = numpy.logical_or(indices, ind)
#    indices = numpy.logical_not(indices) # for indices to include in std calc

    for win in windows:  
        if (len(win) != 2):
            print "Each element in the window list must be a 2-tuple"
        c1, c2 = win
        c1, c2 = sorted((c1, c2))
        ind = numpy.logical_and(x>=c1, x<=c2)
        if (c_loop == 0):
            final_ind = numpy.logical_or(ind,ind)
        else:
            final_ind = numpy.logical_or(final_ind,ind)

    x_windows = x[final_ind]
    for ix in range(lenx):
        for iy in range(leny):
            spectra = data[:,ix,iy]
            spec_windows = spectra[final_ind]

            p = numpy.polyfit(x_windows,spec_windows,order)
            spec_windows -= numpy.polyval(p,len(x_windows))

            sigma[ix, iy] = spec_windows.std()
            if (subtract):
                data[:,ix,iy] -= numpy.polyval(p,lenv)
                
    # this is the original input - with data reduced as needed
    hdu_orig =  pyfits.hdu.image.PrimaryHDU(header = header, data = data)
    
    if returnrms:
        #the following grabs relevant information from the original header
        #and reproduces it in the RMSMap header shifted to account for
        #the different shape of the data
        crpix1 = sxpar(header,"CRPIX1")
        crval1 = sxpar(header,"CRVAL1")
        cdelt1 = sxpar(header,"CDELT1")
        if kms:
            #convert velocity to km/s
            crval1 = crval1/1000.
            cdelt1 = cdelt1/1000.
        ctype1 = sxpar(header,"CTYPE1")
        crpix2 = sxpar(header,"CRPIX2")
        crval2 = sxpar(header,"CRVAL2")
        cdelt2 = sxpar(header,"CDELT2")
        ctype2 = sxpar(header,"CTYPE2")
        crpix3 = sxpar(header,"CRPIX3")
        crval3 = sxpar(header,"CRVAL3")
        cdelt3 = sxpar(header,"CDELT3")
        ctype3 = sxpar(header,"CTYPE3")
        nv = sxpar(header,"NAXIS1")
        nx = sxpar(header,"NAXIS2")
        ny = sxpar(header,"NAXIS3")
        blank = sxpar(header,"BLANK")

        hnew = header.copy()

        sxaddpar(hnew, "CRVAL1", crval2, comment="DEGREES")
        sxaddpar(hnew, "CRPIX1", crpix2)
        sxaddpar(hnew, "CDELT1", cdelt2, comment="DEGREES")
        sxaddpar(hnew, "CTYPE1", ctype2)
        sxaddpar(hnew, "CRVAL2", crval3, comment="DEGREES")
        sxaddpar(hnew, "CRPIX2", crpix3)
        sxaddpar(hnew, "CDELT2", cdelt3, comment="DEGREES")
        sxaddpar(hnew, "CTYPE2", ctype3)
        sxaddpar(hnew, "NAXIS", 2)
        sxaddpar(hnew, "NAXIS1", nx)
        sxaddpar(hnew, "NAXIS2", ny)
        sxaddpar(hnew, "NAXIS3", 1)
        sxaddpar(hnew, "NAXIS4", 1)

        if chan:
            vorc = 'CHANNEL'
        else:
            vorc = 'VELOCITY'

        sxaddhist(hnew, "WINDOW : %s; Window %s LIMITS" % (repr(windows), vorc))
        #sxaddpar(hnew, "BUNIT", units, "Units")
        sxdelpar(hnew, "CRVAL3")
        sxdelpar(hnew, "CRPIX3")
        sxdelpar(hnew, "CDELT3")
        sxdelpar(hnew, "CTYPE3")
        sxdelpar(hnew, "NAXIS3")
        sxdelpar(hnew, "NAXIS4")
        
        hdu_rms =  pyfits.hdu.image.PrimaryHDU(data=sigma, header=hnew)
        return (hdu_orig, hdu_rms)
    else:
        return hdu_orig

if __name__ == '__main__':
    print "test - you've activated the 'if __name__ == '__main__':' clause"
