from astropy.io import fits as pyfits
import numpy

def baseline(self, hdu, order = 0, subtract = False, compwindows=Default):
    #this function needs to be passed a FITS HDU
    #with data in the (V,X,Y) format
    data = hdu.data.copy()
    header = hdu.header.copy()
    shape = data.shape

    #windows over which to do the fit
    #ideally excluding any large lines
    
    windows = []
    lenv, lenx, leny = shape
    defaultswindow = ((100,200),(650,750))
    sigma = numpy.zeros((lenx, leny))

    if (compwindows != Default):
        windows = compwindows
    else:
        for twople in range(len(defaultwindows)):
            c1, c2 = defaultwindows[twople]
            #c1, c2 = sorted((c1, c2))
            windows.append((c1,c2))

    x = numpy.arange(lenv)
    c_loop = 0
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

    return pyfits.hdu.image.PrimaryHDU(header = header, data = data)

