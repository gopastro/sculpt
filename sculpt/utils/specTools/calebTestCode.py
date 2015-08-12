__author__ = 'chw3k5'
import pyfits
import numpy
import getpass
from matplotlib import pyplot as plt
from astropy.modeling import models, fitting

def getSpecDataCube(filename,):
    if getpass.getuser() == 'chw3k5':
        dirname = '/Users/chw3k5/Documents/Grad_School/supercamG30data/'
        filename = dirname+filename
    HDU = pyfits.open(filename)
    fitsData = HDU[0]
    specDataCube = fitsData.data[0,:,:,:] # The 0 is the polarization axis
    return specDataCube, HDU

def findHighIntensitySpec(specDataCube, plotData=False):
    val_Max = 0
    val_galLat = -1
    val_galLon = -1
    mean_Max = 0
    mean_galLat = -1
    mean_galLon = -1
    len_galLon = len(specDataCube[0,0,:])
    len_galLat = len(specDataCube[0,:,0])
    len_spec = len(specDataCube[:,0,0])
    for galLatIndex in range(len_galLat):
        print galLatIndex
        for galLonIndex in range(len_galLon):
            testSpec=specDataCube[:,galLatIndex,galLonIndex]
            test_val = max(testSpec)
            if val_Max < test_val:
                val_Max = test_val
                val_galLat = galLatIndex
                val_galLon = galLonIndex
            test_mean = numpy.mean(testSpec)
            if mean_Max < test_mean:
                mean_Max = test_mean
                mean_galLat = galLatIndex
                mean_galLon = galLonIndex
    if plotData:
        plt.plot(specDataCube[:,val_galLat,val_galLon])
        plt.show()
        plt.plot(specDataCube[:,mean_galLat,mean_galLon])
        plt.show()
    return (val_galLat,mean_galLon), (val_galLat,mean_galLon)




if __name__ == "__main__":
    filename = 'G30_Map_a.fits'
    specDataCube, HDU = getSpecDataCube(filename=filename)
    testSpec = specDataCube[:,342,122]
    len_spec = len(testSpec)
    x = range(len_spec)
    g1 = models.Gaussian1D(1, 270, 4)
    g2 = models.Gaussian1D(2.5, 300, 3)
    gg_init = g1+g2
    fitter = fitting.SLSQPLSQFitter()
    gg_fit = fitter(gg_init,x,testSpec)

    plt.plot(x, testSpec, 'ko')
    plt.plot(x, gg_fit(x), 'r-', lw=2)
    plt.show()
