__author__ = 'chw3k5'
import pyfits
import numpy
import getpass
from matplotlib import pyplot as plt
from astropy.modeling import models, fitting

def getSpecDataCube(filename,):
    if getpass.getuser() == 'chw3k5':
        dirname = '/Users/chw3k5/Documents/Grad_School/supercamG30data/'

    elif getpass.getuser() == 'theodore':
	dirname = '/home/theodore/astroenv/supercamdata/supercamG30data/'
    else:
        dirname = ''
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

def test_for_multiple_peaks(specdatacube, pos_x, pos_y, noise_threshhold, rmsnoise = None, test=False):
    # this function initially draws a single gaussian to fit the peak
    # if the residuals are then anti-symmetric about this new peak
    # we then re-fit the original spectra with two gaussians

    #probably only works with spectra with a single major peak
    if test:
        x = numpy.arange(850)
        peak1 = 6.5 * numpy.exp(-numpy.power(x - 284, 2.) / (2 * numpy.power(10, 2.)))
        peak2 = 10.0 * numpy.exp(-numpy.power(x - 270, 2.) / (2 * numpy.power(4, 2.)))
        spectra = peak1 + peak2 + numpy.random.normal(0.0, size = len(x)) * 1.2
    else:
        spectra = specdatacube[:,pos_x,pos_y]
        x = numpy.arange(len(spectra))

    guess = numpy.where(spectra == spectra.max())[0]
    std = 5.0  #just a guess for now
    init_model = models.Gaussian1D(spectra.max(),guess,std)
    fitter_init = fitting.SLSQPLSQFitter()
    fit_init = fitter_init(init_model,x,spectra)
    residuals_init = spectra - fit_init(x)
    #plt.plot(x, spectra)
    #plt.plot(x, fit_init(x))
    #plt.plot(x, residuals_init)
    if (rmsnoise == None):
        rmsnoise = residuals_init.std()
    before = numpy.mean(residuals_init[guess-6:guess])
    b_sign = numpy.sign(before)
    after = numpy.mean(residuals_init[guess:guess+6])
    a_sign = numpy.sign(after)
    print (rmsnoise,before, after)
    if ((numpy.sign(before) != numpy.sign(after))):
        if ((abs(before) >= rmsnoise) or (abs(after) >= rmsnoise)):
            print "there might be two lines here"
            #time to fit the residuals
            if (b_sign == 1):
                resid_model = models.Gaussian1D(spectra.max() / 10.0, guess - 3, std) -  models.Gaussian1D(spectra.max() / 10.0, guess + 3, std)
            else:
                resid_model = models.Gaussian1D(spectra.max() / 10.0, guess + 3, std) -  models.Gaussian1D(spectra.max() / 10.0, guess - 3, std)
            fitter_resid = fitting.SLSQPLSQFitter()
            fit_resid = fitter_resid(resid_model, x, residuals_init)
            #now that we've fit the residuals, we'll use this as a seed for the final 2-line fit
            final_gauss1 = models.Gaussian1D(spectra.max() / 2.0, fit_resid[0].mean[0][0], fit_resid[0].stddev[0])
            final_gauss2 = models.Gaussian1D(spectra.max() / 2.0, fit_resid[1].mean[0][0], fit_resid[1].stddev[0])
            final_model = final_gauss1 + final_gauss2
            fitter_final = fitting.SLSQPLSQFitter()
            fit_final = fitter_final(final_model, x, spectra)
            plt.figure(1)
            plt.subplot(211)
            plt.plot(x, spectra, label = "Spectra")
            plt.plot(x, fit_init(x), label = "1 Line Fit")
            plt.legend()
            plt.subplot(212)
            plt.plot(x, spectra, label = "Spectra")
            plt.plot(x, fit_final(x), label = "2 Line Fit")
            plt.plot(x, fit_final[0](x), label = "Gaussian 1")
            plt.plot(x, fit_final[1](x), label = "Gaussian 2")
    else:
        plt.plot(x, spectra, label = "Spectra")
        plt.plot(x, fit_init(x), label = "1 Line Fit")
    plt.legend()
if __name__ == "__main__":
    filename = 'G30_Map_a.fits'
    specdatacube, HDU = getSpecDataCube(filename=filename)
    test_for_multiple_peaks(specdatacube, 370, 122, 3.0)
    plt.show()
