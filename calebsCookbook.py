import sys, os, time, numpy
from matplotlib import pyplot as plt
# Taurus data
from astropy.io import fits as pyfits
from sculpt.idealpy.radio import momentcube
from aplpy import FITSFigure
conversionFolder = '/Users/chw3k5/Dissertation/Combined Kappa/code/sculpt/sculpt/scripts/'
sys.path.append(conversionFolder)
from convert_4dim_3dim_fits import convert_gildas_and_transpose
# from sculpt.idealpy.radio.transpose import transpose



parentFolder = '/Users/chw3k5/Dissertation/Combined Kappa/code/sculpt/AmherstData/'
outputFileName = parentFolder + 'blah.fits'

taurasInputFileName = parentFolder + 'taurustest/iras04361+2547_12.fits'

supercamInputFileName = parentFolder + 'supercamG30data/' + 'G30_Map_c.fits'
supercamConvertedFileName = parentFolder + 'supercamG30data/' + 'converted_G30_Map_c.fits'

sReduceInputFileName = '/Users/chw3k5/Dissertation/Combined Kappa/code/sReduce/supercam_testdata/test.fits'

convertSupercamData = False # this only needs to be done once
do_taurasData = False
do_superCamData = False
do_sReduce = True

def lookAtHDU(hdulist, outputFileName, v1=0, v2=10, channelMode=False, coAddSpectra=False,
              lowVelCut = 0, highVelCut = 0, plotInKmPerS = True):
    hdu = hdulist[0]
    workingHeader = hdu.header
    img = momentcube(hdu, v1, v2, moment=0, chan=channelMode)
    workingImageHeader = img.header
    imglist = pyfits.PrimaryHDU(data=img.data, header=img.header)
    # write an output file
    hdulist2 = pyfits.hdu
    hdulist2 = pyfits.HDUList([imglist])
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    hdulist2.writeto(outputFileName)

    if coAddSpectra:
        (axis3Len, axis2Len, specLen) = numpy.shape(hdu.data)
        specLen -= lowVelCut
        specLen -= highVelCut
        meanSpectra = numpy.zeros((specLen))
        for axis3Index in range(axis3Len):
            for axis2Index in range(axis2Len):
                meanSpectra +=  hdu.data[axis3Index, axis2Index, highVelCut:specLen+highVelCut]
        meanSpectra = meanSpectra / float(axis3Len * axis2Len)
        vAxisStart = hdu.header['CRVAL1'] - ((hdu.header['CRPIX1'] - 1.0) * hdu.header['CDELT1'])
        vAxisEnd = hdu.header['CRVAL1'] + (((specLen + 1) - hdu.header['CRPIX1'] ) * hdu.header['CDELT1'])
        vAxisStep = hdu.header['CDELT1']
        velocityAxis = numpy.arange(vAxisStart, vAxisEnd, vAxisStep)
        if plotInKmPerS:
            velocityAxis = velocityAxis / 1000.0
            plt.xlabel("Relative Velocity (km/s)")
        else:
            plt.xlabel("Relative Velocity (m/s)")
        velocityAxis = list(velocityAxis)
        while specLen < len(velocityAxis):
            velocityAxis.pop()
        plt.plot(meanSpectra)
        # plt.plot(velocityAxis, meanSpectra)
        plt.title("max value:" + str(max(meanSpectra)) + "   Min value:" + str(min(meanSpectra)))
        plt.show()


    fig = FITSFigure(img)
    fig.show_colorscale( )
    # time.sleep(20)
    # raw_input('press enter to continue')
    return

# Taurus data
if do_taurasData:
    hdulist = pyfits.open(taurasInputFileName)
    lookAtHDU(hdulist=hdulist, outputFileName=outputFileName,v1=0, v2=10, coAddSpectra=True)


# # Supercam Data
# if do_superCamData:
#     if convertSupercamData:
#         convert_gildas_and_transpose(supercamInputFileName, supercamConvertedFileName)
#     hdulist = pyfits.open(supercamConvertedFileName)
#     lookAtHDU(hdulist=hdulist, outputFileName=outputFileName, v1=100, v2=200,  lowVelCut = -5, highVelCut = 50,
#               coAddSpectra=True)

# sReduce
if do_sReduce:
    hdulist = pyfits.open(sReduceInputFileName)
    lookAtHDU(hdulist=hdulist, outputFileName=outputFileName, v1=520, v2=550,
              lowVelCut = 0, highVelCut = 0,channelMode=True, coAddSpectra=True)


