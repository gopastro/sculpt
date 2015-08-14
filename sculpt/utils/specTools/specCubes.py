__author__ = 'chw3k5'

# This file is currently mostly empty, but will be used to handle cubes of spectral data use the analysis tools
# in spec only functions


if __name__ == "__main__":
    # Getting data here and selecting a test spectrum
    filename = 'G30_Map_a.fits'
    specDataCube, HDU = getSpecDataCube(filename=filename)
    testSpec = specDataCube[:,342,122]

    fractionalRimprovment = 0.005
    showMaximaPlot = False
    maximaConvSigma = 5.
    maximaConvRadius = 20

    showFitterPlots = True