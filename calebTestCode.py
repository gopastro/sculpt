__author__ = 'chw3k5'
import pyfits, numpy, socket, platform, os, getpass
from matplotlib import pyplot as plt

if getpass.getuser() == 'chw3k5':
    dirname = '/Users/chw3k5/Documents/Grad_School/supercamG30data/'
else:
    dirname = ''

filename = dirname+'G30_Map_a.fits'
HDU = pyfits.open(filename)
temp = HDU[0]
data = temp.data
testSpec=data[0,:,0,0]
#
# print testSpec

plt.plot(testSpec)
plt.show()


print dirname

