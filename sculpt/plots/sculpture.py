from .sculpture_fitsfigure import SculptureFITSFigure
import matplotlib.pyplot as mpl
import numpy

class Sculpture(SculptureFITSFigure):
    """
    The Sculpture Class is a top level class 
    for holding FITS images and cubes and provides
    an interactive platform for data visualization
    """
    def __init__(self, data, hdu=0, figure=None,
                 subplot=(1,1,1), downsample=False,
                 north=False, convention=None,
                 dimensions=[0, 1], slices=[],
                 gauss_width=2, extra_hdus=[],
                 extra_hdu_titles=[],
                 blue_windows=[], red_windows=[],
                 sigma_crit=3.0, source_distance=140.0,
                 avgfigure=None, pvangle=None,
                 **kwargs):
        SculptureFITSFigure.__init__(self, data, hdu, figure=figure,
                                     subplot=subplot, downsample=downsample,
                                     north=north, convention=convention,
                                     dimensions=dimensions, slices=slices,
                                     **kwargs)
        self.extra_hdus = extra_hdus
        self.extra_hdu_titles = extra_hdu_titles
        self.blue_windows = blue_windows
        self.red_windows = red_windows
        self.sigma_crit = sigma_crit
        self.avgfigure = avgfigure
        self.pvangle = pvangle
        self.polygon_vert_blue = []
        self.polygon_vert_red = []
        self.blue_polygon_active = False
        self.red_polygon_active = False
        self.source_distance = source_distance
        self.specfigure = mpl.figure(10)
        self.keystrokes = {#key: ('help-text', 'help-prompt')
            'b': ('Start Box-car average',
                  'Click on other edge of box and press e'),
            'e': ('End Box-car average',
                  'First click b to start box'),
            'B': ('Blue-wing Sigma-weighted plot',
                  'Click on the other edge of box and press E'),
            'R': ('Red-wing Sigma-weighted plot',
                  'Click on the other edge of box and press E'),
            'E': ('End blue or red sigma-weighted plot',
                  'First click B or R to start box'),
            'S': ('Popup Spectrum',
                  'Pops up spectrum at this position'),
            'p': ('Position-Velocity Start',
                  'PV Cut: Click on other end of line and press v'),
            'v': ('Position-Velocity End',
                  'First start p to start position-velocity cut'),
            'P': ('Position-Velocity by Angle Start',
                  'PV Cut: Click A to enter angle for PV Cut'),
            'A': ('Position-Velocity by Angle End',
                  'First click P to start position-velocity cut by angle'),
            'd': ('Distance measurement start',
                  'Click D to mark end of line for linear distance'),
            'D': ('Distance measurement end',
                  'First click d to start distance measurement'),
            'G': ('Start definition of blue-shifted polygon',
                  'Keep pressing G to define polygon vertices and press c to close'),
            'g': ('Start definition of red-shifted polygon',
                  'Keep pressing G to define polygon vertices and press c to close'),
            'c': ('Close blue or red-shifted polygon',
                  'First press G or g to start blue or red-shifted polygons')
            }
        self.start_events()
        
        
    def print_keystroke_help(self):
        for key, (help_text, help_prompt) in self.keystrokes.items():
            print "%s: %s" % (key, help_text)

    def on_keypress(self, event):
        if not event.inaxes: return
        if event.key not in self.keystrokes.keys():
            print "Key %s not a recognized keystroke" % event.key
            self.print_keystroke_help()
        print event.xdata, event.ydata
        if event.key == 'b':
            #starts a box for box-averaged plot
            self.box = []
            self.box.append((event.xdata, event.ydata))
            print self.keystrokes[event.key][1]
            self.avgfigure = None
        if event.key == 'e':
            if not self.box:
                print self.keystrokes[event.key][1]
                return
            else:
                self.box.append((event.xdata, event.ydata))
                self.draw_box_and_average()
        if event.key == 'B':
            #starts a box for box-averaged blue sigma-weighted plot
            self.box = []
            self.box.append((event.xdata, event.ydata))
            print self.keystrokes[event.key][1]
            self.blue = True
            self.red = False
        if event.key == 'R':
            #starts a box for box-averaged blue sigma-weighted plot
            self.box = []
            self.box.append((event.xdata, event.ydata))
            print self.keystrokes[event.key][1]
            self.red = True
            self.blue = False
        if event.key == 'E':
            if not self.box:
                print self.keystrokes[event.key][1]
                return
            else:
                self.box.append((event.xdata, event.ydata))
                self.draw_box_and_sigma_average()
        if event.key == 'S':
            #popup a spectrum
            self.draw_spectrum(event.xdata, event.ydata)
        if event.key == 'p':
            #position-velocity cut
            self.pvcut = []
            self.pvcut.append((event.xdata, event.ydata))
            print self.keystrokes[event.key][1]
            self.pvfigure = None
        if event.key == 'v':
            if not self.pvcut:
                print self.keystrokes[event.key][1]
                return
            else:
                self.pvcut.append((event.xdata, event.ydata))
                self.make_posvel_cut()
        if event.key == 'P':
            #print ra, dec for location
            print event.xdata, event.ydata
            r = self._wcs.wcs_pix2sky([[event.xdata,event.ydata]], 1)
            ra = r[0][0]
            dec = r[0][1]
            print "RA = %s, Dec= %s" % (ra, dec)
            print "RA = %s, Dec = %s" % (sixty(ra*24./360.), sixty(dec))
        if event.key == 'A':
            print event.xdata, event.ydata
            if self.pvangle is None:
                angle = get_input_with_default("Enter angle (reckoned from +ve X-axis) in degrees for PV cut",
                                               90., typeconvert=float)
            else:
                angle = self.pvangle
            self.make_posvel_cut_angle(event.xdata, event.ydata,
                                       angle)
        if event.key == 'd':
            #starts distance cut
            self.distance = []
            self.distance.append((event.xdata, event.ydata))
            print self.keystrokes[event.key][1]
        if event.key == "D":
            if not self.distance:
                print self.keystrokes[event.key][1]
                return
            else:
                self.distance.append((event.xdata, event.ydata))
                self.find_dist()
        if event.key == "G":
            #start blue polygon definition
            if not self.blue_polygon_active:
                self.polygon_vert_blue = []
                self.blue_polygon_active = True
            print event.xdata, event.ydata
            print self.keystrokes[event.key][1]
            r = self._wcs.wcs_pix2sky([[event.xdata,event.ydata]], 1)
            ra = r[0][0]
            dec = r[0][1]
            self.polygon_vert_blue.append([ra, dec])
            if len(self.polygon_vert_blue) > 2:
                self.show_lines([numpy.array(self.polygon_vert_blue[-2:]).T, ],
                                color='blue', linewidth=0.5)
            self.blue_polygon = True
            self.red_polygon = False
        if event.key == "g":
            #start red polygon definition
            if not self.red_polygon_active:
                self.polygon_vert_red = []
                self.red_polygon_active = True                
            print event.xdata, event.ydata
            print self.keystrokes[event.key][1]
            r = self._wcs.wcs_pix2sky([[event.xdata,event.ydata]], 1)
            ra = r[0][0]
            dec = r[0][1]
            self.polygon_vert_red.append([ra, dec])
            if len(self.polygon_vert_red) > 2:
                self.show_lines([numpy.array(self.polygon_vert_red[-2:]).T, ],
                                color='red', linewidth=0.5)
            self.blue_polygon = False
            self.red_polygon = True
        if event.key == "c":
            #close polygon definition
            if not self.polygon_vert_blue:
                if not self.polygon_vert_red:
                    print self.keystrokes[event.key][1]
                    return
            if self.blue_polygon:
                self.polygon_vert_blue.append(self.polygon_vert_blue[0])
                self.pvert_blue = [numpy.array(self.polygon_vert_blue),]
                print self.polygon_vert_blue
                print self.pvert_blue
                #self.show_polygons(self.pvert_blue, edgecolor='b')
                self.blue_polygon_active = False
                self.polygon_average(self.pvert_blue, bluered='blue')
            else:
                self.polygon_vert_red.append(self.polygon_vert_red[0])
                self.pvert_red = [numpy.array(self.polygon_vert_red),]
                print self.polygon_vert_red
                print self.pvert_red                
                #self.show_polygons(self.pvert_red, edgecolor='r')
                self.red_polygon_active = False
                self.polygon_average(self.pvert_red, bluered='red')

    def draw_spectrum(self, xdata, ydata):
        if not self.extra_hdus:
            print "Class was not instantiated with any extra HDU cubes"
            return
        self.specfigure.clf()
        specax = self.specfigure.add_subplot(111)
        for i, hdu in enumerate(self.extra_hdus):
            spec = hdu.data[ydata, xdata, :]
            vel = getaxes(hdu.header, 1)
            if not self.extra_hdu_titles:
                title = 'HDU 0'
            else:
                title = self.extra_hdu_titles[i]
            specax.plot(vel, co12, linestyle='steps-middle', label=title)
        specax.set_title('Spectrum at %.1f, %.1f' % (xdata, ydata))
        specax.legend(loc='best')
        self.specfigure.canvas.draw()    


    def start_events(self):
        self.cid = self._figure.canvas.mpl_connect('key_press_event', self.on_keypress)
        
