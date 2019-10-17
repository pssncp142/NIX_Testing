from astropy.io import fits
import numpy as np
import sep
import os
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

class NIX_Base:

    full_path = ''
    f_name = ''
    header = ''

    def __init__(self, path, test_id=None, verbose=False):

        if verbose:
            print "Loading file ", path.split('/')[-1]
        self.f_name = path.split('/')[-1]
        self.full_path = path
        self.header = fits.getheader(path)
        if test_id is not None:
            self.test_id = test_id

    def getImage(self, hdu=None, mask=None, dark=None):

        hdul = fits.open(self.full_path)
        image = hdul[0].data*1. 
        # Handle exception where data format is (1, 2048, 2048)
        if image.shape[0] == 1:
            image = image[0,:,:]

        if dark is not None:
            image = image - dark.getImage()
        if mask is not None:
            image = image * mask.getImage()
        return image

    def plotImage(self, mask=None, dark=None):

        plt.imshow(self.getImage(mask=mask, dark=dark))
        plt.gca().invert_yaxis()
        plt.show()

class NIX_Spectra(NIX_Base):
    
    def _singleLineFit(self, xx, yy, window=10):

        from lmfit.models import GaussianModel, ConstantModel

        model = GaussianModel()+ConstantModel()
        pars = ConstantModel().make_params()
        pars += GaussianModel().guess(yy, x=xx)
        out = model.fit(yy, pars, x=xx)

        return out

    def _polFit(self, xx, yy, order):

        from lmfit.models import PolynomialModel

        model = PolynomialModel(order)
        pars = model.guess(yy, x=xx)
        out = model.fit(yy, pars, x=xx)
        
        return out

    def _continuumFit(self, method='polynomial', order=5):

        import pandas as pd

        line_data = pd.read_csv('data/ThAr_lines.csv', header=None)
        line_data = line_data.sort_values(4, ascending=False)
        line_data = np.array(line_data[2])

        self.spectraFiltered = self.spectra1d
        self.waveFiltered = self.wave

        for i in range(50):
            ndx = np.where(np.abs(self.waveFiltered - line_data[i]) < 5)[0]
            self.spectraFiltered = np.delete(self.spectraFiltered, ndx)
            self.waveFiltered = np.delete(self.waveFiltered, ndx)
            #plt.plot([line_data[i], line_data[i]], [0, 1.2e6], 'r--')

        if method == "polynomial":
            from lmfit.models import PolynomialModel
            model = PolynomialModel(order)
            pars = model.guess(self.spectraFiltered, x=self.waveFiltered)
            out = model.fit(self.spectraFiltered, pars, x=self.waveFiltered)
            self.continuum = np.poly1d([out.best_values['c%d' % i] for i in range(order+1)[::-1]])
        elif method == "spline":
            from scipy.interpolate import UnivariateSpline
            self.continuum = UnivariateSpline(self.waveFiltered, self.spectraFiltered, 
                w=1/self.spectraFiltered, k=order, s=6)

    def getOrderTracing(self, bright_line=None, mask=None, dark=None):

        window = 10

        if bright_line is not None:
            image = self.getImage(mask=mask, dark=dark)
            xys = bright_line['xys']
            centers = []
            plt.figure(figsize=(4*2, 4))
            for i in range(2):
                yy = image[xys[i][1]-window:xys[i][1]+window, xys[i][0]]
                xx = np.arange(xys[i][1]-window,xys[i][1]+window)
                out = self._singleLineFit(xx, yy, window=window) 
                centers.append(out.best_values['center'])

                plt.subplot(1,2,i+1)
                plt.plot(xx, yy, 'k')
                plt.plot(xx, out.best_fit, 'r')

            plt.show()
            
        self.A = (centers[1]-centers[0])/(xys[1][0] - xys[0][0])

    def getSpectra1D(self, bright_line=None, mask=None, dark=None):

        if not hasattr(self, 'A'):
            self.getOrderTracing(bright_line=bright_line, mask=mask, dark=dark)

        image = self.getImage(mask=mask, dark=dark)

        X = np.arange(2048)
        XX, YY = np.meshgrid(X, X)

        B = YY-self.A*XX
        B = np.floor(B).astype(int)

        spec1d = np.zeros(2048)

        for i in range(2048):
            for j in range(600, 1600):
                spec1d[B[i,j]] += image[i,j]

        plt.figure(figsize=(16, 4))

        plt.plot(X, spec1d)
        plt.show()

        self.spectra1d = spec1d

    def subtractContinuum(self, method='polynomial', order=5):

        plt.figure(figsize=(16, 4))
        out = self._continuumFit(method=method, order=order)

        plt.title('Continuum Model')
        plt.plot(self.wave, self.spectra1d)
        plt.plot(self.wave, self.continuum(self.wave))
        plt.show()
 
        plt.figure(figsize=(16, 4))
        plt.title('Continuum Subtracted')
        plt.plot(self.wave, self.spectra1d - self.continuum(self.wave))
        plt.show()

    def calibrate(self, ref_lines=None, mask=None, dark=None, order=2):
        
        window = 10

        if ref_lines is not None:
            image = self.getImage(mask=mask, dark=dark)
            sz = len(ref_lines)
            plt.figure(figsize=(4*sz, 4))
            xs = [ref_line['x'] for ref_line in ref_lines]
            waves = [ref_line['wave'] for ref_line in ref_lines]
            for i in range(sz):
                yy = self.spectra1d[xs[i]-window:xs[i]+window]
                xx = np.arange(xs[i]-window,xs[i]+window)
                out = self._singleLineFit(xx, yy, window=window) 
                xs[i] = out.best_values['center']

                plt.subplot(1,sz,i+1)
                plt.plot(xx, yy, 'k')
                plt.plot(xx, out.best_fit, 'r')

            plt.show()
        
        out = self._polFit(xs, waves, order=order)

        print out.fit_report()

        xx = np.arange(2048)
        
        waveSol = np.poly1d(out.best_values.values())
        self.wave = waveSol(xx)

        plt.figure(figsize=(16, 4))
        plt.plot(self.wave, self.spectra1d)
        plt.show()

    def extract(self, ref_lines, mask=None, dark=None):
        
        spectra2d = self.getImage(mask=mask, dark=dark)
        self.plotImage(mask=mask, dark=dark)




class NIX_Image(NIX_Base):

    def getMedian(self, mask=None):

        if not hasattr(self, 'median'):
            image = self.getImage(mask=mask)
            self.median = np.median(image)
        return self.median

    def getObjects(self, mask=None, search=None):
        
        if not hasattr(self, 'objects'):
            image = self.getImage(mask=mask)
            self.objects = sep.extract(image, 500, minarea=20)

        if search is None:
            return self.objects
        else:
            # Caution. This will return first found object within the search radius
            for object in self.objects:
                if ((object['x'] - search['x'])**2 + (object['y'] - search['y'])**2 < search['r']**2):
                    return object
      
    def plotObjects(self, mask=None):

        fig, ax = plt.subplots()

        ax.imshow(self.getImage(mask=mask))
        ax.invert_yaxis()

        for object in self.objects:
        
            e = Ellipse(xy=(object['x'], object['y']),
                        width=6*object['a'],
                        height=6*object['b'],
                        angle=object['theta'] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        
        plt.ylim([0, 2048])
        plt.xlim([0, 2048])
        plt.show()


class NIX_Image_List:

    def __init__(self, config):

        self.NIX_Files = []
        self.config = config
        self.loadFiles()
        self.Filtered = self.NIX_Files

    def __getitem__(self, ndx):

        return self.Filtered[ndx]

    def __len__(self):

        return len(self.Filtered)

    def loadFiles(self):

        for test_id in self.config['test_ids']:
            full_path = self.config['data_dir'] + '/'+ ('%s-%s' % tuple(test_id.split('-')[:2])) + '/' + test_id + '/'
            files = os.listdir(full_path)
            files.sort()
    
            for f_name in files:
                if f_name.endswith('fits'):
                    self.NIX_Files.append(NIX_Image(full_path + f_name, test_id=test_id))

    def printTable(self, keywords, tbl_fmt):

        table_header =  ['NDX', 'TEST_ID', 'FILENAME']
        table_header.extend([strip_prefix_from_keyword(keyword) for keyword in keywords])
        print ('%4s' + tbl_fmt) % tuple(table_header)
        
        ctr = 0 
        for File in self.NIX_Files:
            out_list = [ctr, File.test_id, File.f_name]
            out_list.extend([File.header[keyword] for keyword in keywords])
            print ('%04d' + tbl_fmt) % tuple(out_list) 
            ctr += 1

    def printFiltered(self, keywords, tbl_fmt):

        table_header =  ['NDX', 'TEST_ID', 'FILENAME']
        table_header.extend([strip_prefix_from_keyword(keyword) for keyword in keywords])
        print ('%4s' + tbl_fmt) % tuple(table_header)
        
        ctr = 0 
        for File in self.Filtered:
            out_list = [ctr, File.test_id, File.f_name]
            out_list.extend([File.header[keyword] for keyword in keywords])
            print ('%04d' + tbl_fmt) % tuple(out_list) 
            ctr += 1


    # revisit later on for filtering data
    def filter(self, val):
        if isinstance(val, slice):
            self.Filtered = self.NIX_Files[val]

    def medianCombine(self, out=None):

        sz = len(self.Filtered)
        ims = np.zeros([2048, 2048, sz])

        for i in range(len(self.Filtered)):
            ims[:,:,i] = self.Filtered[i].getImage()

        if out is not None:
            fits.PrimaryHDU(np.median(ims, axis=2)).writeto(out, overwrite=True)
        
        return np.median(ims, axis=2)

    # very crude BP mask generation to be revisited later on
    def doBadPixelMask(self, thresh=20, out=None):
        
        im = self.medianCombine()

        im[np.where(im < thresh)] = 1
        im[np.where(im >= thresh)] = 0

        if out is not None:
            fits.PrimaryHDU(im).writeto(out, overwrite=True)
        
        return im

    def getMedian(self, mask=None):

        return [NI.getMedian(mask=mask) for NI in self.Filtered]

    def getObjects(self, mask=None, search=None, sepdict=None):
        
        if sepdict is not None:
            return [NI.getObjects(mask=mask, search=search)[sepdict] for NI in self.Filtered]
        else:
            return [NI.getObjects(mask=mask, search=search) for NI in self.Filtered]

    def getHeaderValue(self, keyword):

        return [NI.header['MJD-OBS'] for NI in self.Filtered]

# Some utility functions

def strip_prefix_from_keyword(keyword):

    if keyword.startswith('HIERARCH ESO'):
        keyword = keyword[13:]

    keyword = '_'.join(keyword.split(' '))
    return keyword





