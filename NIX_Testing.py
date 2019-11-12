from astropy.io import fits
from astropy.visualization import ImageNormalize, PercentileInterval
import numpy as np
import sep
import os
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import time
from IPython.display import clear_output
from scipy.stats import linregress
import multiprocessing as mp

class NIX_Base(object):

    __slots__ = ('full_path', 'f_name', 'header', 'test_id')

    def __init__(self, path, test_id=None, verbose=False):

        if verbose:
            print "Loading file ", path.split('/')[-1]
        self.f_name = path.split('/')[-1]
        self.full_path = path
        #self.header = fits.getheader(path)
        if test_id is not None:
            self.test_id = test_id
        else:
            self.test_id = None

    def getImage(self, hdu=None, mask=None, dark=None, linearize=None, gain=None):

        hdul = fits.open(self.full_path)
        image = hdul[0].data*1. 
        hdul.close()
        # Handle exception where data format is (1, 2048, 2048)
        if image.shape[0] == 1:
            image = image[0,:,:]

        if dark is not None:
            image = image - dark.getImage()
        if linearize is not None:
            image = linearize(image)
        if mask is not None:
            image = image * mask.getImage()
        if gain is not None:
            image *= gain
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

    def _multiLineFitFromTable(self, nlines=20):

        import pandas as pd
        from lmfit.models import GaussianModel, LorentzianModel

        line_data = pd.read_csv('data/ThAr_lines.csv', header=None)
        line_data = line_data.sort_values(4, ascending=False)
        line_data = np.array(line_data[2])

        line_data = line_data[np.where((line_data > 3090) & (line_data < 3900))[0]]

        plt.figure(figsize=(16,4))

        model_l = []

        for i in range(nlines):

            model = GaussianModel(prefix="l%02d_" % i)
            model_l.append(model)

            if i == 0:
                pars = model.make_params()
            else:
                pars.update(model.make_params())

            #print pars
            plt.plot([line_data[i], line_data[i]], [0, 1e6], 'k--')

            pars['l%02d_center' % i].set(value=line_data[i], 
                min=(line_data[i]-1.), max=(line_data[i]+1.))
            pars['l%02d_sigma' % i].set(value=1.4, max=2.5, min=1.)
            pars['l%02d_amplitude' % i].set(value=1e6)

        print "%20s %20s %20s" % ("Original(nm)", "Measured(nm)", "FWHM(nm)")
        for i in range(nlines):
            if i == 0:
                mod = model_l[i]
            else:
                mod += model_l[i]
 
        init = mod.eval(pars, x=self.wave)
        out=mod.fit((self.spectra1d-self.continuum(self.wave)), pars, x=self.wave)

        for i in range(nlines):
            print "%20.2f %20.2f %20.5f" % (line_data[i], 
                out.best_values["l%02d_center" % i], out.best_values["l%02d_sigma" % i]*2.35)


        #print out.fit_report()
        plt.plot(self.wave, out.best_fit)
        plt.plot(self.wave, self.spectra1d - self.continuum(self.wave))
        plt.show()

    def _polFit(self, xx, yy, order):

        from lmfit.models import PolynomialModel

        model = PolynomialModel(order)
        pars = model.guess(yy, x=xx)
        out = model.fit(yy, pars, x=xx)
        
        return out

    def _continuumFit(self, method='polynomial', order=5, smooth=6):

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
                w=1/self.spectraFiltered, k=order, s=smooth)

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

    def subtractContinuum(self, method='polynomial', order=5, smooth=6):

        plt.figure(figsize=(16, 4))
        out = self._continuumFit(method=method, order=order, smooth=smooth)

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

    def getMedian(self, dark=None, mask=None):

        if not hasattr(self, 'median'):
            image = self.getImage(mask=mask, dark=dark)
            self.median = np.median(image)
        return self.median

    def getObjects(self, mask=None, search=None, thresh=500, minarea=20, force=False):
        
        if not hasattr(self, 'objects') or force:
            image = self.getImage(mask=mask)
            self.objects = sep.extract(image, thresh, minarea=minarea)

        if search is None:
            return self.objects
        else:
            # Caution. This will return first found object within the search radius
            for object in self.objects:
                if ((object['x'] - search['x'])**2 + (object['y'] - search['y'])**2 < search['r']**2):
                    return object

    def plotDistortionPSF(self, grid=5, perc=95, mask=None, title=None):

        objects = self.getObjects()
        im_sz = 20

        xs = np.array([obj['x'] for obj in objects])
        ndxs = xs.argsort()
        objects = objects[ndxs]

        for i in range(grid):
            ys = np.array([obj['y'] for obj in objects[grid*i:grid*(i+1)]])
            ndxs = ys.argsort()
            objects[grid*i:grid*(i+1)] = objects[grid*i:grid*(i+1)][ndxs]

        fig, axs = plt.subplots(grid, grid, figsize=(12,12))

        im = np.zeros([2048+2*im_sz, 2048+2*im_sz])
        if mask is not None:
            im[im_sz:im_sz+2048, im_sz:im_sz+2048] = self.getImage(mask=mask)
        else:
            im[im_sz:im_sz+2048, im_sz:im_sz+2048] = self.getImage()

        norm = ImageNormalize(im, interval=PercentileInterval(perc))

        for i in range(grid):
            for j in range(grid):
                obj = objects[grid*j+i]
                x_c, y_c = int(obj['x'])+im_sz, int(obj['y'])+im_sz

                crop_im = im[y_c-im_sz:y_c+im_sz, x_c-im_sz:x_c+im_sz]
                axs[grid-1-i][j].imshow(crop_im, norm=norm, cmap='gray')
                axs[grid-1-i][j].invert_yaxis()
                axs[grid-1-i][j].get_yaxis().set_ticks([])
                axs[grid-1-i][j].get_xaxis().set_ticks([])

        if title is not None:
            fig.suptitle(title)

    def plotObjects(self, mask=None, figsize=None):

        if figsize is not None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
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

        if 'prefix' in self.config:
            prefix = True
        else:
            prefix = False

        if 'test_ids' in self.config:

            for test_id in self.config['test_ids']:
                full_path = self.config['data_dir'] + '/'+ ('%s-%s' % tuple(test_id.split('-')[:2])) + '/' + test_id + '/'
                files = os.listdir(full_path)
                files.sort()

                for f_name in files:
                    if prefix:
                        if f_name.startswith(self.config['prefix']) and f_name.endswith('fits'):
                            self.NIX_Files.append(NIX_Image(full_path + f_name, test_id=test_id))                        
                    else:
                        if f_name.endswith('fits'):
                            self.NIX_Files.append(NIX_Image(full_path + f_name, test_id=test_id))

        else:

            full_path = self.config['data_dir'] + '/'
            files = os.listdir(full_path)
            files.sort()

            for f_name in files:
                    if prefix:
                        if f_name.startswith(self.config['prefix']) and f_name.endswith('fits'):
                            self.NIX_Files.append(NIX_Image(full_path + f_name))                        
                    else:
                        if f_name.endswith('fits'):
                            self.NIX_Files.append(NIX_Image(full_path + f_name))
    

    def printTable(self, keywords, tbl_fmt):

        table_header =  ['NDX', 'TEST_ID', 'FILENAME']
        table_header.extend([strip_prefix_from_keyword(keyword) for keyword in keywords])
        print ('%4s' + tbl_fmt) % tuple(table_header)
        
        ctr = 0 
        for File in self.NIX_Files:
            header = fits.getheader(File.full_path)
            out_list = [ctr, File.test_id, File.f_name]
            try:
                out_list.extend([header[keyword] for keyword in keywords])
                print ('%04d' + tbl_fmt) % tuple(out_list)
            except:
                out_list.extend(None for keyword in keywords)
                print ('%04d' + tbl_fmt) % tuple(out_list)
            ctr += 1

    def printFiltered(self, keywords, tbl_fmt):

        table_header =  ['NDX', 'TEST_ID', 'FILENAME']
        table_header.extend([strip_prefix_from_keyword(keyword) for keyword in keywords])
        print ('%4s' + tbl_fmt) % tuple(table_header)
        
        ctr = 0 
        for File in self.Filtered:
            header = fits.getheader(File.full_path)
            out_list = [ctr, File.test_id, File.f_name]
            try:
                out_list.extend([header[keyword] for keyword in keywords])
                print ('%04d' + tbl_fmt) % tuple(out_list) 
            except:
                out_list.extend(None for keyword in keywords)
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

    def doLinearRegress(self, exps, select, cpu=1, dark=None, mask=None):

        sz = exps.shape[0]

        exps = np.tile(exps, 2048*2048)
        exps = exps.reshape([2048, 2048, sz])

        ims = self.getImage(select=select, mask=mask, dark=dark)

        angs, lins, Rs = self._multiLinearRegress(exps, ims, cpu=cpu)

        return angs, lins, Rs

    def getMedian(self, mask=None, dark=None):

        return [NI.getMedian(mask=mask, dark=dark) for NI in self.Filtered]

    def getImage(self, select=None, mask=None, dark=None):

        if select is not None:
            sz = select.shape[0]
        else:
            sz = len(self.Filtered)
        ims = np.zeros([2048, 2048, sz])

        if dark is None:
            darks = [None]*sz
        else:
            darks = [self.Filtered[i] for i in range(sz)]

        if select is not None:
            for i in range(sz):
                ims[:,:,i] = self.Filtered[select[i]].getImage(mask=mask, dark=darks[i])
        else:
            for i in range(sz):
                ims[:,:,i] = self.Filtered[i].getImage(mask=mask, dark=darks[i])

        return ims

    def getPixelVariance(self, mask=None, shift=True, gain=None):

        sz = len(self.Filtered)
        data = np.zeros([2048, 2048, sz])
        diff = np.zeros([2048, 2048, sz])

        for i in range(sz):
            data[:,:,i] = self.Filtered[i].getImage(mask=None, gain=gain)
            diff[:,:,i] = data[:,:,i] - data[:,:,0]
            diff[:,:,i] = np.median(diff[:,:,i])

        if shift:
            data[:,:,:] -= diff

        return np.var(data, axis=2, ddof=1), diff[0,0,:]

    def getObjects(self, mask=None, search=None, sepdict=None, thresh=500, minarea=20, force=False):
        
        if sepdict is not None:
            return [NI.getObjects(mask=mask, search=search, thresh=thresh, minarea=minarea, force=force)[sepdict] for NI in self.Filtered]
        else:
            return [NI.getObjects(mask=mask, search=search, thresh=thresh, minarea=minarea, force=force) for NI in self.Filtered]

    def getHeaderValue(self, keyword):

        return [fits.getheader(NI.full_path)[keyword] for NI in self.Filtered]

    def _multiLinearRegress(self, data1, data2, cpu=1):

        pb = ProgressBar(function="_multiLinearRegress")
        angs = np.zeros([2048, 2048])
        lins = np.zeros([2048, 2048])
        Rs = np.zeros([2048, 2048])
        pool = mp.Pool(cpu)
        for i in range(2048):
            arr1s = [data1[i,j,:] for j in range(2048)]
            arr2s = [data2[i,j,:] for j in range(2048)]
            out = pool.map(_singleLinearRegress, zip(arr1s, arr2s))
            for j in range(2048):
                angs[i,j], lins[i,j], Rs[i,j] = out[j][0], out[j][1], out[j][2]
            
            pb.report(i+1, 2048)

        return angs, lins, Rs


# Some utility functions
def strip_prefix_from_keyword(keyword):

    if keyword.startswith('HIERARCH ESO'):
        keyword = keyword[13:]

    keyword = '_'.join(keyword.split(' '))
    return keyword


def doGridAnalysis(data, grid, window, start, func, factor=1., index=None):

    sz = data.shape[2]
    result = np.zeros([grid*grid*sz])

    for k in range(sz):
        ndx = index if index is not None else k
        for i in range(grid):
            for j in range(grid):
                result[i*grid*sz+j*sz+k] = func(data[start+window*i:start+window*(i+1),
                            start+window*j:start+window*(j+1), ndx])

    return result*factor

def _singleLinearRegress(args):

    arr1, arr2 = args
    ang, lin, R, _, _ = linregress(arr1, arr2)
    return ang, lin, R


class ProgressBar:

    __slots__ = ("t_start", "function")

    def __init__(self, function=None):
        self.t_start = time.time()
        if function is not None:
            self.function = function
        else:
            self.function = ""

    def init(self, function=None):
        self.__init__(function=function)

    def report(self, current, full):
        progress = float(current)/full*100
        bar_done = int(round(progress/5.))
        t_now = time.time()
        t_elapsed = t_now-self.t_start
        t_full = t_elapsed/progress*100

        clear_output(wait=True)
        print "%s ==> [%s%s] %3d%% %s/%s" % (self.function, "#"*bar_done, '-'*(20-bar_done), progress,
            self._formatTime(t_elapsed), self._formatTime(t_full))

    def _formatTime(self, t):

        mins = t / 60.
        secs = t % 60

        return "%02d:%02d" % (mins, secs)


