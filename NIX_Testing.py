from astropy.io import fits
import numpy as np
import sep
import os
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

class NIX_Image:

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

    def getImage(self, hdu=None, mask=None):

        hdul = fits.open(self.full_path)
        image = hdul[0].data*1. 
        # Handle exception where data format is (1, 2048, 2048)
        if image.shape[0] == 1:
            image = image[0,:,:]

        if mask is not None:
            image = image * mask.getImage()
        return image

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


    def plotImage(self, mask=None):

        plt.imshow(self.getImage(mask=mask))
        plt.gca().invert_yaxis()
        plt.show()

        
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





