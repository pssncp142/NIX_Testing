from ctypes import c_int, c_float, c_double, c_size_t, Structure
import numpy as np
import numpy.ctypeslib as npct
from astropy.io import fits
import os

basedir = os.path.abspath(os.path.dirname(__file__)) + '/lib/'
print basedir

hdrl = npct.load_library("libhdrldemo.so", basedir)

arr_1d_dbl = npct.ndpointer(dtype=np.float64, ndim=1, flags='CONTIGUOUS') 
arr_1d_ubyte = npct.ndpointer(dtype=np.ubyte, ndim=1, flags='CONTIGUOUS') 


class hdrl_value(Structure):
    _fields_=[('data', c_double), ('error', c_double)]

class hdrl_strehl_result(Structure):
    _fields_=[('strehl_value', hdrl_value), ('star_x', c_double), ('star_y', c_double), 
            ('star_peak', hdrl_value), ('star_flux', hdrl_value), ('star_background', hdrl_value), 
            ('computed_background_error', c_double), ('nbackground_pixels', c_size_t)]

def compute_strehl(image, wave, m1_rad, m2_rad, pixsc, flux_r, bkg_low, bkg_high):

    hdrl.hdrl_compute_strehl_numpy_float64.restype = hdrl_strehl_result
    hdrl.hdrl_compute_strehl_numpy_float64.argtypes = [arr_1d_dbl, c_int, c_int,
                                                            c_double, c_double, c_double,
                                                            c_double, c_double, 
                                                            c_double, c_double, c_double]

    sz = image.shape[0]
    image = np.ravel(image).astype(np.float64)

    strehl = hdrl.hdrl_compute_strehl_numpy_float64(image, sz, sz, wave, m1_rad, m2_rad, pixsc, pixsc,
                                flux_r, bkg_low, bkg_high)

    return strehl

def bpm_2d_compute(image):

    hdrl.hdrl_bpm_2d_compute_numpy_float64.restype = c_int
    hdrl.hdrl_bpm_2d_compute_numpy_float64.argtypes = [arr_1d_dbl, arr_1d_ubyte, c_int, c_int]

    szx = image.shape[0]
    szy = image.shape[1]
 
    image = np.ravel(image).astype(np.float64)
    mask = np.zeros([szx*szy]).astype(np.ubyte)

    res = hdrl.hdrl_bpm_2d_compute_numpy_float64(image, mask, szx, szy)

    return np.reshape(mask, [szx, szy])

def bpm_3d_compute(images):

    hdrl.hdrl_bpm_3d_compute_numpy_float64.restype = c_int
    hdrl.hdrl_bpm_3d_compute_numpy_float64.argtypes = [arr_1d_dbl, arr_1d_ubyte, c_int, c_int, c_int]

    szx = images.shape[0]
    szy = images.shape[1]
    szz = images.shape[2]

    #ims_in = np.zeros([szz*szx*szy]).astype(np.float64)
    ims_out = np.zeros([szz*szx*szy]).astype(np.ubyte)

    ims_in = np.ravel(np.moveaxis(images, -1, 0)).astype(np.float64)

    res = hdrl.hdrl_bpm_3d_compute_numpy_float64(ims_in, ims_out, szx, szy, szz)

    return np.moveaxis(np.reshape(ims_out, [szz, szx, szy]), 0, -1)

def bpm_fit_compute(images, exptime):

    hdrl.hdrl_bpm_fit_compute_numpy_float64.restype = c_int
    hdrl.hdrl_bpm_fit_compute_numpy_float64.argtypes = [arr_1d_dbl, arr_1d_ubyte, arr_1d_dbl, c_int, c_int, c_int]

    szx = images.shape[0]
    szy = images.shape[1]
    szz = images.shape[2]

    exptime = exptime.astype(np.float64)
    ims_out = np.zeros([szx*szy]).astype(np.ubyte)

    ims_in = np.ravel(np.moveaxis(images, -1, 0)).astype(np.float64)

    res = hdrl.hdrl_bpm_fit_compute_numpy_float64(ims_in, ims_out, exptime, szx, szy, szz)

    return np.reshape(ims_out, [szx, szy])


def collapse_median(images):

    hdrl.hdrl_collapse_median_numpy_float64.restype = c_int
    hdrl.hdrl_collapse_median_numpy_float64.argtypes = [arr_1d_dbl, arr_1d_dbl,
                                                        c_int, c_int, c_int]
    szx = images.shape[0]
    szy = images.shape[1]
    szz = images.shape[2]

    im_out = np.zeros([szx*szy]).astype(np.float64)

    ims_in = np.ravel(np.moveaxis(images, -1, 0)).astype(np.float64)

    res = hdrl.hdrl_collapse_median_numpy_float64(ims_in, im_out, szx, szy, szz)

    return np.reshape(im_out, [szx, szy])







