from ctypes import c_int, c_uint, c_float, c_double, c_bool, c_void_p, c_size_t, POINTER, Structure, byref
import numpy as np
import numpy.ctypeslib as npct
from astropy.io import fits
import pkg_resources
import os

from enum import IntEnum

class CtypesEnum(IntEnum):
    @classmethod
    def from_param(cls, obj):
        return int(obj)

class CPL_BORDER(CtypesEnum):
    FILTER=0
    ZERO=1
    CROP=2
    NOP=3
    COPY=4

class CPL_FILTER(CtypesEnum):
    EROSION=0
    DILATION=1
    OPENING=2
    CLOSING=3
    LINEAR=4
    LINEAR_SCALE=5
    AVERAGE=6
    AVERAGE_FAST=7
    MEDIAN=8
    STDEV=9
    STDEV_FAST=10
    MORPHO=11
    MORPHO_SCALE=12

so_name = pkg_resources.resource_filename(__name__, 'lib/libhdrldemo.so')
hdrl = npct.load_library(so_name, ".")

arr_1d_dbl = npct.ndpointer(dtype=np.float64, ndim=1, flags='CONTIGUOUS') 
arr_1d_ubyte = npct.ndpointer(dtype=np.ubyte, ndim=1, flags='CONTIGUOUS') 

ERROR_IMAGE = c_bool.in_dll(hdrl, 'ERROR_IMAGE')
ERROR_METHOD = c_int.in_dll(hdrl, 'ERROR_METHOD')
RN_ADU = c_double.in_dll(hdrl, 'RN_ADU')
GAIN = c_double.in_dll(hdrl, 'GAIN')

class hdrl_value(Structure):
    _fields_=[('data', c_double), ('error', c_double)]

class hdrl_strehl_result(Structure):
    _fields_=[('strehl_value', hdrl_value), ('star_x', c_double), ('star_y', c_double), 
            ('star_peak', hdrl_value), ('star_flux', hdrl_value), ('star_background', hdrl_value), 
            ('computed_background_error', c_double), ('nbackground_pixels', c_size_t)]

class hdrl_bpm_fit_parameter(Structure):
    _fields_=[('base', c_void_p), ('degree', c_int), ('pval', c_double),
            ('rel_chi_l', c_double), ('rel_chi_h', c_double),
            ('rel_coef_l', c_double), ('rel_chi_h', c_double)]

class hdrl_bpm_2d_parameter(Structure):
    _fields_=[('base', c_void_p), ('filter', c_int), ('border', c_int),
            ('kappa_low', c_double), ('kappa_high', c_double), ('maxiter', c_int),
            ('steps_x', c_int), ('steps_y', c_int),
            ('filter_size_x', c_int), ('filter_size_y', c_int),
            ('order_x', c_int), ('order_y', c_int),
            ('smooth_x', c_int), ('smooth_y', c_int), ('method', c_int)]

def init():

    hdrl.hdrl_init.restype = None
    hdrl.hdrl_init()

def end():

    hdrl.hdrl_end.restype = None
    hdrl.hdrl_end()

def bpm_2d_parameter_create(filter, border,
                            kappa_low, kappa_high, maxiter,
                            steps_x, steps_y,
                            filter_size_x, filter_size_y,
							order_x, order_y,
                            smooth_x, smooth_y):

    hdrl.hdrl_bpm_2d_parameter_create.restype = POINTER(hdrl_bpm_2d_parameter)
    hdrl.hdrl_bpm_2d_parameter_create.argtypes = [c_int, c_int,
                                                    c_double, c_double, c_int,
                                                    c_int, c_int, c_int, c_int,
                                                    c_int, c_int, c_int, c_int]

    params = hdrl.hdrl_bpm_2d_parameter_create(filter, border,
                            kappa_low, kappa_high, maxiter,
                            steps_x, steps_y,
                            filter_size_x, filter_size_y,
							order_x, order_y,
                            smooth_x, smooth_y)

    return params

def bpm_fit_parameter_create(degree, pval, rel_chi_l, rel_chi_h, rel_coef_l, rel_coef_h):

    hdrl.hdrl_bpm_fit_parameter_create.restype = POINTER(hdrl_bpm_fit_parameter)
    hdrl.hdrl_bpm_fit_parameter_create.argtypes = [c_int, c_double,
                                                        c_double, c_double,
                                                        c_double, c_double]

    params = hdrl.hdrl_bpm_fit_parameter_create(degree, pval, rel_chi_l, rel_chi_h, rel_coef_l, rel_coef_h)

    return params

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

def bpm_2d_compute(image, method='filter', kappa_low=5., kappa_high=10., maxiter=10,
                                                steps_x=5, steps_y=5,
                                                filter_size_x=5, filter_size_y=5,
                                                order_x=2, order_y=2,
                                                smooth_x=5, smooth_y=5,
                                                filter=CPL_FILTER.MEDIAN, border=CPL_BORDER.FILTER):

    hdrl.hdrl_bpm_2d_compute_numpy_float64.restype = c_int
    hdrl.hdrl_bpm_2d_compute_numpy_float64.argtypes = [arr_1d_dbl, arr_1d_ubyte, POINTER(hdrl_bpm_2d_parameter),
                                                        c_int, c_int]

    if (method == 'filter'):
        params = bpm_2d_parameter_create(filter, border, kappa_low, kappa_high, maxiter,
                                -1, -1, -1, -1, -1, -1, smooth_x, smooth_y)
    elif (method == 'legendre'):
        params = bpm_2d_parameter_create(filter, border, kappa_low, kappa_high, maxiter,
                                                steps_x, steps_y,
                                                filter_size_x, filter_size_y,
                                                order_x, order_y, -1, -1)
    else:
        print("Undefined method")
        return -1

    szx = image.shape[0]
    szy = image.shape[1]
 
    image = np.ravel(image).astype(np.float64)
    mask = np.zeros([szx*szy]).astype(np.ubyte)

    res = hdrl.hdrl_bpm_2d_compute_numpy_float64(image, mask, params, szx, szy)

    return np.reshape(mask, [szx, szy])

def bpm_3d_compute(images):

    hdrl.hdrl_bpm_3d_compute_numpy_float64.restype = c_int
    hdrl.hdrl_bpm_3d_compute_numpy_float64.argtypes = [arr_1d_dbl, arr_1d_ubyte, c_int, c_int, c_int]

    szx = images.shape[0]
    szy = images.shape[1]
    szz = images.shape[2]

    ims_out = np.zeros([szz*szx*szy]).astype(np.ubyte)

    ims_in = np.ravel(np.moveaxis(images, -1, 0)).astype(np.float64)

    res = hdrl.hdrl_bpm_3d_compute_numpy_float64(ims_in, ims_out, szx, szy, szz)

    return np.moveaxis(np.reshape(ims_out, [szz, szx, szy]), 0, -1)

def bpm_fit_compute(images, exptime, method='pval', degree=1, pval=10.,
                    rel_chi_low=3., rel_chi_high=10., rel_coef_low=3., rel_coef_high=10.):

    hdrl.hdrl_bpm_fit_compute_numpy_float64.restype = c_int
    hdrl.hdrl_bpm_fit_compute_numpy_float64.argtypes = [arr_1d_dbl, arr_1d_ubyte, arr_1d_dbl, POINTER(hdrl_bpm_fit_parameter),
                                                         c_int, c_int, c_int]

    szx = images.shape[0]
    szy = images.shape[1]
    szz = images.shape[2]

    exptime = exptime.astype(np.float64)
    ims_out = np.zeros([szx*szy]).astype(np.ubyte)

    ims_in = np.ravel(np.moveaxis(images, -1, 0)).astype(np.float64)

    if method=='pval':
        params = bpm_fit_parameter_create(degree, pval, -1, -1, -1, -1)
    elif method=='rel_chi':
        params = bpm_fit_parameter_create(degree, -1, rel_chi_low, rel_chi_high, -1, -1)
    elif method=='rel_coef':
        params = bpm_fit_parameter_create(degree, -1, -1, -1, rel_coef_low, rel_coef_high)
    else:
        print("Undefined method")
        return -1

    res = hdrl.hdrl_bpm_fit_compute_numpy_float64(ims_in, ims_out, exptime, params, szx, szy, szz)

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







