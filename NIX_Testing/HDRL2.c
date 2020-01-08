#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include "structmember.h"
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "hdrl.h"
#include "cpl.h"
#include "cpl_image.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// Documentation
PyDoc_STRVAR(mod__doc__, "Some utility functions from HDRL to be used with numpy arrays.\n");

PyDoc_STRVAR(bpm_fit_compute__doc__, 
"bpm_fit_compute(image, exptime, [method, degree, pval, rel_chi_low, rel_chi_high, rel_coef_low, rel_coef_high])\n"
"\n"
"A python wrapper for 'hdrl_bpm_fit_compute' function. 'hdrl_bpm_fit_parameter' is created inside the wrapper "
"based on the keyword values.\n"
"\n"
"Computes bad pixel mask from a set of flat images with a different exposure time. The algorithm performs fit "
"on a stack of images and finds outliers using one of three defined methods.\n"
"\n"
"Parameters\n"
"----------\n"
"\n"
"images : (N,N,N,) ndarray\n"
"\tStack of images - along the third axis.\n"
"exptime : (N,) ndarray\n"
"\tExptime corresponding to images. Should have same length as third axis of **images**.\n"
"method : str, optional \n"
"\tMethod used to determine bad pixels. If given method is none of the below,"
"the function will throw an error. Accepted values are (default='pval'):\n\n"
"\t* 'pval'    : threshold on the p-values.\n"
"\t* 'rel_chi' : sigma cut on the relative chi values.\n"
"\t* 'rel_coef': sigma cut on the relative coefficients of polynomial fit.\n"
"degree : int, optional\n"
"\tDegree of polynomial to fit. (default=1)\n"
"pval : double, optional\n"
"\tThreshold in p-values. Enabled with *pval* (default=10.)\n"
"rel_chi_low, rel_chi_high : double, optional\n"
"\tLow/High threshold in relative chi values. Enabled with *rel_chi* (default=(5., 10.))\n"
"rel_coef_low, rel_coef_high : double, optional\n"
"\tLow/High threshold in relative coefficent values. Enabled with *rel_coef* (default=(5., 10.))\n"
"\n"
"Returns\n"
"-------\n"
"bp_mask : (N,N) ndarray\n"
"\tOutput of *hdrl_bpm_fit_compute* converted to a numpy array. Good pixels are marked "
"with 0 and bad pixels with 1. When *rel_coef* method used, algorithm computes bad "
"pixels for each coefficient of the polynomial. Hence, the resulting bad pixel map is "
"[2^0(0th)+2^1(1st)+2^2(2nd)+..] where bad pixel map derived from Nth order coefficients "
"is given by (Nth). If the function fails, it returns *Py_None* instead of "
"throwing an error."
"\n");

PyDoc_STRVAR(bpm_2d_compute__doc__, 
"bpm_2d_compute(image, [method, kappa_low, kappa_high, maxiter, steps_x, steps_y,"
"filter_size_x, filter_size_y, order_x, order_y, smooth_x, smooth_y)\n"
"A python wrapper for 'hdrl_bpm_2d_compute' function. 'hdrl_bpm_2d_parameter' is created inside the wrapper "
"based on the keyword values.\n"
"\n"
"Parameters\n"
"----------\n"
"\n"
"image : (N,N,) ndarray\n"
"\tAn image.\n"
"kappa_low : double, optional\n"
"\t(default=5.)\n"
"kappa_high : double, optional\n"
"\t(default=10.)\n"
"maxiter : int, optional\n"
"\tNumber of iterations (default=10).\n"
"method : (N,N,) ndarray\n"
"\tMethod used to determine bad pixels. If given method is none of the below,"
"the function will throw an error. Accepted values are (default='filter'):\n\n"
"\t* 'filter'    : threshold on the p-values.\n"
"\t* 'legendre'  : sigma cut on the relative chi values.\n"
"steps_x,steps_y : int, optional\n"
"\t(default=10)\n"
"filter_size_x,filter_size_y : int, optional\n"
"\t(default=10)\n"
"order_x,order_y : int, optional\n"
"\t(default=2)\n"
"smooth_x,smooth_y : int, optional\n"
"\t(default=5)\n"
"\n"
"Returns\n"
"-------\n"
"out : (N,N,) ndarray\n"
"\tReturns an image in which bad pixels are interpolated."
"\n");

PyDoc_STRVAR(bpm_interpolate__doc__, 
"bpm_interpolate(image, mask)\n"
"A python wrapper for CPL function "
"`cpl_detector_interpolate_rejected <https://www.eso.org/sci/software/cpl/reference/group__cpl__detector.html#ga83fd1a1d48eeeb444d6a1a84c6d5e2de>`_ \n"
"\n" 
"Parameters\n"
"----------\n"
"\n"
"image : (N,N,) ndarray\n"
"\tAn image.\n"
"mask : (N,N,) ndarray\n"
"\tBad pixel masks indicating pixels to interpolate.\n"
"\n"
"Returns\n"
"-------\n"
"out : (N,N,) ndarray\n"
"\tReturns an image in which bad pixels are interpolated."
"\n");

PyDoc_STRVAR(compute_strehl__doc__, 
"compute_strehl(image, wave, m1_rad, m2_rad, pixsc_x, pixsc_y, flux_r, bkg_low_r, bkg_high_r)\n"
"A python wrapper for hdrl_compute_strehl function\n"
"\n" 
"Parameters\n"
"----------\n"
"\n"
"image : (N,N,) ndarray\n"
"\tImage to compute Strehl ratio, must contain a point source.\n"
"wave : double\n"
"\tWavelength to compute PSF.\n"
"m1_rad,m2_rad : double\n"
"\tRadius (in m) of m1 and m2. (Or, radius of the aperture and obscuration.)\n"
"pixsc_x,pixsc_y : double\n"
"\tPixel scale in (x,y) arcsec/pix.\n"
"flux_r : double\n"
"\tFlux radius to estimate source count within an aperture. (in arcsec)\n"
"bkg_low_r,bkg_high_r : double\n"
"\tInner and outer radius of the sky annulus to estimate background (in arcsec)\n"
"\n"
"Returns\n"
"-------\n"
"out : HDRL2.strehl_result\n"
"\tPythonized version of hdrl_strehl_result struct. out.strehl_value.data gives the Strehl ratio.."
"\n");

PyDoc_STRVAR(lacosmic_edgedetect__doc__, 
"lacosmic_edgedetect(image, wave, m1_rad, m2_rad, pixsc_x, pixsc_y, flux_r, bkg_low_r, bkg_high_r)\n"
"A python wrapper for hdrl_lacosmic_edgedetect function\n"
"\n" 
"hdrl_lacosmics.c\n"
"----------------\n"
"This routine determines bad-pixels on a single image via edge detection\n"
"following the algorithm (LA-Cosmic) describe in van Dokkum,\n"
"PASP,113,2001,p1420-27. The HDRL implementation does not use use error model\n"
"as described in the paper but the error image passed to the function. Moreover\n"
"we do several iterations and replace the detected bad pixels in each iteration\n"
"by the information of the surrounding pixels.\n"
"\n"
"Parameters\n"
"----------\n"
"\n"
"image : (N,N,) ndarray\n"
"\tImage to compute Strehl ratio, must contain a point source.\n"
"sigma_lim : double\n"
"\tLimiting sigma for detection on the sampling image\n"
"f_lim : double\n"
"\tLimiting f factor for detection on the modified Laplacian image.\n"
"maxiter : int\n"
"\tMaximum number of iterations"
"\n"
"Returns\n"
"-------\n"
"out : (N,N) ndarray\n"
"\tBad Pixel mask."
"\n");

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// Just an example
/*
void HDRL_simplePrint(char * text){
    
    PyObject * sys = PyImport_ImportModule( "sys");
    PyObject * s_out = PyObject_GetAttrString(sys, "stdout"); 
    PyObject * result = PyObject_CallMethod(s_out, "write", "s", strcat(text, "\n")); 
    Py_DECREF(result);
    Py_DECREF(s_out);
    Py_DECREF(sys);

}*/

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//hdrl_value definition

typedef struct {
    PyObject_HEAD
    hdrl_data_t data;
    hdrl_error_t error;
} py_hdrl_value;

static PyObject *
hdrlValue_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    py_hdrl_value *self;
    self = (py_hdrl_value *) type->tp_alloc(type, 0);
    self->data = 0;
    self->error = 0;

    return (PyObject *) self;
}

static int
hdrlValue_init(py_hdrl_value *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"data", "error", NULL};
    double data=0, error=0;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist,
                                     &data, &error))
        return -1;

    self->data = data;
    self->error = error;

    return 0;
}

static void
hdrlValue_dealloc(py_hdrl_value *self)
{
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
hdrlValue_repr(py_hdrl_value *self)
{
    char *out = malloc(sizeof(char)*10);
    sprintf(out, "(%.3f +- %.3f)", self->data, self->error);
    PyObject * res = Py_BuildValue("s", out);
    free(out);
    return res;
}

static PyMemberDef hdrlValue_members[] = {
    {"data", T_DOUBLE, offsetof(py_hdrl_value, data), 0,
     "Stored value of the variable."},
    {"error", T_DOUBLE, offsetof(py_hdrl_value, error), 0,
     "Stored error of the variable."},
    {NULL} 
};

static PyTypeObject hdrlValueType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "HDRL2.hdrl_value",
    .tp_doc = "A minimalistic class imitating hdrl_value.",
    .tp_basicsize = sizeof(py_hdrl_value),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = hdrlValue_new,
    .tp_init = (initproc) hdrlValue_init,
    .tp_dealloc = (destructor) hdrlValue_dealloc,
    .tp_repr = (reprfunc) hdrlValue_repr,
    .tp_members = hdrlValue_members
};

static PyObject *
hdrlValue_initFromC(double data, double error)
{
    py_hdrl_value * out = (py_hdrl_value*) hdrlValue_new(&hdrlValueType, NULL, NULL);
    out->data = data;
    out->error = error;

    return (PyObject*) out;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// strehl_result definition
typedef struct {
    PyObject_HEAD
    PyObject* strehl_value;
    double star_x, star_y;
    PyObject* star_peak;
    PyObject* star_flux;
    PyObject* star_background;
    double computed_background_error;
    Py_ssize_t nbackground_pixels;
} py_strehl_result;

static PyObject *
strehlResult_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    py_strehl_result *self;
    self = (py_strehl_result *) type->tp_alloc(type, 0);
    self->star_x = 0;
    self->star_y = 0;
    self->computed_background_error = 0;
    self->nbackground_pixels = 0;
    self->strehl_value = hdrlValue_initFromC(0,0);
    self->star_peak = hdrlValue_initFromC(0,0);
    self->star_flux = hdrlValue_initFromC(0,0);
    self->star_background = hdrlValue_initFromC(0,0);

    return (PyObject *) self;
}

static void
strehlResult_dealloc(py_strehl_result *self)
{
    Py_XDECREF(self->strehl_value);
    Py_XDECREF(self->star_peak);
    Py_XDECREF(self->star_flux);
    Py_XDECREF(self->star_background);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyMemberDef strehlResult_members[] = {
    {"strehl_value", T_OBJECT_EX, offsetof(py_strehl_result, strehl_value), 0, "strehl value (hdrl_value)"},
    {"star_x", T_DOUBLE, offsetof(py_strehl_result, star_x), 0, "X center in pixel (double)"},
    {"star_y", T_DOUBLE, offsetof(py_strehl_result, star_y), 0, "Y center in pixel (double)"},
    {"star_peak", T_OBJECT_EX, offsetof(py_strehl_result, star_peak), 0, "Star peak (hdrl_value)"},
    {"star_flux", T_OBJECT_EX, offsetof(py_strehl_result, star_flux), 0, "Star flux (hdrl_value)"},
    {"star_background", T_OBJECT_EX, offsetof(py_strehl_result, star_background), 0, "Star Background (hdrl_value)"},
    {"computed_background_err", T_DOUBLE, offsetof(py_strehl_result, computed_background_error), 0, "not sure (double)"},
    {"nbackground_pixels", T_PYSSIZET, offsetof(py_strehl_result, nbackground_pixels), 0, "Number of background puxels (size_t)"},
    {NULL} 
};

static PyTypeObject strehlResultType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "HDRL2.strehl_result",
    .tp_doc = "A minimalistic class imitating hdrl_strehl_result. Some attributes derived from hdrl_value.",
    .tp_basicsize = sizeof(py_strehl_result),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = strehlResult_new,
    .tp_dealloc = (destructor) strehlResult_dealloc,
    .tp_members = strehlResult_members
};

static PyObject * strehlResult_initFromOriginal(hdrl_strehl_result result){

    py_strehl_result * out = (py_strehl_result*) strehlResult_new(&strehlResultType, NULL, NULL);

    out->star_x = result.star_x;
    out->star_y = result.star_y;
    out->computed_background_error = result.computed_background_error;
    out->nbackground_pixels = result.nbackground_pixels;
    out->strehl_value = hdrlValue_initFromC(result.strehl_value.data,
                                            result.strehl_value.error);
    out->star_peak = hdrlValue_initFromC(result.star_peak.data,
                                            result.star_peak.error);
    out->star_flux= hdrlValue_initFromC(result.star_flux.data,
                                            result.star_flux.error);
    out->star_background = hdrlValue_initFromC(result.star_background.data,
                                            result.star_background.error);

    return (PyObject *) out;

    
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// Global variables

static int ERROR_IMAGE = 1;
static int ERROR_METHOD = 0;
static double GAIN = 5.7;
static double RN_ADU = 4.5;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// Data conversion from numpy arrays

cpl_image * _cpl_image_from_numpy2d(PyArrayObject* arr, npy_intp * len){

	cpl_image * c_im = cpl_image_new(len[0], len[1], CPL_TYPE_DOUBLE);
    //double * p_arr = PyArray_DATA(arr);
    //double * p_im = cpl_image_get_data(c_im);
    //for (int i=0; i<len[0]*len[1]; i++){
    //    p_im[i] = p_arr[i];
    //}
    
	memcpy(cpl_image_get_data_double(c_im), 
			PyArray_DATA(arr), 
			sizeof(double)*len[0]*len[1]);
	return c_im;
    
}

cpl_mask * _cpl_mask_from_numpy2d(PyArrayObject* arr, npy_intp * len){

	cpl_mask * c_im = cpl_mask_new(len[0], len[1]);
    memcpy((cpl_binary *) cpl_mask_get_data(c_im), 
			(cpl_binary *) PyArray_DATA(arr), 
			sizeof(cpl_binary)*len[0]*len[1]);
    return c_im;

}

cpl_vector * _cpl_vector_from_numpy1d(PyArrayObject* arr, npy_intp * len){

	cpl_vector * c_vec = cpl_vector_new(len[0]);
	memcpy((double *) cpl_vector_get_data(c_vec), 
			(double *) PyArray_DATA(arr), 
			sizeof(double)*len[0]);
	return c_vec;

}

PyArrayObject * _numpy2d_from_cpl_image(cpl_image * c_im, npy_intp * len){

	PyArrayObject * arr = (PyArrayObject*) PyArray_SimpleNew(2, len, NPY_FLOAT64);
	memcpy((double *) PyArray_DATA(arr),
			(double *) cpl_image_get_data(c_im), 
			sizeof(double)*len[0]*len[1]);
	return arr;

}

PyArrayObject * _numpy2d_from_cpl_image_INT(cpl_image * c_im, npy_intp * len){

	PyArrayObject * arr = (PyArrayObject*) PyArray_SimpleNew(2, len, NPY_INT);
	memcpy((int *) PyArray_DATA(arr),
			(int *) cpl_image_get_data(c_im), 
			sizeof(int)*len[0]*len[1]);
	return arr;

}

PyArrayObject * _numpy2d_from_hdrl_image(hdrl_image * h_im, npy_intp *len){

    cpl_image * c_im = hdrl_image_get_image(h_im);
    return _numpy2d_from_cpl_image(c_im, len);

}

hdrl_image * _hdrl_image_from_numpy2d(PyArrayObject* arr, npy_intp *len){

    cpl_image * c_im = _cpl_image_from_numpy2d(arr, len);
    cpl_image * c_err;

	if (ERROR_IMAGE == true) {

		if (ERROR_METHOD == 0){
			c_err = cpl_image_abs_create(c_im);
			cpl_image_multiply_scalar(c_err, GAIN);
			cpl_image_add_scalar(c_err, RN_ADU*RN_ADU*GAIN*GAIN);
		} else {
			c_err = cpl_image_new(len[0], len[1], CPL_TYPE_DOUBLE);
			cpl_image_add_scalar(c_err, cpl_image_get_median(c_im)*GAIN);
		}
        cpl_image_power(c_err, 0.5);
		cpl_image_divide_scalar(c_err, GAIN);
	} else {
		c_err = cpl_image_new(len[0], len[1], CPL_TYPE_DOUBLE);
	}

	hdrl_image * image_hdrl;
    image_hdrl = hdrl_image_new(2048, 2048);
    hdrl_image_insert(image_hdrl, c_im, c_err, 1, 1);

    cpl_image_delete(c_im);
    cpl_image_delete(c_err);

    return image_hdrl;

}

hdrl_imagelist * _hdrl_imagelist_from_numpy2d(PyArrayObject* arr, npy_intp *len){

    hdrl_image * im_tmp;
    hdrl_imagelist * imlist;
    PyObject * arr_ndx, *arr_tmp;
    int * data_tmp;
    npy_intp dims[NPY_MAXDIMS];
    dims[0] = 1;

    arr_ndx = PyArray_SimpleNew(1, dims, NPY_INT);
    data_tmp = (int*) PyArray_DATA((PyArrayObject*) arr_ndx);

    imlist = hdrl_imagelist_new();

    for (int i=0; i<len[2]; i++){
        data_tmp[0] = i;
        arr_tmp = PyArray_TakeFrom(arr, arr_ndx, 2, NULL, NPY_RAISE);
        //arr_tmp = (PyArrayObject*) PyArray_SimpleNew(3, len, NPY_FLOAT64);
        im_tmp = _hdrl_image_from_numpy2d((PyArrayObject*) arr_tmp, len);
        //im_tmp = hdrl_image_new(2048, 2048);
        hdrl_imagelist_set(imlist, im_tmp, i);
        Py_DECREF(arr_tmp);
    }

    Py_DECREF(arr_ndx);

    return imlist;

}

PyArrayObject * _numpy2d_from_cpl_imagelist(cpl_imagelist * imlist, npy_intp *len){

    PyObject * im2D;
    PyArrayObject * im3D = (PyArrayObject *) PyArray_SimpleNew(3, len, NPY_FLOAT64);
    PyObject * indices;

    cpl_image * im_tmp;

   
    for (int i=0; i<len[2]; i++){
   
        im_tmp = cpl_imagelist_get(imlist, i);
        im2D = (PyObject *) _numpy2d_from_cpl_image(im_tmp, len);
        indices = PyArray_Arange(i, len[0]*len[1]*len[2], len[2], NPY_INTP);
        PyArray_PutTo(im3D, im2D, indices, NPY_RAISE);
        Py_DECREF(im2D);
        Py_DECREF(indices);
    }

    return im3D;

}

PyArrayObject * _numpy2d_from_cpl_imagelist_INT(cpl_imagelist * imlist, npy_intp *len){

    PyObject * im2D;
    PyArrayObject * im3D = (PyArrayObject *) PyArray_SimpleNew(3, len, NPY_INT);
    PyObject * indices;

    cpl_image * im_tmp;

   
    for (int i=0; i<len[2]; i++){
   
        im_tmp = cpl_imagelist_get(imlist, i);
        im2D = (PyObject *) _numpy2d_from_cpl_image_INT(im_tmp, len);
        indices = PyArray_Arange(i, len[0]*len[1]*len[2], len[2], NPY_INT);
        PyArray_PutTo(im3D, im2D, indices, NPY_RAISE);
        Py_DECREF(im2D);
        Py_DECREF(indices);
    }

    return im3D;

}


PyArrayObject * _numpy2d_from_hdrl_imagelist(hdrl_imagelist * imlist, npy_intp *len){

    PyObject * im2D;
    PyArrayObject * im3D = (PyArrayObject*) PyArray_SimpleNew(3, len, NPY_FLOAT64);
    PyObject * indices; 

    hdrl_image * im_tmp;

    printf("1\n");
    for (int i=0; i<len[2]; i++){
    printf("2\n");
        im_tmp = hdrl_imagelist_get(imlist, i);
    printf("3\n");
        im2D = (PyObject*) _numpy2d_from_hdrl_image(im_tmp, len);
    printf("4\n");
        indices = PyArray_Arange(i, len[0]*len[1]*len[2], len[2], NPY_INTP);
    printf("5\n");
        PyArray_PutTo(im3D, im2D, indices, NPY_RAISE);
    printf("6\n");
        Py_DECREF(im2D);
        Py_DECREF(indices);
    }

    return im3D;

}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// Converter functions

cpl_filter_mode _cplFilterConverter(char * filter_str){

    cpl_filter_mode filter; 

    if(strcmp(filter_str, "erosion") == 0){
        filter = CPL_FILTER_EROSION;
    } else if(strcmp(filter_str, "dilation") == 0){
        filter = CPL_FILTER_DILATION;
    } else if(strcmp(filter_str, "opening") == 0){
        filter = CPL_FILTER_OPENING;
    } else if(strcmp(filter_str, "closing") == 0){
        filter = CPL_FILTER_CLOSING;
    } else if(strcmp(filter_str, "linear") == 0){
        filter = CPL_FILTER_LINEAR;
    } else if(strcmp(filter_str, "average") == 0){
        filter = CPL_FILTER_AVERAGE;
    } else if(strcmp(filter_str, "average_fast") == 0){
        filter = CPL_FILTER_AVERAGE_FAST;
    } else if(strcmp(filter_str, "median") == 0){
        filter = CPL_FILTER_MEDIAN;
    } else if(strcmp(filter_str, "stdev") == 0){
        filter = CPL_FILTER_STDEV;
    } else if(strcmp(filter_str, "stdev_fast") == 0){
        filter = CPL_FILTER_STDEV_FAST;
    } else if(strcmp(filter_str, "morpho") == 0){
        filter = CPL_FILTER_MORPHO;
    } else {
        PyErr_Format(PyExc_ValueError, "Undefined filter : %s", filter_str);
        return -1;}

    return filter;

}

cpl_border_mode _cplBorderConverter(char * border_str){

    cpl_border_mode border;
 
    if(strcmp(border_str, "filter") == 0){
        border = CPL_BORDER_FILTER;
    } else if(strcmp(border_str, "zero") == 0){
        border = CPL_BORDER_ZERO;
    } else if(strcmp(border_str, "crop") == 0){
        border = CPL_BORDER_CROP;
    } else if(strcmp(border_str, "nop") == 0){
        border = CPL_BORDER_NOP;
    } else if(strcmp(border_str, "copy") == 0){
        border = CPL_BORDER_COPY;
    } else {
        PyErr_Format(PyExc_ValueError, "Undefined border : %s", border_str);
        return -1;}
    return border;

}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// Module functions

//A function to test stuff. Not for public use.
static PyObject* test(PyObject* self, PyObject* args){

    PyObject * arr1, *arr2;
    //char * f_name;
    cpl_filter_mode filter;


    if (!PyArg_ParseTuple(args, "OOO&", &arr1, &arr2, &_cplFilterConverter, &filter)){
        return NULL;
	}

    //PyObject *mod, *func, *tuple, *res;
    PyObject * res;

    printf("%d\n", filter);

    res = PyNumber_Add(arr1, arr2);

    return res;

}

static PyObject* globals(PyObject* self, PyObject* args, PyObject* keywds){

    int report = 0;

    static char *kwlist[] = {"gain", "rn_adu",
                            "error_method", "error_image",
                            "report",
                            NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ddipp", kwlist, &GAIN, &RN_ADU,
                                                                &ERROR_METHOD, 
                                                                &ERROR_IMAGE, &report))
        return NULL;

    if (report){
        printf("GAIN : %.3f (e-/ADU)\n", GAIN);
        printf("RN_ADU : %.3f (ADU rms)\n", RN_ADU);
        printf("ERROR_METHOD : %d\n", ERROR_METHOD);
        printf("ERROR_IMAGE : %d\n", ERROR_IMAGE);
    }
    Py_RETURN_NONE;
}

static PyObject* bpm_fit_compute(PyObject* self, PyObject* args, PyObject* keywds){

	PyObject * arr_in = NULL, *exptime = NULL;
    PyObject *arr_out = NULL;
    npy_intp *len, *len2;
	cpl_image *c_im = NULL;
    cpl_vector * c_vec = NULL;
    hdrl_imagelist *imlist = NULL;
    hdrl_parameter *params = NULL;
    int res;
    
    int degree=1; 
    double pval=10., rel_chi_low=3., rel_chi_high=10., rel_coef_low=3., rel_coef_high=10.;
    char * method = "pval";

    static char *kwlist[] = {"image", "exptime",
                            "method", "degree", "pval", 
                            "rel_chi_low", "rel_chi_high",
                            "rel_coef_low", "rel_coef_high",
                            NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|siddddd", kwlist, &arr_in, &exptime, &method,
                                                                    &degree, &pval, 
                                                                    &rel_chi_low, &rel_chi_high,
                                                                    &rel_coef_low, &rel_coef_high)){
        return NULL;
    }

	cpl_init(CPL_INIT_DEFAULT);

    if(strcmp(method, "pval") == 0){
        params = hdrl_bpm_fit_parameter_create_pval(degree, pval);
	} else if (strcmp(method, "rel_chi") == 0) {
		params = hdrl_bpm_fit_parameter_create_rel_chi(degree, rel_chi_low, rel_chi_high);
	} else if (strcmp(method, "rel_coef") == 0) {
		params = hdrl_bpm_fit_parameter_create_rel_coef(degree, rel_coef_low, rel_coef_high);
	} else {
        PyErr_Format(PyExc_ValueError, "Undefined method : %s", method);
        goto except;
    }
    if(params == NULL){
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    len = PyArray_SHAPE((PyArrayObject*) arr_in);
    len2 = PyArray_SHAPE((PyArrayObject*) exptime);    
    PyObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    imlist = _hdrl_imagelist_from_numpy2d((PyArrayObject*)tmp, len);
    Py_DECREF(tmp);
    tmp = PyArray_FROM_OTF(exptime, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    c_vec = _cpl_vector_from_numpy1d((PyArrayObject *) tmp, len2);
    Py_DECREF(tmp);
    
    //Py_BEGIN_ALLOW_THREADS
    res = hdrl_bpm_fit_compute(params, imlist, c_vec, &c_im);
    //Py_END_ALLOW_THREADS
    if(res > 0){
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    arr_out = PyArray_SimpleNew(2, len, NPY_INT);
    memcpy((int*) PyArray_DATA((PyArrayObject *) arr_out), (int*) cpl_image_get_data(c_im), sizeof(int)*len[0]*len[1]);

    goto finally;
except:
    arr_out = NULL;
finally:
    if(imlist) hdrl_imagelist_delete(imlist);
    if(params) hdrl_parameter_delete(params);
    if(c_vec) cpl_vector_delete(c_vec);
    if(c_im) cpl_image_delete(c_im);
	cpl_end();
    return PyArray_Return((PyArrayObject *) arr_out);

}


static PyObject* bpm_2d_compute(PyObject* self, PyObject* args, PyObject* keywds){
  
    PyObject *arr_in = NULL;
    PyObject *arr_out = NULL;
    npy_intp * len;

    hdrl_image * im;
    hdrl_parameter * params;
    cpl_mask * mask;

    cpl_filter_mode filter;
    cpl_border_mode border;
    char * filter_str = "median";
    char * border_str = "filter";
    char * method = "filter";
    double kappa_low=5., kappa_high=10.; 
    int maxiter=10;
    int steps_x=10., steps_y=10.;
	int filter_size_x=10, filter_size_y=10;
	int order_x=2, order_y=2;
	int smooth_x=5, smooth_y=5;

    static char *kwlist[] = {"image", "method", "filter", "border",
                            "kappa_low", "kappa_high", "maxiter", 
                            "steps_x", "steps_y",
                            "filter_size_x", "filter_size_y",
                            "order_x", "order_y",
                            "smooth_x", "smooth_y",
                            NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|sssddiiiiiiiii", kwlist, &arr_in, &method,
                                                                    &filter_str, &border_str,
                                                                    &kappa_low, &kappa_high, &maxiter,
                                                                    &steps_x, &steps_y,
                                                                    &filter_size_x, &filter_size_y,
                                                                    &order_x, &order_y,
                                                                    &smooth_x, &smooth_y))
        return NULL;

    len = PyArray_SHAPE((PyArrayObject*)arr_in);

    filter = _cplFilterConverter(filter_str);
    border = _cplBorderConverter(border_str);
    if((filter == -1)|(border == -1)){
        return NULL;
    }

	cpl_init(CPL_INIT_DEFAULT);

    if (strcmp(method, "filter") == 0) {
        params = hdrl_bpm_2d_parameter_create_filtersmooth(kappa_low, kappa_high, maxiter, filter, 
                                                            border, smooth_x, smooth_y);
    } else if (strcmp(method, "legendre") == 0){
        params = hdrl_bpm_2d_parameter_create_legendresmooth(kappa_low, kappa_high, maxiter,
															steps_x, steps_y,
															filter_size_x, filter_size_y,
															order_x, order_y);
    } else {
        PyErr_Format(PyExc_ValueError, "Undefined method : %s", method);
        goto except;
    }
    if(params == NULL){
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    PyObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
	im = _hdrl_image_from_numpy2d((PyArrayObject*)tmp, len);
    Py_DECREF(tmp);

    //Py_BEGIN_ALLOW_THREADS
    mask = hdrl_bpm_2d_compute(im, params);
    //Py_END_ALLOW_THREADS
    if (mask == NULL) {
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    arr_out = PyArray_SimpleNew(2, len, NPY_UBYTE);
    memcpy((cpl_binary*) PyArray_DATA((PyArrayObject *)arr_out), (cpl_binary *) cpl_mask_get_data(mask), sizeof(cpl_binary)*len[0]*len[1]);

    goto finally;

except:
    arr_out = NULL;
finally:
    if(im) {hdrl_image_delete(im);}
    if(params) {hdrl_parameter_delete(params);}
    if(mask) {cpl_mask_delete(mask);}
    cpl_end();
    return PyArray_Return((PyArrayObject*)arr_out);
}

static PyObject* bpm_3d_compute(PyObject* self, PyObject* args, PyObject* keywds){

	PyObject * arr_in = NULL;
    PyObject *arr_out = NULL;
    npy_intp *len;
    cpl_imagelist * result = NULL;
    hdrl_imagelist *imlist = NULL;
    hdrl_parameter *params = NULL;
    int offset = 0;
    
    double kappa_low=3., kappa_high=10.;
    char * method = "relative";

    static char *kwlist[] = {"image", "method", 
                            "kappa_low", "kappa_high", "offset",
                            NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|sddp", kwlist, &arr_in, &method,
                                                                    &kappa_low, &kappa_high, &offset)){
        return NULL;
    }

	cpl_init(CPL_INIT_DEFAULT);

    if (strcmp(method, "relative") == 0){
        params = hdrl_bpm_3d_parameter_create(kappa_low, kappa_high, HDRL_BPM_3D_THRESHOLD_RELATIVE);
    } else if (strcmp(method, "absolute") == 0){
        params = hdrl_bpm_3d_parameter_create(kappa_low, kappa_high, HDRL_BPM_3D_THRESHOLD_ABSOLUTE);
    } else if (strcmp(method, "error") == 0){
        params = hdrl_bpm_3d_parameter_create(kappa_low, kappa_high, HDRL_BPM_3D_THRESHOLD_ERROR);
    } else {
        PyErr_Format(PyExc_ValueError, "Undefined method : %s", method);
        goto except;
    }
    if(params == NULL){
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    

    len = PyArray_SHAPE((PyArrayObject*) arr_in);
    
    PyObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    imlist = _hdrl_imagelist_from_numpy2d((PyArrayObject*)tmp, len);
    Py_DECREF(tmp);

    
    if (offset){
        hdrl_image * tmp_im = hdrl_imagelist_get(imlist, 0);
        hdrl_imagelist * tmp_list = hdrl_imagelist_duplicate(imlist);
        hdrl_imagelist_sub_image(tmp_list, tmp_im);
        cpl_size n_images = hdrl_imagelist_get_size(tmp_list);
        hdrl_value median;
        for (int i=0; i<n_images; i++){
            tmp_im = hdrl_imagelist_get(tmp_list,i);
            median = hdrl_image_get_median(tmp_im);
            tmp_im = hdrl_imagelist_get(imlist, i);
            hdrl_image_sub_scalar(tmp_im, median);
        }
        hdrl_imagelist_delete(tmp_list);
    }

    
    result = hdrl_bpm_3d_compute(imlist, params);
    if (!result){
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    arr_out = _numpy2d_from_cpl_imagelist_INT(result, len);

    goto finally;
except:
    arr_out = NULL;
finally:
    if(result) cpl_imagelist_delete(result);
    if(imlist) hdrl_imagelist_delete(imlist);
    if(params) hdrl_parameter_delete(params);
	cpl_end();
    return PyArray_Return((PyArrayObject *) arr_out);

}

static PyObject* fit_polynomial_imagelist(PyObject* self, PyObject* args, PyObject* keywds){

	PyObject * arr_in = NULL, *exptime = NULL;
    PyObject *n_coef = NULL, *n_chi2 = NULL, *n_dof = NULL;
    npy_intp *len, *len2;
	cpl_image *c_im = NULL;
    cpl_vector * c_vec = NULL;
    hdrl_imagelist *imlist = NULL;
    
    hdrl_imagelist * coef;
    cpl_image * chi2;
    cpl_image * dof;
    
    int degree=1; 
    int res;

    static char *kwlist[] = {"image", "exptime", "degree", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOi", kwlist, &arr_in, &exptime, &degree)){
        return NULL;
    }

	cpl_init(CPL_INIT_DEFAULT);


    len = PyArray_SHAPE((PyArrayObject*) arr_in);
    len2 = PyArray_SHAPE((PyArrayObject*) exptime);    
    PyObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    imlist = _hdrl_imagelist_from_numpy2d((PyArrayObject*)tmp, len);
    Py_DECREF(tmp);
    tmp = PyArray_FROM_OTF(exptime, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    c_vec = _cpl_vector_from_numpy1d((PyArrayObject *) tmp, len2);
    Py_DECREF(tmp);
    
    //Py_BEGIN_ALLOW_THREADS
    res = hdrl_fit_polynomial_imagelist(imlist, c_vec, degree, &coef, &chi2, &dof);
    //Py_END_ALLOW_THREADS
    if(res > 0){
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    n_chi2 = _numpy2d_from_cpl_image(chi2, len);
    n_dof = _numpy2d_from_cpl_image(dof, len);
    npy_intp dims[3] = {len[0], len[1], degree+1};
    n_coef = _numpy2d_from_hdrl_imagelist(coef, dims);

    //arr_out = PyArray_SimpleNew(2, len, NPY_INT);
    //memcpy((int*) PyArray_DATA((PyArrayObject *) arr_out), (int*) cpl_image_get_data(c_im), sizeof(int)*len[0]*len[1]);

    goto finally;
except:
    //arr_out = NULL;
finally:
    if(imlist) hdrl_imagelist_delete(imlist);
    //if(params) hdrl_parameter_delete(params);
    if(c_vec) cpl_vector_delete(c_vec);
    if(coef) hdrl_imagelist_delete(coef);
    if(chi2) cpl_image_delete(chi2);
    if(dof) cpl_image_delete(dof);
	cpl_end();
    //return PyArray_Return((PyArrayObject *) n_coef);
    return Py_BuildValue("OOO", 
            PyArray_Return((PyArrayObject *) n_chi2), 
            PyArray_Return((PyArrayObject *) n_dof),
            PyArray_Return((PyArrayObject *) n_coef));

}



static PyObject* bpm_interpolate(PyObject* self, PyObject* args, PyObject* keywds){
  
    PyObject *arr_in, *arr2_in;
    PyArrayObject *arr_out;
    npy_intp * len, * len2;

    cpl_image * im;
    cpl_mask * mask;

    static char *kwlist[] = {"image", "bp_mask", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO", kwlist, &arr_in, &arr2_in))
    {return NULL;}

	cpl_init(CPL_INIT_DEFAULT);

    len = PyArray_SHAPE((PyArrayObject*)arr_in);
    len2 = PyArray_SHAPE((PyArrayObject*)arr2_in);
    if ((len[0] != len2[0]) || (len[1] != len2[1])){
        PyErr_Format(PyExc_ValueError, "Shape of the image (%ld,%ld) and mask (%ld,%ld) does not match",
                                        len[0], len[1], len2[0], len2[1]);
        goto except;
    }

    PyObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
	im = _cpl_image_from_numpy2d((PyArrayObject*)tmp, len);
    Py_DECREF(tmp);

    //explicit casting crashes without NPY_ARRAY_FORCECAST.
    PyObject *tmp2 = PyArray_FROM_OTF(arr2_in, NPY_UBYTE, NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST);
	mask = _cpl_mask_from_numpy2d((PyArrayObject*)tmp2, len);
    Py_DECREF(tmp2);

    int res=0;

    mask = cpl_image_set_bpm(im, mask);

    //Py_BEGIN_ALLOW_THREADS
    res = cpl_detector_interpolate_rejected(im);
    //Py_END_ALLOW_THREADS
    if (res != 0) {
        PyErr_Format(PyExc_ValueError, "CPL Error (%d): %s", cpl_error_get_code(), cpl_error_get_message());
        goto except;
    }

    arr_out =  _numpy2d_from_cpl_image(im, len);

    goto finally;
except:
    arr_out = NULL;
finally:
    if(im) {cpl_image_delete(im);}
    cpl_end();
    return PyArray_Return(arr_out);
}


static PyObject* lacosmic_edgedetect(PyObject* self, PyObject* args, PyObject* keywds){
  
    PyObject *arr_in;
    PyObject *arr_out;
    npy_intp * len;

    hdrl_image * im;
    hdrl_parameter * params;
    cpl_mask * mask;

    int maxiter=5;
    double sigma_lim=20., f_lim=2.;

    static char *kwlist[] = {"image", "sigma_lim",
                            "f_lim", "maxiter", 
                            NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|ddi", kwlist, &arr_in,
                                            &sigma_lim, &f_lim, &maxiter))
        return NULL;

    len = PyArray_SHAPE((PyArrayObject*)arr_in);

	cpl_init(CPL_INIT_DEFAULT);

    params = hdrl_lacosmic_parameter_create(sigma_lim, f_lim, maxiter);
    if(params == NULL){
        //HDRL_simplePrint(cpl_error_get_message());
        printf("%02d %s\n", cpl_error_get_code(), cpl_error_get_message());
        Py_RETURN_NONE;
    }

    PyObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_C_CONTIGUOUS);
	im = _hdrl_image_from_numpy2d((PyArrayObject*)tmp, len);
    Py_DECREF(tmp);

    //Py_BEGIN_ALLOW_THREADS
    mask = hdrl_lacosmic_edgedetect(im, params);
    //Py_END_ALLOW_THREADS
    if (mask == NULL) {
        //HDRL_simplePrint("Warning: lacosmic_edgedetect failed. Returned None.");
        printf("Warning: lacosmic_edgedetect failed. Returned None.\n");
        Py_RETURN_NONE;
    }

    arr_out = PyArray_SimpleNew(2, len, NPY_UBYTE);
    memcpy((cpl_binary*) PyArray_DATA((PyArrayObject *)arr_out), (cpl_binary*) cpl_mask_get_data(mask), sizeof(cpl_binary)*len[0]*len[1]);

    hdrl_image_delete(im);
    hdrl_parameter_delete(params);
    cpl_mask_delete(mask);

    cpl_end();

    return PyArray_Return((PyArrayObject *)arr_out);
}

PyObject * compute_strehl(PyObject* self, PyObject* args){

    PyObject * arr_in;
    hdrl_image * im;
    hdrl_strehl_result strehl;
    npy_intp * len;

    double wave, m1_rad, m2_rad;
	double pixsc_x, pixsc_y;
	double flux_r, bkg_low, bkg_high;

    if (!PyArg_ParseTuple(args, "Odddddddd", &arr_in, 
                                &wave, &m1_rad, &m2_rad,
                                &pixsc_x, &pixsc_y,
                                &flux_r, &bkg_low, &bkg_high)){
        return NULL;
	}

    len = PyArray_SHAPE((PyArrayObject*)arr_in);
	cpl_init(CPL_INIT_DEFAULT);

    PyObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_C_CONTIGUOUS);
    im = _hdrl_image_from_numpy2d((PyArrayObject*)tmp, len);
    Py_DECREF(tmp);

    hdrl_parameter * params = hdrl_strehl_parameter_create(wave, m1_rad, m2_rad,
							pixsc_x, pixsc_y,
							flux_r, bkg_low, bkg_high);

    //if (params == NULL)
        //HDRL_simplePrint("Oops\n");

    //Py_BEGIN_ALLOW_THREADS
    strehl = hdrl_strehl_compute(im, params);
    //Py_END_ALLOW_THREADS

    hdrl_image_delete(im);
    hdrl_parameter_delete(params);

    cpl_end();

    return strehlResult_initFromOriginal(strehl);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
// Module definitions.

static PyMethodDef hdrlMethods[] =
{
    {"bpm_2d_compute", (PyCFunction) bpm_2d_compute, METH_VARARGS | METH_KEYWORDS,
        bpm_2d_compute__doc__},
    {"bpm_3d_compute", (PyCFunction) bpm_3d_compute, METH_VARARGS | METH_KEYWORDS,
        "asdasd"},
    {"fit_polynomial_imagelist", (PyCFunction) fit_polynomial_imagelist, METH_VARARGS | METH_KEYWORDS,
        "asdasdasda"},
    {"lacosmic_edgedetect", (PyCFunction) lacosmic_edgedetect, METH_VARARGS | METH_KEYWORDS,
        lacosmic_edgedetect__doc__},
    {"bpm_interpolate", (PyCFunction) bpm_interpolate, METH_VARARGS | METH_KEYWORDS,
        bpm_interpolate__doc__},
    {"compute_strehl", compute_strehl, METH_VARARGS,
        compute_strehl__doc__},
    {"bpm_fit_compute", (PyCFunction) bpm_fit_compute, METH_VARARGS | METH_KEYWORDS,
        bpm_fit_compute__doc__},
    {"globals", (PyCFunction) globals, METH_VARARGS | METH_KEYWORDS,
        "globals"},
    {"test", test, METH_VARARGS,
        "evaluate the cosine on a numpy array"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef hdrlModule = {
    PyModuleDef_HEAD_INIT,
    "HDRL2",   /* name of module */
    mod__doc__, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    hdrlMethods
};


PyMODINIT_FUNC PyInit_HDRL2(void) {
    
    PyObject *module;

    module = PyModule_Create(&hdrlModule);
    if(module==NULL) return NULL;
    import_array();

    if (PyType_Ready(&hdrlValueType) < 0)
        return NULL;
    Py_INCREF(&hdrlValueType);
    if (PyModule_AddObject(module, "hdrl_value", (PyObject *) &hdrlValueType) < 0) {
        Py_DECREF(&hdrlValueType);
        Py_DECREF(module);
        return NULL;
    }

    if (PyType_Ready(&strehlResultType) < 0)
        return NULL;
    Py_INCREF(&strehlResultType);
    if (PyModule_AddObject(module, "strehl_result", (PyObject *) &strehlResultType) < 0) {
        Py_DECREF(&strehlResultType);
        Py_DECREF(module);
        return NULL;
    }

    //PyModule_AddObject(module, NAME_ERROR_IMAGE, Py_BuildValue("O", Py_True)); 
    //PyModule_AddObject(module, NAME_ERROR_METHOD, Py_BuildValue("i", ERROR_METHOD)); 
    //PyModule_AddObject(module, NAME_RN_ADU, Py_BuildValue("d", RN_ADU)); 
    //PyModule_AddObject(module, NAME_GAIN, Py_BuildValue("d", GAIN)); 

    if (PyErr_Occurred()) return NULL;
    return module;
}

