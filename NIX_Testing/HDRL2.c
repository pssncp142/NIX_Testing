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
PyDoc_STRVAR(mod_doc, "Some utility functions from HDRL to be used with numpy arrays.\n");
PyDoc_STRVAR(bpm_fit_compute_doc, 
"**bpm_fit_compute(image, exptime, [method, degree, pval, rel_chi_low, rel_chi_high, rel_coef_low, rel_coef_high])**\n"
"\n"
"A python wrapper around 'hdrl_bpm_fit_compute' function. 'hdrl_bpm_fit_parameter' is created inside the wrapper "
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
"method : str\n"
"\tMethod used to determine bad pixels. If given method is none of the below"
"function returns *Py_None*. Accepted values are (default='pval'):\n\n"
"\t* 'pval'    : threshold on the p-values.\n"
"\t* 'rel_chi' : sigma cut on the relative chi values.\n"
"\t* 'rel_coef': sigma cut on the relative coefficients of polynomial fit.\n"
"degree : int, optional\n"
"\tDegree of polynomial to fit. (default=1)\n"
"pval : double, optional\n"
"\tThreshold in p-values. Enabled with *pval* (default=10.)\n"
"rel_chi_low : double, optional\n"
"\tLow threshold in relative chi values. Enabled with *rel_chi* (default=5.)\n"
"rel_chi_low : double, optional\n"
"\tLow threshold in relative chi values. Enabled with *rel_chi* (default=10.)\n"
"rel_coef_low : double, optional\n"
"\tLow threshold in relative coefficients. Enabled with *rel_coef* (default=5.)\n"
"rel_coef_low : double, optional\n"
"\tLow threshold in relative coefficients. Enabled with *rel_coef* (default=10.)\n"
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


/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

void HDRL_simplePrint(char * text){
    
    PyObject * sys = PyImport_ImportModule( "sys");
    PyObject * s_out = PyObject_GetAttrString(sys, "stdout"); 
    PyObject * result = PyObject_CallMethod(s_out, "write", "s", strcat(text, "\n")); 
    Py_DECREF(result);
    Py_DECREF(s_out);
    Py_DECREF(sys);

}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

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
    char *out = malloc(100*sizeof(char));
    sprintf(out, "(%.3f +- %.3f)", self->data, self->error);
    PyObject * res = Py_BuildValue("s", out);
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

cpl_image * _cpl_image_from_numpy2d(PyArrayObject* arr, npy_intp * len){

    printf("%d %d\n", len[0], len[1]);
    //PyArrayObject * tmp = PyArray_GETCONTIGUOUS(arr);
	cpl_image * c_im = cpl_image_new(len[0], len[1], CPL_TYPE_DOUBLE);
	memcpy((double *) cpl_image_get_data(c_im), 
			(double *) PyArray_DATA(arr), 
			sizeof(double)*len[0]*len[1]);
	return c_im;

}

cpl_vector * _cpl_vector_from_numpy1d(PyArrayObject* arr, npy_intp * len){

    //PyArrayObject * tmp = PyArray_GETCONTIGUOUS(arr);
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

hdrl_image * _hdrl_image_from_numpy2d(PyArrayObject* arr, npy_intp *len){

    cpl_image * c_im = _cpl_image_from_numpy2d(arr, len);
    cpl_image * c_err;

    bool ERROR_IMAGE = true;
    int ERROR_METHOD = 0;
    double GAIN = 5.7;
    double RN_ADU = 4.5;

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

	image_hdrl = hdrl_image_create(c_im, c_err);

    cpl_image_delete(c_im);
    cpl_image_delete(c_err);

    return image_hdrl;

}

hdrl_imagelist * _hdrl_imagelist_from_numpy2d(PyArrayObject* arr, npy_intp *len){

    hdrl_image * im_tmp;
    hdrl_imagelist * imlist;
    PyArrayObject * arr_ndx, *arr_tmp;
    int * data_tmp;
    double * data;
    npy_intp test[NPY_MAXDIMS];
    test[0] = 1;

    arr_ndx = PyArray_SimpleNew(1, test, NPY_INT);
    Py_INCREF(arr_ndx);
    data_tmp = (int*) PyArray_DATA(arr_ndx);

    imlist = hdrl_imagelist_new();

    for (int i=0; i<len[2]; i++){
        data_tmp[0] = i;
        arr_tmp = PyArray_TakeFrom(arr, arr_ndx, 2, NULL, NPY_RAISE);
        im_tmp = _hdrl_image_from_numpy2d(arr_tmp, len);
        hdrl_imagelist_set(imlist, im_tmp, i);
        Py_DECREF(arr_tmp);
    }

    Py_DECREF(arr_ndx);
    Py_DECREF(arr_ndx);

    return imlist;

}


hdrl_parameter * hdrl_bpm_2d_parameter_create(cpl_filter_mode filter, cpl_border_mode border,
											double kappa_low, double kappa_high, int maxiter,
											int steps_x, int steps_y,
											int filter_size_x, int filter_size_y,
											int order_x, int order_y,
											int smooth_x, int smooth_y){

	if (order_x == -1){
		return hdrl_bpm_2d_parameter_create_filtersmooth(kappa_low, kappa_high, maxiter, filter, border, smooth_x, smooth_y);
	} else {
		return hdrl_bpm_2d_parameter_create_legendresmooth(kappa_low, kappa_high, maxiter,
															steps_x, steps_y,
															filter_size_x, filter_size_y,
															order_x, order_y);

	}

}

static PyObject* test(PyObject* self, PyObject* args){

    return hdrlValue_initFromC(4, 3);

}

static PyObject* bpm_fit_compute(PyObject* self, PyObject* args, PyObject* keywds){
  
	PyObject * arr_in, *exptime;
    PyArrayObject *arr_out;
    double factor;
    npy_intp *len, *len2;
	cpl_image *c_im;
    cpl_vector * c_vec;
    cpl_mask * mask;
    hdrl_imagelist *imlist;
    hdrl_parameter * params;
    hdrl_image * im;

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
                                                                    &rel_coef_low, &rel_coef_high))

	len = PyArray_SHAPE(arr_in);
    len2 = PyArray_SHAPE(exptime);
	cpl_init(CPL_INIT_DEFAULT);

    if(strcmp(method, "pval") == 0){
        params = hdrl_bpm_fit_parameter_create_pval(degree, pval);
	} else if (strcmp(method, "rel_chi") == 0) {
		params = hdrl_bpm_fit_parameter_create_rel_chi(degree, rel_chi_low, rel_chi_high);
	} else if (strcmp(method, "rel_coef") == 0) {
		params = hdrl_bpm_fit_parameter_create_rel_coef(degree, rel_coef_low, rel_coef_high);
	} else {
        printf("Undefined Method.\n");
        cpl_end();
        Py_RETURN_NONE;
    }
    
    PyArrayObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_C_CONTIGUOUS);
    len = PyArray_SHAPE(tmp);
    len2 = PyArray_SHAPE(exptime);
    imlist = _hdrl_imagelist_from_numpy2d(tmp, len);
    Py_DECREF(tmp);

    tmp = PyArray_FROM_OTF(exptime, NPY_FLOAT64, NPY_ARRAY_C_CONTIGUOUS);
    c_vec = _cpl_vector_from_numpy1d(tmp, len2);
    Py_DECREF(tmp);

    int res;
    Py_BEGIN_ALLOW_THREADS
    res = hdrl_bpm_fit_compute(params, imlist, c_vec, &c_im);
    Py_END_ALLOW_THREADS
    printf("fit_compute res : %d\n", res);

    arr_out = PyArray_SimpleNew(2, len, NPY_INT);
    memcpy((int*) PyArray_DATA(arr_out), (int*) cpl_image_get_data(c_im), sizeof(int)*len[0]*len[1]);

    //cpl_binary * tmp1 = (cpl_binary*) PyArray_DATA(arr_out);
    //int * tmp2 = (int*) cpl_image_get_data(c_im);
    //for (int i=0; i<len[0]*len[1]; i++){
    //    tmp1[i] = (cpl_binary) tmp2[i];
   // }
 
    hdrl_imagelist_delete(imlist);
    hdrl_parameter_delete(params);
    cpl_image_delete(c_im);
 
	cpl_end();

    return PyArray_Return(arr_out);

}


static PyObject* bpm_2d_compute(PyObject* self, PyObject* args, PyObject* keywds){
  
    PyObject *arr_in;
    PyArrayObject *arr_out;
    npy_intp * len;

    hdrl_image * im;
    hdrl_parameter * params;
    cpl_mask * mask;

    cpl_filter_mode filter = CPL_FILTER_MEDIAN;
    cpl_border_mode border = CPL_BORDER_FILTER;
    char * method = "filter";
    double kappa_low=5., kappa_high=10.; 
    int maxiter=10;
    int steps_x=10., steps_y=10.;
	int filter_size_x=10, filter_size_y=10;
	int order_x=2, order_y=2;
	int smooth_x=5, smooth_y=5;

    static char *kwlist[] = {"image", "method",
                            "kappa_low", "kappa_high", "maxiter", 
                            "steps_x", "steps_y",
                            "filter_size_x", "filter_size_y",
                            "order_x", "order_y",
                            "smooth_x", "smooth_y",
                            NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|sddiiiiiiiii", kwlist, &arr_in, &method,
                                                                    &kappa_low, &kappa_high, &maxiter,
                                                                    &steps_x, &steps_y,
                                                                    &filter_size_x, &filter_size_y,
                                                                    &order_x, &order_y,
                                                                    &smooth_x, &smooth_y))
        return NULL;

    len = PyArray_SHAPE(arr_in);
    printf("%d %d\n", len[0], len[1]);

	cpl_init(CPL_INIT_DEFAULT);

    if (strcmp(method, "filter") == 0) 
        order_x = -1;

    params = hdrl_bpm_2d_parameter_create(filter, border,
										kappa_low, kappa_high, maxiter,
										steps_x, steps_y,
										filter_size_x, filter_size_y,
										order_x, order_y,
										smooth_x, smooth_y);

    if(params == NULL){
        HDRL_simplePrint(cpl_error_get_message());
        printf("%02d %s\n", cpl_error_get_code(), cpl_error_get_message());
        Py_RETURN_NONE;
    }

    PyArrayObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_C_CONTIGUOUS);
	im = _hdrl_image_from_numpy2d(tmp, len);
    Py_DECREF(tmp);

    Py_BEGIN_ALLOW_THREADS
    mask = hdrl_bpm_2d_compute(im, params);
    Py_END_ALLOW_THREADS
    if (mask == NULL) {
        HDRL_simplePrint("Warning: bpm_2d_compute failed. Returned None.");
        printf("Warning: bpm_2d_compute failed. Returned None.\n");
        Py_RETURN_NONE;
    }

    arr_out = PyArray_SimpleNew(2, len, NPY_UBYTE);
    memcpy((cpl_binary*) PyArray_DATA(arr_out), (cpl_binary*) cpl_mask_get_data(mask), sizeof(cpl_binary)*len[0]*len[1]);

    hdrl_image_delete(im);
    hdrl_parameter_delete(params);
    cpl_mask_delete(mask);

    cpl_end();

    return PyArray_Return(arr_out);
}

PyObject * compute_strehl(PyObject* self, PyObject* args){

    PyObject * arr_in;
    PyObject * result;
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

    len = PyArray_SHAPE(arr_in);
	cpl_init(CPL_INIT_DEFAULT);

    PyArrayObject * tmp = PyArray_FROM_OTF(arr_in, NPY_FLOAT64, NPY_ARRAY_C_CONTIGUOUS);
    im = _hdrl_image_from_numpy2d(tmp, len);
    Py_DECREF(tmp);

    hdrl_parameter * params = hdrl_strehl_parameter_create(wave, m1_rad, m2_rad,
							pixsc_x, pixsc_y,
							flux_r, bkg_low, bkg_high);

    if (params == NULL)
        HDRL_simplePrint("Oops\n");

    Py_BEGIN_ALLOW_THREADS
    strehl = hdrl_strehl_compute(im, params);
    Py_END_ALLOW_THREADS

    hdrl_image_delete(im);
    hdrl_parameter_delete(params);

    cpl_end();

    return strehlResult_initFromOriginal(strehl);
}


static PyMethodDef hdrlMethods[] =
{
    {"bpm_2d_compute", (PyCFunction) bpm_2d_compute, METH_VARARGS | METH_KEYWORDS,
        "evaluate the cosine on a numpy array"},
    {"compute_strehl", compute_strehl, METH_VARARGS,
        "Calculate Strehl ratio from a given image"},
    {"bpm_fit_compute", (PyCFunction) bpm_fit_compute, METH_VARARGS | METH_KEYWORDS,
        bpm_fit_compute_doc},
    {"test", test, METH_NOARGS,
        "evaluate the cosine on a numpy array"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef hdrlModule = {
    PyModuleDef_HEAD_INIT,
    "HDRL2",   /* name of module */
    mod_doc, /* module documentation, may be NULL */
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

    if (PyErr_Occurred()) return NULL;
    return module;
}

/*
PyMODINIT_FUNC inittest(void) {
    PyObject *module;
    module = Py_InitModule("HDRL2", hdrlMethods);
    if(module==NULL) return;
    import_array();
    return;
}*/
