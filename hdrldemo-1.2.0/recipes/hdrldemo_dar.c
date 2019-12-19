/*
 * This file is part of the HDRLDEMO pipeline
 * Copyright (C) 2017 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                    Includes
 -----------------------------------------------------------------------------*/

#include "hdrl.h"

#include "hdrldemo_utils.h"
#include "hdrldemo_dfs.h"

#include <cpl.h>

/*----------------------------------------------------------------------------*/
/**
  @defgroup hdrldemo_dar - Differential atmospheric refraction (DAR) recipe
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

/* Clean parameters allocated in the recipe */
static void hdrldemo_dar_clean_recipe_parameters(
	cpl_propertylist *header_raw,
	cpl_wcs *wcs, cpl_imagelist *imlistOrig, cpl_imagelist *imlist,
	cpl_image* img_orig_collapsed, cpl_image* img_correct_collapsed,
	cpl_vector *lambdaIn, cpl_vector *xShift, cpl_vector *yShift,
	cpl_vector *xShiftErr, cpl_vector *yShiftErr, cpl_table *tb);

/* Verify parameters in the propertylist */
static cpl_error_code hdrldemo_dar_verify_propertylist_parameters(
	const cpl_propertylist *header);

/* Get exptime, try to found the different ways to found that value in the header */
static double hdrldemo_dar_get_exptime(
	const cpl_propertylist *header);

/* Calculate the correct value of parallactic angle based in p1 and p2 */
static double hdrldemo_dar_parangle(
	double p1, double p2);

/* Create columns and size in the output table */
static cpl_table * hdrldemo_dar_create_output_table(
	const cpl_size size);

/* Add Wavelenght of planes, shift values and errors to the output cpl_table tb */
static cpl_error_code hdrldemo_dar_add_table_wavelenght_and_shift(
	cpl_table        *tb,
	const cpl_vector *lambda,
	const cpl_vector *xShift,
	const cpl_vector *yShift,
	const cpl_vector *xShiftErr,
	const cpl_vector *yShiftErr);

/* Cataloge a imagelist adding to cpl_table tb Isophotal flux, coordX-Y and errors */
static cpl_error_code hdrldemo_dar_catalogue(
	cpl_boolean   orig,
	cpl_imagelist *imlist,
	cpl_wcs       *wcs,
	double		  maxRatioNaNsInImage,
	cpl_table     *tb);

/* Save output data from original and shifted data to disk */
static void hdrldemo_dar_save_disk(
	cpl_propertylist        *header_raw,
	cpl_table 	        	*tb,
	cpl_imagelist           *imlistOrig,
	cpl_imagelist           *imlistCorrect,
	cpl_image               *imOrigCollapsed,
	cpl_image               *imCorrectCollapsed,
	const cpl_parameterlist *parlist,
	const cpl_propertylist 	*app_list,
	cpl_frameset 			*frameset,
	cpl_frame 				*raw_frame);

/* Functions for add Statistical QC values: Min, Max, Mean, Median, StDEV, MAD */
static cpl_error_code hdrldemo_dar_add_QC_parameters(
	cpl_propertylist *qclist,
	const char       *param,
	const cpl_table  *tb,
	const char       *cname);

/* Functions for add Statistical QC value Median absolute deviation (MAD) */
static cpl_error_code hdrldemo_dar_add_QC_parameter_mad(
		cpl_propertylist *qclist,
		const char       *param,
		const cpl_table  *tb,
		const char       *cname);

/* Function for add RMS in the Distance X-Y */
static cpl_error_code hdrldemo_dar_add_QC_RMS_XY(
	cpl_propertylist *qclist,
	const char       *param,
	const cpl_table  *tb,
	const char       *cnameX,
	const char       *cnameY);

/* Check and add the value in the qclist */
static cpl_error_code hdrldemo_dar_add_property(
	cpl_propertylist *qclist,
	const char       *param,
	const char       *qcName,
	const double     value);

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_dar"

static char hdrldemo_dar_description[] =
"                                                                           \n"
"The recipe calculate and apply the differential atmospheric refraction.    \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:                     Explanation:           Required:        \n"
"  RAW                              Data                   Yes              \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                     Explanation:                            \n"
"  HDRLDEMO_DAR_ORIG                Original imagelist set                  \n"
"  HDRLDEMO_DAR_CORRECTED           Shifted  imagelist set                  \n"
"  HDRLDEMO_DAR_ORIG_COLLAPSED      Collapsed original image                \n"
"  HDRLDEMO_DAR_CORRECTED_COLLAPSED Collapsed shifted  image                \n"
"  HDRLDEMO_DAR_RESULT              Shift and error plus more QC parameters \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"                                                                           \n"
"The recipe allows to correct the athmospheric refraction in a imagelist    \n"
"as input (type SINFONI cube datasheet). For that is used the algorithm    \n"
"of Filippenko (1982). In the HDRLDEMO_DAR_RESULT you can find the different\n"
"values like the plane, the wavelengh associate and the shift (that is      \n"
"necessary to apply to each pixel in each axis for this wavelenght).        \n"
"The algorithm calculate the error propagation if you put in the input      \n"
"the error in the input parameters of the observation (i.e. parallactic     \n"
"angle, position angle, pressure, ... See more in the default parameters).  \n"
"The airmass is calculated among the header parameters using the hdrl_utils \n"
"function hdrl_airmass(). This function can use the Hardie (1962),          \n"
"Young & Irvine (1967) or Young approximation (1994). You can select which  \n"
"approximation you want with the parameter --.typeAirmassAprox .            \n"
"Once have apply the shift, the image set is collapse and catalogue adding  \n"
"more extra information to the output (like the isophotal flux and the      \n"
"positon of the first source in the image). In the catalogue you can control\n"
"if it's being apply to be in count the ratio between the good and bad      \n"
"pixels (specified in the parameter --.maxRatioNaNsInImage).                \n"
"Also, in the header of the output file, it's offered to the user           \n"
"QC parameters for each column (min, max, mean, median, stdev, ...).        \n"
"                                                                           \n"
"Please note that part of the code in hdrl_dar is parallelized.             \n"
"In order to optimize use of computing resources you should set the         \n"
"environment variable OMP_NUM_THREADS to a proper value, like (for bash),   \n"
"to use 4 cores export OMP_NUM_THREADS=4                                    \n";

/* Name of columns in the output table for the results */
static char index_cname[]                  = "Plane";
static char orig_bpImg_cname[]             = "orig_bad_pixel_image";
static char shifted_bpImg_cname[]          = "shifted_bad_pixel_image";
static char lambdaIn_cname[]               = "wavelength";
static char x_cname[]		               = "xShift";
static char y_cname[]                      = "yShift";
static char xErr_cname[]                   = "xShiftErr";
static char yErr_cname[]                   = "yShiftErr";
static char orig_cat_flux_cname[]          = "orig_cat_flux";
static char orig_cat_xcoord_cname[]        = "orig_cat_xcoord";
static char orig_cat_ycoord_cname[]        = "orig_cat_ycoord";
static char orig_cat_xcoord_err_cname[]    = "orig_cat_xcoord_err";
static char orig_cat_ycoord_err_cname[]    = "orig_cat_ycoord_err";
static char shifted_cat_flux_cname[]       = "shifted_cat_flux";
static char shifted_cat_xcoord_cname[]     = "shifted_cat_xcoord";
static char shifted_cat_ycoord_cname[]     = "shifted_cat_ycoord";
static char shifted_cat_xcoord_err_cname[] = "shifted_cat_xcoord_err";
static char shifted_cat_ycoord_err_cname[] = "shifted_cat_ycoord_err";


/* Standard CPL recipe definition */
cpl_recipe_define(	hdrldemo_dar, HDRLDEMO_BINARY_VERSION, "HDRL Group",
					PACKAGE_BUGREPORT, "2017", "HDRLDEMO - DAR",
					hdrldemo_dar_description);

/*----------------------------------------------------------------------------*/
/**
 * @brief Function needed by hdrldemo_dar_fill_parameterlist to get input parameters
 *
 * @param  self   parameterlist where you need put parameters
 * @param  param  Name of parameter
 * @param  desc   Description of parameter
 *
 * @return (int)cpl_error_code, for error codes see cpl_table_wrap_double
 *
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_fill_parameter(
		cpl_parameterlist *self, const char* param, cpl_type type,
		const char* desc, double value)
{
	cpl_errorstate prestate = cpl_errorstate_get();

	/* Set parameter */
	cpl_parameter *par = NULL;
	switch (type) {
		case CPL_TYPE_DOUBLE:
			par = cpl_parameter_new_value(param,	type, desc, RECIPE_NAME, value);
			break;
		case CPL_TYPE_INT:
			par = cpl_parameter_new_value(param,	type, desc, RECIPE_NAME, (int)value);
			break;
		default:
			return CPL_ERROR_INCOMPATIBLE_INPUT;
	}

	/* Add parameter to default parameters */
	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, param);
	cpl_parameter_disable(  par, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append(self, par);

	/* Check possible errors */
	if (!cpl_errorstate_is_equal(prestate)) {
		return CPL_ERROR_ILLEGAL_INPUT;
	}

	return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Function needed by cpl_recipe_define to fill the input parameters
 *
 * @param  self   parameterlist where you need put parameters
 *
 * @return (int)cpl_error_code
 *
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_fill_parameterlist(cpl_parameterlist *self)
{                                  
	/* Add the different default parameters to the recipe */
	cpl_errorstate prestate = cpl_errorstate_get();

    // --[VARIABLE]

	hdrldemo_dar_fill_parameter(self, "maxRatioNaNsInImage",   CPL_TYPE_DOUBLE,
			"Max ratio [%] between bad pixels (NaNs) and good pixels", 30. );

	hdrldemo_dar_fill_parameter(self, "lambdaRef-err",    CPL_TYPE_DOUBLE,
			"Error [%] of the reference lambda. If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "parang-err",       CPL_TYPE_DOUBLE,
			"Error [%] of the parallactic angle If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "posang-err",       CPL_TYPE_DOUBLE,
			"Error [%] of the position angle. If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "temp-err",         CPL_TYPE_DOUBLE,
			"Error [%] of the temperature. If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "rhum-err",         CPL_TYPE_DOUBLE,
			"Error [%] of the relative humidity. If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "pres-err",         CPL_TYPE_DOUBLE,
			"Error [%] of the pressure. If the variable is zero, the error is an absolute value",
			0.);

    /* --hdrldemo_dar.typeAirmassAprox */
    char *name = cpl_sprintf("airm.method");
    cpl_parameter *p = cpl_parameter_new_enum(name, CPL_TYPE_STRING,
            "Method of approximation to calculate airmass", name, "HARDIE_62", 3,
			"HARDIE_62", "YOUNG_IRVINE_67", "YOUNG_94");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, name);
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_free(name);
    cpl_parameterlist_append(self, p);

	hdrldemo_dar_fill_parameter(self, "airm.ra-err",      CPL_TYPE_DOUBLE,
			"Error [%] of the right ascension (using for calculate airmass). If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "airm.dec-err",     CPL_TYPE_DOUBLE,
			"Error [%] of the declination (using for calculate airmass). If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "airm.lst-err",     CPL_TYPE_DOUBLE,
			"Error [%] of the local sideral time (using for calculate airmass). If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "airm.exptime-err", CPL_TYPE_DOUBLE,
			"Error [%] of the integration time (using for calculate airmass). If the variable is zero, the error is an absolute value",
			0.);
	hdrldemo_dar_fill_parameter(self, "airm.geolat-err",  CPL_TYPE_DOUBLE,
			"Error [%] of the latitude (using for calculate airmass). If the variable is zero, the error is an absolute value",
			0.);


	/* Check possible errors */
	if (!cpl_errorstate_is_equal(prestate)) {
		  return cpl_error_set_message(cpl_func, cpl_error_get_code(),
				  "hdrldemo_dar_fill_parameterlist");
	}

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Recipe demo for check the calculus of differential atmospheric refraction (DAR)
 *
 * @param  frameset   input set of frames
 * @param  parlist    input recipe parameters
 *
 * @return (int)cpl_error_code
 *
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_dar( cpl_frameset            *frameset,
						 const cpl_parameterlist *parlist)
{
    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /* Load INPUT Data, Get the first Required frame */
    cpl_frame *in_frm = cpl_frameset_find(frameset, HDRLDEMO_RAW);
    if (in_frm == NULL) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
				"Missing RAW file");
    }

    /* Get the image set */
	cpl_msg_info(cpl_func,"Charging the image list ...");
    cpl_imagelist *imlistOrig = cpl_imagelist_load_frameset(frameset, CPL_TYPE_DOUBLE, 0, 0);
    if (imlistOrig == NULL) {
    	return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
      			"hdrldemo_dar: Image list in input frameset is NULL");
    }

	/* Get the number of images */
	cpl_size size = cpl_imagelist_get_size(imlistOrig);
	if (size <= 0) {
    	cpl_imagelist_delete(imlistOrig);
		return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
				"Number of images Zero, NULL of incompatible");
	}

    /* Get the primary and extension */
    const char       *filename   = cpl_frame_get_filename(in_frm);
    cpl_propertylist *header_raw = cpl_propertylist_load(filename, 0);
    if (header_raw == NULL) {
    	cpl_imagelist_delete(imlistOrig);
        return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                "Header is NULL");
	}

    /* Check values in the header for the recipe (header SINFONI DATA) */
    cpl_wcs *wcs = cpl_wcs_new_from_propertylist(header_raw);
    if(hdrldemo_dar_verify_propertylist_parameters(header_raw) != CPL_ERROR_NONE || wcs == NULL){
    	cpl_imagelist_delete(imlistOrig);
    	cpl_propertylist_delete(header_raw);
    	cpl_wcs_delete(wcs);
    	return cpl_error_set_message(cpl_func, cpl_error_get_code(),
    			"The header doesn't contain the expected data");
    }

	/* Duplicate the imglist for shifted */
    cpl_imagelist *imlist = cpl_imagelist_duplicate(imlistOrig);


    /********* Get input parameters in the primary header **************/
	cpl_msg_info(cpl_func,"Charging the values of the parameters inside the header ...");

    hdrl_value airmass_start = {cpl_propertylist_get_double(header_raw, "ESO TEL AIRM START"	),  0.};
    hdrl_value airmass_end   = {cpl_propertylist_get_double(header_raw, "ESO TEL AIRM END"		),  0.};
	hdrl_value airmass       = {(airmass_start.data  + airmass_end.data ) / 2.,                     0.};

    hdrl_value ra            = {cpl_propertylist_get_double(header_raw, "RA"					 ), 0.};
    hdrl_value dec           = {cpl_propertylist_get_double(header_raw, "DEC"					 ), 0.};
    hdrl_value lst           = {cpl_propertylist_get_double(header_raw, "LST"					 ), 0.};
    hdrl_value exptime       = {hdrldemo_dar_get_exptime(header_raw), 0.};
    hdrl_value geolat        = {cpl_propertylist_get_double(header_raw, "ESO TEL GEOLAT"		 ), 0.};

    double parang_start      =  cpl_propertylist_get_double(header_raw, "ESO TEL PARANG START"	 );
	double parang_end        =  cpl_propertylist_get_double(header_raw, "ESO TEL PARANG END"	 );
	hdrl_value parang  		 = {hdrldemo_dar_parangle(parang_start, parang_end), 0.};

	hdrl_value posang        = {cpl_propertylist_get_double(header_raw, "ESO ADA POSANG"		 ), 0.};

	hdrl_value temp          = {cpl_propertylist_get_double(header_raw, "ESO TEL AMBI TEMP"		 ), 0.};
	hdrl_value rhum          = {cpl_propertylist_get_double(header_raw, "ESO TEL AMBI RHUM"		 ), 0.};

	double pres_start        =  cpl_propertylist_get_double(header_raw, "ESO TEL AMBI PRES START");
	double pres_end          =  cpl_propertylist_get_double(header_raw, "ESO TEL AMBI PRES END"	 );
	hdrl_value pres          = {(pres_end + pres_start) / 2., 0.};


	/* Wavelenght reference from micron(um) to Angstrom(A) */
	double   centerLambda    =           cpl_propertylist_get_double(header_raw, "CRVAL3");
	double   planeDist       =           cpl_propertylist_get_double(header_raw, "CDELT3");
	cpl_size planeLambda     = (cpl_size)cpl_propertylist_get_double(header_raw, "CRPIX3");
	hdrl_value lambdaRef     = {centerLambda * 10000., 0.};
	double     step          =  planeDist    * 10000.;


	/*** GET ERROR PARAMETERS FOR THE RECIPE - In [%] respect base Values ***/
	cpl_msg_info(cpl_func,"Get default parameter of the recipe ...");
    const cpl_parameter *par;

	/* Max ratio allow between good and bad pixels in the image to get the catalogue */
    par = cpl_parameterlist_find_const(parlist, "maxRatioNaNsInImage");
	double maxRatioNaNsInImage = cpl_parameter_get_double(par) / 100.;


	/* Error parameters for calculate DAR in percent [frac] of the actual value */

    /*** Note: if the error data values are exactly 0., the error parameter is getting as absolute error ***/

    par = cpl_parameterlist_find_const(parlist, "lambdaRef-err");
    if (lambdaRef.data != 0.) {
    	lambdaRef.error = cpl_parameter_get_double(par) / 100. * fabs(lambdaRef.data);
    } else {
    	lambdaRef.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "parang-err");
    if (parang.data != 0.) {
    	parang.error = cpl_parameter_get_double(par) / 100. * fabs(parang.data);
    } else {
    	parang.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "posang-err");
    if (posang.data != 0.) {
    	posang.error = cpl_parameter_get_double(par) / 100. * fabs(posang.data);
    } else {
    	posang.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "temp-err");
    if (temp.data != 0.) {
    	temp.error = cpl_parameter_get_double(par) / 100. * fabs(temp.data);
    } else {
    	temp.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "rhum-err");
    if (rhum.data != 0.) {
    	rhum.error = cpl_parameter_get_double(par) / 100. * fabs(rhum.data);
    } else {
    	rhum.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "pres-err");
    if (pres.data != 0.) {
    	pres.error = cpl_parameter_get_double(par) / 100. * fabs(pres.data);
    } else {
    	pres.error = cpl_parameter_get_double(par);
    }


	/* Standard methods to calculate airmass: Hardie(62), YoungIrvine(67) or Young(94) */
    hdrl_airmass_approx typeAirmassAprox = 0;
    par = cpl_parameterlist_find_const(parlist, "airm.method");
	const char *str = cpl_parameter_get_string(par);
    if (strcmp(str, "HARDIE_62") == 0) {
    	typeAirmassAprox = HDRL_AIRMASS_APPROX_HARDIE ;
    } else if (strcmp(str, "YOUNG_IRVINE_67") == 0) {
    	typeAirmassAprox = HDRL_AIRMASS_APPROX_YOUNG_IRVINE ;
    } else if (strcmp(str, "YOUNG_94") == 0) {
    	typeAirmassAprox = HDRL_AIRMASS_APPROX_YOUNG ;
    } else {
    	cpl_imagelist_delete(imlistOrig);
    	cpl_imagelist_delete(imlist);
    	cpl_wcs_delete(wcs);
    	cpl_propertylist_delete(header_raw);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                "Approx. airmass mode can only be \"HARDIE_62\", "
                "\"YOUNG_IRVINE_67\" or \"YOUNG_94\" (not %s)", str);
    }

	/* Error parameters for calculate airmass in percent [frac] of the actual value */

    /*** Note: if the error data values are exactly 0., the error parameter is getting as absolute error ***/

    par = cpl_parameterlist_find_const(parlist, "airm.ra-err");
    if (ra.data != 0.) {
    	ra.error = cpl_parameter_get_double(par) / 100. * fabs(ra.data);
    } else {
    	ra.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "airm.dec-err");
    if (dec.data != 0.) {
    	dec.error = cpl_parameter_get_double(par) / 100. * fabs(dec.data);
    } else {
    	dec.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "airm.lst-err");
    if (lst.data != 0.) {
    	lst.error = cpl_parameter_get_double(par) / 100. * fabs(lst.data);
    } else {
    	lst.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "airm.exptime-err");
    if (exptime.data != 0.) {
    	exptime.error = cpl_parameter_get_double(par) / 100. * fabs(exptime.data);
    } else {
    	exptime.error = cpl_parameter_get_double(par);
    }

    par = cpl_parameterlist_find_const(parlist, "airm.geolat-err");
    if (geolat.data != 0.) {
    	geolat.error = cpl_parameter_get_double(par) / 100. * fabs(geolat.data);
    } else {
    	geolat.error = cpl_parameter_get_double(par);
    }


	/********** BUILD INPUT PARAMETERS through LambdaRef ******************/

	/* Initialize the cpl_vector lambdaIni through the size and lambdaRef */
	cpl_vector *lambdaIn  = cpl_vector_new(size);
	double lambda = lambdaRef.data - (planeLambda - 1.) * step;
	for (cpl_size i = 0; i < size; i++) {
		cpl_vector_set(lambdaIn,  i, lambda);
		lambda += step;
	}

    /* Create and initialize the output vectors for the shift and the error */
	cpl_vector *xShift    = cpl_vector_new(size);
	cpl_vector *yShift    = cpl_vector_new(size);
	cpl_vector *xShiftErr = cpl_vector_new(size);
	cpl_vector *yShiftErr = cpl_vector_new(size);
	for (cpl_size i = 0; i < size; i++) {
		cpl_vector_set(xShift,    i, 0.    );
		cpl_vector_set(yShift,    i, 0.    );
		cpl_vector_set(xShiftErr, i, 0.    );
		cpl_vector_set(yShiftErr, i, 0.    );
	}


    cpl_errorstate prestate;


	/************ CALCULE AIRMASS **********************/
	cpl_msg_info(cpl_func,"Calculating airmass ...");
    prestate = cpl_errorstate_get();
	airmass = hdrldemo_utils_airmass( airmass_start, airmass_end,
									  ra, dec, lst, exptime, geolat,
									  typeAirmassAprox);
	if (!cpl_errorstate_is_equal(prestate)) {
		hdrldemo_dar_clean_recipe_parameters(header_raw, wcs, imlistOrig, imlist, NULL, NULL,
				lambdaIn, xShift, yShift, xShiftErr, yShiftErr, NULL);
		return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
				"Airmass compute error");
	}


	/*********** DAR COMPUTATION **************/
	cpl_msg_info(cpl_func,"Calculating differential atmospheric refraction (DAR) ...");
	prestate = cpl_errorstate_get();

	cpl_msg_debug(cpl_func, "lambda  = %g+-[%g]",  lambdaRef.data, lambdaRef.error);
	cpl_msg_debug(cpl_func, "parang  = %g+-[%g]",  parang.data,    parang.error   );
	cpl_msg_debug(cpl_func, "posang  = %g+-[%g]",  posang.data,    posang.error   );
	cpl_msg_debug(cpl_func, "T       = %g+-[%g]",  temp.data, 	   temp.error     );
	cpl_msg_debug(cpl_func, "Hum[%%]  = %g+-[%g]", rhum.data, 	   rhum.error     );
	cpl_msg_debug(cpl_func, "Pre     = %g+-[%g]",  pres.data, 	   pres.error     );
	cpl_msg_debug(cpl_func, "airmass = %g+-[%g]",  airmass.data,   airmass.error  );

	hdrl_parameter *h_par = hdrl_dar_parameter_create(airmass, parang, posang, temp, rhum, pres, wcs);
	hdrl_dar_compute(h_par, lambdaRef, lambdaIn, xShift, yShift, xShiftErr, yShiftErr);
	hdrl_parameter_delete(h_par);

	if (!cpl_errorstate_is_equal(prestate)) {
		hdrldemo_dar_clean_recipe_parameters(header_raw, wcs, imlistOrig, imlist, NULL, NULL,
				lambdaIn, xShift, yShift, xShiftErr, yShiftErr, NULL);
		return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
				"DAR compute error");
	}


	/****** Make the output table for fill data *****/
	cpl_msg_info(cpl_func,"Build table with output results ...");
    prestate = cpl_errorstate_get();

    cpl_table *tb = hdrldemo_dar_create_output_table(size);
    hdrldemo_dar_add_table_wavelenght_and_shift(tb, lambdaIn, xShift, yShift, xShiftErr, yShiftErr);

	if (!cpl_errorstate_is_equal(prestate)) {
		hdrldemo_dar_clean_recipe_parameters(header_raw, wcs, imlistOrig, imlist, NULL, NULL,
				lambdaIn, xShift, yShift, xShiftErr, yShiftErr, tb);
		return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
				"Output table fill error");
	}


	/*********** APPLY SHIFT TO IMAGES ******************/
	cpl_msg_info(cpl_func,"Apply the calculate shift to the image list ...");
	for (cpl_size i = 0; i < size; i++) {

		/* Get original and copy shift images */
		cpl_image *imgOrig = cpl_imagelist_get(imlistOrig, i);
		cpl_image *img     = cpl_imagelist_get(imlist,i);

		/* Apply shift in each plane, round the shift in number of pixels*/
		int xshift = round(cpl_vector_get(xShift, i));
		int yshift = round(cpl_vector_get(yShift, i));
		if (xshift != 0. || yshift != 0.) {
			cpl_image_shift(img, xshift, yshift);
		}

		/* Reject NaN values */
		cpl_image_reject_value(imgOrig, CPL_VALUE_NAN);
		cpl_image_reject_value(img,     CPL_VALUE_NAN);
	}

	/* Create a collapsed image */
	cpl_msg_info(cpl_func,"Collapse original and shifted image list ...");
	cpl_image* img_orig_collapsed    = cpl_imagelist_collapse_create(imlistOrig);
	cpl_image* img_correct_collapsed = cpl_imagelist_collapse_create(imlist);


	/***** CREATE CATALOGUE OF ORIGINAL AND SHIFTED IMAGES *****/

	cpl_msg_info(cpl_func,"Make a catalogue for the original image list ...");
	if (hdrldemo_dar_catalogue(CPL_TRUE, imlistOrig, NULL, maxRatioNaNsInImage, tb)
			!= CPL_ERROR_NONE) {
		hdrldemo_dar_clean_recipe_parameters(header_raw, wcs,
				imlistOrig, imlist, img_orig_collapsed, img_correct_collapsed,
				lambdaIn, xShift, yShift, xShiftErr, yShiftErr, tb);
		return cpl_error_set_message(cpl_func, cpl_error_get_code(),
				"Catalogue imglist original error");
	}

	cpl_msg_info(cpl_func,"Make a catalogue for the shifted image list ...");
	if (hdrldemo_dar_catalogue(CPL_FALSE, imlist, NULL, maxRatioNaNsInImage, tb)
			!= CPL_ERROR_NONE) {
		hdrldemo_dar_clean_recipe_parameters(header_raw, wcs,
				imlistOrig, imlist, img_orig_collapsed, img_correct_collapsed,
				lambdaIn, xShift, yShift, xShiftErr, yShiftErr, tb);
		return cpl_error_set_message(cpl_func, cpl_error_get_code(),
				"Catalogue imglist shifted error");
	}


	/* SAVE DATA: Data shift & catalogue, original images and correct images */
	cpl_msg_info(cpl_func,"Save all data and image to disk and adding QC parameters ...");
	prestate = cpl_errorstate_get();

	cpl_propertylist *props = cpl_propertylist_new();
	cpl_propertylist_update_double(props, "ESO DRS LAMBDAREF", lambdaRef.data);

	hdrldemo_dar_save_disk(header_raw, tb, imlistOrig, imlist, img_orig_collapsed,
			img_correct_collapsed,	parlist, props, frameset, in_frm);

    cpl_propertylist_delete(props);

	if (!cpl_errorstate_is_equal(prestate)) {
		hdrldemo_dar_clean_recipe_parameters(header_raw, wcs,
				imlistOrig, imlist, img_orig_collapsed, img_correct_collapsed,
				lambdaIn, xShift, yShift, xShiftErr, yShiftErr, tb);
		return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO, "Error save data");
	}


	/*********** CLEAN UP **************/
	cpl_msg_info(cpl_func,"Finnish and clean variables ...");
	hdrldemo_dar_clean_recipe_parameters(header_raw, wcs,
			imlistOrig, imlist, img_orig_collapsed, img_correct_collapsed,
			lambdaIn, xShift, yShift, xShiftErr, yShiftErr, tb);

    return (int)cpl_error_get_code();
}

/**
 * @brief  Clean parameters allocated in the recipe
 *
 * @param		header_raw              Primary header of raw input file
 * @param		wcs						World coordinate system
 * @param		imlistOrig				Set with the original images
 * @param		imlist					Set with the shifted images
 * @param		img_orig_collapsed		Original collapsed image
 * @param		img_correct_collapsed	Shifted  collapsed image
 * @param		lambdaIn				Vector with the plane wavelenghts
 * @param		xShift					Vector with the shifts in x-axis
 * @param		yShift					Vector with the shifts in y-axis
 * @param		xShiftErr				Vector with the error shifts in x-axis
 * @param		yShiftErr				Vector with the error shifts in y-axis
 * @param		tb						Table for storage output data
 *
 */
/* ---------------------------------------------------------------------------*/
static void hdrldemo_dar_clean_recipe_parameters(
	cpl_propertylist *header_raw,
	cpl_wcs *wcs, cpl_imagelist *imlistOrig, cpl_imagelist *imlist,
	cpl_image* img_orig_collapsed, cpl_image* img_correct_collapsed,
	cpl_vector *lambdaIn, cpl_vector *xShift, cpl_vector *yShift,
	cpl_vector *xShiftErr, cpl_vector *yShiftErr, cpl_table *tb)
{
	if (header_raw)             cpl_propertylist_delete(header_raw);
	if (wcs)                    cpl_wcs_delete(wcs);

	if (imlistOrig)             cpl_imagelist_delete(imlistOrig);
	if (imlist)                 cpl_imagelist_delete(imlist);

	if (img_orig_collapsed)     cpl_image_delete(img_orig_collapsed);
	if (img_correct_collapsed)  cpl_image_delete(img_correct_collapsed);

	if (lambdaIn)               cpl_vector_delete(lambdaIn);
    if (xShift)                 cpl_vector_delete(xShift);
    if (yShift)                 cpl_vector_delete(yShift);
    if (xShiftErr)              cpl_vector_delete(xShiftErr);
    if (yShiftErr)              cpl_vector_delete(yShiftErr);

    if (tb)                     cpl_table_delete(tb);
}

/**
 * @brief  Verify the content of the property list
 *
 * @param  header			input propertylist
 *
 * @return cpl_error_code,	for error codes see cpl_table_wrap_double
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_verify_propertylist_parameters(
	const cpl_propertylist *header){

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL AIRM START"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL AIRM START");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL AIRM END"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL AIRM END");

	cpl_error_ensure(cpl_propertylist_has(header, "RA"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: RA");

	cpl_error_ensure(cpl_propertylist_has(header, "DEC"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: DEC");

	cpl_error_ensure(cpl_propertylist_has(header, "LST"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: LST");

	cpl_error_ensure(hdrldemo_dar_get_exptime(header) != -1.,
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: EXPTIME");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL GEOLAT"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL GEOLAT");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL GEOLAT"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL GEOLAT");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL PARANG START"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL PARANG START");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL PARANG END"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL PARANG END");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO ADA POSANG"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO ADA POSANG");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL AMBI TEMP"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL AMBI TEMP");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL AMBI RHUM"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL AMBI RHUM");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL AMBI PRES START"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL AMBI PRES START");

	cpl_error_ensure(cpl_propertylist_has(header, "ESO TEL AMBI PRES END"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: ESO TEL AMBI PRES END");

	cpl_error_ensure(cpl_propertylist_has(header, "CRVAL3"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: CRVAL3");

	cpl_error_ensure(cpl_propertylist_has(header, "CDELT3"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: CDELT3");

	cpl_error_ensure(cpl_propertylist_has(header, "CRPIX3"),
			CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_DATA_NOT_FOUND,
			"Doesn't exist in propertylist header: CRPIX3");

	return CPL_ERROR_NONE;
}

/**
 * @brief  Get the exposition time from the header
 *
 * @param  header	input propertylist that you hope that contains EXPTIME
 *
 * @return double,	The time of exposition
 *
 */
/* ---------------------------------------------------------------------------*/
static double hdrldemo_dar_get_exptime(
	const cpl_propertylist *header)
{
    cpl_errorstate prestate = cpl_errorstate_get();

    /* Get value for exptime */
    double value = cpl_propertylist_get_double(header, "EXPTIME");

    if (!cpl_errorstate_is_equal(prestate)) {

    	/* If it doesn't exist reset error, maybe it's storaged in DIT form */
    	cpl_errorstate_set(prestate);

        double DIT = cpl_propertylist_get_double(header, "ESO DET DIT");
         if (cpl_errorstate_is_equal(prestate)) {

        	/* Correct value */
        	value = DIT;

        } else {

        	/* If you don't find it either, return the error value -1 */
        	value = -1.;
        }
    }
    return value;
}

/**
 * @brief  Calcule the paralactic angle
 *
 * @param  p1		Value of parallactic angle when start the observation
 * @param  p2       Value of parallactic angle when end   the observation
 *
 * @return double,	The paralactic angle
 *
 */
/* ---------------------------------------------------------------------------*/
static double hdrldemo_dar_parangle(
	double p1, double p2)
{
    double parang = (p1 + p2) / 2.;
    if (fabs(p1 - p2) < 90.) {
    	/* are not going through the meridian with a 360 deg flip in   *
    	 * PARANG, so it should be safe to just average the two values */
    	return parang;
    }

	/* Flip in PARANG while observing this object, special (pretty          *
	* convoluted) handling needed: Both absolute values should be close to *
	* 180., so compute the absolute distance from 180, the average of the  *
	* distance to 180 with the respective sign, and subtract the average   *
	* of those 180. The final sign is taken from the value whose distance  *
	* to the from respectively signed 180 is larger.                       */

	double d1 = copysign(180. - fabs(p1), p1);
	double d2 = copysign(180. - fabs(p2), p2);

	parang = 180. - fabs((d1 + d2) / 2.);
	if (fabs(d1) > fabs(d2)) {
		parang = copysign(parang, p1);
	} else {
		parang = copysign(parang, p2);
	}

	return parang;
}


/**
 * @brief Create the output table.
 *
 * @param  size        Size of the output table
 *
 * @return cpl_table*, the newly allocated table.
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_table * hdrldemo_dar_create_output_table(cpl_size size)
{
	/* Create table */
    cpl_table *tb = cpl_table_new(size);


    /* Create columns */

    cpl_table_new_column(tb, index_cname, 					CPL_TYPE_INT);

	cpl_table_new_column(tb, orig_bpImg_cname, 				CPL_TYPE_INT);
	cpl_table_new_column(tb, shifted_bpImg_cname, 			CPL_TYPE_INT);

	cpl_table_new_column(tb, lambdaIn_cname, 				CPL_TYPE_DOUBLE);

	cpl_table_new_column(tb, x_cname, 						CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, y_cname, 						CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, xErr_cname, 					CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, yErr_cname, 					CPL_TYPE_DOUBLE);

	cpl_table_new_column(tb, orig_cat_flux_cname, 			CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, orig_cat_xcoord_cname, 		CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, orig_cat_ycoord_cname, 		CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, orig_cat_xcoord_err_cname, 	CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, orig_cat_ycoord_err_cname, 	CPL_TYPE_DOUBLE);

	cpl_table_new_column(tb, shifted_cat_flux_cname, 		CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, shifted_cat_xcoord_cname, 		CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, shifted_cat_ycoord_cname, 		CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, shifted_cat_xcoord_err_cname, 	CPL_TYPE_DOUBLE);
	cpl_table_new_column(tb, shifted_cat_ycoord_err_cname, 	CPL_TYPE_DOUBLE);

	/* Initialize several columns */
	for (cpl_size i = 0; i < size; i++) {

		/* Add Plane number */
		cpl_table_set_int(   tb, index_cname,                   i, i+1				);

		/* Initialize as reject */
		cpl_table_set_int(   tb, orig_bpImg_cname, 				i, (int)CPL_BINARY_1);
		cpl_table_set_double(tb, orig_cat_flux_cname,       	i, 0.    			);
		cpl_table_set_double(tb, orig_cat_xcoord_cname,     	i, 0.    			);
		cpl_table_set_double(tb, orig_cat_ycoord_cname,     	i, 0.    			);
		cpl_table_set_double(tb, orig_cat_xcoord_err_cname, 	i, 0.    			);
		cpl_table_set_double(tb, orig_cat_ycoord_err_cname,		i, 0.    			);

		/* Initialize as reject */
		cpl_table_set_int(   tb, shifted_bpImg_cname, 			i, (int)CPL_BINARY_1);
		cpl_table_set_double(tb, shifted_cat_flux_cname,       	i, 0.    			);
		cpl_table_set_double(tb, shifted_cat_xcoord_cname,     	i, 0.    			);
		cpl_table_set_double(tb, shifted_cat_ycoord_cname,     	i, 0.    			);
		cpl_table_set_double(tb, shifted_cat_xcoord_err_cname, 	i, 0.    			);
		cpl_table_set_double(tb, shifted_cat_ycoord_err_cname, 	i, 0.    			);
	}

	return tb;
}


/**
 * @brief Add wavelenght and shift value and error to the output table.
 *
 * @param tb   		        In/Out parameter: Output table
 * @param lambda   		    Vector with the wavelenght of the planes
 * @param xShift   		    Vector with the shift in x-axis in each plane
 * @param yShift   		    Vector with the shift in y-axis in each plane
 * @param xShiftErr		    Vector with the err in shift in x-axis in each plane
 * @param yShiftErr		    Vector with the err in shift in y-axis in each plane
 *
 * @return cpl_error_code,	for error codes see cpl_table_wrap_double
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_add_table_wavelenght_and_shift(
		cpl_table        *tb,
		const cpl_vector *lambda,
		const cpl_vector *xShift,
		const cpl_vector *yShift,
		const cpl_vector *xShiftErr,
		const cpl_vector *yShiftErr)
{
	/* Check inputs */
	cpl_ensure_code(    tb        != NULL && lambda    != NULL
			         && xShift    != NULL && yShift    != NULL
					 && xShiftErr != NULL && yShiftErr != NULL, CPL_ERROR_NULL_INPUT);

	/* Storage values of wavelenght, shifts and errors in the output table */
	cpl_size size = cpl_table_get_nrow(tb);
	for (cpl_size i = 0; i < size; i++) {

		cpl_table_set_double(tb, lambdaIn_cname, i, cpl_vector_get(lambda,    i));
		cpl_table_set_double(tb, x_cname,     	 i, cpl_vector_get(xShift,    i));
		cpl_table_set_double(tb, y_cname,     	 i, cpl_vector_get(yShift,    i));
		cpl_table_set_double(tb, xErr_cname, 	 i, cpl_vector_get(xShiftErr, i));
		cpl_table_set_double(tb, yErr_cname,	 i, cpl_vector_get(yShiftErr, i));
	}

	return CPL_ERROR_NONE;
}

/**
 * @brief Create an extract values from the catalogue of imlist
 *
 * @param  orig						In: Boolean that decide in which columns write
 * @param  imlist 					In: Set of images
 * @param  wcs      				In: World coordinate system
 * @param  maxRatioOfNaNsInImage    In: Max ratio between NaNs and good prixels
 * @param  tb   					In/Out: Output table for storage results
 *
 * @return cpl_error_code,	for error codes see cpl_table_wrap_double
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_catalogue(
		cpl_boolean   orig,
		cpl_imagelist *imlist,
		cpl_wcs       *wcs,
		double		  maxRatioNaNsInImage,
		cpl_table     *tb)
{

	/* Get the number of images */
	cpl_size size = cpl_imagelist_get_size(imlist);
	if (size <= 0) {
		return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
				"Number of images Zero, NULL of incompatible");
	}

	/* Parameters for catalogue shifted images */
    hdrl_parameter *hCat = hdrl_catalogue_parameter_create(
	    					5,						/* obj_min_pixels    */
							2.5,					/* obj_threshold     */
							CPL_TRUE,				/* obj_deblending    */
							3.0,					/* obj_core_radius   */
							CPL_TRUE,				/* bkg_estimate      */
							16,						/* bkg_mesh_size     */
							2.0,					/* bkg_smooth_fwhm   */
							2.8,					/* det_eff_gain      */
							65200,					/* det_saturation    */
							HDRL_CATALOGUE_ALL);	/* Catalogue options */

    /* Inicialize counters in order to found problems fonts */
    cpl_size nProblems    = 0;
    cpl_size nManySources = 0;

	/* loop all images (planes) */
	for (cpl_size i = 0; i < size; i++) {

		/* Get image */
		cpl_image *img = cpl_imagelist_get(imlist,i);
		if (img == NULL) {
			hdrl_parameter_delete(hCat);
			return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Image NULL");
		}

		cpl_errorstate prestate = cpl_errorstate_get();

		/* Get size of the image */
		cpl_size xPixels = cpl_image_get_size_x(img);
		cpl_size yPixels = cpl_image_get_size_y(img);
		cpl_size nPixels = xPixels * yPixels;

		/* Create Confidence matrix, and accept all pixels */
		cpl_image *cnf = cpl_image_new(xPixels, yPixels, CPL_TYPE_INT);
		cpl_image_add_scalar(cnf, 100.);

		/* Get Bad Pixel mask from the image and modify to the confidence matrix */
		const cpl_mask *bpm = cpl_image_get_bpm_const(img);
		if (bpm != NULL) {
			cpl_image_reject_from_mask(cnf, bpm);
			cpl_image_fill_rejected(cnf, 0);
			cpl_image_accept_all(cnf);
		}

		/* Check the NaNs in the image, and apply if it's allow */
		cpl_size nNaNs = cpl_image_count_rejected(img);
		double ratioNaNsInImage = (double)nNaNs/(double)nPixels;
		if (ratioNaNsInImage < maxRatioNaNsInImage) {

			/* Compute cataloge */
			hdrl_catalogue_result *cat_results = hdrl_catalogue_compute(img, cnf, wcs, hCat);

			/* Get cpl_table in the catalogue for treatement */
			cpl_table *catalogue = cat_results->catalogue;
			cpl_size nRows = cpl_table_get_nrow(catalogue);
			if (nRows > 0) {

				/* Exist at least 1 source in the image */

				if (nRows > 1) {

					/* Too many sources in the image, increase counter */
					nManySources++;
					cpl_msg_debug(cpl_func, "IMAGE=%"CPL_SIZE_FORMAT" (Pixels=%"
						CPL_SIZE_FORMAT", nNans=%"CPL_SIZE_FORMAT", ratioNaNs=%f)"
						" => Number of sources in the picture = %"CPL_SIZE_FORMAT,
						i, nPixels, nNaNs, ratioNaNsInImage, nRows);
				}

				/* Get parameters in the catalogue table */
				int    nl;
				double flux  = cpl_table_get_double(catalogue, "Isophotal_flux",   (cpl_size)0, &nl);
				double x     = cpl_table_get_double(catalogue, "X_coordinate",     (cpl_size)0, &nl);
				double y     = cpl_table_get_double(catalogue, "Y_coordinate",     (cpl_size)0, &nl);
				double x_err = cpl_table_get_double(catalogue, "X_coordinate_err", (cpl_size)0, &nl);
				double y_err = cpl_table_get_double(catalogue, "Y_coordinate_err", (cpl_size)0, &nl);


				/* Depending of the orig CPL_Boolean, put the parameters in the
				 * original columns or in the shifted columns
				 */
				if (orig == CPL_TRUE) {
					cpl_table_set_int(   tb, orig_bpImg_cname, 				i, (int)CPL_BINARY_0);
					cpl_table_set_double(tb, orig_cat_flux_cname,       	i, flux    			);
					cpl_table_set_double(tb, orig_cat_xcoord_cname,     	i, x    			);
					cpl_table_set_double(tb, orig_cat_ycoord_cname,     	i, y    			);
					cpl_table_set_double(tb, orig_cat_xcoord_err_cname, 	i, x_err   			);
					cpl_table_set_double(tb, orig_cat_ycoord_err_cname,		i, y_err   			);
				} else {
					cpl_table_set_int(   tb, shifted_bpImg_cname, 			i, (int)CPL_BINARY_0);
					cpl_table_set_double(tb, shifted_cat_flux_cname,       	i, flux    			);
					cpl_table_set_double(tb, shifted_cat_xcoord_cname,     	i, x    			);
					cpl_table_set_double(tb, shifted_cat_ycoord_cname,     	i, y   				);
					cpl_table_set_double(tb, shifted_cat_xcoord_err_cname, 	i, x_err   			);
					cpl_table_set_double(tb, shifted_cat_ycoord_err_cname, 	i, y_err   			);
				}

			} else {

				/* Doesn't exist anything in the image, add problem */
				nProblems++;
				cpl_msg_debug(cpl_func,
					"IMAGE=%"CPL_SIZE_FORMAT" (Pixels=%"CPL_SIZE_FORMAT", nNans=%"
					CPL_SIZE_FORMAT", ratioNaNs=%f) => It doesn't have Objects! ",
					i, nPixels, nNaNs, ratioNaNsInImage);
				cpl_errorstate_set(prestate);
			}

			/* Delete the catalogue */
			hdrl_catalogue_result_delete(cat_results);

		} else {

			nProblems++;
			cpl_msg_debug(cpl_func,	"IMAGE=%"CPL_SIZE_FORMAT" (Pixels=%"CPL_SIZE_FORMAT
					", nNans=%"CPL_SIZE_FORMAT", ratioNaNs=%f) => "
					" Too many bad pixels: It's not possible to make a catalogue!",
					i, nPixels, nNaNs, ratioNaNsInImage);
		}

		/* Delete the confidence matrix */
		cpl_image_delete(cnf);

		/* Check catalogue of the image */
		if (cpl_error_get_code() != CPL_ERROR_NONE) {
			hdrl_parameter_delete(hCat);
			return cpl_error_set_message(cpl_func, cpl_error_get_code(), "Calcule catalogue from img");
		}
	}

	/* Check problems making the catalogue of the set of images */
	cpl_msg_debug(cpl_func, "Catalogue => nImages = %" CPL_SIZE_FORMAT
		", nManySources = %" CPL_SIZE_FORMAT "[ratio = %f %%], nProblems = %"
		CPL_SIZE_FORMAT " [ratio = %f %%]",
		size, nManySources, 100. * (double)nManySources / (double)size,
			  nProblems,    100. * (double)nProblems    / (double)size);

	hdrl_parameter_delete(hCat);

    return CPL_ERROR_NONE;
}

/**
 * @brief save to disk all of the output data
 *
 * @param  header_raw           Primary header of raw input file cube image
 * @param  tb 					output table with all of data correctness
 * @param  imlistOrig			set images original
 * @param  imlistCorrect	    set images correctness
 * @param  imOrigCollapsed      Collapsed the original set of images
 * @param  imCorrectCollapsed   Collapsed the correctness set of images
 * @param  parlist				input with entry parameters
 * @param  app_list				list of calculate parameters (QC)
 * @param  frameset				input set of frames
 * @param  raw_frame			RAW frame info
 *
 */
/* ---------------------------------------------------------------------------*/
static void hdrldemo_dar_save_disk(
	cpl_propertylist        *header_raw,
	cpl_table 				*tb,
	cpl_imagelist           *imlistOrig,
	cpl_imagelist           *imlistCorrect,
	cpl_image               *imOrigCollapsed,
	cpl_image               *imCorrectCollapsed,
	const cpl_parameterlist *parlist,
	const cpl_propertylist 	*app_list,
	cpl_frameset 			*frameset,
	cpl_frame 				*raw_frame)
{

    /* Propertylist for image FITS file */
    cpl_propertylist *plist_image = cpl_propertylist_new();

	cpl_propertylist_update_double(plist_image,     "CRPIX1", cpl_propertylist_get_double(header_raw, "CRPIX1"));
	cpl_propertylist_update_double(plist_image,     "CRPIX2", cpl_propertylist_get_double(header_raw, "CRPIX2"));

	cpl_propertylist_update_double(plist_image,     "CRVAL1", cpl_propertylist_get_double(header_raw, "CRVAL1"));
	cpl_propertylist_update_double(plist_image,     "CRVAL2", cpl_propertylist_get_double(header_raw, "CRVAL2"));

	cpl_propertylist_update_string(plist_image,     "CTYPE1", cpl_propertylist_get_string(header_raw, "CTYPE1"));
	cpl_propertylist_update_string(plist_image,     "CTYPE2", cpl_propertylist_get_string(header_raw, "CTYPE2"));

	cpl_propertylist_update_string(plist_image,     "CUNIT1", cpl_propertylist_get_string(header_raw, "CUNIT1"));
	cpl_propertylist_update_string(plist_image,     "CUNIT2", cpl_propertylist_get_string(header_raw, "CUNIT2"));

	cpl_propertylist_update_double(plist_image,     "CDELT1", cpl_propertylist_get_double(header_raw, "CDELT1"));
	cpl_propertylist_update_double(plist_image,     "CDELT2", cpl_propertylist_get_double(header_raw, "CDELT2"));

	cpl_propertylist_update_double(plist_image,     "CD1_1",  cpl_propertylist_get_double(header_raw, "CD1_1" ));
	cpl_propertylist_update_double(plist_image,     "CD1_2",  cpl_propertylist_get_double(header_raw, "CD1_2" ));

	cpl_propertylist_update_double(plist_image,     "CD2_1",  cpl_propertylist_get_double(header_raw, "CD2_1" ));
	cpl_propertylist_update_double(plist_image,     "CD2_2",  cpl_propertylist_get_double(header_raw, "CD2_2" ));


    /* Propertylist from imagelist FITS file */
    cpl_propertylist *plist_imagelist = cpl_propertylist_duplicate(plist_image);

	cpl_propertylist_update_double(plist_imagelist, "CRPIX3", cpl_propertylist_get_double(header_raw, "CRPIX3"));
	cpl_propertylist_update_double(plist_imagelist, "CRVAL3", cpl_propertylist_get_double(header_raw, "CRVAL3"));
	cpl_propertylist_update_string(plist_imagelist, "CTYPE3", cpl_propertylist_get_string(header_raw, "CTYPE3"));
	cpl_propertylist_update_string(plist_imagelist, "CUNIT3", cpl_propertylist_get_string(header_raw, "CUNIT3"));
	cpl_propertylist_update_double(plist_imagelist, "CDELT3", cpl_propertylist_get_double(header_raw, "CDELT3"));

	cpl_propertylist_update_double(plist_imagelist, "CD1_3",  cpl_propertylist_get_double(header_raw, "CD1_3" ));
	cpl_propertylist_update_double(plist_imagelist, "CD2_3",  cpl_propertylist_get_double(header_raw, "CD2_3" ));
	cpl_propertylist_update_double(plist_imagelist, "CD3_1",  cpl_propertylist_get_double(header_raw, "CD3_1" ));
	cpl_propertylist_update_double(plist_imagelist, "CD3_2",  cpl_propertylist_get_double(header_raw, "CD3_2" ));
	cpl_propertylist_update_double(plist_imagelist, "CD3_3",  cpl_propertylist_get_double(header_raw, "CD3_3" ));


    /* Save original imagelist set */
	cpl_propertylist_update_string(plist_imagelist, CPL_DFS_PRO_CATG, "HDRLDEMO_DAR_ORIG");
    cpl_dfs_save_imagelist(	frameset, NULL,	parlist, frameset, raw_frame,
							imlistOrig,	CPL_TYPE_DOUBLE, RECIPE_NAME,
							plist_imagelist, NULL, PACKAGE "/" PACKAGE_VERSION,
							"hdrldemo_dar_orig.fits");

    /* Save shifted imagelist set */
	cpl_propertylist_update_string(plist_imagelist, CPL_DFS_PRO_CATG, "HDRLDEMO_DAR_CORRECTED");
    cpl_dfs_save_imagelist(	frameset, NULL,	parlist, frameset, raw_frame,
    						imlistCorrect,	CPL_TYPE_DOUBLE, RECIPE_NAME,
							plist_imagelist, NULL, PACKAGE "/" PACKAGE_VERSION,
							"hdrldemo_dar_corrected.fits");

    /* Save collapsed original image */
	cpl_propertylist_update_string(plist_image, CPL_DFS_PRO_CATG, "HDRLDEMO_DAR_ORIG_COLLAPSED");
    cpl_dfs_save_image(		frameset, NULL, parlist, frameset, raw_frame,
							imOrigCollapsed, CPL_TYPE_DOUBLE, RECIPE_NAME,
							plist_image, NULL, PACKAGE "/" PACKAGE_VERSION,
							"hdrldemo_dar_orig_collapsed.fits");

    /* Save collapsed shifted image */
	cpl_propertylist_update_string(plist_image, CPL_DFS_PRO_CATG, "HDRLDEMO_DAR_CORRECTED_COLLAPSED");
    cpl_dfs_save_image(		frameset, NULL, parlist, frameset, raw_frame,
							imCorrectCollapsed, CPL_TYPE_DOUBLE, RECIPE_NAME,
							plist_image, NULL, PACKAGE "/" PACKAGE_VERSION,
							"hdrldemo_dar_corrected_collapsed.fits");

    /* Add a QC parameter  */
    cpl_propertylist *qclist = cpl_propertylist_new();

	/* Extract rows and calculate QC parameters about the original output table.
	 *  Select columns !reject => !(CPL_BINARY_1) */
	cpl_table_select_all(tb);
    cpl_table_and_selected_int(tb, orig_bpImg_cname,    CPL_NOT_EQUAL_TO, (int)CPL_BINARY_1);
    cpl_table* orig_tb_tmp = cpl_table_extract_selected(tb);

	/* Extract rows and calculate QC parameters about the shifted output table.
	 *  Select columns !reject => !(CPL_BINARY_1) */
	cpl_table_select_all(tb);
    cpl_table_and_selected_int(tb, shifted_bpImg_cname, CPL_NOT_EQUAL_TO, (int)CPL_BINARY_1);
    cpl_table* shifted_tb_tmp = cpl_table_extract_selected(tb);

	/* Extract rows and calculate QC parameters about the original and shifted output table.
	 *  Select columns !reject in both tables => !(CPL_BINARY_1) */
	cpl_table_select_all(tb);
    cpl_table_and_selected_int(tb, orig_bpImg_cname,    CPL_NOT_EQUAL_TO, (int)CPL_BINARY_1);
    cpl_table_and_selected_int(tb, shifted_bpImg_cname, CPL_NOT_EQUAL_TO, (int)CPL_BINARY_1);
    cpl_table* completed_tb_tmp = cpl_table_extract_selected(tb);

    /* Info about the number of columns treated */
    cpl_size nRowsCompletedValid = cpl_table_get_nrow(completed_tb_tmp);
    cpl_size nRowsOrigValid      = cpl_table_get_nrow(orig_tb_tmp);
    cpl_size nRowsShiftedValid   = cpl_table_get_nrow(shifted_tb_tmp);
    if (nRowsCompletedValid != nRowsOrigValid || nRowsCompletedValid != nRowsShiftedValid) {
		cpl_msg_debug(cpl_func, "nRowsOutputTable = %" CPL_SIZE_FORMAT
			", nRowsCompleted = %" CPL_SIZE_FORMAT ", NRowsOrig = %"
			CPL_SIZE_FORMAT ", nRowsShifted = %" CPL_SIZE_FORMAT "",
			cpl_table_get_nrow(tb), nRowsCompletedValid, nRowsOrigValid, nRowsShiftedValid);
    }

    /* Insert in the output table QC statistical parameters */
    cpl_errorstate prestate = cpl_errorstate_get();
 	cpl_propertylist_update_string(qclist, CPL_DFS_PRO_CATG, "HDRLDEMO_DAR_RESULT");

 	hdrldemo_dar_add_QC_parameters(qclist, "SHIFTX",     	 completed_tb_tmp, x_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "SHIFTX ERR", 	 completed_tb_tmp, xErr_cname);

 	hdrldemo_dar_add_QC_parameters(qclist, "SHIFTY",     	 completed_tb_tmp, y_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "SHIFTY ERR", 	 completed_tb_tmp, yErr_cname);

	hdrldemo_dar_add_QC_RMS_XY(    qclist, "SHIFT",          completed_tb_tmp, x_cname,    y_cname);
	hdrldemo_dar_add_QC_RMS_XY(    qclist, "SHIFT ERR", 	 completed_tb_tmp, xErr_cname, yErr_cname);

 	hdrldemo_dar_add_QC_parameters(qclist, "CAT ORIG FLUX",  orig_tb_tmp,      orig_cat_flux_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT ORIGX",      orig_tb_tmp,      orig_cat_xcoord_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT ORIGX ERR",  orig_tb_tmp,      orig_cat_xcoord_err_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT ORIGY",      orig_tb_tmp,      orig_cat_ycoord_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT ORIGY ERR",  orig_tb_tmp,      orig_cat_ycoord_err_cname);

 	hdrldemo_dar_add_QC_RMS_XY(    qclist, "CAT ORIG",       orig_tb_tmp,      orig_cat_xcoord_cname,     orig_cat_ycoord_cname);
 	hdrldemo_dar_add_QC_RMS_XY(    qclist, "CAT ORIG ERR",   orig_tb_tmp,      orig_cat_xcoord_err_cname, orig_cat_ycoord_err_cname);

 	hdrldemo_dar_add_QC_parameters(qclist, "CAT SHIFT FLUX", shifted_tb_tmp,   shifted_cat_flux_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT SHIFTX",     shifted_tb_tmp,   shifted_cat_xcoord_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT SHIFTX ERR", shifted_tb_tmp,   shifted_cat_xcoord_err_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT SHIFTY",     shifted_tb_tmp,   shifted_cat_ycoord_cname);
 	hdrldemo_dar_add_QC_parameters(qclist, "CAT SHIFTY ERR", shifted_tb_tmp,   shifted_cat_ycoord_err_cname);

 	hdrldemo_dar_add_QC_RMS_XY(    qclist, "CAT SHIFT",      shifted_tb_tmp,   shifted_cat_xcoord_cname,     shifted_cat_ycoord_cname);
 	hdrldemo_dar_add_QC_RMS_XY(    qclist, "CAT SHIFT ERR",  shifted_tb_tmp,   shifted_cat_xcoord_err_cname, shifted_cat_ycoord_err_cname);

 	/* Reinitialize all the rows in the table as selected */
	cpl_table_select_all(tb);

	/* Delete the temporary tables */
	cpl_table_delete(completed_tb_tmp);
	cpl_table_delete(orig_tb_tmp     );
	cpl_table_delete(shifted_tb_tmp  );

	/* Check if exist errors */
  	if (!cpl_errorstate_is_equal(prestate)) {
		cpl_msg_warning(cpl_func,"Error make QC parameters: %s", cpl_error_get_message());
		cpl_errorstate_set(prestate);
	}

  	/* Save to disk the output table */
	cpl_dfs_save_table(		frameset, NULL, parlist, frameset, raw_frame,
							tb, app_list, RECIPE_NAME,
							qclist, NULL, PACKAGE "/" PACKAGE_VERSION,
							"hdrldemo_dar_result.fits");

	/* Delete propertylist of QC parameters */
	   cpl_propertylist_delete(plist_image);
	   cpl_propertylist_delete(plist_imagelist);
	   cpl_propertylist_delete(qclist);
}

/**
 * @brief Add to a propertylist QC statistical data: (MAX, MIN, MEAN, MEDIAN, STDEV)
 *
 * @param  qclist	Property list of QC parameters
 * @param  parm		Output parameter
 * @param  tb 		Output table with all of data correctness
 * @param  cname	Name of the column in the table
 *
 * @return cpl_error_code,	for error codes see cpl_table_wrap_double
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_add_QC_parameters(
	cpl_propertylist *qclist,
	const char       *param,
	const cpl_table  *tb,
	const char       *cname)
{
	/* Check inputs */
    cpl_ensure_code(	qclist != NULL && param != NULL && tb != NULL
    		         && cname  != NULL, CPL_ERROR_NULL_INPUT);

    cpl_error_code tmp, e = CPL_ERROR_NONE;

    /* Output value for the different statistical parameters in column param */
    double value;

	value = cpl_table_get_column_min(tb, cname);
	tmp = hdrldemo_dar_add_property(qclist, param, "MIN", value);
    if (tmp != CPL_ERROR_NONE) e = tmp;

    value = cpl_table_get_column_max(tb, cname);
    tmp = hdrldemo_dar_add_property(qclist, param, "MAX", value);
    if (tmp != CPL_ERROR_NONE) e = tmp;

	value = cpl_table_get_column_mean(tb, cname);
	tmp = hdrldemo_dar_add_property(qclist, param, "MEAN", value);
    if (tmp != CPL_ERROR_NONE) e = tmp;

    value = cpl_table_get_column_median(tb, cname);
    tmp = hdrldemo_dar_add_property(qclist, param, "MEDIAN", value);
    if (tmp != CPL_ERROR_NONE) e = tmp;

    value = cpl_table_get_column_stdev(tb, cname);
    tmp = hdrldemo_dar_add_property(qclist, param, "STDEV", value);
    if (tmp != CPL_ERROR_NONE) e = tmp;

    hdrldemo_dar_add_QC_parameter_mad(qclist, param, tb, cname);
    if (tmp != CPL_ERROR_NONE) e = tmp;

	return e;
}

static cpl_error_code hdrldemo_dar_add_QC_parameter_mad(
	cpl_propertylist *qclist,
	const char       *param,
	const cpl_table  *tb,
	const char       *cname){

	/* Check inputs */
    cpl_ensure_code(	qclist != NULL && param != NULL && tb != NULL
    		         && cname  != NULL, CPL_ERROR_NULL_INPUT);

    cpl_error_code tmp, e = CPL_ERROR_NONE;


    /* Get values and put in a image*/
	int   nl;
    cpl_size size = cpl_table_get_nrow(tb);
    cpl_image *img = cpl_image_new(size, 1, CPL_TYPE_DOUBLE);

	for (cpl_size i = 0; i < size; i++) {
		double value = cpl_table_get_double(tb, cname, i, &nl);
		cpl_image_set(img, i+1, 1, value);
	}

	double mad = 0.;
	cpl_image_get_mad(img, &mad);
	double std_mad = CPL_MATH_STD_MAD * mad;

	cpl_image_delete(img);

    tmp = hdrldemo_dar_add_property(qclist, param, "MAD", std_mad);
    if (tmp != CPL_ERROR_NONE) e = tmp;

	return e;
}

/**
 * @brief Add to a propertylist QC RMS for two columns (X,Y)
 *
 * @param  qclist	Property list of QC parameters
 * @param  parm		Output parameter
 * @param  tb 		Output table with all of data correctness
 * @param  cnameX	Name of the column in the table that is the X-axis
 * @param  cnameY	Name of the column in the table that is the Y-axis
 *
 * @return cpl_error_code,	for error codes see cpl_table_wrap_double
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_add_QC_RMS_XY(
	cpl_propertylist *qclist,
	const char       *param,
	const cpl_table  *tb,
	const char       *cnameX,
	const char       *cnameY)
{
	/* Check inputs */
    cpl_ensure_code(	qclist != NULL && param  != NULL && tb != NULL
					 && cnameX != NULL && cnameY != NULL, CPL_ERROR_NULL_INPUT);

	cpl_errorstate prestate = cpl_errorstate_get();

	/* Calculate RMS of the distance between the x and y direction */
	double stdevX = cpl_table_get_column_stdev(tb, cnameX);
	double stdevY = cpl_table_get_column_stdev(tb, cnameY);

 	if (!cpl_errorstate_is_equal(prestate)) {
 		cpl_msg_warning(cpl_func,"QC parameter RMS, error calculate stdev: %s",
			cpl_error_get_message());
		return CPL_ERROR_NULL_INPUT;
	}

	/* RMS = SQRT( (X1**2 + X2**2 + ... + XN**2) / (N - 1) ) */
	double value = sqrt(stdevX * stdevX + stdevY * stdevY);

	/* Add property to the list */
	return hdrldemo_dar_add_property(qclist, param, "RMS", value);
}

/**
 * @brief Add to a propertylist QC RMS for two columns (X,Y)
 *
 * @param  qclist	Property list of QC parameters
 * @param  parm		Output parameter
 * @param  qcName	Kind of parameter
 * @param  value    Value of parameter
 *
 * @return cpl_error_code,	for error codes see cpl_table_wrap_double
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_dar_add_property(
	cpl_propertylist *qclist,
	const char       *param,
	const char       *qcName,
	const double     value){

    cpl_error_code e = CPL_ERROR_NONE;

    /* Make the string to insert in the list */
	char *str = cpl_sprintf("ESO QC %s %s", param, qcName);

	/* Check the correctness of the value -> !NaNs and !Inf */
	if (isfinite(value)) {
		cpl_propertylist_update_double(qclist, str, value);
	} else {
		cpl_propertylist_update_double(qclist, str, -99.);
		cpl_msg_warning(cpl_func,"QC parameter [%s] = %g -> Not added to table", str, value);
		e = CPL_ERROR_ILLEGAL_OUTPUT;
	}

	/* Free memory of the string */
	cpl_free((void*)str);

	return e;
}


/**@}*/
