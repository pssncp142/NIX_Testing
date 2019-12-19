/*
 * This file is part of the HDRLDEMO Toolkit.
 * Copyright (C) 2016 European Southern Observatory
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
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/*-----------------------------------------------------------------------------
                               Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>

#include "hdrldemo_utils.h"
#include "hdrl_parameter.h"

/*-----------------------------------------------------------------------------
                               Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_hdrldemo_correct_overscan(void);
static void test_hdrldemo_get_naxis(void);
static void test_hdrldemo_save_image(void);
static void test_hdrldemo_airmass(void);

/*-----------------------------------------------------------------------------
                               Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    test_hdrldemo_correct_overscan();
    test_hdrldemo_get_naxis();
    test_hdrldemo_save_image();
    test_hdrldemo_airmass();

    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                               Private function codes
 -----------------------------------------------------------------------------*/

static void test_hdrldemo_correct_overscan(void){

	/* parameter parsing smoketest */
	hdrl_parameter *rect_region_def = hdrl_rect_region_parameter_create(1, 1, 20, 20);
	hdrl_parameter *sigclip_def     = hdrl_collapse_sigclip_parameter_create(3., 3., 5);
	hdrl_parameter *minmax_def      = hdrl_collapse_minmax_parameter_create(2., 3.);

	cpl_parameterlist *pos = hdrl_overscan_parameter_create_parlist(
								"RECIPE", "oscan", "alongX", 10, 10., rect_region_def,
								"MINMAX", sigclip_def, minmax_def);
	cpl_test_error(CPL_ERROR_NONE);
	cpl_test_eq(cpl_parameterlist_get_size(pos), 13);
	hdrl_parameter_delete(rect_region_def);
	hdrl_parameter_delete(sigclip_def);
	hdrl_parameter_delete(minmax_def);


    cpl_frameset   *frameset  = cpl_frameset_new();
	int            ext_num    = 0;
	hdrl_parameter *os_params = hdrl_overscan_parameter_parse_parlist(pos, "RECIPE.oscan");
    cpl_parameterlist_delete(pos);

    /* Create image */
    cpl_size nx = 64;
    cpl_size ny = 64;
    hdrl_image *masterbias = hdrl_image_new(nx, ny);

    /* Create output imagelist */
    hdrl_imagelist *corrected = hdrl_imagelist_new();

    /* Create buffer */
    hdrl_buffer *buf = hdrl_buffer_new();

    hdrldemo_os_correct(frameset, ext_num, os_params, masterbias, corrected, buf);
    cpl_test_error(CPL_ERROR_NONE);

    /* Clean up */
    cpl_frameset_delete(frameset);
    hdrl_parameter_destroy(os_params);
    hdrl_image_delete(masterbias);
    hdrl_imagelist_delete(corrected);
    hdrl_buffer_delete(buf);
}

static void test_hdrldemo_get_naxis(void)
{
    cpl_propertylist *props = cpl_propertylist_new();

    int result;

    /* Test error handling if passing NULL pointers. */
    result = hdrldemo_get_naxis1(NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(result, 0);
    result = hdrldemo_get_naxis2(NULL);
    cpl_test_eq(result, 0);
    cpl_test_error(CPL_ERROR_NULL_INPUT);


    /* Test error handling when the keywords are missing from the props list. */
    result = hdrldemo_get_naxis1(props);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_eq(result, 0);
    result = hdrldemo_get_naxis2(props);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_eq(result, 0);

    /* Test successful extraction of the NAXIS* keywords. */
    cpl_propertylist_append_int(props, "NAXIS1", 2);
    cpl_propertylist_append_int(props, "NAXIS2", 3);
    result = hdrldemo_get_naxis1(props);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(result, 2);
    result = hdrldemo_get_naxis2(props);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(result, 3);

    /* Test successful extraction of the ZNAXIS* keywords. */
    cpl_propertylist_delete(props);
    props = cpl_propertylist_new();
    cpl_propertylist_append_int(props, "ZNAXIS1", 4);
    cpl_propertylist_append_int(props, "ZNAXIS2", 5);
    result = hdrldemo_get_naxis1(props);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(result, 4);
    result = hdrldemo_get_naxis2(props);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(result, 5);

    cpl_propertylist_delete(props);
}

static void test_hdrldemo_save_image(void)
{
    cpl_propertylist  *app_list = cpl_propertylist_new();
    cpl_propertylist  *qclist   = cpl_propertylist_new();
    cpl_image         *image    = cpl_image_new(10, 10, CPL_TYPE_FLOAT);
    cpl_parameterlist *parlist  = cpl_parameterlist_new();
    cpl_frameset      *frameset = cpl_frameset_new();

    /* Make sure all files we might produce in these tests do not yet exist. */
    const char * filenames[] = {
		"img_test_dummy_raw.fits",
		"img_test1.fits",
		"img_test2.fits",
		"img_test3.fits",
		"img_test4.fits",
		"img_test5.fits",
		NULL
    };

    /* Prepare QC parameter list */
	cpl_propertylist_update_int(   qclist, "ESO QC TEST1", 1      );
	cpl_propertylist_update_double(qclist, "ESO QC TEST2", 2.0    );
	cpl_propertylist_update_string(qclist, "ESO QC TEST3", "value");

	/* Save image to disk */
    cpl_image_save(image, filenames[0], CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);

    /* Prepare a dummy input raw frame. */
	cpl_frame *frame = cpl_frame_new();
	cpl_frame_set_filename( frame, filenames[0]);
	cpl_frame_set_tag(      frame, "RAW");
    cpl_frame_set_type(     frame, CPL_FRAME_TYPE_IMAGE);
    cpl_frame_set_group(    frame, CPL_FRAME_GROUP_RAW);
    cpl_frame_set_level(    frame, CPL_FRAME_LEVEL_TEMPORARY);

    /* Insert frame in frameset */
    cpl_frameset_insert(frameset, frame);


    /****** TESTS ****/
    cpl_error_code result;

    /* Test error for cases where NULL pointers are given */
    result = hdrldemo_save_image(NULL, NULL, NULL, NULL,
    		NULL, CPL_TYPE_INVALID, NULL, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, NULL, NULL, NULL,
    		NULL, CPL_TYPE_INVALID, NULL, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, qclist, NULL, NULL,
    		NULL, CPL_TYPE_INVALID, NULL, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, qclist, "PROD", NULL,
    		NULL, CPL_TYPE_INVALID, NULL, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, qclist, "PROD", "test_recipe",
    		NULL, CPL_TYPE_INVALID, NULL, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, qclist, "PROD", "test_recipe",
    		filenames[1], CPL_TYPE_INVALID, NULL, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, qclist, "PROD", "test_recipe",
    		filenames[2], CPL_TYPE_FLOAT, NULL, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, qclist, "PROD", "test_recipe",
    		filenames[3], CPL_TYPE_FLOAT, image, NULL, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    result = hdrldemo_save_image(app_list, qclist, "PROD", "test_recipe",
    		filenames[4], CPL_TYPE_FLOAT, image, parlist, NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);
    cpl_test(!hdrldemo_check_file_exist(filenames[4]));
    cpl_test(!hdrldemo_check_image_exists(filenames[4], 0));

    result = hdrldemo_save_image(app_list, qclist, "PROD", "test_recipe",
    		filenames[5], CPL_TYPE_FLOAT, image, parlist, frameset);
    cpl_test_eq_error(result, CPL_ERROR_NONE);
    cpl_test(hdrldemo_check_file_exist(filenames[5]));
    cpl_test(hdrldemo_check_image_exists(filenames[5], 0));


    /********* Check parameters QC list **********/

    /* Not exist File with qc parameters */
    cpl_test(!hdrldemo_check_param( filenames[4], 0, CPL_DFS_PRO_CATG,
    		CPL_TRUE, CPL_TYPE_STRING, 0, 0.0, "PROD"));

    /* Check product category value */
    cpl_test(hdrldemo_check_param( filenames[5], 0, CPL_DFS_PRO_CATG,
    		CPL_TRUE, CPL_TYPE_STRING, 0, 0.0, "PROD"));

    /* Not exits QC parameter */
    cpl_test(!hdrldemo_check_param(filenames[5], 0, "ESO QC TEST",
    		CPL_FALSE, CPL_TYPE_INT, 0, 0.0, ""));

    /* Exist int, check value */
    cpl_test(!hdrldemo_check_param(filenames[5], 0, "ESO QC TEST1",
    		CPL_TRUE, CPL_TYPE_INT, 0, 0.0, ""));
    cpl_test(hdrldemo_check_param( filenames[5], 0, "ESO QC TEST1",
    		CPL_TRUE, CPL_TYPE_INT, 1, 0.0, ""));

    /* Exist double, check value */
    cpl_test(!hdrldemo_check_param(filenames[5], 0, "ESO QC TEST2",
    		CPL_TRUE, CPL_TYPE_DOUBLE, 0, 0.0, ""));
    cpl_test(hdrldemo_check_param( filenames[5], 0, "ESO QC TEST2",
    		CPL_TRUE, CPL_TYPE_DOUBLE, 0, 2.0, ""));

    /* Exist string, check value */
    cpl_test(!hdrldemo_check_param(filenames[5], 0, "ESO QC TEST3",
    		CPL_TRUE, CPL_TYPE_STRING, 0, 0.0, ""));
    cpl_test(hdrldemo_check_param( filenames[5], 0, "ESO QC TEST3",
    		CPL_TRUE, CPL_TYPE_STRING, 0, 0.0, "value"));


    /* Cleanup up */
    cpl_propertylist_delete(app_list);
    cpl_propertylist_delete(qclist);
    cpl_frameset_delete(frameset);
    cpl_parameterlist_delete(parlist);
    cpl_image_delete(image);

    /* Delete files */
    for (const char **filename = filenames; *filename != NULL; ++filename) {
        if (hdrldemo_check_file_exist(*filename)) {
        	remove(*filename);
        }
    }
}

static void test_hdrldemo_airmass(void)
{
    hdrl_value airmass_start = {1.542, 1e-2};
    hdrl_value airmass_end   = {1.545, 1e-2};

    hdrl_value ra            = {272.985722, 1e-2};
    hdrl_value dec           = {-32.64367,  1e-2};
    hdrl_value lst           = {79111.582,  1e-2};
    hdrl_value exptime       = {1.,         1e-2};
    hdrl_value geolat        = {-24.627,    1e-2};

	/* Standard methods to calculate airmass: Hardie(62), YoungIrvine(67) or Young(94) */
    hdrl_airmass_approx type1 = HDRL_AIRMASS_APPROX_HARDIE;
    hdrl_airmass_approx type2 = HDRL_AIRMASS_APPROX_YOUNG_IRVINE;
    hdrl_airmass_approx type3 = HDRL_AIRMASS_APPROX_YOUNG;


    /******** Test errors ************/
    hdrl_value err = {-999., -999.};
    hdrl_airmass_approx type0 = 0;

	/* Check invalid airmass_start */
	hdrldemo_utils_airmass(err,           airmass_end, ra,  dec, lst, exptime, geolat, type1);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

	/* Check invalid airmass_end */
	hdrldemo_utils_airmass(airmass_start, err,         ra,  dec, lst, exptime, geolat, type1);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

	/* Check invalid ra */
	hdrldemo_utils_airmass(airmass_start, airmass_end, err, dec, lst, exptime, geolat, type1);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

	/* Check invalid dec */
	hdrldemo_utils_airmass(airmass_start, airmass_end, ra,  err, lst, exptime, geolat, type1);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

	/* Check invalid lst */
	hdrldemo_utils_airmass(airmass_start, airmass_end, ra,  dec, err, exptime, geolat, type1);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

	/* Check invalid exptime */
	hdrldemo_utils_airmass(airmass_start, airmass_end, ra,  dec, lst, err,     geolat, type1);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

	/* Check invalid geolat */
	hdrldemo_utils_airmass(airmass_start, airmass_end, ra,  dec, lst, exptime, err,    type1);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

	/* Check invalid type */
	hdrldemo_utils_airmass(airmass_start, airmass_end, ra,  dec, lst, exptime, geolat, type0);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);


	/************ Test values ***********/

	hdrl_value airmass1 = hdrldemo_utils_airmass(airmass_start, airmass_end,
									  ra, dec, lst, exptime, geolat, type1);
	cpl_test_abs(airmass1.data,  1.54693,     5e-4);
	cpl_test_abs(airmass1.error, 0.000558673, 1e-6);

	hdrl_value airmass2 = hdrldemo_utils_airmass(airmass_start, airmass_end,
									  ra, dec, lst, exptime, geolat, type2);
	cpl_test_abs(airmass2.data,  1.54633,     5e-4);
	cpl_test_abs(airmass2.error, 0.000551374, 1e-6);


	hdrl_value airmass3 = hdrldemo_utils_airmass(airmass_start, airmass_end,
									  ra, dec, lst, exptime, geolat, type3);
	cpl_test_abs(airmass3.data,  1.54579,     5e-4);
	cpl_test_abs(airmass3.error, 0.000550751, 1e-6);
}
