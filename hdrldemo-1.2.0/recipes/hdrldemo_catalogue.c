/*
 * This file is part of the HDRLDEMO pipeline
 * Copyright (C) 2014 European Southern Observatory
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
  @defgroup hdrldemo_catalogue catalogue computation recipe
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_catalogue_save(
        const hdrl_catalogue_result *   res,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   frameset);

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_catalogue"

static char hdrldemo_catalogue_description[] =
"                                                                           \n"
"The recipe extracts a source catalog from an astronomical FITS image       \n"
"Objects are detected and parameterised using the processed images and      \n"
"confidence maps.                                                           \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:                 Explanation:           Required:            \n"
"  RAW                          Input Data             Yes                  \n"
"  RAW_CONFMAP                  Confidence Map         No                   \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                 Explanation:                                \n"
"  HDRLDEMO_CATALOGUE           Source Catalogue                            \n"
"  HDRLDEMO_CATALOGUE_BACKMAP   Background map                              \n"
"  HDRLDEMO_CATALOGUE_SEGMAP    Segmentation map                            \n"
"                                                                           \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"                                                                           \n"
"In brief, the recipe first estimates the local sky background over the     \n"
"field and track any variations at adequate resolution to eventually        \n"
"remove them (if bkg.estimate=TRUE). Then it detects objects/blends of      \n"
"objects and keep a list of pixels belonging to each blend for further      \n"
"analysis. Finally, it parameter the detected objects, i.e. performs        \n"
"astrometry, photometry and some sort of shape analysis.                    \n"
"                                                                           \n"
"Note that pixels are weighted in the catalogue generation according to     \n"
"their value in the confidence map.  Hence if a pixel is marked as bad,     \n"
"then it is not included in the aperture flux.  The number of bad           \n"
"pixels with an aperture is reported in the 'Error_bit_flag' column of      \n"
"the output table.  The presence of bad pixels will also be reflected       \n"
"in the average confidence for the aperture (column 'Av_conf').             \n"
"                                                                           \n"
"                                                                           \n"
"Recipe parameter to characterize the detector:                             \n"
"                                                                           \n"
"--det.effective-gain: Effective detector gain value to convert             \n"
"                      intensity to electrons.                              \n"
"--det.saturation:     Detector saturation value. All pixel values          \n"
"                      exceeding this threshold before background           \n"
"                      subtraction are marked as bad pixels.                \n"
"                                                                           \n"
"Recipe parameter controlling the background determination and              \n"
"subtraction:                                                               \n"
"                                                                           \n"
"--bkg.estimate:  If set, the recipe estimates and subtracts the            \n"
"                 background from the input image, else, no background      \n"
"                 determination and subtraction is done.                    \n"
"                                                                           \n"
"--bkg.mesh-size: It defines the mesh size in pixels to construct the       \n"
"                 low resolution background by removing the objects. If     \n"
"                 the mesh-size is to small, object residuals may           \n"
"                 distort the background.                                   \n"
"                                                                           \n"
"--bkg.smooth-gauss-fwhm: The FWHM of the Gaussian kernel used to           \n"
"                         convolve the background before doing object       \n"
"                         detection.                                        \n"
"                                                                           \n"
"Recipe parameter controlling the background determination and              \n"
"subtraction:                                                               \n"
"                                                                           \n"
"                                                                           \n"
"--obj.min-pixels:  Minimum pixel area (contiguous) of an object to be      \n"
"                   detected.                                               \n"
"--obj.threshold:   Minimum detection threshold in sigma above sky for a    \n"
"                   pixel to be considered for detection.                   \n"
"--obj.deblending:  If set, the algorithm tries to deblend close-by         \n"
"                   objects.                                                \n"
"--obj.core-radius: Value of the core radius. The core radius, among        \n"
"                   others, determines the fixed aperture flux              \n"
"                   measurements - it uses the core radius multiplied       \n"
"                   by predefined factor.                                   \n"
"                                                                           \n"
"Recipe parameter to control the recipe output:                             \n"
"                                                                           \n"
"--output: Requested output, comma separated: Catalogue (CAT),              \n"
"          Background map (BKG), segmentation map (SEGMAP) or all three     \n"
"          outputs (ALL).                                                   \n"
"                                                                           \n"
"Please note that part of the code is paralelised. In order to optimise     \n"
"use of computing resources you should set the environment variable         \n"
"OMP_NUM_THREADS to a proper value, like (for bash), to use 4 cores         \n"
"export OMP_NUM_THREADS=4                                                   \n"
"                                                                           \n";


/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_catalogue, HDRLDEMO_BINARY_VERSION, "HDRL Group", 
        PACKAGE_BUGREPORT, "2016", "HDRLDEMO - catalogue",
        hdrldemo_catalogue_description);                          

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_catalogue_fill_parameterlist(
        cpl_parameterlist   *   self) 
{                                  
    cpl_parameter   *   par ;

    /* --hdrldemo_catalogue.ext-nb-raw */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw", CPL_TYPE_INT,
            "FITS extension of the RAW", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-r");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    
    /* --hdrldemo_catalogue.ext-nb-raw-cnf */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-cnf", CPL_TYPE_INT,
            "FITS extension of the confidence map", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-c");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    hdrl_parameter * c_def =
        hdrl_catalogue_parameter_create(4, 2.5, CPL_TRUE, 5.0, CPL_TRUE, 64,
                                        2., 3.0, HDRL_SATURATION_INIT, 1.0);

    cpl_parameterlist * s_param =
        hdrl_catalogue_parameter_create_parlist(RECIPE_NAME, "", c_def);
    hdrl_parameter_delete(c_def) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(s_param);
            p != NULL; p = cpl_parameterlist_get_next(s_param))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(s_param);

    par = cpl_parameter_new_value(RECIPE_NAME".output", CPL_TYPE_STRING,
            "Requested output, comma separated: CAT, BKG, SEGMAP or ALL",
            RECIPE_NAME, "ALL");
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "output");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    return CPL_ERROR_NONE;
}

static int hdrldemo_catalogue(
        cpl_frameset            *   frameset, 
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *par       = NULL;
    int                 extnum_raw = 0;
    int                 extnum_cnf = 0;
    const char          *output    = NULL;
    hdrl_iter           *it        = NULL;

    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /* Get parameters*/
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw");
    extnum_raw = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw-cnf");
    extnum_cnf = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".output");
    output = cpl_parameter_get_string(par);
    hdrl_parameter * p = hdrl_catalogue_parameter_parse_parlist(parlist,
                                                                RECIPE_NAME);
    /* catalogue currently is always created */
    hdrl_catalogue_options opt = HDRL_CATALOGUE_CAT_COMPLETE;
    if (strstr(output, "CAT")) {
        opt |= HDRL_CATALOGUE_CAT_COMPLETE;
    }
    if (strstr(output, "BKG")) {
        opt |= HDRL_CATALOGUE_BKG;
    }
    if (strstr(output, "SEGMAP")) {
        opt |= HDRL_CATALOGUE_SEGMAP;
    }
    if (strstr(output, "ALL")) {
        opt |= HDRL_CATALOGUE_ALL;
    }
    hdrl_catalogue_parameter_set_option(p, opt);

    cpl_frameset * img_frameset = cpl_frameset_new();
    cpl_frameset * cnf_frameset = cpl_frameset_new();
    /* Classify the input frames into images, errors and object masks */
    for (cpl_size i = 0; i < cpl_frameset_get_size(frameset); i++) {
        const cpl_frame * cur_frame;
        const char      * cur_tag;

        cur_frame = cpl_frameset_get_position_const(frameset, i);
        cur_tag = cpl_frame_get_tag(cur_frame);

        if (!strcmp(cur_tag, HDRLDEMO_RAW)){
            cpl_frameset_insert(img_frameset,
                                cpl_frame_duplicate(cur_frame));
        }
        else if (!strcmp(cur_tag, HDRLDEMO_RAW_CONFMAP)) {
            cpl_frameset_insert(cnf_frameset,
                                cpl_frame_duplicate(cur_frame));
        }
        else {
            cpl_msg_warning(cpl_func,
                            "%s: Unsupported tag %s",
                            cpl_frame_get_filename(cur_frame),
                            cpl_frame_get_tag(cur_frame));
        }
    }

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        goto cleanup;
    }

    /* catalogue COMPUTATION */

    intptr_t dims[] = {-1, 1};
    intptr_t axes[] = {HDRL_FRAMEITER_AXIS_FRAME, HDRL_FRAMEITER_AXIS_EXT};
    hdrl_iter * subiters[] = {
        hdrl_frameiter_new(img_frameset, HDRL_ITER_OWNS_DATA, 2, axes,
                           (intptr_t[]){0, extnum_raw}, NULL, dims),
        hdrl_frameiter_new(cnf_frameset, HDRL_ITER_OWNS_DATA, 2, axes,
                           (intptr_t[]){0, extnum_cnf}, NULL, NULL),
    };
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        /* TODO multiiter taking ownership and cleanup on previous error might
         * be nicer */
        hdrl_iter_delete(subiters[0]);
        hdrl_iter_delete(subiters[1]);
        goto cleanup;
    }

    it = hdrl_multiiter_new(2, subiters, HDRL_ITER_ALLOW_EMPTY);
    if (hdrl_iter_length(it) == 0) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                              "No correctly tagged input in SOF");
        goto cleanup;
    }
    for (hdrl_frameiter_data ** h = hdrl_iter_next(it); h != NULL;
         h = hdrl_iter_next(it)) {
        if (h[0] == NULL) {
            continue;
        }

        /* Check if file has wcs information */
        cpl_errorstate  prestate = cpl_errorstate_get();
        cpl_wcs *wcs = cpl_wcs_new_from_propertylist(h[0]->plist);
        if (!cpl_errorstate_is_equal(prestate)) {
            cpl_msg_warning(cpl_func, "Input image has now WCS information");
            wcs = NULL;
            /* Reset error code */
            cpl_errorstate_set(prestate);
        }
        hdrl_catalogue_result * res = hdrl_catalogue_compute(
        		h[0]->image, h[1] ? h[1]->image : NULL, wcs, p);
        hdrldemo_catalogue_save(res, parlist, frameset);

        hdrl_catalogue_result_delete(res);
        cpl_wcs_delete(wcs);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            goto cleanup;
        }
    }


cleanup:
    hdrl_iter_delete(it);
    hdrl_parameter_delete(p);
    cpl_frameset_delete(img_frameset);
    cpl_frameset_delete(cnf_frameset);

    return (int)cpl_error_get_code();
}

static cpl_error_code hdrldemo_catalogue_save(
        const hdrl_catalogue_result *   res,
        const cpl_parameterlist     *   parlist,
        cpl_frameset                *   frameset)
{
    cpl_propertylist * applist;
    /* Add a QC parameter  */
    if (res->qclist && res->catalogue) {


		applist = cpl_propertylist_duplicate(res->qclist);
		cpl_propertylist_update_string(applist, CPL_DFS_PRO_CATG, "HDRLDEMO_CATALOGUE");

        cpl_dfs_save_table(frameset, NULL, parlist, frameset, NULL, res->catalogue,
                           res->qclist, "hdrldemo_catalogue", applist, NULL,
                           PACKAGE "/" PACKAGE_VERSION, "hdrldemo_catalogue.fits");

    } else {
        applist = cpl_propertylist_new();
    }

    if (res->segmentation_map) {
        cpl_propertylist_update_string(applist, CPL_DFS_PRO_CATG,
                                       "HDRLDEMO_CATALOGUE_SEGMAP");
        cpl_dfs_save_propertylist(frameset, applist, parlist, frameset, NULL,
                                  "hdrldemo_catalogue", applist, NULL,
                                  PACKAGE "/" PACKAGE_VERSION,
                                  "hdrldemo_catalogue_segmap.fits");
        cpl_image_save(res->segmentation_map, "hdrldemo_catalogue_segmap.fits",
                       CPL_TYPE_INT, NULL, CPL_IO_EXTEND |
                       CPL_IO_COMPRESS_RICE);
    }

    if (res->background) {
        cpl_propertylist_update_string(applist, CPL_DFS_PRO_CATG,
                                       "HDRLDEMO_CATALOGUE_BACKMAP");
        cpl_dfs_save_image(frameset, NULL, parlist, frameset, NULL, res->background,
                           CPL_TYPE_DOUBLE, "hdrldemo_catalogue", applist, NULL,
                           PACKAGE "/" PACKAGE_VERSION,
                           "hdrldemo_catalogue_backmap.fits");
    }

    cpl_propertylist_delete(applist);
    return cpl_error_get_code();
}


/**@}*/

