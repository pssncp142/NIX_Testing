/*
 * This file is part of the HDRLDEMO pipeline
 * Copyright (C) 2013,2014 European Southern Observatory
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
  @defgroup hdrldemo_bpm_2d     Bad Pixel Maps 2D recipe
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_error_code hdrldemo_bpm_2d_save(const char *,
        const char * filename, const cpl_type, const cpl_image *, 
        const cpl_parameterlist *, cpl_frameset *);

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_bpm_2d"

static char hdrldemo_bpm_2d_description[] =
"                                                                           \n"
"The recipe derives bad pixels on individual images.                        \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:                 Explanation:           Required:            \n"
"  RAW                          Data                   Yes                  \n"
"  RAW_BPM                      Bad Pixel Mask         No                   \n"
"  RAW_ERROR                    Associated Error       No                   \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                 Explanation:                                \n"
"  HDRLDEMO_MASTER_BPM          Master bad pixel mask                       \n"
"  HDRLDEMO_MASTER_BPM_FILTERED Grown master bad pixel mask                 \n"
"                                                                           \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"There are currently two methods:                                           \n"
"                                                                           \n"
"(--bpm.method=FILTER): Controlled by                                       \n"
"                  bpm.filter.kappa_low/_high; bpm.filter.maxiter;          \n"
"                  bpm.filter.smooth-x/-y;                                  \n"
"                  bpm.filter; bpm.border;                                  \n"
"                                                                           \n"
"                  The algorithm first smoothes the single image by applying\n"
"                  a cpl-filter (e.g. median), then subtracts the           \n"
"                  smoothed images and derives the bad pixels by            \n"
"                  scaling the rms with kappa_low kappa_high on the         \n"
"                  residual image - for the rms a properly scaled Median    \n"
"                  Absolute Deviation (MAD) is used. The algorithm is       \n"
"                  applied                                                  \n"
"                  iteratively controlled by maxiter. During each           \n"
"                  iteration the newly found bad pixels are ignored.        \n"
"                                                                           \n"
"(--bpm.method=LEGENDRE): Controlled by                                     \n"
"                  bpm.legendre.kappa_low/_high;                            \n"
"                  bpm.legendre.maxiter; bpm.legendre.steps-x/-y;           \n"
"                  bpm.legendre.filter-size-x/-y; bpm.legendre.order-x/-y;  \n"
"                                                                           \n"
"                  The algorithm first smoothes the single image by fitting \n"
"                  a Legendre poly, then subtracts the polynomial and       \n"
"                  derives the bad pixels by scaling the rms with           \n"
"                  kappa_low kappa_high on the residual image - for the     \n"
"                  rms a properly scaled Median Absolute Deviation (MAD)    \n"
"                  is used. The algorithm is applied iteratively controlled \n"
"                  by maxiter. During each iteration the newly found bad    \n"
"                  pixels are ignored.                                      \n"
"                                                                           \n"
"TODO: indicate what data each method is best suited for reduction          \n"
"                                                                           \n"
"The derived bad pixel mask is also filtered by a CLOSING or OPENING        \n"
"filter (see cpl_mask_filter for more details). The filtering process       \n"
"can be controlled by the parameters (--pfx), (--pfy), and (--pfm).         \n"
"                                                                           \n"
"Please note that if no error image is given, the error is estimated        \n"
"with a shot-noise model by using the Gain (--gain) and Ron (--ron)         \n"
"                                                                           \n"
"Please note that part of the code is paralelised. In order to optimise use \n"
"of computing resources you should set the environment variable             \n"
"OMP_NUM_THREADS to a proper value, like (for bash), to use 4 cores         \n"
"export OMP_NUM_THREADS=4                                                   \n";


/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_bpm_2d, HDRLDEMO_BINARY_VERSION, "HDRL Group", 
        PACKAGE_BUGREPORT, "2013", "HDRLDEMO - 2D BPM", 
        hdrldemo_bpm_2d_description);                          

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_bpm_2d_fill_parameterlist(
        cpl_parameterlist   *   self) 
{                                  
    cpl_parameter   *   par ;

    /* --hdrldemo_bpm_2d.ext-nb-raw */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw", CPL_TYPE_INT,
            "FITS extension of the RAW", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-r");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    
    /* --hdrldemo_bpm_2d.ext-nb-raw-err */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-err", CPL_TYPE_INT,
            "FITS extension of the ERROR", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-e");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_2d.ext-nb-raw-bpm */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-bpm", CPL_TYPE_INT,
            "FITS extension or the input BPM", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-b");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_2d.gain */
    par = cpl_parameter_new_value(RECIPE_NAME".gain", CPL_TYPE_DOUBLE,
            "Gain in [e- / ADU]", RECIPE_NAME, 2.5);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_2d.ron */
    par = cpl_parameter_new_value(RECIPE_NAME".ron", CPL_TYPE_DOUBLE,
            "Read-Out Noise", RECIPE_NAME, 1.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ron");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_2d.post-filter-x */
    par = cpl_parameter_new_value(RECIPE_NAME".post-filter-x", CPL_TYPE_INT,
            "X Size of the post filtering kernel", RECIPE_NAME, 3);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfx");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_2d.post-filter-y */
    par = cpl_parameter_new_value(RECIPE_NAME".post-filter-y", CPL_TYPE_INT,
            "Y Size of the post filtering kernel", RECIPE_NAME, 3);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfy");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_2d.post-filter-mode */
    par = cpl_parameter_new_enum(RECIPE_NAME".post-filter-mode", 
            CPL_TYPE_STRING, "Post filtering mode", RECIPE_NAME,
            "closing", 2, "closing", "dilation");
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfm");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_2d.region-llx/lly/urx/ury */
    hdrl_parameter * deflts = hdrl_rect_region_parameter_create(1, 1, 0, 0) ;
    cpl_parameterlist * reg_param = hdrl_rect_region_parameter_create_parlist(
                RECIPE_NAME, "", "region-", deflts);
    hdrl_parameter_delete(deflts) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(reg_param) ; 
            p != NULL; p = cpl_parameterlist_get_next(reg_param)) 
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(reg_param);

    /* Create BPM_2D parameters */
    hdrl_parameter * fil_def =
        hdrl_bpm_2d_parameter_create_filtersmooth(3, 3, 5,
                                                  CPL_FILTER_MEDIAN,
                                                  CPL_BORDER_FILTER,
                                                  3, 3);
    hdrl_parameter * leg_def =
        hdrl_bpm_2d_parameter_create_legendresmooth(3, 3, 5,
                                                    20, 20,
                                                    11, 11, 3, 3);
    cpl_parameterlist * bpm_param = hdrl_bpm_2d_parameter_create_parlist(
                RECIPE_NAME, "bpm", "FILTER", fil_def, leg_def);
    hdrl_parameter_delete(leg_def);
    hdrl_parameter_delete(fil_def);
    for (cpl_parameter * p = cpl_parameterlist_get_first(bpm_param) ; 
            p != NULL; p = cpl_parameterlist_get_next(bpm_param)) 
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(bpm_param);

    return CPL_ERROR_NONE;
}

static int hdrldemo_bpm_2d(
        cpl_frameset            *   frameset, 
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter     *   par = NULL;
    int                         extnum_raw = 0;
    int                         extnum_err = 0;
    int                         extnum_bpm = 0;
    int                         pfx, pfy ;
    const char              *   pfm = NULL;
    cpl_filter_mode             filter_mode = CPL_FILTER_CLOSING ;
    double                      gain = 0.0 ;
    double                      ron = 0.0 ;
    cpl_frame               *   in_frm ;
    cpl_frame               *   err_frm ;
    cpl_frame               *   bpm_frm ;
    hdrl_image              *   hima ;
    hdrl_parameter          *   bpm_params = NULL;
    hdrl_parameter          *   region_params = NULL;
    cpl_mask                *   out_mask ;
    cpl_mask                *   out_mask_filtered = NULL ;
    cpl_image               *   out_image ;

    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /* Get parameters*/
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw");
    extnum_raw = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw-err");
    extnum_err = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw-bpm");
    extnum_bpm = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".gain");
    gain = cpl_parameter_get_double(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ron");
    ron = cpl_parameter_get_double(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".post-filter-x");
    pfx = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".post-filter-y");
    pfy = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".post-filter-mode");
    pfm = cpl_parameter_get_string(par);
    if(!strcmp(pfm, "closing")) {
        filter_mode = CPL_FILTER_CLOSING ;
    } else if(!strcmp(pfm, "dilation")) {
        filter_mode = CPL_FILTER_DILATION ;
    } else {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                "Filter mode can only be \"closing\" or \"dilation\" (not %s)", 
                pfm);
    }

    /* Parse the Region Parameters */
    region_params=hdrl_rect_region_parameter_parse_parlist(parlist, 
            RECIPE_NAME, "region-") ;
    if (region_params == NULL) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                "Parsing of the region parameters failed");
    }

    /* Parse the BPM Parameters */
    bpm_params=hdrl_bpm_2d_parameter_parse_parlist(parlist, RECIPE_NAME".bpm") ;
    if (bpm_params == NULL) {
        if (region_params) hdrl_parameter_delete(region_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                "Parsing of the recipe parameters failed");
    }

    /* Load INPUT Data */
    /* Get the first Required frame */
    if ((in_frm = cpl_frameset_find(frameset, HDRLDEMO_RAW))==NULL) {
        hdrl_parameter_delete(bpm_params) ;
        if (region_params) hdrl_parameter_delete(region_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                "Missing RAW file");
    }
    bpm_frm = cpl_frameset_find(frameset, HDRLDEMO_RAW_BPM) ;
    err_frm = cpl_frameset_find(frameset, HDRLDEMO_RAW_ERROR) ;

    /* Load the image */
    if (hdrldemo_hdrl_image_load(in_frm, extnum_raw, err_frm, extnum_err, 
                bpm_frm, extnum_bpm, region_params, ron, gain, 
                &hima) != CPL_ERROR_NONE) {
        hdrl_parameter_delete(bpm_params) ;
        if (region_params) hdrl_parameter_delete(region_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                "Cannot load RAW image");
    }
    if (region_params) hdrl_parameter_delete(region_params) ;

    /* BPM COMPUTATION */
    if ((out_mask = hdrl_bpm_2d_compute(hima, bpm_params)) == NULL) {
        hdrl_parameter_delete(bpm_params) ;
        hdrl_image_delete(hima);
        return cpl_error_get_code();
    }
    /* Cleanup */
    hdrl_image_delete(hima);
    
    /* Delete the parameters */
    hdrl_parameter_delete(bpm_params) ;

    /* Post Filtering */
    if (pfx > 0 && pfy > 0) {
        out_mask_filtered=hdrl_bpm_filter(out_mask, pfx, pfy, filter_mode);
    }

    /* Save the result */
    /* Convert the mask to image */
    out_image = cpl_image_new_from_mask(out_mask) ;
    cpl_mask_delete(out_mask) ;

    /* Save the image */
    hdrldemo_bpm_2d_save("HDRLDEMO_MASTER_BPM",
            "hdrldemo_bpm_2d_image.fits", CPL_TYPE_INT, out_image, parlist,
            frameset);
    cpl_image_delete(out_image) ;
    if (out_mask_filtered) {

        /* Convert the mask to image */
    	cpl_image *out_image_filtered = cpl_image_new_from_mask(out_mask_filtered) ;
        cpl_mask_delete(out_mask_filtered) ;

        /* Save the filtered image */
        hdrldemo_bpm_2d_save("HDRLDEMO_MASTER_BPM_FILTERED",
							 "hdrldemo_bpm_2d_image_filtered.fits", CPL_TYPE_INT,
							 out_image_filtered, parlist, frameset);

        cpl_image_delete(out_image_filtered) ;
    }
    return (int)cpl_error_get_code();
}

static cpl_error_code hdrldemo_bpm_2d_save(
        const char              *   procatg,
        const char              *   filename,
        const cpl_type              savetype,
        const cpl_image         *   image,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   frameset)
{
    /* Add a QC parameter  */
    cpl_propertylist * qclist = cpl_propertylist_new();

    /* Add the product category and save image */
    cpl_propertylist_update_string(qclist, CPL_DFS_PRO_CATG, procatg);

    cpl_dfs_save_image(frameset, NULL, parlist, frameset, NULL, image, 
            savetype, "hdrldemo_bpm_2d", qclist, NULL, 
            PACKAGE "/" PACKAGE_VERSION, filename);
    cpl_propertylist_delete(qclist);
    return cpl_error_get_code();
}

/**@}*/

