/*
 * This file is part of the HDRLDEMO pipeline
 * Copyright (C) 2013 European Southern Observatory
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

#include "hdrldemo_dfs.h"
#include "hdrldemo_utils.h"
#include <cpl.h>

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrldemo_overscan  Overscan Correction
 * @par Synopsis: TBD
 * @par Input frames: TBD
 * @par Output frames: TBD
 * @code
 * @endcode
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_error_code hdrldemo_overscan_save_image(const char *,
        const char * filename, const cpl_type, const cpl_image *, 
        const cpl_parameterlist *, cpl_frameset *);

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_overscan"

static char hdrldemo_overscan_description[] = "  TBD";

/* ---------------------------------------------------------------------------*/
/**
 * @brief merges bpm while translating codes
 *
 * @param bmap     input bpm to merge with
 * @param bmap_src input bpm to be translated and merged
 * @param src      source code
 * @param dst      destination code
 *
 * performs: bmap[(bmap_src & src) == src] |= dst
 * no error checking is done
 *
 * @return cpl_error_code
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code translate_bpm(cpl_image * bmap, cpl_image * bmap_src,
                                    hdrl_bitmask_t src, hdrl_bitmask_t dst)
{
    int * ibmap = cpl_image_get_data_int(bmap);
    int * ibmap_src = cpl_image_get_data_int(bmap_src);
    const size_t npix = (size_t)(cpl_image_get_size_x(bmap) * cpl_image_get_size_y(bmap));
    for (size_t i = 0; i < npix; i++) {
        if ((ibmap_src[i] & src) == src)
            ibmap[i] |= dst;
    }
    return cpl_error_get_code();
}

/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_overscan, 
                    HDRLDEMO_BINARY_VERSION, 
                    "HDRL Group", 
                    PACKAGE_BUGREPORT, 
                    "2013", 
                    "HDRLDEMO - Overscan correction",
                    hdrldemo_overscan_description);                          

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_overscan_fill_parameterlist(
        cpl_parameterlist   *   self) {                                  
    cpl_parameter   *   par ;

    /* --hdrldemo.overscan.extension-number */
    par = cpl_parameter_new_value(RECIPE_NAME".extension-number", CPL_TYPE_INT,
                                  "FITS extension to load", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo.overscan.flag-no-cor-value */
    par = cpl_parameter_new_value(RECIPE_NAME".flag-no-cor-value",
                                  CPL_TYPE_INT,
                                  "Power of two value used to mark pixels not "
                                  "corrected by overscan correction due being "
                                  "unable to determine a value from the "
                                  "overscan region", RECIPE_NAME, 256);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "flag-no-cor-value");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* Create Overscan Computation parameters */
    hdrl_parameter * rect_region_def =
        hdrl_rect_region_parameter_create(1, 1, 20, 1000);
    hdrl_parameter * sigclip_def =
        hdrl_collapse_sigclip_parameter_create(3., 3., 5);
    hdrl_parameter * minmax_def =
        hdrl_collapse_minmax_parameter_create(1., 1.);
    cpl_parameterlist * os_comp = hdrl_overscan_parameter_create_parlist(
                RECIPE_NAME, "", "alongX", 0, 0., rect_region_def, "MEDIAN",
                sigclip_def, minmax_def);
    hdrl_parameter_delete(rect_region_def);
    hdrl_parameter_delete(sigclip_def);
    hdrl_parameter_delete(minmax_def);
    for (cpl_parameter * p = cpl_parameterlist_get_first(os_comp) ; 
            p != NULL; p = cpl_parameterlist_get_next(os_comp)) 
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(os_comp);

    /* Create Overscan Correction Region parameters */
    /* --hdrldemo.overscan.apply-llx/lly/urx/ury */
    hdrl_parameter * deflts = hdrl_rect_region_parameter_create(1, 1, 0, 0);
    cpl_parameterlist * os_corr_reg = hdrl_rect_region_parameter_create_parlist(
                RECIPE_NAME, "corr", "apply-", deflts);
    hdrl_parameter_delete(deflts) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(os_corr_reg) ;       
            p != NULL; p = cpl_parameterlist_get_next(os_corr_reg)) 
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(os_corr_reg);

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrldemo_overscan  Overscan Correction
 *
 * This module provides functions to correct bias contribute by evaluating the
 * detector's overscan level
 *
 * @par Synopsis:
 *   Functions to evaluate and correct overscan detector's level
 *
 * @code
 * @endcode
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief  
  @param
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_overscan(
        cpl_frameset            *   frameset, 
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter          *par            = NULL;
    int                          extnum          = 0;
    hdrl_bitmask_t               flagnocor;
    cpl_image                    *image          = NULL;
    cpl_image                    *ima_error      = NULL;
    const cpl_frame              *frm_ima_error  = NULL;
    const cpl_frame              *frm_ima        = NULL;
    hdrl_parameter               *os_params      = NULL;
    hdrl_overscan_compute_result *os_computation = NULL;
    hdrl_overscan_correct_result *os_correction  = NULL;

    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /* Get parameters*/
    par = cpl_parameterlist_find_const(parlist,
                                       RECIPE_NAME".extension-number");
    extnum = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist,
                                       RECIPE_NAME".flag-no-cor-value");
    flagnocor = cpl_parameter_get_int(par);
    if (flagnocor & (flagnocor - 1)) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                "--flag-no-cor-value needs to be a power of two");
    }

    /* Parse the Overscan Parameters */
    os_params = hdrl_overscan_parameter_parse_parlist(parlist, RECIPE_NAME);
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        hdrl_parameter_destroy(os_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                        "Parsing of the recipe parameters failed");
    }

    /*Reqired frame*/
    frm_ima = cpl_frameset_find_const(frameset, HDRLDEMO_RAW);
    if(!frm_ima){
        hdrl_parameter_destroy(os_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                        "Cannot find any RAW frame");
    }

    {
        /*Optional frame*/
        const cpl_frame * frm_mask = cpl_frameset_find_const(frameset, 
                HDRLDEMO_MASTER_BPM);
        cpl_mask * mask = NULL;
        if(frm_mask != NULL) {
            mask = cpl_mask_load(cpl_frame_get_filename(frm_mask), 0, extnum);
        }

        image = cpl_image_load(cpl_frame_get_filename(frm_ima), CPL_TYPE_FLOAT,
                        0, extnum);
        if (!image) {
            cpl_mask_delete(mask);
            hdrl_parameter_destroy(os_params) ;
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                            "Cannot load RAW frame");
        }

        /*
        mask=cpl_mask_threshold_image_create(image,-100,270);
        cpl_mask_save(mask,"mask.fits", NULL, CPL_IO_CREATE);
        */

        if(mask != NULL) {
            cpl_image_reject_from_mask(image, mask);
            cpl_mask_delete(mask);
        }
    }

    /* Overscan Computation */
    os_computation = hdrl_overscan_compute(image, os_params);
    if(!os_computation) {
    	hdrl_parameter_destroy(os_params);
        cpl_image_delete(image);
        return cpl_error_get_code();
    }

    cpl_image *rej_low  = NULL;
    cpl_image *rej_high = NULL;
    if (   hdrl_collapse_parameter_is_sigclip(hdrl_overscan_parameter_get_collapse(os_params))
    	|| hdrl_collapse_parameter_is_minmax( hdrl_overscan_parameter_get_collapse(os_params)) )
    {
    	rej_low  = hdrl_overscan_compute_result_get_minmax_reject_low( os_computation);
		rej_high = hdrl_overscan_compute_result_get_minmax_reject_high(os_computation);
    }

    /* Delete the parameters */
    hdrl_parameter_destroy(os_params) ;

    /* Optional errorframe */
    frm_ima_error = cpl_frameset_find_const(frameset, HDRLDEMO_RAW_ERROR);
    if(frm_ima_error != NULL) {
        ima_error = cpl_image_load(cpl_frame_get_filename(frm_ima_error),
                        CPL_TYPE_FLOAT, 0, extnum);
    } else {
        ima_error = cpl_image_duplicate(image);
        cpl_image_multiply_scalar(ima_error, 0.0);
    }

    if(!ima_error) {
        cpl_image_delete(image);
        return cpl_error_get_code();
    }
    hdrl_image * himage = hdrl_image_create(image, ima_error);
    cpl_image_delete(image);
    cpl_image_delete(ima_error);

    {
        hdrl_parameter * os_corr = hdrl_rect_region_parameter_parse_parlist(
                parlist, RECIPE_NAME".corr", "apply-") ;
        hdrl_rect_region_fix_negatives(os_corr, hdrl_image_get_size_x(himage),
                                       hdrl_image_get_size_y(himage));
        os_correction = hdrl_overscan_correct(himage, os_corr, os_computation);
        hdrl_parameter_delete(os_corr);
    }
    hdrl_image_delete(himage);
    if(!os_correction) {
    	hdrl_overscan_compute_result_delete(os_computation);
    	return cpl_error_get_code();
    }

    /*Save the results*/
    cpl_image *rImg;

    hdrldemo_overscan_save_image("HDRLDEMO_OVERSCAN",
                   "hdrldemo_overscan.fits",                 CPL_TYPE_FLOAT,
                   hdrl_image_get_image(
                   hdrl_overscan_compute_result_get_correction(os_computation)),
                   parlist, frameset);

    hdrldemo_overscan_save_image("HDRLDEMO_OVERSCAN_ERROR",
                   "hdrldemo_overscan_error.fits",           CPL_TYPE_FLOAT,
                   hdrl_image_get_error(
                   hdrl_overscan_compute_result_get_correction(os_computation)),
                   parlist, frameset);

    hdrl_image *r = hdrl_overscan_compute_result_unset_correction(os_computation);
    hdrl_image_delete(r);

    hdrldemo_overscan_save_image("HDRLDEMO_OVERSCAN_CONTRIBUTION",
                   "hdrldemo_overscan_contribution.fits",    CPL_TYPE_FLOAT,
                   hdrl_overscan_compute_result_get_contribution(os_computation),
                   parlist, frameset);
    rImg = hdrl_overscan_compute_result_unset_contribution(os_computation);
    cpl_image_delete(rImg);

    hdrldemo_overscan_save_image("HDRLDEMO_OVERSCAN_CHI2",
                   "hdrldemo_overscan_chi2.fits",         CPL_TYPE_FLOAT,
                   hdrl_overscan_compute_result_get_chi2(os_computation),
                   parlist, frameset);
    rImg = hdrl_overscan_compute_result_unset_chi2(os_computation);
    cpl_image_delete(rImg);

    hdrldemo_overscan_save_image("HDRLDEMO_OVERSCAN_REDCHI2",
                   "hdrldemo_overscan_redchi2.fits",         CPL_TYPE_FLOAT,
                   hdrl_overscan_compute_result_get_red_chi2(os_computation),
                   parlist, frameset);
    rImg = hdrl_overscan_compute_result_unset_red_chi2(os_computation);
    cpl_image_delete(rImg);

    if (rej_low) {
        hdrldemo_overscan_save_image("HDRLDEMO_OVERSCAN_REJ_LOW",
           "hdrldemo_overscan_rej_low.fits", CPL_TYPE_FLOAT, rej_low, parlist, frameset);
        rImg = hdrl_overscan_compute_result_unset_minmax_reject_low(os_computation);
        cpl_image_delete(rImg);
    }

    if (rej_high) {
        hdrldemo_overscan_save_image("HDRLDEMO_OVERSCAN_REJ_HIGH",
           "hdrldemo_overscan_rej_high.fits", CPL_TYPE_FLOAT, rej_high, parlist, frameset);
        rImg = hdrl_overscan_compute_result_unset_minmax_reject_high(os_computation);
        cpl_image_delete(rImg);
    }

    hdrldemo_overscan_save_image("HDRLDEMO_OVERSCANCORRECTED",
                   "hdrldemo_overscancorrected.fits",        CPL_TYPE_FLOAT,
                   hdrl_image_get_image(
                    hdrl_overscan_correct_result_get_corrected(os_correction)),
                   parlist, frameset);

    hdrldemo_overscan_save_image("HDRLDEMO_OVERSCANCORRECTED_ERROR",
                   "hdrldemo_overscancorrected_error.fits",  CPL_TYPE_FLOAT,
                   hdrl_image_get_error(
                    hdrl_overscan_correct_result_get_corrected(os_correction)),
                   parlist, frameset);

    {
        cpl_image * bpm =
            hdrl_overscan_correct_result_get_badmask(os_correction);
        /* new empty bpm, would normally be the pipeline bpm */
        cpl_image * tbpm = cpl_image_duplicate(bpm);
        cpl_image_multiply_scalar(tbpm, 0);

        /* translate code 1 to user requested value */
        translate_bpm(tbpm, bpm, 1, flagnocor);
        hdrldemo_overscan_save_image("HDRLDEMO_OVERSCANCORRECTED_BITMASK",
                                     "hdrldemo_overscancorrected_bitmask.fits",
                                     CPL_TYPE_INT, tbpm, parlist, frameset);
        cpl_image_delete(tbpm);

        cpl_image *rImg2 = hdrl_overscan_correct_result_unset_badmask(os_correction);
        cpl_image_delete(rImg2);
    }

    /*Release the memory*/
    hdrl_overscan_compute_result_delete(os_computation);
    hdrl_overscan_correct_result_delete(os_correction);

    return (int)cpl_error_get_code();
}

static cpl_error_code
hdrldemo_overscan_save_image(const char * procatg,
                const char * filename,
                const cpl_type savetype,
                const cpl_image * image,
                const cpl_parameterlist* parlist,
                cpl_frameset* frameset)
{
    /* Add a QC parameter  */
    cpl_propertylist * qclist = cpl_propertylist_new();

    /* Add the product category and save image */
    cpl_propertylist_update_string(qclist, CPL_DFS_PRO_CATG, procatg);

    cpl_dfs_save_image(frameset, NULL, parlist, frameset, NULL,
                    image, savetype,
                    "hdrldemo_overscan", qclist, NULL,
                    PACKAGE "/" PACKAGE_VERSION, filename);

    cpl_propertylist_delete(qclist);

    return cpl_error_get_code();
}

/**@}*/

