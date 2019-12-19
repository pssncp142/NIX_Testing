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
 * @defgroup hdrldemo_bias  master bias creation
 * @par Synopsis: TBD
 * @par Input frames: TBD
 * @par Output frames: TBD
 * @code
 * @endcode
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_bias"

static char hdrldemo_bias_description[] =
"                                                                           \n"
"The recipe combines raw bias or dark frames into a master frame            \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:                      Explanation:           Required:       \n"
"  (BIAS | DARK)                     Data                   Yes             \n"
"  MASTER_BIAS                       Master Bias            No              \n"
"  MASTER_BIAS_ERROR                 Master Bias error      No              \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                      Explanation:                           \n"
"  MASTER_(BIAS | DARK)              Master-frame                           \n"
"  MASTER_(BIAS | DARK)_ERROR        Error of the Master-frame              \n"
"  MASTER_(BIAS | DARK)_CONTRIBUTION Contribution Map of the Master-frame   \n"
"                                                                           \n"
"                                                                           \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"                                                                           \n"
"The recipe allows to combine dark or bias frames into a master-frame       \n"
"by using one of the implemented collapsation methods, i.e. mean,           \n"
"median, weighted mean, kappa-sigma clipping, and min-max                   \n"
"rejection. The recipe optionally first subtracts the overscan region       \n"
"(controlled by the --oscan.XXX parameter and deactivated by setting        \n"
"--oscan.apply=NO) and then combines the images into a single               \n"
"master-frame. The combination methods can be controlled by the             \n"
"--collapse.XXX recipe parameter.                                           \n"
"                                                                           \n"
"Please note, that if the input frames are dark frames and a (optional)     \n"
"MASTER_BIAS frame is passed to the recipe, the latter is subtracted        \n"
"from each dark frame before combining them into a master-dark              \n"
"                                                                           \n";


/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_bias,
                    HDRLDEMO_BINARY_VERSION,
                    "HDRL Group",
                    PACKAGE_BUGREPORT,
                    "2013",
                    "HDRLDEMO - Master Bias creation",
                    hdrldemo_bias_description);

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_bias_fill_parameterlist(
        cpl_parameterlist   *   self) {
    cpl_parameter   *   par ;

    /* --hdrldemo.bias.ExtensionNumber */
    par = cpl_parameter_new_value(RECIPE_NAME".extension-number", CPL_TYPE_INT,
                                  "FITS extension to load", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo.bias.Gain */
    par = cpl_parameter_new_value(RECIPE_NAME".gain", CPL_TYPE_DOUBLE,
                                  "Gain in [e- / ADU]", RECIPE_NAME, 2.5);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo.bias.oscan.apply */
    par = cpl_parameter_new_enum(RECIPE_NAME".oscan.apply", CPL_TYPE_STRING,
                                 "Apply overscan correction", RECIPE_NAME,
                                 "YES", 2, "YES", "NO");
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "oscan.apply");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* Overscan-correction related parameters */
    hdrl_parameter * rect_region_def = 
        hdrl_rect_region_parameter_create(1,1,20,0) ;
    hdrl_parameter * sigclip_def = 
        hdrl_collapse_sigclip_parameter_create(3., 3., 5);
    hdrl_parameter * minmax_def =
        hdrl_collapse_minmax_parameter_create(1., 1.);
    cpl_parameterlist * os_comp = hdrl_overscan_parameter_create_parlist(
                RECIPE_NAME, "oscan", "alongX", 10, 10., rect_region_def,
                "MEDIAN", sigclip_def, minmax_def);
    hdrl_parameter_delete(rect_region_def); 
    for (cpl_parameter * p = cpl_parameterlist_get_first(os_comp) ;
            p != NULL; p = cpl_parameterlist_get_next(os_comp))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(os_comp);

    /* Master-bias related parameters */
    cpl_parameterlist * pbiascollapse = hdrl_collapse_parameter_create_parlist(
                RECIPE_NAME, "collapse", "MEDIAN", sigclip_def, minmax_def) ;
    hdrl_parameter_delete(sigclip_def); 
    hdrl_parameter_delete(minmax_def);
    for (cpl_parameter * p = cpl_parameterlist_get_first(pbiascollapse) ;
            p != NULL; p = cpl_parameterlist_get_next(pbiascollapse))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(pbiascollapse);

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   perform master bias combination
  @param   frameset   input set of frames
  @param   parlist    input recipe parameters
  @return  cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_bias(
        cpl_frameset            *   frameset, 
        const cpl_parameterlist *   parlist)
{
    const cpl_frame     *   frm_ima;
    cpl_frameset        *   frameset_final = NULL;
    const cpl_parameter *   par;
    int                     is_bias = 0;
    int                     do_os = 1;
    hdrl_parameter      *   os_region_params = NULL;
    hdrl_parameter      *   os_params = NULL;
    hdrl_parameter      *   collapse_params = NULL;

    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /* Get parameters*/
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".extension-number");
    int ext_num = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".gain");
    double gain = cpl_parameter_get_double(par);
 
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".oscan.apply");
    if(!strcmp(cpl_parameter_get_string(par), "YES"))   do_os = 1;
    else                                                do_os = 0;

    /* Parse the Overscan Parameters */
    os_params = hdrl_overscan_parameter_parse_parlist(parlist, RECIPE_NAME".oscan");

    /* Collapse parameters */
    collapse_params = hdrl_collapse_parameter_parse_parlist(parlist, RECIPE_NAME".collapse")
        ;

    cpl_size fsize = cpl_frameset_get_size(frameset);
    cpl_frameset * frameset_bias = cpl_frameset_new();
    cpl_frameset * frameset_dark = cpl_frameset_new();
    cpl_image * master_bias = NULL;
    cpl_image * master_bias_err = NULL;

    for (cpl_size i = 0; i < fsize; i++) {
        cpl_frame * cur_frame=cpl_frameset_get_position(frameset, i);

        if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_BIAS)) {
            cpl_frameset_insert(frameset_bias,
                                cpl_frame_duplicate(cur_frame));
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_DARK)) {
            cpl_frameset_insert(frameset_dark,
                                cpl_frame_duplicate(cur_frame));
        } else if (!strcmp(cpl_frame_get_tag(cur_frame),
                         HDRLDEMO_MASTER_BIAS)) {
            master_bias = cpl_image_load(cpl_frame_get_filename(cur_frame),
                            CPL_TYPE_DOUBLE, 0, 0);
        } else if (!strcmp(cpl_frame_get_tag(cur_frame),
                            HDRLDEMO_MASTER_BIAS_ERROR)) {
            master_bias_err = cpl_image_load(cpl_frame_get_filename(cur_frame),
                                           CPL_TYPE_DOUBLE, 0, 0);
        }
    }


    /* verify that we have at least a bias or a dark input frames */
    if (cpl_frameset_get_size(frameset_bias) > 0 &&
        cpl_frameset_get_size(frameset_dark) > 0) {
        cpl_frameset_delete(frameset_bias);
        cpl_frameset_delete(frameset_dark);
        cpl_image_delete(master_bias);
        cpl_image_delete(master_bias_err);
        hdrl_parameter_destroy(os_params) ;
        hdrl_parameter_delete(collapse_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
               "Bias and Dark frames can not be mixed in the SOF");
    }

    hdrl_image * hdrl_master_bias = NULL;
    if( master_bias != NULL ) {
        hdrl_master_bias = hdrl_image_create(master_bias, master_bias_err);
        cpl_image_delete(master_bias);
        cpl_image_delete(master_bias_err);
    }

    /* are we in case master bias or master dark? */
    if (cpl_frameset_get_size(frameset_bias)==0) {
        frameset_final = cpl_frameset_duplicate(frameset_dark);
        cpl_frameset_delete(frameset_bias);
        cpl_frameset_delete(frameset_dark);
    } else {
        frameset_final = cpl_frameset_duplicate(frameset_bias);
        is_bias = 1;
        cpl_frameset_delete(frameset_bias);
        cpl_frameset_delete(frameset_dark);
    }

    /* Required frame */
    frm_ima = (is_bias == 1) ?
                    cpl_frameset_find_const(frameset_final, HDRLDEMO_BIAS)
                    : cpl_frameset_find_const(frameset_final, HDRLDEMO_DARK);
    if (!frm_ima) {
        hdrl_parameter_destroy(os_params) ;
        hdrl_parameter_delete(collapse_params) ;
        cpl_frameset_delete(frameset_final);
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Cannot find any RAW frame");
    }

    cpl_propertylist * plist =
        cpl_propertylist_load(cpl_frame_get_filename(frm_ima), ext_num);
    cpl_size nx = hdrldemo_get_naxis1(plist);
    cpl_size ny = hdrldemo_get_naxis2(plist);
    cpl_propertylist_delete(plist);

    if (nx <= 0 || ny <= 0) {
        hdrl_parameter_destroy(os_params) ;
        hdrl_parameter_delete(collapse_params) ;
        cpl_frameset_delete(frameset_final);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Could not retrieve image dimensions");
    }

    /* wrap negative indices to positive ones in nx/ny */

    os_region_params = hdrl_overscan_parameter_get_rect_region(os_params) ;
    hdrl_rect_region_fix_negatives(os_region_params, nx, ny);
    if (!cpl_error_get_code()) {
        cpl_msg_info(cpl_func, "Computing overscan on %lld,%lld to %lld,%lld",
                hdrl_rect_region_get_llx(os_region_params),
                hdrl_rect_region_get_lly(os_region_params),
                hdrl_rect_region_get_urx(os_region_params),
                hdrl_rect_region_get_ury(os_region_params)) ;
    }

    /* apply overscan correction and bias normalization from images the input
     * provides, allocate output and write scaled images into them */
    hdrl_imagelist * corrected = hdrl_imagelist_new();
    hdrl_buffer * buf = hdrl_buffer_new();
    /* allow 2GiB of anonymous memory */
    hdrl_buffer_set_malloc_threshold(buf, 2048);
    if(do_os==1) {
        hdrldemo_os_correct(frameset_final, ext_num, os_params,
                            hdrl_master_bias, corrected, buf);
    }
    else {
        hdrldemo_hdrl_imagelist_load(frameset_final, ext_num, NULL, -1, NULL, 
                -1, NULL, hdrl_overscan_parameter_get_ccd_ron(os_params), gain,
                corrected, buf);
    }
    hdrl_parameter_destroy(os_params) ;

    if (cpl_error_get_code()) {
        cpl_frameset_delete(frameset_final);
        hdrl_imagelist_delete(corrected);
        hdrl_parameter_delete(collapse_params) ;
        hdrl_buffer_delete(buf);
        return cpl_error_get_code();
    }

    hdrl_image * master;
    cpl_image * contrib_map;
    /* Get the proper collapsing function and perform frames combination */
    hdrl_imagelist_collapse(corrected, collapse_params, &master, &contrib_map);
    hdrl_parameter_delete(collapse_params);
    hdrl_imagelist_delete(corrected);


    /*Save the results*/
    if (is_bias == 1) {
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_BIAS,
        		"hdrldemo_bias", "hdrldemo_masterbias.fits", CPL_TYPE_FLOAT,
                        hdrl_image_get_image(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_BIAS_ERROR,
        		"hdrldemo_bias", "hdrldemo_masterbias_error.fits", CPL_TYPE_FLOAT,
                            hdrl_image_get_error(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_BIAS_CONTRIBUTION,
        		"hdrldemo_bias", "hdrldemo_masterbias_contribution.fits",
                CPL_TYPE_FLOAT, contrib_map, parlist, frameset);
    } else {
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_DARK,
        		"hdrldemo_bias", "hdrldemo_masterdark.fits", CPL_TYPE_FLOAT,
                hdrl_image_get_image(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_DARK_ERROR,
        		"hdrldemo_bias", "hdrldemo_masterdark_error.fits", CPL_TYPE_FLOAT,
                hdrl_image_get_error(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_DARK_CONTRIBUTION,
        		"hdrldemo_bias", "hdrldemo_masterdark_contribution.fits",
                CPL_TYPE_FLOAT, contrib_map, parlist, frameset);
    }

    /* cleanup */
    hdrl_image_delete(master);
    hdrl_image_delete(hdrl_master_bias);
    cpl_image_delete(contrib_map);

    /* hdrldemo_cube_it_delete(input_it); */
    cpl_frameset_delete(frameset_final);

    hdrl_buffer_delete(buf);
    return (int)cpl_error_get_code();
}

/**@}*/

