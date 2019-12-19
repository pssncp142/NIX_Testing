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
#include "hdrl_fringe.h" // TODO remove once in hdrl.h !!

#include "hdrldemo_dfs.h"
#include "hdrldemo_utils.h"

#include <cpl.h>

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrldemo_fringe  master fringe creation
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

#define RECIPE_NAME "hdrldemo_fringe"

static char hdrldemo_fringe_description[] =
"                                                                           \n"
"The recipe derives or subtracts a master-fringe image.                     \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:                      Explanation:           Required:       \n"
"  RAW                               Data                   Yes             \n"
"  RAW_BPM                           Bad pixel masks        NO              \n"
"  OBJ_MASK                          Object masks           NO              \n"
"  STATIC_MASK                       Static mask            NO              \n"
"  MASTER_FRINGE                     Master-fringe          NO              \n"
"  MASTER_FRINGE_ERROR               Error of master-fringe NO              \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                      Explanation:                           \n"
"  MASTER_FRINGE                     Master-fringe                          \n"
"  MASTER_FRINGE_ERROR               Error of master-fringe                 \n"
"  MASTER_FRINGE_CONTRIBUTION        Pixel contribution map                 \n"
"                                                                           \n"
"  or:                                                                      \n"
"                                    Explanation:                           \n"
"  FRINGE_CORRECTED                  Fringe corrected images                \n"
"  FRINGE_CORRECTED_ERROR            Error of the fringe corrected images   \n"
"                                                                           \n"
"                                                                           \n"
"This recipe has polymorphic behaviour, combining two                       \n"
"functionalities. It computes OR subtracts a master-fringe image from a     \n"
"set of dithered input images depending on the presence or absence of a     \n"
"master-fringe image in the SoF. If there is NO master-fringe image         \n"
"(tagged as MASTER_FRINGE) in the SoF a master-fringe image will be         \n"
"calculated from the input images. If there is a master-fringe image        \n"
"the latter will be properly scaled and subtracted from the input           \n"
"images. This allows the user to e.g. compute the master-fringe from        \n"
"one (larger) data-set and apply the result to another data-set.            \n"
"                                                                           \n"
"  Usage of the recipe:                                                     \n"
"  There are currently two modes:                                           \n"
"                                                                           \n"
"                                                                           \n"
"  - Computing a master-fringe image (no MASTER_FRINGE in the SoF):         \n"
"                                                                           \n"
"     The algorithm consists of two main steps: i) Using a Gaussian         \n"
"     mixture model the recipe estimates the background and the fringe      \n"
"     scaling-factor (amplitude) of each image tagged as RAW. The           \n"
"     background is then subtracted and the amplitude used to normalize     \n"
"     the background-subtracted image (multiplicatively) to the same        \n"
"     scale. ii) The resulting images are then collapsed into a             \n"
"     master-fringe image. There are different collapsing methods           \n"
"     available and they are controlled by the -collapse.XXX recipe         \n"
"     parameters. In other words: the amplitude of the fringes is           \n"
"     computed for each input image (RAW) and used to re-scale the          \n"
"     image before stacking. For convenience, the error images              \n"
"     associated to the RAW images are derived by a shot-noise model        \n"
"     using the values passed with the recipe parameters --gain and         \n"
"     --ccd-ron.                                                            \n"
"                                                                           \n"
"     Different masks controlling the algorithm in different stages can     \n"
"     be passed to the recipe:                                              \n"
"                                                                           \n"
"     A) RAW_BPM: One for each RAW image. Bad pixel mask of each single     \n"
"     input image. Each step in the algorithm takes this mask into          \n"
"     account.                                                              \n"
"                                                                           \n"
"     B) OBJ_MASK: One for each RAW image. Images flagging the objects      \n"
"     in the input image. Each step in the algorithm takes this mask        \n"
"     into account.                                                         \n"
"                                                                           \n"
"     C) STATIC_MASK: Single mask. This mask is only taken into account     \n"
"     in step i), i.e. when computing the image background and fringe       \n"
"     amplitude for each RAW image. Step ii), i.e. the collapsing           \n"
"     algorithm ignores this mask.                                          \n"
"                                                                           \n"
"                                                                           \n"
"  - Subtracting a master-fringe image: (MASTER_FRINGE in the SoF):         \n"
"                                                                           \n"
"     The algorithm consists of two main steps: i) The recipe estimates     \n"
"     the background and the fringe scaling-factor (amplitude) for the      \n"
"     master-fringe image (MASTER_FRINGE) of each image tagged as           \n"
"     RAW. The required offset and scaling factor is determined using a     \n"
"     least-square fit of the master-fringe image to the RAW image.         \n"
"     ii) The properly scaled master-fringe is then subtracted from         \n"
"     each individual RAW frame. For convenience, the error images          \n"
"     associated to the RAW images are derived by a shot-noise model        \n"
"     using the values passed with the recipe parameters --gain and         \n"
"     --ccd-ron.                                                            \n"
"                                                                           \n"
"     Different masks controlling the algorithm in different stages can     \n"
"     be passed to the recipe:                                              \n"
"                                                                           \n"
"     A) RAW_BPM: One for each RAW image. Bad pixel mask of each single     \n"
"     input image. Each step in the algorithm takes this mask into          \n"
"     account.                                                              \n"
"                                                                           \n"
"     B) OBJ_MASK: One for each RAW image. Images flagging the objects      \n"
"     in the input image. These masks are only taken into account in        \n"
"     step i), i.e. when computing the scaling factors between              \n"
"     MASTER_FRINGE and RAW, but ignored for step ii).                      \n"
"                                                                           \n"
"     C) STATIC_MASK: Single mask. This mask is only taken into account     \n"
"     in step i), i.e. when computing the scaling factors between           \n"
"     MASTER_FRINGE and RAW, but ignored for step ii).                      \n"
"                                                                           \n"
"Please note that part of the code is paralelised. In order to optimise     \n"
"use of computing resources you should set the environment variable         \n"
"OMP_NUM_THREADS to a proper value, like (for bash), to use 4 cores         \n"
"export OMP_NUM_THREADS=4                                                   \n"
"                                                                           \n";

/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_fringe,
                  HDRLDEMO_BINARY_VERSION,
                  "HDRL Group",
                  PACKAGE_BUGREPORT,
                  "2015",
                  "HDRLDEMO - Master Fringe creation / correction",
                  hdrldemo_fringe_description);

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_fringe_fill_parameterlist(
                cpl_parameterlist   *   self) {
    cpl_parameter   *   par ;


    /* --hdrldemo_fringe.ext-nb-raw */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw", CPL_TYPE_INT,
            "FITS extension of the RAW images", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-r");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    
    /* --hdrldemo_fringe.ext-nb-raw-bpm */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-bpm", CPL_TYPE_INT,
            "FITS extension or the Bad Pixel Masks", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-b");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_fringe.ext-nb-raw-obj */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-obj", CPL_TYPE_INT,
            "FITS extension or the Object Masks", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-obj");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_fringe.ext-nb-raw-stat */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-stat", CPL_TYPE_INT,
            "FITS extension or the Static Mask", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-stat");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo.fringe.Gain */
    par = cpl_parameter_new_value(RECIPE_NAME".gain", CPL_TYPE_DOUBLE,
                    "Gain in [e- / ADU]", RECIPE_NAME, 2.5);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo.fringe.ccd-ron */
    par = cpl_parameter_new_value(RECIPE_NAME".ccd-ron", CPL_TYPE_DOUBLE,
                    "Readout noise in ADU.", RECIPE_NAME, 10.);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ccd-ron");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* Master-fringe related parameters */

    hdrl_parameter * sigclip_def =
                    hdrl_collapse_sigclip_parameter_create(3., 3., 5);
    hdrl_parameter * minmax_def =
                    hdrl_collapse_minmax_parameter_create(1., 1.);

    cpl_parameterlist * pfringecollapse = hdrl_collapse_parameter_create_parlist(
                    RECIPE_NAME, "collapse", "MEDIAN", sigclip_def, minmax_def);
    hdrl_parameter_delete(sigclip_def); 
    hdrl_parameter_delete(minmax_def);
    for (cpl_parameter * p = cpl_parameterlist_get_first(pfringecollapse) ;
                    p != NULL; p = cpl_parameterlist_get_next(pfringecollapse))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(pfringecollapse);

    /* Create fringe parameters */


    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
  @brief   perform master fringe combination
  @param   frameset   input set of frames
  @param   parlist    input recipe parameters
  @return  cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_fringe(
                cpl_frameset            *   frameset,
                const cpl_parameterlist *   parlist)
{
    const cpl_frame     *   frm_ima;
    const cpl_parameter *   par;
    hdrl_parameter      *   collapse_params;

    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /* Get parameters*/

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw");
    int extnum_raw = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw-bpm");
    int extnum_bpm = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw-obj");
    int extnum_obj = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw-stat");
    int extnum_stat = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".gain");
    double gain = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ccd-ron");
    double ron = cpl_parameter_get_double(par);

    /* Collapse parameters */
    collapse_params = hdrl_collapse_parameter_parse_parlist(parlist,
                    RECIPE_NAME".collapse");
    /* Parse the FRINGE Parameters */

    cpl_size fsize = cpl_frameset_get_size(frameset);
    cpl_frameset * frameset_fringe = cpl_frameset_new();
    cpl_frameset * frameset_fringe_bpm = cpl_frameset_new();
    cpl_frameset * frameset_fringe_obj = cpl_frameset_new();
    cpl_mask  * stat_mask = NULL;
    cpl_image  * masterfringe_tmp = NULL;
    cpl_image  * masterfringe_error_tmp = NULL;
    hdrl_image  * masterfringe = NULL;

    for (cpl_size i = 0; i < fsize; i++) {
        cpl_frame * cur_frame=cpl_frameset_get_position(frameset, i);

        if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_RAW)) {
            cpl_frameset_insert(frameset_fringe,
                            cpl_frame_duplicate(cur_frame));
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_RAW_BPM)) {
            cpl_frameset_insert(frameset_fringe_bpm,
                            cpl_frame_duplicate(cur_frame));
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_OBJ_MASK)) {
            cpl_frameset_insert(frameset_fringe_obj,
                            cpl_frame_duplicate(cur_frame));
        } else if (stat_mask == NULL && !strcmp(cpl_frame_get_tag(cur_frame),
                        HDRLDEMO_STATIC_MASK)) {
            cpl_msg_info(cpl_func,"Static mask found in the SOF");
            stat_mask = cpl_mask_load(cpl_frame_get_filename(cur_frame), 0,
                            extnum_stat);
        } else if (masterfringe_tmp == NULL &&
                        !strcmp(cpl_frame_get_tag(cur_frame),
                        HDRLDEMO_MASTER_FRINGE)) {
            cpl_msg_info(cpl_func,"Master Fringe found in the SOF");
            masterfringe_tmp = cpl_image_load(cpl_frame_get_filename(cur_frame),
                            CPL_TYPE_DOUBLE, 0, 0);
        } else if (masterfringe_error_tmp == NULL &&
                        !strcmp(cpl_frame_get_tag(cur_frame),
                        HDRLDEMO_MASTER_FRINGE_ERROR)) {
            cpl_msg_info(cpl_func,"Master Fringe error found in the SOF");
            masterfringe_error_tmp = cpl_image_load(
                            cpl_frame_get_filename(cur_frame), CPL_TYPE_DOUBLE,
                            0, 0);
        }
    }

    if (masterfringe_tmp != NULL) {
        masterfringe = hdrl_image_create(masterfringe_tmp,
                        masterfringe_error_tmp);
    }
    cpl_image_delete(masterfringe_tmp);
    cpl_image_delete(masterfringe_error_tmp);

    /* verify that we have at least a fringe input frames */
    if (cpl_frameset_get_size(frameset_fringe) == 0) {
        cpl_frameset_delete(frameset_fringe);
        cpl_frameset_delete(frameset_fringe_bpm);
        cpl_frameset_delete(frameset_fringe_obj);
        cpl_mask_delete(stat_mask);
        hdrl_image_delete(masterfringe);
        hdrl_parameter_delete(collapse_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                        "Can not find Fringes in the SOF");
    }



    /* Required frame */
    frm_ima = cpl_frameset_find_const(frameset_fringe, HDRLDEMO_RAW);
    if (!frm_ima) {
        cpl_frameset_delete(frameset_fringe);
        cpl_frameset_delete(frameset_fringe_bpm);
        cpl_frameset_delete(frameset_fringe_obj);
        cpl_mask_delete(stat_mask);
        hdrl_image_delete(masterfringe);
        hdrl_parameter_delete(collapse_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                        "Cannot find any RAW frame");
    }



    /* read the fringe images and add bad pixel/object masks */
    hdrl_imagelist * ilist_fringe = hdrl_imagelist_new();
    hdrl_buffer * buf = hdrl_buffer_new();
    /* allow 2GiB of anonymous memory */
    hdrl_buffer_set_malloc_threshold(buf, 2048);

    if(cpl_frameset_get_size(frameset_fringe_bpm) > 0) {
        hdrldemo_hdrl_imagelist_load(frameset_fringe, extnum_raw, NULL, -1,
                        frameset_fringe_bpm, extnum_bpm, NULL, ron, gain,
                        ilist_fringe, buf);
    } else {
        hdrldemo_hdrl_imagelist_load(frameset_fringe, extnum_raw, NULL, -1,
                        NULL, -1, NULL, ron, gain,
                        ilist_fringe, buf);

    }
    /* TODO: AMO: why the following check ?
     * AGA: the buffer allocation could fail - not enough ram, ... */
    if (cpl_error_get_code()) {
        cpl_frameset_delete(frameset_fringe);
        cpl_frameset_delete(frameset_fringe_bpm);
        cpl_frameset_delete(frameset_fringe_obj);
        cpl_mask_delete(stat_mask);
        hdrl_image_delete(masterfringe);
        hdrl_imagelist_delete(ilist_fringe);
        hdrl_parameter_delete(collapse_params);
        hdrl_buffer_delete(buf);
        return cpl_error_get_code();
    }

    /* Read the object mask into an cpl_imaglist
     * TODO use cpl_masklist if implemented in the future */

    cpl_imagelist * ilist_obj = NULL;
    cpl_size fsize_fringe = cpl_frameset_get_size(frameset_fringe);

    if(cpl_frameset_get_size(frameset_fringe_obj) == fsize_fringe){
        ilist_obj = cpl_imagelist_new();
        cpl_msg_info(cpl_func,"Matching number of object masks found -reading");
        for (cpl_size i = 0; i < fsize_fringe; i++) {
            cpl_frame * cur_frame = cpl_frameset_get_position(frameset_fringe_obj, i);
            cpl_image * ima_obj = cpl_image_load(cpl_frame_get_filename(cur_frame),
                            CPL_TYPE_INT, 0, extnum_obj);
            cpl_imagelist_set(ilist_obj, ima_obj, i);
        }
    }


    hdrl_image * master = NULL;
    cpl_image * contrib_map = NULL;

    if(masterfringe == NULL) {
        /* Do the actual fringemap computation */

        hdrl_fringe_compute(ilist_fringe, ilist_obj, stat_mask, collapse_params,
                        &master, &contrib_map, NULL);
        /*Save the results*/
        /* If needed, replace the bad pixels with a value by including the next
         * line of code */
        //cpl_image_fill_rejected(hdrl_image_get_image(master), 0.);

        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_FRINGE,
        	"hdrldemo_fringe", "hdrldemo_masterfringe.fits",
			CPL_TYPE_FLOAT, hdrl_image_get_image(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_FRINGE_ERROR,
        	"hdrldemo_fringe", "hdrldemo_masterfringe_error.fits",
			CPL_TYPE_FLOAT, hdrl_image_get_error(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_FRINGE_CONTRIBUTION,
        	"hdrldemo_fringe", "hdrldemo_masterfringe_contribution.fits",
            CPL_TYPE_FLOAT, contrib_map, parlist, frameset);
    } else {

        hdrl_fringe_correct(ilist_fringe, ilist_obj, stat_mask, masterfringe, NULL);
        /* Save all fringecorrected images in single fits files */
        cpl_size numframes;
        numframes = hdrl_imagelist_get_size(ilist_fringe);

        for (cpl_size i = 0; i < numframes; i++){

        	char *outfile = cpl_sprintf("%s_%03d.fits", "hdrldemo_fringe_corrected",
                            (int)i+1);
            hdrldemo_save_image(NULL, NULL, HDRLDEMO_FRINGE_CORRECTED,
            	"hdrldemo_fringe", outfile, CPL_TYPE_FLOAT,
				hdrl_image_get_image_const(hdrl_imagelist_get_const(ilist_fringe, i)),
				parlist, frameset);
            cpl_free(outfile);
            outfile = cpl_sprintf("%s_%03d_error.fits",
                            "hdrNULL, ldemo_fringe_corrected", (int)i+1);
            hdrldemo_save_image(NULL, NULL, HDRLDEMO_FRINGE_CORRECTED_ERROR,
            	"hdrldemo_fringe", outfile, CPL_TYPE_FLOAT,
                hdrl_image_get_error_const(hdrl_imagelist_get_const(ilist_fringe, i)),
				parlist, frameset);
            cpl_free(outfile);
        }
    }
    /* cleanup */
    hdrl_parameter_delete(collapse_params);
    hdrl_imagelist_delete(ilist_fringe);
    cpl_imagelist_delete(ilist_obj);
    hdrl_image_delete(master);
    cpl_image_delete(contrib_map);
    if(stat_mask != NULL) cpl_mask_delete(stat_mask);
    if(masterfringe != NULL) hdrl_image_delete(masterfringe);
    cpl_frameset_delete(frameset_fringe);
    cpl_frameset_delete(frameset_fringe_bpm);
    cpl_frameset_delete(frameset_fringe_obj);


    hdrl_buffer_delete(buf);
    return (int)cpl_error_get_code();
}

/**@}*/

