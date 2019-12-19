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
 * @defgroup hdrldemo_flat_nir  master flat creation in the NIR regime
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

#define RECIPE_NAME "hdrldemo_flat_nir"

static char hdrldemo_flat_nir_description[] =
"                                                                           \n"
"The recipe derives a master flatfield in the near infrared regime.         \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:                      Explanation:           Required:       \n"
"  FLAT or FLAT_ON                   Data                   Yes             \n"
"  FLAT_OFF                          Data                   NO              \n"
"  STATIC_MASK                       Static mask            NO              \n"
"  MASTER_DARK                       Master Dark            NO              \n"
"  MASTER_DARK_ERROR                 Master Dark error      NO              \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                      Explanation:                           \n"
"  MASTER_FLAT                       Master flatfield                       \n"
"  MASTER_FLAT_ERROR                 Error of the master flatfield          \n"
"  MASTER_FLAT_CONTRIBUTION          Pixel contribution map                 \n"
"                                                                           \n"
"                                                                           \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"There are currently two methods implemented:                               \n"
"                                                                           \n"
"                                                                           \n"
"(--flat.method=low): Controlled by flat.filter-size-x/-y;                  \n"
"                                                                           \n"
"           The algorithm multiplicatively normalizes the input images      \n"
"           by the median of the image to unity. A static mask              \n"
"           (STATIC_MASK) can be provided to the algorithm in order to      \n"
"           define the pixels that should be taken into account when        \n"
"           computing the normalisation factor. This allows the user to     \n"
"           normalize the flatfield e.g.  only by the illuminated           \n"
"           section.  Then all normalized images are collapsed into a       \n"
"           single master flatfield. The collapsing can be done with        \n"
"           all methods currently implemented in hdrl. Finally, the         \n"
"           master flatfield is smoothed by a median filter. The size       \n"
"           of the median-filter is controlled by the above mentioned       \n"
"           recipe parameter flat.filter-size-x/-y. The associated          \n"
"           error of the final masterframe is the error derived via         \n"
"           error propagation of the previous steps, i.e. the smoothing     \n"
"           itself is considered noiseless.                                 \n"
"                                                                           \n"
"                                                                           \n"
"(--flat.method=high): Controlled by flat.filter-size-x/-y;                 \n"
"                                                                           \n"
"           The algorithm first smoothes the input images by a median       \n"
"           filter and divides each input image through the smoothed        \n"
"           image. The size of the medianfilter is controlled by the        \n"
"           above mentioned recipe parameter flat.filter-size-x/-y. The     \n"
"           smoothed images is considered to be noiseless i.e. the          \n"
"           relative error of the resulting images is the same as the       \n"
"           one of the input image. Then all residual images are            \n"
"           collapsed into a single master flatfield. The collapsing        \n"
"           can be done with all methods currently implemented in           \n"
"           hdrl. Moreover, it is also possible to give a static mask       \n"
"           (STATIC_MASK) to the algorithm which e.g. allows the user       \n"
"           to distinguish illuminated and not illuminated regions. In      \n"
"           this case the smoothing procedure is done twice, once for       \n"
"           the illuminated region and once for the blanked                 \n"
"           region. This ensures that the information of one region         \n"
"           does not influence the other regions during the smoothing       \n"
"           process.                                                        \n"
"                                                                           \n"
"Please note that you can have FLAT or FLAT_ON as input frames in the       \n"
"SOF. They behave identical, but if additionaly to the FLAT_ON you also     \n"
"have FLAT_OFF frames in the SOF the FLAT_OFF will be subtrated from        \n"
"the FLAT_ON frames (with error propagation). Moreover, if you also         \n"
"have a MASTER_DARK frame as input this frame will also be subtracted,      \n"
"i.e. (FLAT - MASTER_DARK) or ((FLAT_ON - FLAT_OFF) - MASTER_DARK).         \n"
"Moreover, the errors of the input frames (FLAT, FLAT_ON, FLAT_OFF) are     \n"
"estimated with a shot noise model by using the Gain (--gain) and Ron       \n"
"(--ron).                                                                   \n"
"                                                                           \n"
"Please note that part of the code is paralelised. In order to optimise     \n"
"use of computing resources you should set the environment variable         \n"
"OMP_NUM_THREADS to a proper value, like (for bash), to use 4 cores         \n"
"export OMP_NUM_THREADS=4                                                   \n"
"                                                                           \n";

/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_flat_nir,
                    HDRLDEMO_BINARY_VERSION,
                    "HDRL Group",
                    PACKAGE_BUGREPORT,
                    "2013",
                    "HDRLDEMO - Master Flat creation in the NIR regime",
                    hdrldemo_flat_nir_description);

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_flat_nir_fill_parameterlist(
        cpl_parameterlist   *   self) {
    cpl_parameter   *   par ;

    /* --hdrldemo.flat.ExtensionNumber */
    par = cpl_parameter_new_value(RECIPE_NAME".extension-number", CPL_TYPE_INT,
                                  "FITS extension to load", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo.flat.Gain */
    par = cpl_parameter_new_value(RECIPE_NAME".gain", CPL_TYPE_DOUBLE,
                                  "Gain in [e- / ADU]", RECIPE_NAME, 2.5);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo.flat.ccd-ron */
    par = cpl_parameter_new_value(RECIPE_NAME".ccd-ron", CPL_TYPE_DOUBLE,
                                  "Readout noise in ADU.", RECIPE_NAME, 10.);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ccd-ron");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* Master-flat related parameters */

    hdrl_parameter * sigclip_def =
        hdrl_collapse_sigclip_parameter_create(3., 3., 5);
    hdrl_parameter * minmax_def =
        hdrl_collapse_minmax_parameter_create(1., 1.);

    cpl_parameterlist * pflatcollapse = hdrl_collapse_parameter_create_parlist(
                RECIPE_NAME, "collapse", "MEDIAN", sigclip_def, minmax_def) ;
    hdrl_parameter_delete(sigclip_def); 
    hdrl_parameter_delete(minmax_def);
    for (cpl_parameter * p = cpl_parameterlist_get_first(pflatcollapse) ;
            p != NULL; p = cpl_parameterlist_get_next(pflatcollapse))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(pflatcollapse);

    /* Create flat parameters */
    hdrl_parameter * deflts = hdrl_flat_parameter_create(5, 5,
            HDRL_FLAT_FREQ_LOW);
    cpl_parameterlist * flat_param = hdrl_flat_parameter_create_parlist(
                RECIPE_NAME, "flat", deflts) ;
    hdrl_parameter_delete(deflts) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(flat_param) ;
            p != NULL; p = cpl_parameterlist_get_next(flat_param))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(flat_param);


    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   perform master flat combination
  @param   frameset   input set of frames
  @param   parlist    input recipe parameters
  @return  cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_flat_nir(
        cpl_frameset            *   frameset, 
        const cpl_parameterlist *   parlist)
{
    const cpl_frame     *   frm_ima;
    cpl_frameset        *   frameset_final = NULL;
    const cpl_parameter *   par;
    hdrl_parameter      *   collapse_params;
    hdrl_parameter      *   flat_params;

    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /* Get parameters*/
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".extension-number");
    int ext_num = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".gain");
    double gain = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ccd-ron");
    double ron = cpl_parameter_get_double(par);

    /* Collapse parameters */
    collapse_params = hdrl_collapse_parameter_parse_parlist(parlist, RECIPE_NAME".collapse");

    /* Parse the FLAT Parameters */
    flat_params=hdrl_flat_parameter_parse_parlist(parlist, RECIPE_NAME".flat") ;
    if (flat_params == NULL) {
        if (collapse_params) hdrl_parameter_delete(collapse_params) ;

        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                "Parsing of the recipe parameters failed");
    }

    cpl_size fsize = cpl_frameset_get_size(frameset);
    cpl_frameset * frameset_flat = cpl_frameset_new();
    cpl_frameset * frameset_flat_off = cpl_frameset_new();
    cpl_image * master_dark = NULL;
    cpl_image * master_dark_err = NULL;
    cpl_mask  * stat_mask = NULL;

    for (cpl_size i = 0; i < fsize; i++) {
        cpl_frame * cur_frame=cpl_frameset_get_position(frameset, i);

        if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_FLAT) ||
                        !strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_FLAT_ON)) {
            cpl_frameset_insert(frameset_flat,
                            cpl_frame_duplicate(cur_frame));
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_FLAT_OFF)) {
            cpl_frameset_insert(frameset_flat_off,
                            cpl_frame_duplicate(cur_frame));
        } else if (master_dark == NULL && !strcmp(cpl_frame_get_tag(cur_frame),
                        HDRLDEMO_MASTER_DARK)) {
            cpl_msg_info(cpl_func,"Masterdark found in the SOF");
            master_dark = cpl_image_load(cpl_frame_get_filename(cur_frame),
                            CPL_TYPE_DOUBLE, 0, 0);
        } else if (master_dark_err == NULL && !strcmp(cpl_frame_get_tag(cur_frame),
                        HDRLDEMO_MASTER_DARK_ERROR)) {
            cpl_msg_info(cpl_func,"Masterdark error found in the SOF");
            master_dark_err = cpl_image_load(cpl_frame_get_filename(cur_frame),
                            CPL_TYPE_DOUBLE, 0, 0);
        } else if (stat_mask == NULL && !strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_STATIC_MASK)) {
            cpl_msg_info(cpl_func,"Static mask found in the SOF");
            stat_mask = cpl_mask_load(cpl_frame_get_filename(cur_frame), 0, 0);
        }
    }

    /* verify that we have at least a flat input frames */
    if (cpl_frameset_get_size(frameset_flat) == 0) {
        cpl_frameset_delete(frameset_flat);
        cpl_frameset_delete(frameset_flat_off);
        cpl_image_delete(master_dark);
        cpl_image_delete(master_dark_err);
        cpl_mask_delete(stat_mask);
        hdrl_parameter_delete(collapse_params) ;
        hdrl_parameter_delete(flat_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
               "Can not find Flats in the SOF");
    }

    hdrl_image * hdrl_master_dark = NULL;
    if( master_dark != NULL ) {
        hdrl_master_dark = hdrl_image_create(master_dark, master_dark_err);
        cpl_image_delete(master_dark);
        cpl_image_delete(master_dark_err);
    }

    /* TODO leftover from masterdark - eliminate */
    frameset_final = cpl_frameset_duplicate(frameset_flat);
    cpl_frameset_delete(frameset_flat);


    /* Required frame */
    frm_ima = cpl_frameset_find_const(frameset_final, HDRLDEMO_FLAT);
    if (!frm_ima) {
        frm_ima = cpl_frameset_find_const(frameset_final, HDRLDEMO_FLAT_ON);
    }
    if (!frm_ima) {
        cpl_frameset_delete(frameset_final);
        cpl_frameset_delete(frameset_flat_off);
        cpl_mask_delete(stat_mask);
        hdrl_parameter_delete(collapse_params) ;
        hdrl_parameter_delete(flat_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                        "Cannot find any RAW frame");
    }

    cpl_propertylist * plist =
        cpl_propertylist_load(cpl_frame_get_filename(frm_ima), ext_num);
    cpl_size nx = hdrldemo_get_naxis1(plist);
    cpl_size ny = hdrldemo_get_naxis2(plist);
    cpl_propertylist_delete(plist);

    if (nx <= 0 || ny <= 0) {
        cpl_frameset_delete(frameset_final);
        cpl_frameset_delete(frameset_flat_off);
        cpl_mask_delete(stat_mask);
        hdrl_parameter_delete(collapse_params) ;
        hdrl_parameter_delete(flat_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Could not retrieve image dimensions");
    }


    /* apply overscan correction and dark normalization from images the input
     * provides, allocate output and write scaled images into them */
    hdrl_imagelist * corrected = hdrl_imagelist_new();
    hdrl_buffer * buf = hdrl_buffer_new();
    /* allow 1GiB of anonymous memory */
    hdrl_buffer_set_malloc_threshold(buf, 1024);

    hdrldemo_hdrl_imagelist_load(frameset_final, ext_num, NULL, -1, NULL,
                    -1, NULL, ron, gain,
                    corrected, buf);


    if (cpl_error_get_code()) {
        cpl_frameset_delete(frameset_final);
        cpl_frameset_delete(frameset_flat_off);
        cpl_mask_delete(stat_mask);
        hdrl_imagelist_delete(corrected);
        hdrl_parameter_delete(collapse_params);
        hdrl_parameter_delete(flat_params);
        hdrl_image_delete(hdrl_master_dark);
        hdrl_buffer_delete(buf);
        return cpl_error_get_code();
    }

    /* Only if we have matching OFF frames we do the on-off subtraction */
    if(cpl_frameset_get_size(frameset_final) ==
                    cpl_frameset_get_size(frameset_flat_off)){
        cpl_msg_info(cpl_func,"Reading FLAT_OFF frames and subtracting ... ");
        hdrl_imagelist * corrected_off = hdrl_imagelist_new();

        hdrldemo_hdrl_imagelist_load(frameset_flat_off, ext_num, NULL, -1, NULL,
                                     -1, NULL, ron, gain,
                                     corrected_off, buf);

        hdrl_imagelist_sub_imagelist(corrected, corrected_off);
        hdrl_imagelist_delete(corrected_off);
    }


    /* Subtract the master dark frame */

    if(hdrl_master_dark != NULL){
        hdrl_imagelist_sub_image(corrected, hdrl_master_dark);
        hdrl_image_delete(hdrl_master_dark);
    }

    hdrl_image * master = NULL;
    cpl_image * contrib_map = NULL;

    /* Do the actual flatfield computation */
    hdrl_flat_compute(corrected, stat_mask, collapse_params,
                                     flat_params, &master, &contrib_map);


    hdrl_parameter_delete(collapse_params);
    hdrl_parameter_delete(flat_params);
    hdrl_imagelist_delete(corrected);


    /*Save the results*/
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_FLAT,
        	"hdrldemo_flat_nir", "hdrldemo_masterflat.fits",
			CPL_TYPE_FLOAT, hdrl_image_get_image(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_FLAT_ERROR,
        	"hdrldemo_flat_nir", "hdrldemo_masterflat_error.fits",
			CPL_TYPE_FLOAT, hdrl_image_get_error(master), parlist, frameset);
        hdrldemo_save_image(NULL, NULL, HDRLDEMO_MASTER_FLAT_CONTRIBUTION,
        	"hdrldemo_flat_nir", "hdrldemo_masterflat_contribution.fits",
            CPL_TYPE_FLOAT, contrib_map, parlist, frameset);

    /* cleanup */
    hdrl_image_delete(master);
    cpl_image_delete(contrib_map);

    if(stat_mask != NULL) cpl_mask_delete(stat_mask);

    /* hdrldemo_cube_it_delete(input_it); */
    cpl_frameset_delete(frameset_final);
    cpl_frameset_delete(frameset_flat_off);


    hdrl_buffer_delete(buf);
    return (int)cpl_error_get_code();
}

/**@}*/

