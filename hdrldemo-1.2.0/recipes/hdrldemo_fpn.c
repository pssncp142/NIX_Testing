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

/*----------------------------------------------------------------------------*/
/**
 *                              Includes
 */
/*----------------------------------------------------------------------------*/

#include "hdrl.h"

#include "hdrldemo_utils.h"
#include "hdrldemo_dfs.h"

#include <cpl.h>

/*----------------------------------------------------------------------------*/
/**
 *                              Defines
 */
/*----------------------------------------------------------------------------*/

#define HDRLDEMO_FPN        "HDRLDEMO_FPN"
#define HDRLDEMO_FPN_MASK   "HDRLDEMO_FPN_MASK"

/*----------------------------------------------------------------------------*/
/**
 *                 Typedefs: Structs and enum types
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 *                              Functions prototypes
 */
/*----------------------------------------------------------------------------*/

/* Save to disk a output for one fixed pattern noise execution over a file */
static cpl_error_code hdrldemo_fpn_save(
                cpl_frameset                   *frameset,
                const cpl_frame                *use_frame,
                const cpl_parameterlist        *parlist,
                const cpl_size                 num_raw,
                cpl_image                      *out_img,
                const double                   std,
                const double                   stdmad,
                const char                     *filename);

/*----------------------------------------------------------------------------*/
/**
 *                          Static variables
 */
/*----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_fpn"

static char hdrldemo_fpn_description[] =
"                                                                           \n"
"The recipe detects fixed pattern noise on an image                         \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"   DO category:                 Explanation:                Required:      \n"
"   RAW                          Data                        Yes            \n"
"   POWERSPEC_MASK               Mask for the powerspectrum  No             \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                 Explanation:                                \n"
"  HDRLDEMO_FPN                 Computed power spectrum                     \n"
"  HDRLDEMO_FPN_MASK            Mask used for statistical computations      \n"
"                                                                           \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"                                                                           \n"
"The algorithm first computes the power spectrum of the input image         \n"
"(RAW) using the Fast Fourier Transform (FFT) and then it computes the      \n"
"standard deviation (std) and the mad-based std of the power_spectrum       \n"
"excluding the masked region. For this, the user can provide an             \n"
"optional mask (POWERSPEC_MASK) or use the dc_mask_x and dc_mask_y recipe   \n"
"parameter to create one on the fly. The mask created on the fly will       \n"
"start at pixel (1,1) and extend in both direction up to (dc_mask_x,        \n"
"dc_mask_y).                                                                \n"
"                                                                           \n"
"Please note that the power spectrum contains the DC component (the DC      \n"
"term is the 0 Hz term and is equivalent to the average of all the          \n"
"samples in the window) in pixel (1,1). Moreover, the masked created on     \n"
"the fly and the optional mask are combined and are both taken into         \n"
"account.                                                                   \n";

/* Standard CPL recipe definition */
cpl_recipe_define(  hdrldemo_fpn,
                    HDRLDEMO_BINARY_VERSION,
                    "HDRL Group",
                    PACKAGE_BUGREPORT,
                    "2018",
                    "HDRLDEMO - Fixed pattern noise",
                    hdrldemo_fpn_description);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrldemo_fpn     Fixed pattern noise recipe
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
 *                              Functions code
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief    Interpret the command line options and execute the data processing
 *
 * @param    frameset   the frames list
 * @param    parlist    the parameters list
 *
 * @return   CPL_ERROR_NONE if everything is OK or CPL_ERROR_CODE in other case
 *
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_fpn(
                cpl_frameset            *frameset,
                const cpl_parameterlist *parlist)
{
    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
        return cpl_error_get_code();
    }

    /* Extract and verify data in the parameters of the recipe */
    const cpl_parameter *par;

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-raw");
    int ext_r = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ext-nb-mask");
    int ext_mask = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".dc_mask_x");
    int dc_mask_x = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".dc_mask_y");
    int dc_mask_y = cpl_parameter_get_int(par);

    /* Check the input parameters */
    cpl_error_ensure(ext_r >= 0 && dc_mask_x >= 1 && dc_mask_y >= 1,
                     CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_ILLEGAL_INPUT,
                                     "Error in the input parameters "
                                     "(ext_r, dc_mask_x, dc_mask_y)");


 /* Create the raw and powerspec_mask framesets (only RAW and POWERSPEC_MASK) */
    cpl_size     nframes  = cpl_frameset_get_size(frameset);
    cpl_frameset *fs_raws = cpl_frameset_new();
    cpl_frameset *fs_masks = cpl_frameset_new();
    for (cpl_size i = 0; i < nframes; i++) {

        cpl_frame  *frame    = cpl_frameset_get_position(frameset, i);
        const char *filename = cpl_frame_get_filename(frame);

        if (!strcmp(cpl_frame_get_tag(frame), HDRLDEMO_RAW)) {

            /* Check that the RAW extension exist and insert */
            if (cpl_fits_count_extensions(filename) < ext_r){
                cpl_frameset_delete(fs_raws);
                return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                "RAW extension doesn't exist in file '%s'",
                                filename);
            }
            cpl_frameset_insert(fs_raws, cpl_frame_duplicate(frame));

        } else if (!strcmp(cpl_frame_get_tag(frame), HDRLDEMO_POWERSPEC_MASK)) {

            /* Check that the MASK extension exist and insert */
            if (cpl_fits_count_extensions(filename) < ext_mask){
                cpl_frameset_delete(fs_masks);
                return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                "POWERSPEC_MASK extension doesn't exist in file"
                                " '%s'", filename);
            }
            cpl_frameset_insert(fs_masks, cpl_frame_duplicate(frame));
        }
    }


    /* Loop to compute the pattern noise over the raws input files */
    cpl_size n_raws = cpl_frameset_count_tags(frameset, HDRLDEMO_RAW);
    cpl_size n_masks = cpl_frameset_count_tags(frameset, HDRLDEMO_POWERSPEC_MASK);
    for (cpl_size n_frame = 0; n_frame < n_raws; n_frame++) {

        cpl_frame  *frm_raw      = cpl_frameset_get_position(fs_raws, n_frame);
        const char *filename_raw = cpl_frame_get_filename(frm_raw);

        /* Load the input image */
        cpl_msg_info(cpl_func,"Load image, filename=%s ...", filename_raw);
        cpl_image *img_in = cpl_image_load(cpl_frame_get_filename(frm_raw),
                        CPL_TYPE_DOUBLE, 0, ext_r);
        if (!img_in || cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_frameset_delete(fs_raws);
            if (fs_masks) cpl_frameset_delete(fs_masks);
            if (img_in) cpl_image_delete(img_in);
            return cpl_error_get_code();
        }

        /* Load MASK's */
        cpl_mask *mask_in = NULL;
        if (n_frame < n_masks && n_masks > 0) {

            cpl_frame  *frm_mask = cpl_frameset_get_position(fs_masks, n_frame);
            const char *filename_mask = cpl_frame_get_filename(frm_mask);

            /* Load the input image */
            cpl_msg_info(cpl_func,"Load mask, filename=%s ...", filename_mask);
            mask_in = cpl_mask_load(filename_mask, 0, ext_mask);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_frameset_delete(fs_raws);
                cpl_frameset_delete(fs_masks);
                if (img_in ) cpl_image_delete(img_in);
                if (mask_in) cpl_mask_delete(mask_in);
                return cpl_error_get_code();
            }
        }

        /* Compute Pattern noise */
        cpl_image *power_spectrum  = NULL;
        double std = -1.;
        double std_mad = -1.;
        if (hdrl_fpn_compute(img_in, mask_in, dc_mask_x, dc_mask_y,
                        &power_spectrum, &std, &std_mad) != CPL_ERROR_NONE) {

            cpl_frameset_delete(fs_raws);
            if (fs_masks ) cpl_frameset_delete(fs_masks);
            if (img_in     ) cpl_image_delete(img_in);
            if (mask_in    ) cpl_mask_delete(mask_in);
            if (power_spectrum ) cpl_image_delete(power_spectrum);
            return cpl_error_get_code();
        }

        /* Save the pattern noise results */
        if (hdrldemo_fpn_save( frameset,
                        frm_raw,
                        parlist,
                        n_frame + 1,
                        power_spectrum,
                        std,
                        std_mad,
                        filename_raw) != CPL_ERROR_NONE) {
            cpl_image_delete(img_in);
            if (mask_in) cpl_mask_delete(mask_in);
            cpl_image_delete(power_spectrum);
            cpl_frameset_delete(fs_raws);
            if (fs_masks) cpl_frameset_delete(fs_masks);
            return cpl_error_get_code();
        }

        /* Cleanup */
        cpl_image_delete(img_in);
        if (mask_in) cpl_mask_delete(mask_in);
        cpl_image_delete(power_spectrum);
    }


    /* Cleanup */
    cpl_frameset_delete(fs_raws);
    if (fs_masks) cpl_frameset_delete(fs_masks);


    return CPL_ERROR_NONE;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
 * @brief Function needed by cpl_recipe_define to fill the input parameters
 *
 * @param  parlist   parameterlist where you need put parameters
 *
 * @return cpl_error_code
 *
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_fpn_fill_parameterlist(
                cpl_parameterlist *parlist)
{
    /* Add the different default parameters to the recipe */
    cpl_errorstate prestate = cpl_errorstate_get();

    /* Fill the parameters list */
    cpl_error_code e;
    cpl_boolean    range     = CPL_TRUE;
    const void     *dummyMin = NULL;
    const void     *dummyMax = NULL;


    /* --hdrldemo_fpn.ext-nb-raw */
    int ext_r = 0;
    e = hdrldemo_fill_parameter( RECIPE_NAME, parlist,
                    RECIPE_NAME".ext-nb-raw",
                    "ext-r",
                    !range, dummyMin, dummyMax, CPL_TYPE_INT,
                    &ext_r,
                    "FITS extension of the RAW.");
    if (e != CPL_ERROR_NONE) return (int)e;

    /* --hdrldemo_fpn.ext-nb-mask */
    int ext_mask = 0;
    e = hdrldemo_fill_parameter( RECIPE_NAME, parlist,
                    RECIPE_NAME".ext-nb-mask",
                    "ext-mask",
                    !range, dummyMin, dummyMax, CPL_TYPE_INT,
                    &ext_mask,
                    "FITS extension of the POWERSPEC_MASK.");
    if (e != CPL_ERROR_NONE) return (int)e;

    /* --hdrldemo_fpn.dc_mask_x */
    int dc_mask_x = 1;
    e = hdrldemo_fill_parameter( RECIPE_NAME, parlist,
                    RECIPE_NAME".dc_mask_x",
                    "dc_mask_x",
                    !range, dummyMin, dummyMax, CPL_TYPE_INT,
                    &dc_mask_x,
                    "x-size (pixel) of the mask starting at (x,y) = (1,1)");
    if (e != CPL_ERROR_NONE) return (int)e;

    /* --hdrldemo_fpn.dc_mask_y */
    int dc_mask_y = 1;
    e = hdrldemo_fill_parameter( RECIPE_NAME, parlist,
                    RECIPE_NAME".dc_mask_y",
                    "dc_mask_y",
                    !range, dummyMin, dummyMax, CPL_TYPE_INT,
                    &dc_mask_y,
                    "y-size (pixel) of the mask starting at (x,y) = (1,1)");
    if (e != CPL_ERROR_NONE) return (int)e;


    /* Check possible errors */
    if (!cpl_errorstate_is_equal(prestate)) {
        return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                        "hdrldemo_paternnoise_fill_parameterlist failed!");
    }

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief    Save to disk a new *.fits output file
 *
 * @param  frameset     Input frameset in the recipe.
 * @param  use_frame    Input frame    with the RAW data
 * @param  parlist      Input parameter list in the recipe.
 * @param  num_raw      Number of input raw treatment
 * @param  out_img      cpl_image result of the pattern noise calculation
 * @param  std          STD of the power spectrum
 * @param  stdmad          STDMAD of the power spectrum
 * @param  filename     Name of the output *.fits file
 *
 * @return   cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_fpn_save(
                cpl_frameset                   *frameset,
                const cpl_frame                *use_frame,
                const cpl_parameterlist        *parlist,
                const cpl_size                 num_raw,
                cpl_image                      *out_img,
                const double                   std,
                const double                   stdmad,
                const char                     *filename)
{
    /* Check inputs */
    cpl_error_ensure(frameset && parlist && out_img && filename,
                    CPL_ERROR_NULL_INPUT, return CPL_ERROR_NULL_INPUT,
                                    "Null inputs in save function");

    cpl_errorstate preState = cpl_errorstate_get();

    /* Create output filename */
    char *out_filename      = cpl_sprintf("%s_%03lld.fits",
                    "hdrldemo_fpn_power_spectrum",      num_raw);
    char *out_filename_mask = cpl_sprintf("%s_%03lld.fits",
                    "hdrldemo_fpn_power_spectrum_mask", num_raw);

    /* Add a QC parameter list */
    cpl_propertylist *qclist = cpl_propertylist_new();

    /* Create an output frameset */
    cpl_frameset *usedframes = cpl_frameset_new();
    cpl_frameset_insert(usedframes, cpl_frame_duplicate(use_frame));


    /* Save pattern noise mask image */
    cpl_propertylist_update_string(qclist, CPL_DFS_PRO_CATG, HDRLDEMO_FPN_MASK);
    cpl_image *out_mask_image = cpl_image_new_from_mask(cpl_image_get_bpm(out_img));
    cpl_dfs_save_image(frameset, NULL, parlist, usedframes, use_frame,
                       out_mask_image, CPL_TYPE_INT, RECIPE_NAME, qclist, NULL,
                       PACKAGE "/" PACKAGE_VERSION, out_filename_mask);
    cpl_image_delete(out_mask_image);


    /* Add STD */
    if (!isnan(std)) {
        cpl_propertylist_update_double(qclist, "ESO QC FPN STD", std);
    } else {
        cpl_propertylist_update_double(qclist, "ESO QC FPN STD", -DBL_MAX);
    }

    /* Add Peak value */
    if (!isnan(stdmad)) {
        cpl_propertylist_update_double(qclist, "ESO QC FPN STDMAD", stdmad);
    } else {
        cpl_propertylist_update_double(qclist, "ESO QC FPN STDMAD", -DBL_MAX);
    }

    /* Save pattern noise data image */
    cpl_propertylist_update_string(qclist, CPL_DFS_PRO_CATG, HDRLDEMO_FPN);
    cpl_dfs_save_image(frameset, NULL, parlist, usedframes, use_frame, out_img,
                       CPL_TYPE_DOUBLE, RECIPE_NAME, qclist, NULL,
                       PACKAGE "/" PACKAGE_VERSION, out_filename);


    /* Cleanup */
    cpl_frameset_delete(usedframes);
    cpl_free(out_filename);
    cpl_free(out_filename_mask);
    cpl_propertylist_delete(qclist);

    /* Check possible errors */
    if (!cpl_errorstate_is_equal(preState)) {
        return cpl_error_get_code();
    }

    return CPL_ERROR_NONE;
}
