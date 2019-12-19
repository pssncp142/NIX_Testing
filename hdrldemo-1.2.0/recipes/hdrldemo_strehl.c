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
  @defgroup hdrldemo_strehl Strehl computation recipe
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/
static cpl_error_code hdrldemo_strehl_save(
        const char              *   procatg,
        const char              *   filename,
        const cpl_propertylist  *   header,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   frameset);

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_strehl"

static char hdrldemo_strehl_description[] =
"                                                                           \n"
"The recipe computes the Strehl on individual bias subtracted images.       \n"
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
"  HDRLDEMO_STREHL              Fitsfile containing strehl as QC param      \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"                                                                           \n"
"There is currently one method implemented:                                 \n"
"                                                                           \n"
"The recipe assumes that the input raw image contains a single object       \n"
"and it is already pre-processed. This means that the instrument            \n"
"signature has been removed and the contribute from the sky background      \n"
"subtracted. Additionally, it allows the user to correct a residual         \n"
"local sky background evaluated in an annular region centred on the         \n"
"peak of the object PSF, by setting the values of the corresponding         \n"
"parameters (--bkg-radius-low, and --bkg-radius-high in arcsec units).      \n"
"                                                                           \n"
"The PSF is identified and its integrated flux, whose determination is      \n"
"controlled by the parameter --flux-radius (in arcsec), is normalized       \n"
"to unity. Next the PSF barycentre is computed and used to generate the     \n"
"theoretical normalised PSF. This depends on i) the telescope pupil         \n"
"characteristics (telescope radius and central obstruction parameters,      \n"
"controlled by --m1, and --m2 in meter) ii) the wavelength (parameter       \n"
"--wavelength in meter) at which the image has been obtained and iii)       \n"
"the detector pixel scale on sky in arcsec (--pixel-scale-x, and            \n"
"--pixel-scale-y).                                                          \n"
"                                                                           \n"
"Then the Strehl ratio is obtained by dividing the maximum intensity of     \n"
"the image PSF by the maximum intensity of the ideal PSF.                   \n"
"                                                                           \n"
"Please note that if no error image is provided, the error is estimated     \n"
"by measuring the MAD (scaled to an rms) on the input image, i.e. each      \n"
"pixel has the same associated error.                                       \n"
"                                                                          \n";
/* TODO:
 *  consistency checks on input parameters.
 *  Make sure psf size is at least 1x1 pixels
 *  Test without errors
 *  Add gain and ron input parameters for error model
 */
/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_strehl, HDRLDEMO_BINARY_VERSION, "HDRL Group", 
        PACKAGE_BUGREPORT, "2014", "HDRLDEMO - Strehl",
        hdrldemo_strehl_description);                          

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_strehl_fill_parameterlist(
        cpl_parameterlist   *   self) 
{                                  
    cpl_parameter   *   par ;

    /* --hdrldemo_strehl.ext-nb-raw */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw", CPL_TYPE_INT,
            "FITS extension of the RAW", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-r");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    
    /* --hdrldemo_strehl.ext-nb-raw-err */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-err", CPL_TYPE_INT,
            "FITS extension of the ERROR", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-e");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_strehl.ext-nb-raw-bpm */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-bpm", CPL_TYPE_INT,
            "FITS extension or the input BPM", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-b");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_strehl.region-llx/lly/urx/ury */
    hdrl_parameter * deflts = hdrl_rect_region_parameter_create(1, 1, 0, 0) ;
    cpl_parameterlist * reg_param = hdrl_rect_region_parameter_create_parlist(
                RECIPE_NAME, "", "region-", deflts);
    hdrl_parameter_delete(deflts) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(reg_param) ;
            p != NULL; p = cpl_parameterlist_get_next(reg_param))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(reg_param);

    hdrl_parameter * s_def =
        hdrl_strehl_parameter_create(1.635e-6,
                                     5.08/2, 1.8288/2,
                                     0.0165966, 0.0165966,
                                     1.5, -1 , -1);

    cpl_parameterlist * s_param =
        hdrl_strehl_parameter_create_parlist(RECIPE_NAME, "", s_def);
    hdrl_parameter_delete(s_def) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(s_param);
            p != NULL; p = cpl_parameterlist_get_next(s_param))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(s_param);

    return CPL_ERROR_NONE;
}

static int hdrldemo_strehl(
        cpl_frameset            *   frameset, 
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter     *   par = NULL;
    int                         extnum_raw = 0;
    int                         extnum_err = 0;
    int                         extnum_bpm = 0;
    hdrl_strehl_result          res;

    cpl_frame               *   in_frm ;
    cpl_frame               *   err_frm ;
    cpl_frame               *   bpm_frm ;
    hdrl_image              *   hima ;
    hdrl_parameter          *   region_params = NULL;

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

    /* Parse the Strehl Parameters */
    hdrl_parameter * p = hdrl_strehl_parameter_parse_parlist(parlist,
                                                             RECIPE_NAME);
    if (p == NULL) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                "Parsing of the strehl parameters failed");
    }

    /* Parse the Region Parameters */
    region_params=hdrl_rect_region_parameter_parse_parlist(parlist,
            RECIPE_NAME, "region-") ;
    if (region_params == NULL) {
        hdrl_parameter_delete(p);
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                "Parsing of the region parameters failed");
    }

    /* Load INPUT Data */
    /* Get the first Required frame */
    if ((in_frm = cpl_frameset_find(frameset, HDRLDEMO_RAW))==NULL) {
        if (region_params) hdrl_parameter_delete(region_params) ;
        if (p) hdrl_parameter_delete(p) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                "Missing RAW file");
    }
    bpm_frm = cpl_frameset_find(frameset, HDRLDEMO_RAW_BPM) ;
    err_frm = cpl_frameset_find(frameset, HDRLDEMO_RAW_ERROR) ;


    /* Load the image */
    if (hdrldemo_hdrl_image_load(in_frm, extnum_raw, err_frm, extnum_err,
                bpm_frm, extnum_bpm, region_params, 1, 1,
                &hima) != CPL_ERROR_NONE) {
        if (region_params) hdrl_parameter_delete(region_params) ;
        if (p) hdrl_parameter_delete(p) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                "Cannot load RAW image");
    }
    if (region_params) hdrl_parameter_delete(region_params) ;


    if(err_frm == NULL) {
        /* Replace the error image created by the ron/gain by the mad-scaled
         * rms of the image */
       double dmad = 0.;
       cpl_image * img = cpl_image_load(cpl_frame_get_filename(in_frm),
    		   	   	   	   	   	   	    CPL_TYPE_DOUBLE, 0, extnum_raw);
       cpl_image_get_mad(img, &dmad);
       cpl_image_delete(img);

       cpl_image_multiply_scalar(hdrl_image_get_error(hima), 0.);
       cpl_image_add_scalar(hdrl_image_get_error(hima),
                            (dmad * CPL_MATH_STD_MAD));
    }

    /* Strehl COMPUTATION */
    res = hdrl_strehl_compute(hima, p);
    hdrl_parameter_delete(p);
    cpl_msg_info(cpl_func,"strehl=%g+-%g", res.strehl_value.data,
                 res.strehl_value.error);
    /* expected difference sqrt(pi / 2)  due to median */
    cpl_msg_info(cpl_func, "star peak at %g/%g: %g +- %g", res.star_x,
                 res.star_y, res.star_peak.data, res.star_peak.error);
    cpl_msg_info(cpl_func, "star flux %g+-%g", res.star_flux.data,
                 res.star_flux.error);
    cpl_msg_info(cpl_func,"median estimated background=%g+-%g "
                 "(computed error %g)",
                 res.star_background.data, res.star_background.error,
                 res.computed_background_error);
    cpl_propertylist* header=cpl_propertylist_new();
    cpl_propertylist_update_double(header, "ESO QC STREHL",
                                   res.strehl_value.data);
    cpl_propertylist_update_double(header, "ESO QC STREHL ERROR",
                                   res.strehl_value.error);
    cpl_propertylist_update_string(header, CPL_DFS_PRO_CATG, "HDRLDEMO_STREHL");

    hdrldemo_strehl_save("HDRLDEMO_STREHL","hdrldemo_strehl.fits", header,
                         parlist,frameset);


    /* Cleanup */
    cpl_propertylist_delete(header);
    hdrl_image_delete(hima);
    /* In case of Poisson error model: */

    return (int)cpl_error_get_code();
}

static cpl_error_code hdrldemo_strehl_save(
        const char              *   procatg,
        const char              *   filename,
        const cpl_propertylist  *   header,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   frameset)
{
    /* Add a QC parameter  */
    cpl_propertylist * applist = cpl_propertylist_new();

    /* Add the product category and save image */
    cpl_propertylist_update_string(applist, CPL_DFS_PRO_CATG, procatg);

    cpl_dfs_save_propertylist(frameset, NULL, parlist, frameset, NULL,
            "hdrldemo_strehl", header, NULL, PACKAGE "/" PACKAGE_VERSION,
            filename);
    cpl_propertylist_delete(applist);
    return cpl_error_get_code();
}


/**@}*/

