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

#include "hdrldemo_utils.h"
#include "hdrldemo_dfs.h"
#include <cpl.h>

/*----------------------------------------------------------------------------*/
/**
  @defgroup hdrldemo_bpm_3d     Bad Pixel Maps 3D Recipe
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_error_code hdrldemo_bpm_3d_save_imagelist(const char *,
        const char * filename, const cpl_type, const cpl_imagelist *, 
        const cpl_parameterlist *, cpl_frameset *);

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_bpm_3d"

static char hdrldemo_bpm_3d_description[] =
"                                                                           \n"
"The recipe derives bad pixels on a stack of (identical) images             \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:                      Explanation:           Required:       \n"
"  RAW                               Data                   Yes             \n"
"  RAW_BPM                           Bad Pixel Mask         No              \n"
"  RAW_ERROR                         Associated Error       No              \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                      Explanation:                           \n"
"  HDRLDEMO_MASTER_BPM_LIST          Master bad pixel masks                 \n"
"  HDRLDEMO_MASTER_BPM_LIST_FILTERED Grown master bad pixel masks           \n"
"                                                                           \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"There are currently four methods:                                          \n"
"                                                                           \n"
"                                                                           \n"
"(--bpm.method=absolute): Controlled by bpm.kappa_low/_high;                \n"
"                                                                           \n"
"                  The algorithm first collapses the stack of images by     \n"
"                  using the median. Then it subtracts the collapsed master \n"
"                  image and derives the bad pixels by using kappa_low and  \n"
"                  kappa_high as an absolute threshold on the residual      \n"
"                  images.                                                  \n"
"                                                                           \n"
"(--bpm.method=relative): Controlled by bpm.kappa_low/_high;                \n"
"                                                                           \n"
"                  The algorithm first collapses the stack of images by     \n"
"                  using the median. Then it subtracts the collapsed master \n"
"                  image and derives the bad pixels by scaling the rms      \n"
"                  measured on the residual-image - for the rms a properly  \n"
"                  scaled Median Absolute Deviation (MAD) is used.          \n"
"                                                                           \n"
"(--bpm.method=error): Controlled by bpm.kappa_low/_high;                   \n"
"                                                                           \n"
"                  The algorithm first collapses the stack of images by     \n"
"                  using the median. Then it subtracts the collapsed master \n"
"                  image and derives the bad pixels. It uses kappa_low and  \n"
"                  kappa_high to scale the propagated error of each         \n"
"                  individual pixel and compares it with the measured       \n"
"                  residuals.                                               \n"
"                                                                           \n"
"TODO: indicate what data each method is best suited for reduction          \n"
"                                                                           \n"
"The derived bad pixel masks are also filtered by a CLOSING or OPENING      \n"
"filter (see cpl_mask_filter for more details). The filtering process       \n"
"can be controlled by the parameters (--pfx), (--pfy), and (--pfm).         \n"
"                                                                           \n"
"Please note that if no error images are given, the errors are estimated    \n"
"with a shot noise model by using the Gain (--gain) and Ron (--ron)         \n"
"                                                                           \n"
"Please note that part of the code is paralelised. In order to optimise use \n"
"of computing resources you should set the environment variable             \n"
"OMP_NUM_THREADS to a proper value, like (for bash), to use 4 cores         \n"
"export OMP_NUM_THREADS=4                                                   \n";


/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_bpm_3d, HDRLDEMO_BINARY_VERSION, "HDRL Group", 
        PACKAGE_BUGREPORT, "2013", "HDRLDEMO - 3D BPM", 
        hdrldemo_bpm_3d_description);                          
/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_bpm_3d_fill_parameterlist(
        cpl_parameterlist   *   self) 
{                                  
    cpl_parameter   *   par ;

    /* --hdrldemo_bpm_3d.ext-nb-raw */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw", CPL_TYPE_INT,
            "FITS extension of the RAW", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-r");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    
    /* --hdrldemo_bpm_3d.ext-nb-raw-err */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-err", CPL_TYPE_INT,
            "FITS extension of the ERROR", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-e");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_3d.ext-nb-raw-bpm */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-bpm", CPL_TYPE_INT,
            "FITS extension or the input BPM", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-b");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_3d.gain */
    par = cpl_parameter_new_value(RECIPE_NAME".gain", CPL_TYPE_DOUBLE,
            "Gain in [e- / ADU]", RECIPE_NAME, 2.5);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_3d.ron */
    par = cpl_parameter_new_value(RECIPE_NAME".ron", CPL_TYPE_DOUBLE,
            "Read-Out Noise", RECIPE_NAME, 1.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ron");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_3d.post-filter-x */
    par = cpl_parameter_new_value(RECIPE_NAME".post-filter-x", CPL_TYPE_INT,
            "X Size of the post filtering kernel", RECIPE_NAME, 3);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfx");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_3d.post-filter-y */
    par = cpl_parameter_new_value(RECIPE_NAME".post-filter-y", CPL_TYPE_INT,
            "Y Size of the post filtering kernel", RECIPE_NAME, 3);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfy");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_3d.post-filter-mode */
    par = cpl_parameter_new_enum(RECIPE_NAME".post-filter-mode", 
            CPL_TYPE_STRING, "Post filtering mode", RECIPE_NAME,
            "closing", 2, "closing", "dilation");
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfm");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_3d.region-llx/lly/urx/ury */
    hdrl_parameter * deflts = hdrl_rect_region_parameter_create(1, 1, 0, 0);
    cpl_parameterlist * reg_param = hdrl_rect_region_parameter_create_parlist(
                RECIPE_NAME, "", "region-", deflts);
    hdrl_parameter_delete(deflts) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(reg_param) ;
            p != NULL; p = cpl_parameterlist_get_next(reg_param))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(reg_param);

    /* Create BPM_3d parameters */
    deflts = hdrl_bpm_3d_parameter_create(3., 3.,
            HDRL_BPM_3D_THRESHOLD_RELATIVE);
    cpl_parameterlist * bpm_param = hdrl_bpm_3d_parameter_create_parlist(
                RECIPE_NAME, "bpm", deflts) ;
    hdrl_parameter_delete(deflts) ;
    for (cpl_parameter * p = cpl_parameterlist_get_first(bpm_param) ; 
            p != NULL; p = cpl_parameterlist_get_next(bpm_param)) 
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(bpm_param);

    return CPL_ERROR_NONE;
}

static int hdrldemo_bpm_3d(
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
    cpl_frameset            *   in_set ;
    cpl_frameset            *   err_set ;
    cpl_frameset            *   bpm_set ;
    hdrl_parameter          *   bpm_params = NULL;
    hdrl_parameter          *   region_params = NULL;
    cpl_imagelist           *   out_imlist ;
    cpl_imagelist           *   out_imlist_filtered = NULL ;

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
    bpm_params=hdrl_bpm_3d_parameter_parse_parlist(parlist, RECIPE_NAME".bpm") ;
    if (bpm_params == NULL) {
        if (region_params) hdrl_parameter_delete(region_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                "Parsing of the recipe parameters failed");
    }

    /* Load INPUT Data */
    /* Get the frameset */
    if ((in_set = hdrldemo_extract_frameset(frameset, HDRLDEMO_RAW))==NULL){
        hdrl_parameter_delete(bpm_params) ;
        if (region_params) hdrl_parameter_delete(region_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                "Missing RAW files");
    }
    bpm_set = hdrldemo_extract_frameset(frameset, HDRLDEMO_RAW_BPM) ;
    err_set = hdrldemo_extract_frameset(frameset, HDRLDEMO_RAW_ERROR) ;

    /* Load the image list */
    hdrl_imagelist *himlist = hdrl_imagelist_new();
    hdrl_buffer    *buf     = hdrl_buffer_new();
    hdrl_buffer_set_malloc_threshold(buf, 4096);
    if (hdrldemo_hdrl_imagelist_load(in_set, extnum_raw, err_set, extnum_err, 
                bpm_set, extnum_bpm, region_params, ron, gain, himlist,
                buf) != CPL_ERROR_NONE) {
        hdrl_parameter_delete(bpm_params) ;
        if (region_params) hdrl_parameter_delete(region_params) ;
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                "Cannot load RAW image list");
    }

    cpl_frameset_delete(in_set);

    {
        /* Subtracting the median level to all frames. This makes the recipe
         * more robust to inputframes with slightly different countrates */

        cpl_size numframes;
        hdrl_value median;

        numframes = hdrl_imagelist_get_size(himlist);
        for (cpl_size i = 0; i < numframes; i++){
            median = hdrl_image_get_median(hdrl_imagelist_get(himlist, i));
            hdrl_image_sub_scalar(hdrl_imagelist_get(himlist, i), median);
        }
    }


    if (region_params) hdrl_parameter_delete(region_params) ;

    /* BPM Computation */
    if ((out_imlist = hdrl_bpm_3d_compute(himlist, bpm_params)) == NULL) {
        hdrl_parameter_delete(bpm_params) ;
        hdrl_imagelist_delete(himlist);
        return cpl_error_get_code();
    }
    /* Cleanup */
    hdrl_imagelist_delete(himlist);
    
    /* Delete the parameters */
    hdrl_parameter_delete(bpm_params) ;

    /* Post Filtering */
    if (pfx > 0 && pfy > 0) {
        out_imlist_filtered = hdrl_bpm_filter_list(out_imlist, pfx, pfy, 
                filter_mode);
    }

    /* Save the image list */
    hdrldemo_bpm_3d_save_imagelist("HDRLDEMO_MASTER_BPM",
            "hdrldemo_bpm.fits", CPL_TYPE_SHORT, out_imlist, parlist,
            frameset);
    cpl_imagelist_delete(out_imlist);
    if (out_imlist_filtered) {
        /* Save the filtered image list */
        hdrldemo_bpm_3d_save_imagelist("HDRLDEMO_MASTER_BPM_FILTERED",
                "hdrldemo_bpm_filtered.fits", CPL_TYPE_SHORT,
                out_imlist_filtered, parlist, frameset);
        cpl_imagelist_delete(out_imlist_filtered);
    }
    hdrl_buffer_delete(buf);
    return (int)cpl_error_get_code();
}

static cpl_error_code hdrldemo_bpm_3d_save_imagelist(
                const char              *   procatg,
                const char              *   filename,
                const cpl_type              savetype,
                const cpl_imagelist     *   imlist,
                const cpl_parameterlist *   parlist,
                cpl_frameset            *   frameset)
{
    /* Add a QC parameter  */
    cpl_propertylist * qclist = cpl_propertylist_new();

    /* Add the product category  */
    cpl_propertylist_update_string(qclist, CPL_DFS_PRO_CATG, procatg);

    /* Save all bad pixel masks in single fits files */
    cpl_size numframes;
    numframes = cpl_imagelist_get_size(imlist);
    for (cpl_size i = 0; i < numframes; i++){
    	char *outfile = cpl_sprintf("%04d_%s", (int)i, filename);
        cpl_dfs_save_image(frameset, NULL, parlist, frameset, NULL,
                           cpl_imagelist_get_const(imlist, i), savetype,
                           "hdrldemo_bpm_3d", qclist, NULL,
                           PACKAGE "/" PACKAGE_VERSION, outfile);
        cpl_free(outfile);
    }

    cpl_propertylist_delete(qclist);
    return cpl_error_get_code();
}

/**@}*/

