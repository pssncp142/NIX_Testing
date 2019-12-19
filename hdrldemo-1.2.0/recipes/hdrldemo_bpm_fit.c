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
#include <math.h>

#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------*/
/**
  @defgroup hdrldemo_bpm_fit    BPM
  @par Synopsis: TBD
  @par Input frames: TBD
  @par Output frames: TBD
  @code
  @endcode
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_bpm_fit"

static char hdrldemo_bpm_fit_description[] =
"The recipe derives bad pixels on a sequence of images by fitting a         \n"
"polynomial to the data.                                                    \n"
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
"  HDRLDEMO_MASTER_BPM               Master bad pixel mask                  \n"
"  HDRLDEMO_BPM_FIT_CHI2             Chi-squared of the fit                 \n"
"  HDRLDEMO_BPM_FIT_DOF              Degree of Freedom of the fit           \n"
"  HDRLDEMO_BPM_FIT_RED_CHI2         Reduced Chi-squared of the fit         \n"
"                                                                           \n"
"                                                                           \n"
"Usage of the recipe:                                                       \n"
"                                                                           \n"
"There are currently three methods implemented. The switch between the      \n"
"different methods is done by giving the related parameters meaningful      \n"
"values (currently they have a value of -1 (default) if they are not        \n"
"used). Moreover, the degree of the fit can be controlled with the          \n"
"(--degree) recipe parameter. The different methods are                     \n"
"                                                                           \n"
"                                                                           \n"
"(--rel-chi-low/-high): This method marks pixels as bad by using the        \n"
"                       relative chi (low/high) threshold:                  \n"
"                       pixels with values below/above                      \n"
"                       this threshold times measured-rms are marked        \n"
"                       as bad. Please note that the rms is derived         \n"
"                       on the full reduced-chi-squared-image.              \n"
"                       Moreover, the distribution of the chi-squared       \n"
"                       is not symmetric and the thresholds must            \n"
"                       account for that.                                   \n"
"                                                                           \n"
"(--rel-coef-low/-high): This method marks pixels as bad by using the       \n"
"                        relative coefficient (low/high) threshold:         \n"
"                        pixels with values below/above this threshold      \n"
"                        times measured-rms are marked as bad. Please       \n"
"                        note that the rms is derived on the full           \n"
"                        coefficient-image.                                 \n"
"                        The output image encodes which coefficient         \n"
"                        was not in the threshold as a power of two.        \n"
"                        E.g. a value of 5 means coefficient 0 and 2        \n"
"                        were not within the relative thresholds.           \n"
"                                                                           \n"
"( --pval):             This method uses the p-value (between 0% and        \n"
"                       100%) to discriminate between good and bad          \n"
"                       pixels. The p-value is the integral of the chi2     \n"
"                       distributions probability density function          \n"
"                       from the fit-chi2 to infinity. Fits with a          \n"
"                       p-value below the threshold are considered to       \n"
"                       be bad pixels.                                      \n"
"                                                                           \n"
"TODO: indicate what data each method is best suited for reduction          \n"
"                                                                           \n"
"The derived bad pixel mask is also filtered by a CLOSING or OPENING        \n"
"filter (see cpl_mask_filter for more details). The filtering process       \n"
"can be controlled by the parameters (--pfx), (--pfy), and (--pfm).         \n"
"                                                                           \n"
"Please note, that the RAW files should contain the header keyword          \n"
"EXPTIME which is used as the sampling position of the fit. If the          \n"
"RAW_ERROR frame is not present, the errors are estimated with a            \n"
"shot-noise model assuming Poisson distributed data by using the Gain       \n"
"(--gain) and Ron (--ron) parameter. If gain and ron are smaller than       \n"
"zero the square root of the EXPTIME keyword is used as error. The bad      \n"
"pixel map creation parameters must be given in mutually exclusive          \n"
"groups, e.g. abs-rchi2-low AND abs-rchi2-high OR rel-rchi2-low AND         \n"
"rel-rchi2-high OR pval, etc.                                               \n"
"                                                                           \n"
"Please note that part of the code is parallelised. In order to optimise    \n"
"use of computing resources you should set the environment variable         \n"
"OMP_NUM_THREADS to a proper value, like (for bash), to use 4 cores         \n"
"export OMP_NUM_THREADS=4                                                   \n";




/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_bpm_fit, HDRLDEMO_BINARY_VERSION, "HDRL Group",
                  PACKAGE_BUGREPORT, "2014", "HDRLDEMO - BPM FIT",
                  hdrldemo_bpm_fit_description);


/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_bpm_fit_fill_parameterlist(
        cpl_parameterlist   *   self) 
{                                  
    cpl_parameter   *   par ;

    /* --hdrldemo_bpm_fit.ext-nb-raw */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw", CPL_TYPE_INT,
            "FITS extension of the RAW", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-r");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    
    /* --hdrldemo_bpm_fit.ext-nb-raw-err */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-err", CPL_TYPE_INT,
            "FITS extension of the ERROR", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-e");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm_fit.region-llx/lly/urx/ury */
    hdrl_parameter * deflts = hdrl_rect_region_parameter_create(1, 1, 0, 0);
    cpl_parameterlist * reg_param = hdrl_rect_region_parameter_create_parlist(
                RECIPE_NAME, "", "region-", deflts);
    hdrl_parameter_delete(deflts);
    for (cpl_parameter * p = cpl_parameterlist_get_first(reg_param);
            p != NULL; p = cpl_parameterlist_get_next(reg_param))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(reg_param);

    /* --hdrldemo_bpm_fit.ext-nb-raw-bpm */
    par = cpl_parameter_new_value(RECIPE_NAME".ext-nb-raw-bpm", CPL_TYPE_INT,
            "FITS extension or the input BPM", RECIPE_NAME, 0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ext-b");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    deflts = hdrl_bpm_fit_parameter_create_pval(1, 0.1);
    cpl_parameterlist * fpar =
        hdrl_bpm_fit_parameter_create_parlist(RECIPE_NAME, "", deflts) ;
    hdrl_parameter_delete(deflts);
    for (cpl_parameter * p = cpl_parameterlist_get_first(fpar) ;
            p != NULL; p = cpl_parameterlist_get_next(fpar))
        cpl_parameterlist_append(self, cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(fpar);

    /* --hdrldemo_bpm.gain */
    par = cpl_parameter_new_value(RECIPE_NAME".gain", CPL_TYPE_DOUBLE,
            "Gain in [e- / ADU]", RECIPE_NAME, 1.);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_bpm.ron */
    par = cpl_parameter_new_value(RECIPE_NAME".ron", CPL_TYPE_DOUBLE,
            "Read-Out Noise", RECIPE_NAME, 1.);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "ron");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".post-filter-x", CPL_TYPE_INT,
            "X Size of the post filtering kernel", RECIPE_NAME, 3);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfx");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".post-filter-y", CPL_TYPE_INT,
            "Y Size of the post filtering kernel", RECIPE_NAME, 3);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfy");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_enum(RECIPE_NAME".post-filter-mode",
            CPL_TYPE_STRING, "Post filtering mode", RECIPE_NAME,
            "closing", 2, "closing", "dilation");
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "pfm");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    return CPL_ERROR_NONE;
}

static cpl_error_code hdrldemo_bpm_fit_save(
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
            savetype, RECIPE_NAME, qclist, NULL,
            PACKAGE "/" PACKAGE_VERSION, filename);
    cpl_propertylist_delete(qclist);

    return cpl_error_get_code();
}


static cpl_error_code
filter_mask(const cpl_mask * bpm, const cpl_parameterlist * parlist,
            cpl_mask ** filtered)
{
    const cpl_parameter * par;
    int pfx, pfy;
    const char * pfm;
    cpl_filter_mode filter_mode = CPL_FILTER_CLOSING ;
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".post-filter-x");
    pfx = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".post-filter-y");
    pfy = cpl_parameter_get_int(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".post-filter-mode");
    pfm = cpl_parameter_get_string(par);
    if (!strcmp(pfm, "closing")) {
        filter_mode = CPL_FILTER_CLOSING ;
    }
    else if (!strcmp(pfm, "dilation")) {
        filter_mode = CPL_FILTER_DILATION ;
    }
    else {
        *filtered = NULL;
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Filter mode can only be \"closing\" or "
                                     "\"dilation\" (not %s)", pfm);
    }
    /* Post Filtering */
    if (pfx > 0 && pfy > 0) {
        *filtered = hdrl_bpm_filter(bpm, pfx, pfy, filter_mode);
    }
    else {
        *filtered = NULL;
    }

    return cpl_error_get_code();
}


static int hdrldemo_bpm_fit(
        cpl_frameset            *   frameset, 
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter     *   par = NULL;
    int                         extnum_raw;
    int                         extnum_err;
    int                         extnum_bpm;
    int                         degree;
    double                      ron, gain;
    cpl_error_code err = CPL_ERROR_NONE;
    hdrl_parameter * region_params;
    cpl_vector * exptime = NULL;
    hdrl_imagelist * input = NULL;
    hdrl_parameter * fpar = NULL;

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

    /* Parse the Region Parameters */
    region_params = hdrl_rect_region_parameter_parse_parlist(parlist,
            RECIPE_NAME, "region-");
    if (region_params == NULL) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
                "Parsing of the region parameters failed");
    }

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".gain");
    gain = cpl_parameter_get_double(par);
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".ron");
    ron = cpl_parameter_get_double(par);

    fpar = hdrl_bpm_fit_parameter_parse_parlist(parlist, RECIPE_NAME);
    if (fpar == NULL) {
        if (region_params) hdrl_parameter_delete(region_params) ;
        return cpl_error_get_code();
    }
    degree = hdrl_bpm_fit_parameter_get_degree(fpar);

    cpl_frameset * fs_data = cpl_frameset_new();
    cpl_frameset * fs_errs = cpl_frameset_new();
    cpl_frameset * fs_bpms = cpl_frameset_new();
    for (cpl_size i = 0; i < cpl_frameset_get_size(frameset); i++) {
        cpl_frame * frm =
            cpl_frame_duplicate(cpl_frameset_get_position(frameset, i));
        if (!strcmp(cpl_frame_get_tag(frm), HDRLDEMO_RAW)) {
            cpl_frameset_insert(fs_data, frm);
        }
        else if (!strcmp(cpl_frame_get_tag(frm), HDRLDEMO_RAW_ERROR)) {
            cpl_frameset_insert(fs_errs, frm);
        }
        else if (!strcmp(cpl_frame_get_tag(frm), HDRLDEMO_RAW_BPM)) {
            cpl_frameset_insert(fs_bpms, frm);
        }
        else {
            cpl_msg_error(cpl_func,"The following tag is invalid: %s",
                            cpl_frame_get_tag(frm));
            cpl_frame_delete(frm);
            if (region_params) hdrl_parameter_delete(region_params) ;
            if (fpar) hdrl_parameter_delete(fpar) ;
            cpl_frameset_delete(fs_data);
            cpl_frameset_delete(fs_errs);
            cpl_frameset_delete(fs_bpms);

            return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                         "Invalid frame tag");
        }
    }
    if (cpl_frameset_get_size(fs_errs) == 0) {
        cpl_frameset_delete(fs_errs);
        fs_errs =NULL;
    }
    if (cpl_frameset_get_size(fs_bpms) == 0) {
        cpl_frameset_delete(fs_bpms);
        fs_bpms =NULL;
    }

    input = hdrl_imagelist_new();
    exptime = cpl_vector_new(cpl_frameset_get_size(fs_data));
    for (cpl_size i = 0; i < cpl_frameset_get_size(fs_data); i++) {
        cpl_frame * frm_d = cpl_frameset_get_position(fs_data, i);
        cpl_frame * frm_e = fs_errs ?
            cpl_frameset_get_position(fs_errs, i) : NULL;
        cpl_frame * frm_b = fs_bpms ?
            cpl_frameset_get_position(fs_bpms, i) : NULL;
        hdrl_image * hima;
        if (frm_e) {
            ron = -1.;
            gain = -1.;
        }
        cpl_msg_info(cpl_func, "Loading frameset %d", (int)i);

        if (hdrldemo_hdrl_image_load(frm_d, extnum_raw, frm_e, extnum_err, 
                    frm_b, extnum_bpm, region_params, ron, gain,
                    &hima) != CPL_ERROR_NONE) {
            cpl_vector_delete(exptime);
            hdrl_imagelist_delete(input);
            cpl_frameset_delete(fs_data);
            cpl_frameset_delete(fs_errs);
            cpl_frameset_delete(fs_bpms);
            if (region_params) hdrl_parameter_delete(region_params) ;
            if (fpar) hdrl_parameter_delete(fpar) ;
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                         "Cannot load frame");
        }

        cpl_propertylist * plist =
            cpl_propertylist_load(cpl_frame_get_filename(frm_d), 0);
        if (cpl_propertylist_has(plist, "EXPTIME")) {
            cpl_vector_set(exptime, i, cpl_propertylist_get_double(plist, "EXPTIME"));
        }
        else {
            cpl_msg_warning(cpl_func, "No EXPTIME keyword, using %d", (int)i);
            cpl_vector_set(exptime, i, i);
        }
        cpl_propertylist_delete(plist);
        if (frm_e == NULL && ron < 0. && gain <= 0.) {
            double e = sqrt(cpl_vector_get(exptime, i));
            cpl_msg_warning(cpl_func,
                            "No error or ron/gain given assuming "
                            "sqrt(EXPTIME): %g", e);
            hdrl_image_add_scalar(hima, (hdrl_value){0., e});
        }
        hdrl_imagelist_set(input, hima, i);
    }
    cpl_frameset_delete(fs_data);
    cpl_frameset_delete(fs_errs);
    cpl_frameset_delete(fs_bpms);
    hdrl_parameter_delete(region_params);
    cpl_msg_info(cpl_func, "Beginning polynomial fit");

    cpl_image * out_chi2, * out_dof;
    hdrl_imagelist * out_coef;
    err = hdrl_fit_polynomial_imagelist(input, exptime, degree,
                                        &out_coef, &out_chi2, &out_dof);
    if (err) {
        if (fpar) hdrl_parameter_delete(fpar) ;
        cpl_vector_delete(exptime);
        hdrl_imagelist_delete(input);
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                     "Fit failed");
    }

    cpl_msg_info(cpl_func, "Storing results");
    hdrldemo_bpm_fit_save("HDRLDEMO_BPM_FIT_CHI2",
                          "hdrldemo_bpm_fit_chi2.fits", CPL_TYPE_FLOAT,
                          out_chi2, parlist, frameset);
    hdrldemo_bpm_fit_save("HDRLDEMO_BPM_FIT_DOF",
                          "hdrldemo_bpm_fit_dof.fits", CPL_TYPE_FLOAT,
                          out_dof, parlist, frameset);
    cpl_image_divide(out_chi2, out_dof);
    hdrldemo_bpm_fit_save("HDRLDEMO_BPM_FIT_RED_CHI2",
                          "hdrldemo_bpm_fit_redchi2.fits", CPL_TYPE_FLOAT,
                          out_chi2, parlist, frameset);
    if (cpl_image_count_rejected(out_chi2) ==
        (cpl_image_get_size_x(out_chi2) * cpl_image_get_size_y(out_chi2))) {
        cpl_msg_error(cpl_func, "Too few good pixels to fit polynomial of "
                      "degree %d in all pixels", degree);
        goto end;
    }
    for (cpl_size i = 0; i < hdrl_imagelist_get_size(out_coef); i++) {
        hdrl_image * hd = hdrl_imagelist_get(out_coef, i);
        cpl_image * d = hdrl_image_get_image(hd);
        cpl_image * e = hdrl_image_get_error(hd);
        char fnbuffer[256];
        char tagbuffer[256];
        cpl_msg_info(cpl_func, "Coefficient %d:", (int)i);
        cpl_msg_info(cpl_func, "  Mean: %g",
                     cpl_image_get_mean(d));
        cpl_msg_info(cpl_func, "  Standard deviation: %g",
                     cpl_image_get_stdev(d));
        cpl_msg_info(cpl_func, "  Mean error of fit: %g", cpl_image_get_mean(e));
        sprintf(fnbuffer,"hdrldemo_bpm_fit_coeff%d.fits", (int)i);
        sprintf(tagbuffer,"HDRLDEMO_BPM_FIT_COEF%d", (int)i);
        hdrldemo_bpm_fit_save(tagbuffer, fnbuffer, CPL_TYPE_FLOAT, d,
                              parlist, frameset);
        cpl_image_save(e, fnbuffer, CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
    }
    cpl_msg_info(cpl_func, "Mean of reduced chi2: %g", cpl_image_get_mean(out_chi2));
    cpl_msg_info(cpl_func, "Median of reduced chi2: %g", cpl_image_get_median(out_chi2));
    cpl_msg_info(cpl_func, "Standard deviation of reduced chi2: %g", cpl_image_get_stdev(out_chi2));

    {
        cpl_image * bpm;
        /* fits again but has less output */
        hdrl_bpm_fit_compute(fpar, input, exptime, &bpm);

        cpl_mask * m = cpl_mask_threshold_image_create(bpm, 0, INT_MAX);
        size_t n = cpl_mask_count(m);
        double p = (double)n / (cpl_mask_get_size_x(m) * cpl_mask_get_size_y(m));
        cpl_msg_info(cpl_func, "%zu bad pixels (%g%%)", n, p * 100.);
        hdrldemo_bpm_fit_save("HDRLDEMO_MASTER_BPM",
                              "hdrldemo_bpm_fit_bpm.fits", CPL_TYPE_INT, bpm,
                              parlist, frameset);

        cpl_mask * filtered;
        filter_mask(m, parlist, &filtered);
        if (filtered) {
            cpl_image * imf = cpl_image_new_from_mask(filtered);
            cpl_mask_delete(filtered) ;
            hdrldemo_bpm_fit_save("HDRLDEMO_MASTER_BPM_FILTERED",
                                  "hdrldemo_bpm_fit_bpm_filtered.fits", CPL_TYPE_INT,
                                  imf, parlist, frameset);
            cpl_image_delete(imf) ;
        }

        cpl_mask_delete(m);
        cpl_image_delete(bpm);
    }

end:

    hdrl_imagelist_delete(out_coef);
    cpl_image_delete(out_chi2);
    cpl_image_delete(out_dof);
    cpl_vector_delete(exptime);
    hdrl_imagelist_delete(input);
    hdrl_parameter_delete(fpar);

    return (int)cpl_error_get_code();
}

/**@}*/

