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
#include <string.h>

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrldemo_response response calculation
 * @par Synopsis: The recipe calculates the response from observations of a
 * standard star.
 * @par Input frames:
 *
 *  DO category:                     Explanation:                 Required:
 *  SPECTRUM_1D                      Flux                         Yes
 *  ATM_EXT                          Atmospheric Extinction       Yes
 *  FLUX_CATG                        Std Stars Catalog            Yes
 *  TELLURIC_CATG					 Catalog of Telluric models   No
 *  FIT_AREAS						 Areas used for Telluric 	  No
 *  								 model fitting
 *  QUALITY_AREAS					 Areas where the quality of   No
 *  								 Telluric fitting is
 *  								 calculated
 *  HIGH_ABS_REGIONS				 High absorption regions 	  No
 *  INTERPOL_POINTS_RESPONSE		 Points used for the
 *  								 response interpolation.      Yes
 *
 * @par Output frames:
 *
 *  DO category:                      Explanation:
 *  HDRLDEMO_RESPONSE                 Response
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_response"

const hdrl_spectrum1D_wave_scale
global_scale = hdrl_spectrum1D_wave_scale_linear;

static char hdrldemo_response_description[] =
"                                                                           \n"
"The recipe derives response from observations of the standard star.        \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"DO category:                    Explanation:                Required:      \n"
"SPECTRUM_1D                     Flux                        Yes            \n"
"ATM_EXT                         Atmospheric Extinction      Yes            \n"
"FLUX_CATG                       Std Stars Catalog           Yes            \n"
"TELLURIC_CATG					 Catalog of Telluric models  No             \n"
"FIT_AREAS						 Areas used for Telluric 	 No             \n"
"                                model fitting                              \n"
"QUALITY_AREAS					 Areas where the quality of  No             \n"
"                                Telluric fitting is                        \n"
"                                calculated                                 \n"
"HIGH_ABS_REGIONS				 High absorption regions 	 No             \n"
"INTERPOL_POINTS_RESPONSE		 Points used for the                        \n"
"                                response interpolation.     Yes            \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"DO category:                      Explanation:                             \n"
"HDRLDEMO_RESPONSE                 Response                                 \n"
"                                                                           \n"
"This recipe has polymorphic behaviour, it calculates the response given    \n"
"the observed std star spectrum, the catalog of the models of the std stars.\n"
"The recipe selects the correct std star model from FLUX_CATG using the     \n"
"metadata of the observed std star spectrum.                                \n"
"The recipe can, optionally, perform telluric correction if TELLURIC_CATG,  \n"
"QUALITY_AREAS and FIT_AREAS are provided."
"Lastly, the recipe can also correct the std star flux for doppler shift.	\n"
"Usage of the recipe:                                                       \n"
"                                                                           \n"
"The recipe calculates the response using the following formula:            \n"
"resp = (Tex*G*Istd(l)*10^(0.4*(Ap-Am)*Ex(l)))/(Iobs(l))        			\n"
"where:                                                                     \n"
"Iobs(l): observed std star spectrum, provided by SPECTRUM_1D. This spectrum\n"
"can, optionally, be telluric corrected.									\n"
"Istd(l): model of the std star spectrum, provided by FLUX_CATG. This       \n"
"spectrum can be, optionally, doppler-shift corrected.						\n"
"Ex(l): model of the atmospheric extinction, selected from ATM_EXT          \n"
"Am: air mass, read from the FITS file header                               \n"
"Ap: air mass correction factor, hardcoded to 0.0                           \n"
"G: gain, read from the FITS file header                                    \n"
"Tex: exposure time, read from FITS file header                             \n";

/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_response,
                  HDRLDEMO_BINARY_VERSION,
                  "HDRL Group",
                  PACKAGE_BUGREPORT,
                  "2017",
                  "HDRLDEMO - Response computation",
                  hdrldemo_response_description);


/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_response_fill_parameterlist(
                cpl_parameterlist   *   self) {

    cpl_parameter   *   par ;

    /* hdrldemo_response window for noise calculation */
    par = cpl_parameter_new_value(RECIPE_NAME".der-snr-half-window", CPL_TYPE_INT,
        "Half window used for noise calculation following the DER-SNR approach."
        "\nA good value is 2.5 x the resolving power of a spectral line.",
        RECIPE_NAME, 10);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "der-snr-half-window");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* hdrldemo_response sampling period, to calculate wavelengths */
    par = cpl_parameter_new_value(RECIPE_NAME".sampling-period", CPL_TYPE_DOUBLE,
            "Conversion factor to put the wavelength in [nm]", RECIPE_NAME, 1.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "sampling-period");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* ------------------ Telluric Calculation parameters ------------------ */

    par = cpl_parameter_new_value(RECIPE_NAME".telluric-xcorr-step",
            CPL_TYPE_DOUBLE, "Cross correlation step", RECIPE_NAME,
            10.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "telluric-xcorr-step");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".telluric-xcorr-half-window",
            CPL_TYPE_INT, "Search window for telluric-xcorr", RECIPE_NAME,
            200.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "telluric-xcorr-half-window");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);


    par = cpl_parameter_new_value(RECIPE_NAME".telluric-xcorr-normalize",
            CPL_TYPE_BOOL, "Normalize telluric-xcorr", RECIPE_NAME,
            CPL_TRUE);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "telluric-xcorr-normalize");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);


    par = cpl_parameter_new_value(RECIPE_NAME".telluric-xcorr-l-min",
            CPL_TYPE_DOUBLE, "Cross correlation lmin", RECIPE_NAME,
            7.);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "telluric-xcorr-l-min");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".telluric-xcorr-l-max",
            CPL_TYPE_DOUBLE, "Cross correlation lmax", RECIPE_NAME,
            7.1);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "telluric-xcorr-l-max");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* ------------------ Velocity Calculation parameters ------------------ */
    par = cpl_parameter_new_value(RECIPE_NAME".velocity-wguess",
            CPL_TYPE_DOUBLE,
			"Velocity Calculation: reference line wavelength position",
			RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "velocity-wguess");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".velocity-enable",
			CPL_TYPE_BOOL,
			"Velocity Calculation: enable",
			RECIPE_NAME, CPL_FALSE);
	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "velocity-enable");
	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
	cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".velocity-range-wmin",
            CPL_TYPE_DOUBLE,
			"Velocity Calculation: minimum of wavelength box for line fit ",
			RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "velocity-range-wmin");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".velocity-range-wmax",
            CPL_TYPE_DOUBLE,
			"Velocity Calculation: maximum of wavelength box for line fit ",
			RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "velocity-range-wmax");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".velocity-fit-wmin",
            CPL_TYPE_DOUBLE,
			"Velocity Calculation: minimum wavelength value used to fit line slope",
			RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "velocity-fit-wmin");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".velocity-fit-wmax",
            CPL_TYPE_DOUBLE,
			"Velocity Calculation: maximum wavelength value used to fit line slope",
			RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "velocity-fit-wmax");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".velocity-fit-half-win",
            CPL_TYPE_DOUBLE,
			"Velocity Calculation: half box where polynomial fit is performed",
			RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "velocity-fit-half-win");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    par = cpl_parameter_new_value(RECIPE_NAME".response-wrange-median-final-interpolation",
            CPL_TYPE_DOUBLE,
			"Response calculation: wavelength range to be taken around an interpolation point",
			RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-wrange-median-final-interpolation");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

   	/*-------------------------- Exposure Time -------------------------*/
   	par = cpl_parameter_new_value(RECIPE_NAME".response-exposure-time", CPL_TYPE_DOUBLE,
   			"Exposure time in [s]", RECIPE_NAME, 0.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-exposure-time");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);

   	par = cpl_parameter_new_value(RECIPE_NAME".response-exposure-time-error",
   			CPL_TYPE_DOUBLE, "Exposure time error in [s]", RECIPE_NAME, 0.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-exposure-time-error");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);
   	/*-------------------------- Exposure Time -------------------------*/

   	/*----------------------------- Airmass ----------------------------*/
   	par = cpl_parameter_new_value(RECIPE_NAME".response-airmass", CPL_TYPE_DOUBLE,
   			"Airmass at which the standard star was observed", RECIPE_NAME, 1.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-airmass");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);

   	par = cpl_parameter_new_value(RECIPE_NAME".response-airmass-error", CPL_TYPE_DOUBLE,
   		 "Error on the airmass at which the standard star was observed",
   		 RECIPE_NAME, 0.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-airmass-error");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);
   	/*----------------------------- Airmass ----------------------------*/

   	/*----------------------- Airmass correction -----------------------*/
   	par = cpl_parameter_new_value(RECIPE_NAME".response-airmass-correction",
   			CPL_TYPE_DOUBLE, "Parameter to indicate if the response is "
   			"computed at arimass=0, or at a given non-zero value",
   			RECIPE_NAME, 0.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-airmass-correction");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);

   	par = cpl_parameter_new_value(RECIPE_NAME".response-airmass-correction-error",
   			CPL_TYPE_DOUBLE, "Error on the airmass-correction parameter",
   			RECIPE_NAME, 0.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI,
   			"response-airmass-correction-error");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);
   	/*----------------------- Airmass correction -----------------------*/

   	/*------------------------------ Gain ------------------------------*/
   	par = cpl_parameter_new_value(RECIPE_NAME".response-gain", CPL_TYPE_DOUBLE,
   		   "Detector gain in [e/ADU]", RECIPE_NAME, 1.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-gain");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);

   	par = cpl_parameter_new_value(RECIPE_NAME".response-gain-error",
   			CPL_TYPE_DOUBLE, "Error on the detector gain in [e/ADU]",
			RECIPE_NAME, 0.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "response-gain-error");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);
   	/*------------------------------ Gain ------------------------------*/


   	/*--------------------- Flux scaling factor ------------------------*/
   	par = cpl_parameter_new_value(RECIPE_NAME".flux-scaling-factor",
   			CPL_TYPE_DOUBLE, "Factor to convert the observed flux from "
   		   "actual bin size to Angstrom units", RECIPE_NAME, 1.0);
   	cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "flux-scaling-factor");
   	cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   	cpl_parameterlist_append(self, par);
   	/*--------------------- Flux scaling factor ------------------------*/

    return CPL_ERROR_NONE;
}


static double get_airmass(cpl_propertylist * header){

    const double airmass_start =
            cpl_propertylist_get_double(header, "ESO TEL AIRM START");
    const double airmass_end =
            cpl_propertylist_get_double(header, "ESO TEL AIRM END");

    return 0.5 * (airmass_start + airmass_end);
}

static double
propertylist_get_double_if_available(const cpl_propertylist * l, const char *name){
	if(!cpl_propertylist_has(l, name)) return NAN;
	return cpl_propertylist_get_double(l, name);
}

static void get_gain_exptime(const cpl_propertylist * header,
                double * gain, double * exptime){

    const char* instr_name = cpl_propertylist_get_string(header, "INSTRUME");

    const char * bname = !strcmp(instr_name, "XSHOOTER")
            ? "ESO SEQ ARM" : "ESO INS FILT1 NAME";

    cpl_boolean need_hardcoded_data = !strcmp(instr_name, "SINFONI");

    need_hardcoded_data |= !strcmp(instr_name, "SYNTHETIC");

    need_hardcoded_data |= !strcmp(instr_name, "XSHOOTER")
            && !strcmp(cpl_propertylist_get_string(header, bname), "NIR");

    if (need_hardcoded_data){
        *gain = 2.12;
        *exptime = propertylist_get_double_if_available(header,"ESO DET DIT");
    }
    else{
        *gain = propertylist_get_double_if_available(header, "ESO DET OUT1 GAIN");
        *exptime = propertylist_get_double_if_available(header, "ESO DET WIN1 DIT1");
    }
}

static double
override_default(const cpl_parameterlist * parlist, const double def,
        const char * name){

    const cpl_parameter * par =
            cpl_parameterlist_find_const(parlist, name);

    const int par_flag = cpl_parameter_get_default_flag(par);
    if(!par_flag)
        return def;

    return cpl_parameter_get_double(par);
}

static cpl_error_code get_data_from_header(
    const cpl_frameset * obs_stars, cpl_size * n_ext,
    double * amass, double * gain, double * etime,
    const cpl_parameterlist * parlist){

    const cpl_frame * obs_star =
            cpl_frameset_get_position_const(obs_stars, 0);

    *n_ext = cpl_frame_get_nextensions(obs_star) + 1;

    /*Find instrument*/
    cpl_propertylist * header =
            cpl_propertylist_load(cpl_frame_get_filename(obs_star), 0);

    if(!header)
        return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                       "header is NULL");

    double gain_this = 0.0;
    double etime_this = 0.0;
    double amass_this = get_airmass(header);

    get_gain_exptime(header, &gain_this, &etime_this);

    *gain = override_default(parlist, gain_this, RECIPE_NAME".response-gain");
    *etime = override_default(parlist, etime_this, RECIPE_NAME".response-exposure-time");
    *amass = override_default(parlist, amass_this, RECIPE_NAME".response-airmass");

    cpl_propertylist_delete(header);

    cpl_ensure_code(!isnan(*gain) && !isnan(*etime)
    		&& !isnan(*amass), CPL_ERROR_ILLEGAL_OUTPUT);

    return cpl_error_get_code();
}

static hdrl_spectrum1D * read_atm_ext(const cpl_frameset * frameset_extintion_cat){

    const char * flux_column = "EXTINCTION";
    const char * wavelength_column = "LAMBDA";

    const cpl_frame * atm_ext_frame =
                    cpl_frameset_get_position_const(frameset_extintion_cat, 0);
    const char * fname = cpl_frame_get_filename(atm_ext_frame);

    cpl_table * atm_ext_table = cpl_table_load(fname, 1, 0);

    hdrl_spectrum1D * sp = hdrl_spectrum1D_convert_from_table
            (atm_ext_table, flux_column, wavelength_column, NULL, NULL,
                    global_scale);
    cpl_table_delete(atm_ext_table);

    return sp;
}



static hdrl_spectrum1D *
read_std_star_model(const cpl_frameset * frameset_std_flux_cat,
        const double ra,
        const double dec,
        const double ra_tolerance,
        const double dec_tolerance){

    const char* col_name_extid = "ext_id";
    const char* col_name_ra    = "ra";
    const char* col_name_dec   = "dec";

    const char * flux_column = "FLUX";
    const char * wavelength_column = "LAMBDA";

    const cpl_frame * std_flux_cat_frame =
                    cpl_frameset_get_position_const(frameset_std_flux_cat, 0);
    const char * fname = cpl_frame_get_filename(std_flux_cat_frame);
    cpl_table* table_cat = cpl_table_load(fname, 1, 0);
    int rej = 0;

    cpl_table * model_extension = NULL;

    for(cpl_size i = 0; i < cpl_table_get_nrow(table_cat); ++i){

        const int ext_id =
                cpl_table_get_int(table_cat, col_name_extid, i ,&rej);
        const double curr_ra = cpl_table_get(table_cat, col_name_ra, i,&rej);
        const double curr_dec = cpl_table_get(table_cat, col_name_dec, i,&rej);

        if ((ext_id > 0) && (fabs(curr_ra - ra) < ra_tolerance)
                && (fabs(curr_dec - dec) < dec_tolerance)){
            model_extension = cpl_table_load(fname, ext_id, 0);
            break;
        }
    }
    cpl_table_delete(table_cat);

    if(model_extension == NULL) return NULL;

    hdrl_spectrum1D * spec = hdrl_spectrum1D_convert_from_table
            (model_extension, flux_column, wavelength_column, NULL, NULL,
                    global_scale);

    cpl_table_delete(model_extension);

    return spec;
}

static void
get_ra_dec(const cpl_frameset * obs_star_frameset, double * ra, double * dec){

    const cpl_frame * frm_obs =
                    cpl_frameset_get_position_const(obs_star_frameset, 0);
    const char * name_obs = cpl_frame_get_filename(frm_obs);

    cpl_propertylist * plist = cpl_propertylist_load(name_obs ,0);

    *ra = cpl_propertylist_get_double(plist,"RA");
    *dec = cpl_propertylist_get_double(plist,"DEC");

    cpl_propertylist_delete(plist);
}

static void
save_tab(const cpl_table * spectrum_table, const cpl_parameterlist * parlist,
        const cpl_propertylist * app_list, cpl_frameset * frameset,
        cpl_frame * raw_frame, const char * fname, const char * tag){

    cpl_propertylist * applist = cpl_propertylist_new();

    cpl_propertylist_update_string(applist, CPL_DFS_PRO_CATG,
            tag);

    cpl_dfs_save_table(frameset,
            NULL,
            parlist,
            frameset,
            raw_frame,
            spectrum_table,
            app_list,
            RECIPE_NAME,
            applist,
            NULL,
            PACKAGE "/" PACKAGE_VERSION,
            fname);

    cpl_propertylist_delete(applist);
}


static hdrl_spectrum1D *
read_spectrum(const cpl_frameset * obs_star_fset, const hdrl_data_t sampling_factor,
const cpl_size der_snr_half_window){

    const cpl_frame * obs_star =
            cpl_frameset_get_position_const(obs_star_fset, 0);
    const char * fname = cpl_frame_get_filename(obs_star);

    cpl_table* tab_obs_std_star=cpl_table_load(fname,1,0);
    cpl_table_multiply_scalar(tab_obs_std_star,"wavelength",sampling_factor);

    hdrl_spectrum1D * table_spectrum = hdrl_spectrum1D_convert_from_table(tab_obs_std_star,
            "counts_bkg", "wavelength", NULL, NULL, global_scale);

    const cpl_image * flx =
            hdrl_image_get_image_const(hdrl_spectrum1D_get_flux(table_spectrum));
    const cpl_array * wlen =
                        hdrl_spectrum1D_get_wavelength(table_spectrum).wavelength;
    const hdrl_spectrum1D_wave_scale scale =
                            hdrl_spectrum1D_get_wavelength(table_spectrum).scale;

    hdrl_spectrum1D * to_ret = hdrl_spectrum1D_create_error_DER_SNR(flx,
                        der_snr_half_window, wlen, scale);

    cpl_table_delete(tab_obs_std_star);
    hdrl_spectrum1D_delete(&table_spectrum);
    return to_ret;
}

static hdrl_spectrum1Dlist *
get_telluric_models(const cpl_frameset * set){

	if(cpl_frameset_get_size(set) == 0) return NULL;

    const cpl_frame * telluric_cat =
            cpl_frameset_get_position_const(set, 0);
    const char * cat_name = cpl_frame_get_filename(telluric_cat);
    const cpl_size next = cpl_frame_get_nextensions(telluric_cat);

    hdrl_spectrum1Dlist * list = hdrl_spectrum1Dlist_new();

    for(cpl_size i = 0; i < next; ++i){
        cpl_table * tab = cpl_table_load(cat_name, 1 + i, 1);

        hdrl_spectrum1D * s =
                hdrl_spectrum1D_convert_from_table(tab,
                        "flux",
                        "lam",
                        NULL,
                        NULL,
                        global_scale);
        cpl_table_delete(tab);
        /*um to nm*/
        hdrl_spectrum1D_wavelength_mult_scalar_linear(s, 1e3);
        hdrl_spectrum1Dlist_set(list, s, i);
    }

    return list;
}

static cpl_bivector *
read_wlen_windows(const cpl_frameset * areas_fset){

	if(cpl_frameset_get_size(areas_fset) == 0)
		return NULL;

	cpl_ensure(cpl_frameset_get_size(areas_fset) == 1,
			CPL_ERROR_ILLEGAL_INPUT, NULL);

    const cpl_frame * fit_areas_frm =
            cpl_frameset_get_position_const(areas_fset, 0);
    cpl_table * tab = NULL;

    if(areas_fset !=NULL) {
        tab =
                cpl_table_load(cpl_frame_get_filename(fit_areas_frm),1,0);
    }
    else{
        return NULL;
    }

    const cpl_size nrow = cpl_table_get_nrow(tab);
    double * pwmin = cpl_table_unwrap(tab,"LAMBDA_MIN");
    double * pwmax = cpl_table_unwrap(tab,"LAMBDA_MAX");

    cpl_vector * wmin = cpl_vector_wrap(nrow, pwmin);
    cpl_vector * wmax = cpl_vector_wrap(nrow, pwmax);
    cpl_bivector * to_ret = cpl_bivector_wrap_vectors(wmin, wmax);
    cpl_table_delete(tab);
    return to_ret;
}

static cpl_array *	read_fit_points(const cpl_frameset * fset,
		const double ra, const double dec, const double ra_dec_tolerance){

	const cpl_frame * f = cpl_frameset_get_position_const(fset, 0);
	const char * fname = cpl_frame_get_filename(f);

	cpl_table * index_table = cpl_table_load(fname, 1, CPL_FALSE);
	const cpl_size sz = cpl_table_get_nrow(index_table);

	cpl_size selected_ext = -1;

	for(cpl_size i = 0; i < sz; ++i){
		int rej;
		const int ext_id = cpl_table_get_int(index_table, "ext_id", i ,&rej);
		const double curr_ra = cpl_table_get(index_table, "ra", i, &rej);
		const double curr_dec = cpl_table_get(index_table, "dec", i, &rej);
		if ((ext_id > 0) && (fabs(curr_ra - ra) < ra_dec_tolerance) &&
				(fabs(curr_dec - dec) < ra_dec_tolerance)){
			selected_ext = ext_id;
		}
	}

	cpl_table_delete(index_table);
	cpl_ensure(selected_ext >= 0, CPL_ERROR_ILLEGAL_OUTPUT, NULL);


	cpl_table * tb = cpl_table_load(fname, selected_ext, CPL_FALSE);

	const cpl_size sz_points = cpl_table_get_nrow(tb);
	cpl_array * fit_points = cpl_array_new(sz_points, CPL_TYPE_DOUBLE);

	for(cpl_size i = 0; i < sz_points; ++i){
		int rej = 0;
		const double d = cpl_table_get(tb, "LAMBDA", i, &rej);
		cpl_array_set(fit_points, i, d);
	}

	cpl_table_delete(tb);
	return fit_points;
}


/*----------------------------------------------------------------------------*/
/**
  @brief   perform master fringe combination
  @param   frameset   input set of frames
  @param   parlist    input recipe parameters
  @return  cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_response(
                cpl_frameset            *   frameset,
                const cpl_parameterlist *   parlist)
{
    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    const double ra_dec_tolerance = 0.0166667;

    const cpl_parameter * par;

    /* Get parameters*/
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".der-snr-half-window");
    const cpl_size der_snr_half_window = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".sampling-period");
    const double sampling_factor = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".telluric-xcorr-step");
    const double x_corr_w_step = cpl_parameter_get_double(par);

    /*Get parameters - Cross Correlation for Telluric*/

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".telluric-xcorr-half-window");
    const cpl_size xcorr_half_win = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".telluric-xcorr-normalize");
    const cpl_boolean xcorr_normalize = (cpl_boolean) cpl_parameter_get_bool(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".telluric-xcorr-l-min");
    const hdrl_data_t lmin = (hdrl_data_t)cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".telluric-xcorr-l-max");
    const hdrl_data_t lmax = (hdrl_data_t)cpl_parameter_get_double(par);

    /*Get parameters - Velocity Estimation*/

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".velocity-enable");
	const cpl_boolean velocity_enable = cpl_parameter_get_bool(par);

	par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".velocity-wguess");
	const hdrl_data_t velocity_wguess = cpl_parameter_get_double(par);

	par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".velocity-range-wmin");
	const hdrl_data_t velocity_range_wmin = cpl_parameter_get_double(par);

	par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".velocity-range-wmax");
	const hdrl_data_t velocity_range_wmax = cpl_parameter_get_double(par);

	par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".velocity-fit-wmin");
	const hdrl_data_t velocity_fit_wmin = cpl_parameter_get_double(par);

	par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".velocity-fit-wmax");
	const hdrl_data_t velocity_fit_wmax = cpl_parameter_get_double(par);

	par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".velocity-fit-half-win");
	const hdrl_data_t velocity_fit_half_win = cpl_parameter_get_double(par);

	par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".response-wrange-median-final-interpolation");
	const hdrl_data_t fit_wrange = cpl_parameter_get_double(par);

    /*Get parameters - Response*/
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".response-airmass-correction");
    const double Ap = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".response-airmass-error");
    const double Am_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist,
            RECIPE_NAME".response-airmass-correction-error");
    const double Ap_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist,
    RECIPE_NAME".response-exposure-time-error");
    const double Tex_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".response-gain-error");
    const double G_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".flux-scaling-factor");
	const double flux_scaling = cpl_parameter_get_double(par);

    cpl_size fsize = cpl_frameset_get_size(frameset);
    cpl_frameset * frameset_flux = cpl_frameset_new();
    cpl_frameset * frameset_telluric_cat = cpl_frameset_new();
    cpl_frameset * frameset_fit_areas = cpl_frameset_new();
    cpl_frameset * frameset_quality_areas = cpl_frameset_new();
    cpl_frameset * frameset_high_abs_regions = cpl_frameset_new();
    cpl_frameset * frameset_interpolation_points = cpl_frameset_new();
    cpl_frameset * frameset_extintion_cat = cpl_frameset_new();
    cpl_frameset * frameset_std_flux_cat = cpl_frameset_new();

    cpl_frame * spectrum_frame = NULL;
    for(cpl_size i = 0; i < fsize; ++i){
        cpl_frame * cur_frame = cpl_frameset_get_position(frameset, i);

        if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_SPECTRUM_1D)) {
            spectrum_frame = cpl_frame_duplicate(cur_frame);
            cpl_frameset_insert(frameset_flux,
                    spectrum_frame);
            cpl_msg_info(cpl_func,"Observed spectrum found in the SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_TELL_CATG)) {
            cpl_frameset_insert(frameset_telluric_cat,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Telluric models catalog found in SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame),  HDRLDEMO_FIT_AREAS)) {
            cpl_frameset_insert(frameset_fit_areas,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Fit areas found in SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_QUALITY_AREAS)) {
            cpl_frameset_insert(frameset_quality_areas,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Quality areas found in SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_HIGH_ABS_REGIONS)) {
            cpl_frameset_insert(frameset_high_abs_regions,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"High ABS regions found in SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_INT_POINTS_RESPONSE)) {
            cpl_frameset_insert(frameset_interpolation_points,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Interpolation points found in SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_ATM_EXT)) {
            cpl_frameset_insert(frameset_extintion_cat,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Atm Extinction found in SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_FLUX_CATG)) {
            cpl_frameset_insert(frameset_std_flux_cat,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Std star catalog found in SOF");
        } else{
            cpl_msg_info(cpl_func, "SOF contains unknown tags");
        }
    }

    if(cpl_frameset_get_size(frameset_flux) != 1 ||
			cpl_frameset_get_size(frameset_high_abs_regions) > 1 ||
			cpl_frameset_get_size(frameset_interpolation_points) != 1 ||
			cpl_frameset_get_size(frameset_extintion_cat) != 1 ||
			cpl_frameset_get_size(frameset_std_flux_cat) != 1){

        cpl_frameset_delete(frameset_flux);
        cpl_frameset_delete(frameset_telluric_cat);
        cpl_frameset_delete(frameset_fit_areas);
        cpl_frameset_delete(frameset_quality_areas);
        cpl_frameset_delete(frameset_high_abs_regions);
        cpl_frameset_delete(frameset_interpolation_points);
        cpl_frameset_delete(frameset_extintion_cat);
        cpl_frameset_delete(frameset_std_flux_cat);

        return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                "SOF contains wrong number of frames");
    }

    if(cpl_frameset_get_size(frameset_telluric_cat) != 1 ||
            cpl_frameset_get_size(frameset_fit_areas) != 1 ||
            cpl_frameset_get_size(frameset_quality_areas) != 1){
        cpl_msg_info(cpl_func, "Data needed for telluric correction was not "
        		"provided, telluric correction will be skipped");
    }

	/*Get parameters - Response core calculation*/
    double Am = 0.0;
    double Tex = 0.0;
    double G = 0.0;

    cpl_size n_ext = 0;

    cpl_error_code fail = get_data_from_header
            (frameset_flux, &n_ext, &Am, &G, &Tex, parlist);

    if(fail){
		cpl_frameset_delete(frameset_flux);
        cpl_frameset_delete(frameset_telluric_cat);
        cpl_frameset_delete(frameset_fit_areas);
        cpl_frameset_delete(frameset_quality_areas);
        cpl_frameset_delete(frameset_high_abs_regions);
        cpl_frameset_delete(frameset_interpolation_points);
        cpl_frameset_delete(frameset_extintion_cat);
        cpl_frameset_delete(frameset_std_flux_cat);

        return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                "Unable to extract header data");
    }

    hdrl_spectrum1Dlist * telluric_models =
            get_telluric_models(frameset_telluric_cat);

    hdrl_spectrum1D * obs_spectrum =
            read_spectrum(frameset_flux, sampling_factor, der_snr_half_window);

    hdrl_spectrum1D_div_scalar(obs_spectrum, (hdrl_value){ flux_scaling, 0.0 });

    cpl_bivector * fit_areas = read_wlen_windows(frameset_fit_areas);
    cpl_bivector * quality_areas = read_wlen_windows(frameset_quality_areas);

    cpl_bivector * high_abs_regions = read_wlen_windows(frameset_high_abs_regions);

    /*fit_areas and quality_areas are NOT mandatory*/
    if(!high_abs_regions){
        cpl_msg_info(cpl_func, "No high absorption regions provided, the execution will continue.");
    }

    double ra = 0;
    double dec = 0;
    get_ra_dec(frameset_flux, &ra, &dec);

    hdrl_spectrum1D * ref_spectrum = read_std_star_model(frameset_std_flux_cat,
    		ra, dec, ra_dec_tolerance, ra_dec_tolerance);

    hdrl_spectrum1D * E_x = read_atm_ext(frameset_extintion_cat);

    cpl_msg_debug(cpl_func, "Number of extensions: %lld, sampling period = %f",
            n_ext, sampling_factor);

    cpl_msg_debug(cpl_func, "Airmass = %f +/- %f", Am, Am_e);
    cpl_msg_debug(cpl_func, "Airmass correction = %f +/- %f", Ap, Ap_e);
    cpl_msg_debug(cpl_func, "Gain = %f +/- %f", G, G_e);
    cpl_msg_debug(cpl_func, "Exposure = %f +/- %f", Tex, Tex_e);

    hdrl_value hAp  = {(hdrl_data_t)Ap,   (hdrl_error_t)Ap_e };
    hdrl_value hAm  = {(hdrl_data_t)Am,   (hdrl_error_t)Am_e };
    hdrl_value hG   = {(hdrl_data_t)G,    (hdrl_error_t)G_e  };
    hdrl_value hTex = {(hdrl_data_t)Tex,  (hdrl_error_t)Tex_e};
	hdrl_parameter * resp_calc_par = hdrl_response_parameter_create(hAp, hAm, hG, hTex);

	cpl_array 	* fit_points = read_fit_points(frameset_interpolation_points,
    		ra, dec, ra_dec_tolerance);

    hdrl_parameter * telluric_par = NULL;
    if(telluric_models != NULL && quality_areas != NULL && fit_areas != NULL)
    	telluric_par = hdrl_response_telluric_evaluation_parameter_create(telluric_models,
            x_corr_w_step, xcorr_half_win, xcorr_normalize, CPL_TRUE, quality_areas,
            fit_areas, lmin, lmax);

    hdrl_parameter * shift_par = NULL;
    if(velocity_enable)
    		shift_par = hdrl_spectrum1D_shift_fit_parameter_create(velocity_wguess,
    				velocity_range_wmin, velocity_range_wmax, velocity_fit_wmin,
					velocity_fit_wmax, velocity_fit_half_win);

    hdrl_parameter * resp_par =
    		hdrl_response_fit_parameter_create(11, fit_points, fit_wrange,
    				high_abs_regions);

    hdrl_response_result * response =
            hdrl_response_compute(obs_spectrum, ref_spectrum,
            		E_x,
					telluric_par,
            		shift_par,
					resp_calc_par ,
            		resp_par);

    hdrl_parameter_delete(telluric_par);
    hdrl_parameter_delete(shift_par);
    hdrl_parameter_delete(resp_calc_par);
    hdrl_parameter_delete(resp_par);

    if(response != NULL){
    	const hdrl_spectrum1D * fin_r = hdrl_response_result_get_final_response(response);
        cpl_table * tab_response_final =
                hdrl_spectrum1D_convert_to_table(
                		fin_r,
                        "response", "wavelength", "response_error", "response_bpm");

        cpl_table * tab_response_selected_points =
                hdrl_spectrum1D_convert_to_table(
                		hdrl_response_result_get_selected_response(response),
                        "response", "wavelength", "response_error", "response_bpm");


        cpl_table * tab_response_raw =
                hdrl_spectrum1D_convert_to_table(
                		hdrl_response_result_get_raw_response(response),
                        "response", "wavelength", "response_error", "response_bpm");


        cpl_propertylist * props = cpl_propertylist_new();

        save_tab(tab_response_final, parlist, props, frameset,
                spectrum_frame, "response.fits", "HDRLDEMO_RESPONSE");

        save_tab(tab_response_selected_points, parlist, props, frameset,
                spectrum_frame, "response_selected_points.fits",
                "HDRLDEMO_RESPONSE_SELECTED_POINTS");

        save_tab(tab_response_raw, parlist, props, frameset,
                spectrum_frame, "response_raw.fits", "HDRLDEMO_RESPONSE_RAW");

        cpl_table_delete(tab_response_final);
        cpl_table_delete(tab_response_selected_points);
        cpl_table_delete(tab_response_raw);

        cpl_propertylist_delete(props);
    }

    hdrl_response_result_delete(response);
    hdrl_spectrum1D_delete(&ref_spectrum);
    hdrl_spectrum1D_delete(&E_x);
    hdrl_spectrum1D_delete(&obs_spectrum);
    hdrl_spectrum1Dlist_delete(telluric_models);
    cpl_array_delete(fit_points);
    cpl_bivector_delete(quality_areas);
    cpl_bivector_delete(fit_areas);
    cpl_bivector_delete(high_abs_regions);
    cpl_frameset_delete(frameset_flux);
    cpl_frameset_delete(frameset_telluric_cat);
    cpl_frameset_delete(frameset_fit_areas);
    cpl_frameset_delete(frameset_quality_areas);
    cpl_frameset_delete(frameset_high_abs_regions);
    cpl_frameset_delete(frameset_interpolation_points);
    cpl_frameset_delete(frameset_extintion_cat);
    cpl_frameset_delete(frameset_std_flux_cat);

    return (int)cpl_error_get_code();
}

/**@}*/

