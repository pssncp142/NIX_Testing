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
 * @defgroup hdrldemo_efficiency efficiency calculation
 * @par Synopsis: The recipe calculates the efficiency from observations of a
 * standard star.
 * @par Input frames:
 *
 *  DO category:                     Explanation:                 Required:
 *  SPECTRUM_1D                       Flux                        Yes
 *  ATM_EXT                           Atmospheric Extinction      Yes
 *  FLUX_CATG                         Std Stars Catalog           Yes
 *
 * @par Output frames:
 *
 *  DO category:                      Explanation:
 *  HDRLDEMO_EFFICIENCY               Efficiency
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*-----------------------------------------------------------------------------
                            Static global variables
 -----------------------------------------------------------------------------*/

#define RECIPE_NAME "hdrldemo_efficiency"

static hdrl_spectrum1D_wave_scale global_scale =
        hdrl_spectrum1D_wave_scale_linear;

static char hdrldemo_efficiency_description[] =
"                                                                           \n"
"The recipe derives efficiency from observations of the standard star.      \n"
"                                                                           \n"
"Input files:                                                               \n"
"                                                                           \n"
"  DO category:   Explanation:                              Required:       \n"
"  SPECTRUM_1D    Flux, wavelength in [nm]                  Yes             \n"
"  ATM_EXT        Atm. Extinction, wavelength in [nm]       Yes             \n"
"  FLUX_CATG      Std Stars Catalog, wavelength in [nm]     Yes             \n"
"                                                                           \n"
"Output files:                                                              \n"
"                                                                           \n"
"  DO category:                      Explanation:                           \n"
"  HDRLDEMO_EFFICIENCY               Efficiency                             \n"
"                                                                           \n"
"This recipe has polymorphic behaviour, it calculates the efficiency given  \n"
"the observed std star spectrum, the catalog of the models of the std stars.\n"
"The recipe selects the correct std star model from FLUX_CATG using the     \n"
"metadata of the observed std star spectrum.                                \n"
"                                                                           \n"
"  Usage of the recipe:                                                     \n"
"                                                                           \n"
"The recipe calculates the efficiency using the following formula:          \n"
"eff = (Iobs(l)*10^(0.4*(Ap-Am)*Ex(l))*G*Eph(l))/(Tex*Atel*Istd(l))         \n"
"where:                                                                     \n"
"Iobs(l): observed std star spectrum, provided by SPECTRUM_1D               \n"
"Istd(l): model of the std star spectrum, provided by FLUX_CATG             \n"
"Ex(l): model of the atmospheric extinction, selected from ATM_EXT          \n"
"Eph(l): energy of one photon, calculated internally in the routine         \n"
"Am: air mass, read from the FITS file header                               \n"
"Ap: air mass correction factor, hardcoded to 0.0                           \n"
"G: gain, read from the FITS file header                                    \n"
"Tex: exposure time, read from FITS file header                             \n"
"Atel: telescope collecting area, provided as recipe input parameter        \n"
"--area and --area-error are the parameters related to Atel                 \n";

/* Standard CPL recipe definition */
cpl_recipe_define(hdrldemo_efficiency,
                  HDRLDEMO_BINARY_VERSION,
                  "HDRL Group",
                  PACKAGE_BUGREPORT,
                  "2017",
                  "HDRLDEMO - Efficiency computation",
                  hdrldemo_efficiency_description);

/* Function needed by cpl_recipe_define to fill the input parameters */
static cpl_error_code hdrldemo_efficiency_fill_parameterlist(
                cpl_parameterlist   *   self) {

    cpl_parameter   *   par ;

    /*------------------------- Area telescope -------------------------*/
    /* --hdrldemo_efficiency area parameter */
    par = cpl_parameter_new_value(RECIPE_NAME".area", CPL_TYPE_DOUBLE,
            "Telescope effective collecting area in [cm2]", RECIPE_NAME, 51.2e4);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "area");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_efficiency error on area parameter */
    par = cpl_parameter_new_value(RECIPE_NAME".area-error", CPL_TYPE_DOUBLE,
        "Error on the Telescope effective collecting area in [cm2]",
        RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "area-error");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    /*------------------------- Area telescope -------------------------*/

    /*-------------------------- Exposure Time -------------------------*/
    /* --hdrldemo_efficiency exposure time */
    par = cpl_parameter_new_value(RECIPE_NAME".exposure-time", CPL_TYPE_DOUBLE,
            "Exposure time in [s]", RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "exposure-time");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_efficiency error on exposure time */
    par = cpl_parameter_new_value(RECIPE_NAME".exposure-time-error",
            CPL_TYPE_DOUBLE, "Exposure time error in [s]", RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "exposure-time-error");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    /*-------------------------- Exposure Time -------------------------*/

    /*----------------------------- Airmass ----------------------------*/
    /* --hdrldemo_efficiency airmass */
    par = cpl_parameter_new_value(RECIPE_NAME".airmass", CPL_TYPE_DOUBLE,
            "Airmass at which the standard star was observed", RECIPE_NAME, 1.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "airmass");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_efficiency error on airmass*/
    par = cpl_parameter_new_value(RECIPE_NAME".airmass-error", CPL_TYPE_DOUBLE,
         "Error on the airmass at which the standard star was observed",
         RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "airmass-error");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    /*----------------------------- Airmass ----------------------------*/

    /*----------------------- Airmass correction -----------------------*/
    /* --hdrldemo_efficiency airmass correction parameter*/
    par = cpl_parameter_new_value(RECIPE_NAME".airmass-correction",
            CPL_TYPE_DOUBLE, "Parameter to indicate if the efficiency is "
            "computed at arimass=0, or at a given non-zero value",
            RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "airmass-correction");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* --hdrldemo_efficiency error on airmass correction parameter*/
    par = cpl_parameter_new_value(RECIPE_NAME".airmass-correction-error",
            CPL_TYPE_DOUBLE, "Error on the airmass-correction parameter",
            RECIPE_NAME, 0.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI,
            "airmass-correction-error");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);
    /*----------------------- Airmass correction -----------------------*/

    /*------------------------------ Gain ------------------------------*/
   /* --hdrldemo_efficiency gain */
   par = cpl_parameter_new_value(RECIPE_NAME".gain", CPL_TYPE_DOUBLE,
           "Detector gain in [e/ADU]", RECIPE_NAME, 1.0);
   cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain");
   cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   cpl_parameterlist_append(self, par);

   /* --hdrldemo_efficiency gain*/
   par = cpl_parameter_new_value(RECIPE_NAME".gain-error", CPL_TYPE_DOUBLE,
        "Error on the detector gain in [e/ADU]", RECIPE_NAME, 0.0);
   cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "gain-error");
   cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
   cpl_parameterlist_append(self, par);
   /*------------------------------ Gain ------------------------------*/

    /* hdrldemo_efficiency window for noise calculation */
    par = cpl_parameter_new_value(RECIPE_NAME".half-window", CPL_TYPE_INT,
        "Half window used for noise calculation following the DER-SNR approach."
        "\nA good value is 2.5 x the resolving power of a spectral line.",
        RECIPE_NAME, 10);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "half-window");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    /* hdrldemo_efficiency sampling period */
    par = cpl_parameter_new_value(RECIPE_NAME".sampling-period", CPL_TYPE_DOUBLE,
            "Wavelength sampling period in [nm]", RECIPE_NAME, 1.0);
    cpl_parameter_set_alias(par, CPL_PARAMETER_MODE_CLI, "sampling-period");
    cpl_parameter_disable(par, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(self, par);

    return CPL_ERROR_NONE;
}

static double get_airmass(cpl_propertylist * header){

    const double airmass_start =
            cpl_propertylist_get_double(header, "ESO TEL AIRM START");
    const double airmass_end =
            cpl_propertylist_get_double(header, "ESO TEL AIRM END");

    return 0.5 * (airmass_start + airmass_end);
}

static void get_gain_exptime(const cpl_propertylist * header,
                double * gain, double * exptime){

    const char* instr_name = cpl_propertylist_get_string(header, "INSTRUME");

    const char * bname = !strcmp(instr_name, "XSHOOTER")
            ? "ESO SEQ ARM" : "ESO INS FILT1 NAME";

    const char * band = cpl_propertylist_get_string(header, bname);

    cpl_boolean need_hardcoded_data = !strcmp(instr_name, "SINFONI");
    need_hardcoded_data |= !strcmp(instr_name, "XSHOOTER")
            && !strcmp(band, "NIR");


    if (need_hardcoded_data){
        *gain = 2.12;
        *exptime = cpl_propertylist_get_double(header,"ESO DET DIT");
    }
    else{
        *gain = cpl_propertylist_get_double(header, "ESO DET OUT1 CONAD");
        *exptime = cpl_propertylist_get_double(header, "ESO DET WIN1 DIT1");
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

    *gain = override_default(parlist, gain_this, RECIPE_NAME".gain");
    *etime = override_default(parlist, etime_this, RECIPE_NAME".exposure-time");
    *amass = override_default(parlist, amass_this, RECIPE_NAME".airmass");

    cpl_propertylist_delete(header);

    return cpl_error_get_code();
}

static hdrl_spectrum1D * read_spectrum(const cpl_frame * obs_star,
        cpl_size ext, double sampling_period, cpl_size half_window){

    const char * fname = cpl_frame_get_filename(obs_star);
    cpl_propertylist * plist = cpl_propertylist_load(fname, ext);
    cpl_vector * vec_flux = cpl_vector_load(fname, ext);
    cpl_size sz = cpl_vector_get_size(vec_flux);

    cpl_image * flux = cpl_image_new(sz, 1, HDRL_TYPE_DATA);
    cpl_array * wavelengths = cpl_array_new(sz, HDRL_TYPE_DATA);

    const double crval1 = cpl_propertylist_get_double(plist,"CRVAL1");
    const double cdelt1 = cpl_propertylist_get_double(plist,"CDELT1");
    const double crpix1 = cpl_propertylist_get_double(plist,"CRPIX1");

    for(cpl_size i = 0; i < sz; ++i){
        const double fx = cpl_vector_get(vec_flux, i);
        cpl_image_set(flux, i + 1, 1, fx);

        const double lambda = sampling_period * (crval1 + cdelt1 * (i - crpix1));
        cpl_array_set(wavelengths, i, lambda);
    }

    cpl_propertylist_delete(plist);
    cpl_vector_delete(vec_flux);

    hdrl_spectrum1D * to_ret = hdrl_spectrum1D_create_error_DER_SNR
                              (flux, half_window, wavelengths, global_scale);

    cpl_image_delete(flux);
    cpl_array_delete(wavelengths);

    return to_ret;
}

static hdrl_spectrum1D * read_atm_ext(const cpl_frameset * frameset_extintion_cat){

    const char * flux_column = "EXTINCTION";
    const char * wavelength_column = "LAMBDA";

    const cpl_frame * atm_ext_frame =
                    cpl_frameset_get_position_const(frameset_extintion_cat, 0);
    const char * fname = cpl_frame_get_filename(atm_ext_frame);

    cpl_table * atm_ext_table = cpl_table_load(fname, 1, 0);

    hdrl_spectrum1D * sp = hdrl_spectrum1D_convert_from_table(atm_ext_table,
            flux_column, wavelength_column, NULL, NULL, global_scale);
    cpl_table_delete(atm_ext_table);

    return sp;
}

static hdrl_spectrum1D *
read_std_star_model(const cpl_frameset * frameset_std_flux_cat,
        const double ra,
        const double ra_tolerance,
        const double dec,
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

        const double ext_id =
                cpl_table_get_int(table_cat, col_name_extid, i ,&rej);
        const double curr_ra = cpl_table_get(table_cat, col_name_ra, i,&rej);
        const double curr_dec = cpl_table_get(table_cat, col_name_dec, i,&rej);

        if ((ext_id > 0) && (fabs(curr_ra - ra) < ra_tolerance)
                && (fabs(curr_dec - dec) < dec_tolerance)){
            model_extension = cpl_table_load(fname, ext_id, 0);
            break;
        }
    }

    hdrl_spectrum1D * spec = hdrl_spectrum1D_convert_from_table(model_extension,
            flux_column, wavelength_column, NULL, NULL, global_scale);

    cpl_table_delete(table_cat);
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

static void save_tab(const cpl_table * eff_table,
        const cpl_parameterlist * parlist, cpl_frameset * frameset,
        cpl_frame * raw_frame, const char * fname){

    cpl_propertylist * empty_list = cpl_propertylist_new();
    cpl_propertylist * applist = cpl_propertylist_new();

    cpl_propertylist_update_string(applist, CPL_DFS_PRO_CATG,
            "HDRLDEMO_EFFICIENCY");

    cpl_dfs_save_table(frameset,
            NULL,
            parlist,
            frameset,
            raw_frame,
            eff_table,
            empty_list,
            RECIPE_NAME,
            applist,
            NULL,
            PACKAGE "/" PACKAGE_VERSION,
            fname);

    cpl_propertylist_delete(applist);
    cpl_propertylist_delete(empty_list);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   perform master fringe combination
  @param   frameset   input set of frames
  @param   parlist    input recipe parameters
  @return  cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static int hdrldemo_efficiency(
                cpl_frameset            *   frameset,
                const cpl_parameterlist *   parlist)
{
    /* Check initial Entries */
    if (hdrldemo_check_and_set_groups(frameset) != CPL_ERROR_NONE) {
    	return cpl_error_get_code();
    }

    /*Conversion factor between nm and Angstrom*/
    const double nm2AA = 10.;

    const cpl_parameter *   par;
    const double ra_dec_tolerance = 0.0166667;

    /* Get parameters*/
    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".area");
    const double area = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".half-window");
    const cpl_size half_window = cpl_parameter_get_int(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".sampling-period");
    const double sampling_period = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".airmass-correction");
    const double Ap = cpl_parameter_get_double(par);

    /* ------------------------ Parameters errors ------------------------*/

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".airmass-error");
    const double Am_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist,
            RECIPE_NAME".airmass-correction-error");
    const double Ap_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".area-error");
    const double area_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist,
            RECIPE_NAME".exposure-time-error");
    const double Tex_e = cpl_parameter_get_double(par);

    par = cpl_parameterlist_find_const(parlist, RECIPE_NAME".gain-error");
    const double G_e = cpl_parameter_get_double(par);

    /* ------------------------ Parameters errors ------------------------*/

    cpl_size fsize = cpl_frameset_get_size(frameset);
    cpl_frameset * frameset_flux = cpl_frameset_new();
    cpl_frameset * frameset_std_flux_cat = cpl_frameset_new();
    cpl_frameset * frameset_extintion_cat = cpl_frameset_new();

    cpl_frame * spectrum_frame = NULL;
    for(cpl_size i = 0; i < fsize; ++i){
        cpl_frame * cur_frame = cpl_frameset_get_position(frameset, i);

        if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_SPECTRUM_1D)) {
            spectrum_frame = cpl_frame_duplicate(cur_frame);
            cpl_frameset_insert(frameset_flux,
                    spectrum_frame);
            cpl_msg_info(cpl_func,"Observed spectrum found in the SOF");
        } else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_ATM_EXT)) {
            cpl_frameset_insert(frameset_extintion_cat,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Atm. Extinction found in the SOF");
        }else if (!strcmp(cpl_frame_get_tag(cur_frame), HDRLDEMO_FLUX_CATG)) {
            cpl_frameset_insert(frameset_std_flux_cat,
                            cpl_frame_duplicate(cur_frame));
            cpl_msg_info(cpl_func,"Catalog of std stars found in the SOF");
        }else{
            cpl_msg_info(cpl_func, "SOF contains unknown tags");
        }
    }

    if(cpl_frameset_get_size(frameset_flux) != 1 ||
            cpl_frameset_get_size(frameset_std_flux_cat) != 1 ||
            cpl_frameset_get_size(frameset_extintion_cat) != 1 ){

        cpl_frameset_delete(frameset_flux);
        cpl_frameset_delete(frameset_std_flux_cat);
        cpl_frameset_delete(frameset_extintion_cat);

        return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                "SOF contains wrong number of frames");
    }


    double Am = 0.0;
    double Tex = 0.0;
    double G = 0.0;

    cpl_size n_ext = 0;

    cpl_error_code fail = get_data_from_header
            (frameset_flux, &n_ext, &Am, &G, &Tex, parlist);

    cpl_msg_debug(cpl_func, "Number of extensions: %lld, sampling period = %f",
            n_ext, sampling_period);

    cpl_msg_debug(cpl_func, "Airmass = %f +/- %f", Am, Am_e);
    cpl_msg_debug(cpl_func, "Airmass correction = %f +/- %f", Ap, Ap_e);
    cpl_msg_debug(cpl_func, "Gain = %f +/- %f", G, G_e);
    cpl_msg_debug(cpl_func, "Exposure = %f +/- %f", Tex, Tex_e);
    cpl_msg_debug(cpl_func, "Area = %f +/- %f", area, area_e);

    cpl_msg_debug(cpl_func, "DER_SNR, noise calculation, half-window: %lld",
            half_window);


    if(fail){

        cpl_frameset_delete(frameset_flux);
        cpl_frameset_delete(frameset_std_flux_cat);
        cpl_frameset_delete(frameset_extintion_cat);

        return fail;
    }

    hdrl_spectrum1D * atm_ext = read_atm_ext(frameset_extintion_cat);

    double ra = 0.0, dec = 0.0;

    get_ra_dec(frameset_flux, &ra, &dec);

    cpl_msg_debug(cpl_func, "RA: %f, DEC: %f", ra, dec);

    hdrl_spectrum1D * std_star_model =
            read_std_star_model(frameset_std_flux_cat, ra,
                    ra_dec_tolerance, dec, ra_dec_tolerance);

    const cpl_frame * obs_star =
                    cpl_frameset_get_position_const(frameset_flux, 0);

    cpl_table * eff_table = NULL;
    cpl_size neff_tot = 0;
    int ord = 0;
    for(cpl_size i = 0; i < n_ext; i+=3){

        hdrl_spectrum1D * obs_spectrum =
                read_spectrum(obs_star, i, sampling_period, half_window);

        const double w1 =
                hdrl_spectrum1D_get_wavelength_value(obs_spectrum, 1, NULL);
        const double w0 =
                hdrl_spectrum1D_get_wavelength_value(obs_spectrum, 0, NULL);

        /* Size of one bin (in nm).
         * NOTE: We assume that the input spectrum is sampled uniformly */
        const double delta = (w1 - w0);
        /* Size of one bin (in Angstrom)*/
        const double conversion_factor = delta * nm2AA;

        cpl_msg_debug(cpl_func, "Extension: %lld, conversion factor: %f", i,
                conversion_factor);
        /*
         * Convert observed flux to units per Angstrom, to have the same
         * unit of measure used in the models. This step is instrument-dependent.
        */
        hdrl_spectrum1D_div_scalar(obs_spectrum,
                (hdrl_value){conversion_factor, 0.0});

        hdrl_parameter * pars = hdrl_efficiency_parameter_create(
                (hdrl_value){(hdrl_data_t)Ap,   (hdrl_error_t)Ap_e},
                (hdrl_value){(hdrl_data_t)Am,   (hdrl_error_t)Am_e},
                (hdrl_value){(hdrl_data_t)G,    (hdrl_error_t)G_e},
                (hdrl_value){(hdrl_data_t)Tex,  (hdrl_error_t)Tex_e},
                (hdrl_value){(hdrl_data_t)area, (hdrl_error_t)area_e});

        hdrl_spectrum1D * eff = hdrl_efficiency_compute(
                obs_spectrum, std_star_model, atm_ext, pars);

        hdrl_parameter_delete(pars);

        if(eff == NULL){
            hdrl_spectrum1D_delete(&obs_spectrum);
            break;
        }

        cpl_table * mode = hdrl_spectrum1D_convert_to_table
                (eff, "efficiency","wavelength", "efficiency_error",
                        "efficiency_bpm");

        cpl_table_new_column(mode, "ord", CPL_TYPE_INT);
        cpl_table_fill_column_window(mode, "ord", 0,
                cpl_table_get_nrow(mode), ord);
        ord++;

        if(eff_table){
            neff_tot += cpl_table_get_nrow(eff_table);
            cpl_table_insert(eff_table, mode, neff_tot);
            cpl_table_delete(mode);
        }
        else{
            eff_table = mode;
            mode = NULL;
        }

        hdrl_spectrum1D_delete(&obs_spectrum);
        hdrl_spectrum1D_delete(&eff);
    }

    if(eff_table)
        save_tab(eff_table, parlist, frameset,
                spectrum_frame, "efficiency.fits");
    cpl_table_delete(eff_table);

    hdrl_spectrum1D_delete(&atm_ext);
    hdrl_spectrum1D_delete(&std_star_model);

    cpl_frameset_delete(frameset_flux);
    cpl_frameset_delete(frameset_std_flux_cat);
    cpl_frameset_delete(frameset_extintion_cat);
    return (int)cpl_error_get_code();
}

/**@}*/

