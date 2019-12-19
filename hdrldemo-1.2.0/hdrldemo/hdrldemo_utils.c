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

#include "hdrldemo_utils.h"
#include "hdrldemo_dfs.h"
#include <cpl.h>

/*-----------------------------------------------------------------------------
                            Static Prototypes
 -----------------------------------------------------------------------------*/

static cpl_image * hdrldemo_mef_get(cpl_frame *, cpl_size, hdrl_parameter *);
static cpl_image * hdrldemo_mef_get_set(cpl_frameset *, cpl_size, cpl_size,
                                        hdrl_parameter *);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrldemo_utils  recipe related utilities
 *
 * TBD
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
 * @brief    Add parameter to the recipe parameterlist
 *
 * @param    recipe        Name of the recipe
 * @param    plist         Propertylist in the recipe for put the parameter
 * @param    param         Name  of the parameter
 * @param    alias         Alias of the parameter
 * @param    range         Flag: If the value is unique or it's a range
 * @param    rMin          In case of range: Min value
 * @param    rMax          In case of range: Max value
 * @param    type          Type of the parameter
 * @param    value         Generic pointer to the value/s of the parameter
 * @param    description   Description of the parameter
 *
 * @return   CPL_ERROR_NONE if everything is OK and CPL_ERROR_INVALID_TYPE if
 *           the type isn't correct.
 *
 * Description:
 *    Create a parameter with the input and insert in the recipe parameterlist.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_fill_parameter(
    const char*       recipe,
    cpl_parameterlist *plist,
    const char        *param,
    const char        *alias,
    cpl_boolean       range,
    const void        *rMin,
    const void        *rMax,
    cpl_type          type,
    const void        *value,
    const char        *description)
{
  cpl_parameter *p;

  if (range) {

      if(type == CPL_TYPE_INT) {

          const int val    = *(const int *)value;
          const int valMin = *(const int *)rMin;
          const int valMax = *(const int *)rMax;
          if(val < valMin || val > valMax) {
              return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Illegal int value range");
          }
          p = cpl_parameter_new_range(param, type, description, recipe, val, valMin, valMax);

      } else if( type == CPL_TYPE_DOUBLE) {

          const double val     = *(const double *)value;
          const double valMin  = *(const double *)rMin;
          const double valMax  = *(const double *)rMax;
          if(val < valMin || val > valMax) {
              return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Illegal double value range");
          }
          p = cpl_parameter_new_range(param, type, description, recipe, val, valMin, valMax);

      } else {
          return cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE, "cpl_type is not correct with ranges");
      }

  } else {

      switch(type) {
        case CPL_TYPE_STRING :
          p = cpl_parameter_new_value(param, type, description, recipe,  (const char        *)value); break;
        case CPL_TYPE_BOOL :
          p = cpl_parameter_new_value(param, type, description, recipe, *(const cpl_boolean *)value); break;
        case CPL_TYPE_INT :
          p = cpl_parameter_new_value(param, type, description, recipe, *(const int         *)value); break;
        case CPL_TYPE_DOUBLE :
          p = cpl_parameter_new_value(param, type, description, recipe, *(const double      *)value); break;
        default :
          return cpl_error_set_message(cpl_func, CPL_ERROR_INVALID_TYPE, "cpl_type is not correct");
      }
  }

  cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, alias);
  cpl_parameter_disable(  p, CPL_PARAMETER_MODE_ENV);

  cpl_parameterlist_append(plist, p);

  return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   check the entries in the recipe and classify the frameset with the tags
 *
 * @param   frameset      input set of frames
 *
 * @return  cpl_error_code
 */
/*----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_check_and_set_groups(
    cpl_frameset *frameset)
{
  /* Check size of frameset for to know if the sof file is not empty */
  cpl_size nframes = cpl_frameset_get_size(frameset);
  if(nframes <= 0){

      /* Empty sof file */
      return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                 "There aren't frames in the frameset");
  } else {
      for (cpl_size i = 0; i < nframes; i++) {

          cpl_frame  *frame    = cpl_frameset_get_position(frameset, i);
          const char *filename = cpl_frame_get_filename(frame);

          /* Check if the FITS file exist and have correct data,
           * return 0 if the fits file is valid without extensions */
          if (cpl_fits_count_extensions(filename) < 0){

              return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                         "Problem with the file '%s' (%s --> Code %d)",
                         filename, cpl_error_get_message(), cpl_error_get_code());
          }
      }
  }

  /* Identify the RAW, CONF and CALIB frames in the input frameset */
  if (hdrldemo_dfs_set_groups(frameset)) {

      /* Error classify frames */
      return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                 "Cannot classify RAW and/or CALIB frames");
  } else {

      /* Check classification */
      for (cpl_size i = 0; i < nframes; i++) {

          cpl_frame       *frame = cpl_frameset_get_position(frameset, i);
          const char      *tag   = cpl_frame_get_tag(frame);
          cpl_frame_group group  = cpl_frame_get_group(frame);

          /* The tag is invalid */
          if (group == CPL_FRAME_GROUP_NONE) {
              return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                         "Frame:%lld with tag:%s is invalid", i, tag);
          }
      }
  }

  return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    correct overscan
  @param    frameset        input frameset
  @param    ext_num         extension to load image from
  @param    os_param        overscan parameters
  @param    gain            input gain for error computation
  @param    ima_subtract    input image to be overscan corrected
  @param    ima_subtract_err error associated to ima_subtract
  @param    out_img_it      iterator on output image
  @param    out_err_it      iterator on error on output image
  @param    buf             buffer to allocate memory from or NULL
  @return   0 if everything is ok

  @doc
  Loop over input image frameset, compute and correct the overscan,
  creates corresponding error image eventually corrected by master bias error
  if master bias correction is used
  If an additional master bias is passed to this function, correct the overscan
  corrected image by subtracting this passed image and optionally propagate the
  error
 */
/*----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_os_correct(
        cpl_frameset            *   frameset,
        cpl_size                    ext_num,
        const hdrl_parameter    *   os_param,
        const hdrl_image        *   masterbias,
        hdrl_imagelist          *   out_hdrl_iml,
        hdrl_buffer             *   buf)
{

    cpl_ensure_code(frameset, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(out_hdrl_iml, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(os_param, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(
            hdrl_overscan_parameter_verify(os_param, -1, -1) == CPL_ERROR_NONE,
            CPL_ERROR_ILLEGAL_INPUT);
    cpl_size nx = -1;
    cpl_size ny = -1;
    double ccd_ron = hdrl_overscan_parameter_get_ccd_ron(os_param) ;

    for (cpl_size i = 0; i < cpl_frameset_get_size(frameset); i++) {
        cpl_image * img = hdrldemo_mef_get_set(frameset, i, ext_num, NULL);
        cpl_msg_info(cpl_func,"Doing Overscan correction ...");
        if( nx < 1 ) {
            nx = cpl_image_get_size_x(img);
            ny = cpl_image_get_size_y(img);
        }
        hdrl_overscan_compute_result * os_res;

        /* compute and correct overscan */
        os_res = hdrl_overscan_compute(img, os_param);

        if(!os_res) {
            cpl_image_delete(img);
            return cpl_error_get_code();
        }


        /* correct overscan
         * create an error image with zero error as the shot-noise pixel error
         * can only be calculated after the overscan correction */
        hdrl_overscan_correct_result * os_cor;
        hdrl_image* himg = hdrl_image_create(img, NULL);
        os_cor = hdrl_overscan_correct(himg, NULL, os_res);
        hdrl_image_delete(himg);

        if (masterbias == NULL) {
            cpl_image *ima_err;
            /* we derive a master bias ... */
            /* Compute the shot-noise error for the overscan-subtracted bias
             * frame itself and propagate it */

            cpl_msg_info(cpl_func, "Calculating the shot-noise error for the "
                         "overscan-subtracted bias and propagating it");

            /*  Please note that the shotnoise model function calculates the
             *  error on the full frame passed to the function, i.e. also on the
             *  overscan region, if the latter has not be trimmed! Therefore
             *  the next lines of code add this error to the overscan region as
             *  we dont trim here.
             *  This is not a basic problem as long as one remembers that the
             *  final propagated error in the overscan region itself is wrong */

            hdrl_image * tmp = hdrl_image_new(nx, ny);
            hdrldemo_detector_shotnoise_model_bias(hdrl_image_get_image(
                        hdrl_overscan_correct_result_get_corrected(os_cor)),
                                                   ccd_ron, &ima_err);
            cpl_image_add(hdrl_image_get_error(tmp), ima_err);
            cpl_image_delete(ima_err);
            hdrl_image_add_image(
                    hdrl_overscan_correct_result_get_corrected(os_cor), tmp);
            hdrl_image_delete(tmp);
        }

        /* If an additional master bias is passed to this function, correct
         * the overscan corrected image by subtracting this passed image and
         * propagate the error */

        if (masterbias != NULL && nx == hdrl_image_get_size_x(masterbias)
                        && ny == hdrl_image_get_size_y(masterbias)) {
            cpl_msg_info(cpl_func,"Subtracting masterbias");
            hdrl_image_sub_image(
             hdrl_overscan_correct_result_get_corrected(os_cor), masterbias);
        }

        if (buf) {
            hdrl_image * kimg = hdrl_image_new_from_buffer(nx, ny, buf);
            hdrl_image_copy(kimg,
                    hdrl_overscan_correct_result_get_corrected(os_cor), 1, 1);

            hdrl_imagelist_set(out_hdrl_iml, kimg, i);
        }
        else {
            hdrl_imagelist_set(out_hdrl_iml,
                       hdrl_overscan_correct_result_unset_corrected(os_cor), i);
        }

        /* cleanup */
        cpl_image_delete(img);

        hdrl_overscan_compute_result_delete(os_res);
        hdrl_overscan_correct_result_delete(os_cor);

        if (cpl_error_get_code())
            break;
    }

    return cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Extract the frames with the given tag from a frameset
   @param   in      A non-empty frameset
   @param   tag     The tag of the requested frames   
   @return  The newly created frameset or NULL on error

   The returned frameset must be de allocated with cpl_frameset_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_frameset * hdrldemo_extract_frameset(
        const cpl_frameset  *   in,
        const char          *   tag)
{
    cpl_frameset    *   out;
    cpl_frame       *   loc_frame;
    int                 nbframes;
    int                 i ;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (tag == NULL) return NULL ;

    /* Initialise */
    nbframes = cpl_frameset_get_size(in) ;

    /* Count the frames with the tag */
    if (cpl_frameset_count_tags(in, tag) == 0) return NULL ;

    /* Create the output frameset */
    out = cpl_frameset_new() ;

    /* Loop on the requested frames and store them in out */
    for (i=0 ; i<nbframes ; i++) {
    	const cpl_frame *cur_frame = cpl_frameset_get_position_const(in, i);
        if (!strcmp(cpl_frame_get_tag(cur_frame), tag)) {
            loc_frame = cpl_frame_duplicate(cur_frame) ;
            cpl_frameset_insert(out, loc_frame) ;
        }
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Load an hdrl_image from some frames. 
         If missing, an error is determined using a simple detector model.
  @param in             Input frame
  @param ext_num        Extension to load the image from
  @param in_err         Input frame containing errors (optional)
  @param ext_num_err    Extension to load the error from
  @param in_bpm         Input frame containing BPMs (optional)
  @param ext_num_bpm    Extension to load the BPM from
  @param region_params  Region to Load (if NULL, load all)
  @param ron            Read out noise
  @param gain           Detector gain
  @param out_hdrl_ima   Output hdrl_image object
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_hdrl_image_load(
        cpl_frame       *   in,
        cpl_size            ext_num,
        cpl_frame       *   in_err,
        cpl_size            ext_num_err,
        cpl_frame       *   in_bpm,
        cpl_size            ext_num_bpm,
        hdrl_parameter  *   region_params,
        double              ron,
        double              gain,
        hdrl_image      **  out_hdrl_ima)
{
    cpl_image       *   img ;

    /* Check Entries */
    cpl_ensure_code(in, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(out_hdrl_ima, CPL_ERROR_NULL_INPUT);

    /* Load in image */
    img = hdrldemo_mef_get(in, ext_num, region_params);
    if (img == NULL) return CPL_ERROR_ILLEGAL_INPUT ;

    /* Load BPM */
    if (in_bpm != NULL) {
        cpl_image * bpm = hdrldemo_mef_get(in_bpm, ext_num_bpm, region_params);
        cpl_mask * bpm_mask = cpl_mask_threshold_image_create(bpm,
                -0.5, 0.5) ;
        cpl_image_delete(bpm) ;
        cpl_mask_not(bpm_mask) ;
        cpl_image_reject_from_mask(img, bpm_mask);
        cpl_mask_delete(bpm_mask) ;
    }

    /* Load error */
    cpl_image * err;
    if (in_err != NULL) {
        err = hdrldemo_mef_get(in_err, ext_num_err, region_params);
    } else if (ron < 0. && gain < 0.) {
        /* no error */
        err = NULL;
    } else {
        /* No passed error -> Create uniform error */
        hdrldemo_detector_shotnoise_model(img, gain, ron, &err);
    }

    /* Create out himage */
    *out_hdrl_ima = hdrl_image_create(img, err);

    /* Cleanup */
    cpl_image_delete(img);
    cpl_image_delete(err);
    return cpl_error_get_code();
}
/*----------------------------------------------------------------------------*/
/**
  @brief Load an hdrl_imagelist from some framesets. 
         If missing, an error is determined using a simple detector model.
  @param in             Input frameset
  @param ext_num        Extension to load images from
  @param in_err         Input frameset containing errors (optional)
  @param ext_num_err    Extension to load errors from
  @param in_bpm         Input frameset containing BPMs (optional)
  @param ext_num_bpm    Extension to load BPMs from
  @param region_params  Region to Load (if NULL, load all)
  @param ron            Read out noise
  @param gain           Detector gain
  @param out_hdrl_iml   Output hdrl_imagelist object
  @param buf            Buffer
  @return   0 if everything is ok
  In case the buffer is uѕed, the region_params is ignored.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_hdrl_imagelist_load(
        cpl_frameset    *   in,
        cpl_size            ext_num,
        cpl_frameset    *   in_err,
        cpl_size            ext_num_err,
        cpl_frameset    *   in_bpm,
        cpl_size            ext_num_bpm,
        hdrl_parameter  *   region_params,
        double              ron,
        double              gain,
        hdrl_imagelist  *   out_hdrl_iml,
        hdrl_buffer     *   buf)
{
    cpl_size    nb_frames, i ;
    
    /* Check Entries */
    cpl_ensure_code(in, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(out_hdrl_iml, CPL_ERROR_NULL_INPUT);
    nb_frames = cpl_frameset_get_size(in) ;
    cpl_ensure_code(out_hdrl_iml, CPL_ERROR_NULL_INPUT);
    if (in_err != NULL) {
        cpl_ensure_code(cpl_frameset_get_size(in_err) == nb_frames, 
                CPL_ERROR_ILLEGAL_INPUT);
    }
    if (in_bpm != NULL) {
        cpl_ensure_code(cpl_frameset_get_size(in_bpm) == nb_frames, 
                CPL_ERROR_ILLEGAL_INPUT);
    }

    for (i = 0; i < nb_frames ; i++) {
        /* Load in image */
        cpl_image * img = hdrldemo_mef_get_set(in, i, ext_num, region_params);
        if (img == NULL) {
            break;
        }

        /* Load BPM */
        if (in_bpm != NULL) {
            cpl_image * bpm = hdrldemo_mef_get_set(in_bpm, i, ext_num_bpm,
                    region_params);
            if (bpm == NULL) {
                cpl_image_delete(img);
                break;
            }
            cpl_mask * bpm_mask = cpl_mask_threshold_image_create(bpm,
                    -0.5, 0.5) ;
            cpl_image_delete(bpm) ;
            cpl_mask_not(bpm_mask) ;
            cpl_image_reject_from_mask(img, bpm_mask);
            cpl_mask_delete(bpm_mask) ;
        }

        /* Load error */
        cpl_image * err;
        if (in_err != NULL) {
            err = hdrldemo_mef_get_set(in_err, i, ext_num_err, region_params);
        } else {
            /* No passed error -> Create uniform error */
            hdrldemo_detector_shotnoise_model(img, gain, ron, &err);
        }
        if (err == NULL) {
            cpl_image_delete(img);
            break;
        }

        /* Handle buffer */
        hdrl_image * hdrl_out;
        if (buf) {
            hdrl_out = hdrl_image_new_from_buffer(
                    cpl_image_get_size_x(img), cpl_image_get_size_y(img), buf);
        } else {
            hdrl_out = hdrl_image_create(img, err);
        }
        hdrl_image_insert(hdrl_out, img, err, 1, 1);
        hdrl_imagelist_set(out_hdrl_iml, hdrl_out, i);

        /* Cleanup */
        cpl_image_delete(img);
        cpl_image_delete(err);

        if (cpl_error_get_code()) break;
    }
    return cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   save an image on a FITS file with proper FITS header
 *
 * @param   header      property list with parameters added to the header
 * @param   prevQClist  property list with QC parameters added previously
 * @param   procatg     FITS product category
 * @param   recipe      recipe name
 * @param   filename    product filename
 * @param   save_type   cpl image type
 * @param   image       actual image to be saved
 * @param   parlist     input recipe parameters
 * @param   frameset    product set of frames
 *
 * @return  cpl_error_code
 */
/*----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_save_image( cpl_propertylist 	    *header,
									const cpl_propertylist  *prevQClist,
		                            const char              *procatg,
									const char              *recipe,
									const char              *filename,
									const cpl_type          save_type,
									const cpl_image         *image,
									const cpl_parameterlist *parlist,
									cpl_frameset            *frameset)
{
    /* Property list of QC parameter */
	cpl_propertylist *qclist;
	if (prevQClist) {
		qclist = cpl_propertylist_duplicate(prevQClist);
	} else {
		qclist = cpl_propertylist_new();
	}

    /* Add the product category and save image */
    cpl_propertylist_update_string(qclist, CPL_DFS_PRO_CATG, procatg);

    cpl_dfs_save_image( frameset, header, parlist, frameset, NULL,
						image, save_type, recipe, qclist, NULL,
						PACKAGE "/" PACKAGE_VERSION, filename);

    cpl_propertylist_delete(qclist);

    return cpl_error_get_code();
}


/*----------------------------------------------------------------------------*/
/**
 * @brief   Checks that the FITS image file exist
 *
 * @param   filename      product filename
 *
 * @return  cpl_boolean   Show if the file exist
 */
/*----------------------------------------------------------------------------*/
cpl_boolean hdrldemo_check_file_exist(
	const char *filename)
{
    cpl_image *image = cpl_image_load(filename, CPL_TYPE_FLOAT, 0, 0);

    if (image == NULL) {
        cpl_test_error(CPL_ERROR_FILE_IO);
        return CPL_FALSE;
    } else {
    	cpl_image_delete(image);
    	return CPL_TRUE;
    }
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Checks if the file exists and an image could be loaded from
 *          the given ext
 *
 * @param   filename      product filename
 * @param   extnum        extension number that should contain the image
 *
 * @return  cpl_boolean   Show if the file exist and the extension have an image
 */
/*----------------------------------------------------------------------------*/
cpl_boolean hdrldemo_check_image_exists(
	const char *filename,
	cpl_size   extnum)
{
    cpl_image *image = cpl_image_load(filename, CPL_TYPE_FLOAT, 0, extnum);

    if (image == NULL) {
        cpl_test_error(CPL_ERROR_FILE_IO);
        return CPL_FALSE;
    } else {
    	cpl_image_delete(image);
    	return CPL_TRUE;
    }
}

/*  */
/*----------------------------------------------------------------------------*/
/**
 * @brief   Checks that a QC parameter exists in the FITS file
 *          in the given extension and if it's mark check it value.
 *
 * @param   filename          product filename
 * @param   extnum            extension number that should contain the image
 * @param   keyname           String with the name of the QC parameter
 * @param   check_value       If it's necessary to check it value
 * @param   type              Type of the QC parameter (INT, DOUBLE, STRING)
 * @param   expected_int      Expected value INT
 * @param   expected_double   Expected value DOUBLE
 * @param   expected_str      Expected value STRING
 *
 * @return  cpl_boolean       Show if the QC exist and test the value
 */
/*----------------------------------------------------------------------------*/
cpl_boolean hdrldemo_check_param(
	const char  *filename,
	cpl_size    extnum,
	const char  *keyname,
	cpl_boolean check_value,
	cpl_type    type,
	int         expected_int,
	double      expected_double,
	const char  *expected_str)
{
	cpl_boolean result = CPL_FALSE;
	cpl_propertylist *qclist = cpl_propertylist_load(filename, extnum);
    if (qclist != NULL) {

    	if(cpl_propertylist_has(qclist, keyname)){

    		if(check_value == CPL_FALSE){
    			result = CPL_TRUE;
    		} else {
				switch (type) {
				case CPL_TYPE_INT :

					if(expected_int == cpl_propertylist_get_int(
												qclist, keyname)){
						result = CPL_TRUE;
					}
					break;
				case CPL_TYPE_DOUBLE :
					if(expected_double == cpl_propertylist_get_double(
												qclist, keyname)){
						result = CPL_TRUE;
					}
					break;
				case CPL_TYPE_STRING :
					if(!strcmp(expected_str, cpl_propertylist_get_string(
												qclist, keyname))){
						result = CPL_TRUE;
					}
					break;
				default :
					break;
				}
    		}
    	}

        cpl_propertylist_delete(qclist);
    }
    return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   compute photon count error in [ADU]
  @param   ima_data in [ADU]
  @param   gain detector's gain in [e- / ADU]
  @param   ron  detector's read out noise in [ADU]
  @param   ima_errs output error image in [ADU]
  @return  cpl_error_code
  @note ima_errs need to be deallocated
        ima_data must contain the photon counts with no offsets
        this usually means the image must be overscan and bias corrected
        Then the shot noise can be calculated from the poissonian distribution
        as sqrt(electron-counts). To this (transformed back into ADUs) the
        readout noise is added in quadrature.
  @doc
  error is computed with standard formula

  \f$ err_{ADU} = \sqrt{ \frac{ counts }{ gain } + ron^{ 2 } } \f$

  If an image value is negative the associated error is set to RON
 */
/*----------------------------------------------------------------------------*/
cpl_error_code
hdrldemo_detector_shotnoise_model(const cpl_image* ima_data, const double gain,
                              const double ron, cpl_image ** ima_errs)
{
    cpl_ensure_code(ima_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(ima_errs, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(gain > 0., CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(ron > 0., CPL_ERROR_ILLEGAL_INPUT);

    *ima_errs = cpl_image_duplicate(ima_data);
    /* set negative values (= zero measurable electrons) to read out noise */
    cpl_image_threshold(*ima_errs, 0., INFINITY, ron, ron);

    /* err_ADU = sqrt(counts/gain + ron * ron)*/

    cpl_image_divide_scalar(*ima_errs, gain);
    cpl_image_add_scalar(*ima_errs, ron * ron);
    cpl_image_power(*ima_errs, 0.5);

    return cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief   compute error in [ADU]
  @param   ima_data in [ADU]
  @param   ron  detector's read out noise in [ADU]
  @param   ima_errs output error image in [ADU]
  @return  cpl_error_code
  @note ima_errs need to be deallocated
        The returned error image is filled by the RON in ADUs
 */
/*----------------------------------------------------------------------------*/
cpl_error_code
hdrldemo_detector_shotnoise_model_bias(const cpl_image* ima_data, const double ron,
                                   cpl_image ** ima_errs)
{
    cpl_ensure_code(ima_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(ima_errs, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(ron > 0., CPL_ERROR_ILLEGAL_INPUT);

    /* Create new image initialised to 0 */
    *ima_errs = cpl_image_new(cpl_image_get_size_x(ima_data),
                              cpl_image_get_size_y(ima_data), CPL_TYPE_DOUBLE);

    if (cpl_image_get_bpm_const(ima_data))
        cpl_image_reject_from_mask(*ima_errs,
                                   cpl_image_get_bpm_const(ima_data));

    cpl_image_add_scalar(*ima_errs, ron);

    return cpl_error_get_code();
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief return NAXIS1 entry of plist, uses ZNAXIS1 if present
 *
 * @param plist		Input property list
 *
 * @return 			value of ZNAXIS1 or NAXIS1 stored in plist
 *
 */
/* ---------------------------------------------------------------------------*/
int hdrldemo_get_naxis1(const cpl_propertylist *plist)
{
    if (cpl_propertylist_has(plist, "ZNAXIS1")) {
        return cpl_propertylist_get_int(plist, "ZNAXIS1");
    } else {
        return cpl_propertylist_get_int(plist, "NAXIS1");
    }
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief return NAXIS2 entry of plist, uses ZNAXIS2 if present
 *
 * @param plist		Input property list
 *
 * @return 			value of ZNAXIS2 or NAXIS2 stored in plist
 *
 */
/* ---------------------------------------------------------------------------*/
int hdrldemo_get_naxis2(const cpl_propertylist *plist)
{
    if (cpl_propertylist_has(plist, "ZNAXIS2")) {
        return cpl_propertylist_get_int(plist, "ZNAXIS2");
    } else {
        return cpl_propertylist_get_int(plist, "NAXIS2");
    }
}

/*----------------------------------------------------------------------------*/
/**
 * @brief   Compute and check the effective airmass of an observation.
 * Takes in count the error propagation if you enter the relative error of
 * the input parameters in a hdrl_value structure {data,error}
 *
 * @param airm1    airmass value at start of the observation
 * @param airm2    airmass value at end   of the observation
 * @param ra       right ascension in degrees
 * @param dec      declination     in degrees
 * @param lst      local sideral time seconds
 * @param exptime  integration   time seconds
 * @param geolat   latitude of observatory site in degrees
 * @param type     Kind of approximation
 *
 * @return Computed airmass or (hdrl_value){-1.,0.} on error. The function gives
 *  a warning if the output is not in between the airmass input start-end values
 *   (with a fuzzyness of 0.005).
 *
 **/
/*----------------------------------------------------------------------------*/
hdrl_value hdrldemo_utils_airmass(hdrl_value airm1, hdrl_value airm2,
		hdrl_value ra, hdrl_value dec, hdrl_value lst, hdrl_value exptime,
		hdrl_value geolat, hdrl_airmass_approx type)
{
	hdrl_value err = {-1.,0.};

	cpl_ensure(airm1.data > 0. && airm2.data > 0., CPL_ERROR_ILLEGAL_INPUT, err);

	hdrl_value airmass = hdrl_utils_airmass(ra, dec, lst, exptime, geolat, type);
	if (airmass.data < 0.) {

		/* airmass computation failed, use simple average and airmass.error */
		hdrl_value average = { (airm1.data  + airm2.data) / 2.,
				airm1.error + airm2.error};

		/* comparison doesn't make sense in this case */
		cpl_msg_info(cpl_func, "Airmass computation unsuccessful (%s), using "
				"simple average of start and end values (%f) and error "
				"of entry airmass",	cpl_error_get_message(), average.data);

		return average;

	} else {

		/* Check range of value */
		cpl_boolean check =    airmass.data > fmin(airm1.data, airm2.data) - 0.005
				&& airmass.data < fmax(airm1.data, airm2.data) + 0.005;
		if (!check) {
			cpl_msg_warning(cpl_func, "Computed airmass=%.3f+-[%f] is NOT "
					"in the range recorded in the input [star: %f, end: %f]",
					airmass.data, airmass.error, airm1.data, airm2.data);
		}
	}

	return airmass;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief load header of frame and wrap negative region parameters
 *
 * @param frm    frame to load from
 * @param extnum extension number to load from
 * @param par    region parameter whose ranges should be wrapped by the
 *               image size
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code
hdrldemo_fix_neg_region(const cpl_frame * frm, cpl_size extnum,
                        hdrl_parameter * par)
{
    cpl_propertylist * plist =
        cpl_propertylist_load(cpl_frame_get_filename(frm), extnum);
    if (!plist) {
        return cpl_error_get_code();
    }
    cpl_size nx = hdrldemo_get_naxis1(plist);
    cpl_size ny = hdrldemo_get_naxis2(plist);
    hdrl_rect_region_fix_negatives(par, nx, ny);
    cpl_propertylist_delete(plist);
    return cpl_error_get_code();
}


/**@}*/

static cpl_image * hdrldemo_mef_get(
        cpl_frame       *   frm,
        cpl_size            eidx,
        hdrl_parameter  *   region_params)
{
    cpl_image   *   out = NULL ;

    if (region_params == NULL) {
        out = cpl_image_load(cpl_frame_get_filename(frm),
                CPL_TYPE_DOUBLE, 0, eidx);
    } else {
        hdrldemo_fix_neg_region(frm, eidx, region_params);
        out = cpl_image_load_window(cpl_frame_get_filename(frm),
                CPL_TYPE_DOUBLE, 0, eidx,
                hdrl_rect_region_get_llx(region_params),
                hdrl_rect_region_get_lly(region_params),
                hdrl_rect_region_get_urx(region_params),
                hdrl_rect_region_get_ury(region_params)) ;
    }
    return out ;
}

static cpl_image * hdrldemo_mef_get_set(
        cpl_frameset    *   set,
        cpl_size            fidx,
        cpl_size            eidx,
        hdrl_parameter  *   region_params)
{
    cpl_frame   *   cur_frame = cpl_frameset_get_position(set, fidx);

    return hdrldemo_mef_get(cur_frame, eidx, region_params);
}


