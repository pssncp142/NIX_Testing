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

#ifndef HDRLDEMO_UTILS_H_
#define HDRLDEMO_UTILS_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include "hdrl.h"

#include "hdrldemo_iterator.h"

#include <string.h>
#include <math.h>
#include <cpl.h>

CPL_BEGIN_DECLS

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/* Insert parameter in the input parameterlist */
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
    const char        *description);

/* */
cpl_error_code hdrldemo_check_and_set_groups(
    cpl_frameset            *frameset);

/* */
cpl_error_code hdrldemo_os_correct(
    cpl_frameset            *frameset,
    cpl_size                ext_num,
    const hdrl_parameter    *os_param,
    const hdrl_image        *masterbias,
    hdrl_imagelist          *out_hdrl_iml,
    hdrl_buffer             *buf);

/* */
cpl_frameset * hdrldemo_extract_frameset(
    const cpl_frameset      *in,
    const char              *tag);

/* */
cpl_error_code hdrldemo_hdrl_imagelist_load(
    cpl_frameset            *in,
    cpl_size                ext_num,
    cpl_frameset            *in_err,
    cpl_size                ext_num_err,
    cpl_frameset            *in_bpm,
    cpl_size                ext_num_bpm,
    hdrl_parameter          *region_params,
    double                  ron,
    double                  gain,
    hdrl_imagelist          *out_hdrl_iml,
    hdrl_buffer             *buf);

/* */
cpl_error_code hdrldemo_hdrl_image_load(
    cpl_frame               *in,
    cpl_size                ext_num,
    cpl_frame               *in_err,
    cpl_size                ext_num_err,
    cpl_frame               *in_bpm,
    cpl_size                ext_num_bpm,
    hdrl_parameter          *region_params,
    double                  ron,
    double                  gain,
    hdrl_image              **out_hdrl_ima);

/* */
cpl_error_code hdrldemo_save_image(
    cpl_propertylist        *header,
    const cpl_propertylist  *prevQClist,
    const char              *procatg,
    const char              *recipe,
    const char              *filename,
    const cpl_type          save_type,
    const cpl_image         *image,
    const cpl_parameterlist *parlist,
    cpl_frameset            *frameset);

/* Checks that the FITS image file exist */
cpl_boolean hdrldemo_check_file_exist(
    const char              *filename);

/* Checks if the file exists and an image could be loaded from the given ext */
cpl_boolean hdrldemo_check_image_exists(
    const char              *filename,
    cpl_size                extnum);

/* Checks that a QC parameter exists in the FITS file in the given extension. */
cpl_boolean hdrldemo_check_param(
    const char              *filename,
    cpl_size                extnum,
    const char              *keyname,
    cpl_boolean             check_value,
    cpl_type                type,
    int                     expected_int,
    double                  expected_double,
    const char              *expected_str);

/* */
cpl_error_code hdrldemo_detector_shotnoise_model(
    const cpl_image         *ima_data,
    const double            gain,
    const double            ron,
    cpl_image               **ima_errs);

/* */
cpl_error_code hdrldemo_detector_shotnoise_model_bias(
    const cpl_image         *ima_data,
    const double            ron,
    cpl_image               **ima_errs);

/* */
int hdrldemo_get_naxis1(
    const cpl_propertylist  *plist);

/* */
int hdrldemo_get_naxis2(
    const cpl_propertylist  *plist);

/* */
hdrl_value hdrldemo_utils_airmass(
    hdrl_value              airm1,
    hdrl_value              airm2,
    hdrl_value              ra,
    hdrl_value              dec,
    hdrl_value              lst,
    hdrl_value              exptime,
    hdrl_value              geolat,
    hdrl_airmass_approx     type);

CPL_END_DECLS

#endif /* HDRLDEMO_UTILS_H_ */
