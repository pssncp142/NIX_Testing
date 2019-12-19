/*
 * This file is part of the HDRLDEMO pipeline
 * Copyright (C) 2012,2013 European Southern Observatory
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
 *
 *  Created on: Nov 30, 2012
 *      Author: agabasch
 */

#ifndef HDRLDEMO_DFS_H_
#define HDRLDEMO_DFS_H_

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <string.h>
#include <math.h>

CPL_BEGIN_DECLS

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* General raw-frame */
#define HDRLDEMO_RAW                        "RAW"
#define HDRLDEMO_RAW_CONFMAP                "RAW_CONFMAP"
#define HDRLDEMO_RAW_ERROR                  "RAW_ERROR"
#define HDRLDEMO_RAW_BPM                    "RAW_BPM"
#define HDRLDEMO_RAW_VARIANCE               "RAW_VARIANCE"

/* Specific raw-frame */
#define HDRLDEMO_BIAS                       "BIAS"
#define HDRLDEMO_DARK                       "DARK"
#define HDRLDEMO_FLAT                       "FLAT"
#define HDRLDEMO_FLAT_ON                    "FLAT_ON"
#define HDRLDEMO_FLAT_OFF                   "FLAT_OFF"
#define HDRLDEMO_SPECTRUM_1D                "SPECTRUM_1D"


/* Dedicated product calib-frame */
#define HDRLDEMO_MASTER_BPM                 "MASTER_BPM"
#define HDRLDEMO_MASTER_BPM_LIST            "MASTER_BPM_LIST"
#define HDRLDEMO_MASTER_BPM_FILTERED        "MASTER_BPM_FILTERED"
#define HDRLDEMO_MASTER_BPM_LIST_FILTERED   "MASTER_BPM_LIST_FILTERED"

#define HDRLDEMO_MASTER_LACOSMIC            "MASTER_LACOSMIC"
#define HDRLDEMO_MASTER_LACOSMIC_FILTERED   "MASTER_LACOSMIC_FILTERED"

#define HDRLDEMO_MASTER_BIAS                "MASTER_BIAS"
#define HDRLDEMO_MASTER_BIAS_ERROR          "MASTER_BIAS_ERROR"
#define HDRLDEMO_MASTER_BIAS_CONTRIBUTION   "MASTER_BIAS_CONTRIBUTION"

#define HDRLDEMO_MASTER_DARK                "MASTER_DARK"
#define HDRLDEMO_MASTER_DARK_ERROR          "MASTER_DARK_ERROR"
#define HDRLDEMO_MASTER_DARK_CONTRIBUTION   "MASTER_DARK_CONTRIBUTION"

#define HDRLDEMO_MASTER_FLAT                "MASTER_FLAT"
#define HDRLDEMO_MASTER_FLAT_ERROR          "MASTER_FLAT_ERROR"
#define HDRLDEMO_MASTER_FLAT_CONTRIBUTION   "MASTER_FLAT_CONTRIBUTION"
#define HDRLDEMO_MASTER_FLAT_SHAPE          "MASTER_FLAT_SHAPE"
#define HDRLDEMO_MASTER_FLAT_SHAPE_ERROR    "MASTER_FLAT_SHAPE_ERROR"

#define HDRLDEMO_MASTER_FRINGE              "MASTER_FRINGE"
#define HDRLDEMO_MASTER_FRINGE_ERROR        "MASTER_FRINGE_ERROR"
#define HDRLDEMO_MASTER_FRINGE_CONTRIBUTION "MASTER_FRINGE_CONTRIBUTION"
#define HDRLDEMO_MASTER_FRINGE_BPM          "MASTER_FRINGE_BPM"

#define HDRLDEMO_FRINGE_MASK                "FRINGE_MASK"
#define HDRLDEMO_FRINGE_CORRECTED           "FRINGE_CORRECTED"
#define HDRLDEMO_FRINGE_CORRECTED_ERROR     "FRINGE_CORRECTED_ERROR"
#define HDRLDEMO_FRINGE_CORRECTED_BPM       "FRINGE_CORRECTED_BPM"

#define HDRLDEMO_WAVELENGTHS_DEST_COLLAPSE  "SPECTRUM_COLLAPSE_WLENGTHS"

/* Efficiency and Response */
#define HDRLDEMO_ATM_EXT                    "ATM_EXT"
#define HDRLDEMO_FLUX_CATG                  "FLUX_CATG"

/* Response */
#define HDRLDEMO_TELL_CATG                  "TELLURIC_CATG"
#define HDRLDEMO_FIT_AREAS                  "FIT_AREAS"
#define HDRLDEMO_QUALITY_AREAS              "QUALITY_AREAS"
#define HDRLDEMO_HIGH_ABS_REGIONS			"HIGH_ABS_REGIONS"
#define HDRLDEMO_INT_POINTS_RESPONSE		"INTERPOL_POINTS_RESPONSE"

/* Dedicated static-frames */
#define HDRLDEMO_STATIC_MASK                "STATIC_MASK"
#define HDRLDEMO_OBJ_MASK                   "OBJ_MASK"
#define HDRLDEMO_POWERSPEC_MASK             "POWERSPEC_MASK"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_dfs_set_groups(cpl_frameset * set);

CPL_END_DECLS

#endif /* HDRLDEMO_DFS_H_ */
