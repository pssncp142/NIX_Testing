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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                                                Includes
 -----------------------------------------------------------------------------*/

#include "hdrldemo_dfs.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrldemo_dfs  DFS related functions
 *
 * TBD
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
 * @brief    Set the group as RAW or CALIB in a frameset
 * @param    set     the input frameset
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code hdrldemo_dfs_set_groups(cpl_frameset * set)
{
    /* Check entries */
    cpl_ensure_code(set != NULL, CPL_ERROR_NULL_INPUT);

    /* Initialize */
    cpl_size nframes = cpl_frameset_get_size(set);

    /* Loop on frames */
    for (cpl_size i = 0; i < nframes; i++) {

    	cpl_frame  *cur_frame = cpl_frameset_get_position(set, i);
    	const char *tag       = cpl_frame_get_tag(cur_frame);

        if (       !strcmp(tag, HDRLDEMO_RAW                       ) ||
                   !strcmp(tag, HDRLDEMO_RAW_CONFMAP               ) ||
                   !strcmp(tag, HDRLDEMO_RAW_ERROR                 ) ||
                   !strcmp(tag, HDRLDEMO_RAW_BPM                   ) ||
                   !strcmp(tag, HDRLDEMO_RAW_VARIANCE              ) ||
                   !strcmp(tag, HDRLDEMO_BIAS                      ) ||
                   !strcmp(tag, HDRLDEMO_DARK                      ) ||
                   !strcmp(tag, HDRLDEMO_FLAT                      ) ||
                   !strcmp(tag, HDRLDEMO_FLAT_ON                   ) ||
                   !strcmp(tag, HDRLDEMO_FLAT_OFF                  ) ||
                   !strcmp(tag, HDRLDEMO_SPECTRUM_1D               ) ){

        	/* raw-frames */
            cpl_frame_set_group(cur_frame, CPL_FRAME_GROUP_RAW);

        } else if (!strcmp(tag, HDRLDEMO_MASTER_BPM                ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_BPM_LIST           ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_BPM_FILTERED       ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_BPM_LIST_FILTERED  ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_LACOSMIC           ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_LACOSMIC_FILTERED  ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_BIAS               ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_BIAS_ERROR         ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_BIAS_CONTRIBUTION  ) ||
				   !strcmp(tag, HDRLDEMO_MASTER_DARK               ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_DARK_ERROR         ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_DARK_CONTRIBUTION  ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FLAT               ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FLAT_ERROR         ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FLAT_CONTRIBUTION  ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FLAT_SHAPE         ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FLAT_SHAPE_ERROR   ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FRINGE             ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FRINGE_ERROR       ) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FRINGE_CONTRIBUTION) ||
                   !strcmp(tag, HDRLDEMO_MASTER_FRINGE_BPM         ) ||
				   !strcmp(tag, HDRLDEMO_FRINGE_MASK               ) ||
                   !strcmp(tag, HDRLDEMO_FRINGE_CORRECTED          ) ||
                   !strcmp(tag, HDRLDEMO_FRINGE_CORRECTED_ERROR    ) ||
				   !strcmp(tag, HDRLDEMO_FRINGE_CORRECTED_BPM      ) ||
				   !strcmp(tag, HDRLDEMO_WAVELENGTHS_DEST_COLLAPSE ) ||
                   !strcmp(tag, HDRLDEMO_ATM_EXT                   ) ||
                   !strcmp(tag, HDRLDEMO_FLUX_CATG                 ) ||
                   !strcmp(tag, HDRLDEMO_TELL_CATG                 ) ||
                   !strcmp(tag, HDRLDEMO_FIT_AREAS                 ) ||
                   !strcmp(tag, HDRLDEMO_QUALITY_AREAS             ) ||
                   !strcmp(tag, HDRLDEMO_HIGH_ABS_REGIONS          ) ||
                   !strcmp(tag, HDRLDEMO_INT_POINTS_RESPONSE       ) ||
                   !strcmp(tag, HDRLDEMO_STATIC_MASK               ) ||
                   !strcmp(tag, HDRLDEMO_POWERSPEC_MASK            ) ||
				   !strcmp(tag, HDRLDEMO_OBJ_MASK                  ) ){

        	/* calib-frame */
            cpl_frame_set_group(cur_frame, CPL_FRAME_GROUP_CALIB);

        } else {

        	/* unknown-frame */
            cpl_frame_set_group(cur_frame, CPL_FRAME_GROUP_NONE);
            cpl_msg_warning(cpl_func, "Frame:%lld with tag:%s, unknown!", i, tag);
        }
    }

    return cpl_error_get_code();
}

/**@}*/


