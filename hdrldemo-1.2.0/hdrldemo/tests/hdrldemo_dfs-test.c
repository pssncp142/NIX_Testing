/*
 * This file is part of the HDRLDEMO Toolkit.
 * Copyright (C) 2016 European Southern Observatory
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
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/*-----------------------------------------------------------------------------
                               Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "hdrldemo_dfs.h"

/*-----------------------------------------------------------------------------
                               Defines
 -----------------------------------------------------------------------------*/

/**
 * @brief Macro used to tests that all the frames in the frameset @a frames have
 * the group value indicated by @a result.
 */
/* NOTE: we use a macro rather than a function so that cpl_test produces more
   understandable output in the logs.
 */
#define test_frame_groups_eq(frames, result)                                   \
    for (cpl_size i = 0; i < cpl_frameset_get_size(frames); ++i) {             \
        const cpl_frame * frame = cpl_frameset_get_position_const(frames, i);  \
        cpl_test_eq(cpl_frame_get_group(frame), result);                       \
    }

/*-----------------------------------------------------------------------------
                               Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_frameset * setup_test_frameset(const char **tags);

/*-----------------------------------------------------------------------------
                               Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* List of tags that are not known by HDRLDEMO. */
    const char * unknown_tags[] = {
		"SOMEDATA",
		"DATA_FRAME",
		NULL
	};

    /* List of raw frame tags for HDRLDEMO. */
    const char *raw_tags[] = {
        HDRLDEMO_RAW,
        HDRLDEMO_RAW_CONFMAP,
        HDRLDEMO_RAW_ERROR,
        HDRLDEMO_RAW_BPM,
        HDRLDEMO_RAW_VARIANCE,
        HDRLDEMO_BIAS,
        HDRLDEMO_DARK,
        HDRLDEMO_FLAT,
        HDRLDEMO_FLAT_ON,
        HDRLDEMO_FLAT_OFF,
        HDRLDEMO_SPECTRUM_1D,
		NULL
	};

    /* List of calibration frame tags for HDRLDEMO. */
    const char *calib_tags[] = {
		HDRLDEMO_MASTER_BPM,
		HDRLDEMO_MASTER_BPM_LIST,
		HDRLDEMO_MASTER_BPM_FILTERED,
		HDRLDEMO_MASTER_BPM_LIST_FILTERED,
		HDRLDEMO_MASTER_LACOSMIC,
		HDRLDEMO_MASTER_LACOSMIC_FILTERED,
		HDRLDEMO_MASTER_BIAS,
		HDRLDEMO_MASTER_BIAS_ERROR,
		HDRLDEMO_MASTER_BIAS_CONTRIBUTION,
		HDRLDEMO_MASTER_DARK,
		HDRLDEMO_MASTER_DARK_ERROR,
		HDRLDEMO_MASTER_DARK_CONTRIBUTION,
		HDRLDEMO_MASTER_FLAT,
		HDRLDEMO_MASTER_FLAT_ERROR,
		HDRLDEMO_MASTER_FLAT_CONTRIBUTION,
		HDRLDEMO_MASTER_FLAT_SHAPE,
		HDRLDEMO_MASTER_FLAT_SHAPE_ERROR,
		HDRLDEMO_MASTER_FRINGE,
		HDRLDEMO_MASTER_FRINGE_ERROR,
		HDRLDEMO_MASTER_FRINGE_CONTRIBUTION,
		HDRLDEMO_MASTER_FRINGE_BPM,
		HDRLDEMO_FRINGE_MASK,
		HDRLDEMO_FRINGE_CORRECTED,
		HDRLDEMO_FRINGE_CORRECTED_ERROR,
		HDRLDEMO_FRINGE_CORRECTED_BPM,
		HDRLDEMO_WAVELENGTHS_DEST_COLLAPSE,
		HDRLDEMO_ATM_EXT,
		HDRLDEMO_FLUX_CATG,
		HDRLDEMO_TELL_CATG,
		HDRLDEMO_FIT_AREAS,
		HDRLDEMO_QUALITY_AREAS,
		HDRLDEMO_HIGH_ABS_REGIONS,
		HDRLDEMO_INT_POINTS_RESPONSE,
		HDRLDEMO_STATIC_MASK,
		HDRLDEMO_OBJ_MASK,
		NULL
	};


    /******* TESTS *******/
    cpl_frameset *frameset;

    /* Test error handling when hdrldemo_dfs_set_groups is given a NULL pointer. */
    cpl_error_code result = hdrldemo_dfs_set_groups(NULL);
    cpl_test_eq_error(result, CPL_ERROR_NULL_INPUT);

    /* Test that hdrldemo_dfs_set_groups does not fail with an empty frameset. */
    const char *no_tags[] = {NULL};
    frameset = setup_test_frameset(no_tags);
    result = hdrldemo_dfs_set_groups(frameset);
    cpl_test_eq_error(result, CPL_ERROR_NONE);
    cpl_frameset_delete(frameset);

    /* Test hdrldemo_dfs_set_groups marks a single single frame correctly. */
    const char * single_tag[] = {HDRLDEMO_RAW, NULL};
    frameset = setup_test_frameset(single_tag);
    result = hdrldemo_dfs_set_groups(frameset);
    cpl_test_eq_error(result, CPL_ERROR_NONE);
    test_frame_groups_eq(frameset, CPL_FRAME_GROUP_RAW);
    cpl_frameset_delete(frameset);

    /* Test that hdrldemo_dfs_set_groups leaves the unknown frames unmarked. */
    frameset = setup_test_frameset(unknown_tags);
    result = hdrldemo_dfs_set_groups(frameset);
    cpl_test_eq_error(result, CPL_ERROR_NONE);
    test_frame_groups_eq(frameset, CPL_FRAME_GROUP_NONE);
    cpl_frameset_delete(frameset);

    /* Test that hdrldemo_dfs_set_groups marks the raw frames correctly. */
    frameset = setup_test_frameset(raw_tags);
    result = hdrldemo_dfs_set_groups(frameset);
    cpl_test_eq_error(result, CPL_ERROR_NONE);
    test_frame_groups_eq(frameset, CPL_FRAME_GROUP_RAW);
    cpl_frameset_delete(frameset);

    /* Test that hdrldemo_dfs_set_groups marks the calibration frames correctly. */
    frameset = setup_test_frameset(calib_tags);
    result = hdrldemo_dfs_set_groups(frameset);
    cpl_test_eq_error(result, CPL_ERROR_NONE);
    test_frame_groups_eq(frameset, CPL_FRAME_GROUP_CALIB);
    cpl_frameset_delete(frameset);

    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                               Private function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @brief    Function to setup a new frameset with frames that have tags as
 *           indicated by the NULL terminated list in 'tags'.
 *           Note: the frame's group values are set to CPL_FRAME_GROUP_NONE.
 *
 * @param    tags     the input tag arry with the name of new frames
 *
 * @return   The frameset with the tag names in each frame
 *
 **/
/*----------------------------------------------------------------------------*/
static cpl_frameset * setup_test_frameset(const char **tags)
{
    cpl_frameset *frameset = cpl_frameset_new();

    for (const char **ptag = tags; *ptag != NULL; ++ptag) {

        /* Prepare a dummy input frame. */
     	cpl_frame *frame = cpl_frame_new();
     	cpl_frame_set_tag(frame, *ptag);

        /* Add to the set */
        cpl_frameset_insert(frameset, frame);
    }

    return frameset;
}

