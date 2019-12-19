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


#include "hdrldemo_utils.h"
#include "hdrldemo_iterator.h"
#include "hdrl_iter.h"

/*-----------------------------------------------------------------------------
                               Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                               Private function prototypes
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                               Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Create dummy images and save to disk */
    const char *filename1 = "iter1_dummy_raw.fits";
    const char *filename2 = "iter2_dummy_raw.fits";
    cpl_image  *image1    = cpl_image_new(10, 10, CPL_TYPE_FLOAT);
    cpl_image  *image2    = cpl_image_new(10, 10, CPL_TYPE_FLOAT);
    cpl_image_save(image1, filename1, CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    cpl_image_save(image2, filename2, CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    cpl_test(hdrldemo_check_file_exist(filename1));
    cpl_test(hdrldemo_check_file_exist(filename2));
    cpl_test(hdrldemo_check_image_exists(filename1, 0));
    cpl_test(hdrldemo_check_image_exists(filename2, 0));

    /* Crate a imagelist and save to disk */
    const char *filename = "cube_dummy_raw.fits";
    cpl_imagelist *imglist = cpl_imagelist_new();
    cpl_imagelist_set(imglist, image1, 0);
    cpl_imagelist_set(imglist, image2, 1);
    cpl_imagelist_save(imglist, filename, CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    cpl_test(hdrldemo_check_file_exist(filename));

    /* Create parameter list */
    cpl_parameterlist *parlist = cpl_parameterlist_new();

    /* Create frameset */
    cpl_frameset *frameset = cpl_frameset_new();

    /* Prepare a dummy input raw frame. */
	cpl_frame *frame1 = cpl_frame_new();
	cpl_frame_set_filename( frame1, filename1);
	cpl_frame_set_tag(      frame1, "RAW");
	cpl_frame_set_type(     frame1, CPL_FRAME_TYPE_IMAGE);
	cpl_frame_set_group(    frame1, CPL_FRAME_GROUP_RAW);
	cpl_frame_set_level(    frame1, CPL_FRAME_LEVEL_TEMPORARY);

	/* Prepare a dummy input raw frame. */
	cpl_frame *frame2 = cpl_frame_new();
	cpl_frame_set_filename( frame2, filename2);
	cpl_frame_set_tag(      frame2, "RAW");
	cpl_frame_set_type(     frame2, CPL_FRAME_TYPE_IMAGE);
	cpl_frame_set_group(    frame2, CPL_FRAME_GROUP_RAW);
	cpl_frame_set_level(    frame2, CPL_FRAME_LEVEL_TEMPORARY);

    /* Insert frame in frameset */
    cpl_frameset_insert(frameset, frame1);
    cpl_frameset_insert(frameset, frame2);



    /****** TESTS ****/
	cpl_size     nExt      = 0;

    cpl_size     llx       = 1;
    cpl_size     lly       = 1;
    cpl_size     urx       = 10;
    cpl_size     ury       = 10;

    size_t       nmaxcache = 2;
    cpl_size     nx        = 10;
    cpl_size     ny        = 20;

    cpl_size     nslices   = 2;
    size_t       nslice    = 1;

    cpl_type     type      = CPL_TYPE_DOUBLE;

    //hdrl_iter    *it;

	/* iter with methods: next and delete & without methods: length and reset */
    hdrl_iter *it1 = hdrldemo_mef_itfiles_new(frameset, nExt);
    cpl_test_nonnull(it1);

    //it = hdrl_iter_next(it1);
    //cpl_test_nonnull(it);
    //hdrl_iter_delete(it)

/*
    int d;
    int cnt = 0;
    intptr_t ptr;
    for (hdrl_frameiter_data * h = hdrl_iter_next(it1); h != NULL; h = hdrl_iter_next(it1)) {
        cpl_test_eq(cpl_image_get(h->image, 1, 1, &d), ptr);
        cnt++;
    }
*/


    /* image iterator iterating over each image in NAXIS3 of a cube in ext
     * iter with methods: next and delete & without methods: length and reset */

    hdrl_iter *it2 = hdrldemo_cube_it_new(NULL, nExt, llx, lly, urx, ury);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(it2);

    it2 = hdrldemo_cube_it_new(filename1, nExt, llx, lly, urx, ury);
	cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
	cpl_test_null(it2);

    it2 = hdrldemo_cube_it_new(filename, nExt, llx, lly, urx, ury);
    cpl_test_nonnull(it2);

    //it = hdrl_iter_next(it2);
    //cpl_test_nonnull(it);
    //hdrl_iter_delete(it);
/*
    int d;
    int cnt = 0;
    intptr_t ptr;
    for (hdrl_frameiter_data * h = hdrl_iter_next(it2); h != NULL; h = hdrl_iter_next(it2)) {
        cpl_test_eq(cpl_image_get(h->image, 1, 1, &d), ptr);
        cnt++;
    }
*/

    /* output iterator capable of storing arbitrary large sets of data
     * by swapping to disk if necessary
     * iter with methods: next, length and delete & without method: reset */
    hdrl_iter *it3 = hdrldemo_temp_output_new(0, nx, ny);
    cpl_test_nonnull(it3);
    cpl_test_eq(hdrl_iter_length(it3), 0);
    hdrl_iter_delete(it3);
    it3 = hdrldemo_temp_output_new(nmaxcache, nx, ny);
    cpl_test_nonnull(it3);
    cpl_test_eq(hdrl_iter_length(it3), 0);
    hdrl_iter_next(it3);


    /* imagelist slice iterator coupled with a temporary output iterator */
    hdrl_iter *it4 = hdrldemo_slice_it_new(nslices, 0, it2);
    cpl_test_null(it4);


    /* Y direction image slice output iterator for hdrl_combine_it */
    hdrl_iter *it5 = hdrldemo_img_slice_out_new(nslice, nx, ny, type);
    cpl_test_nonnull(it5);
    hdrl_iter_next(it5);
    hdrl_iter_reset(it5);


    cpl_image *img = hdrldemo_img_slice_out_get_img(it5);
    cpl_test_nonnull(img);


    /* Clean up */

    cpl_imagelist_delete(imglist);

    cpl_frameset_delete(frameset);

    cpl_parameterlist_delete(parlist);

    hdrl_iter_delete(it1);
    hdrl_iter_delete(it2);
    hdrl_iter_delete(it3);
    hdrl_iter_delete(it4);
    hdrl_iter_delete(it5);

    remove(filename);
    remove(filename1);
    remove(filename2);

    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                               Private function codes
 -----------------------------------------------------------------------------*/

