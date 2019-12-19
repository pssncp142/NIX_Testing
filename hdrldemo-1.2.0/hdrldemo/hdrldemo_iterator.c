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


/*-----------------------------------------------------------------------------
                                   Includes
-----------------------------------------------------------------------------*/
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 500 /* posix 2001, fdopen */
#endif

// TMP
#define HDRL_USE_PRIVATE
#include "hdrldemo_iterator.h"
#include "hdrl_utils.h"

#include <cpl.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>


/*-----------------------------------------------------------------------------
                            Signatures and typedefs
-----------------------------------------------------------------------------*/

static void * hdrldemo_mef_itfiles_next(hdrl_iter * it_);
static void hdrldemo_mef_itfiles_delete(void * it);

/* image iterator for OmegaCAM-like multi-extension files (mef) iterating over
 * a given extension over a set of frames */
typedef struct {
    /* frameset with the frames one wants to iterate over*/
    cpl_frameset * frameset;
    /* current file position in the frameset*/
    cpl_size file_pos;
    /* maximum number of files in the frameset*/
    cpl_size file_pos_max;
    /* extension to work on*/
    cpl_size ext_pos;
} hdrldemo_mef_itfiles_t;

/* ------------------------------------------------------------------------- */

static void * hdrldemo_cube_it_next(hdrl_iter * it);
static void hdrldemo_cube_it_delete(void * it);

/* image iterator iterating over each image in NAXIS3 of a cube in ext */
typedef struct {
    /* filename containing cube */
    char * fn;
    /* extension to load cube from */
    cpl_size ext;
    /* current position in cube */
    cpl_size pos;
    /* window to load from cube */
    cpl_size llx;
    cpl_size lly;
    cpl_size urx;
    cpl_size ury;
    /* number of images in cube */
    cpl_size nz;
} hdrldemo_cube_it_t;

/* ------------------------------------------------------------------------- */

static void *   hdrldemo_temp_output_next(     hdrl_iter * it);
static cpl_size hdrldemo_temp_output_get_size( hdrl_iter * it);
static void     hdrldemo_temp_output_delete(   void * it);

/* output iterator capable of storing arbitrary large sets of data by swapping
 * to disk if necessary */
typedef struct {
    /* current position */
    size_t pos;
    /* image dimensions */
    cpl_size nx;
    cpl_size ny;
    /* in memory cache */
    cpl_imagelist * cache;
    /* number of images in cache */
    size_t ncached;
    /* offset of cache to full data set */
    size_t ncached_offset;
    /* maximum size of cache */
    size_t nmaxcache;
    /* backing store stream (mode a+) */
    FILE * fcache;
} hdrldemo_temp_output_t;

/* ------------------------------------------------------------------------- */

static void * hdrldemo_slice_it_next(hdrl_iter * it);

/* imagelist slice iterator coupled with a temporary output iterator */
typedef struct {
    /* current position */
    intptr_t pos;
    /* maximum number of slices(rows) returned in each iteration */
    intptr_t nslices;
    /* y dimension of images */
    intptr_t ny;
    /* number of images in imagelist */
    intptr_t nz;
    /* filled output iterator used as data source */
    hdrldemo_temp_output_t * outit;
} hdrldemo_slice_it;

/* ------------------------------------------------------------------------- */

static void * hdrldemo_img_slice_out_next(hdrl_iter * it_);
static void hdrldemo_img_slice_out_reset(hdrl_iter * it_);
static void hdrldemo_img_slice_out_delete(void * it);

/* Y direction image slice output iterator for drl_combine_it */
typedef struct {
    /* full buffer dimensions */
    cpl_size nx;
    cpl_size ny;
    /* current position */
    intptr_t iy;
    /* number of rows the returned buffer wraps */
    intptr_t sy;
    /* full buffer */
    cpl_image * oimg;
    /* wrapped buffer returned to caller, managed by iterator */
    cpl_image * oimg_wrap;
} hdrldemo_img_slice_out_t;


/*-----------------------------------------------------------------------------
                                   Code
-----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief   creates an iterator for multi-extension files
  @param   frameset   product set of frames
  @param   extname    extension name
  @return  an iterator
 */
/*----------------------------------------------------------------------------*/
hdrl_iter * hdrldemo_mef_itfiles_new(cpl_frameset * frameset, cpl_size extname)
{
    //allocate memory for the iterator
    hdrldemo_mef_itfiles_t * it = cpl_calloc(1, sizeof(*it));

    // set up the first instance of the iterator(structure)
    it->frameset = cpl_frameset_duplicate(frameset);
    it->file_pos_max = cpl_frameset_get_size(frameset);
    it->ext_pos = extname;
    // iterator should start with the first file in the frameset
    it->file_pos = 0;
    return hdrl_iter_init(hdrldemo_mef_itfiles_next,
                          NULL,
						  NULL,
                          hdrldemo_mef_itfiles_delete,
                          HDRL_ITER_OUTPUT | HDRL_ITER_IMAGE, it);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to implement the "next" on a multi-extension iterator
  @param   it_     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void * hdrldemo_mef_itfiles_next(hdrl_iter * it_)
{
    /* define a new iterator (it) and get the state from the passed
     * iterator (it_) - always to be done <- note the underscore! */
    hdrldemo_mef_itfiles_t * it = hdrl_iter_state(it_);

    // check if there is still something to do
    if ((it->file_pos) < (it->file_pos_max)){
        cpl_frame * cur_frame = cpl_frameset_get_position(it->frameset,
                        it->file_pos);
        cpl_image * cur_image = cpl_image_load(cpl_frame_get_filename(cur_frame),
                        CPL_TYPE_DOUBLE, 0, it->ext_pos);
        cpl_msg_debug(cpl_func, "Iterator returns the following file: %s",
                      cpl_frame_get_filename(cur_frame));
        // go to next file for the next iteration
        it->file_pos++;
        return cur_image;
    }
    else {
        return NULL;
    }
}
/*----------------------------------------------------------------------------*/
/**
  @brief   function to delete a multi-extension iterator
  @param   it_     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void hdrldemo_mef_itfiles_delete(void * it)
{
    //delete everything you created in the XXX_new iterator function
	hdrldemo_mef_itfiles_t * state = hdrl_iter_state(it);

    cpl_frameset_delete(state->frameset);

    // delete the state and the iterator - always to be done
    cpl_free(state);
}



/*----------------------------------------------------------------------------*/
/**
  @brief   creates an iterator for cube data files
  @param   fn   input cube file name
  @param   ext    extension id
  @param   llx lower left X image corner
  @param   lly lower left Y image corner
  @param   urx upper right X image corner
  @param   ury upper right Y image corner
  @return  an iterator
 */
/*----------------------------------------------------------------------------*/
hdrl_iter * hdrldemo_cube_it_new(const char * fn,
								 cpl_size ext,
								 cpl_size llx, cpl_size lly,
								 cpl_size urx, cpl_size ury)
{
    cpl_propertylist *plist = cpl_propertylist_load(fn, ext);
    if (plist == NULL) {
        return NULL;
    }

    if (!cpl_propertylist_has(plist, "NAXIS3")) {
    	cpl_propertylist_delete(plist);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Cube must have NAXIS3 keyword");
        return NULL;
    }

    hdrldemo_cube_it_t * it = cpl_calloc(1, sizeof(*it));

    it->nz = cpl_propertylist_get_int(plist, "NAXIS3");
    cpl_propertylist_delete(plist);
    it->fn = cpl_strdup(fn);
    it->pos = 0;
    it->ext = ext;
    it->llx = llx;
    it->lly = lly;
    it->urx = urx;
    it->ury = ury;
    return hdrl_iter_init(hdrldemo_cube_it_next,
                          NULL,
						  NULL,
                          hdrldemo_cube_it_delete,
                          HDRL_ITER_INPUT | HDRL_ITER_IMAGE, it);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to implement the "next" on a cube iterator
  @param   it     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void * hdrldemo_cube_it_next(hdrl_iter * it)
{
    hdrldemo_cube_it_t * state = hdrl_iter_state(it);
    if (state->pos >= state->nz)
        return NULL;
    return cpl_image_load_window(state->fn, CPL_TYPE_DOUBLE, state->pos++,
                                 state->ext, state->llx, state->lly,
                                 state->urx, state->ury);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to delete a cube iterator
  @param   it_     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
void hdrldemo_cube_it_delete(void * it)
{
    hdrldemo_cube_it_t * state = hdrl_iter_state(it);
    cpl_free(state->fn);
    cpl_free(state);
}




/* ---------------------------------------------------------------------------*/
/**
 * @brief
 *  output iterator capable of storing arbitrary large sets of data by swapping
 *  to disk if necessary create ouutput
 *
 * @param nmaxcache  maximum number of images to store in memory
 * @param nx         number of columns the images will have
 * @param ny         number of rows the images will have
 */
/* ---------------------------------------------------------------------------*/
hdrl_iter * hdrldemo_temp_output_new(size_t nmaxcache, cpl_size nx, cpl_size ny)
{
    int tmpfd = hdrl_get_tempfile(NULL, CPL_TRUE);
    if (tmpfd == -1) {
        return NULL;
    }

    hdrldemo_temp_output_t * it = cpl_calloc(1, sizeof(*it));
    /* open in append + read mode
     * one could use a fits file but that complicates resource cleanup in case
     * of crashes and signals */
    it->fcache = fdopen(tmpfd, "a+");
    it->cache = cpl_imagelist_new();
    it->ncached_offset = 0;
    it->ncached = 0;
    it->nmaxcache = nmaxcache == 0 ? 1 : nmaxcache;
    it->pos = 0;
    it->nx = nx;
    it->ny = ny;
    return hdrl_iter_init(hdrldemo_temp_output_next,
    					  NULL,
                          hdrldemo_temp_output_get_size,
                          hdrldemo_temp_output_delete,
                          HDRL_ITER_OUTPUT | HDRL_ITER_IMAGE, it);
}

/*----------------------------------------------------------------------------*/
/**
  @brief append image to a temporary file on disk, must be opened in append mode
  @param   f    file
  @param   image image
  @return
 */
/*----------------------------------------------------------------------------*/
static size_t append_image_to_file(FILE * f, cpl_image * image)
{
    size_t npix = cpl_image_get_size_x(image) * cpl_image_get_size_y(image);
    size_t r = fwrite(cpl_image_get_data(image), sizeof(double), npix, f);
    if (r !=  npix) {
        cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                              "Writing temporary file failed: %s",
                              strerror(errno));
    }
    return r;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to implement the "next" on a temp iterator
  @param   it     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void * hdrldemo_temp_output_next(hdrl_iter * it)
{
    hdrldemo_temp_output_t * state = hdrl_iter_state(it);
    /* number of files to swap to disk when cache is full, not too big to
     * ensure the cache is always well filled, not too small so the IO is
     * efficient */
    const size_t nswap = CX_MIN(state->nmaxcache, 1u + (1u << 24u) /
                                (state->nx * state->ny * sizeof(double)));
    if (state->ncached >= state->nmaxcache) {
        cpl_msg_debug(cpl_func, "Out of memory, swapping %zu images to disk",
                      nswap);
        for (size_t i = 0; i < nswap; i++) {
            cpl_image * img = cpl_imagelist_unset(state->cache, 0);
            append_image_to_file(state->fcache, img);
            cpl_image_delete(img);
        }
        state->ncached_offset += nswap;
        state->ncached -= nswap;
    }

    /* put the new image in the memory cache and return it to be filled */
    {
        cpl_image * img = cpl_image_new(state->nx, state->ny, CPL_TYPE_DOUBLE);
        cpl_imagelist_set(state->cache, img, state->ncached);
        state->ncached += 1;
        return img;
    }
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to get the iterator size
  @param   it     iterator
  @return  iterator size
 */
/*----------------------------------------------------------------------------*/
static cpl_size hdrldemo_temp_output_get_size(hdrl_iter * it)
{
    hdrldemo_temp_output_t * state = hdrl_iter_state(it);
    return state->ncached_offset + state->ncached;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to delete a temp iterator
  @param   it_     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void hdrldemo_temp_output_delete(void * it)
{
    hdrldemo_temp_output_t * state = hdrl_iter_state(it);
    cpl_imagelist_delete(state->cache);
    fclose(state->fcache);
    cpl_free(state);
}



/* ---------------------------------------------------------------------------*/
/**
 * @brief create imagelist iterator iterating from a filled output iterator
 *
 * @param nslices number of rows the returned imagelist should have
 * @param ny      total number of rows
 * @param outit   output iterator used as data source
 */
/* ---------------------------------------------------------------------------*/
hdrl_iter *
hdrldemo_slice_it_new(cpl_size nslices, cpl_size ny, hdrl_iter * outit)
{
    if (!hdrl_iter_check(outit, HDRL_ITER_OUTPUT | HDRL_ITER_IMAGE)) {
        return NULL;
    }
    hdrldemo_slice_it * it = cpl_calloc(1, sizeof(*it));
    it->nz = hdrl_iter_length(outit);
    if (it->nz < 0) {
        cpl_free(it);
        return NULL;
    }
    it->pos = 1;
    it->nslices = nslices;
    it->ny = ny;
    it->outit = hdrl_iter_state(outit);
    return hdrl_iter_init(hdrldemo_slice_it_next,
                          NULL,
						  NULL,
						  NULL,
                          HDRL_ITER_INPUT | HDRL_ITER_IMAGELIST, it);
}

/*----------------------------------------------------------------------------*/
/**
  @brief load image window from temporary file on disk
  @param f    file
  @param pos
  @param llx lower left X image corner
  @param lly lower left Y image corner
  @param urx upper right X image corner
  @param ury upper right Y image corner
  @param nx    number of columns the images will have
  @param ny    number of rows the images will have
  @return
 */
/*----------------------------------------------------------------------------*/
/*  */
static cpl_image * load_image_from_file( FILE * f, size_t pos,
										 cpl_size llx, cpl_size lly,
										 cpl_size urx, cpl_size ury,
										 cpl_size nx, cpl_size ny)
{
    const size_t snx = urx - llx + 1;
    const size_t sny = ury - lly + 1;
    const size_t esize = sizeof(double);
    cpl_image * img = cpl_image_new(snx, sny, CPL_TYPE_DOUBLE);
    int e = 0;
    double * d = cpl_image_get_data_double(img);

    /* read full rows */
    if (snx == (size_t)nx) {
        e = fseek(f, (nx * ny * pos + (lly - 1) * nx) * esize, SEEK_SET) !=  0;
        e = !e && fread(d, esize, nx * sny, f) != nx *sny;
    }
    else {
        off_t seekpos = (nx * ny * pos + (lly - 1) * nx + llx - 1) * esize;
        e = fseek(f, seekpos, SEEK_SET) != 0;
        e = !e && fread(d, esize, snx, f) != snx;
        for (size_t y = 1; e == 0 && y < sny; y++) {
            e = fseek(f, (nx - snx) * esize, SEEK_CUR) != 0;
            e = !e && fread(&d[y * snx], esize, snx, f) != snx;
        }
    }

    if (e) {
        cpl_image_delete(img);
        cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                              "Reading from temporary file failed: %s",
                              strerror(errno));
        return NULL;
    }
    else
        return img;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief   Utility to get from a given iterator file cache the image
           corresponding to the position pos and region (llx,lly,urx,ury)
  @param   it  iterator
  @param   pos position in the file cache
  @param   llx lower left X image corner
  @param   lly lower left Y image corner
  @param   urx upper right X image corner
  @param   ury upper right Y image corner
  @return  cpl_image
 */
/*----------------------------------------------------------------------------*/
static cpl_image *
hdrldemo_temp_output_get(hdrldemo_temp_output_t * it, size_t pos,
                         cpl_size llx, cpl_size lly,
                         cpl_size urx, cpl_size ury)
{
    if (urx <= 0)
        urx = it->nx;
    if (ury <= 0)
        ury = it->ny;

    if (pos < it->ncached_offset) {
        cpl_msg_debug(cpl_func, "getting %zu from disk (%zu)", pos,
                        it->ncached_offset);
        return load_image_from_file(it->fcache, pos, llx, lly, urx, ury,
                                    it->nx, it->ny);
    }
    else if (pos < it->ncached_offset + it->ncached) {
        cpl_msg_debug(cpl_func, "getting %zu from cache", pos);
        return cpl_image_extract(cpl_imagelist_get(it->cache,
                                                   pos - it->ncached_offset),
                                 llx, lly, urx, ury);
    }
    else {
        cpl_msg_debug(cpl_func, "no image for position %zu", pos);
        return NULL;
    }
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to implement the "next" on a slice iterator
  @param   it     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
void * hdrldemo_slice_it_next(hdrl_iter * it)
{
    hdrldemo_slice_it * state = hdrl_iter_state(it);
    if (state->pos >= state->ny) {
        return NULL;
    }

    cpl_imagelist * l = cpl_imagelist_new();
    /* build an imagelist slice from the images in the output iterator */
    cpl_size ury = CX_MIN(state->pos + state->nslices - 1, state->ny);
    for (intptr_t z = 0; z < state->nz; z++) {
        cpl_image * img = hdrldemo_temp_output_get(state->outit, z,
                                                   1, state->pos, -1, ury);
        cpl_imagelist_set(l, img, z);
    }
    state->pos = ury + 1;
    return l;
}



/* ---------------------------------------------------------------------------*/
/**
 * @brief create image slice output iterator
 *
 * @param nslice  number of rows the returned buffer should have
 * @param nx      total number of columns of the image to be filled
 * @param ny      total number of rows the image to be filled
 * @param type    type of the image buffer (int or double)
 */
/* ---------------------------------------------------------------------------*/
hdrl_iter *
hdrldemo_img_slice_out_new(size_t nslice, cpl_size nx, cpl_size ny,
                           cpl_type type)
{
    if (type != CPL_TYPE_DOUBLE && type != CPL_TYPE_INT) {
        cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                              "Only double and integer images supported");
        return NULL;
    }
    hdrldemo_img_slice_out_t * it = cpl_calloc(1, sizeof(*it));
    it->nx = nx;
    it->ny = ny;
    it->iy = 1;
    it->sy = nslice;
    it->oimg_wrap = NULL;
    it->oimg = cpl_image_new(nx, ny, type);

    hdrl_iter *newIt = hdrl_iter_init(hdrldemo_img_slice_out_next,
                          hdrldemo_img_slice_out_reset,
						  NULL,
                          hdrldemo_img_slice_out_delete,
                          HDRL_ITER_OUTPUT | HDRL_ITER_IMAGE, it);
    return newIt;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to implement the "next" on an image-slice iterator
  @param   it_     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void * hdrldemo_img_slice_out_next(hdrl_iter * it_)
{
    hdrldemo_img_slice_out_t * it = hdrl_iter_state(it_);
    /* upper index to load */
    cpl_size upy = CX_MIN(it->iy + it->sy - 1, it->ny);
    /* done? */
    if (it->iy >= it->ny)
        return NULL;
    /* remove last wrap */
    cpl_image_unwrap(it->oimg_wrap);
    /* create new wrap and return it */
    if (cpl_image_get_type(it->oimg) == CPL_TYPE_DOUBLE) {
        double * d = cpl_image_get_data_double(it->oimg);
        it->oimg_wrap = cpl_image_wrap(it->nx, upy - it->iy + 1,
                                       CPL_TYPE_DOUBLE,
                                       &d[(it->iy - 1) * it->nx]);
    } else {
        int * d = cpl_image_get_data_int(it->oimg);
        it->oimg_wrap = cpl_image_wrap(it->nx, upy - it->iy + 1,
                                       CPL_TYPE_INT,
                                       &d[(it->iy - 1) * it->nx]);
    }
    it->iy = upy + 1;
    return it->oimg_wrap;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to reset an image-slice iterator
  @param   it_     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void hdrldemo_img_slice_out_reset(hdrl_iter * it_)
{
    hdrldemo_img_slice_out_t * it = hdrl_iter_state(it_);
    it->iy = 1;
    cpl_image_unwrap(it->oimg_wrap);
    it->oimg_wrap = NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to delete an image-slice iterator
  @param   it     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
static void hdrldemo_img_slice_out_delete(void * it)
{
    hdrldemo_img_slice_out_t * state = hdrl_iter_state(it);
    cpl_image_delete(state->oimg);
    cpl_image_unwrap(state->oimg_wrap);
    cpl_free(state);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   function to get the image to the iterator
  @param   it     iterator
  @return  void
 */
/*----------------------------------------------------------------------------*/
cpl_image * hdrldemo_img_slice_out_get_img(hdrl_iter * it)
{
    hdrldemo_img_slice_out_t * state = hdrl_iter_state(it);
    return state->oimg;
}
