/*
 * This file is part of the HDRL
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

#ifndef HDRLDEMO_ITERATOR_H_
#define HDRLDEMO_ITERATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include "hdrl_iter.h"
#include <cpl.h>

CPL_BEGIN_DECLS

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

/* image iterator for OmegaCAM-like multi-extension files (mef) iterating over
 * a given extension over a set of frames */
hdrl_iter * hdrldemo_mef_itfiles_new(cpl_frameset * frameset,
                                       cpl_size extname);


/* image iterator iterating over each image in NAXIS3 of a cube in ext */
hdrl_iter *
hdrldemo_cube_it_new(const char * fn, cpl_size ext,
                     cpl_size llx, cpl_size lly,
                     cpl_size urx, cpl_size ury);

/* output iterator capable of storing arbitrary large sets of data by swapping
 * to disk if necessary */
hdrl_iter * hdrldemo_temp_output_new(size_t, cpl_size, cpl_size);

/* imagelist slice iterator coupled with a temporary output iterator */
hdrl_iter * hdrldemo_slice_it_new(cpl_size nslices, cpl_size ny,
                                    hdrl_iter * outit);

/* Y direction image slice output iterator for drl_combine_it */
hdrl_iter * hdrldemo_img_slice_out_new(size_t nslice, cpl_size nx,
                                         cpl_size ny, cpl_type type);
cpl_image * hdrldemo_img_slice_out_get_img(hdrl_iter *);

CPL_END_DECLS

#endif /* HDRLDEMO_ITERATOR_H_ */
