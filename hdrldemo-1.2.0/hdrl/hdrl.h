/* $Id: hdrl.h,v 1.4 2013-10-23 09:42:14 jtaylor Exp $
 *
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

/*----------------------------------------------------------------------------*/
/* Documentation here is used in the reference manual Main Page */
/** @mainpage HDRL Introduction
  
    TEXT
 
    \section releases Releases

    List of release notes for diff releases

    \section usage Usage

    Externals Defintions description
    
    \section hdrldemop

    Offered recipes
    
    \section links Links
    CPL page
    Pipelines page
  
 */
/*----------------------------------------------------------------------------*/

#ifndef HDRL_H
#define HDRL_H

  #include "hdrl_image.h"
  #include "hdrl_imagelist.h"
  #include "hdrl_parameter.h"
  #include "hdrl_imagelist_view.h"
  #include "hdrl_overscan.h"
  #include "hdrl_buffer.h"
  #include "hdrl_collapse.h"
  #include "hdrl_lacosmics.h"
  #include "hdrl_bpm_2d.h"
  #include "hdrl_bpm_3d.h"
  #include "hdrl_bpm_fit.h"
  #include "hdrl_bpm_utils.h"
  #include "hdrl_fit.h"
  #include "hdrl_strehl.h"
  #include "hdrl_flat.h"
  #include "hdrl_catalogue.h"
  #include "hdrl_random.h"
  #include "hdrl_iter.h"
  #include "hdrl_frameiter.h"
  #include "hdrl_multiiter.h"
  #include "hdrl_fringe.h"
  #include "hdrl_spectrum.h"
  #include "hdrl_spectrumlist.h"
  #include "hdrl_response.h"
  #include "hdrl_efficiency.h"
  #include "hdrl_spectrum_resample.h"
  #include "hdrl_dar.h"
  #include "hdrl_fpn.h"

#endif  /* HDRL_H */
