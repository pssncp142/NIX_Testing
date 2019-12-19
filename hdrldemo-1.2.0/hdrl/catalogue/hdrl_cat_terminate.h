/*
 * This file is part of the HDRL
 * Copyright (C) 2017 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef HDRL_TERMINATE_H
#define HDRL_TERMINATE_H


#include "hdrl_cat_def.h"


void hdrl_terminate(   ap_t *ap, double gain, cpl_size *nobjects,
                         cpl_table *tab, hdrl_casu_result *res);
void hdrl_apfu(        ap_t *ap);
void hdrl_extract_data(ap_t *ap, cpl_size ip);
void hdrl_restack(     ap_t *ap, cpl_size ip);


#endif /* HDRL_RADII_H */
