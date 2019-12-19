## Process this file with automake to produce Makefile.in

##   This file is part of the HDRLDEMO Pipeline
##   Copyright (C) 2013,2014 European Southern Observatory

##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.

##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

if MAINTAINER_MODE

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in

endif

AM_CPPFLAGS = -DHDRL_USE_EXPERIMENTAL $(all_includes)

EXTRA_DIST = 

noinst_HEADERS = 

plugin_LTLIBRARIES =    hdrldemo_bias.la \
			            hdrldemo_strehl.la \
			            hdrldemo_bpm_2d.la \
			            hdrldemo_fpn.la \
			            hdrldemo_dar.la \
			            hdrldemo_bpm_3d.la \
			            hdrldemo_bpm_fit.la \
			            hdrldemo_bpm_lacosmic.la \
			            hdrldemo_flat_opt.la \
			            hdrldemo_flat_nir.la \
			            hdrldemo_fringe.la \
			            hdrldemo_overscan.la \
			            hdrldemo_efficiency.la \
			            hdrldemo_response.la \
			            hdrldemo_catalogue.la

hdrldemo_bias_la_SOURCES = hdrldemo_bias.c
hdrldemo_bias_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_bias_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_bias_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_strehl_la_SOURCES = hdrldemo_strehl.c
hdrldemo_strehl_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_strehl_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_strehl_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_bpm_2d_la_SOURCES = hdrldemo_bpm_2d.c
hdrldemo_bpm_2d_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_bpm_2d_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_bpm_2d_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_fpn_la_SOURCES = hdrldemo_fpn.c
hdrldemo_fpn_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_fpn_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_fpn_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_bpm_3d_la_SOURCES = hdrldemo_bpm_3d.c
hdrldemo_bpm_3d_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_bpm_3d_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_bpm_3d_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_bpm_fit_la_SOURCES = hdrldemo_bpm_fit.c
hdrldemo_bpm_fit_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_bpm_fit_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_bpm_fit_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_bpm_lacosmic_la_SOURCES = hdrldemo_bpm_lacosmic.c
hdrldemo_bpm_lacosmic_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_bpm_lacosmic_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_bpm_lacosmic_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_flat_opt_la_SOURCES = hdrldemo_flat_opt.c
hdrldemo_flat_opt_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_flat_opt_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_flat_opt_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_flat_nir_la_SOURCES = hdrldemo_flat_nir.c
hdrldemo_flat_nir_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_flat_nir_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_flat_nir_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_fringe_la_SOURCES = hdrldemo_fringe.c
hdrldemo_fringe_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_fringe_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_fringe_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_overscan_la_SOURCES = hdrldemo_overscan.c
hdrldemo_overscan_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_overscan_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_overscan_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_catalogue_la_SOURCES = hdrldemo_catalogue.c
hdrldemo_catalogue_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_catalogue_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_catalogue_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_efficiency_la_SOURCES = hdrldemo_efficiency.c
hdrldemo_efficiency_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_efficiency_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_efficiency_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_response_la_SOURCES = hdrldemo_response.c
hdrldemo_response_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_response_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_response_la_DEPENDENCIES = $(LIBHDRL)

hdrldemo_dar_la_SOURCES = hdrldemo_dar.c
hdrldemo_dar_la_LIBADD = $(HDRL_LIBS) $(LIBHDRLDEMO)
hdrldemo_dar_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version
hdrldemo_dar_la_DEPENDENCIES = $(LIBHDRL)
