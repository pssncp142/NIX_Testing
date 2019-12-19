## Process this file with automake to produce Makefile.in

##   This file is part of the HDRLDEMO Pipeline
##   Copyright (C) 2013,2014 European Southern Observatory

##   This library is free software; you can redistribute it and/or modify
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

SUBDIRS = . tests

if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"HdrldemoLib\" -DHDRL_USE_EXPERIMENTAL $(all_includes)

noinst_HEADERS =  hdrldemo_dfs.h hdrldemo_iterator.h hdrldemo_utils.h


EXTRA_DIST =

pkginclude_HEADERS = 

privatelib_LTLIBRARIES = libhdrldemo.la

libhdrldemo_la_SOURCES = hdrldemo_dfs.c hdrldemo_iterator.c  hdrldemo_utils.c

libhdrldemo_la_LDFLAGS = $(HDRL_LDFLAGS) $(CPL_LDFLAGS) -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libhdrldemo_la_LIBADD = $(HDRL_LIBS)
libhdrldemo_la_DEPENDENCIES = $(LIBHDRL)
