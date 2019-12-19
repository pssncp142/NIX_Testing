## Process this file with automake to produce Makefile.in

##   This file is part of the ESO Common Pipeline Library
##   Copyright (C) 2012,2013 European Southern Observatory

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

AM_CPPFLAGS = $(all_includes) $(HDRL_INCLUDES) -DHDRL_USE_PRIVATE
AM_LDFLAGS = $(CPL_LDFLAGS) $(HDRL_LDFLAGS)

LDADD = $(top_builddir)/hdrl/libhdrl.la ../libhdrldemo.la \
        $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLDRS) $(LIBCPLUI) $(LIBCEXT) \
        -lgsl -lgslcblas -lm
 
check_PROGRAMS = hdrldemo_dfs-test \
    hdrldemo_utils-test \
    hdrldemo_iterator-test

hdrldemo_dfs_test_SOURCES = hdrldemo_dfs-test.c
hdrldemo_utils_test_SOURCES = hdrldemo_utils-test.c
hdrldemo_iterator_test_SOURCES = hdrldemo_iterator-test.c

# Be sure to reexport important environment variables.
TESTS_ENVIRONMENT = MAKE="$(MAKE)" CC="$(CC)" CFLAGS="$(CFLAGS)" \
	CPPFLAGS="$(CPPFLAGS)" LD="$(LD)" LDFLAGS="$(LDFLAGS)" \
	LIBS="$(LIBS)" LN_S="$(LN_S)" NM="$(NM)" RANLIB="$(RANLIB)" \
	OBJEXT="$(OBJEXT)" EXEEXT="$(EXEEXT)" \
	MALLOC_PERTURB_=31 MALLOC_CHECK_=2 \
	GOMP_SPINCOUNT=10000

TESTS = $(check_PROGRAMS) 

# Expected failures
XFAIL_TESTS =

# We need to remove any files that the above tests created.
clean-local:
	$(RM) *.fits *.tfits *.log
