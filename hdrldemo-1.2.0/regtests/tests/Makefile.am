## Process this file with automake to produce Makefile.in

##   This file is part of the AMBER Pipeline
##   Copyright (C) 2002-2004 European Southern Observatory

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


if MAINTAINER_MODE

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in


endif

reg0:
	p_regtest  -force -id 0 ../../${PACKAGE}_reg/sof

reg1:
	p_regtest  -force -id 1 ../../${PACKAGE}_reg/sof

valgrind:
	p_regtest  -valgrind -preload  ../../${PACKAGE}_reg/sof

regcheck:
	@if test ! -d "${SOF_DATA}" ; then \
		echo "ERROR: The SOF_DATA environment variable was not set to point to the regression data." 1>&2 ; \
		exit 1 ; \
	fi
	@if test -d ../../${PACKAGE}_reg/sof/; then esorex --version ; fi
#	if test -d ../../${PACKAGE}_reg; then p_regtest  -force -id 0 ../../${PACKAGE}_reg/sof; fi
#       link to the Reference data under ${SOF_DATA}/${PACKAGE}/ReferenceData/
	if test -d ../../${PACKAGE}_reg; then  for i in ${SOF_DATA}/${PACKAGE}/ReferenceData/*.0; do ln -f -s  $$i ../../${PACKAGE}_reg/sof/ ; done ; fi
	if test -d ../../${PACKAGE}_reg; then p_regtest  -force -id 1 ../../${PACKAGE}_reg/sof; fi 
	@if test "$$?x" = "0x"; then echo PASS: regtests; fi
	@if test "$$?x" != "0x"; then echo FAILED: regtests; fi

check:
	@if test -d /NRI && (echo `hostname` | grep '^vm.*nri$$' > /dev/null 2>&1) && (echo `pwd` | grep '^/disk.*/NRI/DFS' > /dev/null 2>&1) ; then \
		echo "Detected NRI environment, will run regression tests automatically..." ; \
		$(MAKE) $(AM_MAKEFLAGS) regcheck ; \
	else \
		if test x"${FORCE_REGCHECK}" = xyes ; then \
			$(MAKE) $(AM_MAKEFLAGS) regcheck ; \
		else \
			echo "Skipping regression tests. To execute them, either run 'make regcheck' or set the environment variable FORCE_REGCHECK=yes." ; \
		fi ; \
	fi

