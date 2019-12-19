

# HDRLDEMO_SET_PREFIX(PREFIX)
#---------------------------
AC_DEFUN([HDRLDEMO_SET_PREFIX],
[
    unset CDPATH
    # make $PIPE_HOME the default for the installation
    AC_PREFIX_DEFAULT($1)

    if test "x$prefix" = "xNONE"; then
        prefix=$ac_default_prefix
        ac_configure_args="$ac_configure_args --prefix $prefix"
    fi

    if test "x$exec_prefix" = "xNONE"; then
        exec_prefix=$prefix
    fi

])

# HDRLDEMO_SET_VERSION_INFO(VERSION, [CURRENT], [REVISION], [AGE])
#----------------------------------------------------------------
# Setup various version information, especially the libtool versioning
AC_DEFUN([HDRLDEMO_SET_VERSION_INFO],
[
    hdrldemo_version=`echo "$1" | sed -e 's/[[a-z,A-Z]].*$//'`

    hdrldemo_major_version=`echo "$hdrldemo_version" | \
        sed 's/\([[0-9]]*\).\(.*\)/\1/'`
    hdrldemo_minor_version=`echo "$hdrldemo_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\)\(.*\)/\2/'`
    hdrldemo_micro_version=`echo "$hdrldemo_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test -z "$hdrldemo_major_version"; then hdrldemo_major_version=0
    fi

    if test -z "$hdrldemo_minor_version"; then hdrldemo_minor_version=0
    fi

    if test -z "$hdrldemo_micro_version"; then hdrldemo_micro_version=0
    fi

    HDRLDEMO_VERSION="$hdrldemo_version"
    HDRLDEMO_MAJOR_VERSION=$hdrldemo_major_version
    HDRLDEMO_MINOR_VERSION=$hdrldemo_minor_version
    HDRLDEMO_MICRO_VERSION=$hdrldemo_micro_version

    if test -z "$4"; then HDRLDEMO_INTERFACE_AGE=0
    else HDRLDEMO_INTERFACE_AGE="$4"
    fi

    HDRLDEMO_BINARY_AGE=`expr 100 '*' $HDRLDEMO_MINOR_VERSION + $HDRLDEMO_MICRO_VERSION`
    HDRLDEMO_BINARY_VERSION=`expr 10000 '*' $HDRLDEMO_MAJOR_VERSION + \
                          $HDRLDEMO_BINARY_AGE`

    AC_SUBST(HDRLDEMO_VERSION)
    AC_SUBST(HDRLDEMO_MAJOR_VERSION)
    AC_SUBST(HDRLDEMO_MINOR_VERSION)
    AC_SUBST(HDRLDEMO_MICRO_VERSION)
    AC_SUBST(HDRLDEMO_INTERFACE_AGE)
    AC_SUBST(HDRLDEMO_BINARY_VERSION)
    AC_SUBST(HDRLDEMO_BINARY_AGE)

    AC_DEFINE_UNQUOTED(HDRLDEMO_MAJOR_VERSION, $HDRLDEMO_MAJOR_VERSION,
                       [HDRLDEMO major version number])
    AC_DEFINE_UNQUOTED(HDRLDEMO_MINOR_VERSION, $HDRLDEMO_MINOR_VERSION,
                       [HDRLDEMO minor version number])
    AC_DEFINE_UNQUOTED(HDRLDEMO_MICRO_VERSION, $HDRLDEMO_MICRO_VERSION,
                       [HDRLDEMO micro version number])
    AC_DEFINE_UNQUOTED(HDRLDEMO_INTERFACE_AGE, $HDRLDEMO_INTERFACE_AGE,
                       [HDRLDEMO interface age])
    AC_DEFINE_UNQUOTED(HDRLDEMO_BINARY_VERSION, $HDRLDEMO_BINARY_VERSION,
                       [HDRLDEMO binary version number])
    AC_DEFINE_UNQUOTED(HDRLDEMO_BINARY_AGE, $HDRLDEMO_BINARY_AGE,
                       [HDRLDEMO binary age])

    ESO_SET_LIBRARY_VERSION([$2], [$3], [$4])
])


# HDRLDEMO_SET_PATHS
#------------------
# Define auxiliary directories of the installed directory tree.
AC_DEFUN([HDRLDEMO_SET_PATHS],
[

    if test -z "$plugindir"; then
        plugindir='${libdir}/esopipes-plugins/${PACKAGE}-${VERSION}'
    fi

    if test -z "$privatelibdir"; then
        privatelibdir='${libdir}/${PACKAGE}-${VERSION}'
    fi

    if test -z "$apidocdir"; then
        apidocdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}/html'
    fi

    if test -z "$pipedocsdir"; then
        pipedocsdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}'
    fi

    if test -z "$configdir"; then
       configdir='${datadir}/${PACKAGE}/config'
    fi

    if test -z "$wkfextradir"; then
        wkfextradir='${datadir}/esopipes/${PACKAGE}-${VERSION}/reflex'
    fi

    if test -z "$wkfcopydir"; then
        wkfcopydir='${datadir}/reflex/workflows/${PACKAGE}-${VERSION}'
    fi

    AC_SUBST(plugindir)
    AC_SUBST(privatelibdir)
    AC_SUBST(apidocdir)
    AC_SUBST(pipedocsdir)
    AC_SUBST(configdir)
    AC_SUBST(wkfextradir)
    AC_SUBST(wkfcopydir)


    # Define a preprocesor symbol for the plugin search paths

    AC_DEFINE_UNQUOTED(HDRLDEMO_PLUGIN_DIR, "${PACKAGE}/plugins",
                       [Plugin directory tree prefix])

    eval plugin_dir="$plugindir"
    plugin_path=`eval echo $plugin_dir | \
                sed -e "s/\/${PACKAGE}-${VERSION}.*$//"`

    AC_DEFINE_UNQUOTED(HDRLDEMO_PLUGIN_PATH, "$plugin_path",
                       [Absolute path to the plugin directory tree])

])


# HDRLDEMO_CREATE_SYMBOLS
#-----------------------
# Define include and library related makefile symbols
AC_DEFUN([HDRLDEMO_CREATE_SYMBOLS],
[

    # Symbols for package include file and library search paths

    HDRLDEMO_INCLUDES='-I$(top_srcdir)/hdrldemo'
    HDRLDEMO_LDFLAGS='-L$(top_builddir)/hdrldemo'

    # Library aliases

    LIBHDRLDEMO='$(top_builddir)/hdrldemo/libhdrldemo.la'

    # Substitute the defined symbols

    AC_SUBST(HDRLDEMO_INCLUDES)
    AC_SUBST(HDRLDEMO_LDFLAGS)

    AC_SUBST(LIBHDRL)
    AC_SUBST(LIBHDRLDEMO)

    # Check for CPL and user defined libraries
    AC_REQUIRE([ESO_CHECK_CPL])

    all_includes='$(HDRLDEMO_INCLUDES) $(HDRL_INCLUDES) $(CPL_INCLUDES)'
    all_ldflags='$(HDRLDEMO_LDFLAGS) $(HDRL_LDFLAGS) $(CPL_LDFLAGS)'

    AC_SUBST(all_includes)
    AC_SUBST(all_ldflags)
])
