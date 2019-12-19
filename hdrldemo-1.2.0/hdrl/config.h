/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `floor' function. */
#define HAVE_FLOOR 1

/* Define to 1 iff you have GSL */
#define HAVE_GSL 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 iff you have GSL */
#define HAVE_LIBGSL 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define if you have the `strdup' function */
/* #undef HAVE_STRDUP */

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* HDRL binary age */
#define HDRL_BINARY_AGE 200

/* HDRL binary version number */
#define HDRL_BINARY_VERSION 10200

/* HDRL interface age */
#define HDRL_INTERFACE_AGE 0

/* HDRL major version number */
#define HDRL_MAJOR_VERSION 1

/* HDRL micro version number */
#define HDRL_MICRO_VERSION 0

/* HDRL minor version number */
#define HDRL_MINOR_VERSION 2

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "hdrl"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "usd-help@eso.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "HDRL Pipeline"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "HDRL Pipeline 1.2.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "hdrl"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.2.0"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.2.0"

/* Define if using the dmalloc debugging malloc package */
/* #undef WITH_DMALLOC */

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif


#ifndef HAVE_STRDUP
#  define strdup  cx_strdup
#endif
              
