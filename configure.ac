#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(APFELgrid, 0.Xb, nathan.hartland@physics.ox.ac.uk)
AC_CONFIG_MACRO_DIR([m4])

LT_INIT()

AC_CONFIG_SRCDIR([src/APFELgrid.cc])
AC_CONFIG_HEADERS([src/APFELgrid/config.h])
AC_CONFIG_FILES([Makefile src/APFELgrid/fastkernel.h])
AM_INIT_AUTOMAKE([serial-tests subdir-objects])

AC_SUBST(APFELGRID_HAVE_SSE3, ["#define APFELGRID_HAVE_SSE3 0"])
AC_SUBST(APFELGRID_HAVE_AVX, ["#define APFELGRID_HAVE_AVX 0"])
AC_SUBST(APFELGRID_HAVE_OMP, ["#define APFELGRID_HAVE_OMP 0"])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AX_EXT

# Checks for external libs
AC_SEARCH_ROOT
AC_SEARCH_APFEL
AC_SEARCH_APPLGRID
AC_SEARCH_LHAPDF
AC_SEARCH_OPENMP

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday memset mkdir pow sqrt strstr])

AC_OUTPUT
