AC_DEFUN([AC_SEARCH_OPENMP],[

AC_ARG_ENABLE([openmp],
  [AC_HELP_STRING(--enable-openmp, [Enable OpenMP multicore.])],
  [enable_openmp=yes], [enable_openmp=no])
if test x$enable_openmp == xyes; then
  AC_OPENMP
  AC_SUBST(APFELGRID_HAVE_OMP, ["#define APFELGRID_HAVE_OMP 1"])
fi

])
