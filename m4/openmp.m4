AC_DEFUN([AC_SEARCH_OPENMP],[

AC_ARG_ENABLE([openmp],
    AS_HELP_STRING([--enable-openmp], [Enable OpenMP multicore]))

AS_IF([test "x$enable_openmp" = "xyes"], [
  AC_OPENMP
  AC_SUBST(APFELGRID_HAVE_OMP, ["#define APFELGRID_HAVE_OMP 1"])
])

])
