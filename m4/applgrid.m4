# AC_SEARCH_APPLGRID(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_APPLGRID], [

AC_PATH_PROG(APPLGRIDCONFIG, applgrid-config, [], [$PATH])
if test -f "$APPLGRIDCONFIG"; then
  APPLGRID_CPPFLAGS=`$APPLGRIDCONFIG --cxxflags`
  APPLGRID_LDFLAGS=`$APPLGRIDCONFIG --ldflags`
else
  AC_MSG_ERROR([APPLgrid cannot be found!])
  exit 1
fi
AC_SUBST(APPLGRID_CPPFLAGS)
AC_SUBST(APPLGRID_LDFLAGS)
$1
])

