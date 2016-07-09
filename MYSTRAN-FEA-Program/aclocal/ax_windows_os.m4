# SYNOPSIS
#
#   AX_WINDOWS_OS
#
# DESCRIPTION
#
#   Set up conditional compilation based on if target os is Windows.

AC_DEFUN([AX_WINDOWS_OS],[
  AC_CACHE_CHECK([Windows compilation],
                 [ax_cv_windows_os],
                 [AS_CASE([$target_os],
                          [*linux*],[ax_cv_windows_os=false],
                          [*darwin*],[ax_cv_windows_os=false],
                          [ax_cv_windows_os=true])
                 ])
  AM_CONDITIONAL(BUILD_WINDOWS,[test x$ax_cv_windows_os = xtrue])
])
