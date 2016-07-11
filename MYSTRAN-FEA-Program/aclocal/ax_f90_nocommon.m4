# SYNOPSIS
#
#   AX_F90_NOCOMMON
#
# DESCRIPTION
#
#   Determine if compiler needs flags for special treatment of Fortran
#   modules with data but no exectutable functions.

AC_DEFUN([AX_F90_NOCOMMON],[
  AC_CACHE_CHECK([Special treatment of non-executable modules], [ax_cv_no_common],
                 [ax_cv_no_common=''

                  AS_CASE([$FC],
                          [ifort*],[AS_CASE([$target_os],
                                            [*darwin*],
                                            [ax_cv_no_common="-fno-common"])
                                   ]
                         )

                 ])
  AC_SUBST([ax_cv_no_common])
])
