# SYNOPSIS
#
#   AX_F90_FCFLAGS
#
# DESCRIPTION
#
#   Augment compiler flags based on compiler and host type.

AC_DEFUN([AX_F90_FCFLAGS],[
  AC_CACHE_CHECK([Fortran 90 compiler flags augmentation], [ax_cv_f90_flags],
                 [ax_cv_f90_flags=""

                  AS_CASE([$FC],
                          [ifort*],[AS_CASE([$target_os],
                                            [*darwin*],
                                            [ax_cv_f90_flags="-fno-common"])
                                   ]
                         )

                  FCFLAGS="$ax_cv_f90_flags $FCFLAGS"
                 ])
])
