# SYNOPSIS
#
#   AX_INTEL_MKL
#
# DESCRIPTION
#
#   Check for Intel Math Kernel Library.

AC_DEFUN([AX_INTEL_MKL],[

AC_ARG_WITH([intel-threads],
        [AS_HELP_STRING([--with-intel-threads],[Intel OpenMP Runtime Library])],
        [],
        [with_intel_threads='no'])

AC_ARG_WITH([intel-mkl],
        [AS_HELP_STRING([--with-intel-mkl],[enable support for Intel Math Kernel Library])],
        [],
        [with_intel_mkl='no'])

AS_IF([test "x$with_intel_threads" != "xno"],
      [AC_CHECK_FILE([$with_intel_threads/libiomp5.a],
                     [have_threads='yes'],
                     [AC_MSG_FAILURE([could not find Intel OpenMP runtime.])])],
      [have_threads='no'])

AS_IF([test "x$with_intel_mkl" != "xno"], 
      [AC_MSG_CHECKING(for valid MKL)
       AC_LANG_PUSH(Fortran)
       ax_check_mkl_save_LIBS=$LIBS
       AS_IF([test "x$have_threads" != "xno"], 
             [intel_mkl_ldadd="-L$with_intel_mkl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$with_intel_threads -liomp5 -lpthread"],
             [intel_mkl_ldadd="-L$with_intel_mkl -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"])
       LIBS="$LIBS $intel_mkl_ldadd"
       AC_LINK_IFELSE([

      module mkl_dss_private
        type mkl_dss_handle
          integer(kind=8) dummy
        end type mkl_dss_handle
      end module mkl_dss_private

      program main
        use mkl_dss_private

        integer, parameter :: mkl_dss_defaults = 0
        interface
          function dss_create( handle, opt )
            use mkl_dss_private
            type(mkl_dss_handle), intent(out)   :: handle
            integer,              intent(in)    :: opt
            integer                             :: dss_create
          end function dss_create
        end interface

        type(mkl_dss_handle)            :: handle
        integer                         :: intmkl_ier
        intmkl_ier = dss_create ( handle, mkl_dss_defaults )
        if (intmkl_ier .ne. mkl_dss_success) then
          write(*,*) 'Bad dss_create'
        endif
      end program main
                      ],
                      [AC_DEFINE([USE_INTEL_MKL], [1], [Use Intel Math Kernel Library])
                       AC_MSG_RESULT(yes)],
                      [AC_MSG_FAILURE([Failed to compile MKL test program.])])
       LIBS=$ax_check_mkl_save_LIBS
       AC_LANG_POP(Fortran)],
      [AC_MSG_FAILURE([Intel Math Kernel Library is required.])])

AC_SUBST([intel_mkl_ldadd])
])
