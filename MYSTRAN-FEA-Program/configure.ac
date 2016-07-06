# -*- Autoconf -*- to force emacs autoconf mode

# Process this file with autoconf to produce a configure script.

AC_INIT([mystran],
        [1.0],
        [steveDoyle@gsfc.nasa.gov],
        [mystran])

AC_CONFIG_HEADERS(config.h)
AC_PREREQ(2.59)
AC_COPYRIGHT([This material is declared a work of the U.S. Government
and is not subject to copyright protection in the United States.])

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([src/Windows/Modules/CONSTANTS_1.f90])
AC_CONFIG_MACRO_DIR(aclocal)
AC_SUBST([ACLOCAL_AMFLAGS], ["-I ./aclocal \${ACLOCAL_FLAGS}"])

AC_CANONICAL_TARGET

AH_TOP([/* Required Definitions for MYSTRAN */])

AM_INIT_AUTOMAKE([1.9 tar-ustar -Wno-portability])

AC_PROG_RANLIB
AC_HEADER_STDC

AC_LANG(Fortran)

# Check Fortran compiler operation

if test -n "$FCFLAGS"
then
  strip_debug=0
else
  strip_debug=1
fi

AC_PROG_FC(gfortran ifort pgfortran pgf95)

if test $strip_debug -ne 0
then
  FCFLAGS="`echo $FCFLAGS | sed 's/-g //g' | sed 's/-g$//'`"
fi

# Determine module characteristics

AX_F90_MODULE_EXTENSION
AX_F90_MODULE_FLAG

# Check for Darwin ifort compiler bug
CHECK_IFORT_DARWIN

# define macros to call C functions
AC_FC_WRAPPERS

AC_ARG_WITH([intel-iomp],
        [AS_HELP_STRING([--with-intel-iomp],[Intel OpenMP Runtime Library])],
        [],
        [with_intel_iomp='no'])

AC_ARG_WITH([intel-mkl],
        [AS_HELP_STRING([--with-intel-mkl],[enable support for Intel Math Kernel Library])],
        [],
        [with_intel_mkl='no'])
          
AS_IF([test "x$with_intel_mkl" != "xno"], 
      [AS_IF([test "x$with_intel_iomp" != "xno"],
             [AC_CHECK_FILE([$with_intel_iomp/libiomp5.a],
                            [AC_MSG_CHECKING(for valid MKL)
                             AC_LANG_PUSH(Fortran)
                             ax_check_mkl_save_LIBS=$LIBS
                             LIBS="$LIBS -L$with_intel_mkl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$with_intel_iomp -liomp5 -lpthread"
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
                                            [AC_DEFINE([USE_INTEL_MKL], [1],
                                                       [Use Intel Math Kernel Library])
                                             intel_mkl_ldadd="-L$with_intel_mkl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$with_intel_iomp -liomp5"
                                             AC_MSG_RESULT(yes)
                                            ],
                                            [AC_MSG_FAILURE([Failed to compile MKL test program.])]
                                           )
                             LIBS=$ax_check_mkl_save_LIBS
                             AC_LANG_POP(Fortran)],
                            [AC_MSG_FAILURE([could not find Intel OpenMP runtime.])
                           ])
             ],
             [AC_MSG_FAILURE([Intel OpenMP runtime required by MKL.])])
      ],
      [AC_MSG_FAILURE([Intel Math Kernel Library is required.])])

AC_SUBST([intel_mkl_ldadd])

# Done with checks

AC_CONFIG_FILES([ \
	Makefile \
	src/Makefile \
	src/Interfaces/Makefile \
	src/USE_IFs/Makefile \
	src/Linux_differences/Makefile \
	src/Windows/Makefile \
	src/Windows/ARPACK/Makefile \
	src/Windows/BANDIT/Makefile \
	src/Windows/EMG/EMG1/Makefile \
	src/Windows/EMG/EMG2/Makefile \
	src/Windows/EMG/EMG3/Makefile \
	src/Windows/EMG/EMG4/Makefile \
	src/Windows/EMG/EMG5/Makefile \
	src/Windows/EMG/EMG6/Makefile \
	src/Windows/EMG/EMG7/Makefile \
	src/Windows/EMG/EMG8/Makefile \
	src/Windows/EMG/Makefile \
	src/Windows/LK1/L1A/Makefile \
	src/Windows/LK1/L1A-CC/Makefile \
	src/Windows/LK1/L1C/Makefile \
	src/Windows/LK1/L1E/Makefile \
	src/Windows/LK1/LINK1/Makefile \
	src/Windows/LK1/L1A-BD/Makefile \
	src/Windows/LK1/L1B/Makefile \
	src/Windows/LK1/L1D/Makefile \
	src/Windows/LK1/L1U/Makefile \
	src/Windows/LK1/Makefile \
	src/Windows/LK2/Makefile \
	src/Windows/LK3/Makefile \
	src/Windows/LK4/Makefile \
	src/Windows/LK5/Makefile \
	src/Windows/LK6/Makefile \
	src/Windows/LK9/L91/Makefile \
	src/Windows/LK9/L92/Makefile \
	src/Windows/LK9/LINK9/Makefile \
	src/Windows/LK9/Makefile \
	src/Windows/MAIN/Makefile \
	src/Windows/Modules/ARPACK/Makefile \
	src/Windows/Modules/BANDIT/Makefile \
	src/Windows/Modules/LAPACK/Makefile \
	src/Windows/Modules/Makefile \
	src/Windows/UTIL/Makefile \
	])

AC_OUTPUT
