# ===========================================================================
#  extracted from http://git.savannah.gnu.org/gitweb/?p=autoconf.git;a=blob_plain;f=lib/autoconf/fortran.m4;hb=3e3412d3dd1ff7734f0a9ccd1a55252bd7ea8790
# ===========================================================================
#
# SYNOPSIS
#
#   AX_F90_MODULE_FLAG
#
# DESCRIPTION
#
#   Find a flag to include Fortran 90 modules from another directory.
#   If successful, run ACTION-IF-SUCCESS (defaults to nothing), otherwise
#   run ACTION-IF-FAILURE (defaults to failing with an error message).
#   The module flag is cached in the ax_cv_f90_modflag variable.
#   It may contain significant trailing whitespace.
#
#   Known flags:
#   gfortran: -Idir, -I dir (-M dir, -Mdir (deprecated), -Jdir for writing)
#   g95: -I dir (-fmod=dir for writing)
#   SUN: -Mdir, -M dir (-moddir=dir for writing;
#                       -Idir for includes is also searched)
#   HP: -Idir, -I dir (+moddir=dir for writing)
#   IBM: -Idir (-qmoddir=dir for writing)
#   Intel: -Idir -I dir (-mod dir for writing)
#   Absoft: -pdir
#   Lahey: -mod dir
#   Cray: -module dir, -p dir (-J dir for writing)
#         -e m is needed to enable writing .mod files at all
#   Compaq: -Idir
#   NAGWare: -I dir
#   PathScale: -I dir  (but -module dir is looked at first)
#   Portland: -module dir (first -module also names dir for writing)
#   Fujitsu: -Am -Idir (-Mdir for writing is searched first, then '.', then -I)
#                      (-Am indicates how module information is saved)
#
# LICENSE
#
# This file is part of Autoconf.  This program is free
# software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Under Section 7 of GPL version 3, you are granted additional
# permissions described in the Autoconf Configure Script Exception,
# version 3.0, as published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# and a copy of the Autoconf Configure Script Exception along with
# this program; see the files COPYINGv3 and COPYING.EXCEPTION
# respectively.  If not, see <http://www.gnu.org/licenses/>.
#
# Written by David MacKenzie, with help from
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.

AC_DEFUN([AX_F90_MODULE_FLAG],[
AC_CACHE_CHECK([Fortran 90 module inclusion flag], [ax_cv_f90_modflag],
[AC_LANG_PUSH([Fortran])
ax_cv_f90_modflag=unknown
mkdir conftest.dir
cd conftest.dir
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [cd ..
   ac_fc_module_flag_FCFLAGS_save=$FCFLAGS
   # Flag ordering is significant for gfortran and Sun.
   for ac_flag in -I '-I ' -M '-M ' -p '-mod ' '-module ' '-Am -I'; do
     # Add the flag twice to prevent matching an output flag.
     FCFLAGS="$ac_fc_module_flag_FCFLAGS_save ${ac_flag}conftest.dir ${ac_flag}conftest.dir"
     AC_COMPILE_IFELSE([[
      program main
      use conftest_module
      call conftest_routine
      end program]],
       [ax_cv_f90_modflag="$ac_flag"])
     if test "$ax_cv_f90_modflag" != unknown; then
       break
     fi
   done
   FCFLAGS=$ac_fc_module_flag_FCFLAGS_save
])
rm -rf conftest.dir
AC_LANG_POP([Fortran])
])
if test "$ax_cv_f90_modflag" != unknown; then
  FC_MODINC=$ax_cv_f90_modflag
  $1
else
  FC_MODINC=
  m4_default([$2],
    [AC_MSG_ERROR([unable to find compiler flag for module search path])])
fi
AC_SUBST([FC_MODINC])
# Ensure trailing whitespace is preserved in a Makefile.
AC_SUBST([ac_empty], [""])
AC_CONFIG_COMMANDS_PRE([case $FC_MODINC in #(
  *\ ) FC_MODINC=$FC_MODINC'${ac_empty}' ;;
esac])dnl
])
