# ===========================================================================
#  extracted from http://git.savannah.gnu.org/gitweb/?p=autoconf.git;a=blob_plain;f=lib/autoconf/fortran.m4;hb=3e3412d3dd1ff7734f0a9ccd1a55252bd7ea8790
# ===========================================================================
#
# SYNOPSIS
#
#   AX_F90_MODULE_EXTENSION
#
# DESCRIPTION
#
#   Find the Fortran 90 module file extension.  The module extension is stored
#   in the variable FC_MODEXT and empty if it cannot be determined.  The result
#   or "unknown" is cached in the cache variable ax_cv_f90_modext.
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
#

AC_DEFUN([AX_F90_MODULE_EXTENSION],
[AC_CACHE_CHECK([Fortran 90 module extension], [ax_cv_f90_modext],
[AC_LANG_PUSH(Fortran)
mkdir conftest.dir
cd conftest.dir
ax_cv_f90_modext=unknown
AC_COMPILE_IFELSE([[
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'gotcha!'
      end subroutine
      end module]],
  [ax_cv_f90_modext=`ls | sed -n 's,conftest_module\.,,p'`
   if test x$ax_cv_f90_modext = x; then
dnl Some F90 compilers use upper case characters for the module file name.
     ax_cv_f90_modext=`ls | sed -n 's,CONFTEST_MODULE\.,,p'`
   fi])
cd ..
rm -rf conftest.dir
AC_LANG_POP(Fortran)
])
FC_MODEXT=$ax_cv_f90_modext
if test "$FC_MODEXT" = unknown; then
  FC_MODEXT=
fi
AC_SUBST([FC_MODEXT])dnl
])
