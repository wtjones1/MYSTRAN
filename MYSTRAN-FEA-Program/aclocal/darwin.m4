# autoconf macros for checking to see if running Mac OSX intel ifort compiler
# 
# This will turn off vectorization for loops that go bad due to a bug in
# the ifort compiler for at least all versions <= 15.0.3
#
# Assigned AC_DEFINE:
# DARWIN_VEC_OFF
# DARWIN_VEC_WARN
#

AC_DEFUN([CHECK_IFORT_DARWIN],[

compiler_base=`basename $FC`
if test -n "${MPIF90}" -a -n "`echo ${FC}x | grep ${MPIF90}x`"
then
  compiler_base=`$FC -show | grep -v "^ln" | grep -v "^rm" | grep -v "compchk.sh"`
  compiler_base=`echo "$compiler_base" | sed 's/ .*//'`
  compiler_base=`basename "$compiler_base"`
  compiler_base=`echo $compiler_base | sed 's/ .*//'`
else
  compiler_base=`echo $compiler_base | sed 's/ .*//'`
fi
AC_MSG_RESULT([  Check for Mac OSX ifort compiler bug  ])

case $compiler_base in
  ifort)

    if (uname -s | grep Darwin) > /dev/null 2>&1
    then
      current_version=`ifort -v 2>&1 >/dev/null | cut -d' ' -f3 | cut -c1-2`
      vec_bug_version=15
      if test ${current_version} -le ${vec_bug_version}
      then
        AC_DEFINE([DARWIN_VEC_OFF],[1],[Turn off vectorization of problem loop])
      else
        AC_DEFINE([DARWIN_VEC_WARN],[1],[Warn that vectorization is on])
      fi
    fi

    ;;

esac

])
