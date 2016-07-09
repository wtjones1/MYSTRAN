noinst_LIBRARIES = liblinuxdiffs.a

liblinuxdiffs_a_LIBADD =
liblinuxdiffs_a_LINK = $(F90LINK)

MODULE_DIR=@top_builddir@/src/Windows/Modules
ARPACK_DIR=@top_builddir@/src/Windows/Modules/ARPACK
BANDIT_DIR=@top_builddir@/src/Windows/Modules/BANDIT
LAPACK_DIR=@top_builddir@/src/Windows/Modules/LAPACK
INTERF_DIR=@top_builddir@/src/Interfaces
USEIFS_DIR=@top_builddir@/src/USE_IFs

AM_FCFLAGS = \
        $(FC_MODINC)$(MODULE_DIR) \
        $(FC_MODINC)$(ARPACK_DIR) \
        $(FC_MODINC)$(BANDIT_DIR) \
        $(FC_MODINC)$(LAPACK_DIR) \
        $(FC_MODINC)$(INTERF_DIR) \
        $(FC_MODINC)$(USEIFS_DIR)

linuxdiffs_SOURCES = \
	GET_INI_FILNAM.f90 \
	GET_MYSTRAN_DIR.f90 \
	READ_CL.f90 \
	mkl_dtfi.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT)

linuxdiffs_f90s=$(linuxdiffs_SOURCES:.F90=.f90)

lib_MODULES=$(linuxdiffs_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblinuxdiffs_a_SOURCES = $(linuxdiffs_SOURCES)
