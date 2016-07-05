noinst_LIBRARIES = libemg6.a

libemg6_a_LIBADD =
libemg6_a_LINK = $(F90LINK)

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

emg6_SOURCES = \
	B3D_ISOPARAMETRIC.f90 \
	BBDKQ.f90 \
	BBMIN3.f90 \
	BBMIN4.f90 \
	BCHECK_2D.f90 \
	BCHECK_3D.f90 \
	BMQMEM.f90 \
	BSMIN3.f90 \
	BSMIN4.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

emg6_f90s=$(emg6_SOURCES:.F90=.f90)

lib_MODULES=$(emg6_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libemg6_a_SOURCES = $(emg6_SOURCES)
