noinst_LIBRARIES = libemg2.a

libemg2_a_LIBADD =
libemg2_a_LINK = $(F90LINK)

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

emg2_SOURCES = \
	ELMOFF.f90 \
	ELMOUT.f90 \
	ELMTLB.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

emg2_f90s=$(emg2_SOURCES:.F90=.f90)

lib_MODULES=$(emg2_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libemg2_a_SOURCES = $(emg2_SOURCES)
