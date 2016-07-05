noinst_LIBRARIES = libemg3.a

libemg3_a_LIBADD =
libemg3_a_LINK = $(F90LINK)

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

emg3_SOURCES = \
	BAR1.f90 \
	BEAM.f90 \
	BREL1.f90 \
	BUSH.f90 \
	ELAS1.f90 \
	KUSER1.f90 \
	PINFLG.f90 \
	ROD1.f90 \
	USERIN.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

emg3_f90s=$(emg3_SOURCES:.F90=.f90)

lib_MODULES=$(emg3_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libemg3_a_SOURCES = $(emg3_SOURCES)
