noinst_LIBRARIES = libemg7.a

libemg7_a_LIBADD =
libemg7_a_LINK = $(F90LINK)

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

emg7_SOURCES = \
	MIN4SH.f90 \
	ORDER_GAUSS.f90 \
	ORDER_TETRA.f90 \
	ORDER_TRIA.f90 \
	SHP2DQ.f90 \
	SHP3DH.f90 \
	SHP3DP.f90 \
	SHP3DT.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

emg7_f90s=$(emg7_SOURCES:.F90=.f90)

lib_MODULES=$(emg7_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libemg7_a_SOURCES = $(emg7_SOURCES)
