noinst_LIBRARIES = libemg4.a

libemg4_a_LIBADD =
libemg4_a_LINK = $(F90LINK)

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

emg4_SOURCES = \
	CALC_PHI_SQ.f90 \
	GET_ELEM_NUM_PLIES.f90 \
	QDEL1.f90 \
	QMEM1.f90 \
	QPLT1.f90 \
	QPLT2.f90 \
	QPLT3.f90 \
	QSHEAR.f90 \
	TMEM1.f90 \
	TPLT1.f90 \
	TPLT2.f90 \
	TREL1.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

emg4_f90s=$(emg4_SOURCES:.F90=.f90)

lib_MODULES=$(emg4_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libemg4_a_SOURCES = $(emg4_SOURCES)
