noinst_LIBRARIES = libemg1.a

libemg1_a_LIBADD =
libemg1_a_LINK = $(F90LINK)

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

emg1_SOURCES = \
	CHECK_TE_MATRIX.f90 \
	ELMDAT1.f90 \
	ELMDAT2.f90 \
	ELMGM1.f90 \
	ELMGM2.f90 \
	ELMGM3.f90 \
	EMG.f90 \
	GET_ELEM_AGRID_BGRID.f90 \
	GET_ELEM_ONAME.f90 \
	GET_ELGP.f90 \
	GET_MATANGLE_FROM_CID.f90 \
	GET_PCOMP_SECT_PROPS.f90 \
	GRID_ELEM_CONN_TABLE.f90 \
	IS_ELEM_PCOMP_PROPS.f90 \
	ROT_COMP_ELEM_AXES.f90 \
	SHELL_ABD_MATRICES.f90 \
	SOLVE_SHELL_ALP.f90 \
	WRITE_PCOMP_EQUIV.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

emg1_f90s=$(emg1_SOURCES:.F90=.f90)

lib_MODULES=$(emg1_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libemg1_a_SOURCES = $(emg1_SOURCES)
