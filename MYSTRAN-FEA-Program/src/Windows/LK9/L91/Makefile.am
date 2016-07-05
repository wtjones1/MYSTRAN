noinst_LIBRARIES = libl91.a

libl91_a_LIBADD =
libl91_a_LINK = $(F90LINK)

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

l91_SOURCES = \
	BAR_MARGIN.f90 \
	GET_MAX_MIN_ABS_STR.f90 \
	PRINCIPAL_2D.f90 \
	PRINCIPAL_3D.f90 \
	ROD_MARGIN.f90 \
	WRITE_BAR.f90 \
	WRITE_ELEM_ENGR_FORCE.f90 \
	WRITE_ELEM_NODE_FORCE.f90 \
	WRITE_ELEM_STRAINS.f90 \
	WRITE_ELEM_STRESSES.f90 \
	WRITE_FEMAP_ELFO_VECS.f90 \
	WRITE_FEMAP_GRID_VECS.f90 \
	WRITE_FEMAP_STRE_VECS.f90 \
	WRITE_FEMAP_STRN_VECS.f90 \
	WRITE_GRID_OUTPUTS.f90 \
	WRITE_MEFFMASS.f90 \
	WRITE_MPFACTOR.f90 \
	WRITE_PLY_STRAINS.f90 \
	WRITE_PLY_STRESSES.f90 \
	WRITE_ROD.f90 \
	WRT_REAL_TO_CHAR_VAR.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l91_f90s=$(l91_SOURCES:.F90=.f90)

lib_MODULES=$(l91_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl91_a_SOURCES = $(l91_SOURCES)
