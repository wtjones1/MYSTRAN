noinst_LIBRARIES = libl92.a

libl92_a_LIBADD =
libl92_a_LINK = $(F90LINK)

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

l92_SOURCES = \
	CALC_ELEM_NODE_FORCES.f90 \
	CALC_ELEM_STRAINS.f90 \
	CALC_ELEM_STRESSES.f90 \
	ELEM_STRE_STRN_ARRAYS.f90 \
	ELMDIS.f90 \
	ELMDIS_PLY.f90 \
	GET_COMP_SHELL_ALLOWS.f90 \
	GP_FORCE_BALANCE_PROC.f90 \
	INDEP_FAILURE_INDEX.f90 \
	OFP1.f90 \
	OFP2.f90 \
	OFP3.f90 \
	OFP3_ELFE_1D.f90 \
	OFP3_ELFE_2D.f90 \
	OFP3_ELFN.f90 \
	OFP3_STRE_NO_PCOMP.f90 \
	OFP3_STRE_PCOMP.f90 \
	OFP3_STRN_NO_PCOMP.f90 \
	OFP3_STRN_PCOMP.f90 \
	ONE_D_STRAIN_OUTPUTS.f90 \
	ONE_D_STRESS_OUTPUTS.f90 \
	POLYNOM_FIT_STRE_STRN.f90 \
	POLY_FAILURE_INDEX.f90 \
	SHELL_ENGR_FORCE_OGEL.f90 \
	SHELL_STRAIN_OUTPUTS.f90 \
	SHELL_STRESS_OUTPUTS.f90 \
	SOLID_STRAIN_OUTPUTS.f90 \
	SOLID_STRESS_OUTPUTS.f90 \
	STR_TENSOR_TRANSFORM.f90 \
	SUSER1.f90 \
	TRANSFORM_NODE_FORCES.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l92_f90s=$(l92_SOURCES:.F90=.f90)

lib_MODULES=$(l92_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl92_a_SOURCES = $(l92_SOURCES)
