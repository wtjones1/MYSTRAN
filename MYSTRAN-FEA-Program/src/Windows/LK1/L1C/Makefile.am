noinst_LIBRARIES = libl1c.a

libl1c_a_LIBADD =
libl1c_a_LINK = $(F90LINK)

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

l1c_SOURCES = \
	CONM2_PROC_1.f90 \
	CONM2_PROC_2.f90 \
	ELEM_PROP_MATL_IIDS.f90 \
	ELEM_TRANSFORM_LBG.f90 \
	ELESORT.f90 \
	ELSAVE.f90 \
	GPWG.f90 \
	GPWG_PMOI.f90 \
	GPWG_USERIN.f90 \
	RB_DISP_MATRIX_PROC.f90 \
	SUBCASE_PROC.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1c_f90s=$(l1c_SOURCES:.F90=.f90)

lib_MODULES=$(l1c_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1c_a_SOURCES = $(l1c_SOURCES)
