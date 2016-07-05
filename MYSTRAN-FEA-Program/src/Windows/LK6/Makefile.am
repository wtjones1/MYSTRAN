noinst_LIBRARIES = liblk6.a

liblk6_a_LIBADD =
liblk6_a_LINK = $(F90LINK)

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

lk6_SOURCES = \
	ALLOCATE_L6_2.f90 \
	CALC_CB_MEFM_MPF.f90 \
	CALC_KRRcb.f90 \
	CALC_MRN.f90 \
	CALC_MRRcb.f90 \
	CALC_PHIZL.f90 \
	DEALLOCATE_L6_2.f90 \
	INTERFACE_FORCE_LTM.f90 \
	LINK6.f90 \
	MERGE_KXX.f90 \
	MERGE_LTM.f90 \
	MERGE_MXX.f90 \
	MERGE_PHIXA.f90 \
	NET_CG_LOADS_LTM.f90 \
	SOLVE6_SETUP.f90 \
	SOLVE_DLR.f90 \
	SOLVE_PHIZL1.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

lk6_f90s=$(lk6_SOURCES:.F90=.f90)

lib_MODULES=$(lk6_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblk6_a_SOURCES = $(lk6_SOURCES)
