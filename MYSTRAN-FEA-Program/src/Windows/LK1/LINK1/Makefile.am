noinst_LIBRARIES = liblink1.a

liblink1_a_LIBADD =
liblink1_a_LINK = $(F90LINK)

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

link1_SOURCES = \
	ALLOCATE_EMS_ARRAYS.f90 \
	ALLOCATE_L1_MGG.f90 \
	ALLOCATE_STF_ARRAYS.f90 \
	ALLOCATE_TEMPLATE.f90 \
	DEALLOCATE_EMS_ARRAYS.f90 \
	DEALLOCATE_L1_MGG.f90 \
	DEALLOCATE_STF_ARRAYS.f90 \
	DEALLOCATE_TEMPLATE.f90 \
	LINK0.f90 \
	LINK1.f90 \
	LINK1_RESTART_DATA.f90 \
	PRINT_CONSTANTS_1.f90 \
	PRINT_ORDER.f90 \
	WRITE_ENF_TO_L1O.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

link1_f90s=$(link1_SOURCES:.F90=.f90)

lib_MODULES=$(link1_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblink1_a_SOURCES = $(link1_SOURCES)
