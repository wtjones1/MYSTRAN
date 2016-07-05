noinst_LIBRARIES = libl1e.a

libl1e_a_LIBADD =
libl1e_a_LINK = $(F90LINK)

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

l1e_SOURCES = \
	EMP.f90 \
	EMP0.f90 \
	ESP.f90 \
	ESP0.f90 \
	ESP0_FINAL.f90 \
	KGG_SINGULARITY_PROC.f90 \
	MGGC_MASS_MATRIX.f90 \
	MGGS_MASS_MATRIX.f90 \
	SPARSE_KGG.f90 \
	SPARSE_KGGD.f90 \
	SPARSE_MGG.f90 \
	SPARSE_PG.f90 \
	SPARSE_RMG.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1e_f90s=$(l1e_SOURCES:.F90=.f90)

lib_MODULES=$(l1e_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1e_a_SOURCES = $(l1e_SOURCES)
