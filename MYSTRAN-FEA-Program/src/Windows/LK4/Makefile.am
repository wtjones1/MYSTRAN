noinst_LIBRARIES = liblk4.a

liblk4_a_LIBADD =
liblk4_a_LINK = $(F90LINK)

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

lk4_SOURCES = \
	CALC_GEN_MASS.f90 \
	EIG_GIV_MGIV.f90 \
	EIG_INV_PWR.f90 \
	EIG_LANCZOS_ARPACK.f90 \
	EIG_SUMMARY.f90 \
	INVERT_EIGENS.f90 \
	LINK4.f90 \
	RENORM_ON_MASS.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

lk4_f90s=$(lk4_SOURCES:.F90=.f90)

lib_MODULES=$(lk4_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblk4_a_SOURCES = $(lk4_SOURCES)
