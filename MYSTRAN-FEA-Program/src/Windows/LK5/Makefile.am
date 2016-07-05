noinst_LIBRARIES = liblk5.a

liblk5_a_LIBADD =
liblk5_a_LINK = $(F90LINK)

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

lk5_SOURCES = \
	BUILD_A_LR.f90 \
	BUILD_F_AO.f90 \
	BUILD_G_NM.f90 \
	BUILD_N_FS.f90 \
	EXPAND_PHIXA_TO_PHIXG.f90 \
	LINK5.f90 \
	RENORM.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

lk5_f90s=$(lk5_SOURCES:.F90=.f90)

lib_MODULES=$(lk5_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblk5_a_SOURCES = $(lk5_SOURCES)
