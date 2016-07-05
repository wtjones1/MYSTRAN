noinst_LIBRARIES = libbandit.a

libbandit_a_LIBADD =
libbandit_a_LINK = $(F90LINK)

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

bandit_SOURCES = \
	BANDIT_BLOCK_DATA.f90 \
	BANDIT_FILES.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

bandit_f90s=$(bandit_SOURCES:.F90=.f90)

lib_MODULES=$(bandit_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libbandit_a_SOURCES = $(bandit_SOURCES)
