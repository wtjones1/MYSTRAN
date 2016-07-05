noinst_LIBRARIES = liblink9.a

liblink9_a_LIBADD =
liblink9_a_LINK = $(F90LINK)

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

link9_SOURCES = \
	ALLOCATE_FEMAP_DATA.f90 \
	ALLOCATE_LINK9_STUF.f90 \
	DEALLOCATE_FEMAP_DATA.f90 \
	DEALLOCATE_LINK9_STUF.f90 \
	LINK9.f90 \
	LINK9S.f90 \
	MAXREQ_OGEL.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

link9_f90s=$(link9_SOURCES:.F90=.f90)

lib_MODULES=$(link9_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblink9_a_SOURCES = $(link9_SOURCES)
