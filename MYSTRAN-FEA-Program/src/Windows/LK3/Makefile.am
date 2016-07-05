noinst_LIBRARIES = liblk3.a

liblk3_a_LIBADD =
liblk3_a_LINK = $(F90LINK)

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

lk3_SOURCES = \
	EPSCALC.f90 \
	LINK3.f90 \
	VECINORM.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

lk3_f90s=$(lk3_SOURCES:.F90=.f90)

lib_MODULES=$(lk3_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblk3_a_SOURCES = $(lk3_SOURCES)
