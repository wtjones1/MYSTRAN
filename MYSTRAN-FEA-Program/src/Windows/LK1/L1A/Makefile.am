noinst_LIBRARIES = libl1a.a

libl1a_a_LIBADD =
libl1a_a_LINK = $(F90LINK)

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

l1a_SOURCES = \
	CSHIFT.f90 \
	EC_DEBUG.f90 \
	EC_IN4FIL.f90 \
	EC_OUTPUT4.f90 \
	EC_PARTN.f90 \
	ELEPRO.f90 \
	FFIELD.f90 \
	FFIELD2.f90 \
	IS_THIS_A_RESTART.f90 \
	LOADB.f90 \
	LOADB0.f90 \
	LOADB_RESTART.f90 \
	LOADC.f90 \
	LOADC0.f90 \
	LOADE.f90 \
	LOADE0.f90 \
	READ_INCLUDE_FILNAM.f90 \
	REPLACE_TABS_W_BLANKS.f90 \
	RW_INCLUDE_FILES.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1a_f90s=$(l1a_SOURCES:.F90=.f90)

lib_MODULES=$(l1a_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1a_a_SOURCES = $(l1a_SOURCES)
