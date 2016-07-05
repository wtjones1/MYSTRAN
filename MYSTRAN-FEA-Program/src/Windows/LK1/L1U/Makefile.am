noinst_LIBRARIES = libl1u.a

libl1u_a_LIBADD =
libl1u_a_LINK = $(F90LINK)

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

l1u_SOURCES = \
	CHAR_FLD.f90 \
	CHECK_BAR_MOIs.f90 \
	CRDERR.f90 \
	GET_ANSID.f90 \
	GET_SETID.f90 \
	I4FLD.f90 \
	IP6CHK.f90 \
	LEFT_ADJ_BDFLD.f90 \
	MKCARD.f90 \
	MKJCARD.f90 \
	MKJCARD_08.f90 \
	NEXTC.f90 \
	NEXTC0.f90 \
	NEXTC2.f90 \
	NEXTC20.f90 \
	R8FLD.f90 \
	STOKEN.f90 \
	TOKCHK.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1u_f90s=$(l1u_SOURCES:.F90=.f90)

lib_MODULES=$(l1u_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1u_a_SOURCES = $(l1u_SOURCES)
