noinst_LIBRARIES = libl1acc.a

libl1acc_a_LIBADD =
libl1acc_a_LINK = $(F90LINK)

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

l1acc_SOURCES = \
	CC_ACCE.f90 \
	CC_DISP.f90 \
	CC_ECHO.f90 \
	CC_ELDA.f90 \
	CC_ELFO.f90 \
	CC_ENFO.f90 \
	CC_GPFO.f90 \
	CC_LABE.f90 \
	CC_LOAD.f90 \
	CC_METH.f90 \
	CC_MPC.f90 \
	CC_MPCF.f90 \
	CC_NLPARM.f90 \
	CC_OLOA.f90 \
	CC_OUTPUTS.f90 \
	CC_SET.f90 \
	CC_SET0.f90 \
	CC_SPC.f90 \
	CC_SPCF.f90 \
	CC_STRE.f90 \
	CC_STRN.f90 \
	CC_SUBC.f90 \
	CC_SUBT.f90 \
	CC_TEMP.f90 \
	CC_TITL.f90 \
	CHK_CC_CMD_DESCRIBERS.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1acc_f90s=$(l1acc_SOURCES:.F90=.f90)

lib_MODULES=$(l1acc_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1acc_a_SOURCES = $(l1acc_SOURCES)
