noinst_LIBRARIES = liblapack.a

liblapack_a_LIBADD =
liblapack_a_LINK = $(F90LINK)

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

lapack_SOURCES = \
 	$(top_srcdir)/src/Interfaces/OUTA_HERE_Interface.f90 \
 	$(top_srcdir)/src/Interfaces/OURTIM_Interface.f90 \
	LAPACK_BLAS_AUX.f90 \
	LAPACK_LANCZOS_EIG.f90 \
	LAPACK_LIN_EQN_DGB.f90 \
	LAPACK_LIN_EQN_DPB.f90 \
	LAPACK_SYM_MAT_INV.f90 \
	LAPACK_LIN_EQN_DGE.f90 \
	LAPACK_MISCEL.f90 \
	LAPACK_GIV_MGIV_EIG.f90 \
	LAPACK_STD_EIG_1.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

lapack_f90s=$(lapack_SOURCES:.F90=.f90)

lib_MODULES=$(lapack_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblapack_a_SOURCES = $(lapack_SOURCES)
