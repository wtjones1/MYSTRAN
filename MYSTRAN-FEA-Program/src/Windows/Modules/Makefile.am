SUBDIRS = . LAPACK ARPACK BANDIT

noinst_LIBRARIES = libmystran.a libnocommon.a

libmystran_a_LIBADD =
libmystran_a_LINK = $(F90LINK)

MODULE_DIR=@top_builddir@/src/Windows/Modules
ARPACK_DIR=@top_builddir@/src/Windows/Modules/ARPACK
BANDIT_DIR=@top_builddir@/src/Windows/Modules/BANDIT
LAPACK_DIR=@top_builddir@/src/Windows/Modules/LAPACK

AM_FCFLAGS = \
	$(FC_MODINC)$(MODULE_DIR) \
	$(FC_MODINC)$(ARPACK_DIR) \
	$(FC_MODINC)$(BANDIT_DIR) \
	$(FC_MODINC)$(LAPACK_DIR)

mystran_SOURCES = \
	PENTIUM_II_KIND.f90 \
	ARPACK_MATRICES_1.f90 \
	CONSTANTS_1.f90 \
	ALLOCATED_ARRAY_DATA.f90 \
	SCONTR.f90 \
	CC_OUTPUT_DESCRIBERS.f90 \
	COL_VECS.f90 \
	CONSTANTS_GAUSS.f90 \
	DEBUG_PARAMETERS.f90 \
	DERIVED_DATA_TYPES.f90 \
	DOF_TABLES.f90 \
	EIGEN_MATRICES_1.f90 \
	EMS_ARRAYS.f90 \
	FEMAP_ARRAYS.f90 \
	FULL_MATRICES.f90 \
	INPUTT4_MATRICES.f90 \
	IOUNT1.f90 \
	KEY_VALUES.f90 \
	LAPACK_DPB_MATRICES.f90 \
	LINK9_STUFF.f90 \
	LSQ_MYSTRAN.f90 \
	MISC_MATRICES.f90 \
	MKL_DSS_PRIVATE.f90 \
	MKL_DSS.f90 \
	MODEL_STUF.f90 \
	MYSTRAN_Version.f90 \
	NONLINEAR_PARAMS.f90 \
	OUTPUT4_MATRICES.f90 \
	PARAMS.f90 \
	RIGID_BODY_DISP_MATS.f90 \
	SCRATCH_MATRICES.f90 \
	SPARSE_ALG_ARRAYS.f90 \
	SPARSE_MATRICES.f90 \
	STF_ARRAYS.f90 \
	STF_TEMPLATE_ARRAYS.f90 \
	SUBR_BEGEND_LEVELS.f90

libnocommon_a_LIBADD =
libnocommon_a_LINK = $(F90LINK)

libnocommon_a_FCFLAGS = @ax_cv_no_common@

nocommon_SOURCES = \
	MACHINE_PARAMS.f90 \
	SP_ENVIRON_VARS.f90 \
	TIMDAT.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

mystran_f90s=$(mystran_SOURCES:.F90=.f90)
nocommon_f90s=$(nocommon_SOURCES:.F90=.f90)

lib_MODULES=$(mystran_f90s:.f90=.$(FC_MODEXT))
lib_MODULES+=$(nocommon_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libmystran_a_SOURCES = $(mystran_SOURCES)
libnocommon_a_SOURCES = $(nocommon_SOURCES)
