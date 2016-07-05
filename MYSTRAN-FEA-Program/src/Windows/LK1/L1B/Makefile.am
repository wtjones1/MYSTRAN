noinst_LIBRARIES = libl1b.a

libl1b_a_LIBADD =
libl1b_a_LINK = $(F90LINK)

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

l1b_SOURCES = \
	CORD_PROC.f90 \
	DOF_PROC.f90 \
	GRID_PROC.f90 \
	OU4_PARTVEC_PROC.f90 \
	RDOF.f90 \
	SEQ_PROC.f90 \
	TDOF_COL_NUM.f90 \
	TDOF_PROC.f90 \
	TSET_PROC.f90 \
	TSET_PROC_FOR_MPCS.f90 \
	TSET_PROC_FOR_OMITS.f90 \
	TSET_PROC_FOR_RIGELS.f90 \
	TSET_PROC_FOR_SPCS.f90 \
	TSET_PROC_FOR_SUPORTS.f90 \
	USET_PROC.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1b_f90s=$(l1b_SOURCES:.F90=.f90)

lib_MODULES=$(l1b_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1b_a_SOURCES = $(l1b_SOURCES)
