noinst_LIBRARIES = libl1d.a

libl1d_a_LIBADD =
libl1d_a_LINK = $(F90LINK)

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

l1d_SOURCES = \
	EPTL.f90 \
	FORCE_MOM_PROC.f90 \
	GET_GRID_6X6_MASS.f90 \
	GRAV_PROC.f90 \
	MPC_PROC.f90 \
	PRESSURE_DATA_PROC.f90 \
	RBE2_PROC.f90 \
	RBE3_PROC.f90 \
	RFORCE_PROC.f90 \
	RIGID_ELEM_PROC.f90 \
	RSPLINE_PROC.f90 \
	SLOAD_PROC.f90 \
	TEMPERATURE_DATA_PROC.f90 \
	YS_ARRAY.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1d_f90s=$(l1d_SOURCES:.F90=.f90)

lib_MODULES=$(l1d_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1d_a_SOURCES = $(l1d_SOURCES)
