noinst_LIBRARIES = libemg8.a

libemg8_a_LIBADD =
libemg8_a_LINK = $(F90LINK)

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

emg8_SOURCES = \
	JAC2D.f90 \
	JAC3D.f90 \
	MATERIAL_PROPS_2D.f90 \
	MATERIAL_PROPS_3D.f90 \
	MATGET.f90 \
	MATPUT.f90 \
	ROT_AXES_MATL_TO_LOC.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

emg8_f90s=$(emg8_SOURCES:.F90=.f90)

lib_MODULES=$(emg8_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libemg8_a_SOURCES = $(emg8_SOURCES)
