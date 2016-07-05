noinst_LIBRARIES = liblinuxdiffs.a

liblinuxdiffs_a_LIBADD =
liblinuxdiffs_a_LINK = $(F90LINK)

linuxdiffs_SOURCES = \
	GET_INI_FILNAM.f90 \
	GET_MYSTRAN_DIR.f90 \
	READ_CL.f90 \
	mkl_dtfi.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT)

linuxdiffs_f90s=$(linuxdiffs_SOURCES:.F90=.f90)

lib_MODULES=$(linuxdiffs_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblinuxdiffs_a_SOURCES = $(linuxdiffs_SOURCES)
