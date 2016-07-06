bin_PROGRAMS = mystran

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

mystran_LDADD = \
	@top_builddir@/src/Windows/UTIL/libutil.a \
	@top_builddir@/src/Windows/LK9/LINK9/liblink9.a \
	@top_builddir@/src/Windows/LK9/L92/libl92.a \
	@top_builddir@/src/Windows/LK9/L91/libl91.a \
	@top_builddir@/src/Windows/LK6/liblk6.a \
	@top_builddir@/src/Windows/LK5/liblk5.a \
	@top_builddir@/src/Windows/LK4/liblk4.a \
	@top_builddir@/src/Windows/LK3/liblk3.a \
	@top_builddir@/src/Windows/LK2/liblk2.a \
	@top_builddir@/src/Windows/LK1/LINK1/liblink1.a \
	@top_builddir@/src/Windows/LK1/L1U/libl1u.a \
	@top_builddir@/src/Windows/LK1/L1E/libl1e.a \
	@top_builddir@/src/Windows/LK1/L1D/libl1d.a \
	@top_builddir@/src/Windows/LK1/L1C/libl1c.a \
	@top_builddir@/src/Windows/LK1/L1B/libl1b.a \
	@top_builddir@/src/Windows/LK1/L1A-CC/libl1acc.a \
	@top_builddir@/src/Windows/LK1/L1A-BD/libl1abd.a \
	@top_builddir@/src/Windows/LK1/L1A/libl1a.a \
	@top_builddir@/src/Windows/EMG/EMG8/libemg8.a \
	@top_builddir@/src/Windows/EMG/EMG7/libemg7.a \
	@top_builddir@/src/Windows/EMG/EMG6/libemg6.a \
	@top_builddir@/src/Windows/EMG/EMG5/libemg5.a \
	@top_builddir@/src/Windows/EMG/EMG4/libemg4.a \
	@top_builddir@/src/Windows/EMG/EMG3/libemg3.a \
	@top_builddir@/src/Windows/EMG/EMG2/libemg2.a \
	@top_builddir@/src/Windows/EMG/EMG1/libemg1.a \
	@top_builddir@/src/Windows/BANDIT/libbandit.a \
	@top_builddir@/src/Windows/ARPACK/libarpack.a \
	$(USEIFS_DIR)/libuseifs.a \
	$(INTERF_DIR)/libinterfaces.a \
	$(LAPACK_DIR)/liblapack.a \
	$(BANDIT_DIR)/libbandit.a \
	$(ARPACK_DIR)/libarpack.a \
	$(MODULE_DIR)/libmystran.a \
	@intel_mkl_ldadd@ \
	-lpthread

mystran_SOURCES = \
	MYSTRAN.f90 \
	MYSTRAN_FILES.f90 \
	PROCESS_INCLUDE_FILES.f90 \
	READ_INI.f90 \
	READ_INPUT_FILE_NAME.f90

mystran_SOURCES += \
	$(top_srcdir)/src/Linux_differences/GET_INI_FILNAM.f90 \
	$(top_srcdir)/src/Linux_differences/GET_MYSTRAN_DIR.f90

#mystran_SOURCES += \
#	GET_INI_FILNAM.f90 \
#	GET_MYSTRAN_DIR.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

mystran_f90s=$(mystran_SOURCES:.F90=.f90)

lib_MODULES=$(mystran_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules
