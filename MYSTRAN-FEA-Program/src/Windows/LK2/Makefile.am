noinst_LIBRARIES = liblk2.a

liblk2_a_LIBADD =
liblk2_a_LINK = $(F90LINK)

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

lk2_SOURCES = \
	ALLOCATE_L2_GMN_2.f90 \
	ALLOCATE_L2_GOA_2.f90 \
	DEALLOCATE_L2_GMN_2.f90 \
	DEALLOCATE_L2_GOA_2.f90 \
	LINK2.f90 \
	REDUCE_A_LR.f90 \
	REDUCE_F_AO.f90 \
	REDUCE_G_NM.f90 \
	REDUCE_KAAD_TO_KLLD.f90 \
	REDUCE_KAA_TO_KLL.f90 \
	REDUCE_KFFD_TO_KAAD.f90 \
	REDUCE_KFF_TO_KAA.f90 \
	REDUCE_KGGD_TO_KNND.f90 \
	REDUCE_KGG_TO_KNN.f90 \
	REDUCE_KNND_TO_KFFD.f90 \
	REDUCE_KNN_TO_KFF.f90 \
	REDUCE_MAA_TO_MLL.f90 \
	REDUCE_MFF_TO_MAA.f90 \
	REDUCE_MGG_TO_MNN.f90 \
	REDUCE_MNN_TO_MFF.f90 \
	REDUCE_N_FS.f90 \
	REDUCE_PA_TO_PL.f90 \
	REDUCE_PF_TO_PA.f90 \
	REDUCE_PG_TO_PN.f90 \
	REDUCE_PN_TO_PF.f90 \
	SOLVE_GMN.f90 \
	SOLVE_GOA.f90 \
	SOLVE_UO0.f90 \
	STIFF_MAT_EQUIL_CHK.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

lk2_f90s=$(lk2_SOURCES:.F90=.f90)

lib_MODULES=$(lk2_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

liblk2_a_SOURCES = $(lk2_SOURCES)
