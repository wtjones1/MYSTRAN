noinst_LIBRARIES = libl1abd.a

libl1abd_a_LIBADD =
libl1abd_a_LINK = $(F90LINK)

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

l1abd_SOURCES = \
	BD_ASET.f90 \
	BD_ASET1.f90 \
	BD_BAROR.f90 \
	BD_BAROR0.f90 \
	BD_BEAMOR.f90 \
	BD_BEAMOR0.f90 \
	BD_CBAR.f90 \
	BD_CBAR0.f90 \
	BD_CBUSH.f90 \
	BD_CBUSH0.f90 \
	BD_CELAS1.f90 \
	BD_CELAS2.f90 \
	BD_CELAS3.f90 \
	BD_CELAS4.f90 \
	BD_CHEXA.f90 \
	BD_CHEXA0.f90 \
	BD_CMASS1.f90 \
	BD_CMASS2.f90 \
	BD_CMASS3.f90 \
	BD_CMASS4.f90 \
	BD_CONM2.f90 \
	BD_CONROD.f90 \
	BD_CORD.f90 \
	BD_CPENTA.f90 \
	BD_CPENTA0.f90 \
	BD_CQUAD.f90 \
	BD_CQUAD0.f90 \
	BD_CROD.f90 \
	BD_CSHEAR.f90 \
	BD_CTETRA.f90 \
	BD_CTETRA0.f90 \
	BD_CTRIA.f90 \
	BD_CTRIA0.f90 \
	BD_CUSER1.f90 \
	BD_CUSERIN.f90 \
	BD_CUSERIN0.f90 \
	BD_DEBUG.f90 \
	BD_DEBUG0.f90 \
	BD_EIGR.f90 \
	BD_EIGRL.f90 \
	BD_FORMOM.f90 \
	BD_GRAV.f90 \
	BD_GRDSET.f90 \
	BD_GRDSET0.f90 \
	BD_GRID.f90 \
	BD_IMBEDDED_BLANK.f90 \
	BD_LOAD.f90 \
	BD_LOAD0.f90 \
	BD_MAT1.f90 \
	BD_MAT2.f90 \
	BD_MAT8.f90 \
	BD_MAT9.f90 \
	BD_MPC.f90 \
	BD_MPC0.f90 \
	BD_MPCADD.f90 \
	BD_MPCADD0.f90 \
	BD_NLPARM.f90 \
	BD_PARAM.f90 \
	BD_PARAM0.f90 \
	BD_PARVEC.f90 \
	BD_PARVEC1.f90 \
	BD_PBAR.f90 \
	BD_PBARL.f90 \
	BD_PBEAM.f90 \
	BD_PBUSH.f90 \
	BD_PCOMP.f90 \
	BD_PCOMP0.f90 \
	BD_PCOMP1.f90 \
	BD_PCOMP10.f90 \
	BD_PELAS.f90 \
	BD_PLOAD2.f90 \
	BD_PLOAD4.f90 \
	BD_PLOTEL.f90 \
	BD_PMASS.f90 \
	BD_PROD.f90 \
	BD_PSHEAR.f90 \
	BD_PSHEL.f90 \
	BD_PSOLID.f90 \
	BD_PUSER1.f90 \
	BD_PUSERIN.f90 \
	BD_RBAR.f90 \
	BD_RBE1.f90 \
	BD_RBE2.f90 \
	BD_RBE3.f90 \
	BD_RBE30.f90 \
	BD_RFORCE.f90 \
	BD_RSPLINE.f90 \
	BD_RSPLINE0.f90 \
	BD_SEQGP.f90 \
	BD_SLOAD.f90 \
	BD_SLOAD0.f90 \
	BD_SPC.f90 \
	BD_SPC1.f90 \
	BD_SPCADD.f90 \
	BD_SPCADD0.f90 \
	BD_SPOINT.f90 \
	BD_SPOINT0.f90 \
	BD_SUPORT.f90 \
	BD_TEMP.f90 \
	BD_TEMPD.f90 \
	BD_TEMPRP.f90 \
	BD_USET.f90 \
	BD_USET1.f90

# remove *.mod *.interface when mod_suffix is repaired for OS X
CLEANFILES = *.$(FC_MODEXT) *.interface *.d *.fh

l1abd_f90s=$(l1abd_SOURCES:.F90=.f90)

lib_MODULES=$(l1abd_f90s:.f90=.$(FC_MODEXT))

include $(top_srcdir)/make.rules

libl1abd_a_SOURCES = $(l1abd_SOURCES)
