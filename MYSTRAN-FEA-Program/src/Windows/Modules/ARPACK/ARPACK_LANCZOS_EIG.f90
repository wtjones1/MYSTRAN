! ################################################################################################################################## 

      MODULE ARPACK_LANCZOS_EIG

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F04, F06, SC1, WRT_LOG 
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EIG_MSGLVL
      USE SUBR_BEGEND_LEVELS, ONLY    :  ARPACK_BEGEND
      USE ARPACK_UTIL
      USE LAPACK_BLAS_AUX
      USE LAPACK_LANCZOS_EIG
      USE LAPACK_MISCEL                                    ! This contains DSTEQR, used in this module
      USE LAPACK_LIN_EQN_DGB
      USE LAPACK_LIN_EQN_DPB

      USE OURTIM_Interface
      USE MATMULT_SFF_Interface
      USE ARPACK_INFO_MSG_Interface

      CHARACTER(44*BYTE)             :: MODNAM1            ! Name to write to screen to describe module being run.
      CHARACTER(44*BYTE)             :: MODNAM2            ! Name to write to screen to describe module being run.

      INTEGER(LONG), PARAMETER, PRIVATE :: SUBR_BEGEND = ARPACK_BEGEND

! This is the set of ARPACK routines that are used in the Lanczos algorithm. Below are listed the subroutines included in this
! module and the calls to other subroutines in this module

!     dsband: the driver for Lanczos. dsband calls:

!             dsaupd: reverse communication interface, which calls:

!                     dstats: Initialize statistic and timing information

!                     dsaup2: Intermediate level interface called by dsaupd, which calls

!                             dgetv0: Generate a random initial residual vector for the Arnoldi process. Force the residual vector
!                                     to be in the range of the operator OP.  

!                             dsaitr: Reverse comm I/F for applying NP additional steps to a K step symmetric Arnoldi factorization.

!                                     dgetv0: see above

!                             dseigt: 
!                                     dstqrb:
 
!                             dsgets: Given the eigenvalues of the symmetric tridiagonal matrix H, computes the NP shifts AMU that
!                                     are zeros of the polynomial of  degree NP which filters out components of the unwanted
!                                     eigenvectors  corresponding to the AMU's based on some given criteria.

!                                     dsortr: Sort an array
 
!                             dsconv: Convergence testing for the symmetric Arnoldi eigenvalue routine.

!                             dsortr: Sort an array

!                             dsapps: Given the Arnoldi factorization
!                                     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!                                     apply NP shifts implicitly resulting in
!                                     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!                                     where Q is an orthogonal matrix of order KEV+NP. Q is the product of rotations resulting from
!                                     the NP bulge chasing sweeps.  The updated Arnoldi factorization becomes:
!                                     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.


!             dseupd: returns the converged approximations to eigenvalues of A*z = lambda*B*z and (optionally):

!                     dsgets: see above

!                     dsteqr: see above

!                     dsesrt: Sort an array

!                     dsortr: see above

      CONTAINS

! ################################################################################################################################## 

! \BeginDoc
!
! \Name: dsband
!
! \Description:
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal (Lanczos) basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  Matrices A and B are stored in LAPACK-style band form.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  (Lanczos) basis is always computed.  There is an additional storage cost 
!  of n*nev if both are requested (in this case a separate array Z must be 
!  supplied).
!
!  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
!  are called Ritz values and Ritz vectors respectively.  They are referred 
!  to as such in the comments that follow.  The computed orthonormal basis 
!  for the invariant subspace corresponding to these Ritz values is referred 
!  to as a Lanczos basis.
!
!  dsband can be called with one of the following modes:
!
!  Mode 2:  AB*x = lambda*MB*x, AB symmetric, MB symmetric positive definite
!           ===>  OP = inv[MB]*AB
!                 NOTE: ARPACK's implementation of this was intended for:
!                       AB = K (stiffness) and MB = M (mass) which required
!                       M symm pos def which it generally is not

! **********************************************************************************************************************************
! **********************************************************************************************************************************
! **                                                              |                                                               **
! **      N O R M A L   E I G E N V A L U E   P R O B L E M       |                      B U C K L I N G                          **
! **      -------------------------------------------------       |                      ---------------                          **
! ** For mode 2, I changed the roles of AB and MB and added a     | For mode 2, I changed the roles of AB and MB and added a      **
! ** shift so that:                                               | shift so that (KD is differential stiffness matrix):          **
! **                                                              |                                                               **
! **            AB = M                                            |            AB = -KD                                           **
! **            MB = [K - sigma*M]                                |            MB = [K + sigma*KD]                                **
! ** then                                                         | then                                                          **
! **            OP = (inv[K - sigma*M])*M                         |            OP = -(inv[K + sigma*KD])*KD                       **
! **                                                              |                                                               **
! ** The generalized eigenvalue problem is:                       | The generalized eigenvalue problem is:                        **
! **                                                              |                                                               **
! **                [K - lam*M]*x = 0                      (1)    |                [K + lam*KD]*x = 0                      (1)    **
! ** define                                                       | define                                                        **
! **               lam = lam_bar + sigma                          |               lam = lam_bar + sigma                           **
! ** Then (1) is:                                                 | Then (1) is:                                                  **
! **            [K - (lam_bar + sigma)*M]x = 0                    |            [K + (lam_bar + sigma)*KD]x = 0                    **
! ** or                                                           | or                                                            **
! **             [K - sigma*M]x = lam_bar*Mx               (2)    |             [K + sigma*KD]x = -lam_bar*KDx             (2)    **
! ** define                                                       | define                                                        **
! **                 mu = 1/lam_bar                               |                 mu = 1/lam_bar                                **
! ** so that (2) is:                                              | so that (2) is:                                               **
! **             mu*x = (inv[K - sigma*M])*M*x                    |            -mu*x = (inv[K + sigma*KD])*KD*x                   **
! ** or                                                           | or                                                            **
! **                   A*x = mu*x                          (3)    |                   A*x = mu*x                           (3)    **
! ** where                                                        | where                                                         **
! **              A = (inv[K - sigma*M])*M                        |              A = -(inv[K + sigma*KD])*KD                      **
! **                                                              |                                                               **
! ** eqn (3) is the standard eigenvalue problem solved where      | eqn (3) is the standard eigenvalue problem solved where       **
! ** the eigenvalues, lam, of the gen eigen problem are:          | the eigenvalues, lam, of the gen eigen problem are:           **
! **                                                              |                                                               **
! **                lam = sigma + 1/mu                            |                lam = sigma + 1/mu                             **
! **                                                              |                                                               **
! ** The calc of lambda = sigma + 1/mu is done in my subr         | The calc of lambda = sigma + 1/mu is done in my subr          **
! ** INVERT_EIGENS after dsband returns to my subr EIG_LANCZOS    | INVERT_EIGENS after dsband returns to my subr EIG_LANCZOS     **
! **                                                              |                                                               **
! **********************************************************************************************************************************
! **********************************************************************************************************************************


!
!  Mode 3:  AB*x = lambda*MB*x, AB symmetric, MB symmetric positive semi-definite
!           ===>  OP = (inv[K - sigma*M])*M  and  B = M.
!           ===>  Shift-and-Invert mode
!                 NOTE: dseupd does the calc: lambda = sigma + 1/mu   
!
!  The choice of mode must be specified in IPARAM(7) defined below.
!
! \Usage
!     CALL DSBAND ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, 
!    &            N, LDA, RFAC, KL, KU, WHICH, BMAT, NEV, 
!    &            TOL, RESID, NCV, V, LDV, IPARAM, WORKD, WORKL, 
!    &            LWORKL, IWORK, INFO, INFO_LAPACK, DTBSV_MSG, PITERS)
!
! \Arguments
!
!  RVEC    Logical (INPUT)
!          Specifies whether Ritz vectors corresponding to the Ritz value 
!          approximations to the eigenproblem A*z = lambda*B*z are computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute the associated Ritz vectors. 
!
!  HOWMNY  Character*1  (INPUT) 
!          Specifies how many Ritz vectors are wanted and the form of Z
!          the matrix of Ritz vectors. See remark 1 below.
!          = 'A': compute all Ritz vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
!          If HOWMNY = 'A' , SELECT is not referenced.
!
!  D       REAL(DOUBLE) array of dimension NEV.  (OUTPUT)
!          On exit, D contains the Ritz value approximations to the
!          eigenvalues of A*z = lambda*B*z. The values are returned
!          in ascending order. If IPARAM(7) = 3,4,5 then D represents
!          the Ritz values of OP computed by dsaupd transformed to
!          those of the original eigensystem A*z = lambda*B*z. If 
!          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
!          as the those of A*z = lambda*B*z.
!
!  Z       REAL(DOUBLE) N by NEV array if HOWMNY = 'A'.  (OUTPUT)
!          On exit, Z contains the B-orthonormal Ritz vectors of the
!          eigensystem A*z = lambda*B*z corresponding to the Ritz
!          value approximations.
!
!          If  RVEC = .FALSE. then Z is not referenced.
!          NOTE: The array Z may be set equal to first NEV columns of the 
!          Lanczos basis array V computed by DSAUPD.
!
!  LDZ     Integer.  (INPUT) 
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
!
!  SIGMA   REAL(DOUBLE)  (INPUT)
!          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
!          IPARAM(7) = 1 or 2.
! 
!  N       Integer.  (INPUT) 
!          Dimension of the eigenproblem.  
! 
!  AB      REAL(DOUBLE) array of dimension LDA by N. (INPUT)
!          The matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!  MB      REAL(DOUBLE) array of dimension LDA by N. (INPUT)
!          The matrix M in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set. 
!          The j-th column of M is stored in the j-th column of the
!          array AB as follows:
!          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!          Not referenced if IPARAM(7) = 1
!
!  LDA     Integer. (INPUT)
!          Leading dimension of AB, MB, RFAC.
!
!  RFAC    REAL(DOUBLE) array of LDA by N. (WORKSPACE/OUTPUT)
!          RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
!          is invoked.  It is used to store the LU factors of
!          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
!          It is not referenced when IPARAM(7) = 1.
!
!  KL      Integer. (INPUT)
!          Max(number of subdiagonals of A, number of subdiagonals of M)
!
!  KU      Integer. (OUTPUT)
!          Max(number of superdiagonals of A, number of superdiagonals of M)
!
!  WHICH   Character*2.  (INPUT)
!          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
!          the following.
!  
!            'LM' -> want the NEV eigenvalues of largest magnitude.
!            'SM' -> want the NEV eigenvalues of smallest magnitude.
!            'LA' -> want the NEV eigenvalues of largest REAL part.
!            'SA' -> want the NEV eigenvalues of smallest REAL part.
!            'BE' -> Compute NEV eigenvalues, half from each end of the 
!                    spectrum.  When NEV is odd, compute one more from 
!                    the high end than from the low end. 
!
!          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only. 
!          
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
!          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

!  NEV     Integer. (INPUT)
!          Number of eigenvalues of OP to be computed.
!   
!  TOL     REAL(DOUBLE) scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value 
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
!          If TOL .LE. 0. is passed a default is set:
!          DEFAULT = DLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine DLAMCH).
!
!  RESID   REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
!          On INPUT:
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V (less than or equal to N).
!          Represents the dimension of the Lanczos basis constructed
!          by dsaupd for OP.
!
!  V       REAL(DOUBLE) array N by NCV.  (OUTPUT)
!          Upon INPUT: the NCV columns of V contain the Lanczos basis 
!                      vectors as constructed by dsaupd for OP.
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
!                       represent the Ritz vectors that span the desired 
!                       invariant subspace.
!          NOTE: The array Z may be set equal to first NEV columns of the 
!          Lanczos basis vector array V computed by dsaupd. In this case
!          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
!          the desired Ritz vectors. 
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: 
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          It is set to 1 in this subroutine.  The user do not need
!          to set this parameter.
!          ------------------------------------------------------------
!          ISHIFT = 1: exact shifts with respect to the reduced 
!                      tridiagonal matrix T.  This is equivalent to 
!                      restarting the iteration with a starting vector 
!                      that is a linear combination of Ritz vectors 
!                      associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced. 
!
!          IPARAM(3) = MXITER
!          On INPUT:  max number of Arnoldi update iterations allowed.
!          On OUTPUT: actual number of Arnoldi update iterations taken.
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" eigenvalues.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used. 
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4,5; See under \Description of dsband for the 
!          five modes available.
!
!          IPARAM(8) = NP
!          Not referenced.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.
!
! WORKD    REAL(DOUBLE) work array of length at least 3*n. (WORKSPACE)
!
! WORKL    REAL(DOUBLE) work array of length LWORKL.  (WORKSPACE)
!
! LWORKL   Integer.  (INPUT)
!          LWORKL must be at least NCV**2 + 8*NCV.
!
! IWORK    Integer array of dimension at least N. (WORKSPACE)
!          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the 
!          factorization of M or (A-SIGMA*M).
!            
! INFO     Integer.  (INPUT/OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)  
!                returns the number of wanted converged Ritz values.
!          =  3: No shifts could be applied during a cycle of the 
!                Implicitly restarted Arnoldi iteration. One possibility 
!                is to increase the size of NCV relative to NEV. 
!                See remark 4 in DSAUPD.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Informational error from LAPACK routine dsteqr.
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: NEV and WHICH = 'BE' are incompatible.
!          = -13: HOWMNY must be one of 'A' or 'P'
!          = -14: DSAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current
!                   Arnoldi factorization.
!
!  info_lapack (output) INTEGER: Error from LAPACK

!  dtbsv_msg (input) CHARACTER
!           = 'Y', have subr DTBSV print Fwd, Back pass messages

!  piters (input) This is debug(47) to tell whether to write eigen
!                 info at each iteration
! \EndDoc
!
!------------------------------------------------------------------------
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
!     May 1995.
!
!\Routines called:
!     dsaupd  ARPACK reverse communication interface routine.
!     dseupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dgbtrf  LAPACK band matrix factorization routine.
!     dgbtrs  LAPACK band linear system solve routine. 
!     dlacpy  LAPACK matrix copy routine.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     ddot    Level 1 BLAS that computes the dot product of two vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     dgbmv   Level 2 BLAS that computes the band matrix vector product.
!
!\Remarks
!  1. The converged Ritz values are always returned in increasing 
!     (algebraic) order.
!
!  2. Currently only HOWMNY = 'A' is implemented. It is included at this
!     stage for the user who wants to incorporate it.
!
!\Author    
!     Danny Sorensen
!     Richard Lehoucq
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: sband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
!
!\EndLib
!
!---------------------------------------------------------------------
!
      subroutine dsband( rvec, howmny, select, d, z, ldz, sigma,               &
                 n, lda, rfac, kl, ku, which, bmat, nev,                       &
                 tol, resid, ncv, v, ldv, iparam, workd, workl,                &
                 lworkl, iwork, info, info_lapack, dtbsv_msg, piters)

      USE SCONTR, ONLY                :  NTERM_KLLDn, NTERM_MLLn,              &
                                         NTERM_KMSMn, NTERM_ALL
      USE SPARSE_MATRICES, ONLY       :  I_KLLDn, J_KLLDn, KLLDn,              &
                                         I_MLLn , J_MLLn , MLLn,               &
                                         I_KMSMn, J_KMSMn, KMSMn
      USE MODEL_STUF, ONLY            :  EIG_LAP_MAT_TYPE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DSBAND'

      INTEGER(LONG)                   :: SUBR_BEGEND = ARPACK_BEGEND
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
! 
      character        which*2, bmat, howmny
      character*1      dtbsv_msg
      character*80     pagetitle
      integer          n, lda, kl, ku, nev, ncv, ldv,                          &
                       ldz, lworkl, info
      integer          info_lapack, dsaupd_loop_count  
      real(double)                                                             &
                       tol, sigma
      logical          rvec
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer          iparam(*), iwork(*)
      logical          select(*)
      real(double)  :: workd1(n), workd2(n)
      real(double)                                                             &
                       d(*), resid(*), v(ldv,*), z(ldz,*),                     &
                       rfac(lda,*), workd(*), workl(*)
! B 05/24/04 //////////////////////////////////////////////////////////B
      integer nev_user
      real(double) eig_old(n)
! E ///////////////////////////////////////////////////////////////////E

! B 06/08/08 //////////////////////////////////////////////////////////B
      integer iter, iter_old
! E ///////////////////////////////////////////////////////////////////E

! B 11/30/09 //////////////////////////////////////////////////////////B
      integer piters
! E ///////////////////////////////////////////////////////////////////E

!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer          ipntr(14)
      REAL(DOUBLE)     KLLDn_DIAG(N), MLLn_DIAG(N), KMSMn_DIAG(N)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer          ido, i, j, type, imid, itop, ibot, ierr
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                        one, zero
      parameter        (one = 1.0, zero = 0.0)
!
!
!     %-----------------------------%
!     | LAPACK & BLAS routines used |
!     %-----------------------------%
!
!:!   REAL(DOUBLE)
!:!  &                 ddot, dnrm2, dlapy2
!:!   external         ddot, dcopy, dgbmv, dgbtrf, 
!:!  &                 dgbtrs, dnrm2, dlapy2, dlacpy
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!     
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF
      ierr = 0
! **********************************************************************************************************************************
! B 05/24/04 //////////////////////////////////////////////////////////B

      do j=1,n
         eig_old(j) = 0.d0
      enddo

      nev_user = nev
      
! E ///////////////////////////////////////////////////////////////////E

!     %----------------------------------------------------------------%
!     | Set type of the problem to be solved.                          |
!     |                                                                |
!     | type = 3 --> Solving generalized problem in regular mode.      |
!     |              NOTE: see expalnation of my mode 2 above          |
!     | type = 4 --> Solving generalized problem in shift-invert mode. |
!     %----------------------------------------------------------------%
!
      if      ( iparam(7) .eq. 1 ) then
         type = 1                                                            ! ARPACK "mode" = 1
      else if ( iparam(7) .eq. 2 ) then
         type = 3                                                            ! ARPACK "mode" = 2
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
         type = 4                                                            ! ARPACK "mode" = 3
      else                                                                   ! Err msg given in BD_EIGRL
         go to 9000
      end if
!
!     %------------------------%
!     | Initialize the reverse |
!     | communication flag.    |
!     %------------------------%
!
      ido   = 0
!
!     %----------------%
!     | Exact shift is |
!     | used.          | 
!     %----------------%
! 
      iparam(1) = 1
!
!     %-----------------------------------%
!     | Both matrices A and M are stored  |
!     | between rows itop and ibot.  Imid |
!     | is the index of the row that      |
!     | stores the diagonal elements.     |
!     %-----------------------------------%
!
      if      (eig_lap_mat_type(1:3) == 'DGB') then
         itop = kl + 1
         ibot = 2*kl + ku + 1
         imid = kl + ku + 1
      else if (eig_lap_mat_type(1:3) == 'DPB') then
         itop = 1
         ibot = ku + 1
         imid = ibot
      endif

! Factorize matrix [K - sigma*M]

!     if (sol_name(1:8) == 'BUCKLING') then
!        modnam1 = 'LAPACK DGB TRIANG FACTOR OF [KAA+sigma*KLLDn]'
!        modnam2 = 'LAPACK DPB TRIANG FACTOR OF [KAA+sigma*KLLDn]'
!     else
!        modnam1 = 'LAPACK DGB TRIANG FACTOR OF [KAA-sigma*MLLn]'
!        modnam2 = 'LAPACK DPB TRIANG FACTOR OF [KAA-sigma*MLLn]'
!     endif
      if ( type .eq. 3 ) then
!
!        %----------------------------------------------%
!        | Solving generalized eigenvalue problem in    |
!        | regular mode. Copy MB = [K - sigma*M] to rfac|
!        | and Call LAPACK routine dpbtrf to factor rfac|
!        %----------------------------------------------%
!
         if      (eig_lap_mat_type(1:3) == 'DGB') then

            call ourtim
!           write(sc1,4092) linkno,modnam1,hour,minute,sec,sfrac
            call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)

         else if (eig_lap_mat_type(1:3) == 'DPB') then

            call ourtim
!           write(sc1,4092) linkno,modnam2,hour,minute,sec,sfrac
            call dpbtrf ( 'U', n, ku, rfac, ku+1, ierr )
            do i=1,n
               iwork(i) = i                                                  ! Pivot indices (no pivoting in DPBTRF)
            enddo

         endif

         if (ierr .ne. 0) then 
            info_lapack = ierr
            go to 9000 
         end if
!
      else if ( type .eq. 4 ) then
!
!        %-------------------------------------------%
!        | Solving generalized eigenvalue problem in |
!        | shift-invert, Buckling, or Cayley mode.   |
!        %-------------------------------------------%
! 
!        %-------------------------------------%
!        | Construct and factor (A - sigma*M). |
!        %-------------------------------------%
!
         if      (eig_lap_mat_type(1:3) == 'DGB') then

            call ourtim
!           write(sc1,4092) linkno,modnam1,hour,minute,sec,sfrac
            call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
!
         else if (eig_lap_mat_type(1:3) == 'DPB') then

            call ourtim
!           write(sc1,4092) linkno,modnam2,hour,minute,sec,sfrac
            call dpbtrf ( 'U', n, ku, rfac, ku+1, ierr )
            do i=1,n
               iwork(i) = i   ! Pivot indices (no pivoting in DPBTRF)
            enddo

         endif

         if ( ierr .ne. 0 )  then
             info_lapack = ierr
             go to 9000
         end if
! 
      end if 

      IF (EIG_MSGLVL > 0) THEN
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            DO I=1,N
               IF (I_KLLDn(I) == I_KLLDn(I+1)) THEN
                  KLLDn_DIAG(I) = ZERO
               ELSE
                  KLLDn_DIAG(I) = KLLDn(I_MLLn(I))
               ENDIF
               IF (I_KMSMn(I) == I_KMSMn(I+1)) THEN
                  KMSMn_DIAG(I) = ZERO
               ELSE
                  KMSMn_DIAG(I) = KMSMn(I_KMSMn(I))
               ENDIF
            ENDDO
         ELSE
            DO I=1,N
               IF (I_MLLn(I) == I_MLLn(I+1)) THEN
                  MLLn_DIAG(I) = ZERO
               ELSE
                  MLLn_DIAG(I) = MLLn(I_MLLn(I))
               ENDIF
               IF (I_KMSMn(I) == I_KMSMn(I+1)) THEN
                  KMSMn_DIAG(I) = ZERO
               ELSE
                  KMSMn_DIAG(I) = KMSMn(I_KMSMn(I))
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!     %--------------------------------------------%
!     |  M A I N   L O O P (reverse communication) |
!     %--------------------------------------------%
!
      iter_old          = 0
      dsaupd_loop_count = 0
      write(sc1, * )
  90  continue 
 
      IF (EIG_MSGLVL > 0) THEN
         WRITE(F06,99990)
         WRITE(F06,98710) dsaupd_loop_count, ido, type,                        &
                          eig_lap_mat_type(1:3)
      ENDIF
!
! B 04/16/10 //////////////////////////////////////////////////////////B
      iter = 0
! E ///////////////////////////////////////////////////////////////////E
      call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv,                 &
                    v, ldv, iparam, ipntr, workd, workl, lworkl,               &
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &              info )
!    &              info, nev_user, sigma, eig_old )
! E ///////////////////////////////////////////////////////////////////E
! B 06/08/08
! //////////////////////////////////////////////////////////B
                    info, nev_user, sigma, eig_old, iter, piters )
! E ///////////////////////////////////////////////////////////////////E
!
      if (iter > iter_old) then
         dsaupd_loop_count = 1
         iter_old = iter
      else
         dsaupd_loop_count = dsaupd_loop_count + 1
      endif
      write(sc1,12345) iter+1, dsaupd_loop_count, ido
      write(f04, 9876) iter+1, dsaupd_loop_count, ido
      write(f06,*) ' In ARPACK_LANCZOS_EIG: type = ', type

! **********************************************************************
      if (ido .eq. -1) then
!
         if ( type .eq. 3 ) then
!
!           %-------------------------------------%
!           | Perform  y <--- OP*x = inv[MB]*AB*x |
!           | MB = [K - sigma*M] and AB = M       |
!           | to force the starting vector into   | 
!           | the range of OP.                    |
!           %-------------------------------------%
!
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn', n, n, nterm_klldn, 'N',             &
                                   i_klldn, j_klldn, klldn,                    &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn', n, n, nterm_mlln, 'N',               &
                                   i_mlln, j_mlln, mlln,                       &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo

            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)

            if      (eig_lap_mat_type(1:3) == 'DGB') then

               call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,            &
                             iwork, workd(ipntr(2)), n, ierr,                  &
                             dtbsv_msg)

            else if (eig_lap_mat_type(1:3) == 'DPB') then

               call dpbtrs ( 'U', n, ku, 1, rfac, ku+1, workd(ipntr(2)),       &
                              n, ierr, 'N' )
            endif
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

            if (ierr .ne. 0) then
               info_lapack = ierr
               go to 9000
            end if
!
         else if ( type .eq. 4 ) then
!
!           %-----------------------------------------%
!           | Perform y <-- OP*x                      |
!           |           = inv[A-SIGMA*M]*M            |  
!           | to force the starting vector into the   |
!           | range of OP.                            |
!           %-----------------------------------------%
!
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',            &
                                   i_klldn, j_klldn, klldn,                    &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',            &
                                   i_mlln, j_mlln, mlln,                       &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)

            if      (eig_lap_mat_type(1:3) == 'DGB') then

               call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,            &
                            iwork, workd(ipntr(2)), n, ierr,                   &
                             dtbsv_msg)

            else if (eig_lap_mat_type(1:3) == 'DPB') then

               call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,                        &
                              workd(ipntr(2)), n, ierr, 'N' )
            endif
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

            if (ierr .ne. 0) then
               info_lapack = ierr
               go to 9000
            end if

         endif
!
! **********************************************************************
      else if (ido .eq. 1) then
!
         if ( type .eq. 3 ) then
!
!           %----------------------------------------------%
!           | Perform  y <--- OP*x = inv[MB]*AB*x          |
!           |                      = inv[K - sigma*M]*M*x  |
!           %----------------------------------------------%
!
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',            &
                                   i_klldn, j_klldn, klldn,                    &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',            &
                                   i_mlln , j_mlln, mlln,                      &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1) 

            if      (eig_lap_mat_type(1:3) == 'DGB') then

               call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,            &
                             iwork, workd(ipntr(2)), n, ierr,                  &
                             dtbsv_msg)

            else if (eig_lap_mat_type(1:3) == 'DPB') then

               call dpbtrs ( 'U', n, ku, 1, rfac, ku+1, workd(ipntr(2)),       &
                              n, ierr, 'N' )
            endif
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

            if (ierr .ne. 0) then
               info_lapack = ierr
               go to 9000
            end if
!
         else if ( type .eq. 4 ) then
!
!           %----------------------------------------%
!           | Perform y <-- inv[AB-sigma*MB]*(MB*x). |
!           |               inv[K - sigma*M]*m*x     |
!           | (MB*x) has been computed and stored    |
!           | in workd(ipntr(3)).                    |           
!           %----------------------------------------%
!
            if      (eig_lap_mat_type(1:3) == 'DGB') then

               call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
               call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,            &
                             iwork, workd(ipntr(2)), n, ierr,                  &
                             dtbsv_msg)

            else if (eig_lap_mat_type(1:3) == 'DPB') then

               call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
               call dpbtrs ( 'U', n, ku, 1, rfac, ku+1, workd(ipntr(2)),       &
                              n, ierr, 'N' )
            endif
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

            if (ierr .ne. 0) then 
               info_lapack = ierr
               go to 9000
            end if
! 
         end if
!
! **********************************************************************
      else if (ido .eq. 2) then
!
!        %----------------------------------%
!        |        Perform y <-- B*x         | 
!        %----------------------------------%
!
         if (type == 3) then

            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo
            call matmult_sff ( 'KMSMn', n, n, nterm_kmsmn, 'N',                &
                                i_kmsmn, j_kmsmn, kmsmn,                       &
                               'workd1', n, 1, workd1, 'N',                    &
                               'workd2', one,  workd2 )

!           CALL DSBMV ( 'U', N, KU, 1.0D0, RFAC, KU+1, WORKD1, 1,  ! This does not work. RFAC is not K-sigma*M anymore
!    $                   0.D0, WORKD2, 1, 1 )                       ! It was overwritten in subr dpbtrf by its triangular factor
            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0)CALL ARP_DEB(3,N,IDO,IPNTR)

         else if (type == 4) then

            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',            &
                                   i_klldn, j_klldn, klldn,                    &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',            &
                                   i_mlln , j_mlln, mlln,                      &
                                  'workd1', n, 1, workd1, 'N',                 &
                                  'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)

         endif

!
      else 
!
!        %--------------------------------------%
!        | Either we have convergence, or error | 
!        %--------------------------------------%
!
         if ( info .lt. 0) then
!
            call arpack_info_msg ('dsaupd',info,iparam,lworkl,nev,ncv)
            go to 9000
!
         else 
!
            if ((info == 1) .or. (info == 3)) then
               call arpack_info_msg('dsaupd',info,iparam,lworkl,nev,ncv)
            endif
!
            if (iparam(5) .gt. 0) then
!
               call dseupd ( rvec, 'A', select, d, z, ldz, sigma,              &
                        bmat, n, which, nev, tol, resid, ncv, v, ldv,          &
                        iparam, ipntr, workd, workl, lworkl, info )            
!
               if ( info .ne. 0) then
! 
              call arpack_info_msg ('dseupd',info,iparam,lworkl,nev,ncv)
                  go to 9000
! 
               end if
!
            end if
!
         end if
!
         go to 9000
!
      end if
!
!     %----------------------------------------%
!     | L O O P  B A C K to call DSAUPD again. |
!     %----------------------------------------%
!
      go to 90 
!
 9000 continue
!
! B 05/24/04 //////////////////////////////////////////////////////////B
      write(f06,*) ! blank line after eigen iteration results
! E ///////////////////////////////////////////////////////////////////E
! **********************************************************************************************************************************
12345 format("+",5X,'Iteration',i4,' Reverse comm loop',i4,                    &
      ' with IDO =',I3)

 9876 format(7X,'Iteration',i4,' Reverse comm loop ',i4,                       &
      ' with IDO =',I3)

 4907 FORMAT(/,22X,A                                                           &
            ,/,7X,'1',12X,'2',12X,'3',12X,'4',12X,'5',12X,'6',12X,             &
              '7',12X,'8',12X,'9',12X,'10')    

 4908 FORMAT(10(1X,1ES12.5))

 4092 FORMAT(1X,I2,'/',A44,18X,2X,I2,':',I2,':',I2,'.',I3)

98710 FORMAT(' dsaupd loop count = ',I4,' ido = ',i4,', "type" = ',I3,         &
             ', using ',a,' LAPACK matrices',/)

99990 FORMAT('**********************************************************',     &
             '***************************')
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

! ##################################################################################################################################
 
      CONTAINS
 
! ##################################################################################################################################

      SUBROUTINE ARP_DEB ( WHICH, N, IDO, IPNTR )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  PROG_NAME, FATAL_ERR, WARN_ERR
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, WRT_LOG, ERR,       &
                                         F04, F06

      IMPLICIT NONE

      INTEGER, INTENT(IN)            :: IDO                ! 
      INTEGER, INTENT(IN)            :: IPNTR(14)          ! 
      INTEGER, INTENT(IN)            :: N                  ! 
      INTEGER, INTENT(IN)            :: WHICH              ! 
      INTEGER                        :: I,J                ! DO loop

! **********************************************************************************************************************************
      IF (WHICH == 1) THEN

         WRITE(F06,98700) IDO, IPNTR(1), IPNTR(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98711)
               DO J=1,N
                  WRITE(F06,98713) J, KLLDn_DIAG(J), WORKD1(J),WORKD2(J)
               ENDDO
            ELSE
               WRITE(F06,98712)
               DO J=1,N
                  WRITE(F06,98713) J, MLLn_DIAG(J), WORKD1(J), WORKD2(J)
               ENDDO
            ENDIF
            WRITE(F06,*)
         ENDIF

      ELSE IF (WHICH == 2) THEN

         WRITE(F06,98700) IDO, IPNTR(1), IPNTR(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98721)
            ELSE
               WRITE(F06,98722)
            ENDIF
            DO J=1,N
               WRITE(F06,98723) J, RFAC(IMID,J), WORKD1(J), WORKD2(J)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ELSE IF (WHICH == 3) THEN

         WRITE(F06,98700) IDO, IPNTR(1), IPNTR(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98731)
            ELSE
               WRITE(F06,98732)
            ENDIF
            DO J=1,N
               WRITE(F06,98733) J, KMSMn_DIAG(J), WORKD1(J), WORKD2(J)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF

! **********************************************************************************************************************************
98700 FORMAT(' ido, ipntr(1), ipntr(2) = ',I9,3x,2i3,/)

98711 FORMAT('                  Diagonal of matrix:',/,                        &
      ' L-set DOF                KLL                  WORKD1(J) WORKD2(J)',/,  &
      ' ---------        -------------------------------   ------------')

98712 FORMAT('                  Diagonal of matrix:',/,                        &
      ' L-set DOF                MLL                  WORKD1(J) WORKD2(J)',/,  &
      ' ---------        -------------------------------   ------------')

98713 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

98721 FORMAT('                  Diagonal of matrix:',/,                        &
      ' L-set DOF  tridiag factor of [KLL-sig*MLL]    WORKD1(J) WORKD2(J)',/,  &
      ' ---------  -------------------------------------------   ------------')

98722 FORMAT('                  Diagonal of matrix:',/,                        &
      ' L-set DOF  tridiag factor of [KLL-sig*KLL]    WORKD1(J) WORKD2(J)',/,  &
      ' ---------  -------------------------------------------   ------------')

98723 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

98731 FORMAT('                  Diagonal of matrix:',/,                        &
      ' L-set DOF          [KLL-sigma*MLL]            WORKD1(J) WORKD2(J)',/,  &
      ' ---------        -------------------------------   ------------')

98732 FORMAT('                  Diagonal of matrix:',/,                        &
      ' L-set DOF          [KLL-sigma*KLL]            WORKD1(J) WORKD2(J)',/,  &
      ' ---------        -------------------------------   ------------')

98733 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

! **********************************************************************************************************************************

      END SUBROUTINE ARP_DEB

      end subroutine dsband

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsaupd 
!
!\Description: 
!
!  Reverse communication interface for the Implicitly Restarted Arnoldi 
!  Iteration.  For symmetric problems this reduces to a variant of the Lanczos 
!  method.  This method has been designed to compute approximations to a 
!  few eigenpairs of a linear operator OP that is real and symmetric 
!  with respect to a real positive semi-definite symmetric matrix B, 
!  i.e.
!                   
!       B*OP = (OP`)*B.  
!
!  Another way to express this condition is 
!
!       < x,OPy > = < OPx,y >  where < z,w > = z`Bw  .
!  
!  In the standard eigenproblem B is the identity matrix.  
!  ( A` denotes transpose of A)
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  dsaupd  is usually called iteratively to solve one of the 
!  following problems:
!
!  Mode 1:  A*x = lambda*x, A symmetric 
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)
!
!  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
!           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
!           ===> Shift-and-Invert mode
!
!  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite, 
!           KG symmetric indefinite
!           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
!           ===> Buckling mode
!
!  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
!           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
!           ===> Cayley transformed mode
!
!  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
!        should be accomplished either by a direct method
!        using a sparse matrix factorization and solving
!
!           [A - sigma*M]*w = v  or M*w = v,
!
!        or through an iterative method for solving these
!        systems.  If an iterative method is used, the
!        convergence test must be more stringent than
!        the accuracy requirements for the eigenvalue
!        approximations.
!
!\Usage:
!  call dsaupd  
!     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO, PITERS )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first 
!          call to dsaupd .  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          dsaupd  with the result.  The operand is given in
!          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
!          (If Mode = 2 see remark 5 below)
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    In mode 3,4 and 5, the vector B * X is already
!                    available in WORKD(ipntr(3)).  It does not
!                    need to be recomputed in forming OP * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO =  3: compute the IPARAM(8) shifts where
!                    IPNTR(11) is the pointer into WORKL for
!                    placing the shifts. See remark 6 below.
!          IDO = 99: done
!          -------------------------------------------------------------
!             
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  WHICH   Character*2.  (INPUT)
!          Specify which of the Ritz values of OP to compute.
!
!          'LA' - compute the NEV largest (algebraic) eigenvalues.
!          'SA' - compute the NEV smallest (algebraic) eigenvalues.
!          'LM' - compute the NEV largest (in magnitude) eigenvalues.
!          'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
!          'BE' - compute NEV eigenvalues, half from each end of the
!                 spectrum.  When NEV is odd, compute one more from the
!                 high end than from the low end.
!           (see remark 1 below)
!
!  NEV     Integer.  (INPUT)
!          Number of eigenvalues of OP to be computed. 0 < NEV < N.
!
!  TOL     REAL(DOUBLE)  scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value 
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
!          If TOL .LE. 0. is passed a default is set:
!          DEFAULT = DLAMCH ('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine DLAMCH ).
!
!  RESID   REAL(DOUBLE)  array of length N.  (INPUT/OUTPUT)
!          On INPUT: 
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector. 
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V (less than or equal to N).
!          This will indicate how many Lanczos vectors are generated 
!          at each iteration.  After the startup phase in which NEV 
!          Lanczos vectors are generated, the algorithm generates 
!          NCV-NEV Lanczos vectors at each subsequent update iteration.
!          Most of the cost in generating each Lanczos vector is in the 
!          matrix-vector product OP*x. (See remark 4 below).
!
!  V       REAL(DOUBLE)  N by NCV array.  (OUTPUT)
!          The NCV columns of V contain the Lanczos basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are provided by the user via
!                      reverse communication.  The NCV eigenvalues of
!                      the current tridiagonal matrix T are returned in
!                      the part of WORKL array corresponding to RITZ.
!                      See remark 6 below.
!          ISHIFT = 1: exact shifts with respect to the reduced 
!                      tridiagonal matrix T.  This is equivalent to 
!                      restarting the iteration with a starting vector 
!                      that is a linear combination of Ritz vectors 
!                      associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = LEVEC
!          No longer referenced. See remark 2 below.
!
!          IPARAM(3) = MXITER
!          On INPUT:  maximum number of Arnoldi update iterations allowed. 
!          On OUTPUT: actual number of Arnoldi update iterations taken. 
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" Ritz values.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used. 
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4,5; See under \Description of dsaupd  for the 
!          five modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), dsaupd  returns NP, the number
!          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
!          6 below.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.        
!
!  IPNTR   Integer array of length 11.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD and WORKL
!          arrays for matrices/vectors used by the Lanczos iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X in WORKD.
!          IPNTR(2): pointer to the current result vector Y in WORKD.
!          IPNTR(3): pointer to the vector B * X in WORKD when used in 
!                    the shift-and-invert mode.
!          IPNTR(4): pointer to the next available location in WORKL
!                    that is untouched by the program.
!          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
!          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
!          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
!                    with the Ritz values located in RITZ in WORKL.
!          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
!
!          Note: IPNTR(8:10) is only referenced by dseupd . See Remark 2.
!          IPNTR(8): pointer to the NCV RITZ values of the original system.
!          IPNTR(9): pointer to the NCV corresponding error bounds.
!          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
!                     of the tridiagonal matrix T. Only referenced by
!                     dseupd  if RVEC = .TRUE. See Remarks.
!          -------------------------------------------------------------
!          
!  WORKD   REAL(DOUBLE)  work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD 
!          as temporary workspace during the iteration. Upon termination
!          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
!          subroutine dseupd  uses this output.
!          See Data Distribution Note below.  
!
!  WORKL   REAL(DOUBLE)  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least NCV**2 + 8*NCV .
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)  
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the 
!                Implicitly restarted Arnoldi iteration. One possibility 
!                is to increase the size of NCV relative to NEV. 
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV must be greater than NEV and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iterations allowed
!                must be greater than zero.
!          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array WORKL is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Informatinal error from LAPACK routine dsteqr .
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -13: NEV and WHICH = 'BE' are incompatable.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization. The user is advised to check that
!                   enough workspace and array storage has been allocated.
!
!
!\Remarks
!  1. The converged Ritz values are always returned in ascending 
!     algebraic order.  The computed Ritz values are approximate
!     eigenvalues of OP.  The selection of WHICH should be made
!     with this in mind when Mode = 3,4,5.  After convergence, 
!     approximate eigenvalues of the original problem may be obtained 
!     with the ARPACK subroutine dseupd . 
!
!  2. If the Ritz vectors corresponding to the converged Ritz values
!     are needed, the user must call dseupd  immediately following completion
!     of dsaupd . This is new starting with version 2.1 of ARPACK.
!
!  3. If M can be factored into a Cholesky factorization M = LL`
!     then Mode = 2 should not be selected.  Instead one should use
!     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular 
!     linear systems should be solved with L and L` rather
!     than computing inverses.  After convergence, an approximate
!     eigenvector z of the original problem is recovered by solving
!     L`z = x  where x is a Ritz vector of OP.
!
!  4. At present there is no a-priori analysis to guide the selection
!     of NCV relative to NEV.  The only formal requrement is that NCV > NEV.
!     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
!     the same type are to be solved, one should experiment with increasing
!     NCV while keeping NEV fixed for a given test problem.  This will 
!     usually decrease the required number of OP*x operations but it
!     also increases the work and storage required to maintain the orthogonal
!     basis vectors.   The optimal "cross-over" with respect to CPU time
!     is problem dependent and must be determined empirically.
!
!  5. If IPARAM(7) = 2 then in the Reverse commuication interface the user
!     must do the following. When IDO = 1, Y = OP * X is to be computed.
!     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
!     must overwrite X with A*X. Y is then the solution to the linear set
!     of equations B*Y = A*X.
!
!  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
!     NP = IPARAM(8) shifts in locations: 
!     1   WORKL(IPNTR(11))           
!     2   WORKL(IPNTR(11)+1)         
!                        .           
!                        .           
!                        .      
!     NP  WORKL(IPNTR(11)+NP-1). 
!
!     The eigenvalues of the current tridiagonal matrix are located in 
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
!     order defined by WHICH. The associated Ritz estimates are located in
!     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------
!
!\Data Distribution Note:
!
!  Fortran-D syntax:
!  ================
!  REAL       RESID(N), V(LDV,NCV), WORKD(3*N), WORKL(LWORKL)
!  DECOMPOSE  D1(N), D2(N,NCV)
!  ALIGN      RESID(I) with D1(I)
!  ALIGN      V(I,J)   with D2(I,J)
!  ALIGN      WORKD(I) with D1(I)     range (1:N)
!  ALIGN      WORKD(I) with D1(I-N)   range (N+1:2*N)
!  ALIGN      WORKD(I) with D1(I-2*N) range (2*N+1:3*N)
!  DISTRIBUTE D1(BLOCK), D2(BLOCK,:)
!  REPLICATED WORKL(LWORKL)
!
!  Cray MPP syntax:
!  ===============
!  REAL       RESID(N), V(LDV,NCV), WORKD(N,3), WORKL(LWORKL)
!  SHARED     RESID(BLOCK), V(BLOCK,:), WORKD(BLOCK,:)
!  REPLICATED WORKL(LWORKL)
!  
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!  8. R.B. Lehoucq, D.C. Sorensen, "Implementation of Some Spectral
!     Transformations in a k-Step Arnoldi Method". In Preparation.
!
!\Routines called:
!     dsaup2   ARPACK routine that implements the Implicitly Restarted
!             Arnoldi Iteration.
!     dstats   ARPACK routine that initialize timing and other statistics
!             variables.
!     ivout   ARPACK utility routine that prints integers.
!     second  ARPACK utility routine for timing.
!     dvout    ARPACK utility routine that prints vectors.
!     dlamch   LAPACK routine that determines machine constants.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     12/15/93: Version ' 2.4' 
!
!\SCCS Information: @(#) 
! FILE: saupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2 
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsaupd                                                        &
         ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,          &
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &     ipntr, workd, workl, lworkl, info )
!    &     ipntr, workd, workl, lworkl, info, nev_user, sigma, eig_old )
! E ///////////////////////////////////////////////////////////////////E
! B 06/08/08
! //////////////////////////////////////////////////////////B
           ipntr, workd, workl, lworkl, info, nev_user, sigma, eig_old,        &
           iter, piters )
! E ///////////////////////////////////////////////////////////////////E


! B 06/08/08 //////////////////////////////////////////////////////////B
      integer iter     
! E ///////////////////////////////////////////////////////////////////E

! B 05/24/04 //////////////////////////////////////////////////////////B
      integer n, nev_user
      real(double) sigma, eig_old(n)
! E ///////////////////////////////////////////////////////////////////E

! B 11/30/09 //////////////////////////////////////////////////////////B
      integer piters
! E ///////////////////////////////////////////////////////////////////E

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAUPD'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, ncv, nev
      REAL(DOUBLE)                                                             &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(11), ipntr(11)
      REAL(DOUBLE)                                                             &
                 resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    bounds, ierr, ih, iq, ishift, iupd, iw,                       &
                 ldh, ldq, msglvl, mxiter, mode, nb,                           &
                 nev0, next, np, ritz, j
      INTEGER                                                                  &
                 MXITER_ARRAY(1), NP_ARRAY(1)
      save       bounds, ierr, ih, iq, ishift, iupd, iw,                       &
                 ldh, ldq, msglvl, mxiter, mode, nb,                           &
                 nev0, next, np, ritz
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!     external   dsaup2 ,  dvout , ivout, second, dstats 
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!:!   REAL(DOUBLE) 
!:!  &           dlamch 
!:!   external   dlamch 
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
! 
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+1) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      if (ido .eq. 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call dstats 
         call cpu_time (t0)
         msglvl = msaupd
!
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
!         nb     = iparam(4)
         nb     = 1
!
!        %--------------------------------------------%
!        | Revision 2 performs only implicit restart. |
!        %--------------------------------------------%
!
         iupd   = 1
         mode   = iparam(7)
!
!        %----------------%
!        | Error checking |
!        %----------------%
!
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev .or.  ncv .gt. n) then
            ierr = -3
         end if
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        %----------------------------------------------%
!
         np     = ncv - nev
! 
         if (mxiter .le. 0)                     ierr = -4
         if (which .ne. 'LM' .and.                                             &
             which .ne. 'SM' .and.                                             &
             which .ne. 'LA' .and.                                             &
             which .ne. 'SA' .and.                                             &
             which .ne. 'BE')                   ierr = -5
         if (bmat .ne. 'I' .and. bmat .ne. 'G') ierr = -6
!
         if (lworkl .lt. ncv**2 + 8*ncv)        ierr = -7
         if (mode .lt. 1 .or. mode .gt. 5) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         else if (nev .eq. 1 .and. which .eq. 'BE') then
                                                ierr = -13
         end if
! 
!        %------------%
!        | Error Exit |
!        %------------%
!
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
! 
!        %------------------------%
!        | Set default parameters |
!        %------------------------%
!
         if (nb .le. 0)                         nb = 1
         if (tol .le. zero)                     tol = dlamch ('EpsMach')
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        | NEV0 is the local variable designating the   |
!        | size of the invariant subspace desired.      |
!        %----------------------------------------------%
!
         np     = ncv - nev
         nev0   = nev 
! 
!        %-----------------------------%
!        | Zero out internal workspace |
!        %-----------------------------%
!
         do 10 j = 1, ncv**2 + 8*ncv
            workl(j) = zero
 10      continue
! 
!        %-------------------------------------------------------%
!        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
!        | etc... and the remaining workspace.                   |
!        | Also update pointer to be used on output.             |
!        | Memory is laid out as follows:                        |
!        | workl(1:2*ncv) := generated tridiagonal matrix        |
!        | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
!        | workl(3*ncv+1:3*ncv+ncv) := computed error bounds     |
!        | workl(4*ncv+1:4*ncv+ncv*ncv) := rotation matrix Q     |
!        | workl(4*ncv+ncv*ncv+1:7*ncv+ncv*ncv) := workspace     |
!        %-------------------------------------------------------%
!
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + 2*ldh
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ncv**2
         next   = iw     + 3*ncv
!
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritz
         ipntr(7) = bounds
         ipntr(11) = iw
      end if
!
!     %-------------------------------------------------------%
!     | Carry out the Implicitly restarted Lanczos Iteration. |
!     %-------------------------------------------------------%
!
      call dsaup2                                                              &
         ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,              &
           ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz),                &
           workl(bounds), workl(iq), ldq, workl(iw), ipntr, workd,             &
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &     info )
!    &     info, nev_user, sigma, eig_old )
! B 06/08/08
! //////////////////////////////////////////////////////////B
           info, nev_user, sigma, eig_old, iter, piters )
     
! E ///////////////////////////////////////////////////////////////////E
!
!     %--------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication |
!     | to compute operations involving OP or shifts.    |
!     %--------------------------------------------------%
!
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
! 
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
!
!     %------------------------------------%
!     | Exit if there was an informational |
!     | error within dsaup2 .               |
!     %------------------------------------%
!
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
!
      if (msglvl .gt. 0) then
         MXITER_ARRAY(1) = MXITER
         NP_ARRAY(1)     = NP
         call ivout (logfil, 1, MXITER_ARRAY, ndigit,                          &
                     '_saupd: number of update iterations taken')
         call ivout (logfil, 1, NP_ARRAY, ndigit,                              &
                     '_saupd: number of "converged" Ritz values')
         call dvout  (logfil, np, workl(Ritz), ndigit,                         &
                     '_saupd: final Ritz values')
         call dvout  (logfil, np, workl(Bounds), ndigit,                       &
                     '_saupd: corresponding error bounds')
      end if 
!
      call cpu_time (t1)
      tsaupd = t1 - t0
! 
      if (msglvl .gt. 0) then
!
!        %--------------------------------------------------------%
!        | Version Number & Version Date are defined in version.h |
!        %--------------------------------------------------------%
!
         write (logfil,1000)
         write (logfil,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,        &
                        tmvopx, tmvbx, tsaupd, tsaup2, tsaitr, titref,         &
                        tgetv0, tseigt, tsgets, tsapps, tsconv
 1000    format (//,                                                           &
            5x, '==========================================',/                 &
            5x, '= Symmetric implicit Arnoldi update code =',/                 &
            5x, '= Version Number:', ' 2.4' , 19x, ' =',/                      &
            5x, '= Version Date:  ', ' 07/31/96' , 14x, ' =',/                 &
            5x, '==========================================',/                 &
            5x, '= Summary of timing statistics           =',/                 &
            5x, '==========================================',//)
 1100    format (                                                              &
            5x, 'Total number update iterations             = ', i5,/          &
            5x, 'Total number of OP*x operations            = ', i5,/          &
            5x, 'Total number of B*x operations             = ', i5,/          &
            5x, 'Total number of reorthogonalization steps  = ', i5,/          &
            5x, 'Total number of iterative refinement steps = ', i5,/          &
            5x, 'Total number of restart steps              = ', i5,/          &
            5x, 'Total time in user OP*x operation          = ', f12.6,/       &
            5x, 'Total time in user B*x operation           = ', f12.6,/       &
            5x, 'Total time in Arnoldi update routine       = ', f12.6,/       &
            5x, 'Total time in saup2 routine                = ', f12.6,/       &
            5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/       &
            5x, 'Total time in reorthogonalization phase    = ', f12.6,/       &
            5x, 'Total time in (re)start vector generation  = ', f12.6,/       &
            5x, 'Total time in trid eigenvalue subproblem   = ', f12.6,/       &
            5x, 'Total time in getting the shifts           = ', f12.6,/       &
            5x, 'Total time in applying the shifts          = ', f12.6,/       &
            5x, 'Total time in convergence testing          = ', f12.6)
      end if
! 
 9000 continue
! 
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+1) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!
!     %---------------%
!     | End of dsaupd |
!     %---------------%
!
      end subroutine dsaupd

! ################################################################################################################################## 

!\BeginDoc
!
!\Name: dseupd 
!
!\Description: 
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) the corresponding approximate eigenvectors,
!
!      (2) an orthonormal (Lanczos) basis for the associated approximate
!          invariant subspace,
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  (Lanczos) basis is always computed.  There is an additional storage cost 
!  of n*nev if both are requested (in this case a separate array Z must be 
!  supplied).
!
!  These quantities are obtained from the Lanczos factorization computed
!  by DSAUPD  for the linear operator OP prescribed by the MODE selection
!  (see IPARAM(7) in DSAUPD  documentation.)  DSAUPD  must be called before
!  this routine is called. These approximate eigenvalues and vectors are 
!  commonly called Ritz values and Ritz vectors respectively.  They are 
!  referred to as such in the comments that follow.   The computed orthonormal 
!  basis for the invariant subspace corresponding to these Ritz values is 
!  referred to as a Lanczos basis.
!
!  See documentation in the header of the subroutine DSAUPD  for a definition 
!  of OP as well as other terms and the relation of computed Ritz values 
!  and vectors of OP with respect to the given problem  A*z = lambda*B*z.  
!
!  The approximate eigenvalues of the original problem are returned in
!  ascending algebraic order.  The user may elect to call this routine
!  once for each desired Ritz vector and store it peripherally if desired.
!  There is also the option of computing a selected set of these vectors
!  with a single call.
!
!\Usage:
!  call dseupd  
!     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
!       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!  RVEC    LOGICAL  (INPUT) 
!          Specifies whether Ritz vectors corresponding to the Ritz value 
!          approximations to the eigenproblem A*z = lambda*B*z are computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute Ritz vectors.
!
!  HOWMNY  Character*1  (INPUT) 
!          Specifies how many Ritz vectors are wanted and the form of Z
!          the matrix of Ritz vectors. See remark 1 below.
!          = 'A': compute NEV Ritz vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
!          If HOWMNY = 'A' , SELECT is used as a workspace for
!          reordering the Ritz values.
!
!  D       REAL(DOUBLE)  array of dimension NEV.  (OUTPUT)
!          On exit, D contains the Ritz value approximations to the
!          eigenvalues of A*z = lambda*B*z. The values are returned
!          in ascending order. If IPARAM(7) = 3,4,5 then D represents
!          the Ritz values of OP computed by dsaupd  transformed to
!          those of the original eigensystem A*z = lambda*B*z. If 
!          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
!          as the those of A*z = lambda*B*z.
!
!  Z       REAL(DOUBLE)  N by NEV array if HOWMNY = 'A'.  (OUTPUT)
!          On exit, Z contains the B-orthonormal Ritz vectors of the
!          eigensystem A*z = lambda*B*z corresponding to the Ritz
!          value approximations.
!          If  RVEC = .FALSE. then Z is not referenced.
!          NOTE: The array Z may be set equal to first NEV columns of the 
!          Arnoldi/Lanczos basis array V computed by DSAUPD .
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
!
!  SIGMA   REAL(DOUBLE)   (INPUT)
!          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
!          IPARAM(7) = 1 or 2.
!
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to DSAUPD  that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to DSEUPD  following the last call
!         to DSAUPD .  These arguments MUST NOT BE MODIFIED between
!         the the last call to DSAUPD  and the call to DSEUPD .
!
!  Two of these parameters (WORKL, INFO) are also output parameters:
!
!  WORKL   REAL(DOUBLE)  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:4*ncv) contains information obtained in
!          dsaupd .  They are not changed by dseupd .
!          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
!          untransformed Ritz values, the computed error estimates,
!          and the associated eigenvector matrix of H.
!
!          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
!          of the above information computed by dseupd .
!          -------------------------------------------------------------
!          IPNTR(8): pointer to the NCV RITZ values of the original system.
!          IPNTR(9): pointer to the NCV corresponding error bounds.
!          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
!                     of the tridiagonal matrix T. Only referenced by
!                     dseupd  if RVEC = .TRUE. See Remarks.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV must be greater than NEV and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Information error from LAPACK routine dsteqr .
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: NEV and WHICH = 'BE' are incompatible.
!          = -14: DSAUPD  did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
!          = -16: HOWMNY = 'S' not yet implemented
!          = -17: DSEUPD  got a different count of the number of converged
!                 Ritz values than DSAUPD  got.  This indicates the user
!                 probably made an error in passing data from DSAUPD  to
!                 DSEUPD  or that the data was modified before entering 
!                 DSEUPD .
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Remarks
!  1. The converged Ritz values are always returned in increasing 
!     (algebraic) order.
!
!  2. Currently only HOWMNY = 'A' is implemented. It is included at this
!     stage for the user who wants to incorporate it. 
!
!\Routines called:
!     dsesrt   ARPACK routine that sorts an array X, and applies the
!             corresponding permutation to a matrix A.
!     dsortr   dsortr   ARPACK sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     dvout    ARPACK utility routine that prints vectors.
!     dgeqr2   LAPACK routine that computes the QR factorization of
!             a matrix.
!     dlacpy   LAPACK matrix copy routine.
!     dlamch   LAPACK routine that determines machine constants.
!     dorm2r   LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     dsteqr   LAPACK routine that computes eigenvalues and eigenvectors
!             of a tridiagonal matrix.
!     dger     Level 2 BLAS rank one update to a matrix.
!     dcopy    Level 1 BLAS that copies one vector to another .
!     dnrm2    Level 1 BLAS that computes the norm of a vector.
!     dscal    Level 1 BLAS that scales a vector.
!     dswap    Level 1 BLAS that swaps the contents of two vectors.

!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas
!     Dept. of Computational & 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     12/15/93: Version ' 2.1'
!
!\SCCS Information: @(#) 
! FILE: seupd.F   SID: 2.11   DATE OF SID: 04/10/01   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
      subroutine dseupd (rvec  , howmny, select, d    ,                        &
                         z     , ldz   , sigma , bmat ,                        &
                         n     , which , nev   , tol  ,                        &
                         resid , ncv   , v     , ldv  ,                        &
                         iparam, ipntr , workd , workl,                        &
                         lworkl, info )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSEUPD'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      REAL(DOUBLE)                                                             &
                 sigma, tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      REAL(DOUBLE)                                                             &
                 d(nev)     , resid(n)  , v(ldv,ncv),                          &
                 z(ldz, nev), workd(2*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  type*6
      integer    bounds , ierr   , ih    , ihb   , ihd   ,                     &
                 iq     , iw     , j     , k     , ldh   ,                     &
                 ldq    , mode   , msglvl, nconv , next  ,                     &
                 ritz   , irz    , ibd   , np    , ishift,                     &
                 leftptr, rghtptr, numcnv, jj
      INTEGER                                                                  &
                 NUMCNV_ARRAY(1), NCONV_ARRAY(1)
      REAL(DOUBLE)                                                             &
                 bnorm2 , rnorm, temp, temp1, eps23
      REAL(DOUBLE) TEMP_ARRAY(1)
      logical    reord
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   dcopy  , dger   , dgeqr2 , dlacpy , dorm2r , dscal , 
!:!  &           dsesrt , dsteqr , dswap  , dvout  , ivout , dsortr 
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!:!   REAL(DOUBLE) 
!:!  &           dnrm2 , dlamch 
!:!   external   dnrm2 , dlamch 
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+1) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      TEMP = 0.0D0
! 
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
      msglvl = mseupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
!
!     %--------------%
!     | Quick return |
!     %--------------%
!
      if (nconv .eq. 0) go to 9000
      ierr = 0
!
      if (nconv .le. 0)                        ierr = -14 
      if (n .le. 0)                            ierr = -1
      if (nev .le. 0)                          ierr = -2
      if (ncv .le. nev .or.  ncv .gt. n)       ierr = -3
      if (which .ne. 'LM' .and.                                                &
          which .ne. 'SM' .and.                                                &
          which .ne. 'LA' .and.                                                &
          which .ne. 'SA' .and.                                                &
          which .ne. 'BE')                     ierr = -5
      if (bmat .ne. 'I' .and. bmat .ne. 'G')   ierr = -6
      if ( (howmny .ne. 'A' .and.                                              &
                 howmny .ne. 'P' .and.                                         &
                 howmny .ne. 'S') .and. rvec )                                 &
                                               ierr = -15
      if (rvec .and. howmny .eq. 'S')           ierr = -16
!
      if (rvec .and. lworkl .lt. ncv**2+8*ncv) ierr = -7
!     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else if (mode .eq. 4 ) then
         type = 'BUCKLE'
      else if (mode .eq. 5 ) then
         type = 'CAYLEY'
      else 
                                               ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')     ierr = -11
      if (nev .eq. 1 .and. which .eq. 'BE')    ierr = -12
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
!     
!     %-------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
!     | etc... and the remaining workspace.                   |
!     | Also update pointer to be used on output.             |
!     | Memory is laid out as follows:                        |
!     | workl(1:2*ncv) := generated tridiagonal matrix H      |
!     |       The subdiagonal is stored in workl(2:ncv).      |
!     |       The dead spot is workl(1) but upon exiting      |
!     |       dsaupd  stores the B-norm of the last residual   |
!     | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
!     |       The wanted values are in the first NCONV spots. |
!     | workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   |
!     |       The wanted values are in the first NCONV spots. |
!     | NOTE: workl(1:4*ncv) is set by dsaupd  and is not      |
!     |       modified by dseupd .                             |
!     %-------------------------------------------------------%
!
!     %-------------------------------------------------------%
!     | The following is used and set by dseupd .              |
!     | workl(4*ncv+1:4*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the diagonal of H. Upon EXIT contains the NCV   |
!     |       Ritz values of the original system. The first   |
!     |       NCONV spots have the wanted values. If MODE =   |
!     |       1 or 2 then will equal workl(2*ncv+1:3*ncv).    |
!     | workl(5*ncv+1:5*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the subdiagonal of H. Upon EXIT contains the    |
!     |       NCV corresponding Ritz estimates of the         |
!     |       original system. The first NCONV spots have the |
!     |       wanted values. If MODE = 1,2 then will equal    |
!     |       workl(3*ncv+1:4*ncv).                           |
!     | workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  |
!     |       the eigenvector matrix for H as returned by     |
!     |       dsteqr . Not referenced if RVEC = .False.        |
!     |       Ordering follows that of workl(4*ncv+1:5*ncv)   |
!     | workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         |
!     |       Workspace. Needed by dsteqr  and by dseupd .      |
!     | GRAND total of NCV*(NCV+8) locations.                 |
!     %-------------------------------------------------------%
!
!
      ih     = ipntr(5)
      ritz   = ipntr(6)
      bounds = ipntr(7)
      ldh    = ncv
      ldq    = ncv
      ihd    = bounds + ldh
      ihb    = ihd    + ldh
      iq     = ihb    + ldh
      iw     = iq     + ldh*ncv
      next   = iw     + 2*ncv
      ipntr(4)  = next
      ipntr(8)  = ihd
      ipntr(9)  = ihb
      ipntr(10) = iq
!
!     %----------------------------------------%
!     | irz points to the Ritz values computed |
!     |     by _seigt before exiting _saup2.   |
!     | ibd points to the Ritz estimates       |
!     |     computed by _seigt before exiting  |
!     |     _saup2.                            |
!     %----------------------------------------%
!
      irz = ipntr(11)+ncv
      ibd = irz+ncv
!
!
!     %---------------------------------%
!     | Set machine dependent constant. |
!     %---------------------------------%
!
      eps23 = dlamch ('Epsilon-Machine') 
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
!
!     %---------------------------------------%
!     | RNORM is B-norm of the RESID(1:N).    |
!     | BNORM2 is the 2 norm of B*RESID(1:N). |
!     | Upon exit of dsaupd  WORKD(1:N) has    |
!     | B*RESID(1:N).                         |
!     %---------------------------------------%
!
      rnorm = workl(ih)
      if (bmat .eq. 'I') then
         bnorm2 = rnorm
      else if (bmat .eq. 'G') then
         bnorm2 = dnrm2 (n, workd, 1)
      end if
!
      if (msglvl .gt. 2) then
         call dvout (logfil, ncv, workl(irz), ndigit,                          &
         '_seupd: Ritz values passed in from _SAUPD.')
         call dvout (logfil, ncv, workl(ibd), ndigit,                          &
         '_seupd: Ritz estimates passed in from _SAUPD.')
      end if
!
      if (rvec) then
!
         reord = .false.
!
!        %---------------------------------------------------%
!        | Use the temporary bounds array to store indices   |
!        | These will be used to mark the select array later |
!        %---------------------------------------------------%
!
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
!
!        %-------------------------------------%
!        | Select the wanted Ritz values.      |
!        | Sort the Ritz values so that the    |
!        | wanted ones appear at the tailing   |
!        | NEV positions of workl(irr) and     |
!        | workl(iri).  Move the corresponding |
!        | error estimates in workl(bound)     |
!        | accordingly.                        |
!        %-------------------------------------%
!
         np     = ncv - nev
         ishift = 0
         call dsgets (ishift, which       , nev          ,                     &
                      np    , workl(irz)  , workl(bounds),                     &
                      workl)
!
         if (msglvl .gt. 2) then
            call dvout (logfil, ncv, workl(irz), ndigit,                       &
            '_seupd: Ritz values after calling _SGETS.')
            call dvout (logfil, ncv, workl(bounds), ndigit,                    &
            '_seupd: Ritz value indices after calling _SGETS.')
         end if
!
!        %-----------------------------------------------------%
!        | Record indices of the converged wanted Ritz values  |
!        | Mark the select array for possible reordering       |
!        %-----------------------------------------------------%
!
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23, abs(workl(irz+ncv-j)) )
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.                                        &
                workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
!
!        %-----------------------------------------------------------%
!        | Check the count (numcnv) of converged Ritz values with    |
!        | the number (nconv) reported by _saupd.  If these two      |
!        | are different then there has probably been an error       |
!        | caused by incorrect passing of the _saupd data.           |
!        %-----------------------------------------------------------%
!
         NUMCNV_ARRAY(1) = NUMCNV
         NCONV_ARRAY(1)  = NCONV
         if (msglvl .gt. 2) then
             call ivout(logfil, 1, NUMCNV_ARRAY, ndigit,                       &
                  '_seupd: Number of specified eigenvalues')
             call ivout(logfil, 1, NCONV_ARRAY, ndigit,                        &
                  '_seupd: Number of "converged" eigenvalues')
         end if
!
         if (numcnv .ne. nconv) then
            info = -17
            go to 9000
         end if
!
!        %-----------------------------------------------------------%
!        | Call LAPACK routine _steqr to compute the eigenvalues and |
!        | eigenvectors of the final symmetric tridiagonal matrix H. |
!        | Initialize the eigenvector matrix Q to the identity.      |
!        %-----------------------------------------------------------%
!
         call dcopy (ncv-1, workl(ih+1), 1, workl(ihb), 1)
         call dcopy (ncv, workl(ih+ldh), 1, workl(ihd), 1)
!
         call dsteqr ('Identity', ncv, workl(ihd), workl(ihb),                 &
                      workl(iq) , ldq, workl(iw), ierr)
!
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
!
         if (msglvl .gt. 1) then
            call dcopy (ncv, workl(iq+ncv-1), ldq, workl(iw), 1)
            call dvout (logfil, ncv, workl(ihd), ndigit,                       &
                '_seupd: NCV Ritz values of the final H matrix')
            call dvout (logfil, ncv, workl(iw), ndigit,                        &
                 '_seupd: last row of the eigenvector matrix for H')
         end if
!
         if (reord) then
!
!           %---------------------------------------------%
!           | Reordered the eigenvalues and eigenvectors  |
!           | computed by _steqr so that the "converged"  |
!           | eigenvalues appear in the first NCONV       |
!           | positions of workl(ihd), and the associated |
!           | eigenvectors appear in the first NCONV      |
!           | columns.                                    |
!           %---------------------------------------------%
!
            leftptr = 1
            rghtptr = ncv
!
            if (ncv .eq. 1) go to 30
!
 20         if (select(leftptr)) then
!
!              %-------------------------------------------%
!              | Search, from the left, for the first Ritz |
!              | value that has not converged.             |
!              %-------------------------------------------%
!
               leftptr = leftptr + 1
!
            else if ( .not. select(rghtptr)) then
!
!              %----------------------------------------------%
!              | Search, from the right, the first Ritz value |
!              | that has converged.                          |
!              %----------------------------------------------%
!
               rghtptr = rghtptr - 1
!
            else
!
!              %----------------------------------------------%
!              | Swap the Ritz value on the left that has not |
!              | converged with the Ritz value on the right   |
!              | that has converged.  Swap the associated     |
!              | eigenvector of the tridiagonal matrix H as   |
!              | well.                                        |
!              %----------------------------------------------%
!
               temp = workl(ihd+leftptr-1)
               workl(ihd+leftptr-1) = workl(ihd+rghtptr-1)
               workl(ihd+rghtptr-1) = temp
               call dcopy (ncv, workl(iq+ncv*(leftptr-1)), 1,                  &
                          workl(iw), 1)
               call dcopy (ncv, workl(iq+ncv*(rghtptr-1)), 1,                  &
                          workl(iq+ncv*(leftptr-1)), 1)
               call dcopy (ncv, workl(iw), 1,                                  &
                          workl(iq+ncv*(rghtptr-1)), 1)
               leftptr = leftptr + 1
               rghtptr = rghtptr - 1
!
            end if
!
            if (leftptr .lt. rghtptr) go to 20
!
 30      end if
!
         if (msglvl .gt. 2) then
             call dvout  (logfil, ncv, workl(ihd), ndigit,                     &
             '_seupd: The eigenvalues of H--reordered')
         end if
!
!        %----------------------------------------%
!        | Load the converged Ritz values into D. |
!        %----------------------------------------%
!
         call dcopy (nconv, workl(ihd), 1, d, 1)
!
      else
!
!        %-----------------------------------------------------%
!        | Ritz vectors not required. Load Ritz values into D. |
!        %-----------------------------------------------------%
!
         call dcopy (nconv, workl(ritz), 1, d, 1)
         call dcopy (ncv, workl(ritz), 1, workl(ihd), 1)
!
      end if
!
!     %------------------------------------------------------------------%
!     | Transform the Ritz values and possibly vectors and corresponding |
!     | Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values |
!     | (and corresponding data) are returned in ascending order.        |
!     %------------------------------------------------------------------%
!
      if (type .eq. 'REGULR') then
!
!        %---------------------------------------------------------%
!        | Ascending sort of wanted Ritz values, vectors and error |
!        | bounds. Not necessary if only Ritz values are desired.  |
!        %---------------------------------------------------------%
!
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
         end if
!
      else 
! 
!        %-------------------------------------------------------------%
!        | *  Make a copy of all the Ritz values.                      |
!        | *  Transform the Ritz values back to the original system.   |
!        |    For TYPE = 'SHIFTI' the transformation is                |
!        |             lambda = 1/theta + sigma                        |
!        |    For TYPE = 'BUCKLE' the transformation is                |
!        |             lambda = sigma * theta / ( theta - 1 )          |
!        |    For TYPE = 'CAYLEY' the transformation is                |
!        |             lambda = sigma * (theta + 1) / (theta - 1 )     |
!        |    where the theta are the Ritz values returned by dsaupd .  |
!        | NOTES:                                                      |
!        | *The Ritz vectors are not affected by the transformation.   |
!        |  They are only reordered.                                   |
!        %-------------------------------------------------------------%
!
         call dcopy  (ncv, workl(ihd), 1, workl(iw), 1)
         if (type .eq. 'SHIFTI') then 
            do 40 k=1, ncv
               workl(ihd+k-1) = one / workl(ihd+k-1) + sigma
  40        continue
         else if (type .eq. 'BUCKLE') then
            do 50 k=1, ncv
               workl(ihd+k-1) = sigma * workl(ihd+k-1) /                       &
                                (workl(ihd+k-1) - one)
  50        continue
         else if (type .eq. 'CAYLEY') then
            do 60 k=1, ncv
               workl(ihd+k-1) = sigma * (workl(ihd+k-1) + one) /               &
                                (workl(ihd+k-1) - one)
  60        continue
         end if
! 
!        %-------------------------------------------------------------%
!        | *  Store the wanted NCONV lambda values into D.             |
!        | *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   |
!        |    into ascending order and apply sort to the NCONV theta   |
!        |    values in the transformed system. We will need this to   |
!        |    compute Ritz estimates in the original system.           |
!        | *  Finally sort the lambda`s into ascending order and apply |
!        |    to Ritz vectors if wanted. Else just sort lambda`s into  |
!        |    ascending order.                                         |
!        | NOTES:                                                      |
!        | *workl(iw:iw+ncv-1) contain the theta ordered so that they  |
!        |  match the ordering of the lambda. We`ll use them again for |
!        |  Ritz vector purification.                                  |
!        %-------------------------------------------------------------%
!
         call dcopy (nconv, workl(ihd), 1, d, 1)
         call dsortr ('LA', .true., nconv, workl(ihd), workl(iw))
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
            call dscal (ncv, bnorm2/rnorm, workl(ihb), 1)
            call dsortr ('LA', .true., nconv, d, workl(ihb))
         end if
!
      end if 
! 
!     %------------------------------------------------%
!     | Compute the Ritz vectors. Transform the wanted |
!     | eigenvectors of the symmetric tridiagonal H by |
!     | the Lanczos basis matrix V.                    |
!     %------------------------------------------------%
!
      if (rvec .and. howmny .eq. 'A') then
!    
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(iq,ldq).                                |
!        %----------------------------------------------------------%
!     
         call dgeqr2 (ncv, nconv        , workl(iq) ,                          &
                      ldq, workl(iw+ncv), workl(ihb),                          &
                      ierr)
!
!        %--------------------------------------------------------%
!        | * Postmultiply V by Q.                                 |   
!        | * Copy the first NCONV columns of VQ into Z.           |
!        | The N by NCONV matrix Z is now a matrix representation |
!        | of the approximate invariant subspace associated with  |
!        | the Ritz values in workl(ihd).                         |
!        %--------------------------------------------------------%
!     
         call dorm2r ('Right', 'Notranspose', n        ,                       &
                      ncv    , nconv        , workl(iq),                       &
                      ldq    , workl(iw+ncv), v        ,                       &
                      ldv    , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
!
!        %-----------------------------------------------------%
!        | In order to compute the Ritz estimates for the Ritz |
!        | values in both systems, need the last row of the    |
!        | eigenvector matrix. Remember, it`s in factored form |
!        %-----------------------------------------------------%
!
         do 65 j = 1, ncv-1
            workl(ihb+j-1) = zero 
  65     continue
         workl(ihb+ncv-1) = one
         temp_array(1) = temp
         call dorm2r ('Left', 'Transpose'  , ncv       ,                       &
                      1     , nconv        , workl(iq) ,                       &
                      ldq   , workl(iw+ncv), workl(ihb),                       &
                      ncv   , TEMP_ARRAY   , ierr)
!
      else if (rvec .and. howmny .eq. 'S') then
!
!     Not yet implemented. See remark 2 above.
!
      end if
!
      if (type .eq. 'REGULR' .and. rvec) then
!
            do 70 j=1, ncv
               workl(ihb+j-1) = rnorm * abs( workl(ihb+j-1) )
 70         continue
!
      else if (type .ne. 'REGULR' .and. rvec) then
!
!        %-------------------------------------------------%
!        | *  Determine Ritz estimates of the theta.       |
!        |    If RVEC = .true. then compute Ritz estimates |
!        |               of the theta.                     |
!        |    If RVEC = .false. then copy Ritz estimates   |
!        |              as computed by dsaupd .             |
!        | *  Determine Ritz estimates of the lambda.      |
!        %-------------------------------------------------%
!
         call dscal  (ncv, bnorm2, workl(ihb), 1)
         if (type .eq. 'SHIFTI') then 
!
            do 80 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) )                          &
                              / workl(iw+k-1)**2
 80         continue
!
         else if (type .eq. 'BUCKLE') then
!
            do 90 k=1, ncv
               workl(ihb+k-1) = sigma * abs( workl(ihb+k-1) )                  &
                              / (workl(iw+k-1)-one )**2
 90         continue
!
         else if (type .eq. 'CAYLEY') then
!
            do 100 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1)                            &
                              / workl(iw+k-1)*(workl(iw+k-1)-one) )
 100        continue
!
         end if
!
      end if
!
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, d, ndigit,                                 &
                '_seupd: Untransformed converged Ritz values')
         call dvout (logfil, nconv, workl(ihb), ndigit,                        &
           '_seupd: Ritz estimates of the untransformed Ritz values')
      else if (msglvl .gt. 1) then
         call dvout (logfil, nconv, d, ndigit,                                 &
                '_seupd: Converged Ritz values')
         call dvout (logfil, nconv, workl(ihb), ndigit,                        &
           '_seupd: Associated Ritz estimates')
      end if
! 
!     %-------------------------------------------------%
!     | Ritz vector purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 3,4,5. See reference 7               |
!     %-------------------------------------------------%
!
      if (rvec .and. (type .eq. 'SHIFTI' .or. type .eq. 'CAYLEY')) then
!
         do 110 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)                                &
                        / workl(iw+k)
 110     continue
!
      else if (rvec .and. type .eq. 'BUCKLE') then
!
         do 120 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)                                &
                        / (workl(iw+k)-one)
 120     continue
!
      end if 
!
      if (type .ne. 'REGULR')                                                  &
         call dger  (n, nconv, one, resid, 1, workl(iw), 1, z, ldz)
!
 9000 continue
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+1) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!
!     %---------------%
!     | End of dseupd |
!     %---------------%
!
      end subroutine dseupd

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsaup2
!
!\Description: 
!  Intermediate level interface called by dsaupd.
!
!\Usage:
!  call dsaup2 
!     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, 
!       IPNTR, WORKD, INFO, PITERS )
!
!\Arguments
!
!  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dsaupd.
!  MODE, ISHIFT, MXITER: see the definition of IPARAM in dsaupd.
!  
!  NP      Integer.  (INPUT/OUTPUT)
!          Contains the number of implicit shifts to apply during 
!          each Arnoldi/Lanczos iteration.  
!          If ISHIFT=1, NP is adjusted dynamically at each iteration 
!          to accelerate convergence and prevent stagnation.
!          This is also roughly equal to the number of matrix-vector 
!          products (involving the operator OP) per Arnoldi iteration.
!          The logic for adjusting is contained within the current
!          subroutine.
!          If ISHIFT=0, NP is the number of shifts the user needs
!          to provide via reverse comunication. 0 < NP < NCV-NEV.
!          NP may be less than NCV-NEV since a leading block of the current
!          upper Tridiagonal matrix has split off and contains "unwanted"
!          Ritz values.
!          Upon termination of the IRA iteration, NP contains the number 
!          of "converged" wanted Ritz values.
!
!  IUPD    Integer.  (INPUT)
!          IUPD .EQ. 0: use explicit restart instead implicit update.
!          IUPD .NE. 0: use implicit update.
!
!  V       REAL(DOUBLE) N by (NEV+NP) array.  (INPUT/OUTPUT)
!          The Lanczos basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling 
!          program.
!
!  H       REAL(DOUBLE) (NEV+NP) by 2 array.  (OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          The subdiagonal is stored in the first column of H starting 
!          at H(2,1).  The main diagonal is stored in the second column
!          of H starting at H(1,2). If dsaup2 converges store the 
!          B-norm of the final residual vector in H(1,1).
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling 
!          program.
!
!  RITZ    REAL(DOUBLE) array of length NEV+NP.  (OUTPUT)
!          RITZ(1:NEV) contains the computed Ritz values of OP.
!
!  BOUNDS  REAL(DOUBLE) array of length NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
!
!  Q       REAL(DOUBLE) (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
!          Private (replicated) work array used to accumulate the 
!          rotation in the shift application step.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!          
!  WORKL   REAL(DOUBLE) array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  It is used in the computation of the 
!          tridiagonal eigenvalue problem, the calculation and
!          application of the shifts and convergence checking.
!          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
!          of WORKL are used in reverse communication to hold the user 
!          supplied shifts.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD for 
!          vectors used by the Lanczos iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in one of  
!                    the spectral transformation modes.  X is the current
!                    operand.
!          -------------------------------------------------------------
!          
!  WORKD   REAL(DOUBLE) work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Lanczos iteration
!          for reverse communication.  The user should not use WORKD
!          See Data Distribution Note in dsaupd.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =     0: Normal return.
!          =     1: All possible eigenvalues of OP has been found.  
!                   NP returns the size of the invariant subspace
!                   spanning the operator OP. 
!          =     2: No shifts could be applied.
!          =    -8: Error return from trid. eigenvalue calculation;
!                   This should never happen.
!          =    -9: Starting vector is zero.
!          = -9999: Could not build an Lanczos factorization.
!                   Size that was built in returned in NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Routines called:
!     dgetv0  ARPACK initial vector generation routine. 
!     dsaitr  ARPACK Lanczos factorization routine.
!     dsapps  ARPACK application of implicit shifts routine.
!     dsconv  ARPACK convergence of Ritz values routine.
!     dseigt  ARPACK compute Ritz values and error bounds routine.
!     dsgets  ARPACK reorder Ritz values and error bounds routine.
!     dsortr  ARPACK sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     second  ARPACK utility routine for timing.
!     dvout   ARPACK utility routine that prints vectors.
!     dlamch  LAPACK routine that determines machine constants.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     dscal   Level 1 BLAS that scales a vector.
!     dswap   Level 1 BLAS that swaps two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     12/15/93: Version ' 2.4'
!     xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)
!
!\SCCS Information: @(#) 
! FILE: saup2.F   SID: 2.7   DATE OF SID: 5/19/98   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsaup2                                                        &
         ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd,               &
           ishift, mxiter, v, ldv, h, ldh, ritz, bounds,                       &
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &     q, ldq, workl, ipntr, workd, info )
!    &     q, ldq, workl, ipntr, workd, info, nev_user, sigma, eig_old )
! E ///////////////////////////////////////////////////////////////////E
! B 06/08/08
! //////////////////////////////////////////////////////////B
           q, ldq, workl, ipntr, workd, info, nev_user, sigma, eig_old,        &
           iter, piters )   ! 11/30/09: add piters
! E ///////////////////////////////////////////////////////////////////E

! B 05/24/04 //////////////////////////////////////////////////////////B
      
      integer(long) n, numout, num_left, jj, jstart, nev_user
      real(double)  sigma
      real(double) eigout(n)
      real(double) eig_new(n)
      real(double) eig_old(n)
      real(double) eig_pc(n)

! E ///////////////////////////////////////////////////////////////////E

! B 06/08/08 //////////////////////////////////////////////////////////B
      integer iter_old
! E ///////////////////////////////////////////////////////////////////E

! B 11/30/09 //////////////////////////////////////////////////////////B
      integer piters
! E ///////////////////////////////////////////////////////////////////E

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAUP2'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, ldh, ldq, ldv, mxiter,               &
                 mode, nev, nev_array(1), np, np_array(1)
      REAL(DOUBLE)                                                             &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      REAL(DOUBLE)                                                             &
                 bounds(nev+np), h(ldh,2), q(ldq,nev+np), resid(n),            &
                 ritz(nev+np), v(ldv,nev+np), workd(3*n),                      &
                 workl(3*(nev+np))
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift
      integer    ierr, iter, j, kplusp, msglvl, nconv, nevbef, nev0,           &
                 np0, nptemp, nevd2, nevm2, kp(3),                             &
                 iter_array(1), nconv_array(1)
      REAL(DOUBLE)                                                             &
                 rnorm, temp, eps23
      REAL(DOUBLE)                                                             &
                 RNORM_ARRAY(1)
      save       cnorm, getv0, initv, update, ushift,                          &
! B 06/08/08 //////////////////////////////////////////////////////////B
!                iter, kplusp, msglvl, nconv, nev0, np0,                       &
                       kplusp, msglvl, nconv, nev0, np0,                       &
! E
! ///////////////////////////////////////////////////////////////////E
                 rnorm, eps23
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   dcopy, dgetv0, dsaitr, dscal, dsconv, dseigt, dsgets, 
!:!  &           dsapps, dsortr, dvout, ivout, second, dswap
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!:!   REAL(DOUBLE)
!:!  &           ddot, dnrm2, dlamch
!:!   external   ddot, dnrm2, dlamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+2) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      if (ido .eq. 0) then
! 
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call cpu_time (t0)
         msglvl = msaup2
!
!        %---------------------------------%
!        | Set machine dependent constant. |
!        %---------------------------------%
!
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0/3.0D+0)
!
!        %-------------------------------------%
!        | nev0 and np0 are integer variables  |
!        | hold the initial values of NEV & NP |
!        %-------------------------------------%
!
         nev0   = nev
         np0    = np
!
!        %-------------------------------------%
!        | kplusp is the bound on the largest  |
!        |        Lanczos factorization built. |
!        | nconv is the current number of      |
!        |        "converged" eigenvlues.      |
!        | iter is the counter on the current  |
!        |      iteration step.                |
!        %-------------------------------------%
!
         kplusp = nev0 + np0
         nconv  = 0
         iter   = 0
! 
!        %--------------------------------------------%
!        | Set flags for computing the first NEV steps |
!        | of the Lanczos factorization.              |
!        %--------------------------------------------%
!
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
!
         if (info .ne. 0) then
!
!        %--------------------------------------------%
!        | User provides the initial residual vector. |
!        %--------------------------------------------%
!
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
! 
!     %---------------------------------------------%
!     | Get a possibly random starting vector and   |
!     | force it into the range of the operator OP. |
!     %---------------------------------------------%
!
   10 continue
!
      if (getv0) then
         call dgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,         &
                      ipntr, workd, info)
!
         if (ido .ne. 99) go to 9000
!
         if (rnorm .eq. zero) then
!
!           %-----------------------------------------%
!           | The initial vector is zero. Error exit. | 
!           %-----------------------------------------%
!
            info = -9
            go to 1200
         end if
         getv0 = .false.
         ido  = 0
      end if
! 
!     %------------------------------------------------------------%
!     | Back from reverse communication: continue with update step |
!     %------------------------------------------------------------%
!
      if (update) go to 20
!
!     %-------------------------------------------%
!     | Back from computing user specified shifts |
!     %-------------------------------------------%
!
      if (ushift) go to 50
!
!     %-------------------------------------%
!     | Back from computing residual norm   |
!     | at the end of the current iteration |
!     %-------------------------------------%
!
      if (cnorm)  go to 100
! 
!     %----------------------------------------------------------%
!     | Compute the first NEV steps of the Lanczos factorization |
!     %----------------------------------------------------------%
!
      call dsaitr (ido, bmat, n, 0, nev0, mode, resid, rnorm, v, ldv,          &
                   h, ldh, ipntr, workd, info)
! 
!     %---------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication  |
!     | to compute operations involving OP and possibly B |
!     %---------------------------------------------------%
!
      if (ido .ne. 99) go to 9000
!
      if (info .gt. 0) then
!
!        %-----------------------------------------------------%
!        | dsaitr was unable to build an Lanczos factorization |
!        | of length NEV0. INFO is returned with the size of   |
!        | the factorization built. Exit main loop.            |
!        %-----------------------------------------------------%
!
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
! 
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |
!     |           Each iteration implicitly restarts the Lanczos     |
!     |           factorization in place.                            |
!     |                                                              |
!     %--------------------------------------------------------------%
! 
 1000 continue
!
         iter = iter + 1
!
         iter_array(1) = iter
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, iter_array, ndigit,                         &
                 '_saup2: **** Start of major iteration number ****')
         end if
         if (msglvl .gt. 1) then
            nev_array(1) = nev
            np_array(1)  = np
            call ivout (logfil, 1, nev_array, ndigit,                          &
           '_saup2: The length of the current Lanczos factorization')
            call ivout (logfil, 1, np_array, ndigit,                           &
                 '_saup2: Extend the Lanczos factorization by')
         end if
! 
!        %------------------------------------------------------------%
!        | Compute NP additional steps of the Lanczos factorization. |
!        %------------------------------------------------------------%
!
         ido = 0
   20    continue
         update = .true.
!
         call dsaitr (ido, bmat, n, nev, np, mode, resid, rnorm, v,            &
                      ldv, h, ldh, ipntr, workd, info)
! 
!        %---------------------------------------------------%
!        | ido .ne. 99 implies use of reverse communication  |
!        | to compute operations involving OP and possibly B |
!        %---------------------------------------------------%
!
         if (ido .ne. 99) go to 9000
!
         if (info .gt. 0) then
!
!           %-----------------------------------------------------%
!           | dsaitr was unable to build an Lanczos factorization |
!           | of length NEV0+NP0. INFO is returned with the size  |  
!           | of the factorization built. Exit main loop.         |
!           %-----------------------------------------------------%
!
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
!
         if (msglvl .gt. 1) then
            rnorm_array(1) = rnorm
            call dvout (logfil, 1, rnorm_array, ndigit,                        &
                 '_saup2: Current B-norm of residual for factorization')
         end if
! 
!        %--------------------------------------------------------%
!        | Compute the eigenvalues and corresponding error bounds |
!        | of the current symmetric tridiagonal matrix.           |
!        %--------------------------------------------------------%
!
         call dseigt (rnorm, kplusp, h, ldh, ritz, bounds, workl, ierr)
         
! B 05/24/04 //////////////////////////////////////////////////////////B
      if (piters > 0) then  ! 11/30/09: change debug(47) to piters

         do j=1,kplusp
            jj = kplusp - j + 1
            if (dabs(ritz(j)) > dlamch('s')) then
               if (mode == 2) then
                  eigout(jj) = sigma + one/ritz(j)
               else
                  eigout(jj) = ritz(j)
               endif
            else
               if (mode == 2) then
                  if (ritz(j) >= 0.d0) then
                     eigout(jj) = sigma + 1.0d0/dlamch('s')
                  else
                     eigout(jj) = sigma - 1.0d0/dlamch('s')
                  endif
               else
                  eigout(j) = ritz(j)
               endif
            endif
            eig_new(jj) = eigout(jj)
            eig_pc(jj)  = 100.d0*(eig_new(jj) - eig_old(jj))/eig_new(jj)
            eig_old(jj) = eig_new(jj)
         enddo

         numout = min(nev_user,kplusp)
         if (numout <= 10) then
            write(f06,99101)iter,(eigout(j),j=1,numout)
         else
            write(f06,99101)iter,(eigout(j),j=1,10)
            num_left = numout - 10
            jstart = 11
            do j=11,numout
               write(f06,99102) (eigout(jj),jj=jstart,jstart+9)
               jstart   = jstart + 10
               num_left = num_left - 10
               if (num_left < 10) exit
            enddo
            if (num_left > 0) then
               write(f06,99102) (eigout(j),j=jstart,numout)
            endif
            write(f06,*)
         endif

      endif

99101 format(' Eigenvalue estimates in Lanczos iteration ',i5,': ',            &
               10(1es14.6))

99102 format(48x,10(1es14.6))

!9103 format(48x,10(a,1es13.5))

! E////////////////////////////////////////////////////////////////////E
!
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
!
!        %----------------------------------------------------%
!        | Make a copy of eigenvalues and corresponding error |
!        | bounds obtained from _seigt.                       |
!        %----------------------------------------------------%
!
         call dcopy(kplusp, ritz, 1, workl(kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(2*kplusp+1), 1)
!
!        %---------------------------------------------------%
!        | Select the wanted Ritz values and their bounds    |
!        | to be used in the convergence test.               |
!        | The selection is based on the requested number of |
!        | eigenvalues instead of the current NEV and NP to  |
!        | prevent possible misconvergence.                  |
!        | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |
!        | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |
!        %---------------------------------------------------%
!
         nev = nev0
         np = np0
         call dsgets (ishift, which, nev, np, ritz, bounds, workl)
! 
!        %-------------------%
!        | Convergence test. |
!        %-------------------%
!
         call dcopy (nev, bounds(np+1), 1, workl(np+1), 1)
         call dsconv (nev, ritz(np+1), workl(np+1), tol, nconv)
!
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call ivout (logfil, 3, kp, ndigit,                                 &
                        '_saup2: NEV, NP, NCONV are')
            call dvout (logfil, kplusp, ritz, ndigit,                          &
                 '_saup2: The eigenvalues of H')
            call dvout (logfil, kplusp, bounds, ndigit,                        &
                '_saup2: Ritz estimates of the current NCV Ritz values')
         end if
!
!        %---------------------------------------------------------%
!        | Count the number of unwanted Ritz values that have zero |
!        | Ritz estimates. If any Ritz estimates are equal to zero |
!        | then a leading block of H of order equal to at least    |
!        | the number of Ritz values with zero Ritz estimates has  |
!        | split off. None of these Ritz values may be removed by  |
!        | shifting. Decrease NP the number of shifts to apply. If |
!        | no shifts may be applied, then prepare to exit          |
!        %---------------------------------------------------------%
!
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
! 
         if ( (nconv .ge. nev0) .or.                                           &
              (iter .gt. mxiter) .or.                                          &
              (np .eq. 0) ) then
!     
!           %------------------------------------------------%
!           | Prepare to exit. Put the converged Ritz values |
!           | and corresponding bounds in RITZ(1:NCONV) and  |
!           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
!           | careful when NCONV > NP since we don't want to |
!           | swap overlapping locations.                    |
!           %------------------------------------------------%
!
            if (which .eq. 'BE') then
!
!              %-----------------------------------------------------%
!              | Both ends of the spectrum are requested.            |
!              | Sort the eigenvalues into algebraically decreasing  |
!              | order first then swap low end of the spectrum next  |
!              | to high end in appropriate locations.               |
!              | NOTE: when np < floor(nev/2) be careful not to swap |
!              | overlapping locations.                              |
!              %-----------------------------------------------------%
!
               wprime = 'SA'
               call dsortr (wprime, .true., kplusp, ritz, bounds)
               nevd2 = nev0 / 2
               nevm2 = nev0 - nevd2 
               if ( nev .gt. 1 ) then
                  call dswap ( min(nevd2,np), ritz(nevm2+1), 1,                &
                       ritz( max(kplusp-nevd2+1,kplusp-np+1) ), 1)
                  call dswap ( min(nevd2,np), bounds(nevm2+1), 1,              &
                       bounds( max(kplusp-nevd2+1,kplusp-np+1)), 1)
               end if
!
            else
!
!              %--------------------------------------------------%
!              | LM, SM, LA, SA case.                             |
!              | Sort the eigenvalues of H into the an order that |
!              | is opposite to WHICH, and apply the resulting    |
!              | order to BOUNDS.  The eigenvalues are sorted so  |
!              | that the wanted part are always within the first |
!              | NEV locations.                                   |
!              %--------------------------------------------------%
!
               if (which .eq. 'LM') wprime = 'SM'
               if (which .eq. 'SM') wprime = 'LM'
               if (which .eq. 'LA') wprime = 'SA'
               if (which .eq. 'SA') wprime = 'LA'
!
               call dsortr (wprime, .true., kplusp, ritz, bounds)
!
            end if
!
!           %--------------------------------------------------%
!           | Scale the Ritz estimate of each Ritz value       |
!           | by 1 / max(eps23,magnitude of the Ritz value).   |
!           %--------------------------------------------------%
!
            do 35 j = 1, nev0
               temp = max( eps23, abs(ritz(j)) )
               bounds(j) = bounds(j)/temp
 35         continue
!
!           %----------------------------------------------------%
!           | Sort the Ritz values according to the scaled Ritz  |
!           | esitmates.  This will push all the converged ones  |
!           | towards the front of ritzr, ritzi, bounds          |
!           | (in the case when NCONV < NEV.)                    |
!           %----------------------------------------------------%
!
            wprime = 'LA'
            call dsortr(wprime, .true., nev0, bounds, ritz)
!
!           %----------------------------------------------%
!           | Scale the Ritz estimate back to its original |
!           | value.                                       |
!           %----------------------------------------------%
!
            do 40 j = 1, nev0
                temp = max( eps23, abs(ritz(j)) )
                bounds(j) = bounds(j)*temp
 40         continue
!
!           %--------------------------------------------------%
!           | Sort the "converged" Ritz values again so that   |
!           | the "threshold" values and their associated Ritz |
!           | estimates appear at the appropriate position in  |
!           | ritz and bound.                                  |
!           %--------------------------------------------------%
!
            if (which .eq. 'BE') then
!
!              %------------------------------------------------%
!              | Sort the "converged" Ritz values in increasing |
!              | order.  The "threshold" values are in the      |
!              | middle.                                        |
!              %------------------------------------------------%
!
               wprime = 'LA'
               call dsortr(wprime, .true., nconv, ritz, bounds)
!
            else
!
!              %----------------------------------------------%
!              | In LM, SM, LA, SA case, sort the "converged" |
!              | Ritz values according to WHICH so that the   |
!              | "threshold" value appears at the front of    |
!              | ritz.                                        |
!              %----------------------------------------------%

               call dsortr(which, .true., nconv, ritz, bounds)
!
            end if
!
!           %------------------------------------------%
!           |  Use h( 1,1 ) as storage to communicate  |
!           |  rnorm to _seupd if needed               |
!           %------------------------------------------%
!
            h(1,1) = rnorm
!
            if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritz, ndigit,                       &
                  '_saup2: Sorted Ritz values.')
               call dvout (logfil, kplusp, bounds, ndigit,                     &
                  '_saup2: Sorted ritz estimates.')
            end if
!
!           %------------------------------------%
!           | Max iterations have been exceeded. | 
!           %------------------------------------%
!
            if (iter .gt. mxiter .and. nconv .lt. nev) info = 1
!
!           %---------------------%
!           | No shifts to apply. | 
!           %---------------------%
!
            if (np .eq. 0 .and. nconv .lt. nev0) info = 2
!
            np = nconv
            go to 1100
!
         else if (nconv .lt. nev .and. ishift .eq. 1) then
!
!           %---------------------------------------------------%
!           | Do not have all the requested eigenvalues yet.    |
!           | To prevent possible stagnation, adjust the number |
!           | of Ritz values and the shifts.                    |
!           %---------------------------------------------------%
!
            nevbef = nev
            nev = nev + min (nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 2) then
               nev = 2
            end if
            np  = kplusp - nev
!     
!           %---------------------------------------%
!           | If the size of NEV was just increased |
!           | resort the eigenvalues.               |
!           %---------------------------------------%
!     
            if (nevbef .lt. nev)                                               &
               call dsgets (ishift, which, nev, np, ritz, bounds,              &
                    workl)
!
         end if
!
         if (msglvl .gt. 0) then
            nconv_array(1) = nconv
            call ivout (logfil, 1, nconv_array, ndigit,                        &
                 '_saup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit,                              &
                    '_saup2: NEV and NP are')
               call dvout (logfil, nev, ritz(np+1), ndigit,                    &
                    '_saup2: "wanted" Ritz values.')
               call dvout (logfil, nev, bounds(np+1), ndigit,                  &
                    '_saup2: Ritz estimates of the "wanted" values ')
            end if
         end if

! 
         if (ishift .eq. 0) then
!
!           %-----------------------------------------------------%
!           | User specified shifts: reverse communication to     |
!           | compute the shifts. They are returned in the first  |
!           | NP locations of WORKL.                              |
!           %-----------------------------------------------------%
!
            ushift = .true.
            ido = 3
            go to 9000
         end if
!
   50    continue
!
!        %------------------------------------%
!        | Back from reverse communication;   |
!        | User specified shifts are returned |
!        | in WORKL(1:*NP)                   |
!        %------------------------------------%
!
         ushift = .false.
! 
! 
!        %---------------------------------------------------------%
!        | Move the NP shifts to the first NP locations of RITZ to |
!        | free up WORKL.  This is for the non-exact shift case;   |
!        | in the exact shift case, dsgets already handles this.   |
!        %---------------------------------------------------------%
!
         if (ishift .eq. 0) call dcopy (np, workl, 1, ritz, 1)
!
         if (msglvl .gt. 2) then
            np_array(1) = np
            call ivout (logfil, 1, np_array, ndigit,                           &
                        '_saup2: The number of shifts to apply ')
            call dvout (logfil, np, workl, ndigit,                             &
                        '_saup2: shifts selected')
            if (ishift .eq. 1) then
               call dvout (logfil, np, bounds, ndigit,                         &
                        '_saup2: corresponding Ritz estimates')
             end if
         end if
! 
!        %---------------------------------------------------------%
!        | Apply the NP0 implicit shifts by QR bulge chasing.      |
!        | Each shift is applied to the entire tridiagonal matrix. |
!        | The first 2*N locations of WORKD are used as workspace. |
!        | After dsapps is done, we have a Lanczos                 |
!        | factorization of length NEV.                            |
!        %---------------------------------------------------------%
!
         call dsapps (n, nev, np, ritz, v, ldv, h, ldh, resid, q, ldq,         &
              workd)
!
!        %---------------------------------------------%
!        | Compute the B-norm of the updated residual. |
!        | Keep B*RESID in WORKD(1:N) to be used in    |
!        | the first step of the next call to dsaitr.  |
!        %---------------------------------------------%
!
         cnorm = .true.
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
! 
!           %----------------------------------%
!           | Exit in order to compute B*RESID |
!           %----------------------------------%
! 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd, 1)
         end if
! 
  100    continue
! 
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(1:N) := B*RESID            |
!        %----------------------------------%
!
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
! 
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
  130    continue
!
         if (msglvl .gt. 2) then
            rnorm_array(1) = rnorm
            call dvout (logfil, 1, rnorm_array, ndigit,                        &
            '_saup2: B-norm of residual for NEV factorization')
            call dvout (logfil, nev, h(1,2), ndigit,                           &
                 '_saup2: main diagonal of compressed H matrix')
            call dvout (logfil, nev-1, h(2,1), ndigit,                         &
                 '_saup2: subdiagonal of compressed H matrix')
         end if
! 
      go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
! 
 1100 continue
!
      mxiter = iter
      nev = nconv
! 
 1200 continue
      ido = 99
! B 05/24/04 //////////////////////////////////////////////////////////B
                                 ! 11/30/09: change debug(47) to piters
      if ((piters > 0) .and. (mxiter > 1)) then
         write(f06,*)
         if (numout <= 10) then
            write(f06,99201) (eig_pc(j),j=1,numout)
         else
            write(f06,99201 ) (eig_pc(j),j=1,10)
            num_left = numout - 10
            jstart = 11
            do j=11,numout
               write(f06,99202) (eig_pc(jj),jj=jstart,jstart+9)
               jstart   = jstart + 10
               num_left = num_left - 10
               if (num_left < 10) exit
            enddo
            if (num_left > 0) then
               write(f06,99202) (eig_pc(j),j=jstart,numout)
            endif
         endif
      endif

99201 format(' % diffs in eigenvals from last 2 iterations    : ',             &
               10(1es14.6))

99202 format(48x,10(1es14.6))

! E////////////////////////////////////////////////////////////////////E
!
!     %------------%
!     | Error exit |
!     %------------%
!
      call cpu_time (t1)
      tsaup2 = t1 - t0
! 
 9000 continue

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+2) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************
!
!     %---------------%
!     | End of dsaup2 |
!     %---------------%
!
      end subroutine dsaup2

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsesrt
!
!\Description:
!  Sort the array X in the order specified by WHICH and optionally 
!  apply the permutation to the columns of the matrix A.
!
!\Usage:
!  call dsesrt
!     ( WHICH, APPLY, N, X, NA, A, LDA)
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X is sorted into increasing order of magnitude.
!          'SM' -> X is sorted into decreasing order of magnitude.
!          'LA' -> X is sorted into increasing order of algebraic.
!          'SA' -> X is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to A.
!          APPLY = .FALSE. -> do not apply the sorted order to A.
!
!  N       Integer.  (INPUT)
!          Dimension of the array X.
!
!  X      REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  NA      Integer.  (INPUT)
!          Number of rows of the matrix A.
!
!  A      REAL(DOUBLE) array of length NA by N.  (INPUT/OUTPUT)
!         
!  LDA     Integer.  (INPUT)
!          Leading dimension of A.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines
!     dswap  Level 1 BLAS that swaps the contents of two vectors.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     12/15/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO and 
!               the ARPACK code dsortr
!
!\SCCS Information: @(#) 
! FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsesrt (which, apply, n, x, na, a, lda)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSESRT'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      logical    apply
      integer    lda, n, na
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      REAL(DOUBLE)                                                             &
                 x(0:n-1), a(lda, 0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, igap, j
      REAL(DOUBLE)                                                             &
                 temp
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   dswap
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+2) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      igap = n / 2
! 
      if (which .eq. 'SA') then
!
!        X is sorted into decreasing order of algebraic.
!
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if (j.lt.0) go to 30
!
            if (x(j).lt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
!
      else if (which .eq. 'SM') then
!
!        X is sorted into decreasing order of magnitude.
!
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j.lt.0) go to 60
!
            if (abs(x(j)).lt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
!
      else if (which .eq. 'LA') then
!
!        X is sorted into increasing order of algebraic.
!
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j.lt.0) go to 90
!           
            if (x(j).gt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
! 
      else if (which .eq. 'LM') then
!
!        X is sorted into increasing order of magnitude.
!
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j.lt.0) go to 120
!
            if (abs(x(j)).gt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
!
 9000 continue

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+2) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************
!
!     %---------------%
!     | End of dsesrt |
!     %---------------%
!
      end subroutine dsesrt

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsgets
!
!\Description: 
!  Given the eigenvalues of the symmetric tridiagonal matrix H,
!  computes the NP shifts AMU that are zeros of the polynomial of 
!  degree NP which filters out components of the unwanted eigenvectors 
!  corresponding to the AMU's based on some given criteria.
!
!  NOTE: This is called even in the case of user specified shifts in 
!  order to sort the eigenvalues, and error bounds of H for later use.
!
!\Usage:
!  call dsgets
!     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
!
!\Arguments
!  ISHIFT  Integer.  (INPUT)
!          Method for selecting the implicit shifts at each iteration.
!          ISHIFT = 0: user specified shifts
!          ISHIFT = 1: exact shift with respect to the matrix H.
!
!  WHICH   Character*2.  (INPUT)
!          Shift selection criteria.
!          'LM' -> KEV eigenvalues of largest magnitude are retained.
!          'SM' -> KEV eigenvalues of smallest magnitude are retained.
!          'LA' -> KEV eigenvalues of largest value are retained.
!          'SA' -> KEV eigenvalues of smallest value are retained.
!          'BE' -> KEV eigenvalues, half from each end of the spectrum.
!                  If KEV is odd, compute one more from the high end.
!
!  KEV      Integer.  (INPUT)
!          KEV+NP is the size of the matrix H.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be computed.
!
!  RITZ    REAL(DOUBLE) array of length KEV+NP.  (INPUT/OUTPUT)
!          On INPUT, RITZ contains the eigenvalues of H.
!          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues 
!          are in the first NP locations and the wanted part is in 
!          the last KEV locations.  When exact shifts are selected, the
!          unwanted part corresponds to the shifts to be applied.
!
!  BOUNDS  REAL(DOUBLE) array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!  SHIFTS  REAL(DOUBLE) array of length NP.  (INPUT/OUTPUT)
!          On INPUT:  contains the user specified shifts if ISHIFT = 0.
!          On OUTPUT: contains the shifts sorted into decreasing order 
!          of magnitude with respect to the Ritz estimates contained in
!          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     dsortr  ARPACK utility sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     second  ARPACK utility routine for timing.
!     dvout   ARPACK utility routine that prints vectors.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     dswap   Level 1 BLAS that swaps the contents of two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     xx/xx/93: Version ' 2.1'
!
!\SCCS Information: @(#) 
! FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsgets ( ishift, which, kev, np, ritz, bounds, shifts )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSGETS'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      integer    ishift, kev, np
      INTEGER    KEV_ARRAY(1), NP_ARRAY(1)
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      REAL(DOUBLE)                                                             &
                 bounds(kev+np), ritz(kev+np), shifts(np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    kevd2, msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   dswap, dcopy, dsortr, second
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    max, min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
! 
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = msgets
! 
      if (which .eq. 'BE') then
!
!        %-----------------------------------------------------%
!        | Both ends of the spectrum are requested.            |
!        | Sort the eigenvalues into algebraically increasing  |
!        | order first then swap high end of the spectrum next |
!        | to low end in appropriate locations.                |
!        | NOTE: when np < floor(kev/2) be careful not to swap |
!        | overlapping locations.                              |
!        %-----------------------------------------------------%
!
         call dsortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2 
         if ( kev .gt. 1 ) then
            call dswap ( min(kevd2,np), ritz, 1,                               &
                         ritz( max(kevd2,np)+1 ), 1)
            call dswap ( min(kevd2,np), bounds, 1,                             &
                         bounds( max(kevd2,np)+1 ), 1)
         end if
!
      else
!
!        %----------------------------------------------------%
!        | LM, SM, LA, SA case.                               |
!        | Sort the eigenvalues of H into the desired order   |
!        | and apply the resulting order to BOUNDS.           |
!        | The eigenvalues are sorted so that the wanted part |
!        | are always in the last KEV locations.               |
!        %----------------------------------------------------%
!
         call dsortr (which, .true., kev+np, ritz, bounds)
      end if
!
      if (ishift .eq. 1 .and. np .gt. 0) then
!     
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first.       |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when the shifts  |
!        | are applied in subroutine dsapps.                     |
!        %-------------------------------------------------------%
!     
         call dsortr ('SM', .true., np, bounds, ritz)
         call dcopy (np, ritz, 1, shifts, 1)
      end if
! 
      call cpu_time (t1)
      tsgets = tsgets + (t1 - t0)
!
      if (msglvl .gt. 0) then
         kev_array(1) = kev
         np_array(1)  = np
         call ivout (logfil, 1, KEV_ARRAY, ndigit, '_sgets: KEV is')
         call ivout (logfil, 1, NP_ARRAY, ndigit, '_sgets: NP is')
         call dvout (logfil, kev+np, ritz, ndigit,                             &
              '_sgets: Eigenvalues of current H matrix')
         call dvout (logfil, kev+np, bounds, ndigit,                           &
              '_sgets: Associated Ritz estimates')
      end if
! 
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************
!
!     %---------------%
!     | End of dsgets |
!     %---------------%
!
      end subroutine dsgets

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsortr
!
!\Description:
!  Sort the array X1 in the order specified by WHICH and optionally 
!  applies the permutation to the array X2.
!
!\Usage:
!  call dsortr
!     ( WHICH, APPLY, N, X1, X2 )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X1 is sorted into increasing order of magnitude.
!          'SM' -> X1 is sorted into decreasing order of magnitude.
!          'LA' -> X1 is sorted into increasing order of algebraic.
!          'SA' -> X1 is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to X2.
!          APPLY = .FALSE. -> do not apply the sorted order to X2.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  X1      REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  X2      REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
!          Only referenced if APPLY = .TRUE.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     12/16/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#) 
! FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsortr (which, apply, n, x1, x2)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSORTR'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      logical    apply
      integer    n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      REAL(DOUBLE)                                                             &
                 x1(0:n-1), x2(0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, igap, j
      REAL(DOUBLE)                                                             &
                 temp
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      igap = n / 2
! 
      if (which .eq. 'SA') then
!
!        X1 is sorted into decreasing order of algebraic.
!
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if (j.lt.0) go to 30
!
            if (x1(j).lt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
!
      else if (which .eq. 'SM') then
!
!        X1 is sorted into decreasing order of magnitude.
!
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j.lt.0) go to 60
!
            if (abs(x1(j)).lt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
!
      else if (which .eq. 'LA') then
!
!        X1 is sorted into increasing order of algebraic.
!
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j.lt.0) go to 90
!           
            if (x1(j).gt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
! 
      else if (which .eq. 'LM') then
!
!        X1 is sorted into increasing order of magnitude.
!
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j.lt.0) go to 120
!
            if (abs(x1(j)).gt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
!
 9000 continue
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!     %---------------%
!     | End of dsortr |
!     %---------------%
!
      end subroutine dsortr

! ################################################################################################################################## 

!
!\SCCS Information: @(#) 
! FILE: stats.F   SID: 2.1   DATE OF SID: 4/19/96   RELEASE: 2
!     %---------------------------------------------%
!     | Initialize statistic and timing information |
!     | for symmetric Arnoldi code.                 |
!     %---------------------------------------------%
 
      subroutine dstats

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSTATS'

!     %--------------------------------%
!     | See stat.doc for documentation |
!     %--------------------------------%
      include   'stat.h'
 
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+2) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
 
      tsaupd = 0.0D+0
      tsaup2 = 0.0D+0
      tsaitr = 0.0D+0
      tseigt = 0.0D+0
      tsgets = 0.0D+0
      tsapps = 0.0D+0
      tsconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0
 
!     %----------------------------------------------------%
!     | User time including reverse communication overhead |
!     %----------------------------------------------------%
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0
 
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+2) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************
!
!     End of dstats
!
      end subroutine dstats

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dgetv0
!
!\Description: 
!  Generate a random initial residual vector for the Arnoldi process.
!  Force the residual vector to be in the range of the operator OP.  
!
!\Usage:
!  call dgetv0
!     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, 
!       IPNTR, WORKD, IERR )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to dgetv0.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B in the (generalized)
!          eigenvalue problem A*x = lambda*B*x.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  ITRY    Integer.  (INPUT)
!          ITRY counts the number of times that dgetv0 is called.  
!          It should be set to 1 on the initial call to dgetv0.
!
!  INITV   Logical variable.  (INPUT)
!          .TRUE.  => the initial residual vector is given in RESID.
!          .FALSE. => generate a random initial residual vector.
!
!  N       Integer.  (INPUT)
!          Dimension of the problem.
!
!  J       Integer.  (INPUT)
!          Index of the residual vector to be generated, with respect to
!          the Arnoldi process.  J > 1 in case of a "restart".
!
!  V       REAL(DOUBLE) N by J array.  (INPUT)
!          The first J-1 columns of V contain the current Arnoldi basis
!          if this is a "restart".
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling 
!          program.
!
!  RESID   REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
!          Initial residual vector to be generated.  If RESID is 
!          provided, force RESID into the range of the operator OP.
!
!  RNORM   REAL(DOUBLE) scalar.  (OUTPUT)
!          B-norm of the generated residual.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!
!  WORKD   REAL(DOUBLE) work array of length 2*N.  (REVERSE COMMUNICATION).
!          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
!
!  IERR    Integer.  (OUTPUT)
!          =  0: Normal exit.
!          = -1: Cannot generate a nontrivial restarted residual vector
!                in the range of the operator OP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     second  ARPACK utility routine for timing.
!     dvout   ARPACK utility routine for vector output.
!     dlarnv  LAPACK routine for generating a random vector.
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\SCCS Information: @(#) 
! FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dgetv0                                                        &
         ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm,                 &
           ipntr, workd, ierr )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DGETV0'

! 
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      REAL(DOUBLE)                                                             &
                 rnorm
      REAL(DOUBLE)  RNORM_ARRAY(1)
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      REAL(DOUBLE)                                                             &
                 resid(n), v(ldv,j), workd(2*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    first, inits, orth
      integer    idist, iseed(4), iter, msglvl, jj
      REAL(DOUBLE)                                                             &
                 rnorm0
      REAL(DOUBLE)  RNORM0_ARRAY(1)
      save       first, iseed, inits, iter, msglvl, orth, rnorm0
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   dlarnv, dvout, dcopy, dgemv, second
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!:!   REAL(DOUBLE)
!:!  &           ddot, dnrm2
!:!   external   ddot, dnrm2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %-----------------%
!     | Data Statements |
!     %-----------------%
!
      data       inits /.true./
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************
!
!     %-----------------------------------%
!     | Initialize the seed of the LAPACK |
!     | random number generator           |
!     %-----------------------------------%
!
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      end if
!
      if (ido .eq.  0) then
! 
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call cpu_time (t0)
         msglvl = mgetv0
! 
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
!
!        %-----------------------------------------------------%
!        | Possibly generate a random starting vector in RESID |
!        | Use a LAPACK random number generator used by the    |
!        | matrix generation routines.                         |
!        |    idist = 1: uniform (0,1)  distribution;          |
!        |    idist = 2: uniform (-1,1) distribution;          |
!        |    idist = 3: normal  (0,1)  distribution;          |
!        %-----------------------------------------------------%
!
         if (.not.initv) then
            idist = 2
            call dlarnv (idist, iseed, n, resid)
         end if
! 
!        %----------------------------------------------------------%
!        | Force the starting vector into the range of OP to handle |
!        | the generalized problem when B is possibly (singular).   |
!        %----------------------------------------------------------%
!
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call dcopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
! 
!     %-----------------------------------------%
!     | Back from computing OP*(initial-vector) |
!     %-----------------------------------------%
!
      if (first) go to 20
!
!     %-----------------------------------------------%
!     | Back from computing B*(orthogonalized-vector) |
!     %-----------------------------------------------%
!
      if (orth)  go to 40
! 
      if (bmat .eq. 'G') then
         call cpu_time (t3)
         tmvopx = tmvopx + (t3 - t2)
      end if
! 
!     %------------------------------------------------------%
!     | Starting vector is now in the range of OP; r = OP*r; |
!     | Compute B-norm of starting vector.                   |
!     %------------------------------------------------------%
!
      call cpu_time (t2)
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
! 
   20 continue
!
      if (bmat .eq. 'G') then
         call cpu_time (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
! 
      first = .FALSE.
      if (bmat .eq. 'G') then
          rnorm0 = ddot (n, resid, 1, workd, 1)
          rnorm0 = sqrt(abs(rnorm0))
      else if (bmat .eq. 'I') then
           rnorm0 = dnrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
!
!     %---------------------------------------------%
!     | Exit if this is the very first Arnoldi step |
!     %---------------------------------------------%
!
      if (j .eq. 1) go to 50
! 
!     %----------------------------------------------------------------
!     | Otherwise need to B-orthogonalize the starting vector against |
!     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
!     | This is the case where an invariant subspace is encountered   |
!     | in the middle of the Arnoldi factorization.                   |
!     |                                                               |
!     |       s = V^{T}*B*r;   r = r - V*s;                           |
!     |                                                               |
!     | Stopping criteria used for iter. ref. is discussed in         |
!     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
!     %---------------------------------------------------------------%
!
      orth = .TRUE.
   30 continue
!
      call dgemv ('T', n, j-1, one, v, ldv, workd, 1,                          &
                  zero, workd(n+1), 1)
      call dgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1,                    &
                  one, resid, 1)
! 
!     %----------------------------------------------------------%
!     | Compute the B-norm of the orthogonalized starting vector |
!     %----------------------------------------------------------%
!
      call cpu_time (t2)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
! 
   40 continue
!
      if (bmat .eq. 'G') then
         call cpu_time (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
! 
      if (bmat .eq. 'G') then
         rnorm = ddot (n, resid, 1, workd, 1)
         rnorm = sqrt(abs(rnorm))
      else if (bmat .eq. 'I') then
         rnorm = dnrm2(n, resid, 1)
      end if
!
!     %--------------------------------------%
!     | Check for further orthogonalization. |
!     %--------------------------------------%
!
      if (msglvl .gt. 2) then
          rnorm0_array(1) = rnorm0
          rnorm_array(1)  = rnorm
          call dvout (logfil, 1, RNORM0_ARRAY, ndigit,                         &
                      '_getv0: re-orthonalization ; rnorm0 is')
          call dvout (logfil, 1, RNORM_ARRAY, ndigit,                          &
                      '_getv0: re-orthonalization ; rnorm is')
      end if
!
      if (rnorm .gt. 0.717*rnorm0) go to 50
! 
      iter = iter + 1
      if (iter .le. 5) then
!
!        %-----------------------------------%
!        | Perform iterative refinement step |
!        %-----------------------------------%
!
         rnorm0 = rnorm
         go to 30
      else
!
!        %------------------------------------%
!        | Iterative refinement step "failed" |
!        %------------------------------------%
!
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = zero
         ierr = -1
      end if
! 
   50 continue
!
      if (msglvl .gt. 0) then
         rnorm_array(1) = rnorm
         call dvout (logfil, 1, RNORM_ARRAY, ndigit,                           &
              '_getv0: B-norm of initial / restarted starting vector')
      end if
      if (msglvl .gt. 3) then
         call dvout (logfil, n, resid, ndigit,                                 &
              '_getv0: initial / restarted starting vector')
      end if
      ido = 99
! 
      call cpu_time (t1)
      tgetv0 = tgetv0 + (t1 - t0)
! 
 9000 continue
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!     %---------------%
!     | End of dgetv0 |
!     %---------------%
!
      end subroutine dgetv0

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsaitr
!
!\Description: 
!  Reverse communication interface for applying NP additional steps to 
!  a K step symmetric Arnoldi factorization.
!
!  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
!
!          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
!
!  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
!
!          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
!
!  where OP and B are as in dsaupd.  The B-norm of r_{k+p} is also
!  computed and returned.
!
!\Usage:
!  call dsaitr
!     ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH, 
!       IPNTR, WORKD, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!                    This is for the restart phase to force the new
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y,
!                    IPNTR(3) is the pointer into WORK for B * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!          When the routine is used in the "shift-and-invert" mode, the
!          vector B * Q is already available and does not need to be
!          recomputed in forming OP * Q.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of matrix B that defines the
!          semi-inner product for the operator OP.  See dsaupd.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  K       Integer.  (INPUT)
!          Current order of H and the number of columns of V.
!
!  NP      Integer.  (INPUT)
!          Number of additional Arnoldi steps to take.
!
!  MODE    Integer.  (INPUT)
!          Signifies which form for "OP". If MODE=2 then
!          a reduction in the number of B matrix vector multiplies
!          is possible since the B-norm of OP*x is equivalent to
!          the inv(B)-norm of A*x.
!
!  RESID   REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
!          On INPUT:  RESID contains the residual vector r_{k}.
!          On OUTPUT: RESID contains the residual vector r_{k+p}.
!
!  RNORM   REAL(DOUBLE) scalar.  (INPUT/OUTPUT)
!          On INPUT the B-norm of r_{k}.
!          On OUTPUT the B-norm of the updated residual r_{k+p}.
!
!  V       REAL(DOUBLE) N by K+NP array.  (INPUT/OUTPUT)
!          On INPUT:  V contains the Arnoldi vectors in the first K 
!          columns.
!          On OUTPUT: V contains the new NP Arnoldi vectors in the next
!          NP columns.  The first K columns are unchanged.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling 
!          program.
!
!  H       REAL(DOUBLE) (K+NP) by 2 array.  (INPUT/OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          with the subdiagonal in the first column starting at H(2,1)
!          and the main diagonal in the second column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling 
!          program.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORK for 
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the 
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!          
!  WORKD   REAL(DOUBLE) work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The calling program should not 
!          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K step Arnoldi factorization. Used to save some 
!          computation at the first step. 
!          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K+NP step Arnoldi factorization.
!
!  INFO    Integer.  (OUTPUT)
!          = 0: Normal exit.
!          > 0: Size of an invariant subspace of OP is found that is
!               less than K + NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     dgetv0  ARPACK routine to generate the initial vector.
!     ivout   ARPACK utility routine that prints integers.
!     dmout   ARPACK utility routine that prints matrices.
!     dvout   ARPACK utility routine that prints vectors.
!     dlamch  LAPACK routine that determines machine constants.
!     dlascl  LAPACK routine for careful scaling of a matrix.
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     daxpy   Level 1 BLAS that computes a vector triad.
!     dscal   Level 1 BLAS that scales a vector.
!     dcopy   Level 1 BLAS that copies one vector to another .
!     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     xx/xx/93: Version ' 2.4'
!
!\SCCS Information: @(#) 
! FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2
!
!\Remarks
!  The algorithm implemented is:
!  
!  restart = .false.
!  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
!  r_{k} contains the initial residual vector even for k = 0;
!  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
!  computed by the calling program.
!
!  betaj = rnorm ; p_{k+1} = B*r_{k} ;
!  For  j = k+1, ..., k+np  Do
!     1) if ( betaj < tol ) stop or restart depending on j.
!        if ( restart ) generate a new starting vector.
!     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
!        p_{j} = p_{j}/betaj
!     3) r_{j} = OP*v_{j} where OP is defined as in dsaupd
!        For shift-invert mode p_{j} = B*v_{j} is already available.
!        wnorm = || OP*v_{j} ||
!     4) Compute the j-th step residual vector.
!        w_{j} =  V_{j}^T * B * OP * v_{j}
!        r_{j} =  OP*v_{j} - V_{j} * w_{j}
!        alphaj <- j-th component of w_{j}
!        rnorm = || r_{j} ||
!        betaj+1 = rnorm
!        If (rnorm > 0.717*wnorm) accept step and go back to 1)
!     5) Re-orthogonalization step:
!        s = V_{j}'*B*r_{j}
!        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
!        alphaj = alphaj + s_{j};   
!     6) Iterative refinement step:
!        If (rnorm1 > 0.717*rnorm) then
!           rnorm = rnorm1
!           accept step and go back to 1)
!        Else
!           rnorm = rnorm1
!           If this is the first time in step 6), go to 5)
!           Else r_{j} lies in the span of V_{j} numerically.
!              Set r_{j} = 0 and rnorm = 0; go to 1)
!        EndIf 
!  End Do
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsaitr                                                        &
         (ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh,             &
          ipntr, workd, info)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAITR'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, mode, np
      REAL(DOUBLE)                                                             &
                 rnorm
      REAL(DOUBLE) RNORM_ARRAY(1)
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      REAL(DOUBLE)                                                             &
                 h(ldh,2), resid(n), v(ldv,k+np), workd(3*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      logical    first, orth1, orth2, rstart, step3, step4
      integer    i, ierr, ipj, irj, ivj, iter, itry, j, msglvl,                &
                 infol, jj
      INTEGER    J_ARRAY(1)
      REAL(DOUBLE)                                                             &
                 rnorm1, wnorm, safmin, temp1
      save       orth1, orth2, rstart, step3, step4,                           &
                 ierr, ipj, irj, ivj, iter, itry, j, msglvl,                   &
                 rnorm1, safmin, wnorm
!
!     %-----------------------%
!     | Local Array Arguments | 
!     %-----------------------%
!
      REAL(DOUBLE)                                                             &
                 xtemp(2)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   daxpy, dcopy, dscal, dgemv, dgetv0, dvout, dmout,
!:!  &           dlascl, ivout, second
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!:!   REAL(DOUBLE)
!:!  &           ddot, dnrm2, dlamch
!:!   external   ddot, dnrm2, dlamch
!
!     %-----------------%
!     | Data statements |
!     %-----------------%
!
      data      first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      if (first) then
         first = .false.
!
!        %--------------------------------%
!        | safmin = safe minimum is such  |
!        | that 1/sfmin does not overflow |
!        %--------------------------------%
!
         safmin = dlamch('safmin')
      end if
!
      if (ido .eq. 0) then
! 
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call cpu_time (t0)
         msglvl = msaitr
! 
!        %------------------------------%
!        | Initial call to this routine |
!        %------------------------------%
!
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
! 
!        %--------------------------------%
!        | Pointer to the current step of |
!        | the factorization to build     |
!        %--------------------------------%
!
         j      = k + 1
! 
!        %------------------------------------------%
!        | Pointers used for reverse communication  |
!        | when using WORKD.                        |
!        %------------------------------------------%
!
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
! 
!     %-------------------------------------------------%
!     | When in reverse communication mode one of:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | will be .true.                                  |
!     | STEP3: return from computing OP*v_{j}.          |
!     | STEP4: return from computing B-norm of OP*v_{j} |
!     | ORTH1: return from computing B-norm of r_{j+1}  |
!     | ORTH2: return from computing B-norm of          |
!     |        correction to the residual vector.       |
!     | RSTART: return from OP computations needed by   |
!     |         dgetv0.                                 |
!     %-------------------------------------------------%
!
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
!
!     %------------------------------%
!     | Else this is the first step. |
!     %------------------------------%
! 
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%
!
 1000 continue
!
         if (msglvl .gt. 2) then
            j_array(1)     = j
            rnorm_array(1) = rnorm
            call ivout (logfil, 1, J_ARRAY, ndigit,                            &
                        '_saitr: generating Arnoldi vector no.')
            call dvout (logfil, 1, RNORM_ARRAY, ndigit,                        &
                        '_saitr: B-norm of the current residual =')
         end if
! 
!        %---------------------------------------------------------%
!        | Check for exact zero. Equivalent to determing whether a |
!        | j-step Arnoldi factorization is present.                |
!        %---------------------------------------------------------%
!
         if (rnorm .gt. zero) go to 40
!
!           %---------------------------------------------------%
!           | Invariant subspace found, generate a new starting |
!           | vector which is orthogonal to the current Arnoldi |
!           | basis and continue the iteration.                 |
!           %---------------------------------------------------%
!
            if (msglvl .gt. 0) then
               j_array(1) = j
               call ivout (logfil, 1, J_ARRAY, ndigit,                         &
                           '_saitr: ****** restart at step ******')
            end if
! 
!           %---------------------------------------------%
!           | ITRY is the loop variable that controls the |
!           | maximum amount of times that a restart is   |
!           | attempted. NRSTRT is used by stat.h         |
!           %---------------------------------------------%
!
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
!
!           %--------------------------------------%
!           | If in reverse communication mode and |
!           | RSTART = .true. flow returns here.   |
!           %--------------------------------------%
!
            
            call dgetv0 (ido, bmat, itry, .false., n, j, v, ldv,               &
                         resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
!
!              %------------------------------------------------%
!              | Give up after several restart attempts.        |
!              | Set INFO to the size of the invariant subspace |
!              | which spans OP and exit.                       |
!              %------------------------------------------------%
!
               info = j - 1
               call cpu_time (t1)
               tsaitr = tsaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
! 
   40    continue
!
!        %---------------------------------------------------------%
!        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
!        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
!        | when reciprocating a small RNORM, test against lower    |
!        | machine bound.                                          |
!        %---------------------------------------------------------%
!
         call dcopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. safmin) then
             temp1 = one / rnorm
             call dscal (n, temp1, v(1,j), 1)
             call dscal (n, temp1, workd(ipj), 1)
         else
!
!            %-----------------------------------------%
!            | To scale both v_{j} and p_{j} carefully |
!            | use LAPACK routine SLASCL               |
!            %-----------------------------------------%
!
             call dlascl ('General', i, i, rnorm, one, n, 1,                   &
                          v(1,j), n, infol)
             call dlascl ('General', i, i, rnorm, one, n, 1,                   &
                          workd(ipj), n, infol)
         end if
! 
!        %------------------------------------------------------%
!        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
!        | Note that this is not quite yet r_{j}. See STEP 4    |
!        %------------------------------------------------------%
!
         step3 = .true.
         nopx  = nopx + 1
         call cpu_time (t2)
         call dcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
! 
!        %-----------------------------------%
!        | Exit in order to compute OP*v_{j} |
!        %-----------------------------------%
! 
         go to 9000
   50    continue
! 
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
!        %-----------------------------------%
!
         call cpu_time (t3)
         tmvopx = tmvopx + (t3 - t2)
! 
         step3 = .false.
!
!        %------------------------------------------%
!        | Put another copy of OP*v_{j} into RESID. |
!        %------------------------------------------%
!
         call dcopy (n, workd(irj), 1, resid, 1)
! 
!        %-------------------------------------------%
!        | STEP 4:  Finish extending the symmetric   |
!        |          Arnoldi to length j. If MODE = 2 |
!        |          then B*OP = B*inv(B)*A = A and   |
!        |          we don't need to compute B*OP.   |
!        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
!        | assumed to have A*v_{j}.                  |
!        %-------------------------------------------%
!
         if (mode .eq. 2) go to 65
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
! 
!           %-------------------------------------%
!           | Exit in order to compute B*OP*v_{j} |
!           %-------------------------------------%
! 
            go to 9000
         else if (bmat .eq. 'I') then
              call dcopy(n, resid, 1 , workd(ipj), 1)
         end if
   60    continue
! 
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
!        %-----------------------------------%
!
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if 
!
         step4 = .false.
!
!        %-------------------------------------%
!        | The following is needed for STEP 5. |
!        | Compute the B-norm of OP*v_{j}.     |
!        %-------------------------------------%
!
   65    continue
         if (mode .eq. 2) then
!
!           %----------------------------------%
!           | Note that the B-norm of OP*v_{j} |
!           | is the inv(B)-norm of A*v_{j}.   |
!           %----------------------------------%
!
            wnorm = ddot (n, resid, 1, workd(ivj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'G') then         
            wnorm = ddot (n, resid, 1, workd(ipj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = dnrm2(n, resid, 1)
         end if
!
!        %-----------------------------------------%
!        | Compute the j-th residual corresponding |
!        | to the j step factorization.            |
!        | Use Classical Gram Schmidt and compute: |
!        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
!        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
!        %-----------------------------------------%
!
!
!        %------------------------------------------%
!        | Compute the j Fourier coefficients w_{j} |
!        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
!        %------------------------------------------%
!
         if (mode .ne. 2 ) then
            call dgemv('T', n, j, one, v, ldv, workd(ipj), 1, zero,            &
                        workd(irj), 1)
         else if (mode .eq. 2) then
            call dgemv('T', n, j, one, v, ldv, workd(ivj), 1, zero,            &
                        workd(irj), 1)
         end if
!
!        %--------------------------------------%
!        | Orthgonalize r_{j} against V_{j}.    |
!        | RESID contains OP*v_{j}. See STEP 3. | 
!        %--------------------------------------%
!
         call dgemv('N', n, j, -one, v, ldv, workd(irj), 1, one,               &
                     resid, 1)
!
!        %--------------------------------------%
!        | Extend H to have j rows and columns. |
!        %--------------------------------------%
!
         h(j,2) = workd(irj + j - 1)
         if (j .eq. 1  .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         end if
         call cpu_time (t4)
! 
         orth1 = .true.
         iter  = 0
! 
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
! 
!           %----------------------------------%
!           | Exit in order to compute B*r_{j} |
!           %----------------------------------%
! 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
! 
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH1 = .true. |
!        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
!        %---------------------------------------------------%
!
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
! 
         orth1 = .false.
!
!        %------------------------------%
!        | Compute the B-norm of r_{j}. |
!        %------------------------------%
!
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
!
!        %-----------------------------------------------------------%
!        | STEP 5: Re-orthogonalization / Iterative refinement phase |
!        | Maximum NITER_ITREF tries.                                |
!        |                                                           |
!        |          s      = V_{j}^T * B * r_{j}                     |
!        |          r_{j}  = r_{j} - V_{j}*s                         |
!        |          alphaj = alphaj + s_{j}                          |
!        |                                                           |
!        | The stopping criteria used for iterative refinement is    |
!        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
!        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
!        | Determine if we need to correct the residual. The goal is |
!        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
!        %-----------------------------------------------------------%
!
         if (rnorm .gt. 0.717*wnorm) go to 100
         nrorth = nrorth + 1
! 
!        %---------------------------------------------------%
!        | Enter the Iterative refinement phase. If further  |
!        | refinement is necessary, loop back here. The loop |
!        | variable is ITER. Perform a step of Classical     |
!        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
!        %---------------------------------------------------%
!
   80    continue
!
         if (msglvl .gt. 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call dvout (logfil, 2, xtemp, ndigit,                              &
                 '_saitr: re-orthonalization ; wnorm and rnorm are')
         end if
!
!        %----------------------------------------------------%
!        | Compute V_{j}^T * B * r_{j}.                       |
!        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1,                    &
                     zero, workd(irj), 1)
!
!        %----------------------------------------------%
!        | Compute the correction to the residual:      |
!        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
!        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
!        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
!        | H(j,j) is updated.                           |
!        %----------------------------------------------%
!
         call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1,                   &
                     one, resid, 1)
!
         if (j .eq. 1  .or.  rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workd(irj + j - 1)
! 
         orth2 = .true.
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
! 
!           %-----------------------------------%
!           | Exit in order to compute B*r_{j}. |
!           | r_{j} is the corrected residual.  |
!           %-----------------------------------%
! 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH2 = .true. |
!        %---------------------------------------------------%
!
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
!        %-----------------------------------------------------%
!        | Compute the B-norm of the corrected residual r_{j}. |
!        %-----------------------------------------------------%
! 
         if (bmat .eq. 'G') then         
             rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = dnrm2(n, resid, 1)
         end if
!
         if (msglvl .gt. 0 .and. iter .gt. 0) then
            j_array(1) = j
            call ivout (logfil, 1, J_ARRAY, ndigit,                            &
                 '_saitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call dvout (logfil, 2, xtemp, ndigit,                          &
                 '_saitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
! 
!        %-----------------------------------------%
!        | Determine if we need to perform another |
!        | step of re-orthogonalization.           |
!        %-----------------------------------------%
!
         if (rnorm1 .gt. 0.717*rnorm) then
!
!           %--------------------------------%
!           | No need for further refinement |
!           %--------------------------------%
!
            rnorm = rnorm1
! 
         else
!
!           %-------------------------------------------%
!           | Another step of iterative refinement step |
!           | is required. NITREF is used by stat.h     |
!           %-------------------------------------------%
!
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
!
!           %-------------------------------------------------%
!           | Otherwise RESID is numerically in the span of V |
!           %-------------------------------------------------%
!
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
! 
!        %----------------------------------------------%
!        | Branch here directly if iterative refinement |
!        | wasn't necessary or after at most NITER_REF  |
!        | steps of iterative refinement.               |
!        %----------------------------------------------%
!
  100    continue
! 
         rstart = .false.
         orth2  = .false.
! 
         call cpu_time (t5)
         titref = titref + (t5 - t4)
! 
!        %----------------------------------------------------------%
!        | Make sure the last off-diagonal element is non negative  |
!        | If not perform a similarity transformation on H(1:j,1:j) |
!        | and scale v(:,j) by -1.                                  |
!        %----------------------------------------------------------%
!
         if (h(j,1) .lt. zero) then
            h(j,1) = -h(j,1)
            if ( j .lt. k+np) then 
               call dscal(n, -one, v(1,j+1), 1)
            else
               call dscal(n, -one, resid, 1)
            end if
         end if
! 
!        %------------------------------------%
!        | STEP 6: Update  j = j+1;  Continue |
!        %------------------------------------%
!
         j = j + 1
         if (j .gt. k+np) then
            call cpu_time (t1)
            tsaitr = tsaitr + (t1 - t0)
            ido = 99
!
            if (msglvl .gt. 1) then
               call dvout (logfil, k+np, h(1,2), ndigit,                       &
               '_saitr: main diagonal of matrix H of step K+NP.')
               if (k+np .gt. 1) then
               call dvout (logfil, k+np-1, h(2,1), ndigit,                     &
               '_saitr: sub diagonal of matrix H of step K+NP.')
               end if
            end if
!
            go to 9000
         end if
!
!        %--------------------------------------------------------%
!        | Loop back to extend the factorization by another step. |
!        %--------------------------------------------------------%
!
      go to 1000
! 
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 9000 continue
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!     %---------------%
!     | End of dsaitr |
!     %---------------%
!
      end subroutine dsaitr

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsapps
!
!\Description:
!  Given the Arnoldi factorization
!
!     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!
!  apply NP shifts implicitly resulting in
!
!     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!
!  where Q is an orthogonal matrix of order KEV+NP. Q is the product of 
!  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi 
!  factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call dsapps
!     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Problem size, i.e. dimension of matrix A.
!
!  KEV     Integer.  (INPUT)
!          INPUT: KEV+NP is the size of the input matrix H.
!          OUTPUT: KEV is the size of the updated matrix HNEW.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFT   REAL(DOUBLE) array of length NP.  (INPUT)
!          The shifts to be applied.
!
!  V       REAL(DOUBLE) N by (KEV+NP) array.  (INPUT/OUTPUT)
!          INPUT: V contains the current KEV+NP Arnoldi vectors.
!          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
!          are in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       REAL(DOUBLE) (KEV+NP) by 2 array.  (INPUT/OUTPUT)
!          INPUT: H contains the symmetric tridiagonal matrix of the
!          Arnoldi factorization with the subdiagonal in the 1st column
!          starting at H(2,1) and the main diagonal in the 2nd column.
!          OUTPUT: H contains the updated tridiagonal matrix in the 
!          KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   REAL(DOUBLE) array of length (N).  (INPUT/OUTPUT)
!          INPUT: RESID contains the the residual vector r_{k+p}.
!          OUTPUT: RESID is the updated residual vector rnew_{k}.
!
!  Q       REAL(DOUBLE) KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations during the bulge
!          chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKD   REAL(DOUBLE) work array of length 2*N.  (WORKSPACE)
!          Distributed array used in the application of the accumulated
!          orthogonal matrix Q.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers. 
!     second  ARPACK utility routine for timing.
!     dvout   ARPACK utility routine that prints vectors.
!     dlamch  LAPACK routine that determines machine constants.
!     dlartg  LAPACK Givens rotation construction routine.
!     dlacpy  LAPACK matrix copy routine.
!     dlaset  LAPACK matrix initialization routine.
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     daxpy   Level 1 BLAS that computes a vector triad.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     dscal   Level 1 BLAS that scales a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     12/16/93: Version ' 2.4'
!
!\SCCS Information: @(#) 
! FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2
!
!\Remarks
!  1. In this version, each shift is applied to all the subblocks of
!     the tridiagonal matrix H and not just to the submatrix that it 
!     comes from. This routine assumes that the subdiagonal elements 
!     of H that are stored in h(1:kev+np,1) are nonegative upon input
!     and enforce this condition upon output. This version incorporates
!     deflation. See code for documentation.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsapps                                                        &
         ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAPPS'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    kev, ldh, ldq, ldv, n, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      REAL(DOUBLE)                                                             &
                 h(ldh,2), q(ldq,kev+np), resid(n), shift(np),                 &
                 v(ldv,kev+np), workd(2*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, iend, istart, itop, j, jj, kplusp, msglvl
      INTEGER    I_ARRAY(1), JJ_ARRAY(1)
      logical    first
      REAL(DOUBLE)                                                             &
                 a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      save       epsmch, first
!
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   daxpy, dcopy, dscal, dlacpy, dlartg, dlaset, dvout, 
!:!  &           ivout, second, dgemv
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!:!   REAL(DOUBLE)
!:!  &           dlamch
!:!   external   dlamch
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs
!
!     %----------------%
!     | Data statments |
!     %----------------%
!
      data       first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      if (first) then
         epsmch = dlamch('Epsilon-Machine')
         first = .false.
      end if
      itop = 1
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = msapps
! 
      kplusp = kev + np 
! 
!     %----------------------------------------------%
!     | Initialize Q to the identity matrix of order |
!     | kplusp used to accumulate the rotations.     |
!     %----------------------------------------------%
!
      call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)
!
!     %----------------------------------------------%
!     | Quick return if there are no shifts to apply |
!     %----------------------------------------------%
!
      if (np .eq. 0) go to 9000
! 
!     %----------------------------------------------------------%
!     | Apply the np shifts implicitly. Apply each shift to the  |
!     | whole matrix and not just to the submatrix from which it |
!     | comes.                                                   |
!     %----------------------------------------------------------%
!
      do 90 jj = 1, np
! 
         istart = itop
!
!        %----------------------------------------------------------%
!        | Check for splitting and deflation. Currently we consider |
!        | an off-diagonal element h(i+1,1) negligible if           |
!        |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
!        | for i=1:KEV+NP-1.                                        |
!        | If above condition tests true then we set h(i+1,1) = 0.  |
!        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
!        %----------------------------------------------------------%
!
   20    continue
!
!        %------------------------------------------------%
!        | The following loop exits early if we encounter |
!        | a negligible off diagonal element.             |
!        %------------------------------------------------%
!
         do 30 i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) .le. epsmch*big) then
               if (msglvl .gt. 0) then
                  i_array(1)  = i
                  jj_array(1) = jj
                  call ivout (logfil, 1, I_ARRAY, ndigit,                      &
                       '_sapps: deflation at row/column no.')
                  call ivout (logfil, 1, JJ_ARRAY, ndigit,                     &
                       '_sapps: occured before shift number.')
                  call dvout (logfil, 1, h(i+1,1), ndigit,                     &
                       '_sapps: the corresponding off diagonal element')
               end if
               h(i+1,1) = zero
               iend = i
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
!
         if (istart .lt. iend) then
! 
!           %--------------------------------------------------------%
!           | Construct the plane rotation G'(istart,istart+1,theta) |
!           | that attempts to drive h(istart+1,1) to zero.          |
!           %--------------------------------------------------------%
!
             f = h(istart,2) - shift(jj)
             g = h(istart+1,1)
             call dlartg (f, g, c, s, r)
! 
!            %-------------------------------------------------------%
!            | Apply rotation to the left and right of H;            |
!            | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
!            | This will create a "bulge".                           |
!            %-------------------------------------------------------%
!
             a1 = c*h(istart,2)   + s*h(istart+1,1)
             a2 = c*h(istart+1,1) + s*h(istart+1,2)
             a4 = c*h(istart+1,2) - s*h(istart+1,1)
             a3 = c*h(istart+1,1) - s*h(istart,2) 
             h(istart,2)   = c*a1 + s*a2
             h(istart+1,2) = c*a4 - s*a3
             h(istart+1,1) = c*a3 + s*a4
! 
!            %----------------------------------------------------%
!            | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!            %----------------------------------------------------%
!
             do 60 j = 1, min(istart+jj,kplusp)
                a1            =   c*q(j,istart) + s*q(j,istart+1)
                q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
                q(j,istart)   = a1
   60        continue
!
!
!            %----------------------------------------------%
!            | The following loop chases the bulge created. |
!            | Note that the previous rotation may also be  |
!            | done within the following loop. But it is    |
!            | kept separate to make the distinction among  |
!            | the bulge chasing sweeps and the first plane |
!            | rotation designed to drive h(istart+1,1) to  |
!            | zero.                                        |
!            %----------------------------------------------%
!
             do 70 i = istart+1, iend-1
! 
!               %----------------------------------------------%
!               | Construct the plane rotation G'(i,i+1,theta) |
!               | that zeros the i-th bulge that was created   |
!               | by G(i-1,i,theta). g represents the bulge.   |
!               %----------------------------------------------%
!
                f = h(i,1)
                g = s*h(i+1,1)
!
!               %----------------------------------%
!               | Final update with G(i-1,i,theta) |
!               %----------------------------------%
!
                h(i+1,1) = c*h(i+1,1)
                call dlartg (f, g, c, s, r)
!
!               %-------------------------------------------%
!               | The following ensures that h(1:iend-1,1), |
!               | the first iend-2 off diagonal of elements |
!               | H, remain non negative.                   |
!               %-------------------------------------------%
!
                if (r .lt. zero) then
                   r = -r
                   c = -c
                   s = -s
                end if
! 
!               %--------------------------------------------%
!               | Apply rotation to the left and right of H; |
!               | H <- G * H * G',  where G = G(i,i+1,theta) |
!               %--------------------------------------------%
!
                h(i,1) = r
! 
                a1 = c*h(i,2)   + s*h(i+1,1)
                a2 = c*h(i+1,1) + s*h(i+1,2)
                a3 = c*h(i+1,1) - s*h(i,2)
                a4 = c*h(i+1,2) - s*h(i+1,1)
! 
                h(i,2)   = c*a1 + s*a2
                h(i+1,2) = c*a4 - s*a3
                h(i+1,1) = c*a3 + s*a4
! 
!               %----------------------------------------------------%
!               | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!               %----------------------------------------------------%
!
                do 50 j = 1, min( i+jj, kplusp )
                   a1       =   c*q(j,i) + s*q(j,i+1)
                   q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                   q(j,i)   = a1
   50           continue
!
   70        continue
!
         end if
!
!        %--------------------------%
!        | Update the block pointer |
!        %--------------------------%
!
         istart = iend + 1
!
!        %------------------------------------------%
!        | Make sure that h(iend,1) is non-negative |
!        | If not then set h(iend,1) <-- -h(iend,1) |
!        | and negate the last column of Q.         |
!        | We have effectively carried out a        |
!        | similarity on transformation H           |
!        %------------------------------------------%
!
         if (h(iend,1) .lt. zero) then
             h(iend,1) = -h(iend,1)
             call dscal(kplusp, -one, q(1,iend), 1)
         end if
!
!        %--------------------------------------------------------%
!        | Apply the same shift to the next block if there is any |
!        %--------------------------------------------------------%
!
         if (iend .lt. kplusp) go to 20
!
!        %-----------------------------------------------------%
!        | Check if we can increase the the start of the block |
!        %-----------------------------------------------------%
!
         do 80 i = itop, kplusp-1
            if (h(i+1,1) .gt. zero) go to 90
            itop  = itop + 1
   80    continue
!
!        %-----------------------------------%
!        | Finished applying the jj-th shift |
!        %-----------------------------------%
!
   90 continue
!
!     %------------------------------------------%
!     | All shifts have been applied. Check for  |
!     | more possible deflation that might occur |
!     | after the last shift is applied.         |                               
!     %------------------------------------------%
!
      do 100 i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) .le. epsmch*big) then
            if (msglvl .gt. 0) then
               i_array(1) = i
               call ivout (logfil, 1, I_ARRAY, ndigit,                         &
                    '_sapps: deflation at row/column no.')
               call dvout (logfil, 1, h(i+1,1), ndigit,                        &
                    '_sapps: the corresponding off diagonal element')
            end if
            h(i+1,1) = zero
         end if
 100  continue
!
!     %-------------------------------------------------%
!     | Compute the (kev+1)-st column of (V*Q) and      |
!     | temporarily store the result in WORKD(N+1:2*N). |
!     | This is not necessary if h(kev+1,1) = 0.         |
!     %-------------------------------------------------%
!
      if ( h(kev+1,1) .gt. zero )                                              &
         call dgemv ('N', n, kplusp, one, v, ldv,                              &
                      q(1,kev+1), 1, zero, workd(n+1), 1)
! 
!     %-------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order    |
!     | taking advantage that Q is an upper triangular matrix |    
!     | with lower bandwidth np.                              |
!     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
!     %-------------------------------------------------------%
!
      do 130 i = 1, kev
         call dgemv ('N', n, kplusp-i+1, one, v, ldv,                          &
                     q(1,kev-i+1), 1, zero, workd, 1)
         call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  130 continue
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
      call dlacpy ('All', n, kev, v(1,np+1), ldv, v, ldv)
! 
!     %--------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the |
!     | appropriate place if h(kev+1,1) .ne. zero. |
!     %--------------------------------------------%
!
      if ( h(kev+1,1) .gt. zero )                                              &
           call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)
! 
!     %-------------------------------------%
!     | Update the residual vector:         |
!     |    r <- sigmak*r + betak*v(:,kev+1) |
!     | where                               |
!     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
!     |    betak = e_{kev+1}'*H*e_{kev}     |
!     %-------------------------------------%
!
      call dscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,1) .gt. zero)                                                &
         call daxpy (n, h(kev+1,1), v(1,kev+1), 1, resid, 1)
!
      if (msglvl .gt. 1) then
         call dvout (logfil, 1, q(kplusp,kev), ndigit,                         &
            '_sapps: sigmak of the updated residual vector')
         call dvout (logfil, 1, h(kev+1,1), ndigit,                            &
            '_sapps: betak of the updated residual vector')
         call dvout (logfil, kev, h(1,2), ndigit,                              &
            '_sapps: updated main diagonal of H for next iteration')
         if (kev .gt. 1) then
         call dvout (logfil, kev-1, h(2,1), ndigit,                            &
            '_sapps: updated sub diagonal of H for next iteration')
         end if
      end if
!
      call cpu_time (t1)
      tsapps = tsapps + (t1 - t0)
! 
 9000 continue 
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!     %---------------%
!     | End of dsapps |
!     %---------------%
!
      end subroutine dsapps

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsconv
!
!\Description: 
!  Convergence testing for the symmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call dsconv
!     ( N, RITZ, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZ    REAL(DOUBLE) array of length N.  (INPUT)
!          The Ritz values to be checked for convergence.
!
!  BOUNDS  REAL(DOUBLE) array of length N.  (INPUT)
!          Ritz estimates associated with the Ritz values in RITZ.
!
!  TOL     REAL(DOUBLE) scalar.  (INPUT)
!          Desired relative accuracy for a Ritz value to be considered
!          "converged".
!
!  NCONV   Integer scalar.  (OUTPUT)
!          Number of "converged" Ritz values.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines called:
!     second  ARPACK utility routine for timing.
!     dlamch  LAPACK routine that determines machine constants. 
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\SCCS Information: @(#) 
! FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.4, this routine no longer uses the
!        Parlett strategy using the gap conditions. 
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dsconv (n, ritz, bounds, tol, nconv)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSCONV'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    n, nconv
      REAL(DOUBLE)                                                             &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      REAL(DOUBLE)                                                             &
                 ritz(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i
      REAL(DOUBLE)                                                             &
                 temp, eps23
!
!     %-------------------%
!     | External routines |
!     %-------------------%
!
!:!   REAL(DOUBLE)
!:!  &           dlamch
!:!   external   dlamch

!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

      call cpu_time (t0)
!
      eps23 = dlamch('Epsilon-Machine') 
      eps23 = eps23**(2.0D+0 / 3.0D+0)
!
      nconv  = 0
      do 10 i = 1, n
!
!        %-----------------------------------------------------%
!        | The i-th Ritz value is considered "converged"       |
!        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
!        %-----------------------------------------------------%
!
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) then 
            nconv = nconv + 1
         end if
!
   10 continue
! 
      call cpu_time (t1)
      tsconv = tsconv + (t1 - t0)
! 
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!     %---------------%
!     | End of dsconv |
!     %---------------%
!
      end subroutine dsconv

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dseigt
!
!\Description: 
!  Compute the eigenvalues of the current symmetric tridiagonal matrix
!  and the corresponding error bounds given the current residual norm.
!
!\Usage:
!  call dseigt
!     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
!
!\Arguments
!  RNORM   REAL(DOUBLE) scalar.  (INPUT)
!          RNORM contains the residual norm corresponding to the current
!          symmetric tridiagonal matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the symmetric tridiagonal matrix H.
!
!  H       REAL(DOUBLE) N by 2 array.  (INPUT)
!          H contains the symmetric tridiagonal matrix with the 
!          subdiagonal in the first column starting at H(2,1) and the 
!          main diagonal in second column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling 
!          program.
!
!  EIG     REAL(DOUBLE) array of length N.  (OUTPUT)
!          On output, EIG contains the N eigenvalues of H possibly 
!          unsorted.  The BOUNDS arrays are returned in the
!          same sorted order as EIG.
!
!  BOUNDS  REAL(DOUBLE) array of length N.  (OUTPUT)
!          On output, BOUNDS contains the error estimates corresponding
!          to the eigenvalues EIG.  This is equal to RNORM times the
!          last components of the eigenvectors corresponding to the
!          eigenvalues in EIG.
!
!  WORKL   REAL(DOUBLE) work array of length 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from dstqrb.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     dstqrb  ARPACK routine that computes the eigenvalues and the
!             last components of the eigenvectors of a symmetric
!             and tridiagonal matrix.
!     second  ARPACK utility routine for timing.
!     dvout   ARPACK utility routine that prints vectors.
!     dcopy   Level 1 BLAS that copies one vector to another.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!
!\SCCS Information: @(#) 
! FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dseigt                                                        &
         ( rnorm, n, h, ldh, eig, bounds, workl, ierr )


      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSEIGT'
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    ierr, ldh, n
      REAL(DOUBLE)                                                             &
                 rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      REAL(DOUBLE)                                                             &
                 eig(n), bounds(n), h(ldh,2), workl(3*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      REAL(DOUBLE)                                                             &
                 zero
      parameter (zero = 0.0D+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, k, msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
!:!   external   dcopy, dstqrb, dvout, second
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------% 
!
      call cpu_time (t0)
      msglvl = mseigt
!
      if (msglvl .gt. 0) then
         call dvout (logfil, n, h(1,2), ndigit,                                &
                    '_seigt: main diagonal of matrix H')
         if (n .gt. 1) then
         call dvout (logfil, n-1, h(2,1), ndigit,                              &
                    '_seigt: sub diagonal of matrix H')
         end if
      end if
!
      call dcopy  (n, h(1,2), 1, eig, 1)
      call dcopy  (n-1, h(2,1), 1, workl, 1)
      call dstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr .ne. 0) go to 9000
      if (msglvl .gt. 1) then
         call dvout (logfil, n, bounds, ndigit,                                &
                    '_seigt: last row of the eigenvector matrix for H')
      end if
!
!     %-----------------------------------------------%
!     | Finally determine the error bounds associated |
!     | with the n Ritz values of H.                  |
!     %-----------------------------------------------%
!
      do 30 k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
   30 continue
! 
      call cpu_time (t1)
      tseigt = tseigt + (t1 - t0)
!
 9000 continue
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+3) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!     %---------------%
!     | End of dseigt |
!     %---------------%
!
      end subroutine dseigt

! ################################################################################################################################## 

!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dstqrb
!
!\Description:
!  Computes all eigenvalues and the last component of the eigenvectors
!  of a symmetric tridiagonal matrix using the implicit QL or QR method.
!
!  This is mostly a modification of the LAPACK routine dsteqr.
!  See Remarks.
!
!\Usage:
!  call dstqrb
!     ( N, D, E, Z, WORK, INFO )
!
!\Arguments
!  N       Integer.  (INPUT)
!          The number of rows and columns in the matrix.  N >= 0.
!
!  D       REAL(DOUBLE) array, dimension (N).  (INPUT/OUTPUT)
!          On entry, D contains the diagonal elements of the
!          tridiagonal matrix.
!          On exit, D contains the eigenvalues, in ascending order.
!          If an error exit is made, the eigenvalues are correct
!          for indices 1,2,...,INFO-1, but they are unordered and
!          may not be the smallest eigenvalues of the matrix.
!
!  E       REAL(DOUBLE) array, dimension (N-1).  (INPUT/OUTPUT)
!          On entry, E contains the subdiagonal elements of the
!          tridiagonal matrix in positions 1 through N-1.
!          On exit, E has been destroyed.
!
!  Z       REAL(DOUBLE) array, dimension (N).  (OUTPUT)
!          On exit, Z contains the last row of the orthonormal 
!          eigenvector matrix of the symmetric tridiagonal matrix.  
!          If an error exit is made, Z contains the last row of the
!          eigenvector matrix associated with the stored eigenvalues.
!
!  WORK    REAL(DOUBLE) array, dimension (max(1,2*N-2)).  (WORKSPACE)
!          Workspace used in accumulating the transformation for 
!          computing the last components of the eigenvectors.
!
!  INFO    Integer.  (OUTPUT)
!          = 0:  normal return.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = +i, the i-th eigenvalue has not converged
!                              after a total of  30*N  iterations.
!
!\Remarks
!  1. None.
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     daxpy   Level 1 BLAS that computes a vector triad.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     dswap   Level 1 BLAS that swaps the contents of two vectors.
!     lsame   LAPACK character comparison routine.
!     dlae2   LAPACK routine that computes the eigenvalues of a 2-by-2 
!             symmetric matrix.
!     dlaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric 
!             matrix.
!     dlamch  LAPACK routine that determines machine constants.
!     dlanst  LAPACK routine that computes the norm of a matrix.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dlartg  LAPACK Givens rotation construction routine.
!     dlascl  LAPACK routine for careful scaling of a matrix.
!     dlaset  LAPACK matrix initialization routine.
!     dlasr   LAPACK routine that applies an orthogonal transformation to 
!             a matrix.
!     dlasrt  LAPACK sorting routine.
!     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
!             of a symmetric tridiagonal matrix.
!     xerbla  LAPACK error handler routine.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\SCCS Information: @(#) 
! FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.5, this routine is a modified version
!        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
!        only commeted out and new lines inserted.
!        All lines commented out have "c$$$" at the beginning.
!        Note that the LAPACK version 1.0 subroutine SSTEQR contained
!        bugs. 
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine dstqrb ( n, d, e, z, work, info )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSTQRB'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    info, n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      REAL(DOUBLE)                                                             &
                 d( n ), e( n-1 ), z( n ), work( 2*n-2 )
!
!     .. parameters ..
      REAL(DOUBLE)                                                             &
                         zero, one, two, three
      parameter          ( zero = 0.0D+0, one = 1.0D+0,                        &
                           two = 2.0D+0, three = 3.0D+0 )
      integer            maxit
      parameter          ( maxit = 30 )
!     ..
!     .. local scalars ..
      integer            i, icompz, ii, iscale, j, jtot, k, l, l1, lend,       &
                         lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1,         &
                         nm1, nmaxit
      REAL(DOUBLE)                                                             &
                         anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2,         &
                         s, safmax, safmin, ssfmax, ssfmin, tst
!     ..
!     .. external functions ..
!:!   logical            lsame
!:!   REAL(DOUBLE)
!:!  &                   dlamch, dlanst, dlapy2
!:!   external           lsame, dlamch, dlanst, dlapy2
!     ..
!     .. external subroutines ..
!:!   external           dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr,
!:!  &                   dlasrt, dswap, xerbla
!     ..
!     .. intrinsic functions ..
      intrinsic          abs, max, sign, sqrt
!     ..
!     .. executable statements ..
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+4) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! **********************************************************************************************************************************

!     test the input parameters.
!
      info = 0
!
!$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
!$$$         ICOMPZ = 0
!$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
!$$$         ICOMPZ = 1
!$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
!$$$         ICOMPZ = 2
!$$$      ELSE
!$$$         ICOMPZ = -1
!$$$      END IF
!$$$      IF( ICOMPZ.LT.0 ) THEN
!$$$         INFO = -1
!$$$      ELSE IF( N.LT.0 ) THEN
!$$$         INFO = -2
!$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
!$$$     $         N ) ) ) THEN
!$$$         INFO = -6
!$$$      END IF
!$$$      IF( INFO.NE.0 ) THEN
!$$$         CALL XERBLA( 'SSTEQR', -INFO )
!$$$         RETURN
!$$$      END IF
!
!    *** New starting with version 2.5 ***
!
      icompz = 2
!    *************************************
!
!     quick return if possible
!
      if( n.eq.0 )                                                             &
         return
!
      if( n.eq.1 ) then
         if( icompz.eq.2 )  z( 1 ) = one
         return
      end if
!
!     determine the unit roundoff and over/underflow thresholds.
!
      eps = dlamch( 'e' )
      eps2 = eps**2
      safmin = dlamch( 's' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2
!
!     compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
!$$      if( icompz.eq.2 )
!$$$     $   call dlaset( 'full', n, n, zero, one, z, ldz )
!
!     *** New starting with version 2.5 ***
!
      if ( icompz .eq. 2 ) then
         do 5 j = 1, n-1
            z(j) = zero
  5      continue
         z( n ) = one
      end if
!     *************************************
!
      nmaxit = n*maxit
      jtot = 0
!
!     determine where the matrix splits and choose ql or qr iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      l1 = 1
      nm1 = n - 1
!
   10 continue
      if( l1.gt.n )                                                            &
         go to 160
      if( l1.gt.1 )                                                            &
         e( l1-1 ) = zero
      if( l1.le.nm1 ) then
         do 20 m = l1, nm1
            tst = abs( e( m ) )
            if( tst.eq.zero )                                                  &
               go to 30
            if( tst.le.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+                &
                1 ) ) ) )*eps ) then
               e( m ) = zero
               go to 30
            end if
   20    continue
      end if
      m = n
!
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if( lend.eq.l )                                                          &
         go to 10
!
!     scale submatrix in rows and columns l to lend
!
      anorm = dlanst( 'i', lend-l+1, d( l ), e( l ) )
      iscale = 0
      if( anorm.eq.zero )                                                      &
         go to 10
      if( anorm.gt.ssfmax ) then
         iscale = 1
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n,        &
                      info )
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n,          &
                      info )
      else if( anorm.lt.ssfmin ) then
         iscale = 2
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n,        &
                      info )
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n,          &
                      info )
      end if
!
!     choose between ql and qr iteration
!
      if( abs( d( lend ) ).lt.abs( d( l ) ) ) then
         lend = lsv
         l = lendsv
      end if
!
      if( lend.gt.l ) then
!
!        ql iteration
!
!        look for small subdiagonal element.
!
   40    continue
         if( l.ne.lend ) then
            lendm1 = lend - 1
            do 50 m = l, lendm1
               tst = abs( e( m ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+              &
                   safmin )go to 60
   50       continue
         end if
!
         m = lend
!
   60    continue
         if( m.lt.lend )                                                       &
            e( m ) = zero
         p = d( l )
         if( m.eq.l )                                                          &
            go to 80
!
!        if remaining matrix is 2-by-2, use dlae2 or dlaev2
!        to compute its eigensystem.
!
         if( m.eq.l+1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
               work( l ) = c
               work( n-1+l ) = s
!$$$               call dlasr( 'r', 'v', 'b', n, 2, work( l ),
!$$$     $                     work( n-1+l ), z( 1, l ), ldz )
!
!              *** New starting with version 2.5 ***
!
               tst      = z(l+1)
               z(l+1) = c*tst - s*z(l)
               z(l)   = s*tst + c*z(l)
!              *************************************
            else
               call dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
            end if
            d( l ) = rt1
            d( l+1 ) = rt2
            e( l ) = zero
            l = l + 2
            if( l.le.lend )                                                    &
               go to 40
            go to 140
         end if
!
         if( jtot.eq.nmaxit )                                                  &
            go to 140
         jtot = jtot + 1
!
!        form shift.
!
         g = ( d( l+1 )-p ) / ( two*e( l ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
!
         s = one
         c = one
         p = zero
!
!        inner loop
!
         mm1 = m - 1
         do 70 i = mm1, l, -1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m-1 )                                                     &
               e( i+1 ) = r
            g = d( i+1 ) - p
            r = ( d( i )-g )*s + two*c*b
            p = s*r
            d( i+1 ) = g + p
            g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = -s
            end if
!
   70    continue
!
!        if eigenvectors are desired, then apply saved rotations.
!
         if( icompz.gt.0 ) then
            mm = m - l + 1
!$$$            call dlasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ),
!$$$     $                  z( 1, l ), ldz )
!
!             *** New starting with version 2.5 ***
!
              call dlasr( 'r', 'v', 'b', 1, mm, work( l ),                     &
                          work( n-1+l ), z( l ), 1 )
!             *************************************                             
         end if
!
         d( l ) = d( l ) - p
         e( l ) = g
         go to 40
!
!        eigenvalue found.
!
   80    continue
         d( l ) = p
!
         l = l + 1
         if( l.le.lend )                                                       &
            go to 40
         go to 140
!
      else
!
!        qr iteration
!
!        look for small superdiagonal element.
!
   90    continue
         if( l.ne.lend ) then
            lendp1 = lend + 1
            do 100 m = l, lendp1, -1
               tst = abs( e( m-1 ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+              &
                   safmin )go to 110
  100       continue
         end if
!
         m = lend
!
  110    continue
         if( m.gt.lend )                                                       &
            e( m-1 ) = zero
         p = d( l )
         if( m.eq.l )                                                          &
            go to 130
!
!        if remaining matrix is 2-by-2, use dlae2 or dlaev2
!        to compute its eigensystem.
!
         if( m.eq.l-1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
!$$$               work( m ) = c
!$$$               work( n-1+m ) = s
!$$$               call dlasr( 'r', 'v', 'f', n, 2, work( m ),
!$$$     $                     work( n-1+m ), z( 1, l-1 ), ldz )
!
!               *** New starting with version 2.5 ***
!
                tst      = z(l)
                z(l)   = c*tst - s*z(l-1)
                z(l-1) = s*tst + c*z(l-1)
!               ************************************* 
            else
               call dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
            end if
            d( l-1 ) = rt1
            d( l ) = rt2
            e( l-1 ) = zero
            l = l - 2
            if( l.ge.lend )                                                    &
               go to 90
            go to 140
         end if
!
         if( jtot.eq.nmaxit )                                                  &
            go to 140
         jtot = jtot + 1
!
!        form shift.
!
         g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )
!
         s = one
         c = one
         p = zero
!
!        inner loop
!
         lm1 = l - 1
         do 120 i = m, lm1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m )                                                       &
               e( i-1 ) = r
            g = d( i ) - p
            r = ( d( i+1 )-g )*s + two*c*b
            p = s*r
            d( i ) = g + p
            g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = s
            end if
!
  120    continue
!
!        if eigenvectors are desired, then apply saved rotations.
!
         if( icompz.gt.0 ) then
            mm = l - m + 1
!$$$            call dlasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ),
!$$$     $                  z( 1, m ), ldz )
!
!           *** New starting with version 2.5 ***
!
            call dlasr( 'r', 'v', 'f', 1, mm, work( m ), work( n-1+m ),        &
                        z( m ), 1 )
!           *************************************                             
         end if
!
         d( l ) = d( l ) - p
         e( lm1 ) = g
         go to 90
!
!        eigenvalue found.
!
  130    continue
         d( l ) = p
!
         l = l - 1
         if( l.ge.lend )                                                       &
            go to 90
         go to 140
!
      end if
!
!     undo scaling if necessary
!
  140 continue
      if( iscale.eq.1 ) then
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1,               &
                      d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ),       &
                      n, info )
      else if( iscale.eq.2 ) then
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1,               &
                      d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ),       &
                      n, info )
      end if
!
!     check for no convergence to an eigenvalue after a total
!     of n*maxit iterations.
!
      if( jtot.lt.nmaxit )                                                     &
         go to 10
      do 150 i = 1, n - 1
         if( e( i ).ne.zero )                                                  &
            info = info + 1
  150 continue
      go to 190
!
!     order eigenvalues and eigenvectors.
!
  160 continue
      if( icompz.eq.0 ) then
!
!        use quick sort
!
         call dlasrt( 'i', n, d, info )
!
      else
!
!        use selection sort to minimize swaps of eigenvectors
!
         do 180 ii = 2, n
            i = ii - 1
            k = i
            p = d( i )
            do 170 j = ii, n
               if( d( j ).lt.p ) then
                  k = j
                  p = d( j )
               end if
  170       continue
            if( k.ne.i ) then
               d( k ) = d( i )
               d( i ) = p
!$$$               call dswap( n, z( 1, i ), 1, z( 1, k ), 1 )
!           *** New starting with version 2.5 ***
!
               p    = z(k)
               z(k) = z(i)
               z(i) = p
!           *************************************
            end if
  180    continue
      end if
!
  190 continue
!
! **********************************************************************************************************************************
      IF (WRT_LOG >= SUBR_BEGEND+4) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN

! **********************************************************************************************************************************

!     %---------------%
!     | End of dstqrb |
!     %---------------%
!
      end subroutine dstqrb

      END MODULE ARPACK_LANCZOS_EIG
