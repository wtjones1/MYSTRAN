! ################################################################################################################################## 

      MODULE BANDIT_MODULE

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                   :  F04
      USE SCONTR, ONLY                   :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                   :  HOUR, MINUTE, SEC,                 &
                                            SFRAC, TSEC
      USE SUBR_BEGEND_LEVELS, ONLY       :  BANDIT_BEGEND

      INTEGER(LONG), PARAMETER, PRIVATE :: SUBR_BEGEND = BANDIT_BEGEND

! Notes:
! ------

! (1) The correspondence of array IPARAM and the $ directive entries described in Gordon's documentation is as follows:

!          IPARAM( 1) corresponds to $LOOP       (not documented: $LOOP is used for program testing - loop on several input decks)
!          IPARAM( 2) corresponds to $TABLE       
!          IPARAM( 3) corresponds to $METHOD     
!          IPARAM( 4) corresponds to $MPC        
!          IPARAM( 5) corresponds to $SEQUENCE   
!          IPARAM( 6) corresponds to $CRITERION  
!          IPARAM( 7) corresponds to $SPRING     
!          IPARAM( 8) corresponds to $NASTRAN    (not documented - $NASTRAN YES card causes message to be printed to output file )
!          IPARAM( 9) corresponds to $ELEMENTS   
!          IPARAM(10) corresponds to $PRINT      
!          IPARAM(11) corresponds to $FRONTAL

! (2) The statement IPARAM(I) = J has the following meaning for various J:

!          J =  1 means IPARAM(I) is the I-th directive (above) with value SEQGP
!          J =  2 means IPARAM(I) is the I-th directive (above) with value ALL
!          J =  3 means IPARAM(I) is the I-th directive (above) with value NO
!          J =  4 means IPARAM(I) is the I-th directive (above) with value YES
!          J =  5 means IPARAM(I) is the I-th directive (above) with value MIN
!          J =  6 means IPARAM(I) is the I-th directive (above) with value MAX
!          J =  7 means IPARAM(I) is the I-th directive (above) with value CM
!          J =  8 means IPARAM(I) is the I-th directive (above) with value GPS
!          J =  9 means IPARAM(I) is the I-th directive (above) with value BOTH
!          J = 10 means IPARAM(I) is the I-th directive (above) with value BAND
!          J = 11 means IPARAM(I) is the I-th directive (above) with value PROFILE
!          J = 12 means IPARAM(I) is the I-th directive (above) with value RMS
!          J = 13 means IPARAM(I) is the I-th directive (above) with value WAVEFRONT

! For the most part, all changes to Gordon's code (except cosmetic ones) are delineated as shown in the example below where Program
! BANDIT was changed to a subr and where file unit in1 from module iount1 is added for use.

! Example of how mods are delineated:
! B////////////////////////////////////////////////////////////////////B
! E////////////////////////////////////////////////////////////////////E


! Mods to Gordon's original bandit.f file (for use as a set of subrs called by MYSTRAN):
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! Aug 03 - Oct 03 Mods:
! *********************

!    General changes:
!    ----------------
!    ( 1) Add ! ##############... comment lines as separators between each subroutine
!    ( 2) Modify subroutine END statements to append subroutine name to the END statement

!    Main:
!    -----
!    ( 1) Change Program BANDIT to subroutine bandit and add statement: "use iount1, only:  in1" 
!    ( 2) Add arg IER to subr BANDIT so that MYSTRAN can check for error and give message.
!    ( 3) Do not open IOU5, instead set iou5 = in1 (MYSTRAN input file)
!    ( 4) Change the statement: 135 IOU5=5, to 135 iou5 = in1 to reflect that IOU5 is the MYSTRAN input file.
!    ( 5) Change STOP 13 and STOP 5 to writing a message and returning to MYSTRAN
!    ( 6) Change final STOP to RETURN and change END to end subroutine bandit
!    ( 7) Change open status from 'UNKNOWN' to 'REPLACE' for units IOU6,7,8,9,10,14,16,17. This is done to make sure that MYSTRAN
!         does not attempt to use an old file.

!    BLOCK DATA:
!    -----------
!    ( 1) Remove definition of unit IOU5 (since it is now MYSTRAN IN1
!    ( 2) for all other file numbers, add 100 (e.g. IOU6 goes from unit 6 to unit 106, etc.)

!    Subr DOLLAR:
!    ------------
!    ( 1) chge default for IPARAM(6), the variable for $CRITERION, from RMS to BAND

!    Subr ELTYPE:
!    ------------
!    ( 1) Comment out the code that sets fatal error and return if SEQGP cards are in the input deck. MYSTRAN checks
!         this (before Bandit runs) and ignores SEQGP cards if the user has requested Bandit resequencing.

!    Subr FINISH:
!    ------------
!    ( 1) Add: use iount1, only: seq
!    ( 2) Add several lines regarding MYSTRAN file SEQ (used for writing SEQGP card images as well on Bandit file unit IOU6) 

!    Subr SEQGP:
!    -----------
!    ( 1) Add: use iount1, only: seq
!    ( 2) Add WRITE(IOU7,*) just before ENDFILE IOU7 to get 1 blank line written to IOU7 (in case there are no SEQGP card images
!         written to IOU7)
!    ( 3) Add several lines regarding MYSTRAN file SEQ so that the SEQGP card images are written to that file as well as Bandit IOU6 


! 11/30/03 Mods:
! **************

!     Subr TIMER:
!     -----------
!     ( 1) Comment out 2 lines in subroutine TIMER (variables tarray, iwall, time not used and Layhey complains about type specifier)
!     ( 2) Change "call second(t)" to "call cpu_time(t)". Calling second was causing the -NaN problem when BANDIT was imbedded in 
!          MYSTRAN since there is no intrinsic function in Layhey called  second. Not sure what is being used since there is no
!          intrinsic  function in Layhey documentation called second. There is a subr called second in module ARPACK_UTIL that I
!          modified to call cpu_time and, if we put the statement: use arpack_util above then we get no problem with -NaN. Module
!          ARPACK_UTIL calls cpu_time to get t and this is a valid procedure in Layhey. Thus, switch to call cpu_time here.


! Mods to Gordon's second bandit.f file (after he added RBE's):
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! 01/08/04 Mods:
! **************

!     General changes:
!     ----------------
!     ( 1) Change this collection of subroutines to a module and change its name to BANDIT_MODULE. This required changing several
!          scalar variables to appear as if they are arrays to avoid Lahey giving fatal errs regarding shapes of arrays not matching

!     Subr bandit:
!     ------------
!     ( 1) Earlier changes regarding STOP 13 and STOP 5 no longer apply. Gordon has removed these stops and, apparently, is handling
!          these situations in subr FINISH

!     BLOCK DATA:
!     -----------
!     ( 1) Remove this code from the module (it can't be in a module) and create a new procedure called BLOCK DATA BANDIT_BLOCK_DATA

!     Subr ELTYPE:
!     ------------
!     ( 1) Add integer array NCON_array(1) to be equiv to scalar NCON and call subr READIT with NCON_array(1) instead of NCON 

!     Subr FINISH:
!     ------------
!     ( 1) Cosmetic changes (not delineated with ! B////////...) to make subr more readable

!     Subr GIBSTK:
!     ------------
!     ( 1) Remove SORT2 from INTEGER declar. It is an integer function that is now part of this module and can't be redefined here

!     Subr INSERT:
!     ------------
!     ( 1) Add integer array NCARD_array(1) to be equiv to scalar NCARD and call subr READIT with NCARD_array(1) instead of NCARD 

!     Subr REED:
!     ----------
!     ( 1) Add integer array EID_array(1) to be equiv to scalar EID and call subr READIT with EID_array(1) instead of EID 

!     Subr RIGID:
!     -----------
!     ( 1) Add integer array IG_array(1) to be equiv to scalar IG and call subr READIT with IG_array(1) instead of IG (2 locations)
!     ( 2) Call SCAT with integer array IG_array(1) instead of scalar IG

!     Subr TAXI:
!     ----------
!     ( 1) Add integer array NAXIC_array(1) to be equiv to scalar NAXIC and call subr READIT with NAXIC_array(1) instead of NAXIC 

! 01/19/04 Mods:
! **************

!     General changes:
!     ----------------
!     ( 1) Add type declarations for all variables not alre3ady typed so that IMPLICIT NONE switch can be used in the compiler
!          directive.

! 02/01/04 Mods:
! **************

!     Subr BANDIT:
!     ------------
!     ( 1) Add MYSTRAN_NGRID to arg list (INTENT(IN)) in SUBROUTINE BANDIT and change NGRID = KOR/30  to MAX(MYSTRAN_NGRID,KORE/30)
!     ( 2) Just before setting NGRID as stated above, comment out the line KDIM=250 and, just after the line NGRID = MYSTRAN_NGRID
!          set KDIM = MAX(250,NGRID/10)

! 02/21/04 Mods:
! **************

!     Subr RIGID:
!     ------------
!     ( 1) Change dimension of LG(2) to LK(*) per 2/19/04 email from Gordon Everstine

! 03/02/04 Mods:
! **************

!     Subr BANDIT:
!     ------------
!     ( 1) Change to KDIM = MAX(250,MYSTRAN_NGRID) . It apparently is not large enough when there is a large "Nodal degree limit"
!          (as when RBE's connect to many grids). Later in Bandit, KDIM seems to be set to as large as NGRID so use max of the
!          250 that Gordon had here and NGRID.

!     ( 2) Change the 02/01/04 mod to: change original NGRID=KORE/30 to NGRID = MAX(MYSTRAN_NGRID,10).
!          The 10 is a guess to cover small problems since NGRID = MYSTRAN_NGRID doesn't seem to work for small problems and since
!          Bandit uses NGRID/10 later.
!          I had problem with chassis_grav. Bandit wouldn't run because Nodal degree limit of 21 was exceeded. Apparently, since
!          NGRID = MAX(MYSTRAN_NGRID,KORE/30) used KORE/30 = 10,000,000/30 (since KORE/30 > MYSTRAN_NGRID) was using much of the
!          memory for NGRID and then there was not enough for the Nodal degree limit. The Nodal degree limit can be set in Bandit
!          with the $GRID N directive. This sets NGRID to N, so when I input $GRID 7298 (actual number of grids in chassis_grav)
!          Bandit ran successfully. 

! 08/01/06 Mods:
! **************

!     Subr BANDIT:
!     ------------

!     ( 1) Add var NEW_BW to be returned from subr SUMUP to be the new bandwidth found
!     ( 2) Change call to SUMUP to have NEW_BW as an arg
!     ( 3) Add var NEW_BW to subr SUMUP ans set NEW_BW = NBW (from COMMON/D/ which is printed out there as the new bandwidth)

! 10/05/06 Mods:
! **************

!     Subr BANDIT_MODULE:
!     ------------------

!     Add call to BANDIT_FILES (new subr to make sure no files are open in another program and all old files are deleted)

! 07/23/08 Mods:
! **************

!     (1) Add DEN to arg list in CALL SUMUP
!     (2) Add DEN to arg list in subr SUMUP (returned to calling program)
!     (3) Add DEN as return arg to BANDIT subr

! 01/26/09 Mods:
! **************

!     Mod to version 6.00c
!     Change FORMAT #10 in subr SEQGP to write SEQGP card images in large field format. Need this when MYSTRAN SEQGP calls
!     subr BD_SEQGP. This was causing an error when GRID numbers in large field format were used

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################
! B////////////////////////////////////////////////////////////////////B
      SUBROUTINE BANDIT ( MYSTRAN_NGRID, NEW_BW, DEN, IER )
      USE IOUNT1, ONLY :  WRT_LOG, IN1
! E////////////////////////////////////////////////////////////////////E
!
!     Gordon C. Everstine, Gaithersburg, MD, geversti@comcast.net
!     1/8/04
!
!     This is the PC version.  The only machine or compiler dependencies
!     are in Subroutine Timer related to getting the CPU time.  This
!     subroutine must be checked before compiling.
!
!     Initial Bandit version: December 1969
!
!  Recent revisions:
!     4/16/93  Read MPCAX cards.
!     11/8/93  Minor changes needed to suit MS-FORTRAN.
!     4/22/94  Write to file the components and grids (unit 15) if
!              CM is executed.
!     12/5/94  Move connection table output from unit 11 to 16.
!              Make unit 11 a scratch file.
!              Delete selected files at end of job.
!              Add elapsed time report for HP workstations.
!     1/13/95  Add log file, unit 17; increase default KDIM to 250.
!     5/22/96  Speed up element sort in Subroutine FRONT.
!              Change file names from TAPE* to bandit.f*.
!    11/18/96  Change author's mailing address.
!     2/25/98  Add timing call for PC (subroutine timer).
!     6/20/03  Change author's email address.
!    10/27/03  Enable component table with $print max (unit 15).
!    11/20/03  Add additional elastic elements.
!      1/8/04  Add additional rigid elements and some cosmetic changes.
!              Remove extra plus signs, a holdover from BCD/EBCDIC days.
!
!  This is a NASTRAN preprocessor which reads a NASTRAN deck and
!  resequences the grid points to reduce the matrix bandwidth, profile,
!  or wavefront.  Two algorithms are provided,
!  Gibbs-Poole-Stockmeyer (GPS) and Reverse Cuthill-McKee (CM).
!  In addition, the elements can be resequenced for frontal solvers
!  using the Razzaque approach to determine an element sequence,
!  given a grid point sequence determined by GPS or CM.
!
!  The input data deck to BANDIT consists of a standard NASTRAN deck
!  plus one or more $ cards to tell BANDIT what to do.  Output includes
!  a set of SEQGP card images on bandit.f07 and the complete NASTRAN deck
!  (including SEQGP cards) on bandit.f08.  The new element sequence is on
!  bandit.f14.  See "BANDIT User's Guide" by G.C. Everstine.
!
!                    Program Installation
!                    --------------------
!  Set variable MEM (the amount of working storage in real, single-
!  precision words) to the desired value.
!
! B////////////////////////////////////////////////////////////////////B
      INTEGER, PARAMETER :: MEM = 10000000
! E////////////////////////////////////////////////////////////////////E
!
!  To get timing information, add in Subroutine TIMER a call to the
!  appropriate CPU clock routine.
!
!  Integer-packing is not available in this version, but may be
!  enabled by acquiring from the author the packing routines and
!  replacing some code wherever the string "CPACK" appears.
!  With large memories, integer-packing is rarely needed today and
!  inhibits vectorization.
!
!---------------------------------------------------------------------
!
!  I/O FILES - - -
!
!  (SEE BLOCK DATA FOR DEFINITION OF FILES IN COMMON BLOCK /IOUNIT/)
!
!  Unit  File name       Format       I/O     Opn STAT                                  Use
!  ----  -----------   -----------   ------   --------  ---------------------------------------------------------------------------
!    5   INFILE          FORMATTED   input    OLD        MYSTRAN input file
!    6   bandit.out'     FORMATTED   output   REPLACE    Bandit printed output file
!    7   bandit.f07'     FORMATTED   output   REPLACE    SEQGP card images
!    8   bandit.f08'     FORMATTED   output   REPLACE    MYSTRAN input file with SEQGP card images
!    9   bandit.f09'     FORMATTED   output   REPLACE    References in NASNUM, SPRING, and FINISH
!   10   bandit.ins'     FORMATTED   input    UNKNOWN    Insert file (card images)
!   11   bandit.f11'     FORMATTED   output   REPLACE    Scratch file - references in DOLLAR, NASNUM, and NEWIN
!   12   bandit.f12'   UNFORMATTED   output   REPLACE    Scratch file - references in MPC, RESET, RIGID, and TIGER
!   13   bandit.f13    UNFORMATTED   output   SCRATCH    Scratch file - store elems for generating elem ordering for frontal solver
!   14   bandit.f14'     FORMATTED   output   REPLACE    Output new element list for frontal solution
!   15   bandit.f15'     FORMATTED   output   REPLACE    Component and grid list (CM must be run to get this list)
!   16   bandit.f16'     FORMATTED   output   REPLACE    Connection table output ($TABLE YES)
!   17   bandit.log'     FORMATTED   output   REPLACE    Some run-time messages to benefit interactive running
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IADD   ,IB     ,IBYTE  ,IDIM   ,IER    ,IFIR   ,        &
               IFL    ,IIG    ,IGNORE ,IGDEG  ,IH     ,INP    ,IOP    ,        &
               IPARAM ,IPASS  ,ISTA   ,ISTART ,II1    ,II3    ,IWALL1 ,        &
               IWALL2
       
      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20
      
      INTEGER  K1     ,K2     ,K3     ,K4     ,K5     ,K6     ,K7     ,        &
               K8     ,K9     ,KDIM4  ,KDIM   ,KMOD   ,KNEW   ,KORE   ,        &
               KORIG
      
      INTEGER  LINES
      
      INTEGER  MA     ,MB     ,MAXDEG ,MAXGRD ,MINDEG ,MDIM   ,ME     ,        &
               MM
      
      INTEGER  NAXIC  ,NBITIN ,NBYTE  ,NCM    ,NEDGE  ,NEL    ,NELEM  ,        &
               NEQ    ,NEQR   ,NGRID  ,NLINK  ,NN     ,NW     ,NTYPE  ,        &
               NUM    ,NZERO

      REAL     DUMG   ,DUMW   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
      INTEGER   NEW_BW
      REAL      DEN
! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'BANDIT'
      INTEGER(LONG), INTENT(IN)       :: MYSTRAN_NGRID
! E////////////////////////////////////////////////////////////////////E

      COMMON /S/ NN,MM,IH,IB,LINES,NEDGE,IADD,MINDEG,NAXIC
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /B/ IPARAM(20)
      COMMON /BITS/ NBITIN,KORE,IFL,NGRID,IPASS,NW,NBYTE,IBYTE,KDIM
      COMMON /D/ KORIG,KNEW,IOP,INP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /W/ DUMW(6)
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),               &
                    NELEM(160),MDIM
      INTEGER vype,TYPE,WYPE
      COMMON /GRA/ DUMG(3)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      integer KOM(MEM)

! B////////////////////////////////////////////////////////////////////B
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9001) SUBR_NAME,TSEC
 9001    FORMAT(1X,A,' BEGN ',F10.3)
      ENDIF

! E////////////////////////////////////////////////////////////////////E
!
      IPARAM(1)=0
!
!     OPEN FILES (EXCEPT MAYBE 5 AND 6)
!
! B////////////////////////////////////////////////////////////////////B
      IOU5 = IN1              ! IOU5 IS MYSTRAN INPUT FILE. ALREADY OPEN
! E////////////////////////////////////////////////////////////////////E
      CALL BANDIT_FILES ( IOU6 , IOU7 , IOU8 , IOU9 , IOU11, IOU12,            &
                          IOU13, IOU14, IOU15, IOU16, IOU17 )
!
      OPEN(IOU6,FILE='bandit.out',FORM='FORMATTED',STATUS='replace')

      OPEN(IOU7,FILE='bandit.f07',FORM='FORMATTED',STATUS='replace')

      OPEN(IOU8,FILE='bandit.f08',FORM='FORMATTED',STATUS='replace')

      OPEN(IOU9,FILE='bandit.f09',FORM='FORMATTED',STATUS='replace')

!     OPEN(IOU11,FORM='FORMATTED',STATUS='SCRATCH')
!     OPEN(IOU12,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(IOU16,FILE='bandit.f16',FORM='FORMATTED',STATUS='replace')
!     OPEN(IOU17,FILE='bandit.log',FORM='FORMATTED',STATUS='replace')
!     See below for opening of scratch files
!     (placed there to avoid unexplained problem when looping on
!     multiple data sets)
!
!     TOP OF LOOP ON MULTIPLE DATA DECKS (used for code testing).
!     USE '$LOOP YES' IN FIRST DATA DECK TO ENABLE LOOPING AND
!     '$LOOP NO' IN LAST DATA DECK TO DISABLE LOOPING.
!
25    CALL TIMER(TA,IWALL1,0,IOU6)
      IF(IPARAM(1).EQ.4) write(iou6,'(80(1H#))')
      WRITE(IOU6,30)
30    FORMAT('Bandit 1/8/04, G.C. Everstine, Gaithersburg, Maryland,',         &
             ' geversti@comcast.net')
!
!     Open scratch files as named files.
!
      OPEN(IOU11,file='bandit.f11',FORM='FORMATTED',STATUS='replace')
      OPEN(IOU12,file='bandit.f12',FORM='UNFORMATTED',STATUS='replace')
!
!     INITIALIZE SOME VARIABLES.
!
      KORE=MEM
      IFL=KORE
      IER=0
!     NAXIC=NUMBER OF HARMONICS ON AXIC CARD (DEFAULT SET HERE)
      NAXIC=99999
!
!     SET DEFAULT NUMBER OF PRINTED LINES PER PAGE.
!
      LINES=55
!     USE 55 FOR STANDARD 11 IN. PAPER, 40 FOR 8.5 IN. PAPER.
!
!     SET DEFAULT FOR KDIM AND MAX NUMBER OF GRID POINTS.
!
! B////////////////////////////////////////////////////////////////////B
!     KDIM=250
      KDIM=MAX(250,MYSTRAN_NGRID)
! E////////////////////////////////////////////////////////////////////E
      KDIM4=4*KDIM
      KORE=KORE-KDIM4
! B////////////////////////////////////////////////////////////////////B
!     NGRID=KORE/30
      NGRID = MAX(MYSTRAN_NGRID,10)
! E////////////////////////////////////////////////////////////////////E
      KORE=KORE+KDIM4
!
!     READ EXECUTIVE AND CASE CONTROL DECKS.
!
      CALL CASE(IER)
      IF(IER.GT.0) GO TO 135
!
!     PRINT ELEMENT LIBRARY IF REQUESTED.
!
!     ICHAR=MA(3)
      IF(IPARAM(9).EQ.4) then
         WRITE(IOU6,40) (I,vype(i),TYPE(I),WYPE(I),ME(I),I=4,NTYPE)
40       FORMAT(/,'Element Library:'/(I6,3X,A1,2A4,I7))
      end if
!
!     INITIALIZE ELEMENT COUNTERS.
!
      DO I=1,NTYPE
         NELEM(I)=0
      end do
!
      IF(IPARAM(5).EQ.4) GO TO 80
!
!     COPY BULK DATA TO UNIT 8 IF RESEQUENCING IS NOT REQUESTED.
!
      CALL NOSEQ(KOM)
!
      GO TO 120
!
!     COMPUTE MAXGRD AND MAXDEG.
!
80    CONTINUE
      KDIM4=4*KDIM
      KORE=KORE-KDIM4
      CALL GRID(IER,IOU6)
      IF(IER.GT.0) GO TO 135
!
!     PRINT CORE ALLOCATION INFORMATION.
!
      CALL COREKO
!
!     PARTITION BLANK COMMON.
!
      II1=MAXGRD/NW
      II3=2*MAXGRD
      K1=1
      K2=K1+KDIM4
      K3=K2+MAXGRD+1
      K4=K3+2*II3+1
      K5=K4+MAXGRD
      K6=K5+MAXGRD
      K7=K6+MAXGRD
      K8=K7+MAXDEG+1
      K9=K8+II1*MAXDEG
!     (K9 DEFINES BEGINNING OF NEXT ARRAY IF THERE WERE ONE.)
!     IG, THE BIG ARRAY IN NASNUM, IS LOCATED LAST IN OPEN CORE TO KEEP
!     CALLING ADDRESSES SMALLER, A BENEFIT ON UNIVAC WHEN USING OVER 65K.
!
!     READ BULK DATA, SET UP CONNECTION TABLE, RESEQUENCE NODES, AND
!            GENERATE SEQGP CARDS.
!
      CALL NASNUM(KOM(K8),II1,KOM(K3),II3,KOM(K4),KOM(K5),KOM(K2),             &
                  KOM(K6),KOM(K7),KOM(K1),KDIM4,IER)
      IF(IER.GT.0) GO TO 135
!
!     PRINT OUT ELEMENT COUNTERS.
!
      WRITE(IOU6,90)
90    FORMAT(/,'Element Counts for Data Deck:')
      DO I=4,NTYPE
         IF(NELEM(I).gt.0)                                                     &
              WRITE(IOU6,100) vype(i),TYPE(I),WYPE(I),NELEM(I)
100      FORMAT(5X,A1,2A4,I8)
      end do
!
!     PRINT BANDIT SUMMARY.
!
! B////////////////////////////////////////////////////////////////////B
      CALL SUMUP ( NEW_BW, DEN)
! E////////////////////////////////////////////////////////////////////E
!
120   CONTINUE
      ENDFILE IOU8
      REWIND IOU8
      IF(IPASS.GT.0) WRITE(IOU6,130) IPASS
130   FORMAT(/,'Number of calls to the pack/unpack routines',I15)
!
!     WRAP UP JOB
!
!     Reset iou5 in case there is looping on multiple data sets.
! B////////////////////////////////////////////////////////////////////B
!!135 IOU5=5
  135 IOU5 = IN1
! E////////////////////////////////////////////////////////////////////E
!     delete scratch files
      close(iou11,status='delete')
      close(iou12,status='delete')
      CALL TIMER(TB,IWALL2,0,IOU6)
      TB=TB-TA
      IWALL2=IWALL2-IWALL1
!     IF(IWALL2.GT.0) WRITE(IOU6,138) IWALL2
138   FORMAT(/,'Elapsed time =',I6,' seconds.')
      WRITE(IOU6,140) TB
140   FORMAT(/,'End of BANDIT Job.  Total CP time =',F12.3,' seconds.')
      write(iou6,*)
!     IF LOOPING ON MULTIPLE DECKS IS REQUESTED, GO BACK TO BEGINNING.
      IF(IPARAM(1).EQ.4) GO TO 25
      if(iparam(7).eq.3) close(iou9,status='delete')
      if(iparam(2).eq.3) close(iou16,status='delete')
!     close(iou17,status='delete')
!     if(ier.gt.0) stop 13
!     $NASTRAN YES CARD RETURNS CONDITION CODE 5 FOR USE BY IBM.
!     (IBM JCL was written to check Bandit's exit code; now obsolete.
!     Stop 13 above served same purpose.)
!     IF(IPARAM(8).EQ.4) STOP 5
! B////////////////////////////////////////////////////////////////////B
      close(iou6,status="keep")
      close(iou7,status="keep")
      close(iou8,status="keep")
 9000 continue
      IF (WRT_LOG >= SUBR_BEGEND) THEN
         CALL OURTIM
         WRITE(F04,9002) SUBR_NAME,TSEC
 9002    FORMAT(1X,A,' END  ',F10.3)
      ENDIF

      RETURN
! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
      END SUBROUTINE BANDIT
! E////////////////////////////////////////////////////////////////////E

! ##################################################################################################################################
      SUBROUTINE BRIGIT(INV,II3,INT,ILD,IER)
!
!     THIS ROUTINE SORTS THE ORIGINAL GRID NUMBERS AND OUTPUTS THE LIST
!         IN INT, WHERE INT(I)=THE ITH ORIGINAL GRID NUMBER.
!     ALSO OUTPUT IS ILD, WHERE ILD(I)=SORTED INTERNAL NUMBER
!         CORRESPONDING TO THE UNSORTED BANDIT INTERNAL LABEL I.
!
!     INPUT - INV
!     OUTPUT - ILD,INT
!
      INTEGER II3, INV(2,II3),INT(*),ILD(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER I      ,IER    ,IS     ,J      ,KFAC   ,KFM   ,                  &
              KMOD   ,L      ,MAXDEG ,MAXGRD ,MINI   ,NN
              
      REAL     DUMS
! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
!
!     PERFORM A ROUGH SORT OF THE ORIGINAL GRID NUMBERS.
!
      L=0
      KFAC=-1
   20 KFAC=KFAC+1
      MINI=99999999
      KFM=KFAC*KMOD
      DO 50 I=1,KMOD
      IF(INV(1,I).GT.KFM) MINI=MIN(MINI,INV(1,I))
   50 CONTINUE
      KFAC=(MINI-1)/KMOD
      DO 80 I=1,KMOD
      IS=INV(1,I)
      IF(IS.LE.(KFAC*KMOD).OR.IS.GT.(KFAC+1)*KMOD)GO TO 80
      L=L+1
      INT(L)=INV(1,I)
   80 CONTINUE
      IF(L.LT.NN)GO TO 20
!
!     COMPLETE THE SORTING OF THE ORIGINAL GRID NUMBERS.
!
      CALL SORT(INT,NN)
!
!     DETERMINE CORRESPONDENCE (ILD) BETWEEN NORIG AND INT ARRAYS.
!
      DO 40 I=1,NN
      J=INT(I)
      L=INTERN(J,INV,II3,IER)
      IF(IER.GT.0) RETURN
!     J IS DEFINED ABOVE TO ALLOW EXECUTION ON UNIVAC IN OVER 65K.
   40 ILD(L)=I
      RETURN
      END SUBROUTINE BRIGIT

! ##################################################################################################################################
      SUBROUTINE CASE(IER)
!
!     READ EXECUTIVE AND CASE CONTROL DECKS.
!
      INTEGER KA(80)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IER    ,IPARAM ,KT     ,L      ,MA     ,MB     ,        &
               NUM 

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      KT=0
      REWIND IOU8
      WRITE(IOU6,10)
   10 FORMAT(/,'Echo of Data Through BEGIN BULK Card:')
!
!     INITIALIZE PARAMETERS.
!
      CALL DOLLAR(KA,1,KT,IER)
      IF(IER.GT.0) RETURN
!
!     READ CARD.
!
   20 READ(IOU5,30,END=100) KA
   30 FORMAT(80A1)
      KT=KT+1
      L=LENCAS(KA,80)
      WRITE(IOU6,60) KT,(KA(I),I=1,L)
60    FORMAT(I8,'- ',80A1)
      WRITE(IOU8,30) (KA(I),I=1,L)
!
!     IF COLUMN 1 IS $, PROCESS $-CARD.
!
      IF(KA(1).EQ.MB(1)) then
         CALL DOLLAR(KA,2,KT,IER)
         IF(IER.GT.0) RETURN
         GO TO 20
      end if
!
!     IF THE CARD IS NOT BEGIN BULK, GO BACK AND READ ANOTHER CARD.
!
      IF(NBULK(KA).EQ.0) GO TO 20
!
!     CHECK FOR ILLEGAL PARAMETERS.
!
   80 CALL DOLLAR(KA,3,KT,IER)
      IF(IER.GT.0) RETURN
      RETURN
!
!     END-OF-FILE ENCOUNTERED
!
100   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE CASE

! ##################################################################################################################################
      SUBROUTINE COREKO
!
!     PRINT SUMMARY OF CORE ALLOCATION.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  IBYTE  ,IFL    ,KDIM   ,KMOD   ,KORE   ,LL     ,MAXDEG ,        &
               MAXGRD ,NBITIN ,NBYTE  ,NW
     
      REAL     DUM 
! E////////////////////////////////////////////////////////////////////E
      COMMON /BITS/ NBITIN,KORE,IFL,DUM(2),NW,NBYTE,IBYTE,KDIM
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      LL=IFL-KORE
      WRITE(IOU6,10) KORE,MAXGRD,MAXDEG,KDIM
10    FORMAT(/,'Working Storage:'/                                             &
       5X,'Length of Working Storage',I10,' words'/                            &
       5X,'Grid Point Limit         ',I10/                                     &
       5X,'Nodal Degree Limit       ',I10/                                     &
       5X,'$DIMENSION Value         ',I10)
      IF(NW.GT.1) WRITE(IOU6,20) NW
20    FORMAT(5X,'Packing Density',I20,' integers/word')
      RETURN
      END SUBROUTINE COREKO

! ##################################################################################################################################
      SUBROUTINE CUTHIL(NT,NUM,NOM,IO,IP,IG,II1,IC,IDEG,IDIS,IW,               &
        NEW,ICC,ILD,IPP,JUMP,NODESL,KORDIM,IER)
!
! THIS IS THE EXECUTIVE FOR THE CUTHILL-MCKEE GRID POINT RENUMBERING
!      STRATEGY.   THE PRINCIPAL INPUTS ARE THE CONNECTIVITY MATRIX IG
!      AND THE NUMBER OF GRID POINTS (NODES) NN.
!
! INPUT -- NT,NUM,NOM,IO,IP,IG,II1,NN,MAXGRD,ILD,NBITIN,ISTART,ISTA
! OUTPUT -- NEW,ILD,MM,IH0,IHE,KORIG,KNEW,NCM
! SCRATCH -- IC,IDEG,IDIS,IW,ICC,IPP
! SET FOLLOWING DIMENSIONS IN CALLING PROGRAM --
!   IG(II1,M),IC(L),IDEG(L),IDIS(L),IW(L),NEW(L),ICC(L),ILD(L),IP(M)
!        L=MAXGRD EXCEEDS NUMBER OF GRID POINTS
!        II1=MAXGRD/(PACKING DENSITY IN INTEGERS/WORD)
!        M EXCEEDS MAX NODAL DEGREE
!
! NT=MAX NUMBER OF STARTING NODES TO BE CONSIDERED
! NUM AND NOM GIVE THE FRACTION OF THE RANGE FROM MIN DEGREE TO MAX
!        DEGREE TO CONSIDER FOR STARTING NODES
! IO=STRATEGY OPTION - - -
!          1=RMS WAVEFRONT
!          2=BANDWIDTH
!          3=PROFILE
!          4=WAVEFRONT (MAX)
! IP=PRINTING OPTION (0 FOR NO PRINTED OUTPUT FROM CUTHILL)
! IG(I,J) CONTAINS THE GRID POINT LABEL FOR THE JTH NODE ADJACENT TO
!      NODE I  (THE CONNECTIVITY MATRIX).   THE CONNECTION OF A NODE
!      TO ITSELF IS NOT LISTED.
! II1=ROW DIMENSION OF IG
! NN=NUMBER OF GRID POINTS (NODES)
! MM=COLUMN DIMENSION OF IG ON INPUT, MAX NODAL DEGREE ON OUTPUT
! MAXGRD=EFFECTIVE ROW DIMENSION OF IG (NEGLECTING INTEGER PACKING)
! LINE=NUMBER OF PRINTED LINES PER PAGE
! NBITIN=NUMBER OF BITS PER INTEGER (FOR PACKING)
! ISTART(I)=ITH STARTING NODE SELECTED BY USER
! ISTA=NUMBER OF STARTING NODES IN ISTART
! NEW(I)=OLD LABEL FOR GRID POINT NOW LABELLED I
! ILD(I)=NEW LABEL FOR GRID POINT ORIGINALLY LABELLED I
!        ILD AND NEW ARE INVERSES
! ILD MUST BE INPUT TO CUTHILL TO INDICATE AN INITIAL SEQUENCE.
!           NORMALLY, ON INPUT, SET ILD(I)=I FOR ALL I.
! JUMP=1 IF RESEQUENCING ATTEMPTS RESULT IN NO IMPROVEMENT
!     =0 OTHERWISE.
! IH0=ORIG PROFILE
! IHE=NEW PROFILE
! KORIG=ORIG BANDWIDTH
! KNEW=NEW BW
! NCM=NUMBER OF COMPONENTS
!     NODESL IS SCRATCH SPACE.
! IN CALLING PROGRAM, TRY  CALL CUTHIL(80,1,2,2,1, . . . )
!
!    THE FOLLOWING SUBROUTINES WERE WRITTEN BY E. CUTHILL AND
!    J. MCKEE OF NSRDC - - -
!        DEGREE,DIAM,IDIST,KOMPNT,MAXBND,MAXDGR,MINDEG,RELABL,CUTHILL
!     CUTHILL WAS MODIFIED BY G. C. EVERSTINE, DTRC.
!
      INTEGER II1, KORDIM
      DIMENSION IG(II1,1),IC(*),IDEG(*),IDIS(*),IW(*)
      INTEGER NEW(*),ICC(*),ILD(*),IPP(*),NODESL(KORDIM),SUMW
      REAL IM1,IM2
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,         &
              IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,         &
              IOU19  ,IOU20

      INTEGER  I      ,IAJDIM ,IB     ,IC     ,IDEG   ,IDEM   ,IDEM1  ,        &
               IDIM   ,IDIS   ,IER    ,IG     ,IH     ,IH0    ,IHE    ,        &
               ij     ,IO     ,IP     ,IPARAM ,IS     ,ISTA   ,                &
               ISTART ,IW     ,IWALL

      INTEGER  J      ,JMAX   ,JUMP

      INTEGER  K      ,K2     ,KMOD   ,KNEW ,KORIG

      INTEGER  L

      INTEGER  M      ,MA     ,MAD    ,MAXD   ,MAXDEG ,MAXGRD ,MAXLEV ,        &
               MAXW   ,MAXW0  ,MAXW1  ,MEDIAN ,MI     ,MM     ,MODD

      INTEGER  NBITIN ,NC     ,NCM    ,NL     ,NN     ,NNODE  ,NOM    ,        &
                       NT     ,NUM

      REAL     AVERW  ,BRMS   ,BRMS0  ,BRMS1  ,CRIT1  ,CRIT2  ,                &
               DUMBB  ,DUMD   ,DUML   ,DUMO   ,DUMS   ,RMS    ,                &
               RMS0   ,RMS1   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,IH,IB,DUMS(5)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /D/ KORIG,KNEW,IH0,IHE,NCM,DUMD(5)
      COMMON /W/ MAXW0,RMS0,MAXW1,RMS1,BRMS0,BRMS1
      COMMON /B/ IPARAM(20)
      COMMON /DOL/ ISTART(100),DUMO(100)
      COMMON /DOLL/ IDIM,ISTA,DUML(3)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      CALL TIMER(TA,IWALL,0,IOU6)
!     SET UP SCRATCH SPACE NODESL.
      IDEM=KORDIM/4
      K2=IDEM+1
      IAJDIM=3 *IDEM
      JUMP=0
!     DETERMINE THE DEGREE OF EACH NODE.
      CALL DEGREE(IG,II1,IDEG)
!     DETERMINE THE NUMBER OF COMPONENTS, NCM.
      NCM=KOMPNT(IG,II1,IC,IDEG,IW,ICC)
!     DETERMINE THE MAXIMUM DEGREE OF ANY NODE.
      MAXD=MAXDGR(0,IC,IDEG)
      MM=MAXD
! INITIALIZE NEW ARRAY FROM THE ILD ARRAY.
!   ILD MUST BE INPUT TO CUTHILL.
      DO 30 I=1,NN
      K=ILD(I)
   30 NEW(K)=I
!   COMPUTE ORIGINAL BANDWIDTH AND PROFILE.
!   COMPUTE ORIGINAL WAVEFRONT AND ACTIVE COLUMN DATA.
      CALL WAVEY(IG,II1,ILD,NEW,0,IC,IW,IS,MAXW,AVERW,SUMW,RMS,BRMS)
      IH=SUMW
      MAXW0=MAXW
      RMS0=RMS
      BRMS0=BRMS
      KORIG=IS
!    IH0=ORIGINAL PROFILE, IS=ORIGINAL BW.
      IH0=IH
! COMPUTE NODAL DEGREE STATISTICS.
      CALL DIST(IDEG,IPP,IP,MEDIAN,MODD)
!     IF REQUESTED, PRINT INTERNAL NUMBER CONNECTION TABLE.
      IF(IP.ne.0) CALL STABLE(IG,II1,IC,IDEG,ILD,IPP)
      WRITE(IOU6,29)
   29 FORMAT(/,'Before Resequencing:')
      WRITE(IOU6,51) IS,IH,MAXW,AVERW,RMS,BRMS
   51 FORMAT(5X,'Bandwidth',I18/5X,'Profile',I20/                              &
             5X,'Max Wavefront',I14/5X,'Avg Wavefront',F14.3/                  &
             5X,'RMS Wavefront',F14.3/5X,'RMS Bandwidth',F14.3)
      IF(ip.ne.0.and.ISTA.gt.0) then
         WRITE(IOU6,701)
701      FORMAT('Starting nodes supplied by user:')
         WRITE(IOU6,100) (ISTART(I),I=1,ISTA)
      end if
! INITIALIZE ILD AND NEW ARRAYS.
      DO I=1,NN
         NEW(I)=0
         ILD(I)=0
      end do
!
!     GENERATE NUMBERING SCHEME FOR EACH COMPONENT, NC.
!
      DO 500 NC=1,NCM
!     DETERMINE THE RANGE OF DEGREES (MI  TO  MAD) OF NODES OF INTEREST.
      MI=MINDEG(NC,IC,IDEG)
      MAD=MI
      IF(NOM.EQ.0) GO TO 87
   90 MA=MAXDGR(NC,IC,IDEG)
      MAD=MI+((MA-MI)*NUM)/NOM
!     MAKE SURE MAD DOES NOT EXCEED MEDIAN.
      MAD=MIN(MAD,MEDIAN-1)
      MAD=MAX(MAD,MI)
!     DETERMINE BANDWIDTH OR SUM CRITERION FOR EACH NODE MEETING
!        SPECIFIED CONDITION.
   87 CONTINUE
      CALL DIAM(NC,MAD,NL,NODESL,IDEM,MAXLEV,IG,II1,IC,IDEG,IDIS,              &
        IW,ICC)
      IF(IP.EQ.0) GO TO 67
      WRITE(IOU6,39) NC,MAD
   39 FORMAT(/,'Component',I5,', Max Degree Used',I6)
!     COMPUTE THE NUMBER OF NODES IN THIS COMPONENT AND WRITE OUT.
      NNODE=ICC(NC+1)-ICC(NC)
      WRITE(IOU6,41) NNODE
   41 FORMAT('Number of nodes in this component',I8)
      WRITE(IOU6,59) MAXLEV
   59 FORMAT('Starting nodes for minmax number of nodes per level',I6)
      if(NL.gt.0) WRITE(IOU6,100) (NODESL(J),J=1,NL)
  100 FORMAT(10I7)
   67 CONTINUE
      IF(ISTA.LE.0) GO TO 760
      IDEM1=IDEM-1
      M=0
      DO 750 I=1,ISTA
         J=ISTART(I)
         IF(IC(J).NE.NC) GO TO 750
         M=M+1
         DO K=1,IDEM1
            L=IDEM+1-K
            NODESL(L)=NODESL(L-1)
         end do
         NODESL(1)=J
  750 CONTINUE
      NL=MIN(NL+M,IDEM)
      CALL FIXIT(NODESL,NL)
  760 CONTINUE
      IF(IP.EQ.0) GO TO 63
      IF(ISTA.LE.0) GO TO 63
      WRITE(IOU6,730)
  730 FORMAT('Merged list of starting nodes supplied by user',                 &
             ' and by BANDIT -')
      WRITE(IOU6,100) (NODESL(I),I=1,NL)
   63 CONTINUE
      JMAX=MIN(NT,NL)
      JMAX=MAX(JMAX,1)
      IM1=1.E8
      IM2=IM1
!
!  CHECK SEQUENCE FOR EACH STARTING NODE SELECTED.
!
      ij=1
      DO 400 J=1,JMAX
      CALL RELABL(1,NODESL(J),IG,II1,IC,IDEG,IDIS,IW,NEW,ICC,                  &
                  ILD,NODESL(K2),IAJDIM,IER)
      IF(IER.GT.0) RETURN
!  COMPUTE NEW BANDWIDTH,PROFILE,WAVEFRONT DATA.
      CALL WAVEY(IG,II1,ILD,NEW,NC,IC,IW,IB,MAXW,AVERW,SUMW,RMS,BRMS)
      IH=SUMW
!     IB=BANDWIDTH, IH=PROFILE.
      IF(IP.EQ.0) GO TO 70
      WRITE(IOU6,69) NODESL(J),IB,IH,MAXW,RMS
   69 FORMAT('Starting Node',I6,', Band',I6,', Profile',I8,                    &
       ', Max W',I6,', RMS W',F9.3 )
   70 CONTINUE
      GO TO (205,210,215,220), IO
  205 CRIT1=RMS
      CRIT2=IH
      GO TO 71
  210 CRIT1=IB
      CRIT2=IH
      GO TO 71
  215 CRIT1=IH
      CRIT2=IB
      GO TO 71
  220 CRIT1=MAXW
      CRIT2=RMS
      GO TO 71
   71 CONTINUE
      IF(IM1-CRIT1) 400,350,300
  300 IM1=CRIT1
      IM2=CRIT2
      IJ=J
      GO TO 400
  350 IF(IM2.LE.CRIT2) GO TO 400
      IM2=CRIT2
      IJ=J
400   CONTINUE
!
!   RECOMPUTE SEQUENCE FOR STARTING NODE WHICH IS BEST FOR CRITERION
!       SELECTED.
!
      CALL RELABL(1,NODESL(IJ),IG,II1,IC,IDEG,IDIS,IW,NEW,ICC,                 &
                  ILD,NODESL(K2),IAJDIM,IER)
      IF(IER.GT.0) RETURN
!
  500 CONTINUE
!
!
!  DETERMINE NODES OF ZERO DEGREE AND STACK LAST.
      CALL STACK(IDEG,NEW,ILD,IW)
!   COMPUTE BANDWIDTH, PROFILE AND WAVEFRONT DATA.
      CALL WAVEY(IG,II1,ILD,NEW,0,IC,IW,IB,MAXW,AVERW,SUMW,RMS,BRMS)
      IH=SUMW
!
      WRITE(IOU6,705)
  705 FORMAT(/,'After resequencing by Reverse Cuthill-McKee (CM):')
      WRITE(IOU6,51) IB,IH,MAXW,AVERW,RMS,BRMS
!
!   CHECK CM LABELING AGAINST ORIGINAL LABELING TO SEE IF BETTER.
!
      GO TO (255,260,265,270), IO
  255 IM1=RMS0
      IM2=IH0
      CRIT1=RMS
      CRIT2=IH
      GO TO 711
  260 IM1=IS
      IM2=IH0
      CRIT1=IB
      CRIT2=IH
      GO TO 711
! IB = BANDWIDTH,  IH= PROFILE.
  265 IM1=IH0
      IM2=IS
      CRIT1=IH
      CRIT2=IB
      GO TO 711
  270 IM1=MAXW0
      IM2=RMS0
      CRIT1=MAXW
      CRIT2=RMS
      GO TO 711
  711 CONTINUE
      IF(CRIT1-IM1) 715,742,744
  742 IF(CRIT2.LT.IM2) GO TO 715
!   IF NO IMPROVEMENT RETURN TO ORIGINAL SEQUENCE.
  744 IB=IS
      IH=IH0
      MAXW=MAXW0
      RMS=RMS0
      BRMS=BRMS0
      JUMP=1
      DO 713 I=1,NN
      ILD(I)=I
  713 NEW(I)=I
!
  715 CONTINUE
!   SET FINAL VALUES OF B , P , RMS , W .
      KNEW = IB
      IHE=IH
      MAXW1=MAXW
      RMS1=RMS
      BRMS1=BRMS
      CALL TIMER(TB,IWALL,0,IOU6)
      TB=TB-TA
!     IF(TB.GT.1.E-5) WRITE(IOU6,610) TB
      WRITE(IOU6,610) TB
  610 FORMAT(5X,'CP time',F20.3)
      RETURN
      END SUBROUTINE CUTHIL

! ##################################################################################################################################
      SUBROUTINE DEGREE(IG,II1,IDEG)
!
!     SET UP THE IDEG ARRAY CONTAINING THE DEGREE OF EACH NODE STORED
!     IN THE IG ARRAY.
!     IDEG(I) = DEGREE OF NODE I
!
      INTEGER II1
      DIMENSION IG(II1,1),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
      INTEGER  I      ,IDEG   ,IG    ,J      ,KMOD   ,MAXDEG ,                 &
               MAXGRD ,MM     ,NN     ,NBITIN

      REAL     DUMBB  ,DUMS
 
! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      DO 100 I=1,NN
      IDEG(I)=0
      DO 80 J=1,MM
      IF(IG(I,J)) 100,100,50
!PACK IF(IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)) 100,100,50
   50 IDEG(I)=IDEG(I)+1
   80 CONTINUE
  100 CONTINUE
      RETURN
      END SUBROUTINE DEGREE

! ##################################################################################################################################
      SUBROUTINE DGREE(NDSTK,NR,NDEG,IOLD,IBW1,IPF1)
!
!     DGREE COMPUTES THE DEGREE OF EACH NODE IN NDSTK AND STORES
!     IT IN THE ARRAY NDEG.  THE BANDWIDTH AND PROFILE FOR THE ORIGINAL
!     OR INPUT RENUMBERING OF THE GRAPH IS COMPUTED ALSO.
!
!     COMPUTE MAXIMUM DEGREE MM AND STORE IN IDEG.
!
! B////////////////////////////////////////////////////////////////////B
      INTEGER  I      ,IBW1   ,IDEG   ,IDIF   ,IDPTH  ,IOLD   ,IPF1   ,        &
               IRW    ,ITST   ,J      ,KMOD
     
      INTEGER  MAXDEG ,MAXGRD ,MM     ,N      ,NBITIN ,NDEG   ,NDSTK  ,        &
               NN     ,NR

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /GRA/ N,IDPTH,IDEG
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /S/ NN,MM,DUMS(7)
      DIMENSION NDSTK(NR,1),NDEG(*),IOLD(*)
      IBW1=0
      IPF1=0
      IDEG=MM
      MM=0
      DO 100 I=1,N
        NDEG(I)=0
        IRW=0
        DO 80 J=1,IDEG
          ITST=NDSTK(I,J)
!PACK     ITST=IUNPK(NDSTK,MAXGRD*(J-1)+I,NBITIN)
          IF(ITST.EQ.0) GO TO 90
   50     NDEG(I)=NDEG(I)+1
          IDIF=IOLD(I)-IOLD(ITST)
          IF(IRW.LT.IDIF) IRW=IDIF
          MM=MAX(MM,J)
   80   CONTINUE
   90   IPF1=IPF1+IRW
        IF(IRW.GT.IBW1) IBW1=IRW
  100 CONTINUE
      IDEG=MM
!     INCLUDE DIAGONAL TERMS IN BANDWIDTH AND PROFILE
      IBW1=IBW1+1
      IPF1=IPF1+N
      RETURN
      END SUBROUTINE DGREE

! ##################################################################################################################################
      SUBROUTINE DIAM(NC,MAXDG,NL,NODESL,IDEM,MAXLEV,                          &
                      IG,II1,IC,IDEG,IDIS,IW,ICC)
!
!     DETERMINE NL STARTING POINTS AND STORE IN NODESL.
!
      INTEGER II1
      DIMENSION IG(II1,*),IDIS(*),IW(*),ICC(*),IC(*),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
      INTEGER  I      ,IC     ,ICC    ,IDEG   ,IDEM   ,IDIS   ,IG     ,        &
               IW     ,KMOD   ,MAXDEG ,MAXDG  ,MAXGRD ,MAXLEV ,                &
               MD     ,ML     ,NBITIN ,NC     ,NL     ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER NODESL(*)
      NL=0
      MAXLEV=600000
      DO 100 I=1,NN
      IF((NC-IC(I)).NE.0) GO TO 100
40    IF((MAXDG-IDEG(I)).LT.0) GO TO 100
  105 MD=IDIST(I,ML,MAXLEV,IG,II1,IC,IDEG,IDIS,IW,ICC)
      IF(MD) 115,115,56
   56 IF(ML-MAXLEV)58,64,100
   58 MAXLEV=ML
      NL=1
      NODESL(1)=I
      GO TO 100
   64 IF(NL.GE.IDEM) GO TO 100
      NL=NL+1
      NODESL(NL)=I
  100 CONTINUE
  110 RETURN
  115 ML=1
      NODESL(1)=I
      MAXLEV=0
      RETURN
      END SUBROUTINE DIAM

! ##################################################################################################################################
      SUBROUTINE DIST(IDEG,HIST,IP,MEDIAN,MODD)
!
!     COMPUTE AND PRINT THE DISTRIBUTION OF NODAL DEGREES WITH MEDIAN
!        AND MODE.
!
!     IDEG(I) = DEGREE OF NODE I
!     HIST(I) = NUMBER OF NODES OF DEGREE I
!     IP      = PRINT OPTION PARAMETER (IF 0, NO PRINTED OUTPUT)
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IP     ,ISUM   ,K      ,MAXI   ,MEDIAN ,MM     ,        &
               MM1    ,MODD   ,NN     ,NN2

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER IDEG(*),HIST(*)
!
!     COMPUTE HISTOGRAM.
!
      MM1=MM+1
      DO 10 I=1,MM1
   10 HIST(I)=0
      DO 20 I=1,NN
      K=IDEG(I)+1
   20 HIST(K)=HIST(K)+1
!
!     COMPUTE MODE (MODD).
!
      MODD=0
      MAXI=0
      DO 25 I=1,MM1
      K=HIST(I)
      IF(K.LE.MAXI) GO TO 25
      MAXI=K
      MODD=I-1
   25 CONTINUE
      IF(IP.EQ.0) GO TO 60
!
!     PRINT HISTOGRAM.
!
      WRITE(IOU6,30)
   30 FORMAT(/,'Distribution of Nodal Degrees:'//5X,                           &
             'Degree  Number  Cum. Total')
      ISUM=0
      DO 40 I=1,MM1
      ISUM=ISUM+HIST(I)
      K=I-1
   40 WRITE(IOU6,50) K,HIST(I),ISUM
   50 FORMAT(3X,2I8,I12)
!
!     COMPUTE CUMULATIVE DISTRIBUTION.
!
   60 DO 70 I=2,MM1
   70 HIST(I)=HIST(I)+HIST(I-1)
!     COMPUTE MEDIAN.
      NN2=NN/2
      DO 80 I=1,MM1
      IF(HIST(I).GT.NN2) GO TO 90
   80 CONTINUE
   90 MEDIAN=I-1
      IF(IP.NE.0) WRITE(IOU6,100) MEDIAN,MODD
  100 FORMAT(/5X,'Median',I6/5X,'  Mode',I6)
      RETURN
      END SUBROUTINE DIST

! ##################################################################################################################################
      SUBROUTINE DOLLAR(KA,JUMP,KT,IER)
!
!     INTERPRET DOLLAR SIGN CARD.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,         &
              IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,         &
              IOU19  ,IOU20

      INTEGER  I      ,IADD   ,IDIM   ,IER    ,IFIR   ,IFL    ,IFLD   ,        &
               IGDEG  ,IGNORE ,IIG    ,IPARAM ,ISTA   ,ISTART ,ITYPE  ,        &
               J      ,JUMP   ,K      ,KDIM   ,KMOD   ,KORE   ,KT     ,        &
               LINES                                                  ,        &
               MA     ,MAXDEG ,MAXGRD ,MAXI   ,MB     ,MDIM   ,ME     ,        &
               NBITIN ,NCON   ,NEDGE  ,NELEM  ,NGRID  ,NIP    ,NTYPE  ,        &
               NUM

      REAL     DMY    ,DUM    ,DUMS   ,DUMY

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /B/ IPARAM(20)
      COMMON /S/ DMY(4),LINES,NEDGE,IADD,DUMS,DUMY
      COMMON /BITS/ NBITIN,KORE,IFL,NGRID,DUM(4),KDIM
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),               &
                    NELEM(160),MDIM
      INTEGER vype,TYPE,WYPE,KA(80),IP(35)
      DATA MAXI/35/
!
      GO TO (10,30,80), JUMP
!
!     INITIALIZE PARAMETERS.
!
!     IF LOOP PARAMETER IS ALREADY SET ($LOOP YES), INDICATING THAT
!     MULTIPLE DATA SETS ARE TO BE PROCESSED IN ONE EXECUTION, DO NOT
!     RESET IPARAM(1).
!
10    J=1
      IF(IPARAM(1).EQ.4) J=2
      DO 20 I=J,20
20    IPARAM(I)=0
      IDIM=100
      IIG=0
      ISTA=0
      IGDEG=0
      IADD=0
      RETURN
!
!* PARAMETERS - - - -
!
!        KEYWORD 1           KEYWORD 2
!        1 - LOOP            1 - SEQGP
!        2 - TABLE           2 - ALL
!        3 - METHOD          3 - NO
!        4 - MPC             4 - YES
!        5 - SEQUENCE        5 - MIN
!        6 - CRITERION       6 - MAX
!        7 - SPRING          7 - CM
!        8 - NASTRAN         8 - GPS
!        9 - ELEMENTS        9 - BOTH CM AND GPS
!       10 - PRINT          10 - BANDWIDTH
!       11 - FRONTAL        11 - PROFILE
!                           12 - RMS WAVEFRONT
!                           13 - WAVEFRONT (MAX)
!
!     REMAINING PARAMETERS UP TO 20 IN LEFT COLUMN ARE CURRENTLY UNUSED
!
!*   OTHER KEYWORDS - -   IGNORE,GRID,START,DEGREE,INSERT,LINES,
!             DIMENSION,ADD,APPEND
!
!
!     LOOK FOR FIRST KEYWORD.
!
   30 ITYPE=0
!     CHECK COLUMN 2 FOR BLANK.
      IF(KA(2).EQ.MB(2)) RETURN
!     SEE BLOCK DATA FOR ALPHABET KEY.
!    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
!    A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
      IF(KA(2).EQ.MA(12).AND.KA(3).EQ.MA(15)) ITYPE=1      ! LO $LOOP
      IF(KA(2).EQ.MA(20).AND.KA(3).EQ.MA( 1)) ITYPE=2      ! TA $TABLE
      IF(KA(2).EQ.MA(13).AND.KA(3).EQ.MA( 5)) ITYPE=3      ! ME $METHOD
      IF(KA(2).EQ.MA(13).AND.KA(3).EQ.MA(16)) ITYPE=4      ! MP $MPC
      IF(KA(2).EQ.MA(19).AND.KA(3).EQ.MA( 5)) ITYPE=5      ! SE $SEQUENCE
      IF(KA(2).EQ.MA( 3).AND.KA(3).EQ.MA(18)) ITYPE=6      ! CR $RITERION
      IF(KA(2).EQ.MA(19).AND.KA(3).EQ.MA(16)) ITYPE=7      ! SP $SPRING
      IF(KA(2).EQ.MA(14).AND.KA(3).EQ.MA( 1)) ITYPE=8      ! NA $
      IF(KA(2).EQ.MA( 5).AND.KA(3).EQ.MA(12)) ITYPE=9      ! EL $ELEMENTS
      IF(KA(2).EQ.MA(16).AND.KA(3).EQ.MA(18)) ITYPE=10     ! PR $PRINT
      IF(KA(2).EQ.MA( 6).AND.KA(3).EQ.MA(18)) ITYPE=11     ! FR $FRONTAL
      IF(ITYPE.GT.0) GO TO 35
      IF(KA(2).EQ.MA( 7).AND.KA(3).EQ.MA(18)) GO TO 120    ! GR $GRID
      IF(KA(2).EQ.MA( 9).AND.KA(3).EQ.MA( 7)) GO TO 100    ! IG $IGNORE
      IF(KA(2).EQ.MA(19).AND.KA(3).EQ.MA(20)) GO TO 130    ! ST $START
      IF(KA(2).EQ.MA( 4).AND.KA(3).EQ.MA( 5)) GO TO 150    ! DE $DEGREE
      IF(KA(2).EQ.MA( 9).AND.KA(3).EQ.MA(14)) GO TO 160    ! IN $INSERT
      IF(KA(2).EQ.MA(12).AND.KA(3).EQ.MA( 9)) GO TO 180    ! LI $LINES
      IF(KA(2).EQ.MA( 4).AND.KA(3).EQ.MA( 9)) GO TO 190    ! DI $DIMENSION
      IF(KA(2).EQ.MA( 1).AND.KA(3).EQ.MA( 4)) GO TO 230    ! AD $ADD
      IF(KA(2).EQ.MA( 1).AND.KA(3).EQ.MA(16)) GO TO 240    ! AP $APPEND
      RETURN
!
!     LOOK FOR SECOND KEYWORD.
!
   35 DO I=4,70
         IF(KA(I).EQ.MB(2)) GO TO 50
      end do
      RETURN
   50 K=I+1
      DO I=K,71
         IF(KA(I).NE.MB(2)) GO TO 70
      end do
      RETURN
!
70    J=0
      IF(KA(I).EQ.MA(19).AND.KA(I+1).EQ.MA( 5)) J=1        ! SE
      IF(KA(I).EQ.MA( 1).AND.KA(I+1).EQ.MA(12)) J=2        ! AL
      IF(KA(I).EQ.MA(14).AND.KA(I+1).EQ.MA(15)) J=3        ! NO
      IF(KA(I).EQ.MA(25).AND.KA(I+1).EQ.MA( 5)) J=4        ! YE
      IF(KA(I).EQ.MA(13).AND.KA(I+1).EQ.MA( 9)) J=5        ! MI
      IF(KA(I).EQ.MA(13).AND.KA(I+1).EQ.MA( 1)) J=6        ! MA
      IF(KA(I).EQ.MA( 3).AND.KA(I+1).EQ.MA(13)) J=7        ! CM
      IF(KA(I).EQ.MA( 7).AND.KA(I+1).EQ.MA(16)) J=8        ! GP
      IF(KA(I).EQ.MA( 2).AND.KA(I+1).EQ.MA(15)) J=9        ! BO
      IF(KA(I).EQ.MA( 2).AND.KA(I+1).EQ.MA( 1)) J=10       ! BA
      IF(KA(I).EQ.MA(16).AND.KA(I+1).EQ.MA(18)) J=11       ! PR
      IF(KA(I).EQ.MA(18).AND.KA(I+1).EQ.MA(13)) J=12       ! RM
      IF(KA(I).EQ.MA(23).AND.KA(I+1).EQ.MA( 1)) J=13       ! WA
      IF(J.GT.0) IPARAM(ITYPE)=J
      RETURN
!
!     CHECK FOR ILLEGAL PARAMETERS AND SET TO DEFAULTS.
!
80    IF(IPARAM(1).NE.4)                    IPARAM(1)=3
      IF(IPARAM(2).NE.4)                    IPARAM(2)=3
      IF(IPARAM(3).NE.7.AND.IPARAM(3).NE.9) IPARAM(3)=8
      IF(IPARAM(4).NE.4)                    IPARAM(4)=3
      IF(IPARAM(5).NE.3)                    IPARAM(5)=4
! B////////////////////////////////////////////////////////////////////B
! Change default for $CRITERION from RMS to BAND


      IF(IPARAM( 6).NE.11.AND.IPARAM( 6).NE.12.AND.IPARAM(6).NE.13)            &
                                              IPARAM( 6)=10
! E////////////////////////////////////////////////////////////////////E
      IF(IPARAM(7).NE.4)                    IPARAM(7)=3
      IF(IPARAM(8).NE.4)                    IPARAM(8)=3
      IF(IPARAM(9).NE.4)                    IPARAM(9)=3
      IF(IPARAM(10).NE.6)                   IPARAM(10)=5
      IF(IPARAM(11).NE.4)                   IPARAM(11)=3
!     INVOKE CONNECTION TABLE OPTION IF SPRINGS ARE REQUESTED.
      IF(IPARAM(7).EQ.4) IPARAM(2)=4
      KDIM=MIN(KDIM,NGRID)
!     INVOKE RMS CRITERION FOR FRONTAL JOB.
      IF(IPARAM(11).EQ.4) IPARAM(6)=12
!     Force CM to be run if $print max is enabled (to get component list)
      if(iparam(10).eq.6.and.iparam(3).eq.8) iparam(3)=9
!     TO PREVENT USE OF A PARTICULAR $ CARD (SUCH AS $SPRING), ENFORCE
!          THE NO OPTION HERE.
!
      RETURN
!
!     $IGNORE G1,G2, ...     (NODES TO IGNORE)
!
  100 CALL READIT(KA(4),MAXI,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      I=IIG
      IIG=IIG+NIP
      IF(IIG.LE.IDIM) GO TO 105
      WRITE(IOU6,102)
  102 FORMAT(/,'Fatal Error.  Too many points on $-card.')
      call finish(6,IER)
      RETURN
  105 DO 110 J=1,NIP
      K=I+J
  110 IGNORE(K)=IP(J)
      RETURN
!
!     $GRID N    (UPPER BOUND ON NUMBER OF GRID POINTS)
!
  120 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.EQ.0) RETURN
      NGRID=MAX(IP(1),10)
      KDIM=MAX(KDIM,NGRID/10)
      KDIM=MIN(KDIM,NGRID)
      GO TO 195
!
!     $START G1,G2, ...     (USER-SELECTED STARTING NODES FOR CM METHOD)
!
  130 CALL READIT(KA(4),MAXI,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      I=ISTA
      ISTA=ISTA+NIP
      IF(ISTA.LE.IDIM) GO TO 135
      WRITE(IOU6,102)
      call finish(6,IER)
      RETURN
  135 DO 140 J=1,NIP
      K=I+J
  140 ISTART(K)=IP(J)
      RETURN
!
!     $DEGREE N     (TO IGNORE NODES OF DEGREE EXCEEDING N)
!
  150 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.GT.0) IGDEG=IP(1)
      RETURN
!
!     $INSERT       OR      $INSERT N
!
  160 CALL INSERT(KA,KT)
      RETURN
!
!     $LINES N     (NO. OF LINES PER PAGE)
!
  180 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      LINES=MAX(IP(1),10)
      RETURN
!
!     $DIMENSION N    (DIMENSION OF SOME SCRATCH ARRAYS)
!
  190 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.LE.0) RETURN
      KDIM=MAX(IP(1),10)
195   KDIM=MIN(KDIM,KORE/4)
      RETURN
!
!     $ADD N     (TO ADD N TO NEW SEQUENCE NUMBERS)
!
  230 CALL READIT(KA(4),1,69,IP,NIP)
      IF(NIP.GT.0) IADD=IP(1)
      RETURN
!
!     $APPEND     CNAME     NCON     IFLD
!
!     USER-DEFINED CONNECTION CARD, WHERE
!
!     CNAME=NAME OF CONNECTION CARD (E.G., CBAR) LEFT-ADJUSTED STARTING
!          IN COLUMN 9
!     NCON=NUMBER OF CONNECTIONS ON CARD (I.E., NODES IN ELEMENT)
!     IFLD=NASTRAN FIELD NUMBER ON PARENT CARD IN WHICH FIRST
!          CONNECTION APPEARS
!     NCON AND IFLD MAY APPEAR ANYWHERE IN COLUMNS 17 - 32 SEPARATED BY
!          ONE OR MORE BLANKS.
!
!     REMARKS  -
!          1. NO LONG-FIELD CONNECTION CARDS MAY BE DEFINED.
!          2. CONNECTIONS MUST BE LISTED CONSECUTIVELY ON PARENT AND
!             CONTINUATION CARDS, IF ANY.
!          3. EACH TIME A $APPEND CARD APPEARS, A NEW ELEMENT IS
!             DEFINED.
!          4. IF CNAME MATCHES AN EXISTING ELEMENT NAME, THE OLD
!             ELEMENT IS REPLACED BY THE NEW ONE.
!          5. THE NUMBER OF ELEMENT CONNECTIONS IS FIXED RATHER
!             THAN VARIABLE.
!          6. IF "$FRONTAL YES" IS SPECIFIED, FIELD 2 MUST CONTAIN THE
!             ELEMENT NUMBER (SINCE ELEMENT NUMBERS ARE NEEDED FOR
!             RESEQUENCING).
!
  240 CALL READIT(KA(17),2,16,IP,NIP)
      IF(NIP.LT.2) GO TO 260
      NCON=IP(1)
      IFLD=IP(2)
      IF(NCON.LT.1) GO TO 260
      IF(IFLD.LT.2) GO TO 260
      IF(IFLD.GT.9) GO TO 260
      REWIND IOU11
      WRITE(IOU11,245) (KA(I),I=9,16)
  245 FORMAT(8A1)
      rewind iou11
      READ(IOU11,250) (KA(I),I=1,3)
  250 FORMAT(A1,A4,A3)
      REWIND IOU11
      IF(KA(1).NE.MA(3)) GO TO 260
      IF(KA(2).EQ.MB(2)) GO TO 260
!     CHECK IF CNAME MATCHES NAME ALREADY IN LIST.
      DO 251 I=19,NTYPE
         IF(I.GE.70.AND.I.LE.78) GO TO 251
         IF(I.EQ.93) GO TO 251
         IF(I.EQ.96) GO TO 251
         IF(I.EQ.97) GO TO 251
         IF(KA(2).EQ.TYPE(I).AND.KA(3).EQ.WYPE(I)) GO TO 252
251   CONTINUE
      NTYPE=NTYPE+1
      I=NTYPE
      GO TO 254
252   WRITE(IOU6,253)
253   FORMAT(10X,'This element replaces another of the same name.')
254   IF(NTYPE.GT.MDIM) GO TO 260
      vype(i)=ka(1)
      TYPE(I)=KA(2)
      WYPE(I)=KA(3)
      ME(I)=10*NCON+IFLD
!     TURN ON SWITCH TO PRINT OUT ELEMENT LIBRARY.
      IPARAM(9)=4
!     SET KDIM SO ENOUGH SPACE IS AVAILABLE.
      KDIM=MAX(KDIM, NCON/4+1)
      IF(KDIM.GT.(KORE/4)) GO TO 260
      I=MDIM-NTYPE
      WRITE(IOU6,255) I
  255 FORMAT(10X,'Space exists to define',I4,' more elements.')
      RETURN
!     ERROR FOUND ON $APPEND CARD.
  260 CONTINUE
      WRITE(IOU6,265)
  265 FORMAT(/,'Fatal Error.  Illegal data on $APPEND card.')
      call finish(6,IER)
      RETURN
!
      END SUBROUTINE DOLLAR

! ##################################################################################################################################
      SUBROUTINE ELTYPE(KA,ITYPE,NCON,IFLD,LOOP,MAXI,LEN,LESSOK,IER)
!
!     DETERMINE BULK DATA CARD ELEMENT TYPE.
!
!     KA     = CONTENTS OF BULK DATA CARD (A1,A4,A3,64A1,A1,A4,A3)
!              (INPUT)
!     ITYPE  = ELEMENT TYPE NUMBER (OUTPUT)
!     NCON   = NUMBER OF CONNECTIONS FOR ELEMENT (OUTPUT)
!              (FOR SOME ELEMENTS  (E.G., CELAS1), THE CONNECTIONS ARE
!              NOT LISTED CONSECUTIVELY, IN WHICH CASE NCON BECOMES THE
!              NUMBER OF FIELDS.  LATER, THE MIDDLE FIELD IS BLANKED
!              OUT.)
!     IFLD   = FIELD NUMBER OF FIRST CONNECTION (OUTPUT)
!     LOOP   = 2 IF TWO ELEMENTS CAN BE DEFINED ON ONE CARD,
!              OTHERWISE 1. (OUTPUT)
!              THE ONLY CARD TYPES FOR WHICH TWO ELEMENTS CAN BE DEFINED
!              ON ONE CARD ARE CDAMP3, CDAMP4, CELAS3, CELAS4, CMASS3,
!              CMASS4, CDAMP4*,CELAS4*, CMASS4*, CRIGDR.
!              CROD(10), CTUBE(11), CVISC(12) were previously included
!              in this list and removed.
!     MAXI   = MAXIMUM NUMBER OF GRID POINTS ON A CONNECTION CARD
!              (=kdim4 from main) (INPUT)
!     LEN    = 1 FOR SHORT FIELD DATA CARDS,
!              2 FOR LONG  FIELD DATA CARDS (OUTPUT)
!     LESSOK = .TRUE. IF IT IS OK FOR THE ELEMENT TO HAVE FEWER THAN
!              THE MAXIMUM NUMBER OF NODES PRESENT ON CONNECTION CARD
!              (E.G., CELAS1 OR CPENTA), OTHERWISE, .FALSE. (OUTPUT)
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IER    ,IFLD   ,IPARAM ,ITYPE  ,J      ,L      ,        &
               LEN    ,LOOP   ,MA     ,MAXI   ,MB     ,MDIM   ,ME     ,        &
               NCON   ,NELEM  ,NTYPE  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),               &
                    NELEM(160),MDIM
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER vype,TYPE,WYPE,KA(70),EQGP
      LOGICAL LESSOK
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with NCON we will use NCON_array
! instead. Needed so Lahey doesn't complain about shape of NCON being
! different than IP (array) in subr READIT
      INTEGER NCON_array(1)
! E////////////////////////////////////////////////////////////////////E
!
!     SEE BLOCK DATA FOR A LISTING OF THE BANDIT ELEMENT LIBRARY.
!     IT CAN BE LISTED AT EXECUTION TIME WITH  $ELEMENT YES  CARD.
!
! B////////////////////////////////////////////////////////////////////B
      DATA EQGP/4HEQGP/
! E////////////////////////////////////////////////////////////////////E
      LOOP=1
      LEN=1
      ITYPE=0
      LESSOK=.FALSE.
      IF(KA(1).NE.MA(3)) GO TO 30
!
!     LOOK FOR CONNECTION CARD.
!
      DO I=4,NTYPE
         IF(KA(2).EQ.TYPE(I).AND.KA(3).EQ.WYPE(I)) GO TO 20
      end do
      RETURN
!
   20 ITYPE=I
      NCON=ME(I)/10
      IFLD=ME(I)-10*NCON
      IF(NCON.eq.0) then
         WRITE(IOU6,24) TYPE(I),TYPE(I)
24       FORMAT(/,'Fatal Error.  A',A4,' card does not precede first!',        &
                A4,' card.')
         call finish(6,IER)
         RETURN
      end if
!
!     SET SPECIAL PARAMETERS LOOP, LEN, AND LESSOK.
!
!     IF(I.GE.10.AND.I.LE.12) LOOP=2
      IF(I.GE.13.AND.I.LE.18) LOOP=2
      IF(I.GE.70.AND.I.LE.72) LOOP=2
      IF(I.EQ.93) LOOP=2
      IF(I.GE.70.AND.I.LE.78) LEN =2
      IF(I.GE.04.AND.I.LE.09) LESSOK=.TRUE.
      IF(I.GE.13.AND.I.LE.18) LESSOK=.TRUE.
      IF(I.GE.70.AND.I.LE.75) LESSOK=.TRUE.
      IF(I.EQ.96.OR.I.EQ.97) LESSOK=.TRUE.
      IF(I.GE.121.AND.I.LE.127) LESSOK=.TRUE.
      RETURN
!
!     LOOK FOR ENDDATA, MPC, MPC*, MPCAX, OR SEQGP.
!
   30 IF(KA(1).EQ.MA( 5).AND.KA(2).EQ.TYPE(1))   ITYPE=1              ! enddata
      IF(KA(1).EQ.MA(13).AND.KA(2).EQ.TYPE(2))   ITYPE=2              ! mpc
      IF(KA(1).EQ.MA(13).AND.KA(2).EQ.TYPE(3))   ITYPE=3              ! mpc*
      IF(KA(1).EQ.MA(13).AND.KA(2).EQ.TYPE(112)) ITYPE=112            ! mpcax
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(129)) ITYPE=129            ! rbar
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(131)) ITYPE=131            ! rbe1
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(133)) ITYPE=133            ! rbe2
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(135)) ITYPE=135            ! rrod
      IF(KA(1).EQ.MA(18).AND.KA(2).EQ.TYPE(137)) ITYPE=137            ! rtrplt
      IF(KA(1).EQ.MA(19).AND.KA(2).EQ.EQGP) GO TO 35                  ! seqgp
      IF(KA(1).EQ.MA( 1)) GO TO 50                                    ! a
      RETURN
!
!     SEQGP CARD FOUND.   IF RESEQUENCING REQUESTED, ABORT JOB.
!
   35 IF(IPARAM(5).EQ.3) RETURN
! B////////////////////////////////////////////////////////////////////B
! Comment out this fatal error, MYSTRAN disallows SEQGP card images when
! it reads the input deck prior to Bandit being called
! E////////////////////////////////////////////////////////////////////E
      RETURN
!
!     LOOK FOR ADUMI CARD AND EXTRACT NCON, THE NUMBER OF GRID POINTS.
!
   50 DO 60 I=58,67
      IF(KA(2).EQ.TYPE(I)) GO TO 70
   60 CONTINUE
      RETURN
! B////////////////////////////////////////////////////////////////////B
70    CALL READIT(KA(4),1,8,NCON_array(1),J)
      NCON = NCON_array(1)
! E////////////////////////////////////////////////////////////////////E
      L = ME(I) - (ME(I)/10) * 10
      ME(I) = 10 * NCON + L
!     THE DEFAULT NUMBER OF NODES FOR CDUMMY IS GIVEN IN BLOCK DATA.
!          IT CAN BE CHANGED WITH AN ADUMMY CARD.
!     THE DEFAULT NUMBER OF NODES FOR CDUMI IS 0.
      IF(NCON.LE.MAXI) RETURN
      WRITE(IOU6,80) (KA(I),I=1,3),MAXI
80    FORMAT(/,'Fatal Error.  Number of grid connections on ',A1,A4,A3,        &
       'exceeds limit of',I7)
      WRITE(IOU6,85)
85    FORMAT('Use $DIM N card, where 4*N exceeds the maximum number'/          &
             'of grid points for any one element.')
      call finish(6,IER)
      RETURN
      END SUBROUTINE ELTYPE

! ##################################################################################################################################
      SUBROUTINE FINISH(KUMP,IER)
!
! B////////////////////////////////////////////////////////////////////B
      USE IOUNT1, ONLY   :  WRT_LOG, SEQ
! E////////////////////////////////////////////////////////////////////B
!     TERMINATE JOB AFTER FATAL ERROR.
!     THE SUBSEQUENT EXECUTION OF NASTRAN IS PREVENTED BY ERASING UNIT 8.
!
!
!     Kump   Reason
!      1     maxgrd or maxdeg exceeded
!      2     $DIM too small for rigid element
!      3     $DIM exceeded in GPS
!      4     EOF encountered
!      5     $DIM exceeded in CM
!      6     all other reasons
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IER    ,IPARAM ,ITIME  ,KUMP   ,NN     ,KDIM   ,        &
               KA     ,NCARD

      REAL     DUM    ,DUMO   ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /BITS/ DUM(8),KDIM
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      COMMON /DOL/ KA(20),DUMO(180)
      GO TO (10,60,30,40,30,20), KUMP
!
! *********************************************************************
! KUMP = 4: EOF encountered
! --------
!
40    WRITE(IOU6,42)
42    FORMAT(/,'Fatal Error.  End-of-file encountered.')
      GO TO 20
!
! *********************************************************************
! KUMP = 3 or 5: Abort job when scratch dimension exceeded in GPS or CM
! -------------
!
   30 CONTINUE
      KDIM=2*KDIM
      KDIM=MIN(KDIM,NN)
      WRITE(IOU6,32) NN,KDIM
   32 FORMAT(/,'Fatal Error.  Scratch array dimension exceeded.',              &
             '  Resubmit job with'/                                            &
             15X,'$GRID',I8,'  and'/15X,'$DIM ',I8)
      IF(KUMP.EQ.5) GO TO 20
!
!     RECOVER SEQGP CARDS FROM IOU9 IF STOP OCCURRED IN GPS AFTER
!     FINISHING CM.
!
      IF(IPARAM(3).NE.9) GO TO 20
      REWIND IOU7
! B////////////////////////////////////////////////////////////////////B
      REWIND SEQ
      READ (SEQ,'(1X,I11)') ITIME
! E////////////////////////////////////////////////////////////////////E
      REWIND IOU9
!     NCARD=NUMBER OF SEQGP CARDS.
      NCARD=(NN-1)/4 + 1
!     COPY SEQGP CARDS TO UNIT 7.
      DO 37 I=1,NCARD
      READ(IOU9,36) KA
   36 FORMAT(20A4)
!     KA IS SCRATCH SPACE HERE.
37    WRITE(IOU7,36) KA
! B////////////////////////////////////////////////////////////////////B
      WRITE(SEQ,36) KA
! E////////////////////////////////////////////////////////////////////E
      REWIND IOU7
! B////////////////////////////////////////////////////////////////////B
      REWIND SEQ
      READ (SEQ,'(1X,I11)') ITIME
! E////////////////////////////////////////////////////////////////////E
      WRITE(IOU6,38)
   38 FORMAT(/,'SEQGP cards generated by CM have been recovered',              &
             ' and placed on bandit.f07')
      GO TO 20
!
! *********************************************************************
! KUMP = 1: Quit since MAXGRD or MAXDEG too small
! --------
!
   10 WRITE(IOU6,15)
   15 FORMAT('Use $GRID N card and/or increase memory in code.')
      go to 20
!
! *********************************************************************
! KUMP = 6: Quit since $DIM too small for rigid element
! --------
!
60    KDIM=2*KDIM
      WRITE(IOU6,62) KDIM
62    FORMAT('Resubmit job with'/5X,'$DIM',I8)
      go to 20
!
! *********************************************************************
! KUMP = 2: Quit for all other reasons than above
! --------

   20 WRITE(IOU8,'(A)') 'Fatal Error.  bandit.f08 deleted.'
      if(iparam(1).eq.3) close(iou8,status='delete')
!     WRITE(IOU6,'(A)') '0End of BANDIT Job.'
!     STOP 13
      IER=KUMP
      RETURN
      END SUBROUTINE FINISH

! ##################################################################################################################################
      SUBROUTINE FIXIT(LIST,NL)
!
!     THIS ROUTINE COMPRESSES OUT ZEROS AND MULTIPLE ENTRIES IN A LIST
!     ORIGINALLY OF LENGTH NL.  A CORRECTED LENGTH NL IS RETURNED TO
!     THE CALLING PROGRAM.
!
      INTEGER LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,I1     ,J      ,NL     ,NL1

! E////////////////////////////////////////////////////////////////////E
!
!     DELETE ZEROS.
!
      CALL ZERO(LIST,NL)
!
!     DELETE DUPLICATE ENTRIES.
!
      IF(NL.LE.1) RETURN
      NL1=NL-1
      DO 20 I=1,NL1
      I1=I+1
      DO 10 J=I1,NL
      IF(LIST(I).NE.LIST(J)) GO TO 10
      LIST(I)=0
      GO TO 20
10    CONTINUE
20    CONTINUE
!
!     DELETE ZEROS AGAIN.
!
      CALL ZERO(LIST,NL)
      RETURN
      END SUBROUTINE FIXIT

! ##################################################################################################################################
      SUBROUTINE FLIP(LIST,N,INV,II3,IER)
!
!     CONVERT $-ARRAY LIST OF LENGTH N FROM ORIGINAL TO INTERNAL NUMBERS
!
      INTEGER II3, INV(2,II3),LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IER    ,N

! E////////////////////////////////////////////////////////////////////E
!     CHECK FOR DUPLICATE AND ZERO ENTRIES AND REDUCE N IF NECESSARY.
      CALL FIXIT(LIST,N)
      IF(N.LE.0) RETURN
      DO 10 I=1,N
   10 LIST(I)=INTERN(LIST(I),INV,II3,IER)
      IF(IER.GT.0) RETURN
      RETURN
      END SUBROUTINE FLIP

! ##################################################################################################################################
      SUBROUTINE FNDIAM(SND1,SND2,NDSTK,NR,NDEG,LVL,LVLS1,LVLS2,               &
        IWK,IDFLT,NDLST,IDIM,IER)
!
!  FNDIAM IS THE CONTROL PROCEDURE FOR FINDING THE PSEUDO-DIAMETER OF
!  NDSTK AS WELL AS THE LEVEL STRUCTURE FROM EACH END
!
!  SND1-        ON INPUT THIS IS THE NODE NUMBER OF THE FIRST
!               ATTEMPT AT FINDING A DIAMETER.  ON OUTPUT IT
!               CONTAINS THE ACTUAL NUMBER USED.
!  SND2-        ON OUTPUT CONTAINS OTHER END OF DIAMETER
!  LVLS1-       ARRAY CONTAINING LEVEL STRUCTURE WITH SND1 AS ROOT
!  LVLS2-       ARRAY CONTAINING LEVEL STRUCTURE WITH SND2 AS ROOT
!  IDFLT-       FLAG USED IN PICKING FINAL LEVEL STRUCTURE, SET
!               =1 IF WIDTH OF LVLS1 @ WIDTH OF LVLS2, OTHERWISE =2
!  LVL,IWK-     WORKING STORAGE
!
      INTEGER FLAG,NR,SND,SND1,SND2
      COMMON /GRA/ N,IDPTH,DUMG
      DIMENSION NDSTK(NR,*),NDEG(*),LVL(*),LVLS1(*),LVLS2(*),IWK(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IDFLT  ,IDIM   ,IDPTH  ,IER    ,IWK    ,LVL    ,        &
               LVLBOT ,LVLN   ,LVLS1  ,LVLS2  ,LVLWTH ,MAXLW  ,mtw1   ,        &
               MTW2   ,N      ,NDEG   ,NDSTK  ,NDXL   ,ndxn

      REAL    DUMG   

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER NDLST(IDIM)
!     DIMENSION OF NDLST IS THE MAX NUMBER OF NODES IN LAST LEVEL.
      mtw1=0
      ndxn=0
      FLAG=0
      MTW2=N
      SND=SND1
!  ZERO LVL TO INDICATE ALL NODES ARE AVAILABLE TO TREE
   20 DO 25 I=1,N
        LVL(I)=0
   25 CONTINUE
      LVLN=1
!  DROP A TREE FROM SND
      CALL TREE(SND,NDSTK,NR,LVL,IWK,NDEG,LVLWTH,LVLBOT,LVLN,MAXLW,MTW2)
      IF(FLAG.GE.1) GO TO 110
      FLAG=1
   70 IDPTH=LVLN-1
      MTW1=MAXLW
!  COPY LEVEL STRUCTURE INTO LVLS1
      DO 75 I=1,N
        LVLS1(I)=LVL(I)
   75 CONTINUE
      NDXN=1
      NDXL=0
      MTW2=N
!  SORT LAST LEVEL BY DEGREE  AND STORE IN NDLST
      CALL SORTDG(NDLST,IWK(LVLBOT),NDXL,LVLWTH,NDEG)
      IF(NDXL.LE.IDIM) GO TO 100
!     DIMENSION EXCEEDED; STOP JOB.
      CALL FINISH(3,IER)
      RETURN
  100 CONTINUE
      SND=NDLST(1)
      GO TO 20
  110 IF(IDPTH.GE.LVLN-1) GO TO 120
!  START AGAIN WITH NEW STARTING NODE
      SND1=SND
      GO TO 70
  120 IF(MAXLW.GE.MTW2) GO TO 130
      MTW2=MAXLW
      SND2=SND
!  STORE NARROWEST REVERSE LEVEL STRUCTURE IN LVLS2
      DO 125 I=1,N
        LVLS2(I)=LVL(I)
  125 CONTINUE
  130 IF(NDXN.EQ.NDXL) GO TO 140
!  TRY NEXT NODE IN NDLST
      NDXN=NDXN+1
      SND=NDLST(NDXN)
      GO TO 20
  140 IDFLT=1
      IF(MTW2.LE.MTW1) IDFLT=2
      RETURN
      END SUBROUTINE FNDIAM

! ##################################################################################################################################
      SUBROUTINE FRONT(KG,ILD,NN,EL,MEM,IER)
!
!     GENERATE ELEMENT SEQUENCE FOR FRONTAL SOLVERS.
!
!     G.C. EVERSTINE, NSWCCD 204, 12/3/90 (revised 5/22/96)
!
!     KG      = ARRAY FOR STORING THE CONNECTIONS FOR AN ELEMENT
!               (SCRATCH)
!     ILD(I)  = NEW LABEL FOR NODE WITH ORIGINAL INTERNAL LABEL I
!               (INPUT)
!     NN      = NUMBER OF GRID POINTS IN MODEL (INPUT)
!     EL(1,I) = ORIGINAL ELEMENT ID FOR ELEMENT I
!     EL(2,I) = SMALLEST NODE NUMBER IN NEW SEQUENCE IN ELEMENT I
!               (EL IS A SCRATCH ARRAY)
!     MEM     = MEMORY AVAILABLE FOR EL ARRAY (INPUT)
!
      INTEGER KG(*),ILD(*),EL(2,*),EID,EMOD
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IEL    ,IER    ,J      ,J1     ,K      ,KFLAG  ,        &
               KMIN   ,L      ,LOC    ,MEM    ,mflag  ,NBW    ,NCM    ,        &
               NEED   ,NEL    ,NEL1   ,NEQ    ,NEQR   ,NLINK  ,NN     ,        &
               NP     ,NPT    ,NZERO  ,OBW    ,OP     

! E////////////////////////////////////////////////////////////////////E
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!
      REWIND IOU13
!
!     CHECK MEMORY
!
!     I=MAX(NEL,NN)
      I=MAX(NEL,NN,mem/4)
      NEED=4*I
      IF(MEM.LT.NEED) then
         WRITE(IOU6,5) NEED,MEM
5        FORMAT(/,'Fatal Error.  Insufficient memory available',               &
                ' for element sort for frontal solver.'/                       &
                12X,I10,' words needed,',I10,' words available.')
         CALL FINISH(1,IER)
         RETURN
      end if
!
!     INITIALIZE THE SCATTER SORT ARRAY
!
      EMOD=2*I-IFIX(2.3715*SQRT(FLOAT(I)))
      DO I=1,EMOD
         EL(1,I)=0
         EL(2,I)=0
      end do
!
!     READ ELEMENT CONNECTIONS IN ORIGINAL INTERNAL SORT.
!     FIND LOWEST NUMBERED NODE IN NEW SEQUENCE.
!     LOOK FOR DUPLICATE ELEMENT NUMBERS.
!     FILL UP SCATTER SORT ARRAY (EL) SO AS TO PERFORM ROUGH SORT OF
!        ELEMENTS IN ORDER OF LOWEST NUMBERED NODE IN NEW SEQUENCE.
!        (FOR MOST STRUCTURES, THE ELEMENTS WILL BE IN PROPER SORT
!        AFTER THIS STEP.)
!
      DO 50 IEL=1,NEL
         READ(IOU13) EID,NPT,(KG(I),I=1,NPT)
         KMIN=NN
         DO I=1,NPT
            KMIN=MIN(ILD(KG(I)),KMIN)
         end do
         IF(KMIN.EQ.0) GO TO 180
!        Add more space for tet meshes, which have lots of elements
!        for each grid point.  Ideally, the memory available would
!        be at least 40 times the number of grids.
         kmin=10*kmin
         LOC=KMIN-1
30       LOC=MOD(LOC,EMOD)+1
         IF(EL(1,LOC).EQ.0) GO TO 40
         IF(EL(1,LOC).EQ.EID) GO TO 160
         GO TO 30
40       EL(1,LOC)=EID
         EL(2,LOC)=KMIN
50    CONTINUE
!
!     SQUEEZE OUT THE ZEROS IN THE SORT ARRAY.
!
      K=0
      DO 60 I=1,EMOD
         IF(EL(1,I).EQ.0) GO TO 60
         K=K+1
         EL(1,K)=EL(1,I)
         EL(2,K)=EL(2,I)
60    CONTINUE
      IF(K.NE.NEL) then
         WRITE(IOU6,'(A)') ' Logic error in FRONT'
         call finish(6,IER)
         RETURN
      end if
!
!     COMPLETE SORT OF ELEMENTS.  BY NOW, ELEMENTS ARE PROBABLY
!     ALREADY IN SORT.  IF NOT, THIS BUBBLE SORT WILL FINISH IT.
!
      IF(NEL.LE.1) RETURN
      NEL1=NEL-1
      DO 100 I=1,NEL1
         K=NEL-I
         KFLAG=0
         DO 90 J=1,K
            J1=J+1
            IF(EL(2,J).LE.EL(2,J1)) GO TO 90
            KFLAG=1
            L=EL(1,J)
            EL(1,J)=EL(1,J1)
            EL(1,J1)=L
            L=EL(2,J)
            EL(2,J)=EL(2,J1)
            EL(2,J1)=L
90       CONTINUE
         IF(KFLAG.EQ.0) GO TO 110
100   CONTINUE
!
!     WRITE OUT ELEMENT ORDER ON OUTPUT FILE AND UNIT 14.
!        mflag = 0    write to units 6 and 14
!                1    write to unit 14 only
!
110   mflag=1
      if(mflag.eq.0) then
         WRITE(IOU6,120)
120      FORMAT(/,'Element ordering for frontal solution (bandit.f14)')
         DO I=1,NEL,5
            K=MIN(I+4,NEL)
            WRITE(IOU6,140) I,(EL(1,J),J=I,K)
140         FORMAT(I10,' =',5I10)
         end do
      else
         write(iou6,145)
145      format('See bandit.f14 for element ordering for',                     &
                ' frontal solution.')
      end if
      OPEN(IOU14,FILE='bandit.f14',FORM='FORMATTED',STATUS='replace')
      REWIND IOU14
      WRITE(IOU14,150) (EL(1,J),J=1,NEL)
150   FORMAT(5I10)
      close(iou14)
!
      RETURN
!
!     Fatal error messages.
!
160   WRITE(IOU6,170) EID
170   FORMAT(/,'Fatal Error.  Duplicate element ID',i10)
      call finish(6,IER)
      RETURN
!
180   WRITE(IOU6,190) EID
190   FORMAT(/,'Fatal Error.  Element',I9,                                     &
             ' has a zero grid point connection.')
      call finish(6,IER)
      RETURN
!
      END SUBROUTINE FRONT

! ##################################################################################################################################
      SUBROUTINE GIBSTK(NDSTK,NR,IOLD,RENUM,NDEG,LVL,LVLS1,LVLS2,CCSTOR,       &
       IBW2,IPF2,JUMP,ICRIT,NHIGH,NLOW,NACUM,SIZE,STPT,IDIM,IER)
!
!  GIBBSTOCK USES GRAPH THEORETICAL METHODS TO PRODUCE A PERMUTATION
!  OF AN INPUT ARRAY WHICH REDUCES ITS BANDWIDTH.
!  PROGRAMMED BY H.L.CRANE JR.  COLLEGE OF WILLIAM + MARY   3/74.
!
!      MODIFIED BY G. C. EVERSTINE, DTRC.
!
!  THE FOLLOWING INPUT PARAMETERS ARE REQUIRED--NDSTK,NR,N,IDEG,IOLD
!
!  THESE INTEGER ARRAYS MUST BE DIMENSIONED IN THE CALLING PROGRAM--
!  NDSTK(NR,D1),RENUM(D2+1),NDEG(D2),IOLD(D2),LVL(D2),LVLS1(D2),
!  LVLS2(D2),CCSTOR(D2)   WHERE D1 .GE. MAX DEGREE OF ANY NODE AND
!  D2 AND NR ARE .GE. THE TOTAL NUMBER OF NODES IN THE GRAPH.
!
!  EXPLANATION OF PARAMETERS--
!  NDSTK-       ADJACENCY ARRAY REPRESENTING GRAPH TO BE PROCESSED
!               NDSTK(I,J)=NODE NUMBER OF JTH CONNECTION TO NODE
!               NUMBER I.  A CONNECTION OF A NODE TO ITSELF IS NOT
!               LISTED.  EXTRA POSITIONS MUST HAVE ZERO FILL.
!  NR-          ROW DIMENSION ASSIGNED NDSTK IN CALLING PROGRAM
!  IOLD(I)-     RENUMBERING OF ITH NODE BEFORE GIBBSTOCK PROCESSING
!               IF NO RENUMBERING EXISTS THEN ILD(1)=1,ILD(2)=2, ETC.
!  N-           NUMBER OF NODES IN GRAPH BEING PROCESSED
!  IDEG-        MAX DEGREE OF ANY NODE IN GRAPH BEING PROCESSED
!   JUMP IS SET TO 0 IF EITHER CRITERION IS REDUCED.
!     ICRIT=RESEQUENCING CRITERION
!          1=RMS WAVEFRONT
!          2=BANDWIDTH
!          3=PROFILE
!          4=WAVEFRONT (MAX)
!  ON OUTPUT THESE VARIABLES CONTAIN THE FOLLOWING INFORMATION--
!  RENUM(I)-    THE NEW NUMBER FOR THE ITH NODE
!  NDEG(I)-     THE DEGREE OF THE ITH NODE
!  IDPTH-       NUMBER OF LEVELS IN GIBBSTOCK LEVEL STRUCTURE
!  IBW2-        THE BANDWITH AFTER RENUMBERING
!  IPF2-        THE PROFILE AFTER RENUMBERING
!  THE FOLLOWING ONLY HAVE MEANING IF THE GRAPH WAS ALL ONE COMPONENT
!  LVL(I)-      INDEX INTO LVLS1 TO THE FIRST NODE IN LEVEL I
!               LVL(I+1)-LVL(I)= NUMBER OF NODES IN ITH LEVEL
!  LVLS1-       LEVEL STRUCTURE CHOSEN BY GIBBSTOCK
!  LVLS2(I)-    THE LEVEL ASSIGNED TO NODE I BY GIBBSTOCK
!
!    THE FOLLOWING SUBROUTINES WERE WRITTEN BY N. GIBBS, W. POOLE,
!    P. STOCKMEYER, AND H. CRANE OF THE COLLEGE OF WILLIAM AND MARY - -
!     DGREE,FNDIAM,GIBSTK,NUMBER,PIKLVL,RSETUP,SORTDG,SORT2,TREE.
!
! B////////////////////////////////////////////////////////////////////B
! Can't define SORT2 here. It is an integer function that is part of
! this BANDIT_MODULE module. In the version of BANDIT_SUBRS that is not
! a module but a collection of subroutines in a file, then it must be
! defined here (or in an EXTERNAL command)
      INTEGER STNODE,RVNODE,RENUM,XC,STNUM,SBNUM
! E////////////////////////////////////////////////////////////////////E
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IBW1   ,IBW2   ,ICRIT  ,IDEG   ,IDFLT  ,IDIM   ,        &
               IDPTH  ,IER    ,IPARAM ,IPF1   ,IPF2   ,ISDIR  ,IWALL  ,        &
               JUMP   ,                                                        &
               LOWDG  ,LROOT  ,LVL    ,LVLBOT ,LVLN   ,LVLS1  ,LVLS2  ,        &
               LVLWTH ,                                                        &
               MAXB   ,MAXB0  ,MAXLW  ,MAXW0  ,MAXW1  ,MAXWA  ,MAXWB  ,        &
               MM     ,                                                        &
               N      ,NBW    ,NCM    ,NDEG   ,NDSTK  ,NFLG   ,NN     ,        &
               NP     ,NR     ,NUM    ,NZERO  

      REAL     AVERW0 ,AVERWB ,BRMS0  ,BRMS1  ,BRMSA  ,BRMSB  ,                &
               CRIT1  ,CRIT2  ,DUMD   ,DUMS   ,RMS0   ,RMS1   ,                &
               RMSA   ,RMSB   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,DUMD(4)
!     OLD AND NEW MAX AND RMS WAVEFRONT FOR ENTIRE PROBLEM,
!          NOT JUST GIBSTK.
      COMMON /W/ MAXW0,RMS0,MAXW1,RMS1,BRMS0,BRMS1
      COMMON /B/ IPARAM(20)
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      COMMON /GRA/ N,IDPTH,IDEG
      INTEGER NHIGH(IDIM),NLOW(IDIM),NACUM(IDIM),SIZE(*),STPT(*)
!     SIZE AND STPT HAVE DIMENSION IDIM/2 AND SHOULD BE CONTIGUOUS IN
!                  CORE WITH SIZE FIRST.
!     XC=NUMBER OF SUB-COMPONENTS RESULTING AFTER REMOVING DIAMETER
!        FROM ONE COMPONENT OF ORIGINAL GRAPH.
      DIMENSION NDSTK(NR,*),LVL(*),LVLS1(*),LVLS2(*),RENUM(*),NDEG(*)
      INTEGER CCSTOR(*),IOLD(*),OBW,OP,XCMAX,SUMW0,SUMWB
      REAL IM1,IM2
      CALL TIMER(TA,IWALL,0,IOU6)
      XCMAX=IDIM/2
      NCM=0
      N=NN
      IBW2=0
      IPF2=0
!  SET RENUM(I)=0 FOR ALL I TO INDICATE NODE I IS UNNUMBERED
      DO 10 I=1,N
        RENUM(I)=0
   10 CONTINUE
!   COMPUTE DEGREE OF EACH NODE AND ORIGINAL B AND P.
      CALL DGREE(NDSTK,NR,NDEG,IOLD,IBW1,IPF1)
!   COMPUTE ORIGINAL ACTIVE COLUMN DATA IF NOT ALREADY SUPPLIED BY
!        CUTHILL.
      IF(IPARAM(3)-9) 20,15,20
   15 MAXWA=MAXW1
      RMSA=RMS1
      BRMSA=BRMS1
      GO TO 25
   20 CONTINUE
      CALL WAVEY(NDSTK,NR,IOLD,LVL,0,LVLS2,LVLS1,MAXB0,MAXW0,AVERW0,           &
       SUMW0,RMS0,BRMS0)
      MAXWA=MAXW0
      RMSA=RMS0
      BRMSA=BRMS0
      WRITE(IOU6,29)
   29 FORMAT(/,'Before resequencing:')
      WRITE(IOU6,51) MAXB0,SUMW0,MAXW0,AVERW0,RMS0,BRMS0
   51 FORMAT(5X,'Bandwidth',I18/5X,'Profile',I20/                              &
             5X,'Max Wavefront',I14/5X,'Avg Wavefront',F14.3/                  &
             5X,'RMS Wavefront',F14.3/5X,'RMS Bandwidth',F14.3)
   25 CONTINUE
!  SBNUM= LOW END OF AVAILABLE NUMBERS FOR RENUMBERING
!  STNUM= HIGH END OF AVAILABLE NUMBERS FOR RENUMBERING
      SBNUM=1
      STNUM=N
!  NUMBER THE NODES OF DEGREE ZERO
      DO 40 I=1,N
        IF(NDEG(I).GT.0) GO TO 40
        RENUM(I)=STNUM
        STNUM=STNUM-1
   40 CONTINUE
!   NODES OF ZERO DEGREE APPEAR LAST IN NEW SEQUENCE.
      NZERO=N-STNUM
      NCM=NZERO
!  FIND AN UNNUMBERED NODE OF MIN DEGREE TO START ON
   50 LOWDG=IDEG+1
      NCM=NCM + 1
      NFLG=1
      ISDIR=1
      DO 70 I=1,N
        IF(NDEG(I).GE.LOWDG) GO TO 70
        IF(RENUM(I).GT.0) GO TO 70
        LOWDG=NDEG(I)
        STNODE=I
   70 CONTINUE
!  FIND PSEUDO-DIAMETER AND ASSOCIATED LEVEL STRUCTURES.
!  STNODE AND RVNODE ARE THE ENDS OF THE DIAM AND LVLS1 AND LVLS2
!  ARE THE RESPECTIVE LEVEL STRUCTURES.
      CALL FNDIAM(STNODE,RVNODE,NDSTK,NR,NDEG,LVL,LVLS1,LVLS2,CCSTOR,          &
       IDFLT,SIZE,IDIM,IER)
      IF(IER.GT.0) RETURN
      IF(NDEG(STNODE).LE.NDEG(RVNODE)) GO TO 75
!  NFLG INDICATES THE END TO BEGIN NUMBERING ON
      NFLG=-1
      STNODE=RVNODE
   75 CALL RSETUP(LVL,LVLS1,LVLS2,NACUM,IDIM,IER)
      IF(IER.GT.0) RETURN
!  FIND ALL THE CONNECTED COMPONENTS  (XC COUNTS THEM)
      XC=0
      LROOT=1
      LVLN=1
      DO  80 I=1,N
        IF(LVL(I).NE.0) GO TO 80
        XC=XC+1
        IF(XC.LE.XCMAX) GO TO 85
!     DIMENSION EXCEEDED  . . .  STOP JOB.
      CALL FINISH(3,IER)
      RETURN
   85   CONTINUE
        STPT(XC)=LROOT
        CALL TREE(I,NDSTK,NR,LVL,CCSTOR,NDEG,LVLWTH,LVLBOT,LVLN,MAXLW,N)
        SIZE(XC)=LVLBOT+LVLWTH-LROOT
        LROOT=LVLBOT+LVLWTH
        LVLN=LROOT
   80 CONTINUE
      IF(SORT2(XC,SIZE,STPT).EQ.0) GO TO 90
      CALL PIKLVL(LVLS1,LVLS2,CCSTOR,IDFLT,ISDIR,XC,NHIGH,NLOW,                &
        NACUM,SIZE,STPT)
!  ON RETURN FROM PIKLVL, ISDIR INDICATES THE DIRECTION THE LARGEST
!  COMPONENT FELL.  ISDIR IS MODIFIED NOW TO INDICATE THE NUMBERING
!  DIRECTION.  NUM IS SET TO THE PROPER VALUE FOR THIS DIRECTION.
   90 ISDIR=ISDIR*NFLG
      NUM=SBNUM
      IF(ISDIR.LT.0) NUM=STNUM
!
      CALL NUMBR(STNODE,NUM,NDSTK,LVLS2,NDEG,RENUM,LVLS1,LVL,NR,NFLG,          &
        IBW2,IPF2,CCSTOR,ISDIR,NHIGH,NLOW,NACUM,SIZE,IDIM,IER)
      IF(IER.GT.0) RETURN
!
!  UPDATE STNUM OR SBNUM AFTER NUMBERING
      IF(ISDIR.LT.0) STNUM=NUM
      IF(ISDIR.GT.0) SBNUM=NUM
      IF(SBNUM.LE.STNUM) GO TO 50
!
!  COMPUTE THE NEW BANDWIDTH, PROFILE, AND WAVEFRONT.
!
      CALL WAVEY(NDSTK,NR,RENUM,LVL,0,LVLS2,LVLS1,MAXB,MAXWB,AVERWB,           &
         SUMWB,RMSB,BRMSB)
!
      IBW2=MAXB
      IPF2=SUMWB
      WRITE(IOU6,705)
  705 FORMAT(/,'After resequencing by Gibbs-Poole-Stockmeyer (GPS):')
      WRITE(IOU6,51) MAXB,SUMWB,MAXWB,AVERWB,RMSB,BRMSB
!
!     CHECK NEW NUMBERING AGAINST OLD NUMBERING.
!
      GO TO (130,135,140,145), ICRIT
  130 IM1=RMSA
      IM2=IPF1
      CRIT1=RMSB
      CRIT2=IPF2
      GO TO 92
  135 IM1=IBW1
      IM2=IPF1
      CRIT1=IBW2
      CRIT2=IPF2
      GO TO 92
  140 IM1=IPF1
      IM2=IBW1
      CRIT1=IPF2
      CRIT2=IBW2
      GO TO 92
  145 IM1=MAXWA
      IM2=RMSA
      CRIT1=MAXWB
      CRIT2=RMSB
      GO TO 92
   92 CONTINUE
      IF(CRIT1-IM1) 110,94,97
   94 IF(CRIT2.LT.IM2) GO TO 110
!
   97 CONTINUE
!  IF ORIGINAL NUMBERING IS BETTER THAN NEW ONE, SET UP TO RETURN IT
      DO 100 I=1,N
        RENUM(I)=IOLD(I)
  100 CONTINUE
      IBW2=IBW1
      IPF2=IPF1
      MAXWB=MAXWA
      RMSB=RMSA
      BRMSB=BRMSA
      GO TO 112
!
!    EQUATE CORRESPONDING GPS AND BANDIT VARIABLES.
!
  110 CONTINUE
      JUMP=0
  112 CONTINUE
      IF(IPARAM(3).NE.8) GO TO 115
      OBW=IBW1
      OP=IPF1
  115 NBW=IBW2
      NP=IPF2
      MAXW1=MAXWB
      RMS1=RMSB
      BRMS1=BRMSB
      CALL TIMER(TB,IWALL,0,IOU6)
      TB=TB-TA
!     IF(TB.GT.1.E-5) WRITE(IOU6,610) TB
      WRITE(IOU6,610) TB
  610 FORMAT(5X,'CP time',F20.3)
      RETURN
      END SUBROUTINE GIBSTK

! ##################################################################################################################################
      SUBROUTINE GRID(IER,IOU6)
!
!     PARTITION OPEN CORE AND COMPUTE PROBLEM SIZE LIMITS.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IBYTE  ,IER    ,IOU6   ,IPASS  ,KDIM   ,KMOD   ,KOR    ,        &
               MAXDEG ,MAXGRD ,NBITIN ,NBYTE  ,NGRID  ,NW

      REAL     DUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /BITS/ NBITIN,KOR,DUM,NGRID,IPASS,NW,NBYTE,IBYTE,KDIM
      COMMON /A/ MAXGRD,MAXDEG,KMOD
!
      NW=1
!PACK NW=2
!     NW = INTEGER PACKING DENSITY (NUMBER OF INTEGERS PER WORD)
!     NBITIN IS USED ON CDC.  NBYTE AND IBYTE ARE USED ON PACKED CRAY.
      NBITIN=30
      NBYTE=8/NW
      IBYTE=9-NBYTE
      MAXGRD=NGRID+NW-1
      MAXGRD=MAXGRD-MOD(MAXGRD,NW)
      MAXDEG=(KOR-8*MAXGRD-3)/(MAXGRD/NW+1)
      MAXDEG=MIN(MAXDEG,MAXGRD-1)
      IF(MAXDEG.GT.0) RETURN
!     MAXDEG NEGATIVE  - - -   FATAL ERROR.
      CALL COREKO
      WRITE(IOU6,30)
   30 FORMAT(/,'Fatal Error.  Insufficient memory.')
      CALL FINISH(1,IER)
      RETURN
      END SUBROUTINE GRID

! ##################################################################################################################################
      FUNCTION IDIST(NS,ML,MAXLEV,IG,II1,IC,IDEG,IDIS,IW,ICC)
!
!     THIS FUNCTION HAS AS ITS VALUE THE MAXIMUM DISTANCE OF ANY NODE
!          IN COMPONENT IC(NS) FROM THE NODE NS.
!     THE DISTANCE OF EACH NODE IN THIS COMPONENT IS STORED IN THE
!          ARRAY IDIS.
!     THE MAXIMUM NUMBER OF NODES AT THE SAME DISTANCE FROM NS IS
!          STORED IN ML.
!
!     INPUT:  IG,IC,IDEG,ICC,NS,MAXLEV
!     OUTPUT: IDIS,IW,ML
!
      INTEGER II1
      DIMENSION IG(II1,*),IC(*),IDEG(*),IDIS(*),IW(*),ICC(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IA     ,IC     ,ICC    ,ICN    ,IDEG   ,IDIS   ,        &
               IDIST  ,IG     ,II     ,IW     ,                                &
               K      ,KI     ,KMOD   ,KO     ,                                &
               L      ,LL     ,                                                &
               MAXDEG ,MAXGRD ,MAXLEV ,ML     ,                                &
               N      ,NBITIN ,NN     ,NNC    ,NS

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      ICN=IC(NS)
      NNC=ICC(ICN+1)-ICC(ICN)
      DO 50 I=1,NN
      IF(IC(I)-IC(NS)) 50,40,50
40    IDIS(I)=0
   50 CONTINUE
      LL=1
      L=0
      KI=0
      KO=1
      ML=0
      IW(1)=NS
      IDIS(NS)=-1
  130 KI=KI+1
      IF(KI-LL)135,132,135
  132 L=L+1
      LL=KO+1
      K=KO-KI+1
      IF(K-ML) 135,135,133
  133 ML=K
      IF(ML-MAXLEV) 135,135,220
135    II=IW(KI)
      N=IDEG(II)
      IF(N)140,215,140
  140 DO 200 I=1,N
      IA=IG(II,I)
!PACK IA=IUNPK(IG,MAXGRD*(I-1)+II,NBITIN)
      IF(IDIS(IA))200,150,200
150   IDIS(IA)=L
      KO=KO+1
      IW(KO)=IA
  200 CONTINUE
      IF(KO-NNC)130,205,205
  205 IDIST=L
      IDIS(NS)=0
      K=KO-LL+1
      IF(K-ML) 206,206,207
  207 ML=K
  206 CONTINUE
      RETURN
  215 L=0
      GO TO 205
  220 IDIST=1
      RETURN
      END FUNCTION IDIST

! ##################################################################################################################################
      SUBROUTINE IGNOR(IG,II1,INV,II3,LIST,N,IER)
!
!     SET UP LIST OF POINTS TO IGNORE IN LIST ARRAY OF LENGTH N.
!
! UPON ENTRY, LIST ALREADY CONTAINS SET OF MPC DEPENDENT NODES, IF ANY.
      INTEGER II1, II3
      DIMENSION IG(II1,*),INV(2,II3),LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDIM   ,IER    ,IFIR   ,IG     ,IGDEG  ,IGNORE ,        &
               IIG    ,INV    ,ISTA   ,ISTART ,                                &
               J      ,KMOD   ,LIST   ,                                        &
               MAXDEG ,MAXGRD ,MM     ,N      ,NBITIN ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
!
!     ADD IGNORE POINTS TO LIST.
!
      IF(IIG.LE.0) GO TO 20
      CALL FLIP(IGNORE,IIG,INV,II3,IER)
      IF(IER.GT.0) RETURN
      IF(IIG.LE.0) GO TO 20
      DO 10 I=1,IIG
      J=IGNORE(I)
   10 LIST(J)=J
!
!     ADD POINTS OF DEGREE GT IGDEG TO LIST.
!
   20 IF(IGDEG.LE.0) GO TO 40
      IF(IGDEG.GE.MM) GO TO 40
      DO 30 I=1,NN
      IF(IG(I,IGDEG+1).GT.0) LIST(I)=I
!PACK IF(IUNPK(IG,MAXGRD*IGDEG+I,NBITIN).GT.0) LIST(I)=I
   30 CONTINUE
!
!     COMPRESS OUT ZEROES FROM LIST.
!
   40 N=NN
      CALL ZERO(LIST,N)
      RETURN
      END SUBROUTINE IGNOR

! ##################################################################################################################################
      SUBROUTINE INSERT(KA,KT)
!
!     READ AND INSERT CONTENTS OF INSERT FILE (bandit.ins).
!
!     $INSERT CARD MUST APPEAR BEFORE BEGIN BULK TO BE RECOGNIZED BY
!     BANDIT.  FORMAT IS  $INSERT NCARD , WHERE THE PARAMETER NCARD IS
!     THE NUMBER OF CARDS TO BE INSERTED.  IF NCARD IS MISSING OR ZERO,
!     THE ENTIRE FILE 10 WILL BE INSERTED.  ONLY THE FIRST $INSERT CARD
!     IN DECK WILL BE HONORED.
!
!     KT = RUNNING COUNTER ON TOTAL NUMBER OF CASE CONTROL CARDS
!
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with NCARD we will use NCARD_array
! instead. Needed so Lahey doesn't complain about shape of NCARD being
! different than IP (array) in subr READIT
      INTEGER NCARD_array(1)
! E////////////////////////////////////////////////////////////////////E
      INTEGER KA(20)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,II     ,IN     ,KT     ,L      ,NCARD

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      DATA IN/0/
      IF(IN.EQ.1) RETURN
      IN=1
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(4),1,69,NCARD_array(1),I)
      NCARD = NCARD_array(1)
! E////////////////////////////////////////////////////////////////////E
      IF(NCARD.LE.0) NCARD=999999
!
      OPEN(IOU10,FILE='bandit.ins',FORM='FORMATTED',STATUS='UNKNOWN')
      REWIND IOU10
      DO 40 II=1,NCARD
      READ(IOU10,20,END=50) KA
20    FORMAT(20A4)
      KT=KT+1
      L=LENCAS(KA,20)
      WRITE(IOU6,30) KT,(KA(I),I=1,L)
30    FORMAT(I8,'= ',20A4)
      WRITE(IOU8,20) (KA(I),I=1,L)
40    CONTINUE
50    RETURN
      END SUBROUTINE INSERT

! ##################################################################################################################################
      FUNCTION INTERN(IGRID,INV,II3,IER)
!
!     THIS FUNCTION HAS AS ITS VALUE THE INTERNAL NODE LABEL ASSIGNED
!     BY BANDIT TO ORIGINAL GRID POINT IGRID.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  IER    ,IGRID  ,II3    ,INTERN ,KMOD   ,LOC    ,MAXDEG ,        &
               MAXGRD

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER INV(2,II3)
      intern=0
      IF(IGRID.LE.0) GO TO 20
      LOC=IGRID-1
   10 LOC=MOD(LOC,KMOD)+1
      IF(INV(1,LOC).EQ.0) GO TO 20
      IF(INV(1,LOC).NE.IGRID) GO TO 10
      INTERN=INV(2,LOC)
      RETURN
!
!     ABORT JOB DUE TO REFERENCE TO NON-EXISTENT GRID POINT.
!
   20 WRITE(IOU6,30) IGRID
   30 FORMAT(/,'Fatal Error.  Grid point',I9,' on $ card not found.')
      call finish(6,IER)
      RETURN
      END FUNCTION INTERN

! ##################################################################################################################################
      FUNCTION KOMPNT(IG,II1,IC,IDEG,IW,ICC)
!
!     THIS FUNCTION HAS AS ITS VALUE THE NUMBER OF COMPONENTS STORED
!     IN THE CONNECTION ARRAY IG.
!     ALSO, IC AND ICC ARE SET UP.
!     IC(I)=COMPONENT INDEX FOR NODE I
!     ICC(I)=THE STARTING POSITION TO BE USED FOR LABELS IN COMPONENT I
!     THUS, ICC(I+1)-ICC(I)= THE NUMBER OF NODES IN COMPONENT I
!
      INTEGER II1
      DIMENSION IG(II1,*),IC(*),IDEG(*),IW(*),ICC(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IA     ,IC     ,ICC    ,IDEG   ,IG     ,II     ,        &
               IS     ,IW     ,                                                &
               KI     ,KMOD   ,KO     ,KOMPNT ,                                &
               MAXDEG ,MAXGRD ,MM     ,                                        &
               N      ,NBITIN ,NC     ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      kompnt=0
      DO I=1,NN
         ICC(I)=0
         IC(I)=0
      end do
      NC=0
      ICC(1)=1
  105 DO 110 I=1,NN
      IF(IC(I)) 110,120,110
  110 KOMPNT=NC
      RETURN
  120 NC=NC+1
      KI=0
      KO=1
      IW(1)=I
      IC(I)=NC
      IF(NC-1)130,125,125
125   IS=ICC(NC)+1
      ICC(NC+1)=IS
  130 KI=KI+1
      II=IW(KI)
      N=IDEG(II)
      IF(N)140,105,140
  140 DO 200 I=1,N
      IA=IG(II,I)
!PACK IA=IUNPK(IG,MAXGRD*(I-1)+II,NBITIN)
      IF(IC(IA)) 200,150,200
150   IC(IA)=NC
      KO=KO+1
      IW(KO)=IA
      IS=ICC(NC+1)+1
      ICC(NC+1)=IS
  200 CONTINUE
      IF(KO-KI)105,105,130
      END FUNCTION KOMPNT

! ##################################################################################################################################
      SUBROUTINE LEFT(KA,ITYPE)
!
!     LEFT-ADJUST BULK DATA CARD MNEUMONIC AND RETURN IN KA THE NEW
!     CARD IN FORMAT OF A1,A4,A3,64A1,A1,A4,A3.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,ITYPE  ,MA     ,MB     ,NUM
      
! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER KA(70)
      ITYPE=0
!
!     RETURN IF MNEUMONIC IS ALREADY LEFT-ADJUSTED.
!              (SEE IF COL. 1 IS BLANK)
!
      IF(KA(1).NE.MB(2)) RETURN
!
      BACKSPACE IOU5
      READ(IOU5,10) (KA(I),I=1,70)
   10 FORMAT(70A1)
!
!     RETURN IF CARD IS ENDDATA CARD.
!
      ITYPE=NDATA(KA)
      IF(ITYPE.EQ.1) RETURN
      BACKSPACE IOU5
!
!     COUNT NUMBER OF LEADING BLANKS.
!
      DO I=2,8
         IF(KA(I).NE.MB(2)) GO TO 40
      end do
      READ(IOU5,30) (KA(I),I=1,70)
   30 FORMAT(A1,A4,A3,65A1,A4,A3)
      RETURN
!
!     VARIABLE FORMATS FOR LEFT-ADJUSTING.
!
   40 KA(2)=MB(2)
      KA(3)=MB(2)
      I=I-1
      GO TO (50,60,70,80,90,100,110), I
!     I=NUMBER OF LEADING BLANKS
!
!     1 BLANK
!
   50 READ(IOU5,55) (KA(I),I=1,70)
   55 FORMAT(1X,A1,A4,A2,65A1,A4,A3)
      RETURN
!
!     2 BLANKS
!
   60 READ(IOU5,65) (KA(I),I=1,70)
   65 FORMAT(2X,A1,A4,66A1,A4,A3)
      RETURN
!
!     3 BLANKS
!
   70 READ(IOU5,75) KA(1),KA(2),(KA(I),I=4,70)
   75 FORMAT(3X,A1,A4,65A1,A4,A3)
      RETURN
!
!     4 BLANKS
!
   80 READ(IOU5,85) KA(1),KA(2),(KA(I),I=4,70)
   85 FORMAT(4X,A1,A3,65A1,A4,A3)
      RETURN
!
!     5 BLANKS
!
   90 READ(IOU5,95) KA(1),KA(2),(KA(I),I=4,70)
   95 FORMAT(5X,A1,A2,65A1,A4,A3)
      RETURN
!
!     6 BLANKS
!
  100 READ(IOU5,105) KA(1),KA(2),(KA(I),I=4,70)
  105 FORMAT(6X,67A1,A4,A3)
      RETURN
!
!     7 BLANKS
!
  110 READ(IOU5,115) KA(1),(KA(I),I=4,70)
  115 FORMAT(7X,66A1,A4,A3)
      RETURN
!
      END SUBROUTINE LEFT

! ##################################################################################################################################
      FUNCTION LENCAS(KA,N)
!
!     DETERMINE THE LENGTH OF AN EXECUTIVE/CASE CONTROL LINE.
!     N = DIMENSION OF ARRAY KA
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,LENCAS ,N

! E////////////////////////////////////////////////////////////////////E
      INTEGER KA(N),BLANK
! B////////////////////////////////////////////////////////////////////B
      DATA BLANK/1H /
! E////////////////////////////////////////////////////////////////////E
      DO 10 I=N,1,-1
      IF(KA(I).NE.BLANK) GO TO 20
10    CONTINUE
      LENCAS=1
      RETURN
20    LENCAS=I
      RETURN
      END FUNCTION LENCAS

! ##################################################################################################################################
      FUNCTION MAXDGR(NC,IC,IDEG)
!
!     THIS FUNCTION HAS AS ITS VALUE THE MAXIMUM DEGREE OF ANY NODE OF
!     COMPONENT NC IF NC.GT.0
!     IF NC.LE.0, ALL COMPONENTS ARE CONSIDERED.
!
      INTEGER IC(*),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,M      ,MAXDGR ,NC     ,NN

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      M=0
      DO 100 I=1,NN
      IF(NC)40,50,40
40    IF(IC(I)-NC) 100,50,100
50    IF(IDEG(I)-M) 100,100,60
60    M=IDEG(I)
  100 CONTINUE
      MAXDGR=M
      RETURN
      END FUNCTION MAXDGR

! ##################################################################################################################################
      FUNCTION MINDEG(NC,IC,IDEG)
!
!     THIS FUNCTION HAS AS ITS VALUE THE MINIMUM DEGREE OF ANY NODE OF
!     COMPONENT NC IF NC.GT.0
!     IF NC.LE.0, ALL COMPONENTS ARE CONSIDERED.
!
      INTEGER IC(*),IDEG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,M      ,MINDEG ,NC     ,NN

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      M=600000
      DO 100 I=1,NN
      IF(NC)40,50,40
40    IF(IC(I)-NC) 100,50,100
50    IF(M-IDEG(I)) 100,100,60
60    M=IDEG(I)
  100 CONTINUE
      MINDEG=M
      RETURN
      END FUNCTION MINDEG

! ##################################################################################################################################
      SUBROUTINE MORRIS(LIST,NL,IG,II1)
!
!     THIS ROUTINE DELETES ALL REFERENCE IN THE CONNECTION TABLE IG
!     TO THOSE POINTS IN A LIST OF LENGTH NL.
!
!     NEDGE = NUMBER OF UNIQUE EDGES.
!
!     REVISED 12/4/91 (TO AVOID CRAY COMPILER BUG)
!
      INTEGER II1
      DIMENSION IG(II1,*),LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IER    ,IG     ,IJ     ,                                &
               J      ,K      ,KMOD   ,L      ,LIST   ,                        &
               MAXDEG ,MAXGRD ,MM     ,MM1    ,                                &
               N      ,NBITIN ,NEDGE  ,NL     ,NN

      REAL     DUM    ,DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUM(3),NEDGE,DUMS(3)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!
      IF(NL.LE.0) RETURN
      MM1=MM-1
      DO 60 IJ=1,NL
      I=LIST(IJ)
      DO 50 J=1,MM
      L=IG(I,J)
!PACK L=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
      IF(L.EQ.0) GO TO 60
      NEDGE=NEDGE-1
      DO 10 K=1,MM
      IF(IG(L,K).EQ.I) GO TO 15
!PACK IF(IUNPK(IG,MAXGRD*(K-1)+L,NBITIN).EQ.I) GO TO 15
10    CONTINUE
      WRITE(IOU6,'(A)') ' Logic error in MORRIS'
      call finish(6,IER)
15    IF(K.GE.MM) GO TO 40
      DO 30 N=K,MM1
      IG(L,N)=IG(L,N+1)
!PACK IS=IUNPK(IG,MAXGRD*N+L,NBITIN)
!PACK CALL PACK(IG,MAXGRD*(N-1)+L,NBITIN,IS)
30    CONTINUE
40    CONTINUE
      IG(L,MM)=0
      IG(I,J)=0
!PACK CALL PACK(IG,MAXGRD*MM1+L,NBITIN,0)
!PACK CALL PACK(IG,MAXGRD*(J-1)+I,NBITIN,0)
   50 CONTINUE
   60 CONTINUE
      RETURN
      END SUBROUTINE MORRIS

! ##################################################################################################################################
      SUBROUTINE MPC(KA,ITYPE,KG,MAXI,INV,II3,NORIG,IER)
!
!     EXTRACT GRID POINTS FROM MPC EQUATION AND STORE IN KG.
!
!     MAXI=MAXIMUM NUMBER OF GRID POINTS ALLOWED PER ELEMENT.
!     MODIFIED 4/16/93 TO ADD MPCAX CARDS (SHORT FIELD ONLY).
!     ASSUMES FIELD 6 OF FIRST LOGICAL MPCAX CARD IS NOT BLANK.
!
      INTEGER MAXI
      INTEGER KA(70),KG(MAXI),INV(*),NORIG(*),ENDP(2)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IDUM   ,IER    ,II3    ,IPARAM ,ITYPE  ,                &
               J      ,K      ,L      ,                                        &
               MA     ,MB     ,NIP    ,NOUT   ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      NOUT=IOU12
!
!     L=1 OR 2 FOR SHORT OR LONG FIELDS, RESPECTIVELY
      L=ITYPE-1
      IF(ITYPE.EQ.112) L=1
!     I=COUNTER ON THE GRID POINT COMING UP (I.E., THE NEXT ONE).
      I=1
!
   20 J=4+8*L
      IF(ITYPE.EQ.112.AND.I.EQ.1) J=36
      IF(ITYPE.EQ.112.AND.I.GT.1) J=4
      IF(I.GT.MAXI) GO TO 100
      CALL READIT(KA(J),1,8*L,KG(I),NIP)
      I=I+NIP
      GO TO (30,40), L
   30 J=36
      IF(ITYPE.EQ.112.AND.I.EQ.2) GO TO 65
      GO TO 60
   40 ENDP(1)=KA(69)
      ENDP(2)=KA(70)
      READ(IOU5,50,END=120) KA
   50 FORMAT(A1,A4,A3,65A1,A4,A3)
!
!     LEFT-ADJUST FIELD 1.
!
      CALL LEFT(KA,IDUM)
      IF(KA(1).NE.MB(4)) GO TO 90
      IF(KA(2).NE.ENDP(1)) GO TO 90
      IF(KA(3).NE.ENDP(2)) GO TO 90
      WRITE(IOU8,50) KA
      J=4
   60 IF(I.GT.MAXI) GO TO 100
      CALL READIT(KA(J),1,8*L,KG(I),NIP)
      I=I+NIP
65    ENDP(1)=KA(69)
      ENDP(2)=KA(70)
      READ(IOU5,50,END=120) KA
!
!     LEFT-ADJUST FIELD 1.
!
      CALL LEFT(KA,IDUM)
      DO K=3,4
         IF(KA(1).EQ.MB(K)) GO TO 80
      end do
      GO TO 90
   80 IF(KA(2).NE.ENDP(1)) GO TO 90
      IF(KA(3).NE.ENDP(2)) GO TO 90
      WRITE(IOU8,50) KA
      L=1
      IF(K.EQ.5) L=2
      GO TO 20
!
!     END OF LOOP STARTING AT STATEMENT 20.
!
   90 BACKSPACE IOU5
      I=I-1
!
!     CONVERT ORIGINAL GRID NUMBERS TO INTERNAL LABELS.
!
      CALL SCAT(KG,I,INV,II3,NORIG,IER)
      IF(IER.GT.0) RETURN
!
!     DELETE DUPLICATE ENTRIES IN LIST.
!
      DO J=2,I
         IF(KG(J).EQ.KG(1)) KG(J)=0
      end do
      CALL FIXIT(KG,I)
!
!     WRITE OUT LIST OF NODES.
!
      WRITE(NOUT) I,(KG(J),J=1,I)
      RETURN
!
!     FATAL ERROR MESSAGE IF MPC EQUATION HAS MORE THAN MAXI TERMS.
!
  100 WRITE(IOU6,110) MAXI
  110 FORMAT(/,'Fatal Error.  MPC equation has more than',I8,' terms.')
      WRITE(IOU6,115)
  115 FORMAT('Use $DIM N card, where 4*N exceeds the maximum number',          &
       ' of terms in any one MPC equation.')
      call finish(6,IER)
      RETURN
!
!     END-OF-FILE ENCOUNTERED
!
120   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE MPC

! ##################################################################################################################################
      SUBROUTINE NASNUM(IG,II1,INV,II3,INT,ICC,ILD,NORIG,IP,KOR,               &
                        KORDIM,IER)
!
!     READ BULK DATA, SET UP CONNECTION TABLE, RESEQUENCE NODES,
!     AND GENERATE SEQGP CARDS.
!
      DIMENSION IG(*),INV(*),INT(*),ICC(*),ILD(*),NORIG(*),IP(*)
      INTEGER KORDIM ,KOR(KORDIM)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IADD   ,IBW2   ,IBYTE  ,ICC    ,ICRIT  ,IDIM   ,        &
               IER    ,IFIR   ,IG     ,IGDEG  ,IGNORE ,II1    ,II3    ,        &
               IIG    ,ILD    ,INT    ,INV    ,IP     ,IPARAM ,IPASS  ,        &
               IPF2   ,IPR    ,ISTA   ,ISTART ,IWALL  ,                        &
               J      ,JUMP   ,                                                &
               K      ,K1     ,K2     ,K3     ,K4     ,K5     ,KDIM   ,        &
               KMOD   ,KORE   ,L      ,                                        &
               MAXDEG ,MAXGRD ,MEM    ,MINDEG ,MM     ,                        &
               N      ,NBW    ,NBYTE  ,NCM    ,NEDGE  ,NEL    ,NEQ    ,        &
               NEQR   ,NLINK  ,NN     ,NORIG  ,NP     ,NW     ,NZERO  ,        &
               OBW    ,OP

      REAL     DUM    ,DUMM   ,DUMS   ,DUMY   ,TA     ,TB

! E////////////////////////////////////////////////////////////////////E
      COMMON /BITS/ DUM,KORE,DUMM(2),IPASS,NW,NBYTE,IBYTE,KDIM
      COMMON /B/ IPARAM(20)
      COMMON /S/ NN,MM,DUMS(3),NEDGE,IADD,MINDEG,DUMY
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /DOL/ ISTART(100),IGNORE(100)
      COMMON /DOLL/ IDIM,ISTA,IIG,IFIR,IGDEG
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!
      CALL TIMER(TA,IWALL,0,IOU6)
!
!     COPY IOU5 TO IOU11 AND REDEFINE IOU5.
!
      CALL NEWIN(IG,IER)
      IF(IER.GT.0) RETURN
!
!     ZERO OUT THE WORKING STORAGE.  THE ARRAY NAME USED HERE MUST
!        BE THE FIRST ONE.
!
      DO 10 I=1,KORE
10    ILD(I)=0
!
      NN=0
      MM=0
      NEDGE=0
      IPASS=0
      MINDEG=500000
      KMOD=2*MAXGRD-IFIX(2.3715*SQRT(FLOAT(MAXGRD)))
      REWIND IOU9
!
!     READ BULK DATA DECK AND SET UP CONNECTION TABLE IG.
!
      CALL REED(IG,II1,INV,II3,NORIG,KOR,KORDIM,INT,IER)
      IF(IER.GT.0) RETURN
      REWIND IOU11
!
      IF(NN.GT.0) GO TO 16
!
      WRITE(IOU6,15)
   15 FORMAT(/,'Fatal Error.  No grid points found on',                        &
             ' connection cards.')
      call finish(6,IER)
      RETURN
!
   16 IF(MM.GT.0) GO TO 18
      WRITE(IOU6,17)
   17 FORMAT(/,'Fatal Error.  No connections found.')
      call finish(6,IER)
      RETURN
!
   18 I=NEQ+IIG+IGDEG
      IF(I.LE.0) GO TO 19
!
!     MODIFY IG TO ACCOUNT FOR MPC EQUATIONS.
!
      CALL TIGER(IG,II1,ICC,NORIG,KOR,KORDIM,IER)
      IF(IER.GT.0) RETURN
!
!     SET UP LIST OF POINTS TO IGNORE IN ICC ARRAY.
!
      CALL IGNOR(IG,II1,INV,II3,ICC,N,IER)
      IF(IER.GT.0) RETURN
!
!     DELETE POINTS LISTED IN ICC FROM IG.
!
      CALL MORRIS(ICC,N,IG,II1)
!
!     SORT ORIGINAL GRID NUMBERS AND OUTPUT LIST IN INT.   DETERMINE
!     CORRESPONDENCE ILD BETWEEN UNSORTED AND SORTED INTERNAL LABELS.
!
   19 CALL BRIGIT(INV,II3,INT,ILD,IER)
      IF(IER.GT.0) RETURN
!
      CALL TIMER(TB,IWALL,0,IOU6)
      TB=TB-TA
      WRITE(IOU6,195) TB
195   FORMAT('CP time to set up connection table',F12.3,' seconds')
      WRITE(IOU6,196) NN
196   FORMAT('Number of grid points appearing on connection cards',I9/         &
             'Grid cards are not used.')
      WRITE(IOU6,198) MM
  198 FORMAT('Maximum nodal degree before any nodes are ignored',I8)
!
!     WRITE OUT CONNECTION TABLE IF REQUESTED.
!
      IF(IPARAM(2).EQ.4) CALL PUNCON(IG,II1,IP,ILD)
!
      IF(IPARAM(10).EQ.5) GO TO 20
!
!     PRINT TABLES IF MAXIMUM PRINTING REQUESTED.
!
!     CALL TABLE1(INV,II3,INT,IER)
!     The following four lines replace the call to Table1.
      WRITE(IOU6,23)
23    FORMAT(/4(5X,'Grid   Grid BANDIT')/4(7X,'ID  Count  Label'))
      write(iou6,24) (int(i),i,intern(int(i),inv,ii3,ier),i=1,nn)
24    FORMAT(4(I9,2I7))
!
      IF(IER.GT.0) RETURN
!
      CALL TABLE2(IG,II1,NORIG,IP)
20    continue
!
!     CONVERT USER-SELECTED STARTING NODES, IF ANY, TO INTERNAL LABELS.
!
      IF(ISTA.GT.0) CALL FLIP(ISTART,ISTA,INV,II3,IER)
      IF(IER.GT.0) RETURN
!
!     SAVE ORIGINAL ORDERING (ILD) IF NECESSARY.
!
      CALL RESET(1,ILD,INT,JUMP)
!
      IPR=0
      IF(IPARAM(10).EQ.6) IPR=1
!     INITIALIZE JUMP.
      JUMP=1
!     JUMP=1, NO IMPROVEMENT OF CRITERION SELECTED
!         =0, IMPROVEMENT
!
!     CHOOSE CRITERION FOR RESEQUENCING
!
      ICRIT=2
      IF(IPARAM(6).EQ.11) ICRIT=3
      IF(IPARAM(6).EQ.12) ICRIT=1
      IF(IPARAM(6).EQ.13) ICRIT=4
!
      I=MAXGRD+2
      J=I+MAXGRD
      K=J+MAXGRD
      IF(IPARAM(3).EQ.8) GO TO 25
!
!     RESEQUENCE NODES WITH CUTHILL-MCKEE ALGORITHM.
!
      CALL CUTHIL(80,1,2,ICRIT,IPR,IG,II1,INV,INV(I),INV(J),INV(K),            &
        INT,ICC,ILD,IP,JUMP,KOR,KORDIM,IER)
!
!     Write to a file the component and original grid id for each grid
!     (to allow sorting to get list of grids in each component).
!     On exit from subroutine cuthil, inv contains array 'ic'.
!     Method CM must be included to get this list.
!
      if(iparam(10).eq.6) then
         open(iou15,file='bandit.f15',form='formatted',status='replace')
         write(iou15,'(a)') 'Component   Grid'
         write(iou15,'(i5,i11)') (inv(i),norig(i),i=1,nn)
         write(iou6,22)
22       format(/,'See bandit.f15 for component list.')
         close(iou15)
      end if
      IF(IER.GT.0) RETURN
!
      IF(IPARAM(3).EQ.7) GO TO 28
!
   25 CONTINUE
!
!     RESET ILD, IF NECESSARY, AND COPY TO INT.
!
      CALL RESET(2,ILD,INT,JUMP)
!
!     SAVE SEQGP CARDS ON IOU9 AFTER EXECUTING CM.  THEN, IN CASE
!     GPS ABORTS DUE TO EXCEEDING SCRATCH DIMENSION, THE CM RESULTS
!     CAN BE RECOVERED AND WRITTEN TO UNIT 7.  SEE SUBROUTINE FINISH.
!     IADD IS IGNORED HERE.
!
      REWIND IOU9
      WRITE(IOU9,26) (NORIG(L),ILD(L),L=1,NN)
   26 FORMAT('SEQGP   ',8I8)
      REWIND IOU9
!
!
!     RESEQUENCE NODES WITH GPS ALGORITHM.
!
      KDIM=KORDIM/4
      K1=1
      K2=K1+KDIM
      K3=K2+KDIM
      K4=K3+KDIM
      K5=K4+KDIM/2
      CALL GIBSTK(IG,II1,INT,ILD,INV(I),INV,INV(J),INV(K),ICC,IBW2,IPF2,       &
        JUMP,ICRIT,KOR(K1),KOR(K2),KOR(K3),KOR(K4),KOR(K5),KDIM,IER)
      IF(IER.GT.0) RETURN
!
   28 CONTINUE
!
!     GENERATE SEQGP CARDS.
!
      CALL SEQGP(NORIG,ILD,INT,JUMP)
!
!     IF REQUESTED, GENERATE SET OF SCALAR SPRINGS (CELAS3) WITH SAME
!     CONNECTIVITY AS ORIGINAL STRUCTURE ON bandit.f09
!
      IF(IPARAM(7).EQ.4) CALL SPRING(IP)
!
!     GENERATE ELEMENT ORDERING FOR FRONTAL SOLVERS.
!
      IF(IPARAM(11).EQ.4) THEN
!        COMPUTE MEMORY AVAILABLE FOR ELEMENT SORT (INV ARRAY).
         MEM=KORE-(MAXGRD+1)
         CALL FRONT(KOR,ILD,NN,INV,MEM,IER)
         IF(IER.GT.0) RETURN
      END IF
!
      RETURN
      END SUBROUTINE NASNUM

! ##################################################################################################################################
      FUNCTION NBULK(KA)
!
!     THIS FUNCTION RETURNS 1 AS ITS VALUE IF A CARD READ BY 80A1 IS
!     THE BEGIN BULK CARD. IF NOT, 0 IS RETURNED.
!
      INTEGER KA(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,K      ,MA     ,MB     ,NBULK  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      NBULK=0
!
!     LOOK FOR FIRST NON-BLANK.
!
      DO I=1,64
         IF(KA(I).NE.MB(2)) GO TO 20
      end do
      RETURN
!
   20 IF(KA(I  ).NE.MA( 2)) RETURN
      IF(KA(I+1).NE.MA( 5)) RETURN
      IF(KA(I+2).NE.MA( 7)) RETURN
      IF(KA(I+3).NE.MA( 9)) RETURN
      IF(KA(I+4).NE.MA(14)) RETURN
      K=I+5
!
!     LOOK FOR FIRST NON-BLANK AFTER -BEGIN-.
!
      DO I=K,69
         IF(KA(I).NE.MB(2)) GO TO 40
      end do
      RETURN
!
   40 IF(KA(I  ).NE.MA( 2)) RETURN
      IF(KA(I+1).NE.MA(21)) RETURN
      IF(KA(I+2).NE.MA(12)) RETURN
      IF(KA(I+3).NE.MA(11)) RETURN
      NBULK=1
      RETURN
      END FUNCTION NBULK

! ##################################################################################################################################
      FUNCTION NDATA(KA)
!
!     THIS FUNCTION RETURNS 1 AS ITS VALUE IF A CARD READ BY 70A1 IS
!     THE ENDDATA CARD.  OTHERWISE, 0 IS RETURNED.
!
      INTEGER KA(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,MA     ,MB     ,NDATA  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      NDATA=0
!
!     LOOK FOR FIRST NON-BLANK.
!
      DO I=1,64
         IF(KA(I).NE.MB(2)) GO TO 20
      end do
      RETURN
!
!     LOOK FOR -ENDDATA-.
!
   20 IF(KA(I  ).NE.MA( 5)) RETURN
      IF(KA(I+1).NE.MA(14)) RETURN
      IF(KA(I+2).NE.MA( 4)) RETURN
      IF(KA(I+3).NE.MA( 4)) RETURN
      IF(KA(I+4).NE.MA( 1)) RETURN
      IF(KA(I+5).NE.MA(20)) RETURN
      IF(KA(I+6).NE.MA( 1)) RETURN
      NDATA=1
      RETURN
      END FUNCTION NDATA

! ##################################################################################################################################
      SUBROUTINE NEWIN(KA,IER)
!
!     COPY BULK DATA DECK FROM UNIT 5 TO UNIT 11.
!     REDEFINE IOU5 TO BE IOU11.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  IER    ,L

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER KA(80)
!     KA IS SCRATCH SPACE.
      L=IOU11
      REWIND L
    5 READ(IOU5,10,END=30) KA
   10 FORMAT(80A1)
      WRITE(L,10) KA
      IF(NDATA(KA).EQ.0) GO TO 5
      REWIND L
!     REDEFINE IOU5
      IOU5=L
      RETURN
!
30    CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE NEWIN

! ##################################################################################################################################
      SUBROUTINE NOSEQ(KA)
!
!     WRITE BULK DATA DECK IF RESEQUENCING NOT REQUESTED.
!
      INTEGER KA(80)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  IPARAM

! E////////////////////////////////////////////////////////////////////E
!     KA IS SCRATCH SPACE.
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      WRITE(IOU6,10)
   10 FORMAT(/,'Grid point resequencing not requested.' )
   20 READ(IOU5,30,END=40) KA
   30 FORMAT(80A1)
      WRITE(IOU8,30) KA
      IF(NDATA(KA).EQ.0) GO TO 20
40    RETURN
      END SUBROUTINE NOSEQ

! ##################################################################################################################################
      FUNCTION NTHRU(KA,N)
!
!     THIS FUNCTION RETURNS 1 AS ITS VALUE IF A FIELD READ N A1
!     CONTAINS THE CHARACTER STRING "THRU".  IF NOT, 0 IS RETURNED.
!
!     N=NUMBER OF CARD COLUMNS TO SEARCH STARTING AT KA(1)
!
      INTEGER KA(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,MA     ,MB     ,N      ,NA     ,NTHRU  ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      NTHRU=0
!
!     LOOK FOR FIRST NON-BLANK.
!
      NA=N-3
      IF(NA.LT.1) RETURN
      DO I=1,NA
         IF(KA(I).NE.MB(2)) GO TO 20
      end do
      RETURN
!
!     LOOK FOR  -THRU-.
!
   20 IF(KA(I  ).NE.MA(20)) RETURN
      IF(KA(I+1).NE.MA( 8)) RETURN
      IF(KA(I+2).NE.MA(18)) RETURN
      IF(KA(I+3).NE.MA(21)) RETURN
      NTHRU=1
      RETURN
      END FUNCTION NTHRU

! ##################################################################################################################################
      SUBROUTINE NUMBR(SND,NUM,NDSTK,LVLS2,NDEG,RENUM,LVLST,LSTPT,NR,          &
                       NFLG,IBW2,IPF2,IPFA,ISDIR,STKA,STKB,STKC,STKD,          &
                       IDIM,IER)
!
!  NUMBR PRODUCES THE NUMBERING OF THE GRAPH FOR MIN BANDWIDTH
!
!  SND-         ON INPUT THE NODE TO BEGIN NUMBERING ON
!  NUM-         ON INPUT AND OUTPUT, THE NEXT AVAILABLE NUMBER
!  LVLS2-       THE LEVEL STRUCTURE TO BE USED IN NUMBERING
!  RENUM-       THE ARRAY USED TO STORE THE NEW NUMBERING
!  LVLST-       ON OUTPUT CONTAINS LEVEL STRUCTURE
!  LSTPT(I)-    ON OUTPUT, INDEX INTO LVLST TO FIRST NODE IN ITH LVL
!               LSTPT(I+1) - LSTPT(I) = NUMBER OF NODES IN ITH LVL
!  NFLG-        =+1 IF SND IS FORWARD END OF PSEUDO-DIAM
!               =-1 IF SND IS REVERSE END OF PSEUDO-DIAM
!  IBW2-        BANDWIDTH OF NEW NUMBERING COMPUTED BY NUMBER
!  IPF2-        PROFILE OF NEW NUMBERING COMPUTED BY NUMBER
!      IBW2 AND IPF2 HERE DO NOT INCLUDE DIAGONAL TERMS.
!  IPFA-        WORKING STORAGE USED TO COMPUTE PROFILE AND BANDWIDTH
!  ISDIR-       INDICATES STEP DIRECTION USED IN NUMBERING(+1 OR -1)
!
      INTEGER IDIM, SND,XA,XB,XC,XD,CX,END,RENUM,TEST
      COMMON /GRA/ N,IDPTH,IDEG
      INTEGER STKA(IDIM),STKB(IDIM),STKC(IDIM),STKD(IDIM)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER NR, IPFA(*)
      DIMENSION NDSTK(NR,*),LVLS2(*),NDEG(*),RENUM(*),LVLST(*),LSTPT(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IBW2   ,IDEG   ,IDPTH  ,IER    ,INX    ,                &
               IPF2   ,IPRO   ,ISDIR  ,                                        &
               J      ,                                                        &
               KMOD   ,                                                        &
               LND    ,LST    ,LSTPT  ,LVLN   ,LVLS2  ,LVLST  ,                &
               MAXDEG ,MAXGRD ,MAXI   ,                                        &
               N      ,NBITIN ,NBW    ,NDEG   ,NDSTK  ,NFLG   ,                &
               NSTPT  ,NUM

      REAL     DUMBB 

! E////////////////////////////////////////////////////////////////////E
!  SET UP LVLST AND LSTPT FROM LVLS2
      DO 3 I=1,N
        IPFA(I)=0
    3 CONTINUE
      NSTPT=1
      DO 5 I=1,IDPTH
        LSTPT(I)=NSTPT
        DO 5 J=1,N
          IF(LVLS2(J).NE.I) GO TO 5
          LVLST(NSTPT)=J
          NSTPT=NSTPT+1
    5 CONTINUE
      LSTPT(IDPTH+1)=NSTPT
!  THIS ROUTINE USES FOUR STACKS, A,B,C,AND D, WITH POINTERS
!  XA,XB,XC, AND XD.  CX IS A SPECIAL POINTER INTO STKC WHICH
!  INDICATES THE PARTICULAR NODE BEING PROCESSED.
!  LVLN KEEPS TRACK OF THE LEVEL WE ARE WORKING AT.
!  INITIALLY STKC CONTAINS ONLY THE INITIAL NODE, SND.
      LVLN=0
      IF(NFLG.LT.0) LVLN=IDPTH+1
      XC=1
      STKC(XC)=SND
   10 CX=1
      XD=0
      LVLN=LVLN+NFLG
      LST=LSTPT(LVLN)
      LND=LSTPT(LVLN+1)-1
!  BEGIN PROCESSING NODE STKC(CX)
   20 IPRO=STKC(CX)
      RENUM(IPRO)=NUM
      NUM=NUM+ISDIR
      END=NDEG(IPRO)
      XA=0
      XB=0
!  CHECK ALL ADJACENT NODES
      DO 50 I=1,END
        TEST=NDSTK(IPRO,I)
!PACK   TEST=IUNPK(NDSTK,MAXGRD*(I-1)+IPRO,NBITIN)
26      INX=RENUM(TEST)
!  ONLY NODES NOT NUMBERED OR ALREADY ON A STACK ARE ADDED
        IF(INX.EQ.0) GO TO 30
        IF(INX.LT.0) GO TO 50
!  DO PRELIMINARY BANDWIDTH AND PROFILE CALCULATIONS
        NBW=(RENUM(IPRO)-INX)*ISDIR
        IF(ISDIR.GT.0) INX=RENUM(IPRO)
        IF(IPFA(INX).LT.NBW) IPFA(INX)=NBW
        GO TO 50
   30   RENUM(TEST)=-1
!  PUT NODES ON SAME LEVEL ON STKA, ALL OTHERS ON STKB
        IF(LVLS2(TEST).EQ.LVLS2(IPRO)) GO TO 40
        XB=XB+1
      IF(XB.GT.IDIM) GO TO 150
        STKB(XB)=TEST
        GO TO 50
   40   XA=XA+1
      IF(XA.GT.IDIM) GO TO 150
        STKA(XA)=TEST
   50 CONTINUE
!  SORT STKA AND STKB INTO INCREASING DEGREE AND ADD STKA TO STKC
!  AND STKB TO STKD
      IF(XA.EQ.0) GO TO 55
      IF(XA.EQ.1) GO TO 52
      CALL SORTDG(STKC,STKA,XC,XA,NDEG)
      GO TO 55
   52 XC=XC+1
      IF(XC.GT.IDIM) GO TO 150
      STKC(XC)=STKA(XA)
   55 IF(XB.EQ.0) GO TO 65
      IF(XB.EQ.1) GO TO 62
      CALL SORTDG(STKD,STKB,XD,XB,NDEG)
      GO TO 65
   62 XD=XD+1
      IF(XD.GT.IDIM) GO TO 150
      STKD(XD)=STKB(XB)
!  BE SURE TO PROCESS ALL NODES IN STKC
   65 CX=CX+1
      IF(XC.GE.CX) GO TO 20
!  WHEN STKC IS EXHAUSTED LOOK FOR MIN DEGREE NODE IN SAME LEVEL
!  WHICH HAS NOT BEEN PROCESSED
      MAXI=IDEG+1
      SND=N+1
      DO 70 I=LST,LND
        TEST=LVLST(I)
        IF(RENUM(TEST).NE.0) GO TO 70
        IF(NDEG(TEST).GE.MAXI) GO TO 70
        RENUM(SND)=0
        RENUM(TEST)=-1
        MAXI=NDEG(TEST)
        SND=TEST
   70 CONTINUE
      IF(SND.EQ.N+1) GO TO 75
      XC=XC+1
      IF(XC.GT.IDIM) GO TO 150
      STKC(XC)=SND
      GO TO 20
!  IF STKD IS EMPTY WE ARE DONE, OTHERWISE COPY STKD ONTO STKC
!  AND BEGIN PROCESSING NEW STKC
   75 IF(XD.EQ.0) GO TO 100
      DO 80 I=1,XD
        STKC(I)=STKD(I)
   80 CONTINUE
      XC=XD
      GO TO 10
!  DO FINAL BANDWIDTH AND PROFILE CALCULATIONS
  100 DO 120 I=1,N
        IF(IPFA(I).GT.IBW2) IBW2=IPFA(I)
        IPF2=IPF2+IPFA(I)
  120 CONTINUE
      RETURN
!     DIMENSION EXCEEDED  . . .  STOP JOB.
  150 CALL FINISH(3,IER)
      RETURN
      END SUBROUTINE NUMBR

! ##################################################################################################################################
      SUBROUTINE PIKLVL(LVLS1,LVLS2,CCSTOR,IDFLT,ISDIR,XC,NHIGH,NLOW,          &
                        NACUM,SIZE,STPT)
!
!  PIKLVL CHOOSES THE LEVEL STRUCTURE  USED IN NUMBERING GRAPH
!
!  LVLS1-       ON INPUT CONTAINS FORWARD LEVELING INFO
!  LVLS2-       ON INPUT CONTAINS REVERSE LEVELING INFO
!               ON OUTPUT THE FINAL LEVEL STRUCTURE CHOSEN
!  CCSTOR-      ON INPUT CONTAINS CONNECTED COMPONENT INFO
!  IDFLT-       ON INPUT =1 IF WDTH LVLS1@WDTH LVLS2, =2 OTHERWISE
!  NHIGH        KEEPS TRACK OF LEVEL WIDTHS FOR HIGH NUMBERING
!  NLOW-        KEEPS TRACK OF LEVEL WIDTHS FOR LOW NUMBERING
!  NACUM-       KEEPS TRACK OF LEVEL WIDTHS FOR CHOSEN LEVEL STRUCTURE
!  XC-          NUMBER OF CONNECTED COMPONENTS
!  SIZE(I)-     SIZE OF ITH CONNECTED COMPONENT
!  STPT(I)-     INDEX INTO CCSTORE OF 1ST NODE IN ITH CON COMPT
!  ISDIR-       FLAG WHICH INDICATES WHICH WAY THE LARGEST CONNECTED
!               COMPONENT FELL.  =+1 IF LOW AND -1 IF HIGH
!
      INTEGER XC,END
      COMMON /GRA/ N,IDPTH,DUMG
!     DIMENSION OF NHIGH IS MAXIMUM ALLOWABLE NUMBER OF LEVELS.
!     DIMENSION OF SIZE IS MAXIMUM ALLOWABLE NUMBER OF COMPONENTS.
      INTEGER NHIGH(*),NLOW(*),NACUM(*),SIZE(*),STPT(*)
      INTEGER LVLS1(*),LVLS2(*),CCSTOR(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDFLT  ,IDPTH  ,INODE  ,ISDIR  ,IT     ,                &
               J      ,K      ,LVLNH  ,LVLNL  ,                                &
               MAX1   ,MAX2   ,N

      REAL     DUMG

! E////////////////////////////////////////////////////////////////////E
!  FOR EACH CONNECTED COMPONENT DO
      DO 270 I=1,XC
        J=STPT(I)
        END=SIZE(I)+J-1
!  SET NHIGH AND NLOW EQUAL TO NACUM
        DO 205 K=1,IDPTH
          NHIGH(K)=NACUM(K)
          NLOW(K)=NACUM(K)
  205   CONTINUE
!  UPDATE NHIGH AND NLOW FOR EACH NODE IN CONNECTED COMPONENT
        DO 210 K=J,END
          INODE=CCSTOR(K)
          LVLNH=LVLS1(INODE)
          NHIGH(LVLNH)=NHIGH(LVLNH)+1
          LVLNL=LVLS2(INODE)
          NLOW(LVLNL)=NLOW(LVLNL)+1
  210   CONTINUE
        MAX1=0
        MAX2=0
!  SET MAX1=LARGEST NEW NUMBER IN NHIGH
!  SET MAX2=LARGEST NEW NUMBER IN NLOW
        DO 240 K=1,IDPTH
          IF(2*NACUM(K).EQ.NLOW(K)+NHIGH(K)) GO TO 240
          IF(NHIGH(K).GT.MAX1) MAX1=NHIGH(K)
          IF(NLOW(K).GT.MAX2) MAX2=NLOW(K)
  240   CONTINUE
!  SET IT= NUMBER OF LEVEL STRUCTURE TO BE USED
        IT=1
        IF(MAX1.GT.MAX2) IT=2
        IF(MAX1.EQ.MAX2) IT=IDFLT
        IF(IT.EQ.2) GO TO 265
        IF(I.EQ.1) ISDIR=-1
!  COPY LVLS1 INTO LVLS2 FOR EACH NODE IN CONNECTED COMPONENT
        DO 260 K=J,END
          INODE=CCSTOR(K)
          LVLS2(INODE)=LVLS1(INODE)
  260   CONTINUE
!  UPDATE NACUM TO BE THE SAME AS NHIGH
        DO 262 K=1,IDPTH
          NACUM(K)=NHIGH(K)
  262   CONTINUE
        GO TO 270
!  UPDATE NACUM TO BE THE SAME AS NLOW
  265   DO 267 K=1,IDPTH
          NACUM(K)=NLOW(K)
  267   CONTINUE
  270 CONTINUE
      RETURN
      END SUBROUTINE PIKLVL

! ##################################################################################################################################
      SUBROUTINE PUNCON(IG,II1,IP,ILD)
!
!     WRITE CONNECTION TABLE IG ON PERIPHERAL FILE (IOU16).
!
!     ILD(I) = SORTED INTERNAL LABEL FOR NODE WITH UNSORTED INTERNAL
!              LABEL I
!     NN     = NUMBER OF NODES.
!     M      = MAX NODAL DEGREE.
!     IP     = TEMPORARY STORAGE.
!
!     HEADER CARD GIVES NUMBER OF NODES (NN) AND MAX NODAL DEGREE (M)
!     IN 2I5, THEN ONE CARD PER NODE, 24I5.
!     FIELD 1 IS NODE, OTHER FIELDS ARE CONNECTIONS.
!
      INTEGER II1
      DIMENSION IG(II1,*),IP(*),ILD(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IG     ,ILD    ,IP     ,                                &
               J      ,K      ,KMOD   ,                                        &
               M      ,MAXDEG ,MAXGRD ,                                        &
               NBITIN ,NN

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,M,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!
      WRITE(IOU16,30) NN,M
      DO 20 I=1,NN
      DO J=1,M
         IP(J)=0
      end do
      DO J=1,M
         K=IG(I,J)
!PACK    K=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
         IF(K.EQ.0) GO TO 15
         IP(J)=ILD(K)
      end do
   15 CONTINUE
      WRITE(IOU16,30) ILD(I),(IP(J),J=1,M)
   30 FORMAT(24I5)
!
!     24I5 IS CHOSEN SINCE 24*5=120 IS AN INTEGER MULTIPLE OF 10 (CDC),
!     4 (IBM), AND 6 (UNIVAC), MAKING IT CONVENIENT FOR I/O BETWEEN
!     DISSIMILAR COMPUTERS.
!     120 CHARACTER LINES CAN ALSO BE LISTED ON A PRINTER IF DESIRED.
!     IF THIS FORMAT IS CHANGED, SIMILAR CHANGES ARE REQUIRED IN
!     SUBROUTINES SEQGP AND SPRING.
!
   20 CONTINUE
      WRITE(IOU6,40) IOU16
   40 FORMAT(/,'Connection table generated on bandit.f',I2)
      RETURN
      END SUBROUTINE PUNCON

! ##################################################################################################################################
      SUBROUTINE READIT(KA,N,MAXI,IP,NIP)
!
!     INTERPRET NUMERIC DATA READ IN 80A1 FORMAT.  ONLY NON-NEGATIVE
!     INTEGERS ARE FOUND.
!
!     KA(I) = ITH CHARACTER (A1 FORMAT) (INPUT)
!     N     = NUMBER OF INTEGERS SOUGHT (INPUT)
!     MAXI  = MAXIMUM VALUE OF I FOR SEARCH (INPUT)
!     IP(J) = JTH INTEGER FOUND (OUTPUT)
!     NIP   = NUMBER OF INTEGERS FOUND (OUTPUT)
!
!     TO CONVERT READIT TO A STAND-ALONE UTILITY NOT DEPENDENT ON
!     COMMON BLOCK /ALPHA/, WHICH SUPPLIES THE ARRAY NUM, REPLACE
!     THE COMMON STATEMENT BELOW WITH THE DATA STATEMENT FOR NUM
!     WHICH IS COMMENTED OUT.
!
!     G.C. EVERSTINE, DTRC 128, SEPT. 1973 (REVISED JULY 1978)
!
      INTEGER KA(*),IP(*),FLAG,BLANK
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,KOL    ,                                        &
               MA     ,MAXI   ,MB     ,N      ,NIP    ,NUM

! E////////////////////////////////////////////////////////////////////E
!     DATA NUM/1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
! B////////////////////////////////////////////////////////////////////B
      DATA BLANK/1H /
! E////////////////////////////////////////////////////////////////////E
      NIP=0
      KOL=0
      DO 40 I=1,N
      IP(I)=0
      FLAG=0
   10 KOL=KOL+1
      IF(KOL.LE.MAXI) GO TO 15
      IF(FLAG.EQ.0) RETURN
      NIP=NIP+1
      RETURN
   15 CONTINUE
      IF(KA(KOL).EQ.BLANK) GO TO 25
      DO 20 J=1,10
      IF(KA(KOL).EQ.NUM(J)) GO TO 30
   20 CONTINUE
   25 IF(FLAG) 10,10,40
   30 FLAG=1
      IP(I)=10*IP(I)+J-1
      GO TO 10
   40 NIP=NIP+1
      RETURN
      END SUBROUTINE READIT

! ##################################################################################################################################
      SUBROUTINE REED(IG,II1,INV,II3,NORIG,KG,MAXI,IDUM,IER)
!
!     READ BULK DATA DECK AND SET UP CONNECTION TABLE IG.
!
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with EID we will use EID_array
! instead. Needed so Lahey doesn't complain about shape of EID being
! different than IP (array) in subr READIT
      INTEGER EID_array(1)
! E////////////////////////////////////////////////////////////////////E
      INTEGER II1, II3, MAXI
      DIMENSION IG(II1,*),INV(2,II3),NORIG(*)
      INTEGER vype,TYPE,WYPE,KA(70),KG(MAXI),IDUM(*),ENDP(2),EID
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,I1     ,I2     ,IER    ,IFLD   ,IG     ,II     ,        &
               INV    ,IPARAM ,ISAVE  ,ITYPE  ,                                &
               J      ,K      ,L      ,LEN    ,LOOP   ,                        &
               MA     ,MB     ,MDIM   ,ME     ,                                &
               NCARD  ,NCON   ,NEL    ,NELEM  ,NEQ    ,NEQR   ,NLINK  ,        &
               NORIG  ,npt    ,NTYPE  ,NUM

      REAL     DUMD   

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /ELEM/ NTYPE,VYPE(160),TYPE(160),WYPE(160),ME(160),               &
                    NELEM(160),MDIM
      COMMON /B/ IPARAM(20)
      COMMON /D/ DUMD(6),NEL,NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!     MAXI=MAXIMUM NUMBER OF GRID POINTS ALLOWED PER ELEMENT
!          (=kdim4 from main) (input)
      LOGICAL LESSOK,rigide
      IF(IPARAM(11).EQ.4) THEN
         OPEN(IOU13,FORM='UNFORMATTED',STATUS='SCRATCH')
         REWIND IOU13
      END IF
      NEL=0
      NEQ=0
      NEQR=0
      NLINK=0
      NCARD=0
!     NEL   = NUMBER OF ELEMENTS NOT COUNTING RIGID ELEMENTS
!     NEQ   = OVERALL NUMBER OF MPC EQUATIONS, INCLUDING NEQR
!     NEQR  = NUMBER OF MPC EQUATIONS ARISING FROM RIGID LINKS.
!     NLINK = NUMBER OF RIGID LINKS (ELEMENTS)
!     NCARD = NUMBER OF BULK DATA CARDS
!
!     READ BULK DATA CARD.
!
   20 READ(IOU5,30,END=170) KA
   30 FORMAT(A1,A4,A3,65A1,A4,A3)
      NCARD=NCARD+1
!     IF(MOD(NCARD,100).EQ.0) WRITE(IOU17,35) NCARD
35    FORMAT(I9)
!
!     LEFT-ADJUST MNEUMONIC.
!
      CALL LEFT(KA,ITYPE)
!
!     RETURN IF ENDDATA CARD  (A LEFT-ADJUSTED ENDDATA WILL NOT BE
!        FOUND UNTIL LATER).
!
   55 IF(ITYPE.EQ.1) THEN
         WRITE(IOU6,60) NCARD
60       FORMAT(/,'Number of Bulk Data records read',I14)
         RETURN
      END IF
!
!     DETERMINE CARD TYPE.
!
      CALL ELTYPE(KA,ITYPE,NCON,IFLD,LOOP,MAXI,LEN,LESSOK,IER)
      IF(IER.GT.0) RETURN
!
!     Set rigid element flag.
!
      rigide=.false.
      if(itype.ge.91.and.itype.le.93) rigide=.true.
      if(itype.ge.128.and.itype.le.137) rigide=.true.
!
!     IF NOT ENDDATA CARD, WRITE CARD AND CONTINUE.
!
      IF(ITYPE.EQ.1) THEN
         WRITE(IOU6,60) NCARD
         RETURN
      END IF
!
      WRITE(IOU8,30) KA
!
!     IF NO MATCH FOUND, CHECK FOR AXIC, THEN GO BACK AND READ
!             ANOTHER CARD.
!
      IF(ITYPE.EQ.0) THEN
         CALL TAXI(KA,IOU6)
         GO TO 20
      END IF
!
!     CHECK IF MPC CARDS ARE TO BE PROCESSED (MPC, MPC*, MPCAX)
!
      if(itype.le.3.or.itype.eq.112) then
         IF(IPARAM(4).EQ.3) GO TO 20
         NEQ=NEQ+1
!        NEQ=NUMBER OF MPC EQUATIONS
         CALL MPC(KA,ITYPE,KG,MAXI,INV,II3,NORIG,IER)
         IF(IER.GT.0) RETURN
         GO TO 20
      end if
!
!     CHECK FOR AND INTERPRET selected RIGID ELEMENTs
!
      IF(ITYPE.EQ.91.OR.ITYPE.EQ.92.or.                                        &
            (itype.ge.130.and.itype.le.133)) then
         CALL RIGID(KA,ITYPE,KG,MAXI,INV,II3,NORIG,IDUM,IER,npt)
84       format('frigid',8i6)
         IF(IER.GT.0) RETURN
!        Since RBE1 is read like a rigid element, but processed like a
!        regular element, keep going.
         if(itype.eq.130.or.itype.eq.131) go to 87
         NLINK=NLINK+1
         NELEM(ITYPE)=NELEM(ITYPE)+1
         GO TO 20
      END IF
87    continue
!
!     BLANK OUT FIELD 5 FOR SCALAR ELEMENTS WITH CONNECTIONS NOT LISTED
!         CONSECUTIVELY.
!
      IF(ITYPE.LE.9.or.(ITYPE.GE.73.AND.ITYPE.LE.75)) then
!        LEN=1 FOR SHORT FIELD CARD, 2 FOR LONG FIELD CARD
         I1=4+24*LEN
         I2=I1-1+8*LEN
         DO I=I1,I2
            KA(I)=MB(2)
         end do
      end if
!
!     PROCESS CONNECTION CARD.
!
!     LOOP=NUMBER OF ELEMENTS DEFINABLE ON ONE CARD (1 OR 2, USUALLY 1)
!
      DO 140 II=1,LOOP
!
!     GET ELEMENT ID (EID) (needed for frontal ordering)
!     I1=SUBSCRIPT IN KA ARRAY CORRESPONDING TO BEGINNING OF EID FIELD
      I1=4
      IF(II.EQ.2.AND.LEN.EQ.1) I1=36
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(I1),1,8*LEN,EID_array(1),J)
      EID = EID_array(1)
! E////////////////////////////////////////////////////////////////////E
!
!     DETERMINE GRID CONNECTIONS AND STORE IN KG.
!     THE ACTUAL NUMBER OF NON-ZERO NODES IS RETURNED IN NPT.
!     Grid connections have already been found for RBE1 by Subroutine RIGID.
!
      if(itype.eq.130.or.itype.eq.131) go to 103
        CALL SETUP(KA,NCON,IFLD,KG,NPT,LEN,LESSOK,IER)
        IF(IER.GT.0) RETURN
        IF(II.EQ.2.AND.NPT.EQ.0) GO TO 140
        IF(NPT.EQ.NCON.OR.LESSOK) GO TO 103
        WRITE(IOU6,102) KA
102     FORMAT('Warning.  Missing connection(s) on'/A1,A4,A3,65A1,A4,A3)
103   IF(NPT.LE.0) GO TO 130
!
      if(rigide) then
         nlink=nlink+1
      else
         nel=nel+1
      end if
      NELEM(ITYPE)=NELEM(ITYPE)+1
!
!     CONVERT KG FROM ORIGINAL TO INTERNAL NUMBERS.
!
      CALL SCAT(KG,npt,INV,II3,NORIG,IER)
      IF(IER.GT.0) RETURN
!
!     WRITE EID AND INTERNAL CONNECTIONS ON UNIT 13.
!     SKIP FOR "$FRONTAL NO" and for rigid elements.
!
      IF(IPARAM(11).EQ.3.or.eid.eq.0.or.rigide) GO TO 108
      WRITE(IOU13) EID,NPT,(KG(I),I=1,NPT)
108   CONTINUE
!
!     FILL CONNECTION TABLE IG.
!
      IF(npt.le.1) then
         npt=2
         KG(2)=0
      end if
      K=npt-1
      DO I=1,K
         L=I+1
         DO J=L,npt
            CALL SETIG(KG(I),KG(J),IG,II1,NORIG,IER)
            IF(IER.GT.0) RETURN
         end do
      end do
  130 IFLD=IFLD+4
      IF(LEN.EQ.1) GO TO 140
      IF(LOOP.EQ.1) GO TO 140
      IF(II.EQ.2) GO TO 140
      ENDP(1) = KA(69)
      ENDP(2) = KA(70)
      READ(IOU5,30,END=170) KA
      CALL LEFT(KA,ISAVE)
      IF(KA(1).NE.MB(4))   GO TO 150
      IF(KA(2).NE.ENDP(1)) GO TO 150
      IF(KA(3).NE.ENDP(2)) GO TO 150
      WRITE(IOU8,30) KA
      IFLD=4
140   CONTINUE
!
!     GO BACK AND READ ANOTHER CARD.
!
      GO TO 20
!
!     WARNING MESSAGE.
!
  150 WRITE(IOU6,160) KA
  160 FORMAT(/,'Warning.  The continuation to the card preceding'/             &
       A1,A4,A3,65A1,A4,A3/'does not immediately follow its parent',           &
       ' and hence will not be found by BANDIT.')
      ITYPE=ISAVE
      GO TO 55
!
!     END-OF-FILE ENCOUNTERED
!
170   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE REED

! ##################################################################################################################################
      SUBROUTINE RELABL(NS,NODES,IG,II1,IC,IDEG,IDIS,IW,NEW,ICC,               &
                        ILD,IAJ,IDIM,IER)
!
!     GENERATE A RELABELING SCHEME STARTING WITH NS NODES FOR WHICH
!     LABELS HAVE BEEN STORED IN ARRAY NODES.
!     SET UP ILD AND NEW.
!     ILD(OLD)=NEW
!     NEW(NEW)=OLD,   THE INVERSE OF ILD
!
      DIMENSION IG(II1,*),IC(*),IDEG(*),IDIS(*),IW(*),NEW(*),ICC(*)
      INTEGER ILD(*)
      COMMON /S/ NN,DUMS(8)
      INTEGER X
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER IDIM,NODES(*),IAJ(IDIM)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IA     ,IC     ,ICC    ,ICN    ,IDEG   ,                &
               IDIS   ,IER    ,IG     ,II     ,II1    ,IJ     ,IW     ,        &
               IZ     ,                                                        &
               J      ,JJ     ,JT     ,                                        &
               KI     ,KMOD   ,KO     ,                                        &
               L      ,LL     ,                                                &
               MAXDEG ,MAXGRD ,                                                &
               N      ,NBITIN ,NEW    ,NN     ,NNC    ,NS     ,NT

      REAL     DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      I=NODES(1)
      ICN=IC(I)
      NT=ICC(ICN)-1
      DO 50 I=1,NN
      IF(IC(I)-ICN) 50,40,50
40    IDIS(I)=0
50    CONTINUE
      DO 100 J=1,NS
      JJ=NODES(J)
      IDIS(JJ)=-1
      JT=J+NT
      NEW(JT)=JJ
100   ILD(JJ)=JT
      KI=NT
      KO=NS+NT
      LL=KO
      L=1
      J=KO
      NNC=ICC(ICN+1)-1
130   KI=KI+1
      IF(KI-LL)135,132,135
132   L=L+1
      LL=KO+1
135   II=NEW(KI)
      N=IDEG(II)
      IF(N)140,255,140
140   IJ=0
      DO 200 I=1,N
      IA=IG(II,I)
!PACK IA=IUNPK(IG,MAXGRD*(I-1)+II,NBITIN)
      IF(IDIS(IA)) 200,150,200
150   IJ=IJ+1
      IF(IJ.GT.IDIM) GO TO 270
      IDIS(IA)=L
      KO=KO+1
      IAJ(IJ)=IA
      IW(IJ)=IDEG(IA)
200   CONTINUE
      IF(IJ-1)250,210,220
210   J=KO
      IZ=IAJ(1)
      NEW(KO)=IZ
      ILD(IZ)=KO
      GO TO 250
220   X=0
221   DO 230 I=2,IJ
      IF(IW(I)-IW(I-1))224,230,230
224   CONTINUE
      X=IW(I)
      IW(I)=IW(I-1)
      IW(I-1)=X
225   X=IAJ(I)
      IAJ(I)=IAJ(I-1)
      IAJ(I-1)=X
230   CONTINUE
      IF(X)235,235,220
235   DO 240 I=1,IJ
      J=J+1
      IZ=IAJ(I)
      NEW(J)=IZ
      ILD(IZ)=J
240   CONTINUE
250   IF(KO-NNC)130,255,255
255   CONTINUE
!
!     REVERSE SEQUENCE FOR THIS COMPONENT (ICN).
!
!     ICC IS AN ARRAY USED FOR IDENTIFYING COMPONENTS IN THE NEW ARRAY.
!     ICC(I) CONTAINS THE INDEX FOR THE NEW ARRAY AT WHICH COMPONENT I
!        STARTS.
!     NEW(I) = OLD LABEL FOR NODE NOW LABELLED I.
!     I = NUMBER OF NODES IN PREVIOUS COMPONENTS.
!     J = NUMBER OF NODES IN LATER COMPONENTS.
      I=ICC(ICN)-1
      J=NN-ICC(ICN+1) +1
      IF(J.GT.NN) J=0
!
!     GET REVERSE CM SEQUENCE.
!
      CALL REVERS(NEW,ILD,NN,I,J)
!
      RETURN
!
!     DIMENSION EXCEEDED.  STOP JOB.
!
270   CALL FINISH(5,IER)
      RETURN
      END SUBROUTINE RELABL

! ##################################################################################################################################
      SUBROUTINE RESET(IGOTO,ILD,INT,JUMP)
!
!     THIS ROUTINE MAKES A COPY OF THE ORIGINAL GRID POINT SEQUENCE IF
!     BOTH METHODS ARE TO BE USED.  ALSO, IF CUTHILL ACHIEVES NO BW
!     REDUCTION, THIS SEQUENCE IS RETRIEVED BEFORE CALLING GIBSTK.
!
      INTEGER ILD(*),INT(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IGOTO  ,IPARAM ,JUMP   ,L      ,NN

      REAL     DUMS   

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DUMS(8)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      L=IOU12
      GO TO (10,20), IGOTO
!
!     SAVE ORIGINAL ORDERING (ILD) BEFORE CUTHILL IS CALLED.
!
   10 IF(IPARAM(3).NE. 9) RETURN
      REWIND L
      WRITE(L) (ILD(I),I=1,NN)
      RETURN
!
!     RESET ILD, IF NECESSARY, AND COPY TO INT.
!
   20 IF(IPARAM(3).NE. 9) GO TO 30
      REWIND L
!     IF CM MAKES NO IMPROVEMENT IN SEQUENCE, READ ORIGINAL SEQUENCE
!        BACK INTO ILD.
      IF(JUMP.NE.0) READ (L) (ILD(I),I=1,NN)
      REWIND L
!     ENDFILE L
!     REWIND L
!
   30 DO 40 I=1,NN
   40 INT(I)=ILD(I)
      RETURN
      END SUBROUTINE RESET

! ##################################################################################################################################
      SUBROUTINE REVERS(NEW,ILD,NN,N1,N2)
!
!     REVERSE THE NODAL SEQUENCE, OMITTING THE FIRST N1 AND THE LAST
!     N2 POINTS.
!
!     NEW(I) = OLD LABEL FOR NODE NOW LABELLED I (INPUT AND OUTPUT)
!     ILD(I) = NEW LABEL FOR NODE ORIGINALLY LABELED I (OUTPUT)
!     NN     = NUMBER OF NODES (INPUT)
!     N1     = NUMBER OF POINTS AT BEGINNING OF SEQUENCE TO OMIT FROM
!              REVERSAL (INPUT)
!     N2     = NUMBER OF POINTS AT END OF SEQUENCE TO OMIT FROM
!              REVERSAL (INPUT)
!
      INTEGER NEW(*),ILD(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,K      ,L      ,LL     ,M      ,                &
               N1     ,N2     ,NN

! E////////////////////////////////////////////////////////////////////E
!     J = NUMBER OF INTERCHANGES TO MAKE.
      J = (NN - N1 - N2)/2
      IF(J.LE.0) RETURN
      LL = NN - N2 + 1
!     MAKE INTERCHANGES IN NEW ARRAY.
      DO 10 I = 1 , J
      L=LL-I
      K=NEW(L)
      M=N1+I
      NEW(L) = NEW(M)
   10 NEW(M) = K
!     CORRECT ILD, THE INVERSE OF NEW.
      L=1+N1
      M=NN-N2
      DO 20 I=L,M
      K=NEW(I)
   20 ILD(K)=I
      RETURN
      END SUBROUTINE REVERS

! ##################################################################################################################################
      SUBROUTINE RIGID(KA,ITYPE,KG,MAXI,INV,II3,NORIG,LG,IER,npoint)
!
!     Extract grid points from rigid element card, and, in most cases,
!     generate equivalent MPCs.
!
!     CARD CONTENTS ARE STORED IN KA IN THE FORMAT A1,A4,A3,65A1,A4,A3
!
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with IG we will use IG_array
! instead. Needed so Lahey doesn't complain about shape of IG being
! different than IP (array) in subr READIT
      INTEGER IG_array(1)
! E////////////////////////////////////////////////////////////////////E

! B 02/21/04 //////////////////////////////////////////////////////////B
! When running RBE2-09-CBAR-10-A, I got msg that the subscript of array
! LG was out of range (tried to use LG(3) but LG only dimensioned 2.
! Gordon says to change dim to *
!!    INTEGER KA(70),KG(MAXI),INV(*),NORIG(*),ENDP(2),LG(2)
      INTEGER MAXI
      INTEGER KA(70),KG(MAXI),INV(*),NORIG(*),ENDP(2),LG(*)
! E////////////////////////////////////////////////////////////////////E

! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IBEG   ,IDUM   ,IEND   ,IER    ,IG     ,II3    ,        &
               IPARAM ,IT     ,ITYPE  ,                                        &
               K      ,KMOD   ,                                                &
               L      ,LOCDGP ,LOCIGP ,                                        &
               MA     ,MAXDEG ,MAXGRD ,MB     ,                                &
               N      ,NDGP   ,NEQ    ,NEQR   ,NF     ,NFSDGP ,NLINK  ,        &
               npoint ,NUM

      REAL     DUMD

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /B/ IPARAM(20)
      COMMON /D/ DUMD(7),NEQ,NEQR,NLINK
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!
!     DETERMINE RIGID ELEMENT TYPE.
!
!     IT = 1  CRIGD1 ON LEVEL 16 (of Nastran)
!          2  CRIGD2 ON LEVEL 16
!          3  CRIGD1 ON LEVEL 17 AND ABOVE
!          4  CRIGD2 ON LEVEL 17 AND ABOVE
!          5  CRIGD1 ON LEVEL 17 AND ABOVE (THRU OPTION)
!          6  CRBE1 or RBE1 (read here but treated like regular element)
!          7  CRBE2 or RBE2
!
!     This routine processes the elements listed above.
!     CRIGD1, CRIGD2, and RBE2 generate equivalent MPCs.
!     RBE1 is read in this routine because of the complicated syntax of
!     the element, but RBE1 is otherwise treated like a regular element.
!     The following rigid elements are processed and also handled like
!     regular elements:
!     CRIGDR(93), CRBAR(128), RBAR(129), CRROD(134), RROD(135),
!     CRTRPLT(136), RTRPLT(137).
!     No other rigid elements are currently recognizable.
!
!     CHECK FOR BLANK FIELD 3 TO DISTINGUISH BETWEEN LEVEL 16 AND
!     LEVEL 17 FORMATS FOR CRIGD1 AND CRIGD2
!     (WHICH ARE DIFFERENT, BELIEVE IT OR NOT)
!
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(12),1,8,IG_array(1),I)
      IG = IG_array(1)
! E////////////////////////////////////////////////////////////////////E
!
      IT=0
      IF(ITYPE.EQ.91.AND.I.EQ.0) IT=1
      IF(ITYPE.EQ.92.AND.I.EQ.0) IT=2
      IF(ITYPE.EQ.91.AND.I.EQ.1) IT=3
      IF(ITYPE.EQ.92.AND.I.EQ.1) IT=4
      IF(IT.EQ.3.AND.NTHRU(KA(28),8).EQ.1) IT=5
      IF(ITYPE.EQ.130.or.itype.eq.131) IT=6
      IF(ITYPE.EQ.132.or.itype.eq.133) IT=7
      IF(IT.EQ.0) RETURN
!
!     SET PARAMETERS FOR INTERPRETING physical CARD
!     (the first card of the logical card)
!
!     LOCIGP=LOCATION IN KA ARRAY OF INDEPENDENT GRID POINT
!     LOCDGP=LOCATION IN KA ARRAY OF FIRST DEPENDENT GRID POINT ON CARD
!     NDGP  =NUMBER OF DEP. G.P. ON CARD
!     NFSDGP=NUMBER OF FIELDS SEPARATING DEP. G.P. ON CARD
!
      GO TO (40,42,44,46,48,50,52), IT
   40 LOCIGP=20
      LOCDGP=28
      NDGP=5
      NFSDGP=1
      GO TO 75
   42 LOCIGP=20
      LOCDGP=36
      NDGP=2
      NFSDGP=2
      GO TO 75
   44 LOCIGP=12
      LOCDGP=20
      NDGP=6
      NFSDGP=1
      GO TO 75
   46 LOCIGP=12
      LOCDGP=20
      NDGP=3
      NFSDGP=2
      GO TO 75
   48 LOCIGP=12
      LOCDGP=20
      NDGP=2
      NFSDGP=2
      GO TO 75
   50 LOCIGP=12
      LOCDGP=28
      NDGP=2
      NFSDGP=2
      GO TO 75
   52 LOCIGP=12
      LOCDGP=28
      NDGP=5
      NFSDGP=1
      GO TO 75
   75 CONTINUE
!
!     GET INDEPENDENT GRID POINT IG.
!
! B////////////////////////////////////////////////////////////////////B
      CALL READIT(KA(LOCIGP),1,8,IG_array(1),N)
      IG = IG_array(1)
! E////////////////////////////////////////////////////////////////////E
!
!     IS IG THERE?
      IF(IG.eq.0) then
         WRITE(IOU6,20)
20       FORMAT(/,'Fatal Error.  Independent grid point on rigid',             &
                ' element not found.')
         WRITE(IOU6,295) KA
         call finish(6,IER)
         RETURN
      end if
      kg(1)=ig
      npoint=1
!
!     CONVERT IG TO INTERNAL LABEL.
!
! B////////////////////////////////////////////////////////////////////B
! Modify call to SCAT so that 1st arg is an array (like in subr SCAT)
      IG_array(1) = IG
      CALL SCAT(IG_array(1),1,INV,II3,NORIG,IER)
      IG = IG_array(1)
! E////////////////////////////////////////////////////////////////////E
      IF(IER.GT.0) RETURN
!
!     INTERPRET DEPENDENT GRID POINTS ON PARENT OR CONTINUATION CARD.
!
90    CONTINUE
      DO I=1,NDGP
         K=LOCDGP+(I-1)*8*NFSDGP
         CALL READIT(KA(K),1,8,LG(I),NF) ! This is OK. LG is an array
      end do
!     ELIMINATE ZEROS OR DUPLICATE ENTRIES FROM LIST OF DEP. G.P.
      NF=NDGP
      CALL FIXIT(LG,NF)
!     SET UP FOR NEXT CARD.
      GO TO (110,112,114,116,118,124,126), IT
  110 LOCDGP=4
      NDGP=8
      GO TO 150
  112 LOCDGP=4
      NDGP=4
      GO TO 150
  114 LOCDGP=4
      NDGP=8
      GO TO 150
  116 LOCDGP=4
      NDGP=4
      GO TO 150
!     LIST OUT POINTS INCLUDED ON THRU OPTION.
  118 CONTINUE
      IF(NF.LT.2) GO TO 120
      IBEG=LG(1)
      IEND=LG(2)
      N=IEND-IBEG
      IF(N.LT.1) GO TO 120
      IF(N.GT.(MAXGRD-1)) GO TO 120
      DO I=1,N
         LG(I+1)=IBEG+I
      end do
      NF=N+1
      GO TO 150
  120 WRITE(IOU6,121) KA
  121 FORMAT(/,'Fatal Error.  Illegal data on'/A1,A4,A3,65A1,A4,A3)
      call finish(6,IER)
      RETURN
  124 LOCDGP=12
      NDGP=3
      GO TO 150
  126 LOCDGP=4
      NDGP=8
      GO TO 150
150   CONTINUE
!     add points to KG list containing all points on rigid element
      if((npoint+nf).gt.maxi) then
         write(iou6,155)
155      format(/,'Fatal error.  Insufficient space for rigid element.')
         call finish(2,ier)
         return
      end if
      do i=1,nf
         kg(npoint+i)=LG(i)
      end do
      npoint=npoint+nf
      IF(NF.LE.0) GO TO 235
!
!     CONVERT new DEPENDENT POINTS TO INTERNAL LABELS
!     (except those rigid elements which are read in this routine but
!     treated like regular elements: RBE1)
!
      if(it.ne.6) then
         CALL SCAT(LG,NF,INV,II3,NORIG,IER)
         IF(IER.GT.0) RETURN
!
!        WRITE OUT EQUIVALENT MPCs ON IOU12
!
         DO L=1,NF
!           KG(1)=LG(L)
!           KG(2)=IG
!           I=2
!           WRITE(IOU12) I,(KG(J),J=1,I)
            write(iou12) 2,LG(L),IG
!           INCREMENT OVERALL MPC EQUATION COUNTER.
            NEQ=NEQ+1
!           INCREMENT COUNTER OF MPC EQUATIONS ARISING FROM RIGID ELEMENTS.
            NEQR=NEQR+1
         end do
      end if
      IF(IT.EQ.5) RETURN
!
!     READ ANOTHER CARD.
!
!     SAVE OLD FIELD 10 FOR COMPARISON TO NEW FIELD 1.
  235 ENDP(1)=KA(69)
      ENDP(2)=KA(70)
      READ(IOU5,240,END=310) KA
  240 FORMAT(A1,A4,A3,65A1,A4,A3)
!     LEFT-ADJUST FIELD 1.
      CALL LEFT(KA,IDUM)
!     CHECK IF CONTINUATION CARD.
      DO I=3,4
         IF(KA(1).EQ.MB(I)) GO TO 290
      end do
!     NO CONTINUATION.  END OF LOGICAL CARD FOUND.
  280 BACKSPACE IOU5
!     ELIMINATE ZEROS and DUPLICATE ENTRIES FROM LIST OF grids
      call fixit(kg,npoint)
      RETURN
!
!     CONTINUATION FOUND.
  290 CONTINUE
      IF(KA(2).NE.ENDP(1)) GO TO 280
      IF(KA(3).NE.ENDP(2)) GO TO 280
!     IS MATCHING CONTINUATION A LONG-FIELD CARD . . .
      IF(I.NE.5) GO TO 296
      WRITE(IOU6,293)
  293 FORMAT(/,'Fatal Error.  BANDIT cannot read a long field',                &
             ' continuation to a rigid element card.')
      WRITE(IOU6,295) KA
295   FORMAT(A1,A4,A3,65A1,A4,A3)
      call finish(6,IER)
      RETURN
  296 CONTINUE
!     WRITE CARD OUT ON UNIT 8.
      WRITE(IOU8,240)  KA
!     GO BACK AND INTERPRET THIS CARD.
      GO TO 90
!
!     END-OF-FILE ENCOUNTERED
!
310   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE RIGID

! ##################################################################################################################################
      SUBROUTINE RSETUP(LVL,LVLS1,LVLS2,NACUM,IDIM,IER)
!
!     SETUP COMPUTES THE REVERSE LEVELING INFO FROM LVLS2 AND STORES
!     IT INTO LVLS2.  NACUM(I) IS INITIALIZED TO NODES/ITH LEVEL FOR
!     NODES ON THE PSEUDO-DIAMETER OF THE GRAPH.  LVL IS INITIALIZED
!     TO NONZERO FOR NODES ON THE PSEUDO-DIAM AND NODES IN A DIFFERENT
!     COMPONENT OF THE GRAPH.
!
      COMMON /GRA/ N,IDPTH,DUMG
      INTEGER IDIM
      INTEGER NACUM(IDIM),LVL(*),LVLS1(*),LVLS2(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDPTH  ,IER    ,ITEMP  ,N

      REAL     DUMG

! E////////////////////////////////////////////////////////////////////E
!     IDIM=NUMBER OF LEVELS IN A GIVEN COMPONENT.
      IF(IDPTH.LE.IDIM)  GO TO 20
!     DIMENSION EXCEEDED  . . .  STOP JOB.
      CALL FINISH(3,IER)
      RETURN
   20 CONTINUE
      DO 30 I=1,IDPTH
        NACUM(I)=0
   30 CONTINUE
      DO 140 I=1,N
        LVL(I)=1
        LVLS2(I)=IDPTH+1-LVLS2(I)
        ITEMP=LVLS2(I)
        IF(ITEMP.GT.IDPTH) GO TO 140
        IF(ITEMP.NE.LVLS1(I)) GO TO 100
        NACUM(ITEMP)=NACUM(ITEMP)+1
        GO TO 140
  100   LVL(I)=0
  140 CONTINUE
      RETURN
      END SUBROUTINE RSETUP

! ##################################################################################################################################
      SUBROUTINE SCAT(KG,NCON,INV,II3,NORIG,IER)
!
!     THIS ROUTINE USES SCATTER SORT TECHNIQUES FOR EACH GRID POINT
!     ENCOUNTERED TO DETERMINE WHETHER OR NOT THE POINT HAS BEEN
!     SEEN BEFORE.  IF NOT, INV, NORIG, AND NN ARE UPDATED.
!
!     INV(1,I) CONTAINS AN ORIGINAL GRID POINT NUMBER
!     INV(2,I) CONTAINS THE INTERNAL NUMBER ASSIGNED TO IT (BEFORE
!              SORTING)
!
      INTEGER II3,INV(2,II3),NORIG(*),KG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IER    ,KMOD   ,LOC    ,                                &
               MAXDEG ,MAXGRD ,NCON   ,NN     ,NOLD

      REAL     DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /S/ NN,DUMS(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      IF(NCON.LT.1) RETURN
      DO 50 I=1,NCON
      NOLD=KG(I)
      IF(NOLD.EQ.0) GO TO 50
      LOC=NOLD-1
10    LOC=MOD(LOC,KMOD)+1
20    IF(INV(1,LOC).NE.0) GO TO 30
      INV(1,LOC)=NOLD
      NN=NN+1
      IF(NN.GT.MAXGRD) GO TO 60
      NORIG(NN)=NOLD
      INV(2,LOC)=NN
      GO TO 40
30    IF(INV(1,LOC).NE.NOLD) GO TO 10
40    KG(I)=INV(2,LOC)
50    CONTINUE
      RETURN
60    WRITE(IOU6,70) MAXGRD
70    FORMAT(/,'Fatal Error.  Grid point limit of',I9,' exceeded.')
      CALL FINISH(1,IER)
      RETURN
      END SUBROUTINE SCAT

! ##################################################################################################################################
      SUBROUTINE SEQGP(NORIG,ILD,NEW,JUMP)
!
! B////////////////////////////////////////////////////////////////////B
      USE IOUNT1, ONLY                :  WRT_LOG, SEQ
! E////////////////////////////////////////////////////////////////////E
!     WRITE SEQGP BULK DATA CARDS ON 7 and 8.
!
!     NN       = NUMBER OF GRID POINTS
!     NORIG(I) = ORIGINAL GRID POINT CORRESPONDING TO BANDIT INTERNAL
!                LABEL I
!     NEW(I)   = SCRATCH ARRAY
!     ILD(I)   = NEW RESEQUENCED LABEL CORRESPONDING TO BANDIT INTERNAL
!                LABEL I
!     JUMP     = 1 IF NO SEQGP CARDS ARE TO BE GENERATED (E.G., IF NO
!                IMPROVEMENT IN BW HAS RESULTED).
!
      INTEGER NORIG(*),ILD(*),NEW(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IADD   ,ILOOP  ,IPARAM ,ITIME  ,                        &
               J      ,JUMP   ,                                                &
               MAXW   ,MM     ,                                                &
               NAXIC  ,NBW    ,NLOOP  ,NN     ,NP     ,NRMS

      REAL     DMY    ,DUM    ,DUMD   ,DUMS   ,DUMW   ,RMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,DMY(5),IADD,DUMS,NAXIC
      COMMON /B/ IPARAM(20)
      COMMON /D/ OBW,NBW,OP,NP,DUMD(6)
      COMMON /W/ DUM(2),MAXW,RMS,DUMW(2)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER OBW,OP
      REWIND IOU7
! B////////////////////////////////////////////////////////////////////B
      REWIND SEQ
      READ (SEQ,'(1X,I11)') ITIME
! E////////////////////////////////////////////////////////////////////E
      IF(NN.LE.0) RETURN
      WRITE(IOU6,2)
2     FORMAT(/,'Field 10 of first SEQGP card shows new grid',                  &
             ' point rms wavefront.')
      IF(IADD.NE.0) WRITE(IOU6,3) IADD
    3 FORMAT(/,'Integer added to new sequence numbers',I8)
!
      NRMS=RMS+0.5
!     NRMS = NEW RMS WAVEFRONT ROUNDED TO NEAREST INTEGER
      IF(JUMP.NE.1) GO TO 9
!
      WRITE(IOU6,7)
    7 FORMAT(/,'No SEQGP cards generated.')
      REWIND IOU7
      ENDFILE IOU7
! B////////////////////////////////////////////////////////////////////B
      CLOSE (SEQ, STATUS='DELETE')
! E////////////////////////////////////////////////////////////////////E
      GO TO 60
!
    9 CONTINUE
      WRITE(IOU6,13)
13    FORMAT(/,'See bandit.f07 for SEQGP cards generated.')
!
!     IF SPRINGS ARE REQUESTED WITH $SPRING YES, THE SEQGP CARDS
!     RESEQUENCE THE CONSECUTIVELY NUMBERED OLD NUMBERS RATHER THAN
!     THE ORIGINAL GRID NUMBERS ASSIGNED BY THE USER.
!
      IF(IPARAM(7).EQ.3) GO TO 150
      WRITE(IOU6,125)
  125 FORMAT('Since scalar springs were requested with $SPRING YES,',          &
       ' the SEQGP cards resequence the consecutively numbered old'/           &
       'numbers rather than the user-assigned original grid number.'/)
      REWIND iou16
      READ(iou16,130) I,MM
  130 FORMAT(24I5)
!     DETERMINE REPLACEMENT NORIG ARRAY.
      DO I=1,NN
         READ(iou16,130) NORIG(I),(NEW(J),J=1,MM)
!        NEW IS DUMMY SPACE HERE TO INSURE THAT THE RIGHT NUMBER OF
!        CARDS IS READ.
      end do
      REWIND iou16
  150 CONTINUE
!
!     ADD THE NON-NEGATIVE INTEGER IADD TO THE NEW SEQUENCE NUMBERS.
!
      DO I=1,NN
         ILD(I)=ILD(I)+IADD
      end do
!
      NLOOP=0
      IF(NAXIC.NE.99999) NLOOP=IABS(NAXIC)
      DO 30 ILOOP=0,NLOOP
         IF(NAXIC.NE.99999) THEN
            DO I=1,NN
               NORIG(I)=NORIG(I)+1000000
               ILD(I)=ILD(I)+1000000
            end do
            IF(NAXIC.LT.0.AND.ILOOP.NE.NLOOP) GO TO 30
         END IF
         IF(NN.LE.4) THEN
            WRITE(IOU8,10) (NORIG(I),ILD(I),I=1,NN)
            WRITE(IOU7,10) (NORIG(I),ILD(I),I=1,NN)
! B////////////////////////////////////////////////////////////////////B
            WRITE(SEQ, 10) (NORIG(I),ILD(I),I=1,NN)
! E////////////////////////////////////////////////////////////////////E
         ELSE
            WRITE(IOU8,10) (NORIG(I),ILD(I),I=1,4),NRMS,                       &
                           (NORIG(J),ILD(J),J=5,NN)
            WRITE(IOU7,10) (NORIG(I),ILD(I),I=1,4),NRMS,                       &
                           (NORIG(J),ILD(J),J=5,NN)
! B////////////////////////////////////////////////////////////////////B
            WRITE(SEQ, 10) (NORIG(I),ILD(I),I=1,4),NRMS,                       &
                           (NORIG(J),ILD(J),J=5,NN)
! E////////////////////////////////////////////////////////////////////E
10       FORMAT('SEQGP*          ',8I16,I16/('SEQGP*          ',8I16))
         END IF
30    CONTINUE
!
60    WRITE(IOU8,70)
70    FORMAT('ENDDATA')
      REWIND IOU7
      RETURN
      END SUBROUTINE SEQGP

! ##################################################################################################################################
      SUBROUTINE SETIG(KG1,KG2,IG,II1,NORIG,IER)
!
!     THIS ROUTINE SETS IG(KG1,-)=KG2 AND IG(KG2,-)=KG1 IF THIS
!        CONNECTION HAS NOT ALREADY BEEN SET.
!
!     NEDGE = NUMBER OF UNIQUE EDGES.
      INTEGER II1
      DIMENSION IG(II1,*),NORIG(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  IER    ,IG     ,IS     ,                                        &
               K      ,KG1    ,KG2    ,KMOD   ,                                &
               L      ,LOOP   ,                                                &
               M      ,MAXDEG ,MAXGRD ,MM     ,                                &
               NBITIN ,NEDGE  ,NN     ,NORIG

      REAL     DUM    ,DUMBB  ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUM(3),NEDGE,DUMS(3)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!
      IF(KG1.EQ.0 .OR. KG2.EQ.0 .OR. KG1.EQ.KG2) RETURN
      L=KG1
      K=KG2
      DO 50 LOOP=1,2
      IF(LOOP.EQ.1) GO TO 20
      L=KG2
      K=KG1
   20 M=0
   30 M=M+1
      IF(M.GT.MAXDEG) GO TO 60
      IS=IG(L,M)
!PACK IS=IUNPK(IG,MAXGRD*(M-1)+L,NBITIN)
      IF(IS.EQ.0) GO TO 40
      IF(IS.NE.K) GO TO 30
      GO TO 55
40    CONTINUE
      IG(L,M)=K
!PACK CALL PACK(IG,MAXGRD*(M-1)+L,NBITIN,K)
      MM=MAX(MM,M)
      IF(LOOP.EQ.1) NEDGE = NEDGE + 1
   50 CONTINUE
   55 RETURN
!
! B////////////////////////////////////////////////////////////////////B
   60 WRITE(IOU6,70) NORIG(L),M,MAXDEG
70    FORMAT(/,'Fatal Error.  The number of connections for point',I9/        &
       ' = ',I9,' and exceeds the nodal degree limit of',I8)
! E////////////////////////////////////////////////////////////////////E
      CALL FINISH(1,IER)
      RETURN
      END SUBROUTINE SETIG

! ##################################################################################################################################
      SUBROUTINE SETUP(KA,NCON,IFLD,KG,NPT,LEN,LESSOK,IER)
!
!     SET UP ARRAY KG CONTAINING THE LIST OF CONNECTIONS FOR AN ELEMENT.
!
!     DATA CARD ARRAY KA WAS READ A1,A4,A3,64A1,A1,A4,A3
!
!     NCON   = NUMBER OF CONNECTIONS (INPUT)
!     IFLD   = FIELD NUMBER OF FIRST CONNECTION (INPUT)
!     KG     = LIST OF CONNECTIONS (OUTPUT)
!     NPT    = ACTUAL NUMBER OF NON-ZERO CONNECTIONS FOUND (OUTPUT)
!     LEN    = 1 FOR SHORT FIELD DATA CARD, 2 FOR  LONG FIELD DATA CARD
!              (INPUT)
!     LESSOK = .TRUE. IF THE CONNECTION CARD NEED NOT HAVE ALL GRID
!              POINTS PRESENT (E.G., CELAS1 OR CPENTA), OTHERWISE FALSE
!              (INPUT)
!
      INTEGER KA(70),KG(*),ENDP(2)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,      IER    ,IFLD   ,IPARAM ,ITYPE  ,                  &
               J      ,K      ,L      ,LEN    ,LEN8   ,                        &
               M      ,M1     ,MA     ,MB     ,                                &
               N      ,NCON   ,NIP    ,NPT    ,NUM

! E////////////////////////////////////////////////////////////////////E
      COMMON /ALPHA/ MA(26),NUM(10),MB(4)
      COMMON /B/ IPARAM(20)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      LOGICAL LESSOK
!
      DO I=1,NCON
         KG(I)=0
      end do
      N=0
      NPT=0
      LEN8=LEN*8
      J = 4 + 8 * (IFLD- 2 ) * LEN
      L = 14 - IFLD - 4 * LEN
!
!     N=NO. OF CONNECTION FIELDS PROCESSED SO FAR
!     J=SUBSCRIPT IN KA CORRESPONDING TO FIRST CONNECTION FIELD ON CARD
!     L=NO. OF FIELDS FROM FIRST CONNECTION FIELD ON CARD TO END OF CARD
!     K=MAX. NO. OF CONNECTIONS POSSIBLE ON A GIVEN PHYSICAL CARD
!
20    K=MIN(L,NCON-N)
      DO M=1,K
         M1=J+LEN8*(M-1)
         I=NPT+1
         CALL READIT(KA(M1),1,LEN8,KG(I),NIP)
         IF(KG(I).NE.0) NPT=NPT+1
      end do
      N=N+K
      IF(N.GE.NCON) RETURN
!     SAVE ENDPUNCHING (FIELD 10).
      ENDP(1)=KA(69)
      ENDP(2)=KA(70)
!     READ NEXT CARD.
      READ(IOU5,40,END=100) KA
40    FORMAT(A1,A4,A3,65A1,A4,A3)
!     LEFT-ADJUST FIELD 1.
      CALL LEFT(KA,ITYPE)
!     CHECK NEW FIELD 1 AGAINST PREVIOUS FIELD 10 FOR MATCH.
      DO I=3,4
         IF(KA(1).EQ.MB(I)) GO TO 60
      end do
      GO TO 70
60    IF(KA(2).NE.ENDP(1)) GO TO 70
      IF(KA(3).NE.ENDP(2)) GO TO 70
      WRITE(IOU8,40) KA
      LEN=1
      IF(I.EQ.5) LEN=2
      L=8/LEN
      J=4
      GO TO 20
!     END OF LOOP STARTING AT STATEMENT 20.
!
!     ABORT JOB IF CONNECTION CARD CONTINUATIONS DO NOT IMMEDIATELY
!        FOLLOW THEIR PARENTS, UNLESS IT IS FOR AN ELEMENT FOR WHICH
!        IT IS OK IF NOT ALL CONNECTIONS ARE PRESENT (E.G., CELAS1
!        OR CPENTA)
!
70    CONTINUE
      IF(LESSOK) GO TO 90
      WRITE(IOU6,80) KA
80    FORMAT(/,'Fatal Error.  The following data card is out of sort.'/        &
       A1,A4,A3,65A1,A4,A3/                                                    &
       'Since BANDIT does not sort the data, continuations must'/              &
       'immediately follow their parents.')
      call finish(6,IER)
      RETURN
!
90    CONTINUE
      BACKSPACE IOU5
      RETURN
!
!     END-OF-FILE ENCOUNTERED
!
100   CALL FINISH(4,IER)
      RETURN
      END SUBROUTINE SETUP

! ##################################################################################################################################
      SUBROUTINE SORT(LIST,NL)
!
!     SORT A LIST OF INTEGERS OF LENGTH NL.  THIS ROUTINE OPERATES
!        FASTEST FOR THOSE LISTS NOT BADLY OUT OF SORT.
!
      INTEGER LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,J1     ,K      ,KFLAG  ,L      ,                &
               NL     ,NL1

! E////////////////////////////////////////////////////////////////////E
      IF(NL.LE.1) RETURN
      NL1=NL-1
      DO 20 I=1,NL1
      K=NL-I
      KFLAG=0
      DO 10 J=1,K
      J1=J+1
      IF(LIST(J).LE.LIST(J1))  GO TO 10
      KFLAG=1
      L=LIST(J)
      LIST(J)=LIST(J1)
      LIST(J1)=L
10    CONTINUE
      IF(KFLAG.EQ.0) RETURN
20    CONTINUE
      RETURN
      END SUBROUTINE SORT

! ##################################################################################################################################
      SUBROUTINE SORTDG(STK1,STK2,X1,X2,NDEG)
!
!     SORTDG SORTS STK2 BY DEGREE OF THE NODE AND ADDS IT TO THE END
!     OF STK1 IN ORDER OF LOWEST TO HIGHEST DEGREE.  X1 AND X2 ARE THE
!     NUMBER OF NODES IN STK1 AND STK2 RESPECTIVELY.
!
      COMMON /GRA/ N,IDPTH,DUMG
      INTEGER NDEG(*),STK1(*),STK2(*),X1,X2,TEMP
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IDPTH  ,IND    ,ISTK2  ,ITEST  ,                        &
               J      ,JSTK2  ,N

      REAL     DUMG

! E////////////////////////////////////////////////////////////////////E
      IND=X2
   10 ITEST=0
      IND=IND-1
      IF(IND.LT.1) GO TO 40
      DO 30 I=1,IND
        J=I+1
        ISTK2=STK2(I)
        JSTK2=STK2(J)
        IF(NDEG(ISTK2).LE.NDEG(JSTK2)) GO TO 30
        ITEST=1
        TEMP=STK2(I)
        STK2(I)=STK2(J)
        STK2(J)=TEMP
   30 CONTINUE
      IF(ITEST.EQ.1) GO TO 10
   40 DO 50 I=1,X2
        X1=X1+1
        STK1(X1)=STK2(I)
   50 CONTINUE
      RETURN
      END SUBROUTINE SORTDG

! ##################################################################################################################################
      INTEGER FUNCTION SORT2(XC,SIZE,STPT)
!
!     SORT2 SORTS SIZE AND STPT INTO DESENDING ORDER ACCORDING TO
!     VALUES OF SIZE.
!
!     XC=NUMBER OF ENTRIES IN EACH ARRAY
!
      INTEGER SIZE(*),STPT(*),TEMP,XC
! B////////////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IND    ,ITEST  ,J

! E////////////////////////////////////////////////////////////////////E
!     THE DIMENSION OF SIZE AND STPT IS THE MAXIMUM ALLOWABLE NUMBER
!            OF CONNECTED COMPONENTS.
      SORT2=0
      IF(XC.EQ.0) RETURN
      SORT2=1
      IND=XC
   10 ITEST=0
      IND=IND-1
      IF(IND.LT.1) RETURN
      DO 17 I=1,IND
        J=I+1
        IF(SIZE(I).GE.SIZE(J)) GO TO 17
        ITEST=1
        TEMP=SIZE(I)
        SIZE(I)=SIZE(J)
        SIZE(J)=TEMP
        TEMP=STPT(I)
        STPT(I)=STPT(J)
        STPT(J)=TEMP
   17 CONTINUE
      IF(ITEST.EQ.1) GO TO 10
      RETURN
      END FUNCTION SORT2

! ##################################################################################################################################
      SUBROUTINE SPRING(IP)
!
!     GENERATE SCALAR SPRING ELEMENTS CONSISTENT WITH THE CONNECTIVITY
!     MATRIX IG.  CELAS3 ELEMENTS ARE WRITTEN ON UNIT 9.
!     IP = TEMPORARY STORAGE.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,II     ,J      ,MM     ,NN

      REAL     STIFF

! E////////////////////////////////////////////////////////////////////E
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER IP(*),EID,PID
      PID=101
      STIFF=2.71828
      REWIND IOU9
      REWIND iou16
      EID=0
      READ(iou16,10) NN,MM
      IF((NN*MM).EQ.0) GO TO 55
!
!     GENERATE SPRINGS BETWEEN NODES.
!
      DO 25 II=1,NN
         READ(iou16,10) I,(IP(J),J=1,MM)
10       FORMAT(24I5)
         DO 20 J=1,MM
            IF(IP(J).LE.0) GO TO 25
            IF(IP(J).LT.I) GO TO 20
            EID=EID+1
            WRITE(IOU9,15) EID,PID,I,IP(J)
15          FORMAT('CELAS3  ',4I8)
20       CONTINUE
25    CONTINUE
!
!     INSERT BLANK CARD FOR LEVY'S WAVEFRONT PROGRAM.
!
      WRITE(IOU9,'(1x)')
!
!     GENERATE SPRINGS TO GROUND.
!
      PID=PID+1
      DO I=1,NN
         EID=EID+1
         WRITE(IOU9,35) EID,PID,I
35       FORMAT('CELAS3  ',3I8)
      end do
!
!     CREATE PROPERTY CARD.
!
      PID=PID-1
      WRITE(IOU9,50) PID,STIFF,NN,MM
   50 FORMAT('PELAS   ',I8,F8.5,48X,2I4)
      PID=PID+1
      WRITE(IOU9,50) PID,STIFF
!
   55 CONTINUE
      REWIND IOU9
      REWIND iou16
      WRITE(IOU6,60) EID,IOU9
   60 FORMAT(/I8,' CELAS3 spring elements generated on bandit.f',I2)
      RETURN
      END SUBROUTINE SPRING

! ##################################################################################################################################
      SUBROUTINE STABLE(IG,II1,IC,IDEG,ILD,IP)
!
!     PRINT CONNECTION TABLE IN TERMS OF SORTED INTERNAL LABELS.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IC     ,IDEG   ,IG     ,II1    ,ILD    ,IP     ,        &
               J      ,K      ,KMOD   ,                                        &
               MAXDEG ,MAXGRD ,MDIFF  ,MM     ,                                &
               NBITIN ,NN

      REAL     DUM    ,DUMBB

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /S/ NN,MM,DUM(7)
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      DIMENSION IG(II1,*),IC(*),IDEG(*),ILD(*),IP(*)
      IF((NN*MM).EQ.0) RETURN
      WRITE(IOU6,10)
   10 FORMAT(/,'Sorted Internal Label Connection Table:'/                      &
       8X,'Sort       Max'/' Label Label Comp Diff Degr  Connections')
      DO 40 I=1,NN
      DO 20 J=1,MM
   20 IP(J)=0
      MDIFF=0
      DO 30 J=1,MM
      K=IG(I,J)
!PACK K=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
      IF(K.EQ.0) GO TO 40
      MDIFF=MAX(MDIFF,IABS(ILD(I)-ILD(K)))
   30 IP(J)=ILD(K)
   40 WRITE(IOU6,50) I,ILD(I),IC(I),MDIFF,IDEG(I),(IP(J),J=1,MM)
   50 FORMAT(2I6,24I5/(27X,21I5))
      RETURN
      END SUBROUTINE STABLE

! ##################################################################################################################################
      SUBROUTINE STACK(IDEG,NEW,ILD,IW)
!
!     STACK POINTS OF ZERO DEGREE AT END OF THE NUMBERING.
!
      INTEGER IDEG(*),NEW(*),ILD(*),IW(*)
!     IW IS SCRATCH STORAGE.
      COMMON /S/ NN,DUMS(8)
      COMMON /D/ DUM(5),KT,DUMD(4)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,J      ,K      ,KT     ,L      ,NN     ,NN1

      REAL     DUM    ,DUMD   ,DUMS

! E////////////////////////////////////////////////////////////////////E
      KT=0
      NN1=NN-1
!     LIST POINTS OF ZERO DEGREE AND INCREMENT COUNTER KT.
      DO 10 I=1,NN
      IF(IDEG(I).GT.0) GO TO 10
      KT=KT+1
      IW(KT)=ILD(I)
   10 CONTINUE
      IF(KT.LE.0) GO TO 70
!     SORT LIST OF RENUMBERED NUMBERS TO BE STACKED.
      CALL SORT(IW,KT)
!     STACK POINTS OF ZERO DEGREE AT END OF NEW.
      DO 40 L=1,KT
      I=IW(L)-L+1
      K=NEW(I)
      IF(I.GE.NN) GO TO 30
      DO 20 J=I,NN1
   20 NEW(J)=NEW(J+1)
   30 NEW(NN)=K
   40 CONTINUE
!     CORRECT ILD, THE INVERSE OF NEW.
   70 DO 80 I=1,NN
      K=NEW(I)
   80 ILD(K)=I
      RETURN
      END SUBROUTINE STACK

! ##################################################################################################################################
! B////////////////////////////////////////////////////////////////////B
      SUBROUTINE SUMUP ( NEW_BW, DEN )
! E////////////////////////////////////////////////////////////////////E
!
!     PRINT BANDIT SUMMARY.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  IADD   ,IPARAM ,                                                &
               MAXW0  ,MAXW1  ,MINDEG ,MM     ,                                &
               NBW    ,NCM    ,NEDGE  ,NEL    ,NEQ    ,NEQR   ,NLINK  ,        &
               NN     ,NONZ   ,NP     ,NZERO

! B////////////////////////////////////////////////////////////////////B
      INTEGER  NEW_BW
! E////////////////////////////////////////////////////////////////////E
      REAL     AN     ,ANN    ,AV1    ,AV2    ,                                &
               BRMS0  ,BRMS1  ,                                                &
               DEN    ,DUM    ,DUMY   ,                                        &
               rawf   ,rbw    ,RMS0   ,RMS1   ,rmwf   ,rp     ,                &
               rrbw   ,rrwf

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUM(3),NEDGE,IADD,MINDEG,DUMY
      COMMON /B/ IPARAM(20)
      COMMON /D/ OBW,NBW,OP,NP,NCM,NZERO,NEL,NEQ,NEQR,NLINK
      COMMON /W/ MAXW0,RMS0,MAXW1,RMS1,BRMS0,BRMS1
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      INTEGER OBW,OP
! B////////////////////////////////////////////////////////////////////B
      NEW_BW = NBW         ! New bandwidth printed out in this subr
! E////////////////////////////////////////////////////////////////////E
      ANN=FLOAT(NN)
      AV1=FLOAT(OP)/ANN
      AV2=FLOAT(NP)/ANN
      NONZ=2*NEDGE+NN
      AN=ANN*ANN
      DEN=FLOAT(NONZ)*100./AN
      NEQ=NEQ-NEQR
      rbw=float(nbw)/float(obw)
      rp=float(np)/float(op)
      rmwf=float(maxw1)/float(maxw0)
      rawf=av2/av1
      rrwf=rms1/rms0
      rrbw=brms1/brms0
      WRITE(IOU6,20) OBW,NBW,rbw,OP,NP,rp,MAXW0,MAXW1,rmwf,                    &
                     AV1,AV2,rawf,RMS0,RMS1,rrwf,BRMS0,BRMS1,rrbw
   20 FORMAT(/,'BANDIT Summary:'/39X,'Before',7X,'After',5x,'Ratio'/           &
       5X,'Bandwidth (B)',15X,2I12,f10.3/                                      &
       5X,'Profile (P)', 17X,2I12,f10.3/                                       &
       5X,'Maximum Wavefront (C-MAX)',3X,2I12,f10.3/                           &
       5X,'Average Wavefront (C-AVG)',3X,2F12.3,f10.3/                         &
       5X,'RMS Wavefront (C-RMS)',7X,2F12.3,f10.3/                             &
       5X,'RMS Bandwidth (B-RMS)',7X,2F12.3,f10.3/)
      IF(IPARAM(6).EQ.10) WRITE(IOU6,24)
      IF(IPARAM(6).EQ.11) WRITE(IOU6,26)
      IF(IPARAM(6).EQ.12) WRITE(IOU6,28)
      IF(IPARAM(6).EQ.13) WRITE(IOU6,30)
  30  FORMAT(5X,'Criterion',30X,'Max Wavefront')
  24  FORMAT(5X,'Criterion',34X,'Bandwidth')
  26  FORMAT(5X,'Criterion',36X,'Profile')
  28  FORMAT(5X,'Criterion',30X,'RMS Wavefront')
      IF(IPARAM(3).EQ.7) WRITE(IOU6,32)
      IF(IPARAM(3).EQ.8) WRITE(IOU6,34)
      IF(IPARAM(3).EQ.9) WRITE(IOU6,36)
   36 FORMAT(5X,'Method Selected',27X,'CM and GPS')
   32 FORMAT(5X,'Method Selected',35X,'CM')
   34 FORMAT(5X,'Method Selected',34X,'GPS')
      WRITE(IOU6,60) NN,NEL,NLINK,NCM,MM,MINDEG,NEDGE,DEN,NZERO,NEQR,NEQ
   60 FORMAT(5X,'Number of Grid Points (N)',I27/                               &
             5X,'Number of Elements (Non-Rigid)',I22/                          &
             5X,'Number of Rigid Elements Processed',I18/                      &
             5X,'Number of Components',I32/                                    &
             5X,'Maximum Nodal Degree',I32/                                    &
             5X,'Minimum Nodal Degree',I32/                                    &
             5X,'Number of Unique Edges',I30/                                  &
             5X,'Matrix Density in Percent',F27.4/                             &
             5X,'Number of Points of Zero Degree',I21/                         &
             5X,'Number of Rigid Element MPC Equations',I15/                   &
             5X,'Number of MPC Equations Processed',I19)
      WRITE(IOU6,80)
80    FORMAT(/,'All BANDIT statistics use grid point, rather than',            &
       ' DOF, connectivity and'/                                               &
       'include matrix diagonal terms.  Statistics such as',                   &
       ' C-MAX, C-AVG, C-RMS,'/                                                &
       'and N should each be multiplied by the average number',                &
       ' of DOF per grid point'/                                               &
       'before estimating NASTRAN time and core requirements.')
      RETURN
      END SUBROUTINE SUMUP

! ##################################################################################################################################
      SUBROUTINE TABLE2(IG,II1,NORIG,IP)
!
!     PRINT CONNECTION TABLE IN TERMS OF ORIGINAL GRID NUMBERS.
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IG     ,II1    ,IP     ,J      ,K      ,KMOD   ,        &
               MAXDEG ,MAXGRD ,MM     ,NBITIN ,NN     ,NORIG

      REAL     DUM    ,DUMBB 

! E////////////////////////////////////////////////////////////////////E
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /S/ NN,MM,DUM(7)
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
      DIMENSION IG(II1,*),IP(*),NORIG(*)
      IF((NN*MM).EQ.0) RETURN
      WRITE(IOU6,10)
   10 FORMAT(/,'Original Grid Number Connection Table:'/                       &
             ' Label  Grid ID     Connections')
      DO 30 I=1,NN
      DO J=1,MM
         IP(J)=0
      end do
      DO 20 J=1,MM
      K=IG(I,J)
!PACK K=IUNPK(IG,MAXGRD*(J-1)+I,NBITIN)
      IF(K.EQ.0) GO TO 30
   20 IP(J)=NORIG(K)
   30 WRITE(IOU6,40) I,NORIG(I),(IP(J),J=1,MM)
   40 FORMAT(I6,14I9/(15X,13I9))
      RETURN
      END SUBROUTINE TABLE2

! ##################################################################################################################################
      SUBROUTINE TAXI(KA,IOU6)
!
!     CHECK FOR AXIC BULK DATA CARD, AND, IF PRESENT, READ NUMBER OF
!     HARMONICS.  NEGATIVE NUMBER OF HARMONICS IMPLIES SINGLE HARMONIC
!     PROBLEM.  IT IS ASSUMED THAT THE MNEUMONIC IS LEFT-ADJUSTED.
!
! B////////////////////////////////////////////////////////////////////B
! Add this so when READIT is called with NAXIC we will use NAXIC_array
! instead. Needed so Lahey doesn't complain about shape of NAXIC being
! different than IP (array) in subr READIT
      INTEGER NAXIC_array(1)
! E////////////////////////////////////////////////////////////////////E
      INTEGER KA(*),A,XIC,BLANK,MINUS
! B////////////////////////////////////////////////////////////////////B
      DATA A,XIC,BLANK,MINUS/1HA,3HXIC,1H ,1H-/
! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ DUM(8),NAXIC
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IOU6   ,NAXIC  ,NIP

      REAL     DUM

! E////////////////////////////////////////////////////////////////////E
!
!     CHECK FOR AXIC
!
      IF(KA(1).NE.A  ) RETURN
      IF(KA(2).NE.XIC) RETURN
!
!     LOOK FOR FIRST NON-BLANK IN FIELD 2 OF AXIC CARD
!
      DO 10 I=4,11
      IF(KA(I).NE.BLANK) GO TO 20
10    CONTINUE
      RETURN
!
!     READ AXIC INTEGER AND WRITE MESSAGE
!
! B////////////////////////////////////////////////////////////////////B
20    CALL READIT(KA(I),1,12-I,NAXIC_array(1),NIP)
      NAXIC = NAXIC_array(1)
! E////////////////////////////////////////////////////////////////////E
      IF(KA(I).EQ.MINUS) NAXIC=-NAXIC
      WRITE(IOU6,30) NAXIC
30    FORMAT(/,'AXIC',I8,' card read.')
!
      RETURN
      END SUBROUTINE TAXI

! ##################################################################################################################################
      SUBROUTINE TIGER(IG,II1,LIST,NORIG,KG,MAXI,IER)
!
!     THIS ROUTINE MAKES ADDITIONS TO THE CONNECTION TABLE IG TO REFLECT
!     THE PRESENCE OF MPC*S AND STORES THE DEPENDENT POINTS IN LIST.
!
!     NEQ = NUMBER OF MPC EQUATIONS.
!
      INTEGER II1, MAXI
      DIMENSION IG(II1,*),LIST(*),NORIG(*)
      INTEGER KG(MAXI)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOU5   ,IOU6   ,IOU7   ,IOU8   ,IOU9   ,IOU10  ,IOU11  ,        &
               IOU12  ,IOU13  ,IOU14  ,IOU15  ,IOU16  ,IOU17  ,IOU18  ,        &
               IOU19  ,IOU20

      INTEGER  I      ,IER    ,IG     ,IGRID  ,II     ,IOPT   ,                &
               J      ,KMOD   ,L      ,LIST   ,                                &
               MAXDEG ,MAXGRD ,MM     ,                                        &
               NBITIN ,NEQ    ,NN     ,NORIG  ,NTERM

      REAL     DUM    ,DUMBB  ,DUMD   ,DUMS

! E////////////////////////////////////////////////////////////////////E
      COMMON /S/ NN,MM,DUMS(7)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /D/ DUM(7),NEQ,DUMD(2)
      COMMON /BITS/ NBITIN,DUMBB(8)
      COMMON /IOUNIT/ IOU5,IOU6,IOU7,IOU8,IOU9,IOU10,IOU11,IOU12,IOU13,        &
                      IOU14,IOU15,IOU16,IOU17,IOU18,IOU19,IOU20
!     IOPT=1 IF DEPENDENT POINTS ARE TO BE ACCUMULATED IN A LIST FOR
!            LATER REMOVAL FROM THE CONNECTION TABLE IG (OTHERWISE 0)
      DATA IOPT/0/
      IF(NEQ.EQ.0)RETURN
      REWIND IOU12
!     GENERATE NEW CONNECTIONS.
      DO 30 II=1,NEQ
!     READ MPC EQUATION.
      READ(IOU12) NTERM,(KG(I),I=1,NTERM)
!     IGRID=DEPENDENT GRID POINT IN AN MPC EQUATION.
      IGRID=KG(1)
      IF(IOPT.EQ.1) LIST(IGRID)=IGRID
      DO 20 I=1,MAXDEG
         L=IG(IGRID,I)
!PACK    L=IUNPK(IG,MAXGRD*(I-1)+IGRID,NBITIN)
!        L=A GRID POINT THAT IGRID IS CONNECTED TO BEFORE THE MPC IS APPLIED
         IF(L.LE.0) GO TO 30
         IF(NTERM.LT.2) GO TO 20
         DO J=2,NTERM
            CALL SETIG(L,KG(J),IG,II1,NORIG,IER)
            IF(IER.GT.0) RETURN
         end do
20    CONTINUE
   30 CONTINUE
      REWIND IOU12
!     ENDFILE IOU12
!     REWIND IOU12
      RETURN
      END SUBROUTINE TIGER

! ##################################################################################################################################
      SUBROUTINE TIMER(T,IWALL,IOPT,IOU6)
!
!     RETURN IN T THE CURRENT VALUE OF THE CP CLOCK IN SECONDS.
!     RETURN IN IWALL THE ELAPSED TIME IN INTEGER SECONDS.
!     IOPT = 1 IF PRINTOUT OF TIME IS DESIRED, OTHERWISE 0 (INPUT)
!
!     real*4 tarray(2)
!     integer*4 iwall,time
!
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IOPT   ,IOU6   ,IWALL

      REAL     T

! E////////////////////////////////////////////////////////////////////E
      t=0.
      iwall=0
!-----------------------------------------------------------------------
!     Cray:
!     call second(t)
!-----------------------------------------------------------------------
!     HP workstation:
!     etime on the HP requires the +U77 Fortran compile-line option.
!     The etime function value is the sum of user and system time,
!     where tarray(1)=user time, tarray(2)=system time.
!     To use these calls, activate the next 3 lines and the real*4 and
!     integer*4 lines above.
!     t=etime(tarray)
!     t=tarray(1)
!     iwall=time(0.)
!-----------------------------------------------------------------------
!     PC using MS Fortran Powerstation or Compaq Visual Fortran:
!     t=secnds(0.0)
!-----------------------------------------------------------------------
!     PC using G77 compiler:
!     call second(t)
!-----------------------------------------------------------------------
!     PC using Lahey compiler or G77 compiler:
!     G77's 'cpu_time' is an alias for 'second', so either can be used.
      call cpu_time(t)
!-----------------------------------------------------------------------
      IF(IOPT.GT.0) WRITE(IOU6,10) T
10    FORMAT('CP time =',F12.3,' seconds')
      RETURN
      END SUBROUTINE TIMER

! ##################################################################################################################################
      SUBROUTINE TREE(IROOT,NDSTK,NR,LVL,IWK,NDEG,LVLWTH,LVLBOT,               &
                      LVLN,MAXLW,IBORT)
!
!  TREE DROPS A TREE IN NDSTK FROM IROOT
!
!  LVL-         ARRAY INDICATING AVAILABLE NODES IN NDSTK WITH ZERO
!               ENTRIES. TREE ENTERS LEVEL NUMBERS ASSIGNED
!               DURING EXECUTION OF OF THIS PROCEDURE
!  IWK-         ON OUTPUT CONTAINS NODE NUMBERS USED IN TREE
!               ARRANGED BY LEVELS (IWK(LVLN) CONTAINS IROOT
!               AND IWK(LVLBOT+LVLWTH-1) CONTAINS LAST NODE ENTERED)
!  LVLWTH-      ON OUTPUT CONTAINS WIDTH OF LAST LEVEL
!  LVLBOT-      ON OUTPUT CONTAINS INDEX INTO IWK OF FIRST
!               NODE IN LAST LEVEL
!  MAXLW-       ON OUTPUT CONTAINS THE MAXIMUM LEVEL WIDTH
!  LVLN-        ON INPUT THE FIRST AVAILABLE LOCATION IN IWK
!               USUALLY ONE BUT IF IWK IS USED TO STORE PREVIOUS
!               CONNECTED COMPONENTS, LVLN IS NEXT AVAILABLE LOCATION.
!               ON OUTPUT THE TOTAL NUMBER OF LEVELS + 1
!  IBORT-       INPUT PARAM WHICH TRIGGERS EARLY RETURN IF
!               MAXLW BECOMES .GE. IBORT
!
      INTEGER NR
      DIMENSION NDSTK(NR,*),LVL(*),IWK(*),NDEG(*)
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  IBORT  ,INOW   ,IROOT  ,ITEST  ,ITOP   ,IWK    ,IWKNOW ,        &
               J      ,KMOD   ,                                                &
               LVL    ,LVLBOT ,LVLN   ,LVLTOP ,LVLWTH ,                        &
               MAXDEG ,MAXGRD ,MAXLW  ,                                        &
               NBITIN ,NDEG   ,NDROW  ,NDSTK

      REAL     DUMBB 

! E////////////////////////////////////////////////////////////////////E
      MAXLW=0
      ITOP=LVLN
      INOW=LVLN
      LVLBOT=LVLN
      LVLTOP=LVLN+1
      LVLN=1
      LVL(IROOT)=1
      IWK(ITOP)=IROOT
   30 LVLN=LVLN+1
   35 IWKNOW=IWK(INOW)
      NDROW=NDEG(IWKNOW)
      DO 40 J=1,NDROW
        ITEST=NDSTK(IWKNOW,J)
!PACK   ITEST=IUNPK(NDSTK,MAXGRD*(J-1)+IWKNOW,NBITIN)
        IF(LVL(ITEST).NE.0) GO TO 40
        LVL(ITEST)=LVLN
        ITOP=ITOP+1
        IWK(ITOP)=ITEST
   40 CONTINUE
      INOW=INOW+1
      IF(INOW.LT.LVLTOP) GO TO 35
      LVLWTH=LVLTOP-LVLBOT
      IF(MAXLW.LT.LVLWTH) MAXLW=LVLWTH
      IF(MAXLW.GE.IBORT) RETURN
      IF(ITOP.LT.LVLTOP) RETURN
      LVLBOT=INOW
      LVLTOP=ITOP+1
      GO TO 30
      END SUBROUTINE TREE

! ##################################################################################################################################
      SUBROUTINE WAVEY(IG,II1,ILD,NEW,NC,IC,KACT,MAXB,MAXW,AVERW,SUMW,         &
                       RMS,BRMS)
!
!     COMPUTE WAVEFRONT AND ACTIVE COLUMN DATA - -
!     MAXIMUM WAVEFRONT, AVERAGE WAVEFRONT, SUM OF ROW WAVEFRONTS,
!     SUM OF SQUARES OF ROW WAVEFRONTS, RMS WAVEFRONT, BANDWIDTH,
!     RMS BANDWIDTH, AND MINIMUM NODAL DEGREE.
!     DIAGONAL TERMS ARE INCLUDED.
!
!     IG      = CONNECTION TABLE
!     II1     = ROW DIMENSION OF IG
!     ILD(I)  = NEW LABEL FOR NODE WITH ORIGINAL INTERNAL LABEL I
!     NEW(I)  = INTERNAL LABEL CORRESPONDING TO NEW LABEL I
!               NEW AND ILD ARE INVERSES OF EACH OTHER
!     NC      = COMPONENT ID
!               IF NC.LE.0, USE ALL COMPONENTS.
!     IC(I)   = COMPONENT INDEX FOR ORIGINAL NODE I.
!     KACT(I) = LIST OF ACTIVE COLUMN FLAGS (UPDATED FOR EACH ROW)
!                 = 1 IF COL I IS ACTIVE AT GIVEN ROW
!               KACT IS SCRATCH SPACE(TEMPORARY STORAGE)
!     MAXB    = BANDWIDTH
!     MAXW    = MAXIMUM WAVEFRONT
!     AVERW   = AVERAGE WAVEFRONT
!     SUMW    = SUM OF ROW WAVEFRONTS
!     SUMSQ   = SUM OF SQUARES OF ROW WAVEFRONTS
!     BSUMSQ  = SUM OF SQUARES OF ROW BANDWIDTHS
!     RMS     = RMS WAVEFRONT
!     BRMS    = RMS BANDWIDTH
!     MINDEG  = MINIMUM NODAL DEGREE
!     NN      = NUMBER OF NODES
!     MM      = MAX NODAL DEGREE
!     MAXGRD  = EFFECTIVE ROW DIMENSION OF IG
!     NBITIN  = NUMBER OF BITS PER INTEGER(CDC)
!     INPUT   - IG,II1,ILD,NN,MM,MAXGRD,NBITIN,NC,IC
!     OUTPUT  - NEW,KACT,MAXW,AVERW,SUMW,RMS,MAXB,BRMS,MINDEG
!
      COMMON /S/ NN,MM,DUMS(5),MINDEG,DUMY
      COMMON /A/ MAXGRD,MAXDEG,KMOD
      COMMON /BITS/ NBITIN,DUMBB(8)
      INTEGER II1
      DIMENSION IG(II1,*),ILD(*),NEW(*),KACT(*)
      INTEGER IC(*),SUMW
      DOUBLE PRECISION SUMSQ,BSUMSQ
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,IB     ,IB1    ,IG     ,ILD    ,IWAVE  ,                &
               J      ,K      ,KACT   ,KMOD   ,KT     ,L      ,                &
               M      ,MAXB   ,MAXDEG ,MAXGRD ,MAXW   ,MINDEG ,MM     ,        &
               NBITIN ,NC     ,NEW    ,NN

      REAL     ANN    ,AVERW  ,BRMS   ,DUMBB  ,DUMS   ,DUMY   ,                &
               RMS

! E////////////////////////////////////////////////////////////////////E
!
!     INITIALIZE WAVEFRONT DATA.
!
      MAXB=0
      MAXW=0
      SUMW=0
      SUMSQ=0.D0
      BSUMSQ=0.D0
      AVERW=0.
      RMS=0.
      BRMS=0.
      MINDEG=MIN(MINDEG,MM)
      IF((NN*MM).LE.0) RETURN
!
!     INITIALIZE NEW, THE INVERSE OF ILD
!
      IF(NC.GT.0) GO TO 8
      DO 5 I=1,NN
      K=ILD(I)
      IF(K.LE.0) GO TO 5
      NEW(K)=I
    5 CONTINUE
    8 CONTINUE
!
!     INITIALIZE ACTIVE COLUMN FLAGS (1 FOR ACTIVE)
!
      DO 10 I=1,NN
   10 KACT(I)=0
!
!     COMPUTE WAVEFRONT DATA.
!
      IWAVE=1
      KT=0
      DO 40 I=1,NN
!     COMPUTE NUMBER OF ACTIVE COLUMNS FOR ROW I
      K=NEW(I)
      IF(NC) 18,18,15
   15 IF(K.LE.0) GO TO 40
      IF(NC-IC(K)) 40,18,40
18    KT=KT+1
      IB=0
      DO 28 J=1,MM
      L=IG(K,J)
!PACK L=IUNPK(IG,MAXGRD*(J-1)+K,NBITIN)
      IF(L.EQ.0) GO TO 30
      M=ILD(L)
      IB=MAX(IB,I-M)
      IF(M.LE.I) GO TO 28
      IF(KACT(M).EQ.1) GO TO 28
      IWAVE=IWAVE+1
      KACT(M)=1
   28 CONTINUE
      GO TO 35
   30 CONTINUE
      MINDEG=MIN(MINDEG,J-1)
   35 CONTINUE
      IB1=IB+1
!     IB1=ROW BANDWIDTH FOR ROW I (DIAGONAL INCLUDED)
      MAXB=MAX(MAXB,IB1)
      IF(KACT(I).EQ.1) IWAVE=IWAVE-1
!   IWAVE=CURRENT NUMBER OF ACTIVE COLUMNS FOR ROW I (DIAGONAL INCLUDED)
      MAXW=MAX(MAXW,IWAVE)
      SUMW=SUMW+IWAVE
      SUMSQ=SUMSQ+FLOAT(IWAVE)*FLOAT(IWAVE)
      BSUMSQ=BSUMSQ+FLOAT(IB1)*FLOAT(IB1)
   40 CONTINUE
!
      ANN=FLOAT(KT)
      AVERW=FLOAT(SUMW)/ANN
      RMS=SQRT(SNGL(SUMSQ)/ANN)
      BRMS=SQRT(SNGL(BSUMSQ)/ANN)
      RETURN
      END SUBROUTINE WAVEY

! ##################################################################################################################################
      SUBROUTINE ZERO(LIST,NL)
!
!     DELETE ZEROS FROM A LIST OF INTEGERS ORIGINALLY OF LENGTH NL.
!     A CORRECTED LENGTH NL IS RETURNED.
!
      INTEGER LIST(*)
! B////////////////////////////////////////////////////////////////////B
! Add the following so we can use IMPLICIT NONE

      INTEGER  I      ,KT     ,NL

! E////////////////////////////////////////////////////////////////////E
      IF(NL.LE.0) RETURN
      KT=0
      DO 10 I=1,NL
      IF(LIST(I).EQ.0) GO TO 10
      KT=KT+1
      LIST(KT)=LIST(I)
   10 CONTINUE
      NL=KT
      RETURN
      END SUBROUTINE ZERO

! ##################################################################################################################################

      END MODULE BANDIT_MODULE
