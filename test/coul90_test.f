program  coul90_test
    !!---- COULOMB DRIVER    (FULL TEST VERSION)        FILE = COUL90_T.FOR
    !!---- A.R.BARNETT       MANCHESTER                 EASTER 1995 (REVISED)
    !!---- COMPUTES COULOMBS, SPHERICAL BESSELS, AND CYLINDRICAL BESSELS
    !!---- PREPARED FOR 'ATOMIC PHYSICS' BOOK 1995, EDITOR K. BARTSCHAT
    !!----          THE DRIVER TESTS COUL90, SBESJY AND RICBES SUBROUTINES

    use iso_fortran_env,      only: dp => real64, stdout => output_unit
    use libcoul90,            only: coul90, ricbes, sbesjy
    use libcoul90__stede,     only: ERR, NF
    use libcoul90__steed,     only: PACCQ, NFP, NPQ, IEXP, MINL
    use libcoul90__deset,     only: CF1, P, Q, F, GAMM, WRONSK
    use libcoul90__constants, only: ZERO

    implicit none

    integer, parameter :: M = 1200

    real(dp) ::  FC(0:M),  GC(0:M),  FCP(0:M),  GCP(0:M)
    real(dp) ::  XJ(0:M),   Y(0:M),  XJP(0:M),   YP(0:M)
    real(dp) ::  PSI(0:M), CHI(0:M), PSIP(0:M), CHIP(0:M)
    real(dp) ::  X(21), ETA, XM

    logical :: SPHRIC, FIRST

    integer :: LAMBDA(10), LRANGE, NL, NX, KFN
    integer :: funit_in
    integer :: funit_out
    integer :: LMIN, IFAIL, K, N, L

    character(4) :: CHFN
    character(7) :: INP
    character(72) ::  TEXT

    character(16) :: COULTEST_IN = "test/COULTEST.IN"
    character(17) :: COULTEST_OUT = "test/COULTEST.OUT"

    ! common /STEDE/   ERR, NF
    ! common /STEED/   PACCQ, NFP, NPQ, IEXP, MINL
    ! common /DESET/   CF1, P, Q, F, GAMM, WRONSK
    ! ZERO = 0.D0
    OPEN  (newunit = funit_in, FILE = COULTEST_IN)
    OPEN  (newunit = funit_out, FILE = COULTEST_OUT)
    WRITE (funit_out, 1000)
    READ  (funit_in, 1001, END = 50) TEXT

    10   CONTINUE
    NL=0
    NX=0
    READ  (funit_in, 1002, ERR=41, END=50) CHFN, INP, NL, NX, ETA, XM
    IF ( NL.LE.0 )  STOP
    READ  (funit_in, *, ERR=42) (LAMBDA(K), K = 1, NL)  ! number of L-values
    READ  (funit_in, *, ERR=43) (X(K),  K = 1, ABS(NX)) ! number of X-values
    IF ( NX.LT.0 ) THEN
       READ  (funit_in, 1001) TEXT
       WRITE (funit_out, 1003) TEXT                ! annotate output
    END IF
    ! DO 40  K = 1, ABS(NX)                  ! main loop over x-values
    DO K = 1, ABS(NX)                  ! main loop over x-values
       LRANGE = LAMBDA(NL) - LAMBDA(1)     ! additional L-values
       FIRST  = ( K.EQ.1 )
       IF ( CHFN.EQ.'COUL' ) KFN = 0
       IF ( CHFN.EQ.'SBES' ) KFN = 1
       IF ( CHFN.EQ.'BESS' ) KFN = 2
       IF ( CHFN.EQ.'BESI' .AND. FIRST ) THEN
          KFN = 2                          ! Special case
          XM  = 1.D0 / XM                  ! eg for 1/3 Bessels
       END IF
       ! KFN = 0 (COULOMBS), = 1 (SPHERICAL BESSELS),  = 2 (CYLINDER BESSELS)
       CALL COUL90( X(K), ETA, XM, LRANGE, FC, GC, FCP, GCP, KFN, IFAIL)
       IF ( IFAIL.NE.0 )  WRITE (funit_out,1200) '  COUL90 ERROR! IFAIL=',IFAIL
       SPHRIC = ( KFN.NE.2 .AND. ETA.EQ.ZERO )
       IF ( SPHRIC ) THEN
           CALL SBESJY(X(K), LRANGE,  XJ,   Y,  XJP,   YP, IFAIL)
           IF ( IFAIL.NE.0 )  WRITE (2,1200) '  SBESJY ERROR! IFAIL=',IFAIL
           CALL RICBES(X(K), LRANGE, PSI, CHI, PSIP, CHIP, IFAIL)
           IF ( IFAIL.NE.0 )  WRITE (funit_out,1200) '  RICBES ERROR! IFAIL=',IFAIL
       END IF
       IF ( KFN.EQ.0 ) WRITE (funit_out, 1004) ETA, X(K), XM, KFN,' (Coulomb)'
       IF ( KFN.EQ.1 ) WRITE (funit_out, 1005)      X(K), XM, KFN,' (SphBess)'
       IF ( KFN.EQ.funit_out ) WRITE (funit_out, 1006)      X(K), XM, KFN,' (CylBess)'
       IF ( IFAIL.NE.0 ) THEN
          WRITE (funit_out, 1007) IFAIL, NL, KFN, X(K), ETA, XM, LAMBDA(NL)
          GOTO 10                 ! abort & try next parameters
       END IF
       TEXT = ' The 3 sets of results are COUL90(KFN), SBESJY & (1/X) RICBES'
       IF ( SPHRIC )   WRITE (funit_out, 1003) TEXT
       LMIN = LAMBDA(1) - IDINT(XM)
       WRITE (6, 1201) X(K), ETA, XM, LRANGE, LMIN
       !  compile these lines for diagnostics :
       !                        WRITE (funit_out, 1008) NFP, NPQ, PACCQ
       !        IF ( KFN.GT.0 ) WRITE (funit_out, 1009) NF , KFN, ERR
       !                        WRITE (funit_out, 1202) CF1, P, Q, F, GAMM, WRONSK
       ! DO 30  N = 1, NL
       DO N = 1, NL
           L = LMIN + LAMBDA(N)
           WRITE (funit_out, 1010) LAMBDA(N),  FC(L),  GC(L),  FCP(L),  GCP(L)
           IF ( SPHRIC ) THEN
               WRITE (funit_out, 1010) LAMBDA(N),  XJ(L),   Y(L),  XJP(L),   YP(L)
               WRITE (funit_out, 1010) LAMBDA(N), PSI(L), CHI(L), PSIP(L), CHIP(L)
           END IF
       END DO
       ! 30      CONTINUE
       IF ( IEXP.GT.1 ) WRITE (funit_out, 1011) IEXP
    END DO
    ! 40   CONTINUE
    GOTO 10                   ! more data
    41   WRITE(stdout,1002)   CHFN, INP, NL, NX, ETA, XM
    42   WRITE(stdout,1012)   NL, NX, (REAL(LAMBDA(K)), K = 1, NL)
    43   WRITE(stdout,1012)   NL, NX, (X(K), K = 1, ABS(NX))

    50   CLOSE (funit_in)
    CLOSE (funit_out)
    !----
    1000 FORMAT (10X, ' TEST OF THE CONTINUED-FRACTION COULOMB & BESSEL', ' PROGRAM - COUL90'//12X, &
      ' WHEN LAMBDA IS REAL (IN GENERAL) = L + XM '// '  F IS REGULAR AT THE ORIGIN ( X = 0 ) WHILE',/ &
      '  G IS IRREGULAR ( => INFINITY AT X = 0 )'//5X, 'L', 2X, ' F(ETA,X,LAMBDA)', 2X, ' G(ETA,X,LAMBDA)', 1X, &
      ' FCP = D/DX (F)', 2X, ' GCP = D/DX (G)'/)
    1001 FORMAT (A72)
    1002 FORMAT (1X, A4, A7, 2I3, 2F10.2 )
    1003 FORMAT (/4X, A72)
    1004 FORMAT (/1X, 'ETA =', F9.3, 4X, ' X = ', F10.3, 4X, 'XLMIN = ', F6.2, 4X, 'KFN = ', I2, A10)
    1005 FORMAT (/5X, 'sph Bessels', 3X, ' X = ', F10.3, 4X, 'XLMIN = ', F6.2, 4X, 'KFN = ', I2, A10)
    1006 FORMAT (/5X, 'CYL Bessels', 3X, ' X = ', F10.3, 4X, 'XLMIN = ', F6.2, 4X, 'KFN = ', I2, A10)
    1007 FORMAT (1X, 'IFAIL =', I4, ' NL,KFN =', 2I4, 3F12.4, I12)
    ! 1008 FORMAT (1X,': NFP,NPQ,PACCQ = ', 2I5, 1PD9.0, '---')
    ! 1009 FORMAT (1X,': NF ,KFN,ERR   = ', 2I5, 1PD9.0, '---')
    1010 FORMAT (1X, I5, 1P4D17.7 )
    1011 FORMAT (12X, ' **** IEXP = ', I6, '  F,FP *10**(-IEXP)', '  G,GP *10**(+IEXP)')
    1012 FORMAT (1X, 2I3, F10.3, 9F6.2)
    1200 FORMAT (1X, A25, I5)
    1201 FORMAT (3F10.3, 3I10, F10.3)
    ! 1202 FORMAT (6(1PD13.5))

end program coul90_test
