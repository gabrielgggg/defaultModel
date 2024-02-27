MODULE defMod
  USE iso_Fortran_env, ONLY: wp => real64
  USE NL
  USE sim
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: outDir = "./results/"
  CHARACTER(LEN=1), PARAMETER :: TAB = CHAR(9)

  REAL(wp), PARAMETER :: crra = 2.0_wp
  REAL(wp), PARAMETER :: lbd0 = -0.48_wp, lbd1 = 0.525_wp
  REAL(wp), PARAMETER :: rhoY = 0.95_wp
  REAL(wp), PARAMETER :: sigmaY = 0.005_wp
  REAL(wp), PARAMETER :: rf = 1.04_wp**0.25_wp - 1.0_wp
  REAL(wp), PARAMETER :: beta = 0.9775_wp
  REAL(wp), PARAMETER :: gamm = 1.0_wp / (2.0_wp * 4.0_wp)

  REAL(wp), PARAMETER :: macaulay = 5.0_wp * 4.0_wp ! 20 quarters Macaulay duration
  REAL(wp), PARAMETER :: delta = (1 + rf) / macaulay - rf
  REAL(wp), PARAMETER :: kappa = delta + rf

  REAL(wp), PARAMETER :: rhoD = 5.0D-4
  REAL(wp), PARAMETER :: rhoB = 1.0D-5

  INTEGER, PARAMETER :: maxIter = 1000, simSz = 100000
  INTEGER, PARAMETER :: ySz = 31
  INTEGER, PARAMETER :: bSz = 600
  REAL(wp), PARAMETER :: tolErrV = 1.0D-6, tolErrQ = 1.0D-6
  REAL(wp), PARAMETER :: bMin = 0.0_wp, bMax = 0.75_wp
  REAL(wp), PARAMETER :: veryNegative = -1.0D+6

  REAL(wp), DIMENSION(ySz) :: yGrid(ySz), yStat(ySz)
  REAL(wp), DIMENSION(ySz, ySz) :: yPi(ySz, ySz)
  REAL(wp), DIMENSION(bSz) :: bGrid(bSz)

  REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: V0, V1
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: Vd0, Vd1
  REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: Vr, dPol, q0, q1, EbPr
  REAL(wp), ALLOCATABLE, DIMENSION(:, :, :) :: bPol, W

CONTAINS

  SUBROUTINE allocateAll()
    ALLOCATE( V0(ySz, bSz), V1(ySz, bSz) )
    ALLOCATE( Vd0(ySz), Vd1(ySz) )
    ALLOCATE( Vr(ySz, bSz), dPol(ySz, bSz), q0(ySz, bSz), q1(ySz, bSz) )
    ALLOCATE( bPol(ySz, bSz, bSz), W(ySz, bSz, bSz) )
    ALLOCATE( EbPr(ySz, bSz) )
  END SUBROUTINE allocateAll

  SUBROUTINE prepareShocksAndGrids()
    CALL linspace(bGrid, bMin, bMax, bSz)

    CALL discretizeAR1(0.0_wp, rhoY, sigmaY, ySz, 3.0_wp, yGrid, yPi, yStat)
    yGrid = EXP( yGrid - 0.5_wp * sigmaY**2 / (1.0_wp - rhoY**2) )
  END SUBROUTINE prepareShocksAndGrids

  ELEMENTAL FUNCTION uFun(cons)
    REAL(wp), INTENT(IN) :: cons
    REAL(wp) :: uFun

    IF (cons <= 0.0_wp) THEN
      ERROR STOP "Negative consumption passed to u(.)"
    ELSE
      uFun = ( cons**(1.0_wp - crra) - 1.0_wp ) / (1.0_wp - crra)
    END IF
  END FUNCTION uFun

  ELEMENTAL FUNCTION hFun(yVal)
    REAL(wp), INTENT(IN) :: yVal
    REAL(wp) :: hFun

    hFun = yVal - MAX(0.0_wp, lbd0 * yVal + lbd1 * yVal**2)
  END FUNCTION hFun

  SUBROUTINE simulate()
    INTEGER, DIMENSION(:), ALLOCATABLE :: ySimIx, bSimIx, bPrSimIx, dSimIx
    REAL(wp), DIMENSION(:), ALLOCATABLE :: spSim, cSim, tbSim, gdpSim
    INTEGER :: tIx, iunit
    REAL(wp) :: unifDraw

    CALL fixSeed()
    
    ALLOCATE( ySimIx(simSz), bSimIx(simSz), bPrSimIx(simSz), dSimIx(simSz) )
    ALLOCATE( spSim(simSz), cSIm(simSz), tbSim(simSz), gdpSim(simSz) )

    ySimIx(1) = CEILING(REAL(ySz, wp) / 2.0_wp)
    bSimIx(1) = 1
    bPrSimIx(1) = 1
    dSimIx(1) = 0 ! Not in default initially

    DO tIx = 2,simSz
      IF (dSimIx(tIx-1) == 1) THEN ! last period in default/exclusion
        CALL simUniform(unifDraw)
        IF (unifDraw <= gamm) THEN ! Come back to market with zero debt
          dSimIx(tIx) = 0
          bSimIx(tIx) = 1 ! assuming first bGrid element corresponds to zero
        ELSE
          dSimIx(tIx) = 1
          bSimIx(tIx) = bSimIx(tIx-1)
        END IF
      ELSE
        ! If not in default, yesterday's b' is today's b
        bSimIx(tIx) = bPrSimIx(tIx-1)
        dSimIx(tIx) = 0
      END IF

      CALL simMarkov(yPi, ySimIx(tIx), ySimIx(tIx-1))

      IF (dSimIx(tIx) == 0) THEN ! Default today?
        CALL simUniform(unifDraw)
        IF (unifDraw <= dPol(ySimIx(tIx), bSimIx(tIx))) THEN
          dSimIx(tIx) = 1
        END IF
      END IF

      IF (dSimIx(tIx) == 0) THEN ! not in default today
        bPrSimIx(tIx) = simDiscrete( bPol(ySimIx(tIx), bSimIx(tIx), :) )
        spSim(tIx) = kappa * ( 1.0_wp / q1(ySimIx(tIx), bPrSimIx(tIx)) - 1.0_wp )
        gdpSim(tIx) = yGrid(ySimIx(tIx))
        cSim(tIx) = gdpSim(tIx) - kappa * bGrid(bSimIx(tIx)) &
          + q1(ySimIx(tIx), bPrSimIx(tIx)) * (bGrid(bPrSimIx(tIx)) - (1.0_wp - delta) * bGrid(bSimIx(tIx)))
        tbSim(tIx) = gdpSim(tIx) - cSim(tIx)
      ELSE ! in default today
        bPrSimIx(tIx) = bSimIx(tIx)
        spSim(tIx) = veryNegative
        gdpSim(tIx) = hFun(yGrid(ySimIx(tIx)))
        cSim(tIx) = gdpSim(tIx)
        tbSim(tIx) = 0.0_wp
      END IF

      IF (MOD(tIx, 10000) == 0) THEN
        WRITE (*, *) "Simulated ", tIx, " out of ", simSz
      END IF
    END DO

    OPEN(newunit=iunit, file=outDir // "sim.tab")
    DO tIx = 300,simSz
      WRITE (iunit, "(4(I5,A1),4(ES20.5,A1))") ySimIx(tIx), TAB, bSimIx(tIx), TAB, & 
        bPrSimIx(tIx), TAB, dSimIx(tIx), TAB, spSim(tIx), TAB, cSim(tIx), TAB, &
        gdpSim(tIx), TAB, tbSim(tIx), TAB
    END DO
    CLOSE(iunit)
  END SUBROUTINE simulate

  SUBROUTINE saveResults()
    INTEGER :: iunit

    OPEN(newunit=iunit, file=outDir // "parameters.tab")
    WRITE (iunit, "(I20)") ySz
    WRITE (iunit, "(I20)") bSz
    WRITE (iunit, "(E20.10)") crra
    WRITE (iunit, "(E20.10)") rf
    WRITE (iunit, "(E20.10)") delta
    ! TODO fill out later!
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "yGrid.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) yGrid
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "yPi.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) yPi
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "bGrid.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) bGrid
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "V.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) V1
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "Vr.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) Vr
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "Vd.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) Vd1
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "q.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) q1
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "dPol.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) dPol
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "bPol.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) bPol
    CLOSE(iunit)

    OPEN(newunit=iunit, file=outDir // "EbPr.bin", &
      FORM="unformatted", ACCESS="stream", STATUS="unknown")
    WRITE (iunit) EbPr
    CLOSE(iunit)
  END SUBROUTINE saveResults
 
!   SUBROUTINE loadEquilibrium()
!     INTEGER :: iunit
! 
!     OPEN(newunit=iunit, file=outDir // "q.bin", &
!       FORM="unformatted", ACCESS="stream", STATUS="old")
!     READ (iunit) q1
!     CLOSE(iunit)
! 
!   END SUBROUTINE loadEquilibrium

END MODULE defMod
