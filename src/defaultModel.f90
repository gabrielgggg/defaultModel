!
! Sovereign default model with long-term debt
!
PROGRAM defaultModel
  USE iso_Fortran_env, ONLY: wp => real64
  USE defMod
  IMPLICIT NONE

  INTEGER :: iter
  REAL(wp) :: errV, errQ
  REAL(wp), DIMENSION(bSz) :: cc
  INTEGER :: yIx, bIx, yPrIx, bPrIx

  REAL(wp), DIMENSION(ySz) :: hy, uhy

  REAL(wp) :: Wbar, theSum
  REAL(wp), DIMENSION(bSz) :: theExps

  WRITE (*, *) "Start!"
  CALL allocateAll()
  CALL prepareShocksAndGrids()

  hy = hFun(yGrid)
  uhy = uFun(hy)

!  DO yIx = 1,ySz
!    WRITE (*, "(I5,2(A,F10.5))") yIx, CHAR(9), yGrid(yIx), CHAR(9), hy(yIx)
!  END DO

  ! Initial values
  Vd0 = uhy
  q0 = 1.0_wp

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx,bIx)
  DO yIx = 1,ySz
  DO bIx = 1,bSz
    V0(yIx, bIx) = uFun(MAX(yGrid(yIx) - kappa * bGrid(bIx), 0.01_wp))
  END DO
  END DO

  ! Main lopp
  iter = 1
  errV = 1.0_wp
  errQ = 1.0_wp
  DO WHILE ( iter <= maxIter .AND. ( errV > tolErrV .OR. errQ > tolErrQ ))
    ! V0, Vd0 => Vd1
    !
    !$OMP PARALLEL DO PRIVATE(yIx)
    DO yIx = 1,ySz
      Vd1(yIx) = uhy(yIx) + beta * DOT_PRODUCT( yPi(yIx, :), &
        gamm * V0(:, 1) + (1.0_wp - gamm) * Vd0 )
    END DO

    ! V0, q0 => W, Vr, bPol
    ! Vd1, Vr => V1, dPol
    !
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx,bIx,cc,Wbar,theExps,theSum)
    DO yIx = 1,ySz
    DO bIx = 1,bSz
      cc = yGrid(yIx) - kappa * bGrid(bIx) &
        + q0(yIx, :) * ( bGrid - (1.0_wp - delta) * bGrid(bIx) )

      WHERE (cc <= 0.0_wp)
        W(yIx, bIx, :) = veryNegative
      ELSEWHERE
        W(yIx, bIx, :) = uFun(cc) + beta * MATMUL( yPi(yIx, :), V0 )
      END WHERE

      Wbar = MAXVAL(W(yIx, bIx, :))
      theExps = EXP( (W(yIx, bIx, :) - Wbar) / rhoB )
      theSum = SUM(theExps)

      Vr(yIx, bIx) = Wbar + rhoB * LOG( theSum )
      bPol(yIx, bIx, :) = theExps / theSum

      Wbar = MAX( Vr(yIx, bIx), Vd1(yIx) )
      theSum = EXP( (Vd1(yIx) - Wbar) / rhoD ) + EXP( (Vr(yIx, bIx) - Wbar) / rhoD ) 
      dPol(yIx, bIx) = EXP( (Vd1(yIx) - Wbar) / rhoD ) / theSum
      V1(yIx, bIx) = Wbar + rhoD * LOG( theSum )
    END DO
    END DO

    ! q0, dPol, bPol => q1
    !
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx,bPrIx,yPrIx)
    DO yIx = 1,ySz
    DO bPrIx = 1,bSz
      q1(yIx, bPrIx) = 0.0_wp
      DO yPrIx = 1,ySz
        q1(yIx, bPrIx) = q1(yIx, bPrIx) &
          + yPi(yIx, yPrIx) * (1.0_wp - dPol(yPrIx, bPrIx)) &
          * ( kappa + (1.0_wp - delta) * &
              DOT_PRODUCT( bPol(yPrIx, bPrIx, :), q0(yPrIx, :) ) )
      END DO
    END DO
    END DO
    q1 = q1 / (1.0_wp + rf)

    ! Check for convergence and iterate
    errV = MAX( &
      MAXVAL(ABS( V1 - V0  )), &
      MAXVAL(ABS( Vd1 - Vd0 )) )
    errQ = MAXVAL(ABS( q1 - q0 ))
    IF (MOD(iter, 10) == 0) WRITE (*, "(I5,2ES20.5)") iter, errV, errQ
    iter = iter + 1
    V0 = V1
    Vd0 = Vd1
    q0 = q1
  END DO ! end of main loop WHILE


  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx,bIx)
  DO yIx = 1,ySz
  DO bIx = 1,bSz
    EbPr(yIx, bIx) = DOT_PRODUCT( bGrid, bPol(yIx, bIx, :) )
  END DO
  END DO

  WRITE (*, *) "Saving to disk..."
  CALL saveResults()

  WRITE (*, *) "Simulate..."
  CALL simulate()

  WRITE (*, *) "The end! ^_^"

CONTAINS

END PROGRAM defaultModel
