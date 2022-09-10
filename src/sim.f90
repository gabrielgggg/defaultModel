MODULE sim
  USE iso_Fortran_env, ONLY: wp => real64
  IMPLICIT NONE

  REAL(wp) :: piwp = 4.0 * ATAN(1.0)
CONTAINS

  SUBROUTINE fixSeed()
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: seed
    INTEGER :: n
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    seed = 1989
    CALL RANDOM_SEED(put = seed)
  END SUBROUTINE fixSeed

  FUNCTION simDiscrete(pmf) RESULT(val)
    IMPLICIT NONE
    REAL(wp), DIMENSION(:), INTENT(IN) :: pmf
    INTEGER :: val
    INTEGER :: sz, ix
    REAL(wp) :: draw, theSum

    sz = SIZE(pmf)
    CALL simUniform(draw)
    theSum = 0.0_wp
    val = -1
    DO ix = 1,sz
      theSum = theSum + pmf(ix)
      IF (draw <= theSum) THEN
        val = ix
        EXIT
      END IF
    END DO
  END FUNCTION simDiscrete

  SUBROUTINE simMarkov(piMat, newIx, oldIx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: oldIx
    REAL(wp), DIMENSION(:, :), INTENT(IN) :: piMat
    INTEGER, INTENT(OUT) :: newIx

    newIx = simDiscrete( piMat(oldIx, :) )
  END SUBROUTINE simMarkov

  SUBROUTINE simStdNormal(val)
    IMPLICIT NONE
    REAL(wp), INTENT(OUT) :: val
    REAL(wp) :: v1, v2

    CALL RANDOM_NUMBER(v1)
    CALL RANDOM_NUMBER(v2)
    val = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * piwp * v2 )
  END SUBROUTINE simStdNormal

  SUBROUTINE simUniform(val)
    IMPLICIT NONE
    REAL(wp), INTENT(OUT) :: val

    CALL RANDOM_NUMBER(val)
  END SUBROUTINE simUniform

  SUBROUTINE simDiscreteUniform(sz, retIx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sz
    INTEGER, INTENT(OUT) :: retIx
    REAL(wp) :: rno

    CALL RANDOM_NUMBER(rno)
    retIx = FLOOR(sz * rno) + 1
  END SUBROUTINE simDiscreteUniform

END MODULE sim

