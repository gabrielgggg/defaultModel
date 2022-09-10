MODULE NL
  USE, INTRINSIC :: iso_Fortran_env, ONLY: wp => real64
  IMPLICIT NONE

  LOGICAL, PARAMETER :: extrapLin = .TRUE.

CONTAINS

  PURE FUNCTION inv_normal_cdf(p)
    implicit none
  
    REAL(wp), PARAMETER, DIMENSION(8) :: a = (/ &
      3.3871328727963666080D+00, &
      1.3314166789178437745D+02, &
      1.9715909503065514427D+03, &
      1.3731693765509461125D+04, &
      4.5921953931549871457D+04, &
      6.7265770927008700853D+04, &
      3.3430575583588128105D+04, &
      2.5090809287301226727D+03 /)
    REAL(wp), PARAMETER, DIMENSION(8) :: b = (/ &
      1.0D+00, &
      4.2313330701600911252D+01, &
      6.8718700749205790830D+02, &
      5.3941960214247511077D+03, &
      2.1213794301586595867D+04, &
      3.9307895800092710610D+04, &
      2.8729085735721942674D+04, &
      5.2264952788528545610D+03 /)
    real   ( kind = 8 ), PARAMETER, DIMENSION(8) :: c = (/ &
      1.42343711074968357734D+00, &
      4.63033784615654529590D+00, &
      5.76949722146069140550D+00, &
      3.64784832476320460504D+00, &
      1.27045825245236838258D+00, &
      2.41780725177450611770D-01, &
      2.27238449892691845833D-02, &
      7.74545014278341407640D-04 /)
    REAL(wp), PARAMETER :: const1 = 0.180625D+00
    REAL(wp), PARAMETER :: const2 = 1.6D+00
    REAL(wp), PARAMETER, DIMENSION(8) :: d = (/ &
      1.0D+00, &
      2.05319162663775882187D+00, &
      1.67638483018380384940D+00, &
      6.89767334985100004550D-01, &
      1.48103976427480074590D-01, &
      1.51986665636164571966D-02, &
      5.47593808499534494600D-04, &
      1.05075007164441684324D-09 /)
    REAL(wp), PARAMETER, DIMENSION(8) :: e = (/ &
      6.65790464350110377720D+00, &
      5.46378491116411436990D+00, &
      1.78482653991729133580D+00, &
      2.96560571828504891230D-01, &
      2.65321895265761230930D-02, &
      1.24266094738807843860D-03, &
      2.71155556874348757815D-05, &
      2.01033439929228813265D-07 /)
    REAL(wp), PARAMETER, DIMENSION(8) :: f = (/ &
      1.0D+00, &
      5.99832206555887937690D-01, &
      1.36929880922735805310D-01, &
      1.48753612908506148525D-02, &
      7.86869131145613259100D-04, &
      1.84631831751005468180D-05, &
      1.42151175831644588870D-07, &
      2.04426310338993978564D-15 /)
    REAL(wp), INTENT(IN) :: p
    REAL(wp) :: inv_normal_cdf
    REAL(wp) :: q
    REAL(wp) :: r
    REAL(wp), PARAMETER :: split1 = 0.425_wp
    REAL(wp), PARAMETER :: split2 = 5.0_wp
  
    if ( p <= 0.0_wp ) then
      inv_normal_cdf = -HUGE(p)
      RETURN
    end if
  
    if ( 1.0_wp <= p ) then
      inv_normal_cdf = HUGE(p)
      RETURN
    end if
  
    q = p - 0.5_wp
  
    if (abs( q ) <= split1) then
      r = const1 - q * q
      inv_normal_cdf = q * poly_value(8, a, r) / poly_value(8, b, r)
    else
  
      if ( q < 0.0_wp ) then
        r = p
      else
        r = 1.0_wp - p
      end if
  
      if ( r <= 0.0_wp ) then
        inv_normal_cdf = -1.0_wp
      end if
  
      r = sqrt(-log(r))
  
      if (r <= split2) then
        r = r - const2
        inv_normal_cdf = poly_value(8, c, r) / poly_value(8, d, r)
      else
        r = r - split2
        inv_normal_cdf = poly_value (8, e, r) / poly_value(8, f, r)
      end if
  
      if ( q < 0.0D+00 ) then
        inv_normal_cdf = -inv_normal_cdf
      end if
    end if
  END FUNCTION inv_normal_cdf

  PURE FUNCTION poly_value(n, a, x)
    implicit none
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: a(n)
    INTEGER :: i
    REAL(wp) :: poly_value
    REAL(wp), INTENT(IN) :: x
  
    poly_value = 0.0D+00
    do i = n, 1, -1
      poly_value = poly_value * x + a(i)
    end do
  END FUNCTION poly_value

  SUBROUTINE moments(ddata, vvalid, mmean, sstd)
    IMPLICIT NONE
    REAL(wp), INTENT(IN), DIMENSION(:) :: ddata
    LOGICAL, INTENT(IN), DIMENSION(SIZE(ddata)) :: vvalid
    REAL(wp), INTENT(OUT) :: mmean, sstd
    INTEGER :: sz, validSz, ix
    REAL(wp) :: tmpSum, vvar

    sz = SIZE(ddata)

    validSz = 0
    tmpSum = 0.0_wp
    DO ix = 1,sz
      IF (vvalid(ix)) THEN
        validSz = validSz + 1
        tmpSum = tmpSum + ddata(ix)
      END IF
    END DO

    mmean = tmpSum / validSz

    tmpSum = 0.0_wp
    DO ix = 1,sz
      IF (vvalid(ix)) THEN
        tmpSum = tmpSum + ddata(ix) * ddata(ix)
      END IF
    END DO
    vvar = tmpSum / (validSz - 1.0_wp) - mmean * mmean

    sstd = SQRT(vvar)
  END SUBROUTINE moments

  !
  !
  !
  FUNCTION corr(ddata1, ddata2, vvalid) RESULT(ccorr)
    IMPLICIT NONE
    REAL(wp), INTENT(IN), DIMENSION(:) :: ddata1, ddata2
    LOGICAL, INTENT(IN), DIMENSION(SIZE(ddata1)) :: vvalid
    REAL(wp) :: ccorr
    INTEGER :: sz, validSz, ix
    REAL(wp) :: mmean1, mmean2, sstd1, sstd2, tmpSum

    sz = SIZE(ddata1)
    CALL moments(ddata1, vvalid, mmean1, sstd1)
    CALL moments(ddata2, vvalid, mmean2, sstd2)

    ! WRITE (*, *) mmean1, mmean2, sstd1, sstd2

    IF (sstd1 < 1.0D-10 .OR. sstd2 < 1.0D-10) THEN
      ccorr = 0.0_wp
    ELSE
      validSz = 0
      tmpSum = 0.0_wp
      DO ix = 1,sz
        IF (vvalid(ix)) THEN
          validSz = validSz + 1
          tmpSum = tmpSum + ddata1(ix) * ddata2(ix)
        END IF
      END DO
      ccorr = tmpSum / validSz - mmean1 * mmean2
      ccorr = ccorr / sstd1 / sstd2
    END IF
  END FUNCTION corr


  SUBROUTINE discreteNormalPrs(mmean, sstd, ggrid, pprs)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: mmean, sstd
    REAL(wp), INTENT(IN), DIMENSION(:) :: ggrid
    REAL(wp), INTENT(OUT), DIMENSION(SIZE(ggrid)) :: pprs
    INTEGER :: ix, nn
    REAL(wp) :: tmp1, tmp2, tmp3, tmp4

    nn = SIZE(ggrid)
    DO ix = 1,nn
      IF (ix == 1) THEN
        CALL cumnor( ((ggrid(2) + ggrid(1)) * 0.5_wp - mmean) / sstd, tmp1, tmp2)
        pprs(1) = tmp1
      ELSEIF (ix == nn) THEN
        CALL cumnor( ((ggrid(nn) + ggrid(nn-1)) * 0.5_wp - mmean) / sstd, tmp1, tmp2)
        pprs(nn) = tmp2
      ELSE
        CALL cumnor( ((ggrid(ix) + ggrid(ix-1)) * 0.5_wp - mmean) / sstd, tmp1, tmp2)
        CALL cumnor( ((ggrid(ix+1) + ggrid(ix)) * 0.5_wp - mmean) / sstd, tmp3, tmp4)
        pprs(ix) = tmp3-tmp1
      END IF
    END DO
  END SUBROUTINE discreteNormalPrs

  !
  !
  !
  SUBROUTINE discretizeNormal(mmean, sstd, nn, nsds, ggrid, pprs)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: mmean, sstd
    INTEGER, INTENT(IN) :: nn, nsds
    REAL(wp), INTENT(OUT), DIMENSION(nn) :: ggrid, pprs

    IF (nn == 1) THEN
      ggrid(1) = mmean
      pprs(1) = 1.0_wp
      RETURN
    END IF

    CALL linspace(ggrid, mmean-nsds*sstd, mmean+nsds*sstd, nn)
    CALL discreteNormalPrs(mmean, sstd, ggrid, pprs)
  END SUBROUTINE discretizeNormal

  !
  !
  SUBROUTINE discretizeAR1(mmean, rrho, sstd, nn, nsds, ggrid, ttran, stationary)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: mmean, rrho, sstd, nsds
    INTEGER, INTENT(IN) :: nn
    REAL(wp), INTENT(OUT), DIMENSION(nn) :: ggrid, stationary
    REAL(wp), INTENT(OUT), DIMENSION(nn, nn) :: ttran
    INTEGER :: ix
    REAL(wp) :: tmp1, condMean
    REAL(wp), DIMENSION(nn) :: stt0

    IF (nn == 1) THEN
      ggrid(1) = mmean
      ttran(1, 1) = 1.0_wp
      stationary(1) = mmean
      RETURN
    END IF

    CALL linspace(ggrid, &
      mmean-nsds*sstd/SQRT(1.0_wp - rrho*rrho), &
      mmean+nsds*sstd/SQRT(1.0_wp - rrho*rrho), nn)

    DO ix = 1,nn
      condMean = (1.0_wp - rrho) * mmean + rrho * ggrid(ix)
      CALL discreteNormalPrs(condMean, sstd, ggrid, ttran(ix, :))
    END DO

     stationary = 0.0_wp
     stationary(MAX(1, nn/2)) = 1.0_wp
     stt0 = 0.0_wp
     tmp1 = 1.0_wp
     DO WHILE (tmp1 > 1.0D-12)
       stt0 = MATMUL( stationary, ttran )
       tmp1 = MAXVAL(ABS( stt0 - stationary ))
       stationary = stt0
     END DO
  END SUBROUTINE discretizeAR1

  !
  !
  !
  SUBROUTINE degenerateMarkov(mmean, plusminus, nn, ggrid, ttran)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: mmean, plusminus
    INTEGER, INTENT(IN) :: nn
    REAL(wp), INTENT(OUT), DIMENSION(nn) :: ggrid
    REAL(wp), INTENT(OUT), DIMENSION(nn, nn) :: ttran

    CALL linspace(ggrid, mmean - plusminus, mmean + plusminus, nn)
    ttran = 0.0_wp
    ttran(:, nn/2+1) = 1.0_wp
  END SUBROUTINE degenerateMarkov

  !
  !
  !
  SUBROUTINE mesh2Markovs(grid1, tran1, grid2, tran2, bigGrid1, bigGrid2, bigTran, mmap, revmap)
    IMPLICIT NONE
    REAL(wp), DIMENSION(:), INTENT(IN) :: grid1, grid2
    REAL(wp), INTENT(IN), DIMENSION(SIZE(grid1), SIZE(grid1)) :: tran1
    REAL(wp), INTENT(IN), DIMENSION(SIZE(grid2), SIZE(grid2)) :: tran2
    REAL(wp), INTENT(OUT), DIMENSION(SIZE(grid1) * SIZE(grid2)) :: bigGrid1, bigGrid2
    REAL(wp), INTENT(OUT), &
      DIMENSION(SIZE(grid1) * SIZE(grid2), SIZE(grid1) * SIZE(grid2)) :: bigTran
    INTEGER, INTENT(OUT), DIMENSION(SIZE(grid1) * SIZE(grid2), 3) :: mmap
    INTEGER, INTENT(OUT), DIMENSION(SIZE(grid1), SIZE(grid2)) :: revmap
    INTEGER :: iixx, iiPPxx, ix1, ix2, sz1, sz2

    sz1 = SIZE(grid1)
    sz2 = SIZE(grid2)

    iixx = 0
    DO ix1 = 1,sz1
    DO ix2 = 1,sz2
      iixx = iixx + 1
      bigGrid1(iixx) = grid1(ix1)
      bigGrid2(iixx) = grid2(ix2)
      mmap(iixx, :) = (/ iixx, ix1, ix2 /)
      revmap(ix1, ix2) = iixx
    END DO
    END DO

    DO iixx = 1,(sz1*sz2)
    DO iiPPxx = 1,(sz1*sz2)
      bigTran(iixx, iiPPxx) = tran1(mmap(iixx, 2), mmap(iiPPxx, 2)) &
        * tran2(mmap(iixx, 3), mmap(iiPPxx, 3))
    END DO
    END DO
  END SUBROUTINE mesh2Markovs

  !
  !
  !
  SUBROUTINE mesh3Markovs(grid1, tran1, grid2, tran2, grid3, tran3, &
      bigGrid1, bigGrid2, bigGrid3, bigTran, mmap, revmap)
    IMPLICIT NONE
    REAL(wp), DIMENSION(:), INTENT(IN) :: grid1, grid2, grid3
    REAL(wp), INTENT(IN), DIMENSION(SIZE(grid1), SIZE(grid1)) :: tran1
    REAL(wp), INTENT(IN), DIMENSION(SIZE(grid2), SIZE(grid2)) :: tran2
    REAL(wp), INTENT(IN), DIMENSION(SIZE(grid3), SIZE(grid3)) :: tran3
    REAL(wp), INTENT(OUT), DIMENSION(SIZE(grid1) * SIZE(grid2) * SIZE(grid3)) &
      :: bigGrid1, bigGrid2, bigGrid3
    REAL(wp), INTENT(OUT), &
      DIMENSION(SIZE(grid1) * SIZE(grid2) * SIZE(grid3), &
      SIZE(grid1) * SIZE(grid2) * SIZE(grid3)) :: bigTran
    INTEGER, INTENT(OUT), DIMENSION(SIZE(grid1) * SIZE(grid2) * SIZE(grid3), 4) :: mmap
    INTEGER, INTENT(OUT), DIMENSION(SIZE(grid1), SIZE(grid2), SIZE(grid3)) :: revmap
    INTEGER :: iixx, iiPPxx, ix1, ix2, ix3, sz1, sz2, sz3

    sz1 = SIZE(grid1)
    sz2 = SIZE(grid2)
    sz3 = SIZE(grid3)

    iixx = 0
    DO ix1 = 1,sz1
    DO ix2 = 1,sz2
    DO ix3 = 1,sz3
      iixx = iixx + 1
      bigGrid1(iixx) = grid1(ix1)
      bigGrid2(iixx) = grid2(ix2)
      bigGrid3(iixx) = grid3(ix3)
      mmap(iixx, :) = (/ iixx, ix1, ix2, ix3 /)
      revmap(ix1, ix2, ix3) = iixx
    END DO
    END DO
    END DO

    DO iixx = 1,(sz1*sz2*sz3)
      DO iiPPxx = 1,(sz1*sz2*sz3)
        bigTran(iixx, iiPPxx) = tran1(mmap(iixx, 2), mmap(iiPPxx, 2)) &
          * tran2(mmap(iixx, 3), mmap(iiPPxx, 3)) &
          * tran3(mmap(iixx, 4), mmap(iiPPxx, 4))
      END DO
    END DO
  END SUBROUTINE mesh3Markovs

  !
  !
  !
  PURE SUBROUTINE linspace(ddata, startVal, endVal, noEl)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: noEl
    REAL(wp), DIMENSION(noEl), INTENT(OUT) :: ddata
    REAL(wp), INTENT(IN) :: startVal, endVal
    INTEGER :: ix
    REAL(wp) :: dx

    IF (noEl == 1) THEN
      ddata = 0.5_wp * (endVal + startVal)
    ELSE
      dx = (endVal - startVal) / (noEl - 1)
      ddata = [ ( startVal + (ix - 1) * dx, ix = 1,noEl ) ]
    END IF
  END SUBROUTINE linspace

  !
  !
  !
  PURE FUNCTION findNear(xval, yy) RESULT(yloc)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: xval
    REAL(wp), INTENT(in), DIMENSION(:) :: yy
    INTEGER, DIMENSION(SIZE(yy)) :: mloc
    INTEGER :: yloc

    mloc = MINLOC( (yy - xval)**2 )
    yloc = mloc(1)
  END FUNCTION findNear

  !
  !
  !
  PURE FUNCTION linInterp(xx, yy, xval) RESULT(yval)
    IMPLICIT NONE
    REAL(wp) :: yval
    REAL(wp), DIMENSION(:), INTENT(in) :: xx, yy
    REAL(wp), INTENT(in) :: xval
    INTEGER :: nn
    INTEGER :: low, high, mid

    nn = SIZE(xx)

    IF (nn == 1) THEN
      yval = yy(1)
    ELSE

      IF ( xval >= xx(nn)) THEN
        IF (extrapLin) THEN
          yval = yy(nn-1) + (yy(nn) - yy(nn-1)) / (xx(nn) - xx(nn-1)) * (xval - xx(nn-1))
        ELSE
          yval = yy(nn)
        END IF
      ELSEIF (xval <= xx(1)) THEN
        IF (extrapLin) THEN
          yval = yy(1) + (yy(2) - yy(1)) / (xx(2) - xx(1)) * (xval - xx(1))
        ELSE
          yval = yy(1)
        END IF
      ELSE
        low = 1
        high = nn
        DO WHILE (high - low > 1)
          mid = (low + high) / 2
          IF ( xx(mid) >= xval ) THEN
            ! low = low
            high = mid
          ELSE
            low = mid
            ! high = high
          END IF
        END DO

        yval = yy(low) + (yy(high) - yy(low)) / (xx(high) - xx(low)) * (xval - xx(low))
      END IF

    END IF
  END FUNCTION linInterp

  PURE FUNCTION normal_pdf (xx, mmean, ssd)
    IMPLICIT NONE
    REAL(wp) :: normal_pdf
    REAL(wp), INTENT(IN) :: mmean
    REAL(wp), INTENT(IN) :: ssd
    REAL(wp), PARAMETER :: piwp = 4.0_wp * ATAN(1.0_wp)
    REAL(wp), PARAMETER :: denom = SQRT( 2.0_wp * piwp )
    REAL(wp), INTENT(IN) :: xx
    REAL(wp) :: yy

    yy = ( xx - mmean ) / ssd
    normal_pdf = exp ( -0.5_wp * yy * yy )  / ( ssd * denom )
  END FUNCTION normal_pdf
  
  PURE FUNCTION normal_cdf(xx, mmean, ssd)
        IMPLICIT NONE
    REAL(wp) :: normal_cdf
        REAL(wp), INTENT(IN) :: xx
    REAL(wp), INTENT(IN) :: mmean
    REAL(wp), INTENT(IN) :: ssd
        REAL(wp) :: complement 
        
        CALL cumnor( (xx - mmean) / ssd, normal_cdf, complement)
  END FUNCTION normal_cdf

  PURE subroutine r8_swap ( x, y )
    implicit none
    real (wp), INTENT(INOUT) :: x
    real (wp), INTENT(INOUT) :: y
    real (wp) z

    z = x
    x = y
    y = z
  end

  PURE subroutine cumnor ( arg, cum, ccum )
!*****************************************************************************80
!
!! CUMNOR computes the cumulative normal distribution.
!
!  Discussion:
!
!    This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!    This transportable program uses rational functions that
!    theoretically approximate the normal distribution function to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine dependent
!    constants.
!
!  Author:
!
!    William Cody
!    Mathematics and Computer Science Division
!    Argonne National Laboratory
!    Argonne, IL 60439
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev approximations for the error function,
!    Mathematics of Computation,
!    1969, pages 631-637.
!
!    William Cody,
!    Algorithm 715:
!    SPECFUN - A Portable FORTRAN Package of Special Function Routines
!    and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, 1993, pages 22-32.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the upper limit of integration.
!
!    Output, real ( kind = 8 ) CUM, CCUM, the Normal density CDF and
!    complementary CDF.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) EPS, the argument below which anorm(x)
!    may be represented by 0.5 and above which  x*x  will not underflow.
!    A conservative value is the largest machine number X
!    such that   1.0D+00 + X = 1.0D+00   to machine precision.
!
    implicit none

    real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
      2.2352520354606839287D+00, &
      1.6102823106855587881D+02, &
      1.0676894854603709582D+03, &
      1.8154981253343561249D+04, &
      6.5682337918207449113D-02 /)
    real ( kind = 8 ), INTENT(IN) :: arg
    real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
      4.7202581904688241870D+01, &
      9.7609855173777669322D+02, &
      1.0260932208618978205D+04, &
      4.5507789335026729956D+04 /)
    real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
      3.9894151208813466764D-01, &
      8.8831497943883759412D+00, &
      9.3506656132177855979D+01, &
      5.9727027639480026226D+02, &
      2.4945375852903726711D+03, &
      6.8481904505362823326D+03, &
      1.1602651437647350124D+04, &
      9.8427148383839780218D+03, &
      1.0765576773720192317D-08 /)
    real ( kind = 8 ), INTENT(OUT) :: ccum
    real ( kind = 8 ), INTENT(OUT) :: cum
    real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
      2.2266688044328115691D+01, &
      2.3538790178262499861D+02, &
      1.5193775994075548050D+03, &
      6.4855582982667607550D+03, &
      1.8615571640885098091D+04, &
      3.4900952721145977266D+04, &
      3.8912003286093271411D+04, &
      1.9685429676859990727D+04 /)
    real ( kind = 8 ) del
    real ( kind = 8 ) eps
    integer ( kind = 4 ) i
    real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
      2.1589853405795699D-01, &
      1.274011611602473639D-01, &
      2.2235277870649807D-02, &
      1.421619193227893466D-03, &
      2.9112874951168792D-05, &
      2.307344176494017303D-02 /)
    real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
      1.28426009614491121D+00, &
      4.68238212480865118D-01, &
      6.59881378689285515D-02, &
      3.78239633202758244D-03, &
      7.29751555083966205D-05 /)
    real ( kind = 8 ), parameter :: root32 = 5.656854248D+00
    real ( kind = 8 ), parameter :: sixten = 16.0D+00
    real ( kind = 8 ) temp
    real ( kind = 8 ), parameter :: sqrpi = 3.9894228040143267794D-01
    real ( kind = 8 ), parameter :: thrsh = 0.66291D+00
    real ( kind = 8 ) x
    real ( kind = 8 ) xden
    real ( kind = 8 ) xnum
    real ( kind = 8 ) y
    real ( kind = 8 ) xsq
!
!  Machine dependent constants
!
    eps = epsilon ( 1.0D+00 ) * 0.5D+00

    x = arg
    y = abs ( x )

    if ( y <= thrsh ) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
      if ( eps < y ) then
        xsq = x * x
      else
        xsq = 0.0D+00
      end if

      xnum = a(5) * xsq
      xden = xsq
      do i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
      end do
      cum = x * ( xnum + a(4) ) / ( xden + b(4) )
      temp = cum
      cum = 0.5D+00 + temp
      ccum = 0.5D+00 - temp
!
!  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
!
    else if ( y <= root32 ) then

      xnum = c(9) * y
      xden = y
      do i = 1, 7
        xnum = ( xnum + c(i) ) * y
        xden = ( xden + d(i) ) * y
      end do
      cum = ( xnum + c(8) ) / ( xden + d(8) )
      xsq = aint ( y * sixten ) / sixten
      del = ( y - xsq ) * ( y + xsq )
      cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cum
      ccum = 1.0D+00 - cum

      if ( 0.0D+00 < x ) then
        call r8_swap ( cum, ccum )
      end if
!
!  Evaluate ANORM for sqrt(32) < |X|.
!
    else

      cum = 0.0D+00
      xsq = 1.0D+00 / ( x * x )
      xnum = p(6) * xsq
      xden = xsq
      do i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
      end do

      cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
      cum = ( sqrpi - cum ) / y
      xsq = aint ( x * sixten ) / sixten
      del = ( x - xsq ) * ( x + xsq )
      cum = exp ( - xsq * xsq * 0.5D+00 ) &
        * exp ( - del * 0.5D+00 ) * cum
      ccum = 1.0D+00 - cum

      if ( 0.0D+00 < x ) then
        call r8_swap ( cum, ccum )
      end if

    end if

    if ( cum < tiny ( cum ) ) then
      cum = 0.0D+00
    end if

    if ( ccum < tiny ( ccum ) ) then
      ccum = 0.0D+00
    end if

    return
  end

  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  FUNCTION inv(A) RESULT(Ainv)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: A
    REAL(wp), DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv

    REAL(wp), DIMENSION(SIZE(A,1)) :: work  ! work array for LAPACK
    INTEGER, DIMENSION(SIZE(A,1)) :: ipiv   ! pivot indices
    INTEGER :: n, info

    ! External procedures defined in LAPACK
    EXTERNAL DGETRF
    EXTERNAL DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = SIZE(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL DGETRF(n, n, Ainv, n, ipiv, info)

    IF (info /= 0) THEN
      STOP 'Matrix is numerically singular!'
    END IF

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

    IF (info /= 0) THEN
      STOP 'Matrix inversion failed!'
    END IF
  END FUNCTION inv

  PURE FUNCTION eye(n) RESULT(iden)
    INTEGER, INTENT(IN) :: n
    REAL(wp), DIMENSION(n, n) :: iden
    INTEGER :: ix

    iden = 0.0_wp
    DO ix = 1,n
      iden(ix, ix) = 1.0_wp
    END DO
  END FUNCTION eye


END MODULE NL
