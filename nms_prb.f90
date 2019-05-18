program nms_prb
!
!*******************************************************************************
!
!! NMS_PRB calls the NMS tests.
!
  character ( len = 8 ) date
  character ( len = 10 ) time
!
  call date_and_time ( date, time )

  write ( *, * ) ' '
  write ( *, * ) 'NMS_PRB'
  write ( *, * ) '  Run the NMS tests.'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
 
  call test004
  call test005
  call test01
  call test02
  call test03
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09
  call test10
  call test11
  call test12
! call test13
  call test14
  call test15
  call test16
  call test17
  call test18
  call test19
  call test20
  call test21
  call test23
  call test24

  write ( *, * ) ' '
  write ( *, * ) 'NMS_PRB'
  write ( *, * ) '  Normal end of NMS tests.'

  stop
end
subroutine test004
!
!*******************************************************************************
!
!! TEST004 computes an autocorrelation using a direct method.
!
  integer, parameter :: n = 168
!
  real acov(0:n-1)
  real el(0:2*n-1)
  real el_sum
  logical ex
  integer i
  integer j
  integer m
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST004'
  write ( *, * ) '  Compute the autocorrelation of El Nino data'
  write ( *, * ) '  using a direct method.'
!
!  Check to see if the required data file exists.
!
  inquire ( file = 'elnino.dat', exist = ex )

  if ( .not. ex ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST004 - Fatal error!'
    write ( *, * ) '  Cannot find the data file: elnino.dat '
    return
  end if

  open ( unit = 8, file = 'elnino.dat', status = 'old' )
!
!  Read the data, find the mean of the data.
!
  do i = 0, n-1
    read(8,*) el(i)
  end do

  close ( unit = 8 )

  el_sum = sum ( el(0:n-1) )
!
!  Subtract the mean, and append N zeroes.
!
  el(0:n-1) = el(0:n-1) - el_sum / real ( n )

  el(n:2*n-1) = 0.0E+00
!
!  Direct calculation.
!  Only sum as far as there is data.
!  Simple, but slow.
!
  do j = 0, n-1
    acov(j) = 0.0E+00
    do m = 0, n-1-j
      acov(j) = acov(j) + el(m) * el(m+j)
    end do
  end do
!
!  Write the scaled autocorrelation.
!
  write ( *, * ) ' '
  write ( *, * ) '  Autocorrelation by the direct method.'
  write ( *, * ) ' '
  write ( *, '(5e14.6)' ) acov(0:19) / acov(0)

  return
end
subroutine test005
!
!*******************************************************************************
!
!! TEST005 tests CFFTB.
!! TEST005 tests CFFTI.
!! TEST005 tests CFFTF.
!
!  Find the autocorrelation to El Nino data using FFT methods.
!
  integer, parameter :: n = 168
!
  real acov(0:n-1)
  complex cel(0:2*n-1)
  complex corr(0:2*n-1)
  real el(0:2*n-1)
  real el_sum
  logical ex
  integer i
  real wsave(4*(2*n)+15)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST005'
  write ( *, * ) '  For Fourier transforms of complex data,'
  write ( *, * ) ' '
  write ( *, * ) '  CFFTI initializes,'
  write ( *, * ) '  CFFTF forward transforms data,'
  write ( *, * ) '  CFFTB backward transforms coefficient.'
!
!  Check to see if the required data file exists.
!
  inquire ( file = 'elnino.dat', exist = ex )

  if ( .not. ex ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST005 - Fatal error!'
    write ( *, * ) '  Cannot find the data file: elnino.dat '
    return
  end if

  open ( unit = 8, file = 'elnino.dat', status = 'old' )
!
!  Read the data, find the mean of the data.
!
  do i = 0, n-1
    read(8,*) el(i)
  end do

  close ( unit = 8 )

  el_sum = sum ( el(0:n-1) )
!
!  Subtract the mean, and append N zeroes.
!
  el(0:n-1) = el(0:n-1) - el_sum / real ( n )

  el(n:2*n-1) = 0.0E+00
!
!  Make a complex copy of EL.
!
  cel(0:2*n-1) = cmplx ( el(0:2*n-1), 0.0E+00 )
!
!    compute fft of data of length 2n.
!    compute square of magnitude of transform components and place
!       in complex array as real parts.
!    compute inverse transform, throwing away second half and
!       imaginary parts (which are zero), and multiply by length of
!       sequence, 2n.
!
  call cffti ( 2*n, wsave )

  call cfftf ( 2*n, cel, wsave )
!
!  CFFTF returns unscaled transforms.
!  The actual transforms are output divided by 2*N.
!
  corr(0:2*n-1) = abs ( cel(0:2*n-1) / real ( 2 * n ) ) **2
!
!  Since we compute transform times its conjugate, we must divide by
!  (2n) for each, i.e., (2n)**2.
!
  call cfftb ( 2*n, corr, wsave )

  acov(0:n-1) = real ( corr(0:n-1) ) * real ( 2 * n )
!
!  Autocovariance is the inverse transform times the sequence length, 2*N.
!
!  Normally, all the scaling would be done once
!  by dividing by 2*N.  We've broken it up for exposition.
!
  write ( *, * ) ' '
  write ( *, * ) '  Autocorrelation  by the complex FFT method.'
  write ( *, * ) ' '
  write ( *, '(5e14.6)' ) acov(0:19) / acov(0)

  return
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests EZFFTB.
!! TEST01 tests EZFFTI.
!! TEST01 tests EZFFTB.
!
!  Find the autocorrelation to El Nino data using real FFT methods.
!
  integer, parameter :: n = 168
!
  real a(2*n)
  real acovr(0:2*n-1)
  real azero
  real b(2*n)
  real el(0:2*n-1)
  real el_sum
  logical ex
  integer i
  real wsave(4*(2*n)+15)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  For the FFT of a real data sequence:'
  write ( *, * ) '  EZFFTI initializes,'
  write ( *, * ) '  EZFFTF does forward transforms,'
  write ( *, * ) '  EZFFTB does backward transforms.'
!
!  Check to see if the required data file exists.
!
  inquire ( file = 'elnino.dat', exist = ex )

  if ( .not. ex ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST01 - Fatal error!'
    write ( *, * ) '  Cannot find the data file: elnino.dat '
    return
  end if

  open ( unit = 8, file = 'elnino.dat', status = 'old' )
!
!  Read the data, find the mean of the data.
!
  do i = 0, n-1
    read(8,*) el(i)
  end do

  close ( unit = 8 )

  el_sum = sum ( el(0:n-1) )
!
!  Subtract the mean, and append N zeroes.
!
  el(0:n-1) = el(0:n-1) - el_sum / real ( n )

  el(n:2*n-1) = 0.0E+00
!
!  fft approach (real).
!
!  compute fft of data of length 2n.
!  ezfftf produces correctly scaled a's and b's so no extra scaling
!    needed to get transform.

  call ezffti ( 2*n, wsave )

  call ezfftf ( 2*n, el, azero, a, b, wsave )
!
!  compute array of square of each frequency component and place
!    in cosine array (a's) to be back transformed. set b's to 0.
!  there are n a's, and n b's.
!  note that care must be taken to compute magnitude correctly,
!    0.5*(a(i)**2+b(i)**2) for i < n, twice that for i=n.
!
  azero = azero**2

  a(1:n-1) = ( a(1:n-1)**2 + b(1:n-1)**2 ) / 2.0E+00
  a(n) = a(n)**2 + b(n)**2

  b(1:n) = 0.0E+00
!
!  Compute the back transform, throwing away its second half.
!
  call ezfftb ( 2*n, acovr, azero, a, b, wsave )

  write ( *, * ) ' '
  write ( *, * ) 'ex 11.6: autocorrelation (real fft) output reduced.'
  write ( *, * ) ' '
  write ( *, '(5e14.6)' ) acovr(0:19) / acovr(0)

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests CFFTB_2D.
!! TEST02 tests CFFTF_2D.
!
!  Plot the image and transform of an 8 by 8 unit source
!  in a 64 by 64 array.
!
  integer, parameter :: n = 64
  integer, parameter :: lda = n
!
  real a(n,n)
  real dat
  real ermax
  real err
  integer i
  complex image(lda,n)
  complex image2(lda,n)
  integer j
  real wsave(4*n+15)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  For two dimensional complex data:'
  write ( *, * ) '  CFFTF_2D computes the forward FFT transform;'
  write ( *, * ) '  CFFTB_2D computes the backward FFT transform.'
  write ( *, * ) ' '
!
!  Initialize WSAVE.
!
  call cffti ( n, wsave )
!
!  Set up the data.  
!
!    IMAGE is the original data.  
!
!    IMAGE2 is IMAGE scaled by (-1)**(I+J), to place the Fourier coefficients 
!    in the correct place for viewing.
!
  do i = 1, n
    do j = 1, n

      if ( i >=(n/2-4) .and.  i<=(n/2+4) .and. &
           j>=(n/2-4) .and. j<=(n/2+4) ) then
        a(i,j) = 1.0E+00
      else
        a(i,j) = 0.0E+00
      end if

      image(i,j) = cmplx ( a(i,j), 0.0E+00 )
      image2(i,j) = image(i,j) * (-1)**(i-1+j-1)

    end do
  end do
!
!  Compute the forward Fourier transform of IMAGE and IMAGE2.
!
  call cfftf_2d ( lda, n, image, wsave )

  call cfftf_2d ( lda, n, image2, wsave )
!
!  Compute the magnitude of the components of transforms.
!  The actual transforms are unscaled and need to be divided by N*N
!  to be correct.
!
  a(1:n,1:n) = abs ( image(1:n,1:n) )
!
!  Compute the inverse Fourier transform of IMAGE.
!
  call cfftb_2d ( lda, n, image, wsave )
!
!  The transforms need to be divided by N*N to be correct.
!
  image(1:n,1:n) = image(1:n,1:n) / real ( n**2 )
!
!  See if the result agrees with the original data.
!
  ermax = 0.0E+00

  do i = 1, n
    do j = 1, n

      if ( i >=(n/2-4) .and. i <=(n/2+4) .and. &
           j >=(n/2-4) .and. j <=(n/2+4)) then
        dat = 1.0E+00
      else
        dat = 0.0E+00
      end if

      err = abs ( dat - abs ( image(i,j) ) )

      ermax = max ( ermax, err )

    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Maximum error in CFFT2D calculation:'
  write ( *, * ) ' '
  write ( *, * ) ermax

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests CFFTI.
!! TEST03 tests CFFTF.
!
!  Using complex discrete Fourier transform, find the approximate
!  Fourier coefficients to Runge's function on [-1,1] with N=16 and
!  N=17.
!
  complex coeff(0:16)
  real del
  real f
  integer j
  integer n
  real pi
  real runge
  complex sqtm1
  real wsave(150)
  real x0
  real xj
!
  external runge
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  CFFTI initializes the complex FFT routines.'
  write ( *, * ) '  CFFTF does a forward Fourier transform on'
  write ( *, * ) '    complex data.'
  write ( *, * ) ' '

  x0 = -1.0E+00
  sqtm1 = cmplx(0.0E+00,-1.0E+00)

  do n = 16, 17

    call cffti ( n, wsave )
!
!  Function assumed to be periodic on [-1,1], an interval of
!  length 2.
!
    del = 2.0E+00 / real ( n )
    f = 2.0E+00 * pi() / ( real ( n ) * del )
!
!  First sample point at -1, last at 1-del
!
    do j = 0, n-1
      xj = (-1.0E+00)+j*del
      coeff(j) = cmplx(runge(xj),0.0E+00)
    end do

    call cfftf ( n, coeff, wsave )
!
!  Returned coefficients must be divided by N for correct
!  normaliziation.  Also, note repetition after N/2 in original
!  coefficients.  Scaling because X0 not at origin destroys this
!  to some extent.
!
    write ( *, * ) ' '
    write ( *, * ) ' cfftf results for n=' ,n
    write ( *, * ) ' czero=',coeff(0)/n*2
    write ( *, * ) ' j          output from cfftf,             scaled coeffs'
    do j = 1, n-1
       write ( *, '(i5,2e15.6,5x,2e15.6)' ) &
         j, coeff(j), exp(-sqtm1*j*f*x0) * coeff(j)/n *2
    end do

  end do

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests EZFFTI.
!! TEST04 tests EZFFTB.
!! TEST04 tests EZFFTF.
!
!  Using the real discrete fourier transform, find the approximate
!  fourier coefficients to runge's function on [-1,1] with n=16 and
!  n=17.
!
  integer, parameter :: mcoef = 17
!
  real a(mcoef/2)
  real azero
  real b(mcoef/2)
  real c(mcoef/2)
  logical debug
  real del
  real dfta(mcoef/2)
  real dftb(mcoef/2)
  real error
  real f
  integer j
  integer k
  integer m
  integer n
  real pi
  real r(mcoef)
  real runge
  real s(mcoef/2)
  real tn
  real wsave(3*mcoef+15)
  real x
  real x0
  real xj
!
  external runge
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  The "EZ" FFT package:'
  write ( *, * ) '  EZFFTI initializes,'
  write ( *, * ) '  EZFFTF does forward transforms,'
  write ( *, * ) '  EZFFTB does backward transforms.'

  debug = .true.
  x0 = -1.0E+00

  do n = mcoef-1, mcoef

    call ezffti ( n, wsave )
!
! function assumed to be periodic on [-1,1], of length 2.
!
    del = 2.0E+00 / n
    f = 2.0E+00 * pi() / ( real ( n ) * del )
!
! The first sample point is at -1, The last at 1-del
!
    do j = 1, n

      xj = (-1.0E+00) + (j-1) * del
      r(j) = runge ( xj )
!
! compute sines and cosines to adjust output of ezfftf to give
! approximate fourier coefficients.
!
      if ( j <= n/2 ) then
        c(j) = cos ( j * f * x0 )
        s(j) = sin ( j * f * x0 )
      end if

    end do

    call ezfftf ( n, r, azero, a, b, wsave )
!
!  As a convenience this loop can go to n/2. if n is even last b is
! zero.
!
    do j = 1, n/2
      dfta(j) = a(j) * c(j) - b(j) * s(j)
      dftb(j) = a(j) * s(j) + b(j) * c(j)
    end do

    write ( *, * ) ' '
    write ( *, * ) ' ezfftf results for n= ' ,n, ' azero=',azero
    write ( *, * ) ' '
    write ( *, * ) '  j             dfta(j)              dftb(j) '
    write ( *, * ) ' '

    do j = 1, n/2
      write ( *, '(i6,2g14.6)' ) j, dfta(j), dftb(j)
    end do
!
!  Evaluate interpolant at points on [-1,1]
!
    if ( debug ) then

      m = 21

      write ( *, * ) ' '
      write ( *, * ) ' Trigonometric polynomial order n= ', n
      write ( *, * ) ' '
      write ( *, * ) '  X    Interpolant    Runge     Error'
      write ( *, * ) ' '

      do k = 1, m

        x = - 1.0E+00 + 2.0E+00 * ( k - 1.0E+00 ) / real ( m - 1 )

        tn = azero
        do j = 1, n/2
          tn = tn + dfta(j) * cos(j*f*x) + dftb(j) * sin(j*f*x)
        end do

        error = tn - runge ( x )
        write ( *, '(4g14.6)' ) x, tn, runge(x), error

      end do

    end if

  end do

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests FMIN.
!
  real a
  real b
  real fmin
  real fp05
  real fx05
  real tol
  real xstar
!
  external fmin
  external fx05
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  FMIN, function minimizer.'
  write ( *, * ) '  Find a minimizer of F(X) = X^3 - 2 * X - 5.'
  write ( *, * ) ' '

  a = 0.1E+00
  b = 0.9E+00
  tol = 1.0E-06

  xstar = fmin ( a, b, fx05, tol )

  write ( *, * ) 'Results:'
  write ( *, * ) '  X* =     ', xstar
  write ( *, * ) '  F(X*) =  ', fx05 ( xstar)
  write ( *, * ) '  F''(X*) = ', fp05 ( xstar )

  return
end
function fx05 ( x )
!
!*******************************************************************************
!
!! FX05 is a function to be minimized.
!
  real fx05
  real x
!
  fx05 = x * ( x * x - 2.0E+00 ) - 5.0E+00

  return
end
function fp05 ( x )
!
!*******************************************************************************
!
!! FP05 is the derivative of a function to be minimized.
!
  real fp05
  real x
!
  fp05 = 3.0E+00 * x * x - 2.0E+00

  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 tests FZERO.
!
  real ae
  real b
  real c
  real fx06
  integer iflag
  real re
!
  external fx06
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  FZERO, single nonlinear equation solver.'
  write ( *, * ) '  F(X) = X^3 - 2 * X - 5'
  write ( *, * ) ' '

  b = 2.0E+00
  c = 3.0E+00
  ae = 1.0E-06
  re = 1.0E-06

  write ( *, * ) '  Initial interval: ', b, c
  write ( *, * ) '  Absolute error tolerance=',ae
  write ( *, * ) '  Relative error tolerance=',re

  call fzero ( fx06, b, c, c, re, ae, iflag )

  write ( *, * ) ' '
  write ( *, * ) '  FZERO results'
  write ( *, * ) ' '

  if ( iflag /= 1 ) then
    write ( *, * ) '  FZERO returned error code =', iflag
  end if

  write ( *, * ) '  Estimate of zero = ', b
  write ( *, * ) '  Function value=    ', fx06(b)

  return
end
function fx06(x)
!
!*******************************************************************************
!
!! FX06 is a function whose zero is desired.
!
  real fx06
  real x
!
  fx06 = x * ( x * x - 2.0E+00 ) - 5.0E+00

  return
end
subroutine test07
!
!*******************************************************************************
!
!! TEST07 tests RNOR.
!
  integer, parameter :: nbins = 32
  real, parameter :: a = -3.0E+00
  real, parameter :: b = 3.0E+00
  integer, parameter :: nr = 10000
!
  integer h(nbins)
  integer i
  integer inbin
  integer iseed
  integer j
  real r
  real rnor
  real rseed
  real rstart
  real width
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST07'
  write ( *, * ) '  RNOR, normal random number generator.'
  write ( *, * ) ' '

  iseed = 305
  rseed = rstart ( iseed )
  width = ( b - a ) / real  ( nbins - 2 )

  h(1:nbins) = 0

  write ( *, * ) 'Running ', nr, ' normals into ', nbins, ' bins.'

  do i = 1, nr
    r = rnor()
    j = inbin ( r, nbins, a, b, width )
    h(j) = h(j) + 1
  end do

  write ( *, * ) 'Histogram for rnor: number in bin 1,...,32'
  write ( *, * ) '   (-infinity,-3],(-3,-2.8],...,(2.8,3],(3,infinity)'
  write ( *, * ) ' (values are slightly computer dependent)'
  write ( *, '(9i8)' ) h(1:nbins)

  return
end
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests UNCMIN.
!
  integer, parameter :: n = 2
  integer, parameter :: lwork = n*(n+10)
!
  real f
  integer i
  integer ierror
  real work(lwork)
  real x(n)
  real x0(n)
!
  external f8
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST08'
  write ( *, * ) '  UNCMIN, unconstrained minimization code.'
  write ( *, * ) ' '
!
!  Specify an initial estimate of the solution.
!
  x0(1) = 1.0E+00
  x0(2) = 1.0E+00
!
!  Minimize function
!
  call uncmin ( n, x0, f8, x, f, ierror, work, lwork )

  write ( *, * ) ' '
  write ( *, * ) '  UNCMIN return code =', ierror
  write ( *, * ) '  f(x*) =', f
  write ( *, * ) '  x* =', x(1:n)

  return
end
subroutine f8 ( n, x, f )
!
!*******************************************************************************
!
!! F8 is a function to be minimized.
!
  integer, parameter :: m = 4
!
  integer n
!
  real, parameter, dimension ( m ) :: b = &
    (/ 20.0E+00, 9.0E+00, 3.0E+00, 1.0E+00 /)
  real f
  integer j
  real, parameter, dimension ( m ) :: t = &
    (/ 0.0E+00, 1.0E+00, 2.0E+00, 3.0E+00 /)
  real x(n)
!
  f = 0.0E+00
  do j = 1, m
    f = f + ( b(j) - x(1) * exp ( x(2) * t(j) ) )**2
  end do

  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests PCHEV.
!! TEST09 tests PCHEZ.
!! TEST09 tests PCHQA.
!
  integer, parameter :: n = 21
  integer, parameter :: nwk = 42
  integer, parameter :: ne = 101
!
  real a
  real b
  real d(n)
  real error
  real errord
  real f(n)
  real fd(ne)
  real fe(ne)
  integer i
  integer ierr
  real pchqa
  real q
  real runge
  real rungep
  logical spline
  real wk(nwk)
  real x(n)
  real xe(ne)
!
  external runge
  external rungep
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST09'
  write ( *, * ) '  PCHEZ carries out piecewise cubic '
  write ( *, * ) '    Hermite interpolation.'
  write ( *, * ) '  PCHEV evaluates the interpolant.'
  write ( *, * ) '  PCHQA integrates the interpolant.'
  write ( *, * ) ' '
!
!  Compute Runge's function at N points in [-1,1].
!
  do i = 1, n
    x(i) = -1.0E+00 + real ( i - 1 ) / 10.0E+00
  end do

  f(1:n) = runge ( x(1:n) )

  spline = .false.
!
!  Compute cubic hermite interpolant because spline is .false.
!
  call pchez ( n, x, f, d, spline, wk, nwk, ierr )

  if ( ierr < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST09 - Error!'
    write ( *, * ) '  PCHEZ returned error code IERR = ',ierr
    return
  else
    write ( *, * ) ' '
    write ( *, * ) '  PCHEZ has set up the interpolant.'
  end if
!
!  Evaluate interpolant and derivative at NE points from -1 to 0.
!
  do i = 1, ne
    xe(i) = -1.0E+00 + real ( i - 1 ) / real ( ne - 1 )
  end do

  call pchev ( n, x, f, d, ne, xe, fe, fd, ierr )

  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST09 - Error!'
    write ( *, * ) '  PCHEV returned error code IERR = ',ierr
    return
  else
    write ( *, * ) ' '
    write ( *, * ) '  PCHEV has evaluated the interpolant.'
  end if
!
!  Examine the error in the function and derivative.
!
  do i = 1, ne
    error = fe(i) - runge ( xe(i) )
    errord = fd(i) - rungep ( xe(i) )
    write ( *, '(2f8.4,2g12.4)' ) xe(i), fe(i), error, errord
  end do
!
!  Compute the integral over the interval [0,1].
!
  a = 0.0E+00
  b = 1.0E+00
  q = pchqa ( n, x, f, d, a, b, ierr )

  write ( *, * ) ' '
  write ( *, * ) '  PCHQA estimates the integral from A to B.'
  write ( *, * ) '  A = ', a
  write ( *, * ) '  B = ', b
  write ( *, * ) '  Integral estiamte = ', q
  write ( *, * ) '  Return code IERR = ', ierr

  return
end
subroutine test10
!
!*******************************************************************************
!
!! TEST10 tests Q1DA.
!
  real a
  real b
  real e
  real eps
  real f10
  integer iflag
  integer kf
  real r
!
  external f10
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST10'
  write ( *, * ) '  Q1DA, quadrature routine.'
  write ( *, * ) ' '

  a = 0.0E+00
  b = 1.0E+00

  eps = 0.001E+00

  call q1da ( f10, a, b, eps, r, e, kf, iflag )

  write ( *, * ) 'Q1DA results: a, b, eps, r, e, kf, iflag'
  write ( *, '(3f7.4,2e16.8,2i4)' ) a,b,eps,r,e,kf,iflag

  return
end
function f10(x)
!
!*******************************************************************************
!
!! F10 is a function to be integrated.
!
  real f10
  real x
!
  f10 = sin ( 2.0E+00 * x ) - sqrt ( x )

  return
end
subroutine test11
!
!*******************************************************************************
!
!! TEST11 tests QAGI.
!
!  compute integral of exp(-x)*cos(x*x)**2 on [0,infinity)
!  Correct result is 0.70260...
!
  integer, parameter :: limit = 100
  integer, parameter :: lenw = limit * 4
!
  real abserr
  real bound
  real epsabs
  real epsrel
  real f11
  integer ier
  integer inf
  integer iwork(limit)
  integer last
  integer neval
  real result
  real work(lenw)
!
  external f11
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST11'
  write ( *, * ) '  QAGI estimates an integral on a semi-infinite interval.'

  bound = 0.0E+00
  inf = 1
  epsabs = 0.0E+00
  epsrel = 1.0E-05

  call qagi ( f11, bound, inf, epsabs, epsrel, result, abserr, neval, &
    ier, limit, lenw, last, iwork, work )

  write ( *, * ) ' '
  write ( *, * ) '  Estimated integral =   ', result
  write ( *, * ) '  Estimated error =      ', abserr
  write ( *, * ) '  Function evaluations = ', neval
  write ( *, * ) '  Return code IER =      ', ier

  return
end
function f11 ( x )
!
!*******************************************************************************
!
!! F11 is a function to be integrated.
!
  real f11
  real x
!
  f11 = exp ( -x ) * cos ( x**2 )**2

  return
end
subroutine test12
!
!*******************************************************************************
!
!! TEST12 tests QK15.
!
!  Compute erf(1), i.e. integral of 2/sqrt(pi) * exp(-x*x) from 0
!  to 1.0E+00
!
  real a
  real abserr
  real b
  real f12
  real pi
  real resabs
  real resasc
  real result
!
  external f12
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST12'
  write ( *, * ) '  QK15 estimates an integral using '
  write ( *, * ) '  Gauss-Kronrod integration.'
  write ( *, * ) ' '

  a = 0.0E+00
  b = 1.0E+00

  call qk15 ( f12, a, b, result, abserr, resabs, resasc )

  write ( *, * ) ' '
  write ( *, * ) 'qk15 estimate of erf(1) '
  write ( *, * ) '2.0/sqrt(pi())*result,      abserr'
  write ( *, * ) 2.0/sqrt(pi())*result, abserr

  return
end
function f12 ( x )
!
!*******************************************************************************
!
!! F12 is a function to be integrated.
!
  real f12
  real x
!
  f12 = exp ( - x**2 )

  return
end
subroutine test13
!
!*******************************************************************************
!
!! TEST13 ??
!
!  double integral by two one dimensional subroutines
!
  real a1
  real b1
  real eps
  real eps1
  real err
  real err1
  real errmax
  real g13
  integer iflag1
  integer kf1
  real res1
  real xx
!
  external g13
!
  common /comm13/ errmax, xx
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST13'
  write ( *, * ) '  Demonstration of two-dimensional quadrature.'
  write ( *, * ) ' '
  write ( *, * ) '  This test is cancelled for now.'
  write ( *, * ) '  (We forget what routines Q1 and Q2 are...)'
  errmax = 0.0E+00
  eps = 1.0E-04
  eps1 = 9.0E+00 / 10.0E+00 * eps
  a1 = 0.0E+00
  b1 = 1.0E+00

! call q1 ( g13, a1, b1, eps1, res1, err1, kf1, iflag1 )

  err1 = 0.0E+00
  kf1 = 0
  iflag1 = 0

  err = err1 + ( b1 - a1 ) * errmax
  write ( *, * ) res1, err, kf1, iflag1

  return
end
function g13 ( x )
!
!*******************************************************************************
!
!! G13 ??
!
  real a2
  real b2
  real eps2
  real err2
  real errmax
  real f13
  real g13
  integer iflag2
  integer kf2
  real res2
  real x
  real xx
!
  external f13
!
  common /comm13/ errmax, xx
!
  xx = x
  eps2 = 0.1E+00 * 1.0E-04
  a2 = 0.0E+00
  b2 = 2.0E+00
!
! call q2 ( f13, a2, b2, eps2, res2, err2, kf2, iflag2 )
!
  err2 = 0.0E+00
  res2 = 0.0E+00
!
!  test iflag2
!
  if ( iflag2 /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'error in q2, iflag, kf2=', iflag2, kf2
  end if

  errmax = max ( errmax, err2 )

  g13 = res2

  return
end
function f13(y)
!
!*******************************************************************************
!
!! F13 ??
!
  real errmax
  real f13
  real x
  real y
!
  common /comm13/ errmax, x
!
  f13 = exp ( - x**2 * y**2 )

  return
end
subroutine test14
!
!*******************************************************************************
!
!! TEST14 is the reactor shielding problem.
!
  logical absorb
  real azm
  real d
  real dist2c
  real e
  real ea
  real er
  real et
  integer i
  integer iexit
  real mu
  integer na
  integer npart
  integer nr
  integer nt
  integer ntot
  real sa
  real sr
  real st
  real thick
  real x
  real y
  real z
!
  common /comm14/ thick
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST14'
  write ( *, * ) '  The reactor shielding problem'
  write ( *, * ) '  Monte Carlo simulation using the UNI code'
  write ( *, * ) '  for uniform random numbers.'
  write ( *, * ) ' '
!
!  Set total number of particles.
!
  ntot = 10
!
!  Set the slab thickness.
!
  thick = 2.0E+00
!
!  Initialize the particle counts and energy tallies.
!
  call init ( na, nt, nr, npart, ea, et, er, sa, sr, st )
!
!  Main loop over ntot particles
!
    1 continue

  npart = npart + 1
!
!  Finished, compute and print averages, standard deviations
!
  if ( npart == ntot ) then
    call output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )
    return
  end if
!
!  Source generates a new particle with position, direction, energy
!
  call source ( e, mu, azm, x, y, z )

    2 continue
!
!  Compute-distance-to-collision
!
  d = dist2c ( e )
!
!  Update position to decide if particle exits or collides
!
  call update ( mu, azm, d, x, y, z )

  i = iexit ( x, y, z )
!
!  Returns -1, 0 , 1 for 'out on left', 'in', 'out on right'
!
  if ( i < 0 ) then
!
!  Exit on left (reflected), tally and goto source
!
    nr = nr + 1
    er = er + e
    sr = sr + e * e
    go to 1

  else if ( i > 0 ) then
!
!  Exit on right (transmitted), tally and goto source
!
    nt = nt + 1
    et = et + e
    st = st + e * e
    go to 1

  else if ( i == 0 .and. absorb() ) then
!
!  Collision and absorbed, tally and goto source
!
    na = na + 1
    ea = ea + e
    sa = sa + e * e
    go to 1

  else
!
!  Collision and scattered, find scattering angle and energy
!
    call scatt ( mu, azm, e )
!
!  Go to compute-distance-to-collision
!
    go to 2

  end if

  return
end
function absorb()
!
!*******************************************************************************
!
!! ABSORB returns TRUE if the particle is absorbed, which occurs
!  with probability PA.
!
  real, parameter :: pa = 0.1E+00
!
  logical absorb
  real uni
!
  if ( uni() <= pa ) then
    absorb = .true.
  else
    absorb = .false.
  end if

  return
end
function cross ( e )
!
!*******************************************************************************
!
!! CROSS returns cross section (fictional) for energy in range [emin,emax]
!
  real cross
  real e
  real s
  real y
!
  s = abs ( sin(100.0E+00*(exp(e)-1.0E+00)) &
    + sin(18.81E+00*(exp(e)-1.0E+00)) )

  y = max ( 0.02E+00, s )
  cross = 10.0E+00 * exp ( -0.1E+00 / y )

  return
end
function dist2c ( e )
!
!*******************************************************************************
!
!! DIST2C returns distance to collision, with exponential distribution
!  with parameter  `cross section'
!
  real cross
  real dist2c
  real e
  real uni
!
  dist2c = - log ( uni() ) / cross ( e )

  return
end
function energy()
!
!*******************************************************************************
!
!! ENERGY returns energy, with distribution const/sqrt(energy) over [emin,emax]
!  use inverse function approach to compute this
!
!  energy min, max in mev
!
  real, parameter :: emin = 1.0E-03
  real, parameter :: emax = 2.5E+00
!
  real c
  real energy
  real uni
!
  external uni
!
  c = 1.0E+00 / ( 2.0E+00 * ( sqrt ( emax ) - sqrt ( emin ) ) )
  energy = ( uni() / (2.0E+00*c) + sqrt(emin) )**2

  return
end
function iexit ( x, y, z )
!
!*******************************************************************************
!
!! IEXIT ??
!
!  returns -1, 0, +1 as particle is   outside on left, inside,
!  or outside on right.
!
  integer iexit
  real thick
  real x
  real y
  real z
!
  common /comm14/  thick
!
  if ( x > thick ) then
    iexit = 1
  else if ( x < 0.0E+00 )then
    iexit = -1
  else
    iexit = 0
  end if

  return
end
subroutine init ( na, nt, nr, npart, ea, et, er, sa, sr, st )
!
!*******************************************************************************
!
!! INIT initializes data for the reactor shielding problem.
!
  real ea
  real er
  real et
  integer na
  integer npart
  integer nr
  integer nt
  real sa
  real sr
  real st
!
  na = 0
  nt = 0
  nr = 0
  npart = 0
  ea = 0.0E+00
  et = 0.0E+00
  er = 0.0E+00
  sa = 0.0E+00
  sr = 0.0E+00
  st = 0.0E+00

  return
end
subroutine output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )
!
!*******************************************************************************
!
!! OUTPUT ??
!
  real ea
  real er
  real et
  integer na
  integer nr
  integer nt
  integer ntot
  real sa
  real sr
  real st
!
  write ( *, * ) ' '
  write ( *, * ) 'Tallies:'

  if ( na > 0 ) then
    ea = ea / na
    sa = sqrt ( sa / na-ea**2 )
  end if

  write ( *, * ) '% absorbed, energy, sd: ', real(na)/ntot*100, ea, sa

  if ( nr > 0 ) then
    er = er/nr
    sr = sqrt(sr/nr-er*er)
  end if

  write ( *, * ) '% reflected, energy, sd: ', real(nr)/ntot*100, er, sr

  if ( nt > 0 ) then
    et = et/nt
    st = sqrt(st/nt-et*et)
  end if

  write ( *, * ) '% transmitted, energy, sd: ', real(nt)/ntot*100, et, st

  write ( *, * ) ' '

  return
end
subroutine scatt ( mu, azm, e )
!
!*******************************************************************************
!
!! SCATT returns scattering angle and energy.
!
  real azm
  real e
  real mu
  real pi
  real uni
!
!  Isotropic scattering, i.e., uniform on sphere
!
  mu = -1.0E+00 + 2.0E+00 * uni()
  azm = uni() * 2.0E+00 * pi()
!
!  Find scattering energy, uniform in [.3*e,e]
!
  e = ( uni() * 0.7E+00 + 0.3E+00 ) * e

  return
end
subroutine source ( e, mu, azm, x, y, z )
!
!*******************************************************************************
!
!! SOURCE ??
!
!  emitter returns mu uniform on [0,1] as the cosine of an angle,
!  and azm uniform on [0,2*pi] as asimuthal angle
!
  real azm
  real e
  real energy
  real mu
  real pi
  real uni
  real x
  real y
  real z
!
  mu = uni()
  azm = uni() * 2.0E+00 * pi()
  x = 0.0E+00
  y = 0.0E+00
  z = 0.0E+00
!
!  Returns energy from source energy distribution
!
  e = energy()

  return
end
subroutine update ( mu, azm, d, x, y, z )
!
!*******************************************************************************
!
!! UPDATE ??
!
  real azm
  real d
  real mu
  real r
  real x
  real y
  real z
!
  x = x + d * mu
  r = d * sqrt ( 1.0E+00 - mu**2 )
  y = y + r * cos ( azm )
  z = z + r * sin ( azm )

  return
end
subroutine test15
!
!*******************************************************************************
!
!! TEST15 tests RNOR.
!
  integer i
  integer iseed
  real r
  real rnor
  real rseed
  real rstart
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST15'
  write ( *, * ) '  RNOR generates random normal numbers.'
  write ( *, * ) ' '
!
!  set initial seed
!
  iseed = 305
  rseed = rstart(iseed)
!
!  rstart returns floating echo of iseed
!
  write ( *, * ) 'rnor results:'
  write ( *, * ) iseed, rseed

  do i = 1, 3
    r = rnor()
    write ( *, '(g14.6)' )r
  end do

  return
end
subroutine test16
!
!*******************************************************************************
!
!! TEST16 tests SDRIV1.
!
!
!  An example of the use of the ODE solver SDRIV1.
!
!  Here we solve the simple system
!
!    Y1' = Y2
!    Y2' = -Y1
!
!  with initial conditions
!
!    Y1(0) = 0
!    Y2(0) = 1
!
!  with exact solution
!
!    Y1(T) = SIN(T)
!    Y2(T) = COS(T)

!
!  Set the number of equations
!
  integer, parameter :: n = 2
!
!  Set the size of the workspace vector.
!
  integer, parameter :: lenw = n * n + 11 * n + 225
!
  real eps
  integer i
  integer mstate
  integer nstep
  real pi
  real t
  real tout
  real work(lenw)
  real y(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST16'
  write ( *, * ) '  SDRIV1 is a simple interface to the ODE solver.'
  write ( *, * ) ' '
!
!  Set the error tolerance.
!
  eps = 0.00001E+00
!
!  Set the initial time.
!
  t = 0.0E+00
  tout = t
!
!  Set the initial conditions
!
  y(1) = 0.0E+00
  y(2) = 1.0E+00
!
!  Set the number of steps we will take in the DO loop.
!
  nstep = 12
!
!  Tell SDRIV1 that this is the first call for this problem.
!
  mstate = 1
!
!  Print a header for the results.
!
  write ( *, * ) ' '
  write ( *, * ) 'Results'
  write ( *, * ) ' '
  write ( *, * ) '   t        y(1)      y(2)'
  write ( *, * ) '            sin(t)    cos(t)'
  write ( *, * ) '            error     error'
!
!  Call SDRIV1 NSTEP+1 times.
!
  do i = 0, nstep

    tout = real ( 2 * i ) * pi() / real ( nstep ) 

    call sdriv1 ( n, t, y, tout, mstate, eps, work, lenw )

    write ( *, * ) ' '
    write ( *, '(3f11.5)' ) t, y(1), y(2)
    write ( *, '(11x,2f11.5)' ) sin(t), cos(t)
    write ( *, '(11x,2f11.5)' ) y(1)-sin(t), y(2)-cos(t)
!
!  Cancel the computation if we get any output code but 1 or 2.
!
    if ( mstate /= 1 .and. mstate /= 2 ) then
      write ( *, * ) ' '
      write ( *, * ) 'TEST16 - Fatal error!'
      write ( *, * ) 'SDRIV1 returned MSTATE=',mstate
      write ( *, * ) 'The computation is being cancelled.'
      return
    end if

  end do

  return
end
subroutine f ( n, t, y, ydot )
!
!*******************************************************************************
!
!! F evaluates the right hand sides of the ODE's.
!
  integer n
!
  real t
  real y(n)
  real ydot(n)
!
  ydot(1) = y(2)
  ydot(2) = -y(1)

  return
end
subroutine test17
!
!*******************************************************************************
!
!! TEST17 tests SDRIV2.
!
  integer, parameter :: n = 2
  real, parameter :: h = 10.0E+00
  integer, parameter :: nroot = 1
  integer, parameter :: mint = 2
  integer, parameter :: lw = n*n+10*n+2*nroot+204
  integer, parameter :: liw = 23
  real, parameter :: mass = 0.125E+00
!
  real eps
  real ewt
  real gfun
  integer iw(liw)
  integer ms
  real t
  real tout
  real w(lw)
  real y(n+1)
!
  external fsub
  external gfun
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST17'
  write ( *, * ) '  SDRIV2 is an ODE solver.'
  write ( *, * ) ' '
!
  eps = 1.0E-05
!
! set initial point
!
  t = 0.0E+00
  tout = t
!
! set for pure relative error
!
  ewt = 0.0E+00
!
!  Set the initial conditions.
!
  y(1) = h
  y(2) = 0.0E+00
!
!  Set the parameter value.
!
  y(3) = mass
!
!  Set MS to 1, signaling the beginning of the run.
!
  ms = 1

  write ( *, * ) 'sdriv2 results'
  write ( *, * ) '   t,         y(1),      y(2),     ms '

  do

    call sdriv2 ( n, t, y, fsub, tout, ms, nroot, eps, ewt, mint, w, &
      lw, iw, liw, gfun )

    tout = tout + 0.1E+00

    if ( ms == 5 )then
      write ( *, '(3f11.5,i4,a,f11.5)' ) t, y(1), y(2), ms, ' <-- y=0 at t= ', t
      exit
    else
      write ( *, '(3f11.5,i4)' ) t,y(1),y(2),ms
!
!  stop if any output code but 1 or 2.
!
      if ( ms > 2 ) then
        exit
      end if

    end if

  end do

  return
end
subroutine fsub(n,t,y,ydot)
!
!*******************************************************************************
!
!! FSUB ??
!
  integer n
!
  real, parameter :: g = 32.0E+00
!
  real t
  real y(n)
  real ydot(n)
!
  ydot(1) = y(2)
  ydot(2) = - g - y(2) / y(3)

  return
end
function gfun ( n, t, y, iroot )
!
!*******************************************************************************
!
!! GFUN ??
!
  integer n
!
  real gfun
  integer iroot
  real t
  real y(n)

  gfun = y(1)

  return
end
subroutine test18
!
!*******************************************************************************
!
!! TEST18 tests SNSQE.
!
  integer, parameter :: n = 2
  integer, parameter :: lw = 19
!
  real fvec(n)
  integer iflag
  integer info
  integer iopt
  integer nprint
  real tol
  real w(lw)
  real x(n)
!
  external f18
  external j18
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST18'
  write ( *, * ) '  SNSQE, nonlinear equation system solver.'
!
!  Set the parameters for SNSQE.
!
  tol = 1.0E-05

  x(1:2) = (/ 2.0E+00, 3.0E+00 /)

  call f18 ( n, x, fvec, iflag )

  write ( *, * ) ' '
  write ( *, * ) '  Initial solution estimate X0:'
  write ( *, * ) ' '
  write ( *, '(4x,2g14.6)' ) x(1:2)
  write ( *, * ) ' '
  write ( *, * ) '  Function value F(X0):'
  write ( *, * ) ' '
  write ( *, '(4x,2g14.6)' ) fvec(1:2)

  iopt = 2
  nprint = 0
!
!  Solve the nonlinear equations
!
  call snsqe ( f18, j18, iopt, n, x, fvec, tol, nprint, info, w, lw )
!
!  print results
!
  if ( info /= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SNSQE INFO flag = ', info
  end if

  write ( *, * ) ' '
  write ( *, * ) '  SNSQE solution estimate X:'
  write ( *, * ) ' '
  write ( *, '(4x,2g14.6)' ) x(1:2)
  write ( *, * ) ' '
  write ( *, * ) '  Function value F(X):'
  write ( *, * ) ' '
  write ( *, '(4x,2g14.6)' ) fvec(1:2)

  return
end
subroutine f18 ( n, x, fvec, iflag )
!
!*******************************************************************************
!
!! F18 evaluates a set of nonlinear equations whose zero is sought.
!
  integer n
!
  real fvec(n)
  integer iflag
  real x(n)
!
  fvec(1) = x(1) * x(2) - x(2)**3 - 1.0E+00
  fvec(2) = x(1)**2 * x(2) + x(2) - 5.0E+00

  return
end
subroutine j18 ( n, x, fvec, fjac, ldfjac, iflag )
!
!*******************************************************************************
!
!! J18 is a dummy routine.
!
  integer ldfjac
  integer n
!
  real fjac(ldfjac,n)
  real fvec(n)
  integer iflag
  real x(n)
!
  return
end
subroutine test19
!
!*******************************************************************************
!
!! TEST19 tests SQRLS.
!
  integer, parameter :: mm = 5
  integer, parameter :: nn = 3
!
  real a(mm,nn)
  real b(mm)
  integer i
  integer ind
  integer itask
  integer j
  integer jpvt(nn)
  integer kr
  integer m
  integer n
  real qraux(nn)
  real tol
  real work(nn)
  real x(nn)
!
!  set up least-squares problem
!  quadratic model, equally-spaced points
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST19'
  write ( *, * ) '  SQRLS solves linear systems in the least squares sense.'
  write ( *, * ) ' '
!
  m = 5
  n = 3
  do i = 1, m
    a(i,1) = 1.0E+00
    do j = 2, n
      a(i,j) = a(i,j-1) * i
    end do
  end do

  b(1:5) = (/ 1.0E+00, 2.3E+00, 4.6E+00, 3.1E+00, 1.2E+00 /)

  tol = 1.0E-06

  write ( *, * ) ' '
  write ( *, * ) 'Coefficient matrix'
  write ( *, * ) ' '
  do i = 1, m
    write ( *, '(3f12.6)' ) a(i,1:n)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Right-hand side'
  write ( *, * ) ' '
  write ( *, '(5f12.6)' ) b(1:m)
!
!  Solve least-squares problem
!
  itask = 1
  call sqrls ( a, mm, m, n, tol, kr, b, x, b, jpvt, qraux, work, itask, ind )
!
!  Print results
!
  write ( *, * ) ' '
  write ( *, * ) ' error code =', ind
  write ( *, * ) '  Estimated matrix rank =', kr
  write ( *, * ) ' parameters'
  write ( *, '(3f12.6)' ) x(1:n)
  write ( *, * )   ' residuals'
  write ( *, '(5f12.6)' ) b(1:m)

  return
end
subroutine test20
!
!*******************************************************************************
!
!! TEST20 tests SSVDC.
!
  integer, parameter :: ldx = 8
  integer, parameter :: n = 8
  integer, parameter :: p = 3
  integer, parameter :: ldu = n
  integer, parameter :: ldv = p
  integer, parameter :: job = 11
!
  real c(p)
  real e(p)
  integer i
  integer info
  integer j
  integer jb
  real, dimension ( n ) :: pop = (/ 75.994575E+00, 91.972266E+00, &
    105.710620E+00, 122.775046E+00, 131.669275E+00, 150.697361E+00, & 
    179.323175E+00, 203.235298E+00 /)
  real pop80
  real r
  real relerr
  real ri
  real rsq
  real s(p)
  real sum
  real tol
  real u(ldu,ldu)
  real v(ldv,p)
  real w(n)
  real x(ldx,p)
  real y(n)
  real year
!
!  C contains coefficients of polynomial
!
!    c(1)*1+c(2)*t+c(3)*t*t
!
!  where t=year (1900 etc.).
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST20'
  write ( *, * ) '  SSVDC computes the singular value decomposition.'

  do i = 1, 8
    y(i) = 1900.0E+00 + real ( ( i - 1 ) * 10 )
  end do

  x(1:8,1) = 1.0E+00
  x(1:8,2) = y(1:8)
  x(1:8,3) = y(1:8)**2

  call ssvdc ( x, ldx, n, p, s, e, u, ldu, v, ldv, w, job, info )

  write ( *, * ) ' '
  write ( *, * ) '  Computed singular values: '
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) s(1:p)

  c(1:p) = 0.0E+00
!
!  relerr reflects number of accurate digits in data
!  e.g. 6 digits ==> relerr=1.0E-06, etc.
!  making relerr larger increases residuals
!
  relerr = 1.0E-06
  tol = relerr*s(1)
!
!  multiply u-trans * pop, and solve for coefficients c(i)
!
  do j = 1, p

    if ( s(j) > tol ) then

      sum = 0.0E+00
      do i = 1, n
        sum = sum + pop(i) * u(i,j)
      end do
      sum = sum / s(j)

      do i = 1, p
        c(i) = c(i) + sum * v(i,j)
      end do

    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) '  Computed polynomial coefficients:'
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) c(1:p)
!
!  Evaluate the model (horner's rule) and residuals at
!  year =1900,...,1980
!
  write ( *, * ) ' '
  write ( *, * ) '         Model       True'
  write ( *, * ) 'Year  Population  Population  Error'
  write ( *, * ) ' '

  r = 0.0E+00

  do i = 1, 9

    year = 1900.0E+00 + (i-1) * 10.0E+00

    pop80 = 0.0E+00
    do j = p, 1, -1
      pop80 = year * pop80 + c(j)
    end do

    if ( i < 9 ) then
      r = r + ( pop(i) - pop80 )**2
      write ( *, '(i4,3f10.2)' ) int(year), pop80, pop(i), pop(i) - pop80
    else
      write ( *, '(i4,f10.3)' ) int(year), pop80
    end if

  end do

  r = sqrt ( r )

  write ( *, * ) ' '
  write ( *, * ) '  RMS error is ', r

  return
end
subroutine test21
!
!*******************************************************************************
!
!! TEST21 tests UNCMIN.
!
  integer, parameter :: n = 10
  integer, parameter :: lwork = n*(n+10)
!
  real f
  integer i
  integer ierror
  real work(lwork)
  real x(n)
  real x0(n)
!
  external f21
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST21'
  write ( *, * ) '  UNCMIN carries out unconstrained minimization'
  write ( *, * ) '  of a scalar function of several variables.'

  do i = 1, n
    x0(i) = real ( i ) / real ( n + 1 )
  end do
!
!  minimize function
!
  call uncmin ( n, x0, f21, x, f, ierror, work, lwork )

  write ( *, * ) ' '
  write ( *, * ) ' return code =', ierror
  write ( *, * ) ' f(x*) =', f
  write ( *, * ) ' x* ='
  write ( *, '(5f12.6)' ) x(1:n)

  return
end
subroutine f21(n,x,f)
!
!*******************************************************************************
!
!! F21 is a function to be minimized.
!
  integer n
!
  real f
  integer i
  real t1
  real t2
  real x(n)
!
  t1 = 0.0E+00
  t2 = 0.0E+00
  do i = 2, n
    t1 = t1 + ( x(i) - x(i-1)**2 )**2
    t2 = t2 + ( 1.0E+00 - x(i-1) )**2
  end do

  f = 1.0E+00 + 100.0E+00 * t1 + t2

  return
end
subroutine test23
!
!*******************************************************************************
!
!! TEST23 tests UNI.
!
  integer i
  integer iseed
  real u
  real uni
  real ustart
  real useed
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST23'
  write ( *, * ) '  UNI is a uniform random number generator.'
  write ( *, * ) ' '
!
!  Set the initial seed.
!
  iseed = 305
  useed = ustart ( iseed )
!
!  USTART returns floating echo of iseed.
!
  write ( *, * ) ' '
  write ( *, * ) 'The seed value ISEED is ', iseed
  write ( *, * ) 'The starting value is ', useed

  do i = 1, 1000
    u = uni()
  end do

  write ( *, * ) ' '
  write ( *, * ) 'The 1000-th random number generated is ', u

  return
end
subroutine test24
!
!*******************************************************************************
!
!! TEST24 tests SGEFS.
!
  integer, parameter :: lda = 10
!
  real a(lda,lda)
  real b(lda)
  integer i
  integer ind
  integer itask
  integer iwork(lda)
  integer j
  integer n
  real rcond
  real work(lda)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST24'
  write ( *, * ) '  SGEFS solves a system of linear equations.'
  write ( *, * ) ' '
!
!  Set the number of equations.
!
  n = 3

  itask = 1
!
!  Set the coefficient matrix A.
!
  a(1,1) = 10.0E+00
  a(2,1) = -3.0E+00
  a(3,1) =  5.0E+00
  a(1,2) = -7.0E+00
  a(2,2) =  2.0E+00
  a(3,2) = -1.0E+00
  a(1,3) =  0.0E+00
  a(2,3) =  6.0E+00
  a(3,3) =  5.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'Coefficient matrix A:'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3f12.6)' ) a(i,1:n)
  end do
!
!  Set the right hand side vector B.
!
  b(1:3) = (/ 7.0E+00, 4.0E+00, 6.0E+00 /)

  write ( *, * ) ' '
  write ( *, * ) 'Right-hand side B:'
  write ( *, * ) ' '
  write ( *, '(3f12.6)') b(1:n)
!
!  Solve the linear system A*x=b.
!
  call sgefs ( a, lda, n, b, itask, ind, work, iwork, rcond )

  write ( *, * ) ' '
  write ( *, * ) 'SGEFS results:'
  write ( *, * ) ' '

  if ( ind == -10 ) then
    write ( *, * ) '  Error code =',ind
  else if ( ind < 0 ) then
    write ( *, * ) '  Error code =',ind
    return
  else
    write ( *, * ) 'Estimated number of accurate digits =', ind
  end if

  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, * ) ' '
  write ( *, '(3f12.6)' ) b(1:n)

  return
end
function rungep ( x )
!
!*******************************************************************************
!
! RUNGEP evaluates the derivative of Runge's function.
!
!
  real rungep
  real x
!
  rungep = ( -50.0E+00 * x ) / ( 1.0E+00 + 25.0E+00 * x**2 )**2

  return
end
