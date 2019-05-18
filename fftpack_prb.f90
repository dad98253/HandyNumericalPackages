program fftpack_prb
!
!*******************************************************************************
!
!! FFTPACK_PRB calls the FFTPACK test routines.
!
  character ( len = 8 ) date
  character ( len = 10 ) time
!
  call date_and_time ( date, time )

  write ( *, * ) ' '
  write ( *, * ) 'FFTPACK_PRB'
  write ( *, * ) '  A set of tests for FFTPACK.'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
 
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
  call test13
  call test14

  write ( *, * ) ' '
  write ( *, * ) 'FFTPACK_PRB'
  write ( *, * ) '  Normal end of FFTPACK tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests EZFFTB.
!! TEST01 tests EZFFTF.
!! TEST01 tests EZFFTI.
!
  integer, parameter :: n = 4096
!
  real a(n/2)
  real azero
  real b(n/2)
  integer i
  real pimach
  real wsave(3*n+15)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  EZFFTI initializes the EZ FFT routines.'
  write ( *, * ) '  EZFFTF does a forward FFT;'
  write ( *, * ) '  EZFFTB does a backward FFT.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data items is N = ', n
!
!  Set data values of periodic function
!
  do i = 1, n
    x(i) = cos ( 6.0E+00 * pimach ( ) * real ( i - 1 ) / real ( n - 1 ) )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) x(i)
  end do
!
!  Initialize the WSAVE array.
!
  call ezffti ( n, wsave )
!
!  Compute FFT
!
  write ( *, * ) ' '
  write ( *, * ) '  Compute FFT coefficients from data.'

  call ezfftf ( n, x, azero, a, b, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 FFT coefficients:'
  write ( *, * ) ' '
  write ( *, '(g14.6)' ) azero
  do i = 1, 10
    write ( *, '(2g14.6)' ) a(i), b(i)
  end do 
!
!  Now compute inverse FFT of coefficients.  Should get back the
!  original data.  First destroy original data so we're sure
!  that the routine had to recreate them!
!
  x(1:n) = 0.0E+00

  write ( *, * ) ' '
  write ( *, * ) '  Retrieve data from FFT coeficients.'

  call ezfftb ( n, x, azero, a, b, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) x(i)
  end do
  
  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests EZFFTB.
!! TEST02 tests EZFFTF.
!! TEST02 tests EZFFTI.
!
  integer, parameter :: n = 3087
!
  real a(n/2)
  real azero
  real b(n/2)
  integer i
  real pimach
  real wsave(3*n+15)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  EZFFTI initializes the EZ FFT routines.'
  write ( *, * ) '  EZFFTF does a forward FFT;'
  write ( *, * ) '  EZFFTB does a backward FFT.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data items is N = ', n
  write ( *, * ) '  which is not a multiple of two!'
!
!  Set data values of periodic function
!
  do i = 1, n
    x(i) = cos ( 6.0E+00 * pimach() * real ( i - 1 ) / real ( n - 1 ) )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '

  do i = 1, 10
    write ( *, '(g14.6)' ) x(i)
  end do
!
!  Initialize the WSAVE array.
!
  call ezffti ( n, wsave )
!
!  Compute FFT
!
  write ( *, * ) ' '
  write ( *, * ) '  Compute FFT coefficients from data.'

  call ezfftf ( n, x, azero, a, b, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 FFT coefficients:'
  write ( *, * ) ' '
  write ( *, '(g14.6)' ) azero
  do i = 1, 10
    write ( *, '(2g14.6)' ) a(i), b(i)
  end do
!
!  Now compute inverse FFT of coefficients.  Should get back the
!  original data.  First destroy original data so we're sure
!  that the routine had to recreate them!
!
  x(1:n) = 0.0E+00

  write ( *, * ) ' '
  write ( *, * ) '  Retrieve data from FFT coeficients.'

  call ezfftb ( n, x, azero, a, b, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) x(i)
  end do

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests RFFTF.
!
!  The input vector is (1,1,1,...,1) which should produce
!  the output (N,0,0,...,0).
!
  integer, parameter :: n = 36
!
  real error
  integer i
  real wsave(2*n+15)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  RFFTF can compute the Fourier transform of a real'
  write ( *, * ) '  vector of data.  In this case, the vector is'
  write ( *, * ) ' '
  write ( *, * ) '    (1,1,1,...,1)'
  write ( *, * ) ' '
  write ( *, * ) '  and the transform should be'
  write ( *, * ) ' '
  write ( *, * ) '    (N,0,0,...,0), '
  write ( *, * ) ' '
  write ( *, * ) '  where N is the number of entries, N = ', n

  x(1:n) = 1.0E+00
!
!  Initialize the WSAVE array.
!
  call rffti ( n, wsave )
!
!  Compute the Fourier transform.
!
  call rfftf ( n, x, wsave )
!
!  Test results.
!
  error = max ( &
    abs ( real ( n ) - x(1) ), & 
    maxval ( abs ( x(2:n) ) ) )

  write ( *, * ) ' '
  write ( *, * ) '  The maximum error in computation is ', error

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests RFFTB.
!
!  The input vector is (1,0,0,...,0) which should produce
!  the output (1,0,0,...,0).
!
  integer, parameter :: n = 36
!
  real error
  integer i
  real wsave(2*n+15)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  RFFTB can recover a real vector of data from Fourier'
  write ( *, * ) '  coefficients.  In this case, the coefficients are:'
  write ( *, * ) ' '
  write ( *, * ) '    (1,0,0,...,0)'
  write ( *, * ) ' '
  write ( *, * ) '  and the data should be:'
  write ( *, * ) ' '
  write ( *, * ) '    (1,1,1,...,1).'

  x(1) = 1.0E+00
  x(2:n) = 0.0E+00
!
!  Initialize the WSAVE array.
!
  call rffti ( n, wsave )
!
!  Compute the inverse Fourier transform.
!
  call rfftb ( n, x, wsave )
!
!  Compute the maximum error.
!
  error = maxval ( abs ( x(1:n) - 1.0E+00 ) )

  write ( *, * ) ' '
  write ( *, * ) '  The maximum error in computation is ', error

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests EZFFTF.
!
  integer, parameter :: n = 36
!
  real a(n/2)
  real azero
  real b(n/2)
  real error
  integer i
  real wsave(3*n+15)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  EZFFTF can take the Fourier transform of a real vector'
  write ( *, * ) '  of data.  In this case, the vector is'
  write ( *, * ) ' '
  write ( *, * ) '    (1,1,1,...,1)'
  write ( *, * ) ' '
  write ( *, * ) '  and the transform should be'
  write ( *, * ) ' '
  write ( *, * ) '    (N,0,0,...,0),'
  write ( *, * ) ' '
  write ( *, * ) '  where N is the number of entries, ', n

  x(1:n) = 1.0E+00
!
!  Initialize the WSAVE array.
!
  call ezffti ( n, wsave )
!
!  Compute the Fourier transform of the data.
!
  call ezfftf ( n, x, azero, a, b, wsave )
!
!  Compute the maximum error.
!
  error = max ( &
             abs ( 1.0E+00 - azero ), &
    maxval ( abs ( a(1:n/2) ) ), &
    maxval ( abs ( b(1:n/2) ) ) )

  write ( *, * ) ' '
  write ( *, * ) '  The maximum error in computation is ', error

  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 tests EZFFTB.
!
!  The input (1,0,0,...,0) should produce the output (1,0,0,...,0).
!
  integer, parameter :: n = 36
!
  real a(n/2)
  real azero
  real b(n/2)
  real error
  real wsave(3*n+15)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  EZFFTB can be used to recover a real data vector'
  write ( *, * ) '  from a Fourier coefficient vector.'
  write ( *, * ) ' '
  write ( *, * ) '  In this test, the Fourier coefficient vector is:'
  write ( *, * ) ' '
  write ( *, * ) '    (1,0,0,...,0)'
  write ( *, * ) ' '
  write ( *, * ) '  and the recovered data vector should be'
  write ( *, * ) ' '
  write ( *, * ) '    (1,1,1,...,1).'

  azero = 1.0E+00
  a(1:n/2) = 0.0E+00
  b(1:n/2) = 0.0E+00
!
!  Initialize the WSAVE array.
!
  call ezffti ( n, wsave )
!
!  Compute the inverse Fourier transform.
!
  call ezfftb ( n, x, azero, a, b, wsave )
!
!  Compute the maximum error.
!
  error = maxval ( abs ( x(1:n) - 1.0E+00 ) )

  write ( *, * ) ' '
  write ( *, * ) '  The maximum error in the computation was ', error

  return
end
subroutine test07
!
!*******************************************************************************
!
!! TEST07 tests CFFTB_2D.
!! TEST07 tests CFFTF_2D.
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
  write ( *, * ) 'TEST07'
  write ( *, * ) '  For two dimensional complex data:'
  write ( *, * ) '  CFFTF_2D computes the forward FFT transform;'
  write ( *, * ) '  CFFTB_2D computes the backward FFT transform.'
  write ( *, * ) ' '
!
!  Initialize the WSAVE array.
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
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests CFFTB.
!! TEST08 tests CFFTF.
!! TEST08 tests CFFTI.
!
  integer, parameter :: n = 4096
!
  real, parameter :: ahi = 5.0E+00
  real, parameter :: alo = 0.0E+00
  complex c(n)
  integer i
  real wsave(4*n+15)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST08'
  write ( *, * ) '  For Fourier analysis of complex data,'
  write ( *, * ) '  CFFTI initializes the FFT routines.'
  write ( *, * ) '  CFFTF does a forward FFT;'
  write ( *, * ) '  CFFTB does a backward FFT.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  call cvec_random ( alo, ahi, n, c )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(2g14.6)' ) c(i)
  end do
!
!  Initialize the WSAVE array.
!
  call cffti ( n, wsave )
!
!  Compute the FFT coefficients.
!
  write ( *, * ) ' '
  write ( *, * ) '  Compute the FFT coefficients from data.'

  call cfftf ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 FFT coefficients:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(2g14.6)' ) c(i)
  end do
!
!  Now compute inverse FFT of coefficients.  Should get back the
!  original data.

  write ( *, * ) ' '
  write ( *, * ) '  Retrieve data from FFT coeficients.'

  call cfftb ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(2g14.6)' ) c(i) / real ( n )
  end do

  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests COSQB.
!! TEST09 tests COSQF.
!! TEST09 tests COSQI.
!
  integer, parameter :: n = 4096
!
  real, parameter :: ahi = 5.0E+00
  real, parameter :: alo = 0.0E+00
  real c(n)
  integer i
  real wsave(3*n+15)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST09'
  write ( *, * ) '  For cosine analysis of real quarter wave data,'
  write ( *, * ) '  COSQI initializes the FFT routines.'
  write ( *, * ) '  COSQF does a forward FFT;'
  write ( *, * ) '  COSQB does a backward FFT.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  call rvec_random ( alo, ahi, n, c )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Initialize the WSAVE array.
!
  call cosqi ( n, wsave )
!
!  Compute the coefficients.
!
  write ( *, * ) ' '
  write ( *, * ) '  Compute the cosine coefficients from data.'

  call cosqf ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 cosine coefficients:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.

  write ( *, * ) ' '
  write ( *, * ) '  Retrieve data from coeficients.'

  call cosqb ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i) / real ( 4 * n )
  end do

  return
end
subroutine test10
!
!*******************************************************************************
!
!! TEST10 tests EZFFTF.
!! TEST10 tests RSFTF.
!
  integer, parameter :: n = 36
!
  real a_f(n/2)
  real a_s(n/2)
  real, parameter :: ahi = 5.0E+00
  real, parameter :: alo = 0.0E+00
  real azero_f
  real azero_s
  real b_f(n/2)
  real b_s(n/2)
  integer i
  real wsave(3*n+15)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST10'
  write ( *, * ) '  For real data,'
  write ( *, * ) '  EZFFTF takes the fast Fourier transform.'
  write ( *, * ) '  RSFTF computes the "slow" Fourier transform.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data values, N = ', n

  call rvec_random ( alo, ahi, n, x )
!
!  Initialize the WSAVE array.
!
  call ezffti ( n, wsave )
!
!  Compute the fast Fourier transform of the data.
!
  call ezfftf ( n, x, azero_f, a_f, b_f, wsave )
!
!  Compute the slow Fourier transform of the data.
!
  call rsftf ( n, x, azero_s, a_s, b_s )

  write ( *, * ) ' '
  write ( *, * ) '  Fast    Slow'
  write ( *, * ) ' '
  write ( *, * ) '  A coefficients:'
  write ( *, * ) ' '

  write ( *, '(i3,2g14.6)' ) 0, azero_f, azero_s

  do i = 1, n/2
    write ( *, '(i3,2g14.6)' ) i, a_f(i), a_s(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  B coefficients:'
  write ( *, * ) ' '

  do i = 1, n/2
    write ( *, '(i3,2g14.6)' ) i, b_f(i), b_s(i)
  end do

  return
end
subroutine test11
!
!*******************************************************************************
!
!! TEST11 tests RSFTB.
!! TEST11 tests RSFTF.
!
  integer, parameter :: n = 36
!
  real a(n/2)
  real, parameter :: ahi = 5.0E+00
  real, parameter :: alo = 0.0E+00
  real azero
  real b(n/2)
  integer i
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST10'
  write ( *, * ) '  For real data,'
  write ( *, * ) '  RSFTF computes the forward transform.'
  write ( *, * ) '  RSFTB computes the backward transform.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data values, N = ', n

  call rvec_random ( alo, ahi, n, x )

  write ( *, * ) ' '
  write ( *, * ) '  First 10 data values:'
  write ( *, * ) ' '

  do i = 1, 10
    write ( *, '(i3,g14.6)' ) i, x(i)
  end do
!
!  Compute the slow Fourier transform of the data.
!
  call rsftf ( n, x, azero, a, b )

  call rsftb ( n, x, azero, a, b )

  write ( *, * ) ' '
  write ( *, * ) '  First 10 recovered data values:'
  write ( *, * ) ' '

  do i = 1, 10
    write ( *, '(i3,g14.6)' ) i, x(i)
  end do

  return
end
subroutine test12
!
!*******************************************************************************
!
!! TEST12 tests SINQB.
!! TEST12 tests SINQF.
!! TEST12 tests SINQI.
!
  integer, parameter :: n = 4096
!
  real, parameter :: ahi = 5.0E+00
  real, parameter :: alo = 0.0E+00
  real c(n)
  integer i
  real wsave(3*n+15)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST12'
  write ( *, * ) '  For sine analysis of real quarter wave data,'
  write ( *, * ) '  SINQI initializes the FFT routines.'
  write ( *, * ) '  SINQF does a forward FFT;'
  write ( *, * ) '  SINQB does a backward FFT.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  call rvec_random ( alo, ahi, n, c )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Initialize the WSAVE array.
!
  call sinqi ( n, wsave )
!
!  Compute the coefficients.
!
  write ( *, * ) ' '
  write ( *, * ) '  Compute the sine coefficients from data.'

  call sinqf ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 sine coefficients:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.

  write ( *, * ) ' '
  write ( *, * ) '  Retrieve data from coeficients.'

  call sinqb ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i) / real ( 4 * n )
  end do

  return
end
subroutine test13
!
!*******************************************************************************
!
!! TEST13 tests COST.
!! TEST13 tests COSTI.
!
  integer, parameter :: n = 4096
!
  real, parameter :: ahi = 5.0E+00
  real, parameter :: alo = 0.0E+00
  real c(n)
  integer i
  real wsave((5*n+30)/2)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST13'
  write ( *, * ) '  For cosine analysis of real data,'
  write ( *, * ) '  COSTI initializes the FFT routines.'
  write ( *, * ) '  COST does a forward or backward FFT.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  call rvec_random ( alo, ahi, n, c )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Initialize the WSAVE array.
!
  call costi ( n, wsave )
!
!  Compute the coefficients.
!
  write ( *, * ) ' '
  write ( *, * ) '  Compute the cosine coefficients from data.'

  call cost ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 cosine coefficients:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.

  write ( *, * ) ' '
  write ( *, * ) '  Retrieve data from coeficients.'

  call cost ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i) / real ( 2 * ( n - 1 ) )
  end do

  return
end
subroutine test14
!
!*******************************************************************************
!
!! TEST14 tests SINT.
!! TEST14 tests SINTI.
!
  integer, parameter :: n = 4096
!
  real, parameter :: ahi = 5.0E+00
  real, parameter :: alo = 0.0E+00
  real c(n)
  integer i
  real wsave((5*n+30)/2)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST14'
  write ( *, * ) '  For sine analysis of real data,'
  write ( *, * ) '  SINTI initializes the FFT routines.'
  write ( *, * ) '  SINT does a forward or backward FFT.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of data items is N = ', n
!
!  Set the data values.
!
  call rvec_random ( alo, ahi, n, c )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Initialize the WSAVE array.
!
  call sinti ( n, wsave )
!
!  Compute the coefficients.
!
  write ( *, * ) ' '
  write ( *, * ) '  Compute the sine coefficients from data.'

  call sint ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 sine coefficients:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i)
  end do
!
!  Now compute inverse transform of coefficients.  Should get back the
!  original data.

  write ( *, * ) ' '
  write ( *, * ) '  Retrieve data from coeficients.'

  call sint ( n, c, wsave )

  write ( *, * ) ' '
  write ( *, * ) '  The first 10 data values:'
  write ( *, * ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) c(i) / real ( 2 * ( n + 1 ) )
  end do

  return
end
