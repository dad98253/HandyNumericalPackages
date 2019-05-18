!  polpak_prb.f90  17 December 2000
!
program polpak_prb
!
!*******************************************************************************
!
!! POLPAK_PRB calls the POLPAK test routines.
!
  character ( len = 8 ) date
  character ( len = 10 ) time
!
  call date_and_time ( date, time )

  write ( *, * ) ' '
  write ( *, * ) 'POLPAK_PRB'
  write ( *, * ) '  Tests for POLPAK, which computes the values of'
  write ( *, * ) '  certain special functions and polynomials.'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
 
  call test001
  call test002
  call test003
  call test004
  call test005
  call test006
  call test007
  call test008
  call test0085
  call test009
  call test010

  call test08
  call test09
  call test10
  call test105

  call test11
  call test12
  call test13
  call test14
  call test145
  call test146
  call test15
  call test155
  call test16
  call test165
  call test166
  call test26
  call test17
  call test174
  call test175
  call test18
  call test19
  call test20

  call test21
  call test22
  call test23
  call test24
  call test25

  call test265
  call test266
  call test27
  call test28
  call test29
  call test295
  call test30

  write ( *, * ) ' '
  write ( *, * ) 'POLPAK_PRB'
  write ( *, * ) '  Normal end of POLPAK tests.'

  stop
end
subroutine test001
!
!*******************************************************************************
!
!! TEST001 tests ALIGN_ENUM.
!
  integer, parameter :: m_max = 10
  integer, parameter :: n_max = 10
!
  integer align_enum
  integer i
  integer j
  integer table(0:m_max,0:n_max)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST001'
  write ( *, * ) '  ALIGN_ENUM counts the number of possible'
  write ( *, * ) '  alignments of two biological sequences.'

  do i = 0, m_max
    do j = 0, n_max
      table(i,j) = align_enum ( i, j )
    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Alignment enumeration table:'
  write ( *, * ) ' '
  write ( *, '(2x,5i5,6i8)' ) ( j, j = 0, n_max )
  write ( *, * ) ' '
  do i = 0, m_max
    write ( *, '(i2,5i5,6i8)' ) i, table(i,0:n_max)
  end do

  return
end
subroutine test002
!
!*******************************************************************************
!
!! TEST002 tests BELL.
!
  integer, parameter :: n = 10
!
  integer b(0:n)
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST002'
  write ( *, * ) '  BELL computes the Bell numbers.'
  write ( *, * ) ' '
  write ( *, * ) '  I,  BELL(I)'
  write ( *, * ) ' '
  call bell ( b, n )

  do i = 0, n
    write ( *, '(i4,2x,i10)' )  i, b(i)
  end do
 
  return
end
subroutine test003
!
!*******************************************************************************
!
!! TEST003 tests BENFORD.
!
  real benford
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST003'
  write ( *, * ) '  BENFORD(I) is the Benford probability of the'
  write ( *, * ) '  initial digit sequence I.'
  write ( *, * ) ' '
  write ( *, * ) '  I,  BENFORD(I)'
  write ( *, * ) ' '

  do i = 1, 9
    write ( *, '(i4,2x,g14.6)' )  i, benford(i)
  end do
 
  return
end
subroutine test004
!
!*******************************************************************************
!
!! TEST004 tests BERN;
!! TEST004 tests BERN2;
!! TEST004 tests DBERN3.
!
  integer, parameter :: n = 20
!
  real c1(0:n)
  real c2(0:n)
  double precision dbern3
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST004'
  write ( *, * ) '  BERN computes Bernoulli numbers;'
  write ( *, * ) '  BERN2 computes Bernoulli numbers;'
  write ( *, * ) '  DBERN3 computes Bernoulli numbers.'
  write ( *, * ) ' '
  write ( *, * ) '   I      B1               B2                B3'
  write ( *, * ) ' '
 
  call bern ( n, c1 )
  call bern2 ( n, c2 )
 
  do i = 0, n
    write ( *, '(i6,3g18.10)' ) i, c1(i), c2(i), dbern3(i)
  end do
 
  return
end
subroutine test005
!
!*******************************************************************************
!
!! TEST005 tests BERN_POLY;
!! TEST005 tests BERN_POLY2.
!
  integer, parameter :: n = 15
!
  double precision bern_poly2
  real bx
  double precision bx2
  double precision dx
  integer i
  real x
!
  x = 0.2E+00
  dx = dble ( x )
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST005'
  write ( *, * ) '  BERN_POLY evaluates Bernoulli polynomials;'
  write ( *, * ) '  BERN_POLY2 evaluates Bernoulli polynomials. '
  write ( *, * ) ' '
  write ( *, * ) '  X = ', x
  write ( *, * ) ' '
  write ( *, * ) '  I          BX          BX2'
  write ( *, * ) ' '
 
  do i = 1, n
    call bern_poly ( i, x, bx )
    bx2 = bern_poly2 ( i, dx )
    write ( *, '(i2,2x,2g16.8)' ) i, bx, bx2
  end do
 
  return
end
subroutine test006
!
!*******************************************************************************
!
!! TEST006 tests BETA.
!
  integer, parameter :: n = 5
!
  real b(n,n)
  real beta
  integer i
  integer j
  real x(n)
  real y(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST006'
  write ( *, * ) '  BETA evaluates the Beta(X,Y) function.'
  write ( *, * ) ' '

  do i = 1, n
    x(i) = real ( i ) / real ( n )
  end do

  y(1:n) = x(1:n)

  do i = 1, n
    do j = 1, n
      b(i,j) = beta ( x(i), y(j) )
    end do
  end do

  write ( *, * ) ' '
  write ( *, '(7x,5f5.2)' ) y(1:n) 
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(f5.2,2x,5f5.2)' ) x(i), b(i,1:n)
  end do

  return
end
subroutine test007
!
!*******************************************************************************
!
!! TEST007 tests BP01.
!
  integer, parameter :: n = 10
!
  real bern(0:n)
  integer i
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST007'
  write ( *, * ) '  BP01 evaluates Bernstein polynomials.'
  write ( *, * ) ' '
 
  x = 0.3E+00
 
  call bp01 ( n, bern, x )
 
  write ( *, * ) ' '
  write ( *, * ) '  The Bernstein polynomials of degree ',n
  write ( *, * ) '  at X = ',x
  write ( *, * ) ' '
 
  do i = 0, n
    write ( *, '(i4,2x,g14.6)' )  i, bern(i)
  end do
 
  return
end
subroutine test008
!
!*******************************************************************************
!
!! TEST008 tests BPAB.
!
  integer, parameter :: n = 10
!
  real a
  real b
  real bern(0:n)
  integer i
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST008'
  write ( *, * ) '  BPAB evaluates Bernstein polynomials.'
  write ( *, * ) ' '
  x = 0.3E+00
  a = 0.0E+00
  b = 1.0E+00
  call bpab ( n, bern, x, a, b )
 
  write ( *, * ) '  The Bernstein polynomials of degree ',n
  write ( *, * ) '  based on the interval ', a, ' to ', b
  write ( *, * ) '  at X = ', x
  write ( *, * ) ' '
 
  do i = 0, n
    write ( *, '(i4,2x,g14.6)' )  i, bern(i)
  end do
 
  return
end
subroutine test0085
!
!*******************************************************************************
!
!! TEST0085 tests CARDAN.
!! TEST0085 tests CARDAN_COEF.
!
  integer, parameter :: n_max = 10
!
  real c(0:n_max)
  real cx1
  real cx2(0:n_max)
  integer n
  real s
  real x
!
  s = 1.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST0085'
  write ( *, * ) '  CARDAN_COEFF returns the coefficients of a'
  write ( *, * ) '    Cardan polynomial.'
  write ( *, * ) '  CARDAN evaluates a Cardan polynomial directly.'
  write ( *, * ) ' '
  write ( *, * ) '  We use the parameter S = ', s
  write ( *, * ) ' '
  write ( *, * ) '  Table of polynomial coefficients:'
  write ( *, * ) ' '

  do n = 0, n_max
    call cardan_coeff ( n, s, c )
    write ( *, '(i2,11f7.0)' ) n, c(0:n)
  end do

  s = 0.5E+00
  x = 0.25E+00

  write ( *, * ) ' '
  write ( *, * ) '  Compare CARDAN_COEFF + RPOLY_VAL_HORNER versus CARDAN.'
  write ( *, * ) ' '
  write ( *, * ) '  Evaluate polynomials at X = ', x
  write ( *, * ) '  We use the parameter S = ', s
  write ( *, * ) ' '
  write ( *, * ) '  Order, Horner, Direct'
  write ( *, * ) ' '

  call cardan ( n, x, s, cx2 )

  do n = 0, n_max

    call cardan_coeff ( n, s, c )
    call rpoly_val_horner ( n, c, x, cx1 )

    write ( *, '(i2,2g14.6)' ) n, cx1, cx2(n)

  end do

  return
end
subroutine test009
!
!*******************************************************************************
!
!! TEST009 tests CATALAN.
!
  integer, parameter :: n = 10
!
  integer i
  integer c(0:n)
  integer, dimension(0:n), parameter :: c2 = &
    (/ 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 /)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST009'
  write ( *, * ) '  CATALAN computes Catalan numbers.'
  write ( *, * ) ' '
 
  call catalan ( n, c )
 
  write ( *, * ) '  I, computed Cat(I), correct Cat(I)'
  write ( *, * ) ' '
 
  do i = 0, n
    write ( *, '(i3,2i10)' ) i, c(i), c2(i)
  end do
 
  return
end
subroutine test010
!
!*******************************************************************************
!
!! TEST010 tests CATALAN_ROW.
!
  integer, parameter :: n = 10
!
  integer c(0:n)
  integer i
  integer ido
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST010'
  write ( *, * ) '  CATALAN_ROW computes a row of Catalan''s triangle.'
  write ( *, * ) ' '
 
  ido = 0
 
  do i = 0, n
    call catalan_row ( ido, i, c )
    ido = 1
    write ( *, '(i2,2x,11i6)' ) i, c(0:i)
  end do
 
  return
end
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests CHEBY1.
!
  integer, parameter :: n = 10
!
  real c(0:n)
  integer i
  real x
!
  x = 0.2E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST08'
  write ( *, * ) '  CHEBY1 evaluates Cheybyshev polynomials of the '
  write ( *, * ) '  first kind.'
  write ( *, * ) '  Use X = ', x
  write ( *, * ) ' '
 
  call cheby1 ( n, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests CHEBY2.
!
  integer, parameter :: n = 10
!
  real c(0:n)
  integer i
  real x
!
  x = 0.2E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST09'
  write ( *, * ) '  CHEBY2 evaluates Chebyshev polynomials of the '
  write ( *, * ) '  second kind.'
  write ( *, * ) '  Use X = ', x
  write ( *, * ) ' '
 
  call cheby2 ( n, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test10
!
!*******************************************************************************
!
!! TEST10 tests COMBIN.
!
  real cnk
  integer k
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST10'
  write ( *, * ) '  COMBIN evaluates C(N,K).'
  write ( *, * ) ' '
  write ( *, * ) '   N     K    CNK'
  write ( *, * ) ' '
 
  do n = 0, 4
    do k = 0, n
      call combin ( n, k, cnk )
      write ( *, '(i6,i6,g14.6)' ) n, k, cnk
    end do
  end do
 
  return
end
subroutine test105
!
!*******************************************************************************
!
!! TEST105 tests COMBIN2.
!
  integer cnk
  integer k
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST105'
  write ( *, * ) '  COMBIN2 evaluates C(N,K).'
  write ( *, * ) ' '
  write ( *, * ) '   N     K    CNK'
  write ( *, * ) ' '
 
  do n = 0, 4
    do k = 0, n
      call combin2 ( n, k, cnk )
      write ( *, '(i6,i6,g14.6)' ) n, k, cnk
    end do
  end do
 
  return
end
subroutine test11
!
!*******************************************************************************
!
!! TEST11 tests COMB_ROW.
!
  integer, parameter :: n = 10
!
  integer c(0:n)
  integer i
  integer ido
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST11'
  write ( *, * ) '  COMB_ROW computes a row of Pascal''s triangle.'
  write ( *, * ) ' '
 
  ido = 0
 
  do i = 0, n
    call comb_row ( ido, i, c )
    ido = 1
    write ( *, '(i2,2x,11i5)' ) i, c(0:i)
  end do
 
  return
end
subroutine test12
!
!*******************************************************************************
!
!! TEST12 tests EULER;
!! TEST12 tests DEULER2.
!
  integer, parameter :: n = 16
!
  double precision deuler2
  integer i
  integer e(0:n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST12'
  write ( *, * ) '  EULER computes Euler numbers;'
  write ( *, * ) '  DEULER2 computes Euler numbers.'
  write ( *, * ) ' '

  call euler ( n, e )
 
  write ( *, * ) ' '
  write ( *, * ) 'EULER, DEULER2:'
  write ( *, * ) ' '
  do i = 0, n
    write ( *,'(i3,i12,g17.9)') i, e(i), deuler2(i)
  end do

  return
end
subroutine test13
!
!*******************************************************************************
!
!! TEST13 calls EULERIAN for the Eulerian numbers.
!
  integer, parameter :: n = 7
!
  integer e(n,n)
  integer i
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST13'
  write ( *, * ) '  EULERIAN evaluates Eulerian numbers.'
  write ( *, * ) ' '
 
  call eulerian ( e, n )

  do i = 1, n
    write ( *, '(10i6)' )  e(i,1:n)
  end do
 
  return
end
subroutine test14
!
!*******************************************************************************
!
!! TEST14 tests EULER_POLY.
!
  integer, parameter :: n = 15
!
  double precision euler_poly
  integer i
  double precision x
!
  x = 0.5E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST14'
  write ( *, * ) '  EULER_POLY evaluates Euler polynomials.'
  write ( *, * ) '  Evaluate at X = ', x
  write ( *, * ) ' '
 
  do i = 0, n
    write ( *, '(i2,2x,g14.6)' )  i, euler_poly ( i, x )
  end do
 
  return
end
subroutine test145
!
!*******************************************************************************
!
!! TEST145 tests F_HOFSTADTER.
!
  integer f
  integer f_hofstadter
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST145'
  write ( *, * ) '  F_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, * ) '  F function.'

  write ( *, * ) ' '
  write ( *, * ) '     N   F(N)'
  write ( *, * ) ' '

  do i = 0, 30
    f = f_hofstadter ( i )
    write ( *, '(2i6)' ) i, f
  end do

  return
end
subroutine test146
!
!*******************************************************************************
!
!! TEST146 tests DFACTORIAL.
!! TEST146 tests FACTORIAL.
!! TEST146 tests GAMMA.
!! TEST146 tests LOG_FACTORIAL.
!
  double precision dfactorial
  real factorial
  real f1
  real f2
  double precision f3
  real f4
  real gamma
  real log_factorial
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST146'
  write ( *, * ) '  For the factorial function,'
  write ( *, * ) '  DFACTORIAL computes the value in double precision;'
  write ( *, * ) '  FACTORIAL computes the value in single precision;'
  write ( *, * ) '  GAMMA(N+1) computes the value in single precision;'
  write ( *, * ) '  LOG_FACTORIAL computes the logarithm;'
  write ( *, * ) ' '
  write ( *, * ) '  N       F1        F2        F3        F4'
  write ( *, * ) ' '

  do n = 0, 10
    f1 = exp ( log_factorial ( n ) )
    f2 = factorial ( n )
    f3 = dfactorial ( n )
    f4 = gamma ( real ( n + 1 ) )
    write ( *, '(i4,4g16.6)' ) n, f1, f2, f3, f4
  end do

  return
end
subroutine test15
!
!*******************************************************************************
!
!! TEST15 tests FIBONACCI_DIRECT.
!
  integer f
  integer i
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST15'
  write ( *, * ) '  FIBONACCI_DIRECT evalutes a Fibonacci number directly.'
  write ( *, * ) ' '
 
  n = 20
 
  do i = 1, n
    call fibonacci_direct ( i, f )
    write ( *, '(i6,i10)' ) i, f
  end do
 
  return
end
subroutine test155
!
!*******************************************************************************
!
!! TEST155 tests FIBONACCI_FLOOR.
!
  integer f
  integer i
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST155'
  write ( *, * ) '  FIBONACCI_FLOOR computes the largest Fibonacci number'
  write ( *, * ) '  less than or equal to a given positive integer.'
  write ( *, * ) ' '
  write ( *, * ) '     N  Fibonacci  Index'
  write ( *, * ) ' ' 

  do n = 1, 20
    call fibonacci_floor ( n, f, i )
    write ( *, '(i6,2x,i6,2x,i6)' ) n, f, i
  end do
 
  return
end
subroutine test16
!
!*******************************************************************************
!
!! TEST16 tests FIBONACCI_RECURSIVE.
!
  integer, parameter :: n = 20
!
  integer f(n)
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST16'
  write ( *, * ) '  FIBONACCI_RECURSIVE computes the Fibonacci sequence.'
  write ( *, * ) ' '
 
  call fibonacci_recursive ( n, f )
 
  do i = 1, n
    write ( *, '(i6,i10)' ) i, f(i)
  end do
 
  return
end
subroutine test165
!
!*******************************************************************************
!
!! TEST165 tests G_HOFSTADTER.
!
  integer g
  integer g_hofstadter
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST165'
  write ( *, * ) '  G_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, * ) '  G function.'

  write ( *, * ) ' '
  write ( *, * ) '     N   G(N)'
  write ( *, * ) ' '

  do i = 0, 30
    g = g_hofstadter ( i )
    write ( *, '(2i6)' ) i, g
  end do

  return
end
subroutine test166
!
!*******************************************************************************
!
!! TEST166 tests GAMMA.
!
  integer, parameter :: n = 13
!
  real gamma
  real, dimension(n), parameter :: gx = (/ &
    4.590845E+00, 2.218160E+00, 1.489192E+00, 1.164230E+00, 1.000000E+00, &
    0.918169E+00, 0.887264E+00, 0.893515E+00, 0.931384E+00, 1.000000E+00, &
    3.6288000E+05, 1.2164510E+17, 8.8417620E+30 /)
  integer i
  real, dimension(n), parameter :: x = (/ &
    0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00,  1.0E+00, & 
    1.2E+00,  1.4E+00,  1.6E+00,  1.8E+00,  2.0E+00, &
    10.0E+00, 20.0E+00, 30.0E+00 /)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST166'
  write ( *, * ) '  GAMMA evaluates the Gamma(X) function.'
  write ( *, * ) ' '
  write ( *, * ) '       X          Gamma(X)     Exact'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3g14.6)' ) x(i), gamma(x(i)), gx(i)
  end do

  return
end
subroutine test26
!
!*******************************************************************************
!
!! TEST26 tests GAMMA_LOG.
!
  integer, parameter :: n = 13
!
  real gamma_log
  real, dimension(n), parameter :: gx = (/ &
    4.590845E+00, 2.218160E+00, 1.489192E+00, 1.164230E+00, 1.000000E+00, &
    0.918169E+00, 0.887264E+00, 0.893515E+00, 0.931384E+00, 1.000000E+00, &
    3.6288000E+05, 1.2164510E+17, 8.8417620E+30 /)
  integer i
  real, dimension(n), parameter :: x = (/ &
    0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00,  1.0E+00, & 
    1.2E+00,  1.4E+00,  1.6E+00,  1.8E+00,  2.0E+00, &
    10.0E+00, 20.0E+00, 30.0E+00 /)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST26'
  write ( *, * ) '  GAMMA_LOG computes Log(Gamma(X))'
  write ( *, * ) ' '
  write ( *, * ) '  X, Gamma(X), Log(Gamma(X)), GAMMA_LOG(X)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(4g14.6)' ) x(i), gx(i), log ( gx(i) ), gamma_log ( x(i) )
  end do
 
  return
end
subroutine test17
!
!*******************************************************************************
!
!! TEST17 tests GEGENBAUER.
!
  integer, parameter :: n = 10
!
  real alfa
  real c(0:n)
  integer i
  real x
!
  alfa = 0.5E+00
  x = 1.2E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST17'
  write ( *, * ) '  GEGENBAUER evaluates Gegenbauer polynomials.'
  write ( *, * ) '  Use Alfa = ', alfa, ' X = ', x
  write ( *, * ) ' '
 
  call gegenbauer ( n, alfa, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test174
!
!*******************************************************************************
!
!! TEST174 tests HAIL.
!
  integer hail
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST174'
  write ( *, * ) '  HAIL(I) computes the length of the hail sequence'
  write ( *, * ) '  for I, also known as the 3*N+1 sequence.'
  write ( *, * ) ' '
  write ( *, * ) '  I,  HAIL(I)'
  write ( *, * ) ' '

  do i = 1, 20
    write ( *, '(i4,2x,i6)' )  i, hail(i)
  end do
 
  return
end
subroutine test175
!
!*******************************************************************************
!
!! TEST175 tests H_HOFSTADTER.
!
  integer h
  integer h_hofstadter
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST175'
  write ( *, * ) '  H_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, * ) '  H function.'

  write ( *, * ) ' '
  write ( *, * ) '     N   H(N)'
  write ( *, * ) ' '

  do i = 0, 30
    h = h_hofstadter ( i )
    write ( *, '(2i6)' ) i, h
  end do

  return
end
subroutine test18
!
!*******************************************************************************
!
!! TEST18 tests HERMITE.
!
  integer, parameter :: n = 10
!
  real c(0:n)
  integer i
  real x
!
  x = 0.5E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST18' 
  write ( *, * ) '  HERMITE evaluates the Hermite polynomials.'
  write ( *, * ) '  Use X = ', x
  write ( *, * ) ' '

  call hermite ( n, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
  
  return
end
subroutine test19
!
!*******************************************************************************
!
!! TEST19 tests JACOBI.
!
  integer, parameter :: n = 10
!
  real alfa
  real beta
  real c(0:n)
  integer i
  real x
!
  alfa = -0.5E+00
  beta = 0.5E+00
  x = 0.5E+00
  write ( *, * ) ' '
  write ( *, * ) 'TEST19'
  write ( *, * ) '  JACOBI evaluates the Jacobi polynomials.'
  write ( *, * ) '  Use Alfa = ', alfa, ' Beta=', beta, ' X= ', x
  write ( *, * ) ' '
 
  call jacobi ( n, alfa, beta, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test20
!
!*******************************************************************************
!
!! TEST20 tests LAGUERRE_GEN.
!
  integer, parameter :: n = 10
!
  real alfa
  real c(0:n)
  integer i
  real x
!
  x = 0.5E+00
  alfa = 0.1E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST20'
  write ( *, * ) '  LAGUERRE_GEN evaluates the generalized Laguerre '
  write ( *, * ) '  polynomials.'
  write ( *, * ) '  Use ALFA = ', alfa, ' X = ', x
  write ( *, * ) ' '
 
  call laguerre_gen ( n, alfa, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test21
!
!*******************************************************************************
!
!! TEST21 tests LAGUERRE_LNM.
!
  integer, parameter :: n = 5
  integer, parameter :: m = 1
!
  real c(0:n)
  integer i
  real x
!
  x = 0.5E+00
  write ( *, * ) ' '
  write ( *, * ) 'TEST21'
  write ( *, * ) '  LAGUERRE_LNM evaluates the associated Laguerre polynomials,'
  write ( *, * ) '  Use M = ', m, ' X = ', x
  write ( *, * ) ' '
 
  call laguerre_lnm ( n, m, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test22
!
!*******************************************************************************
!
!! TEST22 tests LAGUERRE.
!
  integer, parameter :: n = 10
!
  real c(0:n)
  integer i
  real x
!
  x = 0.5E+00
  write ( *, * ) ' '
  write ( *, * ) 'TEST22'
  write ( *, * ) '  LAGUERRE evaluates the Laguerre polynomials'
  write ( *, * ) '  Use X = ', x
  write ( *, * ) ' '
 
  call laguerre ( n, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test23
!
!*******************************************************************************
!
!! TEST23 tests LEGENDRE_PN.
!
  integer, parameter :: n = 10
!
  real c(0:n)
  real cp(0:n)
  integer i
  real x
!
  x = 0.5E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST23'
  write ( *, * ) '  LEGENDRE_PN evaluates the Legendre polynomials of'
  write ( *, * ) '  the first kind, and derivatives.'
  write ( *, * ) '  Use X =', x
  write ( *, * ) ' '
 
  call legendre_pn ( n, x, c, cp )
 
  do i = 0, n
    write ( *, '(i6,2g14.6)' ) i, c(i), cp(i)
  end do
  
  return
end
subroutine test24
!
!*******************************************************************************
!
!! TEST24 tests LEGENDRE_PNM.
!
  integer, parameter :: n = 10
  integer, parameter :: m = 1
!
  real c(0:n)
  integer i
  real x
!
  x = 0.5E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST24'
  write ( *, * ) '  LEGENDRE_PNM evaluates the associated Legendre'
  write ( *, * ) '  polynomials of the first kind.'
  write ( *, * ) '  Use M = ', m, ' X = ', x
  write ( *, * ) ' '
 
  call legendre_pnm ( n, m, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test25
!
!*******************************************************************************
!
!! TEST25 tests LEGENDRE_QN.
!
  integer, parameter :: n = 4
!
  real c(0:n)
  integer i
  real x
!
  x = 0.5E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST25'
  write ( *, * ) '  LEGENDRE_QN evaluates Legendre polynomials of the second kind.' 
  write ( *, * ) '  Use X = ', x
  write ( *, * ) ' '
 
  call legendre_qn ( n, x, c )
 
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, c(i)
  end do
 
  return
end
subroutine test265
!
!*******************************************************************************
!
!! TEST265 tests PENTAGON_NUM.
!
  integer n
  integer p
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST265'
  write ( *, * ) '  PENTAGON_NUM computes the pentagonal numbers.'
  write ( *, * ) ' '
 
  do n = 1, 10
    call pentagon_num ( n, p )
    write ( *, '(i4,2x,i6)' ) n, p
  end do
 
  return
end
subroutine test266
!
!*******************************************************************************
!
!! TEST266 tests PYRAMID_NUM.
!
  integer n
  integer pyramid_num
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST266'
  write ( *, * ) '  PYRAMID_NUM computes the pyramidal numbers.'
  write ( *, * ) ' '
 
  do n = 1, 10
    write ( *, '(i4,2x,i6)' ) n, pyramid_num ( n )
  end do
 
  return
end
subroutine test27
!
!*******************************************************************************
!
!! TEST27 tests STIRLING1.
!
  integer, parameter :: m = 8
  integer, parameter :: n = m
!
  integer i
  integer j
  integer s1(m,n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST27'
  write ( *, * ) '  STIRLING1: Stirling numbers of first kind.'
  write ( *, * ) '  Get rows 1 through ', m
  write ( *, * ) ' '
 
  call stirling1 ( m, n, s1 )
 
  do i = 1, m
    write ( *, '(i6,8i8)' ) i, s1(i,1:n)
  end do
 
  return
end
subroutine test28
!
!*******************************************************************************
!
!! TEST28 tests STIRLING2.
!
  integer, parameter :: m = 8
  integer, parameter :: n = m
!
  integer i
  integer j
  integer s2(m,n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST28'
  write ( *, * ) '  STIRLING2: Stirling numbers of second kind.'
  write ( *, * ) '  Get rows 1 through ', m
  write ( *, * ) ' '
 
  call stirling2 ( m, n, s2 )
 
  do i = 1, m
    write ( *, '(i6,8i8)' ) i, s2(i,1:n)
  end do
 
  return
end
subroutine test29
!
!*******************************************************************************
!
!! TEST29 tests TRIANGLE_NUM.
!
  integer n
  integer triangle_num
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST29'
  write ( *, * ) '  TRIANGLE_NUM computes the triangular numbers.'
  write ( *, * ) ' '
 
  do n = 1, 10
    write ( *, '(i4,2x,i6)' ) n, triangle_num ( n )
  end do
 
  return
end
subroutine test295
!
!*******************************************************************************
!
!! TEST295 tests VIBONACCI.
!
  integer, parameter :: n = 20
  integer, parameter :: n_time = 3
!
  integer i
  integer j
  integer v(n,n_time)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST295'
  write ( *, * ) '  VIBONACCI computes a Vibonacci sequence.'
 

  write ( *, * ) ' '
  write ( *, * ) '  Compute the series ', n_time, ' times.'

  write ( *, * ) ' '
  do j = 1, n_time
    call vibonacci ( n, v(1,j) ) 
  end do

  do i = 1, n
    write ( *, '(i6,3i6)' ) i, v(i,1:n_time)
  end do
 
  return
end
subroutine test30
!
!*******************************************************************************
!
!! TEST30 tests ZECKENDORF.
!
  integer i
  integer i_list(20)
  integer f_list(20)
  integer f_sum
  integer m
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST30'
  write ( *, * ) '  ZECKENDORF computes the Zeckendorf decomposition of'
  write ( *, * ) '  an integer into nonconsecutive Fibonacci numbers.'
  write ( *, * ) ' '
  write ( *, * ) '   N Sum M Parts'
  write ( *, * ) ' '

  do n = 1, 100

    call zeckendorf ( n, m, i_list, f_list )

    write ( *, '(i4,2x,15i4)' ) n, f_list(1:m)

  end do

  return
end
