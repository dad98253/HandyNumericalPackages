!  minprb.f90  17 June 2000
!
program minprb
!
!*******************************************************************************
!
!! MINPRB runs the MINPACK tests.
!
!
  write ( *, * ) ' '
  write ( *, * ) 'MINPRB'
  write ( *, * ) '  A set of tests for MINPACK.'

  call test01
  call test02
  call test03
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09

  write ( *, * ) ' '
  write ( *, * ) 'MINPRB'
  write ( *, * ) '  Normal end of MINPACK tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests CHKDER.
!
  integer, parameter :: n = 5
  integer, parameter :: m = n
  integer, parameter :: ldfjac = n
!
  real err(m)
  real fjac(ldfjac,n)
  real fvec(m)
  real fvecp(m)
  integer i
  integer ido
  integer iflag
  integer j
  integer mode
  real x(n)
  real xp(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  CHKDER compares a user supplied jacobian'
  write ( *, * ) '  and a finite difference approximation to it'
  write ( *, * ) '  and judges whether the jacobian is correct.'

  do ido = 1, 2

    if ( ido == 2 ) then

       write ( *, * ) ' '
       write ( *, * ) '  Repeat the test, but use a "bad" jacobian'
       write ( *, * ) '  and see if the routine notices!'
       write ( *, * ) ' '

     end if
!
!  set the point at which the test is to be made:
!
    x(1:n) = 0.5

    write ( *, * ) ' '
    write ( *, * ) '  Evaluation point X:'
    write ( *, * ) ' '
    do i = 1, n
      write ( *, '(g14.6)' ) x(i)
    end do
 
    mode = 1
    call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

    iflag = 1

    call f01 ( n, x, fvec, fjac, ldfjac, iflag )
    call f01 ( n, xp, fvecp, fjac, ldfjac, iflag )

    write ( *, * ) ' '
    write ( *, * ) '  Sampled function values F(X) and F(XP)'
    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(i3,2g14.6)' ) i, fvec(i), fvecp(i)
    end do

    iflag = 2
    call f01 ( n, x, fvec, fjac, ldfjac, iflag )
!
!  here's where we put a mistake into the jacobian, on purpose.
!
    if ( ido == 2 ) then
      fjac(1,1) = 0.5 * fjac(1,1)
    end if

    write ( *, * ) ' '
    write ( *, * ) '  Computed jacobian'
    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(5g14.6)' ) fjac(i,1:n)
    end do

    mode = 2
    call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

    write ( *, * ) ' '
    write ( *, * ) '  CHKDER error estimates:'
    write ( *, * ) '     > 0.5, the gradient component is probably correct.'
    write ( *, * ) '     < 0.5, the gradient component is probably incorrect.'
    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(i6,g14.6)' ) i, err(i)
    end do

  end do

  return
end
subroutine f01 ( n, x, fvec, fjac, ldfjac, iflag )
!
!*******************************************************************************
!
!! F01 is a function/jacobian routine.
!
  integer ldfjac
  integer n
!
  real fjac(ldfjac,n)
  real fvec(n)
  integer i
  integer iflag
  integer j
  integer k
  real prod
  real sum
  real x(n)
!
  if ( iflag == 1 ) then

    sum = - real ( n + 1 )
    do i = 1, n
      sum = sum + x(i)
    end do

    do i = 1, n
      fvec(i) = x(i) + sum
    end do

    prod = 1.0
    do i = 1, n
      prod = prod * x(i)
    end do

    fvec(n) = prod - 1.0

  else if ( iflag == 2 ) then

    do j = 1, n
      do i = 1, n-1
        fjac(i,j) = 1.0
      end do
    end do

    do i = 1, n-1
      fjac(i,i) = 2.0
    end do

    do j=1,n
      prod=1.0
      do k=1,n
        if ( k /= j ) then
          prod=x(k)*prod
        end if
      end do
      fjac(n,j)=prod
    end do

  end if

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests HYBRD1.
!
!
!  This is an example of what your main program would look
!  like if you wanted to use minpack to solve n nonlinear equations
!  in n unknowns.  in this version, we avoid computing the jacobian
!  matrix, and request that minpack approximate it for us.
!
!  the set of nonlinear equations is:
!
!  x1*x1-10*x1+x2*x2+8=0
!  x1*x2*x2+x1-10*x2+8=0
!
!  with solution x1=x2=1
!
  integer, parameter :: n = 2
!
  real fvec(n)
  integer info
  real tol
  real x(n)
!
  external f02
!
  x(1) = 3.0
  x(2) = 0.0
  tol = 0.00001

  call hybrd1 ( f02, n, x, fvec, tol, info )

  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  HYBRD1 solves a nonlinear system of equations.'
  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:n)
  write ( *, * ) ' '

  return
end
subroutine f02 ( n, x, fvec, iflag )
!
!*******************************************************************************
!
!! F02 is a function routine.
!
  integer n
!
  real fvec(n)
  integer iflag
  real x(n)
!
  fvec(1) = x(1)**2 - 10.0 * x(1) + x(2)**2 + 8.0
  fvec(2) = x(1) * x(2)**2 + x(1) - 10.0 * x(2) + 8.0

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests HYBRJ1.
!
  integer, parameter :: n = 2
  integer, parameter :: ldfjac = n
!
  real fjac(ldfjac,n)
  real fvec(n)
  integer info
  real tol
  real x(n)
!
  external f03
!
  x(1) = 3.0
  x(2) = 0.0
  tol = 0.00001

  call hybrj1 ( f03, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  HYBRJ1 solves a nonlinear system of equations.'
  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:n)
  write ( *, * ) ' '

  return
end
subroutine f03 ( n, x, fvec, fjac, ldfjac, iflag )
!
!*******************************************************************************
!
!! F03 is a function/jacobian routine.
!
  integer ldfjac
  integer n
!
  real fjac(ldfjac,n)
  real fvec(n)
  integer iflag
  real x(n)
!
  if ( iflag == 1 ) then
    fvec(1) = x(1)**2 - 10.0 * x(1) + x(2)**2 + 8.0
    fvec(2) = x(1) * x(2)**2 + x(1) - 10.0 * x(2) + 8.0
  else if ( iflag == 2 ) then
    fjac(1,1) = 2.0 * x(1) - 10.0
    fjac(1,2) = 2.0 * x(2)
    fjac(2,1) = x(2)**2 + 1.0
    fjac(2,2) = 2.0 * x(1) * x(2) - 10.0
  end if

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests LMDER1.
!
!
!  LMDER1 solves m nonlinear
!  equations in n unknowns, where m is greater than n.
!  generally, you cannot get a solution vector x which will satisfy
!  all the equations.  that is, the vector equation f(x)=0 cannot
!  be solved exactly.  instead, minpack seeks a solution x so that
!  the euclidean norm transpose(f(x))*f(x) is minimized.  the size
!  of the euclidean norm is a measure of how good the solution is.
!
!  in this example, the set of equations is actually linear, but
!  normally they are nonlinear.
!
!  in this problem, we have a set of pairs of data points, and we
!  seek a functional relationship between them.  we assume the
!  relationship is of the form y=a*x+b and we want to know the
!  values of a and b.  therefore, we would like to find numbers
!  a and b which satisfy a set of equations.
!
!  the data points are (2,2), (4,11), (6,28) and (8,40).
!
!  therefore, the equations we want to satisfy are:
!
!  a*2+b-2=0
!  a*4+b-11=0
!  a*6+b-28=0
!  a*8+b-40=0
!
!  the least squares solution of this system is a=6.55, b=-12.5,
!  in other words, the line y=6.55*x-12.5 is the line which "best"
!  models the data in the least squares sense.
!
!  problems with more variables, or higher degree polynomials, would
!  be solved similarly.  for example, suppose we have (x,y,z) data,
!  and we wish to find a relationship of the form f(x,y,z).  we assume
!  that x and y occur linearly, and z quadratically.  then the equation
!  we seek has the form:
!
!  a*x+b*y+c*z + d*z*z + e = 0
!
!  and, supposing that our first two points were (1,2,3), (1,3,8), our set of
!  equations would begin:
!
!  a*1+b*2+c*3 + d*9  + e = 0
!  a*1+b*3+c*8 + d*64 + e = 0
!
!  and so on.
!
!  M is the number of equations, which in this case is the number of
!  (x,y) data values.
!
!  N is the number of variables, which in this case is the number of
!  'free' coefficients in the relationship we are trying to determine.
!
  integer, parameter :: m = 4
  integer, parameter :: n = 2
  integer, parameter :: ldfjac = m
!
  integer fjac(ldfjac,n)
  real fvec(m)
  integer info
  real tol
  real x(n)
!
  external f04
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  LMDER1 minimizes 4 functions in 2 variables.'

  x(1)=0.0
  x(2)=5.0

  tol = 0.00001

  call lmder1 ( f04, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:m)
  write ( *, * ) ' '

  return
end
subroutine f04 ( m, n, x, fvec, fjac, ldfjac, iflag )
!
!*******************************************************************************
!
!! F04 is a function/jacobian routine.
!
  integer ldfjac
  integer m
  integer n
!
  real fjac(ldfjac,n)
  real fvec(m)
  integer iflag
  real x(n)
  real, dimension ( 4 ) :: xdat = (/ 2.0,  4.0,  6.0,  8.0 /)
  real, dimension ( 4 ) :: ydat = (/ 2.0, 11.0, 28.0, 40.0 /)
!
  if ( iflag == 1 ) then

    fvec(1:m) = x(1) * xdat(1:m) + x(2) - ydat(1:m)

  else if ( iflag == 2 ) then

    fjac(1:m,1) = xdat(1:m)
    fjac(1:m,2) = 1.0

  end if

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests LMDER1.
!
!
!  LMDER1 solves m nonlinear equations in n unknowns, where m is greater 
!  than n.  The functional fit is nonlinear this time,
!  of the form y=a+b*x**c, with x and y data, and a, b and c unknown.
!
!  this problem is set up so that the data is exactly fit by by
!  a=1, b=3, c=2.  normally, the data would only be approximately
!  fit by the best possible solution.
!
  integer, parameter :: m = 10
  integer, parameter :: n = 3
  integer, parameter :: ldfjac = m
!
  real fjac(ldfjac,n)
  real fvec(m)
  integer info
  real tol
  real x(n)
!
  external f05
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  LMDER1 minimizes 10 functions in 3 variables.'

  x(1)=0.0
  x(2)=5.0
  x(3)=1.3

  tol = 0.00001

  call lmder1 ( f05, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:m)
  write ( *, * ) ' '

  return
end
subroutine f05 ( m, n, x, fvec, fjac, ldfjac, iflag )
!
!*******************************************************************************
!
!! F05 is a function/jacobian routine.
!
  integer ldfjac
  integer m
  integer n
!
  real fjac(ldfjac,n)
  real fvec(m)
  integer iflag
  real x(n)
  real, dimension ( 10 ) :: xdat = (/ &
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 /)
  real, dimension ( 10 ) :: ydat = (/ &
    4.0, 13.0, 28.0, 49.0, 76.0, 109.0, 148.0, 193.0, 244.0, 301.0 /)
!
  if ( iflag == 1 ) then

    fvec(1:m) = x(1) + x(2) * xdat(1:m)**x(3) - ydat(1:m)

  else if ( iflag == 2 ) then

    fjac(1:m,1) = 1.0
    fjac(1:m,2) = xdat(1:m)**x(3)
    fjac(1:m,3) = x(2) * log ( xdat(1:m) ) * xdat(1:m)**x(3)

  end if

  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 tests LMDIF1.
!
!
!  LMDIF1 solves m nonlinear
!  equations in n unknowns, where m is greater than n.
!  generally, you cannot get a solution vector x which will satisfy
!  all the equations.  that is, the vector equation f(x)=0 cannot
!  be solved exactly.  instead, minpack seeks a solution x so that
!  the euclidean norm transpose(f(x))*f(x) is minimized.  the size
!  of the euclidean norm is a measure of how good the solution is.
!
!  in this example, the set of equations is actually linear, but
!  normally they are nonlinear.
!
!  in this problem, we have a set of pairs of data points, and we
!  seek a functional relationship between them.  we assume the
!  relationship is of the form y=a*x+b and we want to know the
!  values of a and b.  therefore, we would like to find numbers
!  a and b which satisfy a set of equations.
!
!  the data points are (2,2), (4,11), (6,28) and (8,40).
!
!  therefore, the equations we want to satisfy are:
!
!  a*2+b-2=0
!  a*4+b-11=0
!  a*6+b-28=0
!  a*8+b-40=0
!
!  the least squares solution of this system is a=6.55, b=-12.5,
!  in other words, the line y=6.55*x-12.5 is the line which "best"
!  models the data in the least squares sense.
!
!  problems with more variables, or higher degree polynomials, would
!  be solved similarly.  for example, suppose we have (x,y,z) data,
!  and we wish to find a relationship of the form f(x,y,z).  we assume
!  that x and y occur linearly, and z quadratically.  then the equation
!  we seek has the form:
!
!  a*x+b*y+c*z + d*z*z + e = 0
!
!  and, supposing that our first two points were (1,2,3), (1,3,8), our set of
!  equations would begin:
!
!  a*1+b*2+c*3 + d*9  + e = 0
!  a*1+b*3+c*8 + d*64 + e = 0
!
!  and so on.
!
!  M is the number of equations, which in this case is the number of
!  (x,y) data values.
!
!  N is the number of variables, which in this case is the number of
!  'free' coefficients in the relationship we are trying to determine.
!
  integer, parameter :: m = 4
  integer, parameter :: n = 2
!
  real fvec(m)
  integer info
  real tol
  real x(n)
!
  external f06
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  LMDIF1 minimizes 4 functions in 2 variables.'

  x(1)=0.0
  x(2)=5.0

  tol = 0.00001

  call lmdif1 ( f06, m, n, x, fvec, tol, info )

  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:m)
  write ( *, * ) ' '

  return
end
subroutine f06 ( m, n, x, fvec, iflag )
!
!*******************************************************************************
!
!! F06 is a function routine.
!
  integer m
  integer n
!
  real fvec(m)
  integer iflag
  real x(n)
  real, dimension ( 4 ) :: xdat = (/ 2.0,  4.0,  6.0,  8.0 /)
  real, dimension ( 4 ) :: ydat = (/ 2.0, 11.0, 28.0, 40.0 /)
!
  fvec(1:m) = x(1) * xdat(1:m) + x(2) - ydat(1:m)

  return
end
subroutine test07
!
!*******************************************************************************
!
!! TEST07 tests LMDIF1.
!
!
!  LMDIF1 solves m nonlinear
!  equations in n unknowns, where m is greater than n.  it is similar
!  to test02, except that the functional fit is nonlinear this time,
!  of the form y=a+b*x**c, with x and y data, and a, b and c unknown.
!
!  this problem is set up so that the data is exactly fit by by
!  a=1, b=3, c=2.  normally, the data would only be approximately
!  fit by the best possible solution.
!
  integer, parameter :: m = 10
  integer, parameter :: n = 3
!
  real fvec(m)
  integer info
  real tol
  real x(n)
!
  external f07
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST07'
  write ( *, * ) '  LMDIF1 minimizes 10 functions in 3 variables.'

  x(1)=0.0
  x(2)=5.0
  x(3)=1.3

  tol = 0.00001

  call lmdif1 ( f07, m, n, x, fvec, tol, info )

  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:m)
  write ( *, * ) ' '

  return
end
subroutine f07 ( m, n, x, fvec, iflag )
!
!*******************************************************************************
!
!! F07 is a function routine.
!
  integer m
  integer n
!
  real fvec(m)
  integer iflag
  real x(n)
  real, dimension ( 10 ) :: xdat = (/ &
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 /)
  real, dimension ( 10 ) :: ydat = (/ &
    4.0, 13.0, 28.0, 49.0, 76.0, 109.0, 148.0, 193.0, 244.0, 301.0 /)
!
  fvec(1:m) = x(1) + x(2) * xdat(1:m)**x(3) - ydat(1:m)

  return
end
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests LMSTR1.
!
!
!  LMSTR1 solves m nonlinear
!  equations in n unknowns, where m is greater than n.
!  generally, you cannot get a solution vector x which will satisfy
!  all the equations.  that is, the vector equation f(x)=0 cannot
!  be solved exactly.  instead, minpack seeks a solution x so that
!  the euclidean norm transpose(f(x))*f(x) is minimized.  the size
!  of the euclidean norm is a measure of how good the solution is.
!
!  in this example, the set of equations is actually linear, but
!  normally they are nonlinear.
!
!  in this problem, we have a set of pairs of data points, and we
!  seek a functional relationship between them.  we assume the
!  relationship is of the form y=a*x+b and we want to know the
!  values of a and b.  therefore, we would like to find numbers
!  a and b which satisfy a set of equations.
!
!  the data points are (2,2), (4,11), (6,28) and (8,40).
!
!  therefore, the equations we want to satisfy are:
!
!  a*2+b-2=0
!  a*4+b-11=0
!  a*6+b-28=0
!  a*8+b-40=0
!
!  the least squares solution of this system is a=6.55, b=-12.5,
!  in other words, the line y=6.55*x-12.5 is the line which "best"
!  models the data in the least squares sense.
!
!  problems with more variables, or higher degree polynomials, would
!  be solved similarly.  for example, suppose we have (x,y,z) data,
!  and we wish to find a relationship of the form f(x,y,z).  we assume
!  that x and y occur linearly, and z quadratically.  then the equation
!  we seek has the form:
!
!  a*x+b*y+c*z + d*z*z + e = 0
!
!  and, supposing that our first two points were (1,2,3), (1,3,8), our set of
!  equations would begin:
!
!  a*1+b*2+c*3 + d*9  + e = 0
!  a*1+b*3+c*8 + d*64 + e = 0
!
!  and so on.
!
!  M is the number of equations, which in this case is the number of
!  (x,y) data values.
!
!  N is the number of variables, which in this case is the number of
!  'free' coefficients in the relationship we are trying to determine.
!
  integer, parameter :: m = 4
  integer, parameter :: n = 2
  integer, parameter :: ldfjac = m
!
  real fjac(ldfjac,n)
  real fvec(m)
  integer info
  real tol
  real x(n)
!
  external f08
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST08'
  write ( *, * ) '  LMSTR1 minimizes 4 functions in 2 variables.'

  x(1) = 0.0
  x(2) = 5.0

  tol = 0.00001

  call lmstr1 ( f08, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:m)
  write ( *, * ) ' '

  return
end
subroutine f08 ( m, n, x, fvec, fjrow, iflag )
!
!*******************************************************************************
!
!! F08 is a function/jacobian routine.
!
  integer m
  integer n
!
  real fjrow(n)
  real fvec(m)
  integer iflag
  real x(n)
  real, dimension ( 4 ) :: xdat = (/ 2.0,  4.0,  6.0,  8.0 /)
  real, dimension ( 4 ) :: ydat = (/ 2.0, 11.0, 28.0, 40.0 /)
!
  if ( iflag == 1 ) then

    fvec(1:m) = x(1) * xdat(1:m) + x(2) - ydat(1:m)

  else 

    fjrow(1) = xdat(iflag-1)
    fjrow(2) = 1.0

  end if

  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests LMSTR1.
!
!
!  LMSTR1 solves m nonlinear
!  equations in n unknowns, where m is greater than n.  it is similar
!  to test02, except that the functional fit is nonlinear this time,
!  of the form y=a+b*x**c, with x and y data, and a, b and c unknown.
!
!  this problem is set up so that the data is exactly fit by by
!  a=1, b=3, c=2.  normally, the data would only be approximately
!  fit by the best possible solution.
!
  integer, parameter :: m = 10
  integer, parameter :: n = 3
  integer, parameter :: ldfjac = m
!
  real fjac(ldfjac,n)
  real fvec(m)
  integer info
  real tol
  real x(n)
!
  external f09
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST09'
  write ( *, * ) '  LMSTR1 minimizes 10 functions in 3 variables.'

  x(1)=0.0
  x(2)=5.0
  x(3)=1.3

  tol = 0.00001

  call lmstr1 ( f09, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, * ) ' '
  write ( *, * ) '  Returned value of INFO = ',info
  write ( *, * ) ' '
  write ( *, * ) '  X:'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) x(1:n)
  write ( *, * ) ' '
  write ( *, * ) '  F(X):'
  write ( *, * ) ' '
  write ( *, '(5g14.6)' ) fvec(1:m)
  write ( *, * ) ' '

  return
end
subroutine f09 ( m, n, x, fvec, fjrow, iflag )
!
!*******************************************************************************
!
!! F09 is a function/jacobian routine.
!
  integer m
  integer n
!
  real fjrow(n)
  real fvec(m)
  integer iflag
  real x(n)
  real, dimension ( 10 ) :: xdat = (/ &
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 /)
  real, dimension ( 10 ) :: ydat = (/ &
    4.0, 13.0, 28.0, 49.0, 76.0, 109.0, 148.0, 193.0, 244.0, 301.0 /)
!
  if ( iflag == 1 ) then

    fvec(1:m) = x(1) + x(2) * xdat(1:m)**x(3) - ydat(1:m)

  else

    fjrow(1) = 1.0
    fjrow(2) = xdat(iflag-1)**x(3)
    fjrow(3) = x(2) * log ( xdat(iflag-1) ) * xdat(iflag-1)**x(3)

  end if

  return
end
