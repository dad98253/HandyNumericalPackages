!  splprb.f90  11 August 2000
!
program splprb
!
!**********************************************************************
!
!! SPLPRB calls the SPLPAK tests.
!
  write ( *, * ) ' '
  write ( *, * ) 'SPLPRB'
  write ( *, * ) '  Tests for the spline package SPLPAK.'
 
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
  call test15
  call test16
  call test17
  call test18
  call test19
  call test20

  call test21
  call test22
  call test23
  call test24
  call test25

  write ( *, * ) ' '
  write ( *, * ) 'SPLPRB'
  write ( *, * ) '  Normal end of SPLPAK tests.'

  stop
end
subroutine test01
!
!***********************************************************************
!
!! TEST01 is the Runge example with polynomial interpolation.
!
  integer, parameter :: nmax = 20
!
  real algerp
  real aloger
  real d(nmax)
  real decay
  real dx
  real errmax
  integer i
  integer istep
  integer j
  integer k
  integer n
  real pnatx
  real runge
  real tau(nmax)
  real x
!
  external runge
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  The Runge example'
  write ( *, * ) '  Use polynomial interpolation of order N.'
  write ( *, * ) ' '
  write ( *, * ) '  N   Max error   Decay exponent'
  write ( *, * ) ' '
  write ( *, * ) ' '
!
  algerp = 0.0
  decay = 0.0
  istep = 20

  do n = 2, nmax, 2
!
!  Choose interpolation points TAU equally spaced in (-1,1).
!
    do i = 1, n
      tau(i) =  ( -1.0 * ( n - i ) + 1.0 * ( i - 1 ) ) / real(n-1)
      d(i) = runge(tau(i))
    end do
!
!  Calculate the interpolating polynomial, using divided differences.
!
    do k = 1, n-1
      do i = 1, n-k
        d(i) = ( d(i+1) - d(i) ) / ( tau(i+k)-tau(i) )
      end do
    end do
!
!  Estimate the maximum interpolation error on (-1,1).
!
    errmax = 0.0
!
!  Consider the subinterval ( TAU(I-1), TAU(I) ).
!
    do i = 2, n
!
!  Sample ISTEP points in the subinterval.
!
      do j = 1, istep

        x = ( ( istep - j ) * tau(i-1) + ( j - 1 ) * tau(i) ) / real (istep - 1)

        pnatx = d(1)
        do k = 2, n
          pnatx = d(k) + (x-tau(k)) * pnatx
        end do

        errmax = max ( errmax, abs(runge(x)-pnatx) )

      end do
    end do

    aloger = alog(errmax)

    if ( n > 2 ) then
      decay  = ( aloger - algerp ) / log ( real ( n ) / real ( n - 2 ) )
    end if

    algerp = aloger
    write ( *, '(i4,e12.4,f11.2)' ) n, errmax, decay

  end do

  return
end
subroutine test02
!
!***********************************************************************
!
!! TEST02 is the Runge example, with cubic Hermite interpolation.
!
  integer, parameter :: nmax = 20
!
  real aloger
  real algerp
  real c(4,nmax)
  real decay
  real divdf1
  real divdf3
  real dtau
  real dx
  real errmax
  real h
  integer i
  integer istep
  integer j
  integer n
  real pnatx
  real runge
  real rungep
  real step
  real tau(nmax)
  real x
!
  external runge
  external rungep
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  Runge example'
  write ( *, * ) '  Use cubic Hermite interpolation of order N'
  write ( *, * ) ' '
  write ( *, * ) '  N   Max error   Decay exponent'
  write ( *, * ) ' '
  write ( *, * ) ' '

  algerp = 0.0
  decay = 0.0
  istep = 20
  step = 20.0

  do n = 2, nmax, 2
!
!  Choose equally spaced interpolation points.
!
    h = 2.0 / real ( n-1)

    do i = 1, n
      tau(i) =  ( -1.0 * ( n - i ) + 1.0 * ( i - 1 ) ) / real(n-1)
    end do

    do i = 1, n
      c(1,i) = runge(tau(i))
      c(2,i) = rungep(tau(i))
    end do
!
!  Calculate the coefficients of the polynomial pieces.
!
    call spline_hermite_set ( n, tau, c )
!
!  Estimate maximum interpolation error on (-1,1).
!
    errmax = 0.0
    do i = 2, n

      do j = 1, istep

        x = ( ( istep - j ) * tau(i-1) + ( j - 1 ) * tau(i) ) / real (istep - 1)

        call spline_hermite_val ( n, tau, c, x, pnatx )

        errmax = max(errmax,abs(runge(x)-pnatx))

      end do
    end do

    aloger = log(errmax)

    if ( n > 2 ) then
      decay  = (aloger-algerp)/log(real ( n)/real ( n-2))
    end if

    algerp = aloger
    write ( *, '(i3,e12.4,f11.2)' )n,errmax,decay

  end do

  return
end
subroutine test03
!
!***********************************************************************
!
!! TEST03 is the Runge example with cubic spline interpolation.
!
  integer, parameter :: nmax = 20
!
  real algerp
  real aloger
  real c(4,nmax)
  real decay
  real dx
  real errmax
  integer i
  integer ibcbeg
  integer ibcend
  integer ii
  integer istep
  integer j
  integer k
  integer n
  real pnatx
  real runge
  real rungep
  real rungepp
  real tau(nmax)
  real x
!
  external runge
  external rungep
  external rungepp
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  The Runge example'
  write ( *, * ) '  Use cubic spline interpolation of order N.'
  write ( *, * ) ' '

  do ii = 0, 2

    ibcbeg = ii
    ibcend = ii

    write ( *, * ) ' '
    write ( *, * ) 'Boundary conditions are '
    if ( ii == 0 ) then
      write ( *, * ) 'Not-a-knot.'
    else if ( ii == 1 ) then
      write ( *, *) 'Derivatives at endpoints.'
    else if ( ii == 2 ) then
      write ( *, * ) 'Second derivatives at endpoints.'
    end if

    write ( *, * ) ' '
    write ( *, * ) '  N   Max error   Decay exponent'
    write ( *, * ) ' '
    write ( *, * ) ' '

    algerp = 0.0
    decay = 0.0
    istep = 20

    do n = 2, nmax, 2
!
!  Choose interpolation points TAU equally spaced in (-1,1).
!
      do i = 1, n
        tau(i) = ( -1.0 * ( n - i ) + 1.0 * ( i - 1 ) ) / real(n-1)
      end do

      do i = 1, n
        c(1,i) = runge(tau(i))
      end do

      if ( ii == 0 ) then

      else if ( ii == 1 ) then
        c(2,1) = rungep(tau(1))
        c(2,n) = rungep(tau(n))
      else if ( ii == 2 ) then
        c(2,1) = rungepp(tau(1))
        c(2,2) = rungepp(tau(n))
      end if
!
!  Calculate the cubic spline.
!
      call cubspl ( tau, c, n, ibcbeg, ibcend )
!
!  Estimate the maximum interpolation error on (-1,1).
!
      errmax = 0.0
!
!  Consider the subinterval ( TAU(I-1), TAU(I) ).
!
      do i = 2, n
!
!  Sample ISTEP points in the subinterval.
!
        do j = 1, istep

          x = ( ( istep - j ) * tau(i-1) + ( j - 1 ) * tau(i) ) / real (istep - 1)

          call ppvalu ( tau, c, n-1, 4, x, 0, pnatx )

          errmax = max ( errmax, abs(runge(x)-pnatx) )

        end do
      end do

      aloger = log(errmax)

      if ( n > 2 ) then
        decay  = (aloger-algerp) / log ( real ( n)/real ( n-2) )
      end if

      algerp = aloger
      write ( *, '(i4,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test04
!
!***********************************************************************
!
!! TEST04 is the Runge example with simple Hermite interpolation.
!
  integer, parameter :: nmax = 20
!
  real algerp
  real aloger
  real c(4,nmax)
  real decay
  real dx
  real errmax
  integer i
  integer istep
  integer j
  integer k
  integer n
  real pnatx
  real runge
  real rungep
  real tau(nmax)
  real x
!
  external runge
  external rungep
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  The Runge example'
  write ( *, * ) '  Use simple Hermite interpolation of order N.'
  write ( *, * ) ' '
  write ( *, * ) '  N   Max error   Decay exponent'
  write ( *, * ) ' '
  write ( *, * ) ' '
!
  algerp = 0.0
  decay = 0.0
  istep = 20

  do n = 2, nmax, 2
!
!  Choose interpolation points TAU equally spaced in (-1,1).
!
    do i = 1, n
      tau(i) =  ( -1.0 * ( n - i ) + 1.0 * ( i - 1 ) ) / real(n-1)
    end do

    do i = 1, n
      c(1,i) = runge(tau(i))
      c(2,i) = rungep(tau(i))
    end do
!
!  Calculate the interpolant.
!
    call spline_hermite_set ( n, tau, c )
!
!  Estimate the maximum interpolation error on (-1,1).
!
    errmax = 0.0
!
!  Consider the subinterval ( TAU(I-1), TAU(I) ).
!
    do i = 2, n
!
!  Sample ISTEP points in the subinterval.
!
      do j = 1, istep

        x = ( ( istep - j ) * tau(i-1) + ( j - 1 ) * tau(i) ) / real (istep - 1)

        call spline_hermite_val ( n, tau, c, x, pnatx )

        errmax = max ( errmax, abs(runge(x)-pnatx) )

      end do
    end do

    aloger = log(errmax)

    if ( n > 2 ) then
      decay  = (aloger-algerp) / log ( real ( n)/real ( n-2) )
    end if

    algerp = aloger
    write ( *, '(i4,e12.4,f11.2)' ) n, errmax, decay

  end do

  return
end
subroutine test05
!
!***********************************************************************
!
!! TEST05 is the Runge example with simple cubic spline interpolation.
!
  integer, parameter :: nmax = 20
!
  real algerp
  real aloger
  real c(4,nmax)
  real decay
  real dx
  real errmax
  integer i
  integer istep
  integer j
  integer k
  integer n
  real pnatx
  real runge
  real rungep
  real tau(nmax)
  real x
!
  external runge
  external rungep
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  The Runge example'
  write ( *, * ) '  Use simple cubic spline of order N.'
  write ( *, * ) ' '
  write ( *, * ) '  N   Max error   Decay exponent'
  write ( *, * ) ' '
  write ( *, * ) ' '
!
  algerp = 0.0
  decay = 0.0
  istep = 20

  do n = 2, nmax, 2
!
!  Choose interpolation points TAU equally spaced in (-1,1).
!
    do i = 1, n
      tau(i) =  ( -1.0 * ( n - i ) + 1.0 * ( i - 1 ) ) / real(n-1)
    end do

    do i = 1, n
      c(1,i) = runge(tau(i))
    end do

    c(2,1) = rungep(tau(1))
    c(2,n) = rungep(tau(n))
!
!  Calculate the interpolant.
!
    call cubset ( tau, c, n )
!
!  Estimate the maximum interpolation error on (-1,1).
!
    errmax = 0.0
!
!  Consider the subinterval ( TAU(I-1), TAU(I) ).
!
    do i = 2, n
!
!  Sample ISTEP points in the subinterval.
!
      do j = 1, istep

        x = ( ( istep - j ) * tau(i-1) + ( j - 1 ) * tau(i) ) / real (istep - 1)

        call spline_hermite_val ( n, tau, c, x, pnatx )

        errmax = max ( errmax, abs(runge(x)-pnatx) )

      end do
    end do

    aloger = log(errmax)

    if ( n > 2 ) then
      decay  = (aloger-algerp) / log ( real ( n)/real ( n-2) )
    end if

    algerp = aloger
    write ( *, '(i4,e12.4,f11.2)' ) n, errmax, decay

  end do

  return
end
subroutine test06
!
!***********************************************************************
!
!! TEST06 evaluates the piecewise polynomial form of a B spline.
!
  integer, parameter :: k = 4
  integer, parameter :: n = 7
!
  real bcoef(n)
  real break(5)
  real coef(4,4)
  integer i
  integer l
  real scrtch(4,4)
  real t(n+k)
  real value
  real x
!
!  Set the b-coeffs for  b(4,4,t) ....
!
  data bcoef / 3*0.,1.,3*0. /
!
!  Set the knots.
!
  t(1) = 0.0
  t(2) = 0.0
  t(3) = 0.0

  t(4) = 0.0
  t(5) = 1.0
  t(6) = 3.0
  t(7) = 4.0
  t(8) = 6.0

  t(9) = 6.0
  t(10) = 6.0
  t(11) = 6.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  Piecewise representation'
!
!  Construct the piecewise polynomial representation.
!
  call bsplpp(t,bcoef,n,k,scrtch,break,coef,l)
!
!  As a check, evaluate B(4,4,T) from its piecewise polynomial
!  representation on a fine mesh.
!
!  The values should agree with (some of) those generated in 
!  example 2.
!
  do i = 0, 40
    x = real ( i ) * 0.2 - 1.0
    call ppvalu ( break, coef, l, 4, x, 0, value )
    write ( *, '(2g14.6)' ) x, value
  end do

  return
end
subroutine test07
!
!***********************************************************************
!
!! TEST07 constructs a B spline using BVALUE.
!
  real bcoef(1)
  real bvalue
  integer i
  real t(5)
  real value
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST07'
  write ( *, * ) '  Construct a spline via BVALUE'

  t(1) = 0.0
  t(2) = 1.0
  t(3) = 3.0
  t(4) = 4.0
  t(5) = 6.0

  bcoef(1) = 1.0
!
!  Evaluate B(1,4,t) on a fine mesh.  On (0,6), the values should
!  coincide with those obtained in example 3.
!
  do i = 0, 40
    x = real ( i)*0.2 - 1.0
    value = bvalue(t,bcoef,1,4,x,0)
    write ( *, '(2g14.6)' ) x, value
  end do

  return
end
subroutine test08
!
!***********************************************************************
!
!! TEST08 carries out cubic spline interpolation with good knots
!
  integer, parameter :: nmax = 20
!
  real algerp
  real aloger
  real c(4,nmax)
  real decay
  real dx
  real errmax
  real f02
  real h
  integer i
  integer ibcbeg
  integer ibcend
  integer istep
  integer iter
  integer itermx
  integer j
  integer jj
  integer n
  integer nhigh
  integer nlow
  real pnatx
  real scrtch(2,nmax)
  real step
  real tau(nmax)
  real taunew(nmax)
  real x
!
  external f02
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST08'
  write ( *, * ) '  Cubic spline, good knots'

  do jj = 1, 2

    if ( jj == 1 ) then
      itermx = 0
      nlow = 4
      nhigh = 20
    else
      itermx = 3
      nlow = 4
      nhigh = 20
    end if

    write ( *, * ) ' '
    write ( *, * ) 'Run ',jj
    write ( *, * ) ' '
    write ( *, * ) 'ITERMX = ',itermx
    write ( *, * ) 'NLOW = ',nlow
    write ( *, * ) 'NHIGH = ',nhigh
 
    algerp = 0.0
    decay = 0.0
    istep = 20
    step = 20.0

    write ( *, * ) 'Take ',itermx,' cycles through NEWNOT'
    write ( *, * ) ' '
    write ( *, * ) '  N   Max error   Decay exponent'
    write ( *, * ) ' '
!
!  Loop over N = number of data points.
!
    do n = nlow, nhigh, 2
!
!  Knots are initially equispaced.
!
       h = 2.0/real ( n-1)

       do i = 1, n
         tau(i) = real ( i-1)*h - 1.0
       end do

       iter = 1
!
!  Construct cubic spline interpolant. then,itermx  times,
!  determine new knots from it and find a new interpolant.
!
   11       continue

        do i = 1, n
          c(1,i) = f02(tau(i))
        end do

        ibcbeg = 0
        ibcend = 0
        call cubspl ( tau, c, n, ibcbeg, ibcend )

        if ( iter <= itermx ) then

          iter = iter + 1
          call newnot(tau,c,n-1,4,taunew,n-1,scrtch)

          do i = 1, n
            tau(i) = taunew(i)
          enddo

          go to 11

        endif
!
!  Estimate the maximum interpolation error on (-1,1).
!
      errmax = 0.0
      do i = 1, n-1
        dx = (tau(i+1)-tau(i)) / step
        do j = 1, istep
           h = real ( j) * dx
           x = tau(i) + h
!
!  Why not call PPVALU?
!
           pnatx = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
           errmax = max(errmax,abs(f02(x)-pnatx))
        end do
      end do
!
!  Calculate the decay exponent.
!
      aloger = log(errmax)

      if(n > nlow)then
        decay  = (aloger-algerp)/log(real ( n)/real ( n-2))
      end if

      algerp = aloger
      write ( *, '(i3,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test09
!
!***********************************************************************
!
!! TEST09 demonstrates cubic spline interpolation with good knots.
!
  real algerp
  real aloger
  real c(4,20)
  real decay
  real dx
  real errmax
  real g
  real h
  integer i
  integer ibcbeg
  integer ibcend
  integer irate
  integer istep
  integer j
  integer k
  integer n
  integer nhigh
  integer nlow
  real pnatx
  real step
  real tau(20)
  real x
!
  g(x) = sqrt(x+1.0)
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST09'
  write ( *, * ) '  Cubic spline with good knots'

  do k = 1, 3

    if ( k == 1 ) then
      irate = 8
      nlow = 4
      nhigh = 10
    else if ( k == 2 ) then
      irate = 6
      nlow = 4
      nhigh = 20
    else
      irate = 4
      nlow = 4
      nhigh = 20
    end if

    write ( *, * ) ' '
    write ( *, * ) 'Run number ',k
    write ( *, * ) ' '

    algerp = 0.0
    decay = 0.0
    istep = 20
    step = 20.0

    write ( *, * ) '  N   Max error   Decay exponent'
    write ( *, * ) ' '

    do n = nlow, nhigh, 2

      h = 1.0 / real ( n-1)

      do i = 1, n
        tau(i) = 2.0 * (real ( i-1)*h)**irate-1.
        c(1,i) = g(tau(i))
      end do
!
!  Construct cubic spline interpolant.
!
      ibcbeg = 0
      ibcend = 0
      call cubspl(tau,c,n, ibcbeg, ibcend )
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0
      do i = 1, n-1

        dx = (tau(i+1)-tau(i))/step

        do j = 1, istep
          h = real ( j)*dx
          pnatx = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
          x = tau(i)+h
          errmax = max(errmax,abs(g(x)-pnatx))
        end do

      end do

      aloger = log(errmax)

      if ( n > nlow ) then
        decay  = (aloger-algerp)/log(real ( n)/real ( n-2))
      end if

      algerp = aloger
      write(*,'(i3,e12.4,f11.2)')n,errmax,decay

    end do

  end do

  return
end
subroutine test10
!
!***********************************************************************
!
!! TEST10 demonstrates a quasi-interpolant with good knots.
!
  real algerp
  real aloger
  real bcoef(22)
  real break(20)
  real c(4,20)
  real decay
  real dg
  real ddg
  real dtip1
  real dtip2
  real dx
  real errmax
  real g
  real h
  integer i
  integer irate
  integer irun
  integer istep
  integer j
  integer l
  integer m
  integer n
  integer nlow
  integer nhigh
  real x
  real pnatx
  real scrtch(4,4)
  real step
  real t(26)
  real taui
!
!  g is the function to be approximated,
!  dg is its first, and
!  ddg its second derivative.
!
  g(x) = sqrt(x+1.0)
  dg(x) = 0.5/g(x)
  ddg(x) = -0.5*dg(x)/(x+1.0)

  write ( *, * ) ' '
  write ( *, * ) 'TEST10'
  write ( *, * ) '  Quasi-interpolant'

  do irun = 1, 2

    if ( irun == 1 ) then
      irate = 8
      nlow = 4 
      nhigh = 10
    else
      irate = 6
      nlow = 4
      nhigh = 20
    end if

    write ( *, * ) ' '
    write ( *, * ) 'Run ',irun
    write ( *, * ) ' '

    algerp = 0.0
    decay = 0.0
    istep = 20
    step = 20.0

    write ( *, * ) '  N   Max error   Decay exponent'
!
!  Loop over N = dim(spline(4,t))-2.
!  N is chosen as the parameter in order to afford 
!  comparison with examples in which cubic spline 
!  interpolation at N data points was used.
!
    do n = nlow, nhigh, 2

      h = 1.0 / real ( n-1)
      m = n+2
!
!  Interior knots are equidistributed with respect to the
!  function (x+1)**(1/irate).
!
      do i = 1, 4
        t(i) = -1.0
      end do

      do i = 5, m
        t(i) = 2.0*(real ( i-4)*h)**irate-1.0
      end do

      do i = m+1, m+4
        t(i) = 1.0
      end do
!
!  Construct quasi-interpolant.  
!  bcoef(1) = g(-1.) = 0.
!
      bcoef(1) = 0.0
      dtip2 = t(5)-t(4)
      taui = t(5)
!
!  Special choice of tau(2) to avoid infinite
!  derivatives of g at left endpoint.
!
      bcoef(2) = g(taui)-2.*dtip2*dg(taui)/3. + dtip2**2*ddg(taui)/6.

      do i = 3, m
        taui = t(i+2)
        dtip1 = dtip2
        dtip2 = t(i+3)-t(i+2)
!
!  Formula xii(30) of text is used.
!
        bcoef(i) = g(taui) +(dtip2-dtip1)*dg(taui)/3. -dtip1*dtip2*ddg(taui)/6.

      end do
!
!  Convert to piecewise polynomial representation.
!
      call bsplpp(t,bcoef,m,4,scrtch,break,c,l)
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0

      do i = 1, l

        dx = (break(i+1)-break(i))/step

        do j = 1, istep
          h = real ( j)*dx
          pnatx = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
          errmax = max(errmax,abs(g(break(i)+h)-pnatx))
        end do

      end do
!
!  Calculate the decay exponent.
!
      aloger = log(errmax)
      if ( n > nlow ) then
        decay  = (aloger-algerp)/log(real ( n)/real ( n-2) )
      end if
      algerp = aloger
      write(*,'(i3,e12.4,f11.2)')n,errmax,decay

    end do

  end do

  return
end
subroutine test11
!
!***********************************************************************
!
!! TEST11 demonstrates that a large norm amplifies noise.
!
!  An initially uniform data point distribution of N points is
!  changed ITERMX times by moving the JCLOSE-th data point
!  toward its left neighbor, cutting the distance between the two 
!  by a factor of RATE each time.  
!
!  To demonstrate the corresponding increase in the norm of cubic 
!  spline interpolation at these data points, the data are taken 
!  from a cubic polynomial, which would be interpolated exactly,
!  but with noise of size SIZE added.  
!
!  The resulting noise or error in the interpolant, compared to the 
!  cubic, gives the norm or noise amplification factor and is 
!  printed out together with the diminishing distance H between the
!  two data points.
!
  integer, parameter :: nmax = 200
  integer, parameter :: nmaxt4 = nmax*4
  integer, parameter :: nmaxt7 = nmax*7
  integer, parameter :: nmaxp4 = nmax+4
!
  real amax
  real bcoef(nmax)
  real break(nmax)
  real coef(nmaxt4)
  real dx
  real fx
  real g
  real gtau(nmax)
  real h
  integer i
  integer iflag
  integer irun
  integer istep
  integer iter
  integer itermx
  integer j
  integer jclose
  integer l
  integer n
  real rate
  real round
  real scrtch(nmaxt7)
  real size
  real step
  real t(nmaxp4)
  real tau(nmax)
  real x
!
  common /rount/ size
!
!  function to be interpolated.
!
  g(x) = 1.0 + x + x**2 + x**3
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST11'
  write ( *, * ) 'A large norm amplifies noise.'
  write ( *, * ) ' '

  do irun = 1, 2

    if ( irun == 1 ) then
      n = 7
      itermx = 10
      jclose = 4
      size = 0.000001
      rate = 2.0
    else
      n = 7
      itermx = 10
      jclose = 4
      size = 0.0
      rate = 2.0
    end if

    write ( *, * ) ' '
    write ( *, * ) 'Run ',irun
    write ( *, * ) 'Size of noise  = ',size
    write ( *, * ) ' '

    istep = 20
    step = 20.0

    write ( *, * ) ' '
    write ( *, * ) '    H           Max error'
    write ( *, * ) ' '
!
!  Start with uniform data points.
!
    do i = 1, n
      tau(i) = real ( i-1)/real ( n-1)
    end do
!
!  Set up the knot sequence for not-a-knot end condition.
!
    do i = 1, 4
      t(i) = tau(1)
    end do

    do i = 5, n
      t(i) = tau(i-2)
    enddo

    do i = n+1, n+4
      t(i) = tau(n)
    end do

    do iter = 1, itermx

      do i = 1, n
        gtau(i) = round(g(tau(i)))
      end do

      call splint(tau,gtau,t,n,4,scrtch,bcoef,iflag)
 
      if ( iflag == 2 ) then
        write ( *, * ) ' '
        write ( *, * ) 'Error code IFLAG = 2 returned from SPLINT.'
        return
      end if
 
      call bsplpp(t,bcoef,n,4,scrtch,break,coef,l)
!
!  Calculate the maximum interpolation error.
!
      amax = 0.0
      do i = 4, n
        dx = (break(i-2)-break(i-3))/step
        do j = 2, istep
          x = break(i-2)-dx*real ( j-1)
          call ppvalu(break,coef,l,4,x,0,fx)
          amax = max(amax,abs(fx-g(x)))
        end do
      end do

      h = tau(jclose)-tau(jclose-1)

      write ( *, '(e9.2,e15.3)' ) h, amax
!
!  Move TAU(JCLOSE) toward its left neighbor so as to cut
!  their distance by a factor of RATE.
!
      tau(jclose) = (tau(jclose) +(rate-1.)*tau(jclose-1))/rate

    end do

  end do
 
  return
end
function round(x)
!
!***********************************************************************
!
!! ROUND is called to add some noise to data.
!
  real flip
  real round
  real size
  real x
!
  common /rount/ size
!
  save flip
!
  data flip / -1.0 /
!
  flip = -flip
  round = x + flip * size

  return
end
subroutine test12
!
!***********************************************************************
!
!! TEST12 shows a cubic spline interpolant at knot averages with good knots.
!
  integer, parameter :: nmax = 20
  integer, parameter :: nmaxp6 = nmax+6
  integer, parameter :: nmaxp2 = nmax+2
  integer, parameter :: nmaxt7 = nmax*7
!
  real algerp
  real aloger
  real bcoef(nmaxp2)
  real break(nmax)
  real c(4,nmax)
  real decay
  real dx
  real errmax
  real f02
  real gtau(nmax)
  real h
  integer i
  integer iflag
  integer irun
  integer istep
  integer iter
  integer itermx
  integer j
  integer l
  integer n
  integer nhigh
  integer nlow
  real pnatx
  real scrtch(nmaxt7)
  real step
  real t(nmaxp6)
  real tau(nmax)
  real temp
  real tnew(nmax)
  real x
!
  external f02
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST12'
  write ( *, * ) 'Interpolation at knot averages'

  do irun = 1, 3

    if ( irun == 1 ) then
      itermx = 0
      nlow = 4
      nhigh = 20
    else if ( irun == 2 ) then
     itermx = 3
      nlow = 4
      nhigh = 20
    else
      itermx = 6
      nlow = 4
      nhigh = 20
    end if

    write ( *, * ) ' '
    write ( *, * ) 'Run ',irun
    write ( *, * ) 'We will take ',itermx,' cycles through NEWNOT.'
    write ( *, * ) 'NLOW is ',nlow,' and NHIGH is ',nhigh
    write ( *, * ) ' '

    algerp = 0.0
    decay = 0.0
    istep = 20
    step = 20.0
 
    write ( *, * ) ' '
    write ( *, * ) '  N   Max error   Decay exponent'
    write ( *, * ) ' '
!
!  Loop over N = number of data points.
!
    do n = nlow, nhigh, 2

    4     continue

      h = 2.0 / real ( n-3)

      do i = 1, 4
        t(i) = -1.0
      end do

      do i = 5, n
        t(i) = real ( i-4)*h-1.
      enddo

      do i = n+1, n+4
        t(i) = 1.0
      end do

      iter = 1
!
!  Construct cubic spline interpolant.  Then, ITERMX times,
!  determine new knots from it and find a new interpolant.
!
   11     continue

      do i = 1, n
        tau(i) = (t(i+1)+t(i+2)+t(i+3))/3.0
        gtau(i) = f02(tau(i))
      end do

      call splint(tau,gtau,t,n,4,scrtch,bcoef,iflag)
 
      if ( iflag == 2 ) then
        write ( *, * ) ' '
        write ( *, * ) 'Error code IFLAG = 2 from SPLINT.'
        return
      end if
 
      call bsplpp(t,bcoef,n,4,scrtch,break,c,l)

      if ( iter <= itermx ) then

        iter = iter+1
        call newnot(break,c,l,4,tnew,l,scrtch)

        do i = 2, l
          t(3+i) = tnew(i)
        end do

        go to 11

      end if
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0
      do i = 1, l
        dx = (break(i+1)-break(i))/step
        do j = 1, istep
          h = real ( j)*dx
          x = break(i) + h
          pnatx = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
          errmax = max(errmax,abs(f02(x)-pnatx))
        end do
      end do
!
!  Calculate the decay exponent.
!
      aloger = log(errmax)

      if ( n > nlow ) then
        temp = real ( n ) / real ( n - 2 )
        decay = ( aloger - algerp ) / log ( temp )
      end if

      algerp = aloger

      write ( *, '(i3,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test13
!
!***********************************************************************
!
!! TEST13 is a cubic spline interpolant at knot averages with good knots.  
!  modified around label 4.
!
  integer, parameter :: nmax = 20
  integer, parameter :: nmaxp6 = nmax+6
  integer, parameter :: nmaxp2 = nmax+2
  integer, parameter :: nmaxt7 = nmax*7
!
  real algerp
  real aloger
  real bcoef(nmaxp2)
  real break(nmax)
  real c(4,nmax)
  real decay
  real dx
  real errmax
  real f02
  real gtau(nmax)
  real h
  integer i
  integer iflag
  integer irun
  integer istep
  integer iter
  integer itermx
  integer j
  integer l
  integer n
  integer nhigh
  integer nlow
  real pnatx
  real scrtch(nmaxt7)
  real step
  real t(nmaxp6)
  real tau(nmax)
  real tnew(nmax)
  real x
!
  external f02
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST13'
  write ( *, * ) 'Modified example 2'
  write ( *, * ) ' '

  do irun = 1, 3

    write ( *, * ) ' '
    write ( *, * ) 'Run ',irun

    if ( irun == 1 ) then
      itermx = 0
      nlow = 4
      nhigh = 20
    else if ( irun == 2 ) then
      itermx = 1
      nlow = 4
      nhigh = 20
    else
      itermx = 2
      nlow = 4
      nhigh = 20
    end if

    algerp = 0.0
    decay = 0.0
    istep = 20
    step = 20.0
 
    write ( *, * ) 'Take ',itermx,' cycles through NEWNOT'
    write ( *, * ) ' '
    write ( *, * ) '  N   Max error   Decay exponent'
!
!  Loop over number of data points.
!
    do n = nlow, nhigh, 2

      if ( n > nlow ) then

        call newnot(break,c,l,4,tnew,l+2,scrtch)

        l = l+2
        t(5+l) = 1.0
        t(6+l) = 1.0
        iter = 1

        do i = 2, l
          t(i+3) = tnew(i)
        end do

        go to 11

      end if

      h = 2.0/real ( n-3)

      do i = 1, 4
        t(i) = -1.0
      end do

      do i = 5, n
        t(i) = real ( i-4)*h-1.0
      end do

      do i = n+1, n+4
        t(i) = 1.0
      end do

      iter = 1
!
!  Construct cubic spline interpolant. then,itermx  times,
!  determine new knots from it and find a new interpolant.
!
   11     continue

      do i = 1, n
        tau(i) = (t(i+1)+t(i+2)+t(i+3))/3.
        gtau(i) = f02(tau(i))
      end do

      call splint(tau,gtau,t,n,4,scrtch,bcoef,iflag)

      call bsplpp(t,bcoef,n,4,scrtch,break,c,l)

      if ( iter <= itermx ) then

        iter = iter+1
        call newnot(break,c,l,4,tnew,l,scrtch)

        do i = 2, l
          t(3+i) = tnew(i)
        end do

        go to 11

      end if
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0
      do i = 1, l
        dx = (break(i+1)-break(i))/step
        do j = 1, istep
          h = real ( j)*dx
          x = break(i) + h
          pnatx = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
          errmax = max(errmax,abs(f02(x)-pnatx))
        end do
      end do
!
!  Calculate the decay exponent.
!
      aloger = alog(errmax)

      if ( n > nlow ) then
        decay = (aloger-algerp)/alog(real ( n)/real ( n-2))
      end if

      algerp = aloger

      write(*,'(i3,e12.4,f11.2)')n,errmax,decay

    end do

  end do

  return
end
subroutine test14
!
!***********************************************************************
!
!! TEST14 tests optimal spline interpolation.
!
!  lenscr = (n-k)(2k+3)+5k+3 is the length of scrtch required in
!  splopt.
!
  integer, parameter :: n = 12
  integer, parameter :: ntitan = 49
  integer, parameter :: k = 5
  integer, parameter :: npk = n+k
  integer, parameter :: lenscr = (n-k)*(2*k+3)+5*k+3
!
  real a(n)
  real bvalue
  real gtitan(ntitan)
  real gtau(ntitan)
  integer i
  integer iflag
  integer ipick(n)
  integer ipicki
  integer lx
  real scrtch(lenscr)
  real t(npk)
  real tau(n)
  real x(ntitan)
!
  data ipick /1,5,11,21,27,29,31,33,35,40,45,49/
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST14'
  write ( *, * ) 'Optimal spline test'

  call titand(x,gtitan,lx)

  do i = 1, n
    ipicki = ipick(i)
    tau(i) = x(ipicki)
    gtau(i) = gtitan(ipicki)
  end do

  call splopt(tau,n,k,scrtch,t,iflag)

  if ( iflag > 1 ) then
    return
  end if

  call splint(tau,gtau,t,n,k,scrtch,a,iflag)

  if ( iflag > 1 ) then
    return
  end if

  do i = 1, lx
    gtau(i) = bvalue(t,a,n,k,x(i),0)
    scrtch(i) = gtitan(i) - gtau(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  I, data point, data, interpolant, error'
  write ( *, * ) ' '
  do i = 1, lx
    write ( *, 620 ) i, x(i), gtitan(i), gtau(i), scrtch(i)
  end do

  620 format(i3,f8.0,f10.4,f9.4,e11.3)

  write ( *, * ) ' '
  write ( *, * ) 'Optimal knots:'
  write ( *, * ) ' '

  write(*,'(i5,g14.6)')(i,t(k+i),i = 1,n-k)

  return
end
subroutine titand(tau,gtau,n)
!
!***********************************************************************
!
!! TITAND represents a property of titanium as a function of
!  temperature.  They have been used extensively as an example in 
!  spline approximation with variable knots.
!
  integer, parameter :: ntau = 49
!
  integer i
  real gtau(ntau)
  real gtitan(ntau)
  integer n 
  real tau(ntau)
!
  data gtitan / &
    0.644, 0.622, 0.638, 0.649, 0.652, 0.639, 0.646, &
    0.657, 0.652, 0.655, 0.644, 0.663, 0.663, 0.668, &
    0.676, 0.676, 0.686, 0.679, 0.678, 0.683, 0.694, &
    0.699, 0.710, 0.730, 0.763, 0.812, 0.907, 1.044, &
    1.336, 1.881, 2.169, 2.075, 1.598, 1.211, 0.916, &
    0.746, 0.672, 0.627, 0.615, 0.607, 0.606, 0.609, &
    0.603, 0.601, 0.603, 0.601, 0.611, 0.601, 0.608 /

  n = ntau

  do i = 1, ntau
    tau(i) = 585.0 + 10.0*i
    gtau(i) = gtitan(i)
  end do
 
  return
end
subroutine test15
!
!***********************************************************************
!
!! TEST15 demonstrates the cubic smoothing spline.
!
!  Values from a cubic B spline are rounded to NDIGIT places
!  after the decimal point, then smoothed via SMOOTH for
!  various values of the control parameter S.
!
  integer, parameter :: npoint = 61
  integer, parameter :: ns = 7
!
  real a(npoint,4)
  real bcoef(7)
  real break(5)
  real coef(4,4)
  real coefsm(4,60)
  real dely
  real dy(61)
  integer i
  integer is
  integer j
  integer l
  integer ndigit
  real s(ns)
  real scrtch(427)
  real sfp
  real t(11)
  real tenton
  real x(npoint)
  real y(npoint)
!
  equivalence(scrtch,coefsm)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST15'
  write ( *, * ) 'Cubic smoothing spline'

  ndigit = 2

  bcoef(1) = 0.0
  bcoef(2) = 0.0
  bcoef(3) = 0.0
  bcoef(4) = 1.0
  bcoef(5) = 0.0
  bcoef(6) = 0.0
  bcoef(7) = 0.0

  t(1) = 0.0
  t(2) = 0.0
  t(3) = 0.0
  t(4) = 0.0
  t(5) = 1.0
  t(6) = 3.0
  t(7) = 4.0
  t(8) = 6.0
  t(9) = 6.0
  t(10) = 6.0
  t(11) = 7.0

  do i = 1, 7
    s(i) = 6.0 * 10.0**(6-i)
  end do

  call bsplpp(t,bcoef,7,4,scrtch,break,coef,l)

  write ( *, * ) ' '
  write ( *, * ) 'The exact data values were rounded to ',ndigit
  write ( *, * ) 'digits after the decimal point.'

  tenton = 10.0**ndigit
  dely = 0.5/tenton

  do i = 1, npoint
    x(i) = 0.1*real ( i-1)
    call ppvalu(break,coef,l,4,x(i),0,y(i))
    y(i) = real ( int(y(i)*tenton+.5))/tenton
    dy(i) = dely
  end do
 
  do i = 1, npoint, 5
    do j = 1, 4
      call ppvalu(break,coef,l,4,x(i),j-1,a(i,j))
    end do
  end do
 
  write ( *, * ) ' '
  write ( *, * ) 'Value and derivatives of noisefree function ' &
    //'at some points:'
  write ( *, * ) ' '
  write(*,'(i4,4e15.7)')(i,(a(i,j),j = 1,4),i=1,npoint,5)
 
  do is = 1, ns

    call smooth(x,y,dy,npoint,s(is),scrtch,a,sfp)
 
    do i = 1, npoint-1
      do j = 1, 4
        coefsm(j,i) = a(i,j)
      end do
    end do
 
    do i = 1, npoint, 5
      do j = 1, 4
        call ppvalu(x,coefsm,npoint-1,4,x(i),j-1,a(i,j))
      end do
    end do
 
    write ( *, * ) ' '
    write ( *, * ) 'Prescribed S =       ',s(is)
    write ( *, * ) 'S(Smoothing spline)  = ',sfp
    write ( *, * ) ' '
    write ( *, * ) 'Value and derivatives of smoothing spline at ' &
      //'corresponding points:'
    write ( *, * ) ' '
    write(*,'(i4,4e15.7)')(i,(a(i,j),j = 1,4),i=1,npoint,5)
  end do
 
  return
end
subroutine test16
!
!***********************************************************************
!
!! TEST16 demonstrates least-squares approximation by splines.
!
!  The program, though ostensibly written for L2-approximation, is 
!  typical for programs constructing a piecewise polynomial 
!  approximation to a function given in some sense.  The subprogram
!  L2APPR, for instance,could easily be replaced by one carrying 
!  out interpolation or some other form of approximation.
!
!  Input is expected in SETDAT, specifying both the data to
!  be approximated and the order and breakpoint sequence of the
!  piecewise polynomial approximating function to be used.  
!  Further, SETDAT is expected to terminate the run, for lack of 
!  further input or because ICOUNT has reached a critical value.
!
!  The number NTIMES is read in in the main program.  It 
!  specifies the number of passes through the knot improvement 
!  algorithm in NEWNOT to be made.  Also, ADDBRK is read in to 
!  specify that, on the average, ADDBRK knots are to be added per
!  pass through NEWNOT.  For example, ADDBRK = .34 would cause a knot 
!  to be added every third pass, as long as NTIMES.LT.50.
!
!  Printed output
!
!  is governed by the three print control integers
!  iprbco  = 1  gives printout of B spline coeffs. of approxim.
!  iprpco  = 1  gives printout of pp repr. of approximation.
!  iprfun  = 1  gives printout of approximation and error at
!                     every data point.
!  the order  k ,the number of pieces  l,and the interior breakpoints
!  are always printed out as are(in l2err) the mean,mean square,and
!  maximum errors in the approximation.
!
!  ICOUNT provides communication with the data-input-and-
!  termination routine SETDAT. it is initialized to  0  to
!  signal to SETDAT when it is being called for the first time. after
!  that,it is up to SETDAT to use ICOUNT for keeping track of the
!  passes through SETDAT.
!
!  Information about the function to be approximated and order and
!  breakpoint sequence of the approximating pp functions is gathered
!  by calling SETDAT.
!
  integer, parameter :: lpkmax = 100
  integer, parameter :: ntmax = 200
  integer, parameter :: ltkmax = 2000
!
  real addbrk
  real bcoef(lpkmax)
  real break
  real coef
  real gtau
  integer i
  integer icount
  integer ii
  integer inot
  integer iprbco
  integer iprfun
  integer iprpco
  integer irun
  integer j
  integer k
  integer l
  integer lbegin
  integer ll
  integer lnew
  integer n
  integer nt
  integer ntau
  integer ntimes
  real q(ltkmax)
  real scrtch(ntmax)
  real t(ntmax)
  real tau
  real totalw
  real weight
!
  common / data / ntau,tau(ntmax),gtau(ntmax),weight(ntmax),totalw
!
!  common /data/ also occurs in setdat,l2appr and l2err. it is 
!  mentioned here only because it might otherwise become undefined 
!  between calls to those subroutines.
!
  common /approx/ break(lpkmax),coef(ltkmax),l,k
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST16'
  write ( *, * ) '  Least squares approximation by splines.'

  do irun = 1, 3

    write ( *, * ) 'Run number ', irun

    icount = 1

    if ( irun == 1 ) then
      ntimes = 1
      addbrk = 0.0
      iprbco = 1
      iprpco = 1
      iprfun = 1
      inot = 0
      call setd1(icount)
    else if ( irun == 2 ) then
      ntimes = 20
      addbrk = 1.0
      iprbco = 0
      iprpco = 0
      iprfun = 0
      inot = 1
      call setd2(icount)
    else
      ntimes = 4
      addbrk = 2.0
      iprbco = 1
      iprpco = 1
      iprfun = 0
      inot = 0
      call setd3(icount)
    end if
!
!  Breakpoints are translated into knots, and the number N of
!  B splines to be used is obtained by calling L2KNTS.
!
    call l2knts ( break, l, k, t, n )
!
!  The integer NTIMES and the real ADDBRK are requested as well as 
!  the print controls  iprbco ,iprpco  and
!  iprfun.  ntimes  passes  are made through the routine new-
!  not, with an increase of  addbrk  knots for every pass.
!
    lbegin = l
    nt = 1
!
!  The B spline coefficients BCOEF of the L2-approximation are 
!  obtained.
!
   10   continue

    call l2appr ( t, n, k, q, scrtch, bcoef )

    if ( iprbco == 1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'B spline coefficients:'
      write ( *, * ) ' '
      write(*,'(5g14.6)') bcoef(1:n)
    end if
!
!  Convert the B spline representation of the approximation to 
!  piecewise polynomial representation.
!
     call bsplpp(t,bcoef,n,k,q,break,coef,l)

     write ( *, * ) ' '
     write ( *, * ) 'Approximation by splines of order ',k
     write ( *, * ) 'using ',l,' intervals.'
     write ( *, * ) ' '
     write ( *, * ) 'Breakpoints:'
     write ( *, * ) ' '
     write(*,'(5g14.6)') (break(ll),ll = 2,l)

     if ( iprpco == 1 ) then

       write ( *, * ) 'Piecewise polynomial representation:'

       do i = 1, l
         ii = (i-1)*k
         write(*,613)break(i),(coef(ii+j),j = 1,k)
       end do

  613      format(f9.3,4e14.6/(9x,4e14.6))

     end if
!
!  Compute and print out various error norms.
!
     call l2err ( iprfun, scrtch, q )
!
!  If NEWNOT has been applied less than NTIMES times,try
!  it again to obtain,from the current approx. a possibly 
!  improved sequence of breakpoints with  addbrk  more breakpoints
!  (on the average) than the current approximation has.
!  If only an increase in breakpoints is wanted, without the
!  adjustment that newnot provides, a fake newnot routine could be
!  used here which merely returns the breakpoints for LNEW
!  equal intervals.
!
    if ( nt < ntimes ) then

      lnew = lbegin+int(nt*addbrk)

      if ( inot == 0 ) then
        call newnot(break,coef,l,k,scrtch,lnew,t)
      else if ( inot == 1 ) then
        call evnnot(break,coef,l,k,scrtch,lnew,t)
      end if

      call l2knts(scrtch,lnew,k,t,n)
      nt = nt+1
      go to 10

    end if

  end do

  return
end
subroutine setd1(icount)
!
!***********************************************************************
!
!! SETD1 provides data for example 2 in chapter xiv. 
!
!  For a general purpose l2-approximation program,it
!  would have to be replaced by a subroutine reading in
!    ntau,tau(i),gtau(i),i = 1,...,ntau
!  and reading in or setting
!    k,l,break(i),i = 1,...,l+1,and weight(i),i=1,...,ntau,
!  as well as  totalw = sum(weight(i) ,i=1,...,ntau).
!
!  ICOUNT is equal to zero when setd is called in L2MAIN
!  for the first time. after that,it is up to setd to use icount
!  for keeping track of the passes through setd. this is important
!  since l2main relies on setd for termination.
!
  integer, parameter :: lpkmax = 100
  integer, parameter :: ntmax = 200
  integer, parameter :: ltkmax = 2000
!
  real break
  real coef
  real gtau
  integer i
  integer icount
  integer k
  integer l
  integer ntau
  real step
  real tau
  real totalw
  real weight
!
  common / data / ntau,tau(ntmax),gtau(ntmax),weight(ntmax),totalw
  common /approx/ break(lpkmax),coef(ltkmax),l,k

  icount = icount+1
  ntau = 10

  do i = 1, ntau-1
    tau(i) = 1.0-0.5**(i-1)
  end do

  tau(ntau) = 1.0

  do i = 1, ntau
    gtau(i) = tau(i)**2+1.0
  end do

  do i = 1, ntau
    weight(i) = 1.0
  end do

  totalw = ntau
  l = 6
  step = 1.0/real ( l)
  k = 2

  do i = 1, l+1
    break(i) = (i-1)*step
  end do

  return
end
subroutine setd2(icount)
!
!***********************************************************************
!
!! SETD2 provides data for example 3 in chapter xiv.
!
  integer, parameter :: lpkmax = 100
  integer, parameter :: ntmax = 200
  integer, parameter :: ltkmax = 2000
!
  real break
  real coef
  real gtau
  integer i
  integer icount
  integer k
  integer l
  integer ntau
  real roun
  real step
  real tau
  real totalw
  real weight
  real x
!
  common / data / ntau,tau(ntmax),gtau(ntmax),weight(ntmax),totalw

  common /approx/ break(lpkmax),coef(ltkmax),l,k

  roun(x) = real ( ifix(x*100.))/100.

  icount = icount+1
  ntau = 65
  step = 3.0/real ( ntau-1)

  do i = 1, ntau
    tau(i) = i*step
    gtau(i) = roun(exp(tau(i)))
    weight(i) = 1.0
  end do

  totalw = ntau
  l = 1
  break(1) = tau(1)
  break(2) = tau(ntau)
  k = 3

  return
end
subroutine setd3(icount)
!
!***********************************************************************
!
!! SETD3 provides data for example 4 in chapter xiv.
!
  integer, parameter :: lpkmax = 100
  integer, parameter :: n = 9
  integer, parameter :: ntmax = 200
  integer, parameter :: ltkmax = 2000
!
  real break
  real brkpic(n)
  real coef
  real gtau
  integer i
  integer icount
  integer k
  integer l
  integer ntau
  real tau
  real totalw
  real weight
!
  common / data / ntau,tau(ntmax),gtau(ntmax),weight(ntmax),totalw
  common /approx/ break(lpkmax),coef(ltkmax),l,k
!
  brkpic(1) = 595.0
  brkpic(2) = 730.985
  brkpic(3) = 794.414
  brkpic(4) = 844.476
  brkpic(5) = 880.06
  brkpic(6) = 907.814
  brkpic(7) = 938.001
  brkpic(8) = 976.752
  brkpic(9) = 1075.0

  icount = icount+1

  call titand(tau,gtau,ntau)

  do i = 1, ntau
    weight(i) = 1.0
  end do

  totalw = ntau
  l = n-1
  k = 5

  do i = 1, n
    break(i) = brkpic(i)
  end do

  return
end
subroutine test17
!
!***********************************************************************
!
!! TEST17 solves a second order boundary value problem.
!
!  Solution of a second order nonlinear two point boundary value
!  problem  on (0.,1.) by collocation with pp functions having 4
!  pieces of order 6.  
!
!  Two passes through NEWNOT are to be made, without any knots
!  being added. 
!
!  Newton iteration is to be stopped when two iterates agree to 6  
!  decimal places.
!
  real addbrk
  real aleft
  real aright
  integer iorder
  integer lbegin
  integer ntimes
  real relerr
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST17'
  write ( *, * ) 'Solution of second order two point'
  write ( *, * ) 'boundary value problem.'

  aleft = 0.0
  aright = 1.0
  lbegin = 4
  iorder = 6
  ntimes = 2
  addbrk = 0.0
  relerr = 1.0e-6

  call colloc(aleft,aright,lbegin,iorder,ntimes,addbrk,relerr)

  return
end
subroutine test18
!
!***********************************************************************
!
!! TEST18 demonstrates taut spline interpolation.
!
  integer, parameter :: n = 12
  integer, parameter :: npoint = 201
!
  real break(122)
  real coef(4,22)
  real gamma
  real gtau(49)
  real gtitan(49)
  integer i
  integer iflag
  integer ipick(n)
  integer ipicki
  integer k
  integer l
  integer lx
  real plotf(npoint)
  real plotts(npoint)
  real plott(npoint)
  real scrtch(119)
  real step
  real tau(n)
  real x(49)
!
  data ipick /1,5,11,21,27,29,31,33,35,40,45,49/
  data k /4/
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST18'
  write ( *, * ) 'Taut spline interpolation'
  write ( *, * ) ' '

  gamma = 2.5

  call titand(x,gtitan,lx)
 
  do i = 1, n
    ipicki = ipick(i)
    tau(i) = x(ipicki)
    gtau(i) = gtitan(ipicki)
  end do
 
  call tautsp(tau,gtau,n,0.,scrtch,break,coef,l,k,iflag)
 
  if ( iflag > 1 ) then
    return
  end if

  step = (tau(n)-tau(1)) / real(npoint-1)
 
  do i = 1, npoint
    plott(i) = tau(1)+step*real ( i-1)
    call ppvalu(break,coef,l,k,plott(i),0,plotf(i))
  end do
 
  if ( gamma < 0.0 ) then
    return
  end if
 
  call tautsp(tau,gtau,n,gamma,scrtch,break,coef,l,k,iflag)
 
  if ( iflag > 1 ) then
    return
  end if

  do i = 1, npoint
    call ppvalu(break,coef,l,k,plott(i),0,plotts(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Cubic spline versus taut spline with gamma  = ', gamma
  write ( *, * ) ' '

  write(*,'(3g14.6)')(plott(i),plotf(i),plotts(i),i = 1,npoint)

  return
end
subroutine test19
!
!***********************************************************************
!
!! TEST19 demonstrates two parametrizations of some data.
!
  integer, parameter :: k = 4
  integer, parameter :: kpkm1 = k+k-1
  integer, parameter :: n = 8
  integer, parameter :: npk = n+k
  integer, parameter :: npiece = 6
  integer, parameter :: npoint = 21
!
  real bcoef(n)
  real break(npiece)
  real ds
  integer i
  integer icount
  integer iflag
  integer l
  real q(n,kpkm1)
  real s(n)
  real scrtch(k,k)
  real ss
  real t(npk)
  real x(n)
  real xcoef(k,npiece)
  real xx(npoint)
  real y(n)
  real ycoef(k,npiece)
  real yy(npoint)
!
  data x /0.,.1,.2,.3,.301,.4,.5,.6/
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST19'
  write ( *, * ) '  Two parameterizations'
!
!  Compute the Y component and set the natural parametrization.
!
  do i = 1, n
    y(i) = (x(i)-0.3)**2
  end do
!
!  Convert X values to knots.  Note that second and second to
!  last X values are not knots.
!
  do icount = 1, 2
 
    write ( *, * ) ' '

    if ( icount == 1 ) then

      write ( *, * ) 'Using the "natural" parameterization.'

      do i = 1, n
        s(i) = x(i)
      end do

    else

      write ( *, * ) 'Using the "uniform" parameterization.'

      do i = 1, n
        s(i) = i
      end do

    end if

    write ( *, * ) ' '
    write ( *, * ) '      X           Y'
    write ( *, * ) ' '

    do i = 1, k
      t(i) = s(1)
    end do
 
    do i = k+1, n
      t(i) = s(i+2-k)
    end do

    do i = n+1, n+k
      t(i) = s(n)
    end do
!
!  Interpolate to X component.
!
    call splint(s,x,t,n,k,q,bcoef,iflag)

    call bsplpp(t,bcoef,n,k,scrtch,break,xcoef,l)
!
!  Interpolate to Y component.  Since data abscissae and knots are
!  the same for both components, we only need to use 
!  backsubstitution
!
    do i = 1, n
      bcoef(i) = y(i)
    end do
 
    call banslv(q,kpkm1,n,k-1,k-1,bcoef)

    call bsplpp(t,bcoef,n,k,scrtch,break,ycoef,l)
!
!  Evaluate curve at some points near the potential trouble spot,
!  the fourth and fifth data points.
!
    ss = s(3)
    ds = (s(6)-s(3)) / real ( npoint-1)
 
    do i = 1, npoint
      call ppvalu(break,xcoef,l,k,ss,0,xx(i))
      call ppvalu(break,ycoef,l,k,ss,0,yy(i))
      ss = ss+ds
    end do
 
    write(*,'(2g14.6)')(xx(i),yy(i),i = 1,npoint)
 
  end do

  return
end
subroutine test20
!
!***********************************************************************
!
!! TEST20 demonstrates bivariate spline interpolation.
!
!  A function G(X,Y) is interpolated on a rectangular mesh of data
!  points (X(I),Y(J)), I = 1 to 7, J=1 to 6.
!
  integer, parameter :: nx = 7
  integer, parameter :: kx = 3
  integer, parameter :: ny = 6
  integer, parameter :: ky = 4
  integer, parameter :: nq = (2*ky-1)*ny
!
  real bcoef(nx,ny)
  real bvalue
  real f03
  real fdata(nx,ny)
  integer i
  integer iflag
  integer j
  integer jj
  integer lefty
  integer mflag
  real q(nq)
  real taux(nx)
  real tauy(ny)
  real tx(nx+kx)
  real ty(ny+ky)
  real work1(nx,ny)
  real work2(nx)
  real work3(nx,ny)
!
  external f03
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST20'
  write ( *, * ) '  Bivariate interpolation'
  write ( *, * ) '  using B-splines.'
!
!  Set up the X data points TAUX and X spline knots TX.
!
  do i = 1, nx
    taux(i) = real ( i)
  end do
 
  do i = 1, kx
    tx(i) = taux(1)
  end do
 
  do i = kx+1, nx
    tx(i) = (taux(i-kx+1)+taux(i-kx+2))/2.0
  end do

  do i = nx+1, nx+kx
    tx(i) = taux(nx)
  end do
!
!  Set up the Y data points TAUY and Y spline knots TY.
!
  do i = 1, ny
    tauy(i) = real ( i)
  end do

  do i = 1, ky
    ty(i) = tauy(1)
  end do

  do i = ky+1, ny
    ty(i) = tauy(i-ky+2)
  end do

  do i = ny+1, ny+ky
    ty(i) = tauy(ny)
  end do
!
!  Generate the function values at the mesh points.
!
  write ( *, * ) ' '
  write ( *, * ) 'The given data to be interpolated:'
  write ( *, * ) ' '
  write(*,'(6f12.1)')(tauy(i),i = 1,ny)
 
  do i = 1, nx
    do j = 1, ny
      fdata(i,j) = f03(taux(i),tauy(j))
    end do
    write(*,'(f5.1,6e12.5)')taux(i),(fdata(i,j),j = 1,ny)
  end do
!
!  Construct the B spline coefficients of the interpolant.
!
  call spli2d(taux,fdata,tx,nx,kx,ny,work2,q,work1,iflag)

  call spli2d(tauy,work1,ty,ny,ky,nx,work2,q,bcoef,iflag)
!
!  Evaluate the interpolation error at the data points.
!
  do j = 1, ny

    call interv(ty,ny+1,tauy(j),lefty,mflag)

    do i = 1, nx

      do jj = 1, ky
        work2(jj) = bvalue(tx,bcoef(1,lefty-ky+jj),nx,kx,taux(i),0)
      end do

      work3(i,j) = bvalue(ty(lefty-ky+1),work2,ky,ky,tauy(j),0)

    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Interpolating function:'
  write ( *, * ) ' '
  write(*,'(6f12.1)')(tauy(j),j = 1,ny)

  do i = 1, nx
    write(*,'(f5.1,6e12.5)')taux(i),(work3(i,j),j = 1,ny)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Interpolation error:'
  write ( *, * ) ' '
  write(*,'(6f12.1)')(tauy(j),j = 1,ny)

  do i = 1, nx
    write(*,'(f5.1,6e12.5)')taux(i),(fdata(i,j)-work3(i,j),j = 1,ny)
  end do
 
  return
end
subroutine test21
!
!***********************************************************************
!
!! TEST21 demonstrates bivariate spline interpolation.
!
!  This is followed by conversion to piecewise polynomial representation 
!  and evaluation.
!
  integer, parameter :: nx = 7
  integer, parameter :: kx = 3
  integer, parameter :: ny = 6
  integer, parameter :: ky = 4
  integer, parameter :: nq = (2*ky-1)*ny
!
  real bcoef(nx,ny)
  real breakx(6)
  real breaky(4)
  real coef(kx,5,ky,3)
  real f03
  real fdata(nx,ny)
  integer i
  integer iflag
  integer j
  integer jj
  integer lefty
  integer lx
  integer ly
  integer mflag
  real q(nq)
  real taux(nx)
  real tauy(ny)
  real tx(nx+kx)
  real ty(ny+ky)
  real value
  real work1(nx,ny)
  real work2(nx)
  real work3(nx,ny)
  real work4(kx,kx,ny)
  real work5(ny,kx,5)
  real work6(ky,ky,21)
!
  external f03
!
!  Note that, with the above parameters, lx = 5, ly=3
!
  equivalence(work4,work6)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST21'
  write ( *, * ) '  Bivariate interpolation'
  write ( *, * ) '  using piecewise polynomials.'
!
!  Set up data points and knots.
!
!  In X, interpolate between knots by parabolic splines,using
!  not-a-knot end condition
!
  do i = 1, nx
    taux(i) = real ( i)
  end do

  do i = 1, kx
    tx(i) = taux(1)
  end do

  do i = kx+1, nx
    tx(i) = (taux(i-kx+1)+taux(i-kx+2)) / 2.0
  end do

  do i = nx+1, nx+kx
    tx(i) = taux(nx)
  end do
!
!  In Y, interpolate at knots by cubic splines, using not-a-knot
!  end condition.
!
  do i = 1, ny
    tauy(i) = i
  end do
 
  do i = 1, ky
    ty(i) = tauy(1)
  end do
 
  do i = ky+1, ny
    ty(i) = tauy(i-ky+2)
  end do

  do i = ny+1, ny+ky
    ty(i) = tauy(ny)
  end do
!
!  Generate and print out function values.
!
  write ( *, * ) ' '
  write ( *, * ) 'The given data to be interpolated:'
  write ( *, * ) ' '
  write(*,'(6f12.1)')(tauy(i),i = 1,ny)
 
  do i = 1, nx
    do j = 1, ny
      fdata(i,j) = f03(taux(i),tauy(j))
    end do
    write(*,'(f5.1,6e12.5)')taux(i),(fdata(i,j),j = 1,ny)
  end do
!
!  Construct the B spline coefficients of the interpolant.
!
  call spli2d(taux,fdata,tx,nx,kx,ny,work2,q,work1,iflag)
 
  call spli2d(tauy,work1,ty,ny,ky,nx,work2,q,bcoef,iflag)
!
!  Convert to piecewise polynomial representation.
!
  call bspp2d(tx,bcoef,nx,kx,ny,work4,breakx,work5,lx)
 
  call bspp2d(ty,work5,ny,ky,lx*kx,work6,breaky,coef,ly)
!
!  Evaluate the interpolation error at mesh points.
!
  do j = 1, ny

    call interv(breaky,ly,tauy(j),lefty,mflag)

    do i = 1, nx
      do jj = 1, ky
        call ppvalu(breakx,coef(1,1,jj,lefty),lx,kx,taux(i),0,work2(jj))
      end do

      call ppvalu(breaky(lefty),work2,1,ky,tauy(j),0,value)

      work3(i,j) = value

    end do

  end do

  write ( *, * ) ' '
  write ( *, * ) 'Interpolating function:'
  write ( *, * ) ' '
  write(*,'(6f12.1)')(tauy(j),j = 1,ny)

  do i = 1, nx
    write(*,'(f5.1,6e12.5)')taux(i),(work3(i,j),j = 1,ny)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Interpolation error:'
  write ( *, * ) ' '
  write(*,'(6f12.1)')(tauy(j),j = 1,ny)

  do i = 1, nx
    write(*,'(f5.1,6e12.5)')taux(i),(fdata(i,j)-work3(i,j),j = 1,ny)
  end do

  return
end
subroutine test22
!
!***********************************************************************
!
!! TEST22 tests cspint, approximate integration of tabulated data via cubic 
!  splines.
!
  integer, parameter :: maxprb = 5
  integer, parameter :: ntab = 13
!
  real a
  real aleft
  real answer(maxprb)
  real b
  real brite
  real e(ntab)
  real ftab(ntab)
  real func
  integer i
  integer ind
  integer iprob
  integer j
  real result
  real work(ntab)
  real xtab(ntab)
  real y(3,ntab)
!
  external func
!
  common /problm/ iprob
!
  a = 0.0
  b = 1.0
  write ( *, * ) ' '
  write ( *, * ) 'TEST22'
  write ( *, * ) '  CSPINT approximates integrals by cubic splines.'
  write ( *, * ) '  Integrals go from ',a,' to ',b
  write ( *, * ) ' '
  write ( *, * ) '     x           x**2         x**3       ', &
    'exp(x)      1/(1+x*x)'
  write ( *, * ) ' '

  do i = 1, maxprb

    iprob = i
    aleft = a-(b-a)/(ntab-3)
    brite = b+(b-a)/(ntab-3)

    do j = 1, ntab
      xtab(j) = (real(j-1)*brite+real(ntab-j)*aleft)/real(ntab-1)
    end do

    do j = 1, ntab
      ftab(j) = func(xtab(j))
    end do

    call cspint(ftab,xtab,ntab,a,b,y,e,work,result,ind)
    answer(i) = result

  end do

  write(*,'(1x,5g14.6)')(answer(i),i = 1,maxprb)

  return
end
function func(x)
!
!***********************************************************************
!
!! FUNC is used as the function to be integrated.
!
  real func
  integer iprob
  real x
!
  common /problm/ iprob
!
  if ( iprob == 1 ) then
    func = x
  else if ( iprob == 2 ) then
    func = x*x
  else if ( iprob == 3 ) then
    func = x*x*x
  else if ( iprob == 4 ) then
    func = exp(x)
  else if ( iprob == 5 ) then
    func = 1.0/(x*x+1.0)
  else
    func = 0.0
  end if

  return
end
function runge(x)
!
!***********************************************************************
!
!! RUNGE evaluates the Runge function.
!
  real runge
  real x
!
  runge = 1.0 / ( 1.0 + 25.0*x*x )

  return
end
function rungep(x)
!
!***********************************************************************
!
!! RUNGEP evaluates the derivative of the Runge function.
!
  real rungep
  real x
!
  rungep = -50.0 * x / ( 1.0 + 25.0*x*x )**2

  return
end
function rungepp(x)
!
!***********************************************************************
!
!! RUNGEPP evaluates the second derivative of the Runge function.
!
  real rungepp
  real u
  real up
  real v
  real vp
  real x
!
  u = -50.0 * x
  up = -50.0
  v = ( 1.0 + 25.0*x*x )**2
  vp = 2.0 * ( 1.0 + 25.0*x*x ) * ( 50.0 * x )

  rungepp = ( up * v - u * vp ) / v**2

  return
end
function f02(x)
!
!***********************************************************************
!
!! F02 evaluates the function SQRT(X+1).
!
  real f02
  real x
!
  f02 = sqrt(x+1.0)

  return
end
function f03(x,y)
!
!***********************************************************************
!
!! F03 evaluates the function max(x-3.5,0.0)**2 + max(y-3.0,0.0)**3
!
  real f03
  real x
  real y
!
  f03 = max(x-3.5,0.0)**2+max(y-3.0,0.0)**3

  return
end
subroutine test23
!
!***********************************************************************
!
!! TEST23 compares the B-spline representation of a cubic 
!  F with its values at knot averages.
!
  real bcoef(23)
  real d0(4)
  real dtip1
  real dtip2
  real f(23)
  real d(4)
  integer i
  integer id
  integer j
  integer jj
  integer n
  real t(27)
  real tave(23)
  real x
!
!  The Taylor coefficients at 0 for the polynomial f are
!
  d0(1) =  -162.0
  d0(2) =    99.0
  d0(3) =   -18.0
  d0(4) =     1.0
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST23'
  write ( *, * ) 'Compare coefficients to values.'
!
!  Set up the knot sequence in the array  t.
!
  n = 13

  do i = 1, 4
    t(i) = 0.0
  end do

  do i = 5, n
    t(i) = real ( i-4)
  end do

  do i = n+1, n+4
    t(i) = 10.0
  end do

  do i = 1, n
!
!  Use nested multiplication to get Taylor coefficients  d  at
!  t(i+2)  from those at  0.
!
    do j = 1, 4
      d(j) = d0(j)
    end do

    do j = 1, 3
      id = 4
      do jj = j, 3
        id = id-1
        d(id) = d(id)+d(id+1)*t(i+2)
      end do
    end do
!
!  Compute B spline coefficients by formula(9).
!
    dtip1 = t(i+2)-t(i+1)
    dtip2 = t(i+3)-t(i+2)
    bcoef(i) = d(1)+(d(2)*(dtip2-dtip1)-d(3)*dtip1*dtip2)/3.0
!
!  Evaluate F at corresponding knot average.
!
    tave(i) = (t(i+1)+t(i+2)+t(i+3))/3.0
    x = tave(i)
    f(i) = d0(1)+x*(d0(2)+x*(d0(3)+x*d0(4)))

  end do

  write ( *, * ) '  I   tave(i)      F at tave(i)      bcoef(i)'
  write ( *, * ) ' '
  write(*,650)(i,tave(i),f(i),bcoef(i),i = 1,n)
  650 format(i3,f10.5,2f16.5)
  return
end
subroutine test24
!
!***********************************************************************
!
!! TEST24 "plots" some B splines.
!
  integer, parameter :: k = 3
  integer, parameter :: n = 7
!
  real dx
  integer i
  integer j
  integer left
  integer mflag
  integer npoint
  real t(n+k)
  real values(n)
  real x
  real xl
!
!  Knot sequence for spline space 
!
  data t /3*0.,2*1.,3.,4.,3*6./
!
!  B spline values are initialized to 0.
!
  do i = 1, n
    values(i) = 0.0
  end do
!
!  NPOINT is the number of evaluation points.
!
  npoint = 31
  write ( *, * ) ' '
  write ( *, * ) 'TEST24'
  write ( *, * ) 'Plot some B splines'
!
!  Set the leftmost evaluation point XL, and the spacing DX.
!
  xl = t(k)
  dx = (t(n+1)-t(k))/real ( npoint-1)
  write(*,600)(i,i = 1,5)
  600 format('X',8x,5('B',i1,'(X)',7x))
!
!  Locate X with respect to the knot array T.
!
  do i = 1, npoint

     x = xl+real ( i-1)*dx
     call interv(t,n,x,left,mflag)
!
!  Get B(i,k)(x) in  values(i),i = 1,...,n.  
!
!  K of these are supplied by BSPLVB, namely:
!  b(left-k+1,k)(x),...,b(left,k)(x).
!
!  All the others are known to be zero.
!
     call bsplvb(t,k,1,x,left,values(left-k+1))

     write(*,610) x,(values(j),j = 3,7)
  610    format(f7.3,5f12.7)
!
!  Zero out the values just computed in preparation for the next 
!  evaluation point.
!
    do j = 1, k
      values(left-k+j) = 0.0
    end do
 
  end do
 
  return
end
subroutine test25
!
!***********************************************************************
!
!! TEST25 plots the polynomials which make up a B spline.
!
  real biatx(4)
  integer i
  integer left
  real t(11)
  real values(4)
  real x
!
!  Set the sequence of knots.
!
  t(1) = 0.0
  t(2) = 0.0
  t(3) = 0.0

  t(4) = 0.0
  t(5) = 1.0
  t(6) = 3.0
  t(7) = 4.0
  t(8) = 6.0

  t(9) = 6.0
  t(10) = 6.0
  t(11) = 6.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST25'
  write ( *, * ) 'Polynomials that build splines.'
!
!  According to the BSPLVB listing, biatx(.) now contains value
!  at  x  of polynomial which agrees on the interval (t(left),
!  t(left+1))  with the B spline  b(left-4+. ,4,t). hence,
!  biatx(8-left)  now contains value of that polynomial for
!  b(left-4 +(8-left) ,4,t) = b(4,4,t). since this B spline
!  has support (t(4),t(8)),it makes sense to run  left = 4,
! ...,7,storing the resulting values in  values(1),...,
!  values(4)  for later printing.
!
  do i = 1, 40

    x = i*0.2-1.0

    do left = 4, 7
      call bsplvb(t,4,1,x,left,biatx)
      values(left-3) = biatx(8-left)
    end do

    write(*,620) x,values

  end do

  620 format(f10.1,4f20.8)
  return
end
