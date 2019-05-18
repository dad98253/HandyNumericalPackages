!  nspprb.f90  22 June 2000
!
program nspprb
!
!*******************************************************************************
!
!! NSPPRB runs the NSPCG tests.
!
  write ( *, * ) ' '
  write ( *, * ) 'NSPPRB'
  write ( *, * ) '  Tests for NSPCG.'

  call test01
  call test02
  call test03
  call test04
  call test05

  write ( *, * ) ' '
  write ( *, * ) 'NSPPRB'
  write ( *, * ) '  Normal end of NSPCG tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 ...
!
  integer inw
  integer mdim
  integer n
  integer ndim
  integer nw
  integer nx
  integer ny
!
  parameter ( inw = 300 )
  parameter ( mdim = 4 )
  parameter ( nw = 600)
  parameter ( nx = 10 )
  parameter ( ny = 10 )

  parameter ( n = nx * ny )
  parameter ( ndim = n )
!
  real coef(ndim,mdim)
  real h
  integer ier
  integer ip(1)
  integer iparm(30)
  integer isym
  integer iwksp(inw)
  integer jcoef(5)
  integer k
  integer khi
  integer klo
  integer maxnz
  integer p(1)
  real rhs(n)
  real rparm(30)
  real u(n)
  real ubar(1)
  real wksp(nw)
!
  external cg
  external mic2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'

  maxnz = 3
  h = 1.0 / 11.0
!
!  Generate Laplace's equation.
!
  isym = 0
  call matgen ( n, nx, ny, ndim, maxnz, jcoef, coef, rhs, h, isym )
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(2) = 50
  iparm(3) = 3

  rparm(1) = 1.0e-4
!
!  Set an initial value for the solution.
!
  u = 0.0
!
!  Solve the system.
!
  call nspcg ( mic2, cg, ndim, mdim, n, maxnz, coef, jcoef, p, ip, u, ubar, &
    rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, * ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine matgen ( n, nx, ny, ndim, maxnz, jcoef, coef, rhs, h, isym )
!
!*******************************************************************************
!
!! MATGEN generates the data for a discretization of a PDE.
!
!
!  Discussion:
!
!    MATGEN generates the itpack matrices jcoef, coef, and rhs
!    for a 5-point central difference approximation of the
!    self-adjoint pde
!
!      (a*uxx )   +(c*uyy )  +f*u = g
!
!    where the coefficients a, c, f, and g are functions of (x,y)
!    (supplied by the user) and the domain is the unit square
!    (0,1) x (0,1).  dirichlet boundary conditions are imposed
!    upon the boundary in the form of the function ub(x,y) also
!    supplied by the user.
!
!  Parameters:
!
!        n       number of linear equations (output)
!        ndim    row dimension of jcoef and coef in the calling routine
!                 (input)
!        maxnz   maximum number of nonzeros per row (output)
!        jcoef   array of column indices (output)
!        coef    array of equation coefficients (output)
!        rhs     vector of right-hand-side values (output)
!        h       mesh spacing (input)
!
  integer n
  integer ndim
!
  real a
  real ae
  real aw
  real c
  real cc
  real cn
  real coef(ndim,5)
  real cs
  real f
  real fp
  real g
  real gp
  real h
  integer i
  integer isym
  integer j
  integer jcoef(5)
  integer maxnz
  integer neq
  integer nx
  integer ny
  real rhs(n)
  logical symm
  real ub
  real x
  real xx
  real y
  real yy
!
!     statement functions
!
  a(x,y) = 1.0
  c(x,y) = 2.0
  f(x,y) = 0.0
  g(x,y) = 0.0
  ub(x,y) = 1.0e0+x*y
!
  symm = ( isym == 0 )
  maxnz = 5

  if ( symm ) then
    maxnz = 3
  else
    maxnz = 5
  end if
!
!  Loop on equations, from left to right and from down to up.
!
  neq = 0

  do j = 1, ny

    yy = real ( j ) / real ( ny + 1 )

    do i = 1, nx

      xx = real ( i ) / real ( nx + 1 )

      neq = neq + 1
      ae = a(xx + 0.5 * h, yy )
      aw = a(xx - 0.5 * h, yy )
      cn = c(xx,yy+0.5 * h)
      cs = c(xx,yy-0.5 * h)
      fp = f(xx,yy)
      gp = g(xx,yy)

      cc = ae + cn + aw + cs - fp * h**2
!
!  Center point
!
      coef(neq,1) = cc
      rhs(neq) = - h**2 * gp
!
!  East point
!
      if ( i /= nx ) then
        coef(neq,2) = - ae
      else
        coef(neq,2) = 0.0
        rhs(neq) = rhs(neq)+ae*ub( 1.0, yy)
      end if
!
!  North point
!
      if ( j /= ny ) then
        coef(neq,3) = - cn
      else
        coef(neq,3) = 0.0
        rhs(neq) = rhs(neq)+cn*ub( xx, 1.0 )
      end if
!
!  West point
!
      if ( i /= 1 ) then
        if(.not. symm) coef(neq,4) = -aw
      else
        if(.not. symm) coef(neq,4) = 0.0
        rhs(neq) = rhs(neq)+aw*ub(0.0,yy)
      end if
!
!  South point
!
      if ( j /= 1 ) then
        if(.not. symm) coef(neq,5) = -cs
      else
        if ( .not. symm) coef(neq,5) = 0.0
        rhs(neq) = rhs(neq)+cs*ub(xx,0.0)
      end if

    end do
  end do
!
!  Define data structure 2.
!
  jcoef(1) = 0
  jcoef(2) = 1
  jcoef(3) = nx

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 ...
!
  integer mdim
  integer ndim
  integer n
  integer nw
  integer nx
  integer ny
  integer inw

  parameter ( mdim = 3 )
  parameter ( nx = 10 )
  parameter ( ny = 10 )

  parameter ( n = nx * ny )

  parameter ( ndim = n )
  parameter ( nw = 10 * n )
  parameter ( inw = 3 * n )
!
  real coef(ndim,mdim)
  integer i
  integer ier
  integer ip(1)
  integer iparm(30)
  integer iwksp(inw)
  integer j
  integer jcoef(3)
  integer k
  integer khi
  integer klo
  integer maxnz
  integer p(1)
  real rhs(n)
  real rparm(30)
  real u(n)
  real ubar(1)
  real wksp(nw)
  real x
  real y
!
  external bic2
  external cg
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'

  maxnz = 3
!
!  Generate COEF.
!
  do i = 1, n
    coef(i,1) = 6.0
    coef(i,2) = -1.0
    coef(i,3) = -2.0
  end do
!
!  Set the right hand side.
!
  rhs = 0.0
!
!  Make modifications for boundary values.
!
  k = 0

  do j = 1, ny

    y = real ( j ) / real ( ny + 1 )

    do i = 1, nx

      x = real ( i ) / real ( nx + 1 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k) + 2.0
      else if ( j == ny ) then
        rhs(k) = rhs(k) + 2.0 * ( 1.0 + x )
        coef(k,3) = 0.0
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k) + 1.0
      else if ( i == nx ) then
        rhs(k) = rhs(k) + 1.0 + y
        coef(k,2) = 0.0
      end if

    end do
  end do

  jcoef(1) = 0
  jcoef(2) = 1
  jcoef(3) = nx
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  iparm(4) = 6
  iparm(18) = 1
  iparm(19) = nx

  rparm(1) = 1.0e-4
!
!  Set an initial value for the solution.
!
  u = 0.0
!
!  Solve the system.
!
  call nspcg ( bic2, cg, ndim, mdim, n, maxnz, coef, jcoef, p, ip, u, ubar, &
    rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, * ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 ...
!
  integer mdim
  integer inw
  integer n
  integer ndim
  integer nw
  integer nx
  integer ny

  parameter ( mdim = 5 )
  parameter ( inw = 500 )
  parameter ( nw = 800 )
  parameter ( nx = 10 )
  parameter ( ny = 10 )

  parameter ( n = nx * ny )
  parameter ( ndim = n )
!
  real coef(ndim,mdim)
  integer i
  integer ier
  integer ip(n)
  integer iparm(30)
  integer iwksp(inw)
  integer j
  integer jcoef(ndim,5)
  integer k
  integer khi
  integer klo
  integer maxnz
  integer p(n)
  real rhs(n)
  real rparm(30)
  real u(n)
  real ubar(1)
  real wksp(nw)
  real x
  real y
!
  external cg
  external rs6
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
!
!  generate coef, jcoef, and rhs.
!
  maxnz = 5

  do i = 1, n

    coef(i,1) = 6.0
    coef(i,2) = -1.0
    coef(i,3) = -2.0
    coef(i,4) = -1.0
    coef(i,5) = -2.0

    jcoef(i,1) = i
    jcoef(i,2) = i+1
    jcoef(i,3) = i+nx
    jcoef(i,4) = i-1
    jcoef(i,5) = i-nx

  end do

  rhs = 0.0

  k = 0

  do j = 1, ny

    y = real ( j ) / real ( ny + 1 )

    do i = 1, nx

      x = real ( i ) / real ( nx + 1 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k)+2.0
        coef(k,5) = 0.0
        jcoef(k,5) = 0
      else if ( j == ny ) then
        rhs(k) = rhs(k)+2.0*(1.0+x)
        coef(k,3) = 0.0
        jcoef(k,3) = 0
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k)+1.0
        coef(k,4) = 0.0
        jcoef(k,4) = 0
      else if ( i == nx ) then
        rhs(k) = rhs(k)+1.0+y
        coef(k,2) = 0.0
        jcoef(k,2) = 0
      end if

    end do
  end do
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  iparm(12) = 1
  iparm(14) = 1
  iparm(23) = 0

  rparm(1) = 1.0e-4
!
!  Set an initial value for the solution.
!
  u = 0.0
!
!  Call REDBLK to check for property A.
!
  call redblk ( ndim, n, maxnz, coef, jcoef, p, ip, 1, iwksp, ier )
!
!  Solve the system.
!
  call nspcg ( rs6, cg, ndim, mdim, n, maxnz, coef, jcoef, p, ip, &
    u, ubar, rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, * ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 ...
!
  integer inw
  integer mdim
  integer n
  integer ndim
  integer nw
  integer nx
  integer ny
  integer nz

  parameter ( inw = 400 )
  parameter ( mdim = 5 )
  parameter ( nw = 1000 )
  parameter ( nx = 10 )
  parameter ( ny = 10 )
  parameter ( nz = 1 )

  parameter ( n = nx * ny * nz )
  parameter ( ndim = n )
!
  real coef(ndim,mdim)
  integer i
  integer ier
  integer ip(n)
  integer iparm(30)
  integer iwksp(inw)
  integer j
  integer jcoef(5)
  integer k
  integer khi
  integer klo
  integer maxnz
  integer nxp
  integer nyp
  integer nzp
  integer p(n)
  integer patt(2)
  real rhs(n)
  real rparm(30)
  real u(n)
  real ubar(1)
  real wksp(nw)
  real x
  real y
!
  external sor
  external sor7
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
!
!  generate coef, jcoef, and rhs.
!
  maxnz = 3

  do i = 1, n
    coef(i,1) = 6.0
    coef(i,2) = -1.0
    coef(i,3) = -2.0
  end do

  rhs = 0.0

  k = 0
  do j = 1, ny

    y = real ( j ) / real ( ny + 1 )

    do i = 1, nx

      x = real ( i ) / real ( nx + 1 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k) + 2.0
      else if ( j == ny ) then
        rhs(k) = rhs(k) + 2.0 * ( 1.0 + x ) 
        coef(k,3) = 0.0
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k) + 1.0
      else if ( i == nx ) then
        rhs(k) = rhs(k) + 1.0 + y
        coef(k,2) = 0.0
      end if

    end do
  end do

  jcoef(1) = 0
  jcoef(2) = 1
  jcoef(3) = nx
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  iparm(14) = 1
  iparm(19) = nx

  rparm(1) = 1.0e-4
!
!  Set an initial value for the solution.
!
  u = 0.0
!
!  Define a color pattern.
!
  nxp = 1
  nyp = 2
  nzp = 1

  patt(1) = 1
  patt(2) = 2

  call color ( nxp, nyp, nzp, nx, ny, nz, patt, p )
!
!  Solve the system.
!
  call nspcg ( sor7, sor, ndim, mdim, n, maxnz, coef, jcoef, p, ip, &
    u, ubar, rhs, wksp, iwksp, nw, inw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, * ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests SORW.
!
  integer n
  integer nw
  integer nx
  integer ny

  parameter ( nw = 800 )
  parameter ( nx = 10 )
  parameter ( ny = 10 )

  parameter ( n = nx * ny )
!
  real coef(1)
  integer i
  integer ier
  integer iparm(30)
  integer j
  integer jcoef(1)
  integer jwfac(1)
  integer k
  integer khi
  integer klo
  real rhs(n)
  real rparm(30)
  real u(n)
  real ubar(1)
  real wfac(1)
  real wksp(nw)
  real x
  real y
!
  external mult
  external sorpass
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  SORW implements the SOR iteration.'
!
!  Initialize the right hand side to zero, but adjust it for boundary values.
!
  rhs = 0.0

  k = 0
  do j = 1, ny

    y = real ( j ) / real ( ny + 1 )

    do i = 1, nx

      x = real ( i ) / real ( nx + 1 )

      k = k + 1

      if ( j == 1 ) then
        rhs(k) = rhs(k) + 2.0
      else if ( j == ny ) then
        rhs(k) = rhs(k) + 2.0 * ( 1.0 + x )
      end if

      if ( i == 1 ) then
        rhs(k) = rhs(k) + 1.0
      else if ( i == nx ) then
        rhs(k) = rhs(k) + 1.0 + y
      end if

    end do
  end do
!
!  Initialize the parameters to their default values.
!
  call dfault ( iparm, rparm )
!
!  Reset some parameters to more appropriate values.
!
  iparm(3) = 3
  rparm(1) = 1.0e-4
!
!  Set an initial value for the solution.
!
  u = 0.0
!
!  Solve the system.
!
  call sorw ( mult, sorpass, coef, jcoef, wfac, jwfac, n, u, ubar, rhs, &
    wksp, nw, iparm, rparm, ier )
!
!  Print the solution.
!
  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, * ) ' '

  do klo = 1, n, 10
    khi = min ( klo + 9, n )
    write ( *, '(10f8.3)' ) ( u(k), k = klo, khi )
  end do

  return
end
subroutine mult ( coef, jcoef, wfac, jwfac, n, x, y )
!
!*******************************************************************************
!
!! MULT ...
!
  integer n
!
  real coef
  integer i
  integer jcoef
  integer jwfac
  integer nx
  real wfac
  real x(n)
  real y(n)
!
  nx = 10
!
!  Compute the product as if first superdiagonal and first
!  subdiagonal were full.
!
  y(1) = 6.0 * x(1) - x(2) - 2.0 * x(nx+1)
  do i = 2, nx
    y(i) = 6.0 * x(i) - x(i+1) - x(i-1) - 2.0 * x(i+nx)
  end do

  do i = nx+1, n-nx
    y(i) = 6.0 * x(i) - x(i+1) - x(i-1) - 2.0 * ( x(i+nx) + x(i-nx) )
  end do

  do i = n-nx+1, n-1
    y(i) = 6.0 * x(i) - x(i+1) - x(i-1) - 2.0 * x(i-nx)
  end do
  y(n) = 6.0 * x(n) - x(n-1) - 2.0 * x(n-nx)
!
!  Make corrections to y vector for zeros in first superdiagonal
!  and first subdiagonal.
!
  do i = nx, n-nx,nx
    y(i) = y(i) + x(i+1)
  end do

  do i = nx+1, n-nx+1, nx
    y(i) = y(i) + x(i-1)
  end do

  return
end
subroutine sorpass ( coef, jcoef, wfac, jwfac, n, u, rhs, unew )
!
!*******************************************************************************
!
!! SORPASS does an SOR iteration.
!
!        unew=inv((1/w)*d+l)*(((1/w)-1)*d*un+rhs-u*un)
!
!  parameters --
!
!        n      order of system
!        u      current solution vector
!        rhs    right hand side
!        unew   updated solution vector
!
!  specifications for parameters
!
  real alphab
  real betab
  real coef(1)
  real con
  real fff
  integer i
  integer ibgn
  integer iend
  integer j
  integer jcoef(1)
  integer jwfac(1)
  integer n
  integer nx
  real omega
  logical omgadp
  real rhs(1)
  real specr
  real u(1)
  real unew(1)
  real wfac(1)
!
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!  temp=((1/w)-1)*d*un+rhs-u*un
!     unew is used for temp.
!
  nx = 10
  con = 6.0 * ( 1.0 / omega - 1.0 )

  do i = 1, n
    unew(i) = con * u(i)+rhs(i)
  end do

  do i = 1, n-1
    unew(i) = unew(i) + u(i+1)
  end do

  do i = 1, n-nx
    unew(i) = unew(i) + 2.0 * u(i+nx)
  end do

  do i = nx,n-nx, nx
    unew(i) = unew(i) - u(i+1)
  end do
!
!  unew=inv((1/w)*d+l)*temp
!
  con = omega / 6.0

  do j = 1, nx

    ibgn = ( j - 1 ) * nx + 1
    iend = j * nx
    unew(ibgn) = con * unew(ibgn)

    do i = ibgn+1, iend
      unew(i) = con * ( unew(i) + unew(i-1) )
    end do

    if ( j /= nx ) then
      do i = ibgn, iend
        unew(i+nx) = unew(i+nx) + 2.0 * unew(i)
      end do
    end if

  end do

  return
end
