function c_le_l2 ( x, y )
!
!*******************************************************************************
!
!! C_LE_L2 := X <= Y for complex values, and the L2 norm.
!
!
!  Definition:
!
!    The L2 norm can be defined here as:
!
!      C_NORM2(X) = sqrt ( ( real (X) )**2 + ( imag (X) )**2 )
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex X, Y, the values to be compared.
!
!    Output, logical C_LE_L2, is TRUE if X <= Y.
!
  logical c_le_l2
  complex x
  complex y
!
  if ( ( real ( x ) )**2 + ( imag ( x ) )**2 <= &
       ( real ( y ) )**2 + ( imag ( y ) )**2 ) then
    c_le_l2 = .true.
  else
    c_le_l2 = .false.
  end if

  return
end
subroutine c_swap ( x, y )
!
!*******************************************************************************
!
!! C_SWAP swaps two complex values.
!
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, complex X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  complex x
  complex y
  complex z
!
  z = x
  x = y
  y = z

  return
end
subroutine c3_check ( n, ierror )
!
!*******************************************************************************
!
!! C3_CHECK checks the dimensions of a complex tridiagonal matrix.
!
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Output, integer IERROR, error flag.
!    0, no errors detected.
!    1, N was less than 2.
!
  integer ierror
  integer n
!
  ierror = 0

  if ( n < 2 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'C3_CHECK - Fatal error!'
    write ( *, * ) '  N must be at least 2.'
    write ( *, * ) '  The input N was ', n
  end if

  return
end
subroutine c3_cr_fa ( n, subd, diag, supd )
!
!*******************************************************************************
!
!! C3_CR_FA decomposes a complex tridiagonal matrix via cyclic reduction.
!
!
!  Discussion:
!
!    If C3_CR_FA has decomposed a matrix A, then C3_CR_SL may be used to
!    solve linear systems A * x = b.
!
!    On a vector computer, cyclic reduction can be very much faster than
!    standard Gauss elimination techniques, such as SGTSL from LINPACK,
!    which do not vectorize well.
!
!    C3_CR_FA and C3_CR_SL will be slower than the Cray SCILIB routine
!    TRID.  On the other hand, TRID does not provide a factorization, and
!    source code for TRID is not generally available.
!
!    C3_CR_FA does not employ pivoting.  Hence, the results can be more
!    sensitive to ill-conditioning than standard Gauss elimination.  In
!    particular, C3_CR_FA will fail if any diagonal element of the matrix
!    is zero.  Other matrices may also cause C3_CR_FA to fail.
!
!    C3_CR_FA can be guaranteed to work properly if the matrix is strictly
!    diagonally dominant, that is, if the absolute value of the diagonal
!    element is strictly greater than the sum of the absolute values of
!    the offdiagonal elements, for each equation.
!
!    The algorithm may be illustrated by the following figures:
!
!    The initial matrix is given by:
!
!          D1 U1
!          L1 D2 U2
!             L2 D3 U3
!                L3 D4 U4
!                   L4 D5 U5
!                      L5 D6
!
!    Rows and columns are permuted in an odd/even way to yield:
!
!          D1       U1
!             D3    L2 U3
!                D5    L4 U5
!          L1 U2    D2
!             L3 U4    D4
!                L5       D6
!
!    A block LU decomposition is performed to yield:
!
!          D1      | U1
!             D3   | L2 U3
!                D5|    L4 U5
!          --------+---------
!                  | D2'F3
!                  | F1 D4'F4
!                  |    F2 D6'
!
!    Here, D2', D4' and D6' are marked with primes to note that they are
!    altered by the decomposition process.
!
!    For large systems, this reduction is repeated on the lower right hand
!    tridiagonal subsystem until a completely upper triangular system
!    is obtained.  The system is then factored into the product of a
!    lower triangular system and an upper triangular one, and the information
!    defining this factorization may be used by C3_CR_SL to solve linear systems.
!
!  Modified:
!
!    06 December 1998
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, complex SUBD(0:2*N), DIAG(2*N), SUPD(0:2*N).
!
!    On input, SUBD, DIAG and SUPD contain the subdiagonal, diagonal,
!    and superdiagonal entries of the matrix.
!
!    The diagonal entries are stored in DIAG(1) through DIAG(N).
!    The subdiagonal entries are stored in SUBD(1) through SUBD(N-1).
!    The superdiagonal entries are stored in SUPD(1) through SUPD(N-1).
!    The extra entries in DIAG, SUBD, and SUPD need not be initialized
!    by the user.
!
!    On output, SUBD, DIAG, and SUPD contain information defining the
!    factorization of the tridiagonal matrix.  This information
!    will be needed by C3_CR_SL to solve linear systems.
!
!    The extra positions in the arrays, DIAG(N+1) through DIAG(2*N),
!    SUBD(0), SUBD(N) through SUBD(2*N), SUPD(0) and SUPD(N) through
!    SUPD(2*N), are used for workspace and storage.
!
!    In particular, SUBD and SUPD must be declared with an initial,
!    0-th element, or the algorithm will not work.
!
  integer n
!
  complex diag(2*n)
  integer i
  integer iful
  integer ifulp
  integer ihaf
  integer il
  integer ilp
  integer inc
  integer incr
  integer ipnt
  integer ipntp
  complex subd(0:2*n)
  complex supd(0:2*n)
!
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_CR_FA - Fatal error!'
    write ( *, * ) '  Nonpositive N = ', n
    return
  end if

  if ( n == 1 ) then
    diag(1) = 1.0E+00 / diag(1)
    return
  end if
!
!  Zero out the workspace entries.
!
  subd(0) = 0.0E+00
  supd(0) = 0.0E+00
  subd(n) = 0.0E+00
  supd(n) = 0.0E+00

  diag(n+1:2*n) = 0.0E+00
  subd(n+1:2*n) = 0.0E+00
  supd(n+1:2*n) = 0.0E+00

  il = n
  subd(il) = 0.0E+00
  supd(il) = 0.0E+00
  ipntp = 0

  do while ( il > 1 )

    ipnt = ipntp
    ipntp = ipntp + il

    if ( mod(il,2) == 1 ) then
      inc = il + 1
    else
      inc = il
    end if

    incr = inc / 2
    il = il / 2
    ihaf = ipntp + incr + 1
    ifulp = ipnt + inc + 2

!dir$ ivdep
    do ilp = incr, 1, -1
      ifulp = ifulp - 2
      iful = ifulp - 1
      ihaf = ihaf - 1
      diag(iful) = 1.0E+00 / diag(iful)
      subd(iful) = subd(iful) * diag(iful)
      supd(ifulp) = supd(ifulp) * diag(ifulp+1)
      diag(ihaf) = diag(ifulp) - supd(iful) * subd(iful) &
        - supd(ifulp) * subd(ifulp)
      subd(ihaf) = - subd(ifulp) * subd(ifulp+1)
      supd(ihaf) = - supd(ifulp) * supd(ifulp+1)
    end do

  end do

  diag(ipntp+1) = 1.0E+00 / diag(ipntp+1)

  return
end
subroutine c3_cr_sl ( n, subd, diag, supd, rhs )
!
!*******************************************************************************
!
!! C3_CR_SL solves a complex linear system factored by C3_CR_FA.
!
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  C3_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.
!
!    If the factors computed by C3_CR_FA are passed to C3_CR_SL, then
!    one or many linear systems involving the matrix A may be solved.
!
!
!    The form of the equations may be summarized as:
!
!                         DIAG(1) * X(1) + SUPD(1) * X(2)   = RHS(1)
!
!    SUBD(I-1) * X(I-1) + DIAG(I) * X(I) + SUPD(I) * X(I+1) = RHS(I)
!
!    SUBD(N-1) * X(N-1) + DIAG(N) * X(N)                    = RHS(N)
!
!    with the middle form used for equations I = 2 through N-1.
!
!    Here SUBD, DIAG and SUPD are the lower diagonal, diagonal, and upper
!    diagonal coefficients of the tridiagonal system.
!
!
!    Note that C3_CR_FA does not perform pivoting, and so the solution
!    produced by C3_CR_SL may be less accurate than a solution produced
!    by a standard Gauss algorithm.  However, such problems can be guaranteed
!    not to occur if the matrix A is strictly diagonally dominant, that is,
!    if the absolute value of the diagonal coefficient is greater than the
!    sum of the absolute values of the two off diagonal coefficients, for each
!    row of the matrix.
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, complex SUBD(0:2*N), DIAG(2*N), SUPD(0:2*N).  These arrays
!    contain information describing the factorization of the original
!    tridiagonal matrix, as computed by C3_CR_FA.
!
!    Input/output, complex RHS(0:2*N).
!
!    On input, RHS contains the right hand side vector in locations
!    1 thorugh N.  The zero-th entry, and the second N locations are
!    used for workspace.
!
!    On output, the locations 1 through N of RHS contain the solution
!    of the linear system.
!
  integer n
!
  complex diag(2*n)
  integer i
  integer iful
  integer ifulm
  integer ifulp
  integer ihaf
  integer il
  integer ipnt
  integer ipntp
  complex subd(0:2*n)
  integer ndiv
  complex supd(0:2*n)
  complex rhs(0:2*n)
!
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_CR_SL - Fatal error!'
    write ( *, * ) '  Nonpositive N = ', n
    return
  end if

  if ( n == 1 ) then
    rhs(1) = diag(1) * rhs(1)
    return
  end if
!
!  Zero out workspace entries of RHS.
!
  rhs(0) = 0.0E+00
  rhs(n+1:2*n) = 0.0E+00

  subd(0) = 0.0E+00
  supd(0) = 0.0E+00
  il = n
  ndiv = 1
  ipntp = 0

  do while ( il > 1 )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

!dir$ ivdep
    do iful = ipnt+2, ipntp, 2
      ihaf = ihaf + 1
      ifulp = iful + 1
      ifulm = iful - 1
      rhs(ihaf) = rhs(iful) - subd(ifulm) * rhs(ifulm) - supd(iful) * rhs(ifulp)
    end do

  end do

  rhs(ihaf) = rhs(ihaf) * diag(ihaf)
  ipnt = ipntp

  do while ( ipnt > 0 )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

!dir$ ivdep
    do ifulm = ipnt+1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful) = rhs(ihaf)
      rhs(ifulm) = diag(ifulm) * ( rhs(ifulm) - subd(ifulm-1) * rhs(ifulm-1) &
        - supd(ifulm) * rhs(iful) )
    end do

  end do

  return
end
subroutine c3_jac_sl ( n, a1, a2, a3, b, x, maxit, job )
!
!*******************************************************************************
!
!! C3_JAC_SL tries to solve a complex tridiagonal system using Jacobi iteration.
!
!
!  Discussion:
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1), the subdiagonal, diagonal,
!    and superdiagonal of the matrix.
!
!    Input, complex B(N), the right hand side of the linear system.
!
!    Input/output, complex X(N), an approximate solution to the system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  complex b(n)
  integer i
  integer ierror
  integer job
  integer maxit
  integer numit
  complex x(n)
  complex xnew(n)
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
!  Check the dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_JAC_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a2(i) == ZERO ) then
      write ( *, * ) ' '
      write ( *, * ) 'C3_JAC_SL - Fatal error!'
      write ( *, * ) '  Zero diagonal entry, index = ', i
      return
    end if
  end do

  if ( job == 0 ) then

    do numit = 1, maxit

      xnew(1) =   b(1)                  - a3(1) * x(2)
      do i = 2, n - 1
        xnew(i) = b(i) - a1(i) * x(i-1) - a3(i) * x(i+1)
      end do
      xnew(n) =   b(n) - a1(n) * x(n-1)

      xnew(1:n) = xnew(1:n) / a2(1:n)

      x(1:n) = xnew(1:n)

    end do

  else

    do numit = 1, maxit

      xnew(i) =   b(1)                    - a1(2) * x(2)
      do i = 2, n - 1
        xnew(i) = b(i) - a3(i-1) * x(i-1) - a1(i+1) * x(i+1)
      end do
      xnew(n) =   b(n) - a3(n-1) * x(n-1)

      xnew(1:n) = xnew(1:n) / a2(1:n)

      x(1:n) = xnew(1:n)

    end do

  end if

  return
end
subroutine c3_mxv ( n, a1, a2, a3, x, b )
!
!*******************************************************************************
!
!! C3_MXV multiplies a complex tridiagonal matrix times a vector.
!
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1),
!    the nonzero diagonals of the linear system.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A * x.
!
  integer n
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  complex b(n)
  integer i
  integer ierror
  complex x(n)
!
!  Check dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions for the matrix.'
    return
  end if

  b(1:n) = a2(1:n) * x(1:n)

  b(2:n) = b(2:n) + a1(2:n) * x(1:n-1)

  b(1:n-1) = b(1:n-1) + a3(1:n-1) * x(2:n)

  return
end
subroutine c3_np_det ( n, a2, det )
!
!*******************************************************************************
!
!! C3_NP_DET returns the determinant of a complex tridiagonal system factored by C3_NP_FA.
!
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex A2(N), the factor information computed by C3_NP_FA.
!
!    Output, complex DET, the determinant of the matrix.
!
  integer n
!
  complex a2(n)
  complex det
  integer i
  integer ierror
!
!  Check the dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_NP_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  det = product ( a2(1:n) )

  return
end
subroutine c3_np_fa ( n, a1, a2, a3, info )
!
!*******************************************************************************
!
!! C3_NP_FA factors a complex tridiagonal system without pivoting.
!
!
!  Discussion:
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output complex A1(2:N), A2(1:N), A3(1:N-1), the subdiagonal,
!    diagonal, and superdiagonal of the matrix.  On output, these are
!    overwritten by factorization information.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer n
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  integer i
  integer ierror
  integer info
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
!  Check the dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_NP_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do i = 1, n-1

    if ( a2(i) == ZERO ) then
      info = i
      write ( *, * ) ' '
      write ( *, * ) 'C3_NP_FA - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if

    a1(i+1) = a1(i+1) / a2(i)
    a2(i+1) = a2(i+1) - a1(i+1) * a3(i)

  end do

  if ( a2(n) == ZERO ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'C3_NP_FA - Fatal error!'
    write ( *, * ) '  Zero pivot on step ', info
    return
  end if

  return
end
subroutine c3_np_ml ( n, a1, a2, a3, x, b, job )
!
!*******************************************************************************
!
!! C3_NP_ML computes A * x or x * A, where A has been factored by C3_NP_FA.
!
!
!  Modified:
!
!    24 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1), the LU factors from C3_FA.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product.
!
!    Input, integer JOB, specifies the product to find.
!    0, compute A * x.
!    nonzero, compute transpose ( A ) * x.
!
  integer n
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  complex b(n)
  integer i
  integer ierror
  integer job
  complex x(n)
!
!  Check the dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_NP_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Compute X := U * X
!  Compute X: = L * X.
!
  if ( job == 0 ) then

    b(1:n-1) = a2(1:n-1) * x(1:n-1) + a3(1:n-1) * x(2:n)
    b(n) = a2(n) * x(n)

    b(2:n) = b(2:n) + a1(2:n) * b(1:n-1)
!
!  Compute X: = transpose ( L ) * X.
!  Compute X: = transpose ( U ) * X.
!
  else

    b(1:n-1) = x(1:n-1) + a1(2:n) * x(2:n)
    b(n) = x(n)

    b(2:n) = a2(2:n) * b(2:n) + a3(1:n-1) * b(1:n-1)
    b(1) = a2(1) * b(1)

  end if

  return
end
subroutine c3_np_sl ( n, a1, a2, a3, b, job )
!
!*******************************************************************************
!
!! C3_NP_SL solves a tridiagonal system factored by C3_NP_FA.
!
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1), the factor information
!    returned by C3_NP_FA.
!
!    Input/output, complex B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  complex b(n)
  integer i
  integer ierror
  integer job
!
!  Check the dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_NP_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a1(i) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a2(i)
      if ( i > 1 ) then
        b(i-1) = b(i-1) - a3(i-1) * b(i)
      end if
    end do

  else
!
!  Solve tranpose ( U ) * Y = B
!
    do i = 1, n
      b(i) = b(i) / a2(i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a3(i) * b(i)
      end if
    end do
!
!  Solve transpose ( L ) * X = Y.
!
    do i = n-1, 1, -1
      b(i) = b(i) - a1(i+1) * b(i+1)
    end do

  end if

  return
end
subroutine c3_print ( n, a1, a2, a3, title )
!
!*******************************************************************************
!
!! C3_PRINT prints a complex tridiagonal matrix.
!
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of
!    the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  complex a1(2:n)
  complex a2(n)
  complex a3(1:n-1)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call c3_print_some ( n, a1, a2, a3, 1, 1, n, n )

  return
end
subroutine c3_print_some ( n, a1, a2, a3, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! C3_PRINT_SOME prints some of a complex tridiagonal matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of
!    the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column, to be printed.
!
  integer, parameter :: incx = 3
!
  integer n
!
  complex a1(2:n)
  complex a2(n)
  complex a3(1:n-1)
  character ( len = 12 ) citemp(incx)
  character ( len = 12 ) crtemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  real xi
  real xr
!
!  Check the dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_PRINT_SOME - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( crtemp(j2), '(i6,6x)' ) j
      write ( citemp(j2), '(i6,6x)' ) j
    end do

    write ( *, '(''Columns:'',6a12)' ) ( crtemp(j2), citemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 1 )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i-j > 1 .or. j-i > 1 ) then

          crtemp(j2) = ' '
          citemp(j2) = ' '

        else

          if ( j == i-1 ) then
            xr = real ( a1(i) )
            xi = aimag ( a1(i) )
          else if ( j == i ) then
            xr = real ( a2(i) )
            xi = aimag ( a2(i) )
          else if ( j == i+1 ) then
            xr = real ( a3(i) )
            xi = aimag ( a3(i) )
          end if

          if ( xr == 0.0E+00 .and. xi == 0.0E+00 ) then
            crtemp(j2) = '    0.0'
            citemp(j2) = ' '
          else if ( xr == 0.0E+00 .and. xi /= 0.0E+00 ) then
            crtemp(j2) = ' '
            write ( citemp(j2), '(g12.5)' ) xi
          else if ( xr /= 0.0E+00 .and. xi == 0.0E+00 ) then
            write ( crtemp(j2), '(g12.5)' ) xr
            citemp(j2) = ' '
          else
            write ( crtemp(j2), '(g12.5)' ) xr
            write ( citemp(j2), '(g12.5)' ) xi
          end if

        end if

      end do

      write ( *, '(i5,1x,6a12)' ) i, ( crtemp(j2), citemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine c3_random ( n, a1, a2, a3 )
!
!*******************************************************************************
!
!! C3_RANDOM returns a random complex tridiagonal matrix.
!
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Output, complex A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of the
!    matrix.  The entries are all between 0 and 1.
!
  integer n
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  integer i
  integer ierror
  real r1
  real r2
!
!  Check dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions for the matrix.'
    return
  end if
!
  do i = 1, n

    if ( i >= 2 ) then
      call r_random ( 0.0E+00, 1.0E+00, r1 )
      call r_random ( 0.0E+00, 1.0E+00, r2 )
      a1(i) = cmplx ( r1, r2 )
    end if

    call r_random ( 0.0E+00, 1.0E+00, r1 )
    call r_random ( 0.0E+00, 1.0E+00, r2 )

    a2(i) = cmplx ( r1, r2 )

    if ( i <= n - 1 ) then
      call r_random ( 0.0E+00, 1.0E+00, r1 )
      call r_random ( 0.0E+00, 1.0E+00, r2 )
      a3(i) = cmplx ( r1, r2 )
    end if

  end do

  return
end
subroutine c3_to_cge ( lda, n, a1, a2, a3, a )
!
!*******************************************************************************
!
!! C3_TO_CGE copies a complex tridiagonal matrix into a general matrix.
!
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of
!    the matrix.
!
!    Output, complex A(LDA,N), the matrix, stored as a general matrix.
!
  integer lda
  integer n
!
  complex a(lda,n)
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_TO_CGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for input matrix.'
    return
  end if

  call cge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_TO_CGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, n

      if ( j == i-1 ) then
        a(i,j) = a1(i)
      else if ( i == j ) then
        a(i,j) = a2(i)
      else if ( j == i+1 ) then
        a(i,j) = a3(i)
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do

  return
end
subroutine c3_vxm ( n, a1, a2, a3, x, b )
!
!*******************************************************************************
!
!! C3_VXM multiplies the transpose of a complex tridiagonal matrix times a vector.
!
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, complex A1(2:N), A2(1:N), A3(1:N-1),
!    the nonzero diagonals of the linear system.
!
!    Input, complex X(N), the vector to be multiplied by Transpose ( A ).
!
!    Output, complex B(N), the product Transpose ( A ) * x.
!
  integer n
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  complex b(n)
  integer i
  integer ierror
  complex x(n)
!
!  Check dimensions.
!
  call c3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'C3_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions for the matrix.'
    return
  end if

  b(1:n) = a2(1:n) * x(1:n)
  b(2:n) = b(2:n) + a3(1:n-1) * x(1:n-1)
  b(1:n-1) = b(1:n-1) + a1(2:n) * x(2:n)

  return
end
subroutine cci_eval ( n, a, lambda )
!
!*******************************************************************************
!
!! CCI_EVAL returns the eigenvalues of a complex circulant matrix.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis,
!    Circulant Matrices,
!    Wiley, 1979.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(N), the entries of the first row of the circulant matrix.
!
!    Output, complex LAMBDA(N), the eigenvalues.
!
  integer n
!
  complex a(n)
  integer i
  complex lambda(n)
  complex w(n)
!
  call cvec_unity ( n, w )

  lambda(1:n) = a(n)
  do i = n-1, 1, -1
    lambda(1:n) = lambda(1:n) * w(1:n) + a(i)
  end do

  call cvec_sort_a2 ( n, lambda )

  return
end
subroutine cci_mxv ( n, a, x, b )
!
!*******************************************************************************
!
!! CCI_MXV multiplies a complex circulant matrix times a vector.
!
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(N), the entries of the first row of the circulant matrix.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A * x.
!
  integer n
!
  complex a(n)
  complex b(n)
  integer i
  integer j
  complex x(n)
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
  do i = 1, n
    b(i) = ZERO
    do j = 1, i-1
      b(i) = b(i) + a(n+j+1-i) * x(j)
    end do
    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do
  end do

  return
end
subroutine cci_print ( n, a, title )
!
!*******************************************************************************
!
!! CCI_PRINT prints a complex circulant matrix.
!
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(N), the N by N circulant matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  complex a(n)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call cci_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine cci_print_some ( n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! CCI_PRINT_SOME prints some of a complex circulant matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(N), the N by N circulant matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 4
!
  integer n
!
  complex a(n)
  complex aij
  character ( len = 20 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(''Columns:'',4a20)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j >= i ) then
          aij = a(j+1-i)
        else
          aij = a(n+j+1-i)
        end if

        if ( aij == ZERO ) then
          ctemp(j2) = '     0.0            '
        else if ( aimag ( aij ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine cci_random ( n, a )
!
!*******************************************************************************
!
!! CCI_RANDOM randomizes a complex circulant matrix.
!
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, complex A(N), the randomized matrix, with entries between
!    0 and 1.
!
  integer n
!
  complex a(n)
  real ai
  real, parameter :: ahi = 1.0E+00
  real, parameter :: alo = 0.0E+00
  real ar
  integer i
!
  do i = 1, n
    call r_random ( alo, ahi, ar )
    call r_random ( alo, ahi, ai )
    a(i) = cmplx ( ar, ai )
  end do

  return
end
subroutine cci_sl ( n, a, b, x, job )
!
!*******************************************************************************
!
!! CCI_SL solves the complex circulant system A * x = b.
!
!
!  Modified:
!
!    07 March 2001
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(N), the entries of the first row of the circulant matrix.
!
!    Input, complex B(N), the right hand side.
!
!    Output, complex X(N), the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  complex a(n)
  complex b(n)
  integer i
  integer job
  integer nsub
  complex r1
  complex r2
  complex r3
  complex r5
  complex r6
  complex work(2*n-2)
  complex x(n)
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
  if ( job == 0 ) then
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = ZERO
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(n+2-nsub)
      r6 = a(nsub)

      if ( nsub > 2 ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(n+1-i) * work(nsub-i)
          r6 = r6 + a(i+1) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( nsub > 2 ) then

        r6 = work(n)
        work(n-1+nsub-1) = ZERO
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = ZERO
      do i = 1, nsub-1
        r5 = r5 + a(n+1-i) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub-1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  else
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = ZERO
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(nsub)
      r6 = a(n+2-nsub)

      if ( nsub > 2 ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(i+1) * work(nsub-i)
          r6 = r6 + a(n+1-i) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( nsub > 2 ) then

        r6 = work(n)
        work(n-1+nsub-1) = ZERO
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = ZERO
      do i = 1, nsub-1
        r5 = r5 + a(i+1) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub-1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  end if

  return
end
subroutine cci_to_cge ( lda, n, a, a2 )
!
!*******************************************************************************
!
!! CCI_TO_CGE copies a complex circulant matrix into a general matrix.
!
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(N), the circulant matrix.
!
!    Output, complex A2(LDA,N), the circulant matrix, stored as
!    a general matrix.
!
  integer lda
  integer n
!
  complex a(n)
  complex a2(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call cge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CCI_TO_CGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, i-1
      a2(i,j) = a(n+j+1-i)
    end do
    do j = i, n
      a2(i,j) = a(j+1-i)
    end do
  end do

  return
end
subroutine cci_vxm ( n, a, x, b )
!
!*******************************************************************************
!
!! CCI_VXM multiplies a vector times a complex circulant matrix.
!
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(N), the entries of the first row of the circulant matrix.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product Transpose ( A ) * X.
!
  integer n
!
  complex a(n)
  complex b(n)
  integer i
  integer j
  complex x(n)
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
  b(1:n) = ZERO

  do i = 1, n
    do j = 1, i
      b(i) = b(i) + a(i+1-j) * x(j)
    end do
    do j = i+1, n
      b(i) = b(i) + a(n+i+1-j) * x(j)
    end do
  end do

  return
end
subroutine cge_check ( lda, m, n, ierror )
!
!*******************************************************************************
!
!! CGE_CHECK checks the dimensions of a complex general matrix.
!
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  integer ierror
  integer lda
  integer m
  integer n
!
  ierror = 0

  if ( lda < m ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'CGE_CHECK - Illegal LDA = ', lda
  end if

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'CGE_CHECK - Illegal M = ', m
  end if

  if ( n < 1 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'CGE_CHECK - Illegal N = ', n
  end if

  return
end
subroutine cto_mxv ( n, a, x, b )
!
!*******************************************************************************
!
!! CTO_MXV multiplies a complex Toeplitz matrix times a vector.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A * x.
!
  integer n
!
  complex a(2*n-1)
  complex b(n)
  integer i
  integer j
  complex x(n)
!
  do i = 1, n

    b(i) = cmplx ( 0.0E+00, 0.0E+00 )

    do j = 1, i-1
      b(i) = b(i) + a(n+i-j) * x(j)
    end do

    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do

  end do

  return
end
subroutine cto_print ( n, a, title )
!
!*******************************************************************************
!
!! CTO_PRINT prints a complex Toeplitz matrix.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  complex a(2*n-1)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call cto_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine cto_print_some ( n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! CTO_PRINT_SOME prints some of a complex Toeplitz matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, complex A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 4
!
  integer n
!
  complex a(2*n-1)
  complex aij
  character ( len = 20 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(''Columns:'',4a20)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j >= i ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        if ( aij == ZERO ) then
          ctemp(j2) = '    0.0'
        else if ( aimag ( aij ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine cto_random ( n, a )
!
!*******************************************************************************
!
!! CTO_RANDOM randomizes a complex Toeplitz matrix.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, complex A(2*N-1), the randomized matrix, with entries between
!    0 and 1.
!
  integer n
!
  complex a(2*n-1)
  integer i
!
  call cvec_random ( 0.0E+00, 1.0E+00, 2*n-1, a )

  return
end
subroutine cto_sl ( n, a, b, x, job )
!
!***********************************************************************
!
!! CTO_SL solves the complex Toeplitz system A * X = B.
!
!
!  Modified:
!
!    11 March 2001
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(2*N-1), the first row of the Toeplitz matrix, followed
!    by the first column of the Toeplitz matrix, beginning with the second
!    element.
!
!    Input, complex B(N) the right hand side vector.
!
!    Output, complex X(N), the solution vector.  X and B may share the
!    same storage.
!
!    Input, integer JOB,
!    0 to solve A*X=B,
!    nonzero to solve Transpose(A)*X=B.
!
  integer n
!
  complex a(2*n-1)
  complex b(n)
  complex c1(n-1)
  complex c2(n-1)
  integer i
  integer job
  integer nsub
  complex r1
  complex r2
  complex r3
  complex r5
  complex r6
  complex x(n)
  complex, parameter :: ZERO = cmplx ( 0.0E+00, 0.0E+00 )
!
  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = a(n+nsub-1)
      r6 = a(nsub)
    else
      r5 = a(nsub)
      r6 = a(n+nsub-1)
    end if

    if ( nsub > 2 ) then

      c1(nsub-1) = r2

      do i = 1, nsub-2
        if ( job == 0 ) then
          r5 = r5 + a(n+i) * c1(nsub-i)
          r6 = r6 + a(i+1) * c2(i)
        else
          r5 = r5 + a(i+1) * c1(nsub-i)
          r6 = r6 + a(n+i) * c2(i)
        end if
      end do

    end if

    r2 = -r5 / r1
    r3 = -r6 / r1
    r1 = r1 + r5 * r3

    if ( nsub > 2 ) then

      r6 = c2(1)
      c2(nsub-1) = ZERO
      do i = 2, nsub-1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    r5 = ZERO

    do i = 1, nsub-1
      if ( job == 0 ) then
        r5 = r5 + a(n+i) * x(nsub-i)
      else
        r5 = r5 + a(i+1) * x(nsub-i)
      end if
    end do

    r6 = ( b(nsub) - r5 ) / r1

    do i = 1, nsub-1
      x(i) = x(i) + c2(i) * r6
    end do

    x(nsub) = r6

  end do

  return
end
subroutine cto_to_cge ( lda, n, a, a2 )
!
!*******************************************************************************
!
!! CTO_TO_CGE copies a complex Toeplitz matrix into a general matrix.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(2*N-1), the Toeplitz matrix.
!
!    Output, complex A2(LDA,N), the Toeplitz matrix, stored as
!    a general matrix.
!
  integer lda
  integer n
!
  complex a(2*n-1)
  complex a2(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call cge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CTO_TO_CGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, i-1
      a2(i,j) = a(n+i-j)
    end do
    do j = i, n
      a2(i,j) = a(j-i+1)
    end do
  end do

  return
end
subroutine cto_vxm ( n, a, x, b )
!
!*******************************************************************************
!
!! CTO_VXM multiplies a vector times a complex Toeplitz matrix.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, complex A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product Transpose ( A ) * X.
!
  integer n
!
  complex a(2*n-1)
  complex b(n)
  integer i
  integer j
  complex x(n)
!
  do i = 1, n

    b(i) = cmplx ( 0.0E+00, 0.0E+00 )

    do j = 1, i
      b(i) = b(i) + a(i+1-j) * x(j)
    end do

    do j = i+1, n
      b(i) = b(i) + a(n+j-i) * x(j)
    end do

  end do

  return
end
subroutine cvec_identity ( n, a )
!
!*******************************************************************************
!
!! CVEC_IDENTITY sets a complex vector to a sort of identity vector.
!
!
!  Discussion:
!
!    X(1:N) = (0:N-1) * exp ( 2 * PI * I / N )
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, complex A(N), the array to be initialized.
!
  integer n
!
  complex a(n)
  real ai
  real ar
  integer i
  real, parameter :: &
    pi = 3.14159265358979323846264338327950288419716939937510E+00
  real theta
!
  do i = 1, n
    theta = pi * real ( 2 * ( i - 1 ) ) / real ( n )
    ar = real ( i ) * cos ( theta )
    ai = real ( i ) * sin ( theta )
    a(i) = cmplx ( ar, ai )
  end do

  return
end
subroutine cvec_print ( n, a, title )
!
!*******************************************************************************
!
!! CVEC_PRINT prints a complex vector, with an optional title.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, complex A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  complex a(n)
  integer i
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine cvec_print_some ( n, x, max_print )
!
!*******************************************************************************
!
!! CVEC_PRINT_SOME prints some of a complex vector.
!
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    14 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, complex X(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
  integer n
!
  integer i
  integer max_print
  complex x(n)
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,2g14.6)' ) i, x(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,2g14.6)' ) i, x(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i6,2x,2g14.6)' ) i, x(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,2g14.6)' ) i, x(i)
    end do
    i = max_print
    write ( *, '(i6,2x,2g14.6,2x,a)' ) i, x(i), '...more entries...'

  end if

  return
end
subroutine cvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! CVEC_RANDOM returns a random complex vector in a given range.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, complex A(N), the vector of randomly chosen values.
!
  integer n
!
  complex a(n)
  real ahi
  real ai
  real alo
  real ar
  integer i
!
  do i = 1, n

    call r_random ( alo, ahi, ar )
    call r_random ( alo, ahi, ai )

    a(i) = cmplx ( ar, ai )

  end do

  return
end
subroutine cvec_sort_a2 ( n, x )
!
!*******************************************************************************
!
!! CVEC_SORT_A2 ascending sorts a complex array by L2 norm.
!
!
!  Discussion:
!
!    The L2 norm of A+Bi is sqrt ( A**2 + B**2 ).
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input/output, complex X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  integer n
!
  logical c_le_l2
  integer i
  integer indx
  integer isgn
  integer j
  complex x(n)
!
  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( indx > 0 ) then

      call c_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c_le_l2 ( x(i), x(j) ) ) then
        isgn = - 1
      else
        isgn = + 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine cvec_unity ( n, a )
!
!*******************************************************************************
!
!! CVEC_UNITY returns the N roots of unity.
!
!
!  Discussion:
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, complex A(N), the N roots of unity.
!
  integer n
!
  complex a(n)
  integer i
  real, parameter :: &
    pi = 3.14159265358979323846264338327950288419716939937510E+00
  real theta
!
  do i = 1, n
    theta = pi * real ( 2 * ( i - 1 ) ) / real ( n )
    a(i) = cmplx ( cos ( theta ), sin ( theta ) )
  end do

  return
end
subroutine hilb_inv ( lda, n, a )
!
!*******************************************************************************
!
!! HILB_INV returns the inverse of the Hilbert matrix.
!
!
!  Formula:
!
!    A(I,J) =  (-1)**(I+J) * (N+I-1)! * (N+J-1)! /
!           [ (I+J-1) * ((I-1)!*(J-1)!)**2 * (N-I)! * (N-J)! ]
!
!  Example:
!
!    N = 5
!
!       25    -300     1050    -1400     630
!     -300    4800   -18900    26880  -12600
!     1050  -18900    79380  -117600   56700
!    -1400   26880  -117600   179200  -88200
!      630  -12600    56700   -88200   44100
!
!  Properties:
!
!    A is symmetric.
!
!    Because A is symmetric, it is normal, so diagonalizable.
!
!    A is almost impossible to compute accurately by general routines
!    that compute the inverse.
!
!    A is integral.
!
!    The sum of the entries of A is N**2.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of A.
!
!    Output, real A(LDA,N), the inverse Hilbert matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
!
!  Set the (1,1) entry.
!
  a(1,1) = real ( n**2 )
!
!  Define Row 1, Column J by recursion on Row 1 Column J-1
!
  i = 1
  do j = 2, n
    a(i,j) = - a(i,j-1) * real ( ( n + j - 1 ) * ( i + j - 2 ) * &
      ( n + 1 - j ) ) / real ( ( i + j - 1 ) * ( j - 1 )**2 )
  end do
!
!  Define Row I by recursion on row I-1
!
  do i = 2, n
    do j = 1, n

      a(i,j) = - a(i-1,j) * real ( (n+i-1) * (i+j-2) * (n+1-i) ) / &
        real ( (i+j-1) * (i-1)**2 )

    end do
  end do

  return
end
subroutine i_random ( ilo, ihi, i )
!
!*******************************************************************************
!
!! I_RANDOM returns a random integer in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ILO, IHI, the minimum and maximum acceptable values.
!
!    Output, integer I, the randomly chosen integer.
!
  logical, save :: seed = .false.
  integer i
  integer ihi
  integer ilo
  real r
  real rhi
  real rlo
!
  if ( .not. seed ) then
    call random_seed
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = r )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo ) - 0.5E+00
  rhi = real ( ihi ) + 0.5E+00
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0E+00 - r ) * rlo + r * rhi )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
function r_is_int ( r )
!
!*******************************************************************************
!
!! R_IS_INT determines if a real number represents an integer value.
!
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the number to be checked.
!
!    Output, logical R_IS_INT, is TRUE if R is an integer value.
!
  integer i
  real r
  logical r_is_int
!
  if ( r > real ( huge ( i ) ) ) then
    r_is_int = .false.
  else if ( r < - real ( huge ( i ) ) ) then
    r_is_int = .false.
  else if ( r == real ( int ( r ) ) ) then
    r_is_int = .true.
  else
    r_is_int = .false.
  end if

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  logical, save :: seed = .false.
  real r
  real rhi
  real rlo
  real t
!
  if ( .not. seed ) then
    call random_seed
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine rcol_swap ( lda, m, n, a, i, j )
!
!*******************************************************************************
!
!! RCOL_SWAP swaps columns I and J of a real array of column data.
!
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      A = (
!        1.  4.  3.  2.
!        5.  8.  7.  6.
!        9. 12. 11. 10. )
!
!  Modified:
!
!    22 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), the M by N array.
!
!    Input, integer I, J, the columns to be swapped.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer k
  integer m
!
  if ( 1 <= i .and. i <= n .and. 1 <= j .and. j <= n ) then

    do k = 1, m
      call r_swap ( a(k,i), a(k,j) )
    end do

  else

    write ( *, * ) ' '
    write ( *, * ) 'RCOL_SWAP - Fatal error!'
    write ( *, * ) '  I or J is out of bounds.'
    write ( *, * ) '  I =    ', i
    write ( *, * ) '  J =    ', j
    write ( *, * ) '  NCOL = ', n
    stop

  end if

  return
end
subroutine rmat_diag_add_scalar ( lda, n, a, s )
!
!*******************************************************************************
!
!! RMAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of a matrix.
!
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the matrix.
!
!    Input, integer N, the number of rows and columns of the matrix.
!
!    Input/output, real A(LDA,N), the N by N matrix to be modified.
!
!    Input, real S, the value to be added to the diagonal of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  real s
!
  do i = 1, n
    a(i,i) = a(i,i) + s
  end do

  return
end
subroutine rmat_diag_get_vector ( lda, n, a, v )
!
!*******************************************************************************
!
!! RMAT_DIAG_GET_VECTOR gets the value of the diagonal of a matrix.
!
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the matrix.
!
!    Input, integer N, the number of rows and columns of the matrix.
!
!    Input, real A(LDA,N), the N by N matrix.
!
!    Output, real V(N), the diagonal entries of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  real v(n)
!
  do i = 1, n
    v(i) = a(i,i)
  end do

  return
end
subroutine rmat_diag_set_scalar ( lda, n, a, s )
!
!*******************************************************************************
!
!! RMAT_DIAG_SET_SCALAR sets the diagonal of a matrix to a scalar value.
!
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the matrix.
!
!    Input, integer N, the number of rows and columns of the matrix.
!
!    Input/output, real A(LDA,N), the N by N matrix to be modified.
!
!    Input, real S, the value to be assigned to the diagonal of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  real s
!
  do i = 1, n
    a(i,i) = s
  end do

  return
end
subroutine rmat_random ( alo, ahi, lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_RANDOM returns a matrix of uniform random values between AHI and ALO.
!
!
!  Modified:
!
!    01 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the minimum and maximum values that
!    the matrix entries can have.
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, N, the number of rows and columns of A.
!
!    Output, real A(LDA,N), the random matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real ahi
  real alo
  integer i
  integer j
  integer m
!
  do i = 1, m
    do j = 1, n
      call r_random ( alo, ahi, a(i,j) )
    end do
  end do

  return
end
subroutine rrow_swap ( lda, m, n, a, irow1, irow2 )
!
!*******************************************************************************
!
!! RROW_SWAP swaps two rows of a table.
!
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the first dimension of A.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real A(LDA,N), the M by N array.
!
!    Input, integer IROW1, IROW2, the two rows to swap.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer irow1
  integer irow2
  integer j
  integer m
!
  if ( irow1 < 1 .or. irow1 > m ) then
    write ( *, * ) ' '
    write ( *, * ) 'RROW_SWAP - Fatal error!'
    write ( *, * ) '  IROW1 is out of range.'
    stop
  end if

  if ( irow2 < 1 .or. irow2 > m ) then
    write ( *, * ) ' '
    write ( *, * ) 'RROW_SWAP - Fatal error!'
    write ( *, * ) '  IROW2 is out of range.'
    stop
  end if

  if ( irow1 == irow2 ) then
    return
  end if

  do j = 1, n
    call r_swap ( a(irow1,j), a(irow2,j) )
  end do

  return
end
subroutine rvec2_print_some ( n, x1, x2, max_print )
!
!*******************************************************************************
!
!! RVEC2_PRINT_SOME prints some of two real vectors.
!
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vectors, is no more than MAX_PRINT, then
!    the entire vectors are printed, one entry of each per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    10 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vectors.
!
!    Input, real X1(N), X2(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
  integer n
!
  integer i
  integer max_print
  real x1(n)
  real x2(n)
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    write ( *, '(a)' ) '......  ..............  ..............'
    i = n
    write ( *, '(i6,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    i = max_print
    write ( *, '(i6,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), &
      '...more entries...'

  end if

  return
end
subroutine rvec_identity ( n, a )
!
!*******************************************************************************
!
!! RVEC_IDENTITY sets a real vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real A(N), the array to be initialized.
!
  integer n
!
  real a(n)
  integer i
!
  do i = 1, n
    a(i) = real ( i )
  end do

  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real vector, with an optional title.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  real a(n)
  integer i
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine rvec_print_some ( n, x, max_print )
!
!*******************************************************************************
!
!! RVEC_PRINT_SOME prints some of a real vector.
!
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    10 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real X(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
  integer n
!
  integer i
  integer max_print
  real x(n)
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,g14.6)' ) i, x(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,g14.6)' ) i, x(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i6,2x,g14.6)' ) i, x(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,g14.6)' ) i, x(i)
    end do
    i = max_print
    write ( *, '(i6,2x,g14.6,2x,a)' ) i, x(i), '...more entries...'

  end if

  return
end
subroutine rvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_RANDOM returns a random real vector in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(N), the vector of randomly chosen values.
!
  integer n
!
  real a(n)
  real ahi
  real alo
  integer i
!
  do i = 1, n

    call r_random ( alo, ahi, a(i) )

  end do

  return
end
subroutine s3_check ( n, ierror )
!
!*******************************************************************************
!
!! S3_CHECK checks the dimensions of a real tridiagonal matrix.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Output, integer IERROR, error flag.
!    0, no errors detected.
!    1, N was less than 2.
!
  integer ierror
  integer n
!
  ierror = 0

  if ( n < 2 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'S3_CHECK - Fatal error!'
    write ( *, * ) '  N must be at least 2.'
    write ( *, * ) '  The input N was ', n
  end if

  return
end
subroutine s3_cr_fa ( n, subd, diag, supd )
!
!*******************************************************************************
!
!! S3_CR_FA decomposes a real tridiagonal matrix using cyclic reduction.
!
!
!  Discussion:
!
!    Once S3_CR_FA has decomposed a matrix A, then S3_CR_SL may be used to solve
!    linear systems A * x = b.
!
!    On a Cray, S3_CR_FA can be very much faster than standard Gauss
!    elimination techniques, such as SGTSL from LINPACK, which do not
!    vectorize well.
!
!    S3_CR_FA and S3_CR_SL will be slower than the Cray SCILIB routine TRID.
!    On the other hand, TRID does not provide a factorization, and source
!    code for TRID is not generally available.
!
!    S3_CR_FA does not employ pivoting.  Hence, the results can be more
!    sensitive to ill-conditioning than standard Gauss elimination.  In
!    particular, S3_CR_FA will fail if any diagonal element of the matrix
!    is zero.  Other matrices may also cause S3_CR_FA to fail.
!
!    S3_CR_FA can be guaranteed to work properly if the matrix is strictly
!    diagonally dominant, that is, if the absolute value of the diagonal
!    element is strictly greater than the sum of the absolute values of
!    the offdiagonal elements, for each equation.
!
!    The algorithm may be illustrated by the following figures:
!
!    The initial matrix is given by:
!
!          D1 U1
!          L1 D2 U2
!             L2 D3 U3
!                L3 D4 U4
!                   L4 D5 U5
!                      L5 D6
!
!    Rows and columns are permuted in an odd/even way to yield:
!
!          D1       U1
!             D3    L2 U3
!                D5    L4 U5
!          L1 U2    D2
!             L3 U4    D4
!                L5       D6
!
!    A block LU decomposition is performed to yield:
!
!          D1      |U1
!             D3   |L2 U3
!                D5|   L4 U5
!          --------+--------
!                  |D2'F3
!                  |F1 D4'F4
!                  |   F2 D6'
!
!    For large systems, this reduction is repeated on the lower right hand
!    tridiagonal subsystem until a completely upper triangular system
!    is obtained.  The system has now been factored into the product of a
!    lower triangular system and an upper triangular one, and the information
!    defining this factorization may be used by S3_CR_SL to solve linear systems.
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real SUBD(0:2*N), DIAG(2*N), SUPD(0:2*N).
!
!    On input, SUBD, DIAG and SUPD contain the subdiagonal, diagonal,
!    and superdiagonal entries of the matrix.
!
!    The diagonal entries are stored in DIAG(1) through DIAG(N).
!    The subdiagonal entries are stored in SUBD(1) through SUBD(N-1).
!    The superdiagonal entries are stored in SUPD(1) through SUPD(n-1).
!    The extra entries in DIAG, SUBD, and SUPD need not be initialized
!    by the user.
!
!    On output, SUBD, DIAG, and SUPD contain information defining the
!    factorization of the tridiagonal matrix.  This information
!    will be needed by S3_CR_SL to solve linear systems.
!
!    The extra positions in the arrays, DIAG(N+1) through DIAG(2*N),
!    SUBD(0), SUBD(N) through SUBD(2*N), SUPD(0) and SUPD(N) through
!    SUPD(2*N), are used for workspace and storage.
!
!    In particular, SUBD and SUPD must be declared with an initial,
!    0-th element, or the algorithm will not work.
!
  integer n
!
  real diag(2*n)
  integer i
  integer iful
  integer ifulp
  integer ihaf
  integer il
  integer ilp
  integer inc
  integer incr
  integer ipnt
  integer ipntp
  real subd(0:2*n)
  real supd(0:2*n)
!
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_CR_FA - Fatal error!'
    write ( *, * ) '  Nonpositive N = ', n
    return
  end if

  if ( n == 1 ) then
    diag(1) = 1.0E+00 / diag(1)
    return
  end if
!
!  Zero out the workspace entries.
!
  subd(0) = 0.0E+00
  supd(0) = 0.0E+00
  subd(n) = 0.0E+00
  supd(n) = 0.0E+00

  diag(n+1:2*n) = 0.0E+00
  subd(n+1:2*n) = 0.0E+00
  supd(n+1:2*n) = 0.0E+00

  il = n
  subd(il) = 0.0E+00
  supd(il) = 0.0E+00
  ipntp = 0

  do while ( il > 1 )

    ipnt = ipntp
    ipntp = ipntp + il
    if ( mod(il,2) == 1 ) then
      inc = il + 1
    else
      inc = il
    end if

    incr = inc / 2
    il = il / 2
    ihaf = ipntp + incr + 1
    ifulp = ipnt + inc + 2

!dir$ ivdep
    do ilp = incr, 1, -1
      ifulp = ifulp - 2
      iful = ifulp - 1
      ihaf = ihaf - 1
      diag(iful) = 1.0E+00 / diag(iful)
      subd(iful) = subd(iful) * diag(iful)
      supd(ifulp) = supd(ifulp) * diag(ifulp+1)
      diag(ihaf) = diag(ifulp) - supd(iful) * subd(iful) &
        - supd(ifulp) * subd(ifulp)
      subd(ihaf) = - subd(ifulp) * subd(ifulp+1)
      supd(ihaf) = - supd(ifulp) * supd(ifulp+1)
    end do

  end do

  diag(ipntp+1) = 1.0E+00 / diag(ipntp+1)

  return
end
subroutine s3_cr_sl ( n, subd, diag, supd, rhs )
!
!*******************************************************************************
!
!! S3_CR_SL solves a real linear system factored by S3_CR_FA.
!
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  S3_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by S3_CR_FA are passed to S3_CR_SL, then one or many
!    linear systems involving the matrix A may be solved.
!
!
!    The form of the equations may be summarized as:
!
!    DIAG(1)  *X(1)   + SUPD(1)  *X(2)                  = RHS(1)
!
!    SUBD(I-1)*X(I-1) + DIAG(I)  *X(I) + SUPD(I)*X(I+1) = RHS(I)  I=2 to N-1,
!
!                       SUBD(N-1)*X(N-1) +DIAG(N)*X(N)  = RHS(N)
!
!    The names used in this example exactly correspond to the initial
!    storage of information in the arrays SUBD, DIAG, SUPD and RHS.
!
!    Here SUBD, DIAG and SUPD are the lower diagonal, diagonal, and
!    upper diagonal coefficients of the tridiagonal system.
!
!
!    Note that S3_CR_FA does not perform pivoting, and so the solution produced
!    by S3_CR_SL may be less accurate than a solution produced by a standard
!    Gauss algorithm.  However, such problems can be guaranteed not to occur
!    if the matrix A is strictly diagonally dominant, that is, if the
!    absolute value of the diagonal coefficient is greater than the sum of
!    the absolute values of the two off diagonal coefficients, for each
!    row of the matrix.
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real SUBD(0:2*N), DIAG(2*N), SUPD(0:2N).  These arrays contain
!    information describing the factorization of the original tridiagonal
!    matrix, as computed by S3_CR_FA.
!
!    Input/output, real RHS(0:2*N).
!
!    On input, RHS contains the right hand side vector in locations
!    1 thorugh N.  The zero-th entry, and the second N locations are
!    used for workspace.
!
!    On output, the locations 1 through N of RHS contain the solution
!    of the linear system.
!
  integer n
!
  real diag(2*n)
  integer i
  integer iful
  integer ifulm
  integer ifulp
  integer ihaf
  integer il
  integer ipnt
  integer ipntp
  real subd(0:2*n)
  integer ndiv
  real supd(0:2*n)
  real rhs(0:2*n)
!
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_CR_SL - Fatal error!'
    write ( *, * ) '  Nonpositive N = ', n
    return
  end if

  if ( n == 1 ) then
    rhs(1) = diag(1) * rhs(1)
    return
  end if
!
!  Zero out workspace entries of RHS.
!
  rhs(0) = 0.0E+00
  rhs(n+1:2*n) = 0.0E+00
  subd(0) = 0.0E+00
  supd(0) = 0.0E+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( il > 1 )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

!dir$ ivdep
    do iful = ipnt+2, ipntp, 2
      ihaf = ihaf + 1
      ifulp = iful + 1
      ifulm = iful - 1
      rhs(ihaf) = rhs(iful) - subd(ifulm) * rhs(ifulm) - supd(iful) * rhs(ifulp)
    end do

  end do

  rhs(ihaf) = rhs(ihaf) * diag(ihaf)
  ipnt = ipntp

  do while ( ipnt > 0 )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

!dir$ ivdep
    do ifulm = ipnt+1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful) = rhs(ihaf)
      rhs(ifulm) = diag(ifulm) * ( rhs(ifulm) - subd(ifulm-1) * rhs(ifulm-1) &
        - supd(ifulm) * rhs(iful) )
    end do

  end do

  return
end
subroutine s3_gs_sl ( n, a1, a2, a3, b, x, maxit, job )
!
!*******************************************************************************
!
!! S3_GS_SL tries to solve a tridiagonal system using Gauss-Seidel iteration.
!
!
!  Discussion:
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the subdiagonal, diagonal,
!    and superdiagonal of the matrix.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Input/output, real X(N), an approximate solution to the system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer ierror
  integer job
  integer maxit
  integer numit
  real x(n)
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_GS_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a2(i) == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'S3_GS_SL - Fatal error!'
      write ( *, * ) '  Zero diagonal entry, index = ', i
      return
    end if
  end do

  if ( job == 0 ) then

    do numit = 1, maxit

      x(1) =   ( b(1)                  - a3(1) * x(2)   ) / a2(1)
      do i = 2, n - 1
        x(i) = ( b(i) - a1(i) * x(i-1) - a3(i) * x(i+1) ) / a2(i)
      end do
      x(n) =   ( b(n) - a1(n) * x(n-1)                  ) / a2(n)

    end do

  else

    do numit = 1, maxit

      x(1) =   ( b(1)                    - a1(2) * x(2) )     /a2(1)
      do i = 2, n - 1
        x(i) = ( b(i) - a3(i-1) * x(i-1) - a1(i+1) * x(i+1) ) /a2(i)
      end do
      x(n) =   ( b(n) - a3(n-1) * x(n-1) )                    /a2(n)

    end do

  end if

  return
end
subroutine s3_jac_sl ( n, a1, a2, a3, b, x, maxit, job )
!
!*******************************************************************************
!
!! S3_JAC_SL tries to solve a tridiagonal system using Jacobi iteration.
!
!
!  Discussion:
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Modified:
!
!    12 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the subdiagonal, diagonal,
!    and superdiagonal of the matrix.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Input/output, real X(N), an approximate solution to the system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer ierror
  integer job
  integer maxit
  integer numit
  real x(n)
  real xnew(n)
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_JAC_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a2(i) == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'S3_JAC_SL - Fatal error!'
      write ( *, * ) '  Zero diagonal entry, index = ', i
      return
    end if
  end do

  if ( job == 0 ) then

    do numit = 1, maxit

      xnew(1) =   b(1)                  - a3(1) * x(2)
      do i = 2, n - 1
        xnew(i) = b(i) - a1(i) * x(i-1) - a3(i) * x(i+1)
      end do
      xnew(n) =   b(n) - a1(n) * x(n-1)

      xnew(1:n) = xnew(1:n) / a2(1:n)

      x(1:n) = xnew(1:n)

    end do

  else

    do numit = 1, maxit

      xnew(i) =   b(1)                    - a1(2) * x(2)
      do i = 2, n - 1
        xnew(i) = b(i) - a3(i-1) * x(i-1) - a1(i+1) * x(i+1)
      end do
      xnew(n) =   b(n) - a3(n-1) * x(n-1)

      xnew(1:n) = xnew(1:n) / a2(1:n)

      x(1:n) = xnew(1:n)

    end do

  end if

  return
end
subroutine s3_mxv ( n, a1, a2, a3, x, b )
!
!*******************************************************************************
!
!! S3_MXV multiplies a tridiagonal matrix times a vector.
!
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1),
!    the nonzero diagonals of the linear system.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer ierror
  real x(n)
!
!  Check dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions for the matrix.'
    return
  end if

  b(1:n) = a2(1:n) * x(1:n)

  b(2:n) = b(2:n) + a1(2:n) * x(1:n-1)

  b(1:n-1) = b(1:n-1) + a3(1:n-1) * x(2:n)

  return
end
subroutine s3_np_det ( n, a2, det )
!
!*******************************************************************************
!
!! S3_NP_DET returns the determinant of a tridiagonal system factored by S3_NP_FA.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real A2(N), the factor information computed by S3_NP_FA.
!
!    Output, real DET, the determinant of the matrix.
!
  integer n
!
  real a2(n)
  real det
  integer i
  integer ierror
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  det = product ( a2(1:n) )

  return
end
subroutine s3_np_fa ( n, a1, a2, a3, info )
!
!*******************************************************************************
!
!! S3_NP_FA factors a tridiagonal system without pivoting.
!
!
!  Discussion:
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!    S3_NP_FA and S3_NP_SL may be preferable to the corresponding
!    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
!    in one step, and does not save the factorization.
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output real A1(2:N), A2(1:N), A3(1:N-1), the subdiagonal,
!    diagonal, and superdiagonal of the matrix.  On output, these are
!    overwritten by factorization information.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  integer i
  integer ierror
  integer info
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do i = 1, n-1

    if ( a2(i) == 0.0E+00 ) then
      info = i
      write ( *, * ) ' '
      write ( *, * ) 'S3_NP_FA - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if

    a1(i+1) = a1(i+1) / a2(i)
    a2(i+1) = a2(i+1) - a1(i+1) * a3(i)

  end do

  if ( a2(n) == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_FA - Fatal error!'
    write ( *, * ) '  Zero pivot on step ', info
    return
  end if

  return
end
subroutine s3_np_fs ( n, a1, a2, a3, b, x )
!
!*******************************************************************************
!
!! S3_NP_FS factors and solves a tridiagonal linear system.
!
!
!  Note:
!
!    This algorithm requires that each diagonal entry be nonzero.
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real A1(2:N), A2(1:N), A3(1:N-1).
!    On input, the nonzero diagonals of the linear system.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input/output, real B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B has been overwritten by factorization information.
!
!    Output, real X(N), the solution of the linear system.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer ierror
  real x(n)
  real xmult
!
!  Check dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_FS - Fatal error!'
    write ( *, * ) '  Illegal dimensions for input matrix.'
    return
  end if
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a2(i) == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'S3_NP_FS - Fatal error!'
      write ( *, * ) '  A2(', i, ') = 0.'
      return
    end if
  end do

  do i = 2, n-1

    xmult = a1(i) / a2(i-1)
    a2(i) = a2(i) - xmult * a3(i-1)

    b(i) = b(i) - xmult * b(i-1)

  end do

  xmult = a1(n) / a2(n-1)
  a2(n) = a2(n) - xmult * a3(n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
  do i = n-1, 1, -1
    x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
  end do

  return
end
subroutine s3_np_ml ( n, a1, a2, a3, x, b, job )
!
!*******************************************************************************
!
!! S3_NP_ML computes A * x or x * A, where A has been factored by S3_NP_FA.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the LU factors from S3_FA.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product.
!
!    Input, integer JOB, specifies the product to find.
!    0, compute A * x.
!    nonzero, compute transpose ( A ) * x.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer ierror
  integer job
  real x(n)
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute X := U * X
!
    do i = 1, n

      b(i) = a2(i) * b(i)

      if ( i < n ) then
        b(i) = b(i) + a3(i) * b(i+1)
      end if

    end do
!
!  Compute X: = L * X.
!
    do i = n, 2, -1
      b(i) = b(i) + a1(i) * b(i-1)
    end do

  else
!
!  Compute X: = transpose ( L ) * X.
!
    do i = 1, n-1
      b(i) = b(i) + a1(i+1) * b(i+1)
    end do
!
!  Compute X: = transpose ( U ) * X.
!
    do i = n, 1, -1
      b(i) = a2(i) * b(i)
      if ( i > 1 ) then
        b(i) = b(i) + a3(i-1) * b(i-1)
      end if
    end do

  end if

  return
end
subroutine s3_np_sl ( n, a1, a2, a3, b, job )
!
!*******************************************************************************
!
!! S3_NP_SL solves a tridiagonal system factored by S3_NP_FA.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the factor information
!    returned by S3_NP_FA.
!
!    Input/output, real B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer ierror
  integer job
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_NP_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a1(i) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a2(i)
      if ( i > 1 ) then
        b(i-1) = b(i-1) - a3(i-1) * b(i)
      end if
    end do

  else
!
!  Solve tranpose ( U ) * Y = B
!
    do i = 1, n
      b(i) = b(i) / a2(i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a3(i) * b(i)
      end if
    end do
!
!  Solve transpose ( L ) * X = Y.
!
    do i = n-1, 1, -1
      b(i) = b(i) - a1(i+1) * b(i+1)
    end do

  end if

  return
end
subroutine s3_print ( n, a1, a2, a3, title )
!
!*******************************************************************************
!
!! S3_PRINT prints a tridiagonal matrix.
!
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of
!    the matrix.
!
!    Input, character ( len = * ) TITLE, a title to print.
!
  integer n
!
  real a1(2:n)
  real a2(n)
  real a3(1:n-1)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call s3_print_some ( n, a1, a2, a3, 1, 1, n, n )

  return
end
subroutine s3_print_some ( n, a1, a2, a3, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! S3_PRINT_SOME prints some of a tridiagonal matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of
!    the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column, to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a1(2:n)
  real a2(n)
  real a3(1:n-1)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 1 )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i-j > 1 .or. j-i > 1 ) then
          ctemp(j2) = '              '
        else if ( j == i-1 ) then
          if ( r_is_int ( a1(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a1(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a1(i)
          end if
        else if ( j == i ) then
          if ( r_is_int ( a2(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a2(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a2(i)
          end if
        else if ( j == i+1 ) then
          if ( r_is_int ( a3(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a3(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a3(i)
          end if
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine s3_random ( n, a1, a2, a3 )
!
!*******************************************************************************
!
!! S3_RANDOM returns a random tridiagonal matrix.
!
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Output, real A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of the
!    matrix.  The entries are all between 0 and 1.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  integer ierror
!
!  Check dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions for the matrix.'
    return
  end if

  call rvec_random ( 0.0E+00, 1.0E+00, n-1, a1(2:n) )
  call rvec_random ( 0.0E+00, 1.0E+00, n,   a2(1:n) )
  call rvec_random ( 0.0E+00, 1.0E+00, n-1, a3(1:n-1) )

  return
end
subroutine s3_to_sge ( lda, n, a1, a2, a3, a )
!
!*******************************************************************************
!
!! S3_TO_SGE copies a tridiagonal matrix into a general matrix.
!
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1), the nonzero diagonals of
!    the matrix.
!
!    Output, real A(LDA,N), the matrix, stored as a general matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for input matrix.'
    return
  end if

  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, n

      if ( j == i-1 ) then
        a(i,j) = a1(i)
      else if ( i == j ) then
        a(i,j) = a2(i)
      else if ( j == i+1 ) then
        a(i,j) = a3(i)
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do

  return
end
subroutine s3_vxm ( n, a1, a2, a3, x, b )
!
!*******************************************************************************
!
!! S3_VXM multiplies the transpose of a tridiagonal matrix times a vector.
!
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, real A1(2:N), A2(1:N), A3(1:N-1),
!    the nonzero diagonals of the linear system.
!
!    Input, real X(N), the vector to be multiplied by Transpose ( A ).
!
!    Output, real B(N), the product Transpose ( A ) * x.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer ierror
  real x(n)
!
!  Check dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions for the matrix.'
    return
  end if

  b(1:n) = a2(1:n) * x(1:n)
  b(2:n) = b(2:n) + a3(1:n-1) * x(1:n-1)
  b(1:n-1) = b(1:n-1) + a1(2:n) * x(2:n)

  return
end
subroutine s3_zero ( n, a1, a2, a3 )
!
!*******************************************************************************
!
!! S3_ZERO zeroes out a general tridiagonal matrix.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Output, real A1(2:N), A2(1:N), A3(1:N-1), the diagonals of the
!    matrix.
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  integer i
  integer ierror
!
!  Check the dimensions.
!
  call s3_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions for the matrix.'
    return
  end if

  a1(2:n) = 0.0E+00
  a2(1:n) = 0.0E+00
  a3(1:n-1) = 0.0E+00

  return
end
subroutine s3p_check ( n, ierror )
!
!*******************************************************************************
!
!! S3P_CHECK checks the dimensions of a tridiagonal periodic matrix.
!
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Output, integer IERROR, error flag.
!    0, the dimensions are legal.
!    1, N is less than 3.
!
  integer ierror
  integer n
!
  ierror = 0

  if ( n < 3 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'S3P_CHECK - Fatal error!'
    write ( *, * ) '  N must be at least 3.'
    write ( *, * ) '  The input value is N = ', n
  end if

  return
end
subroutine s3p_det ( n, a2, work4, det )
!
!*******************************************************************************
!
!! S3P_DET computes the determinant of a matrix factored by S3P_FA.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A2(N), factored matrix data from S3P_FA.
!
!    Input, real WORK4, factorization information from S3P_FA.
!
!    Output, real DET, the determinant of the matrix.
!
  integer n
!
  real a2(n)
  real det
  integer i
  integer ierror
  real work4
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  det = product ( a2(1:n-1) ) * work4

  return
end
subroutine s3p_fa ( n, a1, a2, a3, info, work2, work3, work4 )
!
!*******************************************************************************
!
!! S3P_FA factors a tridiagonal periodic matrix.
!
!
!  Discussion:
!
!    Once the matrix has been factored by S3P_FA, S3P_SL may be called
!    to solve linear systems involving the matrix.
!
!    The logical matrix has a form which is suggested by this diagram:
!
!      D1 U1          L1
!      L2 D2 U2
!         L3 D3 U3
!            L4 D4 U4
!               L5 D5 U5
!      U6          L6 D6
!
!    The algorithm treats the matrix as a border banded matrix:
!
!      ( A1  A2 )
!      ( A3  A4 )
!
!    where:
!
!      D1 U1          | L1
!      L2 D2 U2       |  0
!         L3 D3 U3    |  0
!            L4 D4 U4 |  0
!               L5 D5 | U5
!      ---------------+---
!      U6  0  0  0 L6 | D6
!
!  Method:
!
!    The algorithm rewrites the system as:
!
!         X1 + inverse(A1) A2 X2 = inverse(A1) B1
!
!      A3 X1 +             A4 X2 = B2
!
!    The first equation can be "solved" for X1 in terms of X2:
!
!         X1 = - inverse(A1) A2 X2 + inverse(A1) B1
!
!    allowing us to rewrite the second equation for X2 explicitly:
!
!      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input/output, real A1(N), A2(N), A3(N).
!    On input, these arrays contain the subdiagonal, diagonal, and
!    superdiagonal entries of the coefficient matrix.  The special
!    cases are that A1(1) is the coefficient of X(N), and A3(N)
!    is the coefficient of X(1).
!
!    On output, the arrays have been modified to hold information
!    defining the border-banded factorization of submatrices A1
!    and A3.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
!    Output, real WORK2(N-1), WORK3(N-1), WORK4, factorization information.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  integer i
  integer ierror
  integer info
  integer job
  real work2(n-1)
  real work3(n-1)
  real work4
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Compute inverse(A1):
!
  call s3_np_fa ( n-1, a1(2), a2, a3, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_FA - Fatal error!'
    write ( *, * ) '  S3_NP_FA returned INFO = ', info
    write ( *, * ) '  Factoring failed for column INFO.'
    write ( *, * ) '  The tridiagonal matrix A1 is singular.'
    write ( *, * ) '  This algorithm cannot continue!'
    write ( *, * ) ' '
    return
  end if
!
!  WORK2 := inverse(A1) * A2.
!
  work2(1) = a1(1)
  work2(2:n-2) = 0.0E+00
  work2(n-1) = a3(n-1)

  job = 0
  call s3_np_sl ( n-1, a1(2), a2, a3, work2, job )
!
!  WORK3 := inverse ( transpose ( A1 ) ) * tranpose ( A3 ).
!
  work3(1) = a3(n)
  work3(2:n-2) = 0.0E+00
  work3(n-1) = a1(n)

  job = 1
  call s3_np_sl ( n-1, a1(2), a2, a3, work3, job )
!
!  A4 := ( A4 - A3 * inverse(A1) * A2 )
!
  work4 = a2(n) - a3(n) * work2(1) - a1(n) * work2(n-1)

  if ( work4 == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'S3P_FA - Fatal error!'
    write ( *, * ) '  The factored A4 submatrix is zero.'
    write ( *, * ) '  This algorithm cannot continue!'
    write ( *, * ) ' '
    return
  end if

  return
end
subroutine s3p_ml ( n, a1, a2, a3, x, b, job )
!
!*******************************************************************************
!
!! S3P_ML computes A * x or x * A, where A has been factored by S3P_FA.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A1(N), A2(N), A3(N), the factors computed by S3P_FA.
!
!    Input, real X(N), the vector to be multiplied by the matrix.
!
!    Output, real B(N), the result of the multiplication.
!
!    Input, integer JOB, indicates what product should be computed.
!    0, compute A * x.
!    nonzero, compute transpose ( A ) * x.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  real b(n)
  integer ierror
  integer job
  real x(n)
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Multiply A(1:N-1,1:N-1) and X(1:N-1).
!
  call s3_np_ml ( n-1, a1(2), a2, a3, x, b, job )
!
!  Add terms from the border.
!
  if ( job == 0 ) then
    b(1) = b(1) + a1(1) * x(n)
    b(n-1) = b(n-1) + a3(n-1) * x(n)
    b(n) = a3(n) * x(1) + a1(n) * x(n-1) + a2(n) * x(n)
  else
    b(1) = b(1) + a3(n) * x(n)
    b(n-1) = b(n-1) + a1(n) * x(n)
    b(n) = a1(1) * x(1) + a3(n-1) * x(n-1) + a2(n) * x(n)
  end if

  return
end
subroutine s3p_mxv ( n, a1, a2, a3, x, b )
!
!*******************************************************************************
!
!! S3P_MXV computes A * x, where A is a tridiagonal periodic matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A1(N), A2(N), A3(N), the subdiagonal, diagonal,
!    and superdiagonal of the matrix.
!    A1(1) is actually the LAST element of the first row, and
!    A3(N) is the FIRST element of the last row.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  real b(n)
  integer i
  integer ierror
  real x(n)
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  b(1) = a1(1) * x(n) + a2(1) * x(1) + a3(1) * x(2)

  do i = 2, n-1
    b(i) = a1(i) * x(i-1) + a2(i) * x(i) + a3(i) * x(i+1)
  end do

  b(n) = a1(n) * x(n-1) + a2(n) * x(n) + a3(n) * x(1)

  return
end
subroutine s3p_print ( n, a1, a2, a3, title )
!
!*******************************************************************************
!
!! S3P_PRINT prints a periodic tridiagonal matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A1(N), A2(N), A3(N), the nonzero "diagonals" of the matrix.
!
!    Input, character ( len = * ) TITLE, a title to print.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call s3p_print_some ( n, a1, a2, a3, 1, 1, n, n )

  return
end
subroutine s3p_print_some ( n, a1, a2, a3, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! S3P_PRINT_SOME prints some of a periodic tridiagonal matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A1(N), A2(N), A3(N), the nonzero "diagonals" of the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column, to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    if ( i2lo > 1 .or. j2hi < n ) then
      i2lo = max ( i2lo, j2lo - 1 )
    end if

    i2hi = min ( ihi, n )

    if ( i2hi < n .or. j2lo > 1 ) then
      i2hi = min ( i2hi, j2hi + 1 )
    end if

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i == 1 .and. j == n ) then
          if ( r_is_int ( a1(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a1(1)
          else
            write ( ctemp(j2), '(g14.6)' ) a1(1)
          end if
        else if ( i == n .and. j == 1 ) then
          if ( r_is_int ( a3(n) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a3(n)
          else
            write ( ctemp(j2), '(g14.6)' ) a3(n)
          end if
        else if ( i-j > 1 .or. j-i > 1 ) then
          ctemp(j2) = '              '
        else if ( j == i-1 ) then
          if ( r_is_int ( a1(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a1(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a1(i)
          end if
        else if ( j == i ) then
          if ( r_is_int ( a2(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a2(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a2(i)
          end if
        else if ( j == i+1 ) then
          if ( r_is_int ( a3(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a3(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a3(i)
          end if
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine s3p_random ( n, a1, a2, a3 )
!
!*******************************************************************************
!
!! S3P_RANDOM randomizes a tridiagonal periodic matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Output, real A1(N), A2(N), A3(N), the subdiagonal, diagonal,
!    and superdiagonal of the matrix.
!    A1(1) is actually the LAST element of the first row, and
!    A3(N) is the FIRST element of the last row.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  integer ierror
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  call rvec_random ( 0.0E+00, 1.0E+00, n, a1 )
  call rvec_random ( 0.0E+00, 1.0E+00, n, a2 )
  call rvec_random ( 0.0E+00, 1.0E+00, n, a3 )

  return
end
subroutine s3p_sl ( n, a1, a2, a3, b, x, job, work2, work3, work4 )
!
!*******************************************************************************
!
!! S3P_SL solves a tridiagonal periodic system factored by S3P_FA.
!
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A1(N), A2(N), A3(N), factor data from S3P_FA.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution to the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
!    Input, real WORK2(N-1), WORK3(N-1), WORK4, factor data from S3P_FA.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  real b(n)
  integer i
  integer ierror
  integer job
  real work2(n-1)
  real work3(n-1)
  real work4
  real x(n)
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  x(1:n) = b(1:n)

  if ( job == 0 ) then
!
!  Solve A1 * X1 = B1.
!
    call s3_np_sl ( n-1, a1(2), a2, a3, x, job )
!
!  X2 = B2 - A3 * X1
!
    x(n) = x(n) - a3(n) * x(1) - a1(n) * x(n-1)
!
!  Solve A4 * X2 = X2
!
    x(n) = x(n) / work4
!
!  X1 := X1 - inverse ( A1 ) * A2 * X2.
!
    x(1:n-1) = x(1:n-1) - work2(1:n-1) * x(n)

  else
!
!  Solve transpose ( A1 ) * X1 = B1.
!
    call s3_np_sl ( n-1, a1(2), a2, a3, x, job )
!
!  X2 := X2 - transpose ( A2 ) * B1
!
    x(n) = x(n) - a1(1) * x(1) - a3(n-1) * x(n-1)
!
!  Solve A4 * X2 = X2.
!
    x(n) = x(n) / work4
!
!  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
!
    x(1:n-1) = x(1:n-1) - work3(1:n-1) * x(n)

  end if

  return
end
subroutine s3p_to_sge ( lda, n, a1, a2, a3, a )
!
!*******************************************************************************
!
!! S3P_TO_SGE copies a tridiagonal periodic matrix into a general matrix.
!
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A1(N), A2(N), A3(N), the periodic tridiagonal matrix.
!
!    Output, real A(LDA,N), the periodic tridiagonal matrix, stored as
!    a general matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real a1(n)
  real a2(n)
  real a3(n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for tridiagonal periodic matrix.'
    return
  end if

  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = a2(i)
      else if ( j == i-1 ) then
        a(i,j) = a1(i)
      else if ( j == i+1 ) then
        a(i,j) = a3(i)
      else if ( i == 1 .and. j == n ) then
        a(i,j) = a1(1)
      else if ( i == n .and. j == 1 ) then
        a(i,j) = a3(n)
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine s3p_vxm ( n, a1, a2, a3, x, b )
!
!*******************************************************************************
!
!! S3P_VXM computes X*A, where A is a tridiagonal periodic matrix.
!
!
!  Modified:
!
!    27 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A1(N), A2(N), A3(N), the subdiagonal, diagonal,
!    and superdiagonal of the matrix.
!    A1(1) is actually the LAST element of the first row, and
!    A3(N) is the FIRST element of the last row.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product X*A.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  real b(n)
  integer i
  integer ierror
  real x(n)
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  b(1) = a3(n) * x(n) + a2(1) * x(1) + a1(2) * x(2)

  do i = 2, n-1
    b(i) = a3(i-1) * x(i-1) + a2(i) * x(i) + a1(i+1) * x(i+1)
  end do

  b(n) = a3(n-1) * x(n-1) + a2(n) * x(n) + a1(1) * x(1)

  return
end
subroutine s3p_zero ( n, a1, a2, a3 )
!
!*******************************************************************************
!
!! S3P_ZERO zeroes out a tridiagonal periodic matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Output, real A1(N), A2(N), A3(N), the subdiagonal, diagonal,
!    and superdiagonal of the matrix.
!    A1(1) is actually the LAST element of the first row, and
!    A3(N) is the FIRST element of the last row.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  integer i
  integer ierror
!
!  Check the dimensions.
!
  call s3p_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S3P_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  a1(1:n) = 0.0E+00
  a2(1:n) = 0.0E+00
  a3(1:n) = 0.0E+00

  return
end
subroutine s5_check ( n, ierror )
!
!*******************************************************************************
!
!! S5_CHECK checks the dimensions of a pentadiagonal matrix.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Output, integer IERROR, error flag.
!    0, no errors detected.
!    1, N was less than 3.
!
  integer ierror
  integer n
!
  ierror = 0

  if ( n < 3 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'S5_CHECK - Fatal error!'
    write ( *, * ) '  N must be at least 3.'
    write ( *, * ) '  The input N was ', n
  end if

  return
end
subroutine s5_fs ( n, a1, a2, a3, a4, a5, b, x )
!
!*******************************************************************************
!
!! S5_FS factors and solves a pentadiagonal linear system.
!
!
!  Note:
!
!    This algorithm requires that each diagonal entry be nonzero.
!
!  Modified:
!
!    05 December 1998
!
!  Reference:
!
!    Cheney and Kincaid,
!    Numerical Mathematics and Computing,
!    1985, pages 233-236.
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A5(1:N-2).
!    On input, the nonzero diagonals of the linear system.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input/output, real B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B has been overwritten by factorization information.
!
!    Output, real X(N), the solution of the linear system.
!
  integer n
!
  real a1(3:n)
  real a2(2:n)
  real a3(1:n)
  real a4(1:n-1)
  real a5(1:n-2)
  real b(n)
  integer i
  integer ierror
  real x(n)
  real xmult
!
!  Check the dimensions.
!
  call s5_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_FS - Fatal error!'
    write ( *, * ) '  Illegal dimensions for pentadiagonal matrix.'
    return
  end if

  do i = 1, n
    if ( a3(i) == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'S5_FS - Fatal error!'
      write ( *, * ) '  A3(', i, ') = 0.'
      return
    end if
  end do

  do i = 2, n-1

    xmult = a2(i) / a3(i-1)
    a3(i) = a3(i) - xmult * a4(i-1)
    a4(i) = a4(i) - xmult * a5(i-1)

    b(i) = b(i) - xmult * b(i-1)

    xmult = a1(i+1) / a3(i-1)
    a2(i+1) = a2(i+1) - xmult * a4(i-1)
    a3(i+1) = a3(i+1) - xmult * a5(i-1)

    b(i+1) = b(i+1) - xmult * b(i-1)

  end do

  xmult = a2(n) / a3(n-1)
  a3(n) = a3(n) - xmult * a4(n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a3(n)
  x(n-1) = ( b(n-1) - a4(n-1) * x(n) ) / a3(n-1)

  do i = n-2, 1, -1
    x(i) = ( b(i) - a4(i) * x(i+1) - a5(i) * x(i+2) ) / a3(i)
  end do

  return
end
subroutine s5_mxv ( n, a1, a2, a3, a4, a5, x, b )
!
!*******************************************************************************
!
!! S5_MXV multiplies a pentadiagonal matrix times a vector.
!
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A5(1:N-2),
!    the nonzero diagonals of the linear system.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer n
!
  real a1(3:n)
  real a2(2:n)
  real a3(1:n)
  real a4(1:n-1)
  real a5(1:n-2)
  real b(n)
  integer ierror
  real x(n)
!
!  Check the dimensions.
!
  call s5_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions for pentadiagonal matrix.'
    return
  end if

  b(1:n) = a3(1:n) * x(1:n)

  b(3:n) = b(3:n) + a1(3:n) * x(1:n-2)
  b(2:n) = b(2:n) + a2(2:n) * x(1:n-1)
  b(1:n-1) = b(1:n-1) + a4(1:n-1) * x(2:n)
  b(1:n-2) = b(1:n-2) + a5(1:n-2) * x(3:n)

  return
end
subroutine s5_print ( n, a1, a2, a3, a4, a5, title )
!
!*******************************************************************************
!
!! S5_PRINT_SOME prints some of a pentadiagonal matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A4(1:N-2), the
!    nonzero diagonals of the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  real a1(3:n)
  real a2(2:n)
  real a3(n)
  real a4(1:n-1)
  real a5(1:n-2)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call s5_print_some ( n, a1, a2, a3, a4, a5, 1, 1, n, n )

  return
end
subroutine s5_print_some ( n, a1, a2, a3, a4, a5, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! S5_PRINT_SOME prints some of a pentadiagonal matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A4(1:N-2), the
!    nonzero diagonals of the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column, to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a1(3:n)
  real a2(2:n)
  real a3(n)
  real a4(1:n-1)
  real a5(1:n-2)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Check the dimensions.
!
  call s5_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 2 )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i-j > 2 .or. j-i > 2 ) then
          ctemp(j2) = '              '
        else if ( j == i-2 ) then
          if ( r_is_int ( a1(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a1(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a1(i)
          end if
        else if ( j == i-1 ) then
          if ( r_is_int ( a2(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a2(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a2(i)
          end if
        else if ( j == i ) then
          if ( r_is_int ( a3(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a3(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a3(i)
          end if
        else if ( j == i+1 ) then
          if ( r_is_int ( a4(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a4(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a4(i)
          end if
        else if ( j == i+2 ) then
          if ( r_is_int ( a5(i) ) ) then
            write ( ctemp(j2), '(f8.0,6x)' ) a5(i)
          else
            write ( ctemp(j2), '(g14.6)' ) a5(i)
          end if
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine s5_random ( n, a1, a2, a3, a4, a5 )
!
!*******************************************************************************
!
!! S5_RANDOM returns a random pentadiagonal matrix.
!
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Output, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A5(1:N-2),
!    the nonzero diagonals of the linear system.  The entries
!    are all between 0 and 1.
!
  integer n
!
  real a1(3:n)
  real a2(2:n)
  real a3(1:n)
  real a4(1:n-1)
  real a5(1:n-2)
  integer i
  integer ierror
!
!  Check the dimensions.
!
  call s5_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions for pentadiagonal matrix.'
    return
  end if

  call rvec_random ( 0.0E+00, 1.0E+00, n-2, a1(3:n) )
  call rvec_random ( 0.0E+00, 1.0E+00, n-1, a2(2:n) )
  call rvec_random ( 0.0E+00, 1.0E+00, n  , a3(1:n) )
  call rvec_random ( 0.0E+00, 1.0E+00, n-1, a4(1:n-1) )
  call rvec_random ( 0.0E+00, 1.0E+00, n-2, a5(1:n-2) )

  return
end
subroutine s5_to_sge ( lda, n, a1, a2, a3, a4, a5, a )
!
!*******************************************************************************
!
!! S5_TO_SGE copies a pentadiagonal matrix into a general matrix.
!
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A5(1:N-2).
!    the nonzero diagonals of the matrix.
!
!    Output, real A(LDA,N), the matrix, stored as a general matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real a1(3:n)
  real a2(2:n)
  real a3(1:n)
  real a4(1:n-1)
  real a5(1:n-2)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call s5_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for pentadiagonal matrix.'
    return
  end if

  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, n

      if ( j == i-2 ) then
        a(i,j) = a1(i)
      else if ( j == i-1 ) then
        a(i,j) = a2(i)
      else if ( i == j ) then
        a(i,j) = a3(i)
      else if ( j == i+1 ) then
        a(i,j) = a4(i)
      else if ( j == i+2 ) then
        a(i,j) = a5(i)
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do

  return
end
subroutine s5_vxm ( n, a1, a2, a3, a4, a5, x, b )
!
!*******************************************************************************
!
!! S5_VXM multiplies the transpose of a pentadiagonal matrix times a vector.
!
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A5(1:N-2),
!    the nonzero diagonals of the linear system.
!
!    Input, real X(N), the vector to be multiplied by Transpose ( A ).
!
!    Output, real B(N), the product Transpose ( A ) * x.
!
  integer n
!
  real a1(3:n)
  real a2(2:n)
  real a3(1:n)
  real a4(1:n-1)
  real a5(1:n-2)
  real b(n)
  integer ierror
  real x(n)
!
!  Check the dimensions.
!
  call s5_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions for pentadiagonal matrix.'
    return
  end if

  b(1:n) = a3(1:n) * x(1:n)
  b(2:n) = b(2:n) + a4(1:n-1) * x(1:n-1)
  b(3:n) = b(3:n) + a5(1:n-2) * x(1:n-2)
  b(1:n-1) = b(1:n-1) + a2(2:n) * x(2:n)
  b(1:n-2) = b(1:n-2) + a1(3:n) * x(3:n)

  return
end
subroutine s5_zero ( n, a1, a2, a3, a4, a5 )
!
!*******************************************************************************
!
!! S5_ZERO zeroes a pentadiagonal matrix.
!
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real A1(3:N), A2(2:N), A3(1:N), A4(1:N-1), A5(1:N-2),
!    the diagonals of the linear system.
!
!    Input, integer N, the order of the linear system.
!
  integer n
!
  real a1(3:n)
  real a2(2:n)
  real a3(1:n)
  real a4(1:n-1)
  real a5(1:n-2)
  integer i
  integer ierror
!
!  Check the dimensions.
!
  call s5_check ( n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S5_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions for pentadiagonal matrix.'
    return
  end if

  a1(3:n) = 0.0E+00
  a2(2:n) = 0.0E+00
  a3(1:n) = 0.0E+00
  a4(1:n-1) = 0.0E+00
  a5(1:n-2) = 0.0E+00

  return
end
subroutine sbb_add ( n1, n2, ml, mu, a, i, j, value )
!
!*******************************************************************************
!
!! SBB_ADD adds a value to an entry in a border banded matrix.
!
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the border banded matrix.
!
!    Input, integer I, J, the row and column of the entry to be incremented.
!    Some combinations of I and J are illegal.
!
!    Input, real VALUE, the value to be added to the (I,J)-th entry.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer j
  real value
!
  if ( value == 0.0E+00 ) then
    return
  end if
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_ADD - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  if ( i <= 0 .or. i > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_ADD - Fatal error!'
    write ( *, * ) '  Illegal input value of row index I = ',i
    stop
  end if

  if ( j <= 0 .or. j > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_ADD - Fatal error!'
    write ( *, * ) '  Illegal input value of row index J = ',j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition (J-I) > MU, but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes (J-I) > MU+ML.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( (j-i) > mu+ml .or. (i-j) > ml ) then
      write ( *, * ) ' '
      write ( *, * ) 'SBB_ADD - warning!'
      write ( *, * ) '  Unable to add to entry (', i, ',', j, ').'
    else
      ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. j > n1 ) then
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( i > n1 ) then
    ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1)
  end if

  a(ij) = a(ij) + value

  return
end
subroutine sbb_check ( n1, n2, ml, mu, ierror )
!
!*******************************************************************************
!
!! SBB_CHECK checks the dimensions of a border banded matrix.
!
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one of
!    N1 and N2 must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1 - 1.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if ML is illegal;
!    IERROR = IERROR + 2 if MU is illegal;
!    IERROR = IERROR + 4 if N1 is illegal;
!    IERROR = IERROR + 8 if N2 is illegal;
!    IERROR = IERROR + 16 if neither N1 nor N2 is positive.
!
  integer ierror
  integer ml
  integer mu
  integer n1
  integer n2
!
  ierror = 0

  if ( ml < 0 .or. ml > max ( ( n1 - 1 ), 0 ) ) then
    ierror = 1
    write ( *, * ) ' '
    write ( *, * ) 'SBB_CHECK - Illegal ML = ', ml
  end if

  if ( mu < 0 .or. mu > max ( ( n1 - 1 ), 0 ) ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SBB_CHECK - Illegal MU = ', mu
  end if

  if ( n1 < 0 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SBB_CHECK - Illegal N1 = ', n1
  end if

  if ( n2 < 0 ) then
    ierror = ierror + 8
    write ( *, * ) ' '
    write ( *, * ) 'SBB_CHECK - Illegal N2 = ', n2
  end if

  if ( n1 + n2 <= 0 ) then
    ierror = ierror + 16
    write ( *, * ) ' '
    write ( *, * ) 'SBB_CHECK - Illegal N1+N2 = ', n1+n2
  end if

  return
end
subroutine sbb_fa ( n1, n2, ml, mu, a, pivot, info )
!
!*******************************************************************************
!
!! SBB_FA factors a border banded matrix.
!
!
!  Discussion:
!
!    Once the matrix has been factored by SBB_FA, SBB_SL may be called
!    to solve linear systems involving the matrix.
!
!    SBB_FA uses LINPACK routines to carry out the factorization.
!
!
!    The linear system must be border banded, of the form:
!
!      ( A1 A2 ) (X1) = (B1)
!      ( A3 A4 ) (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    The algorithm rewrites the system as:
!
!         X1 + inv(A1) A2 X2 = inv(A1) B1
!
!      A3 X1 +         A4 X2 = B2
!
!    and then rewrites the second equation as
!
!      ( A4 - A3 inv(A1) A2 ) X2 = B2 - A3 inv(A1) B1
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input/output, real A( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2 ).
!
!    On input, A contains the border-banded matrix to be factored.
!
!    In particular, A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    On output, A contains information describing a partial factorization
!    of the original coefficient matrix.  This information is required
!    by SBB_SL in order to solve linear systems associated with that
!    matrix.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= n1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    Output, integer PIVOT(N1+N2), contains pivoting information.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer ik
  integer info
  integer pivot(n1+n2)
  integer j
  integer jk
  integer job
  integer k
  integer nband
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  nband = (2*ml+mu+1) * n1
!
!  Factor the A1 band matrix, overwriting A1 by its factors.
!
  if ( n1 > 0 ) then

    call sgb_fa ( 2*ml+mu+1, n1, ml, mu, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SBB_FA - Fatal error!'
      write ( *, * ) '  SGB_FA returned INFO = ', info
      write ( *, * ) '  Factoring failed for column INFO.'
      write ( *, * ) '  The band matrix A1 is singular.'
      write ( *, * ) '  This algorithm cannot continue!'
      write ( *, * ) ' '
      return
    end if

  end if

  if ( n1 > 0 .and. n2 > 0 ) then
!
!  Solve A1 * x = -A2 for x, and overwrite A2 by the results.
!
    do i = nband+1, nband+n1*n2
      a(i) = - a(i)
    end do

    job = 0
    do i = 1, n2
      call sgb_sl ( 2*ml+mu+1, n1, ml, mu, a, pivot, a(nband+(i-1)*n1+1), job )
    end do
!
!  A4 := A4 + A3 * A2.
!
    do i = 1, n2
      do j = 1, n1
        ij = nband + n1*n2 + (j-1)*n2 + i
        do k = 1, n2
          ik = nband + 2*n1*n2 + (k-1)*n2 + i
          jk = nband + (k-1)*n1 + j
          a(ik) = a(ik) + a(ij) * a(jk)
        end do
      end do
    end do

  end if
!
!  Factor A4.
!
  if ( n2 > 0 ) then

    call sge_fa ( n2, n2, a(nband+2*n1*n2+1), pivot(n1+1), info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SBB_FA - Fatal error!'
      write ( *, * ) '  SGE_FA returned INFO = ',info
      write ( *, * ) '  This indicates singularity in column INFO.'
      write ( *, * ) '  of the A4 submatrix, which is column ', n1+info
      write ( *, * ) '  of the full matrix.'
      write ( *, * ) ' '
      write ( *, * ) '  It is possible that the full matrix is '
      write ( *, * ) '  nonsingular, but the algorithm SBB_FA may'
      write ( *, * ) '  not be used for this matrix.'
      return
    end if
  end if

  return
end
subroutine sbb_get ( n1, n2, ml, mu, a, i, j, value )
!
!*******************************************************************************
!
!! SBB_GET returns an entry of a border banded matrix.
!
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the border banded matrix.
!
!    Input, integer I, J, the row and column of the entry to be retrieved.
!
!    Output, real VALUE, the value of the (I,J) entry.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer j
  real value
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_GET - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  if ( i <= 0 .or. i > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_GET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index I = ',i
    stop
  end if

  if ( j <= 0 .or. j > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_GET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index J = ',j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition (J-I) > MU, but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes (J-I) > MU+ML.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( (j-i) > mu+ml .or. (i-j) > ml ) then
      value = 0.0E+00
      return
    else
      ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. j > n1 ) then
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( i > n1 ) then
    ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
  end if

  value = a(ij)

  return
end
subroutine sbb_mxv ( n1, n2, ml, mu, a, x, y )
!
!*******************************************************************************
!
!! SBB_MXV multiplies a border banded matrix times a vector.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the border banded matrix.
!
!    Input, real X(N1+N2), the vector to be multiplied by A.
!
!    Output, real Y(N1+N2), the result of multiplying A by X.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ihi
  integer ij
  integer ilo
  integer j
  real x(n1+n2)
  real y(n1+n2)
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Initialize Y.
!
  y(1:n1+n2) = 0.0E+00
!
!  Multiply by A1.
!
  do j = 1, n1

    ilo = max ( 1, j - mu - ml )
    ihi = min ( n1, j + ml )
    ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1

    y(ilo:ihi) = y(ilo:ihi) + a(ij+ilo:ij+ihi) * x(j)

  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1

    y(1:n1) = y(1:n1) + a(ij+1:ij+n1) * x(j)

  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (2*ml+mu+1)*n1+n1*n2+(j-1)*n2-n1

    y(n1+1:n1+n2) = y(n1+1:n1+n2) + a(ij+n1+1:ij+n1+n2) * x(j)

  end do

  return
end
subroutine sbb_print ( n1, n2, ml, mu, a, title )
!
!*******************************************************************************
!
!! SBB_PRINT prints a border banded matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the N1+N2 by N1+N2
!    border banded matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sbb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2 )

  return
end
subroutine sbb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SBB_PRINT_SOME prints some of a border banded matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the N1+N2 by N1+N2
!    border banded matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ij
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n1+n2 )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n1+n2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0E+00

        if ( i <= n1 .and. j <= n1 ) then
          if ( (j-i) <= mu+ml .and. (i-j) <= ml ) then
            ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
            aij = a(ij)
          end if
        else if ( i <= n1 .and. j > n1 ) then
          ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
          aij = a(ij)
        else if ( i > n1 ) then
          ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
          aij = a(ij)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sbb_random ( n1, n2, ml, mu, a )
!
!*******************************************************************************
!
!! SBB_RANDOM randomizes a border banded matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Output, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the border
!    banded matrix.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ihi
  integer ilo
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if
!
!  Randomize the banded matrix A1.
!
  if ( n1 > 0 ) then
    call sgb_random ( 2*ml+mu+1, n1, n1, ml, mu, a )
  end if
!
!  Randomize the rectangular strip A2.
!
  ilo = (2*ml+mu+1) * n1 + 1
  ihi = (2*ml+mu+1) * n1 + n1*n2
  call rvec_random ( 0.0E+00, 1.0E+00, ihi+1-ilo, a(ilo:ihi) )
!
!  Randomize the rectangular strip A3.
!
  ilo = (2*ml+mu+1) * n1 +   n1*n2 + 1
  ihi = (2*ml+mu+1) * n1 + 2*n1*n2

  call rvec_random ( 0.0E+00, 1.0E+00, ihi+1-ilo, a(ilo:ihi) )
!
!  Randomize the square matrix A4.
!
  ilo = (2*ml+mu+1) * n1 + 2*n1*n2 + 1
  ihi = (2*ml+mu+1) * n1 + 2*n1*n2 + n2**2

  call rvec_random ( 0.0E+00, 1.0E+00, ihi+1-ilo, a(ilo:ihi) )

  return
end
subroutine sbb_set ( n1, n2, ml, mu, a, i, j, value )
!
!*******************************************************************************
!
!! SBB_SET sets an entry of a border banded matrix.
!
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the border banded matrix.
!
!    Input, integer I, J, the row and column of the entry to be set.
!
!    Input, real VALUE, the value to be assigned to the (I,J) entry.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer j
  real value
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_SET - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  if ( i <= 0 .or. i > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_SET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index I = ',i
    stop
  end if

  if ( j <= 0 .or. j > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_SET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index J = ',j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition (J-I) > MU, but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes (J-I) > MU+ML.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( (j-i) > mu+ml .or. (i-j) > ml ) then
      write ( *, * ) ' '
      write ( *, * ) 'SBB_SET - warning!'
      write ( *, * ) '  Unable to set entry (', i, ',', j, ').'
      return
    else
      ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. j > n1 ) then
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( i > n1 ) then
    ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
  end if

  a(ij) = value

  return
end
subroutine sbb_sl ( n1, n2, ml, mu, a, pivot, b )
!
!*******************************************************************************
!
!! SBB_SL solves a border banded linear system factored by SBB_FA.
!
!
!  Discussion:
!
!    The linear system A * x = b is decomposable into the block system:
!
!      ( A1 A2 ) * (X1) = (B1)
!      ( A3 A4 )   (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    All the arguments except B are input quantities only, which are
!    not changed by the routine.  They should have exactly the same values
!    they had on exit from SBB_FA.
!
!    If more than one right hand side is to be solved, with the same matrix,
!    SBB_SL should be called repeatedly.  However, SBB_FA only needs to be called
!    once to create the factorization.
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real A( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
!    the LU factors computed by SBB_FA.
!
!    Input, integer PIVOT(N1+N2), the pivoting information from SBB_FA.
!
!    Input/output, real B(N1+N2).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real b(n1+n2)
  integer i
  integer ierror
  integer ij
  integer pivot(n1+n2)
  integer j
  integer job
  integer nband
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  nband = (2*ml+mu+1)*n1
!
!  Set B1 := inverse(A1) * B1.
!
  if ( n1 > 0 ) then
    job = 0
    call sgb_sl ( 2*ml+mu+1, n1, ml, mu, a, pivot, b, job )
  end if
!
!  Modify the right hand side of the second linear subsystem.
!  Set B2 := B2 - A3*B1.
!
  do i = 1, n2
    do j = 1, n1
      ij = nband + n1*n2 + (j-1)*n2 + i
      b(n1+i) = b(n1+i) - a(ij) * b(j)
    end do
  end do
!
!  Set B2 := inverse(A4) * B2.
!
  if ( n2 > 0 ) then
    job = 0
    call sge_sl ( n2, n2, a(nband+2*n1*n2+1), pivot(n1+1), b(n1+1), job )
  end if
!
!  Modify the first subsolution.
!  Set B1 := B1 + A2*B2.
!
  do i = 1, n1
    do j = 1, n2
      ij = nband + (j-1)*n1 + i
      b(i) = b(i) + a(ij) * b(n1+j)
    end do
  end do

  return
end
subroutine sbb_to_sge ( lda, n1, n2, ml, mu, a, a2 )
!
!*******************************************************************************
!
!! SBB_TO_SGE copies a border banded matrix into a general matrix.
!
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A2.
!    LDA must be at least N1+N2.
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the border banded matrix.
!
!    Output, real A2(LDA,N1+N2), a copy of the matrix, in general storage.
!
  integer lda
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real a2(lda,n1+n2)
  integer i
  integer ierror
  integer ij
  integer j
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for border banded matrix!'
    return
  end if

  call sge_check ( lda, n1+n2, n1+n2, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions general matrix!'
    return
  end if

  do i = 1, n1
    do j = 1, n1

      if ( (j-i) > mu+ml .or. (i-j) > ml ) then
        a2(i,j) = 0.0E+00
      else
        ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
        a2(i,j) = a(ij)
      end if

    end do
  end do

  do i = 1, n1
    do j = n1+1, n2
      ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
      a2(i,j) = a(ij)
    end do
  end do

  do i = n1+1, n2
    do j = 1, n1+n2
      ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
      a2(i,j) = a(ij)
    end do
  end do

  return
end
subroutine sbb_vxm ( n1, n2, ml, mu, a, x, y )
!
!*******************************************************************************
!
!! SBB_VXM multiplies a vector times a border banded matrix.
!
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real A((2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2), the border banded matrix.
!
!    Input, real X(N1+N2), the vector to multiply A.
!
!    Output, real Y(N1+N2), the product X times A.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ihi
  integer ij
  integer ilo
  integer j
  real x(n1+n2)
  real y(n1+n2)
!
!  Check the dimensions.
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Initialize Y.
!
  y(1:n1+n2) = 0.0E+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j - mu - ml )
    ihi = min ( n1, j + ml )
    ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1
    do i = ilo, ihi
      y(j) = y(j) + x(i) * a(ij+i)
    end do
  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (2*ml+mu+1)*n1+(j-n1-1)*n1
    do i = 1, n1
      y(j) = y(j) + x(i) * a(ij+i)
    end do
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (2*ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
    do i = n1+1, n1+n2
      y(j) = y(j) + x(i) * a(ij+i)
    end do
  end do

  return
end
subroutine sbb_zero ( n1, n2, ml, mu, a )
!
!*******************************************************************************
!
!! SBB_ZERO zeroes out a border banded matrix.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Output, real A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the border banded matrix.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SBB_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, (2*ml+mu+1)*n1 + 2*n1*n2 + n2*n2
    a(i) = 0.0E+00
  end do

  return
end
subroutine sbto_mxv ( m, l, a1, a2, x, b )
!
!***********************************************************************
!
!! SBTO_MXV computes the real block Toeplitz matrix product A * X = B.
!
!
!  Discussion:
!
!    The full matrix has order M * L, and can be regarded
!    as comprising L by L blocks.  Each block is of order M.
!
!    Example:
!
!      M = 2, L = 3
!
!      1 2 | 3 4 | 5 6
!      5 5 | 6 6 | 7 7
!      ----+-----+-----
!      7 8 | 1 2 | 3 4
!      8 8 | 5 5 | 6 6
!      ----+-----+-----
!      9 0 | 7 8 | 1 2
!      9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ 91, 134, 73, 125, 97, 129 /)
!
!  Modified:
!
!    20 March 2001
!
!  Parameters:
!
!    Input, integer M, the order of the blocks of the matrix A.
!
!    Input, integer L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, real X(M*L), the vector to be multiplied.
!
!    Output, real B(M*L), the product vector, A * X.
!
  integer l
  integer m
!
  real a1(m,m,l)
  real a2(m,m,l-1)
  real b(m,l)
  integer i
  integer j
  real x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0E+00

    do j = 1, i-1
      b(1:m,i) = b(1:m,i) + matmul ( a2(1:m,1:m,i-j), x(1:m,j) )
    end do

    do j = i, l
      b(1:m,i) = b(1:m,i) + matmul ( a1(1:m,1:m,j+1-i), x(1:m,j) )
    end do

  end do

  return
end
subroutine sbto_print ( m, l, a1, a2, title )
!
!*******************************************************************************
!
!! SBTO_PRINT prints a block Toeplitz matrix.
!
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the order of the blocks of the matrix A.
!
!    Input, integer L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer l
  integer m
!
  real a1(m,m,l)
  real a2(m,m,l-1)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sbto_print_some ( m, l, a1, a2, 1, 1, m*l, m*l )

  return
end
subroutine sbto_print_some ( m, l, a1, a2, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SBTO_PRINT_SOME prints some of a block Toeplitz matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the order of the blocks of the matrix A.
!
!    Input, integer L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer l
  integer m
!
  real a1(m,m,l)
  real a2(m,m,l-1)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i1
  integer i2
  integer i3hi
  integer i3lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j1
  integer j2
  integer j3
  integer j3hi
  integer j3lo
  integer jhi
  integer jlo
  integer n
  logical r_is_int
!
  n = m * l
!
!  Print the columns of the matrix, in strips of 5.
!
  do j3lo = jlo, jhi, incx

    j3hi = j3lo + incx - 1
    j3hi = min ( j3hi, n )
    j3hi = min ( j3hi, jhi )

    inc = j3hi + 1 - j3lo

    write ( *, * ) ' '

    do j = j3lo, j3hi
      j3 = j + 1 - j3lo
      write ( ctemp(j3), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j3), j3 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i3lo = max ( ilo, 1 )
    i3hi = min ( ihi, n )

    do i = i3lo, i3hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j3 = 1, inc

        j = j3lo - 1 + j3
!
!  i = M * ( i1 - 1 ) + i2
!  j = M * ( j1 - 1 ) + j2
!
        i1 = ( i - 1 ) / m + 1
        i2 = i - m * ( i1 - 1 )
        j1 = ( j - 1 ) / m + 1
        j2 = j - m * ( j1 - 1 )

        if ( j1 >= i1 ) then
          aij = a1(i2,j2,j1+1-i1)
        else
          aij = a2(i2,j2,i1-j1)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j3), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j3), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j3), j3 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sbto_sl ( m, l, a1, a2, b, x )
!
!***********************************************************************
!
!! SBTO_SL solves the real block Toeplitz linear system A * X = B.
!
!
!  Discussion:
!
!    The full matrix has order M * L, and can be regarded
!    as comprising L by L blocks.  Each block is of order M.
!
!    Example:
!
!      M = 2, L = 3
!
!      1 2 | 3 4 | 5 6
!      5 5 | 6 6 | 7 7
!      ----+-----+-----
!      7 8 | 1 2 | 3 4
!      8 8 | 5 5 | 6 6
!      ----+-----+-----
!      9 0 | 7 8 | 1 2
!      9 9 | 8 8 | 5 5
!
!  Modified:
!
!    07 March 2001
!
!  Parameters:
!
!    Input, integer M, the order of the blocks of the matrix A.
!
!    Input, integer L, the number of blocks in a row or column of A.
!
!    Input, real A1(M*M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  Each block is represented by columns.
!
!    Input, real A2(M*M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!    Each block is represented by columns.
!
!    Input, real B(M*L), the right hand side vector.
!
!    Output, real X(M*L), the solution vector.  X and B may share storage.
!
  integer l
  integer m
!
  real a1(m*m,l)
  real a2(m*m,l-1)
  real b(m,l)
  real c1(m,m,l-1)
  real c2(m,m,l-1)
  integer i
  integer i1
  integer i2
  integer i3
  integer info
  integer j
  integer n
  real r(m)
  real r1(m,m)
  real r2(m,m)
  real r3(m,m)
  real r5(m,m)
  real r6(m,m)
  real x(m,l)
!
!  Solve the system with the principal minor of order M.
!
  i3 = 1
  do j = 1, m
    do i = 1, m
      c1(i,j,1) = a1(i3,1)
      r1(i,j) = a1(i3,1)
      i3 = i3 + 1
    end do
  end do

  r3(1:m,1:m) = r1(1:m,1:m)
  x(1:m,1) = b(1:m,1)

  call sgefa ( r3, m, m, r, info )

  call sgesl ( r3, m, m, r, x(1,1), 0 )

  if ( l == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system
!  with the block Toeplitz matrix for N = 2 through L.
!
  do n = 2, l
!
!  Compute multiples of the first and last block columns of
!  the inverse of the principal minor of order M*N.
!
    i3 = 1
    do j = 1, m
      do i = 1, m
        r5(i,j) = a2(i3,n-1)
        r6(i,j) = a1(i3,n)
        i3 = i3 + 1
      end do
    end do

    if ( n /= 2 ) then

      c1(1:m,1:m,n-1) = r2(1:m,1:m)

      do i1 = 1, n-2
        i2 = n - i1
        do j = 1, m
          i3 = 1
          do i = 1, m
            call saxpy ( m, c1(i,j,i2), a2(i3,i1), 1, r5(1,j), 1 )
            call saxpy ( m, c2(i,j,i1), a1(i3,i1+1), 1, r6(1,j), 1 )
            i3 = i3 + m
          end do
        end do
      end do

    end if

    do j = 1, m
      r2(1:m,j) = -r5(1:m,j)
      call sgesl ( r3, m, m, r, r2(1,j), 0 )
    end do

    r3(1:m,1:m) = r6(1:m,1:m)
    r6(1:m,1:m) = -c1(1:m,1:m,1)

    do j = 1, m
      do i = 1, m
        call saxpy ( m, r2(i,j), r3(1,i), 1, c1(1,j,1), 1 )
      end do
    end do

    call sgefa ( r6, m, m, r, info )

    do j = 1, m
      call sgesl ( r6, m, m, r, r3(1,j), 0 )
      do i = 1, m
        call saxpy ( m, r3(i,j), r5(1,i), 1, r1(1,j), 1 )
      end do
    end do

    if ( n /= 2 ) then

      r6(1:m,1:m) = c2(1:m,1:m,1)

      do i1 = 2, n-1

        if ( i1 /= n-1 ) then
          r5(1:m,1:m) = c2(1:m,1:m,i1)
        end if

        do j = 1, m
          c2(1:m,j,i1) = r6(1:m,j)
          do i = 1, m
            call saxpy ( m, r3(i,j), c1(1,i,i1), 1, c2(1,j,i1), 1 )
          end do
        end do

        do j = 1, m
          do i = 1, m
            call saxpy ( m, r2(i,j), r6(1,i), 1, c1(1,j,i1), 1 )
          end do
        end do

        r6(1:m,1:m) = r5(1:m,1:m)

      end do

    end if

    c2(1:m,1:m,1) = r3(1:m,1:m)
!
!  Compute the solution of the system with the principal minor of order M*N.
!
    r3(1:m,1:m) = r1(1:m,1:m)
    x(1:m,n) = b(1:m,n)

    do i1 = 1, n-1
      i2 = n - i1
      i3 = 1
      do i = 1, m
        call saxpy ( m, -x(i,i2), a2(i3,i1), 1, x(1,n), 1 )
        i3 = i3 + m
      end do
    end do

    call sgefa ( r3, m, m, r, info )

    call sgesl ( r3, m, m, r, x(1,n), 0 )

    do i1 = 1, n-1
      do i = 1, m
        call saxpy ( m, x(i,n), c2(1,i,i1), 1, x(1,i1), 1 )
      end do
    end do

  end do

  return
end
subroutine sbto_to_sge ( m, l, a1, a2, lda, n, a )
!
!*******************************************************************************
!
!! SBTO_TO_SGE converts a block Toeplitz matrix to a Linpack General matrix.
!
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the order of the blocks of the SBTO matrix.
!
!    Input, integer L, the number of blocks in a row or column of the SBTO matrix.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the SBTO matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the SBTO matrix, beginning with the second block.
!
!    Input, integer LDA, the leading dimension of the GE matrix.
!
!    Output, integer N, the order of the GE matrix.
!
!    Output, real A(LDA,N), the N by N GE matrix.
!
  integer l
  integer lda
  integer m
!
  real a(lda,m*l)
  real a1(m,m,l)
  real a2(m,m,l-1)
  integer i
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer n
!
  n = m * l

  do i = 1, n

    i1 = ( i - 1 ) / m + 1
    i2 = i - m * ( i1 - 1 )

    do j = 1, n


      j1 = ( j - 1 ) / m + 1
      j2 = j - m * ( j1 - 1 )

      if ( j1 >= i1 ) then
        a(i,j) = a1(i2,j2,j1+1-i1)
      else
        a(i,j) = a2(i2,j2,i1-j1)
      end if

    end do

  end do

  return
end
subroutine sbto_vxm ( m, l, a1, a2, x, b )
!
!***********************************************************************
!
!! SBTO_VXM computes the real block Toeplitz matrix product X * A = B.
!
!
!  Discussion:
!
!    The full matrix has order M * L, and can be regarded
!    as comprising L by L blocks.  Each block is of order M.
!
!    Example:
!
!      M = 2, L = 3
!
!      1 2 | 3 4 | 5 6
!      5 5 | 6 6 | 7 7
!      ----+-----+-----
!      7 8 | 1 2 | 3 4
!      8 8 | 5 5 | 6 6
!      ----+-----+-----
!      9 0 | 7 8 | 1 2
!      9 9 | 8 8 | 5 5
!
!    X = (/ 1, 2, 3, 4, 5, 6 /)
!
!    B = (/ ? /)
!
!  Modified:
!
!    20 March 2001
!
!  Parameters:
!
!    Input, integer M, the order of the blocks of the matrix A.
!
!    Input, integer L, the number of blocks in a row or column of A.
!
!    Input, real A1(M,M,L), the M**2 by L matrix containing the first row of
!    blocks of the matrix.  There are L blocks, and each is of order M*M.
!
!    Input, real A2(M,M,L-1), the M**2 by L-1 matrix containing the first
!    column of blocks of the matrix, beginning with the second block.
!
!    Input, real X(M*L), the vector to be multiplied.
!
!    Output, real B(M*L), the product vector, X * A.
!
  integer l
  integer m
!
  real a1(m,m,l)
  real a2(m,m,l-1)
  real b(m,l)
  integer i
  integer j
  real x(m,l)
!
!  Construct the right hand side by blocks.
!
  do i = 1, l

    b(1:m,i) = 0.0E+00

    do j = 1, i
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a1(1:m,1:m,i+1-j) ), x(1:m,j) )
    end do

    do j = i+1, l
      b(1:m,i) = b(1:m,i) + matmul ( transpose ( a2(1:m,1:m,j-i) ), x(1:m,j) )
    end do

  end do

  return
end
subroutine scb_check ( lda, n, ml, mu, ierror )
!
!*******************************************************************************
!
!! SCB_CHECK checks the dimensions of a compact band matrix.
!
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML + MU + 1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if ML is illegal;
!    IERROR = IERROR + 4 if MU is illegal;
!    IERROR = IERROR + 8 if N is illegal.
!
  integer ierror
  integer lda
  integer ml
  integer mu
  integer n
!
  ierror = 0

  if ( lda < ml + mu + 1 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'SCB_CHECK - Illegal LDA = ', lda
  end if

  if ( ml < 0 .or. ml > n - 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SCB_CHECK - Illegal ML = ', ml
  end if

  if ( mu < 0 .or. mu > n - 1 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SCB_CHECK - Illegal MU = ', mu
  end if

  if ( n < 1 ) then
    ierror = ierror + 8
    write ( *, * ) ' '
    write ( *, * ) 'SCB_CHECK - Illegal N = ', n
  end if

  return
end
subroutine scb_det ( lda, n, ml, mu, a, det )
!
!*******************************************************************************
!
!! SCB_DET computes the determinant of a band matrix factored by SCB_NP_FA.
!
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the band matrix, as factored by SCB_NP_FA.
!
!    Output, real DET, the determinant of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  integer ierror
  integer ml
  integer mu
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  det = product ( a(mu+1,1:n) )

  return
end
subroutine scb_ml ( lda, n, ml, mu, a, x, b, job )
!
!*******************************************************************************
!
!! SCB_ML computes A * x or transpose ( A ) * X, using SCB_NP_FA factors.
!
!
!  Discussion:
!
!    It is assumed that SCB_NP_FA has overwritten the original matrix
!    information by LU factors.  SCB_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    SCB_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML + MU + 1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the compact band matrix, factored by SCB_NP_FA.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute transpose ( A ) * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer ihi
  integer ilo
  integer j
  integer jhi
  integer job
  integer ml
  integer mu
  real x(n)
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a(i-j+mu+1,j) * b(j)
      end do
      b(j) = a(j-j+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      ihi = min ( n, j + ml )
      b(j+1:ihi) = b(j+1:ihi) - a(mu+2:ihi-j+mu+1,j) * b(j)

    end do

  else
!
!  Y = transpose ( PL ) * X.
!
    do j = 1, n-1

      jhi = min ( n, j + ml )
      do i = j+1, jhi
        b(j) = b(j) - b(i) * a(i-j+mu+1,j)
      end do

    end do
!
!  B = transpose ( U ) * Y = transpose ( PL * U ) * X = transpose ( A ) * X.
!
    do i = n, 1, -1
      jhi = min ( n, i + mu )
      b(i+1:jhi) = b(i+1:jhi) + b(i) * a(mu:i-jhi+mu+1:-1,j)
      b(i) = b(i) * a(i-i+mu+1,i)
    end do

  end if

  return
end
subroutine scb_mxv ( lda, n, ml, mu, a, x, b )
!
!*******************************************************************************
!
!! SCB_MXV computes A * x, where A is a compact band matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the compact band matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer ml
  integer mu
  real x(n)
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, n
    b(i) = 0.0E+00
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine scb_np_fa ( lda, n, ml, mu, a, info )
!
!*******************************************************************************
!
!! SCB_NP_FA factors a real band matrix by Gaussian elimination.
!
!
!  Discussion:
!
!    SCB_NP_FA is a version of the LINPACK routine SGBFA, but uses no
!    pivoting.  It will fail if the matrix is singular, or if any zero
!    pivot is encountered.
!
!    Because no pivoting is used, SCB_NP_FA uses a compact band matrix
!    storage format that is more compact than the LINPACK general band format.
!
!    If SCB_NP_FA successfully factors the matrix, SCB_NP_SL may be called
!    to solve linear systems involving the matrix.
!
!  Note:
!
!    The matrix is stored in a compact version of LINPACK general
!    band storage, which does not include the fill-in entires.
!    The following program segment will store the entries of a banded
!    matrix in the compact format used by this routine:
!
!      m = mu+1
!      do j = 1, n
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
!        do i = i1, i2
!          k = i-j+m
!          a(k,j) = afull(i,j)
!        end do
!      end do
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real A(LDA,N), the compact band matrix.
!    On input, the coefficient matrix of the linear system.
!    On output, the LU factors of the matrix.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer info
  integer j
  integer ju
  integer k
  integer lm
  integer m
  integer ml
  integer mm
  integer mu
  real t
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_NP_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  The value of M is MU + 1 rather than ML + MU + 1.
!
  m = mu + 1
  info = 0
  ju = 0

  do k = 1, n-1
!
!  If our pivot entry A(MU+1,K) is zero, then we must give up.
!
    if ( a(m,k) == 0.0E+00 ) then
      info = k
      write ( *, * ) ' '
      write ( *, * ) 'SCB_FA - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if
!
!  LM counts the number of nonzero elements that lie below the current
!  diagonal entry, A(K,K).
!
!  Multiply the LM entries below the diagonal by -1/A(K,K), turning
!  them into the appropriate "multiplier" terms in the L matrix.
!
    lm = min ( ml, n-k )
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  MM points to the row in which the next entry of the K-th row is, A(K,J).
!  We then add L(I,K)*A(K,J) to A(I,J) for rows I = K+1 to K+LM.
!
    ju = max ( ju, mu + k )
    ju = min ( ju, n )
    mm = m

    do j = k+1, ju
      mm = mm - 1
      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)
    end do

  end do

  if ( a(m,n) == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'SCB_FA - Fatal error!'
    write ( *, * ) '  Zero pivot on step ', info
  end if

  return
end
subroutine scb_np_sl ( lda, n, ml, mu, a, b, job )
!
!*******************************************************************************
!
!! SCB_NP_SL solves a linear system factored by SCB_NP_FA.
!
!
!  Discussion:
!
!    SCB_NP_SL can also solve the related system transpose ( A ) * x = b.
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the compact band matrix, factored by SCB_NP_FA.
!
!    Input/output, real B(N).
!    On input, B contains the right hand side of the linear system, B.
!    On output, B contains the solution of the linear system, X.
!
!    Input, integer JOB.
!    If JOB is zero, the routine will solve A * x = b.
!    If JOB is nonzero, the routine will solve transpose ( A ) * x = b.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  integer j
  integer job
  integer k
  integer la
  integer lb
  integer lm
  integer m
  integer ml
  integer mu
  real t
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_NP_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  The value of M is ML + 1, rather than MU + ML + 1.
!
  m = mu + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    if ( ml > 0 ) then
      do k = 1, n-1
        lm = min ( ml, n-k )
        b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a(m+1:m+lm,k)
      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a(la:la+lm-1,k)

    end do
!
!  Solve transpose ( A ) * X = B.
!
  else
!
!  Solve transpose ( U ) * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(k) = ( b(k) - dot_product ( a(la:la+lm-1,k), b(lb:lb+lm-1) ) ) &
        / a(m,k)

    end do
!
!  Solve transpose ( PL ) * X = Y.
!
    if ( ml > 0 ) then

      do k = n-1, 1, -1
        lm = min ( ml, n-k )
        b(k) = b(k) + dot_product ( a(m+1:m+lm,k), b(k+1:k+lm) )
      end do

    end if

  end if

  return
end
subroutine scb_print ( lda, n, ml, mu, a, title )
!
!*******************************************************************************
!
!! SCB_PRINT prints a compact banded matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real A(LDA,N), the N by N band matrix, stored in compact
!    band storage mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer ml
  integer mu
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call scb_print_some ( lda, n, ml, mu, a, 1, 1, n, n )

  return
end
subroutine scb_print_some ( lda, n, ml, mu, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SCB_PRINT_SOME prints some of a compact banded matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real A(LDA,N), the N by N band matrix, stored in compact
!    band storage mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
!
  real a(lda,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer ml
  integer mu
  logical r_is_int
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i-j > ml .or. j-i > mu ) then
          ctemp(j2) = '              '
        else if ( r_is_int ( a(i-j+mu+1,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i-j+mu+1,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine scb_random ( lda, n, ml, mu, a )
!
!*******************************************************************************
!
!! SCB_RANDOM randomizes a compact band matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Output, real A(LDA,N), the compact band matrix, set randomly.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer ml
  integer mu
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if
!
!  Zero out the whole array.
!
  call scb_zero ( lda, n, ml, mu, a )
!
!  Set the entries that correspond to matrix elements.
!
  do i = 1, n
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      call r_random ( 0.0E+00, 1.0E+00, a(i-j+mu+1,j) )
    end do
  end do

  return
end
subroutine scb_to_sge ( lda1, lda2, ml, mu, n, a1, a2 )
!
!*******************************************************************************
!
!! SCB_TO_SGE converts a compact band matrix to general matrix format.
!
!
!  Modified:
!
!    27 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA1, the leading dimension of the array A1.
!    LDA1 must be at least ML+MU+1.
!
!    Input, integer LDA2, the leading dimension of the array A2.
!    LDA2 must be at least N.
!
!    Input, integer ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, integer N, the order of the matrices.
!    N must be positive.
!
!    Input, real A1(LDA1,N), the compact band matrix.
!
!    Output, real A2(LDA2,N), the general matrix, which contains the
!    information given in A1.
!
  integer lda1
  integer lda2
  integer n
!
  real a1(lda1,n)
  real a2(lda2,n)
  integer i
  integer ierror
  integer j
  integer ml
  integer mu
!
!  Check the dimensions.
!
  call scb_check ( lda1, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A1.'
    return
  end if

  call sge_check ( lda2, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A2.'
    return
  end if

  do i = 1, n
    do j = 1, n
      if ( j-mu <= i .and. i <= j+ml ) then
        a2(i,j) = a1(mu+1+i-j,j)
      else
        a2(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine scb_vxm ( lda, n, ml, mu, a, x, b )
!
!*******************************************************************************
!
!! SCB_VXM computes X*A, where A is a compact band matrix.
!
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the compact band matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product X*A.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer ml
  integer mu
  real x(n)
!
!  Check the dimensions.
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = 0.0E+00

  do i = 1, n
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(j) = b(j) + x(i) * a(i-j+mu+1,j)
    end do
  end do

  return
end
subroutine scb_zero ( lda, n, ml, mu, a )
!
!*******************************************************************************
!
!! SCB_ZERO zeroes out a compact band matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be nonnegative.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N-1.
!
!    Output, real A(LDA,N), the array holding the band matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer ml
  integer mu
!
  call scb_check ( lda, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  a(1:lda,1:n) = 0.0E+00

  return
end
subroutine scbb_add ( n1, n2, ml, mu, a, i, j, value )
!
!*******************************************************************************
!
!! SCBB_ADD adds a value to an entry of a compact border banded matrix.
!
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the compact border banded matrix.
!
!    Input, integer I, J, the indices of the entry to be incremented.
!
!    Input, real VALUE, the value to be added to the (I,J) entry.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer j
  real value
!
  if ( value == 0.0E+00 ) then
    return
  end if
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_ADD - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
  end if
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. i > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_ADD - Fatal error!'
    write ( *, * ) '  Illegal input value of row index I = ',i
    stop
  end if

  if ( j <= 0 .or. j > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_ADD - Fatal error!'
    write ( *, * ) '  Illegal input value of row index J = ',j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( (j-i) > mu .or. (i-j) > ml ) then
      write ( *, * ) ' '
      write ( *, * ) 'SCBB_VAL - Warning!'
      write ( *, * ) '  Unable to add to entry (', i, ',', j, ').'
      return
    else
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. j > n1 ) then
    ij = (ml+mu+1)*n1+(j-n1-1)*n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( i > n1 ) then
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1)
  end if

  a(ij) = a(ij) + value

  return
end
subroutine scbb_check ( n1, n2, ml, mu, ierror )
!
!*******************************************************************************
!
!! SCBB_CHECK checks the dimensions of a compact border banded matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1 - 1.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if ML is illegal;
!    IERROR = IERROR + 2 if MU is illegal;
!    IERROR = IERROR + 4 if N1 is illegal;
!    IERROR = IERROR + 8 if N2 is illegal;
!    IERROR = IERROR + 16 if neither N1 nor N2 is positive.
!
  integer ierror
  integer ml
  integer mu
  integer n1
  integer n2
!
  ierror = 0

  if ( ml < 0 ) then
    ierror = 1
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_CHECK - Illegal ML = ', ml
    write ( *, * ) '  but ML must be >= 0.'
  else if ( ml > max ( n1 - 1, 0 ) ) then
    ierror = 1
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_CHECK - Illegal ML = ', ml
    write ( *, * ) '  but ML must be <= Max ( N1 - 1, 0 ) = ', max ( n1 - 1, 0 )
  end if

  if ( mu < 0  ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_CHECK - Illegal MU = ', mu
    write ( *, * ) '  but MU must be >= 0.'
  else if ( mu > max ( n1 - 1, 0 ) ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_CHECK - Illegal MU = ', mu
    write ( *, * ) '  but MU must be <= Max ( N1 - 1, 0 ) = ', max ( n1 - 1, 0 )
  end if

  if ( n1 < 0 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_CHECK - Illegal N1 = ', n1
  end if

  if ( n2 < 0 ) then
    ierror = ierror + 8
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_CHECK - Illegal N2 = ', n2
  end if

  if ( n1 + n2 <= 0 ) then
    ierror = ierror + 16
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_CHECK - Illegal N1+N2 = ', n1+n2
  end if

  return
end
subroutine scbb_fa ( n1, n2, ml, mu, a, info )
!
!*******************************************************************************
!
!! SCBB_FA factors a compact border banded matrix.
!
!
!  Discussion:
!
!    Once the matrix has been factored by SCCB_FA, SCCB_SL may be called
!    to solve linear systems involving the matrix.
!
!    SCCB_FA uses special non-pivoting versions of LINPACK routines to
!    carry out the factorization.  The special version of the banded
!    LINPACK solver also results in a space saving, since no entries
!    need be set aside for fill in due to pivoting.
!
!    The linear system must be border banded, of the form:
!
!      ( A1 A2 ) (X1) = (B1)
!      ( A3 A4 ) (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    The algorithm rewrites the system as:
!
!         X1 + inv(A1) A2 X2 = inv(A1) B1
!
!      A3 X1 +         A4 X2 = B2
!
!    and then rewrites the second equation as
!
!      ( A4 - A3 inv(A1) A2 ) X2 = B2 - A3 inv(A1) B1
!
!    The algorithm will certainly fail if the matrix A1 is singular,
!    or requires pivoting.  The algorithm will also fail if the A4 matrix,
!    as modified during the process, is singular, or requires pivoting.
!    All these possibilities are in addition to the failure that will
!    if the total matrix A is singular.
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real A( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
!
!    On input, A contains the compact border-banded coefficient matrix.
!
!    In particular, A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A.
!
!    On output, A contains information describing a partial factorization
!    of the original coefficient matrix.  This information is required
!    by SCBB_SL in order to solve linear systems associated with that
!    matrix.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!      Entries of A1:
!
!        1 <= I <= N1, 1 <= J <= N1,
!        (J-I) <= MU and (I-J) <= ML.
!
!        Store the (I,J) entry into location
!        (I-J+MU+1) + (J-1) * (ML+MU+1).
!
!      Entries of A2:
!
!        1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!        Store the (I,J) entry into location
!        (ML+MU+1)*N1 + (J-N1-1)*N1 + I.
!
!      Entries of A3:
!
!        N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!        Store the (I,J) entry into location
!        (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!      Entries of A4:
!
!        N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!        Store the (I,J) entry into location
!        (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!        (same formula used for A3).
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer ik
  integer info
  integer j
  integer jk
  integer job
  integer k
  integer nband
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
  end if

  nband = (ml+mu+1)*n1
!
!  Factor the A1 band matrix, overwriting A1 by its factors.
!
  if ( n1 > 0 ) then

    call scb_np_fa ( ml+mu+1, n1, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SCBB_FA - Fatal error!'
      write ( *, * ) '  SCB_NP_FA returned INFO = ', info
      write ( *, * ) '  Factoring failed for column INFO.'
      write ( *, * ) '  The band matrix A1 is singular.'
      write ( *, * ) '  This algorithm cannot continue!'
      write ( *, * ) ' '
      return
    end if

  end if

  if ( n1 > 0 .and. n2 > 0 ) then
!
!  Set A2 := - inverse(A1) * A2.
!
    a(nband+1:nband+n1*n2) = - a(nband+1:nband+n1*n2)

    do i = 1, n2
      job = 0
      call scb_np_sl ( ml+mu+1, n1, ml, mu, a, a(nband+(i-1)*n1+1), job )
    end do
!
!  Set A4 := A4 + A3*A2
!
    do i = 1, n2
      do j = 1, n1
        ij = nband + n1*n2 + (j-1)*n2 + i
        do k = 1, n2
          ik = nband + 2*n1*n2 + (k-1)*n2 + i
          jk = nband + (k-1)*n1 + j
          a(ik) = a(ik) + a(ij) * a(jk)
        end do
      end do
    end do

  end if
!
!  Factor A4.
!
  if ( n2 > 0 ) then

    call sge_np_fa ( n2, n2, a(nband+2*n1*n2+1), info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SCBB_FA - Fatal error!'
      write ( *, * ) '  SGE_NP_FA returned INFO = ',info
      write ( *, * ) '  This indicates singularity in column INFO'
      info = n1+info
      write ( *, * ) '  of the A4 submatrix, which is column ',info
      write ( *, * ) '  of the full matrix.'
      write ( *, * ) ' '
      write ( *, * ) '  It is possible that the full matrix is '
      write ( *, * ) '  nonsingular, but the algorithm SCBB_FA may'
      write ( *, * ) '  not be used for this matrix.'
      return
    end if

  end if

  return
end
subroutine scbb_get ( n1, n2, ml, mu, a, i, j, value )
!
!*******************************************************************************
!
!! SCBB_GET returns the value of an entry of a compact border banded matrix.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the compact border banded matrix.
!
!    Input, integer I, J, the row and column of the entry to retrieve.
!
!    Output, real VALUE, the value of the (I,J) entry.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer j
  real value
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_GET - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
  end if
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. i > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_GET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index I = ',i
    stop
  end if

  if ( j <= 0 .or. j > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_GET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index J = ',j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( (j-i) > mu .or. (i-j) > ml ) then
      value = 0.0E+00
      return
    else
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. j > n1 ) then
    ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( i > n1 ) then
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
  end if

  value = a(ij)

  return
end
subroutine scbb_mxv ( n1, n2, ml, mu, a, x, y )
!
!*******************************************************************************
!
!! SCBB_MXV multiplies a compact border banded matrix times a vector.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, real A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the compact border banded matrix.
!
!    Input, real X(N1+N2), the vector to be multiplied by A.
!
!    Output, real Y(N1+N2), the result of multiplying A by X.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ihi
  integer ij
  integer ilo
  integer j
  real x(n1+n2)
  real y(n1+n2)
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
  end if
!
!  Set Y to zero.
!
  y(1:n1+n2) = 0.0E+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j-mu )
    ihi = min ( n1, j+ml )
    ij = (j-1)*(ml+mu+1)-j+mu+1
    y(ilo:ihi) = y(ilo:ihi) + a(ij+ilo:ij+ihi) * x(j)
  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (ml+mu+1)*n1+(j-n1-1)*n1
    y(1:n1) = y(1:n1) + a(ij+1:ij+n1) * x(j)
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
    y(n1+1:n1+n2) = y(n1+1:n1+n2) + a(ij+n1+1:ij+n1+n2) * x(j)
  end do

  return
end
subroutine scbb_print ( n1, n2, ml, mu, a, title )
!
!*******************************************************************************
!
!! SCBB_PRINT prints a compact border banded matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((ML+MU+1)*N1+2*N1*N2+N2*N2), the N1+N2 by N1+N2
!    compact border banded matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call scbb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2 )

  return
end
subroutine scbb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SCBB_PRINT_SOME prints some of a compact border banded matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((ML+MU+1)*N1+2*N1*N2+N2*N2), the N1+N2 by N1+N2
!    compact border banded matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ij
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n1+n2 )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n1+n2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0E+00

        if ( i <= n1 .and. j <= n1 ) then
          if ( (j-i) <= mu+ml .and. (i-j) <= ml ) then
            ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
            aij = a(ij)
          end if
        else if ( i <= n1 .and. j > n1 ) then
          ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
          aij = a(ij)
        else if ( i > n1 ) then
          ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
          aij = a(ij)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine scbb_random ( n1, n2, ml, mu, a )
!
!*******************************************************************************
!
!! SCBB_RANDOM randomizes a compact border banded matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Output, real A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the compact border banded matrix.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ihi
  integer ilo
!
  call sbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if
!
!  Randomize the band matrix A1.
!
  if ( n1 > 0 ) then
    call scb_random ( ml+mu+1, n1, ml, mu, a )
  end if
!
!  Randomize the rectangular strip A2.
!
  ilo = (ml+mu+1) * n1 + 1
  ihi = (ml+mu+1) * n1 + n1*n2

  call rvec_random ( 0.0E+00, 1.0E+00, ihi+1-ilo, a(ilo:ihi) )
!
!  Randomize the rectangular strip A3.
!
  ilo = (ml+mu+1) * n1 +   n1*n2 + 1
  ihi = (ml+mu+1) * n1 + 2*n1*n2

  call rvec_random ( 0.0E+00, 1.0E+00, ihi+1-ilo, a(ilo:ihi) )
!
!  Randomize the square matrix A4.
!
  ilo = (ml+mu+1) * n1 + 2*n1*n2 + 1
  ihi = (ml+mu+1) * n1 + 2*n1*n2 + n2*n2

  call rvec_random ( 0.0E+00, 1.0E+00, ihi+1-ilo, a(ilo:ihi) )

  return
end
subroutine scbb_set ( n1, n2, ml, mu, a, i, j, value )
!
!*******************************************************************************
!
!! SCBB_SET sets the value of an entry in a compact border banded matrix.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the compact border banded matrix.
!
!    Input, integer I, J, the row and column of the entry to set.
!
!    Input, real VALUE, the value to be assigned to the (I,J) entry.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ij
  integer j
  real value
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_SET - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
  end if
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. i > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_SET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index I = ',i
    stop
  end if

  if ( j <= 0 .or. j > n1+n2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_SET - Fatal error!'
    write ( *, * ) '  Illegal input value of row index J = ',j
    stop
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( (j-i) > mu .or. (i-j) > ml ) then
      if ( value /= 0.0E+00 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SCBB_SET - Warning!'
        write ( *, * ) '  Unable to set entry (', i, ',', j, ').'
      end if
      return
    else
      ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. j > n1 ) then
    ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( i > n1 ) then
    ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
  end if

  a(ij) = value

  return
end
subroutine scbb_sl ( n1, n2, ml, mu, a, b )
!
!*******************************************************************************
!
!! SCBB_SL solves a compact border banded system factored by SCBB_FA.
!
!
!  Discussion:
!
!    The linear system A * x = b is decomposable into the block system:
!
!      ( A1 A2 ) * (X1) = (B1)
!      ( A3 A4 )   (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    All the arguments except B are input quantities only, which are
!    not changed by the routine.  They should have exactly the same values
!    they had on exit from SCBB_FA.
!
!    If more than one right hand side is to be solved, with the same
!    matrix, SCBB_SL should be called repeatedly.  However, SCBB_FA only
!    needs to be called once to create the factorization.
!
!    See the documentation of SCBB_FA for details on the matrix storage.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
!    the compact border banded matrix, as factored by SCBB_FA.
!
!    Input/output, real B(N1+N2).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real b(n1+n2)
  integer i
  integer ierror
  integer ij
  integer j
  integer job
  integer nband
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
  end if

  nband = (ml+mu+1)*n1
!
!  Set B1 := inverse(A1) * B1.
!
  if ( n1 > 0 ) then
    job = 0
    call scb_np_sl ( ml+mu+1, n1, ml, mu, a, b, job )
  end if
!
!  Modify the right hand side of the second linear subsystem.
!  Replace B2 by B2-A3*B1.
!
  do i = 1, n2
    do j = 1, n1
      ij = nband + n1*n2 + (j-1)*n2 + i
      b(n1+i) = b(n1+i) - a(ij) * b(j)
    end do
  end do
!
!  Solve A4*B2 = B2.
!
  if ( n2 > 0 ) then
    job = 0
    call sge_np_sl ( n2, n2, a(nband+2*n1*n2+1), b(n1+1), job )
  end if
!
!  Modify the first subsolution.
!  Set B1 = B1+A2*B2.
!
  do i = 1, n1
    do j = 1, n2
      ij = nband + (j-1)*n1 + i
      b(i) = b(i) + a(ij) * b(n1+j)
    end do
  end do

  return
end
subroutine scbb_to_sge ( lda, n1, n2, ml, mu, a, a2 )
!
!*******************************************************************************
!
!! SCBB_TO_SGE copies a compact border banded matrix into a general matrix.
!
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A2.
!    LDA must be at least N1+N2.
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real A((ML+MU+1)*N1+2*N1*N2+N2*N2), the compact border banded matrix.
!
!    Output, real A2(LDA,N1+N2), a copy of the matrix, in general storage.
!
  integer lda
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real a2(lda,n1+n2)
  integer i
  integer ierror
  integer ij
  integer j
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for compact border banded matrix!'
    return
  end if

  call sge_check ( lda, n1+n2, n1+n2, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions general matrix!'
    return
  end if

  do i = 1, n1
    do j = 1, n1

      if ( (j-i) > mu+ml .or. (i-j) > ml ) then
        a2(i,j) = 0.0E+00
      else
        ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
        a2(i,j) = a(ij)
      end if

    end do
  end do

  do i = 1, n1
    do j = n1+1, n2
      ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
      a2(i,j) = a(ij)
    end do
  end do

  do i = n1+1, n2
    do j = 1, n1+n2
      ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
      a2(i,j) = a(ij)
    end do
  end do

  return
end
subroutine scbb_vxm ( n1, n2, ml, mu, a, x, y )
!
!*******************************************************************************
!
!! SCBB_VXM multiplies a vector times a compact border banded matrix.
!
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, real A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the compact border banded matrix.
!
!    Input, real X(N1+N2), the vector to multiply the matrix.
!
!    Output, real Y(N1+N2), the product X * A.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ihi
  integer ij
  integer ilo
  integer j
  real x(n1+n2)
  real y(n1+n2)
!
!  Check the dimensions.
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
  end if
!
!  Set Y to zero.
!
  y(1:n1+n2) = 0.0E+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j-mu )
    ihi = min ( n1, j+ml )
    ij = (j-1)*(ml+mu+1)-j+mu+1
    y(j) = y(j) + dot_product ( x(ilo:ihi), a(ij+ilo:ij+ihi) )
  end do
!
!  Multiply by A2.
!
  do j = n1+1, n1+n2
    ij = (ml+mu+1)*n1+(j-n1-1)*n1
    y(j) = y(j) + dot_product ( x(1:n1), a(ij+1:ij+n1) )
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1+n2
    ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
    do i = n1+1, n1+n2
      y(j) = y(j) + x(i) * a(ij+i)
    end do
  end do

  return
end
subroutine scbb_zero ( n1, n2, ml, mu, a )
!
!*******************************************************************************
!
!! SCBB_ZERO zeroes out a compact border banded matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense blocks.
!    N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Output, real A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the compact border banded matrix.
!
  integer ml
  integer mu
  integer n1
  integer n2
!
  real a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ierror
  integer ihi
!
  call scbb_check ( n1, n2, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCBB_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  ihi = (ml+mu+1)*n1 + 2*n1*n2 + n2*n2

  a(1:ihi) = 0.0E+00

  return
end
subroutine sci_eval ( n, a, lambda )
!
!*******************************************************************************
!
!! SCI_EVAL returns the eigenvalues of a real circulant matrix.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis,
!    Circulant Matrices,
!    Wiley, 1979.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the entries of the first row of the circulant matrix.
!
!    Output, complex LAMBDA(N), the eigenvalues.
!
  integer n
!
  real a(n)
  integer i
  complex lambda(n)
  complex w(n)
!
  call cvec_unity ( n, w )

  lambda(1:n) = cmplx ( a(n), 0.0E+00 )
  do i = n-1, 1, -1
    lambda(1:n) = lambda(1:n) * w(1:n) + cmplx ( a(i), 0.0E+00 )
  end do

  call cvec_sort_a2 ( n, lambda )

  return
end
subroutine sci_mxv ( n, a, x, b )
!
!*******************************************************************************
!
!! SCI_MXV multiplies a circulant matrix times a vector.
!
!
!  Modified:
!
!    07 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the entries of the first row of the circulant matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  integer j
  real x(n)
!
  do i = 1, n
    b(i) = dot_product ( a(n+2-i:n), x(1:i-1) ) &
         + dot_product ( a(1:n+1-i), x(i:n) )
  end do

  return
end
subroutine sci_print ( n, a, title )
!
!*******************************************************************************
!
!! SCI_PRINT prints a circulant matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(N), the N by N circulant matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  real a(n)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sci_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine sci_print_some ( n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SCI_PRINT_SOME prints some of a circulant matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(N), the N by N circulant matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a(n)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j >= i ) then
          aij = a(j+1-i)
        else
          aij = a(n+j+1-i)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sci_random ( n, a )
!
!*******************************************************************************
!
!! SCI_RANDOM randomizes a circulant matrix.
!
!
!  Modified:
!
!    07 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real A(N), the randomized matrix, with entries between
!    0 and 1.
!
  integer n
!
  real a(n)
!
  call rvec_random ( 0.0E+00, 1.0E+00, n, a(1:n) )

  return
end
subroutine sci_sl ( n, a, b, x, job )
!
!*******************************************************************************
!
!! SCI_SL solves the system A * x = b with the circulant matrix A.
!
!
!  Modified:
!
!    16 September 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the entries of the first row of the circulant matrix.
!
!    Input, real B(N), the right hand side.
!
!    Output, real X(N), the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  integer job
  integer nsub
  real r1
  real r2
  real r3
  real r5
  real r6
  real work(2*n-2)
  real x(n)

  if ( job == 0 ) then
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = 0.0E+00
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(n+2-nsub)
      r6 = a(nsub)

      if ( nsub > 2 ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(n+1-i) * work(nsub-i)
          r6 = r6 + a(i+1) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( nsub > 2 ) then

        r6 = work(n)
        work(n-1+nsub-1) = 0.0E+00
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = 0.0E+00
      do i = 1, nsub-1
        r5 = r5 + a(n+1-i) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      x(1:nsub-1) = x(1:nsub-1) + work(n:n+nsub-2) * r6
      x(nsub) = r6

    end do

  else
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = 0.0E+00
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(nsub)
      r6 = a(n+2-nsub)

      if ( nsub > 2 ) then

        work(nsub-1) = r2

        do i = 1, nsub-2
          r5 = r5 + a(i+1) * work(nsub-i)
          r6 = r6 + a(n+1-i) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( nsub > 2 ) then

        r6 = work(n)
        work(n-1+nsub-1) = 0.0E+00
        do i = 2, nsub-1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = 0.0E+00
      do i = 1, nsub-1
        r5 = r5 + a(i+1) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub-1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  end if

  return
end
subroutine sci_to_sge ( lda, n, a, a2 )
!
!*******************************************************************************
!
!! SCI_TO_SGE copies a circulant matrix into a general matrix.
!
!
!  Modified:
!
!    07 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the circulant matrix.
!
!    Output, real A2(LDA,N), the circulant matrix, stored as
!    a general matrix.
!
  integer lda
  integer n
!
  real a(n)
  real a2(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCI_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    a2(i,1:i-1) = a(n+2-i:n+2*1-2*i)
    a2(i,i:n) = a(1:n+1-i)
  end do

  return
end
subroutine sci_vxm ( n, a, x, b )
!
!*******************************************************************************
!
!! SCI_VXM multiplies a vector times a circulant matrix.
!
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the entries of the first row of the circulant matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product Transpose ( A ) * X.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  real x(n)
!
  do i = 1, n
    b(i) = dot_product ( a(i:1:-1), x(1:i) ) &
         + dot_product ( a(n:i+1:-1), x(i+1:n) )
  end do

  return
end
subroutine sgb_check ( lda, m, n, ml, mu, ierror )
!
!*******************************************************************************
!
!! SGB_CHECK checks the dimensions of a general band matrix.
!
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2 * ML + MU + 1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if ML is illegal;
!    IERROR = IERROR + 8 if MU is illegal;
!    IERROR = IERROR + 16 if N is illegal.
!
  integer ierror
  integer lda
  integer m
  integer ml
  integer mu
  integer n
!
  ierror = 0

  if ( lda < 2 * ml + mu + 1 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'SGB_CHECK - Illegal LDA = ', lda
  end if

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SGB_CHECK - Illegal M = ', m
  end if

  if ( ml < 0 .or. ml > min ( m, n ) - 1 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SGB_CHECK - Illegal ML = ', ml
  end if

  if ( mu < 0 .or. mu > min ( m, n ) - 1 ) then
    ierror = ierror + 8
    write ( *, * ) ' '
    write ( *, * ) 'SGB_CHECK - Illegal MU = ', mu
  end if

  if ( n < 1 ) then
    ierror = ierror + 16
    write ( *, * ) ' '
    write ( *, * ) 'SGB_CHECK - Illegal N = ', n
  end if

  return
end
subroutine sgb_det ( lda, n, ml, mu, a, pivot, det )
!
!*******************************************************************************
!
!! SGB_DET computes the determinant of a matrix factored by SGB_FA or SGB_TRF.
!
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the band matrix, as factored by SGB_FA or SGB_TRF.
!
!    Input, integer PIVOT(N), the pivot vector, as computed by SGB_FA
!    or SGB_TRF.
!
!    Output, real DET, the determinant of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  integer i
  integer ierror
  integer pivot(n)
  integer ml
  integer mu
  integer mband
!
!  Check the dimensions.
!
  call sgb_check ( lda, n, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  mband = ml + mu + 1

  det = product ( a(mband,1:n) )

  do i = 1, n
    if ( pivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine sgb_fa ( lda, n, ml, mu, a, pivot, info )
!
!*******************************************************************************
!
!! SGB_FA factors a matrix stored in LINPACK general band storage.
!
!
!  Discussion:
!
!    The matrix is stored in the array using LINPACK general band storage.
!    The following program segment will set up the input.
!
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
!        do i = i1, i2
!          k = i - j + m
!          a(k,j) = afull(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of the array A.
!    In addition, the first ML rows in the array are used for
!    elements generated during the triangularization.
!    The total number of rows needed in A is 2*ML+MU+1.
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real A(LDA,N), the matrix in band storage.  The
!    columns of the matrix are stored in the columns of the array,
!    and the diagonals of the matrix are stored in rows ML+1 through
!    2*ML+MU+1.  On return, A has been overwritten by the LU factors.
!
!    Output, integer PIVOT(N), the pivot vector.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer i0
  integer ierror
  integer info
  integer pivot(n)
  integer j
  integer j0
  integer j1
  integer ju
  integer jz
  integer k
  integer l
  integer lm
  integer m
  integer ml
  integer mm
  integer mu
  real t
!
!  Check the dimensions.
!
  call sgb_check ( lda, n, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  m = ml + mu + 1
  info = 0
!
!  Zero out the initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    a(i0:ml,jz) = 0.0E+00
  end do

  jz = j1
  ju = 0

  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      a(1:ml,jz) = 0.0E+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n-k )

    l = m
    do j = m+1, m+lm
      if ( abs ( a(j,k) ) > abs ( a(l,k) ) ) then
        l = j
      end if
    end do

    pivot(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0E+00 ) then
      info = k
      write ( *, * ) ' '
      write ( *, * ) 'SGB_FA - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange if necessary.
!
    call r_swap ( a(l,k), a(m,k) )
!
!  Compute multipliers.
!
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  Row elimination with column indexing.
!
    ju = max ( ju, mu+pivot(k) )
    ju = min ( ju, n )
    mm = m

    do j = k+1, ju

      l = l - 1
      mm = mm - 1

      if ( l /= mm ) then
        call r_swap ( a(l,j), a(mm,j) )
      end if

      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)

    end do

  end do

  pivot(n) = n
  if ( a(m,n) == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'SGB_FA - Fatal error!'
    write ( *, * ) '  Zero pivot on step ', info
  end if

  return
end
subroutine sgb_ml ( lda, n, ml, mu, a, pivot, x, b, job )
!
!*******************************************************************************
!
!! SGB_ML computes A * x or transpose ( A ) * X, using SGB_FA factors.
!
!
!  Discussion:
!
!    It is assumed that SGB_FA has overwritten the original matrix
!    information by LU factors.  SGB_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    SGB_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(LDA,N), the matrix factors computed by SGB_FA.
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML + MU + 1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, integer PIVOT(N), the pivot vector computed by SGB_FA.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute transpose ( A ) * X.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer ihi
  integer ilo
  integer pivot(n)
  integer j
  integer jhi
  integer job
  integer k
  integer ml
  integer mu
  real temp
  real x(n)
!
!  Check the dimensions.
!
  call sgb_check ( lda, n, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - ml - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a(i-j+ml+mu+1,j) * b(j)
      end do
      b(j) = a(j-j+ml+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      ihi = min ( n, j + ml )
      do i = j+1, ihi
        b(i) = b(i) - a(i-j+ml+mu+1,j) * b(j)
      end do

      k = pivot(j)

      if ( k /= j ) then
        call r_swap ( b(k), b(j) )
      end if

    end do

  else
!
!  Y = transpose ( PL ) * X.
!
    do j = 1, n-1

      k = pivot(j)

      if ( k /= j ) then
        call r_swap ( b(k), b(j) )
      end if

      jhi = min ( n, j + ml )
      do i = j+1, jhi
        b(j) = b(j) - b(i) * a(i-j+ml+mu+1,j)
      end do

    end do
!
!  B = transpose ( U ) * Y = transpose ( PL * U ) * X = transpose ( A ) * X.
!
    do i = n, 1, -1

      jhi = min ( n, i + ml + mu )
      do j = i+1, jhi
        b(j) = b(j) + b(i) * a(i-j+ml+mu+1,j)
      end do
      b(i) = b(i) * a(i-i+ml+mu+1,i)
    end do

  end if

  return
end
subroutine sgb_mu ( lda, n, ml, mu, a, pivot, x, b, job )
!
!*******************************************************************************
!
!! SGB_MU computes A * x or transpose ( A ) * X, using SGB_TRF factors.
!
!
!  Warning:
!
!    This routine must be updated to allow for rectangular matrices.
!
!  Discussion:
!
!    It is assumed that SGB_TRF has overwritten the original matrix
!    information by LU factors.  SGB_MU is able to reconstruct the
!    original matrix from the LU factor data.
!
!    SGB_MU allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML + MU + 1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the matrix factors computed by SGB_TRF.
!
!    Input, integer PIVOT(N), the pivot vector computed by SGB_TRF.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute transpose ( A ) * X.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer ihi
  integer ilo
  integer pivot(n)
  integer j
  integer jhi
  integer job
  integer k
  integer ml
  integer mu
  real temp
  real x(n)
!
!  Check the dimensions.
!
  call sgb_check ( lda, n, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_MU - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - ml - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a(i-j+ml+mu+1,j) * b(j)
      end do
      b(j) = a(j-j+ml+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      ihi = min ( n, j + ml )
      do i = j+1, ihi
        b(i) = b(i) + a(i-j+ml+mu+1,j) * b(j)
      end do

      k = pivot(j)

      if ( k /= j ) then
        call r_swap ( b(k), b(j) )
      end if

    end do

  else
!
!  Y = transpose ( PL ) * X.
!
    do j = 1, n-1

      k = pivot(j)

      if ( k /= j ) then
        call r_swap ( b(k), b(j) )
      end if

      jhi = min ( n, j + ml )
      do i = j+1, jhi
        b(j) = b(j) + b(i) * a(i-j+ml+mu+1,j)
      end do

    end do
!
!  B = transpose ( U ) * Y = transpose ( PL * U ) * X = transpose ( A ) * X.
!
    do i = n, 1, -1

      jhi = min ( n, i + ml + mu )
      do j = i+1, jhi
        b(j) = b(j) + b(i) * a(i-j+ml+mu+1,j)
      end do
      b(i) = b(i) * a(i-i+ml+mu+1,i)
    end do

  end if

  return
end
subroutine sgb_mxv ( lda, m, n, ml, mu, a, x, b )
!
!*******************************************************************************
!
!! SGB_MXV computes A * x, where A is a general band matrix.
!
!
!  Discussion:
!
!    LINPACK and LAPACK storage of general band matrices requires
!    an extra ML upper diagonals for possible fill in entries during
!    Gauss elimination.  This routine does not access any entries
!    in the fill in diagonals, because it assumes that the matrix
!    has NOT had Gauss elimination applied to it.  If the matrix
!    has been Gauss eliminated, then the routine SGB_MU must be
!    used instead.
!
!  Modified:
!
!    19 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real A(LDA,N), the M by N matrix, stored in LINPACK
!    general band matrix storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(M), the product A * x.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(m)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer ml
  integer mu
  real x(n)
!
!  Check the dimensions.
!
  call sgb_check ( lda, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, m
    b(i) = 0.0E+00
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+ml+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine sgb_print ( lda, m, n, ml, mu, a, title )
!
!*******************************************************************************
!
!! SGB_PRINT prints a banded matrix.
!
!
!  Modified:
!
!    25 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real A(LDA,N), the M by N band matrix, stored in LINPACK
!    or LAPACK general band storage mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer ml
  integer mu
  integer m
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sgb_print_some ( lda, m, n, ml, mu, a, 1, 1, m, n )

  return
end
subroutine sgb_print_some ( lda, m, n, ml, mu, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SGB_PRINT_SOME prints some of a banded matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real A(LDA,N), the M by N band matrix, stored in LINPACK
!    or LAPACK general band storage mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
!
  real a(lda,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer m
  integer ml
  integer mu
  logical r_is_int
!
!  Check the dimensions.
!
  call sgb_check ( lda, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )

    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i-j > ml .or. j-i > mu ) then
          ctemp(j2) = '              '
        else if ( r_is_int ( a(i-j+ml+mu+1,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i-j+ml+mu+1,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+ml+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sgb_random ( lda, m, n, ml, mu, a )
!
!*******************************************************************************
!
!! SGB_RANDOM randomizes a general band matrix.
!
!
!  Discussion:
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine assumes it is setting
!    up an unfactored matrix, so it only uses the first MU upper bands,
!    and does not place nonzero values in the fillin bands.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Output, real A(LDA,N), the M by N matrix.  All entries will be
!    between 0 and 1.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer m
  integer ml
  integer mu
!
!  Check the dimensions.
!
  call sgb_check ( lda, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if
!
!  Zero out the whole array.
!
  call sgb_zero ( lda, n, n, ml, mu, a )
!
!  Set the entries that correspond to matrix elements.
!
  do i = 1, m
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      call r_random ( 0.0E+00, 1.0E+00, a(i-j+ml+mu+1,j) )
    end do
  end do

  return
end
subroutine sgb_scan ( lda, m, n, ml, mu, a, nonzer, nzer )
!
!*******************************************************************************
!
!! SGB_SCAN reports the number of zeroes in a general band matrix.
!
!
!  Discussion:
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will examine
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real A(LDA,N), the M by N matrix in general band storage.
!
!    Output, integer NONZER, the number of nonzero entries in A.
!
!    Output, integer NZER, the number of zero entries in A.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer ml
  integer mu
  integer nonzer
  integer nzer
!
!  Check the dimensions.
!
  call sgb_check ( lda, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_SCAN - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  nonzer = 0
  nzer = 0

  do i = 1, m
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu + ml )
    do j = jlo, jhi
      if ( a(i-j+ml+mu+1,j) /= 0.0E+00 ) then
        nonzer = nonzer + 1
      else
        nzer = nzer + 1
      end if
    end do
  end do

  return
end
subroutine sgb_sl ( lda, n, ml, mu, a, pivot, b, job )
!
!*******************************************************************************
!
!! SGB_SL solves a system factored by SGB_FA.
!
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the LU factors from SGB_FA.
!
!    Input, integer PIVOT(N), the pivot vector from SGB_FA.
!
!    Input/output, real B(N).
!    On input, the right hand side vector.
!    On output, the solution.
!
!    Input, integer JOB.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  integer pivot(n)
  integer j
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  integer ml
  integer mu
  real t
!
!  Check the dimensions.
!
  call sgb_check ( lda, n, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  m = mu + ml + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    if ( ml >= 1 ) then

      do k = 1, n-1

        lm = min ( ml, n-k )
        l = pivot(k)

        if ( l /= k ) then
          call r_swap ( b(l), b(k) )
        end if

        b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a(m+1:m+lm,k)

      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a(la:la+lm-1,k)

    end do
!
!  Solve transpose ( A ) * X = B.
!
  else
!
!  Solve transpose(U) * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      b(k) = ( b(k) - dot_product ( a(la:la+lm-1,k), b(lb:lb+lm-1) ) ) &
        / a(m,k)
    end do
!
!  Solve transpose(L) * X = Y.
!
    if ( ml >= 1 ) then

      do k = n-1, 1, -1

        lm = min ( ml, n-k )
        b(k) = b(k) + dot_product ( a(m+1:m+lm,k), b(k+1:k+lm) )
        l = pivot(k)

        if ( l /= k ) then
          call r_swap ( b(l), b(k) )
        end if

      end do

    end if

  end if

  return
end
subroutine sgb_to_sge ( lda1, lda2, m, ml, mu, n, a1, a2 )
!
!*******************************************************************************
!
!! SGB_TO_SGE converts a general band matrix to general matrix format.
!
!
!  Discussion:
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA1, the leading dimension of the array A1.
!    LDA1 must be at least 2*ML+MU+1.
!
!    Input, integer LDA2, the leading dimension of the array A2.
!    LDA2 must be at least M.
!
!    Input, integer M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, integer N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, real A1(LDA1,N), the M by N general band matrix.
!
!    Output, real A2(LDA2,N), the M by N general matrix, which
!    contains the information given in A1.
!
  integer lda1
  integer lda2
  integer n
!
  real a1(lda1,n)
  real a2(lda2,n)
  integer i
  integer ierror
  integer j
  integer m
  integer ml
  integer mu
!
!  Check the dimensions.
!
  call sgb_check ( lda1, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A1.'
    return
  end if

  call sge_check ( lda2, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A2.'
    return
  end if

  do i = 1, m
    do j = 1, n
      if ( i - ml <= j .and. j <= i + mu + ml ) then
        a2(i,j) = a1(ml+mu+1+i-j,j)
      else
        a2(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine sgb_trf ( lda, m, n, ml, mu, a, pivot, info )
!
!*******************************************************************************
!
!! SGB_TRF performs a PLU factorization of an M by N band matrix.
!
!
!  Note:
!
!    SGB_TRF is a simplified, standalone version of the LAPACK
!    routine SGBTRF.
!
!  Modified:
!
!    18 January 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix A.  M >= 0.
!
!    Input, integer N, the number of columns of the matrix A.  N >= 0.
!
!    Input, integer ML, the number of subdiagonals within the band of A.
!    ML >= 0.
!
!    Input, integer MU, the number of superdiagonals within the band of A.
!    MU >= 0.
!
!    Input/output, real A(LDA,N).
!
!    On input, the matrix A in band storage, in rows ML+1 to
!    2*ML+MU+1; rows 1 to ML of the array need not be set.
!    The j-th column of A is stored in the j-th column of the
!    array A as follows:
!    A(ml+mu+1+i-j,j) = A(i,j) for max(1,j-mu)<=i<=min(m,j+ml)
!
!    On exit, details of the factorization: U is stored as an
!    upper triangular band matrix with ML+MU superdiagonals in
!    rows 1 to ML+MU+1, and the multipliers used during the
!    factorization are stored in rows ML+MU+2 to 2*ML+MU+1.
!
!    Output, integer PIVOT(min(M,N)), the pivot indices;
!    for 1 <= i <= min(M,N), row i of the matrix was interchanged with
!    row IPIV(i).
!
!    Output, integer INFO, error flag.
!    = 0: successful exit;
!    < 0: an input argument was illegal;
!    > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer info
  integer pivot(*)
  integer j
  integer jp
  integer ju
  integer k
  integer ml
  integer km
  integer mu
  integer kv
  real piv
  real temp
!
!  Check the dimensions.
!
  call sgb_check ( lda, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_TRF - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    info = - ierror
    return
  end if

  info = 0
!
!  KV is the number of superdiagonals in the factor U, allowing for fill-in.
!
  kv = mu + ml
!
!  Set fill-in elements in columns MU+2 to KV to zero.
!
  do j = mu + 2, min ( kv, n )
    do i = kv - j + 2, ml
      a(i,j) = 0.0E+00
    end do
  end do
!
!  JU is the index of the last column affected by the current stage
!  of the factorization.
!
  ju = 1

  do j = 1, min ( m, n )
!
!  Set the fill-in elements in column J+KV to zero.
!
    if ( j + kv <= n ) then
      a(1:ml,j+kv) = 0.0E+00
    end if
!
!  Find pivot and test for singularity.
!  KM is the number of subdiagonal elements in the current column.
!
    km = min ( ml, m-j )

    piv = abs ( a(kv+1,j) )
    jp = kv+1

    do i = kv + 2, kv + km + 1
      if ( abs ( a(i,j) ) > piv ) then
        piv = abs ( a(i,j ) )
        jp = i
      end if
    end do

    jp = jp - kv

    pivot(j) = jp + j - 1

    if( a(kv+jp,j) /= 0.0E+00 ) then

      ju = max ( ju, min ( j+mu+jp-1, n ) )
!
!  Apply interchange to columns J to JU.
!
      if ( jp /= 1 ) then

        do i = 0, ju - j
          call r_swap ( a(kv+jp-i,j+i), a(kv+1-i,j+i) )
        end do

      end if
!
!  Compute the multipliers.
!
      if ( km > 0 ) then

        a(kv+2:kv+km+1,j) = a(kv+2:kv+km+1,j) / a(kv+1,j)
!
!  Update the trailing submatrix within the band.
!
        if ( ju > j ) then

          do k = 1, ju-j

            if ( a(kv+1-k,j+k) /= 0.0E+00 ) then

              do i = 1, km
                a(kv+i+1-k,j+k) = a(kv+i+1-k,j+k) - a(kv+i+1,j) * a(kv+1-k,j+k)
              end do

            end if

          end do

        end if

      end if

    else
!
!  If pivot is zero, set INFO to the index of the pivot
!  unless a zero pivot has already been found.
!
      if ( info == 0 ) then
        info = j
      end if

    end if

  end do

  return
end
subroutine sgb_trs ( lda, n, ml, mu, nrhs, trans, a, pivot, b, ldb, info )
!
!*******************************************************************************
!
!! SGB_TRS solves a linear system factored by SGB_TRF.
!
!
!  Modified:
!
!    19 January 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer N, the order of the matrix A.
!    N must be positive.
!
!    Input, integer ML, the number of subdiagonals within the band of A.
!    ML must be at least 0, and no greater than N - 1.
!
!    Input. integer MU, the number of superdiagonals within the band of A.
!    MU must be at least 0, and no greater than N - 1.
!
!    Input, integer NRHS, the number of right hand sides and the number of
!    columns of the matrix B.  NRHS must be positive.
!
!    Input, character TRANS, specifies the form of the system.
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real A(LDA,N), contains the LU factorization of the band matrix
!    A, computed by SGB_TRF.  U is stored as an upper triangular band
!    matrix with ML+MU superdiagonals in rows 1 to ML+MU+1, and
!    the multipliers used during the factorization are stored in
!    rows ML+MU+2 to 2*ML+MU+1.
!
!    Input, integer PIVOT(N), the pivot indices; for 1 <= I <= N, row I
!    of the matrix was interchanged with row PIVOT(I).
!
!    Input/output, real B(LDB,NRHS),
!    On entry, the right hand side vectors B for the system of linear equations.
!    On exit, the solution vectors, X.
!
!    Input, integer LDB, the leading dimension of the array B.
!    LDB must be at least N.
!
!    Output, integer INFO, error flag.
!    = 0:  successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!
  integer lda
  integer ldb
  integer n
  integer nrhs
!
  real a(lda,n)
  real b(ldb,nrhs)
  integer i
  integer info
  integer pivot(*)
  integer j
  integer k
  integer kd
  integer l
  integer lm
  integer ml
  integer mu
  real temp
  character trans
!
!  Test the input parameters.
!
  info = 0

  if ( trans /= 'N' .and. trans /= 'n' .and. &
       trans /= 'T' .and. trans /= 't' .and. &
       trans /= 'C' .and. trans /= 'c' ) then
    info = -1
  else if ( n <= 0 ) then
    info = -2
  else if ( ml < 0 ) then
    info = -3
  else if ( mu < 0 ) then
    info = -4
  else if ( nrhs <= 0 ) then
    info = -5
  else if ( lda < ( 2*ml+mu+1 ) ) then
    info = -7
  else if ( ldb < max ( 1, n ) ) then
    info = -10
  end if

  if ( info /= 0 ) then
    return
  end if

  kd = mu + ml + 1
!
!  Solve A * x = b.
!
!  Solve L * x = b, overwriting b with x.
!
!  L is represented as a product of permutations and unit lower
!  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!  where each transformation L(i) is a rank-one modification of
!  the identity matrix.
!
  if ( trans == 'N' .or. trans == 'n' ) then

    if ( ml > 0 ) then

      do j = 1, n - 1

        lm = min ( ml, n-j )
        l = pivot(j)

        do i = 1, nrhs
          call r_swap ( b(l,i), b(j,i) )
        end do

        do k = 1, nrhs
          if ( b(j,k) /= 0.0E+00 ) then
            b(j+1:j+lm,k) = b(j+1:j+lm,k) - a(kd+1:kd+lm,j) * b(j,k)
          end if
        end do

      end do

    end if
!
!  Solve U * x = b, overwriting b with x.
!
    do i = 1, nrhs

      do j = n, 1, -1
        if ( b(j,i) /= 0.0E+00 ) then
          l = ml + mu + 1 - j
          b(j,i) = b(j,i) / a(ml+mu+1,j)
          do k = j - 1, max ( 1, j - ml - mu ), -1
            b(k,i) = b(k,i) - a(l+k,j) * b(j,i)
          end do
        end if
      end do

    end do

  else
!
!  Solve Transpose ( A ) * x = b.
!
!  Solve Transpose ( U ) * x = b, overwriting b with x.
!
    do i = 1, nrhs

      do j = 1, n
        temp = b(j,i)
        l = ml + mu + 1 - j
        do k = max ( 1, j - ml - mu ), j - 1
          temp = temp - a(l+k,j) * b(k,i)
        end do
        temp = temp / a(ml+mu+1,j)
        b(j,i) = temp
      end do

    end do
!
!  Solve Transpose ( L ) * x = b, overwriting b with x.
!
    if ( ml > 0 ) then

      do j = n - 1, 1, -1

        lm = min ( ml, n-j )

        do k = 1, nrhs
          b(j,k) = b(j,k) - dot_product ( b(j+1:j+lm,k), a(kd+1:kd+lm,j) )
        end do

        l = pivot(j)

        do i = 1, nrhs
          call r_swap ( b(l,i), b(j,i) )
        end do

      end do

    end if

  end if

  return
end
subroutine sgb_vxm ( lda, m, n, ml, mu, a, x, b )
!
!*******************************************************************************
!
!! SGB_VXM computes X*A, where A is a general band matrix.
!
!
!  Discussion:
!
!    LINPACK and LAPACK storage of general band matrices requires
!    an extra ML upper diagonals for possible fill in entries during
!    Gauss elimination.  This routine does not access any entries
!    in the fill in diagonals, because it assumes that the matrix
!    has NOT had Gauss elimination applied to it.  If the matrix
!    has been Gauss eliminated, then the routine SGB_MU must be
!    used instead.
!
!  Modified:
!
!    19 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real A(LDA,N), the M by N matrix in LINPACK general
!    band storage.
!
!    Input, real X(M), the vector to be multiplied by A.
!
!    Output, real B(N), the product X*A.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer ml
  integer mu
  real x(m)
!
!  Check the dimensions.
!
  call sgb_check ( lda, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = 0.0E+00

  do i = 1, m
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(j) = b(j) + x(i) * a(i-j+ml+mu+1,j)
    end do
  end do

  return
end
subroutine sgb_zero ( lda, m, n, ml, mu )
!
!*******************************************************************************
!
!! SGB_ZERO zeroes out a general band matrix.
!
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be nonnegative.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be nonnegative.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than min(M,N)-1.
!
!    Output, real A(LDA,N), the array holding the M by N band matrix.
!    The entire LDA by N array is zeroed out, not just the portion
!    representing legal matrix entries.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer m
  integer ml
  integer mu
!
  call sgb_check ( lda, m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  a(1:lda,1:n) = 0.0E+00

  return
end
subroutine sgd_check ( lda, n, ndiag, ierror )
!
!*******************************************************************************
!
!! SGD_CHECK checks the dimensions of a general diagonal matrix.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if N is illegal;
!    IERROR = IERROR + 4 if NDIAG is illegal.
!
  integer ierror
  integer lda
  integer n
  integer ndiag
!
  ierror = 0

  if ( lda < n ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'SGD_CHECK - Illegal LDA = ', lda
  end if

  if ( n < 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SGD_CHECK - Illegal N = ', n
  end if

  if ( ndiag < 1 .or. ndiag > 2 * n - 1 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SGD_CHECK - Illegal NDIAG = ', ndiag
  end if

  return
end
subroutine sgd_mxv ( lda, n, ndiag, offset, a, x, b )
!
!*******************************************************************************
!
!! SGD_MXV computes A * x where A is a general diagonal matrix.
!
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the matrix in general diagonal storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real b(n)
  integer i
  integer ierror
  integer j
  integer jdiag
  integer offset(ndiag)
  real x(n)
!
!  Check the dimensions.
!
  call sgd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGD_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  b(1:n) = 0.0E+00

  do i = 1, n
    do jdiag = 1, ndiag
      j = i + offset(jdiag)
      if ( j >= 1 .and. j <= n ) then
        b(i) = b(i) + a(i,jdiag) * x(j)
      end if
    end do
  end do

  return
end
subroutine sgd_print ( lda, n, ndiag, offset, a, title )
!
!*******************************************************************************
!
!! SGD_PRINT prints a general diagonal matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the N by N general diagonal matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  integer offset(ndiag)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sgd_print_some ( lda, n, ndiag, offset, a, 1, 1, n, n )

  return
end
subroutine sgd_print_some ( lda, n, ndiag, offset, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SGD_PRINT_SOME prints some of a general diagonal matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the N by N general diagonal matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jdiag
  integer jhi
  integer jlo
  integer off
  integer offset(ndiag)
  logical r_is_int
!
!  Check the dimensions.
!
  call sgd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGD_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0E+00
        off = j - i
        do jdiag = 1, ndiag
          if ( off == offset(jdiag) ) then
            aij = a(i,jdiag)
          end if
        end do

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sgd_random ( lda, n, ndiag, offset, a )
!
!*******************************************************************************
!
!! SGD_RANDOM randomizes a general diagonal matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Output, real A(LDA,NDIAG), the matrix in general diagonal storage.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  integer i
  integer ierror
  integer j
  integer jj
  integer offset(ndiag)
!
!  Check the dimensions.
!
  call sgd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGD_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        call r_random ( 0.0E+00, 1.0E+00, a(i,j) )
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine sgd_to_sge ( lda, lda2, n, ndiag, offset, a, a2 )
!
!*******************************************************************************
!
!! SGD_TO_SGE copies a general diagonal matrix into a general matrix.
!
!
!  Modified:
!
!    30 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer LDA2, the leading dimension of the array A2.
!    LDA2 must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the matrix in general diagonal storage.
!
!    Input, real A2(LDA2,N), a copy of the matrix, in general storage.
!
  integer lda
  integer lda2
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real a2(lda2,n)
  integer i
  integer ierror
  integer j
  integer jj
  integer offset(ndiag)
!
!  Check the dimensions.
!
  call sgd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGD_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A.'
    return
  end if

  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGD_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A2.'
    return
  end if

  a2(1:n,1:n) = 0.0E+00

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        a2(i,jj) = a(i,j)
      end if
    end do
  end do

  return
end
subroutine sgd_vxm ( lda, n, ndiag, offset, a, x, b )
!
!*******************************************************************************
!
!! SGD_VXM computes X*A where A is a general diagonal matrix.
!
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the matrix, in general diagonal storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product X*A.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real b(n)
  integer i
  integer ierror
  integer j
  integer jdiag
  integer offset(ndiag)
  real x(n)
!
!  Check the dimensions.
!
  call sgd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGD_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  b(1:n) = 0.0E+00

  do i = 1, n
    do jdiag = 1, ndiag
      j = i + offset(jdiag)
      if ( 1 <= j .and. j <= n ) then
        b(j) = b(j) + x(i) * a(i,jdiag)
      end if
    end do
  end do

  return
end
subroutine sgd_zero ( lda, n, ndiag, a )
!
!*******************************************************************************
!
!! SGD_ZERO zeroes out a general diagonal matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than 2 * N - 1.
!
!    Output, real A(LDA,NDIAG), the matrix in general diagonal storage.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sgd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGD_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  a(1:lda,1:ndiag) = 0.0E+00

  return
end
subroutine sge_check ( lda, m, n, ierror )
!
!*******************************************************************************
!
!! SGE_CHECK checks the dimensions of a general matrix.
!
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  integer ierror
  integer lda
  integer m
  integer n
!
  ierror = 0

  if ( lda < m ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'SGE_CHECK - Illegal LDA = ', lda
  end if

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SGE_CHECK - Illegal M = ', m
  end if

  if ( n < 1 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SGE_CHECK - Illegal N = ', n
  end if

  return
end
subroutine sge_det ( lda, n, a, pivot, det )
!
!*******************************************************************************
!
!! SGE_DET computes the determinant of a matrix factored by SGE_FA or SGE_TRF.
!
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the LU factors computed by SGE_FA or SGE_TRF.
!
!    Input, integer PIVOT(N), as computed by SGE_FA or SGE_TRF.
!
!    Output, real DET, the determinant of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  real diag(n)
  integer i
  integer ierror
  integer pivot(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  call rmat_diag_get_vector ( lda, n, a, diag )

  det = product ( diag(1:n) )

  do i = 1, n
    if ( pivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine sge_dilu ( lda, m, n, a, d )
!
!*******************************************************************************
!
!! SGE_DILU produces the diagonal incomplete LU factors of a real rectangular matrix.
!
!
!  Discussion:
!
!    The D-ILU factors of the M by N matrix A are:
!
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix.
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real A(LDA,N), the M by N matrix to be factored.
!
!    Output, real D(M), the D-ILU factor.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real d(m)
  integer i
  integer j
!
  d(1:n) = 0.0E+00

  do i = 1, min ( m, n )
    d(i) = a(i,i)
  end do

  do i = 1, m
    d(i) = 1.0E+00 / d(i)
    do j = i+1, m
      if ( a(i,j) /= 0.0E+00 .and. a(j,i) /= 0.0E+00 ) then
        d(j) = d(j) - a(j,i) * d(i) * a(i,j)
      end if
    end do
  end do

  return
end
subroutine sge_fa ( lda, n, a, pivot, info )
!
!*******************************************************************************
!
!! SGE_FA factors a general matrix.
!
!
!  Discussion:
!
!    SGE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real A(LDA,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer PIVOT(N), a vector of pivot indices.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer info
  integer pivot(n)
  integer j
  integer k
  integer l
  real t
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(i,k) ) > abs ( a(l,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0E+00 ) then
      info = k
      write ( *, * ) ' '
      write ( *, * ) 'SGE_FA - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r_swap ( a(l,k), a(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n

      if ( l /= k ) then
        call r_swap ( a(l,j), a(k,j) )
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'SGE_FA - Fatal error!'
    write ( *, * ) '  Zero pivot on step ', info
  end if

  return
end
subroutine sge_fs ( lda, n, a, b, info )
!
!*******************************************************************************
!
!! SGE_FS factors and solves a general linear system in one step.
!
!
!  Note:
!
!    SGE_FS does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    SGE_FS uses partial pivoting, but no pivot vector is required.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real A(LDA,N).
!
!    On input, A is the coefficient matrix of the linear system.
!
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real B(N).
!    On input, B is the right hand side of the linear system.
!    On output, B is the solution of the linear system.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer info
  integer ipiv
  integer j
  integer jcol
  integer jj
  real piv
  real temp
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_FS - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol+1, n
      if ( abs ( a(i,jcol) ) > piv ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0E+00 ) then
      info = jcol
      write ( *, * ) ' '
      write ( *, * ) 'SGE_FS - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      do j = 1, n
        call r_swap ( a(jcol,j), a(ipiv,j) )
      end do

      call r_swap ( b(jcol), b(ipiv) )

    end if
!
!  Scale the pivot row.
!
    temp = a(jcol,jcol)
    a(jcol,jcol) = 1.0E+00
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / temp
    b(jcol) = b(jcol) / temp
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol+1, n
      if ( a(i,jcol) /= 0.0E+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0E+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i) = b(i) + temp * b(jcol)
      end if
    end do

  end do
!
!  Back solve.
!
  do jcol = n, 2, -1
    b(1:jcol-1) = b(1:jcol-1) - a(1:jcol-1,jcol) * b(jcol)
  end do

  return
end
subroutine sge_identity ( lda, n, a )
!
!*******************************************************************************
!
!! SGE_IDENTITY sets up the identity matrix in real general storage.
!
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of A.
!
!    Output, real A(LDA,N), the N by N identity matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real, parameter :: one = 1.0E+00
!
  a(1:n,1:n) = 0.0E+00

  call rmat_diag_set_scalar ( lda, n, a, one )

  return
end
subroutine sge_ilu ( lda, m, n, a, l, u )
!
!*******************************************************************************
!
!! SGE_ILU produces the incomplete LU factors of a real rectangular matrix.
!
!
!  Discussion:
!
!    The incomplete LU factors of the M by N matrix A are:
!
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix
!
!    with the property that L and U are computed in the same way as
!    the usual LU factors, except that, whenever an off diagonal element
!    of the original matrix is zero, then the corresponding value of
!    U is forced to be zero.
!
!    This condition means that it is no longer the case that A = L*U.
!
!    On the other hand, L and U will have a simple sparsity structure
!    related to that of A.  The incomplete LU factorization is generally
!    used as a preconditioner in iterative schemes applied to sparse
!    matrices.  It is presented here merely for illustration.
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real A(LDA,N), the M by N matrix to be factored.
!
!    Output, real L(LDA,M), the M by M unit lower triangular factor.
!
!    Output, real U(LDA,N), the M by N upper triangular factor.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer k
  real l(lda,m)
  real u(lda,n)
!
!  Initialize:
!
!    L := M by M Identity
!    U := A
!
  call sge_identity ( lda, m, l )

  u(1:m,1:n) = a(1:m,1:n)

  do j = 1, min ( m-1, n )
!
!  Zero out the entries in column J, from row J+1 to M.
!
    do i = j+1, m

      if ( u(i,j) /= 0.0E+00 ) then

        l(i,j) = u(i,j) / u(j,j)
        u(i,j) = 0.0E+00

        do k = j+1, n
          if ( u(i,k) /= 0.0E+00 ) then
            u(i,k) = u(i,k) - l(i,j) * u(j,k)
          end if
        end do

      end if

    end do

  end do

  return
end
subroutine sge_inv ( lda, n, a, pivot )
!
!*******************************************************************************
!
!! SGE_INV computes the inverse of a matrix factored by SGE_FA.
!
!
!  Note:
!
!    SGE_INV is a simplified standalone version of the LINPACK routine
!    SGEDI.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A,
!    which must be at least N.
!
!    Input, integer N, the order of the matrix A.
!
!    Input/output, real A(LDA,N).
!    On input, the factor information computed by SGE_FA.
!    On output, the inverse matrix.
!
!    Input, integer PIVOT(N), the pivot vector from SGE_FA.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer pivot(n)
  integer j
  integer k
  real temp
  real work(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_INV - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Compute Inverse(U).
!
  do k = 1, n

    a(k,k) = 1.0E+00 / a(k,k)
    a(1:k-1,k) = -a(1:k-1,k) * a(k,k)

    do j = k + 1, n

      temp = a(k,j)
      a(k,j) = 0.0E+00
      a(1:k,j) = a(1:k,j) + temp * a(1:k,k)

    end do

  end do
!
!  Form Inverse(U) * Inverse(L).
!
  do k = n - 1, 1, -1

    work(k+1:n) = a(k+1:n,k)
    a(k+1:n,k) = 0.0E+00

    do j = k + 1, n
      a(1:n,k) = a(1:n,k) + work(j) * a(1:n,j)
    end do

    if ( pivot(k) /= k ) then

      do i = 1, n
        call r_swap ( a(i,k), a(i,pivot(k)) )
      end do

    end if

  end do

  return
end
subroutine sge_ml ( lda, n, a, pivot, x, b, job )
!
!*******************************************************************************
!
!! SGE_ML computes A * x or transpose ( A ) * x, using SGE_FA factors.
!
!
!  Discussion:
!
!    It is assumed that SGE_FA has overwritten the original matrix
!    information by LU factors.  SGE_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    SGE_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix factors computed by SGE_FA.
!
!    Input, integer PIVOT(N), the pivot vector computed by SGE_FA.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute transpose ( A ) * X.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer pivot(n)
  integer j
  integer job
  integer k
  real temp
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      b(1:j-1) = b(1:j-1) + a(1:j-1,j) * b(j)
      b(j) = a(j,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      b(j+1:n) = b(j+1:n) - a(j+1:n,j) * b(j)

      k = pivot(j)

      if ( k /= j ) then
        call r_swap ( b(k), b(j) )
      end if

    end do

  else
!
!  Y = tranpose(PL) * X:
!
    do j = 1, n-1

      k = pivot(j)

      if ( k /= j ) then
        call r_swap ( b(k), b(j) )
      end if

      b(j) = b(j) - dot_product ( b(j+1:n), a(j+1:n,j) )

    end do
!
!  B = transpose(U) * Y = transpose ( PL * U ) * X = transpose ( A ) * X.
!
    do i = n, 1, -1
      b(i+1:n) = b(i+1:n) + b(i) * a(i,i+1:n)
      b(i) = b(i) * a(i,i)
    end do

  end if

  return
end
subroutine sge_mu ( lda, m, n, a, trans, pivot, x, b )
!
!*******************************************************************************
!
!! SGE_MU computes A * x or transpose ( A ) * x, using SGE_TRF factors.
!
!
!  Discussion:
!
!    It is assumed that SGE_TRF has overwritten the original matrix
!    information by PLU factors.  SGE_MU is able to reconstruct the
!    original matrix from the PLU factor data.
!
!    SGE_MU allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Modified:
!
!    14 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows in the matrix.
!
!    Input, integer N, the number of columns in the matrix.
!
!    Input, real A(LDA,N), the matrix factors computed by SGE_TRF.
!
!    Input, character TRANS, specifies the form of the system of equations:
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, integer PIVOT(*), the pivot vector computed by SGE_TRF.
!
!    Input, real X(*), the vector to be multiplied.
!    For the untransposed case, X should have N entries.
!    For the transposed case, X should have M entries.
!
!    Output, real B(*), the result of the multiplication.
!    For the untransposed case, B should have M entries.
!    For the transposed case, B should have N entries.
!
  integer, parameter :: MN_MAX = 100
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(*)
  integer i
  integer ierror
  integer pivot(*)
  integer j
  integer k
  integer npiv
  real temp
  character trans
  real x(*)
  real y(MN_MAX)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_MU - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  npiv = min ( m - 1, n )

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Y[MN] = U[MNxN] * X[N].
!
    y(1:n) = 0.0E+00

    do j = 1, n

      do i = 1, min ( j, m )
        y(i) = y(i) + a(i,j) * x(j)
      end do

    end do
!
!  Z[M] = L[MxMN] * Y[MN] = L[MxMN] * U[MNxN] * X[N].
!
    do i = 1, m

      if ( i <= n ) then
        b(i) = y(i)
      else
        b(i) = 0.0E+00
      end if

    end do

    do j = min ( m-1, n ), 1, -1
      b(j+1:m) = b(j+1:m) + a(j+1:m,j) * y(j)
    end do
!
!  B = P * Z = P * L * Y = P * L * U * X = A * x.
!
    do j = npiv, 1, -1

      k = pivot(j)

      if ( k /= j ) then
        call r_swap ( b(k), b(j) )
      end if

    end do

  else if ( trans == 't' .or. trans == 'T' .or. &
            trans == 'c' .or. trans == 'C' ) then
!
!  Y = tranpose(P) * X:
!
    do i = 1, npiv

      k = pivot(i)

      if ( k /= i ) then
        call r_swap ( x(k), x(i) )
      end if

    end do

    b(1:m) = x(1:m)
    b(m+1:n) = 0.0E+00
!
!  Z = tranpose(L) * Y:
!
    do j = 1, min ( m - 1, n )
      b(j) = b(j) + dot_product ( x(j+1:m), a(j+1:m,j) )
    end do
!
!  B = transpose(U) * Z.
!
    do i = m, 1, -1
      b(i+1:n) = b(i+1:n) + b(i) * a(i,i+1:n)
      if ( i <= n ) then
        b(i) = b(i) * a(i,i)
      end if
    end do
!
!  Now restore X.
!
     do i = npiv, 1, -1

      k = pivot(i)

      if ( k /= i ) then
        call r_swap ( x(k), x(i) )
      end if

    end do

  else

    write ( *, * ) ' '
    write ( *, * ) 'SGE_MU - Fatal error!'
    write ( *, '(a,a)' ) '  Illegal value of TRANS = ', trans
    stop

  end if

  return
end
subroutine sge_mxm ( lda, n, a, b, c )
!
!*******************************************************************************
!
!! SGE_MXM computes A * B = C, where A, B and C are N by N matrices.
!
!
!  Modified:
!
!    12 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the arrays.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrices.
!    N must be positive.
!
!    Input, real A(LDA,N), B(LDA,N), the N by N factor matrices, stored
!    in LINPACK general matrix storage.
!
!    Output, real C(LDA,N), the N by N product matrix, stored in
!    LINPACK general matrix storage.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  integer ierror
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_MXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  return
end
subroutine sge_mxv ( lda, m, n, a, x, b )
!
!*******************************************************************************
!
!! SGE_MXV computes A * x, where A is a general matrix.
!
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N matrix, stored in LINPACK
!    general matrix storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(M), the product A * x.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(m)
  integer ierror
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:m) = matmul ( a(1:m,1:n), x(1:n) )

  return
end
subroutine sge_np_det ( lda, n, a, det )
!
!*******************************************************************************
!
!! SGE_NP_DET computes the determinant of a matrix factored by SGE_NP_FA.
!
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the LU factors computed by SGE_FA.
!
!    Output, real DET, the determinant of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  real diag(n)
  integer ierror
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_NP_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  call rmat_diag_get_vector ( lda, n, a, diag )

  det = product ( diag(1:n) )

  return
end
subroutine sge_np_fa ( lda, n, a, info )
!
!*******************************************************************************
!
!! SGE_NP_FA factors a general matrix by nonpivoting Gaussian elimination.
!
!
!  Discussion:
!
!    SGE_NP_FA is a version of the LINPACK routine SGEFA, but uses no
!    pivoting.  It will fail if the matrix is singular, or if any zero
!    pivot is encountered.
!
!    If SGE_NP_FA successfully factors the matrix, SGE_NP_SL may be called
!    to solve linear systems involving the matrix.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real A(LDA,N).
!
!    On input, A contains the matrix to be factored.
!    On output, A contains information about the factorization,
!    which must be passed unchanged to SGE_NP_SL for solutions.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer info
  integer j
  integer k
  real t
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_NP_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  info = 0

  do k = 1, n-1

    if ( a(k,k) == 0.0E+00 ) then
      info = k
      write ( *, * ) ' '
      write ( *, * ) 'SGE_NP_FA - Fatal error!'
      write ( *, * ) '  Zero pivot on step ', info
      return
    end if

    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
    do j = k+1, n
      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)
    end do

  end do

  if ( a(n,n) == 0.0E+00 ) then
    info = n
    write ( *, * ) ' '
    write ( *, * ) 'SGE_NP_FA - Fatal error!'
    write ( *, * ) '  Zero pivot on step ', info
  end if

  return
end
subroutine sge_np_inv ( lda, n, a )
!
!*******************************************************************************
!
!! SGE_NP_INV computes the inverse of a matrix factored by SGE_NP_FA.
!
!
!  Modified:
!
!    12 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A,
!    which must be at least N.
!
!    Input, integer N, the order of the matrix A.
!
!    Input/output, real A(LDA,N).
!    On input, the factor information computed by SGE_NP_FA.
!    On output, the inverse matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer k
  real temp
  real work(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_NP_INV - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Compute Inverse(U).
!
  do k = 1, n

    a(k,k) = 1.0E+00 / a(k,k)
    a(1:k-1,k) = - a(1:k-1,k) * a(k,k)

    do j = k + 1, n

      temp = a(k,j)
      a(k,j) = 0.0E+00
      a(1:k,j) = a(1:k,j) + temp * a(1:k,k)

    end do

  end do
!
!  Form Inverse(U) * Inverse(L).
!
  do k = n - 1, 1, -1

    work(k+1:n) = a(k+1:n,k)
    a(k+1:n,k) = 0.0E+00

    do j = k + 1, n
      a(1:n,k) = a(1:n,k) + a(1:n,j) * work(j)
    end do

  end do

  return
end
subroutine sge_np_ml ( lda, n, a, x, b, job )
!
!*******************************************************************************
!
!! SGE_NP_ML computes A * x or x * A, for a matrix factored by SGE_NP_FA.
!
!
!  Discussion:
!
!    The matrix A is assumed to have been factored by SGE_NP_FA.
!
!    SGE_NP_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix factors computed by SGE_NP_FA.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real B(N), the result of the multiplication.
!
!    Input, integer JOB, determines the multiplication to
!    be carried out:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute transpose ( A ) * X.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer j
  integer job
  real temp
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_NP_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute U * X = Y:
!
    do i = 1, n
      b(i) = dot_product ( a(i,i:n), b(i:n) )
    end do
!
!  Compute L * Y = B:
!
    do j = n-1, 1, -1
      b(j+1:n) = b(j+1:n) - a(j+1:n,j) * b(j)
    end do

  else
!
!  Compute tranpose(L) * X = Y:
!
    do i = 1, n-1
      b(i) = b(i) - dot_product ( b(i+1:n), a(i+1:n,i) )
    end do
!
!  Compute transpose(U) * Y = B:
!
    do i = n, 1, -1
      b(i) = dot_product ( b(1:i), a(1:i,i) )
    end do

  end if

  return
end
subroutine sge_np_sl ( lda, n, a, b, job )
!
!*******************************************************************************
!
!! SGE_NP_SL solves a system factored by SGE_NP_FA.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix as factored by SGE_NP_FA.
!
!    Input/output, real B(N).
!
!    On input, B contains the right hand side vector B.
!    On output, B contains the solution X.
!
!    Input, integer JOB.
!    If JOB is zero, the routine will solve A * x = b.
!    If JOB is nonzero, the routine will solve transpose ( A ) * x = b.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  integer job
  integer j
  integer k
  real t
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_NP_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if
!
!  Solve A * x = b.
!
  if ( job == 0 ) then

    do k = 1, n-1
      b(k+1:n) = b(k+1:n) + a(k+1:n,k) * b(k)
    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      b(1:k-1) = b(1:k-1) - a(1:k-1,k) * b(k)
    end do
!
!  Solve transpose ( A ) * X = B.
!
  else

    do k = 1, n
      b(k) = ( b(k) - dot_product ( b(1:k-1), a(1:k-1,k) ) ) / a(k,k)
    end do

    do k = n-1, 1, -1
      b(k) = b(k) + dot_product ( b(k+1:n), a(k+1:n,k) )
    end do

  end if

  return
end
subroutine sge_np_trf ( lda, m, n, a, info )
!
!*******************************************************************************
!
!! SGE_NP_TRF computes the LU factorization of a general M by N matrix.
!
!
!  Note:
!
!    SGE_NP_TRF is a nonpivoting version of SGE_TRF, and will fail if
!    a zero element is encountered along the diagonal.
!
!    The factorization has the form
!      A = L * U
!    where L is lower triangular with unit diagonal elements (lower
!    trapezoidal if M > N), and U is upper triangular (upper trapezoidal
!    if M < N).
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= max(1,M).
!
!    Input, integer M, the number of rows of the matrix A.  M >= 0.
!
!    Input, integer N, the number of columns of the matrix A.  N >= 0.
!
!    Input/output, real A(LDA,N).
!    On entry, the M by N matrix to be factored.
!    On exit, the factors L and U from the factorization
!    A = L*U; the unit diagonal elements of L are not stored.
!
!    Output, integer INFO.
!    = 0: successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!    > 0: if INFO = K, U(K,K) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer ii
  integer info
  integer j
  integer jj
  integer m
!
!  Test the input parameters.
!
  info = 0

  if ( m < 0 ) then
    info = - 1
    return
  else if( n < 0 ) then
    info = - 2
    return
  else if ( lda < max ( 1, m ) ) then
    info = - 4
    return
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  do j = 1, min ( m, n )
!
!  Compute elements J+1:M of the J-th column.
!
    if ( a(j,j) /= 0.0E+00 ) then
      a(j+1:m,j) = a(j+1:m,j) / a(j,j)
    else if ( info == 0 ) then
      info = j
    end if
!
!  Update the trailing submatrix.
!
    if ( j < min ( m, n ) ) then

      do ii = j+1, m
        a(ii,j+1:n) = a(ii,j+1:n) - a(ii,j) * a(j,j+1:n)
      end do

    end if

  end do

  return
end
subroutine sge_np_trm ( lda, m, n, a, x, b, job )
!
!*******************************************************************************
!
!! SGE_NP_TRM computes A * x or x * A, for a matrix factored by SGE_NP_TRF.
!
!
!  Discussion:
!
!    The matrix A is assumed to have been factored by SGE_NP_TRF.
!
!    SGE_NP_TRM allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Modified:
!
!    24 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer M, N, the number of rows and columns in the matrix.
!    M and N must be positive.
!
!    Input, real A(LDA,N), the M by N matrix factors computed by SGE_NP_TRF.
!
!    Input, real X(*), the vector to be multiplied.
!    If JOB is 0, X must have dimension N.
!    If JOB is nonzero, X must have dimension M.
!
!    Output, real B(*), the result of the multiplication.
!    If JOB is 0, B must have dimension M.
!    If JOB is nonzero, B must have dimension N.
!
!    Input, integer JOB, determines the multiplication to
!    be carried out:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute transpose ( A ) * X.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(*)
  integer i
  integer ierror
  integer job
  integer m
  real x(*)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_NP_TRM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  if ( job == 0 ) then
!
!  Compute U * X = Y:
!  Compute L * Y = B:
!
    do i = 1, min ( m, n )
      b(i) = dot_product ( a(i,i:n), x(i:n) )
    end do

    if ( n < m ) then
      b(n+1:m) = 0.0E+00
    end if

    do i = m, 2, -1
      b(i) = b(i) + dot_product ( a(i,1:i-1), b(1:i-1) )
    end do
!
!  Compute tranpose(L) * X = Y:
!  Compute transpose(U) * Y = B:
!
  else

    do i = 1, min ( m, n )
      b(i) = x(i) + dot_product ( a(i+1:m,i), x(i+1:m) )
    end do

    if ( m < n ) then
      b(m+1:n) = 0.0E+00
    end if

    do i = m, 1, -1
      b(i) = dot_product ( a(1:i,i), b(1:i) )
    end do

  end if

  return
end
subroutine sge_np_trs ( lda, n, nrhs, trans, a, b, ldb, info )
!
!*******************************************************************************
!
!! SGE_NP_TRS solves a system of linear equations factored by SGE_NP_TRF.
!
!
!  Note:
!
!    SGE_NP_TRS is a nonpivoting version of SGE_TRS.
!
!    SGE_TRS solves a system of linear equations
!      A * x = b  or  A' * X = B
!    with a general N by N matrix A using the LU factorization computed
!    by SGE_NP_TRF.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= max(1,N).
!
!    Input, integer N, the order of the matrix A.  N >= 0.
!
!    Input, integer NRHS, the number of right hand sides.  NRHS >= 0.
!
!    Input, character TRANS, pecifies the form of the system of equations:
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real A(LDA,N), the factors L and U from the factorization
!    A = L*U as computed by SGE_NP_TRF.
!
!    Input/output, real B(LDB,NRHS).
!    On entry, the right hand side matrix B.
!    On exit, the solution matrix X.
!
!    Input, integer LDB, the leading dimension of the array B.
!    LDB >= max(1,N).
!
!    Output, integer INFO
!    = 0:  successful exit
!    < 0:  if INFO = -I, the I-th argument had an illegal value.
!
  integer lda
  integer ldb
  integer n
  integer nrhs
!
  real a(lda,n)
  real b(ldb,nrhs)
  integer i
  integer info
  integer j
  integer k
  character trans
!
  info = 0

  if ( trans /= 'n' .and. trans /= 'N' .and. &
       trans /= 't' .and. trans /= 'T' .and. &
       trans /= 'c' .and. trans /= 'C' ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  else if ( nrhs < 0 ) then
    info = - 3
    return
  else if ( lda < max ( 1, n ) ) then
    info = - 5
    return
  else if ( ldb < max ( 1, n ) ) then
    info = - 8
    return
  end if

  if ( n == 0 .or. nrhs == 0 ) then
    return
  end if

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Solve L * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n - 1
        b(j+1:n,k) = b(j+1:n,k) - a(j+1:n,j) * b(j,k)
      end do
    end do
!
!  Solve U * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 1, -1
        b(j,k) = b(j,k) / a(j,j)
        b(1:j-1,k) = b(1:j-1,k) - a(1:j-1,j) * b(j,k)
      end do
    end do

  else
!
!  Solve Transpose ( U ) * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n
        b(j,k) = b(j,k) / a(j,j)
        do i = j + 1, n
          b(i,k) = b(i,k) - a(j,i) * b(j,k)
        end do
      end do
    end do
!
!  Solve Transpose ( L ) * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 2, -1
        do i = 1, j - 1
          b(i,k) = b(i,k) - a(j,i) * b(j,k)
        end do
      end do
    end do

  end if

  return
end
subroutine sge_plu ( lda, m, n, a, p, l, u )
!
!*******************************************************************************
!
!! SGE_PLU produces the PLU factors of a real rectangular matrix.
!
!
!  Discussion:
!
!    The PLU factors of the M by N matrix A are:
!
!      P, an M by M permutation matrix P,
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix.
!
!  Modified:
!
!    30 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real A(LDA,N), the M by N matrix to be factored.
!
!    Output, real P(LDA,M), the M by M permutation factor.
!
!    Output, real L(LDA,M), the M by M unit lower triangular factor.
!
!    Output, real U(LDA,N), the M by N upper triangular factor.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer k
  real l(lda,m)
  real p(lda,m)
  integer pivot_row
  real pivot_value
  real u(lda,n)
!
!  Initialize:
!
!    P: = M by M Identity
!    L: = M by M Identity
!    U: = A
!
  call sge_identity ( lda, m, l )
  call sge_identity ( lda, m, p )

  u(1:m,1:n) = a(1:m,1:n)
!
!  On step J, find the pivot row and the pivot value.
!
  do j = 1, min ( m-1, n )

    pivot_value = 0.0E+00
    pivot_row = 0

    do i = j, m

      if ( abs ( u(i,j) ) > pivot_value ) then
        pivot_value = abs ( u(i,j) )
        pivot_row = i
      end if

    end do
!
!  If the pivot row is nonzero, swap rows J and PIVOT_ROW.
!
    if ( pivot_row /= 0 ) then

      call rrow_swap ( lda, m, n, u, j, pivot_row )

      call rrow_swap ( lda, m, m, l, j, pivot_row )

      call rcol_swap ( lda, m, m, l, j, pivot_row )

      call rcol_swap ( lda, m, m, p, j, pivot_row )
!
!  Zero out the entries in column J, from row J+1 to M.
!
      do i = j+1, m

        if ( u(i,j) /= 0.0E+00 ) then

          l(i,j) = u(i,j) / u(j,j)
          u(i,j) = 0.0E+00
          u(i,j+1:n) = u(i,j+1:n) - l(i,j) * u(j,j+1:n)

        end if

      end do

    end if

  end do

  return
end
subroutine sge_poly ( lda, n, a, p )
!
!*******************************************************************************
!
!! SGE_POLY computes the characteristic polynomial of a general matrix.
!
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the N by N matrix, stored in LINPACK
!    general matrix storage.
!
!    Output, real P(0:N), the coefficients of the characteristic
!    polynomial of A.  P(I) contains the coefficient of X**I.
!
  integer lda
  integer n
!
  real a(lda,n)
  real diag(n)
  integer i
  integer ierror
  integer iorder
  integer j
  integer k
  real p(0:n)
  real trace
  real work1(n,n)
  real work2(n,n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_POLY - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if
!
!  Initialize WORK1 to the identity matrix.
!
  call sge_identity ( n, n, work1 )

  p(n) = 1.0E+00

  do iorder = n-1, 0, -1
!
!  Work2 = A * WORK1.
!
    work2(1:n,1:n) = matmul ( a(1:n,1:n), work1(1:n,1:n) )
!
!  Take the trace.
!
    call rmat_diag_get_vector ( n, n, work2, diag )

    trace = sum ( diag(1:n) )
!
!  P(IORDER) = - Trace ( WORK2 ) / ( N - IORDER )
!
    p(iorder) = - trace / real ( n - iorder )
!
!  WORK1 := WORK2 + P(IORDER) * Identity.
!
    work1(1:n,1:n) = work2(1:n,1:n)

    call rmat_diag_add_scalar ( n, n, work1, p(iorder) )

  end do

  return
end
subroutine sge_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! SGE_PRINT prints a general matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N matrix, stored in LINPACK
!    or LAPACK general storage mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer m
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  call sge_print_some ( lda, m, n, a, 1, 1, m, n )

  return
end
subroutine sge_print_some ( lda, m, n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SGE_PRINT_SOME prints some of a general matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N matrix, stored in LINPACK
!    or LAPACK general storage mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
!
  real a(lda,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer m
  logical r_is_int
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if (  r_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sge_random ( lda, m, n, a )
!
!*******************************************************************************
!
!! SGE_RANDOM randomizes a general matrix.
!
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, real A(LDA,N), the randomized M by N matrix, with entries
!    between 0 and 1.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  call rmat_random ( 0.0E+00, 1.0E+00, lda, m, n, a )

  return
end
subroutine sge_res ( lda, n, a, b, job, x, r )
!
!*******************************************************************************
!
!! SGE_RES computes the residual vector for a linear system.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the original, UNFACTORED matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input, integer JOB, specifies the linear system being solved:
!    0, A * x = b;
!    nonzero, A' * x = b.
!
!    Input, real X(N), an estimate of the solution the linear system.
!
!    Output, real R(N), the residual vector:
!      b - A * x
!    or
!      b - A' * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  integer job
  real r(n)
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_RES - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  if ( job == 0 ) then
    r(1:n) = b(1:n) - matmul ( a(1:n,1:n), x(1:n) )
  else
    r(1:n) = b(1:n) - matmul ( transpose ( a(1:n,1:n) ), x(1:n) )
  end if

  return
end
subroutine sge_sl ( lda, n, a, pivot, b, job )
!
!*******************************************************************************
!
!! SGE_SL solves a system factored by SGE_FA.
!
!
!  Discussion:
!
!    SGE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the LU factors from SGE_FA.
!
!    Input, integer PIVOT(N), the pivot vector from SGE_FA.
!
!    Input/output, real B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  integer pivot(n)
  integer j
  integer job
  integer k
  integer l
  real t
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n-1

      l = pivot(k)

      if ( l /= k ) then
        call r_swap ( b(l), b(k) )
      end if

      b(k+1:n) = b(k+1:n) + a(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      b(1:k-1) = b(1:k-1) - a(1:k-1,k) * b(k)
    end do
!
!  Solve transpose ( A ) * X = B.
!
  else
!
!  Solve transpose ( U ) * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - dot_product ( b(1:k-1), a(1:k-1,k) ) ) / a(k,k)
    end do
!
!  Solve transpose ( PL ) * X = Y.
!
    do k = n-1, 1, -1

      b(k) = b(k) + dot_product ( b(k+1:n), a(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        call r_swap ( b(l), b(k) )
      end if

    end do

  end if

  return
end
subroutine sge_sl_it ( lda, n, a, alu, pivot, b, job, x, r )
!
!*******************************************************************************
!
!! SGE_SL_IT applies one step of iterative refinement following SGE_SL.
!
!
!  Discussion:
!
!    It is assumed that:
!
!    * the original matrix A has been factored by SGE_FA;
!    * the linear system A * x = b has been solved once by SGE_SL.
!
!    (Actually, it is not necessary to solve the system once using SGE_SL.
!    You may simply supply the initial estimated solution X = 0.)
!
!    Each time this routine is called, it will compute the residual in
!    the linear system, apply one step of iterative refinement, and
!    add the computed correction to the current solution.
!
!  Modified:
!
!    15 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the original, UNFACTORED matrix.
!
!    Input, real ALU(LDA,N), the LU factors of the matrix from SGE_FA.
!
!    Input, integer PIVOT(N), the pivot vector from SGE_FA.
!
!    Input, real B(N), the right hand side vector.
!
!    Input, integer JOB, specifies the operation.
!    0, solve A*X=B.
!    nonzero, solve transpose(A)*X=B.
!
!    Input/output, real X(N), an estimate of the solution of A * x = b.
!    On output, the solution has been improved by one step of iterative
!    refinement.
!
!    Output, real R(N), contains the correction terms added to X.
!
  integer lda
  integer n
!
  real a(lda,n)
  real alu(lda,n)
  real b(n)
  integer i
  integer ierror
  integer pivot(n)
  integer job
  real r(n)
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_SL_IT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Compute the residual vector
!    r = b - A * x
!  or
!    r = b - Transpose ( A ) * x
!
  call sge_res ( lda, n, a, b, job, x, r )
!
!  Solve
!    A * dx = r
!  or
!    Transpose ( A ) * dx = r
!
  call sge_sl ( lda, n, alu, pivot, r, job )
!
!  Add dx to x.
!
  x(1:n) = x(1:n) + r(1:n)

  return
end
subroutine sge_trf ( lda, m, n, a, pivot, info )
!
!*******************************************************************************
!
!! SGE_TRF computes the PLU factorization of a general M by N matrix.
!
!
!  Note:
!
!    SGE_TRF is a standalone version of the LAPACK routine SGETRF.
!
!    The factorization uses partial pivoting with row interchanges,
!    and has the form
!      A = P * L * U
!    where P is a permutation matrix, L is lower triangular with unit
!    diagonal elements (lower trapezoidal if M > N), and U is upper
!    triangular (upper trapezoidal if M < N).
!
!  Modified:
!
!    17 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= max(1,M).
!
!    Input, integer M, the number of rows of the matrix A.  M >= 0.
!
!    Input, integer N, the number of columns of the matrix A.  N >= 0.
!
!    Input/output, real A(LDA,N).
!    On entry, the M by N matrix to be factored.
!    On exit, the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
!
!    Output, integer PIVOT(min(M,N)), the pivot indices;
!    for 1 <= I <= min(M,N), row i of the matrix was interchanged with
!    row PIVOT(I).
!
!    Output, integer INFO.
!    = 0: successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!    > 0: if INFO = K, U(K,K) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ii
  integer info
  integer pivot(*)
  integer j
  integer jj
  integer jp
  integer m
  real temp
!
!  Test the input parameters.
!
  info = 0

  if ( m < 0 ) then
    info = - 1
    return
  else if( n < 0 ) then
    info = - 2
    return
  else if ( lda < max ( 1, m ) ) then
    info = - 4
    return
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  do j = 1, min ( m, n )
!
!  Find the pivot.
!
    temp = abs ( a(j,j) )
    jp = j
    do i = j+1, m
      if ( abs ( a(i,j) ) > temp ) then
        temp = abs ( a(i,j) )
        jp = i
      end if
    end do

    pivot(j) = jp
!
!  Apply the interchange to columns 1:N.
!  Compute elements J+1:M of the J-th column.
!
    if ( a(jp,j) /= 0.0E+00 ) then

      if ( jp /= j ) then
        do jj = 1, n
          call r_swap ( a(j,jj), a(jp,jj) )
        end do
      end if

      if ( j < m ) then
        a(j+1:m,j) = a(j+1:m,j) / a(j,j)
      end if

    else if ( info == 0 ) then

      info = j

    end if
!
!  Update the trailing submatrix.
!
    if ( j < min ( m, n ) ) then

      do ii = j+1, m
        a(ii,j+1:n) = a(ii,j+1:n) - a(ii,j) * a(j,j+1:n)
      end do

    end if

  end do

  return
end
subroutine sge_trs ( lda, n, nrhs, trans, a, pivot, b, ldb, info )
!
!*******************************************************************************
!
!! SGE_TRS solves a system of linear equations factored by SGE_TRF.
!
!
!  Note:
!
!    SGE_TRS is a standalone version of the LAPACK routine SGETRS.
!
!    SGE_TRS solves a system of linear equations
!      A * x = b  or  A' * X = B
!    with a general N by N matrix A using the PLU factorization computed
!    by SGE_TRF.
!
!  Modified:
!
!    17 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA >= max(1,N).
!
!    Input, integer N, the order of the matrix A.  N >= 0.
!
!    Input, integer NRHS, the number of right hand sides.  NRHS >= 0.
!
!    Input, character TRANS, specifies the form of the system of equations:
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real A(LDA,N), the factors L and U from the factorization
!    A = P*L*U as computed by SGE_TRF.
!
!    Input, integer PIVOT(N), the pivot indices from SGE_TRF;
!    for 1<=i<=N, row i of the matrix was interchanged with row PIVOT(I).
!
!    Input/output, real B(LDB,NRHS).
!    On entry, the right hand side matrix B.
!    On exit, the solution matrix X.
!
!    Input, integer LDB, the leading dimension of the array B.
!    LDB >= max(1,N).
!
!    Output, integer INFO
!    = 0:  successful exit
!    < 0:  if INFO = -I, the I-th argument had an illegal value.
!
  integer lda
  integer ldb
  integer n
  integer nrhs
!
  real a(lda,n)
  real b(ldb,nrhs)
  integer i
  integer info
  integer pivot(n)
  integer j
  integer k
  real temp
  character trans
!
  info = 0

  if ( trans /= 'n' .and. trans /= 'N' .and. &
       trans /= 't' .and. trans /= 'T' .and. &
       trans /= 'c' .and. trans /= 'C' ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  else if ( nrhs < 0 ) then
    info = - 3
    return
  else if ( lda < max ( 1, n ) ) then
    info = - 5
    return
  else if ( ldb < max ( 1, n ) ) then
    info = - 8
    return
  end if

  if ( n == 0 .or. nrhs == 0 ) then
    return
  end if

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Apply row interchanges to the right hand sides.
!
    do i = 1, n
      if ( pivot(i) /= i ) then
        do k = 1, nrhs
          call r_swap ( b(i,k), b(pivot(i),k) )
        end do
      end if
    end do
!
!  Solve L * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n - 1
        b(j+1:n,k) = b(j+1:n,k) - a(j+1:n,j) * b(j,k)
      end do
    end do
!
!  Solve U * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 1, -1
        b(j,k) = b(j,k) / a(j,j)
        b(1:j-1,k) = b(1:j-1,k) - a(1:j-1,j) * b(j,k)
      end do
    end do

  else
!
!  Solve Transpose ( U ) * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n
        b(j,k) = b(j,k) / a(j,j)
        b(j+1:n,k) = b(j+1:n,k) - a(j,j+1:n) * b(j,k)
      end do
    end do
!
!  Solve Transpose ( L ) * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 2, -1
        b(1:j-1,k) = b(1:j-1,k) - a(j,1:j-1) * b(j,k)
      end do
    end do
!
!  Apply row interchanges to the solution vectors.
!
    do i = n, 1, -1
      if ( pivot(i) /= i ) then
        do k = 1, nrhs
          call r_swap ( b(i,k), b(pivot(i),k) )
        end do
      end if
    end do

  end if

  return
end
subroutine sge_vxm ( lda, m, n, a, x, b )
!
!*******************************************************************************
!
!! SGE_VXM computes Tranpose(A) * X, where A is a general matrix.
!
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N matrix, stored in LINPACK
!    general matrix storage.
!
!    Input, real X(M), the vector to be multiplied by A.
!
!    Output, real B(N), the product Tranpose ( A ) * x.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  real x(m)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine sge_zero ( lda, m, n, a )
!
!*******************************************************************************
!
!! SGE_ZERO zeroes out a general matrix.
!
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, real A(LDA,N), the M by N matrix.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGE_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  a(1:m,1:n) = 0.0E+00

  return
end
subroutine slt_det ( lda, n, a, det )
!
!*******************************************************************************
!
!! SLT_DET computes the determinant of a lower triangular matrix.
!
!
!  Modified:
!
!    22 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the lower triangular matrix.
!
!    Output, real DET, the determinant of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  real diag(n)
  integer ierror
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SLT_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  call rmat_diag_get_vector ( lda, n, a, diag )

  det = product ( diag(1:n) )

  return
end
subroutine slt_inv ( lda, n, a )
!
!*******************************************************************************
!
!! SLT_INV computes the inverse of a lower triangular matrix.
!
!
!  Reference:
!
!    Combinatorial Algorithms,
!    A Nijenhuis and H Wilf,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6
!
!  Modified:
!
!    22 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(LDA,N).
!
!    On input, the lower triangular matrix to be inverted.
!    On output, the inverse of the lower triangular matrix.
!
  integer n
  integer lda
!
  real a(lda,n)
  integer i
  integer j
  integer k
!
!  Check.
!
  do i = 1, n
    if ( a(i,i) == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SLT_INV - Fatal error!'
      write ( *, * ) '  Zero diagonal element.'
      stop
    end if
  end do

  do j = 1, n

    do i = 1, n

      if ( i < j ) then

        a(i,j) = 0.0E+00

      else if ( i == j ) then

        a(i,j) = 1.0E+00 / a(i,j)

      else if ( i > j ) then

        a(i,j) = - dot_product ( a(i,j:i-1), a(j:i-1,j) ) / a(i,i)

      end if

    end do
  end do

  return
end
subroutine slt_mxv ( lda, m, n, a, x, b )
!
!*******************************************************************************
!
!! SLT_MXV computes A * x, where A is a lower triangular matrix.
!
!
!  Modified:
!
!    22 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N lower triangular matrix, stored
!    in LINPACK general matrix storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(M), the product A * x.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(m)
  integer i
  integer ierror
  integer jmax
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SLT_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, m
    jmax = min ( i, n )
    b(i) = dot_product ( a(i,1:jmax), x(1:jmax) )
  end do

  return
end
subroutine slt_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! SLT_PRINT prints a lower triangular matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N band matrix, stored in LINPACK
!    or LAPACK general band storage mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer m
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call slt_print_some ( lda, m, n, a, 1, 1, m, n )

  return
end
subroutine slt_print_some ( lda, m, n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SLT_PRINT_SOME prints some of a lower triangular matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N band matrix, stored in LINPACK
!    or LAPACK general band storage mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
!
  real a(lda,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer m
  logical r_is_int
!
  if ( jlo > ilo ) then
    return
  end if
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SLT_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo )

    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j > i ) then
          ctemp(j2) = '              '
        else if ( r_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine slt_sl (lda, n, a, b )
!
!*******************************************************************************
!
!! SLT_SL solves a lower triangular system.
!
!
!  Discussion:
!
!    No factorization of the lower triangular matrix is required.
!
!  Modified:
!
!    22 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(LDA,N), the lower triangular matrix.
!
!    Input/output, real B(N).
!
!    On input, the right hand side.
!    On output, the solution vector.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer j
!
  do j = 1, n
    b(j) = b(j) / a(j,j)
    b(j+1:n) = b(j+1:n) - a(j+1:n,j) * b(j)
  end do

  return
end
subroutine slt_vxm ( lda, m, n, a, x, b )
!
!*******************************************************************************
!
!! SLT_VXM computes A' * x, where A is a lower triangular matrix.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N lower triangular matrix, stored
!    in LINPACK general matrix storage.
!
!    Input, real X(M), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  real x(m)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SLT_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, n
    b(i) = dot_product ( x(i:m), a(i:m,i) )
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )
!
!*******************************************************************************
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Modified:
!
!    12 November 2000
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I precedes J, ISGN = +1 if J precedes I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I precedes J;
!    ISGN => 0 means J precedes I.
!
  integer i
  integer indx
  integer isgn
  integer j
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( isgn > 0 ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  INDX > 0, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

  return
end
subroutine spb_cg ( lda, n, mu, a, b, x )
!
!*******************************************************************************
!
!! SPB_CG uses the conjugate gradient method on a symmetric banded system.
!
!
!  Discussion:
!
!    The matrix A must be a positive definite symmetric band matrix.
!    To save storage, A is stored in a compact diagonal format.
!
!    The method is designed to reach the solution after N computational
!    steps.  However, roundoff may introduce unacceptably large errors for
!    some problems.  In such a case, calling the routine again, using
!    the computed solution as the new starting estimate, should improve
!    the results.
!
!  Reference:
!
!    F S Beckman,
!    The Solution of Linear Equations by the Conjugate Gradient Method,
!    Mathematical Methods for Digital Computers, pages 62-72.
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals.
!    MU must be at least 0, and no more than N-1.
!
!    Input, real A(LDA,N), the N by N matrix, stored in LINPACK positive
!    definite symmetric band matrix storage.
!
!    The diagonal is stored in row MU+1 of the array.
!
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.
!
  integer lda
  integer n
!
  real a(lda,n)
  real alpha
  real ap(n)
  real b(n)
  real beta
  integer i
  integer ierror
  integer it
  integer mu
  real p(n)
  real pap
  real pr
  real r(n)
  real rap
  real x(n)
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_CG - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call spb_mxv ( lda, n, mu, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call spb_mxv ( lda, n, mu, a, p, ap )
!
!  Compute the dot products
!    PAP = P*AP,
!    PR  = P*R
!  Set
!    ALPHA = PR / PAP.
!
    pap = dot_product ( p(1:n), ap(1:n) )
    pr = dot_product ( p(1:n), r(1:n) )

    if ( pap == 0.0E+00 ) then
      return
    end if

    alpha = pr / pap
!
!  Set
!    X = X + ALPHA * P
!    R = R - ALPHA * AP.
!
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * ap(1:n)
!
!  Compute the vector dot product
!    RAP = R*AP
!  Set
!    BETA = - RAP / PAP.
!
    rap = dot_product ( r(1:n), ap(1:n) )

    beta = - rap / pap
!
!  Update the perturbation vector
!    P = R + BETA * P.
!
    p(1:n) = r(1:n) + beta * p(1:n)

  end do

  return
end
subroutine spb_check ( lda, n, mu, ierror )
!
!*******************************************************************************
!
!! SPB_CHECK checks the dimensions of a positive definite symmetric band matrix.
!
!
!  Modified:
!
!    14 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the upper bandwidth of the matrix.
!    MU must be at least 0, and no greater than N-1.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if MU is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  integer ierror
  integer lda
  integer mu
  integer n
!
  ierror = 0

  if ( lda < mu + 1 ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'SPB_CHECK - Illegal LDA = ', lda
  end if

  if ( mu < 0 .or. mu > n - 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SPB_CHECK - Illegal MU = ', mu
  end if

  if ( n <= 0 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SPB_CHECK - Illegal N = ', n
    return
  end if

  return
end
subroutine spb_det ( lda, n, mu, a, det )
!
!*******************************************************************************
!
!! SPB_DET computes the determinant of a matrix factored by SPB_FA.
!
!
!  Modified:
!
!    29 October 1998
!
!  Reference:
!
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User's Guide,
!    SIAM, Philadelphia, 1979.
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real A(LDA,N), the matrix, as factored by SPB_FA.
!
!    Output, real DET, the determinant of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  integer i
  integer ierror
  integer mu
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  det = product ( a(mu+1,1:n)**2 )

  return
end
subroutine spb_fa ( lda, n, mu, a, info )
!
!*******************************************************************************
!
!! SPB_FA factors a positive definite symmetric band matrix A.
!
!
!  Discussion:
!
!    The matrix is stored in a compact form.
!
!    Once factored, linear systems A*x=b involving the matrix can be solved
!    by calling SPB_SL.  No pivoting is performed.  Pivoting is not necessary
!    for positive definite symmetric matrices.  If the matrix is not positive
!    definite, the algorithm may behave correctly, but it is also possible
!    that an illegal divide by zero will occur.
!
!  Modified:
!
!    31 October 1998
!
!  Reference:
!
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User's Guide,
!    SIAM, Philadelphia, 1979.
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0, and no more than N-1.
!
!    Input/output, real A(LDA,N), the N by N matrix, stored in LINPACK
!    positive definite symmetric band matrix storage.
!
!    The diagonal is stored in row MU+1 of the array.
!
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    On output, A contains information describing a factored form
!    of the matrix, that can be used to solve linear systems
!    A*x=b, using SPB_SL.
!
!    Output, integer INFO, singularity flag.
!    0, the matrix is nonsingular.
!    nonzero, the matrix is singular.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer ik
  integer info
  integer j
  integer jk
  integer k
  integer mm
  integer mu
  real s
  real t
  real temp
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_FA - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do j = 1, n

    s = 0.0E+00
    ik = mu + 1
    jk = max ( j - mu, 1 )
    mm = max ( mu + 2 - j, 1 )

    s = 0.0E+00

    do k = mm, mu

      a(k,j) = ( a(k,j) - dot_product ( a(ik:ik+k-mm-1,jk), a(mm:k-1,j) ) ) &
        / a(mu+1,jk)

      s = s + a(k,j)**2

      ik = ik - 1
      jk = jk + 1

    end do

    s = a(mu+1,j) - s

    if ( s <= 0.0E+00 ) then
      info = j
      write ( *, * ) ' '
      write ( *, * ) 'SPB_FA - Fatal error!'
      write ( *, * ) '  Nonpositive pivot on step ', info
      return
    end if

    a(mu+1,j) = sqrt ( s )

  end do

  return
end
subroutine spb_ml ( lda, n, mu, a, x, b )
!
!*******************************************************************************
!
!! SPB_ML multiplies a vector times a matrix that was factored by SPB_FA.
!
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real A(LDA,N), the matrix, as factored by SPB_FA.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer ilo
  integer j
  integer jhi
  integer k
  integer mu
  real x(n)
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  b(1:n) = x(1:n)
!
!  Multiply U * X = Y.
!
  do k = 1, n

    ilo = max ( 1, k - mu )
    do i = ilo, k - 1
      b(i) = b(i) + a(mu+1+i-k,k) * b(k)
    end do

    b(k) = a(mu+1,k) * b(k)

  end do
!
!  Multiply L * Y = B.
!
  do k = n, 1, -1

    jhi = min ( k + mu, n )
    do j = k + 1, jhi
      b(j) = b(j) + a(mu+1+k-j,j) * b(k)
    end do

    b(k) = a(mu+1,k) * b(k)

  end do

  return
end
subroutine spb_mxv ( lda, n, mu, a, x, b )
!
!*******************************************************************************
!
!! SPB_MXV multiplies a positive definite symmetric band matrix times a vector.
!
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real A(LDA,N), the matrix, stored in LINPACK positive
!    definite symmetric band storage.
!
!    The diagonal is stored in row MU+1 of the array.
!
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the result vector A * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ieqn
  integer ierror
  integer j
  integer mu
  real x(n)
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Multiply X by the diagonal of the matrix.
!
  b(1:n) = a(mu+1,1:n) * x(1:n)
!
!  Multiply X by the superdiagonals of the matrix.
!
  do i = mu, 1, -1
    do j = mu+2-i, n
      ieqn = i + j - mu - 1
      b(ieqn) = b(ieqn) + a(i,j) * x(j)
      b(j) = b(j) + a(i,j) * x(ieqn)
    end do
  end do

  return
end
subroutine spb_print ( lda, n, mu, a, title )
!
!*******************************************************************************
!
!! SPB_PRINT prints a symmetric banded matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the N by N band matrix, stored in LINPACK
!    or LAPACK symmetric band storage mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer mu
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call spb_print_some ( lda, n, mu, a, 1, 1, n, n )

  return
end
subroutine spb_print_some ( lda, n, mu, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SPB_PRINT_SOME prints some of a symmetric banded matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real A(LDA,N), the N by N band matrix, stored in LINPACK
!    or LAPACK symmetric band storage mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
!
  real a(lda,n)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer mu
  logical r_is_int
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + mu )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j .and. j <= i + mu ) then
          aij = a(mu+1+i-j,j)
        else if ( i - mu <= j .and. j <= i ) then
          aij = a(mu+1+j-i,i)
        else
          aij = 0.0E+00
        end if

        if ( i-j > mu .or. j-i > mu ) then
          ctemp(j2) = '              '
        else if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine spb_random ( lda, n, mu, a )
!
!*******************************************************************************
!
!! SPB_RANDOM randomizes a positive definite symmetric band matrix.
!
!
!  Note:
!
!    The matrix returned will be positive definite, but of limited
!    randomness.  The off diagonal elements are random values between
!    0 and 1, and the diagonal element of each row is selected to
!    ensure strict diagonal dominance.
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Output, real A(LDA,N), the N by N matrix, stored in LINPACK positive
!    definite symmetric band matrix storage.
!
!    The diagonal is stored in row MU+1 of the array.
!
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer jhi
  integer jlo
  integer mu
  real r
  real sum2
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Zero out the "junk" entries.
!
  do j = 1, mu
    a(1:mu+1-j,j) = 0.0E+00
  end do
!
!  Set the off diagonal values.
!
  do i = 1, n
    do j = i+1, min ( i+mu, n )
      call r_random ( 0.0E+00, 1.0E+00, a(mu+1+i-j,j) )
    end do
  end do
!
!  Set the diagonal values.
!
  do i = 1, n

    sum2 = 0.0E+00

    jlo = max ( 1, i - mu )
    do j = jlo, i-1
      sum2 = sum2 + abs ( a(mu+1+j-i,i) )
    end do

    jhi = min ( i + mu, n )
    do j = i+1, jhi
      sum2 = sum2 + abs ( a(mu+1+i-j,j) )
    end do

    call r_random ( 0.0E+00, 1.0E+00, r )

    a(mu+1,i) = ( 1.0E+00 + r ) * ( sum2 + 0.01E+00 )

  end do

  return
end
subroutine spb_sl ( lda, n, mu, a, b )
!
!*******************************************************************************
!
!! SPB_SL solves a linear system A * x = b, factored by SPB_FA.
!
!
!  Modified:
!
!    31 October 1998
!
!  Reference:
!
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User's Guide,
!    SIAM, Philadelphia, 1979.
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real A(LDA,N), the matrix, as factored by SPB_FA.
!
!    Input/output, real B(N).
!
!    On input, B contains the right hand side of the linear system
!    to be solved.
!
!    On output, B contains X, the solution vector.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer ilo
  integer k
  integer mu
  real t
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Solve L * Y = B.
!
  do k = 1, n
    ilo = max ( 1, k - mu )
    b(k) = ( b(k) - dot_product ( b(ilo:k-1), a(mu+1+ilo-k:mu,k) ) ) &
      / a(mu+1,k)
  end do
!
!  Solve U * X = Y.
!
  do k = n, 1, -1

    b(k) = b(k) / a(mu+1,k)

    ilo = max ( 1, k - mu )
    do i = ilo, k - 1
      b(i) = b(i) - b(k) * a(mu+1+i-k,k)
    end do

  end do

  return
end
subroutine spb_sor ( lda, n, mu, a, b, eps, itchk, itknt, itmax, omega, x )
!
!*******************************************************************************
!
!! SPB_SOR uses SOR iteration to solve the PDS band system A*x=b.
!
!
!  Discussion:
!
!    A is a positive definite symmetric band matrix stored in a
!    compact format.  A relaxation factor OMEGA may be used.
!    The iteration will proceed until a convergence test is met,
!    or the iteration limit is reached.
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals in the matrix.
!    MU must be at least 0, and no more than N-1.
!
!    Input, real A(LDA,N), the N by N matrix, stored in LINPACK positive
!    definite symmetric band matrix storage.
!
!    The diagonal is stored in row MU+1 of the array.
!
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    Input, real B(N), the right hand side of the system.
!
!    Input, real EPS, convergence tolerance for the system.  The vector
!    b - A * x is computed every ITCHK iterations, and if the maximum
!    entry of this vector is of norm less than EPS, the program
!    will return.
!
!    Input, integer ITCHK, the interval between convergence checks.  ITCHK steps
!    will be taken before any check is made on whether the iteration
!    has converged.  ITCHK should be at least 1 and no greater
!    than ITMAX.
!
!    Output, integer ITKNT, the number of iterations taken.
!
!    Input, integer ITMAX, the maximum number of iterations allowed.  The
!    program will return to the user if this many iterations are taken
!    without convergence.
!
!    Input, real OMEGA, the relaxation factor.  OMEGA must be strictly between
!    0 and 2.  Use OMEGA = 1 for no relaxation, classical Jacobi iteration.
!
!    Input/output, real X(N).
!
!    On input, a starting vector for the iteration.
!
!    On output, the current approximation to the solution.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  real eps
  real err
  integer i
  integer ierror
  integer it
  integer itchk
  integer itknt
  integer itmax
  integer mu
  real omega
  real x(n)
  real xtemp(n)
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_SOR - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  if ( itchk <= 0 .or. itchk > itmax ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_SOR - Fatal error!'
    write ( *, * ) '  Illegal ITCHK= ', itchk
    return
  end if

  if ( itmax <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_SOR - Fatal error!'
    write ( *, * ) '  Nonpositive ITMAX =', itmax
    return
  end if

  if ( omega <= 0.0E+00 .or. omega >= 2.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_SOR - Fatal error!'
    write ( *, * ) '  Illegal value of OMEGA = ', omega
    return
  end if

  itknt = 0
!
!  Take ITCHK steps of the iteration before doing a convergence check.
!
  do while ( itknt <= itmax )

    do it = 1, itchk
!
!  Compute XTEMP(I) = B(I) + A(I,I) * X(I) - SUM(J=1 to N)A(I,J) * X(J).
!
      call spb_mxv ( lda, n, mu, a, x, xtemp )

      xtemp(1:n) = x(1:n) + ( b(1:n) - xtemp(1:n) ) / a(mu+1,1:n)
!
!  Compute the next iterate as a weighted combination of the
!  old iterate and the just computed standard Jacobi iterate.
!
      if ( omega /= 1.0E+00 ) then
        xtemp(1:n) = ( 1.0E+00 - omega ) * x(1:n) + omega * xtemp(1:n)
      end if

      itknt = itknt + 1
!
!  Copy the new result into the old result vector.
!
      x(1:n) = xtemp(1:n)

    end do
!
!  Compute the maximum residual, the greatest entry in the vector
!  RESID(I) = B(I) - A(I,J) * X(J).
!
    call spb_mxv ( lda, n, mu, a, x, xtemp )

    err = maxval ( abs ( b(1:n) - xtemp(1:n) ) )
!
!  Test to see if
!    we can quit because of convergence,
!
    if ( err <= eps ) then
      return
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'SPB_SOR - Warning!'
  write ( *, * ) '  The iteration did not converge.'

  return
end
subroutine spb_to_sge ( lda1, lda2, n, mu, a1, a2 )
!
!*******************************************************************************
!
!! SPB_TO_SGE converts a positive definite symmetric band matrix to general matrix format.
!
!
!  Modified:
!
!    14 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA1, the leading dimension of the array A1.
!    LDA1 must be at least MU+1.
!
!    Input, integer LDA2, the leading dimension of the array A2.
!    LDA2 must be at least N.
!
!    Input, integer N, the order of the matrices.
!    N must be positive.
!
!    Input, integer MU, the upper bandwidth of A1.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real A1(LDA1,N), the positive definite symmetric band matrix.
!
!    Output, real A2(LDA2,N), the general matrix, which contains the
!    information given in A1.
!
  integer lda1
  integer lda2
  integer n
!
  real a1(lda1,n)
  real a2(lda2,n)
  integer i
  integer ierror
  integer j
  integer mu
!
!  Check the dimensions.
!
  call spb_check ( lda1, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A1.'
    return
  end if

  call sge_check ( lda2, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A2.'
    return
  end if

  do i = 1, n
    do j = 1, n
      if ( i <= j .and. j <= i+mu ) then
        a2(i,j) = a1(mu+1+i-j,j)
      else if ( i-mu <= j .and. j < i ) then
        a2(i,j) = a1(mu+1+j-i,i)
      else
        a2(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine spb_zero ( lda, n, mu, a )
!
!*******************************************************************************
!
!! SPB_ZERO zeroes out a positive definite symmetric band matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Output, real A(LDA,N), the N by N matrix, stored in LINPACK positive
!    definite symmetric band matrix storage.
!
!    The diagonal is stored in row MU+1 of the array.
!
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer mu
!
!  Check the dimensions.
!
  call spb_check ( lda, n, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPB_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  a(1:lda,1:n) = 0.0E+00

  return
end
subroutine spo_det ( lda, n, a, det )
!
!*******************************************************************************
!
!! SPO_DET computes the determinant of an SPD matrix factored by SPO_FA.
!
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(LDA,N), the factor information returned by SPO_FA.
!
!    Output, real DET, the determinant of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  real diag(n)
!
  call rmat_diag_get_vector ( lda, n, a, diag )

  diag = product ( diag(1:n)**2 )

  return
end
subroutine spo_fa ( lda, n, a, info )
!
!*******************************************************************************
!
!! SPO_FA factors a real symmetric positive definite matrix.
!
!
!  Discussion:
!
!    The positive definite symmetric matrix A has a Cholesky factorization
!    of the form:
!
!      A = Transpose ( R ) * R
!
!    where R is an upper triangular matrix with positive elements on
!    its diagonal.  This routine overwrites the matrix A with its
!    factor R.
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(LDA,N).
!
!    On input, the N by N positive definite symmetric matrix.
!    On output, the upper triangular Cholesky factor.
!
!    Output, integer INFO, error flag.
!    0, normal return.
!    K, error condition.  The principal minor of order K is not
!    positive definite, and the factorization was not completed.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer info
  integer j
  integer k
  real s
  real t
!
  do j = 1, n

    s = 0.0E+00

    do k = 1, j - 1

      t = a(k,j)
      do i = 1, k - 1
        t = t - a(i,k) * a(i,j)
      end do
      a(k,j) = t / a(k,k)
      s = s + a(k,j)**2
    end do

    s = a(j,j) - s

    if ( s <= 0.0E+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( s )

  end do
!
!  Zero out the strict lower triangle.
!
  do i = 2, n
    a(i,1:i-1) = 0.0E+00
  end do

  info = 0

  return
end
subroutine spo_inv ( lda, n, a )
!
!*******************************************************************************
!
!! SPO_INV computes the inverse of an SPD matrix factored by SPO_FA.
!
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(LDA,N).
!
!    On input, A contains the factor information returned by SPO_FA.
!    On output, A contains the inverse matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer k
  real t
!
!  Compute Inverse ( R ).
!
  do k = 1, n

    a(k,k) = 1.0E+00 / a(k,k)
    a(1:k-1,k) = - a(1:k-1,k) * a(k,k)

    do j = k + 1, n
      t = a(k,j)
      a(k,j) = 0.0E+00
      a(1:k,j) = a(1:k,j) + t * a(1:k,k)
    end do

  end do
!
!  Compute Inverse ( R ) * Transpose ( Inverse ( R ) ).
!
  do j = 1, n

    do k = 1, j - 1
      t = a(k,j)
      a(1:k,k) = a(1:k,k) + t * a(1:k,j)
    end do

    a(1:j,j) = a(1:j,j) * a(j,j)

  end do
!
!  Copy upper triangle into lower triangle.
!
  do i = 2, n
    a(i,1:i-1) = a(1:i-1,i)
  end do

  return
end
subroutine spo_ml ( lda, n, a, x, b )
!
!*******************************************************************************
!
!! SPO_ML computes A * x = b after A has been factored by SPO_FA.
!
!
!  Modified:
!
!    05 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(LDA,N). the factor information returned by SPO_FA.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer j
  real x(n)
!
!  Compute R * x = y.
!
  do i = 1, n
    b(i) = a(i,i) * x(i) + dot_product ( a(i,i+1:n), x(i+1:n) )
  end do
!
!  Compute Tranpose ( R ) * y = b.
!
  do i = n, 1, -1
    b(i) = a(i,i) * b(i) + dot_product ( b(1:i-1), a(1:i-1,i) )
  end do

  return
end
subroutine spo_random ( lda, n, a )
!
!*******************************************************************************
!
!! SPO_RANDOM randomizes a positive definite symmetric matrix.
!
!
!  Note:
!
!    The matrix is computed by setting a "random" upper triangular
!    Cholesky factor R, and then computing A = transpose(R)*R.
!    The randomness is limited by the fact that all the entries of
!    R will be between 0 and 1.  A truly random R is only required
!    to have positive entries on the diagonal.
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real A(LDA,N), the N by N matrix, stored in LINPACK general
!    storage.  The matrix should be symmetric and positive definite.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  integer k
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPO_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  a(1:n,1:n) = 0.0E+00

  do i = n, 1, -1
!
!  Set row I of R.
!
    do j = i, n
      call r_random ( 0.0E+00, 1.0E+00, a(i,j) )
    end do
!
!  Consider element J of row I, last to first.
!
    do j = n, i, -1
!
!  Add multiples of row I to lower elements of column J.
!
      a(i+1:j,j) = a(i+1:j,j) + a(i,i+1:j) * a(i,j)
!
!  Reset element J.
!
      a(i,j) = a(i,i) * a(i,j)

    end do
  end do
!
!  Now copy the upper triangle to the lower triangle.
!
  do i = 1, n
    do j = 1, i-1
      a(i,j) = a(j,i)
    end do
  end do

  return
end
subroutine spo_sl ( lda, n, a, b )
!
!*******************************************************************************
!
!! SPO_SL solves an SPD system factored by SPO_FA.
!
!
!  Modified:
!
!    04 March 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(LDA,N). the factor information returned by SPO_FA.
!
!    Input/output, real B(N).
!
!    On input, the right hand side.
!    On output, the solution vector.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer k
  real t
!
!   Solve Transpose ( R ) * y = b.
!
  do k = 1, n
    b(k) = ( b(k) - dot_product ( b(1:k-1), a(1:k-1,k) ) ) / a(k,k)
  end do
!
!  Solve R * x = y.
!
  do k = n, 1, -1
    b(k) = b(k) / a(k,k)
    b(1:k-1) = b(1:k-1) - a(1:k-1,k) * b(k)
  end do

  return
end
subroutine spp_print ( n, a, title )
!
!*******************************************************************************
!
!! SPP_PRINT prints a square packed matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(((N+1)*N)/2), the N by N matrix, stored in LINPACK
!    or LAPACK positive definite symmetric packed mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  real a((n*(n+1))/2)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call spp_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine spp_print_some ( n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SPP_PRINT_SOME prints some of a square packed matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(((N+1)*N)/2), the N by N matrix, stored in LINPACK
!    or LAPACK positive definite symmetric packed mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a((n*(n+1))/2)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j >= i ) then
          aij = a(i+(j*(j-1))/2)
        else
          aij = a(j+(i*(i-1))/2)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine spp_random ( n, a )
!
!*******************************************************************************
!
!! SPP_RANDOM randomizes a positive definite symmetric packed matrix.
!
!
!  Note:
!
!    The matrix is computed by setting a "random" upper triangular
!    Cholesky factor R, and then computing A = transpose(R)*R.
!    The randomness is limited by the fact that all the entries of
!    R will be between 0 and 1.  A truly random R is only required
!    to have positive entries on the diagonal.
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real A((N*(N+1))/2), the N by N matrix, stored in LINPACK
!    symmetric packed storage.  The matrix should be positive definite.
!    The entry (I,J) is stored in position I+(J*(J-1))/2
!
  integer n
!
  real a((n*(n+1))/2)
  integer i
  integer ii
  integer ij
  integer ik
  integer j
  integer k
  integer kj
!
  a(1:(n*(n+1))/2) = 0.0E+00

  do i = n, 1, -1
!
!  Set row I of R.
!
    do j = i, n
      ij = i + ( j * ( j - 1 ) ) / 2
      call r_random ( 0.0E+00, 1.0E+00, a(ij) )
    end do
!
!  Consider element J of row I, last to first.
!
    do j = n, i, -1
!
!  Add multiples of row I to lower elements of column J.
!
      ij = i + ( j * ( j - 1 ) ) / 2

      do k = i+1, j
        kj = k + (j*(j-1))/2
        ik = i + (k*(k-1))/2
        a(kj) = a(kj) + a(ik) * a(ij)
      end do
!
!  Reset element J.
!
      ii = i + (i*(i-1))/2
      a(ij) = a(ii) * a(ij)

    end do
  end do

  return
end
subroutine spp_to_sge ( lda, n, a, a2 )
!
!*******************************************************************************
!
!! SPP_TO_SGE copies a packed matrix into a general matrix.
!
!
!  Modified:
!
!    13 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A((N*(N+1))/2), the packed matrix.
!
!    Output, real A2(LDA,N), the matrix, stored as a general matrix.
!
  integer lda
  integer n
!
  real a((n*(n+1))/2)
  real a2(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPP_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, n
      if ( j >= i ) then
        a2(i,j) = a(i+(j*(j-1))/2)
      else
        a2(i,j) = a(j+(i*(i-1))/2)
      end if
    end do
  end do

  return
end
subroutine ssd_cg ( lda, n, ndiag, offset, a, b, x )
!
!*******************************************************************************
!
!! SSD_CG uses conjugate gradient on a symmetric diagonal storage matrix.
!
!
!  Discussion:
!
!    The matrix A must be a positive definite symmetric matrix.
!    Only the nonzero diagonals on or above the main diagonal should be stored.
!
!    The method is designed to reach the solution to the linear system
!      A * x = b
!    after N computational steps.  However, roundoff may introduce
!    unacceptably large errors for some problems.  In such a case,
!    calling the routine a second time, using the current solution estimate
!    as the new starting guess, should result in improved results.
!
!  Reference:
!
!    F S Beckman,
!    The Solution of Linear Equations by the Conjugate Gradient Method,
!    Mathematical Methods for Digital Computers, pages 62-72.
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the matrix in symmetric diagonal storage.
!    Each column of A represents a diagonal of the original matrix.
!    The first entry of the diagonal is stored in the first row of the array.
!    The original location of the diagonal is specified by the OFFSET array.
!    For instance, if column 3 of A stores the main diagonal, OFFSET(3)=0.
!    If column I holds the first superdiagonal, then OFFSET(I) = 1,
!    or if it holds the fifth superdiagonal, then OFFSET(I) = 5.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.  Note that repeated
!    calls to this routine, using the value of X output on the previous
!    call, MAY improve the solution.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real alpha
  real ap(n)
  real b(n)
  real beta
  integer i
  integer ierror
  integer it
  integer offset(ndiag)
  real p(n)
  real pap
  real pr
  real r(n)
  real rap
  real x(n)
!
!  Check the dimensions.
!
  call ssd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSD_CG - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call ssd_mxv ( lda, n, ndiag, offset, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP = A*P.
!
    call ssd_mxv ( lda, n, ndiag, offset, a, p, ap )
!
!  Compute the dot products
!    PAP = P*AP,
!    PR  = P*R
!  Set
!    ALPHA = PR / PAP.
!
    pap = dot_product ( p(1:n), ap(1:n) )
    pr = dot_product ( p(1:n), r(1:n) )

    if ( pap == 0.0E+00 ) then
      return
    end if

    alpha = pr / pap
!
!  Set
!    X = X + ALPHA * P
!    R = R - ALPHA * AP.
!
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * ap(1:n)
!
!  Compute the vector dot product
!    RAP = R*AP
!  Set
!    BETA = - RAP / PAP.
!
    rap = dot_product ( r(1:n), ap(1:n) )

    beta = - rap / pap
!
!  Update the perturbation vector
!    P = R + BETA * P.
!
    p(1:n) = r(1:n) + beta * p(1:n)

  end do

  return
end
subroutine ssd_check ( lda, n, ndiag, ierror )
!
!*******************************************************************************
!
!! SSD_CHECK checks the dimensions of a symmetric diagonal storage matrix.
!
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1, and no more than N.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if N is illegal;
!    IERROR = IERROR + 4 if NDIAG is illegal.
!
  integer ierror
  integer lda
  integer n
  integer ndiag
!
  ierror = 0

  if ( lda < n ) then
    ierror = ierror + 1
    write ( *, * ) ' '
    write ( *, * ) 'SSD_CHECK - Illegal LDA = ', lda
  end if

  if ( n < 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SSD_CHECK - Illegal N = ', n
  end if

  if ( ndiag < 1 .or. ndiag > n ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SSD_CHECK - Illegal NDIAG = ', ndiag
  end if

  return
end
subroutine ssd_mxv ( lda, n, ndiag, offset, a, x, b )
!
!*******************************************************************************
!
!! SSD_MXV computes A * x where A is a symmetric diagonal storage matrix.
!
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the matrix in symmetric diagonal storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real b(n)
  integer i
  integer ierror
  integer j
  integer jdiag
  integer offset(ndiag)
  real x(n)
!
!  Check the dimensions.
!
  call ssd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSD_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  b(1:n) = 0.0E+00

  do i = 1, n
    do jdiag = 1, ndiag
      j = i + offset(jdiag)
      if ( j >= 1 .and. j <= n ) then
        b(i) = b(i) + a(i,jdiag) * x(j)
        if ( offset(jdiag) /= 0 ) then
          b(j) = b(j) + a(i,jdiag) * x(i)
        end if
      end if
    end do
  end do

  return
end
subroutine ssd_print ( lda, n, ndiag, offset, a, title )
!
!*******************************************************************************
!
!! SSD_PRINT prints a symmetric diagonal matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the N by N symmetric diagonal matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  integer offset(ndiag)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call ssd_print_some ( lda, n, ndiag, offset, a, 1, 1, n, n )

  return
end
subroutine ssd_print_some ( lda, n, ndiag, offset, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SSD_PRINT_SOME prints some of a symmetric diagonal matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the N by N symmetric diagonal matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jdiag
  integer jhi
  integer jlo
  integer off
  integer offset(ndiag)
  logical r_is_int
!
!  Check the dimensions.
!
  call ssd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSD_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0E+00
        off = j - i
        do jdiag = 1, ndiag
          if ( off == offset(jdiag) ) then
            aij = a(i,jdiag)
          else if ( off == - offset(jdiag) ) then
            aij = a(j,jdiag)
          end if
        end do

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine ssd_random ( lda, n, ndiag, offset, a )
!
!*******************************************************************************
!
!! SSD_RANDOM randomizes a symmetric diagonal storage matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Output, real A(LDA,NDIAG), the N by N matrix, stored by diagonals.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  integer i
  integer ierror
  integer j
  integer jj
  integer offset(ndiag)
!
!  Check the dimensions.
!
  call ssd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSD_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        call r_random ( 0.0E+00, 1.0E+00, a(i,j) )
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine ssd_to_sge ( lda, lda2, n, ndiag, offset, a, a2 )
!
!*******************************************************************************
!
!! SSD_TO_SGE copies a symmetric diagonal storage matrix to a general matrix.
!
!
!  Modified:
!
!    30 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer LDA2, the leading dimension of the array A2.
!    LDA2 must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal storage.
!
!    Input, real A(LDA,NDIAG), the N by N matrix, stored by diagonals.
!
!    Output, real A2(LDA,N), a copy of the input matrix, as a general matrix.
!
  integer lda
  integer lda2
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  real a2(lda2,n)
  integer i
  integer ierror
  integer j
  integer jj
  integer offset(ndiag)
!
!  Check the dimensions.
!
  call ssd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSD_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A.'
    return
  end if

  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSD_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A2.'
    return
  end if

  a2(1:n,1:n) = 0.0E+00

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        a2(i,jj) = a(i,j)
        if ( i /= jj ) then
          a2(jj,i) = a(i,j)
        end if
      end if
    end do
  end do

  return
end
subroutine ssd_zero ( lda, n, ndiag, a )
!
!*******************************************************************************
!
!! SSD_ZERO zeroes out a symmetric diagonal storage matrix.
!
!
!  Modified:
!
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real A(LDA,NDIAG), the N by N matrix, stored by diagonals.
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
  integer lda
  integer n
  integer ndiag
!
  real a(lda,ndiag)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call ssd_check ( lda, n, ndiag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSD_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  a(1:lda,1:ndiag) = 0.0E+00

  return
end
subroutine ssm_ml ( lda, n, a, u, v, pivot, x, b, job )
!
!*******************************************************************************
!
!! SSM_ML multiplies a factored Sherman Morrison matrix times a vector.
!
!
!  Modified:
!
!    04 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix factors computed by SGE_FA.
!
!    Input, real U(N), V(N), the Sherman Morrison vectors.
!
!    Input, integer PIVOT(N), the pivot vector computed by SGE_FA.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute (A-u*v') * x.
!    JOB nonzero, compute (A-u*v')' * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  integer pivot(n)
  integer job
  real u(n)
  real v(n)
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_ML - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  call sge_ml ( lda, n, a, pivot, x, b, job )

  if ( job == 0 ) then

    b(1:n) = b(1:n) - u(1:n) * dot_product ( v(1:n), x(1:n) )

  else

    b(1:n) = b(1:n) - v(1:n) * dot_product ( u(1:n), x(1:n) )

  end if

  return
end
subroutine ssm_mxv ( lda, n, a, u, v, x, b )
!
!*******************************************************************************
!
!! SSM_MXV multiplies a Sherman-Morrison matrix times a vector.
!
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix A.
!
!    Input, real U(N), V(N), the random vectors U and V that
!    define the Sherman-Morrison matrix (A-u*v').  
!
!    Input, real X(N), the vector to be multiplied by (A-u*v').
!
!    Output, real B(N), the product (A-u*v') * x.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  real u(n)
  real v(n)
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = matmul ( a(1:n,1:n), x(1:n) ) &
    - u(1:n) * dot_product ( v(1:n), x(1:n) )

  return
end
subroutine ssm_print ( lda, n, a, u, v, title )
!
!*******************************************************************************
!
!! SSM_PRINT prints a Sherman Morrison matrix.
!
!
!  Modified:
!
!    24 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the N by N matrix, stored in LINPACK
!    or LAPACK general storage mode.
!
!    Input, real U(N), V(N), the vectors that
!    define the Sherman-Morrison matrix (A-u*v').  
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  character ( len = * ) title
  real u(n)
  real v(n)
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call ssm_print_some ( lda, n, a, u, v, 1, 1, n, n )

  return
end
subroutine ssm_print_some ( lda, n, a, u, v, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SSM_PRINT_SOME prints some of a Sherman Morrison matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the N by N matrix, stored in LINPACK
!    or LAPACK general storage mode.
!
!    Input, real U(N), V(N), the vectors that
!    define the Sherman-Morrison matrix (A-u*v').  
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
!
  real a(lda,n)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
  real u(n)
  real v(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = a(i,j) - u(i) * v(j)

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine ssm_random ( lda, n, a, u, v )
!
!*******************************************************************************
!
!! SSM_RANDOM randomizes a Sherman-Morrison matrix.
!
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real A(LDA,N), the random matrix A.
!
!    Output, real U(N), V(N), the random vectors U and V that
!    define the perturbed matrix (A-u*v').  
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ierror
  integer j
  real u(n)
  real v(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_RANDOM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  call rmat_random ( 0.0E+00, 1.0E+00, lda, n, n, a )

  call rvec_random ( 0.0E+00, 1.0E+00, n, u )
  call rvec_random ( 0.0E+00, 1.0E+00, n, v )

  return
end
subroutine ssm_sl ( lda, n, a, u, v, b, ierror, pivot, job )
!
!*******************************************************************************
!
!! SSM_SL solves a linear system involving a Sherman Morrison matrix.
!
!
!  Discussion:
!
!    The linear system to be solved has the form
!
!      (A-u*v') * x = b.
!
!    The matrix Auv is related to the matrix A by a rank one update:
!
!    It is assumed that A has been decomposed into its LU factors
!    by SGE_FA.  The Sherman Morrison formula allows
!    us to solve linear systems involving (A-u*v') by solving linear
!    systems involving A and adjusting the results.
!
!  Reference:
!
!    Kahaner, Moler, and Nash
!    Numerical Methods and Software,
!    Prentice Hall, 1989
!
!  Modified:
!
!    04 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix A, as factored by SGE_FA.
!
!    Input, real U(N), V(N), the vectors U and V that define the
!    perturbed matrix (A-u*v').
!
!    Input/output, real B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Output, integer IERROR, an error flag.
!    0, no error occurred.  The solution was successfully computed.
!    1, an error occurred.  1 - Transpose(v) * Inverse(A) * u = 0.
!    The solution was not computed.
!
!    Input, integer PIVOT(N), the pivot vector produced by SGE_FA.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve (A-u*v') * X = B.
!    nonzero, solve (A-u*v') * X = B.
!
  integer lda
  integer n
!
  real a(lda,n)
  real alpha
  real b(n)
  real beta
  integer i
  integer ierror
  integer pivot(n)
  integer job
  integer job_local
  real u(n)
  real v(n)
  real w(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_SL - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  ierror = 0

  if ( job == 0 ) then
!
!  Solve A' * w = v.
!
    w(1:n) = v(1:n)

    job_local = 1
    call sge_sl ( lda, n, a, pivot, w, job_local )
!
!  Set beta = w' * b.
!
    beta = dot_product ( w(1:n), b(1:n) )
!
!  Solve A * b = b.
!
    job_local = 0
    call sge_sl ( lda, n, a, pivot, b, job_local )
!
!  Solve A * w = u.
!
    w(1:n) = u(1:n)

    job_local = 0
    call sge_sl ( lda, n, a, pivot, w, job_local )
!
!  Set alpha = 1 / ( 1 - v' * w ).
!
    alpha = 1.0E+00 - dot_product ( v(1:n), w(1:n) )

  else
!
!  Solve A * w = u.
!
    w(1:n) = u(1:n)

    job_local = 0
    call sge_sl ( lda, n, a, pivot, w, job_local )
!
!  Set beta = w' * b.
!
    beta = dot_product ( w(1:n), b(1:n) )
!
!  Solve A' * b = b.
!
    job_local = 1
    call sge_sl ( lda, n, a, pivot, b, job_local )
!
!  Solve A' * w = v.
!
    w(1:n) = v(1:n)

    job_local = 1
    call sge_sl ( lda, n, a, pivot, w, job_local )
!
!  Set alpha = 1 / ( 1 - u' * w ).
!
    alpha = 1.0E+00 - dot_product ( u(1:n), w(1:n) )

  end if

  if ( alpha == 0.0E+00 ) then
    ierror = 1
    write ( *, * ) ' '
    write ( *, * ) 'SSM_SL - Fatal error!'
    write ( *, * ) '  The divisor ALPHA is zero.'
    return
  end if

  alpha = 1.0E+00 / alpha
!
!  Set b = b + alpha * beta * w.
!
  b(1:n) = b(1:n) + alpha * beta * w(1:n)

  return
end
subroutine ssm_to_sge ( lda, lda2, n, a, u, v, a2 )
!
!*******************************************************************************
!
!! SSM_TO_SGE copies a Sherman-Morrison matrix into a general storage matrix.
!
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer LDA2, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix A.
!
!    Input, real U(N), V(N), the Sherman Morrison vectors U and V.
!
!    Output, real A2(LDA2,N), the Sherman Morrison matrix in general storage.
!
  integer lda
  integer lda2
  integer n
!
  real a(lda,n)
  real a2(lda2,n)
  integer i
  integer ierror
  integer j
  real u(n)
  real v(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A!'
    return
  end if

  call sge_check ( lda2, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A2!'
    return
  end if

  do i = 1, n
    a2(i,1:n) = a(i,1:n) - u(i) * v(1:n)
  end do

  return
end
subroutine ssm_vxm ( lda, n, a, u, v, x, b )
!
!*******************************************************************************
!
!! SSM_VXM multiplies a vector times a Sherman-Morrison matrix.
!
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the matrix A.
!
!    Input, real U(N), V(N), the random vectors U and V that
!    define the Sherman-Morrison matrix (A-u*v').  
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real B(N), the product (A-u*v')' * X.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ierror
  real u(n)
  real v(n)
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSM_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = matmul ( transpose ( a(1:n,1:n) ), x(1:n) ) &
    - v(1:n) * dot_product ( u(1:n), x(1:n) )

  return
end
subroutine sss_check ( diag, n, na, ierror )
!
!*******************************************************************************
!
!! SSS_CHECK checks dimensions for symmetric skyline matrix.
!
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIAG(N), the indices in A of the N diagonal elements.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NA, the dimension of the array A.
!    NA must be at least N.
!
!    Output, integer IERROR, error indicator.
!    0, no error.
!    1, N is less than 1.
!    2, NA is less than N.
!    3, DIAG(1) is not 1.
!    4, the elements of DIAG are not strictly increasing.
!    5, DIAG(N) is greater than NA.
!
  integer n
!
  integer diag(n)
  integer i
  integer ierror
  integer na
!
  ierror = 0

  if ( n < 1 ) then
    ierror = 1
    return
  end if

  if ( na < n ) then
    ierror = 2
    return
  end if

  if ( diag(1) /= 1 ) then
    ierror = 3
    return
  end if

  do i = 1, n-1
    if ( diag(i) >= diag(i+1) ) then
      ierror = 4
      return
    end if
  end do

  if ( diag(n) > na ) then
    ierror = 5
    return
  end if

  return
end
subroutine sss_mxv ( diag, n, na, a, x, b )
!
!*******************************************************************************
!
!! SSS_MXV multiplies a symmetric skyline matrix times a vector.
!
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIAG(N), the indices in A of the N diagonal elements.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, real A(NA), the zeroed matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product vector A*x.
!
  integer n
  integer na
!
  real a(na)
  real b(n)
  integer diag(n)
  integer diagold
  integer i
  integer ierror
  integer ilo
  integer j
  integer k
  real x(n)
!
!  Check the dimensions.
!
  call sss_check ( diag, n, na, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSS_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  b(1:n) = 0.0E+00

  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    do i = ilo, j-1
      k = k + 1
      b(i) = b(i) + a(k) * x(j)
      b(j) = b(j) + a(k) * x(i)
    end do

    k = k + 1
    b(j) = b(j) + a(k) * x(j)

    diagold = diag(j)

  end do

  return
end
subroutine sss_print ( n, na, a, diag, title )
!
!*******************************************************************************
!
!! SSS_PRINT prints a symmetric skyline matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NA, the dimension of the array A.
!
!    Input, real A(NA), the N by N matrix, stored in symmetric skyline
!    storage mode.
!
!    Input, integer DIAG(N), the indices in A of the N diagonal elements.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer na
  integer n
!
  real a(na)
  integer diag(n)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sss_print_some ( n, na, a, diag, 1, 1, n, n )

  return
end
subroutine sss_print_some ( n, na, a, diag, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SSS_PRINT_SOME prints some of a symmetric skyline matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIAG(N), the indices in A of the N diagonal elements.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NA, the dimension of the array A.
!
!    Input, real A(NA), the N by N matrix, stored in symmetric skyline
!    storage mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer na
  integer n
!
  real a(na)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer diag(n)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ij
  integer ijm1
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Check the dimensions.
!
  call sss_check ( diag, n, na, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSS_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0E+00

        if ( j < i ) then
          if ( i == 1 ) then
            ijm1 = 0
          else
            ijm1 = diag(i-1)
          end if
          ij = diag(i)
          if ( ij+j-i > ijm1 ) then
            aij = a(ij+j-i)
          end if
        else if ( j == i ) then
          ij = diag(j)
          aij = a(ij)
        else if ( j > i ) then
          if ( j == 1 ) then
            ijm1 = 0
          else
            ijm1 = diag(j-1)
          end if
          ij = diag(j)
          if ( ij+i-j > ijm1 ) then
            aij = a(ij+i-j)
          end if
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sss_random ( n, na, a, diag )
!
!*******************************************************************************
!
!! SSS_RANDOM randomizes a symmetric skyline matrix.
!
!
!  Note:
!
!    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
!    although the actual storage needed will generally be about half of
!    that.
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, integer NA, the dimension of the array A.
!    NA will be at least N and no greater than ( N * ( N + 1 ) ) / 2.
!
!    Output, real A((N*(N+1))/2), the randomized matrix, stored in entries
!    1 through NA.
!
!    Output, integer DIAG(N), the indices in A of the N diagonal elements.
!
  integer n
  integer na
!
  real a((n*(n+1))/2)
  integer diag(n)
  integer diagold
  integer i
  integer ilo
  integer j
  integer k
!
!  Set the values of DIAG.
!
  diag(1) = 1
  na = 1
  do i = 2, n
    call i_random ( 1, i, k )
    diag(i) = diag(i-1) + k
    na = na + k
  end do
!
!  Now set the values of A.
!
  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    do i = ilo, j
      k = k + 1
      call r_random ( 0.0E+00, 1.0E+00, a(k) )
    end do

    diagold = diag(j)

  end do

  return
end
subroutine sss_to_sge ( lda, n, na, a, diag, a2  )
!
!*******************************************************************************
!
!! SSS_TO_SGE copies a symmetric skyline matrix into a general matrix.
!
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NA, the dimension of the array A.
!    NA must be at least N.
!
!    Input, real A(NA), the symmetric skyline matrix.
!
!    Input, integer DIAG(N), the indices in A of the N diagonal elements.
!
!    Output, real A2(LDA,N), a copy of the matrix in general storage.
!
  integer lda
  integer n
  integer na
!
  real a(na)
  real a2(lda,n)
  integer diag(n)
  integer diagold
  integer i
  integer ierror
  integer ilo
  integer j
  integer k
!
!  Check the dimensions.
!
  call sss_check ( diag, n, na, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSS_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for skyline matrix!'
    return
  end if

  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSS_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix!'
    return
  end if

  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)

    a2(1:ilo-1,j) = 0.0E+00
    a2(j,1:ilo-1) = 0.0E+00

    do i = ilo, j-1
      k = k + 1
      a2(i,j) = a(k)
      a2(j,i) = a(k)
    end do

    k = k + 1
    a2(j,j) = a(k)

    diagold = diag(j)

  end do

  return
end
subroutine sss_zero ( n, na, a, diag )
!
!*******************************************************************************
!
!! SSS_ZERO zeroes out a symmetric skyline matrix.
!
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer NA, the dimension of the array A.
!    NA must be at least N.
!
!    Output, real A(NA), the zeroed matrix.
!
!    Input, integer DIAG(N), the indices in A of the N diagonal elements.
!
  integer n
  integer na
!
  real a(na)
  integer diag(n)
  integer diagold
  integer i
  integer ierror
  integer ihi
  integer ilo
  integer j
  integer k
!
!  Check the dimensions.
!
  call sss_check ( diag, n, na, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSS_ZERO - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  diagold = 0
  k = 0

  do j = 1, n

    ilo = j + 1 + diagold - diag(j)
    ihi = j

    do i = ilo, ihi
      k = k + 1
      a(k) = 0.0E+00
    end do

    diagold = diag(j)

  end do

  return
end
subroutine ssto_inv ( n, a, b )
!
!*******************************************************************************
!
!! SSTO_INV computes the inverse of a real symmetric Toeplitz matrix.
!
!
!  Discussion:
!
!    The matrix is also required to be positive definite.
!
!    This implementation of the algorithm assumes that the diagonal element
!    is 1.
!
!    The real symmetric Toeplitz matrix can be described by N numbers, which,
!    for convenience, we will label A(0:N-1).  The inverse matrix is not
!    guaranteed to be a Toeplitz matrix.  It is guaranteed to be symmetric
!    and persymmetric.
!
!  Example:
!
!    To compute the inverse of
!
!     1.0 0.5 0.2
!     0.5 1.0 0.5
!     0.2 0.5 1.0 
!
!    we input:
!
!      N = 3
!      A(0:2) = (/ 1.0, 0.5, 0.2 /)
!
!    with output:
!
!      B(1:3,1:3) = (/ 75, -40,   5,
!                     -40,  96, -40,
!                       5, -40,  75 /) / 56
!
!  Reference:
!
!    Gene Golub and Charles Van Loan,
!    Section 4.7.3, "Computing the Inverse",
!    Matrix Computations, Third Edition,
!    Johns Hopkins, 1996.
!
!  Modified:
!
!    25 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, real A(0:N-1), defines the symmetric Toeplitz matrix.
!
!    Output, real B(N,N), the inverse of the matrix.
!
  integer n
!
  real a(0:n-1)
  real b(n,n)
  integer i
  integer j
  real v(n)
!
  call ssto_yw_sl ( n-1, a(1), v )

  v(n) = 1.0E+00 / ( 1.0E+00 + dot_product ( a(1:n-1), v(1:n-1) ) )
  v(1:n-1) = v(n) * v(n-1:1:-1)

  b(1,1:n) = v(n:1:-1)
  b(n,1:n) = v(1:n)
  b(2:n-1,1) = v(n-1:2:-1)
  b(2:n,n) = v(2:n-1)

  do i = 2, 1+((n-1)/2)
    do j = i, n-i+1
      b(i,j) = b(i-1,j-1) + ( v(n+1-j) * v(n+1-i) - v(i-1) * v(j-1) ) / v(n)
      b(j,i) = b(i,j)
      b(n+1-i,n+1-j) = b(i,j)
      b(n+1-j,n+1-i) = b(i,j)
    end do
  end do

  return
end
subroutine ssto_sl ( n, a, b, x )
!
!*******************************************************************************
!
!! SSTO_SL solves a linear system with a real symmetric Toeplitz matrix.
!
!
!  Discussion:
!
!    The matrix is also required to be positive definite.
!
!    This implementation of the algorithm assumes that the diagonal element
!    is 1.
!
!    The real symmetric Toeplitz matrix can be described by N numbers, which,
!    for convenience, we will label A(0:N-1).
!
!    Note that there is a typographical error in the presentation
!    of this algorithm in the reference, and another in the presentation
!    of a sample problem.  Both involve sign errors.  A minor error
!    makes the algorithm incorrect for the case N = 1.
!
!  Example:
!
!    To solve
!
!     1.0 0.5 0.2    x1    4.0
!     0.5 1.0 0.5 *  x2 = -1.0
!     0.2 0.5 1.0    x3    3.0
!
!    we input:
!
!      N = 3
!      A(0:N-1) = (/ 1.0, 0.5, 0.2 /)
!      B(1:3) = (/ 4.0, -1.0, 3.0 /)
!
!    with output:
!
!      X(1:3) = (/ 355, -376, 285 /) / 56
!             = (/ 6.339, -6.714, 5.089 /)
!
!  Reference:
!
!    Gene Golub and Charles Van Loan,
!    Section 4.7.3, "The General Right Hand Side Problem",
!    Matrix Computations, Third Edition,
!    Johns Hopkins, 1996.
!
!  Modified:
!
!    29 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, real A(0:N-1), the first row of the matrix.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution of the linear system.
!
  integer n
!
  real a(0:n-1)
  real b(n)
  real beta
  integer k
  real x(n)
  real y(n)
!
  k = 0

  beta = 1.0E+00
  x(k+1) = b(k+1) / beta

  if ( k < n-1 ) then
    y(k+1) = -a(k+1) / beta
  end if

  do k = 1, n-1

    beta = ( 1.0E+00 - y(k)**2 ) * beta

    x(k+1) = ( b(k+1) - dot_product ( a(1:k), x(k:1:-1) ) ) / beta

    x(1:k) = x(1:k) + x(k+1) * y(k:1:-1)

    if ( k < n - 1 ) then
      y(k+1) = ( -a(k+1) - dot_product ( a(1:k), y(k:1:-1) ) ) / beta
      y(1:k) = y(1:k) + y(k+1) * y(k:1:-1)
    end if

  end do

  return
end
subroutine ssto_yw_sl ( n, b, x )
!
!*******************************************************************************
!
!! SSTO_YW_SL solves the Yule-Walker equations for a real symmetric Toeplitz matrix.
!
!
!  Discussion:
!
!    The matrix is also required to be positive definite.
!
!    This implementation of the algorithm assumes that the diagonal element
!    is 1.
!
!    The real symmetric Toeplitz matrix can be described by N numbers, which,
!    for convenience, we will label B(0:N-1).  We assume there is one more
!    number, B(N).  If we let A be the symmetric Toeplitz matrix whose first
!    row is B(0:N-1), then the Yule-Walker equations are:
!
!      A * X = -B(1:N)
!
!  Example:
!
!    To solve
!
!     1.0 0.5 0.2    x1   0.5
!     0.5 1.0 0.5 *  x2 = 0.2
!     0.2 0.5 1.0    x3   0.1
!
!    we input:
!
!      N = 3
!      B(1:3) = (/ 0.5, 0.2, 0.1 /)
!
!    with output:
!
!      X(1:3) = (/ -75, 12, -5 /) / 140
!
!  Reference:
!
!    Gene Golub and Charles Van Loan,
!    Section 4.7.2, "Solving the Yule-Walker Equations",
!    Matrix Computations, Third Edition,
!    Johns Hopkins, 1996.
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, real B(N), defines the linear system.  The first row of A is
!    a 1, followed by B(1) through B(N-1).  The right hand side of the 
!    system is -B(1:N).
!
!    Output, real X(N), the solution of the linear system.
!
  integer n
!
  real alpha
  real b(n)
  real beta
  integer i
  real x(n)
!
  x(1) = - b(1)
  beta = 1.0E+00
  alpha = - b(1)

  do i = 1, n-1
    beta = ( 1.0E+00 - alpha**2 ) * beta
    alpha = - ( b(i+1) + dot_product ( b(i:1:-1), x(1:i) ) ) / beta
    x(1:i) = x(1:i) + alpha * x(i:1:-1)
    x(i+1) = alpha
  end do

  return
end
subroutine ssto_mxv ( n, a, x, b )
!
!*******************************************************************************
!
!! SSTO_MXV multiplies a symmetric Toeplitz matrix times a vector.
!
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the entries of the first row of the matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  real x(n)
!
  do i = 1, n
    b(i) = dot_product ( a(i:2:-1), x(1:i-1) ) &
         + dot_product ( a(1:n+1-i), x(i:n) )
  end do

  return
end
subroutine ssto_print ( n, a, title )
!
!*******************************************************************************
!
!! SSTO_PRINT prints a symmetric Toeplitz matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(N), the entries of the first row of the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  real a(n)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call ssto_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine ssto_print_some ( n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! STO_PRINT_SOME prints some of a Toeplitz matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(N), the entries of the first row of the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a(n)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j >= i ) then
          aij = a(1+j-i)
        else
          aij = a(1+i-j)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine ssto_random ( n, a )
!
!*******************************************************************************
!
!! SSTO_RANDOM randomizes a symmetric Toeplitz matrix.
!
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real A(N), the randomized matrix, with entries between 0 and 1.
!
  integer n
!
  real a(n)
  integer i
!
  call rvec_random ( 0.0E+00, 1.0E+00, n, a )

  return
end
subroutine ssto_to_sge ( lda, n, a, a2 )
!
!*******************************************************************************
!
!! SSTO_TO_SGE copies a symmetric Toeplitz matrix into a general matrix.
!
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the first row of the matrix.
!
!    Output, real A2(LDA,N), the matrix stored as a general matrix.
!
  integer lda
  integer n
!
  real a(n)
  real a2(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SSTO_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    a2(i,1:i-1) = a(i:2:-1)
    a2(i,i:n) = a(1:n-i+1)
  end do

  return
end
subroutine sto_mxv ( n, a, x, b )
!
!*******************************************************************************
!
!! STO_MXV multiplies a Toeplitz matrix times a vector.
!
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer n
!
  real a(2*n-1)
  real b(n)
  integer i
  real x(n)
!
  do i = 1, n
    b(i) = dot_product ( a(n+i-1:n+1:-1), x(1:i-1) ) &
         + dot_product ( a(1:n+1-i), x(i:n) )
  end do

  return
end
subroutine sto_print ( n, a, title )
!
!*******************************************************************************
!
!! STO_PRINT prints a Toeplitz matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  real a(2*n-1)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sto_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine sto_print_some ( n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! STO_PRINT_SOME prints some of a Toeplitz matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a(2*n-1)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j >= i ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sto_random ( n, a )
!
!*******************************************************************************
!
!! STO_RANDOM randomizes a Toeplitz matrix.
!
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real A(2*N-1), the randomized matrix, with entries between
!    0 and 1.
!
  integer n
!
  real a(2*n-1)
  integer i
!
  call rvec_random ( 0.0E+00, 1.0E+00, 2*n-1, a )

  return
end
subroutine sto_sl ( n, a, b, x, job )
!
!***********************************************************************
!
!! STO_SL solves the real Toeplitz system A * X = B.
!
!
!  Modified:
!
!    11 March 2001
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(2*N-1), the first row of the Toeplitz matrix, followed by
!    the first column of the Toeplitz matrix beginning with the second element.
!
!    Input, real B(N) the right hand side vector.
!
!    Output, real X(N), the solution vector.  X and B may share the
!    same storage.
!
!    Input, integer JOB,
!    0 to solve A*X=B,
!    nonzero to solve Transpose(A)*X=B.
!
  integer n
!
  real a(2*n-1)
  real b(n)
  real c1(n-1)
  real c2(n-1)
  integer i
  integer job
  integer nsub
  real r1
  real r2
  real r3
  real r5
  real r6
  real x(n)
!
  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = a(n+nsub-1)
      r6 = a(nsub)
    else
      r5 = a(nsub)
      r6 = a(n+nsub-1)
    end if

    if ( nsub > 2 ) then

      c1(nsub-1) = r2

      do i = 1, nsub-2
        if ( job == 0 ) then
          r5 = r5 + a(n+i) * c1(nsub-i)
          r6 = r6 + a(i+1) * c2(i)
        else
          r5 = r5 + a(i+1) * c1(nsub-i)
          r6 = r6 + a(n+i) * c2(i)
        end if
      end do

    end if

    r2 = - r5 / r1
    r3 = - r6 / r1
    r1 = r1 + r5 * r3

    if ( nsub > 2 ) then

      r6 = c2(1)
      c2(nsub-1) = 0.0E+00

      do i = 2, nsub-1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = dot_product ( a(n+1:n+nsub-1), x(nsub-1:1:-1) )
    else
      r5 = dot_product ( a(2:nsub), x(nsub-1:1:-1) )
    end if

    r6 = ( b(nsub) - r5 ) / r1

    x(1:nsub-1) = x(1:nsub-1) + c2(1:nsub-1) * r6
    x(nsub) = r6

  end do

  return
end
subroutine sto_to_sge ( lda, n, a, a2 )
!
!*******************************************************************************
!
!! STO_TO_SGE copies a Toeplitz matrix into a general matrix.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(2*N-1), the Toeplitz matrix.
!
!    Output, real A2(LDA,N), the matrix stored as a general matrix.
!
  integer lda
  integer n
!
  real a(2*n-1)
  real a2(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'STO_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    a2(i,1:i-1) = a(n+i-1:n+1:-1)
    a2(i,i:n) = a(1:n-i+1)
  end do

  return
end
subroutine sto_vxm ( n, a, x, b )
!
!*******************************************************************************
!
!! STO_VXM multiplies a vector times a Toeplitz matrix.
!
!
!  Modified:
!
!    06 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(2*N-1), the entries of the first row of the Toeplitz
!    matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product Transpose ( A ) * X.
!
  integer n
!
  real a(2*n-1)
  real b(n)
  integer i
  real x(n)
!
  do i = 1, n

    b(i) = dot_product ( a(i:1:-1), x(1:i) ) + &
           dot_product ( a(n+1:2*n-i), x(i+1:n) )

  end do

  return
end
subroutine sut_det ( lda, n, a, det )
!
!*******************************************************************************
!
!! SUT_DET computes the determinant of an upper triangular matrix.
!
!
!  Modified:
!
!    22 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the upper triangular matrix.
!
!    Output, real DET, the determinant of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  real det
  real diag(n)
  integer ierror
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SUT_DET - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if

  call rmat_diag_get_vector ( lda, n, a, diag )

  det = product ( diag(1:n) )

  return
end
subroutine sut_inv ( lda, n, a )
!
!*******************************************************************************
!
!! SUT_INV computes the inverse of an upper triangular matrix.
!
!
!  Reference:
!
!    Combinatorial Algorithms,
!    A Nijenhuis and H Wilf,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(LDA,N).
!    On input, the upper triangular matrix to be inverted.
!    On output, the inverse of the upper triangular matrix.
!
  integer n
  integer lda
!
  real a(lda,n)
  integer i
  integer j
  integer k
!
!  Check.
!
  do i = 1, n
    if ( a(i,i) == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SUT_INV - Fatal error!'
      write ( *, * ) '  Zero diagonal element.'
      stop
    end if
  end do

  do j = n, 1, -1

    do i = n, 1, -1

      if ( i > j ) then

        a(i,j) = 0.0E+00

      else if ( i == j ) then

        a(i,j) = 1.0E+00 / a(i,j)

      else if ( i < j ) then

        a(i,j) = - dot_product ( a(i,i+1:j), a(i+1:j,j) ) / a(i,i)

      end if

    end do
  end do

  return
end
subroutine sut_mxv ( lda, m, n, a, x, b )
!
!*******************************************************************************
!
!! SUT_MXV computes A * x, where A is an upper triangular matrix.
!
!
!  Modified:
!
!    05 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N upper triangular matrix, stored
!    in LINPACK general matrix storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(M), the product A * x.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(m)
  integer i
  integer ierror
  integer j
  double precision temp
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SUT_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, m
    b(i) = dot_product ( a(i,i:n), x(i:n) )
  end do

  return
end
subroutine sut_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! SUT_PRINT prints an upper triangular matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N upper triangular matrix, stored in
!    LINPACK or LAPACK general band storage mode.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer m
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call sut_print_some ( lda, m, n, a, 1, 1, m, n )

  return
end
subroutine sut_print_some ( lda, m, n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SUT_PRINT_SOME prints some of an upper triangular matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N upper triangular matrix, stored in
!    LINPACK or LAPACK general band storage mode.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer lda
  integer n
!
  real a(lda,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ierror
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer m
  logical r_is_int
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SUT_PRINT - Fatal error!'
    write ( *, * ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j < i ) then
          ctemp(j2) = '              '
        else if ( r_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine sut_sl ( lda, n, a, b )
!
!*******************************************************************************
!
!! SUT_SL solves an upper triangular system.
!
!
!  Discussion:
!
!    No factorization of the upper triangular matrix is required.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(LDA,N), the upper triangular matrix.
!
!    Input/output, real B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer j
!
  do j = n, 1, -1
    b(j) = b(j) / a(j,j)
    b(1:j-1) = b(1:j-1) - a(1:j-1,j) * b(j)
  end do

  return
end
subroutine sut_vxm ( lda, m, n, a, x, b )
!
!*******************************************************************************
!
!! SUT_VXM computes Transpose ( A ) * x, where A is an upper triangular matrix.
!
!
!  Modified:
!
!    05 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real A(LDA,N), the M by N upper triangular matrix, stored
!    in LINPACK general matrix storage.
!
!    Input, real X(M), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer jhi
  double precision temp
  real x(m)
!
!  Check the dimensions.
!
  call sge_check ( lda, m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SUT_VXM - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, n
    jhi = min ( i, m )
    b(i) = dot_product ( x(1:jhi), a(1:jhi,i) )
  end do

  return
end
subroutine svm_det ( n, a, det )
!
!*******************************************************************************
!
!! SVM_DET computes the determinant of a Vandermonde matrix.
!
!
!  Modified:
!
!    20 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the values that define the Vandermonde matrix.
!
!    Output, real DET, the determinant of the matrix.
!
  integer n
!
  real a(n)
  real det
  integer i
  integer j
!
  det = 1.0E+00
  do j = 1, n
    do i = j+1, n
      det = det * ( a(i) - a(j) )
    end do
  end do

  return
end
subroutine svm_mxv ( n, a, x, b )
!
!*******************************************************************************
!
!! SVM_MXV multiplies a Vandermonde matrix times a vector.
!
!
!  Modified:
!
!    20 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the values that define the Vandermonde matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product A * x.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  integer j
  real x(n)
!
  do i = 1, n
    b(i) = 0.0E+00
    do j = 1, n
      if ( i == 1 ) then
        b(i) = b(i) + x(j)
      else
        b(i) = b(i) + a(j)**(i-1) * x(j)
      end if
    end do
  end do

  return
end
subroutine svm_print ( n, a, title )
!
!*******************************************************************************
!
!! SVM_PRINT prints a Vandermonde matrix.
!
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(N), the entries defining the Vandermonde matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  integer n
!
  real a(n)
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call svm_print_some ( n, a, 1, 1, n, n )

  return
end
subroutine svm_print_some ( n, a, ilo, jlo, ihi, jhi )
!
!*******************************************************************************
!
!! SVM_PRINT_SOME prints some of a Vandermonde matrix.
!
!
!  Discussion:
!
!    Only entries in rows ILO to IHI, columns JLO to JHI are considered.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real A(N), the entries defining the Vandermonde matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  integer, parameter :: incx = 5
!
  integer n
!
  real a(n)
  real aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  logical r_is_int
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns:'',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, * ) '  Row'
    write ( *, * ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i == 1 ) then
          aij = 1.0E+00
        else
          aij = a(j)**(i-1)
        end if

        if ( r_is_int ( aij ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) aij
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine svm_random ( n, a )
!
!*******************************************************************************
!
!! SVM_RANDOM randomizes a Vandermonde matrix.
!
!
!  Modified:
!
!    20 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real A(N), the randomized matrix, with entries between
!    0 and 1.
!
  integer n
!
  real a(n)
  integer i
!
  call rvec_random ( 0.0E+00, 1.0E+00, n, a )

  return
end
subroutine svm_sl ( n, a, b, x, job, info )
!
!*******************************************************************************
!
!! SVM_SL solves the system A * x = b with the Vandermonde matrix A.
!
!
!  Warning:
!
!    Vandermonde systems are very close to singularity.  The singularity
!    gets worse as N increases, and as any pair of values defining
!    the matrix get close.  Even a system as small as N = 10 will
!    involve the 9-th power of the defining values.
!
!  Modified:
!
!    21 November 1998
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the values defining the Vandermonde matrix.
!
!    Input, real B(N), the right hand side.
!
!    Output, real X(N), the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve transpose ( A ) * x = b.
!
!    Output, integer INFO.
!    0, no error.
!    nonzero, at least two of the values in A are equal.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  integer info
  integer j
  integer job
  real x(n)
!
  info = 0
  do j = 1, n - 1
    do i = j+1, n
      if ( a(i) == a(j) ) then
        info = 1
        return
      end if
    end do
  end do

  x(1:n) = b(1:n)

  if ( job == 0 ) then

    do j = 1, n-1
      do i = n, j+1, -1
        x(i) = x(i) - a(j) * x(i-1)
      end do
    end do

    do j = n-1, 1, -1

      do i = j+1, n
        x(i) = x(i) / ( a(i) - a(i-j) )
      end do

      do i = j, n-1
        x(i) = x(i) - x(i+1)
      end do

    end do

  else

    do j = 1, n-1
      do i = n, j+1, -1
        x(i) = ( x(i) - x(i-1) ) / ( a(i) - a(i-j) )
      end do
    end do

    do j = n-1, 1, -1
      do i = j, n-1
        x(i) = x(i) - x(i+1) * a(j)
      end do
    end do

  end if

  return
end
subroutine svm_to_sge ( lda, n, a, a2 )
!
!*******************************************************************************
!
!! SVM_TO_SGE copies a Vandermonde matrix into a general matrix.
!
!
!  Modified:
!
!    21 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A2.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), values defining the Vandermonde matrix.
!
!    Output, real A2(LDA,N), the Vandermonde matrix, stored as
!    a general matrix.
!
  integer lda
  integer n
!
  real a(n)
  real a2(lda,n)
  integer i
  integer ierror
  integer j
!
!  Check the dimensions.
!
  call sge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SVM_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for general matrix.'
    return
  end if

  do i = 1, n
    do j = 1, n
      if ( i == 1 ) then
        a2(i,j) = 1.0E+00
      else
        a2(i,j) = a(j)**(i-1)
      end if
    end do
  end do

  return
end
subroutine svm_vxm ( n, a, x, b )
!
!*******************************************************************************
!
!! SVM_VXM multiplies a vector times a Vandermonde matrix.
!
!
!  Modified:
!
!    20 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N), the values defining the Vandermonde matrix.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(N), the product Transpose ( A ) * x.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  integer j
  real x(n)
!
  do i = 1, n
    b(i) = 0.0E+00
    do j = 1, n
      if ( j == 1 ) then
        b(i) = b(i) + x(j)
      else
        b(i) = b(i) + a(i)**(j-1) * x(j)
      end if
    end do
  end do

  return
end
