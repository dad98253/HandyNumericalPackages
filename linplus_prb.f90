program linplus_prb
!
!*******************************************************************************
!
!! LINPLUS_PRB calls the LINPLUS test routines.
!
  character ( len = 8 ) date
  character ( len = 10 ) time
!
  call date_and_time ( date, time )

  write ( *, * ) ' '
  write ( *, * ) 'LINPLUS_PRB'
  write ( *, * ) '  Problems for LINPLUS.'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
 
  call test01
  call test014
  call test015
  call test016
  call test017
  call test02
  call test03
  call test035
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
  call test155
  call test16
  call test17
  call test18
  call test19
  call test195
  call test20

  call test21
  call test22
  call test23
  call test24
  call test25
  call test26
  call test27
  call test28
  call test29
  call test295
 
  call test30
  call test31
  call test315
  call test32
  call test33
  call test34
  call test35
  call test36
  call test37
  call test38
  call test385
  call test39

  call test40
  call test41
  call test42
  call test43
  call test44
  call test45
  call test46
  call test47
  call test48
  call test49

  call test50
  call test51
  call test52
  call test53
  call test54
  call test55
  call test56
  call test57
  call test58
  call test583
  call test585
  call test587
  call test59

  call test60
  call test61
  call test62
  call test63

  write ( *, * ) ' '
  write ( *, * ) 'LINPLUS_PRB'
  write ( *, * ) '  Normal end of LINPLUS tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests C3_CR_FA;
!! TEST01 tests C3_CR_SL.
!   
  integer, parameter :: n = 10
!
  complex diag(2*n)
  integer i
  complex rhs(0:2*n)
  complex subd(0:2*n)
  complex supd(0:2*n)
  complex x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  C3_CR_FA factors a complex tridiagonal matrix;'
  write ( *, * ) '  C3_CR_SL solves a factored system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) ' '
!
!  Set the matrix values.
!
  do i = 1, n
    diag(i) = cmplx ( 2, 2*i )
  end do

  do i = 1, n-1
    subd(i) = cmplx ( -1, -i )
  end do

  do i = 1, n-1
    supd(i) = cmplx ( -1, -i-1 )
  end do
!
!  Set the desired solution.
!
  do i = 1, n
    x(i) = cmplx ( i, 10 * i )
  end do
!
!  Compute the corresponding right hand side.
!
  i = 1
  rhs(1) = cmplx ( 2, 2*i ) * x(i) + cmplx ( -1, -i-1 ) * x(i+1)

  do i = 2, n-1
    rhs(i) = cmplx (-1, -i+1 ) * x(i-1) + cmplx ( 2, 2*i ) * x(i) &
           + cmplx ( -1, -i-1 ) * x(i+1)
  end do

  rhs(n) = cmplx ( -1, -i+1 ) * x(i-1) + cmplx ( 2, 2*i ) * x(i)
!
!  Factor the matrix.
!
  call c3_cr_fa ( n, subd, diag, supd )
!
!  Solve the linear system.
!
  call c3_cr_sl ( n, subd, diag, supd, rhs )

  write ( *, * ) ' '
  write ( *, * ) '  Solution:'
  write ( *, * ) ' '

  call cvec_print_some ( n, rhs, 10 )

  return
end
subroutine test014
!
!*******************************************************************************
!
!! TEST014 tests C3_NP_FA;
!! TEST014 tests C3_NP_SL.
!! TEST014 tests C3_NP_ML.
!
  integer, parameter :: n = 10
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  complex b(n)
  integer i
  integer info
  integer job
  complex x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST014'
  write ( *, * ) '  For a complex tridiagonal matrix that can be'
  write ( *, * ) '    factored with no pivoting,'
  write ( *, * ) '  C3_NP_FA factors;'
  write ( *, * ) '  C3_NP_SL solves a factored system.'
  write ( *, * ) '  C3_NP_ML multiplies A*X when A has been factored.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c3_random ( n, a1, a2, a3 )

  call c3_print ( n, a1, a2, a3, '  The tridiagonal matrix' )
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( i, 10 * i )
  end do
!
!  Compute the corresponding right hand side.
!
  call c3_mxv ( n, a1, a2, a3, x, b )

  call cvec_print ( n, b, '  The right hand side' )
!
!  Factor the matrix.
!
  call c3_np_fa ( n, a1, a2, a3, info )
!
!  Solve the linear system.
!
  job = 0
  call c3_np_sl ( n, a1, a2, a3, b, job )

  call cvec_print ( n, b, '  The solution' )
!
!  Now set a SECOND desired solution.
!
  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do
!
!  Compute the corresponding right hand side, using the FACTORED matrix.
!
  call c3_np_ml ( n, a1, a2, a3, x, b, job )

  call cvec_print ( n, b, '  The second right hand side' )
!
!  Solve the linear system.
!
  call c3_np_sl ( n, a1, a2, a3, b, job )

  call cvec_print ( n, b, '  The second solution' )

  return
end
subroutine test015
!
!*******************************************************************************
!
!! TEST015 tests C3_NP_FA;
!! TEST015 tests C3_NP_ML.
!! TEST015 tests C3_NP_SL.
!
  integer, parameter :: n = 10
!
  complex a1(2:n)
  complex a2(1:n)
  complex a3(1:n-1)
  complex b(n)
  integer i
  integer info
  integer job
  complex x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST015'
  write ( *, * ) '  For a complex tridiagonal matrix that can be'
  write ( *, * ) '    factored with no pivoting,'
  write ( *, * ) '  C3_NP_FA factors;'
  write ( *, * ) '  C3_NP_SL solves a factored system.'
  write ( *, * ) '  C3_NP_ML multiplies A*X when A has been factored.'
  write ( *, * ) ' '
  write ( *, * ) '  We will look at the TRANSPOSED linear system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c3_random ( n, a1, a2, a3 )

  call c3_print ( n, a1, a2, a3, '  The tridiagonal matrix' )
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( i, 10 * i )
  end do
!
!  Compute the corresponding right hand side.
!
  call c3_vxm ( n, a1, a2, a3, x, b )

  call cvec_print ( n, b, '  The right hand side B1' )
!
!  Factor the matrix.
!
  call c3_np_fa ( n, a1, a2, a3, info )
!
!  Solve the linear system.
!
  job = 1
  call c3_np_sl ( n, a1, a2, a3, b, job )
 
  call cvec_print ( n, b, '  The solution to At * X1 = B1' )
!
!  Set the second solution.
!
  do i = 1, n
    x(i) = cmplx ( 10 * i, i )
  end do
!
!  Compute the corresponding right hand side.
!
  call c3_np_ml ( n, a1, a2, a3, x, b, job )

  call cvec_print ( n, b, '  The second right hand side B2' )
!
!  Solve the linear system.
!
  call c3_np_sl ( n, a1, a2, a3, b, job )
 
  call cvec_print ( n, b, '  Solution to At * X2 = B2' )
 
  return
end
subroutine test016
!
!*******************************************************************************
!
!! TEST016 tests CCI_SL.
!
  integer, parameter :: n = 10
!
  complex a(n)
  complex b(n)
  integer i
  integer job
  complex x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST016'
  write ( *, * ) '  CCI_SL solves a complex circulant system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call cci_random ( n, a )

  call cci_print ( n, a, '  The circulant matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    do i = 1, n
      x(i) = cmplx ( i, 10 * i )
    end do
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call cci_mxv ( n, a, x, b )
    else
      call cci_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call cci_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call cvec_print ( n, x, '  Solution:' )
    else
      call cvec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test017
!
!*******************************************************************************
!
!! TEST017 tests CTO_SL.
!
  integer, parameter :: n = 4
!
  complex a(2*n-1)
  complex b(n)
  integer i
  integer job
  complex x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST017'
  write ( *, * ) '  CTO_SL solves a complex Toeplitz system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call cto_random ( n, a )

  call cto_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call cvec_identity ( n, x )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  Desired solution:'
    else
      write ( *, * ) '  Desired solution to transposed system:'
    end if

    write ( *, * ) ' '
    call cvec_print_some ( n, x, 10 )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call cto_mxv ( n, a, x, b )
    else
      call cto_vxm ( n, a, x, b )
    end if

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  Right Hand Side:'
    else
      write ( *, * ) '  Right Hand Side of transposed system:'
    end if

    write ( *, * ) ' '
    call cvec_print_some ( n, b, 10 )
!
!  Solve the linear system.
!
    call cto_sl ( n, a, b, x, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  Solution:'
    else
      write ( *, * ) '  Solution to transposed system:'
    end if

    write ( *, * ) ' '
    call cvec_print_some ( n, x, 10 )

  end do

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests S3_CR_FA;
!! TEST02 tests S3_CR_SL.
!   
  integer, parameter :: n = 100
!
  real diag(2*n)
  integer i
  integer j
  real rhs(0:2*n)
  real subd(0:2*n)
  real supd(0:2*n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  S3_CR_FA factors a real tridiagonal matrix;'
  write ( *, * ) '  S3_CR_SL solves a factored system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Demonstrate multiple system solution method.'
  write ( *, * ) ' '
!
!  Set the matrix values.
!
  diag(1:n) = 2.0E+00
  subd(1:n-1) = -1.0E+00
  supd(1:n-1) = -1.0E+00
!
!  Factor the matrix once.
!
  call s3_cr_fa ( n, subd, diag, supd )

  do j = 1, 2

    write ( *, * ) ' '
    write ( *, * ) '  Solve linear system number ', J
    rhs(1:n) = 0.0E+00
    if ( j == 1 ) then
      rhs(n) = n + 1
    else
      rhs(1) = 1.0E+00
      rhs(n) = 1.0E+00
    end if
!
!  Solve the linear system.
!
    call s3_cr_sl ( n, subd, diag, supd, rhs )

    write ( *, * ) ' '
    write ( *, * ) '  Solution:'
    write ( *, * ) ' '
    call rvec_print_some ( n, rhs, 10 )

  end do

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests S3_CR_FA;
!! TEST03 tests S3_CR_SL.
!
  integer, parameter :: n = 100
!
  real diag(2*n)
  integer i
  real rhs(0:2*n)
  real subd(0:2*n)
  real supd(0:2*n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  For a real tridiagonal matrix,'
  write ( *, * ) '  using CYCLIC REDUCTION,'
  write ( *, * ) '  S3_CR_FA factors;'
  write ( *, * ) '  S3_CR_SL solves.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ',N
  write ( *, * ) '  The matrix is NOT symmetric.'
  write ( *, * ) ' '
!
!  Set the matrix values.
!
  do i = 1, n
    diag(i) = 4.0E+00 * i
  end do

  subd(0) = 0.0E+00
  do i = 1, n-1
    subd(i) = i
  end do

  do i = 1, n-1
    supd(i) = i
  end do
  supd(n) = 0.0E+00
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  rhs(1:n) = 0.0E+00

  do i = 2, n
    rhs(i) = rhs(i) + subd(i-1) * x(i-1)
  end do

  do i = 1, n
    rhs(i) = rhs(i) + diag(i) * x(i)
  end do

  do i = 1, n-1
    rhs(i) = rhs(i) + supd(i) * x(i+1)
  end do
!
!  Factor the matrix.
!
  call s3_cr_fa ( n, subd, diag, supd )
!
!  Solve the linear system.
!
  call s3_cr_sl ( n, subd, diag, supd, rhs )

  write ( *, * ) ' '
  write ( *, * ) '  The solution:'
  write ( *, * ) ' '
  call rvec_print_some ( n, rhs, 10 )

  return
end
subroutine test035
!
!*******************************************************************************
!
!! TEST035 tests S3_GS_SL.
!   
  integer, parameter :: n = 100
!
  real a1(2:n)
  real a2(n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer job
  integer maxit
  real x(n)
!
  maxit = 1000

  write ( *, * ) ' '
  write ( *, * ) 'TEST035'
  write ( *, * ) '  For a real tridiagonal system,'
  write ( *, * ) '  S3_GS_SL solves a linear system using'
  write ( *, * ) '    Gauss-Seidel iteration'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Iterations per call = ', maxit
  write ( *, * ) ' '
!
!  Set the matrix values.
!
  a1(2:n) = -1.0E+00
  a2(1:n) = 2.0E+00
  a3(1:n-1) = -1.0E+00
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call s3_mxv ( n, a1, a2, a3, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0E+00
!
!  Solve the linear system.
!
  job = 0

  do i = 1, 3

    call s3_gs_sl ( n, a1, a2, a3, b, x, maxit, job )

    write ( *, * ) ' '
    write ( *, * ) '  Solution after call ', i
    write ( *, * ) ' '
    call rvec_print_some ( n, x, 10 )

  end do

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests S3_JAC_SL.
!   
  integer, parameter :: n = 100
!
  real a1(2:n)
  real a2(n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer job
  integer maxit
  real x(n)
!
  maxit = 1000

  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  For a real tridiagonal system,'
  write ( *, * ) '  S3_JAC_SL solves a linear system using'
  write ( *, * ) '    Jacobi iteration'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Iterations per call = ', maxit
  write ( *, * ) ' '
!
!  Set the matrix values.
!
  a1(2:n) = -1.0E+00
  a2(1:n) = 2.0E+00
  a3(1:n-1) = -1.0E+00
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call s3_mxv ( n, a1, a2, a3, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0E+00
!
!  Solve the linear system.
!
  job = 0

  do i = 1, 3

    call s3_jac_sl ( n, a1, a2, a3, b, x, maxit, job )

    write ( *, * ) ' '
    write ( *, * ) '  Solution after call ', i
    write ( *, * ) ' '
    call rvec_print_some ( n, x, 10 )

  end do

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests S3_NP_DET;
!! TEST05 tests S3_NP_FA.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real det
  integer info
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  For a tridiagonal matrix that can be factored'
  write ( *, * ) '    with no pivoting,'
  write ( *, * ) '  S3_NP_FA factors,'
  write ( *, * ) '  S3_NP_DET computes the determinant.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call s3_random ( n, a1, a2, a3 )
!
!  Copy the matrix into general storage.
!
  call s3_to_sge ( lda, n, a1, a2, a3, a )
!
!  Factor the matrix.
!
  call s3_np_fa ( n, a1, a2, a3, info )
!
!  Compute the determinant.
!
  call s3_np_det ( n, a2, det )

  write ( *, * ) ' '
  write ( *, * ) '  S3_NP_DET computes determinant = ', det
!
!  Factor the general matrix.
!
  call sge_np_fa ( lda, n, a, info )
!
!  Compute the determinant.
!
  call sge_np_det ( lda, n, a, det )

  write ( *, * ) '  SGE_DET computes determinant = ', det

  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 tests S3_NP_FA;
!! TEST06 tests S3_NP_SL.
!
  integer, parameter :: n = 10
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  integer info
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  For a tridiagonal matrix that can be factored'
  write ( *, * ) '    with no pivoting,'
  write ( *, * ) '  S3_NP_FA factors;'
  write ( *, * ) '  S3_NP_SL solves a factored system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call s3_random ( n, a1, a2, a3 )

  call s3_print ( n, a1, a2, a3, '  The tridiagonal matrix:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call s3_mxv ( n, a1, a2, a3, x, b )
!
!  Factor the matrix.
!
  call s3_np_fa ( n, a1, a2, a3, info )
!
!  Solve the linear system.
!
  job = 0
  call s3_np_sl ( n, a1, a2, a3, b, job )
 
  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call s3_np_ml ( n, a1, a2, a3, x, b, job )
!
!  Solve the linear system.
!
  job = 1
  call s3_np_sl ( n, a1, a2, a3, b, job )
 
  call rvec_print ( n, b, '  Solution to tranposed system:' )
 
  return
end
subroutine test07
!
!*******************************************************************************
!
!! TEST07 tests S3_NP_FS.
!
  integer, parameter :: n = 10
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST07'
  write ( *, * ) '  S3_NP_FS factors and solves a tridiagonal'
  write ( *, * ) '    linear system.'
!
!  Set the matrix elements.
!
!
!  Compute b = A * x.
!
  call s3_random ( n, a1, a2, a3 )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute b = A * x.
!
  call s3_mxv ( n, a1, a2, a3, x, b )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0E+00
!
!  Solve the system.
!
  call s3_np_fs ( n, a1, a2, a3, b, x )

  call rvec_print ( n, x, '  Solution:' )

  return
end
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests S3_NP_ML.
!
  integer, parameter :: n = 10
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  real b2(n)
  integer info
  integer i
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST08'
  write ( *, * ) '  S3_NP_ML computes A*x or transpose(A)*X'
  write ( *, * ) '    where A has been factored by S3_FA.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
  do job = 0, 1
!
!  Set the matrix.
!
    call s3_random ( n, a1, a2, a3 )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call s3_mxv ( n, a1, a2, a3, x, b )
    else
      call s3_vxm ( n, a1, a2, a3, x, b )
    end if
!
!  Factor the matrix.
!
    call s3_np_fa ( n, a1, a2, a3, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  S3_NP_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call s3_np_ml ( n, a1, a2, a3, x, b2, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  A*x and PLU*x'
    else
      write ( *, * ) '  AT*x and (PLU)T*x'
    end if

    write ( *, * ) ' '
    call rvec2_print_some ( n, b, b2, 10 )

  end do

  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests S3P_DET.
!
  integer, parameter :: n = 12
  integer, parameter :: lda = n
!
  real a(lda,n)
  real a1(n)
  real a2(n)
  real a3(n)
  real det
  integer info
  integer ipivot(n)
  real work2(n-1)
  real work3(n-1)
  real work4
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST09'
  write ( *, * ) '  S3P_DET, determinant of a tridiagonal'
  write ( *, * ) '    periodic matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call s3p_random ( n, a1, a2, a3 )

  call s3p_print ( n, a1, a2, a3, '  The periodic tridiagonal matrix:' )
!
!  Copy the matrix into a general array.
!
  call s3p_to_sge ( lda, n, a1, a2, a3, a )
!
!  Factor the matrix.
!
  call s3p_fa ( n, a1, a2, a3, info, work2, work3, work4 )
!
!  Compute the determinant.
!
  call s3p_det ( n, a2, work4, det )

  write ( *, * ) ' '
  write ( *, * ) '  S3P_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call sge_fa ( lda, n, a, ipivot, info )
!
!  Compute the determinant.
!
  call sge_det ( lda, n, a, ipivot, det )

  write ( *, * ) '  SGE_DET computes the determinant = ', det

  return
end
subroutine test10
!
!*******************************************************************************
!
!! TEST10 tests S3P_FA;
!! TEST10 tests S3P_SL.
!
  integer, parameter :: n = 10
!
  real a1(n)
  real a2(n)
  real a3(n)
  real b(n)
  integer i
  integer info
  integer job
  real work2(n-1)
  real work3(n-1)
  real work4
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST10'
  write ( *, * ) '  S3P_FA factors a tridiagonal periodic system.'
  write ( *, * ) '  S3P_SL solves a tridiagonal periodic system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call s3p_random ( n, a1, a2, a3 )
!
!  Factor the matrix.
!
  call s3p_fa ( n, a1, a2, a3, info, work2, work3, work4 )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  S3P_FA returns INFO = ', info
    return
  end if

  do job = 0, 1
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    call s3p_ml ( n, a1, a2, a3, x, b, job )
!
!  Solve the linear system.
!
    call s3p_sl ( n, a1, a2, a3, b, x, job, work2, work3, work4 )

    if ( job == 0 ) then
      call rvec_print ( n, x, '  Solution:' )
    else
      call rvec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine test11
!
!*******************************************************************************
!
!! TEST11 tests S3P_ML.
!
  integer, parameter :: n = 10
!
  real a1(n)
  real a2(n)
  real a3(n)
  real b(n)
  real b2(n)
  integer info
  integer i
  integer job
  real work2(n-1)
  real work3(n-1)
  real work4
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST11'
  write ( *, * ) '  S3P_ML computes A*x or transpose(A)*X'
  write ( *, * ) '    where A has been factored by S3P_FA.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call s3p_random ( n, a1, a2, a3 )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call s3p_mxv ( n, a1, a2, a3, x, b )
    else
      call s3p_vxm ( n, a1, a2, a3, x, b )
    end if
!
!  Factor the matrix.
!
    call s3p_fa ( n, a1, a2, a3, info, work2, work3, work4 )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  S3P_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call s3p_ml ( n, a1, a2, a3, x, b2, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  A*x and PLU*x'
    else
      write ( *, * ) '  AT*x and (PLU)T*x'
    end if

    write ( *, * ) ' '
    call rvec2_print_some ( n, b, b2, 10 )

  end do

  return
end
subroutine test12
!
!*******************************************************************************
!
!! TEST12 tests S5_FS.
!
  integer, parameter :: n = 10
!
  real a1(3:n)
  real a2(2:n)
  real a3(1:n)
  real a4(1:n-1)
  real a5(1:n-2)
  real b(n)
  integer i
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST12'
  write ( *, * ) '  S5_FS factors and solves a pentadiagonal'
  write ( *, * ) '    linear system.'
  write ( *, * ) ' '
!
!  Set the matrix elements.
!
!
!  Set the matrix to a random value.
!
  call s5_random ( n, a1, a2, a3, a4, a5 )

  call s5_print ( n, a1, a2, a3, a4, a5, '  The pentadiagonal matrix:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute b = A * x.
!
  call s5_mxv ( n, a1, a2, a3, a4, a5, x, b )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0E+00
!
!  Solve the system.
!
  call s5_fs ( n, a1, a2, a3, a4, a5, b, x )

  call rvec_print ( n, x, '  Solution:' )

  return
end
subroutine test13
!
!*******************************************************************************
!
!! TEST13 tests SBB_FA.
!! TEST13 tests SBB_PRINT.
!! TEST13 tests SBB_RANDOM.
!! TEST13 tests SBB_SL.
!
  integer, parameter :: n1 = 8
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2
!
  real a(na)
  real b(n)
  integer i
  integer info
  integer ipivot(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST13'
  write ( *, * ) '  For a border banded matrix:'
  write ( *, * ) '  SBB_FA factors;'
  write ( *, * ) '  SBB_PRINT prints;'
  write ( *, * ) '  SBB_RANDOM randomizes;'
  write ( *, * ) '  SBB_SL solves.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Matrix suborders N1 = ', n1, ' N2 = ', n2
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sbb_random ( n1, n2, ml, mu, a )

  call sbb_print ( n1, n2, ml, mu, a, '  The border-banded matrix:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n1+n2, x )
!
!  Compute the corresponding right hand side.
!
  call sbb_mxv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call sbb_fa ( n1, n2, ml, mu, a, ipivot, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SBB_FA claims the matrix is singular.'
    write ( *, * ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call sbb_sl ( n1, n2, ml, mu, a, ipivot, b )

  call rvec_print ( n, b, '  Solution:' )

  return
end
subroutine test14
!
!*******************************************************************************
!
!! TEST14 tests SBB_FA.
!! TEST14 tests SBB_SL.
!
  integer, parameter :: n1 = 0
  integer, parameter :: n2 = 10
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 0
  integer, parameter :: mu = 0
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2
!
  real a(na)
  real b(n)
  integer i
  integer info
  integer ipivot(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST14'
  write ( *, * ) '  For a border banded matrix:'
  write ( *, * ) '  SBB_FA factors;'
  write ( *, * ) '  SBB_SL solves.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Matrix suborders N1 = ', n1, ' N2 = ', n2
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sbb_random ( n1, n2, ml, mu, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sbb_mxv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call sbb_fa ( n1, n2, ml, mu, a, ipivot, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SBB_FA claims the matrix is singular.'
    write ( *, * ) '  The value of INFO is ',info
    write ( *, * ) ' '
    return
  end if
!
!  Solve the system.
!
  call sbb_sl ( n1, n2, ml, mu, a, ipivot, b )

  call rvec_print ( n, b, '  Solution:' )
 
  return
end
subroutine test15
!
!*******************************************************************************
!
!! TEST15 tests SBB_FA.
!! TEST15 tests SBB_SL.

  integer, parameter :: n1 = 10
  integer, parameter :: n2 = 0
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2
!
  real a(na)
  real b(n)
  integer i
  integer info
  integer ipivot(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST15'
  write ( *, * ) '  For a border banded matrix:'
  write ( *, * ) '  SBB_FA factors;'
  write ( *, * ) '  SBB_SL solves.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Matrix suborders N1 = ', n1, ' N2 = ', n2
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sbb_random ( n1, n2, ml, mu, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sbb_mxv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix
!
  call sbb_fa ( n1, n2, ml, mu, a, ipivot, info )
 
  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SBB_FA claims the matrix is singular.'
    write ( *, * ) '  The value of INFO is ',info
    write ( *, * ) ' '
    return
  end if
!
!  Solve the system.
!
  call sbb_sl ( n1, n2, ml, mu, a, ipivot, b )

  call rvec_print ( n, b, '  Solution:' )
 
  return
end
subroutine test155
!
!***********************************************************************
!
!! TEST155 tests SBTO_MXV.
!! TEST155 tests SBTO_VXM.
!
  integer, parameter :: l = 3
  integer, parameter :: m = 2
!
  real, dimension ( m, m, l ) ::  a1 = reshape ( (/ &
    1.0E+00, 5.0E+00, 2.0E+00, 5.0E+00, &
    3.0E+00, 6.0E+00, 4.0E+00, 6.0E+00, &
    5.0E+00, 7.0E+00, 6.0E+00, 7.0E+00 /), (/ m, m, l /) )

  real, dimension ( m, m, l-1 ) :: a2 = reshape ( (/ &
    7.0E+00, 8.0E+00, 8.0E+00, 8.0E+00, &
    9.0E+00, 9.0E+00, 0.0E+00, 9.0E+00 /), (/ m, m, l-1 /) )
  real b(m,l)
  real x(m,l)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST155'
  write ( *, * ) '  For a real block Toeplitz matrix,'
  write ( *, * ) '  SBTO_MXV computes A * x.'
  write ( *, * ) '  SBTO_VXM computes x * A.'

  call sbto_print ( m, l, a1, a2, '  The block Toeplitz matrix:' )

  call rvec_identity ( m*l, x )

  call rvec_print ( m*l, x, '  The vector x:' )

  call sbto_mxv ( m, l, a1, a2, x, b )

  call rvec_print ( m*l, b, '  The product A*x:' )

  call sbto_vxm ( m, l, a1, a2, x, b )

  call rvec_print ( m*l, b, '  The product x*A:' )

  return
end
subroutine test16
!
!*******************************************************************************
!
!! TEST16 tests SCB_NP_FA.
!! TEST16 tests SCB_DET.
!
  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3
  integer, parameter :: lda = ml + mu + 1
  integer, parameter :: lda2 = n
!
  real a(lda,n)
  real a2(lda2,n)
  real det
  integer info
  integer ipivot(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST16'
  write ( *, * ) '  For a compact band matrix, no pivoting:'
  write ( *, * ) '  SCB_NP_FA factors;'
  write ( *, * ) '  SCB_DET computes the determinant;'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call scb_random ( lda, n, ml, mu, a )

  call scb_print ( lda, n, ml, mu, a, '  The compact band matrix:' )
!
!  Copy the matrix into a general array.
!
  call scb_to_sge ( lda, lda2, ml, mu, n, a, a2 )
!
!  Factor the matrix.
!
  call scb_np_fa ( lda, n, ml, mu, a, info )
!
!  Compute the determinant.
!
  call scb_det ( lda, n, ml, mu, a, det )

  write ( *, * ) ' '
  write ( *, * ) '  SCB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call sge_fa ( lda2, n, a2, ipivot, info )
!
!  Compute the determinant.
!
  call sge_det ( lda2, n, a2, ipivot, det )

  write ( *, * ) '  SGE_DET computes the determinant = ', det

  return
end
subroutine test17
!
!*******************************************************************************
!
!! TEST17 tests SCB_NP_FA.
!! TEST17 tests SCB_NP_SL.
!
  integer, parameter :: n = 10
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2
  integer, parameter :: lda = ml + mu + 1
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST17'
  write ( *, * ) '  For a compact band matrix, no pivoting:'
  write ( *, * ) '  SCB_NP_FA factors;'
  write ( *, * ) '  SCB_NP_SL solves.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call scb_random ( lda, n, ml, mu, a )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the right hand side.
!
    if ( job == 0 ) then
      call scb_mxv ( lda, n, ml, mu, a, x, b )
    else
      call scb_vxm ( lda, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call scb_np_fa ( lda, n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SCB_NP_FA claims the matrix is singular.'
      write ( *, * ) '  The value of info is ',INFO
      write ( *, * ) ' '
      return
    end if
!
!  Solve the system.
!
    call scb_np_sl ( lda, n, ml, mu, a, b, job )

    if ( job == 0 ) then
      call rvec_print ( n, b, '  Solution:' )
    else
      call rvec_print ( n, b, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test18
!
!*******************************************************************************
!
!! TEST18 tests SCB_ML.
!
  integer, parameter :: n = 10
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2
  integer, parameter :: lda = 2 * ml + mu + 1
!
  real a(lda,n)
  real b(n)
  real b2(n)
  integer i
  integer info
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST18'
  write ( *, * ) '  For a compact band matrix:'
  write ( *, * ) '  SCB_ML computes A*x or transpose(A)*X'
  write ( *, * ) '    where A has been factored by SCB_FA.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call scb_random ( lda, n, ml, mu, a )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call scb_mxv ( lda, n, ml, mu, a, x, b )
    else
      call scb_vxm ( lda, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call scb_np_fa ( lda, n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SCB_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call scb_ml ( lda, n, ml, mu, a, x, b2, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  A*x and PLU*x'
    else
      write ( *, * ) '  AT*x and (PLU)T*x'
    end if

    write ( *, * ) ' '
    call rvec2_print_some ( n, b, b2, 10 )

  end do

  return
end
subroutine test19
!
!*******************************************************************************
!
!! TEST19 tests SCBB_FA.
!! TEST19 tests SCBB_PRINT.
!! TEST19 tests SCBB_RANDOM.
!! TEST19 tests SCBB_SL.
!
  integer, parameter :: n1 = 8
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2
!
  real a(na)
  real b(n)
  integer i
  integer info
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST19'
  write ( *, * ) '  For a compressed border banded matrix:'
  write ( *, * ) '  SCBB_RANDOM randomly generates;'
  write ( *, * ) '  SCBB_PRINT prints;'
  write ( *, * ) '  SCBB_FA factors (no pivoting);'
  write ( *, * ) '  SCBB_SL solves.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Matrix suborders N1 = ', n1, ' N2 = ', n2
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call scbb_random ( n1, n2, ml, mu, a )

  call scbb_print ( n1, n2, ml, mu, a, '  The compact border-banded matrix:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call scbb_mxv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix
!
  call scbb_fa ( n1, n2, ml, mu, a, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SCBB_FA claims the matrix is singular.'
    write ( *, * ) '  The value of INFO is ', info
    write ( *, * ) ' '
    return
  end if
!
!  Solve the system.
!
  call scbb_sl ( n1, n2, ml, mu, a, b )

  call rvec_print ( n, b, '  Solution:' )

  return
end
subroutine test195
!
!*******************************************************************************
!
!! TEST195 tests SCI_EVAL.
!
  integer, parameter :: n = 5
!
  real a(n)
  complex lambda(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST195'
  write ( *, * ) '  SCI_EVAL finds the eigenvalues of a real circulant system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sci_random ( n, a )

  call sci_print ( n, a, '  The circulant matrix:' )

  call sci_eval ( n, a, lambda )

  call cvec_print ( n, lambda, '  The eigenvalues:' )

  return
end
subroutine test20
!
!*******************************************************************************
!
!! TEST20 tests SCI_SL.
!
  integer, parameter :: n = 10
!
  real a(n)
  real b(n)
  integer i
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST20'
  write ( *, * ) '  SCI_SL solves a circulant system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sci_random ( n, a )

  call sci_print ( n, a, '  The circulant matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call sci_mxv ( n, a, x, b )
    else
      call sci_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call sci_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call rvec_print ( n, x, '  Solution:' )
    else
      call rvec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine test21
!
!*******************************************************************************
!
!! TEST21 tests SGB_DET.
!
  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 3
  integer, parameter :: mu = 2
  integer, parameter :: lda = 2 * ml + mu + 1
  integer, parameter :: lda2 = n
!
  real a(lda,n)
  real a2(lda2,n)
  real det
  integer info
  integer ipivot(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST21'
  write ( *, * ) '  SGB_DET, determinant of a banded matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Number of rows M = ', m
  write ( *, * ) '  Number of columns N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sgb_random ( lda, m, n, ml, mu, a )
!
!  Copy the matrix into a general array.
!
  call sgb_to_sge ( lda, lda2, m, ml, mu, n, a, a2 )
!
!  Factor the matrix.
!
  call sgb_fa ( lda, n, ml, mu, a, ipivot, info )
!
!  Compute the determinant.
!
  call sgb_det ( lda, n, ml, mu, a, ipivot, det )

  write ( *, * ) ' '
  write ( *, * ) '  SGB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call sge_fa ( lda2, n, a2, ipivot, info )
!
!  Compute the determinant.
!
  call sge_det ( lda2, n, a2, ipivot, det )

  write ( *, * ) '  SGE_DET computes the determinant = ', det

  return
end
subroutine test22
!
!*******************************************************************************
!
!! TEST22 tests SGB_FA.
!! TEST22 tests SGB_SL.
!
  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 1 
  integer, parameter :: mu = 2 
  integer, parameter :: lda = 2 * ml + mu + 1
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  integer ipivot(n)
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST22'
  write ( *, * ) '  SGB_FA factor routine'
  write ( *, * ) '    for general band matrices.'
  write ( *, * ) '  SGB_SL solve routine'
  write ( *, * ) '    for general band matrices.'
  write ( *, * ) ' '
  write ( *, * ) '  Number of matrix rows M = ', m
  write ( *, * ) '  Number of matrix columns N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sgb_random ( lda, m, n, ml, mu, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sgb_mxv ( lda, m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call sgb_fa ( lda, n, ml, mu, a, ipivot, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SGB_FA declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sgb_sl ( lda, n, ml, mu, a, ipivot, b, job )
 
  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call sgb_ml ( lda, n, ml, mu, a, ipivot, x, b, job )
!
!  Solve the linear system.
!
  job = 1
  call sgb_sl ( lda, n, ml, mu, a, ipivot, b, job )
 
  call rvec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine test23
!
!*******************************************************************************
!
!! TEST23 tests SGB_FA.
!! TEST23 tests SGB_TRF.
!
  integer, parameter :: m = 5
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: lda = 2 * ml + mu + 1
!
  real a(lda,n)
  integer info
  integer ipivot(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST23'
  write ( *, * ) '  SGB_FA factors a general band matrix, using'
  write ( *, * ) '    LINPACK conventions;'
  write ( *, * ) '  SGB_TRF factors a general band matrix, using'
  write ( *, * ) '    LAPACK conventions;'
  write ( *, * ) ' '
  write ( *, * ) '  Number of matrix rows M = ', m
  write ( *, * ) '  Number of matrix columns N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sgb_random ( lda, m, n, ml, mu, a )
!
!  Factor the matrix.
!
  call sgb_fa ( lda, n, ml, mu, a, ipivot, info )

  call sgb_print ( lda, m, n, ml, mu, a, '  The SGB_FA factors:' )
!
!  Set the matrix.
!
  call sgb_random ( lda, m, n, ml, mu, a )
!
!  Factor the matrix.
!
  call sgb_trf ( lda, m, n, ml, mu, a, ipivot, info )

  call sgb_print ( lda, m, n, ml, mu, a, '  The SGB_TRF factors:')

  return
end
subroutine test24
!
!*******************************************************************************
!
!! TEST24 tests SGB_ML.
!
  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2
  integer, parameter :: lda = 2 * ml + mu + 1
!
  real a(lda,n)
  real b(n)
  real b2(n)
  integer i
  integer info
  integer ipivot(n)
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST24'
  write ( *, * ) '  SGB_ML computes A*x or transpose(A)*X'
  write ( *, * ) '    where A has been factored by SGB_FA.'
  write ( *, * ) ' '
  write ( *, * ) '  Number of matrix rows M = ', m
  write ( *, * ) '  Number of matrix columns N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call sgb_random ( lda, m, n, ml, mu, a )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call sgb_mxv ( lda, m, n, ml, mu, a, x, b )
    else
      call sgb_vxm ( lda, m, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call sgb_fa ( lda, n, ml, mu, a, ipivot, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SGB_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call sgb_ml ( lda, n, ml, mu, a, ipivot, x, b2, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  A*x and PLU*x'
    else
      write ( *, * ) '  AT*x and (PLU)T*x'
    end if

    write ( *, * ) ' '
    call rvec2_print_some ( n, b, b2, 10 )

  end do

  return
end
subroutine test25
!
!*******************************************************************************
!
!! TEST25 tests SGB_PRINT.
!
  integer, parameter :: m = 10
  integer, parameter :: n = 10
  integer, parameter :: ml = 3
  integer, parameter :: mu = 1
  integer, parameter :: lda = 2 * ml + mu + 1
!
  real a(lda,n)
  integer i
  integer ihi
  integer ilo
  integer j
  integer jhi
  integer jlo
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST25'
  write ( *, * ) '  SGB_PRINT prints out a banded matrix.'

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  do i = 1, m
    do j = 1, n
      if ( i - j <= ml .and. j - i <= mu ) then
        a(i-j+ml+mu+1,j) = real ( 10 * i + j )
      end if
    end do
  end do

  call sgb_print ( lda, m, n, ml, mu, a, '  The banded matrix:' )

  return
end
subroutine test26
!
!*******************************************************************************
!
!! TEST26 tests SGB_SCAN.
!
  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2
  integer, parameter :: lda = 2 * ml + mu + 1
!
  real a(lda,n)
  integer i
  integer j
  integer nonzer
  integer nzer
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST26'
  write ( *, * ) '  SGB_SCAN counts zero/nonzero entries'
  write ( *, * ) '    in a general band matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Number of matrix rows M = ', m
  write ( *, * ) '  Number of matrix columns N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sgb_random ( lda, m, n, ml, mu, a )
!
!  Make some zero entries.
!
  do i = 1, lda
    do j = 1, n
      if ( a(i,j) < 0.3E+00 ) then
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  call sgb_scan ( lda, m, n, ml, mu, a, nonzer, nzer )

  write ( *, * ) ' '
  write ( *, * ) '  Nonzero entries = ', nonzer
  write ( *, * ) '  Zero entries = ', nzer

  return
end
subroutine test27
!
!*******************************************************************************
!
!! TEST27 tests SGB_TRF.
!! TEST27 tests SGB_TRS.
!
  integer, parameter :: n = 10
  integer, parameter :: m = n
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2
  integer, parameter :: lda = 2 * ml + mu + 1
  integer, parameter :: ldb = n
  integer, parameter ::  nrhs = 1
!
  real a(lda,n)
  real b(ldb,nrhs)
  integer i
  integer info
  integer ipivot(n)
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST27'
  write ( *, * ) '  SGB_TRF factor routine'
  write ( *, * ) '    for general band matrices.'
  write ( *, * ) '  SGB_TRS solve routine'
  write ( *, * ) '    for general band matrices.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix rows M = ', m
  write ( *, * ) '  Matrix columns N = ', n
  write ( *, * ) '  Bandwidths are ML = ', ml, ' MU = ', mu
!
!  Set the matrix.
!
  call sgb_random ( lda, m, n, ml, mu, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sgb_mxv ( lda, m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call sgb_trf ( lda, m, n, ml, mu, a, ipivot, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SGB_TRF declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if
!
!  Solve the linear system.
!
  call sgb_trs ( lda, n, ml, mu, nrhs, 'N', a, ipivot, b, ldb, info )

  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call sgb_mu ( lda, n, ml, mu, a, ipivot, x, b, job )
!
!  Solve the linear system.
!
  call sgb_trs ( lda, n, ml, mu, nrhs, 'T', a, ipivot, b, ldb, info )

  call rvec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine test28
!
!*******************************************************************************
!
!! TEST28 tests SGD_MXV.
!! TEST28 tests SGD_PRINT.
!! TEST221 tests SGD_RANDOM.
!! TEST221 tests SGD_VXM.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
  integer, parameter :: ndiag = 4
!
  real a(lda,ndiag)
  real b(n)
  integer i
  integer j
  integer offset(ndiag)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST28'
  write ( *, * ) '  For a general diagonal matrix:'
  write ( *, * ) '  SGD_MXV computes A * x;'
  write ( *, * ) '  SGD_PRINT prints it;'
  write ( *, * ) '  SGD_RANDOM randomly generates one;'
  write ( *, * ) '  SGD_VXM computes Transpose(A)*x;'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1) = -2
  offset(2) = 0
  offset(3) = 1
  offset(4) = n - 1

  call sgd_random ( lda, n, ndiag, offset, a )

  call sge_print ( lda, n, ndiag, a, '  The raw matrix: ' )

  call sgd_print ( lda, n, ndiag, offset, a, '  The general diagonal matrix:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sgd_mxv ( lda, n, ndiag, offset, a, x, b )

  call rvec_print ( n, b, '  A * x:' )
!
!  Compute the corresponding right hand side.
!
  call sgd_vxm ( lda, n, ndiag, offset, a, x, b )

  call rvec_print ( n, b, '  Transpose ( A ) * x:' )

  return
end
subroutine test29
!
!*******************************************************************************
!
!! TEST29 tests SGE_DET.
!
  integer, parameter :: n = 4
  integer, parameter :: lda = n
!
  real a(lda,n)
  real det
  integer i
  integer info
  integer ipivot(n)
  integer j
  real true
  real x
  real y
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST29'
  write ( *, * ) '  SGE_DET, determinant of a general matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) ' '
!
!  Set the matrix.
!
  y = 3.0E+00
  x = 2.0E+00
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = x + y
      else
        a(i,j) = y
      end if
    end do
  end do
!
!  Factor the matrix.
!
  call sge_fa ( lda, n, a, ipivot, info )
!
!  Compute the determinant.
!
  call sge_det ( lda, n, a, ipivot, det )

  write ( *, * ) '  SGE_DET computes the determinant = ',det

  true = x**(n-1) * ( x + n * y )
  write ( *, * ) '  True determinant = ', true

  return
end
subroutine test295
!
!*******************************************************************************
!
!! TEST295 tests SGE_DILU.
!
  integer, parameter :: ncol = 3
  integer, parameter :: nrow = 3
  integer, parameter :: n = nrow * ncol
  integer, parameter :: m = n
  integer, parameter :: lda = n
!
  real a(lda,lda)
  real d(n)
  integer i
  integer j
  integer k
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST295'
  write ( *, * ) '  SGE_DILU returns the DILU factors of a matrix.'

  do i = 1, nrow * ncol
    do j = 1, nrow * ncol

      if ( i == j ) then
        a(i,j) = 4.0E+00
      else if ( i == j + 1 .or. i == j - 1 .or. &
                i == j + nrow .or. i == j - nrow ) then
        a(i,j) = -1.0E+00
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do

  call sge_print ( lda, m, n, a, '  Matrix A:' )
!
!  Compute the incomplete LU factorization.
!
  call sge_dilu ( lda, m, n, a, d )

  call rvec_print ( n, d, '  DILU factor:' )

  return
end
subroutine test30
!
!*******************************************************************************
!
!! TEST30 tests SGE_FA;
!! TEST30 tests SGE_SL.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  integer ipivot(n)
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST30'
  write ( *, * ) '  SGE_FA factors a general linear system,'
  write ( *, * ) '  SGE_SL solves a factored system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sge_random ( lda, n, n, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sge_mxv ( lda, n, n, a, x, b )
!
!  Factor the matrix.
!
  call sge_fa ( lda, n, a, ipivot, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SGE_FA declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sge_sl ( lda, n, a, ipivot, b, job )
 
  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  x(1:n) = 1.0E+00
!
!  Compute the corresponding right hand side.
!
  job = 0
  call sge_ml ( lda, n, a, ipivot, x, b, job )
!
!  Solve the system
!
  job = 0
  call sge_sl ( lda, n, a, ipivot, b, job )

  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call sge_ml ( lda, n, a, ipivot, x, b, job )
!
!  Solve the system
!
  job = 1
  call sge_sl ( lda, n, a, ipivot, b, job )

  call rvec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine test31
!
!*******************************************************************************
!
!! TEST31 tests SGE_FA;
!! TEST31 tests SGE_SL.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  integer ipivot(n)
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST31'
  write ( *, * ) '  SGE_FA factors a general linear system,'
  write ( *, * ) '  SGE_SL solves a factored system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sge_random ( lda, n, n, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sge_mxv ( lda, n, n, a, x, b )
!
!  Factor the matrix.
!
  call sge_fa ( lda, n, a, ipivot, info )
 
  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SGE_FA declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sge_sl ( lda, n, a, ipivot, b, job )

  call rvec_print ( n, b, '  Solution:' )

  return
end
subroutine test315
!
!*******************************************************************************
!
!! TEST315 tests SGE_ILU.
!
  integer, parameter :: ncol = 3
  integer, parameter :: nrow = 3
  integer, parameter :: n = nrow * ncol
  integer, parameter :: m = n
  integer, parameter :: lda = n
!
  real a(lda,lda)
  integer i
  integer j
  integer k
  real l(lda,lda)
  real lu(lda,lda)
  real u(lda,lda)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST315'
  write ( *, * ) '  SGE_ILU returns the ILU factors of a matrix.'

  do i = 1, nrow * ncol
    do j = 1, nrow * ncol

      if ( i == j ) then
        a(i,j) = 4.0E+00
      else if ( i == j + 1 .or. i == j - 1 .or. &
                i == j + nrow .or. i == j - nrow ) then
        a(i,j) = -1.0E+00
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do

  call sge_print ( lda, m, n, a, '  Matrix A:' )
!
!  Compute the incomplete LU factorization.
!
  call sge_ilu ( lda, m, n, a, l, u )

  call sge_print ( lda, m, m, l, '  Factor L:' )

  call sge_print ( lda, m, n, u, '  Factor U:' )

  lu(1:m,1:n) = matmul ( l(1:m,1:m), u(1:m,1:n) )

  call sge_print ( lda, m, n, lu, '  Product L*U:' )

  return
end
subroutine test32
!
!*******************************************************************************
!
!! TEST32 tests SGE_NP_FA;
!! TEST32 tests SGE_NP_SL.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST32'
  write ( *, * ) '  SGE_NP_FA factors without pivoting,'
  write ( *, * ) '  SGE_NP_SL solves factored systems.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sge_random ( lda, n, n, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0E+00
!
!  Compute the corresponding right hand side.
!
  call sge_mxv ( lda, n, n, a, x, b )
!
!  Factor the matrix.
!
  call sge_np_fa ( lda, n, a, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SGE_NP_FA declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call sge_np_sl ( lda, n, a, b, job )
 
  write ( *, * ) ' '
  write ( *, * ) '  Solution:'
  write ( *, * ) ' '
  call rvec_print_some ( n, b, 10 )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call sge_np_ml ( lda, n, a, x, b, job )
!
!  Solve the system
!
  job = 0
  call sge_np_sl ( lda, n, a, b, job )

  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call sge_np_ml ( lda, n, a, x, b, job )
!
!  Solve the system
!
  job = 1
  call sge_np_sl ( lda, n, a, b, job )

  call rvec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine test33
!
!*******************************************************************************
!
!! TEST33 tests SGE_NP_FA;
!! TEST33 tests SGE_NP_INV.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  integer i
  integer info
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST33'
  write ( *, * ) '  SGE_NP_FA factors without pivoting,'
  write ( *, * ) '  SGE_NP_INV computes the inverse.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sge_random ( lda, n, n, a )

  call sge_print ( lda, n, n, a, '  The random matrix:' )
!
!  Factor and invert the matrix.
!
  b(1:n,1:n) = a(1:n,1:n)

  call sge_np_fa ( lda, n, b, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SGE_NP_FA declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if

  call sge_np_inv ( lda, n, b )

  call sge_print ( lda, n, n, b, '  The inverse matrix:' )
!
!  Compute A * B = C.
!
  call sge_mxm ( lda, n, a, b, c )

  call sge_print ( lda, n, n, c, '  The product:' )

  return
end
subroutine test34
!
!*******************************************************************************
!
!! TEST34 tests SGE_FS.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST34'
  write ( *, * ) '  SGE_FS, full storage factor/solve routine.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sge_random ( lda, n, n, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sge_mxv ( lda, n, n, a, x, b )
!
!  Factor and solve the system.
!
  call sge_fs ( lda, n, a, b, info )
  
  call rvec_print ( n, b, '  Solution:' )

  return
end
subroutine test35
!
!*******************************************************************************
!
!! TEST35 tests SGE_INV.
!
  integer, parameter :: n = 4
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  integer info
  integer ipivot(n)
  real, parameter :: x = 2.0E+00
  real, parameter :: y = 3.0E+00
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST35'
  write ( *, * ) '  SGE_INV inverts a general matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  a(1:n,1:n) = y

  call rmat_diag_add_scalar ( lda, n, a, x )

  call sge_print ( lda, n, n, a, '  Matrix A:' )
!
!  Factor and invert the matrix.
!
  b(1:n,1:n) = a(1:n,1:n)

  call sge_fa ( lda, n, b, ipivot, info )

  call sge_inv ( lda, n, b, ipivot )

  call sge_print ( lda, n, n, b, '  Inverse matrix B:' )
!
!  Check.
!
  call sge_mxm ( lda, n, a, b, c )

  call sge_print ( lda, n, n, c, '  Product matrix:' )

  return
end
subroutine test36
!
!*******************************************************************************
!
!! TEST36 tests SGE_ML.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  real b2(n)
  integer info
  integer i
  integer ipivot(n)
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST36'
  write ( *, * ) '  SGE_ML computes A*x or transpose(A)*X'
  write ( *, * ) '    where A has been factored by SGE_FA.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call sge_random ( lda, n, n, a )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call sge_mxv ( lda, n, n, a, x, b )
    else
      call sge_vxm ( lda, n, n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call sge_fa ( lda, n, a, ipivot, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SGE_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call sge_ml ( lda, n, a, ipivot, x, b2, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  A*x and PLU*x'
    else
      write ( *, * ) '  AT*x and (PLU)T*x'
    end if

    write ( *, * ) ' '
    call rvec2_print_some ( n, b, b2, 10 )

  end do

  return
end
subroutine test37
!
!*******************************************************************************
!
!! TEST37 tests SGE_NP_ML.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  real b2(n)
  integer info
  integer i
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST37'
  write ( *, * ) '  SGE_NP_ML computes A*x or transpose(A)*X'
  write ( *, * ) '    where A has been factored by SGE_NP_FA.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call sge_random ( lda, n, n, a )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call sge_mxv ( lda, n, n, a, x, b )
    else
      call sge_vxm ( lda, n, n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call sge_np_fa ( lda, n, a, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SGE_NP_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call sge_np_ml ( lda, n, a, x, b2, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  A*x and PLU*x'
    else
      write ( *, * ) '  AT*x and (PLU)T*x'
    end if

    write ( *, * ) ' '
    call rvec2_print_some ( n, b, b2, 10 )

  end do

  return
end
subroutine test38
!
!*******************************************************************************
!
!! TEST38 tests SGE_MU.
!
  integer, parameter :: m = 5
  integer, parameter :: n = 3
  integer, parameter :: ldam = m
  integer, parameter :: ldan = n
!
  real amn(ldam,n)
  real anm(ldan,m)
  real bm(m)
  real bn(n)
  real cm(m)
  real cn(n)
  integer info
  integer i
  integer ipivot(m+n)
  integer job
  character trans
  real xm(m)
  real xn(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST38'
  write ( *, * ) '  SGE_MU computes A*x or transpose(A)*X'
  write ( *, * ) '    where A has been factored by SGE_TRF.'
  write ( *, * ) ' '
  write ( *, * ) '  Number of matrix rows M = ', m
  write ( *, * ) '  Number of matrix columns N = ', n

  do job = 0, 1

    if ( job == 0 ) then 
      trans = 'N'
    else
      trans = 'T'
    end if
!
!  Set the matrix.
!
    call sge_random ( ldam, m, n, amn )

    if ( job == 0 ) then

      call rvec_identity ( n, xn )

      call sge_mxv ( ldam, m, n, amn, xn, cm )

    else

      call rvec_identity ( m, xm )

      call sge_vxm ( ldam, m, n, amn, xm, cn )

    end if
!
!  Factor the matrix.
!
    call sge_trf ( ldam, m, n, amn, ipivot, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SGE_TRF declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    if ( job == 0 ) then

      call sge_mu ( ldam, m, n, amn, trans, ipivot, xn, bm )

      write ( *, * ) ' '
      write ( *, * ) '  A*x and PLU*x'
      write ( *, * ) ' '
      call rvec2_print_some ( m, cm, bm, 10 )

    else

      call sge_mu ( ldam, m, n, amn, trans, ipivot, xm, bn )

      write ( *, * ) ' '
      write ( *, * ) '  AT*x and (PLU)T*x'
      write ( *, * ) ' '
      call rvec2_print_some ( n, cn, bn, 10 )

    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) '  Matrix is ', n, ' by ', m

  do job = 0, 1

    if ( job == 0 ) then 
      trans = 'N'
    else
      trans = 'T'
    end if
!
!  Set the matrix.
!
    call sge_random ( ldan, n, m, anm )

    if ( job == 0 ) then

      call rvec_identity ( m, xm )

      call sge_mxv ( ldan, n, m, anm, xm, cn )

    else

      call rvec_identity ( n, xn )

      call sge_vxm ( ldan, n, m, anm, xn, cm )

    end if
!
!  Factor the matrix.
!
    call sge_trf ( ldan, n, m, anm, ipivot, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SGE_TRF declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    if ( job == 0 ) then

      call sge_mu ( ldan, n, m, anm, trans, ipivot, xm, bn )

      write ( *, * ) ' '
      write ( *, * ) '  A*x and PLU*x'
      write ( *, * ) ' '
      call rvec2_print_some ( n, cn, bn, 10 )

    else

      call sge_mu ( ldan, n, m, anm, trans, ipivot, xn, bm )

      write ( *, * ) ' '
      write ( *, * ) '  AT*x and (PLU)T*x'
      write ( *, * ) ' '
      call rvec2_print_some ( m, cm, bm, 10 )

    end if

  end do

  return
end
subroutine test385
!
!*******************************************************************************
!
!! TEST385 tests SGE_PLU.
!
  integer, parameter :: lda = 5
!
  real a(lda,lda)
  integer i
  integer j
  integer k
  real l(lda,lda)
  real lu(lda,lda)
  integer m
  integer n
  real p(lda,lda)
  real plu(lda,lda)
  real u(lda,lda)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST385'
  write ( *, * ) '  SGE_PLU returns the PLU factors of a matrix.'

  n = 4

  do m = 3, 5

    k = 0
    do i = 1, m
      do j = 1, n
        k = k + 1
        a(i,j) = k
      end do
    end do

    call sge_print ( lda, m, n, a, '  Matrix A:' )
!
!  Compute the PLU factors.
!
    call sge_plu ( lda, m, n, a, p, l, u )

    call sge_print ( lda, m, m, p, '  Factor P:' )

    call sge_print ( lda, m, m, l, '  Factor L:' )

    call sge_print ( lda, m, n, u, '  Factor U:' )

    do i = 1, m
      do j = 1, n
        lu(i,j) = 0.0E+00
        do k = 1, m
          lu(i,j) = lu(i,j) + l(i,k) * u(k,j)
        end do
      end do
    end do

    do i = 1, m
      do j = 1, n
        plu(i,j) = 0.0E+00
        do k = 1, m
          plu(i,j) = plu(i,j) + p(i,k) * lu(k,j)
        end do
      end do
    end do
        
    call sge_print ( lda, m, n, plu, '  Product P*L*U:')

  end do

  return
end
subroutine test39
!
!*******************************************************************************
!
!! TEST39 tests SGE_POLY.
!
  integer, parameter :: n = 12
  integer, parameter :: lda = n
!
  real a(lda,n)
  integer i
  integer j
  real p(0:n)
  real true(0:n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST39'
  write ( *, * ) '  SGE_POLY computes the characteristic polynomial'
  write ( *, * ) '    of a matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  true(0) =       1.0E+00
  true(1) =    - 23.0E+00
  true(2) =     231.0E+00
  true(3) =  - 1330.0E+00
  true(4) =    4845.0E+00
  true(5) = - 11628.0E+00
  true(6) =   18564.0E+00
  true(7) = - 19448.0E+00
  true(8) =   12870.0E+00
  true(9) =  - 5005.0E+00
  true(10) =   1001.0E+00
  true(11) =   - 78.0E+00
  true(12) =      1.0E+00
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = min ( i, j )
    end do
  end do
!
!  Get the characteristic polynomial.
!
  call sge_poly ( lda, n, a, p )
!
!  Compare.
!
  write ( *, * ) ' '
  write ( *, * ) 'I, P(I), True P(I)'
  write ( *, * ) ' '

  call rvec2_print_some ( n+1, p, true, 10 )

  return
end
subroutine test40
!
!*******************************************************************************
!
!! TEST40 tests SGE_SL_IT.
!
  integer, parameter :: n = 6
  integer, parameter :: lda = n
!
  real a(lda,n)
  real alu(lda,n)
  real b(n)
  integer i
  integer info
  integer ipivot(n)
  integer j
  integer job
  real r(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST40'
  write ( *, * ) '  SGE_SL_IT applies one step of iterative '
  write ( *, * ) '    refinement to an SGE_SL solution.'
!
!  Set the coefficient matrix.
!
  call hilb_inv ( lda, n, a )
!
!  Set the right hand side b.
!
  b(1:n-1) = 0.0E+00
  b(n) = 1.0E+00
!
!  It is necessary to keep both an unfactored and factored copy
!  of the coefficient matrix.
!
  alu = a
!
!  Compute the factored coefficient matrix.
!
  call sge_fa ( lda, n, alu, ipivot, info )
!
!  Solve the system.
!  (Careful!  SGE_SL overwrites the right hand side with the solution!)
!
  x = b

  call sge_sl ( lda, n, alu, ipivot, x, job )
!
!  Compute and print the residual.
!
  call sge_res ( lda, n, a, b, job, x, r )

  write ( *, * ) ' '
  write ( *, * ) '  i, x, b-A*x'
  write ( *, * ) ' '

  call rvec2_print_some ( n, x, r, 10 )
!
!  Take a few steps of iterative refinement.
!
  do j = 1, 5

    write ( *, * ) ' '
    write ( *, * ) 'Iterative refinement step ', j
    write ( *, * ) ' '
!
!  Improve the solution.
!
    call sge_sl_it ( lda, n, a, alu, ipivot, b, job, x, r )

    write ( *, * ) ' '
    write ( *, * ) '  i, dx'
    write ( *, * ) ' '

    call rvec_print_some ( n, r, 10 )
!
!  Compute and print the residual.
!
    call sge_res ( lda, n, a, b, job, x, r )

    write ( *, * ) ' '
    write ( *, * ) '  i, x, b-A*x'
    write ( *, * ) ' '

    call rvec2_print_some ( n, x, r, 10 )

  end do

  return
end
subroutine test41
!
!*******************************************************************************
!
!! TEST41 tests SGE_TRF;
!! TEST41 tests SGE_TRS.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
  integer, parameter :: ldb = n
  integer, parameter :: m = n
  integer, parameter :: nrhs = 1
!
  real a(lda,n)
  real b(ldb,nrhs)
  integer i
  integer info
  integer ipivot(n)
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST41'
  write ( *, * ) '  SGE_TRF factors a general linear system,'
  write ( *, * ) '  SGE_TRS solves a factored system.'
  write ( *, * ) ' '
  write ( *, * ) '  Number of matrix rows M = ', m
  write ( *, * ) '  Number of matrix columns N = ', n

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0E+00
      else if ( i == j - 1 ) then
        a(i,j) = - 1.0E+00
      else if ( i == j + 1 ) then
        a(i,j) = - 1.0E+00
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  call sge_trf ( lda, m, n, a, ipivot, info )

  b(1:n-1,1) = 0.0E+00
  b(n,1) = n + 1

  call sge_trs ( lda, n, nrhs, 'N', a, ipivot, b, ldb, info )

  call rvec_print ( n, b, '  Solution:' )

  b(1:n-1,1) = 0.0E+00
  b(n,1) = n + 1

  call sge_trs ( lda, n, nrhs, 'T', a, ipivot, b, ldb, info )

  call rvec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine test42
!
!*******************************************************************************
!
!! TEST42 tests SGE_NP_TRF;
!! TEST42 tests SGE_NP_TRM.
!! TEST42 tests SGE_NP_TRS.
!
  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: nrhs = 1
  integer, parameter :: lda = m
  integer, parameter :: ldb = m
!
  real a(lda,n)
  real b(ldb,nrhs)
  integer i
  integer info
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST42'
  write ( *, * ) '  Using the LAPACK general matrix format,'
  write ( *, * ) '  SGE_NP_TRF factors without pivoting,'
  write ( *, * ) '  SGE_NP_TRS solves factored systems.'
  write ( *, * ) '  SGE_NP_TRM computes A*X for factored A.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix dimensions:'
  write ( *, * ) '    M = ', m 
  write ( *, * ) '    N = ', n
!
!  Set the matrix.
!
  call sge_random ( lda, m, n, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0E+00
!
!  Compute the corresponding right hand side.
!
  call sge_mxv ( lda, m, n, a, x, b )
!
!  Factor the matrix.
!
  call sge_np_trf ( lda, m, n, a, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SGE_NP_TRF declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if
!
!  Solve the linear system.
!
  call sge_np_trs ( lda, n, nrhs, 'N', a, b, n, info )
 
  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call sge_np_trm ( lda, m, n, a, x, b, job )
!
!  Solve the system
!
  call sge_np_trs ( lda, n, nrhs, 'N', a, b, n, info )

  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call sge_np_trm ( lda, m, n, a, x, b, job )
!
!  Solve the system.
!
  call sge_np_trs ( lda, n, nrhs, 'T', a, b, n, info )

  call rvec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine test43
!
!*******************************************************************************
!
!! TEST43 tests SLT_SL;
!! TEST43 tests SLT_MXV.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer j
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST43'
  write ( *, * ) '  For a lower triangular matrix,'
  write ( *, * ) '  SLT_SL solves systems;'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
  do i = 1, n
    do j = 1, n
      if ( j <= i ) then
        a(i,j) = j
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  call slt_print ( lda, n, n, a, '  The lower triangular matrix:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call slt_mxv ( lda, n, n, a, x, b )
!
!  Solve the linear system.
!
  call slt_sl ( lda, n, a, b )
 
  call rvec_print ( n, b, '  Solution:' )

  return
end
subroutine test44
!
!*******************************************************************************
!
!! TEST44 tests SLT_INV.
!! TEST44 tests SLT_DET.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  real det
  integer i
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST44'
  write ( *, * ) '  For a lower triangular matrix,'
  write ( *, * ) '  SLT_INV computes the inverse.'
  write ( *, * ) '  SLT_DET computes the inverse.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
  do i = 1, n
    do j = 1, n
      if ( j <= i ) then
        a(i,j) = j
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  b(1:n,1:n) = a(1:n,1:n)

  call sge_print ( lda, n, n, a, '  Matrix A:' )
!
!  Compute the determinant.
!
  call slt_det ( lda, n, a, det )

  write ( *, * ) ' '
  write ( *, * ) '  Determinant is ', det
!
!  Compute the inverse matrix.
!
  call slt_inv ( lda, n, b )

  call sge_print ( lda, n, n, b, '  Inverse matrix B:' )
!
!  Check
!
  call sge_mxm ( lda, n, a, b, c )

  call sge_print ( lda, n, n, c, '  Product A * B:' )

  return
end
subroutine test45
!
!*******************************************************************************
!
!! TEST45 tests SPB_CG.
!
  integer, parameter :: n = 50
  integer, parameter :: mu = 1
  integer, parameter :: lda = mu + 1
!
  real a(lda,n)
  real b(n)
  real err
  integer i
  real r(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST45'
  write ( *, * ) '  SPB_CG applies the conjugate gradient method'
  write ( *, * ) '    to a symmetric positive definite banded '
  write ( *, * ) '    linear system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Upper bandwidth is MU = ', mu
!
!  Set the matrix values.
!
  a(2,1:n) = 2.0E+00
  a(1,2:n) = -1.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'The symmetric banded matrix:'
  write ( *, * ) ' '

  call spb_print_some ( lda, n, mu, a, 1, 1, 10, 10 )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the right hand side.
!
  call spb_mxv ( lda, n, mu, a, x, b )
!
!  Set the approximate solution.
!
  x(1:n) = 1.0E+00
!
!  Call the conjugate gradient method.
!
  call spb_cg ( lda, n, mu, a, b, x )
!
!  Compute the residual, A*x-b
!
  call spb_mxv ( lda, n, mu, a, x, r )
 
  err = 0.0E+00
  do i = 1, n
    err = max ( err, abs ( r(i) - b(i) ) )
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Solution:'
  write ( *, * ) ' '
  call rvec_print_some ( n, x, 10 )

  write ( *, * ) ' '
  write ( *, * ) '  Maximum residual=', err
 
  return
end
subroutine test46
!
!*******************************************************************************
!
!! TEST46 tests SPB_DET.
!
  integer, parameter :: n = 10
  integer, parameter :: mu = 3
  integer, parameter :: lda = mu + 1
  integer, parameter :: lda2 = n
!
  real a(lda,n)
  real a2(lda2,n)
  real det
  integer info
  integer ipivot(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST46'
  write ( *, * ) '  SPB_DET, determinant of a positive definite'
  write ( *, * ) '    symmetric banded matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Upper bandwidth is MU = ', mu
!
!  Set the matrix.
!
  call spb_random ( lda, n, mu, a )
!
!  Copy the matrix into a general array.
!
  call spb_to_sge ( lda, lda2, n, mu, a, a2 )
!
!  Factor the matrix.
!
  call spb_fa ( lda, n, mu, a, info )
!
!  Compute the determinant.
!
  call spb_det ( lda, n, mu, a, det )

  write ( *, * ) ' '
  write ( *, * ) '  SPB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call sge_fa ( lda2, n, a2, ipivot, info )
!
!  Compute the determinant.
!
  call sge_det ( lda2, n, a2, ipivot, det )

  write ( *, * ) '  SGE_DET computes the determinant = ', det

  return
end
subroutine test47
!
!*******************************************************************************
!
!! TEST47 tests SPB_FA;
!! TEST47 tests SPB_SL.
!
  integer, parameter :: n = 50
  integer, parameter :: mu = 1
  integer, parameter :: lda = mu + 1
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  integer j
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST47'
  write ( *, * ) '  SPB_FA factors a banded positive definite '
  write ( *, * ) '    symmetric linear system.'
  write ( *, * ) '  SPB_SL solves a banded positive definite '
  write ( *, * ) '    symmetric linear system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Upper bandwidth is MU = ', mu
  write ( *, * ) ' '
!
!  Set the matrix values.
!
  a(mu+1,1:n) = 2.0E+00
  a(mu,2:n) = - 1.0E+00
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the right hand side.
!
  call spb_mxv ( lda, n, mu, a, x, b )
!
!  Factor the matrix.
!
  call spb_fa ( lda, n, mu, a, info )
!
!  Solve the linear system.
!
  call spb_sl ( lda, n, mu, a, b )
 
  write ( *, * ) ' '
  write ( *, * ) '  Solution:'
  write ( *, * ) ' '
  call rvec_print_some ( n, b, 10 )
 
  return
end
subroutine test48
!
!*******************************************************************************
!
!! TEST48 tests SPB_ML.
!
  integer, parameter :: n = 10
  integer, parameter :: mu = 3
  integer, parameter :: lda = mu + 1
!
  real a(lda,n)
  real b(n)
  real b2(n)
  integer i
  integer info
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST48'
  write ( *, * ) '  SPB_ML computes A*x '
  write ( *, * ) '    where A has been factored by SPB_FA.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Upper bandwidth is MU = ', mu
!
!  Set the matrix.
!
  call spb_random ( lda, n, mu, a )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call spb_mxv ( lda, n, mu, a, x, b )
!
!  Factor the matrix.
!
  call spb_fa ( lda, n, mu, a, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SPB_FA declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
  call spb_ml ( lda, n, mu, a, x, b2 )

  write ( *, * ) ' '
  write ( *, * ) '  A*x and PLU*x'
  write ( *, * ) ' '
  call rvec2_print_some ( n, b, b2, 10 )

  return
end
subroutine test49
!
!*******************************************************************************
!
!! TEST49 tests SPB_SOR.
!
  integer, parameter :: n = 50
  integer, parameter :: mu = 1
  real, parameter :: pi = 3.14159265E+00
  integer, parameter :: lda = mu + 1
!
  real a(lda,n)
  real b(n)
  real b2(n)
  real eps
  real err
  integer i
  integer itchk
  integer itknt
  integer itmax
  integer k
  real omega
  real t
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST49'
  write ( *, * ) '  SPB_SOR, SOR routine for iterative'
  write ( *, * ) '    solution of A*x=b.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Upper bandwidth is MU = ', mu

  do k = 1, 3
 
    if ( k == 1 ) then
      omega = 0.25E+00
    else if ( k == 2 ) then
      omega = 0.75E+00
    else
      omega = 1.00E+00
    end if
!
!  Set matrix values.
!
    a(2,1:n) = 2.0E+00

    a(1,1) = 0.0E+00
    a(1,2:n) = -1.0E+00
!
!  Set the desired solution.
!
    do i = 1, n
      t = pi * real ( i - 1 ) / real ( n - 1 )
      x(i) = sin ( t )
    end do
!
!  Compute the right hand side.
!
    call spb_mxv (lda, n, mu, a, x, b ) 
!
!  Set the initial solution estimate.
!
    x(1:n) = 1.0E+00
 
    itchk = 1
    itmax = 8000
    eps = 0.0001E+00

    call spb_sor ( lda, n, mu, a, b, eps, itchk, itknt, itmax, omega, x )
!
!  Compute residual, A*x-b
!
    call spb_mxv ( lda, n, mu, a, x, b2 )
 
    err = 0.0E+00
    do i = 1, n
      err = max ( err, abs ( b2(i) - b(i) ) )
    end do
 
    write ( *, * ) ' '
    write ( *, * ) '  SOR iteration.'
    write ( *, * ) ' '
    write ( *, * ) '  Relaxation factor OMEGA =', omega
    write ( *, * ) '  Iterations taken=', itknt
    write ( *, * ) '  Solution:'
    write ( *, * ) ' '
    call rvec_print_some ( n, x, 10 )
    write ( *, * ) ' '
    write ( *, * ) '  Maximum error = ', err
 
  end do
 
  return
end
subroutine test50
!
!*******************************************************************************
!
!! TEST50 tests SPO_FA;
!! TEST50 tests SPO_SL.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer info
  integer j
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST50'
  write ( *, * ) '  SPO_FA factors a positive definite symmetric'
  write ( *, * ) '    linear system,'
  write ( *, * ) '  SPO_SL solves a factored system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
  do i = 1, n
    do j = 1, n
      a(i,j) = min ( i, j )
    end do
  end do
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sge_mxv ( lda, n, n, a, x, b )
!
!  Factor the matrix.
!
  call spo_fa ( lda, n, a, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Fatal error!'
    write ( *, * ) '  SPO_FA declares the matrix is singular!'
    write ( *, * ) '  The value of INFO is ',info
    return
  end if
!
!  Solve the linear system.
!
  call spo_sl ( lda, n, a, b )
 
  call rvec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  x(1:n) = 1
!
!  Compute the corresponding right hand side.
!
  call spo_ml ( lda, n, a, x, b )
!
!  Solve the linear system.
!
  call spo_sl ( lda, n, a, b )
 
  call rvec_print ( n, b, '  Solution:' )

  return
end
subroutine test51
!
!*******************************************************************************
!
!! TEST51 tests SPO_DET.
!! TEST51 tests SPO_INV.
!
  integer, parameter :: n = 4
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  real det
  integer i
  integer info
  integer j
  integer k
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST51'
  write ( *, * ) '  For a symmetric positive definite matrix'
  write ( *, * ) '    factored by SPO_FA,'
  write ( *, * ) '  SPO_DET computes the determinant;'
  write ( *, * ) '  SPO_INV computes the inverse.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = min ( i,j )
    end do
  end do

  b(1:n,1:n) = a(1:n,1:n)

  call sge_print ( lda, n, n, a, '  Matrix A:' )
!
!  Factor the matrix.
!
  call spo_fa ( lda, n, b, info )
!
!  Compute the determinant.
!
  call spo_det ( lda, n, b, det )
 
  write ( *, * ) ' '
  write ( *, * ) '  Matrix determinant = ', det
!
!  Compute the inverse.
!
  call spo_inv ( lda, n, b )

  call sge_print ( lda, n, n, b, '  Inverse matrix B:' )
!
!  Check.
!
  call sge_mxm ( lda, n, a, b, c )

  call sge_print ( lda, n, n, c, '  Product A * B:' )

  return
end
subroutine test52
!
!*******************************************************************************
!
!! TEST52 tests SPO_RANDOM.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST52'
  write ( *, * ) '  SPO_RANDOM, compute a random positive definite'
  write ( *, * ) '    symmetric matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call spo_random ( lda, n, a )

  call sge_print ( lda, n, n, a, '  Random matrix:' )
 
  return
end
subroutine test53
!
!*******************************************************************************
!
!! TEST53 tests SPP_RANDOM.
!
  integer, parameter :: n = 5
!
  real a((n*(n+1))/2)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST53'
  write ( *, * ) '  SPP_RANDOM, compute a random positive definite'
  write ( *, * ) '    symmetric packed matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call spp_random ( n, a )

  call spp_print ( n, a, '  The matrix:' )
 
  return
end
subroutine test54
!
!*******************************************************************************
!
!! TEST54 tests SSD_CG.
!
!
!  NX and NY are the number of grid points in the X and Y directions.
!  N is the number of unknowns.
!  NDIAG is the number of nonzero diagonals we will store.  We only
!    store the main diagonal, and the superdiagonals.
!  LDA is the leading dimension of the array A, which can be N.
!
  integer, parameter :: ndiag = 3
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10
  integer, parameter :: n = nx * ny
  integer, parameter :: lda = n
!
  real a(lda,ndiag)
  real ap(n)
  real b(n)
  real b2(n)
  real err
  integer i
  integer j
  integer k
  integer offset(ndiag)
  real p(n)
  real r(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST54'
  write ( *, * ) '  SSD_CG applies the conjugate gradient method'
  write ( *, * ) '    to a symmetric positive definite linear'
  write ( *, * ) '    system stored by diagonals.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
  write ( *, * ) '  Number of diagonals is ', ndiag
  write ( *, * ) ' '
!
!  OFFSET tells us where the nonzero diagonals are.  It does this
!  by recording how "high" or to the right the diagonals are from
!  the main diagonal.
!
  offset(1) =   0
  offset(2) =   1
  offset(3) =  nx
!
!  Now we compute the numbers that go into the diagonals.  For this
!  problem, we could simply store a column of 4's, and two columns of
!  -1's, but I wanted to go through the motions of thinking about the
!  value of each entry.  "K" counts the row of the original matrix
!  that we are working on.
!
  k = 0
  do j = 1, ny
    do i = 1, nx

      k = k + 1
!
!  Central
!
      a(k,1) = 4.0E+00
!
!  East ( = West )
!
      if ( i == nx ) then
        a(k,2) = 0.0E+00
      else
        a(k,2) = -1.0E+00
      end if
!
!  North ( = South )
!
      if ( j == ny ) then
        a(k,3) = 0.0E+00
      else
        a(k,3) = -1.0E+00
      end if

    end do
  end do
!
!  Print some of the matrix.
!
  write ( *, * ) ' '
  write ( *, * ) '  First 10 rows and columns of matrix.'
  write ( *, * ) ' '

  call ssd_print_some ( lda, n, ndiag, offset, a, 1, 1, 10, 10 )
!
!  Set the desired solution.
!
  k = 0
  do j = 1, ny
    do i = 1, nx
      k = k + 1
      x(k) = 10 * i + j
    end do
  end do
!
!  Compute the corresponding right hand side.
!
  call ssd_mxv ( lda, n, ndiag, offset, a, x, b )

  write ( *, * ) ' '
  write ( *, * ) '  Right hand side:'
  write ( *, * ) ' '
  call rvec_print_some ( n, b, 10 )
!
!  Set X to zero so no one accuses us of cheating.
!
  x(1:n) = 0.0E+00
!
!  Call the conjugate gradient method.
!
  call ssd_cg ( lda, n, ndiag, offset, a, b, x )
!
!  Compute the residual, A*x-b
!
  call ssd_mxv ( lda, n, ndiag, offset, a, x, b2 )
 
  err = 0.0E+00
  do i = 1, n
    err = max ( err, abs ( b2(i) - b(i) ) )
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Solution:'
  write ( *, * ) ' '
  call rvec_print_some ( n, x, 10 )

  write ( *, * ) ' '
  write ( *, * ) '  Maximum residual=', err
!
!  Note that if we're not satisfied with the solution, we can
!  call again, using the computed X as our starting estimate.
!
!
!  Call the conjugate gradient method AGAIN.
!
  call ssd_cg ( lda, n, ndiag, offset, a, b, x )
!
!  Compute the residual, A*x-b
!
  call ssd_mxv ( lda, n, ndiag, offset, a, x, b2 )
 
  err = 0.0
  do i = 1, n
    err = max ( err, abs ( b2(i) - b(i) ) )
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Second attempt at solution:'
  write ( *, * ) ' '
  call rvec_print_some ( n, x, 10 )

  write ( *, * ) ' '
  write ( *, * ) '  Maximum residual of second attempt =', err

  return
end
subroutine test55
!
!*******************************************************************************
!
!! TEST55 tests SSD_CG.
!
!
!  This is a sample demonstration of how to compute some eigenvalues
!  and corresponding eigenvectors of a matrix.  The matrix is the
!  discretized Laplacian operator, which can be stored by diagonals,
!  and handled by the conjugate gradient method.
!
  integer, parameter :: maxvec = 3
  integer, parameter :: ndiag = 3
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10
  integer, parameter :: n = nx * ny
  integer, parameter :: lda = n
  real, parameter :: pi = 3.141592653589
!
  real a(lda,ndiag)
  real ap(n)
  real del
  real dot
  real eval
  integer i
  integer iter
  integer ivec
  integer j
  integer k
  real lambda
  real lambda_old
  real lamvec(maxvec)
  real norm
  integer nvec
  integer offset(ndiag)
  real p(n)
  real r(n)
  real vec(n,maxvec)
  real x(n)
  real xnew(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST55'
  write ( *, * ) '  SSD_CG is used for linear equation solving'
  write ( *, * ) '    in a demonstration of inverse iteration to'
  write ( *, * ) '    compute eigenvalues and eigenvectors of a'
  write ( *, * ) '    symmetric matrix stored by diagonals.'
  write ( *, * ) ' '
  write ( *, * ) '  Problem size is N = ', n
  write ( *, * ) ' '
  write ( *, * ) '  Here are 25 of the smallest eigenvalues:'
  write ( *, * ) ' '
  write ( *, * ) '  I, J, eigenvalue(I,J):'
  write ( *, * ) ' '

  do i = 1, min ( 5, nx )
    do j = 1, min ( 5, ny )
      eval = 4.0 - 2.0 * cos ( real ( i ) * pi / dble ( nx + 1 ) ) &
                 - 2.0 * cos ( real ( j ) * pi / dble ( ny + 1 ) )
      write ( *, '(2i6,g14.6)' ) i, j, eval
    end do
  end do
!
!  OFFSET tells us where the nonzero diagonals are.  It does this
!  by recording how "high" or to the right the diagonals are from
!  the main diagonal.
!
  offset(1) =   0
  offset(2) =   1
  offset(3) =  nx
!
!  Now we compute the numbers that go into the diagonals.  For this
!  problem, we could simply store a column of 4's, and two columns of
!  -1's, but I wanted to go through the motions of thinking about the
!  value of each entry.  "K" counts the row of the original matrix
!  that we are working on.
!
  k = 0
  do j = 1, ny
    do i = 1, nx

      k = k + 1
!
!  Central
!
      a(k,1) = 4.0E+00
!
!  East ( = West )
!
      if ( i == nx ) then
        a(k,2) = 0.0E+00
      else
        a(k,2) = -1.0E+00
      end if
!
!  North ( = South )
!
      if ( j == ny ) then
        a(k,3) = 0.0E+00
      else
        a(k,3) = -1.0E+00
      end if

    end do
  end do

  nvec = 0
!
!  Set the starting eigenvector and eigenvalue estimates.
!
99    continue

  write ( *, * ) ' '

  lambda = 0.0E+00

  k = 0
  do j = 1, ny
    do i = 1, nx
      k = k + 1
      x(k) = 1.0E+00
    end do
  end do
!
!  Remove any components of previous eigenvectors.
!
  do ivec = 1, nvec

    dot = 0.0E+00
    do i = 1, n
      dot = dot + x(i) * vec(i,ivec)
    end do

    do i = 1, n
      x(i) = x(i) - dot * vec(i,ivec)
    end do

  end do

  xnew(1:n) = x(1:n)
!
!  Iterate
!
  do iter = 1, 40

    norm = sqrt ( sum ( xnew(1:n)**2 ) )

    xnew(1:n) = xnew(1:n) / norm

    lambda_old = lambda
    lambda = 1.0E+00 / norm
!
!  Check for convergence.
!
    if ( iter > 1 ) then
      del = abs ( lambda - lambda_old )
      if ( del < 0.000001E+00 ) then
        write ( *, * ) ' '
        write ( *, * ) 'Lambda estimate = ', lambda
        write ( *, * ) 'Converged on step ', iter
        go to 10
      end if
    end if
!
!  Call the conjugate gradient method, solving
!    A * XNEW = X.
!
    x(1:n) = xnew(1:n)

    call ssd_cg ( lda, n, ndiag, offset, a, x, xnew  )

    do ivec = 1, nvec

      dot = 0.0
      do i = 1, n
        dot = dot + xnew(i) * vec(i,ivec)
      end do

      do i = 1, n
        xnew(i) = xnew(i) - dot * vec(i,ivec)
      end do

    end do

  end do

  write ( *, * ) 'Did not converge, gave up on step ', iter

10    continue

  nvec = nvec + 1
  lamvec(nvec) = lambda
  vec(1:n,nvec) = xnew(1:n)

  if ( ivec < maxvec ) then
    go to 99
  end if

  return
end
subroutine test56
!
!*******************************************************************************
!
!! TEST56 tests SSM_ML.
!
  integer, parameter :: n = 7
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  real b2(n)
  integer info
  integer i
  integer ipivot(n)
  integer job
  real u(n)
  real v(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST56'
  write ( *, * ) '  SSM_ML computes A*x or transpose(A)*X'
  write ( *, * ) '    where A is a Sherman Morrison matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call ssm_random ( lda, n, a, u, v )

    call ssm_print ( lda, n, a, u, v, '  The Sherman Morrison matrix:' )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call ssm_mxv ( lda, n, a, u, v, x, b )
    else
      call ssm_vxm ( lda, n, a, u, v, x, b )
    end if
!
!  Factor the matrix.
!
    call sge_fa ( lda, n, a, ipivot, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SGE_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call ssm_ml ( lda, n, a, u, v, ipivot, x, b2, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  A*x and PLU*x'
    else
      write ( *, * ) '  AT*x and (PLU)T*x'
    end if

    write ( *, * ) ' '
    call rvec2_print_some ( n, b, b2, 10 )

  end do

  return
end
subroutine test57
!
!*******************************************************************************
!
!! TEST57 tests SSM_SL.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer ierror
  integer info
  integer ipivot(n)
  integer job
  real u(n)
  real v(n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST57'
  write ( *, * ) '  SSM_SL implements the Sherman-Morrison method '
  write ( *, * ) '    for solving a perturbed linear system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  do job = 1, 0, -1
!
!  Set the matrix.
!
    call ssm_random ( lda, n, a, u, v )

    call ssm_print ( lda, n, a, u, v, '  The Sherman-Morrison matrix A:' )
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call ssm_mxv ( lda, n, a, u, v, x, b )
    else
      call ssm_vxm ( lda, n, a, u, v, x, b )
    end if

    call rvec_print ( n, b, '  The right hand side vector B:' )
!
!  Factor the matrix.
!
    call sge_fa ( lda, n, a, ipivot, info )

    if ( info /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Fatal error!'
      write ( *, * ) '  SGE_FA declares the matrix is singular!'
      write ( *, * ) '  The value of INFO is ',info
      cycle
    end if
!
!  Solve the linear system.
!
    call ssm_sl ( lda, n, a, u, v, b, ierror, ipivot, job )
 
    if ( job == 0 ) then
      call rvec_print ( n, b, '  Solution to A * X = B:' )
    else
      call rvec_print ( n, b, '  Solution to At * X = B:' )
    end if
 
  end do

  return
end
subroutine test58
!
!*******************************************************************************
!
!! TEST58 tests SSS_MXV.
!! TEST58 tests SSS_PRINT.
!
  integer, parameter :: n = 9
  integer, parameter :: lda = n
!
  real a((n*(n+1))/2)
  real a2(lda,n)
  real b(n)
  real b2(n)
  integer diag(n)
  integer i
  integer ij
  integer ilo
  integer j
  integer na
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST58'
  write ( *, * ) '  For a symmetric skyline storage matrix,'
  write ( *, * ) '  SSS_MXV computes A*x,'
  write ( *, * ) '  SSS_PRINT prints it.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sss_random ( n, na, a, diag )

  write ( *, * ) ' '
  write ( *, * ) '  Number of nonzero entries stored is ', na
  write ( *, * ) ' '
  write ( *, * ) '  Diagonal storage indices:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, * ) i, diag(i)
  end do
!
!  Replace the random entries by marker values.
!
  ij = 0
  do j = 1, n

    if ( j == 1 ) then
      ilo = 1
    else
      ilo = diag(j-1)-diag(j)+j+1
    end if

    do i = ilo, j
      ij = ij + 1
      a(ij) = 10 * i + j
    end do

  end do

  call sss_print ( n, na, a, diag, '  The symmetric skyline storage matrix:' )
!
!  Copy the matrix into a general matrix.
!
  call sss_to_sge ( lda, n, na, a, diag, a2 )
!
!  Set the vector X.
!
  call rvec_identity ( n, x )
!
!  Compute the product.
!
  call sss_mxv ( diag, n, na, a, x, b )
!
!  Compute the product using the general matrix.
!
  call sge_mxv ( lda, n, n, a2, x, b2 )
!
!  Compare the results.
!
  write ( *, * ) ' '
  write ( *, * ) '  SSS_MXV verse SGE_MXV'
  write ( *, * ) ' '
  call rvec2_print_some ( n, b, b2, 10 )

  return
end
subroutine test583
!
!*******************************************************************************
!
!! TEST583 tests SSTO_INV.
!
  integer, parameter :: n = 3
  integer, parameter :: lda = n
!
  real, dimension ( n ) :: a = (/ 1.0E+00, 0.5E+00, 0.2E+00 /)
  real a2(lda,n)
  real b(n,n)
  real c(n,n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST583'
  write ( *, * ) '  SSTO_INV computes the inverse of a positive definite '
  write ( *, * ) '  symmetric Toeplitz matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  call ssto_print ( n, a, '  The symmetric Toeplitz matrix A:' )

  call ssto_inv ( n, a, b )

  call sge_print ( n, n, n, b, '  The inverse matrix B:' )

  call ssto_to_sge ( lda, n, a, a2 )

  c(1:n,1:n) = matmul ( a2(1:n,1:n), b(1:n,1:n) )

  call sge_print ( n, n, n, c, '  The product C = A * B:' )

  return
end
subroutine test585
!
!*******************************************************************************
!
!! TEST585 tests SSTO_MXV.
!! TEST585 tests SSTO_YW_SL.
!
  integer, parameter :: n = 3
!
  real a(n)
  real b(n)
  integer i
  integer job
  real, dimension ( 0:n ) :: r = (/ 1.0E+00, 0.5E+00, 0.2E+00, 0.1E+00 /)
  real x(n)
!
  a(1:n) = r(0:n-1)

  write ( *, * ) ' '
  write ( *, * ) 'TEST585'
  write ( *, * ) '  SSTO_YW_SL solves the Yule-Walker equations for a'
  write ( *, * ) '  symmetric Toeplitz system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  call ssto_print ( n, a, '  The symmetric Toeplitz matrix:' )

  b(1:n) = -r(1:n)
  call rvec_print ( n, b, '  The right hand side, B:' )

  b(1:n) = -b(1:n)
  call ssto_yw_sl ( n, b, x )

  call rvec_print ( n, x, '  The computed solution, X:' )

  call ssto_mxv ( n, a, x, b )

  call rvec_print ( n, b, '  The product A*X:' )

  return
end
subroutine test587
!
!*******************************************************************************
!
!! TEST587 tests SSTO_SL.
!
  integer, parameter :: n = 3
!
  real, dimension ( 0:n-1 ) :: a = (/ 1.0E+00, 0.5E+00, 0.2E+00 /)
  real, dimension ( n ) :: b = (/ 4.0E+00, -1.0E+00, 3.0E+00 /)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST587'
  write ( *, * ) '  SSTO_SL solves a positive definite symmetric Toeplitz system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  call ssto_print ( n, a, '  The symmetric Toeplitz matrix A:' )

  call rvec_print ( n, b, '  The right hand side vector B:' )

  call ssto_sl ( n, a, b, x )

  call rvec_print ( n, x, '  The solution X:' )

  call ssto_mxv ( n, a, x, b )

  call rvec_print ( n, b, '  The product vector  B = A * X:' )

  return
end
subroutine test59
!
!*******************************************************************************
!
!! TEST59 tests STO_SL.
!
  integer, parameter :: n = 10
!
  real a(2*n-1)
  real b(n)
  integer i
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST59'
  write ( *, * ) '  STO_SL solves a Toeplitz system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call sto_random ( n, a )

  call sto_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call sto_mxv ( n, a, x, b )
    else
      call sto_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call sto_sl ( n, a, b, x, job )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  Solution:'
    else
      write ( *, * ) '  Solution to transposed system:'
    end if

    write ( *, * ) ' '
    call rvec_print_some ( n, x, 10 )

  end do
 
  return
end
subroutine test60
!
!*******************************************************************************
!
!! TEST60 tests SUT_SL;
!! TEST60 tests SUT_MXV.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  integer i
  integer j
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST60'
  write ( *, * ) '  For an upper triangular matrix,'
  write ( *, * ) '  SUT_SL solves systems;'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
  do i = 1, n
    do j = 1, n
      if ( j >= i ) then
        a(i,j) = j
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  call sut_print ( lda, n, n, a, '  The upper triangular matrix:' )
!
!  Set the desired solution.
!
  call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
  call sut_mxv ( lda, n, n, a, x, b )
!
!  Solve the linear system.
!
  call sut_sl ( lda, n, a, b )
 
  call rvec_print ( n, b, '  Solution:' )

  return
end
subroutine test61
!
!*******************************************************************************
!
!! TEST61 tests SUT_INV.
!! TEST61 tests SUT_DET.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  real det
  integer i
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST61'
  write ( *, * ) '  For an upper triangular matrix,'
  write ( *, * ) '  SUT_INV computes the inverse.'
  write ( *, * ) '  SUT_DET computes the determinant.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n

  do i = 1, n
    do j = 1, n
      if ( j >= i ) then
        a(i,j) = j
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  call sut_print ( lda, n, n, a, '  The matrix A:' )
!
!  Compute the determinant.
!
  call sut_det ( lda, n, a, det )

  write ( *, * ) ' '
  write ( *, * ) '  Determinant is ', det
!
!  Compute the inverse matrix B.
!
  b(1:n,1:n) = a(1:n,1:n)

  call sut_inv ( lda, n, b )

  call sut_print ( lda, n, n, b, '  The inverse matrix B:' )
!
!  Check
!
  call sge_mxm ( lda, n, a, b, c )

  call sge_print ( lda, n, n, c, '  The product A * B:' )

  return
end
subroutine test62
!
!*******************************************************************************
!
!! TEST62 tests SVM_DET.
!
  integer, parameter :: n = 10
  integer, parameter :: lda = n
!
  real a(n)
  real a2(lda,n)
  real det
  integer info
  integer ipivot(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST62'
  write ( *, * ) '  SVM_DET, determinant of a Vandermonde matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call svm_random ( n, a )

  call svm_print ( n, a, '  The Vandermonde matrix:' )
!
!  Copy the matrix into a general array.
!
  call svm_to_sge ( lda, n, a, a2 )
!
!  Compute the determinant.
!
  call svm_det ( n, a, det )

  write ( *, * ) ' '
  write ( *, * ) '  SVM_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call sge_fa ( lda, n, a2, ipivot, info )
!
!  Compute the determinant.
!
  call sge_det ( lda, n, a2, ipivot, det )

  write ( *, * ) '  SGE_DET computes the determinant = ', det

  return
end
subroutine test63
!
!*******************************************************************************
!
!! TEST63 tests SVM_SL.
!
  integer, parameter :: n = 10
!
  real a(n)
  real b(n)
  integer i
  integer info
  integer job
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST63'
  write ( *, * ) '  SVM_SL solves a Vandermonde system.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call svm_random ( n, a )

  do job = 0, 1
!
!  Set the desired solution.
!
    call rvec_identity ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call svm_mxv ( n, a, x, b )
    else
      call svm_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call svm_sl ( n, a, b, x, job, info )

    write ( *, * ) ' '
    if ( job == 0 ) then
      write ( *, * ) '  Solution:'
    else
      write ( *, * ) '  Solution to transposed system:'
    end if
    write ( *, * ) ' '
    call rvec_print_some ( n, x, 10 )

  end do
 
  return
end
