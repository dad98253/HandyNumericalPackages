!  eisprb.f90  15 June 2000
!
program eisprb
!
!*******************************************************************************
!
!! EISPRB calls the EISPACK sample programs.
!
  write ( *, * ) ' '
  write ( *, * ) 'EISPRB'
  write ( *, * ) '  Sample problems for EISPACK'

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
  call test145
  call test15

  write ( *, * ) ' '
  write ( *, * ) 'EISPRB'
  write ( *, * ) '  Normal end of EISPACK tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests CG.
!
!
!  CG is for the eigenvalues of a complex general matrix.
!
!  eigenvalues and eigenvectors of a complex general matrix
!  note that the eigenvalues of such a matrix are in general complex.
!  however, we will use the same example we used before, namely
!  a hermitian matrix, so the eigenvalues will in fact be real.
!
!  (3     1     0     0+2i)
!  (1     3     0-2i  0   )
!  (0     0+2i  1     1   )
!  (0-2i  0     1     1   )
!
!  The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!  The eigenvector matrix is
!
!  (  1+sqrt(2),  1,                -1,          1)
!  (  1+sqrt(2),  1,                 1,         -1)
!  (     i,       -(1+sqrt(2))*i,    i,          i)
!  (    -i,        (1+sqrt(2))*i,    i,          i)
!
!  Note that the actual eigenvector matrix from EISPACK could
!  be scaled by a real value, or by i, and the columns may
!  appear in any order.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
!
  real ar(n,n)
  real ai(n,n)
  integer i
  integer ierr
  integer j
  integer matz
  real wr(n)
  real wi(n)
  real xr(n,n)
  real xi(n,n)
!
!  Set the values of the matrix.
!
  ar(1,1) = 3.0
  ar(1,2) = 1.0
  ar(1,3) = 0.0
  ar(1,4) = 0.0
 
  ar(2,1) = 1.0
  ar(2,2) = 3.0
  ar(2,3) = 0.0
  ar(2,4) = 0.0
 
  ar(3,1) = 0.0
  ar(3,2) = 0.0
  ar(3,3) = 1.0
  ar(3,4) = 1.0
 
  ar(4,1) = 0.0
  ar(4,2) = 0.0
  ar(4,3) = 1.0
  ar(4,4) = 1.0

  ai(1,1) = 0.0
  ai(1,2) = 0.0
  ai(1,3) = 0.0
  ai(1,4) = 2.0
 
  ai(2,1) = 0.0
  ai(2,2) = 0.0
  ai(2,3) = -2.0
  ai(2,4) = 0.0
 
  ai(3,1) = 0.0
  ai(3,2) = 2.0
  ai(3,3) = 0.0
  ai(3,4) = 0.0
 
  ai(4,1) = -2.0
  ai(4,2) = -0.0
  ai(4,3) = -0.0
  ai(4,4) = 0.0
!
!  matz = 0 for eigenvalues only,
!  matz = 1 for eigenvalues and eigenvectors.
!
  matz = 1
  call cg ( nm, n, ar, ai, wr, wi, matz, xr, xi, ierr )

  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  CG computes the eigenvalues and eigenvectors of '
  write ( *, * ) '  a complex general matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n
  write ( *, * ) ' '
  write ( *, * ) '  Error flag = ',ierr

  call rvec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )
 
  if ( matz /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'The eigenvectors are:'
    do i = 1, n
      write ( *, * ) ' '
      write ( *, * ) '  Eigenvector ',I
      write ( *, * ) ' '
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if
 
  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests CH.
!
!
!  CH is for the eigenvalues of a complex hermitian matrix.
!
!  eigenvalues and eigenvectors of a complex hermitian matrix
!  note that the eigenvalues (though not the eigenvectors) of
!  a hermitian matrix are real.
!
!  (3     1     0     0+2i)
!  (1     3     0-2i  0   )
!  (0     0+2i  1     1   )
!  (0-2i  0     1     1   )
!
!  the eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!  the eigenvector matrix is
!
!  (  1+sqrt(2),  1,                -1,          1)
!  (  1+sqrt(2),  1,                 1,         -1)
!  (     i,       -(1+sqrt(2))*i,    i,          i)
!  (    -i,        (1+sqrt(2))*i,    i,          i)
!
!  Note that the actual eigenvector matrix from EISPACK could
!  be scaled by a real value, or by i, and the columns may
!  appear in any order.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
!
  real ar(n,n)
  real ai(n,n)
  integer i
  integer ierr
  integer j
  integer matz
  real w(n)
  real xr(n,n)
  real xi(n,n)
!
!  Set the values of the matrix.
!
  ar(1,1) = 3.0
  ar(1,2) = 1.0
  ar(1,3) = 0.0
  ar(1,4) = 0.0
 
  ar(2,1) = 1.0
  ar(2,2) = 3.0
  ar(2,3) = 0.0
  ar(2,4) = 0.0
 
  ar(3,1) = 0.0
  ar(3,2) = 0.0
  ar(3,3) = 1.0
  ar(3,4) = 1.0
 
  ar(4,1) = 0.0
  ar(4,2) = 0.0
  ar(4,3) = 1.0
  ar(4,4) = 1.0

  ai(1,1) = 0.0
  ai(1,2) = 0.0
  ai(1,3) = 0.0
  ai(1,4) = 2.0
 
  ai(2,1) = 0.0
  ai(2,2) = 0.0
  ai(2,3) = -2.0
  ai(2,4) = 0.0
 
  ai(3,1) = 0.0
  ai(3,2) = 2.0
  ai(3,3) = 0.0
  ai(3,4) = 0.0
 
  ai(4,1) = -2.0
  ai(4,2) = -0.0
  ai(4,3) = -0.0
  ai(4,4) = 0.0

  matz = 1

  call ch ( nm, n, ar, ai, w, matz, xr, xi, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  CH computes the eigenvalues and eigenvectors of'
  write ( *, * ) '  a complex hermitian matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n
  write ( *, * ) ' '
  write ( *, * ) '  Error flag = ',ierr
  write ( *, * ) ' '

  call rvec_print ( n, w, '  The eigenvalues Lambda:' )

  write ( *, * ) ' '
 
  if ( matz /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Eigenvectors are:'
    do i = 1, n
      write ( *, * ) ' '
      write ( *, * ) '  Eigenvector ',I
      write ( *, * ) ' '
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if
 
  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests MINFIT.
!
  integer, parameter :: m = 5
  integer, parameter :: nb = 1
  integer, parameter :: nm = m
  integer, parameter :: n = 2
!
  real a(m,n)
  real acopy(m,n)
  real b(m,nb)
  integer i
  integer ierr
  integer j
  real r(m)
  real w(n)
  real x(n)
!
  a(1,1) =   1.00
  a(2,1) =   2.05
  a(3,1) =   3.06
  a(4,1) = - 1.02
  a(5,1) =   4.08

  a(1,2) =   1.00
  a(2,2) = - 1.00
  a(3,2) =   1.00
  a(4,2) =   2.00
  a(5,2) = - 1.00

  acopy(1:m,1:n) = a(1:m,1:n)

  b(1,1) = 1.98
  b(2,1) = 0.95
  b(3,1) = 3.98
  b(4,1) = 0.92
  b(5,1) = 2.90

  do i = 1, m
    r(i) = - b(i,1)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  MINFIT solves an overdetermined linear system'
  write ( *, * ) '  using least squares methods.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix rows = ', m
  write ( *, * ) '  Matrix columns = ', n

  call rmat_print ( nm, m, n, a, '  The matrix A:' )

  call rmat_print ( nm, m, nb, b, '  The right hand side B:' )

  call minfit ( nm, m, n, a, w, nb, b, ierr )

  write ( *, * ) ' '
  write ( *, * ) '  MINFIT error code IERR = ', ierr

  call rvec_print ( n, w, '  The singular values:' )
!
!  B now contains U' * B.
!  We need to divide by the singular values, and multiply by V.
!
  b(1:n,1) = b(1:n,1) / w(1:n)

  do i = 1, n
    x(i) = 0.0
    do j = 1, n
      x(i) = x(i) + a(i,j) * b(j,1)
    end do
  end do

  call rvec_print ( n, x, '  The least squares solution X:' )

  do i = 1, m
    do j = 1, n
      r(i) = r(i) + acopy(i,j) * x(j)
    end do
  end do
    
  call rvec_print ( m, r, '  The residual A * X - B:' )

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests RG.
!
!
!  RG is for the eigenvalues of a general real matrix.
!
!  The matrix A is nonsymmetric.  The eigenvalues may therefore be
!  complex numbers.
!
!  ( 33  16  72)
!  (-24 -10 -57)
!  ( -8  -4 -17)
!
!  The eigenvalues of A are (1,2,3)
!
!  The eigenvectors of A are
!
!  (-15 -16 -4)
!  ( 12  13  3)
!  (  4   4  1)
!
  integer, parameter :: n = 3 
  integer, parameter :: nm = n
!
  real a(n,n)
  real acopy(n,n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real sum
  real sum1
  real sum2
  real wi(n)
  real wr(n)
  real x(n,n)
!
!  Set the values of the matrix.
!
  a(1,1) = 33.0
  a(1,2) = 16.0
  a(1,3) = 72.0
 
  a(2,1) = -24.0
  a(2,2) = -10.0
  a(2,3) = -57.0
 
  a(3,1) = -8.0
  a(3,2) = -4.0
  a(3,3) = -17.0
!
!  RG overwrites A with garbage, so save a copy now!
!
  acopy(1:n,1:n) = a(1:n,1:n)
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  RG computes the eigenvalues and eigenvectors of'
  write ( *, * ) '  a real general matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a, '  The matrix A:' )

  matz = 1

  call rg ( nm, n, a, wr, wi, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag was IERR = ',ierr
  end if

  call rvec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )
 
  if ( matz /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  The eigenvectors may be complex:'
    do j = 1, n
      write ( *, * ) ' '
      write ( *, * ) '  Eigenvector ',j
      write ( *, * ) ' '
      do i = 1, n
        if ( wi(j) == 0.0 ) then
          write ( *, '(g14.6)' ) x(i,j)
        else if ( wi(j) > 0.0 ) then
          write ( *, '(2g14.6)' ) x(i,j), x(i,j+1)
        else if ( wi(j) < 0.0 ) then
          write ( *, '(2g14.6)' ) x(i,j-1), -x(i,j)
        end if
      end do
    end do
!
!  Check.
!  First, restore the original values of A.
!
    a(1:n,1:n) = acopy(1:n,1:n)
 
    do k = 1, n
      write ( *, * ) ' '
      write ( *, * ) 'Residuals (A*x-Lambda*x) for eigenvalue ',k
      write ( *, * ) ' '
 
      if ( wi(k)==0.0 ) then
 
        do i = 1, n
          sum = 0.0
          do j = 1, n
            sum = sum+a(i,j)*x(j,k)
          end do
          sum = sum-wr(k)*x(i,k)
          write ( *, '(g14.6)' ) sum
        end do
 
      else if ( wi(k)>0.0 ) then
 
        do i = 1, n
          sum1 = 0.0
          sum2 = 0.0
          do j = 1, n
            sum1 = sum1+a(i,j)*x(j,k)
            sum2 = sum2+a(i,j)*x(j,k+1)
          end do
          sum1 = sum1-wr(k)*x(i,k)+wi(k)*x(i,k+1)
          sum2 = sum2-wi(k)*x(i,k)-wr(k)*x(i,k+1)
          write ( *, '(2g14.6)' ) sum1, sum2
        end do
 
      else if ( wi(k)<0.0 ) then
 
        do i = 1, n
          sum1 = 0.0
          sum2 = 0.0
          do j = 1, n
            sum1 = sum1+a(i,j)*x(j,k-1)
            sum2 = sum2-a(i,j)*x(j,k)
          end do
          sum1 = sum1-wr(k)*x(i,k-1)-wi(k)*x(i,k)
          sum2 = sum2-wi(k)*x(i,k-1)+wr(k)*x(i,k)
          write ( *, '(2g14.6)' ) sum1, sum2
        end do
 
      end if
 
    end do
 
  end if
 
  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests RGG.
!
!
!  RGG is for a real generalized general eigenvalue problem.
!
!  A generalized eigenvalue problem.  Given matrices A and B, find
!  N numbers LAMBDA, and for each LAMBDA a vector X, so that
!
!    A*x = lambda*B*x
!
!  The matrix A is
!
!  ( -7 7  6  6)
!  (-10 8 10  8)
!  ( -8 3 10 11)
!  ( -4 0  4 12)
!
!  The matrix B is
!
!  (2 1 0 0)
!  (1 2 1 0)
!  (0 1 2 1)
!  (0 0 1 2)
!
!  The correct eigenvalues LAMBDA are
!
!  (1,2,3,4)
!
!  The correct eigenvectors X are
!
!  (4 3 2 1)
!  (3 3 2 1)
!  (2 2 2 1)
!  (1 1 1 1)
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
!
  real a(n,n)
  real acopy(n,n)
  real alfi(n)
  real alfr(n)
  real b(n,n)
  real bcopy(n,n)
  real beta(n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real sum
  real sum1
  real sum2
  real x(n,n)
!
!  Set the values in the A matrix.
!
  a(1,1) = -7.0
  a(1,2) = 7.0
  a(1,3) = 6.0
  a(1,4) = 6.0
 
  a(2,1) = -10.0
  a(2,2) = 8.0
  a(2,3) = 10.0
  a(2,4) = 8.0
 
  a(3,1) = -8.0
  a(3,2) = 3.0
  a(3,3) = 10.0
  a(3,4) = 11.0
 
  a(4,1) = -4.0
  a(4,2) = 0.0
  a(4,3) = 4.0
  a(4,4) = 12.0
!
!  Save a copy of A.
!
  acopy(1:n,1:n) = a(1:n,1:n)
!
!  Set the values in the B matrix.
!
  b(1,1) = 2.0
  b(1,2) = 1.0
  b(1,3) = 0.0
  b(1,4) = 0.0
 
  b(2,1) = 1.0
  b(2,2) = 2.0
  b(2,3) = 1.0
  b(2,4) = 0.0
 
  b(3,1) = 0.0
  b(3,2) = 1.0
  b(3,3) = 2.0
  b(3,4) = 1.0
 
  b(4,1) = 0.0
  b(4,2) = 0.0
  b(4,3) = 1.0
  b(4,4) = 2.0
!
!  Save a copy of B.
!
  bcopy(1:n,1:n) = b(1:n,1:n)

  write ( *, * ) ' '
  write ( *, * ) 'TEST05:'
  write ( *, * ) '  RGG for real generalized problem.'
  write ( *, * ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, * ) '    A*X = LAMBDA * B * X'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a, '  The matrix A:' )

  call rmat_print ( nm, n, n, b, '  The matrix B:' )

  matz = 1

  call rgg ( nm, n, a, b, alfr, alfi, beta, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if

  alfr(1:n) = alfr(1:n) / beta(1:n)
  alfi(1:n) = alfi(1:n) / beta(1:n)

  call rvec2_print ( n, alfr, alfi, '  Real and imaginary parts of eigenvalues:' )
 
  if ( matz /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  The eigenvectors are:'
    do i = 1, n
      write ( *, * ) ' '
      write ( *, * ) '  Eigenvector ',I
      write ( *, * ) ' '
      do j = 1, n
        write ( *, '(g14.6)' ) x(i,j)
      end do
    end do
  end if
!
!  Check.
!  First, restore the original values of A and B.
!
  if ( matz /= 0 ) then
 
    a(1:n,1:n) = acopy(1:n,1:n)
    b(1:n,1:n) = bcopy(1:n,1:n)
 
    do k = 1, n
      write ( *, * ) ' '
      write ( *, * ) 'Residuals (A*x-(Alfr+Alfi*I)*B*x) for eigenvalue ',k
      write ( *, * ) ' '
 
      if ( alfi(k)==0.0 ) then
 
        do i = 1, n
 
          sum = 0.0
          do j = 1, n
            sum = sum + a(i,j) * x(j,k)
          end do
 
          do j = 1, n
            sum = sum - alfr(k) * b(i,j) * x(j,k)
          end do
 
          write ( *, '(g14.6)' ) sum
        end do
 
      else if ( alfi(k)>0.0 ) then
 
        do i = 1, n
 
          sum1 = 0.0
          sum2 = 0.0
          do j = 1, n
            sum1 = sum1+a(i,j)*x(j,k)
            sum2 = sum2+a(i,j)*x(j,k+1)
          end do
 
          do j = 1, n
            sum1 = sum1-alfr(k)*b(i,j)*x(j,k)+alfi(k)*b(i,j)*x(j,k+1)
            sum2 = sum2-alfi(k)*b(i,j)*x(j,k)-alfr(k)*b(i,j)*x(j,k+1)
          end do
 
          write ( *, '(2g14.6)' ) sum1, sum2
        end do
 
      else if ( alfi(k)<0.0 ) then
 
        do i = 1, n
 
          sum1 = 0.0
          sum2 = 0.0
          do j = 1, n
            sum1 = sum1+a(i,j)*x(j,k-1)
            sum2 = sum2-a(i,j)*x(j,k)
          end do
 
          do j = 1, n
            sum1 = sum1-alfr(k)*b(i,j)*x(j,k-1)-alfi(k)*b(i,j)*x(j,k)
            sum2 = sum2-alfi(k)*b(i,j)*x(j,k-1)+alfr(k)*b(i,j)*x(j,k)
          end do
 
          write ( *, '(2g14.6)' ) sum1, sum2
        end do
 
      end if
 
    end do
 
  end if
 
  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 tests RS.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
!
  real a(n,n)
  real a2(n,n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real r(n,n)
  real w(n)
  real x(n,n)
!
!  Set the values in the matrix.
!
  a(1,1) = 5.0
  a(1,2) = 4.0
  a(1,3) = 1.0
  a(1,4) = 1.0
 
  a(2,1) = 4.0
  a(2,2) = 5.0
  a(2,3) = 1.0
  a(2,4) = 1.0
 
  a(3,1) = 1.0
  a(3,2) = 1.0
  a(3,3) = 4.0
  a(3,4) = 2.0
 
  a(4,1) = 1.0
  a(4,2) = 1.0
  a(4,3) = 2.0
  a(4,4) = 4.0
!
!  Save a copy of the matrix.
!
  a2(1:n,1:n) = a(1:n,1:n)

  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  RS computes the eigenvalues and eigenvectors'
  write ( *, * ) '  of a real symmetric matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a, '  The matrix A:' )

  matz = 1

  call rs ( nm, n, a, w, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if

  call rvec_print ( n, w, '  The eigenvalues Lambda:' )
 
  if ( matz /= 0 ) then
 
    call rmat_print ( nm, n, n, x, '  The eigenvector matrix:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, n, r, '  The residual (A-Lambda*I)*X:' )
 
  end if
 
  return
end
subroutine test07
!
!*******************************************************************************
!
!! TEST07 tests RSB.
!
  integer, parameter :: n = 5
  integer, parameter :: nm = n
  integer, parameter :: mb = 2
!
  real a(n,mb)
  real a2(n,n)
  integer i
  integer ierr
  integer j
  integer matz
  real r(n,n)
  real w(n)
  real x(n,n)
!
  a(1:nm,1:mb) = 0.0

  a(1:n,mb) = 2.0

  a(2:n,1) = -1.0

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a2(i,j) = 2.0
      else if ( abs ( i - j ) == 1 ) then
        a2(i,j) = - 1.0
      else
        a2(i,j) = 0.0
      end if
    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'TEST07'
  write ( *, * ) '  RSB computes the eigenvalues and eigenvectors'
  write ( *, * ) '  of a real symmetric band matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a2, '  The matrix A:' )

  matz = 1

  call rsb ( nm, n, mb, a, w, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ', ierr
  end if
 
  call rvec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call rmat_print ( nm, n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, n, r, '  The residual (A-Lambda*I)*X:' )
 
  end if
 
  return
end
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests RSG.
!
!
!  RGG is for a real generalized eigenvalue problem of the form
!
!    A*x = lambda*B*x
!
!  with A symmetric and B positive definite symmetric.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
!
  real a(n,n)
  real a2(n,n)
  real b(n,n)
  real b2(n,n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real sum
  real sum1
  real sum2
  real w(n)
  real x(n,n)
!
  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0
      else
        b(i,j) = 0.0
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  write ( *, * ) ' '
  write ( *, * ) 'TEST08:'
  write ( *, * ) '  RSG for real symmetric generalized problem.'
  write ( *, * ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, * ) '    A*X = LAMBDA * B * X'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a, '  The matrix A:' )

  call rmat_print ( nm, n, n, b, '  The matrix B:' )

  matz = 1

  call rsg ( nm, n, a, b, w, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if

  call rvec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call rmat_print ( nm, n, n, x, '  The eigenvector matrix X:' )

    a(1:n,1:n) = a2(1:n,1:n)
    b(1:n,1:n) = b2(1:n,1:n)
 
    do k = 1, n

      write ( *, * ) ' '
      write ( *, * ) 'Residuals (A*x-(w*I)*B*x) for eigenvalue ',k
      write ( *, * ) ' '
  
        do i = 1, n
 
          sum = 0.0
          do j = 1, n
            sum = sum+a(i,j)*x(j,k)
          end do
 
          do j = 1, n
            sum = sum-w(k)*b(i,j)*x(j,k)
          end do
 
          write ( *, '(g14.6)' ) sum
        end do
 
    end do
 
  end if
 
  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests RSGAB.
!
!
!  RGGAB is for a real generalized eigenvalue problem of the form
!
!    A*B*x = lambda*x
!
!  with A symmetric and B positive definite symmetric.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
!
  real a(n,n)
  real a2(n,n)
  real b(n,n)
  real b2(n,n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real r(n,n)
  real w(n)
  real x(n,n)
!
  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0
      else
        b(i,j) = 0.0
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  write ( *, * ) ' '
  write ( *, * ) 'TEST09:'
  write ( *, * ) '  RSGAB for real symmetric generalized problem.'
  write ( *, * ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, * ) '    A*B*X = LAMBDA * X'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a, '  The matrix A:' )

  call rmat_print ( nm, n, n, b, '  The matrix B:' )

  matz = 1

  call rsgab ( nm, n, a, b, w, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if

  call rvec_print ( n, w, '  EThe eigenvalues Lambda:' )
 
  if ( matz /= 0 ) then

    call rmat_print ( nm, n, n, x, '  The eigenvector matrix X:' ) 
 
    r = matmul ( b2, x )

    r = matmul ( a2, r )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, n, r, '  The residual matrix (A*B-Lambda*I)*X:' )
 
  end if
 
  return
end
subroutine test10
!
!*******************************************************************************
!
!! TEST10 tests RSGBA.
!
!
!  RGGBA is for a real generalized eigenvalue problem of the form
!
!    B*A*x = lambda*x
!
!  with A symmetric and B positive definite symmetric.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
!
  real a(n,n)
  real a2(n,n)
  real b(n,n)
  real b2(n,n)
  integer i
  integer ierr
  integer j
  integer matz
  real r(n,n)
  real w(n)
  real x(n,n)
!
  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0
      else
        b(i,j) = 0.0
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  write ( *, * ) ' '
  write ( *, * ) 'TEST10:'
  write ( *, * ) '  RSGBA for real symmetric generalized problem.'
  write ( *, * ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, * ) '    B*A*X = LAMBDA * X'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a, '  The matrix A:' )

  call rmat_print ( nm, n, n, b, '  The matrix B:' )

  matz = 1

  call rsgba ( nm, n, a, b, w, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if

  call rvec_print ( n, w, '  The eigenvalues Lambda:' )
 
  if ( matz /= 0 ) then

    call rmat_print ( nm, n, n, x, '  The eigenvector matrix X:' )
 
    r = matmul ( a2, x )

    r = matmul ( b2, r )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, n, r, '  The residual matrix (B*A-Lambda*I)*X:' )
 
  end if
 
  return
end
subroutine test11
!
!*******************************************************************************
!
!! TEST11 tests RSM.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
  integer, parameter :: m = n
!
  real a(n,n)
  real a2(n,n)
  integer ierr
  integer j
  integer matz
  real r(n,m)
  real w(n)
  real x(n,m)
!
  a(1,1) = 5.0
  a(1,2) = 4.0
  a(1,3) = 1.0
  a(1,4) = 1.0
 
  a(2,1) = 4.0
  a(2,2) = 5.0
  a(2,3) = 1.0
  a(2,4) = 1.0
 
  a(3,1) = 1.0
  a(3,2) = 1.0
  a(3,3) = 4.0
  a(3,4) = 2.0
 
  a(4,1) = 1.0
  a(4,2) = 1.0
  a(4,3) = 2.0
  a(4,4) = 4.0

  a2(1:n,1:n) = a(1:n,1:n)

  write ( *, * ) ' '
  write ( *, * ) 'TEST11'
  write ( *, * ) '  RSM computes some eigenvalues and eigenvectors'
  write ( *, * ) '  of a real symmetric matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n
  write ( *, * ) '  Number of eigenvectors desired = ', m

  call rmat_print ( nm, n, n, a, '  The matrix A:' )

  call rsm ( nm, n, a, w, m, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ',ierr
  end if

  call rvec_print ( n, w, '  The eigenvalues Lambda:' )
 
  if ( m > 0 ) then
 
    call rmat_print ( nm, n, m, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, m
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, m, r, '  The residual (A-Lambda*I)*X:' )
 
  end if
 
  return
end
subroutine test12
!
!*******************************************************************************
!
!! TEST12 tests RSP.
!
!
!  RSP is for the eigenvalues of a real symmetric packed matrix.
!
!  A is symmetric.  Because of this, we know that the eigenvalues
!  of A must be real (rather than complex) numbers.
!
!
!  The entries of A are
!
!  (5 4 1 1)
!  (4 5 1 1)
!  (1 1 4 2)
!  (1 1 2 4)
!
!  The eigenvalues of A are (10, 5, 2, 1)
!
!  One set of eigenvectors of A is:
!
!  ( 2 -1  0 -1)
!  ( 2 -1  0  1)
!  ( 1  2 -1  0)
!  ( 1  2  1  0)
!
!  However, this set is not orthonormal, and EISPACK will compute
!  a different set of values.
!
!  Note that the I-th eigenvector corresponding to the I-th eigenvalue
!  consists of the I-th column of the above matrix of eigenvectors.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
  integer, parameter :: nv = ( n * ( n + 1 ) ) / 2
!
  real a(nv)
  real a2(n,n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real r(n,n)
  real w(n)
  real x(n,n)
!
!  Set the values in the matrix.
!
  a(1) = 5.0
 
  a(2) = 4.0
  a(3) = 5.0
 
  a(4) = 1.0
  a(5) = 1.0
  a(6) = 4.0
 
  a(7) = 1.0
  a(8) = 1.0
  a(9) = 2.0
  a(10) = 4.0

  a2(1,1) = 5.0
  a2(1,2) = 4.0
  a2(1,3) = 1.0
  a2(1,4) = 1.0
 
  a2(2,1) = 4.0
  a2(2,2) = 5.0
  a2(2,3) = 1.0
  a2(2,4) = 1.0
 
  a2(3,1) = 1.0
  a2(3,2) = 1.0
  a2(3,3) = 4.0
  a2(3,4) = 2.0
 
  a2(4,1) = 1.0
  a2(4,2) = 1.0
  a2(4,3) = 2.0
  a2(4,4) = 4.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST12'
  write ( *, * ) '  RSP computes the eigenvalues and eigenvectors'
  write ( *, * ) '  of a real symmetric packed matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a2, '  The matrix A:' )

  matz = 1

  call rsp ( nm, n, nv, a, w, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag was IERR = ',ierr
  end if

  call rvec_print ( n, w, '  The eigenvalues Lambda:' )
 
  if ( matz /= 0 ) then
 
    call rmat_print ( nm, n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if
 
  return
end
subroutine test13
!
!*******************************************************************************
!
!! TEST13 tests RSPP.
!
!
!  RSPP is for some eigenvalues of a real symmetric packed matrix.
!
!  A is symmetric.  Because of this, we know that the eigenvalues
!  of A must be real (rather than complex) numbers.
!
!
!  The entries of A are
!
!  (5 4 1 1)
!  (4 5 1 1)
!  (1 1 4 2)
!  (1 1 2 4)
!
!  The eigenvalues of A are (10, 5, 2, 1)
!
!  One set of eigenvectors of A is:
!
!  ( 2 -1  0 -1)
!  ( 2 -1  0  1)
!  ( 1  2 -1  0)
!  ( 1  2  1  0)
!
!  However, this set is not orthonormal, and EISPACK will compute
!  a different set of values.
!
!  Note that the I-th eigenvector corresponding to the I-th eigenvalue
!  consists of the I-th column of the above matrix of eigenvectors.
!
  integer, parameter :: n = 4
  integer, parameter :: nm = n
  integer, parameter :: m = n
  integer, parameter :: nv = ( n * ( n + 1 ) ) / 2
!
  real a(nv)
  real a2(n,n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real r(n,m)
  logical type
  real w(m)
  real x(n,m)
!
!  Set the values in the matrix.
!
  a(1) = 5.0
 
  a(2) = 4.0
  a(3) = 5.0
 
  a(4) = 1.0
  a(5) = 1.0
  a(6) = 4.0
 
  a(7) = 1.0
  a(8) = 1.0
  a(9) = 2.0
  a(10) = 4.0

  a2(1,1) = 5.0
  a2(1,2) = 4.0
  a2(1,3) = 1.0
  a2(1,4) = 1.0
 
  a2(2,1) = 4.0
  a2(2,2) = 5.0
  a2(2,3) = 1.0
  a2(2,4) = 1.0
 
  a2(3,1) = 1.0
  a2(3,2) = 1.0
  a2(3,3) = 4.0
  a2(3,4) = 2.0
 
  a2(4,1) = 1.0
  a2(4,2) = 1.0
  a2(4,3) = 2.0
  a2(4,4) = 4.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST13'
  write ( *, * ) '  RSPP finds some eigenvalues and eigenvectors of'
  write ( *, * ) '  a real symmetric packed matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a2, '  The matrix A:' )
!
!  Set MATZ = 0 for no eigenvectors, 1 for eigenvectors.
!
  matz = 1
!
!  TYPE = TRUE to find smallest eigenvalues, FALSE for largest.
!
  type = .true.

  call rspp ( nm, n, nv, a, w, matz, x, ierr, m, type )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag was IERR = ',ierr
  end if

  call rvec_print ( m, w, '  The eigenvalues Lambda:' )
 
  if ( matz /= 0 ) then 

    call rmat_print ( nm, n, m, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, m
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, m, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test14
!
!*******************************************************************************
!
!! TEST14 tests RST.
!
  integer, parameter :: n = 5
  integer, parameter :: nm = n
!
  real a(n,n)
  real e(n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real r(n,n)
  real w(n)
  real x(n,n)
!
!  Here is where the matrix is defined.
!
  w(1:n) = 2.0

  e(1) = 0.0
  e(2:n) = -1.0
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0
      else
        a(i,j) = 0.0
      end if
    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'TEST14'
  write ( *, * ) '  RST computes the eigenvalues and eigenvectors'
  write ( *, * ) '  of a real symmetric tridiagonal matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a, '  The matrix A:' )
 
  matz = 1

  call rst ( nm, n, w, e, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ', ierr
  end if
 
  call rvec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call rmat_print ( nm, n, n, x, '  The eigenvector matrix X:' )
 
    r = matmul ( a, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test145
!
!*******************************************************************************
!
!! TEST145 tests RT.
!
  integer, parameter :: n = 5
  integer, parameter :: nm = n
!
  real a(n,3)
  real a2(n,n)
  real e(n)
  integer i
  integer ierr
  integer j
  integer k
  integer matz
  real r(n,n)
  real w(n)
  real x(n,n)
!
!  Here is where the matrix is defined.
!
  a(2:n,  1) = - 1.0
  a(1:n,  2) =   2.0
  a(1:n-1,3) = - 1.0
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a2(i,j) = 2.0
      else if ( abs ( i - j ) == 1 ) then
        a2(i,j) = - 1.0
      else
        a2(i,j) = 0.0
      end if
    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'TEST145'
  write ( *, * ) '  RT computes the eigenvalues and eigenvectors'
  write ( *, * ) '  of a real sign-symmetric tridiagonal matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, n, n, a2, '  The matrix A:' )
 
  matz = 1

  call rt ( nm, n, a, w, matz, x, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ', ierr
  end if
 
  call rvec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz /= 0 ) then

    call rmat_print ( nm, n, n, x, '  The eigenvector matrix X:' )
 
    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call rmat_print ( nm, n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine test15
!
!*******************************************************************************
!
!! TEST15 tests SVD.
!
!
!  In our special example, the matrix is square and symmetric.
!
  integer, parameter :: n = 4
  integer, parameter :: m = n
  integer, parameter :: nm = n
!
  real a(m,n)
  integer i
  integer ierr
  integer j
  integer k
  logical matu
  logical matv
  real r(m,n)
  real u(m,n)
  real v(n,n)
  real w(n)
!
!  Set the values of the matrix.
!
  a(1,1) = 0.9900
  a(1,2) = 0.0020
  a(1,3) = 0.0060
  a(1,4) = 0.0020

  a(2,1) = 0.0020
  a(2,2) = 0.9900
  a(2,3) = 0.0020
  a(2,4) = 0.0060

  a(3,1) = 0.0060
  a(3,2) = 0.0020
  a(3,3) = 0.9900
  a(3,4) = 0.0020

  a(4,1) = 0.0020
  a(4,2) = 0.0060
  a(4,3) = 0.0020
  a(4,4) = 0.9900
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST15'
  write ( *, * ) '  SVD computes the singular value decomposition'
  write ( *, * ) '  of a real general matrix.'
  write ( *, * ) ' '
  write ( *, * ) '  Matrix order = ', n

  call rmat_print ( nm, m, n, a, '  The matrix A:' )

  matu = .TRUE.
  matv = .TRUE.

  call svd ( nm, m, n, a, w, matu, u, matv, v, ierr )
 
  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Warning!'
    write ( *, * ) '  The error return flag IERR = ', ierr
  end if

  call rvec_print ( n, w, '  The singular values S' )
 
  call rmat_print ( nm, m, n, u, '  The U matrix:' )

  call rmat_print ( nm, n, n, v, '  The V matrix:' )

  do j = 1, n
    v(1:n,j) = w(j) * v(1:n,j)
  end do

  r = matmul ( u, transpose ( v ) )

  call rmat_print ( nm, m, n, r, '  The product U * S * Transpose(V):' )

  return
end
