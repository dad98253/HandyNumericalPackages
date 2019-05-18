!  eispack.f90  21 June 2000
!
subroutine bakvec ( nm, n, t, e, m, z, ierr )
!
!*******************************************************************************
!
!! BAKVEC determines eigenvectors by reversing the FIGI transformation.
!
!
!  Discussion:
!
!    This subroutine forms the eigenvectors of a nonsymmetric tridiagonal 
!    matrix by back transforming those of the corresponding symmetric 
!    matrix determined by FIGI.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        t contains the nonsymmetric matrix.  its subdiagonal is
!          stored in the last n-1 positions of the first column,
!          its diagonal in the n positions of the second column,
!          and its superdiagonal in the first n-1 positions of
!          the third column.  t(1,1) and t(n,3) are arbitrary.
!
!        e contains the subdiagonal elements of the symmetric
!          matrix in its last n-1 positions.  e(1) is arbitrary.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        t is unaltered.
!
!        e is destroyed.
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!        ierr is set to
!          zero       for normal return,
!          2*n+i      if e(i) is zero with t(i,1) or t(i-1,3) non-zero.
!                     in this case, the symmetric matrix is not similar
!                     to the original matrix, and the eigenvectors
!                     cannot be found by this program.
!
  integer m
  integer n
  integer nm
!
  real e(n)
  integer i
  integer ierr
  integer j
  real t(nm,3)
  real z(nm,m)
!
  ierr = 0

  if ( m == 0 ) then
    return
  end if

  e(1) = 1.0
  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    if ( e(i) == 0.0 ) then
      if (t(i,1) /= 0.0 .or. t(i-1,3) /= 0.0 ) then
        ierr = 2 * n + i
        return
      end if
      e(i) = 1.0
    else
      e(i) = e(i-1) * e(i) / t(i-1,3)
    end if
  end do

  do j = 1, m
    z(2:n,j) = z(2:n,j) * e(2:n)
  end do

  return
end
subroutine balanc ( nm, n, a, low, igh, scale )
!
!*******************************************************************************
!
!! BALANC balances a real matrix before eigenvalue calculations.
!
!
!  Discussion:
!
!    This subroutine balances a real matrix and isolates eigenvalues 
!    whenever possible.
!
!    Suppose that the principal submatrix in rows LOW through IGH
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!
!      SCALE(J) = P(J),    J = 1,...,LOW-1,
!               = D(J,J),  J = LOW,...,IGH,
!               = P(J)     J = IGH+1,...,N.
!
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is returned for LOW if IGH is zero formally.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the input matrix to be balanced.
!
!     on output
!
!        a contains the balanced matrix.
!
!        low and igh are two integers such that a(i,j)
!          is equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
  integer nm
  integer n
!
  real a(nm,n)
  real b2
  real c
  real f
  real g
  integer i
  integer iexc
  integer igh
  integer j
  integer k
  integer l
  integer low
  integer m
  logical noconv
  real r
  real radix
  real s
  real scale(n)
!
  radix = 16.0

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

   20 continue

  scale(m) = j

  if ( j /= m ) then

    do i = 1, l
      call r_swap ( a(i,j), a(i,m) )
    end do

    do i = k, n
      call r_swap ( a(j,i), a(m,i) )
    end do

  end if

50    continue

  if(iexc==2)go to 130
!
!  Search for rows isolating an eigenvalue and push them down.
!
   80 if (l == 1) go to 280
  l = l - 1

  100 continue

  do j = l, 1, -1

     do i = 1, l
       if (i /= j) then
         if (a(j,i) /= 0.0) go to 120
       end if
     end do

     m = l
     iexc = 1
     go to 20
120      continue
  end do

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
  130 k = k + 1

  140 continue

  do j = k, l

    do i = k, l
      if (i /= j) then
        if (a(i,j) /= 0.0) go to 170
      end if
    end do

    m = k
    iexc = 2
    go to 20

  170   continue

  end do
!
!  Balance the submatrix in rows k to l.
!
  scale(k:l) = 1.0
!
!  Iterative loop for norm reduction.
!
  190 noconv = .false.

  do i = k, l

     c = 0.0
     r = 0.0

     do j = k, l
       if (j /= i) then
         c = c + abs(a(j,i))
         r = r + abs(a(i,j))
       end if
     end do
!
!  Guard against zero c or r due to underflow.
!
     if ( c == 0.0 .or. r == 0.0 ) go to 270

     g = r / radix
     f = 1.0
     s = c + r

     do while ( c < g ) 
       f = f * radix
       c = c * b2
     end do

     g = r * radix

     do while ( c >= g ) 
       f = f / radix
       c = c / b2
     end do
!
!  Balance.
!
     if ( (c + r) / f < 0.95 * s ) then

       g = 1.0 / f
       scale(i) = scale(i) * f
       noconv = .true.

       a(i,k:n) = a(i,k:n) * g
       a(1:l,i) = a(1:l,i) * f

      end if

  270 continue

  end do

  if (noconv) go to 190

  280 continue

  low = k
  igh = l
  return
end
subroutine balbak ( nm, n, low, igh, scale, m, z )
!
!*******************************************************************************
!
!! BALBAK determines eigenvectors by undoing the BALANC transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  balanc.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  balanc.
!
!        scale contains information determining the permutations
!          and scaling factors used by  balanc.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
  integer n
  integer nm
!
  integer i
  integer igh
  integer ii
  integer j
  integer k
  integer low
  integer m
  real s
  real scale(n)
  real z(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  if ( igh /= low ) then

    do i = low, igh
      z(i,1:m) = z(i,1:m) * scale(i)
    end do

  end if

   do ii = 1, n

     i = ii

     if ( i < low .or. i > igh ) then

       if ( i < low ) then
         i = low - ii
       end if

       k = int(scale(i))

       if ( k /= i ) then

         do j = 1, m
           call r_swap ( z(i,j), z(k,j) )
         end do

        end if

      end if

  end do

  return
end
subroutine bandr ( nm, n, mb, a, d, e, e2, matz, z )
!
!*******************************************************************************
!
!! BANDR reduces a symmetric band matrix to tridiagonal form.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure bandrd,
!     num. math. 12, 231-241(1968) by schwarz.
!     handbook for auto. comp., vol.ii-linear algebra, 273-283(1971).
!
!     this subroutine reduces a real symmetric band matrix
!     to a symmetric tridiagonal matrix using and optionally
!     accumulating orthogonal similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        mb is the (half) band width of the matrix, defined as the
!          number of adjacent diagonals, including the principal
!          diagonal, required to specify the non-zero portion of the
!          lower triangle of the matrix.
!
!        a contains the lower triangle of the symmetric band input
!          matrix stored as an n by mb array.  its lowest subdiagonal
!          is stored in the last n+1-mb positions of the first column,
!          its next subdiagonal in the last n+2-mb positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the n positions of the last column.
!          contents of storages not part of the matrix are arbitrary.
!
!        matz should be set to .true. if the transformation matrix is
!          to be accumulated, and to .false. otherwise.
!
!     on output
!
!        a has been destroyed, except for its last two columns which
!          contain a copy of the tridiagonal matrix.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        z contains the orthogonal transformation matrix produced in
!          the reduction if matz has been set to .true.  otherwise, z
!          is not referenced.
!
  integer mb
  integer n
  integer nm
!
  real a(nm,mb)
  real b1
  real b2
  real c2
  real d(n)
  real dmin
  real dminrt
  real e(n)
  real e2(n)
  real f1
  real f2
  real g
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer k
  integer kr
  integer l
  integer m1
  logical matz
  integer maxl
  integer maxr
  integer mr
  integer r
  integer r1
  real s2
  real u
  integer ugl
  real z(nm,n)
!
  dmin = 2.0**(-64)
  dminrt = 2.0**(-32)
!
!  Initialize the diagonal scaling matrix.
!
  d(1:n) = 1.0

  if ( matz ) then
    call rmat_ident ( nm, n, z )
  end if

  m1 = mb - 1
  if (m1 - 1) 900, 800, 70
   70 continue

  do k = 1, n - 2

     maxr = min ( m1, n-k )

     do r1 = 2, maxr

        r = maxr + 2 - r1
        kr = k + r
        mr = mb - r
        g = a(kr,mr)
        a(kr-1,1) = a(kr-1,mr+1)
        ugl = k

        do j = kr, n, m1

           j1 = j - 1
           j2 = j1 - 1

           if (g == 0.0) go to 600

           b1 = a(j1,1) / g
           b2 = b1 * d(j1) / d(j)
           s2 = 1.0 / (1.0 + b1 * b2)

           if (s2 >= 0.5 ) go to 450

           b1 = g / a(j1,1)
           b2 = b1 * d(j) / d(j1)
           c2 = 1.0 - s2
           d(j1) = c2 * d(j1)
           d(j) = c2 * d(j)
           f1 = 2.0 * a(j,m1)
           f2 = b1 * a(j1,mb)
           a(j,m1) = -b2 * (b1 * a(j,m1) - a(j,mb)) - f2 + a(j,m1)
           a(j1,mb) = b2 * (b2 * a(j,mb) + f1) + a(j1,mb)
           a(j,mb) = b1 * (f2 - f1) + a(j,mb)

           do l = ugl, j2
              i2 = mb - j + l
              u = a(j1,i2+1) + b2 * a(j,i2)
              a(j,i2) = -b1 * a(j1,i2+1) + a(j,i2)
              a(j1,i2+1) = u
           end do

           ugl = j
           a(j1,1) = a(j1,1) + b2 * g

           if ( j /= n ) then

             maxl = min(m1,n-j1)

             do l = 2, maxl
               i1 = j1 + l
               i2 = mb - l
               u = a(i1,i2) + b2 * a(i1,i2+1)
               a(i1,i2+1) = -b1 * a(i1,i2) + a(i1,i2+1)
               a(i1,i2) = u
             end do

             i1 = j + m1

             if ( i1 <= n ) then
               g = b2 * a(i1,1)
             end if

           end if

           if ( matz ) then

             do l = 1, n
               u = z(l,j1) + b2 * z(l,j)
               z(l,j) = -b1 * z(l,j1) + z(l,j)
               z(l,j1) = u
             end do

           end if

           go to 500

  450      continue

           u = d(j1)
           d(j1) = s2 * d(j)
           d(j) = s2 * u
           f1 = 2.0 * a(j,m1)
           f2 = b1 * a(j,mb)
           u = b1 * (f2 - f1) + a(j1,mb)
           a(j,m1) = b2 * (b1 * a(j,m1) - a(j1,mb)) + f2 - a(j,m1)
           a(j1,mb) = b2 * (b2 * a(j1,mb) + f1) + a(j,mb)
           a(j,mb) = u

           do l = ugl, j2
              i2 = mb - j + l
              u = b2 * a(j1,i2+1) + a(j,i2)
              a(j,i2) = -a(j1,i2+1) + b1 * a(j,i2)
              a(j1,i2+1) = u
           end do

           ugl = j
           a(j1,1) = b2 * a(j1,1) + g

           if ( j /= n ) then

             maxl = min(m1,n-j1)

             do l = 2, maxl
               i1 = j1 + l
               i2 = mb - l
               u = b2 * a(i1,i2) + a(i1,i2+1)
               a(i1,i2+1) = -a(i1,i2) + b1 * a(i1,i2+1)
               a(i1,i2) = u
             end do

             i1 = j + m1

             if ( i1 <= n ) then
               g = a(i1,1)
               a(i1,1) = b1 * a(i1,1)
             end if

           end if

           if ( matz ) then

             do l = 1, n
               u = b2 * z(l,j1) + z(l,j)
               z(l,j) = -z(l,j1) + b1 * z(l,j)
               z(l,j1) = u
             end do

           end if

  500       continue

           end do

  600    continue

     end do

     if (mod(k,64) /= 0) go to 700
!
!  Rescale to avoid underflow or overflow
!
     do j = k, n

       if ( d(j) < dmin ) then

        maxl = max(1,mb+1-j)

        a(j,maxl:m1) = dminrt * a(j,maxl:m1)

        if ( j /= n ) then

        maxl = min(m1,n-j)

        do l = 1, maxl
          i1 = j + l
          i2 = mb - l
          a(i1,i2) = dminrt * a(i1,i2)
        end do

        end if

        if ( matz ) then
          z(1:n,j) = dminrt * z(1:n,j)
        end if

        a(j,mb) = dmin * a(j,mb)
        d(j) = d(j) / dmin

        end if

      end do

  700 continue

    end do
!
!   Form square root of scaling matrix.
!
  800 continue

  e(2:n) = sqrt(d(2:n))

  if ( matz ) then

    do k = 2, n
      z(1:n,k) = z(1:n,k) * e(k)
    end do

  end if

  u = 1.0

  do j = 2, n
    a(j,m1) = u * e(j) * a(j,m1)
    u = e(j)
    e2(j) = a(j,m1) ** 2
    a(j,mb) = d(j) * a(j,mb)
    d(j) = a(j,mb)
    e(j) = a(j,m1)
  end do

  d(1) = a(1,mb)
  e(1) = 0.0
  e2(1) = 0.0

  return

  900 continue

   d(1:n) = a(1:n,mb)
   e(1:n) = 0.0
   e2(1:n) = 0.0

  return
end
subroutine bandv ( nm, n, mbw, a, e21, m, w, z, ierr, nv, rv )
!
!*******************************************************************************
!
!! BANDV finds eigenvectors from eigenvalues, for a real symmetric band matrix.
!
!
!  Discussion:
!
!     this subroutine finds those eigenvectors of a real symmetric
!     band matrix corresponding to specified eigenvalues, using inverse
!     iteration.  the subroutine may also be used to solve systems
!     of linear equations with a symmetric or non-symmetric band
!     coefficient matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        mbw is the number of columns of the array a used to store the
!          band matrix.  if the matrix is symmetric, mbw is its (half)
!          band width, denoted mb and defined as the number of adjacent
!          diagonals, including the principal diagonal, required to
!          specify the non-zero portion of the lower triangle of the
!          matrix.  if the subroutine is being used to solve systems
!          of linear equations and the coefficient matrix is not
!          symmetric, it must however have the same number of adjacent
!          diagonals above the main diagonal as below, and in this
!          case, mbw=2*mb-1.
!
!        a contains the lower triangle of the symmetric band input
!          matrix stored as an n by mb array.  its lowest subdiagonal
!          is stored in the last n+1-mb positions of the first column,
!          its next subdiagonal in the last n+2-mb positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the n positions of column mb.
!          if the subroutine is being used to solve systems of linear
!          equations and the coefficient matrix is not symmetric, a is
!          n by 2*mb-1 instead with lower triangle as above and with
!          its first superdiagonal stored in the first n-1 positions of
!          column mb+1, its second superdiagonal in the first n-2
!          positions of column mb+2, further superdiagonals similarly,
!          and finally its highest superdiagonal in the first n+1-mb
!          positions of the last column.
!          contents of storages not part of the matrix are arbitrary.
!
!        e21 specifies the ordering of the eigenvalues and contains
!            0.0 if the eigenvalues are in ascending order, or
!            2.0 if the eigenvalues are in descending order.
!          if the subroutine is being used to solve systems of linear
!          equations, e21 should be set to 1.0 if the coefficient
!          matrix is symmetric and to -1.0 if not.
!
!        m is the number of specified eigenvalues or the number of
!          systems of linear equations.
!
!        w contains the m eigenvalues in ascending or descending order.
!          if the subroutine is being used to solve systems of linear
!          equations (a-w(r)*i)*x(r)=b(r), where i is the identity
!          matrix, w(r) should be set accordingly, for r=1,2,...,m.
!
!        z contains the constant matrix columns (b(r),r=1,2,...,m), if
!          the subroutine is used to solve systems of linear equations.
!
!        nv must be set to the dimension of the array parameter rv
!          as declared in the calling program dimension statement.
!
!     on output
!
!        a and w are unaltered.
!
!        z contains the associated set of orthogonal eigenvectors.
!          any vector which fails to converge is set to zero.  if the
!          routine is used to solve systems of linear equations,
!          z contains the solution matrix columns (x(r),r=1,2,...,m).
!
!        ierr is set to
!          zero       for normal return,
!          -r         if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge, or if the r-th
!                     system of linear equations is nearly singular.
!
!        rv is a temporary storage arrays.  note that rv is
!          of dimension at least n*(2*mb-1).  if the subroutine
!          is being used to solve systems of linear equations, the
!          determinant (up to sign) of a-w(m)*i is available, upon
!          return, as the product of the first n elements of rv.
!
  integer mbw
  integer n
  integer nm
  integer nv
!
  real a(nm,mbw)
  real e21
  real eps2
  real eps3
  real eps4
  integer group
  integer i
  integer ierr
  integer ii
  integer ij
  integer ij1
  integer its
  integer j
  integer jj
  integer k
  integer kj
  integer kj1
  integer m
  integer m1
  integer m21
  integer maxj
  integer maxk
  integer mb
  real norm
  real order
  real pythag
  integer r
  real rv(nv)
  real rv6(n)
  real u
  real uk
  real v
  real w(m)
  real x0
  real x1
  real xu
  real z(nm,m)
!
  ierr = 0

  if ( m == 0 ) then
    return
  end if

  x0 = 0.0
  mb = mbw
  if (e21 < 0.0) mb = (mbw + 1) / 2
  m1 = mb - 1
  m21 = m1 + mb
  order = 1.0 - abs ( e21 )
!
!  Find vectors by inverse iteration.
!
  do r = 1, m

     its = 1
     x1 = w(r)
     if (r /= 1) go to 100
!
!  Compute norm of matrix
!
     norm = 0.0

     do j = 1, mb

        jj = mb + 1 - j
        kj = jj + m1
        ij = 1
        v = 0.0

        do i = jj, n

          v = v + abs(a(i,j))

          if ( e21 < 0.0 ) then
            v = v + abs ( a(ij,kj) )
            ij = ij + 1
          end if

        end do

        norm = max(norm,v)

     end do

     if (e21 < 0.0) norm = 0.5 * norm
!
!  eps2 is the criterion for grouping,
!  eps3 replaces zero pivots and equal roots are modified by eps3,
!  eps4 is taken very small to avoid overflow
!
     if (norm == 0.0) norm = 1.0
     eps2 = 1.0e-3 * norm * abs(order)
     eps3 = abs ( norm ) * epsilon ( 1.0 )
     uk = n
     uk = sqrt(uk)
     eps4 = uk * eps3
   80    group = 0
     go to 120
!
!  Look for close or coincident roots.
!
  100    continue

     if (abs(x1-x0) >= eps2) go to 80
     group = group + 1

     if ( order * (x1 - x0) <= 0.0 ) then
       x1 = x0 + order * eps3
     end if
!
!  Expand matrix, subtract eigenvalue, and initialize vector.
!
  120   continue

     do i = 1, n

        ij = i + min(0,i-m1) * n
        kj = ij + mb * n
        ij1 = kj + m1 * n

        if (m1 == 0) go to 180

        do j = 1, m1

          if ( ij <= m1 ) then
            if ( ij <= 0 ) then
              rv(ij1) = 0.0
              ij1 = ij1 + n
            end if
          else
            rv(ij) = a(i,j)
          end if

          ij = ij + n
          ii = i + j

          if ( ii <= n ) then

            jj = mb - j

            if ( e21 < 0.0 ) then
              ii = i
              jj = mb + j
            end if

            rv(kj) = a(ii,jj)
            kj = kj + n

          end if

        end do

  180   continue

        rv(ij) = a(i,mb) - x1
        rv6(i) = eps4
        if ( order == 0.0 ) rv6(i) = z(i,r)

     end do

     if ( m1 /= 0 ) then
!
!  Elimination with interchanges.
!
     do i = 1, n

        ii = i + 1
        maxk = min(i+m1-1,n)
        maxj = min(n-i,m21-2) * n

        do k = i, maxk

           kj1 = k
           j = kj1 + n
           jj = j + maxj

           do kj = j, jj, n
              rv(kj1) = rv(kj)
              kj1 = kj
           end do

           rv(kj1) = 0.0

        end do

        if ( i /= n ) then

        u = 0.0
        maxk = min(i+m1,n)
        maxj = min(n-ii,m21-2) * n

        do j = i, maxk
          if ( abs(rv(j)) >= abs(u) ) then
            u = rv(j)
            k = j
          end if
        end do

        j = i + n
        jj = j + maxj

        if ( k /= i ) then

          kj = k

          do ij = i, jj, n
            call r_swap ( rv(ij), rv(kj) )
            kj = kj + n
          end do

          if ( order == 0.0 ) then
            call r_swap ( rv6(i), rv6(k) )
          end if

        end if

        if ( u /= 0.0 ) then

        do k = ii, maxk

           v = rv(k) / u
           kj = k

           do ij = j, jj, n
              kj = kj + n
              rv(kj) = rv(kj) - v * rv(ij)
           end do

           if ( order == 0.0 ) then
             rv6(k) = rv6(k) - v * rv6(i)
           end if

        end do

       end if

      end if

      end do

     end if
!
!  Back substitution
!
600  continue

     do ii = 1, n

        i = n + 1 - ii
        maxj = min(ii,m21)

        if ( maxj /= 1 ) then

          ij1 = i
          j = ij1 + n
          jj = j + (maxj - 2) * n

          do ij = j, jj, n
            ij1 = ij1 + 1
            rv6(i) = rv6(i) - rv(ij) * rv6(ij1)
          end do

        end if

        v = rv(i)
!
!  Error: nearly singular linear system
!
        if ( abs ( v ) < eps3) then
          if ( order == 0.0 ) then
            ierr = -r
          end if
          v = sign(eps3,v)
        end if

        rv6(i) = rv6(i) / v

     end do

     xu = 1.0

     if (order == 0.0) go to 870
!
!  Orthogonalize with respect to previous members of group.
!
     do jj = 1, group

        j = r - group - 1 + jj
        xu = 0.0

        do i = 1, n
          xu = xu + rv6(i) * z(i,j)
        end do

        rv6(1:n) = rv6(1:n) - xu * z(1:n,j)

     end do

     norm = 0.0
     do i = 1, n
       norm = norm + abs ( rv6(i) )
     end do
!
!  Choose a new starting vector
!
     if ( norm < 0.1 ) then

       if ( its < n ) then
         its = its + 1
         xu = eps4 / (uk + 1.0)
         rv6(1) = eps4
         rv6(2:n) = xu
         rv6(its) = rv6(its) - eps4 * uk
         go to 600
       else
         ierr = -r
         xu = 0.0
         go to 870
       end if

     end if
!
!   Normalize so that sum of squares is 1 and expand to full order.
!
     u = 0.0
     do i = 1, n
       u = pythag ( u, rv6(i) )
     end do

     xu = 1.0 / u

  870    continue

     z(1:n,r) = rv6(1:n) * xu

     x0 = x1

  end do

  return
end
subroutine bisect ( n, eps1, d, e, e2, lb, ub, mm, m, w, ind, ierr )
!
!*******************************************************************************
!
!! BISECT computes some eigenvalues of a real symmetric band matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the bisection technique
!     in the algol procedure tristurm by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvalues of a tridiagonal
!     symmetric matrix which lie in a specified interval,
!     using bisection.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        n is the order of the matrix.
!
!        eps1 is an absolute error tolerance for the computed
!          eigenvalues.  if the input eps1 is non-positive,
!          it is reset for each submatrix to a default value,
!          namely, minus the product of the relative machine
!          precision and the 1-norm of the submatrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        lb and ub define the interval to be searched for eigenvalues.
!          if lb is not less than ub, no eigenvalues will be found.
!
!        mm should be set to an upper bound for the number of
!          eigenvalues in the interval.  warning. if more than
!          mm eigenvalues are determined to lie in the interval,
!          an error return is made with no eigenvalues found.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        m is the number of eigenvalues determined to lie in (lb,ub).
!
!        w contains the m eigenvalues in ascending order.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w:
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc..
!
!        ierr is set to
!          zero       for normal return,
!          3*n+1      if m exceeds mm.
!
  integer mm
  integer n
!
  real d(n)
  real e(n)
  real e2(n)
  real eps1
  integer i
  integer ierr
  integer ii
  integer ind(mm)
  integer isturm
  integer j
  integer k
  integer l
  real lb
  integer m
  integer m1
  integer m2
  integer p
  integer q
  integer r
  real rv4(n)
  real rv5(n)
  integer s
  real t1
  real t2
  integer tag
  real tst1
  real tst2
  real u
  real ub
  real v
  real w(mm)
  real x0
  real x1
  real xu
!
  ierr = 0
  s = 0
  tag = 0
  t1 = lb
  t2 = ub
!
!  Look for small sub-diagonal entries
!
  do i = 1, n
     if (i == 1) then
       e2(i) = 0.0
     else
       tst1 = abs(d(i)) + abs(d(i-1))
       tst2 = tst1 + abs(e(i))
       if ( tst2 <= tst1) then
         e2(i) = 0.0
       end if
    end if
  end do
!
!  Determine the number of eigenvalues in the interval.
!
  p = 1
  q = n
  x1 = ub
  isturm = 1
  go to 320

60 continue

  m = s
  x1 = lb
  isturm = 2
  go to 320

80 continue

  m = m - s
  if (m > mm) go to 980
  q = 0
  r = 0
!
!  Establish and process next submatrix, refining
!  interval by the gerschgorin bounds
!
100 continue

  if (r == m) go to 1001
  tag = tag + 1
  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = 0.0

  do q = p, n

    x1 = u
    u = 0.0
    v = 0.0

    if (q /= n) then
      u = abs(e(q+1))
      v = e2(q+1)
    end if

    xu = min (d(q)-(x1+u), xu )
    x0 = max (d(q)+(x1+u), x0 )

    if ( v == 0.0 ) then
      exit
    end if

  end do

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( 1.0 )
  if (eps1 <= 0.0) eps1 = -x1
  if (p /= q) go to 180
!
!  Check for isolated root within interval
!
  if (t1 > d(p) .or. d(p) >= t2) go to 940
  m1 = p
  m2 = p
  rv5(p) = d(p)
  go to 900

  180 continue

  x1 = x1 * (q - p + 1)
  lb = max(t1,xu-x1)
  ub = min(t2,x0+x1)
  x1 = lb
  isturm = 3
  go to 320

  200 continue

  m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320

  220 continue

  m2 = s
  if (m1 > m2) go to 940
!
!  Find roots by bisection.
!
  x0 = ub
  isturm = 5
  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for k-th eigenvalue
!
  k = m2

  250 continue

     xu = lb

     do ii = m1, k
        i = m1 + k - ii
        if (xu < rv4(i)) then
          xu = rv4(i)
          go to 280
        end if
     end do

  280 continue

   if (x0 > rv5(k)) x0 = rv5(k)
!
!  next bisection step
!
  300    continue

     x1 = (xu + x0) * 0.5

     if ((x0 - xu) <= abs(eps1)) go to 420

     tst1 = 2.0 * (abs(xu) + abs(x0))
     tst2 = tst1 + (x0 - xu)
     if (tst2 == tst1) go to 420
!
!  Sturm sequence
!
  320    s = p - 1
     u = 1.0

     do i = p, q

        if ( u == 0.0 ) then
          v = abs(e(i)) / epsilon ( 1.0 )
          if (e2(i) == 0.0) v = 0.0
        else
          v = e2(i) / u
        end if

        u = d(i) - x1 - v
        if (u < 0.0) s = s + 1

     end do

     go to (60,80,200,220,360), isturm
!
!  refine intervals
!
  360    if (s >= k) go to 400
     xu = x1
     if (s >= m1) go to 380
     rv4(m1) = x1
     go to 300
  380    rv4(s+1) = x1
     if (rv5(s) > x1) rv5(s) = x1
     go to 300
  400    x0 = x1
     go to 300
!
!    k-th eigenvalue found
!
  420    rv5(k) = x1
  k = k - 1
  if (k >= m1) go to 250
!
!  order eigenvalues tagged with their submatrix associations
!
  900 s = r
  r = r + m2 - m1 + 1
  j = 1
  k = m1

  do l = 1, r

     if (j > s) go to 910
     if (k > m2) go to 940
     if (rv5(k) >= w(l)) go to 915

     do ii = j, s
        i = l + s - ii
        w(i+1) = w(i)
        ind(i+1) = ind(i)
     end do

  910    w(l) = rv5(k)
     ind(l) = tag
     k = k + 1
     go to 920
  915    j = j + 1

  920 continue

  end do

  940 if (q < n) go to 100
  go to 1001
!
!  set error: underestimate of number of eigenvalues in interval
!
  980 ierr = 3 * n + 1

 1001 continue

  lb = t1
  ub = t2
  return
end
subroutine bqr ( nm, n, mb, a, t, r, ierr, nv, rv )
!
!*******************************************************************************
!
!! BQR finds the eigenvalue of smallest magnitude, for a real symmetric band matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure bqr,
!     num. math. 16, 85-92(1970) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol ii-linear algebra, 266-272(1971).
!
!     this subroutine finds the eigenvalue of smallest (usually)
!     magnitude of a real symmetric band matrix using the
!     qr algorithm with shifts of origin.  consecutive calls
!     can be made to find further eigenvalues.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        mb is the (half) band width of the matrix, defined as the
!          number of adjacent diagonals, including the principal
!          diagonal, required to specify the non-zero portion of the
!          lower triangle of the matrix.
!
!        a contains the lower triangle of the symmetric band input
!          matrix stored as an n by mb array.  its lowest subdiagonal
!          is stored in the last n+1-mb positions of the first column,
!          its next subdiagonal in the last n+2-mb positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the n positions of the last column.
!          contents of storages not part of the matrix are arbitrary.
!          on a subsequent call, its output contents from the previous
!          call should be passed.
!
!        t specifies the shift (of eigenvalues) applied to the diagonal
!          of a in forming the input matrix. what is actually determined
!          is the eigenvalue of a+ti (i is the identity matrix) nearest
!          to t.  on a subsequent call, the output value of t from the
!          previous call should be passed if the next nearest eigenvalue
!          is sought.
!
!        r should be specified as zero on the first call, and as its
!          output value from the previous call on a subsequent call.
!          it is used to determine when the last row and column of
!          the transformed band matrix can be regarded as negligible.
!
!        nv must be set to the dimension of the array parameter rv
!          as declared in the calling program dimension statement.
!
!     on output
!
!        a contains the transformed band matrix.  the matrix a+ti
!          derived from the output parameters is similar to the
!          input a+ti to within rounding errors.  its last row and
!          column are null (if ierr is zero).
!
!        t contains the computed eigenvalue of a+ti (if ierr is zero).
!
!        r contains the maximum of its input value and the norm of the
!          last column of the input matrix a.
!
!        ierr is set to
!          zero       for normal return,
!          n          if the eigenvalue has not been
!                     determined after 30 iterations.
!
!        rv is a temporary storage array of dimension at least
!          (2*mb**2+4*mb-3).  the first (3*mb-2) locations correspond
!          to the algol array b, the next (2*mb-1) locations correspond
!          to the algol array h, and the final (2*mb**2-mb) locations
!          correspond to the mb by (2*mb-1) algol array u.
!
!     note. for a subsequent call, n should be replaced by n-1, but
!     mb should not be altered even when it exceeds the current n.
!
  integer mb
  integer nm
  integer nv
!
  real a(nm,mb)
  real f
  real g
  integer i
  integer ierr
  integer ii
  integer ik
  integer imult
  integer its
  integer j
  integer jk
  integer jm
  integer k
  integer kj
  integer kj1
  integer kk
  integer km
  integer l
  integer ll
  integer m
  integer m1
  integer m2
  integer m21
  integer m3
  integer m31
  integer m4
  integer mk
  integer mn
  integer mz
  integer n
  integer ni
  real pythag
  real q
  real r
  real rv(nv)
  real s
  real scale
  real t
  real tst1
  real tst2
!
  ierr = 0
  m1 = min ( mb, n )
  m = m1 - 1
  m2 = m + m
  m21 = m2 + 1
  m3 = m21 + m
  m31 = m3 + 1
  m4 = m31 + m2
  mn = m + n
  mz = mb - m1
  its = 0
!
!  Test for convergence
!
   40 continue

  g = a(n,mb)
  if (m == 0) go to 360

  f = 0.0
  do k = 1, m
    mk = k + mz
    f = f + abs ( a(n,mk) )
  end do

  if (its == 0 .and. f > r) r = f
  tst1 = r
  tst2 = tst1 + f
  if (tst2 <= tst1) go to 360

  if ( its >= 30 ) then
    ierr = n
    return
  end if

  its = its + 1
!
!  Form shift from bottom 2 by 2 minor
!
  if ( f <= 0.25 * r .or. its >= 5 ) then

    f = a(n,mb-1)

    if ( f /= 0.0 ) then
      q = (a(n-1,mb) - g) / (2.0 * f)
      s = pythag(q,1.0)
      g = g - f / (q + sign(s,q))
    end if

    t = t + g

    a(1:n,mb) = a(1:n,mb) - g

  end if

  rv(m31:m4) = 0.0

  do ii = 1, mn

     i = ii - m
     ni = n - ii
     if (ni < 0) go to 230
!
!  Form column of shifted matrix a-g*i
!
     l = max(1,2-i)

     rv(1:m3) = 0.0

     do k = l, m1
        km = k + m
        mk = k + mz
        rv(km) = a(ii,mk)
     end do

     ll = min(m,ni)
     if (ll == 0) go to 135

     do k = 1, ll
        km = k + m21
        ik = ii + k
        mk = mb - k
        rv(km) = a(ik,mk)
     end do
!
!  Pre-multiply with Householder reflections
!
  135    ll = m2
     imult = 0
!
!  Multiplication procedure
!
  140    kj = m4 - m1

     do j = 1, ll

        kj = kj + m1
        jm = j + m3
        if (rv(jm) == 0.0) go to 170
        f = 0.0

        do k = 1, m1
           kj = kj + 1
           jk = j + k - 1
           f = f + rv(kj) * rv(jk)
        end do

        f = f / rv(jm)
        kj = kj - m1

        do k = 1, m1
           kj = kj + 1
           jk = j + k - 1
           rv(jk) = rv(jk) - rv(kj) * f
        end do

        kj = kj - m1

  170    continue

     end do

     if (imult /= 0) go to 280
!
!  Householder reflection
!
     f = rv(m21)
     s = 0.0
     rv(m4) = 0.0
     scale = 0.0

     do k = m21, m3
       scale = scale + abs(rv(k))
     end do

     if (scale == 0.0) go to 210

     do k = m21, m3
       s = s + (rv(k)/scale)**2
     end do

     s = scale * scale * s
     g = -sign(sqrt(s),f)
     rv(m21) = g
     rv(m4) = s - f * g
     kj = m4 + m2 * m1 + 1
     rv(kj) = f - g

     do k = 2, m1
        kj = kj + 1
        km = k + m2
        rv(kj) = rv(km)
     end do
!
!  Save column of triangular factor r
!
  210    continue

     do k = l, m1
        km = k + m
        mk = k + mz
        a(ii,mk) = rv(km)
     end do

  230    l = max(1,m1+1-i)
     if (i <= 0) go to 300
!
!  Perform additional steps
!
     rv(1:m21) = 0.0
     ll = min(m1,ni+m1)
!
!  Get row of triangular factor r
!
     do kk = 1, ll
        k = kk - 1
        km = k + m1
        ik = i + k
        mk = mb - k
        rv(km) = a(ik,mk)
     end do
!
!  Post-multiply with Householder reflections
!
     ll = m1
     imult = 1
     go to 140
!
!  Store column of new a matrix
!
  280    continue

     do k = l, m1
        mk = k + mz
        a(i,mk) = rv(k)
     end do
!
!  Update Householder reflections.
!
  300    if (l > 1) l = l - 1
     kj1 = m4 + l * m1

     do j = l, m2

        jm = j + m3
        rv(jm) = rv(jm+1)

        do k = 1, m1
           kj1 = kj1 + 1
           kj = kj1 - m1
           rv(kj) = rv(kj1)
        end do

     end do

  end do

  go to 40
!
!  Convergence
!
  360 t = t + g

  a(1:n,mb) = a(1:n,mb) - g

  do k = 1, m1
     mk = k + mz
     a(n,mk) = 0.0
  end do

  return
end
subroutine cbabk2 ( nm, n, low, igh, scale, m, zr, zi )
!
!*******************************************************************************
!
!! CBABK2 finds eigenvectors by undoing the CBAL transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure
!     cbabk2, which is a complex version of balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  cbal.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  cbal.
!
!        scale contains information determining the permutations
!          and scaling factors used by  cbal.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
  integer m
  integer n
  integer nm
!
  integer i
  integer igh
  integer ii
  integer j
  integer k
  integer low
  real s
  real scale(n)
  real zi(nm,m)
  real zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  if ( igh /= low ) then

    do i = low, igh

      s = scale(i)

      do j = 1, m
        zr(i,j) = zr(i,j) * s
        zi(i,j) = zi(i,j) * s
      end do

    end do

  end if

  do ii = 1, n

    i = ii

    if (i < low .or. i > igh) then

      if (i < low) i = low - ii
      k = scale(i)

      if ( k /= i ) then

        do j = 1, m
          call r_swap ( zr(i,j), zr(k,j) )
          call r_swap ( zi(i,j), zi(k,j) )
        end do

      end if

    end if

  end do

  return
end
subroutine cbal ( nm, n, ar, ai, low, igh, scale )
!
!*******************************************************************************
!
!! CBAL balances a complex matrix before eigenvalue calculations.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure
!     cbalance, which is a complex version of balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a complex matrix and isolates
!     eigenvalues whenever possible.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j)       j = low,...,igh
!                 = p(j)         j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex matrix to be balanced.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the balanced matrix.
!
!        low and igh are two integers such that ar(i,j) and ai(i,j)
!          are equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  real b2
  real c
  real f
  real g
  integer i
  integer iexc
  integer igh
  integer j
  integer jj
  integer k
  integer l
  integer low
  integer m
  logical noconv
  real r
  real radix
  real s
  real scale(n)
!
  radix = 16.0

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

20 continue

  scale(m) = j

  if ( j /= m ) then

    do i = 1, l
      call r_swap ( ar(i,j), ar(i,m) )
      call r_swap ( ai(i,j), ai(i,m) )
    end do

    do i = k, n
      call r_swap ( ar(j,i), ar(m,i) )
      call r_swap ( ai(j,i), ai(m,i) )
    end do

  end if

  if ( iexc == 2 ) then
    go to 130
  end if
!
!  Search for rows isolating an eigenvalue and push them down
!
   80 continue

  if (l == 1) go to 280
  l = l - 1

  100 continue

  do jj = 1, l

     j = l + 1 - jj

     do i = 1, l
        if (i == j) go to 110
        if (ar(j,i) /= 0.0 .or. ai(j,i) /= 0.0) go to 120
  110    continue
     end do

     m = l
     iexc = 1
     go to 20

120      continue

  end do

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left
!
130 continue

  k = k + 1

140 continue

   do j = k, l

     do i = k, l
        if (i == j) go to 150
        if (ar(i,j) /= 0.0 .or. ai(i,j) /= 0.0) go to 170
  150    continue
     end do

     m = k
     iexc = 2
     go to 20
  170 continue

  end do
!
!  Now balance the submatrix in rows k to l
!
  scale(k:l) = 1.0
!
!  Iterative loop for norm reduction
!
  190 continue

  noconv = .false.

  do i = k, l

    c = 0.0
    r = 0.0

    do j = k, l
      if (j /= i) then
        c = c + abs(ar(j,i)) + abs(ai(j,i))
        r = r + abs(ar(i,j)) + abs(ai(i,j))
      end if
    end do
!
!  Guard against zero c or r due to underflow
!
     if ( c == 0.0 .or. r == 0.0 ) go to 270

     g = r / radix
     f = 1.0
     s = c + r

     do while ( c < g )
       f = f * radix
       c = c * b2
     end do

     g = r * radix

     do while  ( c >= g )
       f = f / radix
       c = c / b2
     end do
!
!  Now balance
!
     if ((c + r) / f >= 0.95 * s) go to 270

     g = 1.0 / f
     scale(i) = scale(i) * f
     noconv = .true.

     do j = k, n
        ar(i,j) = ar(i,j) * g
        ai(i,j) = ai(i,j) * g
     end do

     do j = 1, l
        ar(j,i) = ar(j,i) * f
        ai(j,i) = ai(j,i) * f
     end do

  270 continue

  end do

  if (noconv) go to 190

  280 continue

  low = k
  igh = l
  return
end
subroutine cdiv ( ar, ai, br, bi, cr, ci )
!
!*******************************************************************************
!
!! CDIV emulates complex division, using real arithmetic.
!
!
!  Discussion:
!
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
  real ai
  real ais
  real ar
  real ars
  real bi
  real bis
  real br
  real brs
  real ci
  real cr
  real s
!
  s = abs ( br ) + abs ( bi )

  ars = ar / s
  ais = ai / s
  brs = br / s
  bis = bi / s

  s = brs**2 + bis**2
  cr = ( ars * brs + ais * bis ) / s
  ci = ( ais * brs - ars * bis ) / s

  return
end
subroutine cg ( nm, n, ar, ai, wr, wi, matz, zr, zi, ierr )
!
!*******************************************************************************
!
!! CG gets eigenvalues and eigenvectors of a complex general matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex general matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for comqr
!           and COMQR2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  real fv1(n)
  real fv2(n)
  real fv3(n)
  integer ierr
  integer is1
  integer is2
  integer matz
  real wi(n)
  real wr(n)
  real zi(nm,n)
  real zr(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call cbal ( nm, n, ar, ai, is1, is2, fv1 )

  call corth ( nm, n, is1, is2, ar, ai, fv2, fv3 )

  if ( matz == 0 ) then

    call comqr ( nm, n, is1, is2, ar, ai, wr, wi, ierr )

    if ( ierr /= 0 ) then
      return
    end if
  
  else

    call comqr2 ( nm, n, is1, is2, fv2, fv3, ar, ai, wr, wi, zr, zi, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CG - Fatal error!'
      write ( *, * ) '  Nonzero error return from COMQR2.'
      return
    end if

    call cbabk2 ( nm, n, is1, is2, fv1, n, zr, zi )

  end if

  return
end
subroutine ch ( nm, n, ar, ai, w, matz, zr, zi, ierr )
!
!*******************************************************************************
!
!! CH gets eigenvalues and eigenvectors of a complex Hermitian matrix. .
!
!
!  Discussion:
!
!    This subroutine calls the recommended sequence of subroutines from the 
!    EISPACK eigensystem package to find the eigenvalues and eigenvectors 
!    of a complex hermitian matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex hermitian matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and TQL2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  real fm1(2,n)
  real fv1(n)
  real fv2(n)
  integer ierr
  integer matz
  real w(n)
  real zi(nm,n)
  real zr(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call htridi ( nm, n, ar, ai, w, fv1, fv2, fm1 )

  if ( matz == 0 ) then

    call tqlrat ( n, w, fv2, ierr )

  else

    call rmat_ident ( nm, n, zr )

    call tql2 ( nm, n, w, fv1, zr, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call htribk ( nm, n, ar, ai, fm1, n, zr, zi )

  end if

  return
end
subroutine cinvit ( nm, n, ar, ai, wr, wi, select, mm, m, zr, zi, ierr, &
  rm1, rm2, rv1, rv2 )
!
!*******************************************************************************
!
!! CINVIT gets eigenvectors from eigenvalues, for a complex Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure cx invit
!     by peters and wilkinson.
!     handbook for auto. comp. vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvectors of a complex upper
!     Hessenberg matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the Hessenberg matrix.
!
!        wr and wi contain the real and imaginary parts, respectively,
!          of the eigenvalues of the matrix.  the eigenvalues must be
!          stored in a manner identical to that of subroutine  comlr,
!          which recognizes possible splitting of the matrix.
!
!        select specifies the eigenvectors to be found.  the
!          eigenvector corresponding to the j-th eigenvalue is
!          specified by setting select(j) to .true..
!
!        mm should be set to an upper bound for the number of
!          eigenvectors to be found.
!
!     on output
!
!        ar, ai, wi, and select are unaltered.
!
!        wr may have been altered since close eigenvalues are perturbed
!          slightly in searching for independent eigenvectors.
!
!        m is the number of eigenvectors actually found.
!
!        zr and zi contain the real and imaginary parts, respectively,
!          of the eigenvectors.  the eigenvectors are normalized
!          so that the component of largest magnitude is 1.
!          any vector which fails the acceptance test is set to zero.
!
!        ierr is set to
!          zero       for normal return,
!          -(2*n+1)   if more than mm eigenvectors have been specified,
!          -k         if the iteration corresponding to the k-th
!                     value fails,
!          -(n+k)     if both error situations occur.
!
!        rm1, rm2, rv1, and rv2 are temporary storage arrays.
!
  integer mm
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  real eps3
  real growto
  integer i
  integer ierr
  integer ii
  real ilambd
  integer its
  integer j
  integer k
  integer km1
  integer m
  integer mp
  real norm
  real normv
  real pythag
  real rlambd
  real rm1(n,n)
  real rm2(n,n)
  real rv1(n)
  real rv2(n)
  integer s
  logical select(n)
  integer uk
  real ukroot
  real wi(n)
  real wr(n)
  real x
  real y
  real zi(nm,mm)
  real zr(nm,mm)
!
  ierr = 0
  uk = 0
  s = 1

  do k = 1, n

     if (.not. select(k)) go to 980

     if (s > mm) go to 1000
     if (uk >= k) go to 200
!
!  Check for possible splitting
!
     do uk = k, n - 1

       if (ar(uk+1,uk) == 0.0 .and. ai(uk+1,uk) == 0.0) then
         exit
       end if

     end do
!
!  Compute infinity norm of leading uk by uk (Hessenberg) matrix
!   
     norm = 0.0
     mp = 1

     do i = 1, uk

        x = 0.0
        do j = mp, uk
          x = x + pythag(ar(i,j),ai(i,j))
        end do

        norm = max ( norm, x )
        mp = i

     end do
!
!  eps3 replaces zero pivot in decomposition
!                and close roots are modified by eps3
!
     if (norm == 0.0) norm = 1.0
     eps3 = abs ( norm ) * epsilon ( 1.0 )
!
!  GROWTO is the criterion for growth.
!
     ukroot = uk
     ukroot = sqrt(ukroot)
     growto = 0.1 / ukroot
  200    rlambd = wr(k)
     ilambd = wi(k)
     if (k == 1) go to 280
     km1 = k - 1
     go to 240
!
!  Perturb eigenvalue if it is close to any previous eigenvalue
!
  220    rlambd = rlambd + eps3

  240    continue

     do ii = 1, km1
        i = k - ii
        if (select(i) .and. abs(wr(i)-rlambd) < eps3 .and. &
            abs(wi(i)-ilambd) < eps3) then
          go to 220
        end if
     end do

     wr(k) = rlambd
!
!  Form upper Hessenberg (ar,ai)-(rlambd,ilambd)*i
!  and initial complex vector
!
  280    mp = 1

     do i = 1, uk

        do j = mp, uk
          rm1(i,j) = ar(i,j)
          rm2(i,j) = ai(i,j)
        end do

        rm1(i,i) = rm1(i,i) - rlambd
        rm2(i,i) = rm2(i,i) - ilambd
        mp = i
        rv1(i) = eps3

     end do
!
!  Triangular decomposition with interchanges, replacing zero pivots by eps3
!
     do i = 2, uk

        mp = i - 1

        if ( pythag(rm1(i,mp),rm2(i,mp)) > pythag(rm1(mp,mp),rm2(mp,mp)) ) then

          do j = mp, uk
            call r_swap ( rm1(i,j), rm1(mp,j) )
            call r_swap ( rm2(i,j), rm2(mp,j) )
          end do

        end if

        if (rm1(mp,mp) == 0.0 .and. rm2(mp,mp) == 0.0) then
          rm1(mp,mp) = eps3
        end if

        call cdiv(rm1(i,mp),rm2(i,mp),rm1(mp,mp),rm2(mp,mp),x,y)

        if ( x /= 0.0 .or. y /= 0.0 ) then

          do j = i, uk
            rm1(i,j) = rm1(i,j) - x * rm1(mp,j) + y * rm2(mp,j)
            rm2(i,j) = rm2(i,j) - x * rm2(mp,j) - y * rm1(mp,j)
          end do

        end if

     end do

     if (rm1(uk,uk) == 0.0 .and. rm2(uk,uk) == 0.0) then
       rm1(uk,uk) = eps3
     end if

     its = 0
!
!  Back substitution
!
  660   continue

    do ii = 1, uk

        i = uk + 1 - ii
        x = rv1(i)
        y = 0.0

        do j = i+1, uk
          x = x - rm1(i,j) * rv1(j) + rm2(i,j) * rv2(j)
          y = y - rm1(i,j) * rv2(j) - rm2(i,j) * rv1(j)
        end do

        call cdiv(x,y,rm1(i,i),rm2(i,i),rv1(i),rv2(i))

     end do
!
!  Acceptance test for eigenvector and normalization
!
     its = its + 1
     norm = 0.0
     normv = 0.0

     do i = 1, uk
        x = pythag(rv1(i),rv2(i))
        if ( normv < x ) then
          normv = x
          j = i
        end if
        norm = norm + x
     end do

     if (norm < growto) go to 840
!
!  Accept vector
!
     x = rv1(j)
     y = rv2(j)

     do i = 1, uk
       call cdiv(rv1(i),rv2(i),x,y,zr(i,s),zi(i,s))
     end do

     if (uk == n) go to 940
     j = uk + 1
     go to 900
!
!  Choose a new starting vector
!
  840    continue

     if ( its < uk ) then

       x = ukroot
       y = eps3 / (x + 1.0)

       rv1(1) = eps3
       rv1(2:uk) = y

       j = uk - its + 1
       rv1(j) = rv1(j) - eps3 * x
       go to 660

     end if
!
!  Error: unaccepted eigenvector
!
  880    j = 1
     ierr = -k
!
!  Set remaining vector components to zero
!
  900    continue

       zr(j:n,s) = 0.0
       zi(j:n,s) = 0.0

  940    s = s + 1
  980 continue

  end do

  go to 1001
!
!  Set error: underestimate of eigenvector space required
!
 1000 if (ierr /= 0) ierr = ierr - n
  if (ierr == 0) ierr = -(2 * n + 1)
 1001 continue
  m = s - 1
  return
end
subroutine combak ( nm, low, igh, ar, ai, int, m, zr, zi )
!
!*******************************************************************************
!
!! COMBAK determines eigenvectors by undoing the COMHES transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure combak,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     upper Hessenberg matrix determined by  comhes.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        ar and ai contain the multipliers which were used in the
!          reduction by  comhes  in their lower triangles
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  comhes.
!          only elements low through igh are used.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
  integer igh
  integer m
  integer nm
!
  real ai(nm,igh)
  real ar(nm,igh)
  integer i
  integer int(igh)
  integer j
  integer la
  integer low
  integer mm
  integer mp
  real xi
  real xr
  real zi(nm,m)
  real zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  la = igh - 1

  if ( igh - 1 < low + 1 ) then
    return
  end if

  do mm = low + 1, la

     mp = low + igh - mm

     do i = mp+1, igh

        xr = ar(i,mp-1)
        xi = ai(i,mp-1)

        if ( xr /= 0.0 .or. xi /= 0.0 ) then

          do j = 1, m
            zr(i,j) = zr(i,j) + xr * zr(mp,j) - xi * zi(mp,j)
            zi(i,j) = zi(i,j) + xr * zi(mp,j) + xi * zr(mp,j)
          end do

       end if

     end do

     i = int(mp)

     if ( i /= mp ) then

       do j = 1, m
         call r_swap ( zr(i,j), zr(mp,j) )
         call r_swap ( zi(i,j), zi(mp,j) )
       end do

     end if

  end do

  return
end
subroutine comhes ( nm, n, low, igh, ar, ai, int )
!
!*******************************************************************************
!
!! COMHES transforms a complex general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure comhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a complex general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper Hessenberg form by
!     stabilized elementary similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex input matrix.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the Hessenberg matrix.  the
!          multipliers which were used in the reduction
!          are stored in the remaining triangles under the
!          Hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
  integer igh
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  integer i
  integer int(igh)
  integer j
  integer la
  integer low
  integer m
  integer mm1
  real xi
  real xr
  real yi
  real yr
!
  la = igh - 1

  do m = low + 1, la

     mm1 = m - 1
     xr = 0.0
     xi = 0.0
     i = m

     do j = m, igh

        if (abs(ar(j,mm1)) + abs(ai(j,mm1)) > abs(xr) + abs(xi)) then
          xr = ar(j,mm1)
          xi = ai(j,mm1)
          i = j
        end if

     end do

     int(m) = i
!
!  interchange rows and columns of ar and ai
!
     if ( i /= m ) then

       do j = mm1, n
         call r_swap ( ar(i,j), ar(m,j) )
         call r_swap ( ai(i,j), ai(m,j) )
       end do

       do j = 1, igh
         call r_swap ( ar(j,i), ar(j,m) )
         call r_swap ( ai(j,i), ai(j,m) )
       end do

     end if

    if ( xr /= 0.0 .or. xi /= 0.0 ) then

      do i = m+1, igh

        yr = ar(i,mm1)
        yi = ai(i,mm1)

        if ( yr /= 0.0 .or. yi /= 0.0 ) then

          call cdiv(yr,yi,xr,xi,yr,yi)
          ar(i,mm1) = yr
          ai(i,mm1) = yi

          do j = m, n
            ar(i,j) = ar(i,j) - yr * ar(m,j) + yi * ai(m,j)
            ai(i,j) = ai(i,j) - yr * ai(m,j) - yi * ar(m,j)
          end do

          do j = 1, igh
            ar(j,m) = ar(j,m) + yr * ar(j,i) - yi * ai(j,i)
            ai(j,m) = ai(j,m) + yr * ai(j,i) + yi * ar(j,i)
          end do

        end if

      end do

    end if

  end do

  return
end
subroutine comlr ( nm, n, low, igh, hr, hi, wr, wi, ierr )
!
!*******************************************************************************
!
!! COMLR gets all eigenvalues of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure comlr,
!     num. math. 12, 369-376(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
!
!     this subroutine finds the eigenvalues of a complex
!     upper Hessenberg matrix by the modified lr method.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper Hessenberg matrix.
!          their lower triangles below the subdiagonal contain the
!          multipliers which were used in the reduction by  comhes,
!          if performed.
!
!     on output
!
!        the upper Hessenberg portions of hr and hi have been
!          destroyed.  therefore, they must be saved before
!          calling  comlr  if subsequent calculation of
!          eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
  integer n
  integer nm
!
  integer en
  integer enm1
  real hi(nm,n)
  real hr(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer l
  integer ll
  integer low
  integer m
  integer mm
  real si
  real sr
  real ti
  real tr
  real tst1
  real tst2
  real wi(n)
  real wr(n)
  real xi
  real xr
  real yi
  real yr
  real zzi
  real zzr
!
  ierr = 0
!
!  Store roots isolated by CBAL.
!
  do i = 1, n
     if (i < low .or. i > igh ) then
       wr(i) = hr(i,i)
       wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = 0.0
  ti = 0.0
  itn = 30*n
!
!  search for next eigenvalue
!
  220 continue

  if (en < low) then
    return
  end if

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element
!
  240 continue

  do ll = low, en
     l = en + low - ll
     if (l == low) go to 300
     tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
     tst2 = tst1 + abs(hr(l,l-1)) + abs(hi(l,l-1))
     if (tst2 == tst1) go to 300
  end do
!
!  Form shift
!
  300 if (l == en) go to 660

  if (itn == 0) then
    ierr = en
    return
  end if

  if (its == 10 .or. its == 20) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
  xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
  if (xr == 0.0 .and. xi == 0.0) go to 340
  yr = (hr(enm1,enm1) - sr) / 2.0
  yi = (hi(enm1,enm1) - si) / 2.0
  call csroot(yr**2-yi**2+xr,2.0*yr*yi+xi,zzr,zzi)
  if (yr * zzr + yi * zzi >= 0.0) go to 310
  zzr = -zzr
  zzi = -zzi

  310 continue

  call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift
!
  320 continue

  sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
  si = abs(hi(en,enm1)) + abs(hi(enm1,en-2))

  340 continue

  do i = low, en
     hr(i,i) = hr(i,i) - sr
     hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements
!
  xr = abs(hr(enm1,enm1)) + abs(hi(enm1,enm1))
  yr = abs(hr(en,enm1)) + abs(hi(en,enm1))
  zzr = abs(hr(en,en)) + abs(hi(en,en))

  do mm = l, enm1
     m = enm1 + l - mm
     if (m == l) go to 420
     yi = yr
     yr = abs(hr(m,m-1)) + abs(hi(m,m-1))
     xi = zzr
     zzr = xr
     xr = abs(hr(m-1,m-1)) + abs(hi(m-1,m-1))
     tst1 = zzr / yi * (zzr + xr + xi)
     tst2 = tst1 + yr
     if (tst2 == tst1) go to 420
  end do
!
!  Triangular decomposition h=l*r
!
  420 continue

  do i = m+1, en

     xr = hr(i-1,i-1)
     xi = hi(i-1,i-1)
     yr = hr(i,i-1)
     yi = hi(i,i-1)
     if (abs(xr) + abs(xi) >= abs(yr) + abs(yi)) go to 460
!
!    interchange rows of hr and hi
!
     do j = i-1, en
       call r_swap ( hr(i-1,j), hr(i,j) )
       call r_swap ( hi(i-1,j), hi(i,j) )
     end do

     call cdiv(xr,xi,yr,yi,zzr,zzi)
     wr(i) = 1.0
     go to 480

460 continue

     call cdiv(yr,yi,xr,xi,zzr,zzi)
     wr(i) = -1.0

480  continue

     hr(i,i-1) = zzr
     hi(i,i-1) = zzi

     do j = i, en
        hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
        hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
     end do

  end do
!
!  Composition r*l=h
!
  do j = m+1, en

    xr = hr(j,j-1)
    xi = hi(j,j-1)
    hr(j,j-1) = 0.0
    hi(j,j-1) = 0.0
!
!  Interchange columns of hr and hi, if necessary
!
    if ( wr(j) > 0.0 ) then

      do i = l, j
        call r_swap ( hr(i,j-1), hr(i,j) )
        call r_swap ( hi(i,j-1), hi(i,j) )
      end do

    end if

    do i = l, j
      hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
      hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
    end do

  end do

  go to 240
!
!    a root found
!
  660 continue

  wr(en) = hr(en,en) + tr
  wi(en) = hi(en,en) + ti
  en = enm1
  go to 220
end
subroutine comlr2 ( nm, n, low, igh, int, hr, hi, wr, wi, zr,  zi, ierr )
!
!*******************************************************************************
!
!! COMLR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure comlr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a complex upper Hessenberg matrix by the modified lr
!     method.  the eigenvectors of a complex general matrix
!     can also be found if  comhes  has been used to reduce
!     this general matrix to Hessenberg form.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        int contains information on the rows and columns interchanged
!          in the reduction by  comhes, if performed.  only elements
!          low through igh are used.  if the eigenvectors of the hessen-
!          berg matrix are desired, set int(j)=j for these elements.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper Hessenberg matrix.
!          their lower triangles below the subdiagonal contain the
!          multipliers which were used in the reduction by  comhes,
!          if performed.  if the eigenvectors of the Hessenberg
!          matrix are desired, these elements must be set to zero.
!
!     on output
!
!        the upper Hessenberg portions of hr and hi have been
!          destroyed, but the location hr(1,1) contains the norm
!          of the triangularized matrix.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors.  the eigenvectors
!          are unnormalized.  if an error exit is made, none of
!          the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
  integer n
  integer nm
!
  integer en
  integer enm1
  real hi(nm,n)
  real hr(nm,n)
  integer i
  integer iend
  integer ierr
  integer igh
  integer ii
  integer int(igh)
  integer itn
  integer its
  integer j
  integer jj
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer nn
  real norm
  real si
  real sr
  real ti
  real tr
  real tst1
  real tst2
  real wi(n)
  real wr(n)
  real xi
  real xr
  real yi
  real yr
  real zi(nm,n)
  real zr(nm,n)
  real zzi
  real zzr
!
  ierr = 0
!
!  Initialize the eigenvector matrix.
!
  call rmat_ident ( nm, n, zr )

  zi(1:n,1:n) = 0.0
!
!  Form the matrix of accumulated transformations from the information left 
!  by COMHES.
!
  iend = igh - low - 1

  do ii = 1, iend

    i = igh - ii

    do k = i+1, igh
      zr(k,i) = hr(k,i-1)
      zi(k,i) = hi(k,i-1)
    end do

    j = int(i)

    if ( i /= j ) then

      do k = i, igh
        zr(i,k) = zr(j,k)
        zi(i,k) = zi(j,k)
        zr(j,k) = 0.0
        zi(j,k) = 0.0
      end do

      zr(j,i) = 1.0

    end if

  end do
!
!  Store roots isolated by CBAL.
!
  do i = 1, n
    if ( i < low .or. i > igh ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = 0.0
  ti = 0.0
  itn = 30*n
!
!  Search for next eigenvalue
!
  220 continue

  if (en < low) go to 680
  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element
!
  240 continue

  do ll = low, en
     l = en + low - ll
     if (l == low) go to 300
     tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
     tst2 = tst1 + abs(hr(l,l-1)) + abs(hi(l,l-1))
     if (tst2 == tst1) go to 300
  end do
!
!  Form shift.
!
  300 continue

  if (l == en) go to 660
  if (itn == 0) go to 1000
  if (its == 10 .or. its == 20) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
  xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
  if (xr == 0.0 .and. xi == 0.0) go to 340
  yr = (hr(enm1,enm1) - sr) / 2.0
  yi = (hi(enm1,enm1) - si) / 2.0
  call csroot(yr**2-yi**2+xr,2.0*yr*yi+xi,zzr,zzi)
  if (yr * zzr + yi * zzi >= 0.0) go to 310
  zzr = -zzr
  zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
  sr = sr - xr
  si = si - xi
  go to 340
!    form exceptional shift
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
  si = abs(hi(en,enm1)) + abs(hi(enm1,en-2))

  340 continue

  do i = low, en
     hr(i,i) = hr(i,i) - sr
     hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  xr = abs(hr(enm1,enm1)) + abs(hi(enm1,enm1))
  yr = abs(hr(en,enm1)) + abs(hi(en,enm1))
  zzr = abs(hr(en,en)) + abs(hi(en,en))

  do mm = l, enm1
     m = enm1 + l - mm
     if (m == l) go to 420
     yi = yr
     yr = abs(hr(m,m-1)) + abs(hi(m,m-1))
     xi = zzr
     zzr = xr
     xr = abs(hr(m-1,m-1)) + abs(hi(m-1,m-1))
     tst1 = zzr / yi * (zzr + xr + xi)
     tst2 = tst1 + yr
     if (tst2 == tst1) go to 420
  end do
!
!  Triangular decomposition h=l*r
!
  420 continue

  do i = m+1, en

     xr = hr(i-1,i-1)
     xi = hi(i-1,i-1)
     yr = hr(i,i-1)
     yi = hi(i,i-1)
     if (abs(xr) + abs(xi) >= abs(yr) + abs(yi)) go to 460
!
!  Interchange rows of hr and hi
!
     do j = i-1, n
       call r_swap ( hr(i-1,j), hr(i,j) )
       call r_swap ( hi(i-1,j), hi(i,j) )
    end do

     call cdiv(xr,xi,yr,yi,zzr,zzi)
     wr(i) = 1.0
     go to 480
  460    call cdiv(yr,yi,xr,xi,zzr,zzi)
     wr(i) = -1.0
  480    hr(i,i-1) = zzr
     hi(i,i-1) = zzi

     do j = i, n
        hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
        hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
     end do

  end do
!
!  Composition r*l=h
!
  do j = m+1, en

     xr = hr(j,j-1)
     xi = hi(j,j-1)
     hr(j,j-1) = 0.0
     hi(j,j-1) = 0.0
!
!  Interchange columns of hr, hi, zr, and zi.
!
     if ( wr(j) > 0.0 ) then

       do i = 1, j
         call r_swap ( hr(i,j-1), hr(i,j) )
         call r_swap ( hi(i,j-1), hi(i,j) )
       end do

       do i = low, igh
         call r_swap ( zr(i,j-1), zr(i,j) )
         call r_swap ( zi(i,j-1), zi(i,j) )
       end do

    end if

    do i = 1, j
      hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
      hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
    end do
!
!  Accumulate transformations
!
    do i = low, igh
      zr(i,j-1) = zr(i,j-1) + xr * zr(i,j) - xi * zi(i,j)
      zi(i,j-1) = zi(i,j-1) + xr * zi(i,j) + xi * zr(i,j)
    end do

  end do

  go to 240
!
!  A root found
!
  660 continue

  hr(en,en) = hr(en,en) + tr
  wr(en) = hr(en,en)
  hi(en,en) = hi(en,en) + ti
  wi(en) = hi(en,en)
  en = enm1
  go to 220
!
!  All roots found.  
!  Backsubstitute to find vectors of upper triangular form
!
  680 continue

  norm = 0.0

  do i = 1, n
    do j = i, n
      tr = abs(hr(i,j)) + abs(hi(i,j))
      if (tr > norm) norm = tr
    end do
  end do

  hr(1,1) = norm
  if (n == 1 ) then
    return
  end if

  if ( norm == 0.0 ) then
    return
  end if

  do nn = 2, n

     en = n + 2 - nn
     xr = wr(en)
     xi = wi(en)
     hr(en,en) = 1.0
     hi(en,en) = 0.0
     enm1 = en - 1

     do ii = 1, enm1

        i = en - ii
        zzr = 0.0
        zzi = 0.0

        do j = i+1, en
           zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
           zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
        end do

        yr = xr - wr(i)
        yi = xi - wi(i)
        if (yr /= 0.0 .or. yi /= 0.0) go to 765
           tst1 = norm
           yr = tst1
  760          yr = 0.01 * yr
           tst2 = norm + yr
           if (tst2 > tst1) go to 760
  765       continue
        call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
!
!  Overflow control
!
        tr = abs(hr(i,en)) + abs(hi(i,en))
        if (tr == 0.0) go to 780
        tst1 = tr
        tst2 = tst1 + 1.0/tst1
        if (tst2 > tst1) go to 780

        do j = i, en
           hr(j,en) = hr(j,en)/tr
           hi(j,en) = hi(j,en)/tr
        end do

  780    continue

      end do

  end do
!
!  End backsubstitution
!
  enm1 = n - 1
!
!  Vectors of isolated roots
!
  do i = 1, n - 1

    if (i < low .or. i > igh) then

      do j = i+1, n
        zr(i,j) = hr(i,j)
        zi(i,j) = hi(i,j)
      end do

    end if

  end do
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  do jj = low, n - 1

     j = n + low - jj
     m = min(j,igh)

     do i = low, igh
        zzr = 0.0
        zzi = 0.0
        do k = low, m
          zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
          zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
        end do
        zr(i,j) = zzr
        zi(i,j) = zzi
      end do
  end do

  return
!
!  Set error: all eigenvalues have not converged after 30*n iterations
!
 1000 ierr = en
  return
end
subroutine comqr ( nm, n, low, igh, hr, hi, wr, wi, ierr )
!
!*******************************************************************************
!
!! COMQR gets eigenvalues of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of a unitary analogue of the
!     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
!     and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
!     the unitary analogue substitutes the qr algorithm of francis
!     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
!
!     this subroutine finds the eigenvalues of a complex
!     upper Hessenberg matrix by the qr method.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper Hessenberg matrix.
!          their lower triangles below the subdiagonal contain
!          information about the unitary transformations used in
!          the reduction by  corth, if performed.
!
!     on output
!
!        the upper Hessenberg portions of hr and hi have been
!          destroyed.  therefore, they must be saved before
!          calling  comqr  if subsequent calculation of
!          eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
  integer n
  integer nm
!
  integer en
  integer enm1
  real hi(nm,n)
  real hr(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer l
  integer ll
  integer low
  integer lp1
  real norm
  real pythag
  real si
  real sr
  real ti
  real tr
  real tst1
  real tst2
  real wi(n)
  real wr(n)
  real xi
  real xr
  real yi
  real yr
  real zzi
  real zzr
!
  ierr = 0
!
!  Create real subdiagonal elements
!
  l = low + 1

  do i = l, igh

     ll = min(i+1,igh)

     if (hi(i,i-1) /= 0.0) then

     norm = pythag(hr(i,i-1),hi(i,i-1))
     yr = hr(i,i-1) / norm
     yi = hi(i,i-1) / norm
     hr(i,i-1) = norm
     hi(i,i-1) = 0.0

     do j = i, igh
        si = yr * hi(i,j) - yi * hr(i,j)
        hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
        hi(i,j) = si
     end do

     do j = low, ll
        si = yr * hi(j,i) + yi * hr(j,i)
        hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
        hi(j,i) = si
     end do

    end if

  end do
!
!  Store roots isolated by CBAL
!
  do i = 1, n
    if (i < low .or. i > igh) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = 0.0
  ti = 0.0
  itn = 30*n
!
!  Search for next eigenvalue
!
  220 continue

  if (en < low) then
    return
  end if

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element
!
  240 continue

  do ll = low, en
     l = en + low - ll
     if (l == low) go to 300
     tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
     tst2 = tst1 + abs(hr(l,l-1))
     if (tst2 == tst1) go to 300
  end do
!
!  Form shift
!
  300 continue

  if (l == en) go to 660
  if (itn == 0) go to 1000
  if (its == 10 .or. its == 20) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1)
  xi = hi(enm1,en) * hr(en,enm1)
  if (xr == 0.0 .and. xi == 0.0) go to 340
  yr = (hr(enm1,enm1) - sr) / 2.0
  yi = (hi(enm1,enm1) - si) / 2.0
  call csroot(yr**2-yi**2+xr,2.0*yr*yi+xi,zzr,zzi)
  if (yr * zzr + yi * zzi >= 0.0) go to 310
  zzr = -zzr
  zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift
!
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
  si = 0.0

  340 continue

  do i = low, en
     hr(i,i) = hr(i,i) - sr
     hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Reduce to triangle (rows)
!
  lp1 = l + 1

  do i = l+1, en

     sr = hr(i,i-1)
     hr(i,i-1) = 0.0
     norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
     xr = hr(i-1,i-1) / norm
     wr(i-1) = xr
     xi = hi(i-1,i-1) / norm
     wi(i-1) = xi
     hr(i-1,i-1) = norm
     hi(i-1,i-1) = 0.0
     hi(i,i-1) = sr / norm

     do j = i, en
        yr = hr(i-1,j)
        yi = hi(i-1,j)
        zzr = hr(i,j)
        zzi = hi(i,j)
        hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
        hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
        hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
        hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
    end do

  end do

  si = hi(en,en)
  if (si == 0.0) go to 540
  norm = pythag(hr(en,en),si)
  sr = hr(en,en) / norm
  si = si / norm
  hr(en,en) = norm
  hi(en,en) = 0.0
!
!  Inverse operation (columns)
!
  540 continue

  do j = lp1, en

     xr = wr(j-1)
     xi = wi(j-1)

     do i = l, j

        yr = hr(i,j-1)
        yi = 0.0
        zzr = hr(i,j)
        zzi = hi(i,j)
        if ( i /= j ) then
          yi = hi(i,j-1)
          hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
        end if
        hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
        hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
        hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

     end do

  end do

  if (si == 0.0) go to 240

  do i = l, en
     yr = hr(i,en)
     yi = hi(i,en)
     hr(i,en) = sr * yr - si * yi
     hi(i,en) = sr * yi + si * yr
  end do

  go to 240
!
!  A root found
!
  660 wr(en) = hr(en,en) + tr
  wi(en) = hi(en,en) + ti
  en = enm1
  go to 220
!
!  set error: all eigenvalues have not converged after 30*n iterations
!
 1000 ierr = en
  return
end
subroutine comqr2 ( nm, n, low, igh, ortr, orti, hr, hi, wr, wi, zr, zi, ierr )
!
!*******************************************************************************
!
!! COMQR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of a unitary analogue of the
!     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
!     and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!     the unitary analogue substitutes the qr algorithm of francis
!     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a complex upper Hessenberg matrix by the qr
!     method.  the eigenvectors of a complex general matrix
!     can also be found if  corth  has been used to reduce
!     this general matrix to Hessenberg form.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        ortr and orti contain information about the unitary trans-
!          formations used in the reduction by  corth, if performed.
!          only elements low through igh are used.  if the eigenvectors
!          of the Hessenberg matrix are desired, set ortr(j) and
!          orti(j) to 0.0 for these elements.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper Hessenberg matrix.
!          their lower triangles below the subdiagonal contain further
!          information about the transformations which were used in the
!          reduction by  corth, if performed.  if the eigenvectors of
!          the Hessenberg matrix are desired, these elements may be
!          arbitrary.
!
!     on output
!
!        ortr, orti, and the upper Hessenberg portions of hr and hi
!          have been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors.  the eigenvectors
!          are unnormalized.  if an error exit is made, none of
!          the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
  integer igh
  integer n
  integer nm
!
  integer en
  integer enm1
  real hi(nm,n)
  real hr(nm,n)
  integer i
  integer iend
  integer ierr
  integer ii
  integer itn
  integer its
  integer j
  integer jj
  integer k
  integer l
  integer ll
  integer low
  integer lp1
  integer m
  integer nn
  real norm
  real orti(igh)
  real ortr(igh)
  real pythag
  real si
  real sr
  real ti
  real tr
  real tst1
  real tst2
  real wi(n)
  real wr(n)
  real xi
  real xr
  real yi
  real yr
  real zi(nm,n)
  real zr(nm,n)
  real zzi 
  real zzr
!
  ierr = 0
!
!  Initialize eigenvector matrix
!
  call rmat_ident ( nm, n, zr )

  zi(1:n,1:n) = 0.0
!
!  Form the matrix of accumulated transformations from the information 
!  left by corth
!
  iend = igh - low - 1
  if (iend) 180, 150, 105

  105 continue

  do ii = 1, iend

     i = igh - ii
     if (ortr(i) == 0.0 .and. orti(i) == 0.0) go to 140
     if (hr(i,i-1) == 0.0 .and. hi(i,i-1) == 0.0) go to 140
!
!  Norm below is negative of h formed in corth
!
     norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)

     do k = i+1, igh
       ortr(k) = hr(k,i-1)
       orti(k) = hi(k,i-1)
     end do

     do j = i, igh

        sr = 0.0
        si = 0.0

        do k = i, igh
           sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
           si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
        end do

        sr = sr / norm
        si = si / norm

        do k = i, igh
           zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
           zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
        end do

      end do

  140 continue

  end do
!
!  Create real subdiagonal elements
!
  150 continue

  l = low + 1

  do i = l, igh

     ll = min(i+1,igh)

     if (hi(i,i-1) == 0.0) go to 170

     norm = pythag(hr(i,i-1),hi(i,i-1))
     yr = hr(i,i-1) / norm
     yi = hi(i,i-1) / norm
     hr(i,i-1) = norm
     hi(i,i-1) = 0.0

     do j = i, n
        si = yr * hi(i,j) - yi * hr(i,j)
        hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
        hi(i,j) = si
     end do

     do j = 1, ll
        si = yr * hi(j,i) + yi * hr(j,i)
        hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
        hi(j,i) = si
     end do

     do j = low, igh
        si = yr * zi(j,i) + yi * zr(j,i)
        zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
        zi(j,i) = si
     end do

  170 continue

  end do
!
!  Store roots isolated by cbal
!
  180 continue

  do i = 1, n
    if ( i < low .or. i > igh) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    end if
  end do

  en = igh
  tr = 0.0
  ti = 0.0
  itn = 30*n
!
!  Search for next eigenvalue
!
  220 if (en < low) go to 680
  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element
!
  240 continue

  do ll = low, en
     l = en + low - ll
     if (l == low) go to 300
     tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
     tst2 = tst1 + abs(hr(l,l-1))
     if (tst2 == tst1) go to 300
  end do
!
!  Form shift
!
  300 if (l == en) go to 660
  if (itn == 0) go to 1000
  if (its == 10 .or. its == 20) go to 320
  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1)
  xi = hi(enm1,en) * hr(en,enm1)
  if (xr == 0.0 .and. xi == 0.0) go to 340
  yr = (hr(enm1,enm1) - sr) / 2.0
  yi = (hi(enm1,enm1) - si) / 2.0
  call csroot(yr**2-yi**2+xr,2.0*yr*yi+xi,zzr,zzi)
  if (yr * zzr + yi * zzi >= 0.0) go to 310
  zzr = -zzr
  zzi = -zzi

  310 continue

  call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift
!
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
  si = 0.0

  340 continue

  do i = low, en
     hr(i,i) = hr(i,i) - sr
     hi(i,i) = hi(i,i) - si
  end do

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Reduce to triangle (rows)
!
  lp1 = l + 1

  do i = lp1, en

     sr = hr(i,i-1)
     hr(i,i-1) = 0.0
     norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
     xr = hr(i-1,i-1) / norm
     wr(i-1) = xr
     xi = hi(i-1,i-1) / norm
     wi(i-1) = xi
     hr(i-1,i-1) = norm
     hi(i-1,i-1) = 0.0
     hi(i,i-1) = sr / norm

     do j = i, n
        yr = hr(i-1,j)
        yi = hi(i-1,j)
        zzr = hr(i,j)
        zzi = hi(i,j)
        hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
        hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
        hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
        hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
     end do

  end do

  si = hi(en,en)

  if ( si /= 0.0 ) then

    norm = pythag(hr(en,en),si)
    sr = hr(en,en) / norm
    si = si / norm
    hr(en,en) = norm
    hi(en,en) = 0.0

    do j = en+1, n
      yr = hr(en,j)
      yi = hi(en,j)
      hr(en,j) = sr * yr + si * yi
      hi(en,j) = sr * yi - si * yr
    end do

  end if
!
!  Inverse operation (columns)
!
  do j = lp1, en

     xr = wr(j-1)
     xi = wi(j-1)

     do i = 1, j

        yr = hr(i,j-1)
        yi = 0.0
        zzr = hr(i,j)
        zzi = hi(i,j)

        if ( i /= j ) then
          yi = hi(i,j-1)
          hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
        end if

        hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
        hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
        hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

     end do

     do i = low, igh
        yr = zr(i,j-1)
        yi = zi(i,j-1)
        zzr = zr(i,j)
        zzi = zi(i,j)
        zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
        zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
        zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
        zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
     end do

  end do

  if ( si /= 0.0 ) then

    do i = 1, en
      yr = hr(i,en)
      yi = hi(i,en)
      hr(i,en) = sr * yr - si * yi
      hi(i,en) = sr * yi + si * yr
    end do

    do i = low, igh
      yr = zr(i,en)
      yi = zi(i,en)
      zr(i,en) = sr * yr - si * yi
      zi(i,en) = sr * yi + si * yr
    end do

  end if

  go to 240
!
!  A root found
!
  660 continue

  hr(en,en) = hr(en,en) + tr
  wr(en) = hr(en,en)
  hi(en,en) = hi(en,en) + ti
  wi(en) = hi(en,en)
  en = enm1
  go to 220
!
!  All roots found.  
!  Backsubstitute to find vectors of upper triangular form
!
  680 norm = 0.0

  do i = 1, n
    do j = i, n
      tr = abs(hr(i,j)) + abs(hi(i,j))
      if (tr > norm) norm = tr
    end do
  end do

  if ( n == 1 ) then
    return
  end if

  if ( norm == 0.0 ) then
    return
  end if

  do nn = 2, n

     en = n + 2 - nn
     xr = wr(en)
     xi = wi(en)
     hr(en,en) = 1.0
     hi(en,en) = 0.0
     enm1 = en - 1

     do ii = 1, enm1

        i = en - ii
        zzr = 0.0
        zzi = 0.0

        do j = i+1, en
           zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
           zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
        end do

        yr = xr - wr(i)
        yi = xi - wi(i)

        if ( yr == 0.0 .and. yi == 0.0 ) then
           tst1 = norm
           yr = tst1
  760      continue
           yr = 0.01 * yr
           tst2 = norm + yr
           if (tst2 > tst1) go to 760
        end if

        call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
!
!  Overflow control
!
        tr = abs(hr(i,en)) + abs(hi(i,en))

        if ( tr /= 0.0 ) then

          tst1 = tr
          tst2 = tst1 + 1.0/tst1

          if ( tst2 <= tst1 ) then

            do j = i, en
              hr(j,en) = hr(j,en)/tr
              hi(j,en) = hi(j,en)/tr
            end do

          end if

       end if

     end do

  end do
!
!  End backsubstitution
!
  enm1 = n - 1
!
!  Vectors of isolated roots
!
  do i = 1, n - 1

    if (i < low .or. i > igh ) then

      do j = i+1, n
        zr(i,j) = hr(i,j)
        zi(i,j) = hi(i,j)
      end do

    end if

  end do
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  do jj = low, n - 1

     j = n + low - jj
     m = min(j,igh)

     do i = low, igh

        zzr = 0.0
        zzi = 0.0
        do k = low, m
           zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
           zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
        end do

        zr(i,j) = zzr
        zi(i,j) = zzi

      end do

  end do

  return
!
!  Set error: all eigenvalues have not converged after 30*n iterations
!
 1000 ierr = en
  return
end
subroutine cortb ( nm, low, igh, ar, ai, ortr, orti, m, zr, zi )
!
!*******************************************************************************
!
!! CORTB determines eigenvectors by undoing the CORTH transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure ortbak, num. math. 12, 349-368(1968)
!     by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     upper Hessenberg matrix determined by  corth.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        ar and ai contain information about the unitary
!          transformations used in the reduction by  corth
!          in their strict lower triangles.
!
!        ortr and orti contain further information about the
!          transformations used in the reduction by  corth.
!          only elements low through igh are used.
!
!        m is the number of columns of zr and zi to be back transformed.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!        ortr and orti have been altered.
!
  integer igh
  integer m
  integer nm
!
  real ai(nm,igh)
  real ar(nm,igh)
  real gi
  real gr
  real h
  integer i
  integer j
  integer la
  integer low
  integer mm
  integer mp
  real orti(igh)
  real ortr(igh)
  real zi(nm,m)
  real zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  la = igh - 1

  if ( igh - 1 < low + 1 ) then
    return
  end if

  do mm = low + 1, la

    mp = low + igh - mm

    if ( ar(mp,mp-1) /= 0.0 .or. ai(mp,mp-1) /= 0.0 ) then

      h = ar(mp,mp-1) * ortr(mp) + ai(mp,mp-1) * orti(mp)

      do i = mp+1, igh
        ortr(i) = ar(i,mp-1)
        orti(i) = ai(i,mp-1)
      end do

      do j = 1, m

        gr = 0.0
        gi = 0.0

        do i = mp, igh
           gr = gr + ortr(i) * zr(i,j) + orti(i) * zi(i,j)
           gi = gi + ortr(i) * zi(i,j) - orti(i) * zr(i,j)
        end do

        gr = gr / h
        gi = gi / h

        do i = mp, igh
          zr(i,j) = zr(i,j) + gr * ortr(i) - gi * orti(i)
          zi(i,j) = zi(i,j) + gr * orti(i) + gi * ortr(i)
        end do

      end do

    end if

  end do

  return
end
subroutine corth ( nm, n, low, igh, ar, ai, ortr, orti )
!
!*******************************************************************************
!
!! CORTH transforms a complex general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure orthes, num. math. 12, 349-368(1968)
!     by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a complex general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper Hessenberg form by
!     unitary similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex input matrix.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the Hessenberg matrix.  information
!          about the unitary transformations used in the reduction
!          is stored in the remaining triangles under the
!          Hessenberg matrix.
!
!        ortr and orti contain further information about the
!          transformations.  only elements low through igh are used.
!
  integer igh
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  real f
  real fi
  real fr
  real g
  real h
  integer i
  integer ii
  integer j
  integer jj
  integer la
  integer m,mp,low
  real orti(igh)
  real ortr(igh)
  real pythag
  real scale
!
  la = igh - 1

  if ( igh - 1 < low + 1 ) then
    return
  end if

  do m = low + 1, la

     h = 0.0
     ortr(m) = 0.0
     orti(m) = 0.0
     scale = 0.0
!
!  Scale column.
!
     do i = m, igh
       scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
     end do

     if (scale == 0.0) go to 180
     mp = m + igh

     do ii = m, igh
        i = mp - ii
        ortr(i) = ar(i,m-1) / scale
        orti(i) = ai(i,m-1) / scale
        h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
     end do

     g = sqrt(h)
     f = pythag(ortr(m),orti(m))

     if (f == 0.0) go to 103
     h = h + f * g
     g = g / f
     ortr(m) = (1.0 + g) * ortr(m)
     orti(m) = (1.0 + g) * orti(m)
     go to 105

  103    ortr(m) = g
     ar(m,m-1) = scale
!
!  Form (i-(u*ut)/h) * a
!
  105    continue

     do j = m, n

        fr = 0.0
        fi = 0.0

        do ii = m, igh
           i = mp - ii
           fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
           fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
        end do

        fr = fr / h
        fi = fi / h

        do i = m, igh
           ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
           ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
        end do

    end do
!
!  Form (i-(u*ut)/h)*a*(i-(u*ut)/h)
!
     do i = 1, igh

        fr = 0.0
        fi = 0.0

        do jj = m, igh
           j = mp - jj
           fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
           fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
        end do

        fr = fr / h
        fi = fi / h

        do j = m, igh
           ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
           ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
        end do

     end do

     ortr(m) = scale * ortr(m)
     orti(m) = scale * orti(m)
     ar(m,m-1) = -g * ar(m,m-1)
     ai(m,m-1) = -g * ai(m,m-1)

  180 continue

  end do

  return
end
subroutine csroot ( xr, xi, yr, yi )
!
!*******************************************************************************
!
!! CSROOT computes the complex square root of a complex quantity.
!
!
!  Discussion:
!
!     branch chosen so that yr >= 0.0 and sign(yi) == sign(xi)
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
  real pythag
  real s
  real ti
  real tr
  real xi
  real xr
  real yi
  real yr
!
  tr = xr
  ti = xi
  s = sqrt ( 0.5 * ( pythag(tr,ti) + abs(tr)) )
  if (tr >= 0.0) yr = s
  if (ti < 0.0) s = -s
  if (tr <= 0.0) yi = s
  if (tr < 0.0) yr = 0.5*(ti/yi)
  if (tr > 0.0) yi = 0.5*(ti/yr)

  return
end
subroutine elmbak ( nm, low, igh, a, int, m, z )
!
!*******************************************************************************
!
!! ELMBAK determines eigenvectors by undoing the ELMHES transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure elmbak,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     upper Hessenberg matrix determined by  elmhes.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        a contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
  integer igh
  integer m
  integer nm
!
  real a(nm,igh)
  integer i
  integer int(igh)
  integer j
  integer la
  integer low
  integer mm
  integer mp
  real x
  real z(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  la = igh - 1

  if ( la < low + 1 ) then
    return
  end if

  do mm = low + 1, la

     mp = low + igh - mm

     do i = mp+1, igh

       x = a(i,mp-1)
       if (x /= 0.0 ) then
         do j = 1, m
           z(i,j) = z(i,j) + x * z(mp,j)
         end do
       end if

     end do

     i = int(mp)

     if (i /= mp ) then

       do j = 1, m
         call r_swap ( z(i,j), z(mp,j) )
       end do

     end if

  end do

  return
end
subroutine elmhes ( nm, n, low, igh, a, int )
!
!*******************************************************************************
!
!! ELMHES transforms a real general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure elmhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper Hessenberg form by
!     stabilized elementary similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the Hessenberg matrix.  the multipliers
!          which were used in the reduction are stored in the
!          remaining triangle under the Hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
  integer igh
  integer n
  integer nm
!
  real a(nm,n)
  integer i
  integer int(igh)
  integer j
  integer la
  integer low
  integer m
  integer mm1
  real x
  real y
!
  la = igh - 1

  do m = low + 1, la

    mm1 = m - 1
    x = 0.0
    i = m

    do j = m, igh
      if ( abs ( a(j,mm1) ) > abs ( x ) ) then
        x = a(j,mm1)
        i = j
      end if
    end do

    int(m) = i
!
!  Interchange rows and columns of the matrix.
!
    if ( i /= m ) then

      do j = mm1, n
        call r_swap ( a(i,j), a(m,j) )
      end do

      do j = 1, igh
        call r_swap ( a(j,i), a(j,m) )
      end do

    end if

    if ( x /= 0.0 ) then

      do i = m+1, igh

        y = a(i,mm1)

        if ( y /= 0.0 ) then

          y = y / x
          a(i,mm1) = y

          do j = m, n
            a(i,j) = a(i,j) - y * a(m,j)
          end do

          do j = 1, igh
            a(j,m) = a(j,m) + y * a(j,i)
          end do

        end if

      end do

    end if

  end do

  return
end
subroutine eltran ( nm, n, low, igh, a, int, z )
!
!*******************************************************************************
!
!! ELTRAN accumulates similarity transformations used by ELMHES.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure elmtrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the stabilized elementary
!     similarity transformations used in the reduction of a
!     real general matrix to upper Hessenberg form by  elmhes.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!     on output
!
!        z contains the transformation matrix produced in the
!          reduction by  elmhes.
!
  integer igh
  integer n
  integer nm
!
  real a(nm,igh)
  integer i
  integer int(igh)
  integer j
  integer kl
  integer low
  integer mm
  integer mp
  real z(nm,n)
!
!  Initialize Z to the identity matrix.
!
  call rmat_ident ( nm, n, z )

  kl = igh - low - 1

  if ( kl < 1 ) then
    return
  end if

  do mm = 1, kl

     mp = igh - mm

     do i = mp+1, igh
       z(i,mp) = a(i,mp-1)
     end do

     i = int(mp)

     if ( i /= mp ) then

       do j = mp, igh
         z(mp,j) = z(i,j)
         z(i,j) = 0.0
       end do

       z(i,mp) = 1.0

     end if

  end do

  return
end
subroutine figi ( nm, n, t, d, e, e2, ierr )
!
!*******************************************************************************
!
!! FIGI transforms a real nonsymmetric tridiagonal matrix to symmetric form.
!
!
!  Discussion:
!
!     given a nonsymmetric tridiagonal matrix such that the products
!     of corresponding pairs of off-diagonal elements are all
!     non-negative, this subroutine reduces it to a symmetric
!     tridiagonal matrix with the same eigenvalues.  if, further,
!     a zero product only occurs when both factors are zero,
!     the reduced matrix is similar to the original matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        t contains the input matrix.  its subdiagonal is
!          stored in the last n-1 positions of the first column,
!          its diagonal in the n positions of the second column,
!          and its superdiagonal in the first n-1 positions of
!          the third column.  t(1,1) and t(n,3) are arbitrary.
!
!     on output
!
!        t is unaltered.
!
!        d contains the diagonal elements of the symmetric matrix.
!
!        e contains the subdiagonal elements of the symmetric
!          matrix in its last n-1 positions.  e(1) is not set.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        ierr is set to
!          zero       for normal return,
!          n+i        if t(i,1)*t(i-1,3) is negative,
!          -(3*n+i)   if t(i,1)*t(i-1,3) is zero with one factor
!                     non-zero.  in this case, the eigenvectors of
!                     the symmetric matrix are not simply related
!                     to those of  t  and should not be sought.
!
  integer n
  integer nm
!
  real d(n)
  real e(n)
  real e2(n)
  integer i
  integer ierr
  real t(nm,3)
!
  ierr = 0

  do i = 1, n

    if ( i >= 1 ) then

      e2(i) = t(i,1) * t(i-1,3)

      if ( e2(i) < 0.0 ) then

        ierr = n + i
        return

      else if ( e2(i) == 0.0 ) then

        if ( t(i,1) /= 0.0 .or. t(i-1,3) /= 0.0 ) then
          ierr = - 3 * n - i
          return
        end if

        e(i) = 0.0

      else

        e(i) = sqrt ( e2(i) )

      end if

    end if

    d(i) = t(i,2)

  end do

  return
end
subroutine figi2 ( nm, n, t, d, e, z, ierr )
!
!*******************************************************************************
!
!! FIGI2 transforms a real nonsymmetric tridiagonal matrix to symmetric form.
!
!
!  Discussion:
!
!     given a nonsymmetric tridiagonal matrix such that the products
!     of corresponding pairs of off-diagonal elements are all
!     non-negative, and zero only when both factors are zero, this
!     subroutine reduces it to a symmetric tridiagonal matrix
!     using and accumulating diagonal similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        t contains the input matrix.  its subdiagonal is
!          stored in the last n-1 positions of the first column,
!          its diagonal in the n positions of the second column,
!          and its superdiagonal in the first n-1 positions of
!          the third column.  t(1,1) and t(n,3) are arbitrary.
!
!     on output
!
!        t is unaltered.
!
!        d contains the diagonal elements of the symmetric matrix.
!
!        e contains the subdiagonal elements of the symmetric
!          matrix in its last n-1 positions.  e(1) is not set.
!
!        z contains the transformation matrix produced in
!          the reduction.
!
!        ierr is set to
!          zero       for normal return,
!          n+i        if t(i,1)*t(i-1,3) is negative,
!          2*n+i      if t(i,1)*t(i-1,3) is zero with
!                     one factor non-zero.
!
  integer n
  integer nm
!
  real d(n)
  real e(n)
  real h
  integer i
  integer ierr
  integer j
  real t(nm,3)
  real z(nm,n)
!
  ierr = 0

  do i = 1, n

    z(i,1:n) = 0.0

    if ( i == 1 ) then

      z(i,i) = 1.0

    else

      h = t(i,1) * t(i-1,3)

      if ( h < 0 ) then

        ierr = n + i
        return

      else if ( h == 0 ) then

        if ( t(i,1) /= 0.0 .or. t(i-1,3) /= 0.0 ) then
          ierr = 2 * n + i
          return
        end if

        e(i) = 0.0
        z(i,i) = 1.0

      else

        e(i) = sqrt ( h )
        z(i,i) = z(i-1,i-1) * e(i) / t(i-1,3)

      end if

    end if

    d(i) = t(i,2)

  end do

  return
end
subroutine hqr ( nm, n, low, igh, h, wr, wi, ierr )
!
!*******************************************************************************
!
!! HQR computes all eigenvalues of a real upper Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     this subroutine finds the eigenvalues of a real
!     upper Hessenberg matrix by the qr method.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper Hessenberg matrix.  information about
!          the transformations used in the reduction to Hessenberg
!          form by  elmhes  or  orthes, if performed, is stored
!          in the remaining triangle under the Hessenberg matrix.
!
!     on output
!
!        h has been destroyed.  therefore, it must be saved
!          before calling  hqr  if subsequent calculation and
!          back transformation of eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
  integer n
  integer nm
!
  integer en
  integer enm2
  real h(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer na
  real norm
  logical notlas
  real p
  real q
  real r
  real s
  real t
  real tst1
  real tst2
  real w
  real wi(n)
  real wr(n)
  real x
  real y
  real zz
!
  ierr = 0
  norm = 0.0
  k = 1
!
!  Store roots isolated by balanc and compute matrix norm.
!
  do i = 1, n

     do j = k, n
       norm = norm + abs ( h(i,j) )
     end do

     k = i
     if (i < low .or. i > igh) then
       wr(i) = h(i,i)
       wi(i) = 0.0
     end if

  end do

  en = igh
  t = 0.0
  itn = 30*n
!
!  Search for next eigenvalues
!
   60 continue

  if (en < low) then
    return
  end if

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for single small sub-diagonal element
!
   70 continue

  do ll = low, en
     l = en + low - ll
     if (l == low) go to 100
     s = abs(h(l-1,l-1)) + abs(h(l,l))
     if (s == 0.0) s = norm
     tst1 = s
     tst2 = tst1 + abs(h(l,l-1))
     if (tst2 == tst1) go to 100
  end do
!
!  Form shift
!
  100 continue

  x = h(en,en)
  if (l == en) go to 270
  y = h(na,na)
  w = h(en,na) * h(na,en)
  if (l == na) go to 280

  if (itn == 0) then
    ierr = en
    return
  end if

  if (its /= 10 .and. its /= 20) go to 130
!
!  Form exceptional shift
!
  t = t + x

  do i = low, en
    h(i,i) = h(i,i) - x
  end do

  s = abs(h(en,na)) + abs(h(na,enm2))
  x = 0.75 * s
  y = x
  w = -0.4375 * s * s
  130 its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  do mm = l, enm2
     m = enm2 + l - mm
     zz = h(m,m)
     r = x - zz
     s = y - zz
     p = (r * s - w) / h(m+1,m) + h(m,m+1)
     q = h(m+1,m+1) - zz - r - s
     r = h(m+2,m+1)
     s = abs(p) + abs(q) + abs(r)
     p = p / s
     q = q / s
     r = r / s
     if (m == l) go to 150
     tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
     tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
     if (tst2 == tst1) go to 150
  end do

  150 continue

  do i = m+2, en
     h(i,i-2) = 0.0
     if ( i /= m+2 ) then
       h(i,i-3) = 0.0
     end if
  end do
!
!  Double qr step involving rows l to en and columns m to en
!
  do k = m, na

     notlas = k /= na
     if (k == m) go to 170
     p = h(k,k-1)
     q = h(k+1,k-1)
     r = 0.0
     if (notlas) r = h(k+2,k-1)
     x = abs(p) + abs(q) + abs(r)
     if (x == 0.0) go to 260
     p = p / x
     q = q / x
     r = r / x

170  continue

     s = sign(sqrt(p*p+q*q+r*r),p)
     if (k == m) go to 180
     h(k,k-1) = -s * x

     go to 190

180  continue

     if (l /= m) h(k,k-1) = -h(k,k-1)

190  continue

     p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p
     if (notlas) go to 225
!
!  Row modification
!
     do j = k, n
        p = h(k,j) + q * h(k+1,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
     end do

     j = min(en,k+3)
!
!  Column modification
!
     do i = 1, j
        p = x * h(i,k) + y * h(i,k+1)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
     end do

     go to 255

  225    continue
!
!  Row modification
!
     do j = k, n
        p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
        h(k+2,j) = h(k+2,j) - p * zz
     end do

     j = min(en,k+3)
!
!  Column modification
!
     do i = 1, j
        p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
        h(i,k+2) = h(i,k+2) - p * r
     end do

  255    continue

  260 continue

  end do

  go to 70
!
!  One root found
!
270 continue

  wr(en) = x + t
  wi(en) = 0.0
  en = na
  go to 60
!
!  Two roots found
!
  280 continue

  p = (y - x) / 2.0
  q = p * p + w
  zz = sqrt(abs(q))
  x = x + t
  if (q < 0.0) go to 320
!
!  Real pair
!
  zz = p + sign(zz,p)
  wr(na) = x + zz
  wr(en) = wr(na)
  if (zz /= 0.0) wr(en) = x - w / zz
  wi(na) = 0.0
  wi(en) = 0.0
  go to 330
!
!  Complex pair
!
  320 wr(na) = x + p
  wr(en) = x + p
  wi(na) = zz
  wi(en) = -zz
  330 en = enm2
  go to 60
end
subroutine hqr2 ( nm, n, low, igh, h, wr, wi, z, ierr )
!
!*******************************************************************************
!
!! HQR2 computes eigenvalues and eigenvectors of a real upper Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper Hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to Hessenberg form
!     and to accumulate the similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper Hessenberg matrix.
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the Hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output
!
!        h has been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
  integer nm
  integer n
!
  integer en
  integer enm2
  real h(nm,n)
  real hnorm
  integer i
  integer ierr
  integer igh
  integer ii
  integer itn
  integer its
  integer j
  integer jj
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer na
  integer nn
  real norm
  logical notlas
  real p
  real q
  real r
  real ra
  real s
  real sa
  real t
  real temp
  real tst1
  real tst2
  real vi
  real vr
  real w
  real wi(n)
  real wr(n)
  real x
  real y
  real z(nm,n)
  real zz
!
  ierr = 0
  norm = 0.0
  k = 1
!
!  Store roots isolated by BALANC and compute the matrix norm.
!
  do i = 1, n

     do j = k, n
       norm = norm + abs ( h(i,j) )
     end do

     k = i
     if (i < low .or. i > igh) then
       wr(i) = h(i,i)
       wi(i) = 0.0
     end if

  end do

  en = igh
  t = 0.0
  itn = 30*n
!
!  Search for next eigenvalues.
!
   60 continue

  if (en < low) go to 340
  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for single small sub-diagonal element.
!
   70 continue

  do ll = low, en
     l = en + low - ll
     if (l == low) go to 100
     s = abs(h(l-1,l-1)) + abs(h(l,l))
     if (s == 0.0) s = norm
     tst1 = s
     tst2 = tst1 + abs(h(l,l-1))
     if (tst2 == tst1) go to 100
  end do
!
!  Form shift
!
  100 continue

  x = h(en,en)
  if (l == en) go to 270
  y = h(na,na)
  w = h(en,na) * h(na,en)
  if (l == na) go to 280

  if (itn == 0) then
    ierr = en
    return
  end if

  if (its /= 10 .and. its /= 20) go to 130
!
!  Form exceptional shift.
!
  t = t + x

  do i = low, en
    h(i,i) = h(i,i) - x
  end do

  s = abs(h(en,na)) + abs(h(na,enm2))
  x = 0.75 * s
  y = x
  w = -0.4375 * s * s

  130 continue

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  do mm = l, enm2
     m = enm2 + l - mm
     zz = h(m,m)
     r = x - zz
     s = y - zz
     p = (r * s - w) / h(m+1,m) + h(m,m+1)
     q = h(m+1,m+1) - zz - r - s
     r = h(m+2,m+1)
     s = abs(p) + abs(q) + abs(r)
     p = p / s
     q = q / s
     r = r / s
     if (m == l) go to 150
     tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
     tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
     if (tst2 == tst1) go to 150
  end do

  150 continue

  do i = m+2, en
     h(i,i-2) = 0.0
     if (i /= m+2) then
       h(i,i-3) = 0.0
     end if
  end do
!
!  Double QR step involving rows L to EN and columns M to EN.
!
  do k = m, na

     notlas = k /= na

     if ( k /= m ) then

       p = h(k,k-1)
       q = h(k+1,k-1)
       r = 0.0
       if (notlas) r = h(k+2,k-1)
       x = abs(p) + abs(q) + abs(r)
       if (x == 0.0) go to 260
       p = p / x
       q = q / x
       r = r / x

     end if

     s = sign(sqrt(p*p+q*q+r*r),p)

     if (k /= m) then
       h(k,k-1) = -s * x
     else if (l /= m) then
       h(k,k-1) = -h(k,k-1)
     end if

     p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p
     if (notlas) go to 225
!
!  Row modification
!
     do j = k, n
        p = h(k,j) + q * h(k+1,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
     end do

     j = min(en,k+3)
!
!  Column modification
!
     do i = 1, j
        p = x * h(i,k) + y * h(i,k+1)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
     end do
!
!  Accumulate transformations
!
     do i = low, igh
        p = x * z(i,k) + y * z(i,k+1)
        z(i,k) = z(i,k) - p
        z(i,k+1) = z(i,k+1) - p * q
     end do

     go to 255
  225    continue
!
!  Row modification
!
     do j = k, n
        p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
        h(k+2,j) = h(k+2,j) - p * zz
     end do

     j = min(en,k+3)
!
!  Column modification
!
     do i = 1, j
        p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
        h(i,k+2) = h(i,k+2) - p * r
     end do
!
!  Accumulate transformations
!
     do i = low, igh
        p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
        z(i,k) = z(i,k) - p
        z(i,k+1) = z(i,k+1) - p * q
        z(i,k+2) = z(i,k+2) - p * r
     end do

  255    continue

  260 continue

  end do

  go to 70
!
!  One root found
!
  270 continue
  h(en,en) = x + t
  wr(en) = h(en,en)
  wi(en) = 0.0
  en = na
  go to 60
!
!  Two roots found
!
  280 continue
  p = (y - x) / 2.0
  q = p * p + w
  zz = sqrt(abs(q))
  h(en,en) = x + t
  x = h(en,en)
  h(na,na) = y + t
  if (q < 0.0) go to 320
!
!  Real pair
!
  zz = p + sign(zz,p)
  wr(na) = x + zz
  wr(en) = wr(na)
  if (zz /= 0.0) wr(en) = x - w / zz
  wi(na) = 0.0
  wi(en) = 0.0
  x = h(en,na)
  s = abs(x) + abs(zz)
  p = x / s
  q = zz / s
  r = sqrt(p*p+q*q)
  p = p / r
  q = q / r
!
!  Row modification
!
  do j = na, n
     zz = h(na,j)
     h(na,j) = q * zz + p * h(en,j)
     h(en,j) = q * h(en,j) - p * zz
  end do
!
!  Column modification
!
  do i = 1, en
     zz = h(i,na)
     h(i,na) = q * zz + p * h(i,en)
     h(i,en) = q * h(i,en) - p * zz
  end do
!
!  Accumulate transformations
!
  do i = low, igh
     zz = z(i,na)
     z(i,na) = q * zz + p * z(i,en)
     z(i,en) = q * z(i,en) - p * zz
  end do

  go to 330
!
!  Complex pair
!
  320 continue
  wr(na) = x + p
  wr(en) = x + p
  wi(na) = zz
  wi(en) = -zz
  330 continue
  en = enm2
  go to 60
!
!  All roots found.  
!  Backsubstitute to find vectors of upper triangular form.
!
  340 continue

  if ( norm == 0.0 ) then
    return
  end if

  do nn = 1, n

     en = n + 1 - nn
     p = wr(en)
     q = wi(en)
     na = en - 1
     if (q) 710, 600, 800
!
!  Real vector
!
  600    m = en
     h(en,en) = 1.0
     if (na == 0) go to 800

     do ii = 1, na

        i = en - ii
        w = h(i,i) - p
        r = 0.0

        do j = m, en
          r = r + h(i,j) * h(j,en)
        end do

        if (wi(i) >= 0.0) go to 630

        zz = w
        s = r
        go to 700
  630       m = i
        if (wi(i) /= 0.0) go to 640
        t = w
        if (t /= 0.0) go to 635
           tst1 = norm
           t = tst1
  632          t = 0.01 * t
           tst2 = norm + t
           if (tst2 > tst1) go to 632
  635       h(i,en) = -r / t
        go to 680
!
!  Solve real equations
!
  640       x = h(i,i+1)
        y = h(i+1,i)
        q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
        t = (x * s - zz * r) / q
        h(i,en) = t
        if (abs(x) <= abs(zz)) go to 650
        h(i+1,en) = (-r - w * t) / x
        go to 680
  650       h(i+1,en) = (-s - y * t) / zz
!
!  Overflow control.
!
  680       t = abs(h(i,en))
        if (t == 0.0) go to 700
        tst1 = t
        tst2 = tst1 + 1.0/tst1

        if ( tst2 <= tst1 ) then
          h(i:en,en) = h(i:en,en) / t
        end if

  700    continue
    end do
!
!  End real vector
!
     go to 800
!
!  Complex vector
!
  710    m = na
!
!  Last vector component chosen imaginary so that eigenvector matrix is triangular
!
     if (abs(h(en,na)) <= abs(h(na,en))) go to 720
     h(na,na) = q / h(en,na)
     h(na,en) = -(h(en,en) - p) / h(en,na)
     go to 730
  720    call cdiv(0.0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0
     h(en,en) = 1.0
     enm2 = na - 1

     do ii = 1, enm2

        i = na - ii
        w = h(i,i) - p
        ra = 0.0
        sa = 0.0

        do j = m, en
           ra = ra + h(i,j) * h(j,na)
           sa = sa + h(i,j) * h(j,en)
        end do

        if (wi(i) < 0.0) then
          zz = w
          r = ra
          s = sa
        end if

         m = i

        if (wi(i) == 0.0) then
          call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
          go to 790
        end if
!
!  Solve complex equations
!
        x = h(i,i+1)
        y = h(i+1,i)
        vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
        vi = (wr(i) - p) * 2.0 * q

        if (vr == 0.0 .and. vi == 0.0) then
           tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
           vr = tst1
  783      continue
           vr = 0.01 * vr
           tst2 = tst1 + vr
           if (tst2 > tst1) then
             go to 783
           end if
        end if

        call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,h(i,na),h(i,en))

        if (abs(x) > abs(zz) + abs(q)) then
          h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
          h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
        else
           call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,h(i+1,na),h(i+1,en))
        end if
!
!  Overflow control
!
  790   continue

        t = max(abs(h(i,na)), abs(h(i,en)))

        if (t /= 0.0) then
          tst1 = t
          tst2 = tst1 + 1.0/tst1
          if (tst2 <= tst1) then
            h(i:en,na) = h(i:en,na) / t
            h(i:en,en) = h(i:en,en) / t
          end if
        end if

  795   continue

      end do
!
!  End complex vector
!
  800 continue

  end do
!
!  End back substitution.
!
!  Vectors of isolated roots.
!
  do i = 1, n

    if ( i < low .or. i > igh ) then
      z(i,i:n) = h(i,i:n)
    end if

  end do
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  do jj = low, n

     j = n + low - jj
     m = min(j,igh)

     do i = low, igh

        zz = 0.0
        do k = low, m
          zz = zz + z(i,k) * h(k,j)
        end do
        z(i,j) = zz

      end do

  end do

  return
end
subroutine htrib3 ( nm, n, a, tau, m, zr, zi )
!
!*******************************************************************************
!
!! HTRIB3 determines eigenvectors by undoing the HTRID3 transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure trbak3, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a complex hermitian
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by  htrid3.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains information about the unitary transformations
!          used in the reduction by  htrid3.
!
!        tau contains further information about the transformations.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
  integer m
  integer n
  integer nm
!
  real a(nm,n)
  real h
  integer i
  integer j
  integer k
  integer l
  real s
  real si
  real tau(2,n)
  real zi(nm,m)
  real zr(nm,m)
!
  if ( m == 0 ) then
    return
  end if
!
!  Transform the eigenvectors of the real symmetric tridiagonal matrix 
!  to those of the hermitian tridiagonal matrix.
!
  do k = 1, n
    do j = 1, m
      zi(k,j) = -zr(k,j) * tau(2,k)
      zr(k,j) = zr(k,j) * tau(1,k)
    end do
  end do
!
!  Recover and apply the Householder matrices.
!
  do i = 2, n

    l = i - 1
    h = a(i,i)

    if ( h /= 0.0 ) then

      do j = 1, m

        s = 0.0
        si = 0.0

        do k = 1, l
          s = s + a(i,k) * zr(k,j) - a(k,i) * zi(k,j)
          si = si + a(i,k) * zi(k,j) + a(k,i) * zr(k,j)
        end do

        s = ( s / h ) / h
        si = ( si / h ) / h

        do k = 1, l
          zr(k,j) = zr(k,j) - s * a(i,k) - si * a(k,i)
          zi(k,j) = zi(k,j) - si * a(i,k) + s * a(k,i)
        end do

      end do

    end if

  end do

  return
end
subroutine htribk ( nm, n, ar, ai, tau, m, zr, zi )
!
!*******************************************************************************
!
!! HTRIBK determines eigenvectors by undoing the HTRIDI transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure trbak1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a complex hermitian
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by  htridi.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction by  htridi  in their
!          full lower triangles except for the diagonal of ar.
!
!        tau contains further information about the transformations.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
  integer m
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  real h
  integer i
  integer j
  integer k
  integer l
  real s
  real si
  real tau(2,n)
  real zi(nm,m)
  real zr(nm,m)
!
  if (m == 0) then
    return
  end if
!
!  Transform the eigenvectors of the real symmetric tridiagonal matrix to 
!  those of the hermitian tridiagonal matrix.
!
  do k = 1, n
    do j = 1, m
      zi(k,j) = -zr(k,j) * tau(2,k)
      zr(k,j) = zr(k,j) * tau(1,k)
    end do
  end do
!
!  Recover and apply the Householder matrices
!
  do i = 2, n

    l = i - 1
    h = ai(i,i)

    if ( h /= 0.0 ) then

      do j = 1, m

        s = 0.0
        si = 0.0
        do k = 1, l
          s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
          si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
        end do

        s = ( s / h ) / h
        si = ( si / h ) / h

        do k = 1, l
          zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
          zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
        end do

      end do

    end if

  end do

  return
end
subroutine htrid3 ( nm, n, a, d, e, e2, tau )
!
!*******************************************************************************
!
!! HTRID3 tridiagonalizes a complex hermitian packed matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure tred3, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a complex hermitian matrix, stored as
!     a single square array, to a real symmetric tridiagonal matrix
!     using unitary similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the lower triangle of the complex hermitian input
!          matrix.  the real parts of the matrix elements are stored
!          in the full lower triangle of a, and the imaginary parts
!          are stored in the transposed positions of the strict upper
!          triangle of a.  no storage is required for the zero
!          imaginary parts of the diagonal elements.
!
!     on output
!
!        a contains information about the unitary transformations
!          used in the reduction.
!
!        d contains the diagonal elements of the the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        tau contains further information about the transformations.
!
  integer n
  integer nm
!
  real a(nm,n)
  real d(n)
  real e(n)
  real e2(n)
  real f
  real fi
  real g
  real gi
  real h
  real hh
  integer i
  integer ii
  integer j
  integer k
  integer l
  real pythag
  real scale
  real si
  real tau(2,n)
!
  tau(1,n) = 1.0
  tau(2,n) = 0.0

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0
    scale = 0.0

    if ( l < 1 ) then
      e(i) = 0.0
      e2(i) = 0.0
      go to 290
    end if
!
!  Scale row.
!
     do k = 1, l
       scale = scale + abs(a(i,k)) + abs(a(k,i))
     end do

     if ( scale == 0.0 ) then
       tau(1,l) = 1.0
       tau(2,l) = 0.0
       e(i) = 0.0
       e2(i) = 0.0
       go to 290
     end if
   
      do k = 1, l
        a(i,k) = a(i,k) / scale
        a(k,i) = a(k,i) / scale
        h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i)
     end do

     e2(i) = scale * scale * h
     g = sqrt(h)
     e(i) = scale * g
     f = pythag(a(i,l),a(l,i))
!
!  Form next diagonal element of matrix t
!
     if ( f /= 0.0 ) then

       tau(1,l) = (a(l,i) * tau(2,i) - a(i,l) * tau(1,i)) / f
       si = (a(i,l) * tau(2,i) + a(l,i) * tau(1,i)) / f
       h = h + f * g
       g = 1.0 + g / f
       a(i,l) = g * a(i,l)
       a(l,i) = g * a(l,i)

       if (l == 1) go to 270

     else

       tau(1,l) = -tau(1,i)
       si = tau(2,i)
       a(i,l) = g

     end if

     f = 0.0

     do j = 1, l

        g = 0.0
        gi = 0.0
!
!  Form element of A*U.
!
        do k = 1, j-1
           g = g + a(j,k) * a(i,k) + a(k,j) * a(k,i)
           gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k)
        end do

        g = g + a(j,j) * a(i,j)
        gi = gi - a(j,j) * a(j,i)

        do k = j+1, l
           g = g + a(k,j) * a(i,k) - a(j,k) * a(k,i)
           gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k)
        end do
!
!  Form element of P.
!
        e(j) = g / h
        tau(2,j) = gi / h
        f = f + e(j) * a(i,j) - tau(2,j) * a(j,i)

     end do

     hh = f / (h + h)
!
!  Form reduced A.
!
     do j = 1, l

        f = a(i,j)
        g = e(j) - hh * f
        e(j) = g
        fi = -a(j,i)
        gi = tau(2,j) - hh * fi
        tau(2,j) = -gi
        a(j,j) = a(j,j) - 2.0 * (f * g + fi * gi)

        do k = 1, j-1
           a(j,k) = a(j,k) - f * e(k) - g * a(i,k) + fi * tau(2,k) + gi * a(k,i)
           a(k,j) = a(k,j) - f * tau(2,k) - g * a(k,i) - fi * e(k) - gi * a(i,k)
        end do

     end do

  270   continue

     a(i,1:l) = scale * a(i,1:l)
     a(1:l,i) = scale * a(1:l,i)
     tau(2,l) = -si

  290 continue

     d(i) = a(i,i)
     a(i,i) = scale * sqrt(h)

  end do

  return
end
subroutine htridi ( nm, n, ar, ai, d, e, e2, tau )
!
!*******************************************************************************
!
!! HTRIDI tridiagonalizes a complex hermitian matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure tred1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a complex hermitian matrix
!     to a real symmetric tridiagonal matrix using
!     unitary similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex hermitian input matrix.
!          only the lower triangle of the matrix need be supplied.
!
!     on output
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction in their full lower
!          triangles.  their strict upper triangles and the
!          diagonal of ar are unaltered.
!
!        d contains the diagonal elements of the the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        tau contains further information about the transformations.
!
  integer n
  integer nm
!
  real ai(nm,n)
  real ar(nm,n)
  real d(n)
  real e(n)
  real e2(n)
  real f
  real fi
  real g
  real gi
  real h
  real hh
  integer i
  integer ii
  integer j
  integer k
  integer l
  real pythag
  real scale
  real si
  real tau(2,n)
!
  tau(1,n) = 1.0
  tau(2,n) = 0.0

  do i = 1, n
    d(i) = ar(i,i)
  end do

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0
    scale = 0.0

    if ( l < 1 ) then
      e(i) = 0.0
      e2(i) = 0.0
      go to 290
    end if
!
!  Scale row.
!
    do k = 1, l
      scale = scale + abs(ar(i,k)) + abs(ai(i,k))
    end do

    if ( scale == 0.0 ) then
      tau(1,l) = 1.0
      tau(2,l) = 0.0
      e(i) = 0.0
      e2(i) = 0.0
      go to 290
    end if

    do k = 1, l
        ar(i,k) = ar(i,k) / scale
        ai(i,k) = ai(i,k) / scale
        h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
     end do

     e2(i) = scale * scale * h
     g = sqrt(h)
     e(i) = scale * g
     f = pythag(ar(i,l),ai(i,l))
!
!  Form next diagonal element of matrix T.
!
     if (f == 0.0) go to 160
     tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
     si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
     h = h + f * g
     g = 1.0 + g / f
     ar(i,l) = g * ar(i,l)
     ai(i,l) = g * ai(i,l)
     if (l == 1) go to 270
     go to 170
  160    tau(1,l) = -tau(1,i)
     si = tau(2,i)
     ar(i,l) = g
  170    f = 0.0

     do j = 1, l

        g = 0.0
        gi = 0.0
!
!  Form element of A*U.
!
        do k = 1, j
           g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
           gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
        end do

        do k = j+1, l
           g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
           gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
        end do
!
!  Form element of P.
!
  220   e(j) = g / h
        tau(2,j) = gi / h
        f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)

     end do

     hh = f / (h + h)
!
!  Form reduced A.
!
     do j = 1, l

        f = ar(i,j)
        g = e(j) - hh * f
        e(j) = g
        fi = -ai(i,j)
        gi = tau(2,j) - hh * fi
        tau(2,j) = -gi

        do k = 1, j
          ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) + fi * tau(2,k) &
            + gi * ai(i,k)
          ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) - fi * e(k) &
            - gi * ar(i,k)
        end do

    end do

  270   do k = 1, l
          ar(i,k) = scale * ar(i,k)
          ai(i,k) = scale * ai(i,k)
        end do

     tau(2,l) = -si

290  continue

     hh = d(i)
     d(i) = ar(i,i)
     ar(i,i) = hh
     ai(i,i) = scale * sqrt(h)
  
  end do

  return
end
subroutine imtql1 ( n, d, e, ierr )
!
!*******************************************************************************
!
!! IMTQL1 computes all eigenvalues of a symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure imtql1,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the implicit ql method.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
  integer n
!
  real b
  real c
  real d(n)
  real e(n)
  real f
  real g
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer m
  integer mml
  real p
  real pythag
  real r
  real s
  real tst1
  real tst2
!
  ierr = 0
  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do
  e(n) = 0.0

  do l = 1, n

     j = 0
!
!  Look for small sub-diagonal element
!
  105    continue

    do m = l, n

      if (m == n) then
        exit
      end if

      tst1 = abs(d(m)) + abs(d(m+1))
      tst2 = tst1 + abs(e(m))

      if (tst2 == tst1) then
        exit
      end if

    end do

     p = d(l)
     if (m == l) go to 215

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift
!
     g = (d(l+1) - p) / (2.0 * e(l))
     r = pythag(g,1.0)
     g = d(m) - p + e(l) / (g + sign(r,g))
     s = 1.0
     c = 1.0
     p = 0.0
     mml = m - l

     do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)
        r = pythag(f,g)
        e(i+1) = r
        if (r == 0.0) go to 210
        s = f / r
        c = g / r
        g = d(i+1) - p
        r = (d(i) - g) * s + 2.0 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b

     end do

     d(l) = d(l) - p
     e(l) = g
     e(m) = 0.0
     go to 105
!
!  Recover from underflow
!
210  continue

     d(i+1) = d(i+1) - p
     e(m) = 0.0
     go to 105
!
!  Order eigenvalues
!
  215  continue

     do ii = 2, l
        i = l + 2 - ii
        if (p >= d(i-1)) go to 270
        d(i) = d(i-1)
     end do

     i = 1
  270    d(i) = p

  end do

  return
end
subroutine imtql2 ( nm, n, d, e, z, ierr )
!
!*******************************************************************************
!
!! IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure IMTQL2,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the implicit ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
  integer n
  integer nm
!
  real b
  real c
  real d(n)
  real e(n)
  real f
  real g
  integer i
  integer ierr
  integer ii
  integer j
  integer k
  integer l
  integer m
  integer mml
  real p
  real pythag
  real r
  real s
  real tst1
  real tst2
  real z(nm,n)
!
  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do
  e(n) = 0.0

  do l = 1, n

     j = 0
!
!  Look for small sub-diagonal element
!
  105    continue

     do m = l, n

       if (m == n) then
         exit
       end if

       tst1 = abs(d(m)) + abs(d(m+1))
       tst2 = tst1 + abs(e(m))

       if (tst2 == tst1) then
         exit
       end if

     end do

     p = d(l)
     if (m == l) go to 240

     if (j >= 30) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift
!
     g = (d(l+1) - p) / (2.0 * e(l))
     r = pythag(g,1.0)
     g = d(m) - p + e(l) / (g + sign(r,g))
     s = 1.0
     c = 1.0
     p = 0.0
     mml = m - l

     do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)
        r = pythag(f,g)
        e(i+1) = r
        if (r == 0.0) go to 210
        s = f / r
        c = g / r
        g = d(i+1) - p
        r = (d(i) - g) * s + 2.0 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
!
!  Form vector
!
        do k = 1, n
           f = z(k,i+1)
           z(k,i+1) = s * z(k,i) + c * f
           z(k,i) = c * z(k,i) - s * f
        end do

     end do

     d(l) = d(l) - p
     e(l) = g
     e(m) = 0.0
     go to 105
!
!  Recover from underflow
!
  210    d(i+1) = d(i+1) - p
     e(m) = 0.0
     go to 105
  240 continue

  end do
!
!  Order eigenvalues and eigenvectors
!
  do ii = 2, n

     i = ii - 1
     k = i
     p = d(i)

     do j = ii, n
        if (d(j) < p) then
          k = j
          p = d(j)
        end if
     end do

     if (k == i) go to 300
     d(k) = d(i)
     d(i) = p

     do j = 1, n
       call r_swap ( z(j,i), z(j,k) )
     end do

  300 continue

  end do

  return
end
subroutine imtqlv ( n, d, e, e2, w, ind, ierr, rv1 )
!
!*******************************************************************************
!
!! IMTQLV computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine is a variant of  imtql1  which is a translation of
!     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
!     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues of a symmetric tridiagonal
!     matrix by the implicit ql method and associates with them
!     their corresponding submatrix indices.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!     on output
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        w contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        ind contains the submatrix indices associated with the
!          corresponding eigenvalues in w: 1 for eigenvalues
!          belonging to the first submatrix from the top,
!          2 for those belonging to the second submatrix, etc..
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
  integer n
!
  real b
  real c
  real d(n)
  real e(n)
  real e2(n)
  real f
  real g
  integer i
  integer ierr
  integer ii
  integer ind(n)
  integer j
  integer k
  integer l
  integer m
  integer mml
  real p
  real pythag
  real r
  real rv1(n)
  real s
  integer tag
  real tst1
  real tst2
  real w(n)
!
  ierr = 0
  k = 0
  tag = 0
  w(1:n) = d(1:n)
  e2(1) = 0.0
  rv1(1:n-1) = e(2:n)
  rv1(n) = 0.0

  do l = 1, n

     j = 0
!
!  Look for small sub-diagonal element
!
  105    continue

     do m = l, n
        if (m == n) go to 120
        tst1 = abs(w(m)) + abs(w(m+1))
        tst2 = tst1 + abs(rv1(m))
        if (tst2 == tst1) go to 120
!
!  Guard against underflowed element of e2
!
        if (e2(m+1) == 0.0) go to 125
     end do

120  continue

     if (m <= k) go to 130
     if (m /= n) e2(m+1) = 0.0

125  continue

     k = m
     tag = tag + 1

130  continue

    p = w(l)
     if (m == l) go to 215

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift
!
     g = (w(l+1) - p) / (2.0 * rv1(l))
     r = pythag(g,1.0)
     g = w(m) - p + rv1(l) / (g + sign(r,g))
     s = 1.0
     c = 1.0
     p = 0.0
     mml = m - l

     do ii = 1, mml
        i = m - ii
        f = s * rv1(i)
        b = c * rv1(i)
        r = pythag(f,g)
        rv1(i+1) = r
        if (r == 0.0) go to 210
        s = f / r
        c = g / r
        g = w(i+1) - p
        r = (w(i) - g) * s + 2.0 * c * b
        p = s * r
        w(i+1) = g + p
        g = c * r - b
     end do

     w(l) = w(l) - p
     rv1(l) = g
     rv1(m) = 0.0
     go to 105
!
!  Recover from underflow
!
  210    w(i+1) = w(i+1) - p
     rv1(m) = 0.0
     go to 105
!
!  Order eigenvalues
!
  215    continue

     do ii = 2, l
        i = l + 2 - ii
        if (p >= w(i-1)) go to 270
        w(i) = w(i-1)
        ind(i) = ind(i-1)
     end do

     i = 1
  270    w(i) = p
     ind(i) = tag

  end do

  return
end
subroutine invit ( nm, n, a, wr, wi, select, mm, m, z, ierr, rm1, rv1, rv2 )
!
!*******************************************************************************
!
!! INVIT computes eigenvectors given eigenvalues, for a real upper Hessenberg matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure invit
!     by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvectors of a real upper
!     Hessenberg matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the Hessenberg matrix.
!
!        wr and wi contain the real and imaginary parts, respectively,
!          of the eigenvalues of the matrix.  the eigenvalues must be
!          stored in a manner identical to that of subroutine  hqr,
!          which recognizes possible splitting of the matrix.
!
!        select specifies the eigenvectors to be found. the
!          eigenvector corresponding to the j-th eigenvalue is
!          specified by setting select(j) to .true..
!
!        mm should be set to an upper bound for the number of
!          columns required to store the eigenvectors to be found.
!          note that two columns are required to store the
!          eigenvector corresponding to a complex eigenvalue.
!
!     on output
!
!        a and wi are unaltered.
!
!        wr may have been altered since close eigenvalues are perturbed
!          slightly in searching for independent eigenvectors.
!
!        select may have been altered.  if the elements corresponding
!          to a pair of conjugate complex eigenvalues were each
!          initially set to .true., the program resets the second of
!          the two elements to .false..
!
!        m is the number of columns actually used to store
!          the eigenvectors.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the next selected eigenvalue is real, the next column
!          of z contains its eigenvector.  if the eigenvalue is
!          complex, the next two columns of z contain the real and
!          imaginary parts of its eigenvector.  the eigenvectors are
!          normalized so that the component of largest magnitude is 1.
!          any vector which fails the acceptance test is set to zero.
!
!        ierr is set to
!          zero       for normal return,
!          -(2*n+1)   if more than mm columns of z are necessary
!                     to store the eigenvectors corresponding to
!                     the specified eigenvalues.
!          -k         if the iteration corresponding to the k-th
!                     value fails,
!          -(n+k)     if both error situations occur.
!
!        rm1, rv1, and rv2 are temporary storage arrays.  note that rm1
!          is square of dimension n by n and, augmented by two columns
!          of z, is the transpose of the corresponding algol b array.
!
  integer nm
  integer n
!
  real a(nm,n)
  real eps3
  real growto
  integer i
  integer ierr
  integer ii
  real ilambd
  integer ip
  integer its
  integer j
  integer k
  integer km1
  integer l
  integer m
  integer mm
  integer mp
  integer n1
  real norm
  real normv
  integer ns
  real pythag
  real rlambd
  real rm1(n,n)
  real rv1(n)
  real rv2(n)
  integer s
  logical select(n)
  real t
  integer uk
  real ukroot
  real w
  real wi(n)
  real wr(n)
  real x
  real y
  real z(nm,mm)
!
  ierr = 0
  uk = 0
  s = 1
!
!  ip = 0, real eigenvalue
!       1, first of conjugate complex pair
!      -1, second of conjugate complex pair
!
  ip = 0
  n1 = n - 1

  do k = 1, n

     if (wi(k) == 0.0 .or. ip < 0) go to 100
     ip = 1
     if (select(k) .and. select(k+1)) select(k+1) = .false.
  100    if (.not. select(k)) go to 960
     if (wi(k) /= 0.0) s = s + 1
     if (s > mm) go to 1000
     if (uk >= k) go to 200
!
!  Check for possible splitting
!
     do uk = k, n
        if (uk == n) go to 140
        if (a(uk+1,uk) == 0.0) go to 140
     end do
!
!  Compute infinity norm of leading uk by uk (Hessenberg) matrix
!
  140    norm = 0.0
     mp = 1

     do i = 1, uk

        x = 0.0
        do j = mp, uk
          x = x + abs(a(i,j))
        end do

        if (x > norm) norm = x
        mp = i

     end do
!
!  EPS3 replaces zero pivot in decomposition and close roots are modified 
!  by eps3
!
     if (norm == 0.0) norm = 1.0
     eps3 = abs ( norm ) * epsilon ( 1.0 )
!
!  GROWTO is the criterion for the growth
!
     ukroot = uk
     ukroot = sqrt(ukroot)
     growto = 0.1 / ukroot
  200    rlambd = wr(k)
     ilambd = wi(k)
     if (k == 1) go to 280
     km1 = k - 1
     go to 240
!
!  Perturb eigenvalue if it is close to any previous eigenvalue
!
  220    rlambd = rlambd + eps3

  240    continue

     do ii = 1, km1
        i = k - ii
        if (select(i) .and. abs(wr(i)-rlambd) < eps3 .and. &
            abs(wi(i)-ilambd) < eps3) then
          go to 220
        end if
     end do

     wr(k) = rlambd
!
!  Perturb conjugate eigenvalue to match
!
     wr(k+ip) = rlambd
!
!  Form upper Hessenberg a-rlambd*i (transposed) and initial real vector
!
  280    mp = 1

     do i = 1, uk

        do j = mp, uk
          rm1(j,i) = a(i,j)
        end do

        rm1(i,i) = rm1(i,i) - rlambd
        mp = i
        rv1(i) = eps3

     end do

     its = 0
     if (ilambd /= 0.0) go to 520
!
!  Real eigenvalue.
!  Triangular decomposition with interchanges, replacing zero pivots by eps3
!
     do i = 2, uk

        mp = i - 1

        if ( abs ( rm1(mp,i) ) > abs ( rm1(mp,mp) ) ) then

          do j = mp, uk
            call r_swap ( rm1(j,i), rm1(j,mp) )
          end do

        end if

        if (rm1(mp,mp) == 0.0) then
          rm1(mp,mp) = eps3
        end if

        x = rm1(mp,i) / rm1(mp,mp)

        if ( x /= 0.0 ) then
          rm1(i:uk,i) = rm1(i:uk,i) - x * rm1(i:uk,mp)
        end if

      end do

      if (rm1(uk,uk) == 0.0) then
        rm1(uk,uk) = eps3
      end if
!
!  Back substitution for real vector
!
  440 continue

      do ii = 1, uk

        i = uk + 1 - ii
        y = rv1(i)

        do j = i+1, uk
          y = y - rm1(j,i) * rv1(j)
        end do

        rv1(i) = y / rm1(i,i)

     end do

     go to 740
!
!  Complex eigenvalue.
!  triangular decomposition with interchanges,
!  replacing zero pivots by eps3.  store imaginary
!  parts in upper triangle starting at (1,3)
!
520  continue

     ns = n - s
     z(1,s-1) = -ilambd
     z(1,s) = 0.0
     if (n == 2) go to 550
     rm1(1,3) = -ilambd
     z(1,s-1) = 0.0

     rm1(1,4:n) = 0.0

  550 continue

     do i = 2, uk

        mp = i - 1
        w = rm1(mp,i)
        if (i < n) t = rm1(mp,i+1)
        if (i == n) t = z(mp,s-1)
        x = rm1(mp,mp) * rm1(mp,mp) + t * t
        if (w * w <= x) go to 580
        x = rm1(mp,mp) / w
        y = t / w
        rm1(mp,mp) = w
        if (i < n) rm1(mp,i+1) = 0.0
        if (i == n) z(mp,s-1) = 0.0

        do j = i, uk

           w = rm1(j,i)
           rm1(j,i) = rm1(j,mp) - x * w
           rm1(j,mp) = w
           if (j < n1) go to 555
           l = j - ns
           z(i,l) = z(mp,l) - y * w
           z(mp,l) = 0.0
           go to 560
  555          rm1(i,j+2) = rm1(mp,j+2) - y * w
           rm1(mp,j+2) = 0.0
  560       continue

        end do

        rm1(i,i) = rm1(i,i) - y * ilambd
        if (i < n1) go to 570
        l = i - ns
        z(mp,l) = -ilambd
        z(i,l) = z(i,l) + x * ilambd
        go to 640
  570       rm1(mp,i+2) = -ilambd
        rm1(i,i+2) = rm1(i,i+2) + x * ilambd
        go to 640
  580       if (x /= 0.0) go to 600
        rm1(mp,mp) = eps3
        if (i < n) rm1(mp,i+1) = 0.0
        if (i == n) z(mp,s-1) = 0.0
        t = 0.0
        x = eps3 * eps3
  600       w = w / x
        x = rm1(mp,mp) * w
        y = -t * w

        do j = i, uk
           if (j < n1) go to 610
           l = j - ns
           t = z(mp,l)
           z(i,l) = -x * t - y * rm1(j,mp)
           go to 615
  610          t = rm1(mp,j+2)
           rm1(i,j+2) = -x * t - y * rm1(j,mp)
  615          rm1(j,i) = rm1(j,i) - x * rm1(j,mp) + y * t
        end do

        if (i < n1) go to 630
        l = i - ns
        z(i,l) = z(i,l) - ilambd
        go to 640
  630       rm1(i,i+2) = rm1(i,i+2) - ilambd
  640    continue

     end do

     if (uk < n1) go to 650
     l = uk - ns
     t = z(uk,l)
     go to 655
  650    t = rm1(uk,uk+2)
  655    if (rm1(uk,uk) == 0.0 .and. t == 0.0) rm1(uk,uk) = eps3
!
!  Back substitution for complex vector
!
  660   continue

     do ii = 1, uk

        i = uk + 1 - ii
        x = rv1(i)
        y = 0.0

        do j = i+1, uk

           if ( j >= n1 ) then
             t = z(i,j-ns)
           else
             t = rm1(i,j+2)
           end if

           x = x - rm1(j,i) * rv1(j) + t * rv2(j)
           y = y - rm1(j,i) * rv2(j) - t * rv1(j)

        end do

        if ( i >= n1 ) then
          t = z(i,i-ns)
        else
          t = rm1(i,i+2)
        end if

       call cdiv(x,y,rm1(i,i),t,rv1(i),rv2(i))

     end do
!
!  Acceptance test for real or complex eigenvector and normalization
!
  740    its = its + 1
     norm = 0.0
     normv = 0.0

     do i = 1, uk
       if (ilambd == 0.0) x = abs(rv1(i))
       if (ilambd /= 0.0) x = pythag(rv1(i),rv2(i))
       if ( normv < x)  then
         normv = x
         j = i
       end if
       norm = norm + x
     end do

     if (norm < growto) go to 840
!
!  Accept vector
!
     x = rv1(j)
     if (ilambd == 0.0) x = 1.0 / x
     if (ilambd /= 0.0) y = rv2(j)

     do i = 1, uk
        if (ilambd /= 0.0) go to 800
        z(i,s) = rv1(i) * x
        go to 820
  800       call cdiv(rv1(i),rv2(i),x,y,z(i,s-1),z(i,s))
  820    continue
     end do

     if (uk == n) go to 940
     j = uk + 1
     go to 900
!
!  Choose a new starting vector.
!
  840    if (its >= uk) go to 880
     x = ukroot
     y = eps3 / (x + 1.0)

     rv1(1) = eps3
     rv1(2:uk) = y
     
     j = uk - its + 1
     rv1(j) = rv1(j) - eps3 * x
     if (ilambd == 0.0) go to 440
     go to 660
!
!  Set error: unaccepted eigenvector
!
  880    j = 1
     ierr = -k
!
!  Set remaining vector components to zero
!
  900    continue

     do i = j, n
        z(i,s) = 0.0
        if (ilambd /= 0.0) z(i,s-1) = 0.0
     end do

  940    s = s + 1
  960    if (ip == (-1)) ip = 0
     if (ip == 1) ip = -1

  end do

  go to 1001
!
!  Set error: underestimate of eigenvector space required
!
 1000 if (ierr /= 0) ierr = ierr - n
  if (ierr == 0) ierr = -(2 * n + 1)
 1001 m = s - 1 - iabs(ip)
  return
end
subroutine minfit ( nm, m, n, a, w, ip, b, ierr )
!
!*******************************************************************************
!
!! MINFIT solves the least squares problem, for a real overdetermined linear system.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure minfit,
!     num. math. 14, 403-420(1970) by golub and reinsch.
!     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
!
!     this subroutine determines, towards the solution of the linear
!                                                        t
!     system ax=b, the singular value decomposition a=usv  of a real
!                                         t
!     m by n rectangular matrix, forming u b rather than u.  Householder
!     bidiagonalization and a variant of the qr algorithm are used.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.  note that nm must be at least
!          as large as the maximum of m and n.
!
!        m is the number of rows of a and b.
!
!        n is the number of columns of a and the order of v.
!
!        a contains the rectangular coefficient matrix of the system.
!
!        ip is the number of columns of b.  ip can be zero.
!
!        b contains the constant column matrix of the system
!          if ip is not zero.  otherwise b is not referenced.
!
!     on output
!
!        a has been overwritten by the matrix v (orthogonal) of the
!          decomposition in its first n rows and columns.  if an
!          error exit is made, the columns of v corresponding to
!          indices of correct singular values should be correct.
!
!        w contains the n (non-negative) singular values of a (the
!          diagonal elements of s).  they are unordered.  if an
!          error exit is made, the singular values should be correct
!          for indices ierr+1,ierr+2,...,n.
!
!                                   t
!        b has been overwritten by u b.  if an error exit is made,
!                       t
!          the rows of u b corresponding to indices of correct
!          singular values should be correct.
!
!        ierr is set to
!          zero       for normal return,
!          k          if the k-th singular value has not been
!                     determined after 30 iterations.
!
  integer ip
  integer n
  integer nm
!
  real a(nm,n)
  real b(nm,ip)
  real c
  real f
  real g
  real h
  integer i
  integer i1
  integer ierr
  integer ii
  integer its
  integer j
  integer k
  integer k1
  integer kk
  integer l
  integer l1
  integer ll
  integer m
  integer m1
  real pythag
  real rv1(n)
  real s
  real scale
  real tst1
  real tst2
  real w(n)
  real x
  real y
  real z
!
  ierr = 0
!
!  Householder reduction to bidiagonal form
!
  g = 0.0
  scale = 0.0
  x = 0.0

  do i = 1, n

    l = i + 1
    rv1(i) = scale * g
    g = 0.0
    s = 0.0
    scale = 0.0

    if ( i <= m ) then

      do k = i, m
        scale = scale + abs(a(k,i))
      end do

      if ( scale /= 0.0 ) then

        do k = i, m
          a(k,i) = a(k,i) / scale
          s = s + a(k,i)**2
        end do

        f = a(i,i)
        g = - sign ( sqrt(s), f )
        h = f * g - s
        a(i,i) = f - g

        do j = l, n

          s = 0.0
          do k = i, m
            s = s + a(k,i) * a(k,j)
          end do

          f = s / h
          a(i:m,j) = a(i:m,j) + f * a(i:m,i)

        end do

        do j = 1, ip

          s = 0.0
          do k = i, m
            s = s + a(k,i) * b(k,j)
          end do

          b(i:m,j) = b(i:m,j) + s * a(i:m,i) / h

        end do

        a(i:m,i) = scale * a(i:m,i)

      end if

    end if

    w(i) = scale * g
    g = 0.0
    s = 0.0
    scale = 0.0

    if ( i <= m .and. i /= n ) then

      do k = l, n
        scale = scale + abs(a(i,k))
      end do

      if ( scale /= 0.0 ) then

        a(i,l:n) = a(i,l:n) / scale

        do k = l, n
          s = s + a(i,k)**2
        end do

        f = a(i,l)
        g = -sign(sqrt(s),f)
        h = f * g - s
        a(i,l) = f - g
        rv1(l:n) = a(i,l:n) / h

        do j = l, m

          s = 0.0
          do k = l, n
            s = s + a(j,k) * a(i,k)
          end do

          a(j,l:n) = a(j,l:n) + s * rv1(l:n)

        end do

        a(i,l:n) = scale * a(i,l:n)

      end if

    end if

    x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

  end do
!
!  Accumulation of right-hand transformations.
!
  do ii = 1, n

    i = n + 1 - ii

    if ( i /= n ) then

      if ( g /= 0.0 ) then

        a(l:n,i) = ( a(i,l:n) / a(i,l) ) / g

        do j = l, n

          s = 0.0
          do k = l, n
            s = s + a(i,k) * a(k,j)
          end do

          a(l:n,j) = a(l:n,j) + s * a(l:n,i)

        end do

      end if

      a(i,l:n) = 0.0
      a(l:n,i) = 0.0

    end if

    a(i,i) = 1.0
    g = rv1(i)
    l = i

  end do

  if ( m < n .and. ip /= 0 ) then
    m1 = m + 1
    b(m+1:n,1:ip) = 0.0
  end if
!
!  Diagonalization of the bidiagonal form.
!
  tst1 = x

  do kk = 1, n

     k1 = n - kk
     k = k1 + 1
     its = 0
!
!  Test for splitting.
!
  520   continue

     do ll = 1, k

        l1 = k - ll
        l = l1 + 1
        tst2 = tst1 + abs(rv1(l))

        if (tst2 == tst1) go to 565

        tst2 = tst1 + abs(w(k-ll))

        if ( tst2 == tst1 ) then
          exit
        end if

     end do
!
!  Cancellation of rv1(l) if l greater than 1
!
  540    continue

     c = 0.0
     s = 1.0

     do i = l, k

        f = s * rv1(i)
        rv1(i) = c * rv1(i)
        tst2 = tst1 + abs(f)

        if (tst2 == tst1) then
          exit
        end if

        g = w(i)
        h = pythag(f,g)
        w(i) = h
        c = g / h
        s = -f / h

        do j = 1, ip
           y = b(l1,j)
           z = b(i,j)
           b(l1,j) = y * c + z * s
           b(i,j) = -y * s + z * c
        end do

    end do
!
!  Test for convergence
!
  565    z = w(k)

     if (l == k) go to 650
!
!  Shift from bottom 2 by 2 minor
!
     if ( its >= 30 ) then
       ierr = k
       return
     end if

     its = its + 1
     x = w(l)
     y = w(k1)
     g = rv1(k1)
     h = rv1(k)
     f = 0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
     g = pythag(f,1.0)
     f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h)
!
!  Next qr transformation
!
     c = 1.0
     s = 1.0

     do i1 = l, k1

        i = i1 + 1
        g = rv1(i)
        y = w(i)
        h = s * g
        g = c * g
        z = pythag(f,h)
        rv1(i1) = z
        c = f / z
        s = h / z
        f = x * c + g * s
        g = -x * s + g * c
        h = y * s
        y = y * c

        do j = 1, n
           x = a(j,i1)
           z = a(j,i)
           a(j,i1) = x * c + z * s
           a(j,i) = -x * s + z * c
        end do

        z = pythag(f,h)
        w(i1) = z

        if ( z /= 0.0 ) then
          c = f / z
          s = h / z
        end if

        f = c * g + s * y
        x = -s * g + c * y

        do j = 1, ip
           y = b(i1,j)
           z = b(i,j)
           b(i1,j) = y * c + z * s
           b(i,j) = -y * s + z * c
        end do

     end do

     rv1(l) = 0.0
     rv1(k) = f
     w(k) = x
     go to 520
!
!  Convergence
!
  650 continue

    if ( z < 0.0 ) then
      w(k) = - z
      a(1:n,k) = - a(1:n,k)
    end if

  end do

  return
end
subroutine ortbak ( nm, low, igh, a, ort, m, z )
!
!*******************************************************************************
!
!! ORTBAK determines eigenvectors by undoing the ORTHES transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure ortbak,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     upper Hessenberg matrix determined by  orthes.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  orthes
!          in its strict lower triangle.
!
!        ort contains further information about the trans-
!          formations used in the reduction by  orthes.
!          only elements low through igh are used.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!        ort has been altered.
!
  integer igh
  integer m
  integer nm
!
  real a(nm,igh)
  real g
  integer i
  integer j
  integer low
  integer mp
  real ort(igh)
  real z(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  do mp = igh - 1, low + 1, -1

    if ( a(mp,mp-1) /= 0.0 ) then

      do i = mp+1, igh
        ort(i) = a(i,mp-1)
      end do

      do j = 1, m

        g = 0.0
        do i = mp, igh
          g = g + ort(i) * z(i,j)
        end do

        g = ( g / ort(mp) ) / a(mp,mp-1)

        do i = mp, igh
          z(i,j) = z(i,j) + g * ort(i)
        end do

      end do

    end if

  end do

  return
end
subroutine orthes ( nm, n, low, igh, a, ort )
!
!*******************************************************************************
!
!! ORTHES transforms a real general matrix to upper Hessenberg form.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure orthes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper Hessenberg form by
!     orthogonal similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the Hessenberg matrix.  information about
!          the orthogonal transformations used in the reduction
!          is stored in the remaining triangle under the
!          Hessenberg matrix.
!
!        ort contains further information about the transformations.
!          only elements low through igh are used.
!
  integer igh
  integer n
  integer nm
!
  real a(nm,n)
  real f
  real g
  real h
  integer i
  integer ii
  integer j
  integer jj
  integer la
  integer low
  integer m
  integer mp
  real ort(igh)
  real scale
!
  la = igh - 1

  do m = low + 1, la

     h = 0.0
     ort(m) = 0.0
     scale = 0.0
!
!  Scale column.
!
     do i = m, igh
       scale = scale + abs(a(i,m-1))
     end do

     if (scale /= 0.0 ) then

     mp = m + igh

     do ii = m, igh
        i = mp - ii
        ort(i) = a(i,m-1) / scale
        h = h + ort(i) * ort(i)
     end do

     g = -sign(sqrt(h),ort(m))
     h = h - ort(m) * g
     ort(m) = ort(m) - g
!
!  Form (i-(u*ut)/h) * a
!
     do j = m, n

        f = 0.0

        do ii = m, igh
           i = mp - ii
           f = f + ort(i) * a(i,j)
        end do

        f = f / h

        do i = m, igh
          a(i,j) = a(i,j) - f * ort(i)
        end do

     end do
!
!  Form (i-(u*ut)/h)*a*(i-(u*ut)/h)
!
     do i = 1, igh

        f = 0.0
        do jj = m, igh
          j = mp - jj
          f = f + ort(j) * a(i,j)
        end do

        f = f / h

        do j = m, igh
          a(i,j) = a(i,j) - f * ort(j)
        end do

     end do

     ort(m) = scale * ort(m)
     a(m,m-1) = scale * g

    end if

  end do

  return
end
subroutine ortran ( nm, n, low, igh, a, ort, z )
!
!*******************************************************************************
!
!! ORTRAN accumulates similarity transformations generated by ORTHES.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure ortrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the orthogonal similarity
!     transformations used in the reduction of a real general
!     matrix to upper Hessenberg form by  orthes.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  orthes
!          in its strict lower triangle.
!
!        ort contains further information about the trans-
!          formations used in the reduction by  orthes.
!          only elements low through igh are used.
!
!     on output
!
!        z contains the transformation matrix produced in the
!          reduction by  orthes.
!
!        ort has been altered.
!
  integer igh
  integer n
  integer nm
!
  real a(nm,igh)
  real g
  integer i
  integer j
  integer kl
  integer low
  integer mm
  integer mp
  real ort(igh)
  real z(nm,n)
!
!  Initialize Z to the identity matrix.
!
  call rmat_ident ( nm, n, z )

  kl = igh - low - 1

  if (kl < 1) then
    return
  end if

  do mm = 1, kl

    mp = igh - mm

    if ( a(mp,mp-1) /= 0.0 ) then

      ort(mp+1:igh) = a(mp+1:igh,mp-1)

      do j = mp, igh

        g = 0.0
        do i = mp, igh
          g = g + ort(i) * z(i,j)
        end do

        g = ( g / ort(mp) ) / a(mp,mp-1)

        z(mp:igh,j) = z(mp:igh,j) + g * ort(mp:igh)

      end do

    end if

  end do

  return
end
function pythag ( a, b )
!
!*******************************************************************************
!
!! PYTHAG computes SQRT ( A**2 + B**2 ) carefully.
!
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A**2 + B**2 )
!
!    is reasonably accurate, but the formula can actually fail if
!    for example, A**2 is larger than the machine overflow.  The
!    formula can lose most of its accuracy if the sum of the squares
!    is very large or very small.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Modified:
!
!    02 March 2000
!
!  Parameters:
!
!    Input, real A, B, the two legs of a right triangle.
!
!    Output, real PYTHAG, the length of the hypotenuse.
!
  real a
  real b
  real p
  real pythag
  real r
  real s
  real t
  real u
!
  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

   10   continue

    t = 4.0 + r

    if ( t /= 4.0 ) then
      s = r / t
      u = 1.0 + 2.0 * s
      p = u * p
      r = ( s / u )**2 * r
      go to 10
    end if

  end if

  pythag = p

  return
end
subroutine qzhes ( nm, n, a, b, matz, z )
!
!*******************************************************************************
!
!! QZHES carries out transformations for a generalized eigenvalue problem.
!
!
!  Discussion:
!
!     this subroutine is the first step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real general matrices and
!     reduces one of them to upper Hessenberg form and the other
!     to upper triangular form using orthogonal transformations.
!     it is usually followed by QZIT,  qzval  and, possibly,  qzvec.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real general matrix.
!
!        b contains a real general matrix.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!     on output
!
!        a has been reduced to upper Hessenberg form.  the elements
!          below the first subdiagonal have been set to zero.
!
!        b has been reduced to upper triangular form.  the elements
!          below the main diagonal have been set to zero.
!
!        z contains the product of the right hand transformations if
!          matz has been set to .true.  otherwise, z is not referenced.
!
  integer n
  integer nm
!
  real a(nm,n)
  real b(nm,n)
  integer i
  integer j
  integer k
  integer l
  integer l1
  integer lb
  logical matz
  integer nk1
  integer nm1
  real r
  real rho
  real s
  real t
  real u1
  real u2
  real v1
  real v2
  real z(nm,n)
!
!  Set Z to the identity matrix.
!
  if ( matz ) then
    call rmat_ident ( nm, n, z )
  end if
!
!  Reduce b to upper triangular form
!
  if ( n <= 1 ) then
    return
  end if

  nm1 = n - 1

  do l = 1, n-1

    l1 = l + 1

    s = 0.0
    do i = l+1, n
      s = s + abs(b(i,l))
    end do

    if ( s /= 0.0 ) then

      s = s + abs(b(l,l))
      r = 0.0

      do i = l, n
        b(i,l) = b(i,l) / s
        r = r + b(i,l)**2
      end do

      r = sign(sqrt(r),b(l,l))
      b(l,l) = b(l,l) + r
      rho = r * b(l,l)

      do j = l+1, n

        t = 0.0
        do i = l, n
          t = t + b(i,l) * b(i,j)
        end do

        b(l:n,j) = b(l:n,j) - t * b(l:n,l) / rho

      end do

      do j = 1, n

        t = 0.0
        do i = l, n
          t = t + b(i,l) * a(i,j)
        end do

        a(l:n,j) = a(l:n,j) - t * b(l:n,l) / rho

      end do

      b(l,l) = -s * r
      b(l+1:n,l) = 0.0

    end if

  end do
!
!  Reduce A to upper Hessenberg form, while keeping b triangular
!
  if ( n == 2 ) then
    return
  end if

  do k = 1, n-2

     nk1 = nm1 - k

     do lb = 1, nk1

        l = n - lb
        l1 = l + 1
!
!  Zero a(l+1,k)
!
        s = abs(a(l,k)) + abs(a(l1,k))

        if ( s /= 0.0 ) then

        u1 = a(l,k) / s
        u2 = a(l1,k) / s
        r = sign ( sqrt ( u1*u1 + u2*u2 ), u1 )
        v1 =  -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1

        do j = k, n
           t = a(l,j) + u2 * a(l1,j)
           a(l,j) = a(l,j) + t * v1
           a(l1,j) = a(l1,j) + t * v2
        end do

        a(l1,k) = 0.0

        do j = l, n
           t = b(l,j) + u2 * b(l1,j)
           b(l,j) = b(l,j) + t * v1
           b(l1,j) = b(l1,j) + t * v2
        end do
!
!  Zero b(l+1,l)
!
        s = abs(b(l1,l1)) + abs(b(l1,l))

        if (s /= 0.0) then

        u1 = b(l1,l1) / s
        u2 = b(l1,l) / s
        r = sign(sqrt(u1*u1+u2*u2),u1)
        v1 =  -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1

        do i = 1, l1
           t = b(i,l1) + u2 * b(i,l)
           b(i,l1) = b(i,l1) + t * v1
           b(i,l) = b(i,l) + t * v2
        end do

        b(l1,l) = 0.0

        do i = 1, n
           t = a(i,l1) + u2 * a(i,l)
           a(i,l1) = a(i,l1) + t * v1
           a(i,l) = a(i,l) + t * v2
        end do

        if ( matz ) then

          do i = 1, n
            t = z(i,l1) + u2 * z(i,l)
            z(i,l1) = z(i,l1) + t * v1
            z(i,l) = z(i,l) + t * v2
          end do

        end if

        end if

      end if

    end do

  end do

  return
end
subroutine qzit ( nm, n, a, b, eps1, matz, z, ierr )
!
!*******************************************************************************
!
!! QZIT carries out iterations to solve a generalized eigenvalue problem.
!
!
!  Discussion:
!
!     this subroutine is the second step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
!     as modified in technical note nasa tn d-7305(1973) by ward.
!
!     this subroutine accepts a pair of real matrices, one of them
!     in upper Hessenberg form and the other in upper triangular form.
!     it reduces the Hessenberg matrix to quasi-triangular form using
!     orthogonal transformations while maintaining the triangular form
!     of the other matrix.  it is usually preceded by  qzhes  and
!     followed by  qzval  and, possibly,  qzvec.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper Hessenberg matrix.
!
!        b contains a real upper triangular matrix.
!
!        eps1 is a tolerance used to determine negligible elements.
!          eps1 = 0.0 (or negative) may be input, in which case an
!          element will be neglected only if it is less than roundoff
!          error times the norm of its matrix.  if the input eps1 is
!          positive, then an element will be considered negligible
!          if it is less than eps1 times the norm of its matrix.  a
!          positive value of eps1 may result in faster execution,
!          but less accurate results.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!        z contains, if matz has been set to .true., the
!          transformation matrix produced in the reduction
!          by  qzhes, if performed, or else the identity matrix.
!          if matz has been set to .false., z is not referenced.
!
!     on output
!
!        a has been reduced to quasi-triangular form.  the elements
!          below the first subdiagonal are still zero and no two
!          consecutive subdiagonal elements are nonzero.
!
!        b is still in upper triangular form, although its elements
!          have been altered.  the location b(n,1) is used to store
!          eps1 times the norm of b for later use by  qzval  and  qzvec.
!
!        z contains the product of the right hand transformations
!          (for both steps) if matz has been set to .true..
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
  integer n
  integer nm
!
  real a(nm,n)
  real a1
  real a11
  real a12
  real a2
  real a21
  real a22
  real a3
  real a33
  real a34
  real a43
  real a44
  real ani
  real anorm
  real b(nm,n)
  real b11
  real b12
  real b22
  real b33
  real b34
  real b44
  real bni
  real bnorm
  integer en
  integer enm2
  integer enorn
  real ep
  real eps1
  real epsa
  real epsb
  integer i
  integer ierr
  integer ish
  integer itn
  integer its
  integer j
  integer k
  integer k1
  integer k2
  integer km1
  integer l
  integer l1
  integer ld
  integer ll
  integer lm1
  integer lor1
  logical matz
  integer na
  logical notlas
  real r
  real s
  real sh
  real t
  real u1
  real u2
  real u3
  real v1
  real v2
  real v3
  real z(nm,n)
!
  ierr = 0
!
!  Compute epsa, epsb
!
  anorm = 0.0
  bnorm = 0.0

  do i = 1, n

     ani = 0.0
     if (i /= 1) ani = abs(a(i,i-1))
     bni = 0.0

     do j = i, n
        ani = ani + abs(a(i,j))
        bni = bni + abs(b(i,j))
     end do

     if (ani > anorm) anorm = ani
     if (bni > bnorm) bnorm = bni

  end do

  if (anorm == 0.0) anorm = 1.0
  if (bnorm == 0.0) bnorm = 1.0
  ep = eps1
  if (ep > 0.0) go to 50
!
!  Use roundoff level if eps1 is 0.
!
  ep = epsilon ( 1.0 )

   50 epsa = ep * anorm
  epsb = ep * bnorm
!
!  Reduce a to quasi-triangular form, while keeping b triangular.
!
  lor1 = 1
  enorn = n
  en = n
  itn = 30*n
!
!  Begin qz step
!
   60 if (en <= 2) go to 1001
  if (.not. matz) enorn = en
  its = 0
  na = en - 1
  enm2 = na - 1
   70 ish = 2
!
!  Check for convergence or reducibility.
!
  do ll = 1, en
     lm1 = en - ll
     l = lm1 + 1
     if (l == 1) go to 95
     if (abs(a(l,lm1)) <= epsa) go to 90
  end do

   90 continue

  a(l,lm1) = 0.0
  if (l < na) go to 95
!
!  1-by-1 or 2-by-2 block isolated
!
  en = lm1
  go to 60
!
!  Check for small top of b
!
   95 ld = l
  100 l1 = l + 1
  b11 = b(l,l)
  if (abs(b11) > epsb) go to 120
  b(l,l) = 0.0
  s = abs(a(l,l)) + abs(a(l1,l))
  u1 = a(l,l) / s
  u2 = a(l1,l) / s
  r = sign(sqrt(u1*u1+u2*u2),u1)
  v1 = -(u1 + r) / r
  v2 = -u2 / r
  u2 = v2 / v1

  do j = l, enorn
     t = a(l,j) + u2 * a(l1,j)
     a(l,j) = a(l,j) + t * v1
     a(l1,j) = a(l1,j) + t * v2
     t = b(l,j) + u2 * b(l1,j)
     b(l,j) = b(l,j) + t * v1
     b(l1,j) = b(l1,j) + t * v2
  end do

  if (l /= 1) a(l,lm1) = -a(l,lm1)
  lm1 = l
  l = l1
  go to 90
  120 a11 = a(l,l) / b11
  a21 = a(l1,l) / b11
  if (ish == 1) go to 140
!
!  Iteration strategy
!
  if (itn == 0) go to 1000
  if (its == 10) go to 155
!
!  Determine type of shift
!
  b22 = b(l1,l1)
  if (abs(b22) < epsb) b22 = epsb
  b33 = b(na,na)
  if (abs(b33) < epsb) b33 = epsb
  b44 = b(en,en)
  if (abs(b44) < epsb) b44 = epsb
  a33 = a(na,na) / b33
  a34 = a(na,en) / b44
  a43 = a(en,na) / b33
  a44 = a(en,en) / b44
  b34 = b(na,en) / b44
  t = 0.5 * (a43 * b34 - a33 - a44)
  r = t * t + a34 * a43 - a33 * a44
  if (r < 0.0) go to 150
!
!  Determine single shift zeroth column of a
!
  ish = 1
  r = sqrt(r)
  sh = -t + r
  s = -t - r
  if (abs(s-a44) < abs(sh-a44)) sh = s
!
!  Look for two consecutive small sub-diagonal elements of a.
! 
  do ll = ld, enm2
     l = enm2 + ld - ll
     if (l == ld) go to 140
     lm1 = l - 1
     l1 = l + 1
     t = a(l,l)
     if (abs(b(l,l)) > epsb) t = t - sh * b(l,l)
     if (abs(a(l,lm1)) <= abs(t/a(l1,l)) * epsa) go to 100
  end do

  140 continue

  a1 = a11 - sh
  a2 = a21
  if (l /= ld) a(l,lm1) = -a(l,lm1)
  go to 160
!
!  Determine double shift zeroth column of A.
!
  150 a12 = a(l,l1) / b22
  a22 = a(l1,l1) / b22
  b12 = b(l,l1) / b22
  a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) &
    / a21 + a12 - a11 * b12
  a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34
  a3 = a(l1+1,l1) / b22
  go to 160
!
!  Ad hoc shift
!
  155 a1 = 0.0
  a2 = 1.0
  a3 = 1.1605
  160 its = its + 1
  itn = itn - 1
  if (.not. matz) lor1 = ld
!
!  Main loop
!
  do k = l, na

     notlas = k /= na .and. ish == 2
     k1 = k + 1
     k2 = k + 2
     km1 = max(k-1,l)
     ll = min(en,k1+ish)
     if (notlas) go to 190
!
!  Zero a(k+1,k-1)
!
     if ( k /= l ) then
       a1 = a(k,km1)
       a2 = a(k1,km1)
     end if

     s = abs(a1) + abs(a2)

     if (s == 0.0) go to 70

     u1 = a1 / s
     u2 = a2 / s
     r = sign(sqrt(u1*u1+u2*u2),u1)
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     u2 = v2 / v1

     do j = km1, enorn
        t = a(k,j) + u2 * a(k1,j)
        a(k,j) = a(k,j) + t * v1
        a(k1,j) = a(k1,j) + t * v2
        t = b(k,j) + u2 * b(k1,j)
        b(k,j) = b(k,j) + t * v1
        b(k1,j) = b(k1,j) + t * v2
     end do

     if (k /= l) a(k1,km1) = 0.0
     go to 240
!
!  Zero a(k+1,k-1) and a(k+2,k-1)
!
  190    continue

     if ( k /= l ) then
       a1 = a(k,km1)
       a2 = a(k1,km1)
       a3 = a(k2,km1)
     end if

     s = abs(a1) + abs(a2) + abs(a3)

     if (s == 0.0) go to 260

     u1 = a1 / s
     u2 = a2 / s
     u3 = a3 / s
     r = sign(sqrt(u1*u1+u2*u2+u3*u3),u1)
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     v3 = -u3 / r
     u2 = v2 / v1
     u3 = v3 / v1

     do j = km1, enorn

        t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
        a(k,j) = a(k,j) + t * v1
        a(k1,j) = a(k1,j) + t * v2
        a(k2,j) = a(k2,j) + t * v3
        t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
        b(k,j) = b(k,j) + t * v1
        b(k1,j) = b(k1,j) + t * v2
        b(k2,j) = b(k2,j) + t * v3
     end do

     if ( k /= l ) then
       a(k1,km1) = 0.0
       a(k2,km1) = 0.0
     end if
!
!    zero b(k+2,k+1) and b(k+2,k)
!
     s = abs(b(k2,k2)) + abs(b(k2,k1)) + abs(b(k2,k))
     if (s == 0.0) go to 240
     u1 = b(k2,k2) / s
     u2 = b(k2,k1) / s
     u3 = b(k2,k) / s
     r = sign(sqrt(u1*u1+u2*u2+u3*u3),u1)
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     v3 = -u3 / r
     u2 = v2 / v1
     u3 = v3 / v1

     do i = lor1, ll
        t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
        a(i,k2) = a(i,k2) + t * v1
        a(i,k1) = a(i,k1) + t * v2
        a(i,k) = a(i,k) + t * v3
        t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
        b(i,k2) = b(i,k2) + t * v1
        b(i,k1) = b(i,k1) + t * v2
        b(i,k) = b(i,k) + t * v3
     end do

     b(k2,k) = 0.0
     b(k2,k1) = 0.0

     if ( matz ) then

       do i = 1, n
         t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
         z(i,k2) = z(i,k2) + t * v1
         z(i,k1) = z(i,k1) + t * v2
         z(i,k) = z(i,k) + t * v3
       end do

     end if
!
!    zero b(k+1,k)
!
240  continue

     s = abs(b(k1,k1)) + abs(b(k1,k))

     if ( s /= 0.0 ) then

     u1 = b(k1,k1) / s
     u2 = b(k1,k) / s
     r = sign(sqrt(u1*u1+u2*u2),u1)
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     u2 = v2 / v1

     do i = lor1, ll
        t = a(i,k1) + u2 * a(i,k)
        a(i,k1) = a(i,k1) + t * v1
        a(i,k) = a(i,k) + t * v2
        t = b(i,k1) + u2 * b(i,k)
        b(i,k1) = b(i,k1) + t * v1
        b(i,k) = b(i,k) + t * v2
     end do

     b(k1,k) = 0.0

     if ( matz ) then

      do i = 1, n
         t = z(i,k1) + u2 * z(i,k)
         z(i,k1) = z(i,k1) + t * v1
         z(i,k) = z(i,k) + t * v2
       end do

     end if

     end if

260  continue

  end do

  go to 70
!
!  Wet error: all eigenvalues have not converged after 30*n iterations
!
 1000 ierr = en
!
!  Save epsb for use by qzval and qzvec
!
 1001 if (n > 1) b(n,1) = epsb
  return
end
subroutine qzval ( nm, n, a, b, alfr, alfi, beta, matz, z )
!
!*******************************************************************************
!
!! QZVAL computes eigenvalues for a generalized eigenvalue problem.
!
!
!  Discussion:
!
!     this subroutine is the third step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real matrices, one of them
!     in quasi-triangular form and the other in upper triangular form.
!     it reduces the quasi-triangular matrix further, so that any
!     remaining 2-by-2 blocks correspond to pairs of complex
!     eigenvalues, and returns quantities whose ratios give the
!     generalized eigenvalues.  it is usually preceded by  qzhes
!     and  qzit  and may be followed by  qzvec.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper quasi-triangular matrix.
!
!        b contains a real upper triangular matrix.  in addition,
!          location b(n,1) contains the tolerance quantity (epsb)
!          computed and saved in  qzit.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!        z contains, if matz has been set to .true., the
!          transformation matrix produced in the reductions by qzhes
!          and qzit, if performed, or else the identity matrix.
!          if matz has been set to .false., z is not referenced.
!
!     on output
!
!        a has been reduced further to a quasi-triangular matrix
!          in which all nonzero subdiagonal elements correspond to
!          pairs of complex eigenvalues.
!
!        b is still in upper triangular form, although its elements
!          have been altered.  b(n,1) is unaltered.
!
!        alfr and alfi contain the real and imaginary parts of the
!          diagonal elements of the triangular matrix that would be
!          obtained if a were reduced completely to triangular form
!          by unitary transformations.  non-zero values of alfi occur
!          in pairs, the first member positive and the second negative.
!
!        beta contains the diagonal elements of the corresponding b,
!          normalized to be real and non-negative.  the generalized
!          eigenvalues are then the ratios ((alfr+i*alfi)/beta).
!
!        z contains the product of the right hand transformations
!          (for all three steps) if matz has been set to .true.
!
  integer n
  integer nm
!
  real a(nm,n)
  real a1
  real a11
  real a11i
  real a11r
  real a12
  real a12i
  real a12r
  real a1i
  real a2
  real a21
  real a22
  real a22i
  real a22r
  real a2i
  real an
  real alfi(n)
  real alfr(n)
  real b(nm,n)
  real b11
  real b12
  real b22
  real beta(n)
  real bn
  real c
  real cq
  real cz
  real d
  real di
  real dr
  real e
  real ei
  integer en
  real epsb
  integer i
  integer isw
  integer j
  logical matz
  integer na
  integer nn
  real r
  real s
  real sqi
  real sqr
  real ssi
  real ssr
  real szi
  real szr
  real t
  real ti
  real tr
  real u1
  real u2
  real v1
  real v2
  real z(nm,n)
!
  epsb = b(n,1)
  isw = 1
!
!  Find eigenvalues of quasi-triangular matrices.
!
  do nn = 1, n

     en = n + 1 - nn
     na = en - 1
     if (isw == 2) go to 505
     if (en == 1) go to 410
     if (a(en,na) /= 0.0) go to 420
!
!  1-by-1 block, one real root
!
410  continue

     alfr(en) = a(en,en)
     if (b(en,en) < 0.0) alfr(en) = -alfr(en)
     beta(en) = abs(b(en,en))
     alfi(en) = 0.0
     go to 510
!
!  2-by-2 block
!
420  continue

     if (abs(b(na,na)) <= epsb) go to 455
     if (abs(b(en,en)) > epsb) go to 430
     a1 = a(en,en)
     a2 = a(en,na)
     bn = 0.0
     go to 435

430  continue

     an = abs(a(na,na)) + abs(a(na,en)) + abs(a(en,na)) + abs(a(en,en))
     bn = abs(b(na,na)) + abs(b(na,en)) + abs(b(en,en))
     a11 = a(na,na) / an
     a12 = a(na,en) / an
     a21 = a(en,na) / an
     a22 = a(en,en) / an
     b11 = b(na,na) / bn
     b12 = b(na,en) / bn
     b22 = b(en,en) / bn
     e = a11 / b11
     ei = a22 / b22
     s = a21 / (b11 * b22)
     t = (a22 - e * b22) / b22
     if (abs(e) <= abs(ei)) go to 431
     e = ei
     t = (a11 - e * b11) / b11
431  continue

     c = 0.5 * (t - s * b12)
     d = c * c + s * (a12 - e * b12)
     if (d < 0.0) go to 480
!
!  Two real roots.
!  Zero both a(en,na) and b(en,na)
!
     e = e + (c + sign(sqrt(d),c))
     a11 = a11 - e * b11
     a12 = a12 - e * b12
     a22 = a22 - e * b22
     if (abs(a11) + abs(a12) < abs(a21) + abs(a22)) then
       go to 432
     end if
     a1 = a12
     a2 = a11
     go to 435
  432    a1 = a22
     a2 = a21
!
!  Choose and apply real z
!
  435    s = abs(a1) + abs(a2)
     u1 = a1 / s
     u2 = a2 / s
     r = sign(sqrt(u1*u1+u2*u2),u1)
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     u2 = v2 / v1

     do i = 1, en
        t = a(i,en) + u2 * a(i,na)
        a(i,en) = a(i,en) + t * v1
        a(i,na) = a(i,na) + t * v2
        t = b(i,en) + u2 * b(i,na)
        b(i,en) = b(i,en) + t * v1
        b(i,na) = b(i,na) + t * v2
     end do

     if ( matz ) then

       do i = 1, n
         t = z(i,en) + u2 * z(i,na)
         z(i,en) = z(i,en) + t * v1
         z(i,na) = z(i,na) + t * v2
       end do

     end if

  450    if (bn == 0.0) go to 475
     if (an < abs(e) * bn) go to 455
     a1 = b(na,na)
     a2 = b(en,na)
     go to 460
  455    a1 = a(na,na)
     a2 = a(en,na)
!
!  Choose and apply real q
!
  460    s = abs(a1) + abs(a2)
     if (s == 0.0) go to 475
     u1 = a1 / s
     u2 = a2 / s
     r = sign(sqrt(u1*u1+u2*u2),u1)
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     u2 = v2 / v1

     do j = na, n
        t = a(na,j) + u2 * a(en,j)
        a(na,j) = a(na,j) + t * v1
        a(en,j) = a(en,j) + t * v2
        t = b(na,j) + u2 * b(en,j)
        b(na,j) = b(na,j) + t * v1
        b(en,j) = b(en,j) + t * v2
     end do

  475    a(en,na) = 0.0
     b(en,na) = 0.0
     alfr(na) = a(na,na)
     alfr(en) = a(en,en)
     if (b(na,na) < 0.0) alfr(na) = -alfr(na)
     if (b(en,en) < 0.0) alfr(en) = -alfr(en)
     beta(na) = abs(b(na,na))
     beta(en) = abs(b(en,en))
     alfi(en) = 0.0
     alfi(na) = 0.0
     go to 505
!
!  Two complex roots
!
  480    e = e + c
     ei = sqrt(-d)
     a11r = a11 - e * b11
     a11i = ei * b11
     a12r = a12 - e * b12
     a12i = ei * b12
     a22r = a22 - e * b22
     a22i = ei * b22
     if (abs(a11r) + abs(a11i) + abs(a12r) + abs(a12i) < &
            abs(a21) + abs(a22r) + abs(a22i)) then
       go to 482
     end if

     a1 = a12r
     a1i = a12i
     a2 = -a11r
     a2i = -a11i
     go to 485
  482    a1 = a22r
     a1i = a22i
     a2 = -a21
     a2i = 0.0
!
!  Choose complex z
!
  485    cz = sqrt(a1*a1+a1i*a1i)
     if (cz == 0.0) go to 487
     szr = (a1 * a2 + a1i * a2i) / cz
     szi = (a1 * a2i - a1i * a2) / cz
     r = sqrt(cz*cz+szr*szr+szi*szi)
     cz = cz / r
     szr = szr / r
     szi = szi / r
     go to 490
  487    szr = 1.0
     szi = 0.0
  490    if (an < (abs(e) + ei) * bn) go to 492
     a1 = cz * b11 + szr * b12
     a1i = szi * b12
     a2 = szr * b22
     a2i = szi * b22
     go to 495
  492    a1 = cz * a11 + szr * a12
     a1i = szi * a12
     a2 = cz * a21 + szr * a22
     a2i = szi * a22
!
!  Choose complex q
!
  495    cq = sqrt(a1*a1+a1i*a1i)
     if (cq == 0.0) go to 497
     sqr = (a1 * a2 + a1i * a2i) / cq
     sqi = (a1 * a2i - a1i * a2) / cq
     r = sqrt(cq*cq+sqr*sqr+sqi*sqi)
     cq = cq / r
     sqr = sqr / r
     sqi = sqi / r
     go to 500
  497    sqr = 1.0
     sqi = 0.0
!
!  Compute diagonal elements that would result if transformations were applied
!
  500    ssr = sqr * szr + sqi * szi
     ssi = sqr * szi - sqi * szr
     i = 1
     tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22
     ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
     dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
     di = cq * szi * b12 + ssi * b22
     go to 503
  502    i = 2
     tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22
     ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
     dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
     di = -ssi * b11 - sqi * cz * b12
  503    t = ti * dr - tr * di
     j = na
     if (t < 0.0) j = en
     r = sqrt(dr*dr+di*di)
     beta(j) = bn * r
     alfr(j) = an * (tr * dr + ti * di) / r
     alfi(j) = an * t / r
     if (i == 1) go to 502
  505    isw = 3 - isw

  510 continue

  end do

  b(n,1) = epsb

  return
end
subroutine qzvec ( nm, n, a, b, alfr, alfi, beta, z )
!
!*******************************************************************************
!
!! QZVEC computes eigenvectors for a generalized eigenvalue problem.
!
!
!  Discussion:
!
!     this subroutine is the optional fourth step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real matrices, one of them in
!     quasi-triangular form (in which each 2-by-2 block corresponds to
!     a pair of complex eigenvalues) and the other in upper triangular
!     form.  it computes the eigenvectors of the triangular problem and
!     transforms the results back to the original coordinate system.
!     it is usually preceded by  qzhes,  qzit, and  qzval.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper quasi-triangular matrix.
!
!        b contains a real upper triangular matrix.  in addition,
!          location b(n,1) contains the tolerance quantity (epsb)
!          computed and saved in  qzit.
!
!        alfr, alfi, and beta  are vectors with components whose
!          ratios ((alfr+i*alfi)/beta) are the generalized
!          eigenvalues.  they are usually obtained from  qzval.
!
!        z contains the transformation matrix produced in the
!          reductions by  qzhes,  qzit, and  qzval, if performed.
!          if the eigenvectors of the triangular problem are
!          desired, z must contain the identity matrix.
!
!     on output
!
!        a is unaltered.  its subdiagonal elements provide information
!           about the storage of the complex eigenvectors.
!
!        b has been destroyed.
!
!        alfr, alfi, and beta are unaltered.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if alfi(i) == 0.0, the i-th eigenvalue is real and
!            the i-th column of z contains its eigenvector.
!          if alfi(i) /= 0.0, the i-th eigenvalue is complex.
!            if alfi(i) > 0.0, the eigenvalue is the first of
!              a complex pair and the i-th and (i+1)-th columns
!              of z contain its eigenvector.
!            if alfi(i) < 0.0, the eigenvalue is the second of
!              a complex pair and the (i-1)-th and i-th columns
!              of z contain the conjugate of its eigenvector.
!          each eigenvector is normalized so that the modulus
!          of its largest component is 1.0 .
!
  integer n
  integer nm
!
  real a(nm,n)
  real alfi(n)
  real alfm
  real alfr(n)
  real almi
  real almr
  real b(nm,n)
  real beta(n)
  real betm
  real d
  real di
  real dr
  integer en
  integer enm2
  real epsb
  integer i
  integer ii
  integer isw
  integer j
  integer jj
  integer k
  integer m
  integer na
  integer nn
  real q
  real r
  real ra
  real rr
  real s
  real sa
  real t
  real t1
  real t2
  real ti
  real tr
  real w
  real w1
  real x
  real x1
  real y
  real z(nm,n)
  real z1
  real zz
!
  epsb = b(n,1)
  isw = 1

  do nn = 1, n

     en = n + 1 - nn
     na = en - 1
     if (isw == 2) go to 795
     if (alfi(en) /= 0.0) go to 710
!
!  Real vector
!
     m = en
     b(en,en) = 1.0
     if (na == 0) go to 800
     alfm = alfr(m)
     betm = beta(m)

     do ii = 1, na

        i = en - ii
        w = betm * a(i,i) - alfm * b(i,i)
        r = 0.0

        do j = m, en
          r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)
        end do

        if (i == 1 .or. isw == 2) go to 630
        if (betm * a(i,i-1) == 0.0) go to 630
        zz = w
        s = r
        go to 690
  630       m = i
        if (isw == 2) go to 640
!
!  Real 1-by-1 block
!
        t = w
        if (w == 0.0) t = epsb
        b(i,en) = -r / t
        go to 700
!
!  Real 2-by-2 block
!
  640       x = betm * a(i,i+1) - alfm * b(i,i+1)
        y = betm * a(i+1,i)
        q = w * zz - x * y
        t = (x * s - zz * r) / q
        b(i,en) = t
        if (abs(x) <= abs(zz)) go to 650
        b(i+1,en) = (-r - w * t) / x
        go to 690
  650       b(i+1,en) = (-s - y * t) / zz
  690       isw = 3 - isw
  700    continue

     end do
!
!  end real vector
!
     go to 800
!
!  Complex vector
!
  710    m = na
     almr = alfr(m)
     almi = alfi(m)
     betm = beta(m)
!
!  last vector component chosen imaginary so eigenvector matrix is triangular
!
     y = betm * a(en,na)
     b(na,na) = -almi * b(en,en) / y
     b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y
     b(en,na) = 0.0
     b(en,en) = 1.0
     enm2 = na - 1

     do ii = 1, enm2

        i = na - ii
        w = betm * a(i,i) - almr * b(i,i)
        w1 = -almi * b(i,i)
        ra = 0.0
        sa = 0.0

        do j = m, en
           x = betm * a(i,j) - almr * b(i,j)
           x1 = -almi * b(i,j)
           ra = ra + x * b(j,na) - x1 * b(j,en)
           sa = sa + x * b(j,en) + x1 * b(j,na)
        end do

        if (i == 1 .or. isw == 2) go to 770
        if (betm * a(i,i-1) == 0.0) go to 770
        zz = w
        z1 = w1
        r = ra
        s = sa
        isw = 2
        go to 790
  770       m = i
        if (isw == 2) go to 780
!
!  Complex 1-by-1 block
!
        tr = -ra
        ti = -sa
  773       dr = w
        di = w1
!
!  Complex divide (t1,t2) = (tr,ti) / (dr,di)
!
  775       if (abs(di) > abs(dr)) go to 777
        rr = di / dr
        d = dr + di * rr
        t1 = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        go to (787,782), isw
  777       rr = dr / di
        d = dr * rr + di
        t1 = (tr * rr + ti) / d
        t2 = (ti * rr - tr) / d
        go to (787,782), isw
!
!  Complex 2-by-2 block
!
  780       x = betm * a(i,i+1) - almr * b(i,i+1)
        x1 = -almi * b(i,i+1)
        y = betm * a(i+1,i)
        tr = y * ra - w * r + w1 * s
        ti = y * sa - w * s - w1 * r
        dr = w * zz - w1 * z1 - x * y
        di = w * z1 + w1 * zz - x1 * y
        if (dr == 0.0 .and. di == 0.0) dr = epsb
        go to 775
  782       b(i+1,na) = t1
        b(i+1,en) = t2
        isw = 1
        if (abs(y) > abs(w) + abs(w1)) go to 785
        tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
        ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
        go to 773
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y
        t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y
  787       b(i,na) = t1
        b(i,en) = t2
  790    continue

     end do
!
!  End complex vector
!
  795    isw = 3 - isw
  800 continue

  end do
!
!  End back substitution.
!  Transform to original coordinate system.
!
  do jj = 1, n

     j = n + 1 - jj

     do i = 1, n

        zz = 0.0

        do k = 1, j
          zz = zz + z(i,k) * b(k,j)
        end do

        z(i,j) = zz

      end do

  end do
!
!  Normalize so that modulus of largest component of each vector is 1.
!  (isw is 1 initially from before)
!
  do j = 1, n

     d = 0.0
     if (isw == 2) go to 920
     if (alfi(j) /= 0.0) go to 945

     do i = 1, n
       d = max ( d, abs(z(i,j)) )
     end do

     z(1:n,j) = z(1:n,j) / d

     go to 950

  920    continue

     do i = 1, n
       r = abs(z(i,j-1)) + abs(z(i,j))
       if (r /= 0.0) then
         r = r * sqrt((z(i,j-1)/r)**2 + (z(i,j)/r)**2)
       end if
       if (r > d) d = r
     end do

     z(1:n,j-1) = z(1:n,j-1) / d
     z(1:n,j) = z(1:n,j) / d

  945    isw = 3 - isw
  950 continue

  end do

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
subroutine ratqr ( n, eps1, d, e, e2, m, w, ind, bd, type, idef, ierr )
!
!*******************************************************************************
!
!! RATQR computes selected eigenvalues of a real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the algebraically smallest or largest
!    eigenvalues of a symmetric tridiagonal matrix by the
!    rational QR method with Newton corrections.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        n is the order of the matrix.
!
!        eps1 is a theoretical absolute error tolerance for the
!          computed eigenvalues.  if the input eps1 is non-positive,
!          or indeed smaller than its default value, it is reset
!          at each iteration to the respective default value,
!          namely, the product of the relative machine precision
!          and the magnitude of the current eigenvalue iterate.
!          the theoretical absolute error in the k-th eigenvalue
!          is usually not greater than k times eps1.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        m is the number of eigenvalues to be found.
!
!        idef should be set to 1 if the input matrix is known to be
!          positive definite, to -1 if the input matrix is known to
!          be negative definite, and to 0 otherwise.
!
!        type should be set to .true. if the smallest eigenvalues
!          are to be found, and to .false. if the largest eigenvalues
!          are to be found.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered (unless w overwrites d).
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is set to 0.0 if the smallest eigenvalues have been
!          found, and to 2.0 if the largest eigenvalues have been
!          found.  e2 is otherwise unaltered (unless overwritten by bd).
!
!        w contains the m algebraically smallest eigenvalues in
!          ascending order, or the m largest eigenvalues in
!          descending order.  if an error exit is made because of
!          an incorrect specification of idef, no eigenvalues
!          are found.  if the newton iterates for a particular
!          eigenvalue are not monotone, the best estimate obtained
!          is returned and ierr is set.  w may coincide with d.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w:
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc..
!
!        bd contains refined bounds for the theoretical errors of the
!          corresponding eigenvalues in w.  these bounds are usually
!          within the tolerance specified by eps1.  bd may coincide
!          with e2.
!
!        ierr is set to
!          zero       for normal return,
!          6*n+1      if  idef  is set to 1 and  type  to .true.
!                     when the matrix is not positive definite, or
!                     if  idef  is set to -1 and  type  to .false.
!                     when the matrix is not negative definite,
!          5*n+k      if successive iterates to the k-th eigenvalue
!                     are not monotone increasing, where k refers
!                     to the last such occurrence.
!
  integer n
!
  real bd(n)
  real d(n)
  real delta
  real e(n)
  real e2(n)
  real ep
  real eps1
  real err
  real f
  integer i
  integer idef
  integer ierr
  integer ii
  integer ind(n)
  integer j
  integer jdef
  integer jj
  integer k
  integer k1
  integer m
  real p
  real q
  real qp
  real r
  real s
  real tot
  logical type
  real w(n)
!
  ierr = 0
  jdef = idef
  w(1:n) = d(1:n)

  if ( .not. type ) then
    j = 1
    go to 400
  end if

40 continue

  err = 0.0
  s = 0.0
!
!  Look for small sub-diagonal entries and define
!  initial shift from lower gerschgorin bound.
!                copy e2 array into bd
!
  tot = w(1)
  q = 0.0
  j = 0

  do i = 1, n

     p = q
     if (i == 1) go to 60
     if ( p > ( abs ( d(i) ) + abs( d(i-1) ) ) * epsilon ( 1.0 ) ) then
       go to 80
     end if

60   continue

     e2(i) = 0.0

80   continue

     bd(i) = e2(i)
!
!    count also if element of e2 has underflowed
!
     if (e2(i) == 0.0) j = j + 1
     ind(i) = j
     q = 0.0
     if (i /= n) q = abs(e(i+1))
     tot = min(w(i)-p-q,tot)
  end do

  if (jdef == 1 .and. tot < 0.0) then
    go to 140
  end if

  w(1:n) = w(1:n) - tot

  go to 160
  140 tot = 0.0

  160 continue

  do k = 1, m
!
!    next qr transformation
!
  180    tot = tot + s
     delta = w(n) - s
     i = n
     f = abs ( tot ) * epsilon ( 1.0 )
     if (eps1 < f) eps1 = f
     if (delta > eps1) go to 190
     if (delta < (-eps1)) go to 1000
     go to 300
!
!  Replace small sub-diagonal squares by zero to reduce the incidence of 
!  underflows
!
  190    if (k == n) go to 210
     k1 = k + 1

     do j = k+1, n
       if ( bd(j) <= ( abs( w(j) + w(j-1) ) * epsilon ( 1.0 ) ) ** 2 ) then
         bd(j) = 0.0
       end if
     end do

  210    f = bd(n) / delta
     qp = delta + f
     p = 1.0
     if (k == n) go to 260
     k1 = n - k

     do ii = 1, k1

        i = n - ii
        q = w(i) - s - f
        r = q / qp
        p = p * r + 1.0
        ep = f * r
        w(i+1) = qp + ep
        delta = q - ep
        if (delta > eps1) go to 220
        if (delta < (-eps1)) go to 1000
        go to 300
  220   f = bd(i) / q
        qp = delta + f
        bd(i+1) = qp * ep
    end do

  260    w(k) = qp
     s = qp / p
     if (tot + s > tot) go to 180
!
!  Set error: irregular end of iteration.
!  deflate minimum diagonal element
!
     ierr = 5 * n + k
     s = 0.0
     delta = qp

     do j = k, n
       if (w(j) <= delta) then
         i = j
         delta = w(j)
       end if
     end do
!
!  Convergence
!
  300    continue

     if (i < n) bd(i+1) = bd(i) * f / qp
     ii = ind(i)

     k1 = i - k

     do jj = 1, k1
        j = i - jj
        w(j+1) = w(j) - s
        bd(j+1) = bd(j)
        ind(j+1) = ind(j)
     end do

     w(k) = tot
     err = err + abs(delta)
     bd(k) = err
     ind(k) = ii

  end do

  if ( type ) then
    return
  end if

  f = bd(1)
  e2(1) = 2.0
  bd(1) = f
  j = 2
!
!  Negate elements of w for largest values
!
  400 continue

  w(1:n) = - w(1:n)
  jdef = -jdef

  if ( j == 1 ) then
    go to 40
  end if

  return
!
!  Set error: idef specified incorrectly
!
 1000 ierr = 6 * n + 1
  return
end
subroutine rebak ( nm, n, b, dl, m, z )
!
!*******************************************************************************
!
!! REBAK determines eigenvectors by undoing the REDUC transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure rebaka,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine forms the eigenvectors of a generalized
!     symmetric eigensystem by back transforming those of the
!     derived symmetric matrix determined by  reduc.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix system.
!
!        b contains information about the similarity transformation
!          (cholesky decomposition) used in the reduction by  reduc
!          in its strict lower triangle.
!
!        dl contains further information about the transformation.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
  integer m
  integer nm
  integer n
!
  real b(nm,n)
  real dl(n)
  integer i
  integer j
  integer k
  real x
  real z(nm,m)
!
  do j = 1, m

     do i = n, 1, -1

        x = z(i,j)

        do k = i+1, n
          x = x - b(k,i) * z(k,j)
        end do

        z(i,j) = x / dl(i)

    end do
  end do

  return
end
subroutine rebakb ( nm, n, b, dl, m, z )
!
!*******************************************************************************
!
!! REBAKB determines eigenvectors by undoing the REDUC2 transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure rebakb,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine forms the eigenvectors of a generalized
!     symmetric eigensystem by back transforming those of the
!     derived symmetric matrix determined by  reduc2.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix system.
!
!        b contains information about the similarity transformation
!          (cholesky decomposition) used in the reduction by  reduc2
!          in its strict lower triangle.
!
!        dl contains further information about the transformation.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
  integer m
  integer nm
  integer n
!
  real b(nm,n)
  real dl(n)
  integer i
  integer j
  integer k
  real x
  real z(nm,m)
!
  do j = 1, m

    do i = n, 1, -1

      x = dl(i) * z(i,j)

      do k = 1, i-1
        x = x + b(i,k) * z(k,j)
      end do

      z(i,j) = x

    end do

  end do

  return
end
subroutine reduc ( nm, n, a, b, dl, ierr )
!
!*******************************************************************************
!
!! REDUC reduces the eigenvalue problem A*x=lambda*B*x to A*x=lambda*x.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure reduc1,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine reduces the generalized symmetric eigenproblem
!     ax=(lambda)bx, where b is positive definite, to the standard
!     symmetric eigenproblem using the cholesky factorization of b.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices a and b.  if the cholesky
!          factor l of b is already available, n should be prefixed
!          with a minus sign.
!
!        a and b contain the real symmetric input matrices.  only the
!          full upper triangles of the matrices need be supplied.  if
!          n is negative, the strict lower triangle of b contains,
!          instead, the strict lower triangle of its cholesky factor l.
!
!        dl contains, if n is negative, the diagonal elements of l.
!
!     on output
!
!        a contains in its full lower triangle the full lower triangle
!          of the symmetric matrix derived from the reduction to the
!          standard form.  the strict upper triangle of a is unaltered.
!
!        b contains in its strict lower triangle the strict lower
!          triangle of its cholesky factor l.  the full upper
!          triangle of b is unaltered.
!
!        dl contains the diagonal elements of l.
!
!        ierr is set to
!          zero       for normal return,
!          7*n+1      if b is not positive definite.
!
  integer n
  integer nm
!
  real a(nm,n)
  real b(nm,n)
  real dl(n)
  integer i
  integer ierr
  integer j
  integer k
  integer nn
  real x
  real y
!
  ierr = 0
  nn = abs ( n )
!
!  form l in the arrays b and dl
!
  do i = 1, n

     do j = i, n

        x = b(i,j)

        do k = 1, i - 1
          x = x - b(i,k) * b(j,k)
        end do

        if ( j == i ) then

          if ( x <= 0.0 ) then
            write ( *, * ) ' '
            write ( *, * ) 'REDUC - Fatal error!'
            write ( *, * ) '  The matrix is not positive definite.'
            ierr = 7 * n + 1
            return
          end if

          y = sqrt(x)
          dl(i) = y
        else
          b(j,i) = x / y
        end if

    end do

  end do
!
!  form the transpose of the upper triangle of inv(l)*a
!  in the lower triangle of the array a
!
  do i = 1, nn

     y = dl(i)

     do j = i, nn

        x = a(i,j)

        do k = 1, i - 1
          x = x - b(i,k) * a(j,k)
        end do

        a(j,i) = x / y

      end do

  end do
!
!  pre-multiply by inv(l) and overwrite
!
  do j = 1, nn

     do i = j, nn

        x = a(i,j)

        do k = j, i-1
          x = x - a(k,j) * b(i,k)
        end do

        do k = 1, j-1
          x = x - a(j,k) * b(i,k)
        end do

        a(i,j) = x / dl(i)

    end do

  end do

  return
end
subroutine reduc2 ( nm, n, a, b, dl, ierr )
!
!*******************************************************************************
!
!! REDUC2 reduces the eigenvalue problem A*B*x=lamdba*x to A*x=lambda*x.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure reduc2,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine reduces the generalized symmetric eigenproblems
!     abx=(lambda)x or bay=(lambda)y, where b is positive definite,
!     to the standard symmetric eigenproblem using the cholesky
!     factorization of b.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices a and b.  if the cholesky
!          factor l of b is already available, n should be prefixed
!          with a minus sign.
!
!        a and b contain the real symmetric input matrices.  only the
!          full upper triangles of the matrices need be supplied.  if
!          n is negative, the strict lower triangle of b contains,
!          instead, the strict lower triangle of its cholesky factor l.
!
!        dl contains, if n is negative, the diagonal elements of l.
!
!     on output
!
!        a contains in its full lower triangle the full lower triangle
!          of the symmetric matrix derived from the reduction to the
!          standard form.  the strict upper triangle of a is unaltered.
!
!        b contains in its strict lower triangle the strict lower
!          triangle of its cholesky factor l.  the full upper
!          triangle of b is unaltered.
!
!        dl contains the diagonal elements of l.
!
!        ierr is set to
!          zero       for normal return,
!          7*n+1      if b is not positive definite.
!
  integer n
  integer nm
!
  real a(nm,n)
  real b(nm,n)
  real dl(n)
  integer i
  integer ierr
  integer j
  integer k
  integer nn
  real x
  real y
!
  ierr = 0
  nn = abs ( n )
!
!  form l in the arrays b and dl
!
  do i = 1, n

     do j = i, n

        x = b(i,j)

        do k = 1, i - 1
          x = x - b(i,k) * b(j,k)
        end do

        if ( j == i ) then

          if ( x <= 0.0 ) then
            write ( *, * ) ' '
            write ( *, * ) 'REDUC2 - Fatal error!'
            write ( *, * ) '  The matrix is not positive definite.'
            ierr = 7 * n + 1
            return
          end if

          y = sqrt(x)
          dl(i) = y

        else

          b(j,i) = x / y

        end if

    end do

  end do
!
!  form the lower triangle of a*l
!  in the lower triangle of the array a
!
  do i = 1, nn

     do j = 1, i

        x = a(j,i) * dl(j)

        do k = j+1, i
          x = x + a(k,i) * b(k,j)
        end do

        do k = i+1, nn
          x = x + a(i,k) * b(k,j)
        end do

        a(i,j) = x

     end do

  end do
!
!  Pre-multiply by transpose(l) and overwrite
!
  do i = 1, nn

    y = dl(i)

    do j = 1, i

      x = y * a(i,j)

      do k = i+1, nn
        x = x + a(k,j) * b(k,i)
      end do

      a(i,j) = x

    end do

  end do

  return
end
subroutine rg ( nm, n, a, wr, wi, matz, z, ierr )
!
!*******************************************************************************
!
!! RG computes eigenvalues and eigenvectors of a real general matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real general matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.  complex conjugate
!        pairs of eigenvalues appear consecutively with the
!        eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for hqr
!           and hqr2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real a(nm,n)
  real fv1(n)
  integer ierr
  integer is1
  integer is2
  integer iv1(n)
  integer matz
  real wi(n)
  real wr(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call balanc ( nm, n, a, is1, is2, fv1 )

  call elmhes ( nm, n, is1, is2, a, iv1 )

  if ( matz == 0 ) then

    call hqr ( nm, n, is1, is2, a, wr, wi, ierr )

    if ( ierr /= 0 ) then
      return
    end if

  else

    call eltran ( nm, n, is1, is2, a, iv1, z )

    call hqr2 ( nm, n, is1, is2, a, wr, wi, z, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call balbak ( nm, n, is1, is2, fv1, n, z )

  end if

  return
end
subroutine rgg ( nm, n, a, b, alfr, alfi, beta, matz, z, ierr )
!
!*******************************************************************************
!
!! RGG computes eigenvalues/vectors for the generalized problem A*x = lambda*B*x.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real general generalized eigenproblem  ax = (lambda)bx.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real general matrix.
!
!        b  contains a real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        alfr  and  alfi  contain the real and imaginary parts,
!        respectively, of the numerators of the eigenvalues.
!
!        beta  contains the denominators of the eigenvalues,
!        which are thus given by the ratios  (alfr+i*alfi)/beta.
!        complex conjugate pairs of eigenvalues appear consecutively
!        with the eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for QZIT.
!           the normal completion code is zero.
!
  integer n
  integer nm
!
  real a(nm,n)
  real alfi(n)
  real alfr(n)
  real b(nm,n)
  real beta(n)
  real eps1
  integer ierr
  integer matz
  logical tf
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  eps1 = 0.0

  if ( matz == 0 ) then
    tf = .false.
  else
    tf = .true.
  end if

  call qzhes ( nm, n, a, b, tf, z )

  call qzit ( nm, n, a, b, eps1, tf, z, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  call qzval ( nm, n, a, b, alfr, alfi, beta, tf, z )

  if ( matz /= 0 ) then
    call qzvec ( nm, n, a, b, alfr, alfi, beta, z )
  end if

  return
end
subroutine rmat_ident ( lda, n, a )
!
!*******************************************************************************
!
!! RMAT_IDENT sets the square matrix A to the identity.
!
!
!  Modified:
!
!    24 March 2000
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
!    Output, real A(LDA,N), the matrix which has been
!    set to the identity.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 1.0
      else
        a(i,j) = 0.0
      end if
    end do
  end do

  return
end
subroutine rmat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! RMAT_PRINT prints a real matrix, with an optional title.
!
!
!  Modified:
!
!    24 March 2000
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
!    Input, real A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer jhi
  integer jlo
  integer m
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 5
    jhi = min ( jlo + 4, n )
    write ( *, * ) ' '
    write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )
    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine rs ( nm, n, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RS computes eigenvalues and eigenvectors of real symmetric matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and TQL2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real a(nm,n)
  real fv1(n)
  real fv2(n)
  integer ierr
  integer matz
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

  end if

  return
end
subroutine rsb ( nm, n, mb, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSB computes eigenvalues and eigenvectors of a real symmetric band matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric band matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        mb  is the half band width of the matrix, defined as the
!        number of adjacent diagonals, including the principal
!        diagonal, required to specify the non-zero portion of the
!        lower triangle of the matrix.
!
!        a  contains the lower triangle of the real symmetric
!        band matrix.  its lowest subdiagonal is stored in the
!        last  n+1-mb  positions of the first column, its next
!        subdiagonal in the last  n+2-mb  positions of the
!        second column, further subdiagonals similarly, and
!        finally its principal diagonal in the  n  positions
!        of the last column.  contents of storages not part
!        of the matrix are arbitrary.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and TQL2.  the normal completion code is zero.
!
  integer mb
  integer n
  integer nm
!
  real a(nm,mb)
  real fv1(n)
  real fv2(n)
  integer ierr
  integer matz
  logical tf
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( mb <= 0 ) then
    ierr = 12 * n
    return
  end if

  if ( mb > n ) then
    ierr = 12 * n
    return
  end if

  if ( matz == 0 ) then
 
    tf = .false.

    call bandr ( nm, n, mb, a, w, fv1, fv2, tf, z )

    call tqlrat ( n, w, fv2, ierr )

  else

    tf = .true.

    call bandr ( nm, n, mb, a, w, fv1, fv1, tf, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

  end if

  return
end
subroutine rsg ( nm, n, a, b, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSG computes eigenvalues/vectors, A*x=lambda*B*x, A symmetric, B pos-def.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Modified:
!
!    12 June 2000
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and TQL2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real a(nm,n)
  real b(nm,n)
  real fv1(n)
  real fv2(n)
  integer ierr
  integer matz
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    write ( *, * ) ' '
    write ( *, * ) 'RSG - Fatal error!'
    write ( *, * ) '  N > NM.'
    return
  end if

  call reduc ( nm, n, a, b, fv2, ierr )

  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RSG - Fatal error!'
    write ( *, * ) '  Error return from REDUC.'
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RSG - Warning!'
      write ( *, * ) '  Error return from TQLRAT!'
      return
    end if

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RSG - Fatal error!'
      write ( *, * ) '  Error return from TQL2!'
      return
    end if

    call rebak ( nm, n, b, fv2, n, z )

  end if

  return
end
subroutine rsgab ( nm, n, a, b, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSGAB computes eigenvalues/vectors, A*B*x=lambda*x, A symmetric, B pos-def.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  abx = (lambda)x.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and TQL2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real a(nm,n)
  real b(nm,n)
  real fv1(n)
  real fv2(n)
  integer ierr
  integer matz
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call reduc2 ( nm, n, a, b, fv2, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call rebak ( nm, n, b, fv2, n, z )

  end if

  return
end
subroutine rsgba ( nm, n, a, b, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSGBA computes eigenvalues/vectors, B*A*x=lambda*x, A symmetric, B pos-def.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  bax = (lambda)x.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and TQL2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real a(nm,n)
  real b(nm,n)
  real fv1(n)
  real fv2(n)
  integer ierr
  integer matz
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  call reduc2 ( nm, n, a, b, fv2, ierr )

  if (ierr /= 0) then
    return
  end if

  if ( matz == 0 ) then

    call tred1 ( nm, n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( nm, n, a, w, fv1, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      return
    end if

    call rebakb ( nm, n, b, fv2, n, z )

  end if

  return
end
subroutine rsm ( nm, n, a, w, m, z, ierr )
!
!*******************************************************************************
!
!! RSM computes eigenvalues, some eigenvectors, real symmetric matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find all of the eigenvalues and some of the eigenvectors
!     of a real symmetric matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        m  the eigenvectors corresponding to the first m eigenvalues
!           are to be computed.
!           if m = 0 then no eigenvectors are computed.
!           if m = n then all of the eigenvectors are computed.
!
!     on output
!
!        w  contains all n eigenvalues in ascending order.
!
!        z  contains the orthonormal eigenvectors associated with
!           the first m eigenvalues.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat,
!           imtqlv and TINVIT.  the normal completion code is zero.
!
  integer n
  integer nm
  integer m
!
  real a(nm,n)
  real fwork1(n)
  real fwork2(n)
  real fwork3(n)
  real fwork4(n)
  integer ierr
  integer iwork(n)
  integer k1
  integer k2
  integer k3
  integer k4
  real w(n)
  real z(nm,m)
!
  if ( n > nm .or. m > nm ) then
    ierr = 10 * n
    return
  end if

  k1 = 1
  k2 = k1 + n
  k3 = k2 + n
  k4 = k3 + n

  if ( m <= 0 ) then

    call tred1 ( nm, n, a, w, fwork1, fwork2 )

    call tqlrat ( n, w, fwork2, ierr )

  else

    call tred1 ( nm, n, a, fwork1, fwork2, fwork3 )

    call imtqlv ( n, fwork1, fwork2, fwork3, w, iwork, ierr, fwork4 )

    call tinvit ( nm, n, fwork1, fwork2, fwork3, m, w, iwork, z, ierr )

    call trbak1 ( nm, n, a, fwork2, m, z )

  end if

  return
end
subroutine rsp ( nm, n, nv, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RSP computes eigenvalues and eigenvectors of real symmetric packed matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric packed matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        nv  is an integer variable set equal to the
!        dimension of the array  a  as specified for
!        a  in the calling program.  nv  must not be
!        less than  n*(n+1)/2.
!
!        a  contains the lower triangle of the real symmetric
!        packed matrix stored row-wise.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and TQL2.  the normal completion code is zero.
!
  integer n
  integer nm
  integer nv
!
  real a(nv)
  real fv1(n)
  real fv2(n)
  integer ierr
  integer matz
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( ( n * ( n + 1 ) ) / 2 > nv ) then
    ierr = 20 * n
    return
  end if

  call tred3 ( n, nv, a, w, fv1, fv2 )

  if ( matz == 0 ) then

    call tqlrat ( n, w, fv2, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RSP - Fatal error!'
      write ( *, * ) '  Error return from TQLRAT.'
      return
    end if

  else

    call rmat_ident ( nm, n, z )

    call tql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RSP - Fatal error!'
      write ( *, * ) '  Error return from TQL2.'
      return
    end if

    call  trbak3 ( nm, n, nv, a, n, z )

  end if

  return
end
subroutine rspp ( nm, n, nv, a, w, matz, z, ierr, m, type )
!
!*******************************************************************************
!
!! RSPP computes some eigenvalues/vectors, real symmetric packed matrix.
!
!
!  Discussion:
!
!    This routine calls the appropriate routines for the following problem:
!
!    given a symmetric matrix a, which is stored in a packed mode, find
!    the m smallest or largest eigenvalues, and corresponding eigenvectors.
!    the standard eispack routine rsp returns all eigenvalues and eigenvectors.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the eigenvector array z.
!
!    Input, integer N, the order of a, the number of rows and columns in the
!    original matrix.
!
!        nv  is an integer variable set equal to the
!        dimension of the array  a  as specified for
!        a  in the calling program.  nv  must not be
!        less than  n*(n+1)/2.
!
!    Input, real a(*), on input the lower triangle of the
!    real symmetric matrix, stored row-wise in the vector,
!    in the order a(1,1), / a(2,1), a(2,2), / a(3,1), a(3,2), a(3,3)/
!    and so on.  a will require (n*(n+1))/2 entries of storage.
!
!    Output, real W(M), the eigenvalues requested.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!    Output, real Z(NM,M), the eigenvectors.
!
!    Output, integer IERR, error flag from ratqr.  ierr=0 on normal return.
!    ierr nonzero, in this case, means that the algorithm broke
!    down while computing an eigenvalue.
!
!    Input, integer M, the number of eigenvalues to be found.
!
!    Input, logical TYPE, set to .true. if the smallest eigenvalues
!    are to be found, or .false. if the largest ones are sought.
!
  integer m
  integer n
  integer nm
  integer nv
!
  real a(nv)
  real bd(n)
  real eps1
  integer idef
  integer ierr
  integer iwork(n)
  integer matz
  integer nv
  logical type
  real w(m)
  real work1(n)
  real work2(n)
  real work3(n)
  real z(nm,m)
!
!  IDEF =
!    -1 if the matrix is known to be negative definite, 
!    +1 if the matrix is known to be positive definite, or
!    0 otherwise.
!
  idef = 0
!
!  Reduce to symmetric tridiagonal form.
!
  call tred3 ( n, nv, a, work1, work2, work3 )
!
!  Find the eigenvalues.
!
  eps1 = 0.0

  call ratqr ( n, eps1, work1, work2, work3, m, w, iwork, &
    bd, type, idef, ierr )

  if ( ierr /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RSPP - Fatal error!'
    write ( *, * ) '  Error return from RATQR.'
    return
  end if
!
!  Find eigenvectors for the first M eigenvalues.
!
  if ( matz /= 0 ) then

    call tinvit ( nm, n, work1, work2, work3, m, w, iwork, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RSPP - Fatal error!'
      write ( *, * ) '  Error return from TINVIT.'
      return
    end if
!
!  Reverse the transformation.
!
    call trbak3 ( nm, n, nv, a, m, z )

  end if

  return
end
subroutine rst ( nm, n, w, e, matz, z, ierr )
!
!*******************************************************************************
!
!! RST computes eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric tridiagonal matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix.
!
!        w  contains the diagonal elements of the real
!        symmetric tridiagonal matrix.
!
!        e  contains the subdiagonal elements of the matrix in
!        its last n-1 positions.  e(1) is arbitrary.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for imtql1
!           and IMTQL2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real e(n)
  integer ierr
  integer matz
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    ierr = 10 * n
    return
  end if

  if ( matz == 0 ) then

    call imtql1 ( n, w, e, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RST - Fatal error!'
      write ( *, * ) '  Error return from IMTQL1.'
      return
    end if

  else

    call rmat_ident ( nm, n, z )

    call imtql2 ( nm, n, w, e, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RST - Fatal error!'
      write ( *, * ) '  Error return from IMTQL2.'
      return
    end if

  end if

  return
end
subroutine rt ( nm, n, a, w, matz, z, ierr )
!
!*******************************************************************************
!
!! RT computes eigenvalues/vectors, real sign-symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a special real tridiagonal matrix.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the special real tridiagonal matrix in its
!        first three columns.  the subdiagonal elements are stored
!        in the last  n-1  positions of the first column, the
!        diagonal elements in the second column, and the superdiagonal
!        elements in the first  n-1  positions of the third column.
!        elements  a(1,1)  and  a(n,3)  are arbitrary.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for imtql1
!           and IMTQL2.  the normal completion code is zero.
!
  integer n
  integer nm
!
  real fv1(n)
  integer ierr
  integer matz
  real a(nm,3)
  real w(n)
  real z(nm,n)
!
  if ( n > nm ) then
    write ( *, * ) ' '
    write ( *, * ) 'RT - Fatal error!'
    write ( *, * ) '  N greater than NM.'
    ierr = 10 * n
    return
  end if

  if ( matz == 0 ) then

    call figi ( nm, n, a, w, fv1, fv1, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RT - Fatal error!'
      write ( *, * ) '  Error return from FIGI.'
      return
    end if

    call imtql1 ( n, w, fv1, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RT - Fatal error!'
      write ( *, * ) '  Error return from IMTQL1.'
      return
    end if

  else

    call figi2 ( nm, n, a, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RT - Fatal error!'
      write ( *, * ) '  Error return from FIGI2.'
      return
    end if

    call imtql2 ( nm, n, w, fv1, z, ierr )

    if ( ierr /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RT - Fatal error!'
      write ( *, * ) '  Error return from IMTQL2.'
      return
    end if

  end if

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
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine rvec2_print ( n, a1, a2, title )
!
!*******************************************************************************
!
!! RVEC2_PRINT prints a pair of real vectors, with an optional title.
!
!
!  Modified:
!
!    14 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  real a1(n)
  real a2(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,2g14.6)' ) i, a1(i), a2(i)
  end do

  return
end

subroutine svd ( nm, m, n, a, w, matu, u, matv, v, ierr )
!
!*******************************************************************************
!
!! SVD computes the singular value decomposition for a real matrix.
!
!
!  Discussion:
!
!    This subroutine is a translation of the ALGOL procedure SVD.
!
!  Reference:
!
!    Golub and Reinsch,
!    Numerische Mathematik,
!    Volume 14, 1970, pages 403-420.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!     this subroutine determines the singular value decomposition
!          t
!     a=usv  of a real m by n rectangular matrix.  Householder
!     bidiagonalization and a variant of the qr algorithm are used.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.  note that nm must be at least
!          as large as the maximum of m and n.
!
!        m is the number of rows of a (and u).
!
!        n is the number of columns of a (and u) and the order of v.
!
!        a contains the rectangular input matrix to be decomposed.
!
!        matu should be set to .true. if the u matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!        matv should be set to .true. if the v matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!     on output
!
!        a is unaltered (unless overwritten by u or v).
!
!        w contains the n (non-negative) singular values of a (the
!          diagonal elements of s).  they are unordered.  if an
!          error exit is made, the singular values should be correct
!          for indices ierr+1,ierr+2,...,n.
!
!        u contains the matrix u (orthogonal column vectors) of the
!          decomposition if matu has been set to .true.  otherwise
!          u is used as a temporary array.  u may coincide with a.
!          if an error exit is made, the columns of u corresponding
!          to indices of correct singular values should be correct.
!
!        v contains the matrix v (orthogonal) of the decomposition if
!          matv has been set to .true.  otherwise v is not referenced.
!          v may also coincide with a if u is not needed.  if an error
!          exit is made, the columns of v corresponding to indices of
!          correct singular values should be correct.
!
!        ierr is set to
!          zero       for normal return,
!          k          if the k-th singular value has not been
!                     determined after 30 iterations.
!
  integer n
  integer nm
!
  real a(nm,n)
  real c
  real f
  real g
  real h
  integer i
  integer ierr
  integer its
  integer i1
  integer j
  integer k
  integer kk
  integer k1
  integer l
  integer ll
  integer l1
  integer m
  logical matu
  logical matv
  integer mn
  real pythag
  real rv1(n)
  real s
  real scale
  real tst1
  real tst2
  real u(nm,n)
  real v(nm,n)
  real w(n)
  real x
  real y
  real z
!
  ierr = 0
  u(1:m,1:n) = a(1:m,1:n)
!
!  Householder reduction to bidiagonal form.
!
  g = 0.0
  scale = 0.0
  x = 0.0

  do i = 1, n

    l = i + 1
    rv1(i) = scale * g
    g = 0.0
    s = 0.0
    scale = 0.0

    if ( i <= m ) then

      do k = i, m
        scale = scale + abs ( u(k,i) )
      end do

      if ( scale /= 0.0 ) then

        do k = i, m
          u(k,i) = u(k,i) / scale
          s = s + u(k,i)**2
        end do

        f = u(i,i)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        u(i,i) = f - g

        if ( i /= n ) then

          do j = l, n

            s = 0.0
            do k = i, m
              s = s + u(k,i) * u(k,j)
            end do

            f = s / h

            u(i:m,j) = u(i:m,j) + f * u(i:m,i)

          end do

        end if

        u(i:m,i) = scale * u(i:m,i)

      end if

    end if

    w(i) = scale * g
    g = 0.0
    s = 0.0
    scale = 0.0

    if ( i <= m .and. i /= n ) then

    do k = l, n
      scale = scale + abs ( u(i,k) )
    end do

    if ( scale /= 0.0 ) then

      do k = l, n
        u(i,k) = u(i,k) / scale
        s = s + u(i,k)**2
      end do

      f = u(i,l)
      g = -sign(sqrt(s),f)
      h = f * g - s
      u(i,l) = f - g
      rv1(l:n) = u(i,l:n) / h

      if ( i /= m ) then

        do j = l, m

          s = 0.0
          do k = l, n
            s = s + u(j,k) * u(i,k)
          end do

          u(j,l:n) = u(j,l:n) + s * rv1(l:n)

        end do

      end if

      u(i,l:n) = scale * u(i,l:n)

    end if

    end if

    x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

  end do
!
!  Accumulation of right-hand transformations.
!
  if ( matv ) then

    do i = n, 1, -1

      if ( i /= n ) then

         if ( g /= 0.0 ) then

          v(l:n,i) = ( u(i,l:n) / u(i,l) ) / g

          do j = l, n

            s = 0.0
            do k = l, n
              s = s + u(i,k) * v(k,j)
            end do

            v(l:n,j) = v(l:n,j) + s * v(l:n,i)

          end do

        end if

        v(i,l:n) = 0.0
        v(l:n,i) = 0.0

      end if

      v(i,i) = 1.0
      g = rv1(i)
      l = i

    end do

  end if
!
!  Accumulation of left-hand transformations.
!
  if ( matu ) then

    mn = min ( m, n )

    do i = min ( m, n ), 1, -1

      l = i + 1
      g = w(i)

      if ( i /= n ) then
        u(i,l:n) = 0.0
      end if

      if ( g /= 0.0 ) then

        if ( i /= mn ) then

          do j = l, n

            s = 0.0
            do k = l, m
              s = s + u(k,i) * u(k,j)
            end do

            f = (s / u(i,i)) / g

            u(i:m,j) = u(i:m,j) + f * u(i:m,i)

          end do

        end if

        u(i:m,i) = u(i:m,i) / g

      else

        u(i:m,i) = 0.0

      end if

      u(i,i) = u(i,i) + 1.0

    end do

  end if
!
!  Diagonalization of the bidiagonal form.
!
  tst1 = x

  do kk = 1, n

     k1 = n - kk
     k = k1 + 1
     its = 0
!
!  Test for splitting.
!
  520    continue

     do ll = 1, k

        l1 = k - ll
        l = l1 + 1
        tst2 = tst1 + abs(rv1(l))
        if (tst2 == tst1) go to 565
        tst2 = tst1 + abs(w(l1))
        if (tst2 == tst1) go to 540

     end do
!
!  Cancellation of rv1(l) if l greater than 1
!
  540    continue

     c = 0.0
     s = 1.0

     do i = l, k

       f = s * rv1(i)
       rv1(i) = c * rv1(i)
       tst2 = tst1 + abs(f)

       if ( tst2 == tst1 ) then
         go to 565
       end if

       g = w(i)
       h = pythag(f,g)
       w(i) = h
       c = g / h
       s = -f / h

       if ( matu ) then

         do j = 1, m
           y = u(j,l1)
           z = u(j,i)
           u(j,l1) = y * c + z * s
           u(j,i) = -y * s + z * c
         end do

       end if

    end do
!
!  Test for convergence.
!
  565   continue
 
    z = w(k)

     if (l == k) go to 650
!
!  Shift from bottom 2 by 2 minor.
!
     if ( its >= 30 ) then
       ierr = k
       return
     end if

     its = its + 1
     x = w(l)
     y = w(k1)
     g = rv1(k1)
     h = rv1(k)
     f = 0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
     g = pythag(f,1.0)
     f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h)
!
!  Next qr transformation.
!
     c = 1.0
     s = 1.0

     do i1 = l, k1

        i = i1 + 1
        g = rv1(i)
        y = w(i)
        h = s * g
        g = c * g
        z = pythag(f,h)
        rv1(i1) = z
        c = f / z
        s = h / z
        f = x * c + g * s
        g = -x * s + g * c
        h = y * s
        y = y * c

        if ( matv ) then

          do j = 1, n
            x = v(j,i1)
            z = v(j,i)
            v(j,i1) = x * c + z * s
            v(j,i) = -x * s + z * c
          end do

        end if

        z = pythag(f,h)
        w(i1) = z
!
!  Rotation can be arbitrary if z is zero.
!
        if ( z /= 0.0 ) then
          c = f / z
          s = h / z
        end if

  580       continue

        f = c * g + s * y
        x = -s * g + c * y

        if ( matu ) then

          do j = 1, m
            y = u(j,i1)
            z = u(j,i)
            u(j,i1) = y * c + z * s
            u(j,i) = -y * s + z * c
          end do

        end if

     end do

     rv1(l) = 0.0
     rv1(k) = f
     w(k) = x
     go to 520
!
!  Convergence
!
  650    continue

    if ( z <= 0.0 ) then

       w(k) = - z

       if ( matv ) then
         v(1:n,k) = - v(1:n,k)
       end if

     end if

  end do

  return
end
subroutine tinvit ( nm, n, d, e, e2, m, w, ind, z, ierr )
!
!*******************************************************************************
!
!! TINVIT computes eigenvectors from eigenvalues, real tridiagonal symmetric.
!
!
!  Discussion:
!
!     this subroutine is a translation of the inverse iteration tech-
!     nique in the algol procedure tristurm by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvectors of a tridiagonal
!     symmetric matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e,
!          with zeros corresponding to negligible elements of e.
!          e(i) is considered negligible if it is not larger than
!          the product of the relative machine precision and the sum
!          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
!          0.0 if the eigenvalues are in ascending order, or 2.0
!          if the eigenvalues are in descending order.  if  bisect,
!          TRIDIB, or  imtqlv  has been used to find the eigenvalues,
!          their output e2 array is exactly what is expected here.
!
!        m is the number of specified eigenvalues.
!
!        w contains the m eigenvalues in ascending or descending order.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w:
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc.
!
!     on output
!
!        all input arrays are unaltered.
!
!        z contains the associated set of orthonormal eigenvectors.
!          any vector which fails to converge is set to zero.
!
!        ierr is set to
!          zero       for normal return,
!          -r         if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge in 5 iterations.
!
  integer m
  integer n
  integer nm
!
  real d(n)
  real e(n)
  real e2(n)
  real eps2
  real eps3
  real eps4
  integer group
  integer i
  integer ierr
  integer ii
  integer ind(m)
  integer ip
  integer its
  integer j
  integer jj
  real norm
  real order
  integer p
  real pythag
  integer q
  integer r
  real rv1(n)
  real rv2(n)
  real rv3(n)
  real rv4(n)
  real rv6(n)
  integer s
  integer tag
  real u
  real uk
  real v
  real w(m)
  real x0
  real x1
  real xu
  real z(nm,m)
!
  ierr = 0

  if ( m == 0 ) then
    return
  end if

  u = 0.0
  x0 = 0.0

  tag = 0
  order = 1.0 - e2(1)
  q = 0
!
!  Establish and process next submatrix
!
  100 p = q + 1

  do q = p, n
     if (q == n) go to 140
     if (e2(q+1) == 0.0) go to 140
  end do
!
!  Find vectors by inverse iteration
!
  140 tag = tag + 1
  s = 0

  do r = 1, m

     if (ind(r) /= tag) go to 920
     its = 1
     x1 = w(r)
     if (s /= 0) go to 510
!
!  Check for isolated root
!
     xu = 1.0
     if (p /= q) go to 490
     rv6(p) = 1.0
     go to 870
  490    norm = abs(d(p))
     ip = p + 1

     do i = p+1, q
       norm = max(norm, abs(d(i))+abs(e(i)))
     end do
!
!  EPS2 is the criterion for grouping,
!  eps3 replaces zero pivots and equal roots are modified by eps3,
!  eps4 is taken very small to avoid overflow
!
     eps2 = 1.0e-3 * norm
     eps3 = abs ( norm ) * epsilon ( 1.0 )
     uk = q - p + 1
     eps4 = uk * eps3
     uk = eps4 / sqrt(uk)
     s = p
  505    group = 0
     go to 520
!
!  Look for close or coincident roots
!
  510    if (abs(x1-x0) >= eps2) go to 505
     group = group + 1
     if (order * (x1 - x0) <= 0.0) x1 = x0 + order * eps3
!
!  Elimination with interchanges and initialization of vector
!
  520    v = 0.0

     do i = p, q

        rv6(i) = uk

        if (i == p) go to 560

        if (abs(e(i)) < abs(u)) go to 540

        xu = u / e(i)
        rv4(i) = xu
        rv1(i-1) = e(i)
        rv2(i-1) = d(i) - x1
        rv3(i-1) = 0.0
        if (i /= q) rv3(i-1) = e(i+1)
        u = v - xu * rv2(i-1)
        v = -xu * rv3(i-1)
        go to 580

  540   continue

        xu = e(i) / u
        rv4(i) = xu
        rv1(i-1) = u
        rv2(i-1) = v
        rv3(i-1) = 0.0

  560   continue

        u = d(i) - x1 - xu * v
        if (i /= q) v = e(i+1)

  580   continue

     end do

     if (u == 0.0) u = eps3
     rv1(q) = u
     rv2(q) = 0.0
     rv3(q) = 0.0
!
!  Back substitution
!
  600   continue

  do ii = p, q
    i = p + q - ii
    rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
    v = u
    u = rv6(i)
  end do
!
!  Orthogonalize with respect to previous members of group
!
     j = r

     do jj = 1, group

  630       j = j - 1
        if (ind(j) /= tag) go to 630
        xu = 0.0

        do i = p, q
          xu = xu + rv6(i) * z(i,j)
        end do

        rv6(p:q) = rv6(p:q) - xu * z(p:q,j)

     end do

     norm = 0.0
     do i = p, q
       norm = norm + abs(rv6(i))
     end do

     if (norm >= 1.0) go to 840
!
!  Forward substitution
!
     if (its == 5) go to 830
     if (norm /= 0.0) go to 740
     rv6(s) = eps4
     s = s + 1
     if (s > q) s = p
     go to 780
  740    xu = eps4 / norm

     rv6(p:q) = rv6(p:q) * xu
!
!  Elimination operations on next vector iterate
!
  780 continue
!
!  If rv1(i-1) == e(i), a row interchange was performed earlier in the
!  triangularization process.
!
     do i = ip, q

       u = rv6(i)

       if (rv1(i-1) == e(i)) then
         u = rv6(i-1)
         rv6(i-1) = rv6(i)
       end if

       rv6(i) = u - rv4(i) * rv6(i-1)

     end do

     its = its + 1
     go to 600
!
!  Set error: non-converged eigenvector
!
  830    ierr = -r
     xu = 0.0
     go to 870
!
!  Normalize so that sum of squares is 1 and expand to full order
!
  840   continue

     u = 0.0
     do i = p, q
       u = pythag(u,rv6(i))
     end do

     xu = 1.0 / u

  870    continue

     z(1:n,r) = 0.0
     z(p:q,r) = rv6(p:q) * xu

     x0 = x1

  920 continue

  end do

  if (q < n) go to 100

  return
end
subroutine tql1 ( n, d, e, ierr )
!
!*******************************************************************************
!
!! TQL1 computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric tridiagonal 
!    matrix by the QL method.
!
!  References:
!
!    Bowdler, Martin, Reinsch, Wilkinson,
!    Numerische Mathematik,
!    Volume 11, 1968, pages 293-306.
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, is the order of the matrix.
!
!    Input/output, real D(N).
!    On input, the diagonal elements of the matrix.
!    On output, the eigenvalues in ascending order.
!    If an error exit is made, the eigenvalues are correct and
!    ordered for indices 1, 2,... IERR-1, but may not be
!    the smallest eigenvalues.
!
!    Input/output, real E(N).  On input, E(2:N) contains the subdiagonal 
!    elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Output, integer IERR, error flag.
!    0, normal return,
!    J, if the J-th eigenvalue has not been determined after 
!    30 iterations.
!
  integer n
!
  real c
  real c2
  real c3
  real d(n)
  real dl1
  real e(n)
  real el1
  real f
  real g
  real h
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer l1
  integer l2
  integer m
  integer mml
  real p
  real pythag
  real r
  real s
  real s2
  real tst1
  real tst2
!
  ierr = 0
  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do

  f = 0.0
  tst1 = 0.0
  e(n) = 0.0

  do l = 1, n

    j = 0
    h = abs ( d(l) ) + abs ( e(l) )
    tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
    do m = l, n

      tst2 = tst1 + abs ( e(m) )

      if ( tst2 == tst1 ) then
        exit
      end if

    end do

    if ( m == l ) go to 210

  130   continue

    if ( j >= 30 ) then
      ierr = l
      return
    end if

    j = j + 1
!
!  Form the shift.
!
    l1 = l + 1
    l2 = l1 + 1
    g = d(l)
    p = ( d(l1) - g ) / ( 2.0 * e(l) )
    r = pythag ( p, 1.0 )
    d(l) = e(l) / ( p + sign ( r, p ) )
    d(l1) = e(l) * ( p + sign ( r, p ) )
    dl1 = d(l1)
    h = g - d(l)

    d(l2:n) = d(l2:n) - h

    f = f + h
!
!  QL transformation.
!
    p = d(m)
    c = 1.0
    c2 = c
    el1 = e(l1)
    s = 0.0
    mml = m - l

    do ii = 1, mml
      c3 = c2
      c2 = c
      s2 = s
      i = m - ii
      g = c * e(i)
      h = c * p
      r = pythag ( p, e(i) )
      e(i+1) = s * r
      s = e(i) / r
      c = p / r
      p = c * d(i) - s * g
      d(i+1) = h + s * (c * g + s * d(i))
    end do

    p = - s * s2 * c3 * el1 * e(l) / dl1
    e(l) = s * p
    d(l) = c * p
    tst2 = tst1 + abs ( e(l) )
    if ( tst2 > tst1 ) go to 130

  210   continue

    p = d(l) + f
!
!  Order the eigenvalues.
!
    do ii = 2, l
      i = l + 2 - ii
      if ( p >= d(i-1) ) go to 270
      d(i) = d(i-1)
    end do

    i = 1

  270   continue

    d(i) = p

  end do

  return
end
subroutine tql2 ( nm, n, d, e, z, ierr )
!
!*******************************************************************************
!
!! TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
  integer n
  integer nm
!
  real c
  real c2
  real c3
  real d(n)
  real dl1
  real e(n)
  real el1
  real f
  real g
  real h
  integer i
  integer ierr
  integer ii
  integer j
  integer k
  integer l
  integer l1
  integer l2
  integer m
  integer mml
  real p
  real pythag
  real r
  real s
  real s2
  real tst1
  real tst2
  real z(nm,n)
!
  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do

  f = 0.0
  tst1 = 0.0
  e(n) = 0.0

  do l = 1, n

     j = 0
     h = abs(d(l)) + abs(e(l))
     if (tst1 < h) tst1 = h
!
!  look for small sub-diagonal element
!
     do m = l, n
       tst2 = tst1 + abs(e(m))
       if ( tst2 == tst1 ) then
         exit
       end if
     end do

     if (m == l) go to 220

  130 continue

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift
!
     l1 = l + 1
     l2 = l1 + 1
     g = d(l)
     p = (d(l1) - g) / (2.0 * e(l))
     r = pythag(p,1.0)
     d(l) = e(l) / (p + sign(r,p))
     d(l1) = e(l) * (p + sign(r,p))
     dl1 = d(l1)
     h = g - d(l)
     d(l2:n) = d(l2:n) - h
     f = f + h
!
!  QL transformation
!
     p = d(m)
     c = 1.0
     c2 = c
     el1 = e(l1)
     s = 0.0
     mml = m - l

     do ii = 1, mml

        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = pythag(p,e(i))
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * (c * g + s * d(i))
!
!  Form vector
!
        do k = 1, n
           h = z(k,i+1)
           z(k,i+1) = s * z(k,i) + c * h
           z(k,i) = c * z(k,i) - s * h
        end do

     end do

     p = -s * s2 * c3 * el1 * e(l) / dl1
     e(l) = s * p
     d(l) = c * p
     tst2 = tst1 + abs(e(l))

     if (tst2 > tst1) then 
       go to 130
     end if

220  continue

     d(l) = d(l) + f

  end do
!
!  Order eigenvalues and eigenvectors
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n

      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if

    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      do j = 1, n
        call r_swap ( z(j,i), z(j,k) )
      end do

    end if

  end do

  return
end
subroutine tqlrat ( n, d, e2, ierr )
!
!*******************************************************************************
!
!! TQLRAT compute all eigenvalues of a real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure tqlrat,
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e2 contains the squares of the subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e2 has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
  integer n
!
  real b
  real c
  real d(n)
  real e2(n)
  real f
  real g
  real h
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer l1
  integer m
  integer mml
  real p
  real pythag
  real r
  real s
  real t
!
  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e2(i-1) = e2(i)
  end do

  f = 0.0
  t = 0.0
  e2(n) = 0.0

  do l = 1, n

     j = 0
     h = abs(d(l)) + sqrt(e2(l))
     if (t > h) go to 105
     t = h
     b = abs ( t ) * epsilon ( 1.0 )
     c = b * b
!
!  Look for small squared sub-diagonal element
!
105    continue

     do m = l, n
       if ( e2(m) <= c ) then
         exit
       end if
     end do

     if (m == l) go to 210

130    continue

     if ( j >= 30 ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift
!
     l1 = l + 1
     s = sqrt(e2(l))
     g = d(l)
     p = (d(l1) - g) / (2.0 * s)
     r = pythag(p,1.0)
     d(l) = s / (p + sign(r,p))
     h = g - d(l)

     do i = l1, n
       d(i) = d(i) - h
     end do

     f = f + h
!
!  Rational ql transformation
!
     g = d(m)
     if (g == 0.0) g = b
     h = g
     s = 0.0
     mml = m - l

     do ii = 1, mml
        i = m - ii
        p = g * h
        r = p + e2(i)
        e2(i+1) = s * r
        s = e2(i) / r
        d(i+1) = h + s * (h + d(i))
        g = d(i) - e2(i) / g
        if (g == 0.0) g = b
        h = g * p / r
     end do

     e2(l) = s * g
     d(l) = h
!
!  Guard against underflow in convergence test
!
     if (h == 0.0) go to 210
     if (abs(e2(l)) <= abs(c/h)) go to 210
     e2(l) = h * e2(l)
     if (e2(l) /= 0.0) go to 130
210  continue

     p = d(l) + f
!
!  Order eigenvalues
!
     if (l == 1) go to 250

     do ii = 2, l
        i = l + 2 - ii
        if (p >= d(i-1)) go to 270
        d(i) = d(i-1)
     end do

  250    i = 1
  270    d(i) = p

  290 continue

  end do

  return
end
subroutine trbak1 ( nm, n, a, e, m, z ) 
!
!*******************************************************************************
!
!! TRBAK1 determines eigenvectors by undoing the TRED1 transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure trbak1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a real symmetric
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  tred1.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  tred1
!          in its strict lower triangle.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is arbitrary.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
  integer m
  integer n
  integer nm
!
  real a(nm,n)
  real e(n)
  integer i
  integer j
  integer k
  integer l
  real s
  real z(nm,m)
!
  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if

  do i = 2, n

     l = i - 1
     if (e(i) == 0.0) go to 140

     do j = 1, m

        s = 0.0
        do k = 1, l
          s = s + a(i,k) * z(k,j)
        end do

        s = ( s / a(i,l) ) / e(i)

        z(1:l,j) = z(1:l,j) + s * a(i,1:l)

     end do

  140 continue

  end do

  continue

  return
end
subroutine trbak3 ( nm, n, nv, a, m, z )
!
!*******************************************************************************
!
!! TRBAK3 determines eigenvectors by undoing the TRED3 transformation.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure trbak3,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a real symmetric
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  tred3.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        nv must be set to the dimension of the array parameter a
!          as declared in the calling program dimension statement.
!
!        a contains information about the orthogonal transformations
!          used in the reduction by  tred3  in its first
!          n*(n+1)/2 positions.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
  integer m
  integer nm
  integer nv
!
  real a(nv)
  real h
  integer i
  integer ik
  integer iz
  integer j
  integer k
  integer l
  integer n
  real s
  real z(nm,m)
!
  if ( m == 0 ) then
    return
  end if

  do i = 2, n

    l = i - 1
    iz = (i * l) / 2
    ik = iz + i
    h = a(ik)

    if ( h /= 0.0 ) then

      do j = 1, m

        s = 0.0
        ik = iz

        do k = 1, l
          ik = ik + 1
          s = s + a(ik) * z(k,j)
        end do

        s = (s / h) / h
        ik = iz

        do k = 1, l
          ik = ik + 1
          z(k,j) = z(k,j) - s * a(ik)
        end do

      end do

    end if

  end do

  return
end
subroutine tred1 ( nm, n, a, d, e, e2 )
!
!*******************************************************************************
!
!! TRED1 transforms a real symmetric matrix to tridiagonal form.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure tred1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
  integer n
  integer nm
!
  real a(nm,n)
  real d(n)
  real e(n)
  real e2(n)
  real f
  real g
  real h
  integer i
  integer ii
  integer j
  integer k
  integer l
  real scale
!
  d(1:n) = a(n,1:n)

  do i = 1, n
    a(n,i) = a(i,i)
  end do

  do ii = 1, n

     i = n + 1 - ii
     l = i - 1
     h = 0.0
     scale = 0.0

     if ( l < 1 ) go to 130
!
!  Scale row
!
     do k = 1, l
       scale = scale + abs(d(k))
     end do

     if ( scale /= 0.0 ) go to 140

     do j = 1, l
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = 0.0
     end do

  130    e(i) = 0.0
     e2(i) = 0.0
     go to 300

  140 continue

     do k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
     end do

     e2(i) = scale * scale * h
     f = d(l)
     g = -sign(sqrt(h),f)
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g

     if (l == 1) go to 285
!
!  Form a*u
!
     e(1:l) = 0.0

     do j = 1, l

        f = d(j)
        g = e(j) + a(j,j) * f

        do k = j+1, l
           g = g + a(k,j) * d(k)
           e(k) = e(k) + a(k,j) * f
        end do

        e(j) = g

     end do
!
!    form p
!
     f = 0.0

     do j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
     end do

     h = f / (h + h)
!
!  Form Q.
!
     e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
     do j = 1, l

        f = d(j)
        g = e(j)

        do k = j, l
          a(k,j) = a(k,j) - f * e(k) - g * d(k)
        end do

     end do

  285 continue

     do j = 1, l
        f = d(j)
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = f * scale
     end do

  300 continue

  end do

  return
end
subroutine tred2 ( nm, n, a, d, e, z )
!
!*******************************************************************************
!
!! TRED2 transforms a real symmetric matrix to tridiagonal form.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
  integer n
  integer nm
!
  real a(nm,n)
  real d(n)
  real e(n)
  real f
  real g
  real h
  real hh
  integer i
  integer ii
  integer j
  integer k
  integer l
  real scale
  real z(nm,n)
!
  do i = 1, n
    z(i:n,i) = a(i:n,i)
  end do

  d(1:n) = a(n,1:n)

  if (n == 1) go to 510

  do ii = 2, n

     i = n + 2 - ii
     l = i - 1
     h = 0.0
     scale = 0.0
     if (l < 2) go to 130
!
!  Scale row.
!
     do k = 1, l
       scale = scale + abs(d(k))
     end do

     if (scale /= 0.0) go to 140
  130    e(i) = d(l)

     do j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0
        z(j,i) = 0.0
     end do

     go to 290

  140    continue

     do k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
     end do

     f = d(l)
     g = -sign(sqrt(h),f)
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
!
!  Form a*u
!
     e(1:l) = 0.0

     do j = 1, l

        f = d(j)
        z(j,i) = f
        g = e(j) + z(j,j) * f

        do k = j+1, l
           g = g + z(k,j) * d(k)
           e(k) = e(k) + z(k,j) * f
        end do

        e(j) = g
  end do
!
!  Form p
!
     f = 0.0

     do j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
     end do

     hh = f / (h + h)
!
!  Form Q.
!
     e(1:l) = e(1:l) - hh * d(1:l)
!
!    form reduced a
!
     do j = 1, l

        f = d(j)
        g = e(j)

        do k = j, l
          z(k,j) = z(k,j) - f * e(k) - g * d(k)
        end do

        d(j) = z(l,j)
        z(i,j) = 0.0

     end do

  290    d(i) = h
  300 continue

  end do
!
!  Accumulation of transformation matrices
!
  do i = 2, n

     l = i - 1
     z(n,l) = z(l,l)
     z(l,l) = 1.0
     h = d(i)
     if (h == 0.0) go to 380

     d(1:l) = z(1:l,i) / h

     do j = 1, l

        g = 0.0

        do k = 1, l
          g = g + z(k,i) * z(k,j)
        end do

        do k = 1, l
           z(k,j) = z(k,j) - g * d(k)
        end do

     end do

  380 continue

     z(1:l,i) = 0.0

  500 continue

  end do

  510 continue

  d(1:n) = z(n,1:n)

  z(n,1:n-1) = 0.0
  z(n,n) = 1.0

  e(1) = 0.0

  return
end
subroutine tred3 ( n, nv, a, d, e, e2 )
!
!*******************************************************************************
!
!! TRED3 transforms a real symmetric packed matrix to tridiagonal form.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure tred3,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix, stored as
!     a one-dimensional array, to a symmetric tridiagonal matrix
!     using orthogonal similarity transformations.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        n is the order of the matrix.
!
!        nv must be set to the dimension of the array parameter a
!          as declared in the calling program dimension statement.
!
!        a contains the lower triangle of the real symmetric
!          input matrix, stored row-wise as a one-dimensional
!          array, in its first n*(n+1)/2 positions.
!
!     on output
!
!        a contains information about the orthogonal
!          transformations used in the reduction.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
  integer n
  integer nv
!
  real a(nv)
  real d(n)
  real e(n)
  real e2(n)
  real f
  real g
  real h
  real hh
  integer i
  integer ii
  integer iz
  integer j
  integer jk
  integer k
  integer l
  real scale
!
  do ii = 1, n

     i = n + 1 - ii
     l = i - 1
     iz = (i * l) / 2
     h = 0.0
     scale = 0.0
     if (l < 1) go to 130
!
!  Scale row.
!
     do k = 1, l
        iz = iz + 1
        d(k) = a(iz)
        scale = scale + abs(d(k))
     end do

     if (scale /= 0.0) go to 140
  130    e(i) = 0.0
     e2(i) = 0.0
     go to 290

  140    continue

     do k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
     end do

     e2(i) = scale * scale * h
     f = d(l)
     g = -sign(sqrt(h),f)
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
     a(iz) = scale * d(l)

     if ( l == 1 ) go to 290

     jk = 1

     do j = 1, l

        f = d(j)
        g = 0.0

        do k = 1, j-1
           g = g + a(jk) * d(k)
           e(k) = e(k) + a(jk) * f
           jk = jk + 1
        end do

        e(j) = g + a(jk) * f
        jk = jk + 1

     end do
!
!  Form p
!
     f = 0.0

     do j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
     end do

     hh = f / (h + h)
!
!  Form q
!
     do j = 1, l
       e(j) = e(j) - hh * d(j)
     end do

     jk = 1
!
!  Form reduced a
!
     do j = 1, l
        f = d(j)
        g = e(j)
        do k = 1, j
           a(jk) = a(jk) - f * e(k) - g * d(k)
           jk = jk + 1
        end do
     end do

  290    d(i) = a(iz+1)
     a(iz+1) = scale * sqrt(h)
  300 continue

  end do

  return
end
subroutine tridib ( n, eps1, d, e, e2, lb, ub, m11, m, w, ind, ierr )
!
!*******************************************************************************
!
!! TRIDIB computes some eigenvalues of a real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure bisect,
!     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
!
!     this subroutine finds those eigenvalues of a tridiagonal
!     symmetric matrix between specified boundary indices,
!     using bisection.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        n is the order of the matrix.
!
!        eps1 is an absolute error tolerance for the computed
!          eigenvalues.  if the input eps1 is non-positive,
!          it is reset for each submatrix to a default value,
!          namely, minus the product of the relative machine
!          precision and the 1-norm of the submatrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        m11 specifies the lower boundary index for the desired
!          eigenvalues.
!
!        m specifies the number of eigenvalues desired.  the upper
!          boundary index m22 is then obtained as m22=m11+m-1.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        lb and ub define an interval containing exactly the desired
!          eigenvalues.
!
!        w contains, in its first m positions, the eigenvalues
!          between indices m11 and m22 in ascending order.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w:
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc..
!
!        ierr is set to
!          zero       for normal return,
!          3*n+1      if multiple eigenvalues at index m11 make
!                     unique selection impossible,
!          3*n+2      if multiple eigenvalues at index m22 make
!                     unique selection impossible.
!
  integer m
  integer n
!
  real d(n)
  real e(n)
  real e2(n)
  real eps1
  integer i
  integer ierr
  integer ii
  integer ind(m)
  integer isturm
  integer j
  integer k
  integer l
  real lb
  integer m1
  integer m11
  integer m2
  integer m22
  integer p
  integer q
  integer r
  real rv4(n)
  real rv5(n)
  integer s
  real t1
  real t2
  integer tag
  real tst1
  real tst2
  real u
  real ub
  real v
  real w(m)
  real x0
  real x1
  real xu
!
  ierr = 0
  tag = 0
  xu = d(1)
  x0 = d(1)
  s = 0
  u = 0.0
!
!  Look for small sub-diagonal entries and determine an
!  interval containing all the eigenvalues
!
  do i = 1, n

     x1 = u

     if ( i == n ) then
       u = 0.0
     else
       u = abs(e(i+1))
     end if

     xu = min(d(i)-(x1+u),xu)
     x0 = max(d(i)+(x1+u),x0)

     if ( i >= 1 ) then
       tst1 = abs(d(i)) + abs(d(i-1))
       tst2 = tst1 + abs(e(i))
       if ( tst2 <= tst1 ) then
         e2(i) = 0.0
       end if
     else
       e2(i) = 0.0
     end if

  end do

  x1 = n
  x1 = x1 * max ( abs ( xu ), abs ( x0 ) ) * epsilon ( 1.0 )
  xu = xu - x1
  t1 = xu
  x0 = x0 + x1
  t2 = x0
!
!  Determine an interval containing exactly the desired eigenvalues
!
  p = 1
  q = n
  m1 = m11 - 1
  if (m1 == 0) go to 75
  isturm = 1
   50 v = x1
  x1 = xu + (x0 - xu) * 0.5
  if (x1 == v) go to 980
  go to 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
  go to 50
   70 x0 = x1
  go to 50
   73 xu = x1
  t1 = x1
   75 m22 = m1 + m
  if (m22 == n) go to 90
  x0 = t2
  isturm = 2
  go to 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
  r = 0
!
!  Establish and process next submatrix, refining interval by the 
!  gerschgorin bounds
!
  100 continue

  if (r == m) then
    go to 1001
  end if

  tag = tag + 1
  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = 0.0

  do q = p, n
    x1 = u
    u = 0.0
    v = 0.0
    if (q == n) go to 110
    u = abs(e(q+1))
    v = e2(q+1)
  110    xu = min(d(q)-(x1+u),xu)
    x0 = max(d(q)+(x1+u),x0)
    if (v == 0.0) go to 140
  end do

  140 continue

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( 1.0 )
  if (eps1 <= 0.0) eps1 = -x1
  if (p /= q) go to 180
!
!  Check for isolated root within interval
!
  if (t1 > d(p) .or. d(p) >= t2) go to 940
  m1 = p
  m2 = p
  rv5(p) = d(p)
  go to 900
  180 x1 = x1 * (q - p + 1)
  lb = max(t1,xu-x1)
  ub = min(t2,x0+x1)
  x1 = lb
  isturm = 3
  go to 320
  200 m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320
  220 m2 = s
  if (m1 > m2) go to 940
!
!  Find roots by bisection
!
  x0 = ub
  isturm = 5

  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for k-th eigenvalue
!
  k = m2
  250    xu = lb

     do ii = m1, k

        i = m1 + k - ii
        if (xu >= rv4(i)) go to 260
        xu = rv4(i)
        go to 280

  260    continue

      end do

  280    if (x0 > rv5(k)) x0 = rv5(k)
!
!  Next bisection step
!
  300    x1 = (xu + x0) * 0.5
     if ((x0 - xu) <= abs(eps1)) go to 420
     tst1 = 2.0 * (abs(xu) + abs(x0))
     tst2 = tst1 + (x0 - xu)
     if (tst2 == tst1) go to 420
!
!  Sturm sequence
!
  320    s = p - 1
     u = 1.0

     do i = p, q
        if (u /= 0.0) go to 325
        v = abs(e(i)) / epsilon ( 1.0 )
        if (e2(i) == 0.0) v = 0.0
        go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
        if (u < 0.0) s = s + 1
     end do

     go to (60,80,200,220,360), isturm
!
!  Refine intervals
!
  360    if (s >= k) go to 400
     xu = x1
     if (s >= m1) go to 380
     rv4(m1) = x1
     go to 300
  380    rv4(s+1) = x1
     if (rv5(s) > x1) rv5(s) = x1
     go to 300
  400    x0 = x1
     go to 300
!
!  K-th eigenvalue found
!
  420    rv5(k) = x1
  k = k - 1
  if (k >= m1) go to 250
!
!  Order eigenvalues tagged with their submatrix associations
!
  900 s = r
  r = r + m2 - m1 + 1
  j = 1
  k = m1

  do l = 1, r

     if (j > s) go to 910
     if (k > m2) go to 940
     if (rv5(k) >= w(l)) go to 915

     do ii = j, s
        i = l + s - ii
        w(i+1) = w(i)
        ind(i+1) = ind(i)
     end do

  910    w(l) = rv5(k)
     ind(l) = tag
     k = k + 1
     go to 920
  915    j = j + 1

  920 continue

  end do

  940 if (q < n) go to 100
  go to 1001
!
!  Set error: interval cannot be found containing exactly the 
!  desired eigenvalues
!
  980 ierr = 3 * n + isturm
 1001 lb = t1
  ub = t2
  return
end
subroutine tsturm ( nm, n, eps1, d, e, e2, lb, ub, mm, m, w, z, ierr )
!
!*******************************************************************************
!
!! TSTURM computes some eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!     this subroutine is a translation of the algol procedure tristurm
!     by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvalues of a tridiagonal
!     symmetric matrix which lie in a specified interval and their
!     associated eigenvectors, using bisection and inverse iteration.
!
!  Reference:
!
!    J H Wilkinson and C Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        eps1 is an absolute error tolerance for the computed
!          eigenvalues.  it should be chosen commensurate with
!          relative perturbations in the matrix elements of the
!          order of the relative machine precision.  if the
!          input eps1 is non-positive, it is reset for each
!          submatrix to a default value, namely, minus the
!          product of the relative machine precision and the
!          1-norm of the submatrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        lb and ub define the interval to be searched for eigenvalues.
!          if lb is not less than ub, no eigenvalues will be found.
!
!        mm should be set to an upper bound for the number of
!          eigenvalues in the interval.  warning. if more than
!          mm eigenvalues are determined to lie in the interval,
!          an error return is made with no values or vectors found.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        m is the number of eigenvalues determined to lie in (lb,ub).
!
!        w contains the m eigenvalues in ascending order if the matrix
!          does not split.  if the matrix splits, the eigenvalues are
!          in ascending order for each submatrix.  if a vector error
!          exit is made, w contains those values already found.
!
!        z contains the associated set of orthonormal eigenvectors.
!          if an error exit is made, z contains those vectors
!          already found.
!
!        ierr is set to
!          zero       for normal return,
!          3*n+1      if m exceeds mm.
!          4*n+r      if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge in 5 iterations.
!
  integer mm
  integer n
  integer nm
!
  real d(n)
  real e(n)
  real e2(n)
  real eps1
  real eps2
  real eps3
  real eps4
  integer group
  integer i
  integer ierr
  integer ii
  integer ip
  integer isturm
  integer its
  integer j
  integer jj
  integer k
  real lb
  integer m
  integer m1
  integer m2
  real norm
  integer p
  real pythag
  integer q
  integer r
  real rv1(n)
  real rv2(n)
  real rv3(n)
  real rv4(n)
  real rv5(n)
  real rv6(n)
  integer s
  real t1
  real t2
  real tst1
  real tst2
  real u
  real ub
  real uk
  real v
  real w(mm)
  real x0
  real x1
  real xu
  real z(nm,mm)
!
  ierr = 0
  s = 0
  t1 = lb
  t2 = ub
!
!  Look for small sub-diagonal entries
!
  do i = 1, n

    if ( i == 1 ) then

      e2(i) = 0.0

    else

      tst1 = abs(d(i)) + abs(d(i-1))
      tst2 = tst1 + abs(e(i))
      if ( tst2 <= tst1 ) then
        e2(i) = 0.0
     end if

    end if


  end do
!
!  Determine the number of eigenvalues in the interval.
!
  p = 1
  q = n
  x1 = ub
  isturm = 1
  go to 320

60 continue

  m = s
  x1 = lb
  isturm = 2
  go to 320

80 continue

  m = m - s

  if (m > mm ) go to 980

  q = 0
  r = 0
!
!  Establish and process next submatrix, refining interval by the 
!  gerschgorin bounds
!
  100 continue

  if (r == m) go to 1001
  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = 0.0

  do q = p, n

     x1 = u
     u = 0.0
     v = 0.0

     if ( q /= n ) then
       u = abs(e(q+1))
       v = e2(q+1)
     end if

     xu = min(d(q)-(x1+u),xu)
     x0 = max(d(q)+(x1+u),x0)

     if ( v == 0.0 ) then
       exit
     end if

  end do

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( 1.0 )
  if (eps1 <= 0.0) eps1 = -x1
  if (p /= q) go to 180
!
!  Check for isolated root within interval
!
  if (t1 > d(p) .or. d(p) >= t2) go to 940
  r = r + 1

  z(1:n,r) = 0.0

  w(r) = d(p)
  z(p,r) = 1.0
  go to 940
  180 u = q-p+1
  x1 = u * x1
  lb = max(t1,xu-x1)
  ub = min(t2,x0+x1)
  x1 = lb
  isturm = 3
  go to 320
  200 m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320
  220 m2 = s
  if (m1 > m2) go to 940
!
!  Find roots by bisection
!
  x0 = ub
  isturm = 5

  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for k-th eigenvalue
!
  k = m2
  250    continue
  xu = lb

  do ii = m1, k

    i = m1 + k - ii

    if ( xu < rv4(i) ) then
      xu = rv4(i)
      go to 280
    end if

  end do

  280    if (x0 > rv5(k)) x0 = rv5(k)
!
!  Next bisection step
!
  300    x1 = (xu + x0) * 0.5
     if ((x0 - xu) <= abs(eps1)) go to 420
     tst1 = 2.0 * (abs(xu) + abs(x0))
     tst2 = tst1 + (x0 - xu)
     if (tst2 == tst1) go to 420
!
!  Sturm sequence
!
  320    s = p - 1
     u = 1.0

     do i = p, q

        if (u /= 0.0) go to 325
        v = abs(e(i)) / epsilon ( 1.0 )
        if (e2(i) == 0.0) v = 0.0
        go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
        if (u < 0.0) s = s + 1

     end do

     go to (60,80,200,220,360), isturm
!
!  Refine intervals
!
  360    if (s >= k) go to 400
     xu = x1
     if (s >= m1) go to 380
     rv4(m1) = x1
     go to 300
  380    rv4(s+1) = x1
     if (rv5(s) > x1) rv5(s) = x1
     go to 300
  400    x0 = x1
     go to 300
!
!  K-th eigenvalue found
!
  420    rv5(k) = x1
  k = k - 1
  if (k >= m1) go to 250
!
!  Find vectors by inverse iteration
!
  norm = abs(d(p))
  ip = p + 1

  do i = ip, q
    norm = max(norm, abs(d(i)) + abs(e(i)))
  end do
!
!  EPS2 is the criterion for grouping,
!  EPS3 replaces zero pivots and equal roots are modified by eps3,
!  EPS4 is taken very small to avoid overflow.
!
  eps2 = 1.0e-3 * norm
  eps3 = abs ( norm ) * epsilon ( 1.0 )
  uk = q - p + 1
  eps4 = uk * eps3
  uk = eps4 / sqrt(uk)
  group = 0
  s = p

  do k = m1, m2

     r = r + 1
     its = 1
     w(r) = rv5(k)
     x1 = rv5(k)
!
!  Look for close or coincident roots
!
     if (k == m1) go to 520
     if (x1 - x0 >= eps2) group = -1
     group = group + 1
     if (x1 <= x0) x1 = x0 + eps3
!
!  Elimination with interchanges and initialization of vector.
!
  520    v = 0.0

     do i = p, q

        rv6(i) = uk
        if (i == p) go to 560
        if (abs(e(i)) < abs(u)) go to 540
        xu = u / e(i)
        rv4(i) = xu
        rv1(i-1) = e(i)
        rv2(i-1) = d(i) - x1
        rv3(i-1) = 0.0
        if (i /= q) rv3(i-1) = e(i+1)
        u = v - xu * rv2(i-1)
        v = -xu * rv3(i-1)
        go to 580
  540       xu = e(i) / u
        rv4(i) = xu
        rv1(i-1) = u
        rv2(i-1) = v
        rv3(i-1) = 0.0
  560       u = d(i) - x1 - xu * v
        if (i /= q) v = e(i+1)

  580    continue

     end do

     if (u == 0.0) u = eps3
     rv1(q) = u
     rv2(q) = 0.0
     rv3(q) = 0.0
!
!  Back substitution.
!
  600    continue

     do ii = p, q
        i = p + q - ii
        rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
        v = u
        u = rv6(i)
     end do
!
!  Orthogonalize with respect to previous members of group.
!
     do jj = 1, group

        j = r - group - 1 + jj
        xu = 0.0

        do i = p, q
          xu = xu + rv6(i) * z(i,j)
        end do

        rv6(p:q) = rv6(p:q) - xu * z(p:q,j)

     end do

  700    norm = 0.0

     do i = p, q
       norm = norm + abs(rv6(i))
     end do

     if (norm >= 1.0) go to 840
!
!  Forward substitution
!
     if (its == 5) go to 960
     if (norm /= 0.0) go to 740
     rv6(s) = eps4
     s = s + 1
     if (s > q) s = p
     go to 780
  740    xu = eps4 / norm

     rv6(p:q) = rv6(p:q) * xu
!
!  Elimination operations on next vector iterate
!
  780    continue
!
!  If rv1(i-1) == e(i), a row interchange was performed earlier in the
!  triangularization process
!
     do i = p, q

       u = rv6(i)

       if ( rv1(i-1) == e(i) ) then
         u = rv6(i-1)
         rv6(i-1) = rv6(i)
       end if

       rv6(i) = u - rv4(i) * rv6(i-1)

     end do

     its = its + 1
     go to 600
!
!  Normalize so that sum of squares is 1 and expand to full order.
!
  840    u = 0.0

     do i = p, q
       u = pythag(u,rv6(i))
     end do

     xu = 1.0 / u

     z(1:n,r) = 0.0
     z(p:q,r) = rv6(p:q) * xu

     x0 = x1
  
  end do

  940 if (q < n) go to 100
  go to 1001
!
!  Set error: non-converged eigenvector
!
  960 ierr = 4 * n + r
  go to 1001
!
!  Set error: underestimate of number of eigenvalues in interval
!
  980 ierr = 3 * n + 1
 1001 lb = t1
  ub = t2
  return
end
