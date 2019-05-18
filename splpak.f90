!  splpak.f90  11 August 2000
!
subroutine banfac ( w, nroww, nrow, nbandl, nbandu, iflag )
!
!*************************************************************************
!
!! BANFAC factors a banded matrix without pivoting.
!
!
!  Discussion:
!
!    BANFAC returns in W the LU-factorization, without pivoting, of 
!    the banded matrix A of order NROW with (NBANDL+1+NBANDU) bands 
!    or diagonals in the work array W.
! 
!    Gauss elimination without pivoting is used.  The routine is 
!    intended for use with matrices A which do not require row 
!    interchanges during factorization, especially for the totally 
!    positive matrices which occur in spline calculations.
!
!  Parameters:
! 
!    Input/output, real W(NROWW,NROW).
!    On input, W contains the "interesting" part of a banded 
!    matrix A, with the diagonals or bands of A stored in the
!    rows of W, while columns of A correspond to columns of W. 
!
!    This is the storage mode used in LINPACK and results in efficient 
!    innermost loops.
! 
!    Explicitly, A has 
! 
!      NBANDL bands below the diagonal
!      1     main diagonal
!      NBANDU bands above the diagonal
!
!    and thus, with MIDDLE=NBANDU+1,
!    A(I+J,J) is in W(I+MIDDLE,J) for I=-NBANDU,...,NBANDL, J=1,...,NROW.
!
!    For example, the interesting entries of a banded matrix
!    matrix of order 9, with NBANDL=1, NBANDU=2:
!
!      11 12 13  0  0  0  0  0  0
!      21 22 23 24  0  0  0  0  0
!       0 32 33 34 35  0  0  0  0
!       0  0 43 44 45 46  0  0  0
!       0  0  0 54 55 56 57  0  0
!       0  0  0  0 65 66 67 68  0
!       0  0  0  0  0 76 77 78 79
!       0  0  0  0  0  0 87 88 89
!       0  0  0  0  0  0  0 98 99
!
!    would appear in the first 1+1+2=4 rows of W as follows:
!
!       0  0 13 24 35 46 57 68 79
!       0 12 23 34 45 56 67 78 89
!      11 22 33 44 55 66 77 88 99
!      21 32 43 54 65 76 87 98  0
! 
!    All other entries of W not identified in this way with an
!    entry of A are never referenced.
! 
!    On output, W contains the LU-factorization of A into a unit 
!    lower triangular matrix L and an upper triangular matrix U 
!    (both banded) and stored in customary fashion over the 
!    corresponding entries of A.  
!
!    This makes it possible to solve any particular linear system A*X=B 
!    for X by the call
!
!      call banslv(w,nroww,nrow,nbandl,nbandu,b)
!
!    with the solution X contained in B on return.
! 
!    If IFLAG=2, then one of NROW-1, NBANDL, NBANDU failed to be nonnegative, 
!    or else one of the potential pivots was found to be zero 
!    indicating that A does not have an LU-factorization.  This 
!    implies that A is singular in case it is totally positive.
! 
!    Input, integer NROWW, the row dimension of the work array W.
!    NROWW must be at least NBANDL+1 + NBANDU.
! 
!    Input, integer NROW, the number of rows in A.
!
!    Input, integer NBANDL, the number of bands of A below the main diagonal.
! 
!    Input, integer NBANDU, the number of bands of A above the main diagonal.
! 
!    Output, integer IFLAG, error flag.
!    1, success.
!    2, failure, the matrix was not factored.
! 
  integer nrow
  integer nroww
!
  real factor
  integer i
  integer iflag
  integer j
  integer k
  integer middle
  integer nbandl
  integer nbandu
  real pivot
  real w(nroww,nrow)
!
  iflag = 1
  middle = nbandu+1
!
!  W(MIDDLE,*) contains the main diagonal of A.
!
  if ( nrow < 1 ) then
    iflag = 2
    return
  end if
  
  if ( nrow == 1 ) then
    if ( w(middle,nrow) == 0.0 ) then
      iflag = 2
    end if
    return
  end if
!
!  A is upper triangular.  Check that the diagonal is nonzero.
!
  if ( nbandl <= 0 ) then

    do i = 1, nrow-1
      if ( w(middle,i) == 0.0 ) then
        iflag = 2
        return
      end if
    end do

    if ( w(middle,nrow) == 0.0 ) then
      iflag = 2
    end if

    return
!
!  A is lower triangular.  Check that the diagonal is nonzero and
!  divide each column by its diagonal.
!
  else if ( nbandu <= 0 ) then

    do i = 1, nrow-1

      pivot = w(middle,i)

      if ( pivot == 0.0 ) then
        iflag = 2
        return
      end if

      do j = 1, min(nbandl,nrow-i)
        w(middle+j,i) = w(middle+j,i)/pivot
      end do

    end do

    return

  end if
!
!  A is not just a triangular matrix.  
!  Construct the LU factorization.
!
  do i = 1, nrow-1
!
!  W(MIDDLE,I) is the pivot for the I-th step.
!
    if ( w(middle,i) == 0.0 ) then
      iflag = 2
      write(*,*)' '
      write(*,*)'BanFac - Fatal error!'
      write(*,*)'  Zero pivot encountered in column ',i
      stop
    end if
!
!  Divide each entry in column I below the diagonal by PIVOT.
!
    do j = 1, min(nbandl,nrow-i)
      w(middle+j,i) = w(middle+j,i)/w(middle,i)
    end do
!
!  Subtract A(I,I+K)*(I-th column) from (I+K)-th column (below
!  row I).
!
    do k = 1, min(nbandu,nrow-i)
      factor = w(middle-k,i+k)
      do j = 1, min(nbandl,nrow-i)
        w(middle-k+j,i+k) = w(middle-k+j,i+k)-w(middle+j,i)*factor
      end do
    end do
 
  end do
!
!  Check the last diagonal entry.
!
  if ( w(middle,nrow) == 0.0 ) then
    iflag = 2
  end if

  return
end
subroutine banslv ( w, nroww, nrow, nbandl, nbandu, b )
!
!*************************************************************************
!
!! BANSLV solves a banded linear system X * X = B factored by BANFAC.
!
!
!  Parameters:
!
!    Input, real W(NROWW,NROW).  W contains the banded matrix,
!    after it has been factored by BANFAC.
!
!    Input, integer NROWW, the row dimension of the work array W.
!    NROWW must be at least NBANDL+1 + NBANDU.
! 
!    Input, integer NROW, the number of rows in A.
!
!    Input, integer NBANDL, the number of bands of A below the 
!    main diagonal.
! 
!    Input, integer NBANDU, the number of bands of A above the 
!    main diagonal.
! 
!    Input/output, real B(NROW).
!    On input, B contains the right hand side of the system to be solved.
!    On output, B contains the solution, X.
!
  integer nrow
  integer nroww
!
  real b(nrow)
  integer i
  integer j
  integer jmax
  integer middle
  integer nbandl
  integer nbandu
  real w(nroww,nrow)
!
  middle = nbandu + 1

  if ( nrow == 1 ) then
    b(1) = b(1) / w(middle,1)
    return
  end if
!
!  Forward pass
!
!  For I=1,2,...,NROW-1, subtract RHS(I)*(I-th column of L) 
!  from the right side, below the I-th row.
!
  if ( nbandl > 0 ) then
    do i = 1, nrow-1
      jmax = min ( nbandl, nrow-i )
      do j = 1, jmax
        b(i+j) = b(i+j) - b(i) * w(middle+j,i)
      end do
    end do
  end if
!
!  Backward pass
!
!  For I=NROW, NROW-1,...,1, divide RHS(I) by 
!  the I-th diagonal entry of U, then subtract 
!  RHS(I)*(I-th column of U) from right side, above the I-th row.
!
  do i = nrow, 2, -1
   
    b(i) = b(i) / w(middle,i)

    do j = 1, min ( nbandu, i-1 )
      b(i-j) = b(i-j) - b(i) * w(middle-j,i)
    end do

  end do

  b(1) = b(1) / w(middle,1)

  return
end
subroutine bchfac ( w, nbands, nrow, diag )
!
!*************************************************************************
!
!! BCHFAC constructs a Cholesky factorization of a matrix.
!
!
!  Discussion:
!
!    The factorization has the form
!
!      C = L * D * L-transpose 
!  
!    with L unit lower triangular and D diagonal, for a given matrix C of 
!    order NROW, where C is symmetric positive semidefinite and banded, 
!    having NBANDS diagonals at and below the main diagonal.
! 
!    Gauss elimination is used, adapted to the symmetry and bandedness of C.
! 
!    Near zero pivots are handled in a special way.  The diagonal 
!    element C(N,N)=W(1,N) is saved initially in DIAG(N), all N. 
! 
!    At the N-th elimination step, the current pivot element, W(1,N), 
!    is compared with its original value, DIAG(N).  If, as the result 
!    of prior elimination steps, this element has been reduced by about 
!    a word length, (i.e., if W(1,N)+DIAG(N) <= DIAG(N)), then the pivot 
!    is declared to be zero, and the entire N-th row is declared to
!    be linearly dependent on the preceding rows.  This has the effect 
!    of producing X(N) = 0 when solving C*X = B for X, regardless of B.
! 
!    Justification for this is as follows.  In contemplated applications 
!    of this program, the given equations are the normal equations for 
!    some least-squares approximation problem, DIAG(N) = C(N,N) gives 
!    the norm-square of the N-th basis function, and, at this point, 
!    W(1,N) contains the norm-square of the error in the least-squares 
!    approximation to the N-th basis function by linear combinations 
!    of the first N-1.  
!
!    Having W(1,N)+DIAG(N) <= DIAG(N) signifies that the N-th function 
!    is linearly dependent to machine accuracy on the first N-1 
!    functions, therefore can safely be left out from the basis of 
!    approximating functions.
!
!    The solution of a linear system C*X=B is effected by the 
!    succession of the following two calls:
! 
!      CALL BCHFAC(W,NBANDS,NROW,DIAG)
!
!      CALL BCHSLV(W,NBANDS,NROW,B,X)
!
!  Parameters:
!
!  W      Input/output, real W(NBANDS,NROW).
!
!         On input, W contains the NBANDS diagonals in its rows, 
!         with the main diagonal in row 1.  Precisely, W(I,J) 
!         contains C(I+J-1,J), I=1,...,NBANDS, J=1,...,NROW.
!
!         For example, the interesting entries of a seven diagonal
!         symmetric matrix C of order 9 would be stored in W as
! 
!           11 22 33 44 55 66 77 88 99
!           21 32 43 54 65 76 87 98
!           31 42 53 64 75 86 97
!           41 52 63 74 85 96
! 
!
!
!  ???  ANOTHER LINE LOST ???
!
!         entry of C are never referenced.
!
!         On output, W contains the Cholesky factorization 
!         C = L*D*L-transp, with W(1,I) containing 1/D(I,I) and W(I,J) 
!         containing L(I-1+J,J), I=2,...,NBANDS.
!
!  NBANDS Input, integer NBANDS, indicates the bandwidth of the
!         matrix C, i.e., C(I,J) = 0 for ABS(I-J) > NBANDS.
! 
!  NROW   Input, integer NROW, is the order of the matrix C.
! 
!  DIAG   Work array, real DIAG(NROW).
! 
  integer nbands
  integer nrow
!
  real diag(nrow)
  integer i
  integer imax
  integer j
  integer jmax
  integer n
  real ratio
  real w(nbands,nrow)
!
  if ( nrow <= 1 ) then
    if ( w(1,1) > 0.0 ) then
      w(1,1) = 1.0 / w(1,1)
    end if
    return
  end if
!
!  Store the diagonal.
!
  do i = 1, nrow
    diag(i) = w(1,i)
  end do
!
!  Factorization.
!
  do n = 1, nrow
 
    if ( w(1,n)+diag(n) <= diag(n) ) then
 
      do j = 1, nbands
        w(j,n) = 0.0
      end do

    else
 
      w(1,n) = 1.0 / w(1,n)
 
      imax = min(nbands-1,nrow-n)
 
      jmax = imax
 
      do i = 1, imax
 
        ratio = w(i+1,n)*w(1,n)
 
        do j = 1, jmax
          w(j,n+i) = w(j,n+i)-w(j+i,n)*ratio
        end do
 
        jmax = jmax-1
        w(i+1,n) = ratio
 
      end do
 
    end if
 
  end do
 
  return
end
subroutine bchslv ( w, nbands, nrow, b )
!
!*************************************************************************
!
!! BCHSLV solves a banded symmetric positive definite system.
!
!
!  Discussion:
!
!    The system is of the form:
!
!      C*X = B 
!  
!    and the Cholesky factorization of C has been constructed 
!    by BCHFAC.
! 
!    With the factorization C = L*D*L-transpose available, where 
!    L is unit lower triangular and D is diagonal, the triangular 
!    system L*Y = B is solved for Y (forward substitution), Y is stored 
!    in  B, the vector  D**(-1)*Y is computed and stored in B, then the 
!    triangular system L-transpose*X = D**(-1)*Y is solved for X 
!    (backsubstitution).
! 
!  Parameters:
!
!    Input, real W(NBANDS,NROW) contains the Cholesky factorization for C, 
!    as computed by BCHFAC.
! 
!    Input, integer NBANDS, the bandwidth of C.
!
!    Input, integer NROW, the order of the matrix C.
! 
!    Input/output, real B(NROW).
!    On input, the right hand side.
!    On output, the solution.
!
  integer nbands
  integer nrow
!
  real b(nrow)
  integer j
  integer n
  real w(nbands,nrow)
!
  if ( nrow <= 1 ) then
    b(1) = b(1) * w(1,1)
    return
  end if
!
!  Forward substitution. 
!  Solve L*Y=B.
!
  do n = 1, nrow

    do j = 1, min(nbands-1,nrow-n)
      b(j+n) = b(j+n)-w(j+1,n)*b(n)
    end do

  end do
!
!  Backsubstitution. 
!  Solve L-transp*X=D**(-1)*Y.
!
  do n = nrow, 1, -1

    b(n) = b(n)*w(1,n)

    do j = 1, min(nbands-1,nrow-n)
      b(n) = b(n)-w(j+1,n)*b(j+n)
    end do

  end do

  return
end
subroutine bsplpp ( t, bcoef, n, k, scrtch, break, coef, l )
!
!*************************************************************************
!
!! BSPLPP converts from B-spline to piecewise polynomial form.
!
!
!  Discussion:
!
!    The B-spline representation of a spline is ( T, BCOEF, N, K ),
!    while the piecewise polynomial representation is 
!    ( BREAK, COEF, L, K ).
!
!    For each breakpoint interval, the K relevant B-spline coefficients 
!    of the spline are found and then differenced repeatedly to get the 
!    B-spline coefficients of all the derivatives of the spline on that 
!    interval. 
!
!    The spline and its first K-1 derivatives are then evaluated at the 
!    left end point of that interval, using BSPLVB repeatedly to obtain 
!    the values of all B-splines of the appropriate order at that point.
! 
!  Parameters:
!
!    Input, real T(N+K), the knot sequence.
! 
!    Input, real BCOEF(N), the B spline coefficient sequence.
! 
!    Input, integer N, the number of B spline coefficients.
! 
!    Input, integer K, the order of the spline.
! 
!    WARNING: the restriction K <= KMAX (= 20) is imposed
!    by the arbitrary dimension statement for BIATX, but
!    is nowhere checked for.
! 
!    Work array, real SCRTCH(K,K).
! 
!    Output, real BREAK(L+1), the piecewise polynomial breakpoint 
!    sequence.  BREAK contains the distinct points in the 
!    sequence T(K),...,T(N+1)
! 
!    Output, real COEF(K,N), with COEF(I,J) = (I-1)st derivative 
!    of the spline at BREAK(J) from the right.
! 
!    Output, integer L, the number of polynomial pieces which 
!    make up the spline in the interval (T(K),T(N+1)).
!
  integer, parameter :: kmax = 20
!
  integer k
  integer l
  integer n
!
  real bcoef(n)
  real biatx(kmax)
  real break(*)
  real coef(k,n)
  real diff
  integer i
  integer j
  integer jp1
  integer left
  integer lsofar
  real scrtch(k,k)      
  real sum
  real t(n+k)
!
  lsofar = 0
  break(1) = t(k)
  
  do left = k, n
!
!  Find the next nontrivial knot interval.
!
    if ( t(left+1) == t(left) ) then
      cycle
    end if

    lsofar = lsofar+1
    break(lsofar+1) = t(left+1)

    if ( k <= 1 ) then
      coef(1,lsofar) = bcoef(left)
      cycle
    end if
!
!  Store the K B-spline coefficients relevant to current knot 
!  interval in SCRTCH(*,1).
!
    do i = 1, k
      scrtch(i,1) = bcoef(left-k+i)
    end do
!
!  For j=1,...,k-1, compute the  k-j  b-spline coeff.s relevant to
!  current knot interval for the j-th derivative by differencing
!  those for the (j-1)st derivative, and store in scrtch(.,j+1) .
!
    do jp1 = 2, k
      j = jp1-1
      do i = 1, k-j
        diff = t(left+i)-t(left+i-(k-j))
        if ( diff > 0.0 ) then
          scrtch(i,jp1)=((scrtch(i+1,j)-scrtch(i,j)) /diff)*real(k-j)
        end if
      end do
    end do
!
!  For J=0, ..., K-1, find the values at  t(left)  of the  j+1
!  B-splines of order J+1 whose support contains the current
!  knot interval from those of order J (in  biatx ), then comb-
!  ine with the B-spline coeff.s (in scrtch(.,k-j) ) found earlier
!  to compute the (k-j-1)st derivative at  t(left)  of the given
!  spline.
!
!  Note. if the repeated calls to  bsplvb  are thought to gene-
!  rate too much overhead, then replace the first call by
!    biatx(1)=1.
!  and the subsequent call by the statement
!    j=jp1-1
!  followed by a direct copy of the lines
!    deltar(j)=t(left+j)-x
!    ......
!    biatx(j+1)=saved
!  from  bsplvb . deltal(kmax)  and  deltar(kmax)  would have to
!  appear in a dimension statement, of course.
!
    call bsplvb(t,1,1,t(left),left,biatx)

    coef(k,lsofar) = scrtch(1,k)
    
    do jp1 = 2, k
    
      call bsplvb(t,jp1,2,t(left),left,biatx)

      sum = 0.0
      do i = 1, jp1
        sum = biatx(i)*scrtch(i,k+1-jp1)+sum
      end do

      coef(k+1-jp1,lsofar) = sum
      
    end do

  end do
   
  l = lsofar

  return
end
subroutine bsplvb ( t, jhigh, index, x, left, biatx )
!
!*************************************************************************
!
!! BSPLVB calculates all possibly nonzero B-splines at X of order
!
!    JOUT = MAX( JHIGH, (J+1)*(INDEX-1) ) 
!  
!  with knot sequence T.
! 
!  The recurrence relation
! 
!                   X - T(I)              T(I+J+1) - X
!  B(I,J+1)(X) = -----------B(I,J)(X) + ---------------B(I+1,J)(X)
!                 T(I+J)-T(I)            T(I+J+1)-T(I+1)
! 
!  is used (repeatedly) to generate the (J+1)-vector  
!
!    B(LEFT-J,J+1)(X), ..., B(LEFT,J+1)(X)  
!
!  from the J-vector  
!
!    B(LEFT-J+1,J)(X), ..., B(LEFT,J)(X), 
!
!  storing the new values in BIATX over the old. 
!
!  The facts that 
!
!    B(I,1)=1  if  T(I) <= X  < T(I+1)
!
!  and that 
!
!    B(I,J)(X) = 0  unless  T(I) <= X  < T(I+J)
!
!  are used. 
!
!  The particular organization of the calculations follows 
!  algorithm (8) in chapter X of the text.
!
!  Parameters:
!
!    Input, real T(LEFT+JOUT), the knot sequence.
!    T is assumed to be nondecreasing, and also, T(LEFT) must
!    be strictly less than T(LEFT + 1).
! 
!    Input, integer JHIGH, INDEX, determine the order 
!    JOUT = MAX(JHIGH,(J+1)*(INDEX-1))  
!    of the B-splines whose values at X are to be returned.  
!    INDEX is used to avoid recalculations when several 
!    columns of the triangular array of B-spline values are
!    needed (e.g., in  BVALUE  or in  BSPLVD ).
! 
!    If INDEX= 1, the calculation starts from scratch and the entire 
!    triangular array of B-spline values of orders
!    1,2,...,JHIGH  is generated order by order, i.e., 
!    column by column.
!   
!    If INDEX = 2, only the B-spline values of order J+1, J+2, ..., JOUT  
!    are generated, the assumption being that BIATX, J, 
!    DELTAL, DELTAR are, on entry, as they were on exit 
!    at the previous call.  In particular, if JHIGH = 0, 
!    then JOUT = J+1, i.e., just the next column of B-spline 
!    values is generated.
! 
!    WARNING: the restriction  JOUT <= JMAX (= 20) is
!    imposed arbitrarily by the dimension statement for DELTAL
!    and DELTAR, but is nowhere checked for.
! 
!    Input, real X, the point at which the B-splines are to be evaluated.
! 
!    Input, integer LEFT, an integer chosen (usually) so that 
!    T(LEFT) <= X.le.T(LEFT+1).
! 
!    Output, real BIATX(JOUT), with BIATX(I) containing the
!    value at X of the polynomial of order JOUT which agrees 
!    with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval 
!    (T(LEFT),T(LEFT+1)).
!
  integer, parameter :: jmax = 20
!
  integer jhigh
!
  real biatx(jhigh)
  real deltal(jmax)
  real deltar(jmax)
  integer i
  integer index
  integer j
  integer jp1
  integer left
  real saved
  real t(left+jhigh)
  real term
  real x
!
  save deltal
  save deltar
  save j
!
  data j / 1 /
!
  if ( index == 1 ) then 
    j = 1
    biatx(1) = 1.0
    if ( j >= jhigh ) then
      return
    end if
  end if

10    continue
   
  jp1 = j+1
  deltar(j) = t(left+j)-x
  deltal(j) = x-t(left+1-j)

  saved = 0.0
  do i = 1, j
    term = biatx(i)/(deltar(i)+deltal(jp1-i))
    biatx(i) = saved+deltar(i)*term
    saved = deltal(jp1-i)*term
  end do

  biatx(jp1) = saved
  j = jp1
  if ( j < jhigh ) go to 10

  return
end
subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )
!
!*************************************************************************
!
!! BSPLVD calculates the nonvanishing B-splines and derivatives at X.
! 
!
!  Discussion:
!
!    Values at X of all the relevant B-splines of order K, K-1,..., K+1-NDERIV 
!    are generated via BSPLVB and stored temporarily in DBIATX.  
!
!    Then, the B-spline coefficients of the required derivatives 
!    of the B-splines of interest are generated by differencing, 
!    each from the preceding one of lower order, and combined with 
!    the values of B-splines of corresponding order in DBIATX 
!    to produce the desired values.
!
!  Parameters:
!
!    Input, real T(LEFT+K), the knot sequence.  It is assumed that 
!    T(LEFT) < T(LEFT+1).  Also, the output is correct only if 
!    T(LEFT) <= X <= T(LEFT+1) .
!
!    Input, integer K, the order of the B-splines to be evaluated.
! 
!    Input, real X, the point at which these values are sought.
! 
!    Input, integer LEFT, indicates the left endpoint of the interval of 
!    interest.  The K B-splines whose support contains the interval 
!    (T(LEFT), T(LEFT+1)) are to be considered.
! 
!    Workspace, real A(K,K).
! 
!    Output, real DBIATX(K,NDERIV).  DBIATX(I,M) contains the value of the 
!    (M-1)st derivative of the (LEFT-K+I)-th B-spline of order K for knot
!    sequence T, I=M,...,K, M=1,...,NDERIV.
!
!    Input, integer NDERIV, indicates that values of 
!    B-splines and their derivatives up to but not
!    including the NDERIV-th are asked for. 
!
  integer k
  integer left
  integer nderiv
!
  real a(k,k)
  real dbiatx(k,nderiv)
  real factor
  real fkp1mm
  integer i
  integer ideriv
  integer il
  integer j
  integer jlow
  integer jp1mid
  integer ldummy
  integer m
  integer mhigh
  real sum
  real t(left+k)
  real x
!
  mhigh = max ( min ( nderiv, k ), 1 )
!
!  MHIGH is usually equal to nderiv.
!
  call bsplvb(t,k+1-mhigh,1,x,left,dbiatx)
  
  if ( mhigh == 1 ) then
    return
  end if
!
!  The first column of DBIATX always contains the B-spline values
!  for the current order.  These are stored in column K+1-current
!  order  before BSPLVB is called to put values for the next
!  higher order on top of it.
!
  ideriv = mhigh
  do m = 2, mhigh
    jp1mid = 1
    do j = ideriv, k
      dbiatx(j,ideriv) = dbiatx(jp1mid,1)
      jp1mid = jp1mid+1
    end do
    ideriv = ideriv-1
    call bsplvb(t,k+1-ideriv,2,x,left,dbiatx)
  end do
!
!  At this point,  b(left-k+i, k+1-j)(x) is in  dbiatx(i,j) for
!  i=j,...,k and j=1,...,mhigh ('=' nderiv). in particular, the
!  first column of  dbiatx  is already in final form. to obtain cor-
!  ???  LOST A LINE ???
!  rate their b-repr. by differencing, then evaluate at  x.
!
  jlow = 1
  do i = 1, k
    do j = jlow,k
      a(j,i) = 0.0
    end do
    jlow = i
    a(i,i) = 1.0
  end do
!
!  At this point, a(.,j) contains the b-coeffs for the j-th of the
!  k  b-splines of interest here.
!
  do m = 2, mhigh
  
    fkp1mm = real ( k+1-m)
    il = left
    i = k
!
!  For j=1,...,k, construct b-coeffs of  (m-1)st  derivative of
!  b-splines from those for preceding derivative by differencing
!  and store again in  a(.,j) . the fact that  a(i,j)=0  for
!  i < j  is used.
!
    do ldummy = 1, k+1-m
    
      factor = fkp1mm/(t(il+k+1-m)-t(il))
!
!  The assumption that t(left) < t(left+1) makes denominator
!  in  factor  nonzero.
!
      do j = 1, i
        a(i,j) = (a(i,j)-a(i-1,j))*factor
      end do

      il = il-1
      i = i-1
      
    end do
!
!  For i=1,...,k, combine b-coeffs a(.,i) with b-spline values
!  stored in dbiatx(.,m) to get value of  (m-1)st  derivative of
!  i-th b-spline (of interest here) at  x , and store in
!  dbiatx(i,m). storage of this value over the value of a b-spline
!  of order m there is safe since the remaining b-spline derivat-
!  ives of the same order do not use this value due to the fact
!  that  a(j,i)=0  for j < i .
!
    do i = 1, k
 
      sum = 0.0
      jlow = max(i,m)
      do j = jlow,k
        sum = a(j,i)*dbiatx(j,m)+sum
      end do
 
      dbiatx(i,m) = sum
    end do
    
  end do
 
  return
end
subroutine bspp2d ( t, bcoef, n, k, m, scrtch, break, coef, l )
!
!*************************************************************************
!
!! BSPP2D converts a spline from B-spline to piecewise polynomial representation.
!
!
!  Discussion:
!
!    The B-spline representation
! 
!      T, BCOEF(.,J), N, K 
!
!    is converted to its piecewise polynomial representation 
!
!      BREAK, COEF(J,.,.), L, K, J=1, ..., M.
!
!    This is an extended version of BSPLPP for use with tensor products.
!
!    For each breakpoint interval, the K relevant B-spline 
!    coefficients of the spline are found and then differenced 
!    repeatedly to get the B-spline coefficients of all the 
!    derivatives of the spline on that interval. 
!
!    The spline and its first K-1 derivatives are then evaluated 
!    at the left endpoint of that interval, using BSPLVB 
!    repeatedly to obtain the values of all B-splines of the 
!    appropriate order at that point.
! 
!  Parameters:
!
!    Input, real T(N+K), the knot sequence.
! 
!    Input, real BCOEF(N,M).  For each J, B(*,J) is the B-spline coefficient 
!    sequence, of length N.
! 
!    Input, integer N, the length of  BCOEF.
! 
!    Input, integer K, the order of the spline.
! 
!    Input, INTEGE M, the number of data sets.
! 
!    Work array, real SCRTCH(K,K,M).
! 
!    Output, real BREAK(L+1), the breakpoint sequence, of length L+1, 
!    containing the distinct points in the sequence T(K),...,T(N+1)
! 
!    Output, real COEF(M,K,N), with COEF(MM,I,J) = the (I-1)st derivative of 
!    the MM-th spline at BREAK(J) from the right, MM=1, ..., M.
! 
!    Output, integer L, the number of polynomial pieces which make up the 
!    spline in the interval (T(K), T(N+1)).
!
  integer, parameter :: kmax = 20
!
  integer k
  integer m
  integer n
!
  real bcoef(n,m)
  real biatx(kmax)
  real break(*)
  real coef(m,k,*)
  real diff
  real fkmj
  integer i
  integer j
  integer jp1
  integer kmj
  integer l
  integer left
  integer lsofar
  integer mm
  real scrtch(k,k,m)
  real sum
  real t(n+k)
!
!     dimension coef(k,l)
!
  lsofar = 0
  break(1) = t(k)
  
  do left = k, n
!
!  Find the next nontrivial knot interval.
!
    if ( t(left+1) == t(left) ) then
      cycle
    end if

    lsofar = lsofar+1
    break(lsofar+1) = t(left+1)

    if ( k <= 1 ) then
    
      do mm = 1, m
        coef(mm,1,lsofar) = bcoef(left,mm)
      end do

      cycle

    end if
!
!  Store the K b-spline coeff.s relevant to current knot interval
!  in  scrtch(.,1) .
!
    do i = 1, k
      do mm = 1, m
        scrtch(i,1,mm) = bcoef(left-k+i,mm)
      end do
    end do
!
!  for j=1,...,k-1, compute the  k-j  b-spline coeff.s relevant to
!  current knot interval for the j-th derivative by differencing
!  those for the (j-1)st derivative, and store in scrtch(.,j+1) .
!
    do jp1 = 2, k
    
      j = jp1-1
      kmj = k-j
      fkmj = real ( k-j)
      
      do i = 1, k-j
      
        diff = (t(left+i)-t(left+i-kmj))/fkmj
        
        if ( diff > 0.0 ) then
        
          do mm = 1, m
            scrtch(i,jp1,mm)=(scrtch(i+1,j,mm)-scrtch(i,j,mm))/diff
          end do
          
        end if
        
      end do
        
    end do
!
!  For  j=0, ..., k-1, find the values at  t(left)  of the  j+1
!  b-splines of order  j+1  whose support contains the current
!  knot interval from those of order  j  (in  biatx ), then comb-
!  ine with the b-spline coeff.s (in scrtch(.,k-j) ) found earlier
!  to compute the (k-j-1)st derivative at  t(left)  of the given
!  spline.
!
!  Note. if the repeated calls to  bsplvb  are thought to gene-
!  rate too much overhead, then replace the first call by
!    biatx(1)=1.
!  and the subsequent call by the statement
!    j=jp1-1
!  followed by a direct copy of the lines
!    deltar(j)=t(left+j)-x
!    ...
!    biatx(j+1)=saved
!  from  bsplvb . deltal(kmax)  and  deltar(kmax)  would have to
!  appear in a dimension statement, of course.
!
    call bsplvb(t,1,1,t(left),left,biatx)
    
    do mm = 1, m
      coef(mm,k,lsofar) = scrtch(1,k,mm)
    end do
    
    do jp1 = 2, k
    
      call bsplvb (t,jp1,2,t(left),left,biatx)
      kmj = k+1-jp1
      
      do mm = 1, m
      
        sum = 0.0
        do i = 1, jp1
          sum = biatx(i)*scrtch(i,kmj,mm)+sum
        end do
        
        coef(mm,kmj,lsofar) = sum
        
      end do
      
    end do
   
  end do
   
  l = lsofar
  
  return
end
function bvalue ( t, bcoef, n, k, x, jderiv )
!
!*************************************************************************
!
!! BVALUE evaluates a derivative of a spline from its B-spline representation.  
!
!
!  The spline is taken to be continuous from the right.
! 
!  The nontrivial knot interval (T(I),T(I+1)) containing X is 
!  located with the aid of INTERV.  The K B-spline coefficients 
!  of F relevant for this interval are then obtained from BCOEF, 
!  or are taken to be zero if not explicitly available, and are 
!  then differenced JDERIV times to obtain the B-spline 
!  coefficients of (D**JDERIV)F relevant for that interval.  
!
!  Precisely, with J = JDERIV, we have from X.(12) of the text that:
! 
!    (D**J)F = sum ( BCOEF(.,J)*B(.,K-J,T) )
! 
!  where
!                   / BCOEF(.),                    ,  J == 0
!                   /
!    BCOEF(.,J) = / BCOEF(.,J-1) - BCOEF(.-1,J-1)
!                   / -----------------------------,  J > 0
!                   /    (T(.+K-J) - T(.))/(K-J)
! 
!  Then, we use repeatedly the fact that
! 
!    sum ( A(.)*B(.,M,T)(X) ) = sum ( A(.,X)*B(.,M-1,T)(X) )
! 
!  with
!                 (X - T(.))*A(.) + (T(.+M-1) - X)*A(.-1)
!    A(.,X) =   ---------------------------------------
!                 (X - T(.))      + (T(.+M-1) - X)
! 
!  to write (D**J)F(X) eventually as a linear combination of 
!  B-splines of order 1, and the coefficient for B(I,1,T)(X) 
!  must then be the desired number (D**J)F(X).
!  See x.(17)-(19) of text.
!
!  Parameters:
!
!    Input, real T(N+K), the knot sequence.  T is assumed to be nondecreasing.
!
!    Input, real BCOEF(N), B-spline coefficient sequence.
! 
!    Input, integer N, the length of BCOEF.
! 
!    Input, integer K, the order of the spline.
!    WARNING: the restriction K <= KMAX (=20)  is imposed
!    arbitrarily by the dimension statement for AJ, DL, DR
!    but is nowhere checked for.
! 
!    Input, real X, the point at which to evaluate.
! 
!    Input, integer JDERIV, the order of the derivative to 
!    be evaluated.  JDERIV is assumed to be zero or positive.
! 
!    Output, real BVALUE, the value of the (JDERIV)-th 
!    derivative of the spline at X.
! 
  integer, parameter :: kmax = 20
!
  integer k
  integer n
!
  real aj(kmax)
  real bcoef(n)
  real bvalue
  real dl(kmax)
  real dr(kmax)
  integer i
  integer ilo
  integer j
  integer jc
  integer jcmax
  integer jcmin
  integer jderiv
  integer jj
  integer mflag
  real t(n+k)
  real x
!
  bvalue = 0.0
  
  if ( jderiv >= k ) then
    return
  end if
!
!  Find  i   s.t.   1 <= i  < n+k   and   t(i) < t(i+1)   and
!  t(i) <= x  < t(i+1) . if no such i can be found,  x  lies
!  outside the support of  the spline  f  and  bvalue=0.
!  (the asymmetry in this choice of  i  makes  f  rightcontinuous)
!
  call interv(t,n+k,x,i,mflag)
  
  if ( mflag /= 0 ) then
    return
  end if
!
!  If K=1 (and jderiv = 0), bvalue = bcoef(i).
!
  if ( k <= 1 ) then
    bvalue = bcoef(i)
    return
  end if
!
!  Store the k b-spline coefficients relevant for the knot interval
!  (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j)=x-t(i+1-j),
!  dr(j)=t(i+j)-x, j=1,...,k-1 . set any of the aj not obtainable
!  from input to zero. set any t.s not obtainable equal to t(1) or
!  to t(n+k) appropriately.
!
  jcmin = 1
  
  if ( i >= k ) then
  
    do j = 1, k-1
      dl(j) = x-t(i+1-j)
    end do
    
  else
  
    jcmin = 1-(i-k)
 
    do j = 1, i
      dl(j) = x-t(i+1-j)
    end do
 
    do j = i, k-1
      aj(k-j) = 0.0
      dl(j) = dl(i)
    end do
  
  end if
 
  jcmax = k

  if ( n >= i ) then
    go to 90
  end if
 
  jcmax = k+n-i
  do j = 1, k+n-i
    dr(j) = t(i+j)-x
  end do
 
  do j = k+n-i, k-1
    aj(j+1) = 0.0
    dr(j) = dr(k+n-i)
  end do
 
  go to 110
 
   90 continue
 
  do j = 1, k-1
    dr(j) = t(i+j)-x
  end do
 
  110 continue
 
  do jc = jcmin, jcmax
    aj(jc) = bcoef(i-k+jc)
  end do
!
!  Difference the coefficients JDERIV times.
!
  do j = 1, jderiv
 
    ilo = k-j
    do jj = 1, k-j
      aj(jj) = ((aj(jj+1)-aj(jj))/(dl(ilo)+dr(jj)))*real(k-j)
      ilo = ilo-1
    end do
 
  end do
!
!  Compute value at X in (t(i),t(i+1)) of jderiv-th derivative,
!  given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
!
  do j = jderiv+1, k-1
    ilo = k-j
    do jj = 1, k-j
      aj(jj) = (aj(jj+1)*dl(ilo)+aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
      ilo = ilo-1
    end do
  end do
  
  bvalue = aj(1)
 
  return
end
subroutine chol1d ( p, v, qty, npoint, u, qu )
!
!*************************************************************************
!
!! CHOL1D constructs the upper three diagonals of
!
!     V(I,J), I=2 to NPOINT-1, J=1,3, 
!
!  of the matrix 
!
!    6*(1-P)*Q-transpose*(D**2)*Q + P*R.
! 
!  It then computes its L*L-transpose decomposition and stores it also
!  in V, then applies forward and backsubstitution to the right side 
!
!    Q-transpose*Y 
!
!  in QTY to obtain the solution in U.
! 
  integer npoint
!
  integer i
  real p
  real qty(npoint)
  real qu(npoint)
  real u(npoint)
  real v(npoint,7)
  real prev
  real ratio
  real six1mp
  real twop
!
!  Construct 6*(1-p)*q-transp.*(d**2)*q + p*r
!
  six1mp = 6.0*(1.0-p)
  twop = 2.0*p
  
  do i = 2, npoint-1
    v(i,1) = six1mp*v(i,5)+twop*(v(i-1,4)+v(i,4))
    v(i,2) = six1mp*v(i,6)+p*v(i,4)
    v(i,3) = six1mp*v(i,7)
  end do
 
  if ( npoint < 4 ) then
    u(1) = 0.0
    u(2) = qty(2)/v(2,1)
    u(3) = 0.0
!
!  Factorization
!
  else

    do i = 2, npoint-2
      ratio = v(i,2)/v(i,1)
      v(i+1,1) = v(i+1,1)-ratio*v(i,2)
      v(i+1,2) = v(i+1,2)-ratio*v(i,3)
      v(i,2) = ratio
      ratio = v(i,3)/v(i,1)
      v(i+2,1) = v(i+2,1)-ratio*v(i,3)
      v(i,3) = ratio
    end do
!
!  Forward substitution
!
    u(1) = 0.0
    v(1,3) = 0.0
    u(2) = qty(2)
    do i = 2, npoint-2
      u(i+1) = qty(i+1)-v(i,2)*u(i)-v(i-1,3)*u(i-1)
    end do
!
!  Back substitution.
!
    u(npoint) = 0.0
    u(npoint-1) = u(npoint-1) / v(npoint-1,1)

    do i = npoint-2, 2, -1
      u(i) = u(i)/v(i,1)-u(i+1)*v(i,2)-u(i+2)*v(i,3)
    end do

  end if
!
!  Construct Q*U.
!
  prev = 0.0
  do i = 2, npoint
    qu(i) = (u(i)-u(i-1))/v(i-1,4)
    qu(i-1) = qu(i)-prev
    prev = qu(i)
  end do

  qu(npoint) = - qu(npoint)
    
  return
end
subroutine colloc ( aleft, aright, lbegin, iorder, ntimes, addbrk, &
  relerr )
!
!*************************************************************************
!
!! COLLOC solves an ordinary differential equation by collocation.
!
!
!  Method:
!
!    The M-th order ordinary differential equation with M side 
!    conditions, to be specified in subroutine DIFEQU, is solved 
!    approximately by collocation.
!
!    The approximation F to the solution G is piecewise polynomial of order 
!    k+m with l pieces and  m-1 continuous derivatives.   F is determined by 
!    the requirement that it satisfy the differential equation at K points 
!    per interval (to be specified in COLPNT ) and the M side conditions.
!
!    This usually nonlinear system of equations for  f  is solved by
!    newton's method. the resulting linear system for the b-coeffs of an
!    iterate is constructed appropriately in  e q b l o k  and then solved
!    in  s l v b l k , a program designed to solve  a l m o s t  b l o c k
!    d i a g o n a l  linear systems efficiently.
!
!    There is an opportunity to attempt improvement of the breakpoint
!    sequence (both in number and location) through use of  n e w n o t .
!
!    Printed output consists of the pp-representation of the approximate 
!    solution, and of the error at selected points.
!
!  Parameters: 
!
!  aleft, 
!  aright endpoints of interval of approximation
!
!  lbegin   initial number of polynomial pieces in the approximation.
!           a uniform breakpoint sequence is chosen.
!
!  iorder   order of polynomial pieces in the approximation
!
!  ntimes   number of passes through  n e w n o t  to be made
!
!  addbrk   the number (possibly fractional) of breaks to be added per
!           pass through newnot. e.g., if addbrk=.33334, then a break-
!           point will be added at every third pass through newnot.
!
!  relerr   a tolerance. newton iteration is stopped if the difference
!           between the b-coeffs of two successive iterates is no more
!           than  relerr*(absol.largest b-coefficient).
!
  integer, parameter :: npiece = 100
  integer, parameter :: ndim = 200
  integer, parameter :: ncoef = 2000
  integer, parameter :: lenblk = 2000
!
  real a(ndim)
  real addbrk
  real aleft
  real amax
  real aright
  real asave(ndim)
  real b(ndim)
  real bloks(lenblk)
  real break
  real coef
  real dx
  real err
  integer i
  integer iflag
  integer ii
  integer integs(3,npiece)
  integer iorder
  integer iside
  integer itemps(ndim)
  integer iter
  integer itermx
  integer j
  integer k
  integer kpm
  integer l
  integer lbegin
  integer lnew
  integer m
  integer n
  integer nbloks
  integer nt
  integer ntimes
  real relerr
  real rho
  real t(ndim)
  real templ(lenblk)
  real temps(ndim)
  real xside
!
  equivalence (bloks,templ)
!
  common /approx/ break(npiece),coef(ncoef),l,kpm
  common /side/ m,iside,xside(10)
  common /other/ itermx,k,rho(19)
!
  kpm = iorder

  if ( lbegin*kpm > ncoef ) then
    go to 120
  end if
!
!  Set the various parameters concerning the particular dif.equ.
!  including a first approx. in case the de is to be solved by
!  iteration ( itermx > 0).
!
  call difequ ( 1, temps(1), temps )
!
!  Obtain the  k  collocation points for the standard interval.
!
  k = kpm-m
  call colpnt(k,rho)
!
!  The following five statements could be replaced by a read in or-
!  der to obtain a specific (nonuniform) spacing of the breakpnts.
!
  dx = (aright-aleft)/real ( lbegin)
 
  temps(1) = aleft
  do i = 2, lbegin
    temps(i) = temps(i-1)+dx
  end do
  temps(lbegin+1) = aright
!
!  Generate the required knots t(1),...,t(n+kpm).
!
  call knots ( temps, lbegin, kpm, t, n )
  nt = 1
!
!  Generate the almost block diagonal coefficient matrix  bloks  and
!  right side  b  from collocation equations and side conditions.
!  then solve via  slvblk , obtaining the b-representation of the ap-
!  proximation in  t , a ,  n  , kpm  .
!
20    continue

  call eqblok ( t, n, kpm, temps, a, bloks, lenblk, integs, nbloks, b )

  call slvblk(bloks,integs,nbloks,b,itemps,a,iflag)

  iter = 1
  if ( itermx <= 1 ) then
    go to 60
  end if
!
!  Save b-spline coeff. of current approx. in  asave , then get new
!  approx. and compare with old. if coeff. are more than  relerr
!  apart (relatively) or if no. of iterations is less than  itermx ,
!  continue iterating.
!
   30 continue

  call bsplpp(t,a,n,kpm,templ,break,coef,l)
 
  do i = 1, n
    asave(i) = a(i)
  end do
 
  call eqblok ( t, n, kpm, temps, a, bloks, lenblk, integs, nbloks, b )

  call slvblk(bloks,integs,nbloks,b,itemps,a,iflag)
 
  err = 0.0
  amax = 0.0
  do i = 1, n
    amax = max ( amax, abs ( a(i) ) )
    err = max ( err, abs ( a(i)-asave(i) ) )
  end do
 
  if ( err <= relerr*amax ) then
    go to 60
  end if

  iter = iter+1

  if ( iter < itermx ) then
    go to 30
  end if
!
!  Iteration (if any) completed. print out approx. based on current
!  breakpoint sequence, then try to improve the sequence.
!
   60 continue

  write(*,70)kpm,l,n,(break(i),i=2,l)
   70 format (' approximation from a space of splines of order',i3, &
     ' on ',i3,' intervals,'/' of dimension',i4,'.  breakpoints -'/ &
     (5e20.10))

  if ( itermx > 0 ) then
    write(*,*)' '
    write(*,*)'Results on interation ',iter
  end if

  call bsplpp(t,a,n,kpm,templ,break,coef,l)
  
  write ( *, * ) ' '
  write ( *, * ) 'The piecewise polynomial representation of the approximation:'
  write ( *, * ) ' '
   
  do i = 1, l
    ii = (i-1)*kpm
    write(*,'(f9.3,e13.6,10e11.3)')break(i),(coef(ii+j),j=1,kpm)
  end do
!
!  The following call is provided here for possible further analysis
!  of the approximation specific to the problem being solved.
!  it is, of course, easily omitted.
!
  call difequ(4,temps(1),temps)
 
  if ( nt > ntimes ) then
    return
  end if
!
!  From the pp-rep. of the current approx., obtain in  newnot  a new
!  (and possibly better) sequence of breakpoints, adding (on the 
!  average) ADDBRK breakpoints per pass through NEWNOT.
!
  lnew = lbegin+int(real ( nt)*addbrk)

  if ( lnew*kpm > ncoef ) then
    go to 120
  end if

  call newnot(break,coef,l,kpm,temps,lnew,templ)

  call knots ( temps, lnew, kpm, t, n )
  nt = nt+1
  go to 20
  
  120 continue
  write(*,*)' '
  write(*,*)'COLLOC - Fatal error!'
  write(*,*)'  The assigned dimension for COEF is ',ncoef
  write(*,*)'  but this is too small.'
  stop
end
subroutine colpnt ( k, rho )
!
!*************************************************************************
!
!! COLPNT supplies collocation points.
!
!
!  Discussion:
!
!    The collocation points are for the standard interval (-1,1) as the 
!    zeros of the legendre polynomial of degree K, provided K <= 8.  
!
!    Otherwise, uniformly spaced points are given.
!
!  Parameters:
!
!    Input, integer K, the number of collocation points desired.
!
!    Output, real RHO(K), the collocation points.
!
  integer k
!
  integer j
  real rho(k)
!
  if ( k == 1 ) then
    rho(1) = 0.0
  else if ( k == 2 ) then
    rho(1) = -0.577350269189626
    rho(2) = 0.577350269189626
  else if ( k == 3 ) then
    rho(1) = -0.774596669241483
    rho(2) = 0.0
    rho(3) = 0.774596669241483
  else if ( k == 4 ) then
    rho(1) = -0.861136311594053
    rho(2) = -0.339981043584856
    rho(3) = 0.339981043584856
    rho(4) = 0.861136311594053
  else if ( k == 5 ) then
    rho(1) = -0.906179845938664
    rho(2) = -0.538469310105683
    rho(3) = 0.0
    rho(4) = 0.538469310105683
    rho(5) = 0.906179845938664
  else if ( k == 6 ) then
    rho(4) = 0.238619186083197
    rho(3) = -rho(4)
    rho(5) = 0.661209386466265
    rho(2) = -rho(5)
    rho(6) = 0.932469514203152
    rho(1) = -rho(6)
  else if ( k == 7 ) then
    rho(5) = 0.405845151377397
    rho(3) = -rho(5)
    rho(6) = 0.741531185599394
    rho(2) = -rho(6)
    rho(7) = 0.949107912342759
    rho(1) = -rho(7)
    rho(4) = 0.0
  else if ( k == 8 ) then
    rho(5) = 0.183434642495650
    rho(4) = -rho(5)
    rho(6) = 0.525532409916329
    rho(3) = -rho(6)
    rho(7) = 0.796666477413627
    rho(2) = -rho(7)
    rho(8) = 0.960289856497536
    rho(1) = -rho(8)
  else

    write(*,*)' '
    write(*,*)'ColPnt - Warning!'
    write(*,*)'  Equispaced collocation points will be used,'
    write(*,*)'  because K =',k,' which is greater than 8.'
 
    do j = 1, k
      rho(j) = 2.0*real ( j-1)/real ( k-1)-1.0
    end do
    
  end if
 
  return
end
subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result, ind )
!
!*************************************************************************
!
!! CSPINT estimates an integral using a spline interpolant.
!
!
!  Discussion:
!
!    CSPINT estimates the integral from A to B of F(X) by
!    computing the natural spline S(X) that interpolates to F 
!    and integrating that exactly.  
!
!    F is supplied to the routine in the form of tabulated data.
! 
!    Other output from the program includes the definite integral
!    from X(1) to X(I) of the spline, and the coefficients
!    necessary for the user to evaluate the spline outside of
!    this routine.
! 
!  Parameters:
!
!  FTAB   Input, real FTAB(NTAB), contains the tabulated values 
!         of the functions, FTAB(I)=F(XTAB(I)).
! 
!  XTAB   Input, real XTAB(NTAB), contains the points at 
!         which the function was evaluated.  The XTAB's must be 
!         distinct and in ascending order.
! 
!  NTAB   Input, integer NTAB, the number of entries in FTAB 
!         and XTAB.  NTAB must be at least 3.
! 
!  A      Input, real A, lower limit of integration.
!  
!  B      Input, real B, upper limit of integration.
! 
!  Y      Output, real Y(3,NTAB), will contain the coefficients
!         of the interpolating natural spline over each 
!         subinterval.
!
!         For XTAB(I) < = X  <= XTAB(I+1),
!
!           S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I)) + Y(2,I)*(X-XTAB(I))**2
!                      + Y(3,I)*(X-XTAB(I))**3
! 
!  E      Output, real E(NTAB), E(I)=the definite integral 
!         from XTAB(1) to XTAB(I) of S(X).
! 
!  WORK   Workspace, real WORK(NTAB).
! 
!  RESULT Output, real RESULT, the estimated value of the integral.
! 
!  IND    Output, integer IND, error flag.
!         IND=0 if NTAB < 3 or the XTAB's are not distinct and in
!         ascending order.
!         IND=1 otherwise.
!
  integer ntab
!
  real a
  real b
  real e(ntab)
  real ftab(ntab)
  integer i
  integer ind
  integer j
  real r
  real result
  real s
  real term
  real u
  real work(ntab)
  real xtab(ntab)
  real y(3,ntab)
!
  ind = 0
 
  if ( ntab < 3 ) then
    write(*,*)' '
    write(*,*)'CSPINT - Fatal error!'
    write(*,*)'  NTAB must be at least 3,'
    write(*,*)'  but your value was NTAB = ',ntab
    stop
  end if
 
  do i = 1, ntab-1
 
    if ( xtab(i+1) <= xtab(i) ) then
      write(*,*)' '
      write(*,*)'CSPINT - Fatal error!'
      write(*,*)'  Interval ',i,' is illegal.'
      write(*,*)'  XTAB(I) =',xtab(i)
      write(*,*)'  XTAB(I+1)=',xtab(i+1)
      stop
    end if
 
  end do
 
  s = 0.0
  do i = 1, ntab-1
    r = (ftab(i+1)-ftab(i))/(xtab(i+1)-xtab(i))
    y(2,i) = r - s
    s = r
  end do
 
  result = 0.0
  s = 0.0
  r = 0.0
  y(2,1) = 0.0
  y(2,ntab) = 0.0
 
  do i = 2, ntab-1
    y(2,i) = y(2,i)+r*y(2,i-1)
    work(i) = 2.0*(xtab(i-1)-xtab(i+1))-r*s
    s = xtab(i+1)-xtab(i)
    r = s / work(i)
  end do
 
  do j = 2, ntab-1
    i = ntab+1-j
    y(2,i) = ((xtab(i+1)-xtab(i))*y(2,i+1)-y(2,i))/work(i)
  end do
 
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    r = y(2,i+1)-y(2,i)
    y(3,i) = r/s
    y(2,i) = 3.0*y(2,i)
    y(1,i) = (ftab(i+1)-ftab(i)) / s - (y(2,i)+r)*s
  end do
 
  e(1) = 0.0
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    term = (((y(3,i)*.25*s+y(2,i)/3.0)*s+y(1,i)*.5)*s+ftab(i))*s
    e(i+1) = e(i)+term
  end do
!
!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!
  r = a
  u = 1.0
  
  do j = 1, 2
  
    if ( r <= xtab(1) ) then
      result = result-u*((r-xtab(1))*y(1,1)*.5+ftab(1))*(r-xtab(1))
    else if ( r >= xtab(ntab) ) then
      result = result-u*(e(ntab)+(r-xtab(ntab))*(ftab(ntab)+.5* &
        (ftab(ntab-1)+(xtab(ntab)-xtab(ntab-1))*y(1,ntab-1))*(r- &
        xtab(ntab))))
    else
      do i = 1, ntab-1
    
        if ( r <= xtab(i+1) ) then
          r = r - xtab(i)
          result = result-u*(e(i)+(((y(3,i)*.25*r+y(2,i)/3.)*r &
            +y(1,i)*.5)*r+ftab(i))*r)
          go to 100
        end if
      
      end do
    
    end if
   
  100   continue

    u = -1.0
    r = b
    
  end do
  
  ind = 1

  return
end
subroutine cubset ( tau, c, n, ibcbeg, ibcend )
!
!*******************************************************************************
!
!! CUBSET sets up a simple cubic spline interpolant.
!
!
!  WARNING: IBCBEG and IBCEND are not set up yet.
!
!  A tridiagonal linear system for the unknown slopes S(I) of
!  F at TAU(I), I=1,..., N, is generated and then solved by Gauss
!  elimination, with S(I) ending up in C(2,I), for all I.
!
!  Parameters:
!
!  TAU    Input, real TAU(N), the abscissas or X values of
!         the data points.  The entries of TAU are assumed to be
!         strictly increasing.
!
!  N      Input, integer N, the number of data points.  N is
!         assumed to be at least 2.
!
!  C      Input/output, real C(4,N).
!
!         On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
!         or C(2,N) should have been set to the desired derivative
!         values, as described further under IBCBEG and IBCEND.
!
!         On output, C contains the polynomial coefficients of
!         the cubic interpolating spline with interior knots
!         TAU(2) through TAU(N-1).
!
!         In the interval interval (TAU(I), TAU(I+1)), the spline
!         F is given by
!
!           F(X) = 
!             C(1,I) + 
!             C(2,I) * ( X - TAU(I) ) +
!             C(3,I) * ( X - TAU(I) )**2 +
!             C(4,I) * ( X - TAU(I) )**3
!
!  IBCBEG,
!  IBCEND Input, integer IBCBEG, IBCEND, boundary condition
!         indicators.
!
!         IBCBEG=0 means no boundary condition at TAU(1) is given.
!         In this case, the "not-a-knot condition" is used.  That
!         is, the jump in the third derivative across TAU(2) is
!         forced to zero.  Thus the first and the second cubic
!         polynomial pieces are made to coincide.
!
!         IBCBEG=1 means that the slope at TAU(1) is to equal the
!         input value C(2,1).
!
!         IBCBEG=2 means that the second derivative at TAU(1) is
!         to equal C(2,1).
!
!         IBCEND=0, 1, or 2 has analogous meaning concerning the
!         boundary condition at TAU(N), with the additional
!         information taken from C(2,N).
!
  integer n
!
  real c(4,n)
  integer ibcbeg
  integer ibcend
  real tau(n)
!
!  Solve for the slopes at internal nodes.
!
  call cubslo ( tau, c, n ) 
!
!  Now compute the quadratic and cubic coefficients used in the 
!  piecewise polynomial representation.
!
  call spline_hermite_set ( n, tau, c )

  return
end
subroutine cubslo ( tau, c, n )
!
!*******************************************************************************
!
!! CUBSLO solves for slopes defining a cubic spline.
!
!
!  A tridiagonal linear system for the unknown slopes S(I) of
!  F at TAU(I), I=1,..., N, is generated and then solved by Gauss
!  elimination, with S(I) ending up in C(2,I), for all I.
!
!  Parameters:
!
!  TAU    Input, real TAU(N), the abscissas or X values of
!         the data points.  The entries of TAU are assumed to be
!         strictly increasing.
!
!  N      Input, integer N, the number of data points.  N is
!         assumed to be at least 2.
!
!  C      Input/output, real C(4,N).
!
!         On input, C(1,I) contains the function value at TAU(I),
!         for I = 1 to N.  
!         C(2,1) contains the slope at TAU(1) and C(2,N) contains
!         the slope at TAU(N).
!
!         On output, the intermediate slopes at TAU(I) have been
!         stored in C(2,I), for I = 2 to N-1.
!
  integer n
!
  real c(4,n)
  real g
  integer i
  integer ibcbeg
  integer ibcend
  real tau(n)
!
!  Set up the right hand side of the linear system.
!  C(2,1) and C(2,N) are presumably already set.
!
  do i = 2, n-1
    c(2,i) = 3.0 * ( &
      ( tau(i) - tau(i-1) ) * ( c(1,i+1) - c(1,i) ) / ( tau(i+1) - tau(i) ) + &
      ( tau(i+1) - tau(i) ) * ( c(1,i) - c(1,i-1) ) / ( tau(i) - tau(i-1) ) )
  end do
!
!  Set the diagonal coefficients.
!
  c(4,1) = 1.0
  do i = 2, n-1
    c(4,i) = 2.0 * ( tau(i+1) - tau(i-1) )
  end do
  c(4,n) = 1.0
!
!  Set the off-diagonal coefficients.
!
  c(3,1) = 0.0
  do i = 2, n
    c(3,i) = tau(i) - tau(i-1)
  end do
!
!  Forward elimination.
!
  do i = 2, n-1
    g = -c(3,i+1) / c(4,i-1)
    c(4,i) = c(4,i) + g * c(3,i-1)
    c(2,i) = c(2,i) + g * c(2,i-1)
  end do
!
!  Back substitution for the interior slopes.
!
  do i = n-1, 2, -1
    c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
  end do

  return
end
subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
!
!*******************************************************************************
!
!! CUBSPL defines an interpolatory cubic spline.
!
!
!  Discussion:
!
!    A tridiagonal linear system for the unknown slopes S(I) of
!    F at TAU(I), I=1,..., N, is generated and then solved by Gauss
!    elimination, with S(I) ending up in C(2,I), for all I.
!
!  Parameters:
!
!  TAU    Input, real TAU(N), the abscissas or X values of
!         the data points.  The entries of TAU are assumed to be
!         strictly increasing.
!
!  N      Input, integer N, the number of data points.  N is
!         assumed to be at least 2.
!
!  C      Input/output, real C(4,N).
!
!         On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
!         or C(2,N) should have been set to the desired derivative
!         values, as described further under IBCBEG and IBCEND.
!
!         On output, C contains the polynomial coefficients of
!         the cubic interpolating spline with interior knots
!         TAU(2) through TAU(N-1).
!
!         In the interval interval (TAU(I), TAU(I+1)), the spline
!         F is given by
!
!           F(X) = 
!             C(1,I) + 
!             C(2,I) * H +
!             C(3,I) * H**2 / 2 + 
!             C(4,I) * H**3 / 6.
!
!         where H=X-TAU(I).  The routine PPVALU may be used to
!         evaluate F or its derivatives from TAU, C, L=N-1,
!         and K=4.
!
!  IBCBEG,
!  IBCEND Input, integer IBCBEG, IBCEND, boundary condition
!         indicators.
!
!         IBCBEG=0 means no boundary condition at TAU(1) is given.
!         In this case, the "not-a-knot condition" is used.  That
!         is, the jump in the third derivative across TAU(2) is
!         forced to zero.  Thus the first and the second cubic
!         polynomial pieces are made to coincide.
!
!         IBCBEG=1 means that the slope at TAU(1) is to equal the
!         input value C(2,1).
!
!         IBCBEG=2 means that the second derivative at TAU(1) is
!         to equal C(2,1).
!
!         IBCEND=0, 1, or 2 has analogous meaning concerning the
!         boundary condition at TAU(N), with the additional
!         information taken from C(2,N).
!
  integer n
!
  real c(4,n)
  real divdf1
  real divdf3
  real dtau
  real g
  integer i
  integer ibcbeg
  integer ibcend
  real tau(n)
!
!  C(3,*) and C(4,*) are used initially for temporary storage.
!
!  Store first differences of the TAU sequence in C(3,*).
!
!  Store first divided difference of data in C(4,*).
!
  do i = 2, n
    c(3,i) = tau(i) - tau(i-1)
  end do

  do i = 2, n 
    c(4,i) = ( c(1,i) - c(1,i-1) ) / ( tau(i) - tau(i-1) )
  end do
!
!  Construct the first equation from the boundary condition
!  at the left endpoint, of the form:
!
!    C(4,1)*S(1) + C(3,1)*S(2) = C(2,1)
!
!  IBCBEG = 0: Not-a-knot
!
  if ( ibcbeg == 0 ) then

    if ( n <= 2 ) then
      c(4,1) = 1.0
      c(3,1) = 1.0
      c(2,1) = 2.0*c(4,2)
      go to 120
    end if

    c(4,1) = c(3,3)
    c(3,1) = c(3,2)+c(3,3)
    c(2,1) = ((c(3,2)+2.0d0*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3)) / c(3,1)
!
!  IBCBEG = 1: derivative specified.
!
  else if ( ibcbeg == 1 ) then

    c(4,1) = 1.0
    c(3,1) = 0.0

    if ( n == 2 ) then
      go to 120
    end if
!
!  Second derivative prescribed at left end.
!
  else

    c(4,1) = 2.0
    c(3,1) = 1.0
    c(2,1) = 3.0*c(4,2) - c(3,2)/2.0*c(2,1)

    if ( n == 2 ) then
      go to 120
    end if

  end if
!
!  If there are interior knots, generate the corresponding
!  equations and carry out the forward pass of Gauss elimination,
!  after which the I-th equation reads:
!
!    C(4,I) * S(I) + C(3,I) * S(I+1) = C(2,I).
!
  do i = 2, n-1
    g = -c(3,i+1) / c(4,i-1)
    c(2,i) = g*c(2,i-1)+3.0*(c(3,i)*c(4,i+1)+c(3,i+1)*c(4,i))
    c(4,i) = g*c(3,i-1)+2.0*(c(3,i)+c(3,i+1))
  end do
!
!  Construct the last equation from the second boundary condition, of
!  the form
!
!    -G * C(4,N-1) * S(N-1) + C(4,N) * S(N) = C(2,N)
!
!  If slope is prescribed at right end, one can go directly to
!  back-substitution, since the C array happens to be set up just
!  right for it at this point.
!
  if ( ibcend == 1 ) then
    go to 160
  end if

  if ( ibcend > 1 ) then
    go to 110
  end if
 
90    continue
!
!  Not-a-knot and N >= 3, and either N > 3 or also not-a-knot
!  at left end point.
!
  if ( n /= 3 .or. ibcbeg /= 0 ) then
    g = c(3,n-1)+c(3,n)
    c(2,n) = ((c(3,n)+2.0*g)*c(4,n)*c(3,n-1)+c(3,n)**2 &
      *(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
    g = -g/c(4,n-1)
    c(4,n) = c(3,n-1)
    c(4,n) = c(4,n) + g*c(3,n-1)
    c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
    go to 160
  end if
!
!  N=3 and not-a-knot also at left.
!
100   continue
 
  c(2,n) = 2.0*c(4,n)
  c(4,n) = 1.0
  g = -1.0/c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
  go to 160
!
!  IBCEND = 2: Second derivative prescribed at right endpoint.
!
110   continue
 
  c(2,n) = 3.0*c(4,n)+c(3,n)/2.0*c(2,n)
  c(4,n) = 2.0
  g = -1.0/c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
  go to 160
!
!  N = 2.
!
120   continue
  
  if ( ibcend == 2  ) then

    c(2,n) = 3.0*c(4,n)+c(3,n)/2.0*c(2,n)
    c(4,n) = 2.0
    g = -1.0/c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
 
  else if ( ibcend == 0 .and. ibcbeg /= 0 ) then

    c(2,n) = 2.0*c(4,n)
    c(4,n) = 1.0
    g = -1.0/c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)

  else if ( ibcend == 0 .and. ibcbeg == 0 ) then

    c(2,n) = c(4,n)

  end if
!
!  Back solve the upper triangular system 
!    C(4,I) * S(I) + C(3,I) * S(I+1) = B(I)
!  for the slopes C(2,I), given that S(N) is already known.
!
160   continue
 
  do i = n-1, 1, -1
    c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
  end do
!
!  Generate cubic coefficients in each interval, that is, the
!  derivatives at its left endpoint, from value and slope at its
!endpoints.
!
  do i = 2, n
    dtau = c(3,i)
    divdf1 = ( c(1,i) - c(1,i-1) ) / dtau
    divdf3 = c(2,i-1) + c(2,i) - 2.0 * divdf1
    c(3,i-1) = 2.0 * ( divdf1 - c(2,i-1) - divdf3 ) / dtau
    c(4,i-1) = 6.0 * divdf3 / dtau**2
  end do
 
  return
end
subroutine cwidth ( w, b, nequ, ncols, integs, nbloks, d, x, iflag )
!
!*************************************************************************
!
!! CWIDTH is a variation of the theme in the algorithm bandet1
!  by martin and wilkinson (numer.math. 9(1976)279-307). it solves
!  the linear system
!                           a*x = b
!  of  nequ  equations in case  a  is almost block diagonal with all
!  blocks having  ncols  columns using no more storage than it takes to
!  store the interesting part of  a . such systems occur in the determ-
!  ination of the b-spline coefficients of a spline approximation.
!
!  Parameters:
!
!  w     on input, a two-dimensional array of size (nequ,ncols) contain-
!        ing the interesting part of the almost block diagonal coeffici-
!        ent matrix  a (see description and example below). the array
!        integs  describes the storage scheme.
!        on output, w  contains the upper triangular factor  u  of the
!        lu factorization of a possibly permuted version of  a . in par-
!        ticular, the determinant of  a  could now be found as
!            iflag*w(1,1)*w(2,1)* ... * w(nequ,1)  .
!
!  b     on input, the right side of the linear system, of length  nequ.
!        the contents of  b  are changed during execution.
!
!  nequ  number of equations in system
!
!  ncols  block width, i.e., number of columns in each block.
!
!  integs integer array, of size (2,nequ), describing the block 
!         structure of  a .
!         integs(1,i)=no. of rows in block i              = nrow
!         integs(2,i)=no. of elimination steps in block i
!                     =overhang over next block             = last
!  nbloks number of blocks
!
!  d      work array, to contain row sizes . if storage is scarce, the
!         array  x  could be used in the calling sequence for  d .
!
!  x      on output, contains computed solution (if iflag /= 0), of
!         length  nequ .
!
!  iflag  on output, integer
!        =(-1)**(no.of interchanges during elimination)
!                if  a  is invertible
!        = 0   if  a  is singular
!
!  block structure of  a  
!
!  the interesting part of  a  is taken to consist of  nbloks  con-
!  secutive blocks, with the i-th block made up of  nrowi=integs(1,i)
!  consecutive rows and  ncols  consecutive columns of  a , and with
!  the first  lasti=integs(2,i) columns to the left of the next block.
!  these blocks are stored consecutively in the workarray  w .
!
!  for example, here is an 11th order matrix and its arrangement in
!  the workarray  w . (the interesting entries of  a  are indicated by
!  their row and column index modulo 10.)
!
!                  ---   a   ---                          ---   w   ---
!
!                     nrow1=3
!          11 12 13 14                                     11 12 13 14
!          21 22 23 24                                     21 22 23 24
!          31 32 33 34      nrow2=2                        31 32 33 34
!   last1=2      43 44 45 46                               43 44 45 46
!                53 54 55 56         nrow3=3               53 54 55 56
!         last2=3         66 67 68 69                      66 67 68 69
!                         76 77 78 79                      76 77 78 79
!                         86 87 88 89   nrow4=1            86 87 88 89
!                  last3=1   97 98 99 90   nrow5=2         97 98 99 90
!                     last4=1   08 09 00 01                08 09 00 01
!                               18 19 10 11                18 19 10 11
!                        last5=4
!
!  for this interpretation of  a  as an almost block diagonal matrix,
!  we have  nbloks=5 , and the integs array is
!
!                        i= 1   2   3   4   5
!                  k=
!  integs(k,i)=      1      3   2   3   1   2
!                     2      2   3   1   1   4
!
!
!  Method:
!
!  gauss elimination with scaled partial pivoting is used, but mult-
!  ipliers are  n o t  s a v e d  in order to save storage. rather, the
!  right side is operated on during elimination.  the two parameters
!                  i p v t e q   and  l a s t e q
!  are used to keep track of the action.  ipvteq is the index of the
!  variable to be eliminated next, from equations  ipvteq+1,...,lasteq,
!  using equation  ipvteq (possibly after an interchange) as the pivot
!  equation. the entries in the pivot column are  a l w a y s  in column
!  1 of  w . this is accomplished by putting the entries in rows
!  ipvteq+1,...,lasteq  revised by the elimination of the  ipvteq-th
!  variable one to the left in  w . in this way, the columns of the
!  equations in a given block (as stored in  w ) will be aligned with
!  those of the next block at the moment when these next equations be-
!  come involved in the elimination process.
!
!  thus, for the above example, the first elimination steps proceed
!  as follows.
!
!  *11 12 13 14    11 12 13 14    11 12 13 14    11 12 13 14
!  *21 22 23 24   *22 23 24       22 23 24       22 23 24
!  *31 32 33 34   *32 33 34      *33 34          33 34
!   43 44 45 46    43 44 45 46   *43 44 45 46   *44 45 46        etc.
!   53 54 55 56    53 54 55 56   *53 54 55 56   *54 55 56
!   66 67 68 69    66 67 68 69    66 67 68 69    66 67 68 69
!        .              .              .              .
!
!  In all other respects, the procedure is standard, including the
!  scaled partial pivoting.
!
  integer nbloks
  integer ncols
  integer nequ
!
  real awi1od
  real b(nequ)
  real colmax
  real d(nequ)
  integer i
  integer icount
  integer iflag
  integer ii
  integer integs(2,nbloks)
  integer ipvteq
  integer ipvtp1
  integer istar
  integer j
  integer jmax
  integer lastcl
  integer lasteq
  integer lasti
  integer nexteq
  integer nrowad
  real ratio
  real rowmax
  real sum
  real temp
  real w(nequ,ncols)
  real x(nequ)
!
  iflag = 1
  ipvteq = 0
  lasteq = 0
!
!  The I loop runs over the blocks.
!
  do i = 1, nbloks
!
!  The equations for the current block are added to those current-
!  ly involved in the elimination process, by increasing  lasteq
!  by  integs(1,i) after the rowsize of these equations has been
!  recorded in the array D.
!
    nrowad = integs(1,i)
    
    do icount = 1, nrowad

      nexteq = lasteq+icount
      
      rowmax = 0.0
      do j = 1, ncols
        rowmax = max ( rowmax, abs ( w(nexteq,j) ) )
      end do
      
      if ( rowmax == 0.0 ) then
        go to 150
      end if

      d(nexteq) = rowmax

    end do
   
    lasteq = lasteq+nrowad
!
!  There will be  lasti=integs(2,i)  elimination steps before
!  the equations in the next block become involved. further,
!  l a s t c l  records the number of columns involved in the cur-
!  rent elimination step. it starts equal to  ncols  when a block
!  first becomes involved and then drops by one after each elim-
!  ination step.
!
    lastcl = ncols
    lasti = integs(2,i)
    
    do icount = 1, lasti
    
      ipvteq = ipvteq+1

      if ( ipvteq < lasteq ) then
        go to 30
      end if

      if ( abs ( w(ipvteq,1))+d(ipvteq) > d(ipvteq) ) then
        go to 100
      end if

      go to 150
!
!  Determine the smallest  i s t a r  in  (ipvteq,lasteq)  for
!  which  abs(w(istar,1))/d(istar)  is as large as possible, and
!  interchange equations  ipvteq  and  istar  in case  ipvteq
!  < istar .
!
   30     continue

      colmax = abs(w(ipvteq,1)) / d(ipvteq)
      istar = ipvteq
      ipvtp1 = ipvteq+1
      
      do ii = ipvtp1, lasteq
        awi1od = abs(w(ii,1))/d(ii)
        if ( awi1od > colmax ) then
          colmax = awi1od
          istar = ii
        end if
      end do
      
      if ( abs(w(istar,1))+d(istar) == d(istar) ) then
        go to 150
      end if

      if ( istar == ipvteq ) then
        go to 60
      end if

      iflag = -iflag

      temp = d(istar)
      d(istar) = d(ipvteq)
      d(ipvteq) = temp

      temp = b(istar)
      b(istar) = b(ipvteq)
      b(ipvteq) = temp
      
      do j = 1, lastcl
        temp = w(istar,j)
        w(istar,j) = w(ipvteq,j)
        w(ipvteq,j) = temp
      end do
!
!  Subtract the appropriate multiple of equation  ipvteq  from
!  equations  ipvteq+1,...,lasteq to make the coefficient of the
!  ipvteq-th unknown (presently in column 1 of  w ) zero, but
!  store the new coefficients in  w  one to the left from the old.
!
   60     continue
   
      do ii = ipvtp1, lasteq
      
        ratio = w(ii,1)/w(ipvteq,1)
        do j = 2, lastcl
          w(ii,j-1) = w(ii,j)-ratio*w(ipvteq,j)
        end do
        w(ii,lastcl) = 0.0
        b(ii) = b(ii)-ratio*b(ipvteq)
        
      end do
   
      lastcl = lastcl-1
      
    end do
   
100     continue

  end do
!
!  At this point, W and B contain an upper triangular linear system
!  equivalent to the original one, with  w(i,j) containing entry
!  (i, i-1+j ) of the coefficient matrix. solve this system by backsub-
!  stitution, taking into account its block structure.
!
!  i-loop over the blocks, in reverse order
!
  i = nbloks

  110 continue

  lasti = integs(2,i)
  jmax = ncols-lasti
  
  do icount = 1, lasti
  
    sum = 0.0
    do j = 1, jmax
      sum = sum+x(ipvteq+j)*w(ipvteq,j+1)
    end do
  
    x(ipvteq) = (b(ipvteq)-sum)/w(ipvteq,1)
    jmax = jmax+1
    ipvteq = ipvteq-1
    
  end do
  
  i = i-1
  if ( i > 0 ) then
    go to 110
  end if

  return
  
  150 continue

  iflag = 0
  return
end
subroutine difequ ( mode, xx, v )
!
!*************************************************************************
!
!! DIFEQU returns information about a differential equation.
!
!
!  Parameters:
!
!    Input, integer MODE, an integer indicating the task to be performed.
!    1, initialization
!    2, evaluate  de  at  xx
!    3, specify the next side condition
!    4, analyze the approximation
!
!    Input, real XX, a point at which information is wanted
!
!    Output, real V, depends on the  mode  . see comments below
!
  integer, parameter :: npiece = 100
  integer, parameter :: ncoef = 2000
!
  real break
  real coef
  real eps
  real ep1
  real ep2
  real error
  real factor
  integer i
  integer iside
  integer itermx
  integer k
  integer kpm
  integer l
  integer m
  integer mode
  real rho
  real s2ovep
  real solutn
  real un
  real v(20)
  real value
  real x
  real xside
  real xx
!
  common /approx/ break(npiece),coef(ncoef),l,kpm
  common /side/ m,iside,xside(10)
  common /other/ itermx,k,rho(19)
!
!  this sample of  difequ  is for the example in chapter xv. it is a
!  nonlinear second order two point boundary value problem.
!
  go to (10,50,60,110), mode
!
!  initialize everything
!  i.e. set the order  m  of the dif.equ., the nondecreasing sequence
!  xside(i),i=1,...,m, of points at which side cond.s are given and
!  anything else necessary.
!
   10 continue

  m = 2
  xside(1) = 0.0
  xside(2) = 1.0
!
!  print out heading
!
  write(*,*)' '
  write(*,*)'Carrier''s nonlinear perturb. problem'
  write(*,*)' '
  
  eps = 0.005
  write(*,*)'EPS = ',eps
!
!  set constants used in formula for solution below.
!
  factor = (sqrt(2.0)+sqrt(3.0))**2
  s2ovep = sqrt(2.0/eps)
!
!  Initial guess for newton iteration. un(x)=x*x-1.
!
  l = 1
  break(1) = 0.0
  do i = 1, kpm
    coef(i) = 0.0
  end do
  coef(1) = -1.0
  coef(3) = 2.0
  itermx = 10
  return
!
!  Provide value of left side coeff.s and right side at  xx .
!  specifically, at  xx  the dif.equ. reads:
!
!    v(m+1)d**m+v(m)d**(m-1) + ... + v(1)d**0 = v(m+2)
!
!  in terms of the quantities v(i),i=1,...,m+2, to be computed here.
!
   50 continue

  v(3) = eps
  v(2) = 0.0
  call ppvalu(break,coef,l,kpm,xx,0,un)
  v(1) = 2.0*un
  v(4) = un**2+1.0
  return
!
!  provide the  m  side conditions. these conditions are of the form
!        v(m+1)d**m+v(m)d**(m-1) + ... + v(1)d**0 = v(m+2)
!  in terms of the quantities v(i),i=1,...,m+2, to be specified here.
!  note that v(m+1)=0  for customary side conditions.
!
   60 continue

  v(m+1) = 0.0
  if ( iside == 1 ) then
    v(2) = 1.0
    v(1) = 0.0
    v(4) = 0.0
    iside = iside+1
  else if ( iside == 2 ) then
    v(2) = 0.0
    v(1) = 1.0
    v(4) = 0.0
    iside = iside+1
  end if
  
  return
!
!  calculate the error near the boundary layer at  1.
!
  110 continue

  write(*,*)' '
  write(*,*)' X, G(X) and G(X)-F(X) at selected points:'
  write(*,*)' '

  x = 0.75
  
  do i = 1, 9
    ep1 = exp(s2ovep*(1.-x))*factor
    ep2 = exp(s2ovep*(1.+x))*factor
    solutn = 12./(1.+ep1)**2*ep1+12./(1.+ep2)**2*ep2-1.
    call ppvalu(break,coef,l,kpm,x,0,value)
    error = solutn-value
    write ( *, '(1x,3g14.6)' ) x, solutn, error
    x = x+0.03125
  end do
  
  return
end
subroutine dtblok ( bloks, integs, nbloks, ipivot, iflag, detsgn, detlog )
!
!*************************************************************************
!
!! DTBLOK gets the determinant of an almost block diagonal matrix.
!
!
!  The matrix's PLU factorization must have been obtained 
!  previously by FCBLOK.
!
!  The logarithm of the determinant is computed instead of the
!  determinant itself to avoid the danger of overflow or underflow
!  inherent in this calculation.
!
!  Parameters:
!
!  bloks, integs, nbloks, ipivot, iflag  are as on return from fcblok.
!            in particular, iflag=(-1)**(number of interchanges dur-
!            ing factorization) if successful, otherwise iflag=0.
!
!  detsgn  on output, contains the sign of the determinant.
!
!  detlog  on output, contains the natural logarithm of the determi-
!            nant if determinant is not zero. otherwise contains 0.
!
  integer nbloks
!
  real bloks(1)
  real detlog
  real detsgn
  integer i
  integer iflag
  integer index
  integer indexp
  integer integs(3,nbloks)
  integer ip
  integer ipivot(1)
  integer k
  integer last
  integer nrow
!
  detsgn = iflag
  detlog = 0.0

  if ( iflag == 0 ) then
    return
  end if

  index = 0
  indexp = 0
  
  do i = 1, nbloks
  
    nrow = integs(1,i)
    last = integs(3,i)
    
    do k = 1, last
      ip = index+nrow*(k-1)+ipivot(indexp+k)
      detlog = detlog+alog(abs(bloks(ip)))
      detsgn = detsgn*sign(1.0,bloks(ip))
    end do
   
    index = nrow*integs(2,i)+index
    indexp = indexp+nrow
    
  end do
   
  return
end
subroutine eqblok ( t, n, kpm, work1, work2, bloks, lenblk, integs, &
  nbloks, b )
!
!*************************************************************************
!
!! EQBLOK is to be called in COLLOC.
!
!
!  Method:
!
!    each breakpoint interval gives rise to a block in the linear system.
!    this block is determined by the  k  colloc.equations in the interval
!    with the side conditions (if any) in the interval interspersed ap-
!    propriately, and involves the  kpm  b-splines having the interval in
!    their support. correspondingly, such a block has  nrow=k+isidel
!    rows, with  isidel=number of side conditions in this and the prev-
!    ious intervals, and  ncol=kpm  columns.
!
!    further, because the interior knots have multiplicity  k, we can
!    carry out (in slvblk)  k  elimination steps in a block before pivot-
!    ing might involve an equation from the next block. in the last block,
!    of course, all kpm elimination steps will be carried out (in slvblk).
!
!    see the detailed comments in the solveblok package for further in-
!    formation about the almost block diagonal form used here.
!
!  Parameters:
!
!  input 
!
!    Input, real T(N+KPM), the knot sequence.
!
!  n   the dimension of the approximating spline space, i.e., the order
!      of the linear system to be constructed.
!
!  kpm=k+m, the order of the approximating spline
!
!  lenblk   the maximum length of the array  bloks  as allowed by the
!           dimension statement in  colloc .
!
!  work  a r e a s 
!
!  work1    used in  putit, of size (kpm,kpm)
!  work2    used in  putit, of size (kpm,m+1)
!
!  output 
!
!  bloks    the coefficient matrix of the linear system, stored in al-
!           most block diagonal form, of size
!              kpm*sum(integs(1,i) , i=1,...,nbloks)
!
!  integs   an integer array, of size (3,nbloks), describing the block
!           structure.
!           integs(1,i) = number of rows in block  i
!           integs(2,i) = number of columns in block  i
!           integs(3,i) = number of elimination steps which can be
!                       carried out in block  i  before pivoting might
!                       bring in an equation from the next block.
!
!  nbloks number of blocks, equals number of polynomial pieces
!
!  b      the right side of the linear system, stored corresponding to the
!         almost block diagonal form, of size sum(integs(1,i) , i=1,...,
!         nbloks).
!
  integer kpm
  integer n
!
  real b(*)
  real bloks(*)
  integer i
  integer index
  integer indexb
  integer integs(3,*)
  integer iside
  integer isidel
  integer itermx
  integer k
  integer left
  integer lenblk
  integer m
  integer nbloks
  integer nrow
  real rho
  real t(n+kpm)
  real work1(kpm,kpm)
  real work2(kpm,*)
  real xside
!
  common /side/ m,iside,xside(10)
  common /other/ itermx,k,rho(19)
!
  index = 1
  indexb = 1
  i = 0
  iside = 1
  
  do left = kpm, n, k
  
    i = i+1
!
!  determine integs(.,i)
!
    integs(2,i) = kpm
    
    if ( left >= n ) then
      integs(3,i) = kpm
      isidel = m
      go to 30
    end if
    
    integs(3,i) = k
!
!  At this point,  iside-1  gives the number of side conditions
!  incorporated so far. adding to this the side conditions in the
!  current interval gives the number  isidel .
!
    isidel = iside-1

   20   continue

    if ( isidel == m ) then
      go to 30
    end if

    if ( xside(isidel+1) >= t(left+1) ) then
      go to 30
    end if

    isidel = isidel+1
    go to 20

   30   continue

    nrow = k+isidel
    integs(1,i) = nrow
!
!  the detailed equations for this block are generated and put
!  together in PUTIT.
!
    if ( lenblk < index+nrow*kpm-1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'EQBLOK - Fatal error!'
      write ( *, * ) '  The dimension of BLOKS is too small.'
      write ( *, * ) '  LENBLK = ', lenblk
      stop
    end if

    call putit(t,kpm,left,work1,work2,bloks(index),nrow,b(indexb))
    index = index+nrow*kpm
    indexb = indexb+nrow
    
  end do
 
  nbloks = i

  return
end
subroutine evnnot ( break, coef, l, k, brknew, lnew, coefg )
!
!*************************************************************************
!
!! EVNNOT is a "fake" version of NEWNOT.
!
!
!  Discussion:
!
!    EVNNOT returns  lnew+1  knots in  brknew  which are 
!    equidistributed on (a,b) =(break(1),break(l+1)) .
!
!  Parameters:
!
!    Input, real BREAK(L+1), coef, l, k.....contains the pp-representation 
!    of a certain function F of order K.  Specifically,
!    d**(k-1)f(x)=coef(k,i)  for  break(i) <= x < break(i+1)
!
!    Input, integer LNEW, the number of subintervals into which the interval 
!    (a,b) is to be sectioned by the new breakpoint sequence  brknew .
!
!    Output, brknew array of length  lnew+1  containing the new breakpoint 
!    sequence.
!
!    Output, coefg  the coefficient part of the pp-repr.  break, coefg, l, 2
!    for the monotone p.linear function  g  wrto which  brknew  will
!    be equidistributed.
!
  integer k
  integer l
!
  real break(l+1)
  real brknew(lnew+1)
  real coef(k,l)
  real coefg(2,l)
  integer i
  integer lnew
  real step
!
  brknew(1) = break(1)
  brknew(lnew+1) = break(l+1)
  step = (break(l+1)-break(1))/real ( lnew)
 
  do i = 2, lnew
    brknew(i) = break(1)+real ( i-1)*step
  end do
 
  return
end
subroutine factrb ( w, ipivot, d, nrow, ncol, last, iflag )
!
!*************************************************************************
!
!! FACTRB constructs a partial PLU factorization.
!
!
!  Discussion:
!
!    This factorization corresponds to steps 1 through LAST in Gauss 
!    elimination for the matrix W of order ( NROW, NCOL ), using 
!    pivoting of scaled rows.
!
!  Parameters:
!
!  w       contains the (nrow,ncol) matrix to be partially factored
!          on input, and the partial factorization on output.
!
!  ipivot  an integer array of length nrow containing a record of the
!          pivoting strategy used; row ipivot(i) is used during the
!          i-th elimination step, i=1,...,last.
!
!  d       a work array of length nrow used to store row sizes
!          temporarily.
!
!  nrow    number of rows of w.
!
!  ncol    number of columns of w.
!
!  last    number of elimination steps to be carried out.
!
!  iflag   on output, equals iflag on input times (-1)**(number of
!          row interchanges during the factorization process), in
!          case no zero pivot was encountered.
!          otherwise, iflag=0 on output.
!
  integer ncol
  integer nrow
!
  real awikdi
  real colmax
  real d(nrow)
  integer i
  integer iflag
  integer ipivi
  integer ipivk
  integer ipivot(nrow)
  integer j
  integer k
  integer kp1
  integer last
  real ratio
  real rowmax
  real w(nrow,ncol)
!
!  Initialize IPIVOT and D.
!
  do i = 1, nrow
    ipivot(i) = i
  end do

  do i = 1, nrow
    
    rowmax = 0.0
    do j = 1, ncol
      rowmax = max ( rowmax, abs ( w(i,j) ) )
    end do
    
    if ( rowmax == 0.0 ) then
      iflag = 0
      return
    end if

    d(i) = rowmax
    
  end do
!
!  Gauss elimination with pivoting of scaled rows, loop over k=1,.,last
!
  k = 1
!
!  As pivot row for k-th step, pick among the rows not yet used,
!  i.e., from rows ipivot(k),...,ipivot(nrow), the one whose k-th
!  entry (compared to the row size) is largest. then, if this row
!  does not turn out to be row ipivot(k), redefine ipivot(k) ap-
!  propriately and record this interchange by changing the sign
!  of IFLAG.
!
   30 continue

  ipivk = ipivot(k)

  if ( k == nrow ) then
    if ( abs(w(ipivk,nrow))+d(ipivk) <= d(ipivk) ) then
      iflag = 0
    end if
    return
  end if

  j = k
  kp1 = k+1
  colmax = abs(w(ipivk,k))/d(ipivk)
!
!  Find the largest pivot
!
  do i = kp1, nrow
    ipivi = ipivot(i)
    awikdi = abs(w(ipivi,k))/d(ipivi)
    if ( awikdi > colmax ) then
      colmax = awikdi
      j = i
    end if
  end do
  
  if ( j /= k ) then
    ipivk = ipivot(j)
    ipivot(j) = ipivot(k)
    ipivot(k) = ipivk
    iflag = -iflag
  end if
!
!  If pivot element is too small in absolute value, declare
!  matrix to be noninvertible and quit.
!
  if ( abs(w(ipivk,k))+d(ipivk) <= d(ipivk) ) then
    iflag = 0
    return
  end if
!
!  Otherwise, subtract the appropriate multiple of the pivot
!  row from remaining rows, i.e., the rows ipivot(k+1),...,
!  ipivot(nrow), to make k-th entry zero. save the multiplier in
!  its place.
!
  do i = kp1, nrow
  
    ipivi = ipivot(i)
    w(ipivi,k) = w(ipivi,k)/w(ipivk,k)
    
    ratio = -w(ipivi,k)
    do j = kp1, ncol
      w(ipivi,j) = ratio*w(ipivk,j)+w(ipivi,j)
    end do
    
  end do
   
  k = kp1
!
!  Check for having reached the next block.
!
  if ( k <= last ) then
    go to 30
  end if

  return
end
subroutine fcblok ( bloks, integs, nbloks, ipivot, scrtch, iflag )
!
!*************************************************************************
!
!! FCBLOK supervises the PLU factorization with pivoting of
!  scaled rows of the almost block diagonal matrix stored in the arrays
!  b l o k s  and  i n t e g s .
!
!  The FACTRB routine carries out steps 1,...,last of gauss
!  elimination (with pivoting) for an individual block.
!
!  The SHIFTB routine shifts the remaining rows to the top of
!            the next block
!
!  Parameters:
!
!  bloks   an array that initially contains the almost block diagonal
!            matrix  a  to be factored, and on return contains the com-
!            puted factorization of  a .
!
!  integs  an integer array describing the block structure of  a .
!
!  nbloks  the number of blocks in  a .
!
!  ipivot  an integer array of dimension  sum (integs(1,n) ; n=1,
!            ...,nbloks) which, on return, contains the pivoting stra-
!            tegy used.
!
!  scrtch  work area required, of length  max (integs(1,n) ; n=1,
!            ...,nbloks).
!
!  iflag   output parameter;
!          =0  in case matrix was found to be singular.
!            otherwise,
!          =(-1)**(number of row interchanges during factorization)
!
  integer nbloks
!
  real bloks(*)
  integer i
  integer iflag
  integer index
  integer indexb
  integer indexn
  integer integs(3,nbloks)
  integer ipivot(*)
  integer last
  integer ncol
  integer nrow
  real scrtch(*)
!
  iflag = 1
  indexb = 1
  indexn = 1
  i = 1 
!
!  Loop over the blocks.  i  is loop index
!
   10 continue

  index = indexn
  nrow = integs(1,i)
  ncol = integs(2,i)
  last = integs(3,i)
!
!  Carry out elimination on the i-th block until next block
!  enters, i.e., for columns 1,...,last  of i-th block.
!
  call factrb(bloks(index),ipivot(indexb),scrtch,nrow,ncol,last, iflag)
!
!  Check for having reached a singular block or the last block
!
  if ( iflag == 0 .or. i == nbloks ) then
    return
  end if

  i = i+1
  indexn = nrow*ncol+index
!
!  Put the rest of the i-th block onto the next block.
!
  call shiftb (bloks(index),ipivot(indexb),nrow,ncol,last, &
    bloks(indexn),integs(1,i),integs(2,i))
     
  indexb = indexb + nrow
  go to 10
end
subroutine spline_hermite_set ( ndata, tdata, c )
!
!*************************************************************************
!
!! SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
!
!
!  Reference:
!
!    Conte and de Boor,
!    Algorithm CALCCF,
!    Elementary Numerical Analysis, 
!    1973, page 235.
!
!  Modified:
!
!    06 April 1999
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.  
!    NDATA must be at least 2.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.  
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input/output, real C(4,NDATA).
!
!    On input, C(1,I) and C(2,I) should contain the value of the
!    function and its derivative at TDATA(I), for I = 1 to NDATA.
!    These values will not be changed by this routine.
!
!    On output, C(3,I) and C(4,I) contain the quadratic
!    and cubic coefficients of the Hermite polynomial
!    in the interval (TDATA(I), TDATA(I+1)), for I=1 to NDATA-1.
!    C(3,NDATA) and C(4,NDATA) are set to 0.
!
!    In the interval (TDATA(I), TDATA(I+1)), the interpolating Hermite 
!    polynomial is given by
!
!    SVAL(TVAL) =                C(1,I)
!       + ( TVAL - TDATA(I) ) * ( C(2,I)
!       + ( TVAL - TDATA(I) ) * ( C(3,I)
!       + ( TVAL - TDATA(I) ) *   C(4,I) ) )
!
  integer ndata
!
  real c(4,ndata)
  real divdif1
  real divdif3
  real dt
  integer i
  real tdata(ndata)
!
  do i = 1, ndata-1
    dt = tdata(i+1) - tdata(i)
    divdif1 = ( c(1,i+1) - c(1,i) ) / dt
    divdif3 = c(2,i) + c(2,i+1) - 2.0 * divdif1
    c(3,i) = ( divdif1 - c(2,i) - divdif3 ) / dt
    c(4,i) = divdif3 / dt**2
  end do

  c(3,ndata) = 0.0
  c(4,ndata) = 0.0

  return
end
subroutine spline_hermite_val ( ndata, tdata, c, tval, sval )
!
!*************************************************************************
!
!! SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
!
!
!  Discussion:
!
!    SPLINE_HERMITE_SET must be called first, to set up the
!    spline data from the raw function and derivative data.
!
!  Reference:
!
!    Conte and de Boor,
!    Algorithm PCUBIC,
!    Elementary Numerical Analysis, 
!    1973, page 234.
!
!  Modified:
!
!    06 April 1999
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.  
!    NDATA is assumed to be at least 2.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.  
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input, real C(4,NDATA), contains the data computed by 
!    SPLINE_HERMITE_SET.
!
!    Input, real TVAL, the point where the interpolant is to 
!    be evaluated.  
!
!    Output, real SVAL, the value of the interpolant at TVAL.
!
  integer ndata
!
  real c(4,ndata)
  real dt
  integer i
  integer j
  real sval
  real tdata(ndata)
  real tval
!
!  Find the interval J = [ TDATA(J), TDATA(J+1) ] that contains 
!  or is nearest to TVAL.
!
  j = ndata - 1

  do i = 1, ndata-2

    if ( tval < tdata(i+1) ) then
      j = i
      exit
    end if

  end do
!
!  Evaluate the cubic polynomial.
!
  dt = tval - tdata(j)

  sval = c(1,j) + dt * ( c(2,j) + dt * ( c(3,j) + dt * c(4,j) ) )

  return
end
subroutine interv ( xt, lxt, x, left, mflag )
!
!*******************************************************************************
!
!! INTERV computes LEFT, the maximum value of I so that
!
!    1 < = I  <= XT
!
!  and
!
!    XT(I) < = X.
!
!  The routine is designed to be efficient in the common situation
!  that it is called repeatedly, with X taken from an increasing
!  or decreasing sequence.
!
!  This will happen when a piecewise polynomial is to be graphed.
!  The first guess for LEFT is therefore taken to be the value
!  returned at the previous call and stored in the local variable
!  ILO.
!
!  A first check ascertains that ILO.LT.LXT.  This is necessary
!  since the present call may have nothing to do with the previous
!  call.  Then, if XT(ILO) < = X < XT(ILO+1), we set LEFT=ILO
!  and are done after just three comparisons.
!
!  Otherwise, we repeatedly double the difference ISTEP=IHI-ILO
!  while also moving ILO and IHI in the direction of X, until
!    XT(ILO) < = X < XT(IHI)
!  after which we use bisection to get, in addition, ILO+1=IHI.
!  LEFT=ILO is then returned.
!
!  Parameters:
!
!    Input, real XT(LXT), a nondecreasing sequence of values.
!
!    Input, integer LXT, the dimension of XT.
!
!    Input, real X, the point whose location with 
!    respect to the sequence XT is to be determined.
!
!    Output, integer LEFT, integer MFLAG, whose value is
!
!      1     -1      if               X <  XT(1)
!      I      0      if   XT(I)  <= X  < XT(I+1)
!      LXT    1      if  XT(LXT) <= X
!
!    In particular, MFLAG=0 is the 'usual' case.  MFLAG /=0
!    indicates that X lies outside the half open interval
!    XT(1) <= Y < XT(LXT).  The asymmetric treatment of the
!    interval is due to the decision to make all piecewise
!    polynomials continuous from the right.
!
  integer lxt
!
  integer left
  integer mflag
  integer ihi
  integer ilo
  integer istep
  integer middle
  real x
  real xt(lxt)
!
  save ilo
!
  data ilo / 1 /
!
  ihi = ilo+1

  if ( ihi >= lxt ) then

    if ( x >= xt(lxt) ) then
      go to 110
    end if

    if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
    end if

    ilo = lxt-1
    ihi = lxt

  end if

  if ( x >= xt(ihi) ) then
    go to 40
  end if

  if ( x >= xt(ilo) ) then
    mflag = 0
    left = ilo
    return
  end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
  istep = 1

   31 continue

  ihi = ilo
  ilo = ihi - istep

  if ( ilo > 1 ) then
    if (x >= xt(ilo) ) then
      go to 50
    end if
    istep = istep*2
    go to 31
  end if

  ilo = 1

  if ( x < xt(1) ) then
    mflag = -1
    left = 1
    return
  end if

  go to 50
!
!  Now X => XT(IHI).  Increase IHI to capture X.
!
   40 continue

  istep = 1

   41 continue

  ilo = ihi
  ihi = ilo + istep

  if ( ihi < lxt ) then
    if ( x < xt(ihi) ) then
      go to 50
    end if
    istep=istep * 2
    go to 41
  end if

  if (x >= xt(lxt) ) then
    go to 110
  end if

  ihi = lxt
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
   50 continue

  middle = (ilo + ihi)/2

  if ( middle == ilo ) then
    mflag = 0
    left = ilo
    return
  end if
!
!  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
!
  if ( x >= xt(middle) ) then
    ilo = middle
  else
    ihi = middle
  end if

  go to 50
!
!  Set output and return.
!
  110 continue

  mflag = 1

  if ( x == xt(lxt) ) then
    mflag = 0
  end if

  do left = lxt, 1, -1
    if ( xt(left) < xt(lxt) ) then
      return
    end if
  end do

  return
end
subroutine knots ( break, l, kpm, t, n )
!
!*************************************************************************
!
!! KNOTS is to be called in COLLOC.
!
!
!  Discussion:
!
!    The routine constructs from the given breakpoint sequence BREAK the knot
!    sequence T so that
!
!      spline(k+m,t)=pp(k+m,break) 
!
!    with  m-1  continuous derivatives.  This means that
!
!      t(1),...,t(n+kpm) = break(1) kpm times, then break(2),...,
!        break(l) each  k  times, then, finally, break(l+1) kpm times.
!
!  Parameters:
!
!    Input, real BREAK(L+1), the breakpoint sequence.
!
!    Input, integer L, the number of intervals or pieces.
!
!    Input, kpm =k+m, order of the pp function or spline
!
!    Output, real T(N+KPM), the knot sequence.
!
!    Output, n   =l*k+m  = dimension of  spline(k+m,t).
!
  integer kpm
  integer l
  integer n
!
  real break(l+1)
  integer iside
  integer j
  integer jj
  integer jjj
  integer k
  integer kpm
  integer ll
  integer m
  real t(*)
  real xside
!
  common /side/ m,iside,xside(10)
!
  k = kpm-m
  n = l*k+m
  jj = n+kpm
  jjj = l+1
  
  do ll = 1, kpm
    t(jj) = break(jjj)
    jj = jj-1
  end do
  
  do j = 1, l
    jjj = jjj-1
    do ll = 1, k
      t(jj) = break(jjj)
      jj = jj-1
    end do
  end do
   
  do ll = 1, kpm
    t(ll) = break(1)
  end do
   
  return
end
subroutine l2appr ( t, n, k, q, diag, bcoef )
!
!*************************************************************************
!
!! L2APPR constructs the (weighted discrete) l2-approximation by 
!  splines of order k  with knot sequence  t(1), ..., t(n+k)  to 
!  given data points
!  ( tau(i), gtau(i) ), i=1,...,ntau. the b-spline coefficients
!  b c o e f   of the approximating spline are determined from the
!  normal equations using cholesky's method.
!
!  Method:
!
!  the b-spline coefficients of the l2-appr. are determined as the sol-
!  ution of the normal equations
!     sum ( (b(i),b(j))*bcoef(j) : j=1,...,n) =(b(i),g),
!                                               i=1, ..., n .
!  here,  b(i)  denotes the i-th b-spline,  g  denotes the function to
!  be approximated, and the  i n n e r   p r o d u c t  of two funct-
!  ions  f  and  g  is given by
!      (f,g)  := sum ( f(tau(i))*g(tau(i))*weight(i) : i=1,...,ntau) .
!  the arrays  t a u  and  w e i g h t  are given in common block
!   d a t a , as is the array  g t a u  containing the sequence
!  g(tau(i)), i=1,...,ntau.
!  the relevant function values of the b-splines  b(i), i=1,...,n, are
!  supplied by the subprogram  b s p l v b .
!     the coeff.matrix  c , with
!           c(i,j)  := (b(i), b(j)), i,j=1,...,n,
!  of the normal equations is symmetric and (2*k-1)-banded, therefore
!  can be specified by giving its k bands at or below the diagonal. for
!  i=1,...,n,  we store
!   (b(i),b(j)) = c(i,j)  in  q(i-j+1,j), j=i,...,min(i+k-1,n)
!  and the right side
!   (b(i), g )  in  bcoef(i) .
!  since b-spline values are most efficiently generated by finding sim-
!  ultaneously the value of  e v e r y  nonzero b-spline at one point,
!  the entries of  c  (i.e., of  q ), are generated by computing, for
!  each ll, all the terms involving  tau(ll)  simultaneously and adding
!  them to all relevant entries.
!
!  input
!
!    t(1), ..., t(n+k)  the knot sequence
!    n.....the dimension of the space of splines of order k with knots t.
!    k.....the order
!
!    Warning: The restriction   k <= kmax (= 20) is imposed 
!    by the arbitrary dimension statement for  biatx  below, but
!    is  n o w h e r e   c h e c k e d   for.
!
!  work arrays
!
!  q....a work array of size (at least) k*n. its first  k  rows are used
!       for the  k  lower diagonals of the gramian matrix  c .
!  diag.....a work array of length  n  used in bchfac .
!
!  input  via  c o m m o n  /data/ 
!  ntau.....number of data points
!  (tau(i),gtau(i)), i=1,...,ntau     are the  ntau  data points to be
!        fitted .
!  weight(i), i=1,...,ntau    are the corresponding weights .
!
!  output 
!  bcoef(1), ..., bcoef(n)  the b-spline coeffs. of the l2-appr.
!
  integer, parameter :: kmax = 20
  integer, parameter :: ntmax = 200
!
  real bcoef(n)
  real biatx(kmax)
  real diag(n)
  real dw
  real gtau
  integer i
  integer j
  integer jj
  integer k
  integer left
  integer leftmk
  integer ll
  integer mm
  integer n
  integer ntau
  real q(k,n)
  real t(n+k)
  real tau
  real totalw
  real weight
!
  common /data/ ntau,tau(ntmax),gtau(ntmax),weight(ntmax),totalw
!
  do j = 1, n
    bcoef(j) = 0.0
    do i = 1, k
      q(i,j) = 0.0
    end do
  end do

  left = k
  leftmk = 0
  
  do ll = 1, ntau
!
!  locate LEFT s.t. tau(ll) in (t(left),t(left+1))
!
   20   continue

    if ( left == n ) then
      go to 30
    end if

    if ( tau(ll) < t(left+1) ) then
      go to 30
    end if

    left = left+1
    leftmk = leftmk+1
    go to 20

   30   continue

    call bsplvb (t,k,1,tau(ll),left,biatx)
!
!  biatx(mm) contains the value of b(left-k+mm) at tau(ll).
!  hence, with  dw := biatx(mm)*weight(ll), the number dw*gtau(ll)
!  is a summand in the inner product
!     (b(left-k+mm), g)  which goes into  bcoef(left-k+mm)
!  and the number biatx(jj)*dw is a summand in the inner product
!     (b(left-k+jj), b(left-k+mm)), into  q(jj-mm+1,left-k+mm)
!  since  (left-k+jj)-(left-k+mm)+1 = jj - mm + 1 .
!
    do mm = 1, k
    
      dw = biatx(mm)*weight(ll)
      j = leftmk+mm
      bcoef(j) = dw*gtau(ll)+bcoef(j)
      i = 1
      
      do jj = mm, k
        q(i,j) = biatx(jj)*dw+q(i,j)
        i = i+1
      end do
      
    end do
    
  end do
!
!  Construct the Cholesky factorization for  c  in  q , then 
!  use it to solve the normal equations
!    c*x = bcoef
!  for x , and store  x  in  bcoef .
!
  call bchfac(q,k,n,diag)
  
  call bchslv(q,k,n,bcoef)
  
  return
end
subroutine l2err ( iprfun, ftau, error )
!
!*************************************************************************
!
!! L2ERR computes errors of an L2 approximation.
!
!
!  This routine computes various errors of the current l2-approxi-
!  mation , whose pp-repr. is contained in common block  approx  ,
!  to the given data contained in common block  data . it prints out
!  the average error  e r r l 1 , the l2-error  e r r l 2,  and the
!  maximum error  e r r m a x .
!
!  Parameters: 
!
!    Input, iprfun  an integer.  if iprfun= 1, the routine prints out
!          the value of the approximation as well as its error at
!          every data point.
!
!    Output, ftau(1), ..., ftau(ntau),  with  ftau(i)  the approximation  f at
!          tau(i), all i.
!
!    Output, error(1), ..., error(ntau),  with  error(i)=scale*(g-f)
!          at tau(i), all i. here,  s c a l e  equals  1. in case
!          iprfun /= 1 , or the abs.error is greater than 100 some-
!          where. otherwise, s c a l e  is such that the maximum of
!          abs(error))  over all  i  lies between  10  and  100. this
!          makes the printed output more illustrative.
!
  integer, parameter :: lpkmax = 100
  integer, parameter :: ntmax = 200
  integer, parameter :: ltkmax = 2000
!
  integer ntau
!
  real break
  real coef
  real err
  real errl1
  real errl2
  real errmax
  real error(ntau)
  real ftau(ntau)
  real gtau
  integer ie
  integer iprfun
  integer k
  integer l
  integer ll
  real scale
  real tau
  real totalw
  real weight
!
  common /data/ ntau,tau(ntmax),gtau(ntmax),weight(ntmax),totalw
  common /approx/ break(lpkmax),coef(ltkmax),l,k
!
  errl1 = 0.0
  errl2 = 0.0
  errmax = 0.0

  do ll = 1, ntau
    call ppvalu(break,coef,l,k,tau(ll),0,ftau(ll))
    error(ll) = gtau(ll)-ftau(ll)
    err = abs(error(ll))
    if ( errmax < err ) then
      errmax = err
    end if
    errl1 = errl1+err * weight(ll)
    errl2 = errl2+err**2 * weight(ll)
  end do

  errl1 = errl1/totalw
  errl2 = sqrt(errl2/totalw)

  write(*,*)' '
  write(*,*)' Least square error =',errl2
  write(*,*)' Average error     =',errl1
  write(*,*)' Maximum error     =',errmax
  write(*,*)' '
  
  if ( iprfun /= 1 ) then
    return
  end if
!
!  Scale error curve and print
!
  ie = 0
  scale = 1.0

  if ( errmax < 10.0 ) then
  
    do ie = 1, 9
      scale = scale*10.0
      if ( errmax*scale >= 10.0 ) then
        exit
      end if
    end do

  end if  

  do ll = 1, ntau
    error(ll) = error(ll)*scale
  end do
  
  write(*,60)ie,(ll,tau(ll),ftau(ll),error(ll),ll=1,ntau)
   60 format (///14x,'approximation and scaled error curve'/ &
       7x,'data point',7x,'approximation',3x,'deviation x 10**',i1/ &
       (i4, f16.8,f16.8,f17.6))

  return
end
subroutine l2knts ( break, l, k, t, n )
!
!*************************************************************************
!
!! L2KNTS converts breakpoints to knots.
!
!
!  Discussion:
!
!    The breakpoint sequence BREAK is converted into a corresponding 
!    knot sequence T to allow the representation of a piecewise
!    polynomial function of order  k  with
!    k-2 continuous derivatives as a spline of order  k  with knot
!    sequence  t . this means that
!    t(1), ..., t(n+k)= break(1) k times, then break(i), i=2,...,l, each
!    once, then break(l+1) k times.  Therefore,  n=k-1+l.
!
!  Parameters:
!
!    Input, k, the order.
!
!    Input, l, the number of polynomial pieces.
!
!    Input, break(1), ...,break(l+1), the breakpoint sequence.
!
!    Output, t(1),...,t(n+k), the knot sequence.
!
!    Output, n, the dimension of the corresp. spline space of order k.
!
  integer k
  integer l
  integer n
!
  real break(l+1)
  integer i
  real t(k-1+l+k)
!
  do i = 1, k-1
    t(i) = break(1)
  end do
 
  do i = 1, l
    t(k-1+i) = break(i)
  end do
 
  n = k-1+l
 
  do i = 1, k
    t(n+i) = break(l+1)
  end do
 
  return
end
subroutine newnot ( break, coef, l, k, brknew, lnew, coefg )
!
!*************************************************************************
!
!! NEWNOT returns  lnew+1  knots which are equidistributed on (a,b).
!
!
!   (a,b) =(break(1),break(l+1)) wrto a certain monotone fctn  g  related to
!  the k-th root of the k-th derivative of the pp function  f  whose pp-
!  representation is contained in  break, coef, l, k .
!
!  method 
!
!  the k-th derivative of the given pp function  f  does not exist
!  (except perhaps as a linear combination of delta functions). never-
!  theless, we construct a p.constant function  h  with breakpoint se-
!  quence  break  which is approximately proportional to abs(d**k(f)).
!  specifically, on  (break(i), break(i+1)),
!
!     abs(jump at break(i) of pc)    abs(jump at break(i+1) of pc)
!  h=-------------- + ----------------------------
!       break(i+1)-break(i-1)         break(i+2) - break(i)
!
!  with  pc  the p.constant (k-1)st derivative of  f .
!      then, the p.linear function  g  is constructed as
!
!    g(x) = integral of  h(y)**(1/k)  for  y  from  a  to  x
!
!  and its pp coeffs. stored in  coefg .
!
!  then  brknew  is determined by
!
!        brknew(i) = a+g**(-1)((i-1)*step) , i=1,...,lnew+1
!
!  where  step=g(b)/lnew  and  (a,b) = (break(1),break(l+1)) .
!  In the event that  pc=d**(k-1)(f) is constant in  (a,b)  and
!  therefore  h=0 identically,  brknew  is chosen uniformly spaced.
!
!  optional  p r i n t e d  output 
!  coefg.....the pp coeffs of  g  are printed out if  iprint  is set
!        > 0  in data statement below.
!
!  input
!
!  break, coef, l, k.....contains the pp-representation of a certain
!        function  f  of order  k . specifically,
!        d**(k-1)f(x)=coef(k,i)  for  break(i) <= x < break(i+1)
!  lnew.....number of intervals into which the interval (a,b) is to be
!        sectioned by the new breakpoint sequence  brknew .
!
!  output 
!  brknew.....array of length  lnew+1  containing the new breakpoint se-
!        quence
!  coefg.....the coefficient part of the pp-repr.  break, coefg, l, 2
!        for the monotone p.linear function  g  wrto which  brknew  will
!        be equidistributed.
!
  integer k
  integer l
!
  real break(l+1)
  real brknew(lnew+1)
  real coef(k,l)
  real coefg(2,l)
  real dif
  real difprv
  integer i
  integer iprint
  integer j
  integer lnew
  real oneovk
  real step
  real stepi
!
  data iprint / 0 /
!
  brknew(1) = break(1)
  brknew(lnew+1) = break(l+1)
!
!  If G is constant,  brknew  is uniform.
!
  if ( l <= 1) then
    go to 70
  end if
!
!  Construct the continuous p.linear function  g .
!
  oneovk = 1.0/real ( k)
  coefg(1,1) = 0.0
  difprv = abs(coef(k,2)-coef(k,1))/(break(3)-break(1))
  
  do i = 2, l
    dif = abs(coef(k,i)-coef(k,i-1))/(break(i+1)-break(i-1))
    coefg(2,i-1) = (dif+difprv)**oneovk
    coefg(1,i) = coefg(1,i-1)+coefg(2,i-1)*(break(i)-break(i-1))
    difprv = dif
  end do
   
  coefg(2,l) = (2.0*difprv)**oneovk
!
!  step = g(b)/lnew
!
  step=(coefg(1,l)+coefg(2,l)*(break(l+1)-break(l)))/real ( lnew)

  if ( iprint > 0 ) then
    write(*,20)step,(i,coefg(1,i),coefg(2,i),i=1,l)
  end if

   20 format (' step =',e16.7/(i5,2e16.5))
!
!  if  g  is constant,  brknew  is uniform .
!
  if ( step <= 0.0 ) then
    go to 70
  end if
!
!  for i=2,...,lnew, construct  brknew(i)=a+g**(-1)(stepi),
!  with  stepi=(i-1)*step .  this requires inversion of the p.lin-
!  ear function  g .  for this,  j  is found so that
!    g(break(j)) <= stepi .le. g(break(j+1))
!  and then
!    brknew(i) = break(j)+(stepi-g(break(j)))/dg(break(j)) .
!  the midpoint is chosen if  dg(break(j))=0 .
!
  j = 1
  
  do i = 2, lnew
  
    stepi = real ( i-1)*step

   30   continue

    if ( j == l ) then
      go to 40
    end if

    if ( stepi <= coefg(1,j+1) ) then
      go to 40
    end if

    j = j+1
    go to 30
    
   40   continue
   
    if ( coefg(2,j) /= 0.0 ) then
      brknew(i) = break(j)+(stepi-coefg(1,j))/coefg(2,j)
    else
      brknew(i) = (break(j)+break(j+1))/2.
    end if
    
  end do
  
  return
!
!  If G is constant, BRKNEW is uniform.
!
   70 step = (break(l+1)-break(1))/real ( lnew)
 
  do i = 2, lnew
    brknew(i) = break(1)+real ( i-1)*step
  end do
 
  return
end
subroutine ppvalu ( break, coef, l, k, x, jderiv, value )
!
!*******************************************************************************
!
!! PPVALU calculates the value at X of the JDERIV-th derivative of
!  the piecewise polynomial function F from its piecewise
!  polynomial representation.
!
!  The interval index I, appropriate for X, is found through a
!  call to INTERV.  The formula above for the JDERIV-th derivative
!  of F is then evaluated by nested multiplication.
!
!  The J-th derivative of F is given by:
!
!    (d**j)f(x) = 
!      coef(j+1,i) + h * (
!      coef(j+2,i) + h * (
!      ...
!      coef(k-1,i) + h * (
!      coef(k,i) / (k-j-1) ) / (k-j-2) ... ) / 2 ) / 1
!
!  with
!
!    H=X-BREAK(I)
!
!  and
!
!    i = max( 1 , max( j ,  break(j) <= x , 1 .le. j .le. l ) ).
!
!  Parameters:
!
!    Input, real BREAK(L+1), real COEF(*), integer L, for
!    piecewise polynomial representation of the function F to
!    be evaluated.
!
!    Input, integer K, the order of the polynomial pieces
!    that make up the function F.  The most usual value for
!    K is 4, signifying a piecewise cubic polynomial.
!
!    Input, real X, the point at which to evaluate F or
!    of its derivatives.
!
!    Input, integer JDERIV, the order of the derivative to be
!    evaluated.  If JDERIV is 0, then F itself is evaluated,
!    which is actually the most common case.  It is assumed
!    that JDERIV is zero or positive.
!
!    Output, real VALUE, the value of the JDERIV-th
!    derivative of F at X.
!
  integer k
  integer l
!
  real break(l+1)
  real coef(k,l)
  real fmmjdr
  real h
  integer i
  integer jderiv
  integer m
  integer ndummy
  real value
  real x
!
  value = 0.0
 
  fmmjdr = k - jderiv
!
!  Derivatives of order K or higher are identically zero.
!
  if ( k <= jderiv ) then
    return
  end if
!
!  Find the index I of the largest breakpoint to the left of X.
!
  call interv(break,l+1,x,i,ndummy)
!
!  Evaluate the JDERIV-th derivative of the I-th polynomial piece
!  at X.
!
  h = x-break(i)
  m = k
 
10    continue
 
  value = (value/fmmjdr)*h + coef(m,i)
  m = m-1
  fmmjdr = fmmjdr-1.0
  if ( fmmjdr > 0.0 ) then
    go to 10
  end if
 
  return
end
subroutine putit ( t, kpm, left, scrtch, dbiatx, q, nrow, b )
!
!*************************************************************************
!
!! PUTIT puts together one block of the collocation equation system
!
!
!  Method:
!
!    the  k  collocation equations for the interval  (t(left),t(left+1))
!    are constructed with the aid of the subroutine DIFEQU( 2, .,
!    . ) and interspersed (in order) with the side conditions (if any) in
!    this interval, using  d i f e q u ( 3, ., . )  for the information.
!
!    the block  q  has  kpm  columns, corresponding to the  kpm  b-
!    splines of order  kpm  which have the interval (t(left),t(left+1))
!    in their support. the block's diagonal is part of the diagonal of the
!    total system. the first equation in this block not overlapped by the
!    preceding block is therefore equation  l o w r o w , with lowrow =
!    number of side conditions in preceding intervals (or blocks).
!
!  input 
!
!  t     knot sequence, of size  left+kpm (at least)
!  kpm   order of spline
!  left  integer indicating interval of interest, viz the interval
!           (t(left), t(left+1))
!  nrow  number of rows in block to be put together
!
!  work a r e a 
!  scrtch   used in bsplvd, of size (kpm,kpm)
!  dbiatx   used to contain derivatives of b-splines, of size (kpm,m+1)
!           with dbiatx(j,i+1) containing the i-th derivative of the
!           j-th b-spline of interest
!
!  output 
!
!  q  the block, of size (nrow,kpm)
!  b  the corresponding piece of the right side, of size (nrow)
!
  integer kpm
  integer nrow
!
  real b(*)
  real dbiatx(kpm,*)
  real dx
  integer i
  integer irow
  integer iside
  integer itermx
  integer j
  integer k
  integer left
  integer ll
  integer lowrow
  integer m
  integer mode
  integer mp1
  real q(nrow,kpm)
  real rho
  real scrtch(*)
  real sum
  real t(*)
  real v(20)
  real xm
  real xside
  real xx
!
  common /side/ m,iside,xside(10)
  common /other/ itermx,k,rho(19)
!
  mp1 = m+1
  
  do j = 1, kpm
    do i = 1, nrow
      q(i,j) = 0.0
    end do
  end do
  
  xm = (t(left+1)+t(left))/2.0
  dx = (t(left+1)-t(left))/2.0

  ll = 1
  lowrow = iside

  do irow = lowrow, nrow

    if ( ll > k ) then
      go to 20
    end if

    mode = 2
!
!  next collocation point is ...
!
    xx = xm+dx*rho(ll)
    ll = ll+1
!
!  The corresp.collocation equation is next unless the next side
!  condition occurs at a point at, or to the left of, the next
!  collocation point.
!
    if ( iside > m ) then
      go to 30
    end if

    if ( xside(iside) > xx ) then
      go to 30
    end if

    ll = ll-1

   20   continue

    mode = 3
    xx = xside(iside)

   30   continue

    call difequ(mode,xx,v)
!
!  the next equation, a collocation equation (mode=2) or a side
!  condition (mode=3), reads
!    (*)   (v(m+1)*d**m+v(m)*d**(m-1) +...+ v(1)*d**0)f(xx)=v(m+2)
!  in terms of the info supplied by  difequ . the corresponding
!  equation for the b-coeffs of  f  therefore has the left side of
!  (*), evaluated at each of the  kpm  b-splines having  xx  in
!  their support, as its  kpm  possibly nonzero coefficients.
!
    call bsplvd(t,kpm,xx,left,scrtch,dbiatx,mp1)
    
    do j = 1, kpm
    
      sum = 0.0
      do i = 1, mp1
        sum = v(i)*dbiatx(j,i)+sum
      end do
      
      q(irow,j) = sum
      
    end do
   
    b(irow) = v(m+2)
    
  end do
   
  return
end
subroutine sbblok ( bloks, integs, nbloks, ipivot, b, x )
!
!*************************************************************************
!
!! SBBLOK supervises the solution (by forward and backward 
!  substitution) of the linear system  a*x=b  for x, with the 
!  plu factorization of  a
!  already generated in  f c b l o k .  individual blocks of equations
!  are solved via  s u b f o r  and  s u b b a k .
!
!  Parameters:
!
!    bloks, integs, nbloks, ipivot    are as on return from fcblok.
!
!    b       the right side, stored corresponding to the storage of
!            the equations. see comments in  s l v b l k  for details.
!
!    x       solution vector
!
  integer nbloks
!
  real b(*)
  real bloks(*)
  integer i
  integer index
  integer indexb
  integer indexx
  integer integs(3,nbloks)
  integer ipivot(*)
  integer j
  integer last
  integer nbp1
  integer ncol
  integer nrow
  real x(*)
!
!  Forward substitution pass:
!
  index = 1
  indexb = 1
  indexx = 1
  do i = 1, nbloks
    nrow = integs(1,i)
    last = integs(3,i)
    call subfor(bloks(index),ipivot(indexb),nrow,last,b(indexb),x(indexx))
    index = nrow*integs(2,i)+index
    indexb = indexb+nrow
    indexx = indexx+last
  end do
!
!  Back substitution pass.
!
  nbp1 = nbloks+1
  do j = 1, nbloks
    i = nbp1-j
    nrow = integs(1,i)
    ncol = integs(2,i)
    last = integs(3,i)
    index = index-nrow*ncol
    indexb = indexb-nrow
    indexx = indexx-last
    call subbak(bloks(index),ipivot(indexb),nrow,ncol,last,x(indexx))
  end do
   
  return
end
subroutine setupq ( x, dx, y, npoint, v, qty )
!
!*************************************************************************
!
!! SETUPQ is to be called in SMOOTH.
!
!  put  delx=x(.+1)-x(.)  into  v(.,4),
!  put  the three bands of  q-transp*d  into  v(.,1-3), and
!  put the three bands of  (d*q)-transp*(d*q)  at and above the diagonal
!  into  v(.,5-7) .
!
!  here,  q is  the tridiagonal matrix of order (npoint-2,npoint)
!  with general row  1/delx(i) , -1/delx(i)-1/delx(i+1) , 1/delx(i+1)
!  and   d  is the diagonal matrix  with general row  dx(i) .
!
  integer npoint
!
  real diff
  real dx(npoint)
  integer i
  real prev
  real qty(npoint)
  real v(npoint,7)
  real x(npoint)
  real y(npoint)
!
  v(1,4) = x(2)-x(1)
  
  do i = 2, npoint-1
    v(i,4) = x(i+1)-x(i)
    v(i,1) = dx(i-1)/v(i-1,4)
    v(i,2) = -dx(i)/v(i,4)-dx(i)/v(i-1,4)
    v(i,3) = dx(i+1)/v(i,4)
  end do
   
  v(npoint,1) = 0.0
  do i = 2, npoint-1
    v(i,5) = v(i,1)**2 + v(i,2)**2 + v(i,3)**2
  end do
   
  do i = 3, npoint-1
    v(i-1,6) = v(i-1,2)*v(i,1)+v(i-1,3)*v(i,2)
  end do
   
  v(npoint-1,6) = 0.0

  do i = 4, npoint-1
    v(i-2,7) = v(i-2,3)*v(i,1)
  end do
   
  v(npoint-2,7) = 0.0
  v(npoint-1,7) = 0.0
!
!  Construct  q-transp. * y  in  qty.
!
  prev = (y(2)-y(1))/v(1,4)
  do i = 2, npoint-1
    diff = (y(i+1)-y(i))/v(i,4)
    qty(i) = diff-prev
    prev = diff
  end do
  
  return
end
subroutine shiftb ( ai, ipivot, nrowi, ncoli, last, ai1, nrowi1, ncoli1 )
!
!*************************************************************************
!
!! SHIFTB shifts the rows in current block, ai, not used as pivot 
!  rows, if any, i.e., rows ipivot(last+1),...,ipivot(nrowi), 
!  onto the first mmax=nrow-last rows of the next block, ai1, 
!  with column last+j of ai  going to column j , 
!  for j=1,...,jmax=ncoli-last. the remaining columns of these 
!  rows of ai1 are zeroed out.
!
!                             picture
!
!       original situation after         results in a new block i+1
!       last=2 columns have been       created and ready to be
!       done in factrb (assuming no      factored by next factrb call.
!       interchanges of rows)
!                   1
!              x  x 1x  x  x           x  x  x  x  x
!                   1
!              0  x 1x  x  x           0  x  x  x  x
!  block i          1                       ---
!  nrowi=4     0  0 1x  x  x           0  0 1x  x  x  0  01
!  ncoli=5          1                       1             1
!  last=2      0  0 1x  x  x           0  0 1x  x  x  0  01
!              -------------------          1             1   new
!                   1x  x  x  x  x          1x  x  x  x  x1  block
!                   1                       1             1   i+1
!  block i+1        1x  x  x  x  x          1x  x  x  x  x1
!  nrowi1= 5        1                       1             1
!  ncoli1= 5        1x  x  x  x  x          1x  x  x  x  x1
!              -------------------          1-------------1
!                   1
!
  integer ncoli
  integer ncoli1
  integer nrowi1
  integer nrowi
!
  real ai(nrowi,ncoli)
  real ai1(nrowi1,ncoli1)
  integer ip
  integer ipivot(nrowi)
  integer j
  integer last
  integer m
!
  if ( nrowi-last < 1 ) then
    return
  end if

  if ( ncoli-last < 1 ) then
    return
  end if
!
!  Put the remainder of block I into AI1.
!
  do m = 1, nrowi-last
    ip = ipivot(last+m)
    do j = 1, ncoli-last
      ai1(m,j) = ai(ip,last+j)
    end do
  end do
!
!  Zero out the upper right corner of ai1.
!
  do j = ncoli+1-last, ncoli1
    do m = 1, nrowi-last
      ai1(m,j) = 0.0
    end do
  end do
  
  return
end
subroutine slvblk ( bloks, integs, nbloks, b, ipivot, x, iflag )
!
!*************************************************************************
!
!! SLVBLK solves the almost block diagonal linear system A*x=b.  
!
!
!  such almost block diagonal matrices
!  arise naturally in piecewise polynomial interpolation or approx-
!  imation and in finite element methods for two-point boundary value
!  problems.  the plu factorization method is implemented here to take
!  advantage of the special structure of such systems for savings in
!  computing time and storage requirements.
!
!  Parameters:
!
!  bloks  a one-dimenional array, of length
!                   sum( integs(1,i)*integs(2,i) ; i=1,nbloks )
!         on input, contains the blocks of the almost block diagonal
!         matrix  a  .  the array integs (see below and the example)
!         describes the block structure.
!         on output, contains correspondingly the plu factorization
!         of  a  (if iflag /= 0).  certain of the entries into bloks
!         are arbitrary (where the blocks overlap).
!
!  integs integer array description of the block structure of  a .
!         integs(1,i)=no. of rows of block i       = nrow
!         integs(2,i)=no. of colums of block i     = ncol
!         integs(3,i)=no. of elim. steps in block i = last
!                          i =1,2,...,nbloks
!         the linear system is of order
!         n = sum ( integs(3,i) , i=1,...,nbloks ),
!         but the total number of rows in the blocks is
!         nbrows=sum( integs(1,i) ; i = 1,...,nbloks)
!
!  nbloks number of blocks
!  b       right side of the linear system, array of length nbrows.
!          certain of the entries are arbitrary, corresponding to
!          rows of the blocks which overlap (see block structure and
!          the example below).
!  ipivot  on output, integer array containing the pivoting sequence
!          used. length is nbrows
!  x       on output, contains the computed solution (if iflag /= 0)
!          length is n.
!  iflag   on output, integer
!          =(-1)**(no. of interchanges during factorization)
!                   if  a  is invertible
!          =0    if  a  is singular
!
!                   auxiliary programs
!  fcblok (bloks,integs,nbloks,ipivot,scrtch,iflag)  factors the matrix
!           a , and is used for this purpose in slvblk. its arguments
!          are as in slvblk, except for
!              scrtch=a work array of length max(integs(1,i)).
!
!  sbblok (bloks,integs,nbloks,ipivot,b,x)  solves the system a*x=b
!          once  a  is factored. this is done automatically by slvblk
!          for one right side b, but subsequent solutions may be
!          obtained for additional b-vectors. the arguments are all
!          as in slvblk.
!
!  dtblok (bloks,integs,nbloks,ipivot,iflag,detsgn,detlog) computes the
!          determinant of  a  once slvblk or fcblok has done the fact-
!          orization.the first five arguments are as in slvblk.
!              detsgn =sign of the determinant
!              detlog =natural log of the determinant
!
!              block structure of  a  
!  the nbloks blocks are stored consecutively in the array  bloks .
!  the first block has its (1,1)-entry at bloks(1), and, if the i-th
!  block has its (1,1)-entry at bloks(index(i)), then
!         index(i+1)=index(i) + nrow(i)*ncol(i) .
!    the blocks are pieced together to give the interesting part of  a
!  as follows.  for i=1,2,...,nbloks-1, the (1,1)-entry of the next
!  block (the (i+1)st block ) corresponds to the (last+1,last+1)-entry
!  of the current i-th block.  recall last=integs(3,i) and note that
!  this means that
!      a. every block starts on the diagonal of  a .
!      b. the blocks overlap (usually). the rows of the (i+1)st block
!         which are overlapped by the i-th block may be arbitrarily de-
!         fined initially. they are overwritten during elimination.
!    the right side for the equations in the i-th block are stored cor-
!  respondingly as the last entries of a piece of  b  of length  nrow
!  (= integs(1,i)) and following immediately in  b  the corresponding
!  piece for the right side of the preceding block, with the right side
!  for the first block starting at  b(1) . in this, the right side for
!  an equation need only be specified once on input, in the first block
!  in which the equation appears.
!
!              example and test driver 
!    the test driver for this package contains an example, a linear
!  system of order 11, whose nonzero entries are indicated in the fol-
!  lowing schema by their row and column index modulo 10. next to it
!  are the contents of the  integs  arrray when the matrix is taken to
!  be almost block diagonal with  nbloks=5, and below it are the five
!  blocks.
!
!                      nrow1=3, ncol1 = 4
!           11 12 13 14
!           21 22 23 24   nrow2=3, ncol2 = 3
!           31 32 33 34
!  last1=2      43 44 45
!                 53 54 55            nrow3=3, ncol3 = 4
!        last2=3         66 67 68 69   nrow4 = 3, ncol4 = 4
!                          76 77 78 79      nrow5=4, ncol5 = 4
!                          86 87 88 89
!                 last3=1   97 98 99 90
!                    last4=1   08 09 00 01
!                                18 19 10 11
!                       last5=4
!
!         actual input to bloks shown by rows of blocks of  a .
!      (the ** items are arbitrary, this storage is used by slvblk)
!
!  11 12 13 14  / ** ** **  / 66 67 68 69  / ** ** ** **  / ** ** ** **
!  21 22 23 24 /  43 44 45 /  76 77 78 79 /  ** ** ** ** /  ** ** ** **
!  31 32 33 34/   53 54 55/   86 87 88 89/   97 98 99 90/   08 09 00 01
!                                                           18 19 10 11
!
!  index=1      index = 13  index = 22     index = 34     index = 46
!
!         actual right side values with ** for arbitrary values
!  b1 b2 b3 ** b4 b5 b6 b7 b8 ** ** b9 ** ** b10 b11
!
!  (it would have been more efficient to combine block 3 with block 4)
!
  integer nbloks
!
  real b(*)
  real bloks(*)
  integer iflag
  integer integs(3,nbloks)
  integer ipivot(*)
  real x(*)
!
!  In the call to FCBLOK, X is used for temporary storage.
!
  call fcblok ( bloks, integs, nbloks, ipivot, x, iflag )
  
  if ( iflag == 0 ) then
    return
  end if
  
  call sbblok(bloks,integs,nbloks,ipivot,b,x)
  
  return
end
subroutine smooth ( x, y, dy, npoint, s, v, a, sfp )
!
!*************************************************************************
!
!! SMOOTH constructs the cubic smoothing spline F to given data  
!
!    (x(i),y(i)), i=1,...,npoint, 
!
!  which has as small a second derivative as possible while
!
!    s(f)=sum( ((y(i)-f(x(i)))/dy(i))**2 , i=1,...,npoint ) <= s .
!
!  input 
!
!  x(1),...,x(npoint)   data abscissae,  assumed  to be strictly
!        increasing .
!  y(1),...,y(npoint)     corresponding data ordinates .
!  dy(1),...,dy(npoint)     estimate of uncertainty in data,  a s s u m-
!        e d  to be positive .
!  npoint.....number of data points,  assumed  > 1
!  s.....upper bound on the discrete weighted mean square distance of
!        the approximation  f  from the data .
!
!  work arrays:
!
!  v      of size (npoint,7)
!  a      of size (npoint,4)
!
!  output
!
!  a(.,1).....contains the sequence of smoothed ordinates .
!  a(i,j)=d**(j-1)f(x(i)), j=2,3,4, i=1,...,npoint-1 ,  i.e., the
!        first three derivatives of the smoothing spline  f  at the
!        left end of each of the data intervals .
!     Warning . . .   a  would have to be transposed before it
!        could be used in  ppvalu .
!
!  Method: 
!
!  the matrices  q-transp*d  and  q-transp*d**2*q  are constructed in
!   s e t u p q  from  x  and  dy , as is the vector  qty=q-transp*y .
!  then, for given  p , the vector  u  is determined in  c h o l 1 d  as
!  the solution of the linear system
!               (6(1-p)q-transp*d**2*q+p*r)u =qty  .
!  from  u , the smoothing spline  f  (for this choice of smoothing par-
!  ameter  p ) is obtained in the sense that
!                        f(x(.)) = y-6(1-p)d**2*q*u        and
!                  (d**2)f(x(.)) = 6*p*u                      .
!
!  the smoothing parameter  p  is found (if possible) so that
!                sf(p) = s ,
!  with  sf(p)=s(f) , where  f  is the smoothing spline as it depends
!  on  p .  if  s=0, then p = 1 . if  sf(0) <= s , then p = 0 .
!  otherwise, the secant method is used to locate an appropriate  p  in
!  the open interval  (0,1) . specifically,
!                p(0)=0,  p(1) = (s-sf(0))/dsf
!  with  dsf=-24*u-transp*r*u  a good approximation to  d(sf(0)) = dsf
!  +60*(d*q*u)-transp*(d*q*u) , and  u  as obtained for  p=0 .
!  after that, for n=1,2,...  until sf(p(n)) <= 1.01*s, do....
!  determine  p(n+1)  as the point at which the secant to  sf  at the
!  points  p(n)  and  p(n-1)  takes on the value  s .
!  if  p(n+1) >= 1 , choose instead  p(n+1)  as the point at which
!  the parabola  sf(p(n))*((1-.)/(1-p(n)))**2  takes on the value  s.
!
!  Note that, in exact arithmetic, always  p(n+1) < p(n) , hence
!  sf(p(n+1)) < sf(p(n)) . therefore, also stop the iteration,
!  with final  p=1 , in case  sf(p(n+1)) >= sf(p(n)) .
!
  integer npoint
!
  real a(npoint,4)
  real change
  real dy(npoint)
  integer i
  real p
  real prevp
  real prevsf
  real s
  real sfp
  real utru
  real v(npoint,7)
  real x(npoint)
  real y(npoint)
!
  call setupq(x,dy,y,npoint,v,a(1,4))
 
  if ( s > 0.0 ) then
    go to 20
  end if

10    continue

  p = 1.0
  call chol1d (p,v,a(1,4),npoint,a(1,3),a(1,1))
  sfp = 0.0
  go to 70

20    continue
 
  p = 0.0
  call chol1d (p,v,a(1,4),npoint,a(1,3),a(1,1))
  
  sfp = 0.0
  do i = 1, npoint
    sfp = sfp + (a(i,1)*dy(i))**2
  end do
  sfp = sfp * 36.0

  if ( sfp <= s ) then
    go to 70
  end if
  
  prevp = 0.0
  prevsf = sfp

  utru = 0.0
  do i = 2, npoint
    utru = utru+v(i-1,4)*(a(i-1,3)*(a(i-1,3)+a(i,3))+a(i,3)**2)
  end do

  p = (sfp-s) / (24.0*utru)
!
!  Secant iteration for the determination of p starts here.
!
   50 continue

  call chol1d(p,v,a(1,4),npoint,a(1,3),a(1,1))

  sfp = 0.0
  do i = 1, npoint
    sfp = sfp+(a(i,1)*dy(i))**2
  end do
  sfp = sfp*36.0*(1.-p)**2

  if ( sfp <= 1.01*s ) then
    go to 70
  end if

  if ( sfp >= prevsf ) then
    go to 10
  end if

  change = (p-prevp) / (sfp-prevsf)*(sfp-s)
  prevp = p
  p = p-change
  prevsf = sfp

  if ( p >= 1.0 ) then
    p = 1.0-sqrt(s/prevsf) * (1.0-prevp)
  end if

  go to 50
!
!  The correct value of p has been found.
!  compute pol.coefficients from  q*u (in a(.,1)).
!
   70 continue

  do i = 1, npoint
    a(i,1) = y(i)-6*(1-p)*dy(i)**2 * a(i,1)
  end do

  do i = 1, npoint
    a(i,3) = 6.0*a(i,3)*p
  end do

  do i = 1, npoint-1
    a(i,4) = (a(i+1,3)-a(i,3))/v(i,4)
    a(i,2) = (a(i+1,1)-a(i,1))/v(i,4)-(a(i,3)+a(i,4)/3.*v(i,4))/2.*v(i,4)
  end do
 
  return
end
subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )
!
!*************************************************************************
!
!! SPLI2D is an extended version of SPLINT, for use in tensor 
!  product spline interpolation.
! 
!  SPLI2D produces the B-spline coefficients BCOEF(J,.) of the 
!  spline of order K with knots T(I), I=1,..., N+K, which takes on 
!  the value GTAU(I,J) at TAU(I), I=1,..., N, J=1,...,M.
! 
!  The I-th equation of the linear system 
!
!    A*BCOEF=B  
!  
!  for the B-spline coefficients of the interpolant enforces 
!  interpolation at TAU(I), I=1,...,N.  Hence,  B(I)=GTAU(I), 
!  all I, and A is a band matrix with 2K-1 bands, if it is 
!  invertible.
! 
!  The matrix A is generated row by row and stored, diagonal by
!  diagonal, in the rows of the array Q, with the main diagonal
!  going into row K.
! 
!  The banded system is then solved by a call to BANFAC, which 
!  constructs the triangular factorization for A and stores it 
!  again in Q, followed by a call to BANSLV, which then obtains 
!  the solution BCOEF by substitution.
! 
!  Parameters:
!
!  TAU    Input, real TAU(N), contains the data point abscissas.
!         TAU must be strictly increasing
! 
!  GTAU   Input, real GTAU(N), contains the data point ordinates, 
!         J=1,...,M.
! 
!  T      Input, real T(N+K), the knot sequence.
! 
!  N      Input, integer N,  the number of data points and the 
!         dimension of the spline space SPLINE(K,T)
! 
!  K      Input, integer K, the order of the spline.
! 
!  M      Input, integer M, the number of data sets.
!
!  WORK   Work space, real WORK(N).
! 
!  Q      Output, real Q(2*K-1)*N, containing the triangular 
!         factorization of the coefficient matrix of the linear 
!         system for the B-spline coefficients of the spline interpolant.
!       
!         The B-spline coefficients for the interpolant of an additional 
!         data set (TAU(I),HTAU(I)), I=1,...,N  with the same data 
!         abscissae can be obtained without going through all the 
!         calculations in this routine, simply by loading HTAU into 
!         BCOEF and then using the statement
!       
!           CALL BANSLV(Q,2*K-1,N,K-1,K-1,BCOEF)
! 
!  BCOEF  Output, real BCOEF(N), the B-spline coefficients of 
!         the interpolant.
! 
!  IFLAG  Output, integer IFLAG, error indicator.
!
!         1, no error.
!       
!         2, an error occurred, which may have been caused by 
!         singularity of the linear system.
!
!         The linear system to be solved is theoretically invertible if
!         and only if 
!       
!           T(I) < TAU(I) < TAU(I+K), for all I.
!         
!        Violation of this condition is certain to lead to IFLAG=2.
!
  integer m
  integer n
!
  real bcoef(m,n)
  real gtau(n,m)
  integer i
  integer iflag
  integer ilp1mx
  integer j
  integer jj
  integer k
  integer left
  real q((2*k-1)*n)
  real t(n+k)
  real tau(n)
  real taui
  real work(n)
!
  left = k
  
  do i = 1, (2*k-1)*n
    q(i) = 0.0
  end do
!
!  Construct the N interpolation equations.
!
  do i = 1, n
  
    taui = tau(i)
    ilp1mx = min(i+k,n+1)
!
!  Find the index LEFT in the closed interval (I,I+K-1) such 
!  that:
!
!    T(LEFT) < = TAU(I) < T(LEFT+1)
!
!  The matrix will be singular if this is not possible.
!
    left = max(left,i)
    
    if ( taui < t(left) ) then
      iflag = 2
      write(*,*)' '
      write(*,*)'SPLI2D - Fatal error!'
      write(*,*)'  The TAU array is not strictly increasing.'
      stop
    end if
    
   20   continue
   
    if ( taui >= t(left+1) ) then
   
      left = left+1
      if ( left < ilp1mx ) then
        go to 20
      end if
    
      left = left-1
    
      if ( taui > t(left+1) ) then
        iflag = 2
        write(*,*)' '
        write(*,*)'SPLI2D - Fatal error!'
        write(*,*)'  The TAU array is not strictly increasing.'
        stop
      end if
 
    end if
!
!  The I-th equation enforces interpolation at TAUI, hence
!
!    A(I,J)=B(J,K,T)(TAUI), for all J. 
!
!  Only the K entries with J=LEFT-K+1, ..., LEFT actually might be 
!  nonzero.  These K numbers are returned, in WORK (used for 
!  temporary storage here), by the following call:
!
    call bsplvb(t,k,1,taui,left,work)
!
!  We therefore want  
!
!    WORK(J)=B(LEFT-K+J)(TAUI) 
!
!  to go into
!
!        a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
!        a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
!        as a two-dim. array , with  2*k-1  rows (see comments in
!        banfac). in the present program, we treat  q  as an equivalent
!        one-dimensional array (because of fortran restrictions on
!  ??  LOST LINE ??
!        entry
!            i -(left+j)+2*k + ((left+j)-k-1)*(2*k-1)
!                 = i-left+1+(left -k)*(2*k-1) + (2*k-2)*j
!  of  q .
!
    jj = i-left+1+(left-k)*(k+k-1)
    
    do j = 1, k
      jj = jj+k+k-2
      q(jj) = work(j)
    end do
    
  end do
!
!  Factor A, stored again in Q.
!
  call banfac(q,k+k-1,n,k-1,k-1,iflag)
  
  if ( iflag == 2 ) then
    write(*,*)' '
    write(*,*)'SPLI2D - Fatal error!'
    write(*,*)'  BANFAC reports that the matrix is singular.'
    stop
  end if
!
!  Solve A*BCOEF=GTAU by backsubstitution.
!
  do j = 1, m
   
    do i = 1, n
      work(i) = gtau(i,j)
    end do
   
    call banslv(q,k+k-1,n,k-1,k-1,work)
    
    do i = 1, n
      bcoef(j,i) = work(i)
    end do
    
  end do
   
  return
end
subroutine splint ( tau, gtau, t, n, k, q, bcoef, iflag )
!
!*************************************************************************
!
!! SPLINT produces the B-spline coefficients BCOEF of the spline of 
!  order k  with knots  t(i), i=1,..., n+k , which takes on the 
!  value gtau(i) at  tau(i), i=1,..., n .
!
!  Method:
!
!    the i-th equation of the linear system  a*bcoef=b  for the b-co-
!    effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
!    hence,  b(i)=gtau(i), all i, and  a  is a band matrix with  2k-1
!    bands (if it is invertible).
!
!    the matrix  a  is generated row by row and stored, diagonal by di-
!    agonal, in the  r o w s  of the array  q , with the main diagonal go-
!    ing into row  k .  see comments in the program below.
!    the banded system is then solved by a call to  banfac (which con-
!    structs the triangular factorization for  a  and stores it again in
!    q ), followed by a call to  banslv (which then obtains the solution
!    bcoef  by substitution).
!
!    banfac  does no pivoting, since the total positivity of the matrix
!    a  makes this unnecessary.
!
!  input 
!
!  tau    array of length  n , containing data point abscissae.
!         a s s u m p t i o n . . .  tau  is strictly increasing
!
!  gtau   corresponding array of length  n , containing data point or-
!         dinates
!
!  t      knot sequence, of length  n+k
!
!  n      number of data points and dimension of spline space  s(k,t)
!
!  k      order of spline
!
!  output 
!
!  q, array of size  (2*k-1)*n , containing the triangular factoriz-
!        ation of the coefficient matrix of the linear system for the b-
!        coefficients of the spline interpolant.
!        the b-coeffs for the interpolant of an additional data set can
!        be obtained without going through all the calculations in this
!        routine, simply by loading  htau  into  bcoef  and then execut-
!        ing the    call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
!
!  bcoef, the b-coefficients of the interpolant, of length  n.
!
!  iflag, an integer indicating success (= 1)  or failure (= 2)
!        the linear system to be solved is (theoretically) invertible if
!        and only if
!              t(i) < tau(i) < tau(i+k),    all i.
!        violation of this condition is certain to lead to  iflag=2 .
!
  integer n
!
  real bcoef(n)
  real gtau(n)
  integer i
  integer iflag
  integer ilp1mx
  integer j
  integer jj
  integer k
  integer kpkm2
  integer left
  real q((2*k-1)*n)
  real t(n+k)
  real tau(n)
  real taui
!
  kpkm2 = 2*(k-1)
  left = k
!
  do i = 1, (2*k-1)*n
    q(i) = 0.0
  end do
!
!  loop over i to construct the  n  interpolation equations
!
  do i = 1, n
  
    taui = tau(i)
    ilp1mx = min(i+k,n+1)
!
!  find  left  in the closed interval (i,i+k-1) such that
!    t(left) <= tau(i)  < t(left+1)
!  matrix is singular if this is not possible
!
    left = max(left,i)

    if ( taui < t(left)) then
      go to 70
    end if

   20   continue

    if ( taui < t(left+1)) then
      go to 30
    end if

    left = left+1
    if ( left < ilp1mx) then
      go to 20
    end if

    left = left-1
    if ( taui > t(left+1)) then
      go to 70
    end if
!
!  the i-th equation enforces interpolation at taui, hence
!  a(i,j)=b(j,k,t)(taui), all j. only the  k  entries with  j =
!  left-k+1,...,left actually might be nonzero. these  k  numbers
!  are returned, in  bcoef (used for temp.storage here), by the
!  following
!
   30   continue

    call bsplvb(t,k,1,taui,left,bcoef)
!
!  we therefore want  bcoef(j)=b(left-k+j)(taui) to go into
!  a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
!  a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
!  as a two-dim. array , with  2*k-1  rows (see comments in
!  banfac). in the present program, we treat  q  as an equivalent
!  one-dimensional array (because of fortran restrictions on
!  dimension statements) . we therefore want  bcoef(j) to go into
!  entry
!    i -(left+j)+2*k + ((left+j)-k-1)*(2*k-1)
!   = i-left+1+(left -k)*(2*k-1) + (2*k-2)*j
!  of  q .
!
    jj = i-left+1+(left-k)*(k+k-1)
    do j = 1,k
      jj = jj+kpkm2
      q(jj) = bcoef(j)
    end do
    
  end do
!
!  Obtain factorization of A, stored again in Q.
!
  call banfac(q,k+k-1,n,k-1,k-1,iflag)
  
  if ( iflag == 2 ) then
    write(*,*)' '
    write(*,*)'SPLINT - Fatal Error!'
    write(*,*)'  The linear system is not invertible!'
    return
  end if
!
!  Solve  a*bcoef=gtau  by backsubstitution
!
  do i = 1, n
    bcoef(i) = gtau(i)
  end do

  call banslv(q,k+k-1,n,k-1,k-1,bcoef)
  return
  
   70 iflag=2
 
  write(*,*)' '
  write(*,*)'SPLINT - Fatal Error!'
  write(*,*)'  The linear system is not invertible!'
 
  return
end
subroutine splopt ( tau, n, k, scrtch, t, iflag )
!
!*************************************************************************
!
!! SPLOPT computes the knots T for the optimal recovery scheme 
!  of order  k for data at  tau(i), i=1,...,n .
!
!  method 
!     the (interior) knots  t(k+1), ..., t(n)  are determined by newtons
!  method in such a way that the signum function which changes sign at
!   t(k+1), ..., t(n)  and nowhere else in  (tau(1),tau(n)) is orthogon-
!  al to the spline space  spline( k , tau )  on that interval .
!     let  xi(j)  be the current guess for  t(k+j), j=1,...,n-k. then
!  the next newton iterate is of the form
!              xi(j) + (-)**(n-k-j)*x(j)  ,  j=1,...,n-k,
!  with  x  the solution of the linear system
!                        c*x = d  .
!  here,  c(i,j)=b(i)(xi(j)), all j, with  b(i)  the i-th b-spline of
!  order  k  for the knot sequence  tau , all i, and  d  is the vector
!  given by  d(i)=sum( -a(j) , j=i,...,n )*(tau(i+k)-tau(i))/k, all i,
!  with  a(i)=sum ( (-)**(n-k-j)*b(i,k+1,tau)(xi(j)) , j=1,...,n-k )
!  for i=1,...,n-1, and  a(n)=-.5 .
!     (see chapter  xiii  of text and references there for a derivation)
!     the first guess for  t(k+j)  is  (tau(j+1)+...+tau(j+k-1))/(k-1) .
!     iteration terminates if  max(abs(x(j))) < t o l  , with
!                 t o l = t o l r t e *(tau(n)-tau(1))/(n-k) ,
!  or else after  n e w t m x  iterations , currently,
!                 newtmx, tolrte / 10, .000001
!
!  input 
!
!  tau    array of length  n , containing the interpolation points.
!         assumed to be nondecreasing, with tau(i) < tau(i+k),all i.
!  n      number of data points .
!  k      order of the optimal recovery scheme to be used .
!
!  work array
!
!  scrtch, array of length  (n-k)(2k+3)+5k + 3 . the various
!        contents are specified in the text below .
!
!  output 
!  iflag, integer indicating success (=1) or failure (=2) .
!     if iflag=1, then
!  t, array of length  n+k  containing the optimal knots ready for
!        use in optimal recovery. specifically,  t(1)=... = t(k) =
!        tau(1)  and  t(n+1)=... = t(n+k) = tau(n) , while the  n-k
!        interior knots  t(k+1), ..., t(n)  are calculated as described
!        below under  *method* .
!     if iflag=2, then
!        k < 3, or n < k, or a certain linear system was found to
!        be singular.
!
!  p r i n t e d  output 
!  a comment will be printed in case  ilfag=2  or newton iterations
!  failed to converge in  n e w t m x  iterations .
!
  integer k
  integer n
!
  real del
  real delmax
  real floatk
  integer i
  integer id
  integer iflag
  integer index
  integer j
  integer kp1
  integer kpkm1
  integer kpn
  integer l
  integer left
  integer leftmk
  integer lenw
  integer ll
  integer llmax
  integer llmin
  integer na
  integer nb
  integer nc
  integer nd
  integer newtmx
  integer newton
  integer nmk
  integer nx
  real scrtch((n-k)*(2*k+3)+5*k+3)
  real t(n+k)
  real tau(n)
  real sign
  real signst
  real sum
  real tol
  real tolrte
  real xij
!
  data newtmx / 10 /
  data tolrte / 0.000001 /
!
  nmk = n-k
  
  if ( n < k ) then
    write(*,*)' '
    write(*,*)'SPLOPT - Fatal error!'
    write(*,*)'  N < K, N = ',n,' K = ',k
    iflag=2
    return
  end if
  
  if ( n == k ) then
    do i = 1, k
      t(i) = tau(1)
      t(n+i) = tau(n)
    end do
    return
  end if
  
  if ( k <= 2 ) then
    write(*,*)' '
    write(*,*)'SPLOPT - Fatal error!'
    write(*,*)'  K < 2, K = ',k
    iflag = 2
    stop
  end if
 
  floatk = k
  kp1 = k+1
  kpkm1 = k+k-1
  kpn = k+n

  signst = -1.0
  if ( nmk > (nmk/2) * 2 ) then
    signst = 1.0
  end if
!
!  scrtch(i)=tau-extended(i), i=1,...,n+k+k
!
  nx = n+k+k+1
!
!  scrtch(i+nx)=xi(i),i=0,...,n-k+1
!
  na = nx+nmk+1
!
!  scrtch(i+na)=-a(i), i=1,...,n
!
  nd = na+n
!
!  scrtch(i+nd)=x(i) or d(i), i=1,...,n-k
!
  nb = nd+nmk
!
!  scrtch(i+nb)=biatx(i),i=1,...,k+1
!
  nc = nb+kp1
!
!  scrtch(i+(j-1)*(2k-1)+nc)=w(i,j) = c(i-k+j,j), i=j-k,...,j+k,
!                                                     j=1,...,n-k.
!
  lenw = kpkm1*nmk
!
!  Extend TAU to a knot sequence and store in scrtch.
!
  do j = 1, k
    scrtch(j) = tau(1)
    scrtch(kpn+j) = tau(n)
  end do
 
  do j = 1, n
    scrtch(k+j) = tau(j)
  end do
!
!  first guess for  scrtch (.+nx) = xi .
!
  scrtch(nx) = tau(1)
  scrtch(nmk+1+nx) = tau(n)
 
  do j = 1, nmk
 
    sum = 0.0
    do l = 1, k-1
      sum = sum+tau(j+l)
    end do
 
    scrtch(j+nx) = sum/real ( k-1)
 
  end do
!
!  last entry of  scrtch (.+na) =-a  is always ...
!
  scrtch(n+na) = 0.5
!
!  Start newton iteration.
!
  newton = 1
  tol = tolrte*(tau(n)-tau(1))/real ( nmk)
!
!  start newton step
!  compute the 2k-1 bands of the matrix c and store in scrtch(.+nc),
!  and compute the vector  scrtch(.+na)=-a.
!
  100 continue

  do i = 1, lenw
    scrtch(i+nc) = 0.0
  end do
  
  do i = 2, n
    scrtch(i-1+na) = 0.0
  end do
  
  sign = signst
  left = kp1
  
  do j = 1, nmk
  
    xij = scrtch(j+nx)

  130   continue

    if ( xij < scrtch(left+1) ) then
      go to 140
    end if

    left = left+1
    if ( left < kpn ) then
      go to 130
    end if
    left = left-1

  140   continue

    call bsplvb(scrtch,k,1,xij,left,scrtch(1+nb))
!
!  The tau sequence in scrtch is preceded by  k  additional knots
!  therefore,  scrtch(ll+nb)  now contains  b(left-2k+ll)(xij)
!  which is destined for  c(left-2k+ll,j), and therefore for
!    w(left-k-j+ll,j)= scrtch(left-k-j+ll+(j-1)*kpkm1 + nc)
!  since we store the 2k-1 bands of  c  in the 2k-1  r o w s  of
!  the work array w, and  w  in turn is stored in  s c r t c h ,
!  with  w(1,1)=scrtch(1+nc).
!
!  also, c  being of order  n-k, we would want  
!    1 <= left-2k+ll .le. n-k  or
!    llmin=2k-left  <=  ll  .le.  n-left+k = llmax .
!
    leftmk = left-k
    index = leftmk-j+(j-1)*kpkm1+nc
    llmin = max(1,k-leftmk)
    llmax = min(k,n-leftmk)
    do ll = llmin, llmax
      scrtch(ll+index)=scrtch(ll+nb)
    end do
    
    call bsplvb (scrtch,kp1,2,xij,left,scrtch(1+nb))
    id=max(0,leftmk-kp1)
    llmin=1-min(0,leftmk-kp1)
    do ll=llmin, kp1
      id=id+1
      scrtch(id+na)=scrtch(id+na)-sign*scrtch(ll+nb)
    end do
    
    sign=-sign
    
  end do
  
  call banfac(scrtch(1+nc),kpkm1,nmk,k-1,k-1,iflag)
  
  if ( iflag == 2 ) then
    write(*,*)' '
    write(*,*)'SPLOPT - Fatal error!'
    write(*,*)'  Matrix C is not invertible.'
    stop
  end if
!
!  compute  scrtch (.+nd)= d  from  scrtch (.+na) =-a .
!
  do i=n,2,-1
    scrtch(i-1+na)=scrtch(i-1+na)+scrtch(i+na)
  end do
  
  do i=1,nmk
    scrtch(i+nd)=scrtch(i+na)*(tau(i+k)-tau(i))/floatk
  end do
!
!  Compute  scrtch (.+nd)= x .
!
  call banslv(scrtch(1+nc),kpkm1,nmk,k-1,k-1,scrtch(1+nd))
!
!  Compute  scrtch (.+nd)=change in  xi . modify, if necessary, to
!  prevent new  xi  from moving more than 1/3 of the way to its
!  neighbors. then add to  xi  to obtain new  xi  in scrtch(.+nx).
!
  delmax=0.0
  sign=signst
  do i = 1, nmk
    del = sign*scrtch(i+nd)
    delmax = max (delmax,abs(del))
    if ( del > 0.0 ) then
      go to 230
    end if
    del = max (del,(scrtch(i-1+nx)-scrtch(i+nx))/3.)
    go to 240
  230   del = min (del,(scrtch(i+1+nx)-scrtch(i+nx))/3.)
  240   sign = -sign
    scrtch(i+nx) = scrtch(i+nx)+del
  end do
!
!  Call it a day in case change in  xi  was small enough or too many
!  steps were taken.
!
  if ( delmax < tol ) then
    go to 270
  end if

  newton = newton+1
  if ( newton <= newtmx) then
    go to 100
  end if

  write(*,260)newtmx
  260 format (' no convergence in  splopt  after',i3,' newton steps.')
  
  270 do i=1, nmk
    t(k+i)=scrtch(i+nx)
  end do
  
  290 continue

  do i=1,k
    t(i)=tau(1)
    t(n+i)=tau(n)
  end do
  
  return
  
! 310 iflag=2
!
!     return
end
subroutine subbak ( w, ipivot, nrow, ncol, last, x )
!
!*************************************************************************
!
!! SUBBAK carries out backsubstitution for the current block.
!
!
!  Parameters:
!
!    w, ipivot, nrow, ncol, last  are as on return from factrb.
!
!    x(1),...,x(ncol)  contains, on input, the right side for the
!            equations in this block after backsubstitution has been
!            carried up to but not including equation ipivot(last).
!            means that x(j) contains the right side of equation ipi-
!            vot(j) as modified during elimination, j=1,...,last, while
!            for j > last, x(j) is already a component of the solut-
!            ion vector.
!
!    x(1),...,x(ncol) contains, on output, the components of the solut-
!            ion corresponding to the present block.
!
  integer ncol
  integer nrow
!
  integer ip
  integer ipivot(nrow)
  integer j
  integer k
  integer last
  real sum
  real w(nrow,ncol)
  real x(ncol)
!
  do k = last, 1, -1

    ip = ipivot(k)

    sum = 0.0
    do j = k+1, ncol
      sum = w(ip,j) * x(j) + sum
    end do
   
    x(k) = ( x(k) - sum ) / w(ip,k)

  end do

end
subroutine subfor ( w, ipivot, nrow, last, b, x )
!
!*************************************************************************
!
!! SUBFOR carries out the forward pass of substitution for the
!  current block,
!  i.e., the action on the right side corresponding to the elimination
!  carried out in  f a c t r b  for this block.
!     at the end, x(j) contains the right side of the transformed
!  ipivot(j)-th equation in this block, j=1,...,nrow. then, since
!  for i=1,...,nrow-last, b(nrow+i) is going to be used as the right
!  side of equation  i  in the next block (shifted over there from
!  this block during factorization), it is set equal to x(last+i) here.
!
!  Parameters:
!
!    w, ipivot, nrow, last  are as on return from factrb.
!
!    b(j)   is expected to contain, on input, the right side of j-th
!           equation for this block, j=1,...,nrow.
!    b(nrow+j)   contains, on output, the appropriately modified right
!           side for equation j in next block, j=1,...,nrow-last.
!
!    x(j)   contains, on output, the appropriately modified right
!           side of equation ipivot(j) in this block, j=1,...,last (and
!           even for j=last+1,...,nrow).
!
  integer last
  integer nrow
!
  real b(nrow+nrow-last)
  integer ip
  integer ipivot(nrow)
  integer j
  integer k
  real sum
  real w(nrow,last)
  real x(nrow)
!
  ip = ipivot(1)
  x(1) = b(ip)
  
  do k = 2, nrow
  
    ip = ipivot(k)
    
    sum = 0.0
    do j = 1, min(k-1,last)
      sum = w(ip,j)*x(j)+sum
    end do
    
    x(k) = b(ip)-sum

  end do
!
!  Transfer modified right sides of equations ipivot(last+1),...,
!  ipivot(nrow) to next block.
!
  do k = last+1, nrow
    b(nrow-last+k) = x(k)
  end do
  
  return
end
subroutine tautsp ( tau, gtau, ntau, gamma, s, break, coef, l, k, iflag )
!
!*************************************************************************
!
!! TAUTSP constructs a cubic spline interpolant to given data.
!
!
!  Discussion:
!
!    if  gamma > 0., additional knots are introduced where needed to
!    make the interpolant more flexible locally. this avoids extraneous
!    inflection points typical of cubic spline interpolation at knots to
!    rapidly changing data.
!
!  Method:  
!
!  on the i-th interval, (tau(i), tau(i+1)), the interpolant is of the
!  form
!  (*)  f(u(x))=a+b*u + c*h(u,z) + d*h(1-u,1-z) ,
!  with  u=u(x) = (x-tau(i))/dtau(i). here,
!       z=z(i) = addg(i+1)/(addg(i)+addg(i+1))
!  (= .5, in case the denominator vanishes). with
!       addg(j)=abs(ddg(j)), ddg(j) = dg(j+1)-dg(j),
!       dg(j)=divdif ( j) = (gtau(j+1)-gtau(j))/dtau(j)
!  and
!       h(u,z)=alpha*u**3+(1-alpha)*(max(((u-zeta)/(1-zeta)),0)**3
!  with
!       alpha(z)=(1-gamma/3)/zeta
!       zeta(z)=1-gamma*min((1 - z), 1/3)
!  thus, for 1/3 <= z .le. 2/3,  f  is just a cubic polynomial on
!  the interval i. otherwise, it has one additional knot, at
!         tau(i)+zeta*dtau(i) .
!  as  z  approaches  1, h(.,z) has an increasingly sharp bend  near 1,
!  thus allowing  f  to turn rapidly near the additional knot.
!     in terms of f(j)=gtau(j) and
!       fsecnd(j)=2.derivative of  f  at  tau(j),
!  the coefficients for (*) are given as
!       a=f(i)-d
!       b=(f(i+1)-f(i)) - (c - d)
!       c=fsecnd(i+1)*dtau(i)**2/hsecnd(1,z)
!       d=fsecnd(i)*dtau(i)**2/hsecnd(1,1-z)
!  hence can be computed once fsecnd(i),i=1,...,ntau, is fixed.
!
!  F is automatically continuous and has a continuous second derivat-
!  ive (except when z=0 or 1 for some i). we determine fscnd(.) from
!  the requirement that also the first derivative of  f  be contin-
!  uous. in addition, we require that the third derivative be continuous
!  across  tau(2) and across  tau(ntau-1) . this leads to a strictly
!  diagonally dominant tridiagonal linear system for the fsecnd(i)
!  which we solve by gauss elimination without pivoting.
!
!  Parameters:
!
!  input
!
!  tau    sequence of data points. must be strictly increasing.
!
!  gtau   corresponding sequence of function values.
!
!  ntau   number of data points. must be at least  4 .
!
!  gamma  indicates whether additional flexibility is desired.
!        =0., no additional knots
!         in (0.,3.), under certain conditions on the given data at
!                points i-1, i, i+1, and i+2, a knot is added in the
!                i-th interval, i=2,...,ntau-2. see description of meth-
!                od below. the interpolant gets rounded with increasing
!                gamma. a value of  2.5  for gamma is typical.
!          in (3.,6.), same , except that knots might also be added in
!                intervals in which an inflection point would be permit-
!                ted.  a value of  5.5  for gamma is typical.
!            output
!  break, coef, l, k  give the pp-representation of the interpolant.
!          specifically, for break(i) <= x .le. break(i+1), the
!        interpolant has the form
!  f(x)=coef(1,i) +dx(coef(2,i) +(dx/2)(coef(3,i) +(dx/3)coef(4,i)))
!        with  dx=x-break(i) and i=1,...,l .
!  iflag =1, ok
!        =2, input was incorrect. a printout specifying the mistake
!            was made.
!            workspace
!  s     is required, of size (ntau,6). the individual columns of this
!        array contain the following quantities mentioned in the write-
!        up and below.
!     s(.,1)=dtau = tau(.+1)-tau
!     s(.,2)=diag = diagonal in linear system
!     s(.,3)=u = upper diagonal in linear system
!     s(.,4)=r = right side for linear system (initially)
!          = fsecnd = solution of linear system , namely the second
!                       derivatives of interpolant at  tau
!     s(.,5)=z = indicator of additional knots
!     s(.,6)=1/hsecnd(1,x) with x = z or = 1-z. see below.
!
  integer ntau
!
  real alph
  real alpha
  real break(*)
  real c
  real coef(4,*)
  real d
  real del
  real denom
  real divdif
  real entry
  real entry3
  real factor
  real factr2
  real gam
  real gamma
  real gtau(ntau)
  integer i
  integer iflag
  integer k
  integer l
  integer method
  real onemg3
  real onemzt
  real ratio
  real s(ntau,6)
  real sixth
  real tau(ntau)
  real temp
  real x
  real z
  real zeta
  real zt2
!
  alph(x) = min (1.0,onemg3/x)
!
!  There must be at least 4 interpolation points.
!
  if ( ntau < 4 ) then
    write(*,*)' '
    write(*,*)'TAUTSP - Fatal error!'
    write(*,*)'  NTAU < 4, NTAU=',ntau
    iflag = 2
    return
  end if
!
!  Construct delta tau and first and second (divided) differences of data
!
  do i = 1, ntau-1

    s(i,1) = tau(i+1)-tau(i)
    
    if ( s(i,1) <= 0.0 ) then
      write(*,30)i,tau(i),tau(i+1)
   30     format (' point ',i3,' and the next',2e15.6,' are disordered')
      iflag=2
      return
    end if
    
    s(i+1,4)=(gtau(i+1)-gtau(i))/s(i,1)
  end do
   
  do i = 2, ntau-1
    s(i,4) = s(i+1,4)-s(i,4)
  end do
!
!  construct system of equations for second derivatives at  tau. at each
!  interior data point, there is one continuity equation, at the first
!  and the last interior data point there is an additional one for a
!  total of NTAU equations in  ntau  unknowns.
!
  i = 2
  s(2,2) = s(1,1)/3.0
  sixth = 1.0/6.0
  method = 2
  gam = gamma

  if ( gam <= 0.0 ) then
    method=1
  end if

  if ( gam > 3.0 ) then
    method = 3
    gam = gam-3.0
  end if

  onemg3=1.-gam/3.0
!
!  loop over i
!
   70 continue
!
!  Construct z(i) and zeta(i)
!
  z=0.5

  if ( method == 1) then
    go to 100
  end if

  if ( method == 3) then
    go to 90
  end if

  if ( s(i,4)*s(i+1,4) < 0.0 ) then
    go to 100
  end if

   90 temp=abs(s(i+1,4))
  denom=abs(s(i,4))+temp
  
  if ( denom /= 0.0 ) then
    z = temp/denom
    if (abs(z-0.5) <= sixth) then
      z=0.5
    end if
  end if
  
  100 s(i,5) = z
!
!  Set up part of the i-th equation which depends on the i-th interval.
!
  if ( z < 0.5 ) then

    zeta = gam*z
    onemzt = 1.-zeta
    zt2 = zeta**2
    alpha = alph(onemzt)
    factor = zeta/(alpha*(zt2-1.)+1.)
    s(i,6) = zeta*factor/6.
    s(i,2) = s(i,2)+s(i,1)*((1.-alpha*onemzt)*factor/2.-s(i,6))
!
!  if z=0 and the previous z = 1, then d(i) = 0. since then
!  also u(i-1)=l(i+1) = 0, its value does not matter. reset
!  d(i)=1 to insure nonzero pivot in elimination.
!
    if ( s(i,2) <= 0.0 ) then
      s(i,2)=1.0
    end if

    s(i,3)=s(i,1)/6.0
  else if ( z - 0.5 == 0.0 ) then
    s(i,2)=s(i,2)+s(i,1)/3.
    s(i,3)=s(i,1)/6.
  else if ( z - 0.5 > 0.0 ) then
    onemzt = gam*(1.-z)
    zeta = 1.-onemzt
    alpha = alph(zeta)
    factor = onemzt/(1.-alpha*zeta*(1.+onemzt))
    s(i,6) = onemzt*factor/6.
    s(i,2) = s(i,2)+s(i,1)/3.
    s(i,3) = s(i,6)*s(i,1)
  end if

  if ( i > 2 ) then
    go to 190
  end if

  s(1,5) = 0.5
!
!  the first two equations enforce continuity of the first and of
!  the third derivative across tau(2).
!
  s(1,2)=s(1,1)/6.
  s(1,3)=s(2,2)
  entry3=s(2,3)
  if ( z-.5) 150, 160, 170
  150 factr2=zeta*(alpha*(zt2-1.)+1.)/(alpha*(zeta*zt2-1.)+1.)
  ratio=factr2*s(2,1)/s(1,2)
  s(2,2)=factr2*s(2,1)+s(1,1)
  s(2,3)=-factr2*s(1,1)
  go to 180
  
  160 ratio=s(2,1)/s(1,2)
  s(2,2)=s(2,1)+s(1,1)
  s(2,3)=-s(1,1)
  go to 180
  
  170 ratio=s(2,1)/s(1,2)
  s(2,2)=s(2,1)+s(1,1)
  s(2,3)=-s(1,1)*6.*alpha*s(2,6)
!
!  At this point, the first two equations read
!              diag(1)*x1+u(1)*x2 + entry3*x3=r(2)
!       -ratio*diag(1)*x1+diag(2)*x2 + u(2)*x3=0.
!  Eliminate first unknown from second equation
!
  180 s(2,2)=ratio*s(1,3)+s(2,2)
  s(2,3)=ratio*entry3+s(2,3)
  s(1,4)=s(2,4)
  s(2,4)=ratio*s(1,4)
  go to 200
  
  190 continue
!
!  the i-th equation enforces continuity of the first derivative
!  across tau(i). it has been set up in statements 35 up to 40
!  and 21 up to 25 and reads now
!    -ratio*diag(i-1)*xi-1+diag(i)*xi + u(i)*xi+1=r(i) .
!  eliminate (i-1)st unknown from this equation
!
  s(i,2)=ratio*s(i-1,3)+s(i,2)
  s(i,4)=ratio*s(i-1,4)+s(i,4)
!
!  set up the part of the next equation which depends on the
!  i-th interval.
!
  200 if ( z-.5) 210, 220, 230
  210 ratio=-s(i,6)*s(i,1)/s(i,2)
  s(i+1,2)=s(i,1)/3.
  go to 240
  220 ratio=-(s(i,1)/6.)/s(i,2)
  s(i+1,2)=s(i,1)/3.
  go to 240
  230 ratio=-(s(i,1)/6.)/s(i,2)
  s(i+1,2)=s(i,1)*((1.-zeta*alpha)*factor/2.-s(i,6))
!
!  end of i loop
!
  240 i=i+1
  if ( i < ntau-1) then
    go to 70
  end if

  s(i,5)=.5
!
!  the last two equations enforce continuity of third derivative and
!  of first derivative across  tau(ntau-1).
!
  entry=ratio*s(i-1,3)+s(i,2)+s(i,1)/3.
  s(i+1,2)=s(i,1)/6.
  s(i+1,4)=ratio*s(i-1,4)+s(i,4)
  if ( z-.5) 250, 260, 270
  250 ratio=s(i,1)*6.*s(i-1,6)*alpha/s(i-1,2)
  s(i,2)=ratio*s(i-1,3)+s(i,1)+s(i-1,1)
  s(i,3)=-s(i-1,1)
  go to 280
  
  260 ratio=s(i,1)/s(i-1,2)
  s(i,2)=ratio*s(i-1,3)+s(i,1)+s(i-1,1)
  s(i,3)=-s(i-1,1)
  go to 280
  
  270 continue

  factr2=onemzt*(alpha*(onemzt**2-1.)+1.)/(alpha*(onemzt**3-1.)+1.)
  ratio=factr2*s(i,1)/s(i-1,2)
  s(i,2)=ratio*s(i-1,3)+factr2*s(i-1,1)+s(i,1)
  s(i,3)=-factr2*s(i-1,1)
!
!  At this point, the last two equations read:
!
!           diag(i)*xi+     u(i)*xi+1=r(i)
!    -ratio*diag(i)*xi+diag(i+1)*xi+1=r(i+1)
!
!  Eliminate XI from the last equation.
!
  280 s(i,4)=ratio*s(i-1,4)
  ratio=-entry/s(i,2)
  s(i+1,2)=ratio*s(i,3)+s(i+1,2)
  s(i+1,4)=ratio*s(i,4)+s(i+1,4)
!
!  Back substitution.
!
  s(ntau,4)=s(ntau,4)/s(ntau,2)

  290 continue

  s(i,4)=(s(i,4)-s(i,3)*s(i+1,4))/s(i,2)
  i=i-1
  if ( i > 1) then
    go to 290
  end if

  s(1,4)=(s(1,4)-s(1,3)*s(2,4)-entry3*s(3,4))/s(1,2)
!
!  Construct polynomial pieces. 
!
  break(1)=tau(1)
  l=1

  do i=1, ntau-1
    coef(1,l)=gtau(i)
    coef(3,l)=s(i,4)
    divdif=(gtau(i+1)-gtau(i))/s(i,1)
    z=s(i,5)
    if ( z-.5) 300, 310, 320

  300   continue

    if ( z == 0.) go to 330
    zeta=gam*z
    onemzt=1.0-zeta
    c=s(i+1,4)/6.0
    d=s(i,4)*s(i,6)
    l=l+1
    del=zeta*s(i,1)
    break(l)=tau(i)+del
    zt2=zeta**2
    alpha=alph(onemzt)
    factor=onemzt**2*alpha
    coef(1,l)=gtau(i)+divdif*del+s(i,1)**2*(d*onemzt*(factor-1.) &
      +c*zeta*(zt2-1.))
    coef(2,l)=divdif+s(i,1)*(d*(1.-3.*factor)+c*(3.*zt2-1.))
    coef(3,l)=6.*(d*alpha*onemzt+c*zeta)
    coef(4,l)=6.*(c-d*alpha)/s(i,1)
    coef(4,l-1)=coef(4,l)-6.*d*(1.-alpha)/(del*zt2)
    coef(2,l-1)=coef(2,l)-del*(coef(3,l)-(del/2.)*coef(4,l-1))
    go to 340

  310   coef(2,l)=divdif-s(i,1)*(2.*s(i,4)+s(i+1,4))/6.
    coef(4,l)=(s(i+1,4)-s(i,4))/s(i,1)
    go to 340

  320   onemzt=gam*(1.-z)
    if ( onemzt == 0.) go to 330
    zeta=1.-onemzt
    alpha=alph(zeta)
    c=s(i+1,4)*s(i,6)
    d=s(i,4)/6.
    del=zeta*s(i,1)
    break(l+1)=tau(i)+del
    coef(2,l)=divdif-s(i,1)*(2.*d+c)
    coef(4,l)=6.*(c*alpha-d)/s(i,1)
    l=l+1
    coef(4,l)=coef(4,l-1)+6.*(1.-alpha)*c/(s(i,1)*onemzt**3)
    coef(3,l)=coef(3,l-1)+del*coef(4,l-1)
    coef(2,l)=coef(2,l-1)+del*(coef(3,l-1)+(del/2.)*coef(4,l-1))
    coef(1,l)=coef(1,l-1)+del*(coef(2,l-1)+(del/2.)*(coef(3,l-1) &
      +(del/3.)*coef(4,l-1)))
    go to 340

  330   continue

    coef(2,l) = divdif
    coef(3,l) = 0.0
    coef(4,l) = 0.0

  340   continue

    l = l+1
    break(l) = tau(i+1)
    
  end do
  
  l = l-1
  k = 4
  iflag = 1
  return
  
! 360 iflag = 2
!     return
end
