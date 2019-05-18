function align_enum ( m, n )
!
!*******************************************************************************
!
!! ALIGN_ENUM counts the number of alignments of two sequences of M and N elements.
!
!
!  Discussion:
!
!    We assume that we have sequences A and B of M and N characters each.
!    An alignment of the two sequences is a rule matching corresponding
!    elements of one sequence to another.  Some elements of either sequence
!    can be matched to a null element.  If A(I1) and A(I2) are matched
!    to B(J1) and B(J2), and I1 < I2, then it must be the case that J1 < J2.
!
!    The 5 alignments of a sequence of 1 to a sequence of 2 are:
!
!          _1_   _2_   __3__   __4__   __5__
!
!      A:  1 -   - 1   - 1 -   - - 1   1 - -
!      B:  1 2   1 2   1 - 2   1 2 -   - 1 2
!
!    The formula is:
!
!      F(0,0) = 1
!      F(1,0) = 1
!      F(0,1) = 1
!      F(M,N) = F(M-1,N) + F(M-1,N-1) + F(M,N-1)
!
!    To compute F(M,N), it is not necessary to keep an M+1 by N+1
!    array in memory.  A vector of length N will do.
!
!    F(N,N) is approximately ( 1 + Sqrt(2) )**(2*N+1) / sqrt ( N )
!
!  Example:
!
!    The initial portion of the table is:
!
!  
!  M/N   0    1    2    3    4       5       6       7       8       9      10
!  
!   0    1    1    1    1    1       1       1       1       1       1       1
!   1    1    3    5    7    9      11      13      15      17      19      21
!   2    1    5   13   25   41      61      85     113     145     181     221
!   3    1    7   25   63  129     231     377     575     833    1159    1561
!   4    1    9   41  129  321     681    1289    2241    3649    5641    8361
!   5    1   11   61  231  681    1683    3653    7183   13073   22363   36365
!   6    1   13   85  377 1289    3653    8989   19825   40081   75517  134245
!   7    1   15  113  575 2241    7183   19825   48639  108545  224143  433905
!   8    1   17  145  833 3649   13073   40081  108545  265729  598417 1256465
!   9    1   19  181 1159 5641   22363   75517  224143  598417 1462563 3317445
!  10    1   21  221 1561 8361   36365  134245  433905 1256465 3317445 8097453
!
!  Reference:
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995, pages 186-190.
!
!  Modified:
!
!    24 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of elements of the two sequences.
!
!    Output, integer ALIGN_ENUM, the number of possible alignments of the
!    sequences.
!
  integer n
!
  integer align_enum
  integer fi(0:n)
  integer fim1j
  integer fim1jm1
  integer i
  integer j
  integer m
!
  if ( m < 0 ) then
    align_enum = 0
    return
  else if ( n < 0 ) then
    align_enum = 0
    return
  else if ( m == 0 ) then
    align_enum = 1
    return
  else if ( n == 0 ) then
    align_enum = 1
    return
  end if

  fi = 1

  do i = 1, m

    fim1jm1 = 1

    do j = 1, n

      fim1j = fi(j)

      fi(j) = fi(j) + fi(j-1) + fim1jm1

      fim1jm1 = fim1j

    end do
  end do

  align_enum = fi(n)

  return
end
subroutine bell ( b, n )
!
!*******************************************************************************
!
!! BELL returns the Bell numbers from 0 to N.
!
!
!  Discussion:
!
!    The Bell number B(N) is the number of restricted growth functions on N.
!
!    Note that the Stirling numbers of the second kind, S^m_n, count the
!    number of partitions of N objects into M classes, and so it is
!    true that
!
!      B(N) = S^1_N + S^2_N + ... + S^N_N.
!
!    The Bell numbers were named for Eric Temple Bell.
!
!  Definition:
!
!    The Bell number B(N) is defined as the number of partitions (of
!    any size) of a set of N distinguishable objects.  
!
!    A partition of a set is a division of the objects of the set into 
!    subsets.
!
!  Examples:
!
!    There are 15 partitions of a set of 4 objects:
!
!      (1234), 
!      (123) (4), 
!      (124) (3), 
!      (12) (34), 
!      (12) (3) (4), 
!      (134) (2), 
!      (13) (24), 
!      (13) (2) (4), 
!      (14) (23), 
!      (1) (234),
!      (1) (23) (4), 
!      (14) (2) (3), 
!      (1) (24) (3), 
!      (1) (2) (34), 
!      (1) (2) (3) (4).
!
!    and so B(4) = 15.
!
!  First values:
!
!     N         B(N)
!     0           1
!     1           1
!     2           2
!     3           5
!     4          15
!     5          52
!     6         203
!     7         877
!     8        4140
!     9       21147
!    10      115975
!
!  Recursion:
!
!    B(I) = SUM ( J = 1 to I ) Binomial ( I-1, J-1 ) * B(I-J)
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer B(0:N), the Bell numbers from 0 to N.
!
!    Input, integer N, the number of Bell numbers desired.
!
  integer n
!
  integer b(0:n)
  integer combo
  integer i
  integer j
!
  if ( n < 0 ) then
    return
  end if

  b(0) = 1

  do i = 1, n
    b(i) = 0
    do j = 1, i
      call combin2 ( i-1, j-1, combo )
      b(i) = b(i) + combo * b(i-j)
    end do
  end do

  return
end
function benford ( ival )
!
!*******************************************************************************
!
!! BENFORD returns the Benford probability of one or more significant digits.
!
!
!  Discussion:
!
!    Benford's law is an empirical formula explaining the observed
!    distribution of initial digits in lists culled from newspapers,
!    tax forms, stock market prices, and so on.  It predicts the observed
!    high frequency of the initial digit 1, for instance.
!
!    Note that the probabilities of digits 1 through 9 are guaranteed
!    to add up to 1, since
!      LOG10 ( 2/1 ) + LOG10 ( 3/2) + LOG10 ( 4/3 ) + ... + LOG10 ( 10/9 )
!      = LOG10 ( 2/1 * 3/2 * 4/3 * ... * 10/9 ) = LOG10 ( 10 ) = 1.
!
!  Formula:
!
!    Prob ( First significant digits are IVAL ) =
!      LOG10 ( ( IVAL + 1 ) / IVAL ).
!
!  Reference:
!
!    F Benford,
!    The Law of Anomalous Numbers,
!    Proceedings of the American Philosophical Society,
!    Volume 78, pages 551-572, 1938.
!
!    T P Hill,
!    The First Digit Phenomenon,
!    American Scientist, 
!    Volume 86, July/August 1998, pages 358 - 363.
!
!    R Raimi,
!    The Peculiar Distribution of First Digits,
!    Scientific American, 
!    December 1969, pages 109-119.
!
!  Modified:
!
!    13 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, the string of significant digits to be checked.
!    If IVAL is 1, then we are asking for the Benford probability that
!    a value will have first digit 1.  If IVAL is 123, we are asking for
!    the probability that the first three digits will be 123, and so on.
!
!    Note that IVAL must not be 0 or negative.
!
!    Output, real BENFORD, the Benford probability that an item taken
!    from a real world distribution will have the initial digits IVAL.
!
  real benford
  integer ival
!
  if ( ival <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BENFORD - Fatal error!'
    write ( *, * ) '  The input argument must be positive.'
    write ( *, * ) '  Your value was ', ival
    stop
  end if

  benford = log10 ( real ( ival + 1 ) / real ( ival ) )

  return
end
subroutine bern ( n, b )
!
!*******************************************************************************
!
!! BERN computes the value of the Bernoulli numbers B(0) through B(N).
!
!
!  First values:
!
!   B0  1                   =         1.00000000000
!   B1 -1/2                 =        -0.50000000000
!   B2  1/6                 =         1.66666666666
!   B3  0                   =         0
!   B4 -1/30                =        -0.03333333333
!   B5  0                   =         0
!   B6  1/42                =         0.02380952380
!   B7  0                   =         0
!   B8 -1/30                =        -0.03333333333
!   B9  0                   =         0
!  B10  5/66                =         0.07575757575
!  B11  0                   =         0
!  B12 -691/2730            =        -0.25311355311
!  B13  0                   =         0
!  B14  7/6                 =         1.16666666666
!  B15  0                   =         0
!  B16 -3617/510            =        -7.09215686274
!  B17  0                   =         0
!  B18  43867/798           =        54.97117794486
!  B19  0                   =         0
!  B20 -174611/330          =      -529.12424242424
!  B21  0                   =         0
!  B22  854,513/138         =      6192.123
!  B23  0                   =         0
!  B24 -236364091/2730      =    -86580.257
!  B25  0                   =         0
!  B26  8553103/6           =   1425517.16666
!  B27  0                   =         0
!  B28 -23749461029/870     = -27298231.0678
!  B29  0                   =         0
!  B30  8615841276005/14322 = 601580873.901
!
!  Recursion:
!
!    With C(N+1,K) denoting the standard binomial coefficient,
!
!    B(0) = 1.0
!    B(N) = - ( Sum(K=0 to N-1) C(N+1,K)*B(K) ) / C(N+1,N)
!
!  Warning:
!
!    This recursion, which is used in this routine, rapidly results
!    in significant errors.
!
!  Special Values:
!
!    Except for B(1), all Bernoulli numbers of odd index are 0.
!
!  Modified:
!
!    24 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the highest Bernoulli number to compute.
!
!    Output, real B(0:N), B(I) contains the I-th Bernoulli number.
!
  integer n
!
  real b(0:n)
  real b_sum
  integer i
  integer ido
  integer iwork(0:n+1)
  integer j
!
  if ( n < 0 ) then
    return
  end if

  b(0) = 1.0E+00

  if ( n < 1 ) then
    return
  end if

  b(1) = - 0.5E+00

  ido = 0
 
  do i = 2, n

    call comb_row ( ido, i+1, iwork )
    ido = 1
 
    if ( mod ( i, 2 ) == 1 ) then
 
      b(i) = 0.0E+00
 
    else
 
      b_sum = 0.0E+00
      do j = 0, i-1
        b_sum = b_sum + b(j) * real ( iwork(j) )
      end do
 
      b(i) = - b_sum / real ( iwork(i) )
 
    end if

  end do
 
  return
end
subroutine bern2 ( n, b )
!
!*******************************************************************************
!
!! BERN2 evaluates the Bernoulli numbers.
!
!
!  Discussion:
!
!    Note that the Bernoulli numbers grow rapidly.  Bernoulli number
!    62 is probably the last that can be computed on the VAX without
!    overflow.
!
!    A different method than that used in BERN is employed.
!
!  First values:
!
!   B0  1                   =         1.00000000000
!   B1 -1/2                 =        -0.50000000000
!   B2  1/6                 =         1.66666666666
!   B3  0                   =         0
!   B4 -1/30                =        -0.03333333333
!   B5  0                   =         0
!   B6  1/42                =         0.02380952380
!   B7  0                   =         0
!   B8 -1/30                =        -0.03333333333
!   B9  0                   =         0
!  B10  5/66                =         0.07575757575
!  B11  0                   =         0
!  B12 -691/2730            =        -0.25311355311
!  B13  0                   =         0
!  B14  7/6                 =         1.16666666666
!  B15  0                   =         0
!  B16 -3617/510            =        -7.09215686274
!  B17  0                   =         0
!  B18  43867/798           =        54.97117794486
!  B19  0                   =         0
!  B20 -174611/330          =      -529.12424242424
!  B21  0                   =         0
!  B22  854,513/138         =      6192.123
!  B23  0                   =         0
!  B24 -236364091/2730      =    -86580.257
!  B25  0                   =         0
!  B26  8553103/6           =   1425517.16666
!  B27  0                   =         0
!  B28 -23749461029/870     = -27298231.0678
!  B29  0                   =         0
!  B30  8615841276005/14322 = 601580873.901
!
!  Recursion:
!
!    With C(N+1,K) denoting the standard binomial coefficient,
!
!    B(0) = 1.0
!    B(N) = - ( Sum(K=0 to N-1) C(N+1,K)*B(K) ) / C(N+1,N)
!
!  Special Values:
!
!    Except for B(1), all Bernoulli numbers of odd index are 0.
!
!  Modified:
!
!    14 April 1999
!
!  Parameters:
!
!    Input, integer N, the highest order Bernoulli number to compute.
!
!    Output, real B(0:N), the requested Bernoulli numbers.
!
  integer, parameter :: kmax = 400
  real, parameter :: TOL = 1.0E+06
!
  integer n
!
  real altpi
  real b(0:n)
  integer i
  integer k
  real pi
  real sgn
  real sum2
  real t
  real term
!
  if ( n < 0 ) then
    return
  end if

  b(0) = 1.0E+00

  if ( n < 1 ) then
    return
  end if

  b(1) = - 0.5E+00

  if ( n < 2 ) then
    return
  end if

  altpi = log ( 2.0E+00 * pi ( ) )
!
!  Initial estimates for B(I), I = 2 to N
!
  b(2) = log ( 2.0E+00 )
  do i = 3, n
    if ( mod ( i, 2 ) == 1 ) then
      b(i) = 0.0E+00
    else
      b(i) = log ( real ( i * ( i - 1 ) ) ) + b(i-2)
    end if
  end do

  b(2) = 1.0E+00 / 6.0E+00

  if ( n <= 3 ) then
    return
  end if

  b(4) = - 1.0E+00 / 30.0E+00

  sgn = - 1.0E+00
 
  do i = 6, n, 2
 
    sgn = - sgn
    t = 2.0E+00 * sgn * exp ( b(i) - real ( i ) * altpi )
 
    sum2 = 1.0E+00

    do k = 2, kmax

      term = real ( k )**(-i)
      sum2 = sum2 + term

      if ( term <= TOL * sum2 ) then
        exit
      end if

    end do
 
    b(i) = t * sum2
 
  end do
 
  return
end
subroutine bern_poly ( n, x, bx )
!
!*******************************************************************************
!
!! BERN_POLY evaluates the Bernoulli polynomial of order N at X.
!
!
!  Special values:
!
!    B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
!
!    B'(N,X) = N * B(N-1,X)
!
!    B(N,X+1) - B(N,X) = N * X**(N-1)
!    B(N,X) = (-1)**N * B(N,1-X)
!
!  Formula:
!
!    B(N,X) = Sum (K=1 to N) B(K) * C(N,K) * X**(N-K)
!
!  First values:
!
!    B(0,X)  1
!    B(1,X)  X    - 1/2
!    B(2,X)  X**2 -   X      +  1/6
!    B(3,X)  X**3 - 3/2*X**2 +  1/2*X
!    B(4,X)  X**4 - 2*X**3   +      X**2 - 1/30
!    B(5,X)  X**5 - 5/2*X**4 +  5/3*X**3 - 1/6*X
!    B(6,X)  X**6 - 3*X**5   +  5/2*X**4 - 1/2*X**2 + 1/42
!    B(7,X)  X**7 - 7/2*X**6 +  7/2*X**5 - 7/6*X**3 + 1/6*X
!    B(8,X)  X**8 - 4*X**7   + 14/3*X**6 - 7/3*X**4 + 2/3*X**2 - 1/30
!
!  Modified:
!
!    24 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the Bernoulli polynomial to
!    be evaluated.  N must be 0 or greater.
!
!    Input, real X, the value of X at which the polynomial is to
!    be evaluated.
!
!    Output, real BX, the value of B(N,X).
!
  integer n
!
  real bx
  integer i
  integer ido
  integer iwork(0:n)
  real work(0:n)
  real x
!
  call bern ( n, work )
 
  ido = 0
  call comb_row ( ido, n, iwork )
 
  bx = 1.0E+00
  do i = 1, n
    bx = bx * x + work(i) * real ( iwork(i) )
  end do
 
  return
end
function bern_poly2 ( n, x )
!
!*******************************************************************************
!
!! BERN_POLY2 evaluates the N-th Bernoulli polynomial at X.
!
!
!  Special values:
!
!    BERN(N,0) = BERN(N,1) = B(N), the N-th Bernoulli number.
!
!    B'(N,X) = N*B(N-1,X).
!
!    B(N,X+1) - B(N,X) = N*X**(N-1)
!    B(N,X) = (-1)**N * B(N,1-X)
!
!  Formula:
!
!    B(N,X) = Sum (K=1 to N) B(K)*C(N,K)*X**(N-K)
!
!  First values:
!
!    B(0,X)  1
!    B(1,X)  X    - 1/2
!    B(2,X)  X**2 -   X      +  1/6
!    B(3,X)  X**3 - 3*X**2/2 +    X/2
!    B(4,X)  X**4 - 2*X**3   +    X**2   - 1/30
!    B(5,X)  X**5 - 5*X**4/2 +  5*X**3/3 -   X/6
!    B(6,X)  X**6 - 3*X**5   +  5*X**4/2 -   X**2/2 + 1/42
!    B(7,X)  X**7 - 7*X**6/2 +  7*X**5/2 - 7*X**3/6 +   X/6
!    B(8,X)  X**8 - 4*X**7   + 14*X**6/3 - 7*X**4/3 + 2*X**2/3 - 1/30
!
!  Modified:
!
!    14 April 1999
!
!  Parameters:
!
!    Input, integer N, the order of the Bernoulli polynomial to
!    be evaluated.  N must be 0 or greater.
!
!    Input, double precision X, the value at which the polynomial is to
!    be evaluated.
!
!    Output, double precision BERN_POLY2, the value of B(N,X).
!
  double precision dbern3
  double precision bern_poly2
  double precision fact
  integer i
  integer n
  double precision sum
  double precision x
!
  fact = 1.0D+00
  sum = dbern3 ( 0 )

  do i = 1, n
    fact = fact * dble ( n + 1 - i ) / dble ( i )
    sum = sum * x + fact * dbern3 ( i )
  end do

  bern_poly2 = sum

  return
end
function beta ( x, y )
!
!*******************************************************************************
!
!! BETA returns the value of the Beta function.
!
!
!  Formula:
!
!    BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
!
!  Restrictions:
!
!    Both X and Y must be greater than 0.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Properties:
!
!    BETA(X,Y) = BETA(Y,X).
!    BETA(X,Y) = Integral ( 0 <= T <= 1 ) T**(X-1) (1-T)**(Y-1) dT.
!
!  Parameters:
!
!    Input, real X, Y, the two parameters that define the Beta function.
!    X and Y must be greater than 0.
!
!    Output, real BETA, the value of the Beta function.
!
  real beta
  real gamma_log
  real x
  real y
!
  if ( x <= 0.0E+00 .or. y <= 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA - Fatal error!'
    write ( *, * ) '  Both X and Y must be greater than 0.'
    stop
  end if

  beta = exp ( gamma_log ( x ) + gamma_log ( y ) - gamma_log ( x + y ) )

  return
end
subroutine bp01 ( n, bern, x )
!
!*******************************************************************************
!
!! BP01 computes the values of the Bernstein polynomials at a point X.
!
!
!  Discussion:
!
!    The Bernstein polynomials are assumed to be based on [0,1].
!
!  Formula:
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
!
!  First values:
!
!    B(0,0,X) = 1
!
!    B(1,0,X) =      1-X
!    B(1,1,X) =                X
!
!    B(2,0,X) =     (1-X)**2
!    B(2,1,X) = 2 * (1-X)    * X
!    B(2,2,X) =                X**2
!
!    B(3,0,X) =     (1-X)**3
!    B(3,1,X) = 3 * (1-X)**2 * X
!    B(3,2,X) = 3 * (1-X)    * X**2
!    B(3,3,X) =                X**3
!
!    B(4,0,X) =     (1-X)**4
!    B(4,1,X) = 4 * (1-X)**3 * X
!    B(4,2,X) = 6 * (1-X)**2 * X**2
!    B(4,3,X) = 4 * (1-X)    * X**3
!    B(4,4,X) =                X**4
!
!  Special values:
!
!    B(N,I,1/2) = C(N,K) / 2**N
!
!  Modified:
!
!    14 April 1999
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein polynomials to be
!    used.  For any N, there is a set of N+1 Bernstein polynomials,
!    each of degree N, which form a basis for polynomials on [0,1].
!
!    Output, real BERN(0:N), the values of the N+1 Bernstein polynomials at X.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
  integer n
!
  real bern(0:n)
  integer i
  integer j
  real x
!
  if ( n == 0 ) then
 
    bern(0) = 1.0E+00
 
  else if ( n > 0 ) then
 
    bern(0) = 1.0E+00 - x
    bern(1) = x
 
    do i = 2, n
      bern(i) = x * bern(i-1)
      do j = i-1, 1, -1
        bern(j) = x * bern(j-1) + ( 1.0E+00 - x ) * bern(j)
      end do
      bern(0) = ( 1.0E+00 - x ) * bern(0)
    end do
 
  end if
 
  return
end
subroutine bpab ( n, bern, x, a, b )
!
!*******************************************************************************
!
!! BPAB evaluates at X the Bernstein polynomials based in [A,B].
!
!
!  Formula:
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)**(N-I) * (X-A)**I / (B-A)**N
!
!  First values:
!
!    B(0,0,X) =   1
!
!    B(1,0,X) = (      B-X                ) / (B-A)
!    B(1,1,X) = (                 X-A     ) / (B-A)
!
!    B(2,0,X) = (     (B-X)**2            ) / (B-A)**2
!    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)**2
!    B(2,2,X) = (                (X-A)**2 ) / (B-A)**2
!
!    B(3,0,X) = (     (B-X)**3            ) / (B-A)**3
!    B(3,1,X) = ( 3 * (B-X)**2 * (X-A)    ) / (B-A)**3
!    B(3,2,X) = ( 3 * (B-X)    * (X-A)**2 ) / (B-A)**3
!    B(3,3,X) = (                (X-A)**3 ) / (B-A)**3
!
!    B(4,0,X) = (     (B-X)**4            ) / (B-A)**4
!    B(4,1,X) = ( 4 * (B-X)**3 * (X-A)    ) / (B-A)**4
!    B(4,2,X) = ( 6 * (B-X)**2 * (X-A)**2 ) / (B-A)**4
!    B(4,3,X) = ( 4 * (B-X)    * (X-A)**3 ) / (B-A)**4
!    B(4,4,X) = (                (X-A)**4 ) / (B-A)**4
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein polynomials to be used.
!    For any N, there is a set of N+1 Bernstein polynomials, each of
!    degree N, which form a basis for polynomials on [A,B].
!
!    Output, real BERN(0:N), the values of the N+1 Bernstein polynomials at X.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Input, real A, B, the endpoints of the interval on which the
!    polynomials are to be based.  A and B should not be equal.
!
  integer n
!
  real a
  real b
  real bern(0:n)
  integer i
  integer j
  real x
!
  if ( b == a ) then
    write ( *, * ) ' '
    write ( *, * ) 'BPAB - Fatal error!'
    write ( *, * ) '  A = B = ', a
    stop
  end if

  if ( n == 0 ) then
 
    bern(0) = 1.0E+00
 
  else if ( n > 0 ) then
 
    bern(0) = ( b - x ) / ( b - a )
    bern(1) = ( x - a ) / ( b - a )
 
    do i = 2, n
      bern(i) = ( x - a ) * bern(i-1) / ( b - a )
      do j = i-1, 1, -1
        bern(j) = ( ( b - x ) * bern(j) + ( x - a ) * bern(j-1) ) / ( b - a )
      end do
      bern(0) = ( b - x ) * bern(0) / ( b - a )
    end do
 
  end if
 
  return
end
subroutine cardan ( n, x, s, cx )
!
!*******************************************************************************
!
!! CARDAN evaluates the Cardan polynomials.
!
!
!  First terms:
!
!    C( 0,S,X) = 2
!    C( 1,S,X) = X
!    C( 2,S,X) = X**2  -  2 S
!    C( 3,S,X) = X**3  -  3 S X
!    C( 4,S,X) = X**4  -  4 S X**2 +  2 S**2
!    C( 5,S,X) = X**5  -  5 S X**3 +  5 S**2 X
!    C( 6,S,X) = X**6  -  6 S X**4 +  9 S**2 X**2 -  2 S**3
!    C( 7,S,X) = X**7  -  7 S X**5 + 14 S**2 X**3 -  7 S**3 X
!    C( 8,S,X) = X**8  -  8 S X**6 + 20 S**2 X**4 - 16 S**3 X**2 +  2 S**4
!    C( 9,S,X) = X**9  -  9 S X**7 + 27 S**2 X**5 - 30 S**3 X**3 +  9 S**4 X
!    C(10,S,X) = X**10 - 10 S X**8 + 35 S**2 X**6 - 50 S**3 X**4 + 25 S**4 X**2 -  2 S**5
!    C(11,S,X) = X**11 - 11 S X**9 + 44 S**2 X**7 - 77 S**3 X**5 + 55 S**4 X**3 - 11 S**5 X
!
!  Recursion:
!
!    Writing the N-th polynomial in terms of its coefficients:
!
!      C(N,S,X) = Sum ( 0 <= I <= N ) D(N,I) * S**(N-I)/2 * X**I
!
!    then
!
!    D(0,0) = 1
!
!    D(1,1) = 1
!    D(1,0) = 0
!
!    D(N,N) = 1
!    D(N,K) = D(N-1,K-1) - D(N-2,K)
!
!  Reference:
!
!    Thomas Osler,
!    Cardan Polynomials and the Reduction of Radicals,
!    Mathematics Magazine, 
!    Volume 74, Number 1, February 2001, pages 26-32.
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
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real X, the point at which the polynomials are to be computed.
!
!    Input, real S, the value of the parameter, which must be positive.
!
!    Output, real CX(0:N), the values of the Cardan polynomials at X.
!
  integer n
!
  real cx(0:n)
  real fact
  integer i
  real s
  real s2
  real x
  real x2
!
  s2 = sqrt ( s )
  x2 = 0.5E+00 * x / s2

  call cheby1 ( n, x2, cx )

  fact = 1.0E+00

  do i = 0, n
    cx(i) = 2.0E+00 * fact * cx(i)
    fact = fact * s2
  end do
 
  return
end
subroutine cardan_coeff ( n, s, c )
!
!*******************************************************************************
!
!! CARDAN_COEFF computes the coefficients of the N-th Cardan polynomial.
!
!
!  First terms:
!
!    2
!    0       1
!   -2 S     0       1
!    0      -3 S     0       1
!    2 S**2  0      -4 S     0       1
!    0       5 S**2  0      -5 S     0       1
!   -2 S**3  0       9 S**2  0      -6 S     0       1
!    0       7 S**3  0      14 S**2  0      -7 S     0       1
!    2 S**4  0     -16 S**3  0      20 S**2  0      -8 S     0        1
!    0       9 S**4  0     -30 S**3  0      27 S**2  0      -9 S      0     1
!   -2 S**5  0      25 S**4  0     -50 S**3  0      35 S**2  0      -10 S   0   1
!    0     -11 S**5  0      55 S**4  0     -77 S**3  0     +44 S**2   0   -11 S 0 1
!
!  Reference:
!
!    Thomas Osler,
!    Cardan Polynomials and the Reduction of Radicals,
!    Mathematics Magazine, 
!    Volume 74, Number 1, February 2001, pages 26-32.
!
!  Modified:
!
!    01 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial
!
!    Input, real S, the value of the parameter, which must be positive.
!
!    Output, real C(0:N), the coefficients.  C(0) is the constant term,
!    and C(N) is the coefficient of X**N.
!
  integer n
!
  real c(0:n)
  real cm1(0:n)
  real cm2(0:n)
  integer i
  real s
!
  if ( n < 0 ) then
    return
  end if

  c(0) = 2.0E+00
  c(1:n) = 0.0

  if ( n == 0 ) then
    return
  end if

  cm1(0:n) = c(0:n)

  c(0) = 0.0E+00
  c(1) = 1.0E+00
  c(2:n) = 0.0

  do i = 2, n

    cm2(0:i-2) = cm1(0:i-2)
    cm1(0:i-1) = c(0:i-1)

    c(0) = 0.0E+00
    c(1:i) = cm1(0:i-1)
    c(0:i-2) = c(0:i-2) - s * cm2(0:i-2)

  end do

  return
end
subroutine catalan ( n, c )
!
!*******************************************************************************
!
!! CATALAN computes the Catalan numbers, from C(0) to C(N).
!
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) ) 
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = SUM ( I = 1 to N-1 ) C(I) * C(N-I)
!
!  Discussion:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which 
!       satisfy f(i) <= i for all i;
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Modified:
!
!    14 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of Catalan numbers desired.
!
!    Output, integer C(0:N), the Catalan numbers from C(0) to C(N).
!
  integer n
!
  integer c(0:n)
  integer i
!
  if ( n < 0 ) then
    return
  end if

  c(0) = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do
 
  return
end
subroutine catalan_row ( ido, n, irow )
!
!*******************************************************************************
!
!! CATALAN_ROW computes row N of Catalan's triangle.
!
!
!  Example:
!
!    I\J 0   1   2   3   4   5   6
!
!    0   1
!    1   1   1
!    2   1   2   2
!    3   1   3   5   5
!    4   1   4   9  14  14
!    5   1   5  14  28  42  42
!    6   1   6  20  48  90 132 132
!
!  Recursion:
!
!    C(0,0) = 1
!    C(I,0) = 1
!    C(I,J) = 0 for J > I
!    C(I,J) = C(I,J-1) + C(I-1,J)
!    C(I,I) is the I-th Catalan number.
!
!  Modified:
!
!    08 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IDO, indicates whether this is a call for
!    the 'next' row of the triangle.
!
!    0 means this is a startup call.  Row N is desired, but
!    presumably this is a first call, or row N-1 was not computed
!    on the previous call.
!
!    1 means this is not the first call, and row N-1 was computed
!    on the previous call.  In this case, much work can be saved
!    by using the information from the previous values of IROW
!    to build the next values.
!
!    Input, integer N, the row of the triangle desired.  
!
!    Output, integer IROW(0:N), the row of coefficients.
!
  integer n
!
  integer i
  integer ido
  integer irow(0:n)
  integer j
!
  if ( n <= 0 ) then
    return
  end if
 
  if ( ido == 1 ) then
 
    irow(0) = 1
    do j = 1, n
      irow(j) = irow(j) + irow(j-1)
    end do
 
  else
 
    irow(0) = 1
    irow(1:n) = 0
 
    do i = 1, n
      do j = 1, n
        irow(j) = irow(j) + irow(j-1)
      end do
    end do
 
  end if
 
  return
end
subroutine cheby1 ( n, x, cx )
!
!*******************************************************************************
!
!! CHEBY1 evaluates the Chebyshev polynomials of the first kind.
!
!
!  Discussion:
!
!    Chebyshev polynomials are useful as a basis for representing the
!    approximation of functions since they are well conditioned, in the sense
!    that in the interval [-1,1] they each have maximum absolute value 1.
!    Hence an error in the value of a coefficient of the approximation, of
!    size epsilon, is exactly reflected in an error of size epsilon between
!    the computed approximation and the theoretical approximation.
!
!    Typical usage is as follows, where we assume for the moment
!    that the interval of approximation is [-1,1].  The value
!    of N is chosen, the highest polynomial to be used in the
!    approximation.  Then the function to be approximated is
!    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
!    Chebyshev polynomial.  Let these values be denoted by F(XJ).
!
!    The coefficients of the approximation are now defined by
!      C(I) = 2/(N+1)*Sum(J=1 to N+1) F(XJ) T(I,XJ)
!    except that C(0) is given a value which is half that assigned
!    to it by the above formula,
!
!    and the representation is
!
!    F(X) approximated by Sum(J = 0 to N) C(J) T(J,X)
!
!    Now note that, again because of the fact that the Chebyshev polynomials
!    have maximum absolute value 1, if the higher order terms of the
!    coefficients C are small, then we have the option of truncating
!    the approximation by dropping these terms, and we will have an
!    exact value for maximum perturbation to the approximation that
!    this will cause.
!
!    It should be noted that typically the error in approximation
!    is dominated by the first neglected basis function (some multiple of
!    T(N+1,X) in the example above).  If this term were the exact error,
!    then we would have found the minimax polynomial, the approximating
!    polynomial of smallest maximum deviation from the original function.
!    The minimax polynomial is hard to compute, and another important
!    feature of the Chebyshev approximation is that it tends to behave
!    like the minimax polynomial while being easy to compute.
!
!    To evaluate a sum like Sum(J = 0 to N) C(J) T(J,X), Clenshaw's
!    recurrence formula is recommended instead of computing the
!    polynomial values, forming the products and summing.
!
!    Assuming that the coefficients C(J) have been computed
!    for J = 0 to N, then the coefficients of the representation of the
!    indefinite integral of the function may be computed by
!    B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1, with C(N+1)=0
!    and with B(0) arbitrary.  Also, the coefficients of the representation
!    of the derivative of the function may be computed by
!    D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0, with
!    D(N+1) = D(N)=0.
!
!    Some of the above may have to adjusted because of the irregularity of C(0).
!
!  Differential equation:
!
!    (1-X*X) Y'' - X Y' + N N Y = 0
!
!  Formula:
!
!    T(N,X) = COS(N*ARCCOS(X))
!
!  First terms:
!
!    T(0,X) =  1
!    T(1,X) =  1 X
!    T(2,X) =  2 X**2 -   1
!    T(3,X) =  4 X**3 -   3 X
!    T(4,X) =  8 X**4 -   8 X**2 +  1
!    T(5,X) = 16 X**5 -  20 X**3 +  5 X
!    T(6,X) = 32 X**6 -  48 X**4 + 18 X**2 - 1
!    T(7,X) = 64 X**7 - 112 X**5 + 56 X**3 - 7 X
!
!  Inequality:
!
!    ABS(T(N,X))< = 1 FOR -1<=X<=1
!
!  Orthogonality:
!
!    For integration over [-1,1] with weight
!
!    1/SQRT(1-X*X), <T(I,X),T(J,X)> = 0 if I.NE.J
!                                     = PI/2 if I.EQ.J.NE.0
!                                     = PI if I.EQ.J.EQ.0
!
!    A discrete orthogonality relation is also satisfied at each of
!    the N zeroes of T(N,X):  Sum(K = 1 to N) T(I,X) * T(J,X)
!                              = 0 if I.NE.J
!                              = N/2 if I.EQ.J.NE.0
!                              = N if I.EQ.J.EQ.0
!
!  Recursion:
!
!    T(0,X) = 1,
!    T(1,X) = X,
!    T(N,X) = 2*X*T(N-1,X) - T(N-2,X)
!
!    T'(N,X) = N*(-X*T(N,X)+T(N-1,X)) / ((1+X)*(1-X))
!
!  Special values:
!
!    T(N,1) = 1
!    T(N,-1) = (-1)**N
!    T(2N,0) = (-1)**N
!    T(2N+1,0) = 0
!    T(N,X) = (-1)**N * T(N,-X)
!
!  Zeroes and extrema:
!
!    M-th zero of T(N,X) is COS((2*M-1)*PI/(2*N)), M = 1 to N
!    M-th extremum of T(N,X) is COS(PI*M/N), M = 0 to N
!
!  Modified:
!
!    01 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real X, the point at which the polynomials are to be computed.
!
!    Output, real CX(0:N), the values of the N+1 Chebyshev polynomials.
!
  integer n
!
  real cx(0:n)
  integer i
  real x
!
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00

  if ( n < 1 ) then
    return
  end if

  cx(1) = x
 
  do i = 2, n
    cx(i) = 2.0E+00 * x * cx(i-1) - cx(i-2)
  end do
 
  return
end
subroutine cheby2 ( n, x, cx )
!
!*******************************************************************************
!
!! CHEBY2 evaluates the Chebyshev polynomials of the second kind.
!
!
!  Differential equation:
!
!    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
!
!  First terms:
!
!    U(0,X) =   1
!    U(1,X) =   2 X
!    U(2,X) =   4 X**2 -   1
!    U(3,X) =   8 X**3 -   4 X
!    U(4,X) =  16 X**4 -  12 X**2 +  1
!    U(5,X) =  32 X**5 -  32 X**3 +  6 X
!    U(6,X) =  64 X**6 -  80 X**4 + 24 X**2 - 1
!    U(7,X) = 128 X**7 - 192 X**5 + 80 X**3 - 8X
!
!  Recursion:
!
!    U(0,X) = 1,
!    U(1,X) = 2 * X,
!    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
!
!  Special values:
!
!    If ALFA = 1, the Gegenbauer polynomials reduce to the Chebyshev
!    polynomials of the second kind.
!
!  Modified:
!
!    01 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real X, the point at which the polynomials are to be computed.
!
!    Output, real CX(0:N), the values of the N+1 Chebyshev polynomials.
!
  integer n
!
  real cx(0:n)
  integer i
  real x
!
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00

  if ( n < 1 ) then
    return
  end if

  cx(1) = 2.0E+00 * x

  do i = 2, n
    cx(i) = 2.0E+00 * x * cx(i-1) - cx(i-2)
  end do
 
  return
end
subroutine combin ( n, k, cnk )
!
!*******************************************************************************
!
!! COMBIN computes the combinatorial coefficient C(N,K).
!
!
!  Method:
!
!    Real arithmetic is used, and C(N,K) is computed directly, via
!    Gamma functions, rather than recursively.
!
!  Definition:
!
!    C(N,K) is the number of distinct combinations of K objects
!    chosen from a set of N distinct objects.  A combination is
!    like a set, in that order does not matter.
!
!  Examples:
!
!    The number of combinations of 2 things chosen from 5 is 10.
!
!    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
!
!    The actual combinations may be represented as:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3), 
!      (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Formula:
!
!    C(N,K) = N! / ( (N-K)! * K! )
!
!  Modified:
!
!    17 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the value of N.
!
!    Input, integer K, the value of K.
!
!    Output, real CNK, the value of C(N,K)
!
  real arg
  real cnk
  real fack
  real facn
  real facnmk
  real gamma_log
  integer k
  integer n
!
  if ( n < 0 ) then

    cnk = 0.0E+00

  else if ( k == 0 ) then
 
    cnk = 1.0E+00
 
  else if ( k == 1 ) then
 
    cnk = real ( n )
 
  else if ( k > 1 .and. k < n-1 ) then
 
    arg = real ( n + 1 )
    facn = gamma_log ( arg )
 
    arg = real ( k + 1 )
    fack = gamma_log ( arg )
 
    arg = real ( n - k + 1 )
    facnmk = gamma_log ( arg )
 
    cnk = anint ( exp ( facn - fack - facnmk ) )
 
  else if ( k == n-1 ) then
 
    cnk = real ( n )
 
  else if ( k == n ) then
 
    cnk = 1.0E+00
 
  else
 
    cnk = 0.0
 
  end if
 
  return
end
subroutine combin2 ( n, k, icnk )
!
!*******************************************************************************
!
!! COMBIN2 computes the binomial coefficient C(N,K).
!
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!  Formula:
!
!    ICNK = C(N,K) = N! / ( K! * (N-K)! )
!
!  Reference:
!
!    M L Wolfson and H V Wright,
!    Combinatorial of M Things Taken N at a Time,
!    ACM algorithm 160,
!    Communications of the ACM,
!    April, 1963.
!
!  Modified:
!
!    17 January 1999
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, integer ICNK, the number of combinations of N
!    things taken K at a time.
!
  integer i
  integer icnk
  integer k
  integer mn
  integer mx
  integer n
!
  mn = min ( k, n-k )

  if ( mn < 0 ) then

    icnk = 0

  else if ( mn == 0 ) then

    icnk = 1

  else

    mx = max ( k, n-k )
    icnk = mx + 1

    do i = 2, mn
      icnk = ( icnk * ( mx + i ) ) / i
    end do

  end if

  return
end
subroutine comb_row ( ido, n, irow )
!
!*******************************************************************************
!
!! COMB_ROW computes row N of Pascal's triangle.
!
!
!  Discussion:
!
!    Row N contains the combinatorial coefficients
!
!      C(N,0), C(N,1), C(N,2), ... C(N,N)
!
!  Discussion:
!
!    The sum of the elements of row N is equal to 2**N.
!
!  Formula:
!
!    C(N,K) = N! / ( K! * (N-K)! )
!
!  First terms:
!
!     N K:0  1   2   3   4   5   6   7  8  9 10
!
!     0   1
!     1   1  1
!     2   1  2   1
!     3   1  3   3   1
!     4   1  4   6   4   1
!     5   1  5  10  10   5   1
!     6   1  6  15  20  15   6   1
!     7   1  7  21  35  35  21   7   1
!     8   1  8  28  56  70  56  28   8  1
!     9   1  9  36  84 126 126  84  36  9  1
!    10   1 10  45 120 210 252 210 120 45 10  1
!
!  Recursion:
!
!    C(N,K) = C(N-1,K-1)+C(N-1,K)
!
!  Special values:
!
!    C(N,0) = C(N,N) = 1
!    C(N,1) = C(N,N-1) = N
!    C(N,N-2) = Sum of integers from 1 to N.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IDO, indicates whether this is a call for
!    the 'next' row of the triangle.
!
!    0 means this is a startup call.  Row N is desired, but
!    presumably this is a first call, or row N-1 was not computed
!    on the previous call.
!
!    1 means this is not the first call, and row N-1 was computed
!    on the previous call.  In this case, much work can be saved
!    by using the information from the previous values of IROW
!    to build the next values.
!
!    Input, integer N, the row of the triangle desired.  The triangle
!    begins with row 0.
!
!    Output, integer IROW(N+1), the row of coefficients.
!    IROW(I) = C(N,I-1).
!
  integer n
!
  integer i
  integer ido
  integer irow(n+1)
  integer j
  integer k
!
  if ( n < 0 ) then
    return
  end if
 
  if ( ido == 1 ) then
 
    do i = 3, n+1
      j = n + 3 - i
      irow(j) = irow(j) + irow(j-1)
    end do
 
    irow(n+1) = 1
 
  else
 
    irow(1) = 1
    irow(2:n+1) = 0
 
    do k = 1, n
      do i = 2, k+1
        j = k + 3 - i
        irow(j) = irow(j) + irow(j-1)
      end do
    end do
 
  end if
 
  return
end
subroutine commul ( iarray, n, nfact, ncomb )
!
!*******************************************************************************
!
!! COMMUL computes a multinomial combinatorial coefficient.
!
!
!  Definition:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where IARRAY(1) objects are indistinguishable of type 1,
!    ... and IARRAY(K) are indistinguishable of type NFACT.
!
!  Formula:
!
!    NPERM = N! / ( IARRAY(1)! IARRAY(2)! ... IARRAY(NFACT)! )
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IARRAY(NFACT).
!    IARRAY contains the NFACT values used in the denominator.
!    Note that the sum of these entries should be N,
!    and that all entries should be nonnegative.
!
!    Input, integer N, determines the numerator.
!
!    Input, integer NFACT, the number of factors in the numerator.
!
!    Output, integer NCOMB, the value of the multinomial coefficient.
!
  integer nfact
!
  real arg
  real fack
  real facn
  real gamma_log
  integer i
  integer iarray(nfact)
  integer isum
  integer n
  integer ncomb
!
  do i = 1, nfact

    if ( iarray(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'COMMUL - Fatal error'
      write ( *, * ) '  Entry ', I, ' of IARRAY = ', iarray(i)
      write ( *, * ) '  But this value must be nonnegative.'
      stop
    end if

  end do
 
  isum = sum ( iarray(1:nfact) )

  if ( isum /= n ) then
    write ( *, * ) ' '
    write ( *, * ) 'COMMUL - Fatal error!'
    write ( *, * ) '  The sum of the IARRAY entries is ', isum
    write ( *, * ) '  But it must equal N = ', n
    stop
  end if
 
  arg = real ( n + 1 )
  facn = gamma_log ( arg )
 
  do i = 1, nfact
 
    arg = real ( iarray(i) + 1 )
    fack = gamma_log ( arg )
    facn = facn - fack
 
  end do
 
  ncomb = nint ( exp ( facn ) )
 
  return
end
function dbern3 ( n )
!
!*******************************************************************************
!
!! DBERN3 computes the value of the Bernoulli number B(N).
!
!
!  First values:
!
!     B0  1                   =         1.00000000000
!     B1 -1/2                 =        -0.50000000000
!     B2  1/6                 =         1.66666666666
!     B3  0                   =         0
!     B4 -1/30                =        -0.03333333333
!     B5  0                   =         0
!     B6  1/42                =         0.02380952380
!     B7  0                   =         0
!     B8 -1/30                =        -0.03333333333
!     B9  0                   =         0
!    B10  5/66                =         0.07575757575
!    B11  0                   =         0
!    B12 -691/2730            =        -0.25311355311
!    B13  0                   =         0
!    B14  7/6                 =         1.16666666666
!    B15  0                   =         0
!    B16 -3617/510            =        -7.09215686274
!    B17  0                   =         0
!    B18  43867/798           =        54.97117794486
!    B19  0                   =         0
!    B20 -174611/330          =      -529.12424242424
!    B21  0                   =         0
!    B22  854513/138          =      6192.123
!    B23  0                   =         0
!    B24 -236364091/2730      =    -86580.257
!    B25  0                   =         0
!    B26  8553103/6           =   1425517.16666
!    B27  0                   =         0
!    B28 -23749461029/870     = -27298231.0678
!    B29  0                   =         0
!    B30  8615841276005/14322 = 601580873.901
!
!  Recursion:
!
!    With C(N+1,K) denoting the standard binomial coefficient,
!
!    B(0) = 1.0
!    B(N) = - ( Sum(K=0 to N-1) C(N+1,K)*B(K) ) / C(N+1,N)
!
!  Special Values:
!
!    Except for B(1), all Bernoulli numbers of odd index are 0.
!
!  Modified:
!
!    14 April 1999
!
!  Parameters:
!
!    Input, integer N, the order of the Bernoulli number to compute.
!
!    Output, real DBERN3, the desired Bernoulli number.
!
  integer, parameter :: itmax = 1000
  double precision, parameter :: TOL = 5.0D-07
!
  double precision dbern3 
  double precision dfactorial
  double precision dpi
  integer i
  integer n
  double precision sum
  double precision term
!
  if ( n < 0 ) then
    dbern3 = 0.0D+00
    return
  end if

  if ( n == 0 ) then
    dbern3 = 1.0D+00
    return
  else if ( n == 1 ) then
    dbern3 = -0.5D+00
    return
  else if ( n == 2 ) then
    dbern3 = 1.0D+00 / 6.0D+00
    return
  else if ( mod ( n, 2 ) == 1 ) then
    dbern3 = 0.0D+00
    return
  end if

  sum = 0.0D+00
  do i = 1, itmax

    term = 1.0D+00 / dble ( i**n )
    sum = sum + term

    if ( abs ( term ) < TOL .or. abs ( term / sum ) < TOL ) then
      exit
    end if

  end do

  dbern3 = 2.0D+00 * sum * dfactorial ( n ) / ( 2.0D+00 * dpi ( ) )**n

  if ( mod ( n, 4 ) == 0 ) then
    dbern3 = - dbern3
  end if

  return
end
function deuler2 ( n )
!
!*******************************************************************************
!
!! DEULER2 computes the Euler numbers.
!
!
!  First terms:
!
!    E0  = 1
!    E1  = 0
!    E2  = -1
!    E3  = 0
!    E4  = 5
!    E5  = 0
!    E6  = -61
!    E7  = 0
!    E8  = 1385
!    E9  = 0
!    E10 = -50521
!    E11 = 0
!    E12 = 2702765
!    E13 = 0
!    E14 = -199360981
!    E15 = 0
!    E16 = 19391512145
!    E17 = 0
!    E18 = -2404879675441
!    E19 = 0
!    E20 = 370371188237525
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the Euler number to compute.
!
!    Output, double precision DEULER2, the value of E(N).
!
  integer, parameter :: itmax = 1000
!
  double precision deuler2
  double precision dfactorial
  double precision dpi
  double precision, save, dimension ( 0:6 ) :: e = &
    (/ 1.0D+00, -1.0D+00, 5.0D+00, -61.0D+00, 1385.0D+00, &
       -50521.0D+00, 2702765.0D+00 /)
  integer i
  integer n
  double precision sum1
  double precision term
!
  if ( n < 0 ) then
    deuler2 = 0.0D+00
    return
  end if

  if ( n == 0 ) then
    deuler2 = e(0)
    return
  end if

  if ( mod ( n, 2 ) == 1 ) then
    deuler2 = 0.0D+00
    return
  end if

  if ( n <= 12 ) then
    deuler2 = e(n/2)
    return
  end if

  sum1 = 0.0D+00
  do i = 1, itmax

    term = 1.0D+00 / dble ( ( 2 * i - 1 )**( n + 1 ) )

    if ( mod ( i, 2 ) == 1 ) then
      sum1 = sum1 + term
    else
      sum1 = sum1 - term
    end if

    if ( abs ( term ) < 1.0D-10 ) then
      exit
    else if ( abs ( term ) < 1.0D-8 * abs ( sum1 ) ) then
      exit
    end if

  end do

  deuler2 = 2.0D+00**( n + 2 ) * sum1 * dfactorial ( n ) / dpi ( )**( n + 1 )

  if ( mod ( n, 4 ) /= 0 ) then
    deuler2 = - deuler2
  end if

  return
end 
function dfactorial ( n )
!
!*******************************************************************************
!
!! DFACTORIAL computes the factorial N!
!
!
!  Formula:
!
!    DFACTORIAL( N ) = PRODUCT ( 1 <= I <= N ) I
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, DFACTORIAL is returned as 1.
!
!    Output, double precision DFACTORIAL, the factorial of N.
!
  double precision dfactorial
  integer i
  integer n
!
  dfactorial = 1.0D+00

  do i = 1, n
    dfactorial = dfactorial * dble ( i )
  end do

  return
end
function dpi ( )
!
!*******************************************************************************
!
!! DPI returns the value of Pi.
!
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision DPI, the value of Pi.
!
  double precision dpi
!
  dpi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
function e_constant ( )
!
!*******************************************************************************
!
!! E_CONSTANT returns the value of the base of the natural logarithm system.
!
!
!  Definition:
!
!    E = Limit ( N -> Infinity ) ( 1 + 1 / N )**N
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real E_CONSTANT, the base of the natural logarithm system.
!
  real e_constant
!
  e_constant = 2.718281828459045235360287E+00
 
  return
end
subroutine euler ( n, e )
!
!*******************************************************************************
!
!! EULER computes the Euler numbers.
!
!
!  First terms:
!
!    E0  = 1
!    E1  = 0
!    E2  = -1
!    E3  = 0
!    E4  = 5
!    E5  = 0
!    E6  = -61
!    E7  = 0
!    E8  = 1385
!    E9  = 0
!    E10 = -50521
!    E11 = 0
!    E12 = 2702765
!    E13 = 0
!    E14 = -199360981
!    E15 = 0
!    E16 = 19391512145
!    E17 = 0
!    E18 = -2404879675441
!    E19 = 0
!    E20 = 370371188237525
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the last Euler number to compute.
!
!    Output, integer E(0:N), the Euler numbers.
!
  integer n
!
  real cnk
  integer e(0:n)
  integer i
  integer j
  integer sgn
!
  if ( n < 0 ) then
    return
  end if

  e(0) = 1

  if ( n == 0 ) then
    return
  end if

  e(1) = 1
 
  do i = 2, n

    e(i) = 0
    sgn = +1
    do j = 1, i
      call combin ( 2*i, 2*j, cnk )
      e(i) = e(i) + sgn * nint ( cnk ) * e(i-j)
      sgn = - sgn
    end do
  end do
 
  return
end
function euler_constant ( )
!
!*******************************************************************************
!
!! EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
!
!
!  Discussion:
!
!    The Euler-Mascheroni constant is often denoted by a lower-case
!    Gamma.  Gamma is defined as
!
!      Gamma = limit ( M -> Infinity )
!        ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real EULER_CONSTANT, the value of the Euler-Mascheroni constant.
!
  real euler_constant
!
  euler_constant = 0.577215664901532860606512090082402431042E+00

  return
end
subroutine eulerian ( e, n )
!
!*******************************************************************************
!
!! EULERIAN computes the Eulerian number E(N,K).
!
!
!  Definition:
!
!    A run in a permutation is a sequence of consecutive ascending values.
!
!    E(N,K) is the number of permutations of N objects which contain
!    exactly K runs.
!
!  Examples:
!
!     N = 7
!
!     1     0     0     0     0     0     0
!     1     1     0     0     0     0     0
!     1     4     1     0     0     0     0
!     1    11    11     1     0     0     0
!     1    26    66    26     1     0     0
!     1    57   302   302    57     1     0
!     1   120  1191  2416  1191   120     1
!
!  Recursion:
!
!    E(N,K) = K * E(N-1,K) + (N-K+1) * E(N-1,K-1).
!
!  Properties:
!
!    E(N,1) = E(N,N) = 1.
!    E(N,K) = 0 if K <= 0 or K > N.
!    Sum ( K = 1 to N ) E(N,K) = N!.
!    X**N = Sum ( K = 0 to N ) COMB(X+K-1, N ) E(N,K)
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, 1986
!
!  Modified:
!
!    23 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, E(N,N), the first N rows of Eulerian numbers.
!
!    Input, integer N, the number of rows desired.
!
  integer n
!
  integer e(n,n)
  integer i
  integer j
!
  if ( n < 1 ) then
    return
  end if
!
!  Construct rows 1, 2, ..., N of the Eulerian triangle.
!
  e(1,1) = 1
  e(1,2:n) = 0

  do i = 2, n
    e(i,1) = 1
    do j = 2, n
      e(i,j) = j * e(i-1,j) + ( i - j + 1 ) * e(i-1,j-1)
    end do
  end do

  return
end
function euler_poly ( n, x )
!
!*******************************************************************************
!
!! EULER_POLY evaluates the N-th Euler polynomial at X.
!
!
!  First values:
!
!    E(0,X) = 1
!    E(1,X) = X - 1/2
!    E(2,X) = X**2 - X 
!    E(3,X) = X**3 - 3/2 X**2 + 1/4
!    E(4,X) = X**4 - 2*X**3 + X
!    E(5,X) = X**5 - 5/2 X**4 + 5/2 X**2 - 1/2
!    E(6,X) = X**6 - 3 X**5 + 5 X**3 - 3 X
!    E(7,X) = X**7 - 7/2 X**6 + 35/4 X**4 - 21/2 X**2 + 17/8
!    E(8,X) = X**8 - 4 X**7 + 14 X**5 - 28 X**3 + 17 X
!
!  Special values:
!
!    E'(N,X) = N * E(N-1,X)
!
!    E(N,1/2) = E(N) / 2**N, where E(N) is the N-th Euler number.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the Euler polynomial to
!    be evaluated.  N must be 0 or greater.
!
!    Input, double precision X, the value at which the polynomial is to
!    be evaluated.
!
!    Output, double precision EULER_POLY, the value of E(N,X).
!
  double precision bern_poly2
  double precision euler_poly
  integer n
  double precision x
!
  euler_poly = 2.0D+00 / dble ( n + 1 ) * ( bern_poly2 ( n+1, x ) &
    - bern_poly2 ( n+1, 0.5D+00*x ) * 2.0D+00**( n + 1 ) )

  return
end
recursive function f_hofstadter ( n ) result ( value )
!
!*******************************************************************************
!
!! F_HOFSTADTER computes the Hofstadter F sequence.
!
!
!  Discussion:
!
!    F(N) = 0                if N = 0
!         = N - F ( N - 1 ), otherwise.
!
!    F(N) is defined for all nonnegative integers, and turns out
!    to be equal to int ( ( N + 1 ) / 2 ).
!
!  Table:
!
!     N  F(N)
!    --  ----
!
!     0     0
!     1     1
!     2     1
!     3     2
!     4     2
!     5     3
!     6     3
!     7     4
!     8     4
!     9     5
!    10     5
!    11     6
!    12     6
!    13     7
!    14     7
!    15     8
!    16     8
!    17     9
!    18     9
!    19    10
!    20    10
!    21    11
!    22    11
!    23    12
!    24    12
!    25    13
!    26    13
!    27    14
!    28    14
!    29    15
!    30    15
!
!  Reference:
!
!    Douglas Hofstadter,
!    Goedel, Escher, Bach,
!    Basic Books, 1979.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the function.
!
!    Output, integer F_HOFSTADTER, the value of the function.
!
  integer n
  integer value
!
  if ( n <= 0 ) then
    value = 0
  else
    value = n - f_hofstadter ( n-1 )
  end if

  return
end
function factorial ( n )
!
!*******************************************************************************
!
!! FACTORIAL computes the factorial N!
!
!
!  Formula:
!
!    FACTORIAL( N ) = PRODUCT ( 1 <= I <= N ) I
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, FACTORIAL is returned as 1.
!
!    Output, real FACTORIAL, the factorial of N.
!
  real factorial
  integer i
  integer n
!
  factorial = 1.0E+00

  do i = 1, n
    factorial = factorial * real ( i )
  end do

  return
end
subroutine fibonacci_direct ( n, f )
!
!*******************************************************************************
!
!! FIBONACCI_DIRECT computes the N-th Fibonacci number directly.
!
!
!  Formula:
!
!      F(N) = ( PHIP**N - PHIM**N ) / SQRT(5)
!
!    where 
!
!      PHIP = ( 1 + SQRT(5) ) / 2, 
!      PHIM = ( 1 - SQRT(5) ) / 2.
!
!  Example:
!
!     N   F
!    --  --
!     0   0
!     1   1
!     2   1
!     3   2
!     4   3
!     5   5
!     6   8
!     7  13
!     8  21
!     9  34
!    10  55
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the Fibonacci number to compute.
!    N should be nonnegative.
!
!    Output, integer F, the value of the N-th Fibonacci number.
!
  integer f
  integer n
  real, parameter :: sqrt5 = 2.236068E+00
  real, parameter :: phim = ( 1.0E+00 - sqrt5 ) / 2.0E+00
  real, parameter :: phip = ( 1.0E+00 + sqrt5 ) / 2.0E+00
!
  if ( n < 0 ) then
    f = 0
  else
    f = nint ( ( phip**n - phim**n ) / sqrt ( 5.0E+00 ) )
  end if
 
  return
end
subroutine fibonacci_floor ( n, f, i )
!
!*******************************************************************************
!
!! FIBONACCI_FLOOR returns the largest Fibonacci number less than or equal to N.
!
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the positive integer whose Fibonacci "floor" is desired.
!
!    Output, integer F, the largest Fibonacci number less than or equal to N.
!
!    Output, integer I, the index of the F.
!
  integer f
  integer i
  integer n
!
  if ( n <= 0 ) then

    i = 0
    f = 0

  else

    i = int ( &
        log ( 0.5E+00 * real ( 2 * n + 1 ) * sqrt ( 5.0E+00 ) ) &
      / log ( 0.5E+00 * ( 1.0E+00 + sqrt ( 5.0E+00 ) ) ) )

    call fibonacci_direct ( i, f )

    if ( f > n ) then
      i = i - 1
      call fibonacci_direct ( i, f )
    end if

  end if

  return
end
subroutine fibonacci_recursive ( n, f )
!
!*******************************************************************************
!
!! FIBONACCI_RECURSIVE computes the first N Fibonacci numbers.
!
!
!  Algebraic equation:
!
!    The 'golden ratio' PHI = (1+SQRT(5))/2 satisfies the equation X*X-X-1=0.
!
!  Formula:
!
!    Let
!
!      PHIP = (1+SQRT(5))/2
!      PHIM = (1+SQRT(5))/2.
!
!    Then
!
!      F(N) = (PHIP**N+PHIM**N)/SQRT(5)
!
!    Moreover, F(N) can be computed by computing PHIP**N/SQRT(5) and rounding
!    to the nearest whole number.
!
!  First terms:
!
!      1
!      1
!      2
!      3
!      5
!      8
!     13
!     21
!     34
!     55
!     89
!    144
!
!    The 40th number is                  102,334,155.
!    The 50th number is               12,586,269,025.
!    The 100th number is 354,224,848,179,261,915,075.
!
!  Recursion:
!
!    F(1) = 1
!    F(2) = 1
!
!    F(N) = F(N-1)+F(N-2)
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest Fibonacci number to compute.
!
!    Output, integer F(N), the first N Fibonacci numbers.
!
  integer n
!
  integer f(n)
  integer i
!
  if ( n <= 0 ) then
    return
  end if

  f(1) = 1

  if ( n <= 1 ) then
    return
  end if

  f(2) = 1

  do i = 3, n
    f(i) = f(i-1) + f(i-2)
  end do
 
  return
end
recursive function g_hofstadter ( n ) result ( value )
!
!*******************************************************************************
!
!! G_HOFSTADTER computes the Hofstadter G sequence.
!
!
!  Discussion:
!
!    G(N) = 0                      if N = 0
!         = N - G ( G ( N - 1 ) ), otherwise.
!
!    G(N) is defined for all nonnegative integers.
!
!    The value of G(N) turns out to be related to the Zeckendorf 
!    representation of N as a sum of non-consecutive Fibonacci numbers.  
!    To compute G(N), determine the Zeckendorf representation:
!
!      N = Sum ( 1 <= I <= M ) F(I)
!
!    and reduce the index of each Fibonacci number by 1:
!
!      G(N) = Sum ( 1 <= I <= M ) F(I-1)
! 
!    However, this is NOT how the computation is done in this routine.
!    Instead, a straightforward recursive function call is defined
!    to correspond to the definition of the mathematical function.
!
!  Table:
!
!     N  G(N)  Zeckendorf   Decremented
!    --  ----  ----------   -----------
!
!     1   1    1            1
!     2   1    2            1             
!     3   2    3            2
!     4   3    3 + 1        2 + 1
!     5   3    5            3
!     6   4    5 + 1        3 + 1
!     7   4    5 + 2        3 + 1
!     8   5    8            5
!     9   6    8 + 1        5 + 1
!    10   6    8 + 2        5 + 1
!    11   7    8 + 3        5 + 2
!    12   8    8 + 3 + 1    5 + 2 + 1
!    13   8    13           8
!    14   9    13 + 1       8 + 1
!    15   9    13 + 2       8 + 1
!    16  10    13 + 3       8 + 2
!    17  11    13 + 3 + 1   8 + 2 + 1
!    18  11    13 + 5       8 + 3
!    19  12    13 + 5 + 1   8 + 3 + 1
!    20  12    13 + 5 + 2   8 + 3 + 1
!    21  13    21           13
!
!  Reference:
!
!    Douglas Hofstadter,
!    Goedel, Escher, Bach,
!    Basic Books, 1979.
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the function.
!
!    Output, integer G_HOFSTADTER, the value of the function.
!
  integer n
  integer value
!
  if ( n <= 0 ) then
    value = 0
  else
    value = n - g_hofstadter ( g_hofstadter ( n-1 ) )
  end if

  return
end
function gamma ( x )
!
!*******************************************************************************
!
!! GAMMA returns the value of the Gamma function at X.
!
!
!  Definition:
!
!    GAMMA(Z) = Integral(0 to Infinity) T**(Z-1) EXP(-T) DT
!
!  Recursion:
!
!    GAMMA(X+1) = X*GAMMA(X)
!
!  Restrictions:
!
!    X > 0 ( a software restriction).
!
!  Special values:
!
!    GAMMA(0.5) = SQRT(PI)
!    N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the point at which the Gamma function is desired.
!
!    Output, real GAMMA, the Gamma function of X.
!
  real gamma
  real gamma_log
  real x
!
  gamma = exp ( gamma_log ( x ) ) 
 
  return
end
function gamma_log ( x )
!
!*******************************************************************************
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for X > 12 is from reference 3, while approximations
!    for X < 12.0 are similar to those in reference 1, but are unpublished.
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    16 June 1999
!
!  Authors:
!
!    W. J. Cody and L. Stoltz
!    Argonne National Laboratory
!
!  References:
!
!    # 1)
!    W. J. Cody and K. E. Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Math. Comp.
!    Volume 21, 1967, pages 198-203.
!
!    # 2)
!    K. E. Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    # 3)
!    Hart, Et. Al.,
!    Computer Approximations,
!    Wiley and sons, New York, 1968.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.  X must be positive.
!
!    Output, real GAMMA_LOG, the logarithm of the Gamma function of X.
!    If X <= 0.0, or if overflow would occur, the program returns the
!    value XINF, the largest representable floating point number.
!
!*******************************************************************************
!
!  Explanation of machine-dependent constants
!
!  BETA   - radix for the floating-point representation.
!
!  MAXEXP - the smallest positive power of BETA that overflows.
!
!  XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!  XINF   - largest machine representable floating-point number;
!           approximately BETA**MAXEXP.
!
!  EPS    - The smallest positive floating-point number such that
!           1.0+EPS .GT. 1.0
!
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62E+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72E+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08E+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            XINF        EPS        FRTBIG
!
!  CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
!  Cyber 180/855
!    under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
!  IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
!  VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
!  VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
!
  real, parameter :: d1 = - 5.772156649015328605195174E-01
  real, parameter :: d2 =   4.227843350984671393993777E-01
  real, parameter :: d4 =   1.791759469228055000094023E+00
  real, parameter :: EPS = 1.19E-07
  real, parameter :: FRTBIG = 1.42E+09
  real, parameter :: PNT68 = 0.6796875E+00
  real, parameter :: SQRTPI = 0.9189385332046727417803297E+00
  real, parameter :: XBIG = 4.08E+36
  real, parameter :: XINF = 3.401E+38
!
  real, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728E-03, &
     8.4171387781295E-04, &
    -5.952379913043012E-04, &
     7.93650793500350248E-04, &
    -2.777777777777681622553E-03, &
     8.333333333333333331554247E-02, &
     5.7083835261E-03 /)
  real corr
  integer i
  real gamma_log
  real, parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888e+00, &
    2.018112620856775083915565e+02, &
    2.290838373831346393026739e+03, &
    1.131967205903380828685045e+04, &
    2.855724635671635335736389e+04, &
    3.848496228443793359990269e+04, &
    2.637748787624195437963534e+04, &
    7.225813979700288197698961e+03 /)
  real, parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064e+00, &
    5.424138599891070494101986e+02, &
    1.550693864978364947665077e+04, &
    1.847932904445632425417223e+05, &
    1.088204769468828767498470e+06, &
    3.338152967987029735917223e+06, &
    5.106661678927352456275255e+06, &
    3.074109054850539556250927e+06 /)
  real, parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062e+04, &
    2.426813369486704502836312e+06, &
    1.214755574045093227939592e+08, &
    2.663432449630976949898078e+09, &
    2.940378956634553899906876e+10, &
    1.702665737765398868392998e+11, &
    4.926125793377430887588120e+11, &
    5.606251856223951465078242e+11 /)
  real, parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036e+01, &
    1.113332393857199323513008e+03, &
    7.738757056935398733233834e+03, &
    2.763987074403340708898585e+04, &
    5.499310206226157329794414e+04, &
    6.161122180066002127833352e+04, &
    3.635127591501940507276287e+04, &
    8.785536302431013170870835e+03 /)
  real, parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942e+02, &
    7.765049321445005871323047e+03, &
    1.331903827966074194402448e+05, &
    1.136705821321969608938755e+06, &
    5.267964117437946917577538e+06, &
    1.346701454311101692290052e+07, &
    1.782736530353274213975932e+07, &
    9.533095591844353613395747e+06 /)
  real, parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843e+03, &
    6.393885654300092398984238e+05, &
    4.135599930241388052042842e+07, &
    1.120872109616147941376570e+09, &
    1.488613728678813811542398e+10, &
    1.016803586272438228077304e+11, &
    3.417476345507377132798597e+11, &
    4.463158187419713286462081e+11 /)
  real res
  real x
  real xden
  real xm1
  real xm2
  real xm4
  real xnum
  real xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0E+00 .or. x > XBIG ) then
    gamma_log = XINF
    return
  end if

  if ( x <= EPS ) then

    res = - log ( x )

  else if ( x <= 1.5E+00 ) then

    if ( x < PNT68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0E+00
      xm1 = ( x - 0.5E+00 ) - 0.5E+00
    end if

    if ( x <= 0.5E+00 .or. x >= PNT68 ) then

      xden = 1.0E+00
      xnum = 0.0E+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5E+00 ) - 0.5E+00
      xden = 1.0E+00
      xnum = 0.0E+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0E+00 ) then

    xm2 = x - 2.0E+00
    xden = 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0E+00 ) then

    xm4 = x - 4.0E+00
    xden = - 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0E+00

    if ( x <= FRTBIG ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + SQRTPI - 0.5E+00 * corr
    res = res + x * ( corr - 1.0E+00 )

  end if

  gamma_log = res

  return
end
subroutine gegenbauer ( n, alfa, x, cx )
!
!*******************************************************************************
!
!! GEGENBAUER computes the Gegenbauer polynomials C(I,ALFA)(X) for I = 1 to N.
!
!
!  Differential equation:
!
!    (1-X*X) Y'' - (2 ALFA + 1) X Y' + N (N + 2 ALFA) Y = 0
!
!  Recursion:
!
!    C(0,ALFA) = 1,
!    C(1,ALFA) = 2*ALFA*X
!    C(N,ALFA) = ( 2*(N-1+ALFA)*X*C(N-1,ALFA) - (N-2+2*ALFA)*C(N-2,ALFA) )/N
!
!  Restrictions:
!
!    ALFA must be greater than -0.5, and may not be equal to 0.
!
!  Special values:
!
!    If ALFA = 1, the Gegenbauer polynomials reduce to the Chebyshev
!    polynomials of the second kind.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ALFA, a parameter which is part of the definition of
!    the Gegenbauer polynomials.  It must be greater than -0.5.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 Gegenbauer
!    polynomials at the point X.  
!
  integer n
!
  real alfa
  real cx(0:n)
  integer i
  real x
!
  if ( alfa <= -0.5E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GEGENBAUER - Fatal error!'
    write ( *, * ) '  Illegal value of ALFA = ', alfa
    write ( *, * ) '  but ALFA must be greater than -0.5.'
    return
  end if
 
  if ( alfa == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GEGENBAUER - Fatal error!'
    write ( *, * ) '  Illegal value of ALFA = ', alfa
    write ( *, * ) '  but ALFA must be nonzero.'
    return
  end if

  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 2.0E+00 * alfa * x

  do i = 2, n
    cx(i) = ( 2.0E+00 * ( real ( i ) - 1.0E+00 + alfa ) * x * cx(i-1) &
      - ( real ( i ) - 2.0E+00 + 2.0E+00 * alfa ) * cx(i-2) ) / real ( i )
  end do
 
  return
end
function hail ( n )
!
!*******************************************************************************
!
!! HAIL computes the hail function.
!
!
!  Discussion:
!
!    Starting with a positive integer N, we divide it by 2 if it is
!    even, or triple and add 1 if odd, and repeat this process until
!    reaching the value 1.  The number of times the process is carried
!    out is the value of the hail function for the given starting value.
!
!    Actually, HAIL is not well defined, since it is not known if
!    the above process actually terminates at 1 for every starting value N.
!
!  Example:
!
!     N  Sequence                                                  Hail
!
!     1                                                               0
!     2   1                                                           1
!     3  10,  5, 16,  8,  4,  2,  1                                   7
!     4   2   1                                                       2
!     5  16,  8,  4,  2,  1                                           5
!     6   3, 10,  5, 16,  8,  4,  2,  1                               8
!     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   16
!     8   4,  2,  1                                                   3
!     9  28, 14,  7, ...                                             19
!    10   5, 16,  8,  4,  2,  1                                       6
!    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          14
!    12   6,  3, 10,  5, 16,  8,  4,  2,  1                           9
!
!  Modified:
!
!    28 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the starting value for the hail sequence.
!
!    Output, integer HAIL, the number of steps before the hail sequence 
!    reached 1.
!
  integer hail
  integer k
  integer m
  integer n
!
  k = 0
  m = n

  if ( n > 0 ) then

    do while ( m /= 1 )
      k = k + 1
      if ( mod ( m, 2 ) == 0 ) then
        m = m / 2
      else
        m = 3 * m + 1
      end if
    end do

  end if

  hail = k

  return
end
recursive function h_hofstadter ( n ) result ( value )
!
!*******************************************************************************
!
!! H_HOFSTADTER computes the Hofstadter H sequence.
!
!
!  Discussion:
!
!    H(N) = 0                          if N = 0
!         = N - H ( H ( H ( N - 1 ) ), otherwise.
!
!    H(N) is defined for all nonnegative integers.
!
!  Table:
!
!     N  H(N)
!    --  ----
!
!     0     0
!     1     1
!     2     1
!     3     2
!     4     3
!     5     4
!     6     4
!     7     5
!     8     5
!     9     6
!    10     7
!    11     7
!    12     8
!    13     9
!    14    10
!    15    10
!    16    11
!    17    12
!    18    13
!    19    13
!    20    14
!    21    14
!    22    15
!    23    16
!    24    17
!    25    17
!    26    18
!    27    18
!    28    19
!    29    20
!    30    20
!
!  Reference:
!
!    Douglas Hofstadter,
!    Goedel, Escher, Bach,
!    Basic Books, 1979.
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the function.
!
!    Output, integer H_HOFSTADTER, the value of the function.
!
  integer n
  integer value
!
  if ( n <= 0 ) then
    value = 0
  else
    value = n - h_hofstadter ( h_hofstadter ( h_hofstadter ( n-1 ) ) )
  end if

  return
end
subroutine hermite ( n, x, cx )
!
!*******************************************************************************
!
!! HERMITE evaluates the Hermite polynomials at X.
!
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X**2     - 2
!      8 X**3     - 12 X
!     16 X**4     - 48 X**2     + 12
!     32 X**5    - 160 X**3    + 120 X
!     64 X**6    - 480 X**4    + 720 X**2    - 120
!    128 X**7   - 1344 X**5   + 3360 X**3   - 1680 X
!    256 X**8   - 3584 X**6  + 13440 X**4  - 13440 X**2   + 1680
!    512 X**9   - 9216 X**7  + 48384 X**5  - 80640 X**3  + 30240 X
!   1024 X**10 - 23040 X**8 + 161280 X**6 - 403200 X**4 + 302400 X**2 - 30240
!
!  Recursion:
!
!    H(0) = 1,
!    H(1) = 2*X,
!    H(N) = 2*X*H(N-1)-2*(N-1)*H(N-2)
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 Hermite
!    polynomials at the point X.
!
  integer n
!
  real cx(0:n)
  integer i
  real x
! 
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 2.0E+00 * x
 
  do i = 2, n
    cx(i) = 2.0E+00 * x * cx(i-1) - 2.0E+00 * real ( i - 1 ) * cx(i-2)
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
!    23 September 2000
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
subroutine jacobi ( n, alfa, beta, x, cx )
!
!*******************************************************************************
!
!! JACOBI evaluates the Jacobi polynomials at X.
!
!
!  Differential equation:
!
!    (1-X*X) Y'' + (BETA-ALFA-(ALFA+BETA+2) X) Y' + N (N+ALFA+BETA+1) Y = 0
!
!  Recursion:
!
!    P(0,ALFA,BETA) = 1,
!
!    P(1,ALFA,BETA) = (1+0.5*(ALFA+BETA))*X + 0.5*(ALFA-BETA)
!
!    2*N*(N+ALFA+BETA)*(2*N-2+ALFA+BETA) * P(N,ALFA,BETA)  = 
!    (2*N+ALFA+BETA-1)*
!    ((ALFA**2-BETA**2)+(2*N+ALFA+BETA)*(2*N+ALFA+BETA-2)*X)
!    * P(N-1,ALFA,BETA)
!    -2*(N-1+ALFA)*(N-1+BETA)*(2*I+ALFA+BETA) * P(N-2,ALFA,BETA)
!
!  Restrictions:
!
!    ALFA > -1
!    BETA > -1
!
!  Special values:
!
!    P(N,ALFA,BETA)(1) = (N+ALFA)!/(N!*ALFA!) for integer ALFA.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.  Note
!    that polynomials 0 through N will be computed.
!
!    Input, real ALFA, one of the parameters defining the Jacobi
!    polynomials, ALFA must be greater than -1.
!
!    Input, real BETA, the second parameter defining the Jacobi
!    polynomials, BETA must be greater than -1.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 Jacobi
!    polynomials at the point X.
!
  integer n
!
  real alfa
  real beta
  real cx(0:n)
  real c1
  real c2
  real c3
  real c4
  integer i
  real x
!
  if ( alfa <= -1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Jacobi - Fatal error!'
    write ( *, * ) '  Illegal input value of ALFA = ', alfa
    write ( *, * ) '  But ALFA must be greater than -1.'
    stop
  end if
 
  if ( beta <= -1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Jacobi - Fatal error!'
    write ( *, * ) '  Illegal input value of BETA = ', beta
    write ( *, * ) '  But BETA must be greater than -1.'
    stop
  end if
  
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = ( 1.0E+00 + 0.5E+00 * ( alfa + beta ) ) * x + 0.5E+00 * ( alfa - beta )
 
  do i = 2, n

    c1 = 2.0E+00 * real ( i ) * ( real ( i ) + alfa + beta ) &
      * ( real ( 2 * i - 2 ) + alfa + beta )

    c2 = ( real ( 2 * i - 1 ) + alfa + beta ) &
      * ( real ( 2 * i ) + alfa + beta ) &
      * ( real ( 2 * i - 2 ) + alfa + beta )

    c3 = ( real ( 2 * i - 1 ) + alfa + beta ) &
      * ( alfa + beta ) * ( alfa - beta )

    c4 = - real ( 2 ) * ( real ( i - 1 ) + alfa ) &
      * ( real ( i - 1 ) + beta )  * ( real ( 2 * i ) + alfa + beta )

    cx(i) = ( ( c3 + c2 * x ) * cx(i-1) + c4 * cx(i-2) ) / c1

  end do

  return
end
subroutine laguerre ( n, x, cx )
!
!*******************************************************************************
!
!! LAGUERRE evaluates the Laguerre polynomials at X.
!
!
!  Differential equation:
!
!    X * Y'' + (1-X) * Y' + N * Y = 0
!
!  First terms:
!
!     1
!    -X    +  1
!     X**2 -  4 X     +  2
!    -X**3 +  9 X**2 -  18 X    +    6
!     X**4 - 16 X**3 +  72 X**2 -   96 X +      24
!    -X**5 + 25 X**4 - 200 X**3 +  600 X**2 -  600 x    +  120
!     X**6 - 36 X**5 + 450 X**4 - 2400 X**3 + 5400 X**2 - 4320 X + 720
!    -X**7 + 49 X**6 - 882 X**5 + 7350 X**4 - 29400 X**3 
!      + 52920 X**2 - 35280 X + 5040
!
!  Recursion:
!
!    L(0) = 1,
!    L(1) = 1-X,
!    L(N) = (2*N-1-X)*L(N-1)-(N-1)*(N-1)*L(N-2)
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 Laguerre
!    polynomials at the point X.
!
  integer n
!
  real cx(0:n)
  integer i
  real x
!
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 1.0E+00 - x
 
  do i = 2, n

    cx(i) = ( real ( 2 * i - 1 ) - x ) * cx(i-1) &
      - real ( ( i - 1 )**2 ) * cx(i-2) 

  end do

  return
end
subroutine laguerre_gen ( n, alfa, x, cx )
!
!*******************************************************************************
!
!! LAGUERRE_GEN evaluates the generalized Laguerre polynomials at X.
!
!
!  Differential equation:
!
!    X * Y'' + (ALFA+1-X) * Y' + N * Y = 0
!
!  Recursion:
!
!    L(0,ALFA) = 1
!    L(1,ALFA) = 1+ALFA-X
!
!    L(N,ALFA) = ( (2*N-1+ALFA-X)*L(N-1,ALFA)-(N-1+ALFA)*L(N-2,ALFA) ) / N
!
!  Restrictions:
!
!    ALFA > -1
!
!  Special values:
!
!    For ALFA = 0, the generalized Laguerre polynomial becomes the
!    Laguerre polynomial.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ALFA, a parameter which is part of the definition of
!    the generalized Laguerre polynomials.  It must be greater than -1.
!
!    Input, real X, the point at which the polynomials are to be
!    evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 generalized
!    Laguerre polynomials at the point X.
!
  integer n
!
  real alfa
  real cx(0:n)
  integer i
  real x
!
  if ( alfa <= -1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LAGUERRE_GEN - Fatal error!'
    write ( *, * ) '  Input value of ALFA is ', alfa
    write ( *, * ) '  but ALFA must be greater than -1.'
    stop
  end if
 
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00

  if ( n == 0 ) then
    return
  end if

  cx(1) = 1.0E+00 + alfa - x

  do i = 2, n
    cx(i) = ( ( real ( 2 * i - 1 ) + alfa - x ) * cx(i-1) &
      - ( real ( i - 1 ) + alfa ) * cx(i-2) ) / real ( i )
  end do

  return
end
subroutine laguerre_lnm ( n, m, x, cx )
!
!*******************************************************************************
!
!! LAGUERRE_LNM evaluates the associated Laguerre polynomials Lnm at X.
!
!
!  Differential equation:
!
!    X Y'' + (M+1-X) Y' + (N-M) Y = 0
!
!  First terms:
!
!    M = 0
!
!    L(0,0) =   1
!    L(1,0) =  -X    +  1
!    L(2,0) =   X**2 -  4 X     +  2
!    L(3,0) =  -X**3 +  9 X**2 -  18 X    +    6
!    L(4,0) =   X**4 - 16 X**3 +  72 X**2 -   96 X +      24
!    L(5,0) =  -X**5 + 25 X**4 - 200 X**3 +  600 X**2 -  600 x    +  120
!    L(6,0) =   X**6 - 36 X**5 + 450 X**4 - 2400 X**3 + 5400 X**2 - 4320 X + 720
!
!    M = 1
!
!    L(0,1) =    0
!    L(1,1) =   -1,
!    L(2,1) =    2 X-4,
!    L(3,1) =   -3 X**2 + 18 X - 18,
!    L(4,1) =    4 X**3 - 48 X**2 + 144 X - 96
!
!    M = 2
!
!    L(0,2) =    0
!    L(1,2) =    0,
!    L(2,2) =    2,
!    L(3,2) =   -6 X + 18,
!    L(4,2) =   12 X**2 - 96 X + 144
!
!    M = 3
!
!    L(0,3) =    0
!    L(1,3) =    0,
!    L(2,3) =    0,
!    L(3,3) =   -6,
!    L(4,3) =   24 X - 96
!
!    M = 4
!
!    L(0,4) =    0
!    L(1,4) =    0
!    L(2,4) =    0
!    L(3,4) =    0
!    L(4,4) =   24
!
!  Recursion:
!
!    if N < M:
!      L(N,M)   = 0
!    if N = M:
!      L(N,M)   = (-1)**M * M! 
!    if N = M+1:
!      L(N,M) = (M+1)*(M+1-X)*L(M,M)
!    if N = M+2 or greater:
!      L(N,M)  = -( (X+M-2*N+1)*L(N-1,M) + (N-1)*(N-1)*L(N-2,M) )*N / (N-M)
!
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials are equal to the
!    Laguerre polynomials.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, integer M, the parameter.  M must be nonnegative.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 Laguerre
!    polynomials at the point X.
!
  integer n
!
  real cx(0:n)
  integer i
  integer ifact
  integer m
  real x
!
  if ( m < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LAGUERRE_LNM - Fatal error!'
    write ( *, * ) '  Input value of M = ', m
    write ( *, * ) '  but M must be nonnegative.'
    stop
  end if
 
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LAGUERRE_LNM - Fatal error!'
    write ( *, * ) '  Input value of N = ', n
    write ( *, * ) '  but N must be positive.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(0:m-1) = 0.0E+00

  ifact = 1
  do i = 1, m
    ifact = - ifact * i
  end do
 
  cx(m) = real ( ifact )
  cx(m+1) = real ( m + 1 ) * ( real ( m + 1 ) - x ) * cx(m)

  do i = m+2, n
    cx(i) = - ( real ( i ) * &
      ( x + real ( m - 2 * i + 1 ) ) * cx(i-1) &
      + real ( ( i - 1 )**2 ) * cx(i-2) ) / real ( i - m )
  end do

  return
end
subroutine legendre_pn ( n, x, cx, cpx )
!
!*******************************************************************************
!
!! LEGENDRE_PN evaluates the Legendre polynomials P(N)(X) at X.
!
!
!  Differential equation:
!
!    (1-X*X) * P(N)(X)'' - 2 * X * P(N)(X)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0)(X) =       1
!    P( 1)(X) =       1 X
!    P( 2)(X) =  (    3 X**2 -       1)/2
!    P( 3)(X) =  (    5 X**3 -     3 X)/2
!    P( 4)(X) =  (   35 X**4 -    30 X**2 +     3)/8
!    P( 5)(X) =  (   63 X**5 -    70 X**3 +    15 X)/8
!    P( 6)(X) =  (  231 X**6 -   315 X**4 +   105 X**2 -     5)/16
!    P( 7)(X) =  (  429 X**7 -   693 X**5 +   315 X**3 -    35 X)/16
!    P( 8)(X) =  ( 6435 X**8 - 12012 X**6 +  6930 X**4 -  1260 X**2 +   35)/128
!    P( 9)(X) =  (12155 X**9 - 25740 X**7 + 18018 X**5 -  4620 X**3 +  315 X)/128
!    P(10)(X) =  (46189 X**10-109395 X**8 + 90090 X**6 - 30030 X**4 + 3465 X**2
!                 -63 ) /256
!
!  Recursion:
!
!    P(0)(X) = 1
!    P(1)(X) = X
!    P(N)(X) = ( (2*N-1)*X*P(N-1)(X)-(N-1)*P(N-2)(X) ) / N
!
!    P'(0)(X) = 0
!    P'(1)(X) = 1
!    P'(N)(X) = ( (2*N-1)*(P(N-1)(X)+X*P'(N-1)(X)-(N-1)*P'(N-2)(X) ) / N
!
!  Formula:
!
!    P(N)(X) = (1/2**N) * Sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X**(N-2*M)
!
!  Orthogonality:
!
!    Integral ( -1 <= X <= 1 ) P(I)(X) * P(J)(X) dX 
!      = 0 if I =/= J
!      = 2 / ( 2I+1 ) if I = J.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!
!      C0*P(0)(X) + C1*P(1)(X) + ... + CN*P(N)(X)
!
!    where
!
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I)(X) dx.
!
!  Special values:
!
!    P(N)(1) = 1.
!    P(N)(-1) = (-1)**N.
!    | P(N)(X) | <= 1 in [-1,1].
!
!    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
!    function of the first kind and order N equals the Legendre polynomial
!    of the first kind and order N.
!
!    The N zeroes of P(N)(X) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Output, real CX(0:N), the values of the Legendre polynomials 
!    of order 0 through N at the point X.
!
!    Output, real CPX(0:N), the values of the derivatives of the
!    Legendre polynomials of order 0 through N at the point X.
!
  integer n
!
  real cx(0:n)
  real cpx(0:n)
  integer i
  real x
!
  if ( n < 0 ) then
    return
  end if

  cx(0) = 1.0E+00
  cpx(0) = 0.0E+00

  if ( n < 1 ) then
    return
  end if

  cx(1) = x
  cpx(1) = 1.0E+00
 
  do i = 2, n
 
    cx(i) = ( real ( 2 * i - 1 ) * x * cx(i-1) &
      - real ( i - 1 ) * cx(i-2) ) / real ( i )
 
    cpx(i) = ( real ( 2 * i - 1 ) * ( cx(i-1) + x * cpx(i-1) ) &
      - real ( i - 1 ) * cpx(i-2) ) / real ( i )
 
  end do
 
  return
end
subroutine legendre_pnm ( n, m, x, cx )
!
!*******************************************************************************
!
!! LEGENDRE_PNM evaluates the associated Legendre polynomial Pnm(x) of the first kind.
!
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N)(X) )
!
!    P00 =    1
!    P10 =    1 X
!    P20 = (  3 X**2 -   1)/2
!    P30 = (  5 X**3 -   3 X)/2
!    P40 = ( 35 X**4 -  30 X**2 +   3)/8
!    P50 = ( 63 X**5 -  70 X**3 +  15 X)/8
!    P60 = (231 X**6 - 315 X**4 + 105 X**2 -  5)/16
!    P70 = (429 X**7 - 693 X**5 + 315 X**3 - 35 X)/16
!
!    M = 1
!
!    P01 =   0
!    P11 =   1 * SQRT(1-X*X)
!    P21 =   3 * SQRT(1-X*X) * X
!    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
!    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
!
!    M = 2
!
!    P02 =   0
!    P12 =   0
!    P22 =   3 * (1-X*X)
!    P32 =  15 * (1-X*X) * X
!    P42 = 7.5 * (1-X*X) * (7*X*X-1)
!
!    M = 3
!
!    P03 =   0
!    P13 =   0
!    P23 =   0
!    P33 =  15 * (1-X*X)**1.5
!    P43 = 105 * (1-X*X)**1.5 * X
!
!    M = 4
!
!    P04 =   0
!    P14 =   0
!    P24 =   0
!    P34 =   0
!    P44 = 105 * (1-X*X)**2
!
!  Recursion:
!
!    if N < M:
!      P(N,M) = 0
!    if N = M:
!      P(N,M) = (2*M-1)!! * (1-X*X)**(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      P(N,M) = X*(2*M+1)*P(M,M)
!    if N > M+1:
!      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
!
!  Restrictions:
!
!    -1 <= X <= 1
!     0 <= M <= N
!
!  Special values:
!
!    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
!    polynomial of the first kind equals the Legendre polynomial of the
!    first kind.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the maximum first index of the Legendre
!    polynomial, which must be at least 0.
!
!    Input, integer M, the second index of the Legendre polynomial,
!    which must be at least 0, and no greater than N.
!
!    Input, real X, the point at which the polynomial is to be
!    evaluated.  X must satisfy -1 <= X <= 1.
!
!    Output, real CX(0:N), the values of the first N+1 polynomials.
!
  integer n
!
  real cx(0:n)
  real fact
  integer i
  integer m
  real somx2
  real x
!
  if ( m < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LEGENDRE_PNM - Fatal error!'
    write ( *, * ) '  Input value of M is ', m
    write ( *, * ) '  but M must be nonnegative.'
    stop
  end if
 
  if ( m > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'LEGENDRE_PNM - Fatal error!'
    write ( *, * ) '  Input values of M, N = ', m, n
    write ( *, * ) '  but M must be less than or equal to N.'
    stop
  end if
 
  if ( x < -1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LEGENDRE_PNM - Fatal error!'
    write ( *, * ) '  Input value of X = ', x
    write ( *, * ) '  but X must be no less than -1.'
    stop
  end if
 
  if ( x > 1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LEGENDRE_PNM - Fatal error!'
    write ( *, * ) '  Input value of X = ', x
    write ( *, * ) '  but X must be no more than 1.'
    stop
  end if
  
  cx(0:m-1) = 0.0E+00

  cx(m) = 1.0E+00
  somx2 = sqrt ( 1.0E+00 - x**2 )
 
  fact = 1.0E+00
  do i = 1, m
    cx(m) = cx(m) * fact * somx2
    fact = fact + 2.0E+00
  end do
 
  cx(m+1) = x * real ( 2 * m + 1 ) * cx(m)

  do i = m+2, n
    cx(i) = ( real ( 2 * i - 1 ) * x * cx(i-1) &
      - real ( i + m - 1 ) * cx(i-2) ) / real ( i - m )
  end do

  return
end
subroutine legendre_qn ( n, x, cx )
!
!*******************************************************************************
!
!! LEGENDRE_QN evaluates the Legendre polynomials Qn(X).
!
!
!  Differential equation:
!
!    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0
!
!  First terms:
!
!    Q(0)(X) = 0.5 * LOG((1+X)/(1-X))
!    Q(1)(X) = Q(0)(X)*X - 1 
!    Q(2)(X) = Q(0)(X)*(3*X*X-1)/4 - 1.5*X
!    Q(3)(X) = Q(0)(X)*(5*X*X*X-3*X)/4 - 2.5*X**2 + 2/3
!    Q(4)(X) = Q(0)(X)*(35*X**4-30*X**2+3)/16 - 35/8 * X**3 + 55/24 * X
!    Q(5)(X) = Q(0)(X)*(63*X**5-70*X**3+15*X)/16 - 63/8*X**4 + 49/8*X**2 - 8/15
!
!  Recursion:
!
!    Q(0) = 0.5 * Log ( (1+X) / (1-X) )
!    Q(1) = 0.5 * X * Log ( (1+X) / (1-X) ) - 1.0
!
!    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
!
!  Restrictions:
!
!    -1 < X < 1
!
!  Special values:
!
!    Note that the Legendre polynomial Q(N)(X) is equal to the
!    associated Legendre function of the second kind,
!    Q(N,M)(X) with M = 0.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to evaluate.
!
!    Input, real X, the point at which the polynomials are to be
!    evaluated.  X must satisfy -1 < X < 1.
!
!    Output, real CX(0:N), the values of the first N+1 Legendre
!    polynomials at the point X.
!
  integer n
!
  real cx(0:n)
  integer i
  real x
!
!  Check the value of X.
!
  if ( x <= -1.0E+00 .or. x >= 1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LEGENDRE_QN - Fatal error!'
    write ( *, * ) '  Illegal input value of X = ', x
    write ( *, * ) '  But X must be between -1 and 1.'
    stop
  end if
 
  if ( n < 0 ) then
    return
  end if

  cx(0) = 0.5E+00 * log ( ( 1.0E+00 + x ) / ( 1.0E+00 - x ) )

  if ( n == 0 ) then
    return
  end if

  cx(1) = 0.5E+00 * x * log ( ( 1.0E+00 + x ) / ( 1.0E+00 - x ) ) - 1.0E+00

  do i = 2, n
    cx(i) = ( real ( 2 * i - 1 ) * cx(i-1) &
      - real ( i - 1 ) * cx(i-2) ) / real ( i )
  end do
 
  return
end
function log_factorial ( n )
!
!*******************************************************************************
!
!! LOG_FACTORIAL computes the natural logarithm of the factorial N!
!
!
!  Formula:
!
!    LOG ( FACTORIAL ( N ) ) 
!      = LOG ( PRODUCT ( 1 <= I <= N ) I )
!      = SUM ( ( 1 <= I <= N ) LOG ( I ) )
!
!  Modified:
!
!    17 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, LOG_FACTORIAL is returned as 0.
!
!    Output, real LOG_FACTORIAL, the logarithm of the factorial of N.
!
  integer i
  real log_factorial
  integer n
!
  log_factorial = 0.0E+00

  do i = 1, n
    log_factorial = log_factorial + log ( real ( i ) )
  end do

  return
end
subroutine pentagon_num ( n, p )
!
!*******************************************************************************
!
!! PENTAGON_NUM computes the N-th pentagonal number.
!
!
!  Definition:
!
!    The pentagonal number P(N) counts the number of dots in a figure of
!    N nested pentagons.  The pentagonal numbers are defined for both
!    positive and negative N.
!
!  First values:
!
!    N   P
!
!   -5  40
!   -4  26
!   -3  15
!   -2   7
!   -1   2
!    0   0
!    1   1
!    2   5
!    3  12
!    4  22
!    5  35
!
!  Formula:
!
!    P(N) = ( N * ( 3 * N - 1 ) ) / 2
!
!  Modified:
!
!    22 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the pentagonal number desired.
!
!    Output, integer P, the value of the N-th pentagonal number.
!
  integer n
  integer p
!
  p = ( n * ( 3 * n - 1 ) ) / 2

  return
end
function pi ( )
!
!*******************************************************************************
!
!! PI returns the value of pi.
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
!    Output, real PI, the value of pi.
!
  real pi
!
  pi = 3.14159265358979323846264338327950288419716939937510E+00

  return
end
function pyramid_num ( n )
!
!*******************************************************************************
!
!! PYRAMID_NUM returns the N-th pyramidal number.
!
!
!  Definition:
!
!    The N-th pyramidal number P(N) is formed by the sum of the first
!    N triangular numbers T(J):
!
!      T(J) = Sum ( J = 1 to N ) J
!
!      P(N) = Sum ( I = 1 to N ) T(I)
!
!    By convention, T(0) = 0.
!
!  Formula:
!
!    P(N) = ( (N+1)**3 - (N+1) ) / 6
!
!  First Values:
!
!      0
!      1
!      4
!     10
!     20
!     35
!     56
!     84
!    120
!    165
!
!  Modified:
!
!    11 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the desired number, which must be
!    at least 0.
!
!    Output, integer PYRAMID_NUM, the N-th pyramidal number.
!
  integer n
  integer pyramid_num
!
  pyramid_num = ( ( n + 1 )**3 - ( n + 1 ) ) / 6

  return
end
subroutine rpoly_val_horner ( n, c, x, cx )
!
!*******************************************************************************
!
!! RPOLY_VAL_HORNER evaluates a polynomial using Horner's method.
!
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the polynomial.
!
!    Input, real C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X**I.
!
!    Input, real X, the point at which the polynomial is to be evaluated.
!
!    Output, real PVAL, the value of the polynomial at X.
!
  integer n
!
  real c(0:n)
  real cx
  integer i
  real x
!
  cx = c(n)
  do i = n - 1, 0, -1
    cx = cx * x + c(i)
  end do

  return
end
subroutine stirling1 ( m, n, s1 )
!
!*******************************************************************************
!
!! STIRLING1 computes the Stirling numbers of the first kind.
!
!
!  Discussion:
!
!    The absolute value of the Stirling number S1(M,N) gives the number
!    of permutations on M objects having exactly N cycles, while the
!    sign of the Stirling number records the sign (odd or even) of
!    the permutations.  For example, there are six permutations on 3 objects:
!
!      A B C   3 cycles (A) (B) (C)
!      A C B   2 cycles (A) (BC)
!      B A C   2 cycles (AB) (C)
!      B C A   1 cycle  (ABC)
!      C A B   1 cycle  (ABC)
!      C B A   2 cycles (AC) (B)
!
!    There are 
!
!      2 permutations with 1 cycle, and S1(3,1) = 2
!      3 permutations with 2 cycles, and S1(3,2) = -3,
!      1 permutation with 3 cycles, and S1(3,3) = 1.
!
!    Since there are M! permutations of M objects, the sum of the absolute 
!    values of the Stirling numbers in a given row, 
!
!      Sum ( I = 1 to M ) ABS ( S1(M,I) ) = M!
!
!  First terms:
!
!    M/N:  1     2      3     4     5    6    7    8
!
!    1     1     0      0     0     0    0    0    0
!    2    -1     1      0     0     0    0    0    0
!    3     2    -3      1     0     0    0    0    0
!    4    -6    11     -6     1     0    0    0    0
!    5    24   -50     35   -10     1    0    0    0
!    6  -120   274   -225    85   -15    1    0    0
!    7   720 -1764   1624  -735   175  -21    1    0
!    8 -5040 13068 -13132  6769 -1960  322  -28    1
!
!  Recursion:
!
!    S1(M,1) = (-1)**(M-1) * (M-1)! for all M.
!    S1(I,I) = 1 for all I.
!    S1(I,J) = 0 if J>I.
!
!    S1(M,N) = S1(M-1,N-1) - (M-1) * S1(M-1,N)
!
!  Properties:
!
!    Sum ( K = 1 to N ) S2(I,K) * S1(K,J) = Delta(I,J)
!    X_N = Sum ( K = 0 to N ) S1(N,K) X**K
!    where X_N is the falling factorial function.
!
!  Modified:
!
!    24 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the table.
!
!    Input, integer N, the number of columns of the table.
!
!    Output, integer S1(M,N), the Stirling numbers of the first kind.
!
  integer m
  integer n
!
  integer i
  integer j
  integer s1(m,n)
!
  if ( m <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  s1(1,1) = 1
  s1(1,2:n) = 0
 
  do i = 2, m

    s1(i,1) = - ( i - 1 ) * s1(i-1,1)

    do j = 2, n
      s1(i,j) = s1(i-1,j-1) - ( i - 1 ) * s1(i-1,j)
    end do

  end do
 
  return
end
subroutine stirling2 ( m, n, s2 )
!
!*******************************************************************************
!
!! STIRLING2 computes the Stirling numbers of the second kind.
!
!
!  Discussion:
!
!    S2(M,N) represents the number of distinct partitions of M elements
!    into N nonempty sets.  For a fixed M, the sum of the Stirling
!    numbers S2(M,N) is represented by B(M), called "Bell's number",
!    and represents the number of distinct partitions of M elements.
!
!    For example, with 4 objects, there are:
!
!    1 partition into 1 set:
!
!      (A,B,C,D)
!
!    7 partitions into 2 sets:
!
!      (A,B,C) (D)
!      (A,B,D) (C)
!      (A,C,D) (B)
!      (A) (B,C,D)
!      (A,B) (C,D)
!      (A,C) (B,D)
!      (A,D) (B,C)
!
!    6 partitions into 3 sets:
!
!      (A,B) (C) (D)
!      (A) (B,C) (D)
!      (A) (B) (C,D)
!      (A,C) (B) (D)
!      (A,D) (B) (C)
!      (A) (B,D) (C)
!
!    1 partition into 4 sets:
!
!      (A) (B) (C) (D)
!
!    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
!
!
!  First terms:
!
!    M/N: 1    2    3    4    5    6    7    8
!
!    1    1    0    0    0    0    0    0    0
!    2    1    1    0    0    0    0    0    0
!    3    1    3    1    0    0    0    0    0
!    4    1    7    6    1    0    0    0    0
!    5    1   15   25   10    1    0    0    0
!    6    1   31   90   65   15    1    0    0
!    7    1   63  301  350  140   21    1    0
!    8    1  127  966 1701 1050  266   28    1
!
!  Recursion:
!
!    S2(M,1) = 1 for all M.
!    S2(I,I) = 1 for all I.
!    S2(I,J) = 0 if J>I.
!
!    S2(M,N) = N * S2(M-1,N) + S2(M-1,N-1)
!
!  Properties:
!
!    Sum ( K = 1 to N ) S2(I,K) * S1(K,J) = Delta(I,J)
!    X**N = Sum ( K = 0 to N ) S2(N,K) X_K
!    where X_K is the falling factorial function.
!
!  Modified:
!
!    24 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the table.
!
!    Input, integer N, the number of columns of the table.
!
!    Output, integer S2(M,N), the Stirling numbers of the second kind.
!
  integer m
  integer n
!
  integer i
  integer j
  integer s2(m,n)
!
  if ( m <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  s2(1,1) = 1
  s2(1,2:n) = 0
 
  do i = 2, m

    s2(i,1) = 1

    do j = 2, n
      s2(i,j) = j * s2(i-1,j) + s2(i-1,j-1)
    end do

  end do
 
  return
end
function triangle_num ( n )
!
!*******************************************************************************
!
!! TRIANGLE_NUM returns the N-th triangular number.
!
!
!  Definition:
!
!    The N-th triangular number T(N) is formed by the sum of the first
!    N integers:
!
!      T(N) = Sum ( I = 1 to N ) I
!
!    By convention, T(0) = 0.
!
!  Formula:
!
!    T(N) = ( N * ( N + 1 ) ) / 2
!
!  First Values:
!
!     0
!     1
!     3
!     6
!    10
!    15
!    21
!    28
!    36
!    45
!    55
!
!  Modified:
!
!    11 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the desired number, which must be
!    at least 0.
!
!    Output, integer TRIANGLE_NUM, the N-th triangular number.
!
  integer n
  integer triangle_num
!
  triangle_num = ( n * ( n + 1 ) ) / 2

  return
end
subroutine vibonacci ( n, v )
!
!*******************************************************************************
!
!! VIBONACCI computes the first N Vibonacci numbers.
!
!
!  Discussion:
!
!    The "Vibonacci numbers" are a generalization of the Fibonacci numbers:
!      V(N+1) = +/- V(N) +/- V(N-1)
!    where the signs are chosen randomly.
!
!  Reference:
!
!    Brian Hayes,
!    The Vibonacci Numbers,
!    American Scientist,
!    July-August 1999, Volume 87, Number 4.
!
!    Divakar Viswanath,
!    Random Fibonacci sequences and the number 1.13198824,
!    Mathematics of Computation, 1998.
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
!    Input, integer N, the highest number to compute.
!
!    Output, integer V(N), the first N Vibonacci numbers.  By convention,
!    V(1) and V(2) are taken to be 1.
!
  integer n
!
  integer i
  integer j
  integer s1
  integer s2
  integer v(n)
!
  if ( n <= 0 ) then
    return
  end if

  v(1) = 1

  if ( n <= 1 ) then
    return
  end if

  v(2) = 1

  do i = 3, n
    
    call i_random ( 0, 1, j )

    if ( j == 0 ) then
      s1 = -1
    else
      s1 = +1
    end if

    call i_random ( 0, 1, j )

    if ( j == 0 ) then
      s2 = -1
    else
      s2 = +1
    end if

    v(i) = s1 * v(i-1) + s2 * v(i-2)

  end do
 
  return
end
subroutine zeckendorf ( n, m, i_list, f_list )
!
!*******************************************************************************
!
!! ZECKENDORF produces the Zeckendorf decomposition of a positive integer.
!
!
!  Discussion:
!
!    Zeckendorf proved that every positive integer can be represented
!    uniquely as the sum of non-consecutive Fibonacci numbers.
!
!    N = Sum ( 1 <= I <= M ) F_LIST(I)
!
!  Example:
!
!     N    Decomposition
!
!    50    34 + 13 + 3
!    51    34 + 13 + 3 + 1
!    52    34 + 13 + 5
!    53    34 + 13 + 5 + 1
!    54    34 + 13 + 5 + 2
!    55    55
!    56    55 + 1
!    57    55 + 2
!    58    55 + 3
!    59    55 + 3 + 1
!    60    55 + 5
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the positive integer to be decomposed.
!
!    Output, integer M, the number of parts in the decomposition.
!
!    Output, integer I_LIST(M), the index of the Fibonacci numbers
!    in the decomposition.
!
!    Output, integer F_LIST(M), the value of the Fibonacci numbers 
!    in the decomposition.
!
  integer f
  integer f_list(*)
  integer i
  integer i_list(*)
  integer m
  integer m_max
  integer n
  integer n_copy
!
!  Determine how much space we have in F_LIST and I_LIST.
!
!  ...Well, THAT didn't work!
!
! m_max = min ( size ( f_list ), size ( i_list ) )
!
!  FIX THIS SOMEDAY SOON!
!
  m_max = 100

  m = 0

  if ( m >= m_max ) then
    return
  end if

  n_copy = n
!
!  Extract a sequence of Fibonacci numbers.
!
  do while ( n_copy > 0 .and. m < m_max ) 
    call fibonacci_floor ( n_copy, f, i )
    m = m + 1
    i_list(m) = i
    n_copy = n_copy - f
  end do
!
!  Replace any pair of consecutive indices ( I, I-1 ) by I+1.
!
  do i = m, 2, -1

    if ( i_list(i-1) == i_list(i) + 1 ) then
      i_list(i-1) = i_list(i-1) + 1
      i_list(i:m-1) = i_list(i+1:m)
      i_list(m) = 0
      m = m - 1
    end if

  end do
!
!  Fill in the actual values of the Fibonacci numbers.
!
  do i = 1, m
    call fibonacci_direct ( i_list(i), f_list(i) )
  end do

  return
end
