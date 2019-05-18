function acosh2 ( x )
!
!*******************************************************************************
!
!! ACOSH2 returns the inverse hyperbolic cosine of a number.
!
!
!  Definition:
!
!    Applying the inverse function
!
!      Y = ACOSH2(X)
!
!    implies that
!
!      X = COSH(Y) = 0.5 * ( EXP(Y) + EXP(-Y) ).
!
!    For every X greater than or equal to 1, there are two possible
!    choices Y such that X = COSH(Y), differing only in sign.  It
!    is usual to resolve this choice by taking the value of ACOSH2(X)
!    to be nonnegative.
!
!  Discussion:
!
!    Since a library function ACOSH may be available on some systems,
!    this routine is named ACOSH2 to avoid name conflicts.
!
!  Method:
!
!    One formula is:
!
!      ACOSH2 = LOG ( X + SQRT ( X**2 - 1.0 ) )
!
!    but this formula suffers from roundoff and overflow problems.
!    The formula used here was recommended by W Kahan, as discussed
!    by Moler.
!
!  Reference:
!
!    Cleve Moler,
!    Trigonometry is a Complex Subject,
!    MATLAB News and Notes,
!    Summer 1998.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose inverse hyperbolic cosine is desired.
!    X should be greater than or equal to 1.
!
!    Output, real ACOSH2, the inverse hyperbolic cosine of X.  The
!    principal value (that is, the positive value of the two ) is returned.
!
  real acosh2
  real x
!
  if ( x < 1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ACOSH2 - Fatal error!'
    write ( *, * ) '  Argument X must be >= 1.'
    write ( *, * ) '  The input X = ', x
    stop
  end if

  acosh2 = 2.0E+00 * log ( &
    sqrt ( 0.5E+00 * ( x + 1.0E+00 ) ) + sqrt ( 0.5E+00 * ( x - 1.0E+00 ) ) )

  return
end
function agud ( gamma )
!
!*******************************************************************************
!
!! AGUD evaluates the inverse Gudermannian function.
!
!
!  Definition:
!
!    The Gudermannian function relates the hyperbolic and trigonomentric
!    functions.  For any argument X, there is a corresponding value
!    GAMMA so that
!
!      SINH(X) = TAN(GAMMA).
!
!    This value GAMMA(X) is called the Gudermannian of X.  The inverse
!    Gudermannian function is given as input a value GAMMA and computes
!    the corresponding value X.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real GAMMA, the value of the Gudermannian.
!
!    Output, real AGUD, the argument of the Gudermannian.
!
  real agud
  real gamma
  real pi
!
  agud = log ( tan ( 0.25E+00 * pi ( ) + 0.5E+00 * gamma ) )

  return
end
function arc_cosine ( c )
!
!*******************************************************************************
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
!    This routine truncates arguments outside the range.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real C, the argument.
!
!    Output, real ARC_COSINE, an angle whose cosine is C.
!
  real arc_cosine
  real c
  real c2
!
  c2 = c
  c2 = max ( c2, -1.0E+00 )
  c2 = min ( c2, +1.0E+00 )

  arc_cosine = acos ( c2 )

  return
end
function asinh2 ( x )
!
!*******************************************************************************
!
!! ASINH2 returns the inverse hyperbolic sine of a number.
!
!
!  Definition:
!
!    Y = ASINH2(X) implies that
!    X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).
!
!  Discussion:
!
!    Since a library function ASINH may be available on some systems,
!    this routine is named ASINH2 to avoid name conflicts.
!
!  Modified:
!
!    13 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose inverse hyperbolic sine is desired.
!
!    Output, real ASINH2, the inverse hyperbolic sine of X.
!
  real asinh2
  real x
!
  asinh2 = log ( x + sqrt ( x * x + 1.0E+00 ) )

  return
end
function atan4 ( y, x )
!
!*******************************************************************************
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
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
!    Input, real Y, X, two quantities which represent the tangent of
!    an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ATAN4, an angle between 0 and 2 * PI, whose tangent is
!    (Y/X), and which lies in the appropriate quadrant so that the signs
!    of its cosine and sine match those of X and Y.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
!
  real abs_x
  real abs_y
  real atan4
  real theta
  real theta_0
  real x
  real y
!
!  Special cases:
!
  if ( x == 0.0E+00 ) then

    if ( y > 0.0E+00 ) then
      theta = PI / 2.0E+00
    else if ( y < 0.0E+00 ) then
      theta = 3.0E+00 * PI / 2.0E+00
    else if ( y == 0.0E+00 ) then
      theta = 0.0E+00
    end if

  else if ( y == 0.0E+00 ) then

    if ( x > 0.0E+00 ) then
      theta = 0.0E+00
    else if ( x < 0.0E+00 ) then
      theta = PI
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( x > 0.0E+00 .and. y > 0.0E+00 ) then
      theta = theta_0
    else if ( x < 0.0E+00 .and. y > 0.0E+00 ) then
      theta = PI - theta_0
    else if ( x < 0.0E+00 .and. y < 0.0E+00 ) then
      theta = PI + theta_0
    else if ( x > 0.0E+00 .and. y < 0.0E+00 ) then
      theta = 2.0E+00 * PI - theta_0
    end if

  end if

  atan4 = theta

  return
end
function atanh2 ( x )
!
!*******************************************************************************
!
!! ATANH2 returns the inverse hyperbolic tangent of a number.
!
!
!  Definition:
!
!    Y = ATANH2(X) implies that
!    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )
!
!  Discussion:
!
!    Since a library function ATANH may be available on some systems,
!    this routine is named ATANH2 to avoid name conflicts.
!
!  Modified:
!
!    13 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose inverse hyperbolic tangent is desired.
!    The absolute value of X should be less than or equal to 1.
!
!    Output, real ATANH2, the inverse hyperbolic tangent of X.
!
  real atanh2
  real x
!
  if ( abs ( x ) >= 1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ATANH2 - Fatal error!'
    write ( *, * ) '  ABS(X) must be < 1.'
    write ( *, * ) '  Your input is X = ', x
    stop
  end if

  atanh2 = 0.5E+00 * log ( ( 1.0E+00 + x ) / ( 1.0E+00 - x ) )

  return
end
subroutine axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )
!
!*******************************************************************************
!
!! AXIS_LIMITS returns "nice" axis limits for a plot.
!
!
!  Discussion:
!
!    The routine is given information about the range of a variable, and
!    the number of divisions desired.  It returns suggestions for
!    labeling a plotting axis for the variable, including the
!    starting and ending points, the length of a single division,
!    and a suggested tick marking for the axis.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XMIN, XMAX, the lower and upper values that must be
!    included on the axis.  XMIN must be less than XMAX.
!
!    Input, integer NDIVS, the number of divisions desired along
!    the axis.
!
!    Output, real PXMIN, PXMAX, the recommended lower and upper axis
!    bounds.  It will be the case that PXMIN <= XMIN < XMAX <= PXMAX.
!
!    Output, real PXDIV, the recommended size of a single division.
!
!    Output, integer NTICKS, a suggested number of ticks to use,
!    if subdividing each of the NDIVS divisions of the axis.
!
  integer, parameter :: nsteps = 5
!
  real best
  real good
  integer i
  integer i_log_10
  integer ihi
  integer ilo
  integer intlog
  integer iticks(5)
  integer ival
  integer j
  integer ndivs
  integer nticks
  real pxmax
  real pxmax2
  real pxmin
  real pxmin2
  real pxdiv
  real pxdiv2
  real reldif
  real steps(nsteps)
  real temp
  real xmax
  real xmin
!
  if ( xmin == xmax ) then
    xmin = xmin - 0.5E+00
    xmax = xmax + 0.5E+00
  else if ( xmin > xmax ) then
    temp = xmin
    xmin = xmax
    xmax = temp
  end if

  if ( ndivs <= 0 ) then
    ndivs = 5
  end if
!
  steps(1) =  1.0E+00
  steps(2) =  2.0E+00
  steps(3) =  4.0E+00
  steps(4) =  5.0E+00
  steps(5) = 10.0E+00

  iticks(1) = 5
  iticks(2) = 4
  iticks(3) = 4
  iticks(4) = 5
  iticks(5) = 5
!
!  Set RELDIF, the size of the X interval divided by the largest X.
!
  if ( xmax /= xmin ) then
    reldif = ( xmax - xmin ) / max ( abs ( xmax ), abs ( xmin ) )
  else
    reldif = 0.0E+00
  end if
!
!  If RELDIF tells us that XMIN and XMAX are extremely close,
!  do some simple things.
!
  if ( reldif < 0.00001E+00 ) then

    if ( xmax == 0.0E+00 ) then

      pxdiv = 1.0E+00

    else

      intlog = i_log_10 ( xmax )

      if ( intlog < 0 ) then
        intlog = intlog - 1
      end if

      pxdiv = 10.0E+00**intlog

      if ( pxdiv > 1.0E+00 ) then
        pxdiv = 1.0E+00
      end if

    end if

    nticks = 5
    pxmin = xmax - real ( ndivs / 2 ) * pxdiv
    pxmax = xmax + real ( ndivs - ( ndivs / 2 ) ) * pxdiv
!
!  But now handle the more general case, when XMIN and XMAX
!  are relatively far apart.
!
  else

    best = - 999.0E+00
!
!  On second loop, increase INTLOG by 1.
!
    do j = 1, 2
!
!  Compute INTLOG, roughly the logarithm base 10 of the range
!  divided by the number of divisions.
!
      intlog = i_log_10 ( ( xmax - xmin ) / real ( ndivs ) ) + ( j - 1 )

      if ( xmax - xmin  < real ( ndivs ) ) then
        intlog = intlog - 1
      end if
!
!  Now consider taking 1, 2, 4, 5 or 10 steps of size 10**INTLOG:
!
      do i = 1, nsteps
!
!  Compute the size of each step.
!
        pxdiv2 = steps(i) * 10.0E+00**intlog
!
!  Make sure NDIVS steps can reach from XMIN to XMAX, at least.
!
        if ( xmin + ndivs * pxdiv2 >= xmax ) then
!
!  Now decide where to start the axis.
!  Start the axis at PXMIN2, to the left of XMIN, and
!  representing a whole number of steps of size PXDIV2.
!
          if ( xmin >= 0.0E+00 ) then
            ival = int ( xmin / pxdiv2 )
          else
            ival = int ( xmin / pxdiv2 ) - 1
          end if

          pxmin2 = ival * pxdiv2
!
!  PXMAX2 is, of course, NDIVS steps above PXMIN2.
!
          pxmax2 = pxmin2 + ndivs * pxdiv2
!
!  Only consider going on if PXMAX2 is at least XMAX.
!
          if ( pxmax2 >= xmax ) then
!
!  Now judge this grid by the relative amount of wasted axis length.
!
            good = ( xmax - xmin ) / ( pxmax2 - pxmin2 )

            if ( good > best ) then
              best = good
              pxmax = pxmax2
              pxmin = pxmin2
              pxdiv = pxdiv2
              nticks = iticks(i)
            end if

          end if

        end if

      end do

    end do

  end if
!
!  If necessary, adjust the locations of PXMIN and PXMAX so that the
!  interval is more symmetric in containing XMIN through XMAX.
!
  do

    ilo = int ( xmin - pxmin ) / pxdiv
    ihi = int ( pxmax - xmax ) / pxdiv

    if ( ihi < ilo + 2 ) then
      exit
    end if

    pxmin = pxmin - pxdiv
    pxmax = pxmax - pxdiv

  end do

  return
end
subroutine bar_check ( digit )
!
!*******************************************************************************
!
!! BAR_CHECK computes the check digit for a barcode.
!
!
!  Formula:
!
!    CHECK1 = SUM ( I = 1, 11, by 2's ) DIGIT(I)
!       + 3 * SUM ( I = 2, 10, by 2'2 ) DIGIT(I)
!
!    CHECK = MOD ( 10 - MOD ( CHECK1, 10 ), 10 )
!
!  Modified:
!
!    19 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer DIGIT(12).
!
!    On input, entries 1 through 11 of DIGIT contain the digits of the
!    bar code.  Each entry must be between 0 and 9.
!
!    On output, entry 12 of DIGIT contains the check digit.
!
  integer digit(12)
!
  digit(12) = sum ( digit(1:11:2) ) + 3 * sum ( digit(2:10:2) )

  digit(12) = mod ( 10 - mod ( digit(12), 10 ), 10 )

  return
end
subroutine bar_code ( digit, bar )
!
!*******************************************************************************
!
!! BAR_CODE constructs the 113 character barcode from 11 digits.
!
!
!  Modified:
!
!    10 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer DIGIT(12).
!
!    On input, the first 11 entries of DIGIT contain a code to be
!    turned into a barcode.
!
!    On output, the 12-th entry of DIGIT contains a check digit.
!
!    Output, character( len = 113 ) BAR, the bar code corresponding to the
!    digit information.
!
  character ( len = 113 ) bar
  character ( len = 7 ) codel
  character ( len = 7 ) coder
  integer digit(12)
  integer i
!
!  9 character quiet zone.
!
  bar(1:9) = '000000000'
!
!  3 character guard pattern.
!
  bar(10:12) = '101'
!
!  7 character product category.
!
  call bar_digit_code ( digit(1), codel, coder )
  bar(13:19) = codel
!
!  35 characters contain the 5 digit manufacturer code.
!
  do i = 1, 5
    call bar_digit_code ( digit(i+1), codel, coder )
    bar(20+(i-1)*7:20+(i-1)*7+6) = codel
  end do
!
!  Center guard pattern.
!
  bar(55:59) = '01010'
!
!  35 characters contain the 5 digit product code.
!
  do i = 1, 5
    call bar_digit_code ( digit(i+6), codel, coder )
    bar(60+(i-1)*7:60+(i-1)*7+6) = coder
  end do
!
!  Check digit.
!
  call bar_check ( digit )
  call bar_digit_code ( digit(12), codel, coder )
  bar(95:101) = coder
!
!  Guard pattern.
!
  bar(102:104) = '101'
!
!  Quiet zone.
!
  bar(105:113) = '000000000'

  return
end
subroutine bar_digit_code ( digit, codel, coder )
!
!*******************************************************************************
!
!! BAR_DIGIT_CODE returns the 7 character right and left bar codes for a digit.
!
!
!  Example:
!
!    DIGIT = 3
!    CODEL = '0111101'
!    CODER = '1000010'
!
!  Modified:
!
!    10 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit, between 0 and 9.
!
!    Output, character ( len = 7 ) CODEL, CODER; the left and right
!    codes for the digit.
!
  character ( len = 7 ) codel
  character ( len = 7 ) coder
  integer digit
  integer i
!
  if ( digit == 0 ) then
    codel = '0001101'
  else if ( digit == 1 ) then
    codel = '0011001'
  else if ( digit == 2 ) then
    codel = '0010011'
  else if ( digit == 3 ) then
    codel = '0111101'
  else if ( digit == 4 ) then
    codel = '0100011'
  else if ( digit == 5 ) then
    codel = '0110001'
  else if ( digit == 6 ) then
    codel = '0101111'
  else if ( digit == 7 ) then
    codel = '0111011'
  else if ( digit == 8 ) then
    codel = '0110111'
  else if ( digit == 9 ) then
    codel = '0001011'
  else
    codel = '???????'
  end if

  do i = 1, 7
    if ( codel(i:i) == '0' ) then
      coder(i:i) = '1'
    else if ( codel(i:i) == '1' ) then
      coder(i:i) = '0'
    else
      coder(i:i) = '?'
    end if
  end do

  return
end
subroutine bin_search_one_2d ( bin, nset, pset, nbin, bin_start, bin_next, &
  ptest, found_a_neighbor, i_min, d_min_sq, compares )
!
!*******************************************************************************
!
!! BIN_SEARCH_ONE_2D searches one cell in a 2D array of bins.
!
!
!  Discussion:
!
!    The bins are presumed to have been set up by successive calls to:
!
!      R2VEC_BIN_EVEN2,
!      R2VEC_BINNED_REORDER, and
!      R2VEC_BINNED_SORT_A.
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
!    Input, integer BIN(2), the indices of the cell to be examined.
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN(2), the number of cells in the horizontal and
!    vertical directions.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, real PTEST(2), the coordinates of the test point.
!
!    Input/output, logical FOUND_A_NEIGHBOR, is set to TRUE if at least
!    one point of PSET is found in the current bin.  Otherwise, it retains its
!    input value.
!
!    Input/output, integer I_MIN, the index of the nearest neighbor in
!    PSET to PTEST, if at least one neighbor has been found.
!
!    Input/output, real D_MIN_SQ, the square of the distance from the nearest
!    neighbor in PSET to PTEST, if at least one neighbor has been found.
!
!    Input/output, integer COMPARES, the number of elements of PSET whose
!    distance to PTEST has been computed.
!
  integer, parameter :: ndim = 2
!
  integer nbin(ndim)
  integer nset
!
  integer bin(ndim)
  integer bin_next(nset)
  integer bin_start(nbin(1),nbin(2))
  integer compares
  real d_min_sq
  real d_sq
  logical found_a_neighbor
  integer i_min
  integer node
  real pset(ndim,nset)
  real ptest(ndim)
!
  node = bin_start(bin(1),bin(2))

  do while ( node > 0 )

    found_a_neighbor = .true.

    d_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )
    compares = compares + 1

    if ( d_sq < d_min_sq ) then
      d_min_sq = d_sq
      i_min = node
    end if

    node = bin_next(node)

  end do

  return
end
subroutine bin_to_r2_even ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R2_EVEN returns the limits for a given R2 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
!    An initial bin takes everything less than A, and a final bin takes
!    everything greater than B.
!
!  Example:
!
!    NBIN = 7, A(1) = 5, B(1) = 15
!              A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1     -HUGE -HUGE   5     0
!    2, 2       5     0     7     4
!    3, 3       7     4     9     8
!    4, 4       9     8    11    12
!    5, 5      11    12    13    16
!    6, 6      13    16    15    20
!    7, 7      15    20    HUGE HUGE
!
!  Modified:
!
!    02 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.  NBIN is normally at least 3.
!    If NBIN is 1 or 2, then everything is assigned to bin 1.
!
!    Input, integer BIN(2), the index of the bin to be considered.
!    If BIN(I) is less than 1, or greater than NBIN, the user will get what
!    the user deserves.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(2), CMAX(2), the minimum and maximum limits on the bin.
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call bin_to_r_even ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r2_even2 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R2_EVEN2 returns the limits for a given R2 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A(1) = 5, B(1) = 15
!              A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1       5     0     7     4
!    2, 2       7     4     9     8
!    3, 3       9     8    11    12
!    4, 4      11    12    13    16
!    5, 5      13    16    15    20
!
!  Modified:
!
!    07 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, integer BIN(2), the index of the bin to be considered.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(2), CMAX(2), the minimum and maximum limits on the bin.
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call bin_to_r_even2 ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r2_even3 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R2_EVEN3 returns the limits for a given R2 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A(I) to B(I) is divided into NBIN(I) equal
!    subintervals or bins.
!
!  Example:
!
!    NBIN = (/ 4, 5, /)
!
!    A(1) = 5, B(1) = 15
!    A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1       5     0     7     4
!    2, 2       7     4     9     8
!    3, 3       9     8    11    12
!    4, 4      11    12    13    16
!    5, 5      13    16    15    20
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!
!    Input, integer BIN(2), the index of the bin to be considered.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(2), CMAX(2), the minimum and maximum limits on the bin.
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin(ndim)
!
  do i = 1, ndim
    call bin_to_r_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r3_even2 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R3_EVEN2 returns the limits for a given R3 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
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
!    Input, integer NBIN, the number of bins.
!
!    Input, integer BIN(3), the index of the bin to be considered.
!
!    Input, real A(3), B(3), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(3), CMAX(3), the minimum and maximum limits on the bin.
!
  integer, parameter :: ndim = 3
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call bin_to_r_even2 ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r_even ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R_EVEN returns the limits for a given "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
!    An initial bin takes everything less than A, and a final bin takes
!    everything greater than B.
!
!  Example:
!
!    NBIN = 7, A = 10, B = 20
!
!    BIN      CMIN  CMAX
!
!    1         -HUGE 10
!    2         10    12
!    3         12    14
!    4         14    16
!    5         16    18
!    6         18    20
!    7         20    HUGE
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.  NBIN is normally at least 3.
!    If NBIN is 1 or 2, then everything is assigned to bin 1.
!
!    Input, integer BIN, the index of the bin to be considered.
!    If BIN is less than 1, or greater than NBIN, the user will get what
!    the user deserves.
!
!    Input, real A, B, the lower and upper limits of the bin interval.
!    While A is expected to be less than B, the code should return useful
!    results if A is actually greater than B.
!
!    Output, real CMIN, CMAX, the minimum and maximum limits on the bin.
!
  real a
  real b
  integer bin
  real cmax
  real cmin
  integer nbin
!
!  Take care of special cases.
!
  if ( nbin <= 2 ) then
    cmin = - huge ( cmin )
    cmax =   huge ( cmax )
    return
  end if

  if ( b == a ) then
    cmin = - huge ( cmin )
    cmax =   huge ( cmax )
    return
  end if
!
!  Compute the bin limits.
!
  if ( bin == 1 ) then
    cmin = - huge ( cmin )
    cmax = a
  else if ( bin < nbin ) then
    cmin = ( real ( nbin - bin ) * a + real ( bin - 2 ) * b ) &
      / real ( nbin - 2 )
    cmax = ( real ( nbin - bin - 1 ) * a + real ( bin - 1 ) * b ) &
      / real ( nbin - 2 )
  else if ( bin == nbin ) then
    cmin = b
    cmax = huge ( cmax )
  else
    cmin = - huge ( cmin )
    cmax =   huge ( cmax )
  end if

  return
end
subroutine bin_to_r_even2 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R_EVEN2 returns the limits for a given "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A = 10, B = 20
!
!    BIN      CMIN  CMAX
!
!    1         10    12
!    2         12    14
!    3         14    16
!    4         16    18
!    5         18    20
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, integer BIN, the index of the bin to be considered.
!    If BIN is less than 1, or greater than NBIN, the user will get what
!    the user deserves.
!
!    Input, real A, B, the lower and upper limits of the bin interval.
!    While A is expected to be less than B, the code should return useful
!    results if A is actually greater than B.
!
!    Output, real CMIN, CMAX, the minimum and maximum limits on the bin.
!
  real a
  real b
  integer bin
  real cmax
  real cmin
  integer nbin
!
!  Compute the bin limits.
!
  if ( bin < 1 ) then
    cmin = - huge ( cmin )
    cmax = a
  else if ( bin <= nbin ) then
    cmin = ( real ( nbin - bin + 1 ) * a + real ( bin - 1 ) * b ) &
      / real ( nbin )
    cmax = ( real ( nbin - bin ) * a + real ( bin ) * b ) &
      / real ( nbin )
  else if ( bin > nbin ) then
    cmin = b
    cmax = huge ( cmax )
  end if

  return
end
function bmi_english ( w_lb, h_ft, h_in )
!
!*******************************************************************************
!
!! BMI_ENGLISH computes the body mass index given English measurements.
!
!
!  Modified:
!
!    20 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real W_LB, the body weight in pounds.
!
!    Input, real H_FT, H_IN, the body height in feet and inches
!
!    Output, real BMI_ENGLISH, the body mass index.
!
  real bmi_english
  real bmi_metric
  real feet_to_meters
  real h_ft
  real h_in
  real h_m
  real pounds_to_kilograms
  real w_kg
  real w_lb
!
  w_kg = pounds_to_kilograms ( w_lb )

  h_m = feet_to_meters ( h_ft + ( h_in / 12.0E+00 ) )

  bmi_english = bmi_metric ( w_kg, h_m )

  return
end
function bmi_metric ( w_kg, h_m )
!
!*******************************************************************************
!
!! BMI_METRIC computes the body mass index given metric measurements.
!
!
!  Modified:
!
!    20 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real W_KG, the body weight in kilograms.
!
!    Input, real H_M, the body height in meters.
!
!    Output, real BMI_METRIC, the body mass index.
!
  real bmi_metric
  real h_m
  real w_kg
!
  bmi_metric = w_kg / h_m**2

  return
end
function c_cube_root ( x )
!
!*******************************************************************************
!
!! C_CUBE_ROOT returns the principal cube root of a complex number.
!
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex X, the number whose cube root is desired.
!
!    Output, complex C_CUBE_ROOT, the cube root of X.
!
  real a
  real atan4
  real b
  complex c_cube_root
  real mag
  real theta
  complex x
!
  a = real ( x )
  b = imag ( x )
  mag = sqrt ( a * a + b * b )

  if ( mag == 0.0E+00 ) then

    c_cube_root = cmplx ( 0.0E+00, 0.0E+00 )

  else

    theta = atan4 ( b, a )
    c_cube_root = mag**( 1.0E+00 / 3.0E+00 ) &
      * cmplx ( cos ( theta / 3.0E+00 ), sin ( theta / 3.0E+00 ) )

  end if

  return
end
function c_le_l1 ( x, y )
!
!*******************************************************************************
!
!! C_LE_L1 := X <= Y for complex values, and the L1 norm.
!
!
!  Definition:
!
!    The L1 norm can be defined here as:
!
!      C_NORM1(X) = abs ( real (X) ) + abs ( imag (X) )
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
!    Output, logical C_LE_L1, is TRUE if X <= Y.
!
  logical c_le_l1
  complex x
  complex y
!
  if ( abs ( real ( x ) ) + abs ( imag ( x ) ) <= &
       abs ( real ( y ) ) + abs ( imag ( y ) ) ) then
    c_le_l1 = .true.
  else
    c_le_l1 = .false.
  end if

  return
end
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
function c_le_linf ( x, y )
!
!*******************************************************************************
!
!! C_LE_LINF := X <= Y for complex values, and the L Infinity norm.
!
!
!  Definition:
!
!    The L Infinity norm can be defined here as:
!
!      C_NORM_Infinity(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
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
!    Output, logical C_LE_LINF, is TRUE if X <= Y.
!
  logical c_le_linf
  complex x
  complex y
!
  if ( max ( abs ( real ( x ) ), abs ( imag ( x ) ) ) <= &
       max ( abs ( real ( y ) ), abs ( imag ( y ) ) ) ) then
    c_le_linf = .true.
  else
    c_le_linf = .false.
  end if

  return
end
function c_norm1 ( x )
!
!*******************************************************************************
!
!! C_NORM1 evaluates the L1 norm of a complex number.
!
!
!  Discussion:
!
!    Numbers of equal norm lie along diamonds centered at (0,0).
!
!  Definition:
!
!    The L1 norm can be defined here as:
!
!      C_NORM1(X) = abs ( real (X) ) + abs ( imag (X) )
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
!    Input, complex X, the value whose norm is desired.
!
!    Output, real C_NORM1, the norm of X.
!
  real c_norm1
  complex x
!
  c_norm1 = abs ( real ( x ) ) + abs ( imag ( x ) )

  return
end
function c_norm2 ( x )
!
!*******************************************************************************
!
!! C_NORM2 evaluates the L2 norm of a complex number.
!
!
!  Discussion:
!
!    Numbers of equal norm lie on circles centered at (0,0).
!
!  Definition:
!
!    The L2 norm can be defined here as:
!
!      C_NORM2(X) = sqrt ( ( real (X) )**2 + ( imag ( X ) )**2 )
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
!    Input, complex X, the value whose norm is desired.
!
!    Output, real C_NORM2, the 2-norm of X.
!
  real c_norm2
  complex x
!
  c_norm2 = sqrt ( ( real ( x ) )**2 + ( imag ( x ) )**2 )

  return
end
function c_normi ( x )
!
!*******************************************************************************
!
!! C_NORMI evaluates the L-infinity norm of a complex number.
!
!
!  Discussion:
!
!    Numbers of equal norm lie along squares whose centers are at (0,0).
!
!  Definition:
!
!    The L-infinity norm can be defined here as:
!
!      C_NORMI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
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
!    Input, complex X, the value whose norm is desired.
!
!    Output, real C_NORMI, the infinity norm of X.
!
  real c_normi
  complex x
!
  c_normi = max ( abs ( real ( x ) ), abs ( imag ( x ) ) )

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
subroutine ch_cap ( c )
!
!*******************************************************************************
!
!! CH_CAP capitalizes a single character.
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
!    Input/output, character C, the character to capitalize.
!
  character c
  integer itemp
!
  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_is_digit ( c )
!
!*******************************************************************************
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  character c
  logical ch_is_digit
!
  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
subroutine ch_random ( clo, chi, c )
!
!*******************************************************************************
!
!! CH_RANDOM returns a random character in a given range.
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
!    Input, character CLO, CHI, the minimum and maximum acceptable characters.
!
!    Output, character C, the randomly chosen character.
!
  character c
  character chi
  character clo
  integer i
  integer ihi
  integer ilo
!
  ilo = ichar ( clo )
  ihi = ichar ( chi )

  call i_random ( ilo, ihi, i )

  c = char ( i )

  return
end
subroutine ch_swap ( c1, c2 )
!
!*******************************************************************************
!
!! CH_SWAP swaps two character values.
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
!    Input/output, character C1, C2.  On output, the values of C1 and
!    C2 have been interchanged.
!
  character c1
  character c2
  character c3
!
  c3 = c1
  c1 = c2
  c2 = c3

  return
end
subroutine chvec2_print ( m, a, n, b, title )
!
!*******************************************************************************
!
!! CHVEC2_PRINT prints two vectors of characters.
!
!
!  Modified:
!
!    15 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the length of the first sequence.
!
!    Input, character A(M), the first sequence.
!
!    Input, integer N, the length of the second sequence.
!
!    Input, character B(N), the second sequence.
!
!    Input, character ( len = * ) TITLE, a title.
!
  integer m
  integer n
!
  character a(m)
  character ai
  character b(n)
  character bi
  integer i
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '

  do i = 1, max ( m, n )

    if ( i <= m ) then
      ai = a(i)
    else
      ai = ' '
    end if

    if ( i <= n ) then
      bi = b(i)
    else
      bi = ' '
    end if

    write ( *, '(i3,2x,a1,2x,a1)' ) i, ai, bi

  end do

  return
end
subroutine chvec_permute ( n, a, p )
!
!*******************************************************************************
!
!! CHVEC_PERMUTE permutes a character vector in place.
!
!
!  Note:
!
!    This routine permutes an array of character "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (  'B', 'D', 'E', 'A', 'C' )
!
!    Output:
!
!      A    = ( 'A', 'B', 'C', 'D', 'E' ).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input/output, character A(N), the array to be permuted.
!
!    Input, integer P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  integer n
!
  character a(n)
  character a_temp
  integer i
  integer ierror
  integer iget
  integer iput
  integer istart
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHVEC_PERMUTE - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. iget > n ) then
          write ( *, * ) ' '
          write ( *, * ) 'CHVEC_PERMUTE - Fatal error!'
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine chvec_print ( n, a, title )
!
!*******************************************************************************
!
!! CHVEC_PRINT prints a character vector.
!
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, character A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  character a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,2x,a1)' ) i, a(i)
  end do

  return
end
subroutine chvec_reverse ( n, x )
!
!*******************************************************************************
!
!! CHVEC_REVERSE reverses the elements of a character vector.
!
!
!  Example:
!
!    Input:
!
!      N = 4, X = ( 'L', 'I', 'V', 'E' ).
!
!    Output:
!
!      X = ( 'E', 'V', 'I', 'L' ).
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, character X(N), the array to be reversed.
!
  integer n
!
  integer i
  character x(n)
!
  do i = 1, n/2
    call ch_swap ( x(i), x(n+1-i) )
  end do

  return
end
subroutine cmat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! CMAT_PRINT prints a complex matrix.
!
!
!  Modified:
!
!    28 March 2000
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
!    Input, complex A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer lda
  integer n
!
  complex a(lda,n)
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

  do jlo = 1, n, 3
    jhi = min ( jlo + 2, n )
    write ( *, * ) ' '
    write ( *, '(6x,3(10x,i8,10x))' ) ( j, j = jlo, jhi )
    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(i6,6g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
function cot ( angle )
!
!*******************************************************************************
!
!! COT returns the cotangent of an angle.
!
!
!  Definition:
!
!    COT ( THETA ) = COS ( THETA ) / SIN ( THETA )
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle, in radians.
!
!    Output, real COT, the cotangent of the angle.
!
  real angle
  real cot
!
  cot  = cos ( angle ) / sin ( angle )

  return
end
function cotd ( angle )
!
!*******************************************************************************
!
!! COTD returns the cotangent of an angle given in degrees.
!
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle, in degrees.
!
!    Output, real COTD, the cotangent of the angle.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
  real, parameter :: degrees_to_radians = PI / 180.0E+00
!
  real angle
  real cotd
!
  cotd  =     cos ( degrees_to_radians * angle ) &
            / sin ( degrees_to_radians * angle )

  return
end
function csc ( theta )
!
!*******************************************************************************
!
!! CSC returns the cosecant of X.
!
!
!  Definition:
!
!    CSC ( THETA ) = 1.0 / SIN ( THETA )
!
!  Discussion:
!
!    CSC is not a built-in function in FORTRAN, and occasionally it
!    is handier, or more concise, to be able to refer to it directly
!    rather than through its definition in terms of the sine function.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real THETA, the angle, in radians, whose cosecant is desired.
!    It must be the case that SIN ( THETA ) is not zero.
!
!    Output, real CSC, the cosecant of THETA.
!
  real csc
  real theta
!
  csc = sin ( theta )

  if ( csc == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CSC - Fatal error!'
    write ( *, * ) '  Cosecant undefined for THETA = ', theta
    stop
  end if

  csc = 1.0E+00 / csc

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
subroutine cvec_sort_a1 ( n, x )
!
!*******************************************************************************
!
!! CVEC_SORT_A1 ascending sorts a complex array by L1 norm.
!
!
!  Discussion:
!
!    The L1 norm of A+Bi is abs(A) + abs(B).
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
  logical c_le_l1
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

      if ( c_le_l1 ( x(i), x(j) ) ) then
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
subroutine cvec_sort_ainf ( n, x )
!
!*******************************************************************************
!
!! CVEC_SORT_AINF ascending sorts a complex array by L infinity norm.
!
!
!  Discussion:
!
!    The L infinity norm of A+Bi is max ( abs ( A ) , abs ( B ) ).
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
  logical c_le_linf
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

      if ( c_le_linf ( x(i), x(j) ) ) then
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
subroutine d_swap ( x, y )
!
!*******************************************************************************
!
!! D_SWAP swaps two double precision values.
!
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, double precision X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  double precision x
  double precision y
  double precision z
!
  z = x
  x = y
  y = z

  return
end
function dpi ( )
!
!*******************************************************************************
!
!! DPI returns the value of pi as a double precision quantity.
!
!
!  Modified:
!
!    10 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision DPI, the value of pi.
!
  double precision dpi
!
  dpi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
subroutine dvec_print ( n, a, title )
!
!*******************************************************************************
!
!! DVEC_PRINT prints a double precision vector.
!
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, double precision A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  integer n
!
  double precision a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g16.8)' ) i, a(i)
  end do

  return
end
function e_constant ( )
!
!*******************************************************************************
!
!! E_CONSTANT returns the value of E.
!
!
!  Discussion:
!
!    "E" was named in honor of Euler, but is known as Napier's constant.
!
!  Modified:
!
!    29 April 1999
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
  e_constant = 2.71828182845904523536028747135266249775724709369995E+00

  return
end
subroutine erf ( x, erfx )
!
!*******************************************************************************
!
!! ERF evaluates the error function ERF(X).
!
!
!  Reference:
!
!    W J Cody,
!    "Rational Chebyshev approximations for the error function",
!    Math. Comp., 1969, pages 631-638.
!
!  Author:
!
!    W J Cody,
!    Mathematics and Computer Science Division,
!    Argonne National Laboratory,
!    Argonne, Illinois, 60439.
!
!  Parameters:
!
!    Input, double precision X, the argument of ERF.
!
!    Output, double precision ERFX, the value of ERF(X).
!
  double precision, parameter :: SQRPI = 0.56418958354775628695D+00
  double precision, parameter :: THRESH = 0.46875D+00
  double precision, parameter :: XBIG = 26.543D+00
  double precision, parameter :: XSMALL = 1.11D-16
!
  double precision, save, dimension ( 5 ) :: a = (/ &
    3.16112374387056560D+00, &
    1.13864154151050156D+02, &
    3.77485237685302021D+02, &
    3.20937758913846947D+03, &
    1.85777706184603153D-01 /)
  double precision, save, dimension ( 4 ) :: b = (/ &
    2.36012909523441209D+01, &
    2.44024637934444173D+02, &
    1.28261652607737228D+03, &
    2.84423683343917062D+03 /)
  double precision, save, dimension ( 9 ) :: c = (/ &
    5.64188496988670089D-01, &
    8.88314979438837594D+00, &
    6.61191906371416295D+01, &
    2.98635138197400131D+02, &
    8.81952221241769090D+02, &
    1.71204761263407058D+03, &
    2.05107837782607147D+03, &
    1.23033935479799725D+03, &
    2.15311535474403846D-08 /)
  double precision, save, dimension ( 8 ) :: d = (/ &
    1.57449261107098347D+01, &
    1.17693950891312499D+02, &
    5.37181101862009858D+02, &
    1.62138957456669019D+03, &
    3.29079923573345963D+03, &
    4.36261909014324716D+03, &
    3.43936767414372164D+03, &
    1.23033935480374942D+03 /)
  double precision del
  double precision erfx
  integer i
  double precision, save, dimension ( 6 ) :: p = (/ &
    3.05326634961232344D-01, &
    3.60344899949804439D-01, &
    1.25781726111229246D-01, &
    1.60837851487422766D-02, &
    6.58749161529837803D-04, &
    1.63153871373020978D-02 /)
  double precision, save, dimension ( 5 ) :: q = (/ &
    2.56852019228982242D+00, &
    1.87295284992346047D+00, &
    5.27905102951428412D-01, &
    6.05183413124413191D-02, &
    2.33520497626869185D-03 /)
  double precision x
  double precision xabs
  double precision xden
  double precision xnum
  double precision xsq
!
  xabs = abs ( x )
!
!  Evaluate ERF(X) for |X| <= 0.46875.
!
  if ( xabs <= THRESH ) then

    if ( xabs > XSMALL ) then
      xsq = xabs * xabs
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do

    erfx = x * ( xnum + a(4) ) / ( xden + b(4) )
!
!  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
!
  else if ( xabs <= 4.0D+00 ) then

    xnum = c(9) * xabs
    xden = xabs
    do i = 1, 7
      xnum = ( xnum + c(i) ) * xabs
      xden = ( xden + d(i) ) * xabs
    end do

    erfx = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
    del = ( xabs - xsq ) * ( xabs + xsq )
    erfx = exp ( - xsq * xsq ) * exp ( - del ) * erfx

    erfx = ( 0.5D+00 - erfx ) + 0.5D+00

    if ( x < 0.0D+00 ) then
      erfx = - erfx
    end if
!
!  Evaluate ERFC(X) for |X| > 4.0.
!
  else

    if ( xabs >= XBIG ) then

      if ( x > 0.0D+00 ) then
        erfx = 1.0D+00
      else
        erfx = - 1.0D+00
      end if

    else

      xsq = 1.0D+00 / ( xabs * xabs )

      xnum = p(6) * xsq
      xden = xsq
      do i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
      end do

      erfx = xsq * ( xnum + p(5) ) / ( xden + q(5) )
      erfx = ( SQRPI -  erfx ) / xabs
      xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
      del = ( xabs - xsq ) * ( xabs + xsq )
      erfx = exp ( - xsq * xsq ) * exp ( - del ) * erfx

      erfx = ( 0.5D+00 - erfx ) + 0.5D+00
      if ( x < 0.0D+00 ) then
        erfx = - erfx
      end if

    end if

  end if

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
!      Gamma = limit ( M -> Infinity ) ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
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
subroutine fac_div ( nprime, npower1, npower2, npower3 )
!
!*******************************************************************************
!
!! FAC_DIV divides two quantities represented as prime factors.
!
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPRIME, the index of the highest prime number
!    used in the representations.
!
!    Input, integer NPOWER1(NPRIME), contains the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer NPOWER2(NPRIME), contains the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer NPOWER3(NPRIME), contains the powers of primes
!    in the representation of the quotient.
!
  integer nprime
!
  integer npower1(nprime)
  integer npower2(nprime)
  integer npower3(nprime)
!
  npower3(1:nprime) = npower1(1:nprime) - npower2(1:nprime)

  return
end
subroutine fac_gcd ( nprime, npower1, npower2, npower3 )
!
!*******************************************************************************
!
!! FAC_GCD finds the GCD of two products of prime factors.
!
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPRIME, the index of the highest prime number
!    used in the representations.
!
!    Input, integer NPOWER1(NPRIME), contains the powers of primes
!    in the representation of the first quantity.  All the powers
!    must be nonnegative.
!
!    Input, integer NPOWER2(NPRIME), contains the powers of primes
!    in the representation of the second quantity.  All the powers
!    must be nonnegative.
!
!    Output, integer NPOWER3(NPRIME), contains the powers of primes
!    in the representation of the GCD.
!
  integer nprime
!
  integer i
  integer npower1(nprime)
  integer npower2(nprime)
  integer npower3(nprime)
!
  do i = 1, nprime

    if ( npower1(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'FAC_GCD - Fatal error!'
      write ( *, * ) '  One of the powers is negative!'
      stop
    end if

    if ( npower2(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'FAC_GCD - Fatal error!'
      write ( *, * ) '  One of the powers is negative!'
      stop
    end if

    npower3(i) = min ( npower1(i), npower2(i) )

  end do

  return
end
subroutine fac_lcm ( nprime, npower1, npower2, npower3 )
!
!*******************************************************************************
!
!! FAC_LCM finds the LCM of two products of prime factors.
!
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPRIME, the index of the highest prime number
!    used in the representations.
!
!    Input, integer NPOWER1(NPRIME), contains the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer NPOWER2(NPRIME), contains the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer NPOWER3(NPRIME), contains the powers of primes
!    in the representation of the LCM.
!
  integer nprime
!
  integer i
  integer npower1(nprime)
  integer npower2(nprime)
  integer npower3(nprime)
!
  do i = 1, nprime

    if ( npower1(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'FAC_LCM - Fatal error!'
      write ( *, * ) '  One of the powers is negative!'
      stop
    end if

    if ( npower2(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'FAC_LCM - Fatal error!'
      write ( *, * ) '  One of the powers is negative!'
      stop
    end if

    npower3(i) = max ( npower1(i), npower2(i) )

  end do

  return
end
subroutine fac_mul ( nprime, npower1, npower2, npower3 )
!
!*******************************************************************************
!
!! FAC_MUL multiplies two quantities represented as prime factors.
!
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPRIME, the index of the highest prime number
!    used in the representations.
!
!    Input, integer NPOWER1(NPRIME), contains the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer NPOWER2(NPRIME), contains the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer NPOWER3(NPRIME), contains the powers of primes
!    in the representation of the product.
!
  integer nprime
!
  integer i
  integer npower1(nprime)
  integer npower2(nprime)
  integer npower3(nprime)
!
  npower3(1:nprime) = npower1(1:nprime) + npower2(1:nprime)

  return
end
subroutine fac_print ( nprime, npower )
!
!*******************************************************************************
!
!! FAC_PRINT prints a product of prime factors.
!
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPRIME, the index of the highest prime number
!    used in the representations.
!
!    Input, integer NPOWER(NPRIME), contains the powers of primes
!    in the representation of the quantity.
!
  integer nprime
!
  integer i
  integer npower(nprime)
  integer prime
!
  write ( *, * ) ' '
  write ( *, * ) '   Prime     Power'
  write ( *, * ) ' '
  do i = 1, nprime
    if ( npower(i) /= 0 ) then
      write ( *, '(i8,2x,i8)' ) prime(i), npower(i)
    end if
  end do

  return
end
subroutine fac_to_i ( nprime, npower, intval )
!
!*******************************************************************************
!
!! FAC_TO_I converts a product of prime factors into an integer.
!
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPRIME, the index of the highest prime number
!    used in the representations.
!
!    Input, integer NPOWER(NPRIME), contains the powers of primes
!    in the representation of the quantity.  If any of these powers
!    are negative, then INTVAL will be set to 0.
!
!    Output, integer INTVAL, the integer represented by the product of the
!    prime factors.
!
  integer nprime
!
  integer i
  integer intval
  integer npower(nprime)
  integer prime
!
  intval = 1
  do i = 1, nprime

    if ( npower(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'FAC_TO_I - Fatal error!'
      write ( *, * ) '  One of the powers is negative!'
      stop
    end if

    intval = intval * prime(i)**npower(i)

  end do

  return
end
subroutine fac_to_rat ( nprime, npower, itop, ibot )
!
!*******************************************************************************
!
!! FAC_TO_RAT converts a prime factorization into a rational value.
!
!
!  Example:
!
!    Start with the prime factorization representation:
!
!      40/9 = 2**3 * 3**(-2) * 5
!
!    Input:
!
!      NPOWER = ( 3, -2, 1 )
!
!    Output:
!
!      ITOP = 40 ( = 2**3 * 5**1 = PRIME(1)**3                 * PRIME(3)**1 )
!      IBOT = 9  ( = 3**2        =               PRIME(2)**2 )
!
!  Modified:
!
!    16 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPRIME, the index of the highest prime number
!    used in the representations.
!
!    Input, integer NPOWER(NPRIME).  NPOWER(I) is the power of
!    the I-th prime in the prime factorization.  NPOWER(I) may
!    be positive or negative.
!
!    Output, integer ITOP, IBOT, the top and bottom of a rational value.
!
  integer nprime
!
  integer i
  integer ibot
  integer itop
  integer npower(nprime)
  integer prime
!
  itop = 1
  ibot = 1
  do i = 1, nprime
    if ( npower(i) > 0 ) then
      itop = itop * prime(i)**npower(i)
    else if ( npower(i) < 0 ) then
      ibot = ibot * prime(i)**(-npower(i))
    end if
  end do

  return
end
function feet_to_meters ( ft )
!
!*******************************************************************************
!
!! FEET_TO_METERS converts a measurement in feet to meters.
!
!
!  Modified:
!
!    20 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real FT, the length in feet.
!
!    Output, real FEET_TO_METERS, the corresponding length in meters.
!
  real feet_to_meters
  real ft
!
  feet_to_meters = 0.0254E+00 * 12.0E+00 * ft

  return
end
function fgauss ( ngauss, amplitude, center, width, x )
!
!*******************************************************************************
!
!! FGAUSS evaluates a function that is the sum of Gaussians.
!
!
!  Definition:
!
!    FGAUSS(X) = Sum ( 1 <= I <= NGAUSS )
!    AMPLITUDE(I) * exp ( - ( ( X - CENTER(I) ) / WIDTH(I) )**2 )
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NGAUSS, the number of Gaussians that are summed.
!
!    Input, real AMPLITUDE(NGAUSS), CENTER(NGAUSS), WIDTH(NGAUSS),
!    the amplitude, center and width used for each Gaussian function.
!
!    Input, real X, the point at which Y is to be evaluated.
!
!    Output, real FGAUSS, the value of the function.
!
  integer ngauss
!
  real amplitude(ngauss)
  real arg
  real center(ngauss)
  real fgauss
  integer i
  real width(ngauss)
  real x
!
  fgauss = 0.0E+00

  do i = 1, ngauss

    arg = ( x - center(i) ) / width(i)

    fgauss = fgauss + amplitude(i) * exp ( - arg**2 )

  end do

  return
end
subroutine get_seed ( iseed )
!
!*******************************************************************************
!
!! GET_SEED returns a seed for the random number generator.
!
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ISEED, a pseudorandom seed value.
!
  integer, parameter :: I_MAX = 2147483647
!
  integer iseed
  double precision temp
  integer values(8)
!
  character ( len = 10 ) time
  character ( len = 8 ) today
  character ( len = 5 ) zone
!
  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + dble ( values(2) - 1 ) / 11.0D+00
  temp = temp + dble ( values(3) - 1 ) / 30.0D+00
  temp = temp + dble ( values(5) ) / 23.0D+00
  temp = temp + dble ( values(6) ) / 59.0D+00
  temp = temp + dble ( values(7) ) / 59.0D+00
  temp = temp + dble ( values(8) ) / 999.0D+00
  temp = temp / 6.0D+00

  if ( temp <= 0.0D+00 ) then
    temp = 1.0D+00 / 3.0D+00
  else if ( temp >= 1.0D+00 ) then
    temp = 2.0D+00 / 3.0D+00
  end if

  iseed = int ( dble ( I_MAX ) * temp )
!
!  Never use a seed of 0 or I_MAX.
!
  if ( iseed == 0 ) then
    iseed = 1
  end if

  if ( iseed == I_MAX ) then
    iseed = I_MAX - 1
  end if

  return
end
subroutine grid1 ( ndim, nstep, x, x1, x2 )
!
!*******************************************************************************
!
!! GRID1 finds grid points between X1 and X2 in N dimensions.
!
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of the points X1 and X2.
!
!    Input, integer NSTEP, the number of points to be generated.
!    NSTEP must be at least 2.
!
!    Output, real X(NDIM,NSTEP), the set of equally spaced points.  Each
!    column of X represents one point, with X(*,1) = X1 and X(*,NSTEP) = X2.
!
!    Input, real X1(NDIM), X2(NDIM), the first and last
!    points, between which the equally spaced points are
!    to be computed.
!
  integer ndim
  integer nstep
!
  integer i
  integer j
  real x(ndim,nstep)
  real x1(ndim)
  real x2(ndim)
!
  if ( nstep <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID1 - Fatal error!'
    write ( *, * ) '  NSTEP <= 1.'
    write ( *, * ) '  NSTEP = ', nstep
    stop
  end if

  do i = 1, nstep
    x(1:ndim,i) = ( real ( nstep - i ) * x1(1:ndim) &
        + real ( i - 1 ) * x2(1:ndim) ) / real ( nstep - 1 )
  end do

  return
end
subroutine grid1n ( i, ndim, nstep, x, x1, x2 )
!
!*******************************************************************************
!
!! GRID1N finds the I-th grid point between X1 and X2 in N dimensions.
!
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number of the desired point.
!    Normally I would be between 1 and NSTEP, but that is
!    not necessary.  Note that i = 1 returns X1 and i = NSTEP
!    returns X2.
!
!    Input, integer NDIM, the dimension of the points X, X1 and X2.
!
!    Input, integer NSTEP, this is the number of equally
!    spaced points that are between X1 and X2.  NSTEP must
!    be at least 2, because X1 and X2 are always included
!    in the set of points.
!
!    Output, real X(NDIM), the I-th grid point between X1 and X2.
!
!    Input, real X1(NDIM), X2(NDIM), the first and last
!    points, between which the equally spaced points lie.
!
  integer ndim
  integer nstep
!
  integer i
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
!
  if ( nstep <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID1N - Fatal error!'
    write ( *, * ) '  NSTEP <= 1.'
    write ( *, * ) '  NSTEP = ', nstep
    stop
  end if

  x(1:ndim) = ( real ( nstep - i ) * x1(1:ndim) &
              + real ( i - 1 )     * x2(1:ndim) ) / real ( nstep - 1 )

  return
end
subroutine grid2 ( i1, i2, ndim, nstep, x, x1, x2 )
!
!*******************************************************************************
!
!! GRID2 computes grid points between X1 and X2 in N dimensions.
!
!
!  Discussion:
!
!    However, X1 need not be the first point computed, nor X2 the last.
!    The user must specify the steps on which X1 and X2 are passed
!    through.  These steps may even be outside the range of 1 through NSTEP.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I1, I2.  I1 specifies the step on which
!    X1 would be computed, and similarly for I2.  That is,
!    we assume that a set of equally spaced points have
!    been drawn on the line through X1 and X2, and that
!    they have been numbered, with X1 labeled I1 and X2
!    labeled I2.  I1 or I2 may be between 1 and NSTEP,
!    in which case X1 or X2 will actually be returned in the
!    X array, but there is no requirement that I1 or I2
!    satisfy this condition.
!
!    I1 and I2 must be distinct.
!
!    Input, integer NDIM, the dimension of the points X1 and X2.
!
!    Input, integer NSTEP, this is the number of equally
!    spaced points that are to be generated.
!
!    NSTEP should be at least 1.
!
!    Output, real X(NDIM,NSTEP), the set of equally spaced points.
!    Each column of X represents one point.
!    If 1 <= I1 <= NSTEP, then X(*,I1) = X1, and similarly for I2.
!
!    Input, real X1(NDIM), X2(NDIM), the points that define the line along
!    which the equally spaced points are generated, and which may or may
!    not be included in the set of computed points.
!
  integer ndim
  integer nstep
!
  integer i
  integer i1
  integer i2
  integer j
  real x(ndim,nstep)
  real x1(ndim)
  real x2(ndim)
!
  if ( i1 == i2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID2 - Fatal error!'
    write ( *, * ) '  I1 = I2, leading to zero denominator.'
    write ( *, * ) '  I1 = ', i1, ' I2 = ', i2
    stop
  end if

  do i = 1, nstep
    do j = 1, ndim
      x(j,i) = ( real ( i2 - i ) * x1(j) + real ( i - i1 ) * x2(j) )  &
        / real ( i2 - i1 )
    end do
  end do

  return
end
subroutine grid2n ( i, i1, i2, ndim, x, x1, x2 )
!
!*******************************************************************************
!
!! GRID2N computes one grid point between X1 and X2 in N dimensions.
!
!
!  Discussion:
!
!    However, X1 need not be the first point computed, nor X2 the last.
!    The user must specify the steps on which X1 and X2 are passed through.
!
!  Modified:
!
!   02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the I coordinate of the desired point.
!    Note that if I = I1, X will be returned as X1, and if
!    I = I2, X will be returned as X2.
!
!    Input, integer I1, I2.  I1 specifies the step on which
!    X1 would be computed, and similarly for I2.  That is,
!    we assume that a set of equally spaced points have
!    been drawn on the line through X1 and X2, and that
!    they have been numbered, with X1 labeled I1 and X2
!    labeled I2.
!
!    I1 and I2 must be distinct.
!
!    Input, integer NDIM, the dimension of the points X1 and X2.
!
!    Output, real X(NDIM).  X(I) is the I-th point from the
!    set of equally spaced points.
!
!    Input, real X1(NDIM), X2(NDIM), the points that define
!    the line along which the equally spaced points are
!    generated, and which may or may not be included in the
!    set of computed points.
!
  integer ndim
!
  integer i
  integer i1
  integer i2
  integer j
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
!
  if ( i1 == i2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID2N - Fatal error!'
    write ( *, * ) '  I1 = I2, leading to zero denominator.'
    write ( *, * ) '  I1 = ', i1, ' I2 = ', i2
    stop
  end if

  do j = 1, ndim
    x(j) = ( real ( i2 - i ) * x1(j) + real ( i - i1 ) * x2(j) ) &
      / real ( i2 - i1 )
  end do

  return
end
subroutine grid3 ( ndim, nstep1, nstep2, x, x1, x2, x3 )
!
!*******************************************************************************
!
!! GRID3 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of the points X1, X2 and X3.
!
!    Input, integer NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.
!
!    That is, the line between X1 and X2 will have NSTEP1
!    points generated along it, and the line between X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
!    X3 are always included in the set of points.
!
!    Output, real X(NDIM,NSTEP1,NSTEP2), the set of equally
!    spaced points.  Fixing the second and third indices
!    of X represents one point, with the following special
!    values:
!
!      X(*,1,1)      = X1
!      X(*,NSTEP1,1) = X2
!      X(*,1,NSTEP2) = X3.
!
!    Input, real X1(NDIM), X2(NDIM), X3(NDIM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
  integer ndim
  integer nstep1
  integer nstep2
!
  integer i
  integer j
  integer k
  real psi1
  real psi2
  real psi3
  real x(ndim,nstep1,nstep2)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  if ( nstep1 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID3 - Fatal error!'
    write ( *, * ) '  NSTEP1 <= 1.'
    write ( *, * ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID3 - Fatal error!'
    write ( *, * ) '  NSTEP2 <= 1.'
    write ( *, * ) '  NSTEP2 = ', nstep2
    stop
  end if

  do i = 1, nstep1
    do j = 1, nstep2

      psi2 = real ( i - 1 ) / real ( nstep1 - 1 )
      psi3 = real ( j - 1 ) / real ( nstep2 - 1 )
      psi1 = 1.0E+00 - psi2 - psi3

      do k = 1, ndim
        x(k,i,j) = psi1 * x1(k) + psi2 * x2(k) + psi3 * x3(k)
      end do

    end do
  end do

  return
end
subroutine grid3n ( i, j, ndim, nstep1, nstep2, x, x1, x2, x3 )
!
!*******************************************************************************
!
!! GRID3N computes a parallelogram grid on 3 points in N dimensions.
!
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the parallelogram coordinates
!    of the point.  I measures steps from X1 to X2, and
!    J measures steps from X1 to X3.  Normally, I would
!    be between 1 and NSTEP1, J between 1 and NSTEP2,
!    but this is not necessary.
!
!    Input, integer NDIM, the dimension of the points X1, X2 and X3.
!
!    Input, integer NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.
!
!    That is, the line between X1 and X2 will have NSTEP1
!    points generated along it, and the line between X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
!    X3 are always included in the set of points.
!
!    Output, real X(NDIM), the point with coordinates (I,J)
!    from the the set of equally  spaced points.  The
!    following special values are:
!
!       I       J         X
!
!       1       1         X1
!       NSTEP1  1         X2
!       1       NSTEP2    X3
!
!    Input, real X1(NDIM), X2(NDIM), X3(NDIM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
  integer ndim
  integer nstep1
  integer nstep2
!
  integer i
  integer j
  real psi1
  real psi2
  real psi3
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  if ( nstep1 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID3N - Fatal error!'
    write ( *, * ) '  NSTEP1 <= 1.'
    write ( *, * ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID3N - Fatal error!'
    write ( *, * ) '  NSTEP2 <= 1.'
    write ( *, * ) '  NSTEP2 = ', nstep2
    stop
  end if

  psi2 = real ( i - 1 ) / real ( nstep1 - 1 )
  psi3 = real ( j - 1 ) / real ( nstep2 - 1 )
  psi1 = 1.0E+00 - psi2 - psi3

  x(1:ndim) = psi1 * x1(1:ndim) + psi2 * x2(1:ndim) + psi3 * x3(1:ndim)

  return
end
subroutine grid4 ( i1, i2, j1, j2, ndim, nstep1, nstep2, x, x1, x2, x3 )
!
!*******************************************************************************
!
!! GRID4 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!
!
!  Discussion:
!
!    Unlike GRID3, GRID4 does not necessarily place X1 at the
!    "origin" of the parallelogram, with X2 and X3 set at the
!    extreme I and J coordinates.  Instead, the user is free
!    to specify the I and J coordinates of the points, although
!    they are required to lie on a subparallelogram of the
!    larger one.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I1, I2, J1, J2.  If we imagine that the
!    main parallelogram is drawn first, with coordinate
!    ranges 1 <= I <= NSTEP1 and 1 <= J <= NSTEP2, then
!    these indices determine the (I,J) coordinates of the
!    three points, namely:
!
!      X1 : (I1,J1)
!      X2 : (I2,J1)
!      X3 : (I1,J2)
!
!    Of course, we actually start with the points X1, X2,
!    and X3, and they define a parallelogram and an (I,J)
!    coordinate system over the plane containing them.  We
!    then are free to consider the parallelogram defined
!    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
!    which may or may not contain any of the points X1, X2
!    and X3.
!
!    Input, integer NDIM, the dimension of the points X1, X2 and X3.
!
!    Input, integer NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.
!
!    That is, the line through X1 and X2 will have NSTEP1
!    points generated along it, and the line through X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    NSTEP1 and NSTEP2 should be at least 1.
!
!    Output, real X(NDIM,NSTEP1,NSTEP2), the set of equally
!    spaced points.  Fixing the second and third indices
!    of X represents one point.
!
!    Assuming that the indices I1, I2, J1 and J2 are "within
!    bounds", the following special values will be computed:
!
!      X(*,I1,J1) = X1
!      X(*,I2,J1) = X2
!      X(*,I1,J2) = X3.
!
!    Input, real X1(NDIM), X2(NDIM), X3(NDIM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
  integer ndim
  integer nstep1
  integer nstep2
!
  integer i
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer k
  real psi1
  real psi2
  real psi3
  real x(ndim,nstep1,nstep2)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  if ( nstep1 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4 - Fatal error!'
    write ( *, * ) '  NSTEP1 <= 1.'
    write ( *, * ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4 - Fatal error!'
    write ( *, * ) '  NSTEP2 <= 1.'
    write ( *, * ) '  NSTEP2 = ', nstep2
    stop
  end if

  if ( i1 == i2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4 - Fatal error!'
    write ( *, * ) '  I1 = I2'
    write ( *, * ) '  I1 = ', i1, ' I2 = ', i2
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4 - Fatal error!'
    write ( *, * ) '  J1 = J2'
    write ( *, * ) '  J1 = ', j1, ' J2 = ', j2
    stop
  end if

  do i = 1, nstep1
    do j = 1, nstep2

      psi2 = real ( i - i1 ) / real ( i2 - i1 )
      psi3 = real ( j - j1 ) / real ( j2 - j1 )
      psi1 = 1.0E+00 - psi2 - psi3

      do k = 1, ndim
        x(k,i,j) = psi1 * x1(k) + psi2 * x2(k) + psi3 * x3(k)
      end do

    end do
  end do

  return
end
subroutine grid4n ( i, i1, i2, j, j1, j2, ndim, nstep1, nstep2, x, x1, x2, x3)
!
!*******************************************************************************
!
!! GRID4N computes a single point on a parallelogram grid in N space.
!
!
!  Discussion:
!
!    The computation is identical to that of GRID4, except that
!    only one point at a time is computed.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the I coordinate of the point X.
!
!    Input, integer I1, I2.  See discussion of J1, J2.
!
!    Input, integer J, the J coordinate of the point X.
!
!    Input, integer J1, J2.  If we imagine that the
!    main parallelogram is drawn first, with coordinate
!    ranges 1 <= I <= NSTEP1 and 1 <= J <= NSTEP2, then
!    these indices determine the (I,J) coordinates of the
!    three points X1, X2, and X3, namely:
!
!      X1 : (I1,J1)
!      X2 : (I2,J1)
!      X3 : (I1,J2)
!
!    Of course, we actually start with the points X1, X2,
!    and X3, and they define a parallelogram and an (I,J)
!    coordinate system over the plane containing them.  We
!    then are free to consider the parallelogram defined
!    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
!    which may or may not contain any of the points X1, X2
!    and X3.
!
!    Input, integer NDIM, the dimension of the points X, X1, X2 and X3.
!
!    Input, integer NSTEP1, NSTEP2.  These are the number of
!    equally spaced points generated in the first and second
!    directions.
!
!    That is, the line through X1 and X2 will have NSTEP1
!    points generated along it, and the line through X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    NSTEP1 and NSTEP2 should be at least 1.
!
!    Output, real X(NDIM), the point whose parallelogram
!    coordinates are (I,J).
!
!    The following special values will be computed:
!
!      I  J  X
!
!      I1 J1 X1
!      I2 J2 X2
!      I1 J2 X3
!
!    Input, real X1(NDIM), X2(NDIM), X3(NDIM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
  integer ndim
  integer nstep1
  integer nstep2
!
  integer i
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer k
  real psi1
  real psi2
  real psi3
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  if ( nstep1 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4N - Fatal error!'
    write ( *, * ) '  NSTEP1 <= 1.'
    write ( *, * ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4N - Fatal error!'
    write ( *, * ) '  NSTEP2 <= 1.'
    write ( *, * ) '  NSTEP2 = ', nstep2
    stop
  end if

  if ( i1 == i2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4N - Fatal error!'
    write ( *, * ) '  I1 = I2'
    write ( *, * ) '  I1 = ', i1, ' I2 = ', i2
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GRID4N - Fatal error!'
    write ( *, * ) '  J1 = J2'
    write ( *, * ) '  J1 = ', j1,' J2 = ', j2
    stop
  end if

  psi2 = real ( i - i1 ) / real ( i2 - i1 )
  psi3 = real ( j - j1 ) / real ( j2 - j1 )
  psi1 = 1.0E+00 - psi2 - psi3

  do k = 1, ndim
    x(k) = psi1 * x1(k) + psi2 * x2(k) + psi3 * x3(k)
  end do

  return
end
function gud ( x )
!
!*******************************************************************************
!
!! GUD evaluates the Gudermannian function.
!
!
!  Definition:
!
!    The Gudermannian function relates the hyperbolic and trigonometric
!    functions.  For any argument X, there is a corresponding value
!    GAMMA so that
!
!      sinh(x) = tan(gamma).
!
!    The value GAMMA is called the Gudermannian of X.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the Gudermannian.
!
!    Output, real GUD, the value of the Gudermannian.
!
  real gud
  real x
!
  gud = 2.0E+00 * atan ( tanh ( 0.5E+00 * x ) )

  return
end
subroutine hexcol ( angle, r, g, b )
!
!*******************************************************************************
!
!! HEXCOL returns a color on the perimeter of the color hexagon.
!
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle in the color hexagon.  The sextants are
!    defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real R, G, B, RGB specifications for the color that lies
!    at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  real angle
  real b
  real g
  real r
!
  angle = mod ( angle, 360.0E+00 )

  if ( angle < 0.0E+00 ) then
    angle = angle + 360.0E+00
  end if

  if ( angle <= 60.0E+00 ) then
    r = 1.0E+00
    g = angle / 60.0E+00
    b = 0.0E+00
  else if ( angle <= 120.0E+00 ) then
    r = ( 120.0E+00 - angle ) / 60.0E+00
    g = 1.0E+00
    b = 0.0E+00
  else if ( angle <= 180.0E+00 ) then
    r = 0.0E+00
    g = 1.0E+00
    b = ( angle - 120.0E+00 ) / 60.0E+00
  else if ( angle <= 240.0E+00 ) then
    r = 0.0E+00
    g = ( 240.0E+00 - angle ) / 60.0E+00
    b = 1.0E+00
  else if ( angle <= 300.0E+00 ) then
    r = ( angle - 240.0E+00 ) / 60.0E+00
    g = 0.0E+00
    b = 1.0E+00
  else if ( angle <= 360.0E+00 ) then
    r = 1.0E+00
    g = 0.0E+00
    b = ( 360.0E+00 - angle ) / 60.0E+00
  end if

  return
end
subroutine i2vec_print ( n, a, title )
!
!*******************************************************************************
!
!! I2VEC_PRINT prints an I2 vector.
!
!
!  Modified:
!
!    03 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(2,N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  integer a(2,n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,2i10)' ) i, a(1,i), a(2,i)
  end do

  return
end
function i_divp ( i, j )
!
!*******************************************************************************
!
!! I_DIVP returns the smallest multiple of J greater than or equal to I.
!
!
!  Examples:
!
!    I  J  I_DIVP(I,J)
!
!    0  4    0
!    1  4    1
!    2  4    1
!    3  4    1
!    4  4    1
!    5  4    2
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
!    Input, integer I, the number to be analyzed.
!
!    Input, integer J, the number, multiples of which will
!    be compared against I.  J may not be zero.
!
!    Output, integer I_DIVP, the smallest multiple of J that
!    is greater than or equal to I.
!
  integer i
  integer i_divp
  integer j
!
  if ( j /= 0 ) then
    i_divp = 1 + ( i - 1 ) / j
  else
    i_divp = 0
    write ( *, * ) ' '
    write ( *, * ) 'I_DIVP - Fatal error!'
    write ( *, * ) '  The input value of J was zero!'
    stop
  end if

  return
end
subroutine i_factor ( n, maxfactor, nfactor, factor, power, nleft )
!
!*******************************************************************************
!
!! I_FACTOR factors an integer into prime factors.
!
!
!  Formula:
!
!    N = NLEFT * Product ( I = 1 to NFACTOR ) FACTOR(I)**POWER(I).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer to be factored.  N may be positive,
!    negative, or 0.
!
!    Input, integer MAXFACTOR, the maximum number of prime factors for
!    which storage has been allocated.
!
!    Output, integer NFACTOR, the number of prime factors of N discovered
!    by the routine.
!
!    Output, integer FACTOR(MAXFACTOR), the prime factors of N.
!
!    Output, integer POWER(MAXFACTOR).  POWER(I) is the power of
!    the FACTOR(I) in the representation of N.
!
!    Output, integer NLEFT, the factor of N that the routine could not
!    divide out.  If NLEFT is 1, then N has been completely factored.
!    Otherwise, NLEFT represents factors of N involving large primes.
!
  integer maxfactor
!
  integer factor(maxfactor)
  integer i
  integer maxprime
  integer n
  integer nleft
  integer nfactor
  integer p
  integer power(maxfactor)
  integer prime
!
  nfactor = 0

  factor(1:maxfactor) = 0
  power(1:maxfactor) = 0

  nleft = n

  if ( n == 0 ) then
    return
  end if

  if ( abs ( n ) == 1 ) then
    nfactor = 1
    factor(1) = 1
    power(1) = 1
    return
  end if
!
!  Find out how many primes we stored.
!
  maxprime = prime ( -1 )
!
!  Try dividing the remainder by each prime.
!
  do i = 1, maxprime

    p = prime ( i )

    if ( mod ( abs ( nleft ), p ) == 0 ) then

      if ( nfactor < maxfactor ) then

        nfactor = nfactor + 1
        factor(nfactor) = p

        do

          power(nfactor) = power(nfactor) + 1
          nleft = nleft / p

          if ( mod ( abs ( nleft ), p ) /= 0 ) then
            exit
          end if

        end do

        if ( abs ( nleft ) == 1 ) then
          exit
        end if

      end if

    end if

  end do

  return
end
function i_gcd ( i, j )
!
!*******************************************************************************
!
!! I_GCD finds the greatest common divisor of I and J.
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
!    Input, integer I, J, two numbers whose greatest common divisor
!    is desired.
!
!    Output, integer I_GCD, the greatest common divisor of I and J.
!
!    Note that only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I_GCD is the
!    largest common factor of I and J.
!
  integer i
  integer i_gcd
  integer ip
  integer iq
  integer ir
  integer j
!
  i_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set IP to the larger of I and J, IQ to the smaller.
!  This way, we can alter IP and IQ as we go.
!
  ip = max ( abs ( i ), abs ( j ) )
  iq = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    ir = mod ( ip, iq )

    if ( ir == 0 ) then
      exit
    end if

    ip = iq
    iq = ir

  end do

  i_gcd = iq

  return
end
function i_gcdb ( i, j, k )
!
!*******************************************************************************
!
!! I_GCDB finds the greatest common divisor of the form K**N of two numbers.
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
!    Input, integer I, J, two numbers whose greatest common divisor K**N
!    is desired.
!
!    Input, integer K, the possible divisor of I and J.
!
!    Output, integer I_GCDB, the greatest common divisor of
!    the form K**N shared by I and J.
!
!    Note that if J is negative, I_GCDB will also be negative.
!    This is because it is likely that the caller is forming
!    the fraction I/J, and so any minus sign should be
!    factored out of J.
!
!    If I and J are both zero, I_GCDB is returned as 1.
!
!    If I is zero and J is not, I_GCDB is returned as J,
!    and vice versa.
!
!    If I and J are nonzero, and have no common divisor of the
!    form K**N, I_GCDB is returned as 1.
!
!    Otherwise, I_GCDB is returned as the largest common divisor
!    of the form K**N shared by I and J.
!
  integer i
  integer icopy
  integer i_gcdb
  integer j
  integer jcopy
  integer k
!
  i_gcdb = 1
!
!  If both I and J are zero, I_GCDB is 1.
!
  if ( i == 0 .and. j == 0 ) then
    i_gcdb = 1
    return
  end if
!
!  If just one of I and J is zero, I_GCDB is the other one.
!
  if ( i == 0 ) then
    i_gcdb = j
    return
  else if ( j == 0 ) then
    i_gcdb = i
    return
  end if
!
!  Divide out K as long as you can.
!
  if ( j > 0 ) then
    i_gcdb = 1
  else
    i_gcdb = - 1
  end if

  icopy = i
  jcopy = j

  do

    if ( mod ( icopy, k ) /= 0 .or. mod ( jcopy, k ) /= 0 ) then
      exit
    end if

    i_gcdb = i_gcdb * k
    icopy = icopy / k
    jcopy = jcopy / k

  end do

  return
end
function i_is_prime ( n )
!
!*******************************************************************************
!
!! I_IS_PRIME reports whether an integer is prime.
!
!
!  Method:
!
!    A simple, unoptimized sieve of Erasthosthenes is used to
!    check whether N can be divided by any integer between 2
!    and SQRT(N).
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer to be tested.
!
!    Output, logical I_IS_PRIME, is TRUE if N is prime, and FALSE
!    otherwise.  Note that negative numbers and 0 are not
!    considered prime.
!
  integer i
  logical i_is_prime
  integer n
  integer nhi
!
  if ( n <= 0 ) then
    i_is_prime = .false.
    return
  end if

  if ( n <= 3 ) then
    i_is_prime = .true.
    return
  end if

  nhi = int ( sqrt ( real ( n ) ) )

  do i = 2, nhi
    if ( mod ( n, i ) == 0 ) then
      i_is_prime = .false.
      return
    end if
  end do

  i_is_prime = .true.

  return
end
subroutine i_jacobi_symbol ( q, p, j )
!
!*******************************************************************************
!
!! I_JACOBI_SYMBOL evaluates the Jacobi symbol (Q/P).
!
!
!  Definition:
!
!    If P is prime, then
!      Jacobi Symbol (Q/P) = Legendre Symbol (Q/P)
!    Else let P have the prime factorization
!      P = Product ( I = 1 to N ) P(I)**E(I)
!    Then
!      Jacobi Symbol (Q/P) =
!    Product ( I = 1 to N ) Legendre Symbol (Q/P(I))**E(I)
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 86-87.
!
!  Parameters:
!
!    Input, integer Q, an integer whose Jacobi symbol with
!    respect to P is desired.
!
!    Input, integer P, the number with respect to which the Jacobi
!    symbol of Q is desired.  P should be 2 or greater.
!
!    Output, integer L, the Jacobi symbol (Q/P).
!    Ordinarily, L will be -1, 0 or 1.
!    -2, not enough factorization space.
!    -3, an error during Legendre symbol calculation.
!
  integer, parameter :: maxfactor = 20
!
  integer factor(maxfactor)
  integer i
  integer j
  integer l
  integer nfactor
  integer nleft
  integer p
  integer power(maxfactor)
  integer pp
  integer q
  integer qq
!
!  P must be greater than 1.
!
  if ( p <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_JACOBI_SYMBOL - Fatal error!'
    write ( *, * ) '  P must be greater than 1.'
    l = -2
    return
  end if
!
!  Decompose P into factors of prime powers.
!
  call i_factor ( p, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_JACOBI_SYMBOL - Fatal error!'
    write ( *, * ) '  Not enough factorization space.'
    j = -2
    return
  end if
!
!  Force Q to be nonnegative.
!
  qq = q

  do while ( qq < 0 )
    qq = qq + p
  end do
!
!  For each prime factor, compute the Legendre symbol, and
!  multiply the Jacobi symbol by the appropriate factor.
!
  j = 1
  do i = 1, nfactor
    pp = factor(i)
    call i_legendre_symbol ( qq, pp, l )
    if ( l < -1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'I_JACOBI_SYMBOL - Fatal error!'
      write ( *, * ) '  Error during Legendre symbol calculation.'
      j = -3
      return
    end if
    j = j * l**power(i)
  end do

  return
end
function i_lcm ( i, j )
!
!*******************************************************************************
!
!! I_LCM computes the least common multiple of two integers.
!
!
!  Definition:
!
!    The least common multiple may be defined as
!
!      LCM(I,J) = ABS( I * J ) / GCD(I,J)
!
!    where GCD(I,J) is the greatest common divisor of I and J.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the integers whose I_LCM is desired.
!
!    Output, integer I_LCM, the least common multiple of I and J.
!    I_LCM is never negative.  I_LCM is 0 if either I or J is zero.
!
  integer i
  integer i_gcd
  integer j
  integer i_lcm
!
  i_lcm = abs ( i * ( j / i_gcd ( i, j ) ) )

  return
end
subroutine i_legendre_symbol ( q, p, l )
!
!*******************************************************************************
!
!! I_LEGENDRE_SYMBOL evaluates the Legendre symbol (Q/P).
!
!
!  Definition:
!
!    Let P be an odd prime.  Q is a QUADRATIC RESIDUE modulo P
!    if there is an integer R such that R**2 = Q ( mod P ).
!    The Legendre symbol ( Q / P ) is defined to be:
!
!      + 1 if Q ( mod P ) /= 0 and Q is a quadratic residue modulo P,
!      - 1 if Q ( mod P ) /= 0 and Q is not a quadratic residue modulo P,
!        0 if Q ( mod P ) == 0.
!
!    We can also define ( Q / P ) for P = 2 by:
!
!      + 1 if Q ( mod P ) /= 0
!        0 if Q ( mod P ) == 0
!
!  Example:
!
!    (0/7) =   0
!    (1/7) = + 1  ( 1**2 = 1 mod 7 )
!    (2/7) = + 1  ( 3**2 = 2 mod 7 )
!    (3/7) = - 1
!    (4/7) = + 1  ( 2**2 = 4 mod 7 )
!    (5/7) = - 1
!    (6/7) = - 1
!
!  Note:
!
!    For any prime P, exactly half of the integers from 1 to P-1
!    are quadratic residues.
!
!    ( 0 / P ) = 0.
!
!    ( Q / P ) = ( mod ( Q, P ) / P ).
!
!    ( Q / P ) = ( Q1 / P ) * ( Q2 / P ) if Q = Q1 * Q2.
!
!    If Q is prime, and P is prime and greater than 2, then:
!
!      if ( Q == 1 ) then
!
!        ( Q / P ) = 1
!
!      else if ( Q == 2 ) then
!
!        ( Q / P ) = + 1 if mod ( P, 8 ) = 1 or mod ( P, 8 ) = 7,
!        ( Q / P ) = - 1 if mod ( P, 8 ) = 3 or mod ( P, 8 ) = 5.
!
!      else
!
!        ( Q / P ) = - ( P / Q ) if Q = 3 ( mod 4 ) and P = 3 ( mod 4 ),
!                  =   ( P / Q ) otherwise.
!
!  Modified:
!
!    28 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Pinter,
!    A Book of Abstract Algebra,
!    McGraw Hill, 1982, pages 236-237.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 86-87.
!
!  Parameters:
!
!    Input, integer Q, an integer whose Legendre symbol with
!    respect to P is desired.
!
!    Input, integer P, a prime number, greater than 1, with respect
!    to which the Legendre symbol of Q is desired.
!
!    Output, integer L, the Legendre symbol (Q/P).
!    Ordinarily, L will be -1, 0 or 1.
!    L = -2, P is less than or equal to 1.
!    L = -3, P is not prime.
!    L = -4, the internal stack of factors overflowed.
!    L = -5, not enough factorization space.
!
  integer, parameter :: maxfactor = 20
  integer, parameter :: maxstack = 50
!
  integer factor(maxfactor)
  integer i
  logical i_is_prime
  integer l
  integer nfactor
  integer nleft
  integer nmore
  integer nstack
  integer p
  integer power(maxfactor)
  integer pp
  integer pstack(maxstack)
  integer q
  integer qq
  integer qstack(maxstack)
!
!  P must be greater than 1.
!
  if ( p <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_LEGENDRE_SYMBOL - Fatal error!'
    write ( *, * ) '  P must be greater than 1.'
    l = -2
    return
  end if
!
!  P must be prime.
!
  if ( .not. i_is_prime ( p ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_LEGENDRE_SYMBOL - Fatal error!'
    write ( *, * ) '  P is not prime.'
    l = -3
    return
  end if
!
!  ( k*P / P ) = 0.
!
  if ( mod ( q, p ) == 0 ) then
    l = 0
    return
  end if
!
!  For the special case P = 2, (Q/P) = 1 for all odd numbers.
!
  if ( p == 2 ) then
    l = 1
    return
  end if
!
!  Make a copy of Q, and force it to be nonnegative.
!
  qq = q

  do while ( qq < 0 )
    qq = qq + p
  end do

  nstack = 0
  pp = p
  l = 1

  do

    qq = mod ( qq, pp )
!
!  Decompose QQ into factors of prime powers.
!
    call i_factor ( qq, maxfactor, nfactor, factor, power, nleft )

    if ( nleft /= 1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'I_LEGENDRE_SYMBOL - Fatal error!'
      write ( *, * ) '  Not enough factorization space.'
      l = - 5
      return
    end if
!
!  Each factor which is an odd power is added to the stack.
!
    nmore = 0

    do i = 1, nfactor

      if ( mod ( power(i), 2 ) == 1 ) then

        nmore = nmore + 1
        nstack = nstack + 1

        if ( nstack > maxstack ) then
          write ( *, * ) ' '
          write ( *, * ) 'I_LEGENDRE_SYMBOL - Fatal error!'
          write ( *, * ) '  Stack overflow!'
          l = - 4
          return
        end if

        pstack(nstack) = pp
        qstack(nstack) = factor(i)

      end if

    end do

    if ( nmore /= 0 ) then

      qq = qstack(nstack)
      nstack = nstack - 1
!
!  Check for a QQ of 1 or 2.
!
      if ( qq == 1 ) then

        l = + 1 * l

      else if ( qq == 2 .and. &
              ( mod ( pp, 8 ) == 1 .or. mod ( pp, 8 ) == 7 ) ) then

        l = + 1 * l

      else if ( qq == 2 .and. &
              ( mod ( pp, 8 ) == 3 .or. mod ( pp, 8 ) == 5 ) ) then

        l = - 1 * l

      else

        if ( mod ( pp, 4 ) == 3 .and. mod ( qq, 4 ) == 3 ) then
          l = - 1 * l
        end if

        call i_swap ( pp, qq )

        cycle

      end if

    end if
!
!  If the stack is empty, we're done.
!
    if ( nstack == 0 ) then
      exit
    end if
!
!  Otherwise, get the last P and Q from the stack, and process them.
!
    pp = pstack(nstack)
    qq = qstack(nstack)
    nstack = nstack - 1

  end do

  return
end
function i_log_10 ( x )
!
!*******************************************************************************
!
!! I_LOG_10 returns the integer part of the logarithm base 10 of ABS(X).
!
!
!  Example:
!
!    if 0.01 <  X <= 0.1   I_LOG_10(X) = -1,
!    if 0.1  <  X <= 1     I_LOG_10(X) = 0,
!    if 1    <= X <  10    I_LOG_10(X) = 0,
!    if 10   <= X <  100   I_LOG_10(X) = 1,
!    if 100  <= X <  1000  I_LOG_10(X) = 2.
!
!  Remark:
!
!    To get the exact logarithm, base 10, of a positive number, use
!    the intrinsic LOG10 function, as in Y = LOG10 ( X )
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number for which the whole number part of the
!    logarithm base 10 of the absolute value is desired.
!
!    Output, integer I_LOG_10, the whole number part of the
!    logarithm base 10 of the absolute value of X.
!
!    For positive I_LOG_10(X), it should be true that
!
!      10**I_LOG_10(X) <= ABS(X) < 10**(I_LOG_10(X)+1).
!
!    However, if X = 0, then I_LOG_10 = 0.
!
  integer i_log_10
  real temp
  real x
!
  i_log_10 = 0

  if ( x == 0.0E+00 ) then
    return
  end if

  if ( x == 1.0E+00 .or. x == -1.0E+00 ) then
    return
  end if

  temp = abs ( x )

  do while ( temp >= 10.0E+00 )
    temp = temp / 10.0E+00
    i_log_10 = i_log_10 + 1
  end do

  do while ( 10.0E+00 * temp <= 1.0E+00 )
    temp = temp * 10.0E+00
    i_log_10 = i_log_10 - 1
  end do

  return
end
function i_log_2 ( x )
!
!*******************************************************************************
!
!! I_LOG_2 returns the integer part of the logarithm base 2 of ABS(X).
!
!
!  Example:
!
!    if 1/4  <  X <= 1/2  I_LOG_2(X) = -1,
!    if 1/2  <  X <= 1    I_LOG_2(X) = 0,
!    if 1    <= X <  2    I_LOG_2(X) = 0,
!    if 2    <= X <  4    I_LOG_2(X) = 1,
!    if 4    <= X <  8    I_LOG_2(X) = 2.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose logarithm base 2 is desired.
!
!    Output, integer I_LOG_2, the integer part of the logarithm base 2 of
!    the absolute value of X.  For positive I_LOG_2(X), it should be true that
!      2**I_LOG_2(X) <= ABS(X) < 2**(I_LOG_2(X)+1).
!    The special case of I_LOG_2(0.0) returns 0.
!
  integer i_log_2
  real temp
  real x
!
  i_log_2 = 0

  if ( x == 0.0E+00 ) then
    return
  end if

  temp = abs ( x )

  if ( 1.0E+00 <= temp .and. temp < 2.0E+00 ) then
    return
  end if

  do while ( temp >= 2.0E+00 )
    temp = temp / 2.0E+00
    i_log_2 = i_log_2 + 1
  end do

  do while ( temp <= 0.5E+00 )
    temp = 2.0E+00 * temp
    i_log_2 = i_log_2 - 1
  end do

  return
end
function i_log_b ( x, b )
!
!*******************************************************************************
!
!! I_LOG_B returns the integer part of the logarithm base ABS(B) of ABS(X).
!
!
!  Example:
!
!    If B is greater than 1, and X is positive:
!
!    if 1/B**2 <  X <= 1/B   I_LOG_B(X) = -1,
!    if 1/B    <  X <= 1     I_LOG_B(X) = 0,
!    if 1  <= X <  B,    I_LOG_B(X) = 0,
!    if B  <= X <  B**2  I_LOG_B(X) = 1,
!    if B**2   <= X <  B**3  I_LOG_B(X) = 2.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose logarithm base B is desired.
!    If X is 0.0, then I_LOG_B is returned as 0, no matter what B is.
!
!    Input, real B, the absolute value of the base of the
!    logarithms.  B must not be -1, 0, or 1.
!
!    Output, integer I_LOG_B, the integer part of the logarithm
!    base abs(B) of X.
!
!    For positive I_LOG_B(X), it should be true that
!
!      ABS(B)**I_LOG_B(X) <= ABS(X) < ABS(B)**(I_LOG_B(X)+1).
!
  real b
  real base
  integer i_log_b
  real temp
  real x
!
  i_log_b = 0

  if ( x == 0.0E+00 ) then
    return
  end if

  if ( abs ( b ) == 1.0E+00 ) then
    return
  end if

  if ( b == 0.0E+00 ) then
    return
  end if

  temp = abs ( x )
  base = abs ( b )

  if ( base < 1.0E+00 ) then
    base = 1.0E+00 / base
  end if

  if ( 1.0E+00 <= temp .and. temp < base ) then
    return
  end if

  do while ( temp >= base )
    temp = temp / base
    i_log_b = i_log_b + 1
  end do

  do while ( base * temp <= 1.0E+00 )
    temp = base * temp
    i_log_b = i_log_b - 1
  end do

  if ( abs ( b ) < 1.0E+00 ) then
    i_log_b = - i_log_b
  end if

  return
end
subroutine i_mant ( x, is, j, k, l )
!
!*******************************************************************************
!
!! I_MANT computes the "mantissa" of a real number.
!
!
!  Discussion:
!
!    I_MANT computes the "mantissa" or "fraction part" of a real
!    number X, which it stores as a pair of integers, (J/K).
!
!    It also computes the sign, and the integer part of the logarithm
!    (base 2) of X.
!
!    On return from I_MANT, IS, J, K and L satisfy:
!
!      X = IS * (J/K) * 2**L
!
!    where IS is +1 or -1, K is a power of 2, (J/K) is greater than
!    or equal to 1 and strictly less than 2, and L is an integer.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the real number to be decomposed.
!
!    Output, integer IS, the "sign" of the number.
!    IS will be -1 if X is less than 0, and +1 if X is greater
!    than or equal to zero.
!
!    Output, integer J, the top part of the mantissa fraction.
!
!    Output, integer K, the bottom part of the mantissa
!    fraction.  K is a power of 2.
!
!    Output, integer L, the integer part of the logarithm (base 2) of X.
!
  integer is
  integer j
  integer k
  integer l
  real x
  real xtemp
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0E+00 ) then
    is = 1
    j = 0
    k = 1
    l = 0
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( x > 0.0E+00 ) then
    is = 1
    xtemp = x
  else
    is = - 1
    xtemp = - x
  end if
!
!  3: Force XTEMP to lie between 1 and 2, and compute the logarithm L.
!
  l = 0

  do while ( xtemp >= 2.0E+00 )
    xtemp = xtemp / 2.0E+00
    l = l + 1
  end do

  do while ( xtemp < 1.0E+00 )
    xtemp = xtemp * 2.0E+00
    l = l - 1
  end do
!
!  4: Now strip out the mantissa as J/K.
!
  j = 0
  k = 1

  do

    j = 2 * j

    if ( xtemp >= 1.0E+00 ) then
      j = j + 1
      xtemp = xtemp - 1.0E+00
    end if

    if ( xtemp == 0.0E+00 ) then
      exit
    end if

    k = 2 * k
    xtemp = xtemp * 2.0E+00

  end do

  return
end
subroutine i_memory ( action, name, ival )
!
!*******************************************************************************
!
!! I_MEMORY manages a set of runtime integer variables.
!
!
!  Discussion:
!
!    I_MEMORY allows the user to define the name of an integer variable,
!    set it, increment it, push a new value onto a stack of up to
!    five values, pop a value off the stack, or get the current value.
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, desired action.
!
!    'GET',  return value of NAME in IVAL.
!    'INC',  increment variable NAME by IVAL.
!    'INIT', reset all values to zero, wipe out all names.
!    'NAME', add a variable of the given name.
!    'POP',  pop the stack, retrieving the last pushed value.
!    'PRINT', print the value of variable NAME, and return it in IVAL.
!            NAME = '*' prints all variables, and returns IVAL = 0.
!    'PUSH', push the previous value down the stack, insert a new one.
!    'SET',  set variable NAME to IVAL.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    NAME may be blank for the 'INIT' command, but never any other time.
!
!    Input/output, integer IVAL.
!
!    For the 'INC', 'NAME', 'PUSH', and 'SET' commands, IVAL must
!    contain the increment, initial value, push value, or set value
!    of the variable on input.
!
!    For the 'GET' and 'POP' commands, IVAL will contain the current
!    value or popped value of the named variable on output.
!
  integer, parameter :: maxnam = 100
  integer, parameter :: maxcol = 5
!
  character ( len = * ) action
  integer i
  integer, save, dimension ( maxnam ) :: ipoint = (/ ( 0, i = 1, maxnam ) /)
  integer ival
  integer, save, dimension ( maxnam, maxcol ) :: ivals
  integer j
  character ( len = * ) name
  character ( len = 20 ), save, dimension ( maxnam ) :: names = &
    (/ ( ' ', i = 1, maxnam ) /)
  integer, save :: numnam = 0
  logical s_eqi
!
  data ( ( ivals(i,j), i = 1, maxnam ), j = 1, maxcol ) / 500 * 0 /
!
  if ( name == ' ' .and. .not. s_eqi ( action, 'INIT' ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  The value of NAME cannot be blank!'
    stop
  end if
!
!  GET: Get the current value of a variable.
!
  if ( s_eqi ( action, 'GET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        ival = ivals(i,ipoint(i))
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to get value of unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  INC: Increment something.
!
  else if ( s_eqi ( action, 'INC' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        ivals(i,ipoint(i)) = ivals(i,ipoint(i)) + ival
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to increment unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  INIT: Initialize everything.
!
  else if ( s_eqi ( action, 'INIT' ) ) then

    numnam = 0
    ivals(1:maxnam,1:maxcol) = 0
    names(1:maxnam) = ' '
    ipoint(1:maxnam) = 0
!
!  NAME: Declare the name of something and set it.
!
  else if ( s_eqi ( action, 'NAME' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        write ( *, * ) ' '
        write ( *, * ) 'I_MEMORY - Warning!'
        write ( *, '(''  There is ALREADY a variable '', a )' ) trim ( name )
        write ( *, * ) '  The new value has been stored.'
        ipoint(i) = 1
        ivals(i,ipoint(i)) = ival
        return
      end if

    end do

    if ( numnam < maxnam ) then
      numnam = numnam + 1
      i = numnam
      names(i) = name
      ipoint(i) = 1
      ivals(i,ipoint(i)) = ival
    else
      write ( *, * ) ' '
      write ( *, * ) 'I_MEMORY - Fatal error!'
      write ( *, * ) '  We have reached the name limit of ',maxnam
      stop
    end if
!
!  POP: "Pop" a value, decrement pointer.
!
  else if ( s_eqi ( action, 'POP' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) > 1 ) then
          ipoint(i) = ipoint(i) - 1
          ival = ivals(i,ipoint(i))
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'I_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to pop the stack to 0.'
          write ( *, '(''  Variable name is '',a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to pop an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  PRINT: "Print" the value, and return in IVAL.
!
  else if ( s_eqi ( action, 'PRINT' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) .or. name == '*' ) then
        ival = ivals(i,ipoint(i))
        write ( *, '(a,a,a,i8)' ) &
          'I_MEMORY - Value of ', trim ( names(i) ), ' is ', ival
      end if

      if ( s_eqi ( name, names(i) ) ) then
        return
      end if

    end do

    if ( name == '*' ) then
      ival = 0
      return
    end if

    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to print an unknown variable.'
    write ( *, '(''  Variable name is '', a )' ) trim ( name )
    stop
!
!  PUSH: "Push" a value, increment the pointer.
!
  else if ( s_eqi ( action, 'PUSH' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) < maxcol ) then
          ipoint(i) = ipoint(i) + 1
          ivals(i,ipoint(i)) = ival
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'I_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to push the stack past ', maxcol
          write ( *, '(''  Variable name is '',a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to push an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  SET: Set something.
!
  else if ( s_eqi ( action, 'SET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        ivals(i,ipoint(i)) = ival
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  Unrecognized action.
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'I_MEMORY - Fatal error!'
    write ( *, * ) '  Unrecognized action:'
    write ( *, '(a)' ) trim ( action )
    stop

  end if

  return
end
subroutine i_moddiv ( number, ndivid, nmult, nrem )
!
!*******************************************************************************
!
!! I_MODDIV breaks a number into a multiple of a divisor and remainder.
!
!
!  Formula:
!
!    NUMBER = NMULT * NDIVID + NREM
!
!    0 <= || NREM || < || NDIVID ||
!
!    NREM has the sign of NUMBER.
!
!  Example:
!
!    NUMBER   NDIVID  NMULT   NREM
!
!   107       50      2      7
!   107      -50     -2      7
!  -107       50     -2     -7
!  -107      -50      2     -7
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NUMBER, the number to be decomposed.
!
!    Input, integer NDIVID, the divisor.  NDIVID may not be zero.
!
!    Output, integer NMULT, the number of times NUMBER
!    is evenly divided by NDIVID.
!
!    Output, integer NREM, a remainder.
!
  integer ndivid
  integer nmult
  integer nrem
  integer number
!
  if ( ndivid == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_MODDIV - Fatal error!'
    write ( *, * ) '  Input divisor NDIVID = 0'
    stop
  end if

  nmult = number / ndivid
  nrem = number - ndivid * nmult

  return
end
function i_modp ( i, j )
!
!*******************************************************************************
!
!! I_MODP returns the nonnegative remainder of integer division.
!
!
!  Formula:
!
!    If
!      NREM = I_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I_MODP(A,360) is between 0 and 360, always.
!
!  Examples:
!
!        I     J     MOD  I_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I_MODP, the nonnegative remainder when I is
!    divided by J.
!
  integer i
  integer i_modp
  integer j
!
  if ( j == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_MODP - Fatal error!'
    write ( *, * ) '  I_MODP ( I, J ) called with J = ', j
    stop
  end if

  i_modp = mod ( i, j )

  if ( i_modp < 0 ) then
    i_modp = i_modp + abs ( j )
  end if

  return
end
subroutine i_moebius ( n, mu )
!
!*******************************************************************************
!
!! I_MOEBIUS returns the value of MU(N), the Moebius function of N.
!
!
!  Definition:
!
!    MU(N) is defined as follows:
!
!      MU(N) = 1 if N = 1;
!              0 if N is divisible by the square of a prime;
!              (-1)**K, if N is the product of K distinct primes.
!
!  First values:
!
!     N  MU(N)
!
!     1    1
!     2   -1
!     3   -1
!     4    0
!     5   -1
!     6    1
!     7   -1
!     8    0
!     9    0
!    10    1
!    11   -1
!    12    0
!    13   -1
!    14    1
!    15    1
!    16    0
!    17   -1
!    18    0
!    19   -1
!    20    0
!
!  Note:
!
!    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
!    if N is a square, cube, etc.
!
!  Formula:
!
!    The Moebius function is related to Euler's totient function:
!
!  PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the value to be analyzed.
!
!    Output, integer MU, the value of MU(N).
!    If N is less than or equal to 0, MU will be returned as -2.
!    If there was not enough internal space for factoring, MU
!    is returned as -3.
!
  integer, parameter :: maxfactor = 20
!
  integer factor(maxfactor)
  integer i
  integer mu
  integer n
  integer nfactor
  integer nleft
  integer power(maxfactor)
!
  if ( n <= 0 ) then
    mu = - 2
    return
  end if

  if ( n == 1 ) then
    mu = 1
    return
  end if
!
!  Factor N.
!
  call i_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_MOEBIUS - Fatal error!'
    write ( *, * ) '  Not enough factorization space.'
    mu = - 3
    return
  end if

  mu = 1

  do i = 1, nfactor

    mu = - mu

    if ( power(i) > 1 ) then
      mu = 0
      return
    end if

  end do

  return
end
subroutine i_omega ( n, omega )
!
!*******************************************************************************
!
!! I_OMEGA returns OMEGA(N), the number of distinct prime divisors of N.
!
!
!  First values:
!
!     N   OMEGA(N)
!
!     1    1
!     2    1
!     3    1
!     4    1
!     5    1
!     6    2
!     7    1
!     8    1
!     9    1
!    10    2
!    11    1
!    12    2
!    13    1
!    14    2
!    15    2
!    16    1
!    17    1
!    18    2
!    19    1
!    20    2
!
!  Formula:
!
!    If N = 1, then
!
!      OMEGA(N) = 1
!
!    else if the prime factorization of N is
!
!      N = P1**E1 * P2**E2 * ... * PM**EM,
!
!    then
!
!      OMEGA(N) = M
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the value to be analyzed.  N must be 1 or
!    greater.
!
!    Output, integer OMEGA, the value of OMEGA(N).  But if N is 0 or
!    less, OMEGA is returned as 0, a nonsense value.  If there is
!    not enough room for factoring, OMEGA is returned as -1.
!
  integer, parameter :: maxfactor = 20
!
  integer factor(maxfactor)
  integer n
  integer nfactor
  integer nleft
  integer omega
  integer power(maxfactor)
!
  if ( n <= 0 ) then
    omega = 0
    return
  end if

  if ( n == 1 ) then
    omega = 1
    return
  end if
!
!  Factor N.
!
  call i_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_OMEGA - Fatal error!'
    write ( *, * ) '  Not enough factorization space.'
    omega = - 1
    return
  end if

  omega = nfactor

  return
end
subroutine i_phi ( n, phin )
!
!*******************************************************************************
!
!! I_PHI returns the value of PHI(N), the number of relatively prime predecessors.
!
!
!  Definition:
!
!    PHI(N) is the number of integers between 1 and N which are
!    relatively prime to N.  I and J are relatively prime if they
!    have no common factors.  The function PHI(N) is known as
!    "Euler's totient function".
!
!    By convention, 1 and N are relatively prime.
!
!  First values:
!
!     N  PHI(N)
!
!     1    1
!     2    1
!     3    2
!     4    2
!     5    4
!     6    2
!     7    6
!     8    4
!     9    6
!    10    4
!    11   10
!    12    4
!    13   12
!    14    6
!    15    8
!    16    8
!    17   16
!    18    6
!    19   18
!    20    8
!
!  Formula:
!
!    PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
!
!    PHI(P**K) = P**(K-1) * ( P - 1 ) if P is prime.
!
!    PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
!
!    N = Sum ( D divides N ) PHI(D).
!
!  Modified:
!
!    01 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the value to be analyzed.
!
!    Output, integer PHIN, the value of PHI(N).  If N is less than
!    or equal to 0, PHI will be returned as 0.  If there is not enough
!    room for full factoring of N, PHI will be returned as -1.
!
  integer, parameter :: maxfactor = 20
!
  integer factor(maxfactor)
  integer i
  integer n
  integer nfactor
  integer nleft
  integer phin
  integer power(maxfactor)
!
  if ( n <= 0 ) then
    phin = 0
    return
  end if

  if ( n == 1 ) then
    phin = 1
    return
  end if
!
!  Factor N.
!
  call i_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_PHI - Fatal error!'
    write ( *, * ) '  Not enough factorization space.'
    phin = - 1
    return
  end if

  phin = 1
  do i = 1, nfactor
    phin = phin * factor(i)**( power(i) - 1 ) * ( factor(i) - 1 )
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
  integer i
  integer ihi
  integer ilo
  real r
  real, parameter :: rhi = 1.0E+00
  real, parameter :: rlo = 0.0E+00
!
  call r_random ( rlo, rhi, r )

  i = ilo + int ( r * real ( ihi + 1 - ilo ) )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
function i_roundup ( r )
!
!*******************************************************************************
!
!! I_ROUNDUP rounds a real value "up" to the nearest integer.
!
!
!  Examples:
!
!    R     I_ROUNDUP
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the real value to be rounded up.
!
!    Output, integer I_ROUNDUP, the rounded value.
!
  integer i_roundup
  integer itemp
  real r
!
  itemp = int ( r )
  if ( real ( itemp ) < r ) then
    itemp = itemp + 1
  end if

  i_roundup = itemp

  return
end
subroutine i_sigma ( n, sigma_n )
!
!*******************************************************************************
!
!! I_SIGMA returns the value of SIGMA(N), the divisor sum.
!
!
!  Definition:
!
!    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
!
!  First values:
!
!     N  SIGMA(N)
!
!     1    1
!     2    3
!     3    4
!     4    7
!     5    6
!     6   12
!     7    8
!     8   15
!     9   13
!    10   18
!    11   12
!    12   28
!    13   14
!    14   24
!    15   24
!    16   31
!    17   18
!    18   39
!    19   20
!    20   42
!
!  Formula:
!
!    SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
!
!    SIGMA(P**K) = ( P**(K+1) - 1 ) / ( P - 1 ) if P is prime.
!
!  Modified:
!
!    01 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the value to be analyzed.
!
!    Output, integer SIGMA_N, the value of SIGMA(N).  If N is less than
!    or equal to 0, SIGMA_N will be returned as 0.  If there is not
!    enough room for factoring N, SIGMA_N is returned as -1.
!
  integer, parameter :: maxfactor = 20
!
  integer factor(maxfactor)
  integer i
  integer n
  integer nfactor
  integer nleft
  integer power(maxfactor)
  integer sigma_n
!
  if ( n <= 0 ) then
    sigma_n = 0
    return
  end if

  if ( n == 1 ) then
    sigma_n = 1
    return
  end if
!
!  Factor N.
!
  call i_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_SIGMA - Fatal error!'
    write ( *, * ) '  Not enough factorization space.'
    sigma_n = - 1
    return
  end if

  sigma_n = 1
  do i = 1, nfactor
    sigma_n = ( sigma_n * ( factor(i)**( power(i) + 1 ) - 1 ) ) &
      / ( factor(i) - 1 )
  end do

  return
end
function i_sign ( x )
!
!*******************************************************************************
!
!! I_SIGN evaluates the sign of an integer.
!
!
!  Discussion:
!
!    This function differs from the intrinsic SIGN function, because
!    it returns a value of 0 if the input argument is 0.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the number whose sign is desired.
!
!    Output, integer I_SIGN, a result based on the sign of X:
!
!    -1, if X < 0.
!     0, if X = 0.
!    +1, if X > 0.
!
  integer i_sign
  integer x
!
  if ( x < 0 ) then
    i_sign = - 1
  else if ( x == 0 ) then
    i_sign = 0
  else if ( x > 0 ) then
    i_sign = + 1
  end if

  return
end
subroutine i_swap ( i, j )
!
!*******************************************************************************
!
!! I_SWAP swaps two integer values.
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
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  integer i
  integer j
  integer k
!
  k = i
  i = j
  j = k

  return
end
subroutine i_swap3 ( i, j, k )
!
!*******************************************************************************
!
!! I_SWAP3 swaps three integer values.
!
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J, K.  On output, the values of I, J, and K
!    have been interchanged.
!
  integer i
  integer j
  integer k
  integer l
!
  l = i
  i = j
  j = k
  k = l

  return
end
subroutine i_tau ( n, taun )
!
!*******************************************************************************
!
!! I_TAU returns the value of TAU(N), the number of distinct divisors of N.
!
!
!  Definition:
!
!    TAU(N) is the number of divisors of N, including 1 and N.
!
!  First values:
!
!     N   TAU(N)
!
!     1    1
!     2    2
!     3    2
!     4    3
!     5    2
!     6    4
!     7    2
!     8    4
!     9    3
!    10    4
!    11    2
!    12    6
!    13    2
!    14    4
!    15    4
!    16    5
!    17    2
!    18    6
!    19    2
!    20    6
!
!  Formula:
!
!    If the prime factorization of N is
!
!      N = P1**E1 * P2**E2 * ... * PM**EM,
!
!    then
!
!      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
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
!    Input, integer N, the value to be analyzed.  N must be 1 or
!    greater.
!
!    Output, integer TAUN, the value of TAU(N).  But if N is 0 or
!    less, TAUN is returned as 0, a nonsense value.  If there is
!    not enough room for factoring, TAUN is returned as -1.
!
  integer, parameter :: maxfactor = 20
!
  integer factor(maxfactor)
  integer i
  integer n
  integer nfactor
  integer nleft
  integer power(maxfactor)
  integer taun
!
  if ( n <= 0 ) then
    taun = 0
    return
  end if

  if ( n == 1 ) then
    taun = 1
    return
  end if
!
!  Factor N.
!
  call i_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_TAU - Fatal error!'
    write ( *, * ) '  Not enough factorization space.'
    taun = - 1
    return
  end if

  taun = 1
  do i = 1, nfactor
    taun = taun * ( power(i) + 1 )
  end do

  return
end
subroutine i_to_fac ( intval, nprime, npower )
!
!*******************************************************************************
!
!! I_TO_FAC converts an integer into a product of prime factors.
!
!
!  Discussion:
!
!    This routine will fail if the input integer is not positive,
!    or if NPRIME is too small to account for the factors of the integer.
!
!  Formula:
!
!    INTVAL = Product ( I = 1 to NPRIME ) PRIME(I)**NPOWER(I).
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INTVAL, the integer to be factored.
!
!    Input, integer NPRIME, the number of prime factors for
!    which storage has been allocated.
!
!    Output, integer NPOWER(NPRIME), the powers of the primes.
!
  integer nprime
!
  integer i
  integer intcopy
  integer intval
  integer npower(nprime)
  integer p
  integer prime
!
  if ( intval <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_TO_FAC - Fatal error!'
    write ( *, * ) '  Input integer is not positive.'
    stop
  end if
!
!  Try dividing the remainder by each prime.
!
  intcopy = intval

  do i = 1, nprime

    npower(i) = 0

    p = prime ( i )

    do while ( mod ( intcopy, p ) == 0 )
      npower(i) = npower(i) + 1
      intcopy = intcopy / p
    end do

  end do

  return
end
function i_to_isbn ( i )
!
!*******************************************************************************
!
!! I_TO_ISBN converts an integer to an ISBN digit.
!
!
!  Discussion:
!
!    Only the integers 0 through 10 can be input.  The representation
!    of 10 is 'X'.
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, an integer between 0 and 10.
!
!    Output, character I_TO_ISBN, the ISBN character code of the integer.
!    If I is illegal, then I_TO_ISBN is set to '?'.
!
  integer i
  character i_to_isbn
!
       if ( i == 0 ) then
    i_to_isbn = '0'
  else if ( i == 1 ) then
    i_to_isbn = '1'
  else if ( i == 2 ) then
    i_to_isbn = '2'
  else if ( i == 3 ) then
    i_to_isbn = '3'
  else if ( i == 4 ) then
    i_to_isbn = '4'
  else if ( i == 5 ) then
    i_to_isbn = '5'
  else if ( i == 6 ) then
    i_to_isbn = '6'
  else if ( i == 7 ) then
    i_to_isbn = '7'
  else if ( i == 8 ) then
    i_to_isbn = '8'
  else if ( i == 9 ) then
    i_to_isbn = '9'
  else if ( i == 10 ) then
    i_to_isbn = 'X'
  else
    i_to_isbn = '?'
  end if

  return
end
subroutine i_to_ivec_binary ( i, ivec )
!
!*******************************************************************************
!
!! I_TO_IVEC_BINARY makes a vector binary representation of an integer.
!
!
!  Example:
!
!     I       IVEC
!    --  ----------------
!     1  '0,0,...0,0,0,1'
!     2  '0,0,...0,0,1,0'
!     3  '0,0,...0,0,1,1'
!     4  '0,0,...0,1,0,0'
!    -9  '0,0,...1,0,0,1'
!
!  Discussion:
!
!    The vector is assumed to have dimension 32.  Only the absolute value
!    of the input integer is considered.
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
!    Input, integer I, an integer to be represented.
!
!    Output, integer IVEC(32), the binary representation.
!
  integer i
  integer i_copy
  integer ivec(32)
  integer j
!
  i_copy = abs ( i )

  do j = 32, 1, -1

    ivec(j) = mod ( i_copy, 2 )

    i_copy = i_copy / 2

  end do

  return
end
subroutine i_unswap3 ( i, j, k )
!
!*******************************************************************************
!
!! I_UNSWAP3 unswaps three integer values.
!
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J, K.  On output, the values of I, J, and K
!    have been interchanged.
!
  integer i
  integer j
  integer k
  integer l
!
  l = k
  k = j
  j = i
  i = l

  return
end
function i_wrap ( ival, ilo, ihi )
!
!*******************************************************************************
!
!! I_WRAP forces an integer to lie between given limits by wrapping.
!
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer I_WRAP, a "wrapped" version of IVAL.
!
  integer i_modp
  integer i_wrap
  integer ihi
  integer ilo
  integer ival
  integer wide
!
  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i_wrap = ilo
  else
    i_wrap = ilo + i_modp ( ival-ilo, wide )
  end if

  return
end
subroutine icol_compare ( lda, m, n, a, i, j, isgn )
!
!*******************************************************************************
!
!! ICOL_COMPARE compares columns I and J of a integer array.
!
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array, which must
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), an array of N columns of vectors of length M.
!
!    Input, integer I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column I > column J.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer isgn
  integer j
  integer k
  integer m
!
!  Check.
!
  if ( i < 1 .or. i > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'ICOL_COMPARE - Fatal error!'
    write ( *, * ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. j > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'ICOL_COMPARE - Fatal error!'
    write ( *, * ) '  Column index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = - 1
      return
    else if ( a(k,i) > a(k,j) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine icol_find ( lda, m, n, a, ivec, icol )
!
!*******************************************************************************
!
!! ICOL_FIND seeks a table column equal to an integer vector.
!
!
!  Example:
!
!    M = 3, N = 4,
!
!    A = (
!      1  2  3  4
!      5  6  7  8
!      9 10 11 12 )
!
!    IVEC = ( 3, 7, 11 )
!
!    ICOL = 3
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns in
!    the table.  M is also the length of IVEC.
!
!    Input, integer A(LDA,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer IVEC(M), a vector to be matched with the data
!    in the array.
!
!    Output, integer ICOL, the index of the first column of the table
!    which exactly matches every entry of IVEC, or 0 if no match
!    could be found.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer i
  integer icol
  integer ivec(m)
  integer j
!
  if ( m <= 0 ) then
    icol = 0
    return
  end if

  do j = 1, n

    i = 1

    do while ( ivec(i) == a(i,j) )

      if ( i == m ) then
        icol = j
        return
      end if

      i = i + 1

    end do

  end do

  icol = 0

  return
end
subroutine icol_find_item ( lda, m, n, a, item, irow, icol )
!
!*******************************************************************************
!
!! ICOL_FIND_ITEM searches a table by columns for a given value.
!
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns in
!    the table.  M is also the length of IVEC.
!
!    Input, integer A(LDA,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ITEM, the value to search for.
!
!    Output, integer IROW, ICOL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by rows.  If the item is not found, then
!    IROW = ICOL = 0.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer icol
  integer irow
  integer item
  integer j
  integer m
!
  do j = 1, n
    do i = 1, m
      if ( a(i,j) == item ) then
        irow = i
        icol = j
        return
      end if
    end do
  end do

  irow = 0
  icol = 0

  return
end
subroutine icol_find_pair_wrap ( lda, m, n, a, item1, item2, irow, icol )
!
!*******************************************************************************
!
!! ICOL_FIND_PAIR_WRAP searches a table by columns for a pair of items.
!
!
!  Discussion:
!
!    The items must occur consecutively, with ITEM1 occurring
!    first.  However, wrapping is allowed.  That is, if ITEM1
!    occurs in the last row, and ITEM2 in the first, this
!    is also regarded as a match.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A(LDA,N), the array to search.
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, integer ITEM1, ITEM2, the values to search for.
!
!    Output, integer IROW, ICOL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.  The search is conducted by columns.  If the pair of
!    items is not found, then IROW = ICOL = 0.  If IROW = M,
!    the ITEM1 occurs in row M and ITEM2 occurs in row 1.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer icol
  integer ip1
  integer irow
  integer item1
  integer item2
  integer j
  integer m
!
  do j = 1, n
    do i = 1, m

      if ( a(i,j) == item1 ) then

        if ( i < m ) then
          ip1 = i + 1
        else
          ip1 = 1
        end if

        if ( a(ip1,j) == item2 ) then
          irow = i
          icol = j
          return
        end if

      end if

    end do
  end do

  irow = 0
  icol = 0

  return
end
subroutine icol_sort_a ( lda, m, n, a )
!
!*******************************************************************************
!
!! ICOL_SORT_A ascending sorts an integer array of columns.
!
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(LDA,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer indx
  integer isgn
  integer j
  integer m
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call icol_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call icol_compare ( lda, m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine icol_sort_d ( lda, m, n, a )
!
!*******************************************************************************
!
!! ICOL_SORT_D descending sorts an integer array of columns.
!
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(LDA,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in descending
!    lexicographic order.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer indx
  integer isgn
  integer j
  integer m
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call icol_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call icol_compare ( lda, m, n, a, i, j, isgn )
      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine icol_swap ( lda, m, n, a, i, j )
!
!*******************************************************************************
!
!! ICOL_SWAP swaps columns I and J of a integer array of column data.
!
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Modified:
!
!    19 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, integer A(LDA,N), an array of N columns of length M.
!
!    Input, integer I, J, the columns to be swapped.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer j
  integer k
  integer m
!
  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, * ) ' '
    write ( *, * ) 'ICOL_SWAP - Fatal error!'
    write ( *, * ) '  I or J is out of bounds.'
    write ( *, * ) '  I =    ', i
    write ( *, * ) '  J =    ', j
    write ( *, * ) '  N =    ', n
    stop

  end if

  do k = 1, m
    call i_swap ( a(k,i), a(k,j) )
  end do

  return
end
subroutine icol_uniq ( lda, m, n, a, nuniq )
!
!*******************************************************************************
!
!! ICOL_UNIQ keeps the unique elements in a sorted integer array of columns.
!
!
!  Discussion:
!
!    The array can be sorted into ascending or descending order.
!    The important point is that identical elements must be stored
!    in adjacent positions.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array, which
!    must be at least M.
!
!    Input, integer M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, real A(LDA,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of NUNIQ columns of M-vectors.
!
!    Output, integer NUNIQ, the number of unique columns of A.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer isgn
  integer itest
  integer m
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    call icol_compare ( lda, m, n, a, itest, nuniq, isgn )

    if ( isgn /= 0 ) then
      nuniq = nuniq + 1
      a(1:m,nuniq) = a(1:m,itest)
    end if

  end do

  return
end
subroutine iint_to_rint ( imin, imax, i, rmin, rmax, r )
!
!*******************************************************************************
!
!! IINT_TO_RINT maps an integer interval to a real interval.
!
!
!  Formula:
!
!    R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IMIN, IMAX, the range.
!
!    Input, integer I, the integer to be converted.
!
!    Input, real RMIN, RMAX, the real range.
!
!    Output, real R, the corresponding value in [RMIN,RMAX].
!
  integer i
  integer imax
  integer imin
  real r
  real rmax
  real rmin
!
  if ( imax == imin ) then

    r = 0.5E+00 * ( rmin + rmax )

  else

    r = ( real ( imax - i ) * rmin + real ( i - imin ) * rmax ) &
      / real ( imax - imin )

  end if

  return
end
subroutine ij_next ( i, j, n )
!
!*******************************************************************************
!
!! IJ_NEXT returns the next matrix index.
!
!
!  Discussion:
!
!    For N = 3, the sequence of indices returned is:
!
!      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
!
!    Note that once the value (N,N) is returned, the next value returned
!    will be (0,0).
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On input, the current pair of indices.
!    On output, the next pair of indices.  If either index is illegal on
!    input, the output value of (I,J) will be (1,1).
!
!    Input, integer N, the maximum value for I and J.
!
  integer i
  integer j
  integer n
!
  if ( n < 1 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. i > n .or. j < 1 .or. j > n ) then
    i = 1
    j = 1
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n ) then
    i = i + 1
    j = 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine ij_next_gt ( i, j, n )
!
!*******************************************************************************
!
!! IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
!
!
!  Discussion:
!
!    For N = 5, the sequence of indices returned is:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On input, the current pair of indices.
!    On output, the next pair of indices.  If either index is illegal on
!    input, the output value of (I,J) will be (1,2).
!
!    Input, integer N, the maximum value for I and J.
!    A value of N less than 2 is nonsense.
!
  integer i
  integer j
  integer n
!
  if ( n < 2 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. i > n .or. j < 1 .or. j > n .or. j <= i ) then
    i = 1
    j = 2
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n - 1 ) then
    i = i + 1
    j = i + 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine imat_copy ( lda, m, n, a1, a2 )
!
!*******************************************************************************
!
!! IMAT_COPY copies a rectangular integer matrix.
!
!
!  Modified:
!
!    17 November 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading first dimension of A.
!    LDA must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A1(LDA,N), the M by N matrix to be copied.
!
!    Input, integer A2(LDA,N), a copy of A1.
!
  integer lda
  integer n
!
  integer a1(lda,n)
  integer a2(lda,n)
  integer m
!
  a2(1:m,1:n) = a1(1:m,1:n)

  return
end
subroutine imat_elim ( lda, m, n, a )
!
!*******************************************************************************
!
!! IMAT_ELIM carries out exact Gauss elimination on an integer matrix.
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
!    Input, integer LDA, the leading dimension of A, which must be
!    at least M.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input/output, integer A(LDA,N).  On input, the M by N matrix to
!    be Gauss eliminated.  On output, the Gauss-eliminated matrix.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer amax
  integer i
  integer icol(n)
  integer ifact
  integer i_gcd
  integer imax
  integer imult
  integer irow(m)
  integer iswap
  integer j
  integer jcol
  integer jmult
!
!  Initialize the swap parity counter.
!
  iswap = 1
!
!  For each column JCOL...
!
  do jcol = 1, min ( m, n )
!
!  Find the maximum element in rows JCOL through M.
!
    amax = abs ( a(jcol,jcol) )
    imax = jcol

    do i = jcol+1, m
      if ( abs ( a(i,jcol) ) > amax ) then
        amax = abs ( a(i,jcol) )
        imax = i
      end if
    end do
!
!  If the maximum entry is nonzero, then...
!
    if ( amax /= 0 ) then
!
!  If the maximum entry does not occur in row JCOL, then swap rows.
!
      if ( imax /= jcol ) then
        iswap = - iswap
        do j = 1, n
          call i_swap ( a(jcol,j), a(imax,j) )
        end do
      end if
!
!  Eliminate all nonzero entries in column JCOL, below the diagonal entry.
!
      do i = jcol+1, m

        if ( a(i,jcol) /= 0 ) then

          jmult = a(i,jcol)
          imult = a(jcol,jcol)
          ifact = i_gcd ( imult, jmult )
          imult = imult / ifact
          jmult = jmult / ifact

          do j = jcol, n
            a(i,j) = jmult * a(jcol,j) - imult * a(i,j)
          end do

        end if

      end do
!
!  Remove any row or column factors.
!
      call imat_red ( lda, m, n, a, irow, icol )

    end if

  end do

  return
end
subroutine imat_imax ( lda, m, n, a, i, j )
!
!*******************************************************************************
!
!! IMAT_IMAX returns the location of the maximum of an integer M by N matrix.
!
!
!  Modified:
!
!    17 July 2000
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
!    Input, integer A(LDA,N), the M by N matrix.
!
!    Output, integer I, J, the indices of the maximum entry of A.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer ii
  integer j
  integer jj
  integer m
!
  i = 0
  j = 0

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(ii,jj) > a(i,j) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
subroutine imat_imin ( lda, m, n, a, i, j )
!
!*******************************************************************************
!
!! IMAT_IMIN returns the location of the minimum of an integer M by N matrix.
!
!
!  Modified:
!
!    17 July 2000
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
!    Input, integer A(LDA,N), the M by N matrix.
!
!    Output, integer I, J, the indices of the minimum entry of A.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer ii
  integer j
  integer jj
  integer m
!
  i = 0
  j = 0

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(ii,jj) < a(i,j) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
subroutine imat_l1_inverse ( lda, n, a, b )
!
!*******************************************************************************
!
!! IMAT_L1_INVERSE inverts an integer unit lower triangular matrix.
!
!
!  Discussion:
!
!    A unit lower triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's above the main diagonal.
!
!    The inverse of an integer unit lower triangular matrix is also
!    an integer unit lower triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call imat_l1_inverse ( lda, n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    23 March 2000
!
!  Parameters:
!
!    Input, integer LDA, the declared first dimension of the arrays.
!    LDA must be at least N.
!
!    Input, integer N, number of rows and columns in the matrix.
!
!    Input, integer A(LDA,N), the unit lower triangular matrix.
!
!    Output, integer B(LDA,N), the inverse matrix.
!
  integer n
  integer lda
!
  integer a(lda,n)
  integer ab
  integer b(lda,n)
  integer i
  integer j
  integer k
!
  do i = 1, n

    do j = 1, n

      if ( j > i ) then
        b(i,j) = 0
      else if ( j == i ) then
        b(i,j) = 1
      else
        ab = 0
        do k = 1, i-1
          ab = ab + a(i,k) * b(k,j)
        end do
        b(i,j) = - ab
      end if

    end do
  end do

  return
end
function imat_max ( lda, m, n, a )
!
!*******************************************************************************
!
!! IMAT_MAX returns the maximum entry of an integer M by N matrix.
!
!
!  Modified:
!
!    17 January 2000
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
!    Input, integer A(LDA,N), the M by N matrix.
!
!    Output, integer IMAT_MAX, the maximum entry of A.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer imat_max
  integer m
!
  imat_max = maxval ( a(1:m,1:n) )

  return
end
function imat_min ( lda, m, n, a )
!
!*******************************************************************************
!
!! IMAT_MIN returns the minimum entry of an integer M by N matrix.
!
!
!  Modified:
!
!    17 January 2000
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
!    Input, integer A(LDA,N), the M by N matrix.
!
!    Output, integer IMAT_MIN, the minimum entry of A.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer imat_min
  integer m
!
  imat_min = minval ( a(1:m,1:n) )

  return
end
subroutine imat_perm ( lda, n, a, p )
!
!*******************************************************************************
!
!! IMAT_PERM permutes the rows and columns of a square integer matrix.
!
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    27 July 2000
!
!  Parameters:
!
!    Input, integer LDA, the leading first dimension of A.
!    LDA must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, integer A(LDA,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer P(N), the permutation.  P(I) is the new number of
!    row and column I.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer i1
  integer ierror
  integer is
  integer it
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer nc
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IMAT_PERM - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = - p(i)

    if ( i1 > 0 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine imat_perm2 ( lda, m, n, a, p, q )
!
!*******************************************************************************
!
!! IMAT_PERM2 permutes the rows and columns of a rectangular integer matrix.
!
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    28 October 1999
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!    LDA must be at least M.
!
!    Input, integer M, number of rows in the matrix.
!
!    Input, integer N, number of columns in the matrix.
!
!    Input/output, integer A(LDA,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer P(M), the row permutation.  P(I) is the new number of row I.
!
!    Input, integer Q(N), the column permutation.  Q(I) is the new number of column I.
!    Note that this routine allows you to pass a single array as both
!    P and Q.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer i
  integer i1
  integer ierror
  integer is
  integer it
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer nc
  integer p(m)
  integer q(n)
!
  call perm_check ( m, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IMAT_PERM2 - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
  call perm_check ( n, q, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IMAT_PERM2 - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
  call perm_cycle ( m, p, is, nc, 1 )

  if ( q(1) > 0 ) then
    call perm_cycle ( n, q, is, nc, 1 )
  end if

  do i = 1, m

    i1 = - p(i)

    if ( i1 > 0 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( q(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( q(j1) )

              call i_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                cycle
              end if

            end do

            j2 = abs ( q(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:m) = abs ( p(1:m) )

  if ( q(1) <= 0 ) then
    q(1:n) = abs ( q(1:n) )
  end if

  return
end
subroutine imat_perm2_random ( lda, m, n, a )
!
!*******************************************************************************
!
!! IMAT_PERM2_RANDOM selects a random permutation of an integer matrix.
!
!
!  Discussion:
!
!    The matrix may be rectangular.  Separate permutations are
!    applied to the rows and columns.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer LDA, the leading dimension of the array, which
!    must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, integer A(LDA,N), the M by N array to be permuted.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer i2
  integer j
  integer j2
  integer m
  real r
!
!  Permute the rows.
!
  do i = 1, m

    call i_random ( i, m, i2 )

    do j = 1, n
      call i_swap ( a(i2,j), a(i,j) )
    end do

  end do
!
!  Permute the columns.
!
  do j = 1, n

    call i_random ( j, n, j2 )

    do i = 1, m
      call i_swap ( a(i,j2), a(i,j) )
    end do

  end do

  return
end
subroutine imat_perm_random ( lda, n, a )
!
!*******************************************************************************
!
!! IMAT_PERM_RANDOM selects a random permutation of an integer matrix.
!
!
!  Discussion:
!
!    The matrix is assumed to be square.  A single permutation is
!    applied to both rows and columns.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least N.
!
!    Input, integer N, the number of rows and columns in the array.
!
!    Input/output, integer A(LDA,N), the N by N array to be permuted.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer i2
  integer j
  real r
!
!  Permute the rows and columns together.
!
  do i = 1, n

    call i_random ( i, n, i2 )

    i2 = i + int ( r * ( n + 1 - i ) )

    do j = 1, n
      call i_swap ( a(i2,j), a(i,j) )
    end do

    do j = 1, n
      call i_swap ( a(j,i2), a(j,i) )
    end do

  end do

  return
end
subroutine imat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! IMAT_PRINT prints an integer matrix.
!
!
!  Modified:
!
!    08 May 2000
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
!    Input, integer A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer lda
  integer n
!
  integer a(lda,n)
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

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, * ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(i6,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine imat_random ( alo, ahi, lda, m, n, a )
!
!*******************************************************************************
!
!! IMAT_RANDOM returns a matrix of uniform random values between AHI and ALO.
!
!
!  Formula:
!
!    A(I,J) = ALO + (AHI-ALO) * Random_Number
!
!  Modified:
!
!    23 September 2000
!
!  Parameters:
!
!    Input, integer ALO, AHI, the minimum and maximum values that
!    the matrix entries can have.
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, N, the number of rows and columns of A.
!
!    Output, integer A(LDA,N), the random matrix.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer ahi
  integer alo
  integer i
  integer j
  integer m
!
  do i = 1, m
    do j = 1, n
      call i_random ( alo, ahi, a(i,j) )
    end do
  end do

  return
end
subroutine imat_red ( lda, m, n, a, irow, icol )
!
!*******************************************************************************
!
!! IMAT_RED divides out common factors in a row or column of a matrix.
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
!    Input, integer LDA, the leading dimension of the matrix.
!    LDA must be at least M.
!
!    Input, integer M, the number of rows in the matrix.
!
!    Input, integer N, the number of columns in the matrix.
!
!    Input/output, integer A(LDA,N), on input, the M by N matrix to be reduced.
!    On output, A has been reduced.  The greatest common factor in any
!    row or column is 1.
!
!    Output, integer IROW(M), the row factors that were divided out.
!
!    Output, integer ICOL(N), the column factors that were divided out.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer i
  integer icol(n)
  integer ifact
  integer incx
  integer irow(m)
  integer j
!
  if ( m <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IMAT_RED - Warning!'
    write ( *, * ) '  Called with LDA, M, N = ', lda, m, n
    return
  end if

  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IMAT_RED - Warning!'
    write ( *, * ) '  Called with LDA, M, N = ', lda, m, n
    return
  end if

  if ( lda < m ) then
    write ( *, * ) ' '
    write ( *, * ) 'IMAT_RED - Warning!'
    write ( *, * ) '  Called with LDA, M, N = ', lda, m, n
    return
  end if
!
!  Remove factors common to a column.
!
  incx = 1
  do j = 1, n
    call ivec_red ( m, a(1,j), incx, ifact )
    icol(j) = ifact
  end do
!
!  Remove factors common to a row.
!
  incx = lda
  do i = 1, m
    call ivec_red ( n, a(i,1), incx, ifact )
    irow(i) = ifact
  end do

  return
end
subroutine imat_row_compare ( lda, m, n, a1, a2, result )
!
!*******************************************************************************
!
!! IMAT_ROW_COMPARE compares two arrays of row vectors.
!
!
!  Discussion:
!
!    The arrays are compared by comparing the rows.
!    The rows are compared in the lexicographic sense.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the arrays.
!    LDA must be at least M.
!
!    Input, integer M, number of rows in the arrays.
!
!    Input, integer N, number of columns in the arrays.
!
!    Input, integer A1(LDA,N), A2(LDA,N), the arrays.
!
!    Output, integer RESULT:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  integer lda
  integer n
!
  integer a1(lda,n)
  integer a2(lda,n)
  integer i
  integer j
  integer m
  integer result
!
  result = 0

  do i = 1, m

    do j = 1, n

      if ( a1(i,j) < a2(i,j) ) then
        result = - 1
        return
      else if ( a1(i,j) > a2(i,j) ) then
        result = + 1
        return
      end if

    end do

  end do

  return
end
subroutine imat_u1_inverse ( lda, n, a, b )
!
!*******************************************************************************
!
!! IMAT_U1_INVERSE inverts an integer unit upper triangular matrix.
!
!
!  Discussion:
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.
!
!    The inverse of an integer unit upper triangular matrix is also
!    an integer unit upper triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call imat_u1_inverse ( lda, n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    23 March 2000
!
!  Parameters:
!
!    Input, integer LDA, the declared first dimension of the arrays.
!    LDA must be at least N.
!
!    Input, integer N, number of rows and columns in the matrix.
!
!    Input, integer A(LDA,N), the unit upper triangular matrix.
!
!    Output, integer B(LDA,N), the inverse matrix.
!
  integer n
  integer lda
!
  integer a(lda,n)
  integer b(lda,n)
  integer i
  integer j
  integer k
  integer sum
!
  do j = n, 1, -1

    do i = n, 1, -1

      if ( i > j ) then
        b(i,j) = 0.0
      else if ( i == j ) then
        b(i,j) = 1
      else
        sum = 0
        do k = i + 1, j
          sum = sum + a(i,k) * b(k,j)
        end do
        b(i,j) = - sum
      end if

    end do
  end do

  return
end
subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )
!
!*******************************************************************************
!
!! INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
!
!
!  Discussion:
!
!    The box is has center at (IC,JC), and has half-widths N1 and N2.
!    The indices are exactly those which are between (IC-N1,JC-N1) and
!    (IC+N1,JC+N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the half-widths of the box, that is, the
!    maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer IC, JC, the central cell of the box.
!
!    Input/output, integer I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  integer i
  integer ic
  integer j
  integer jc
  logical more
  integer n1
  integer n2
!
  if ( more == .false. ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( j > jc + n2 ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. ( i == ic - n1 .or. i == ic + n1 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more )
!
!*******************************************************************************
!
!! INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
!
!
!  Discussion:
!
!    The box has a central cell of (IC,JC,KC), with a half widths of
!    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
!    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
!    is an "extreme" value.
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, the "half widths" of the box, that is, the
!    maximum distances from the central cell allowed for I, J and K.
!
!    Input, integer IC, JC, KC, the central cell of the box.
!
!    Input/output, integer I, J, K.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  integer i
  integer ic
  integer j
  integer jc
  integer k
  integer kc
  logical more
  integer n1
  integer n2
  integer n3
!
  if ( more == .false. ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    k = kc - n3
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 .and. k == kc + n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( k > kc + n3 ) then
    k = kc - n3
    j = j + 1
  else if ( k < kc + n3 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. j == jc - n2 .or. j == jc + n2 ) ) then
    return
  else
    k = kc + n3
    return
  end if
!
!  Check J.
!
  if ( j > jc + n2 ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. k == kc - n3 .or. k == kc + n3 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine irow_compare ( lda, m, n, a, i, j, isgn )
!
!*******************************************************************************
!
!! IROW_COMPARE compares two rows of a integer array.
!
!
!  Example:
!
!    Input:
!
!  M = 3, N = 4, I = 2, J = 3
!
!  A = (
!    1  2  3  4
!    5  6  7  8
!    9 10 11 12 )
!
!    Output:
!
!  ISGN = -1
!
!  Modified:
!
!    14 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array, which must
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), an array of M rows of vectors of length N.
!
!    Input, integer I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row I > row J.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer isgn
  integer j
  integer k
  integer m
!
!  Check.
!
  if ( i < 1 .or. i > m ) then
    write ( *, * ) ' '
    write ( *, * ) 'IROW_COMPARE - Fatal error!'
    write ( *, * ) '  Row index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. j > m ) then
    write ( *, * ) ' '
    write ( *, * ) 'IROW_COMPARE - Fatal error!'
    write ( *, * ) '  Row index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = - 1
      return
    else if ( a(i,k) > a(j,k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine irow_find_item ( lda, m, n, a, item, irow, icol )
!
!*******************************************************************************
!
!! IROW_FIND_ITEM searches the rows of an integer array for a given value.
!
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), the table to search.
!
!    Input, integer ITEM, the value to search for.
!
!    Output, integer IROW, ICOL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by rows.  If the item is not found, then
!    IROW = ICOL = 0.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer icol
  integer irow
  integer item
  integer j
  integer m
!
  do i = 1, m
    do j = 1, n
      if ( a(i,j) == item ) then
        irow = i
        icol = j
        return
      end if
    end do
  end do

  irow = 0
  icol = 0

  return
end
subroutine irow_find_pair_wrap ( lda, m, n, a, item1, item2, irow, icol )
!
!*******************************************************************************
!
!! IROW_FIND_PAIR_WRAP searches the rows of an integer array for a pair of items.
!
!
!  Discussion:
!
!    The items must occur consecutively, with ITEM1 occurring
!    first.  However, wrapping is allowed.  That is, if ITEM1
!    occurs in the last column, and ITEM2 in the first, this
!    is also regarded as a match.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), the table to search.
!
!    Input, integer ITEM1, ITEM2, the values to search for.
!
!    Output, integer IROW, ICOL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.  The search is conducted by rows.  If the pair of
!    items is not found, then IROW = ICOL = 0.  If ICOL = N,
!    the ITEM1 occurs in column N and ITEM2 occurs in column 1.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer icol
  integer irow
  integer item1
  integer item2
  integer j
  integer jp1
  integer m
!
  do i = 1, m
    do j = 1, n

      if ( a(i,j) == item1 ) then

        if ( j < n ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if

        if ( a(i,jp1) == item2 ) then
          irow = i
          icol = j
          return
        end if

      end if

    end do
  end do

  irow = 0
  icol = 0

  return
end
subroutine irow_max ( lda, m, n, a, imax, amax )
!
!*******************************************************************************
!
!! IROW_MAX returns the maximums of the rows of an integer array.
!
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), the array to be examined.
!
!    Output, integer IMAX(M); IMAX(I) is the column of the array in which
!    the maximum for row I occurs.
!
!    Output, integer AMAX(M), the maximums of the rows of the array.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer amax(m)
  integer i
  integer imax(m)
  integer j
!
  do i = 1, m

    imax(i) = 1
    amax(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) > amax(i) ) then
        imax(i) = j
        amax(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine irow_mean ( lda, m, n, a, mean )
!
!*******************************************************************************
!
!! IROW_MEAN returns the means of the rows of an integer array.
!
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), an array of data.
!
!    Output, real MEAN(M), the mean or average value of each row of TABLE.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer i
  integer j
  real mean(m)
!
  do i = 1, m
    mean(i) = 0.0E+00
    do j = 1, n
      mean(i) = mean(i) + a(i,j)
    end do
    mean(i) = mean(i) / real ( n )
  end do

  return
end
subroutine irow_min ( lda, m, n, a, imin, amin )
!
!*******************************************************************************
!
!! IROW_MIN returns the minimums of the rows of an integer array.
!
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), the array to be examined.
!
!    Output, integer IMIN(M); IMIN(I) is the column in which
!    the minimum for row I occurs.
!
!    Output, integer AMIN(M), the minimums of the rows.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer amin(m)
  integer i
  integer imin(m)
  integer j
!
  do i = 1, m

    imin(i) = 1
    amin(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) < amin(i) ) then
        imin(i) = j
        amin(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine irow_sort2_d ( lda, m, n, a )
!
!*******************************************************************************
!
!! IROW_SORT2_D descending sorts the elements of each row of an integer array.
!
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer A(LDA,N).
!    On input, the array of M rows of N-vectors.
!    On output, the elements of each row of A have been sorted in descending
!    order.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer indx
  integer irow
  integer isgn
  integer j
  integer m
!
!  Initialize.
!
  do irow = 1, m

    i = 0
    indx = 0
    isgn = 0
    j = 0
!
!  Call the external heap sorter.
!
    do

      call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
      if ( indx > 0 ) then

        call i_swap ( a(irow,i), a(irow,j) )
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(irow,i) < a(irow,j) ) then
          isgn = + 1
        else
          isgn = - 1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine irow_sort_a ( lda, m, n, a )
!
!*******************************************************************************
!
!! IROW_SORT_A ascending sorts the rows of an integer array.
!
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer A(LDA,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in ascending
!    lexicographic order.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer indx
  integer isgn
  integer j
  integer m
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call irow_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call irow_compare ( lda, m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine irow_sort_d ( lda, m, n, a )
!
!*******************************************************************************
!
!! IROW_SORT_D descending sorts the rows of an integer array.
!
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, the number of rows and columns of A.
!
!    Input/output, integer A(LDA,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in descending
!    lexicographic order.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer indx
  integer isgn
  integer j
  integer m
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call irow_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call irow_compare ( lda, m, n, a, i, j, isgn )
      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine irow_sum ( lda, m, n, a, sum )
!
!*******************************************************************************
!
!! IROW_SUM returns the sums of the rows of an integer array.
!
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), an array of data.
!
!    Output, integer SUM(M), the sum of the entries of each row of TABLE.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer i
  integer j
  integer sum(m)
!
  do i = 1, m
    sum(i) = 0
    do j = 1, n
      sum(i) = sum(i) + a(i,j)
    end do
  end do

  return
end
subroutine irow_swap ( lda, m, n, a, irow1, irow2 )
!
!*******************************************************************************
!
!! IROW_SWAP swaps two rows of an integer array.
!
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, integer A(LDA,N), an array of data.
!
!    Input, integer IROW1, IROW2, the two rows to swap.
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer m
  integer irow1
  integer irow2
  integer j
!
!  Check.
!
  if ( irow1 < 1 .or. irow1 > m ) then
    write ( *, * ) ' '
    write ( *, * ) 'IROW_SWAP - Fatal error!'
    write ( *, * ) '  IROW1 is out of range.'
    stop
  end if

  if ( irow2 < 1 .or. irow2 > m ) then
    write ( *, * ) ' '
    write ( *, * ) 'IROW_SWAP - Fatal error!'
    write ( *, * ) '  IROW2 is out of range.'
    stop
  end if

  if ( irow1 == irow2 ) then
    return
  end if

  do j = 1, n
    call i_swap ( a(irow1,j), a(irow2,j) )
  end do

  return
end
subroutine irow_variance ( lda, m, n, a, variance )
!
!*******************************************************************************
!
!! IROW_VARIANCE returns the variances of the rows of an integer array.
!
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), the array of data.
!
!    Output, real VARIANCE(M), the variance of each row.
!
  integer lda
  integer m
  integer n
!
  integer a(lda,n)
  integer i
  integer j
  real mean
  real variance(m)
!
  call irow_mean ( lda, m, n, a, variance )

  do i = 1, m

    mean = variance(i)

    variance(i) = 0.0E+00
    do j = 1, n
      variance(i) = variance(i) + ( real ( a(i,j) ) - mean )**2
    end do

    if ( n > 1 ) then
      variance(i) = variance(i) / real ( n - 1 )
    else
      variance(i) = 0.0E+00
    end if

  end do

  return
end
subroutine isbn_check ( isbn, check )
!
!*******************************************************************************
!
!! ISBN_CHECK checks an ISBN code.
!
!
!  Definition:
!
!    ISBN stands for International Standard Book Number.  A unique ISBN
!    is assigned to each new book.  The ISBN includes 10 digits.  There is
!    an initial digit, then a dash, then a set of digits which are a
!    code for the publisher, another digit, and then the check digit:
!
!      initial-publisher-book-check
!
!    The number of digits used for the publisher and book codes can vary,
!    but the check digit is always one digit, and the total number of
!    digits is always 10.
!
!    The check digit is interesting, because it is a way of trying to
!    make sure that an ISBN has not been incorrectly copied.  Specifically,
!    if the ISBN is correct, then its ten digits will satisfy
!
!       10 * A + 9 * B + 8 * C + 7 * D + 6 * E
!      + 5 * F * 4 * G * 3 * H + 2 * I +     J  = 0 mod 11.
!
!    Here, we've taken 'A' to represent the first digit and 'J' the
!    last (which is the check digit).  In order for the code to work,
!    the value of J must be allowed to be anything from 0 to 10.  In
!    cases where J works out to be 10, the special digit 'X' is used.
!    An 'X' digit can only occur in this last check-digit position
!    on an ISBN.
!
!  Example:
!
!    0-8493-9640-9
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ISBN, an ISBN code.
!
!    Output, integer CHECK, the value of the ISBN check sum.
!    If CHECK is zero, the ISBN code is legitimate.
!    If CHECK is -1, then the ISBN code is not legitimate because it does
!    not contain exactly 10 digits.  If CHECK is between 1 and 10, then
!    the ISBN code has the right number of digits, but at least one of
!    the digits is incorrect.
!
  character c
  logical ch_is_digit
  integer check
  integer digit(10)
  integer i
  character ( len = * ) isbn
  integer isbn_to_i
  integer lenc
  integer num_digit
!
!  Determine how many digits have been supplied.
!
  lenc = len_trim ( isbn )

  i = 0
  num_digit = 0

  do

    i = i + 1

    if ( i > lenc ) then
      exit
    end if

    c = isbn(i:i)

    if ( ch_is_digit ( c ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i ( c )

    else if ( ( num_digit == 9 .and. isbn(i:i) == 'X' ) .or. &
              ( num_digit == 9 .and. isbn(i:i) == 'x' ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i ( c )

    end if

    if ( num_digit >= 10 ) then
      exit
    end if

  end do
!
!  If we didn't get exactly 10 digits, return with an error.
!
  if ( num_digit /= 10 ) then
    check = -1
    return
  end if
!
!  Compute the checksum.
!
  check = 0
  do i = 1, 10
    check = check + ( 11 - i ) * digit(i)
  end do

  check = mod ( check, 11 )

  return
end
subroutine isbn_fill ( isbn )
!
!*******************************************************************************
!
!! ISBN_FILL fills in a missing digit in an ISBN code.
!
!
!  Example:
!
!    Input:
!
!      0-8493-9?40-9
!
!    Output:
!
!      0-8493-9640-9
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) ISBN, a partial ISBN code.  On input,
!    a single digit has been replaced by the character '?', signifying
!    that that digit is missing.  The routine replaces the question
!    mark by the correct digit.
!
  character c
  logical ch_is_digit
  integer check
  integer digit(10)
  integer digit_pos
  integer i
  character i_to_isbn
  character ( len = * ) isbn
  integer isbn_pos
  integer isbn_to_i
  integer j
  integer k
  integer lenc

  integer num_digit
!
  lenc = len_trim ( isbn )

  i = 0
  isbn_pos = -1
  digit_pos = -1
  num_digit = 0

  do

    i = i + 1

    if ( i > lenc ) then
      exit
    end if

    c = isbn(i:i)

    if ( ch_is_digit ( c ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i ( c )

    else if ( ( num_digit == 9 .and. isbn(i:i) == 'X' ) .or. &
              ( num_digit == 9 .and. isbn(i:i) == 'x' ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i ( c )

    else if ( c == '?' ) then

      if ( isbn_pos == -1 ) then

        num_digit = num_digit + 1
        digit(num_digit) = 0
        digit_pos = num_digit
        isbn_pos = i

      else
        write ( *, * ) ' '
        write ( *, * ) 'ISBN_FILL - Fatal error!'
        write ( *, * ) '  Only one question mark is allowed!'
        stop
      end if

    end if

    if ( num_digit >= 10 ) then
      exit
    end if

  end do

  if ( num_digit /= 10 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ISBN_FILL - Fatal error!'
    write ( *, * ) '  The input ISBN code did not have 10 digits.'
    stop
  end if

  if ( isbn_pos == -1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ISBN_FILL - Fatal error!'
    write ( *, * ) '  A question mark is required!'
    stop
  end if

  check = 0
  do i = 1, 10
    check = check + ( 11 - i ) * digit(i)
  end do

  check = mod ( check, 11 )

  if ( check == 0 ) then

    k = 0
!
!  Need to solve the modular equation:
!
!    A * X = B mod C
!
!  Below is a stupid way.  One day I will come back and fix this up.
!
  else

    do i = 1, 10
      j = ( 11 - digit_pos ) * i + check
      if ( mod ( j, 11 ) == 0 ) then
        k = i
      end if
    end do

  end if

  isbn(isbn_pos:isbn_pos) = i_to_isbn ( k )

  return
end
function isbn_to_i ( c )
!
!*******************************************************************************
!
!! ISBN_TO_I converts an ISBN character into an integer.
!
!
!  Discussion:
!
!    The characters '0' through '9' stand for themselves, but
!    the character 'X' or 'x' stands for 10.
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the ISBN character code to be converted.
!
!    Output, integer ISBN_TO_I, the numeric value of the character
!    code, between 0 and 10.  This value is returned as -1 if C is
!    not a valid character code.
!
  character c
  integer isbn_to_i
!
       if ( c == '0' ) then
    isbn_to_i = 0
  else if ( c == '1' ) then
    isbn_to_i = 1
  else if ( c == '2' ) then
    isbn_to_i = 2
  else if ( c == '3' ) then
    isbn_to_i = 3
  else if ( c == '4' ) then
    isbn_to_i = 4
  else if ( c == '5' ) then
    isbn_to_i = 5
  else if ( c == '6' ) then
    isbn_to_i = 6
  else if ( c == '7' ) then
    isbn_to_i = 7
  else if ( c == '8' ) then
    isbn_to_i = 8
  else if ( c == '9' ) then
    isbn_to_i = 9
  else if ( c == 'X' .or. c == 'x' ) then
    isbn_to_i = 10
  else
    isbn_to_i = -1
  end if

  return
end
function iset2_compare ( x1, y1, x2, y2 )
!
!*******************************************************************************
!
!! ISET2_COMPARE compares two I2 sets.
!
!
!  Discussion:
!
!    The I2 set (X1,Y1) < (X2,Y2) if
!
!      min ( X1, Y1 ) < min ( X2, Y2 ) or
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) < max ( X2, Y2 )
!
!    The I2 set (X1,Y1) = (X2,Y2) if
!
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) = max ( X2, Y2 )
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X1, Y1, the first I2 set.
!
!    Input, integer X2, Y2, the second I2 set.
!
!    Output, character ISET2_COMPARE: '<', '>' or '=' if the first I2 set
!    is less, greater or equal to the second.
!
  integer a1
  integer a2
  integer b1
  integer b2
  character c
  character iset2_compare
  integer x1
  integer x2
  integer y1
  integer y2
!
  a1 = min ( x1, y1 )
  b1 = max ( x1, y1 )

  a2 = min ( x2, y2 )
  b2 = max ( x2, y2 )

  if ( a1 < a2 ) then
    c = '<'
  else if ( a1 > a2 ) then
    c = '>'
  else if ( b1 < b2 ) then
    c = '<'
  else if ( b1 > b2 ) then
    c = '>'
  else
    c = '='
  end if

  iset2_compare = c

  return
end
subroutine iset2_index_insert_unique ( maxn, n, x, y, indx, &
  xval, yval, ival, ierror )
!
!*******************************************************************************
!
!! ISET2_INDEX_INSERT_UNIQUE inserts a unique I2 set value in an indexed sorted list.
!
!
!  Discussion:
!
!    If the input value does not occur in the list, then N, X, Y and INDX
!    are updated.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input/output, integer N, the size of the list.
!
!    Input/output, integer X(N), Y(N), the list of I2 sets.
!
!    Input/output, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer IERROR, 0 for no error, 1 if an error occurred.
!
  integer maxn
!
  integer equal
  integer ierror
  integer indx(maxn)
  integer ival
  integer less
  integer more
  integer n
  integer x(maxn)
  integer xval
  integer y(maxn)
  integer yval
!
  ierror = 0

  if ( n <= 0 ) then

    if ( maxn <= 0 ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, * ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = min ( xval, yval )
    y(1) = max ( xval, yval )
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL ) already occur in the list?
!
  call iset2_index_search ( maxn, n, x, y, indx, xval, yval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( n >= maxn ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, * ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = min ( xval, yval )
    y(n+1) = max ( xval, yval )
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine iset2_index_search ( maxn, n, x, y, indx, xval, yval, &
  less, equal, more )
!
!*******************************************************************************
!
!! ISET2_INDEX_SEARCH searches for an I2 set value in an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input, integer N, the size of the current list.
!
!    Input, integer X(N), Y(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, YVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    list entries that are just less than, equal to, and just greater
!    than the test value.  If the test value does not occur in the list,
!    then EQUAL is zero.  If the test value is the minimum entry of the
!    list, then LESS is 0.  If the test value is the greatest entry of
!    the list, then MORE is N+1.
!
  integer maxn
!
  character c
  integer equal
  integer hi
  integer indx(maxn)
  integer less
  integer lo
  integer mid
  integer more
  integer n
  character iset2_compare
  integer x(maxn)
  integer xhi
  integer xlo
  integer xmid
  integer xval
  integer y(maxn)
  integer yhi
  integer ylo
  integer ymid
  integer yval
!
  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))

  c = iset2_compare ( xval, yval, xlo, ylo )

  if ( c == '<' ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( c == '=' ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  c = iset2_compare ( xval, yval, xhi, yhi )

  if ( c == '>' ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( c == '=' ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))

    c = iset2_compare ( xval, yval, xmid, ymid )

    if ( c == '=' ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( c == '<' ) then
      hi = mid
    else if ( c == '>' ) then
      lo = mid
    end if

  end do

  return
end
subroutine ivec2_compare ( n, a1, a2, i, j, isgn )
!
!*******************************************************************************
!
!! IVEC2_COMPARE compares pairs of integers stored in two vectors.
!
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data items.
!
!    Input, integer A1(N), A2(N), contain the two components of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer i
  integer isgn
  integer j
!
  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(i) > a2(j) ) then
      isgn = +1
    end if

  else if ( a1(i) > a1(j) ) then

    isgn = +1

  end if

  return
end
subroutine ivec2_print ( n, a, b, title )
!
!*******************************************************************************
!
!! IVEC2_PRINT prints a pair of integer vectors.
!
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), B(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  integer a(n)
  integer b(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,2i10)' ) i, a(i), b(i)
  end do

  return
end
subroutine ivec2_sort_a ( n, a1, a2 )
!
!*******************************************************************************
!
!! IVEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), A2(N), the data to be sorted..
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call i_swap ( a1(i), a1(j) )
      call i_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call ivec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine ivec2_sort_d ( n, a1, a2 )
!
!*******************************************************************************
!
!! IVEC2_SORT_D descending sorts a vector of pairs of integers.
!
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), A2(N), the data to be sorted..
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call i_swap ( a1(i), a1(j) )
      call i_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call ivec2_compare ( n, a1, a2, i, j, isgn )
      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine ivec2_uniq ( n, a1, a2, nuniq )
!
!*******************************************************************************
!
!! IVEC2_UNIQ keeps the unique elements in a array of pairs of integers.
!
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items.
!
!    Input/output, integer A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer NUNIQ, the number of unique items.
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer itest
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a1(itest) /= a1(nuniq) .or. a2(itest) /= a2(nuniq) ) then

      nuniq = nuniq + 1

      a1(nuniq) = a1(itest)
      a2(nuniq) = a2(itest)

    end if

  end do

  return
end
subroutine ivec_amax ( n, a, aamax )
!
!*******************************************************************************
!
!! IVEC_AMAX returns the largest magnitude in an integer vector.
!
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector to be searched.
!
!    Output, integer AAMAX, the value of the entry of largest magnitude.
!
  integer n
!
  integer a(n)
  integer aamax
  integer i
!
  if ( n <= 0 ) then

    aamax = 0

  else

    aamax = abs ( a(1) )

    do i = 2, n
      aamax = max ( aamax, abs ( a(i) ) )
    end do

  end if

  return
end
subroutine ivec_amin ( n, a, aamin )
!
!*******************************************************************************
!
!! IVEC_AMIN returns the smallest magnitude in an integer vector.
!
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries to be checked.
!
!    Input, integer A(N), the vector to be checked.
!
!    Output, integer AAMIN, the value of the smallest magnitude.
!
  integer n
!
  integer a(n)
  integer aamin
  integer i
!
  if ( n <= 0 ) then

    aamin = 0

  else

    aamin = abs ( a(1) )

    do i = 2, n
      aamin = min ( aamin, abs ( a(i) ) )
    end do

  end if

  return
end
subroutine ivec_aminz ( n, a, aminz )
!
!*******************************************************************************
!
!! IVEC_AMINZ returns the smallest nonzero magnitude in an integer vector.
!
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries to be checked.
!
!    Input, integer A(N), the vector to be checked.
!
!    Output, integer AMINZ, the value of the smallest nonzero magnitude.
!    If all entries are zero, AMINZ is 0.
!
  integer n
!
  integer a(n)
  integer aminz
  integer i
  integer iset
!
  aminz = 0
  iset = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( iset == 0 ) then
        aminz = abs ( a(i) )
        iset = 1
      else
        aminz = min ( aminz, abs ( a(i) ) )
      end if

    end if

  end do

  return
end
subroutine ivec_axpy ( n, ia, x, incx, y, incy )
!
!*******************************************************************************
!
!! IVEC_AXPY:  Y(I) := Y(I) + A * X(I).
!
!
!  Discussion:
!
!    If X and Y are simple vectors, then IAXPY is equivalent to:
!
!      DO I = 1, N
!        Y(I) = Y(I) + IA * X(I)
!      END DO
!
!    However, by using the increments correctly, IAXPY can also be used
!    to manipulate rows or columns of matrices.
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
!    Input, integer N, the number of entries of X and Y.
!
!    Input, integer IA, the scalar value by which each entry
!    of X is multiplied before being added to Y.
!
!    Input, integer X(*), the vector, a multiple of which is to be
!    added to Y.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, integer Y(*).
!
!    On output, each entry of Y has been increased by
!    IA times the corresponding entry of X.
!
!    Input, integer INCY, the increment between successive entries of Y.
!
  integer i
  integer ia
  integer incx
  integer incy
  integer indx
  integer indy
  integer n
  integer x(*)
  integer y(*)
!
  indx = 1
  indy = 1

  do i = 1, n

    y(indy) = y(indy) + ia * x(indx)

    indx = indx + incx
    indy = indy + incy

  end do

  return
end
subroutine ivec_bracket ( n, a, xval, left, right )
!
!*******************************************************************************
!
!! IVEC_BRACKET searches a sorted array for successive brackets of a value.
!
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the number line, then this routine searches for the interval
!    containing the given value.
!
!  Modified:
!
!    25 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input, integer A(N), an array that has been sorted into ascending order.
!
!    Input, integer XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
!    and A(LEFT) <= XVAL <= A(RIGHT).
!
!    Special cases:
!      Value is less than all data values:
!        LEFT = -1, RIGHT = 1, and XVAL < A(RIGHT).
!      Value is greater than all data values:
!        LEFT = N, RIGHT = -1, and A(LEFT) < XVAL.
!      Value is equal to a data value:
!        LEFT = RIGHT, and A(LEFT) = A(RIGHT) = XVAL.
!
  integer n
!
  integer a(n)
  integer high
  integer left
  integer low
  integer mid
  integer right
  integer xval
!
!  XVAL < A(1).
!
  if ( xval < a(1) ) then
    left = -1
    right = 1
!
!  A(N) < XVAL.
!
  else if ( xval > a(n) ) then
    left = n
    right = -1
!
!  N = 1
!
  else if ( n == 1 ) then
    left = 1
    right = 1
!
!  A(1) <= XVAL <= A(N).
!
  else

    low = 1
    high = n - 1

    do

      mid = ( low + high ) / 2

      if ( low > high ) then
        write ( *, * ) ' '
        write ( *, * ) 'IVEC_BRACKET - Fatal error!'
        write ( *, * ) '  Algorithm or data failure.'
        stop
      end if

      if ( a(mid) == xval ) then
        left = mid
        right = mid
        exit
      else if ( a(mid+1) == xval ) then
        left = mid + 1
        right = mid + 1
        exit
      else if ( a(mid) < xval .and. xval < a(mid+1) ) then
        left = mid
        right = mid + 1
        exit
      else if ( a(mid+1) < xval ) then
        low = mid + 1
      else if ( a(mid) > xval ) then
        high = mid - 1
      end if

    end do

  end if

  return
end
subroutine ivec_compare ( n, a1, a2, isgn )
!
!*******************************************************************************
!
!! IVEC_COMPARE compares two integer vectors.
!
!
!  Discussion:
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2, 6, 2 )
!      A2 = ( 2, 8, 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, integer A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer isgn
  integer k
!
  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = - 1
      return
    else if ( a1(k) > a2(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine ivec_cum ( n, a, a_cum )
!
!*******************************************************************************
!
!! IVEC_CUM computes the cumulutive sum of the entries of a vector.
!
!
!  Example:
!
!    Input:
!
!      A = ( 1, 2, 3, 4 )
!
!    Output:
!
!      A_CUM = ( 0, 1, 3, 6, 10 )
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector to be summed.
!
!    Output, integer A_CUM(N+1), the cumulative sum of the entries of A.
!
  integer n
!
  integer a(n)
  integer a_cum(n+1)
  integer i
  integer sum
!
  sum = 0

  do i = 1, n
    a_cum(i) = sum
    sum = sum + a(i)
  end do

  a_cum(n+1) = sum

  return
end
function ivec_descends ( n, x )
!
!*******************************************************************************
!
!! IVEC_DESCENDS determines if an integer vector is decreasing.
!
!
!  Example:
!
!    X = ( 9, 7, 7, 3, 2, 1, -8 )
!
!    IVEC_DESCENDS = TRUE
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the array.
!
!    Input, integer X(N), the array to be examined.
!
!    Output, logical IVEC_DESCENDS, is TRUE if the entries of X descend.
!
  integer n
!
  integer i
  logical ivec_descends
  integer x(n)
!
  ivec_descends = .false.

  do i = 1, n-1
    if ( x(i) < x(i+1) ) then
      return
    end if
  end do

  ivec_descends = .true.

  return
end
subroutine ivec_frac ( n, a, k, iafrac )
!
!*******************************************************************************
!
!! IVEC_FRAC searches for the K-th smallest element in an N-vector.
!
!
!  Discussion:
!
!    Hoare's algorithm is used.
!
!  Modified:
!
!    17 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, integer A(N), array to search.  On output,
!    the elements of A have been somewhat rearranged.
!
!    Input, integer K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.
!    Other values of K search for the entry which is K-th in size.
!    K must be at least 1, and no greater than N.
!
!    Output, integer IAFRAC, the value of the K-th fractile of A.
!
  integer n
!
  integer i
  integer a(n)
  integer iafrac
  integer iryt
  integer iw
  integer ix
  integer j
  integer k
  integer left
!
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_FRAC  - Fatal error!'
    write ( *, * ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_FRAC  - Fatal error!'
    write ( *, * ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( k > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_FRAC  - Fatal error!'
    write ( *, * ) '  Illegal K > N, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( left >= iryt ) then
      iafrac = a(k)
      exit
    end if

    ix = a(k)
    i = left
    j = iryt

    do

      if ( i > j ) then

        if ( j < k ) then
          left = i
        end if

        if ( k < i ) then
          iryt = j
        end if

        exit

      end if
!
!  Find I so that IX <= A(I).
!
      do while ( a(i) < ix )
        i = i + 1
      end do
!
!  Find J so that A(J) <= IX.
!
      do while ( ix < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call i_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
function ivec_gcd ( n, a )
!
!*******************************************************************************
!
!! IVEC_GCD finds the greatest common divisor of an integer vector.
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
!    Input, integer N, the dimension of the vector.
!
!    Input, integer A(N), the vector to be checked.
!
!    Output, integer IVEC_GCD, the greatest common divisor of all entries of A.
!
  integer n
!
  integer a(n)
  integer i
  integer i_gcd
  integer ivec_gcd
!
  ivec_gcd = maxval ( abs ( a(1:n) ) )

  do i = 1, n
    ivec_gcd = i_gcd ( ivec_gcd, a(i) )  
  end do

  return
end
subroutine ivec_heap_a ( n, a )
!
!*******************************************************************************
!
!! IVEC_HEAP_A reorders an array of integers into an ascending heap.
!
!
!  Definition:
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the size of the input array.
!
!    Input/output, integer A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  integer n
!
  integer a(n)
  integer i
  integer ifree
  integer key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( m > n ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) >= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine ivec_heap_d ( n, a )
!
!*******************************************************************************
!
!! IVEC_HEAP_D reorders an array of integers into an descending heap.
!
!
!  Definition:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the size of the input array.
!
!    Input/output, integer A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  integer n
!
  integer a(n)
  integer i
  integer ifree
  integer key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( m > n ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m+1) > a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine ivec_heap_d_extract ( n, a, val )
!
!*******************************************************************************
!
!! IVEC_HEAP_D_EXTRACT extracts the maximum value from a descending heap.
!
!
!  Discussion:
!
!    In other words, the routine finds the maximum value in the
!    heap, returns that value to the user, deletes that value from
!    the heap, and restores the heap to its proper form.
!
!    This is one of three functions needed to model a priority queue.
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N, the number of items in the heap.
!
!    Input/output, integer A(N), the heap.
!
!    Output, integer VAL, the item of maximum value, which has been
!    removed from the heap.
!
  integer a(*)
  integer n
  integer val
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_HEAP_D_EXTRACT - Fatal error!'
    write ( *, * ) '  The heap is empty.'
    stop
  end if
!
!  Get the maximum value.
!
  val = a(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last value down.
!
  a(1) = a(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call ivec_sort_heap_d ( n, a )

  return
end
subroutine ivec_heap_d_insert ( n, a, val )
!
!*******************************************************************************
!
!! IVEC_HEAP_D_INSERT inserts a new value into a descending heap.
!
!
!  Discussion:
!
!    This is one of three functions needed to model a priority queue.
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N, the number of items in the heap.
!
!    Input/output, integer A(N), the heap.
!
!    Input, integer VAL, the value to be inserted.
!
  integer a(*)
  integer i
  integer n
  integer parent
  integer val
!
  n = n + 1
  i = n

  do while ( i > 1 )

    parent = i / 2

    if ( a(parent) >= val ) then
      exit
    end if

    a(i) = a(parent)
    i = parent

  end do

  a(i) = val

  return
end
subroutine ivec_heap_d_max ( n, a, val_max )
!
!*******************************************************************************
!
!! IVEC_HEAP_D_MAX returns the maximum value in a descending heap of integers.
!
!
!  Discussion:
!
!    This is one of three functions needed to model a priority queue.
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items in the heap.
!
!    Input, integer A(N), the heap.
!
!    Output, integer VAL_MAX, the maximum value in the heap.
!
  integer n
!
  integer a(n)
  integer val_max
!
  val_max = a(1)

  return
end
subroutine ivec_iamax ( n, a, iamax )
!
!*******************************************************************************
!
!! IVEC_IAMAX returns the index of the largest magnitude in an integer vector.
!
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector to be searched.
!
!    Output, integer IAMAX, the index of the entry of largest magnitude.
!
  integer n
!
  integer a(n)
  integer aamax
  integer i
  integer iamax
!
  if ( n <= 0 ) then

    iamax = 0

  else

    aamax = abs ( a(1) )
    iamax = 1

    do i = 2, n

      if ( abs ( a(i) ) > aamax ) then
        aamax = abs ( a(i) )
        iamax = i
      end if

    end do

  end if

  return
end
subroutine ivec_iamin ( n, a, iamin )
!
!*******************************************************************************
!
!! IVEC_IAMIN returns the index of the smallest magnitude in an integer vector.
!
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries to be checked.
!
!    Input, integer A(N), the vector to be checked.
!
!    Output, integer IAMIN, the entry of the smallest magnitude.
!
  integer n
!
  integer a(n)
  integer aamin
  integer i
  integer iamin
!
  if ( n <= 0 ) then

    iamin = 0

  else

    aamin = a(1)
    iamin = 1

    do i = 2, n

      if ( abs ( a(i) ) < aamin ) then
        aamin = abs ( a(i) )
        iamin = i
      end if

    end do

  end if

  return
end
subroutine ivec_iaminz ( n, a, iaminz )
!
!*******************************************************************************
!
!! IVEC_IAMINZ returns the smallest nonzero magnitude in an integer vector.
!
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries to be checked.
!
!    Input, integer A(N), the vector to be checked.
!
!    Output, integer IAMINZ, the entry of the smallest nonzero magnitude.
!    If all entries are zero, IAMINZ is 0.
!
  integer n
!
  integer a(n)
  integer aaminz
  integer i
  integer iaminz
!
  aaminz = 0
  iaminz = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( iaminz == 0 .or. abs ( a(i) ) < aaminz ) then
        aaminz = abs ( a(i) )
        iaminz = i
      end if

    end if

  end do

  return
end
subroutine ivec_identity ( n, a )
!
!*******************************************************************************
!
!! IVEC_IDENTITY sets an integer vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  integer n
!
  integer a(n)
  integer i
!
  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine ivec_imax ( n, a, imax )
!
!*******************************************************************************
!
!! IVEC_IMAX computes the index of a maximum element of an integer array.
!
!
!  Discussion:
!
!    If more than one element has the maximum value, this routine returns
!    the index of the first such element.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer IMAX, the index of the largest entry.
!
  integer n
!
  integer a(n)
  integer amax
  integer i
  integer imax
!
  if ( n <= 0 ) then

    imax = 0

  else

    amax = a(1)
    imax = 1

    do i = 2, n

      if ( a(i) > amax ) then
        amax = a(i)
        imax = i
      end if

    end do

  end if

  return
end
function ivec_imax_last ( n, x )
!
!*******************************************************************************
!
!! IVEC_MAX_LAST returns the index of the last maximal integer vector element.
!
!
!  Example:
!
!    X = ( 5, 1, 2, 5, 0, 5, 3 )
!
!    IVEC_IMAX_LAST = 6
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the array.
!
!    Input, integer X(N), the array to be examined.
!
!    Output, integer IVEC_IMAX_LAST, the index of the last element of
!    X of maximal value.
!
  integer n
!
  integer i
  integer ivec_imax_last
  integer max_last
  integer x(n)
!
  ivec_imax_last = 0

  do i = 1, n
    if ( i == 1 ) then
      ivec_imax_last = 1
      max_last = x(1)
    else if ( x(i) >= max_last ) then
      ivec_imax_last = i
      max_last = x(i)
    end if
  end do

  return
end
subroutine ivec_imin ( n, a, imin )
!
!*******************************************************************************
!
!! IVEC_IMIN computes the index of the minimum element of an integer array.
!
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer IMIN, the index of the smallest entry.
!
  integer n
!
  integer a(n)
  integer amin
  integer i
  integer imin
!
  if ( n <= 0 ) then

    imin = 0

  else

    amin = a(1)
    imin = 1

    do i = 2, n

      if ( a(i) < amin ) then
        amin = a(i)
        imin = i
      end if

    end do

  end if

  return
end
function ivec_index ( n, a, aval )
!
!*******************************************************************************
!
!! IVEC_INDEX returns the location of the first occurrence of a given value.
!
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector to be searched.
!
!    Input, integer AVAL, the value to be indexed.
!
!    Output, integer IVEC_INDEX, the first location in A which has the
!    value AVAL, or 0 if no such index exists.
!
  integer n
!
  integer a(n)
  integer aval
  integer i
  integer ivec_index
!
  do i = 1, n
    if ( a(i) == aval ) then
      ivec_index = i
      return
    end if
  end do

  ivec_index = 0

  return
end
subroutine ivec_index_delete_all ( n, x, indx, xval )
!
!*******************************************************************************
!
!! IVEC_INDEX_DELETE_ALL deletes all occurrences of an integer value from an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, integer X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer equal1
  integer equal2
  integer get
  integer i
  integer indx(*)
  integer j
  integer less
  integer more
  integer put
  integer x(*)
  integer xval
!
  if ( n < 1 ) then
    n = 0
    return
  end if

  call ivec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal /= 0 ) then

    equal1 = equal

    do

      if ( equal1 <= 1 ) then
        exit
      end if

      if ( x(indx(equal1-1)) /= xval ) then
        exit
      end if

      equal1 = equal1 - 1

    end do

    equal2 = equal

    do

      if ( equal2 >= n ) then
        exit
      end if

      if ( x(indx(equal2+1)) /= xval ) then
        exit
      end if

      equal2 = equal2 + 1

    end do
!
!  Discard certain X values.
!
    put = 0

    do get = 1, n

      if ( x(get) /= xval ) then
        put = put + 1
        x(put) = x(get)
      end if

    end do

    x(put+1:n) = 0
!
!  Adjust the INDX values.
!
    do equal = equal1, equal2
      do i = 1, n
        if ( indx(i) > indx(equal) ) then
          indx(i) = indx(i) - 1
        end if
      end do
    end do
!
!  Discard certain INDX values.
!
    indx(equal1:n+equal1-equal2-1) = indx(equal2+1:n)
    indx(n+equal1-equal2:n) = 0
!
!  Adjust N.
!
    n = put

  end if

  return
end
subroutine ivec_index_delete_dupes ( n, x, indx, n2 )
!
!*******************************************************************************
!
!! IVEC_INDEX_DELETE_DUPES deletes duplicate integer values from an indexed sorted list.
!
!
!  Discussion:
!
!    If any value occurs more than once in the input list, all duplicate
!    values are removed.
!
!    On output, the list has been sorted, and the index array has
!    been set to the identity.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input/output, integer X(N), the list.  On output, the list
!    has only unique entries, and they have been sorted explicitly.
!
!    Input/output, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, the value to be sought.
!
!    Output, integer N2, the number of unique entries in X.
!
  integer n
!
  integer i
  integer indx(*)
  integer n2
  integer x(n)
  integer y(n)
!
  i = 0
  n2 = 0

  do

    i = i + 1

    if ( i > n ) then
      exit
    end if

    if ( i > 1 ) then
      if ( x(indx(i)) == y(n2) ) then
        cycle
      end if
    end if

    n2 = n2 + 1
    y(n2) = x(indx(i))

  end do

  x(1:n2) = y(1:n2)
  call ivec_identity ( n2, indx )

  return
end
subroutine ivec_index_delete_one ( n, x, indx, xval )
!
!*******************************************************************************
!
!! IVEC_INDEX_DELETE_ONE deletes one copy of an integer value from an indexed sorted list.
!
!
!  Discussion:
!
!    If the value occurs in the list more than once, only one copy is deleted.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, integer X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer i
  integer indx(*)
  integer j
  integer less
  integer more
  integer x(*)
  integer xval
!
  if ( n < 1 ) then
    n = 0
    return
  end if

  call ivec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal /= 0 ) then
    j = indx(equal)
    x(j:n-1) = x(j+1:n)
    indx(equal:n-1) = indx(equal+1:n)
    do i = 1, n-1
      if ( indx(i) > j ) then
        indx(i) = indx(i) - 1
      end if
    end do
    n = n - 1
  end if

  return
end
subroutine ivec_index_insert ( n, x, indx, xval )
!
!*******************************************************************************
!
!! IVEC_INDEX_INSERT inserts an integer value in an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, integer X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer indx(*)
  integer less
  integer more
  integer x(*)
  integer xval
!
  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if

  call ivec_index_search ( n, x, indx, xval, less, equal, more )

  x(n+1) = xval
  indx(n+1:more+1:-1) = indx(n:more:-1)
  indx(more) = n + 1
  n = n + 1

  return
end
subroutine ivec_index_insert_unique ( n, x, indx, xval )
!
!*******************************************************************************
!
!! IVEC_INDEX_INSERT_UNIQUE inserts a unique integer value in an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N, the size of the current list.
!    If the input value XVAL does not already occur in X, then N is increased.
!
!    Input/output, integer X(N), the list.
!    If the input value XVAL does not already occur in X, then it is added
!    to X.
!
!    Input/output, integer INDX(N), the sort index of the list.
!    If the input value XVAL does not already occur in X, then INDX is updated.
!
!    Input, integer XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer indx(*)
  integer less
  integer more
  integer x(*)
  integer xval
!
  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if
!
!  Does XVAL already occur in X?
!
  call ivec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    x(n+1) = xval
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1
  end if

  return
end
subroutine ivec_index_order ( n, x, indx )
!
!*******************************************************************************
!
!! IVEC_INDEX_ORDER sorts an integer vector using an index vector.
!
!
!  Discussion:
!
!    The index vector itself is not modified.  Therefore, the pair
!    (X,INDX) no longer represents an index sorted vector.  If this
!    relationship is to be preserved, then simply set INDX(1:N)=(1:N).
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input/output, integer X(N), the list.  On output, the list
!    has been sorted.
!
!    Input, integer INDX(N), the sort index of the list.
!
  integer n
!
  integer i
  integer indx(n)
  integer x(n)
  integer y(n)
!
  y(1:n) = x(indx(1:n))
  x(1:n) = y(1:n)

  return
end
subroutine ivec_index_search ( n, x, indx, xval, less, equal, more )
!
!*******************************************************************************
!
!! IVEC_INDEX_SEARCH searches for an integer value in an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, integer X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  integer n
!
  integer equal
  integer hi
  integer indx(n)
  integer less
  integer lo
  integer mid
  integer more
  integer x(n)
  integer xhi
  integer xlo
  integer xmid
  integer xval
!
  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n
  xlo = x(indx(lo))
  xhi = x(indx(hi))

  if ( xval < xlo ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( xval == xlo ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  if ( xval > xhi ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( xval == xhi ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))

    if ( xval == xmid ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( xval < xmid ) then
      hi = mid
    else if ( xval > xmid ) then
      lo = mid
    end if

  end do

  return
end
subroutine ivec_index_sort_unique ( n, x, indx, n2 )
!
!*******************************************************************************
!
!! IVEC_INDEX_SORT_UNIQUE creates a sort index for an integer vector.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input/output, integer X(N), the list.  On output, X contains only
!    unique elements.
!
!    Output, integer INDX(N), the sort index of the list.
!
!    Output, integer N2, the number of unique elements in X.
!
  integer n
!
  integer i
  integer indx(n)
  integer n2
  integer x(n)
  integer y(n)
!
  n2 = 0

  do i = 1, n
    call ivec_index_insert_unique ( n2, y, indx, x(i) )
  end do

  x(1:n2) = y(1:n2)

  x(n2+1:n) = 0
  indx(n2+1:n) = 0

  return
end
subroutine ivec_insert ( n, a, pos, value )
!
!*******************************************************************************
!
!! IVEC_INSERT inserts a value into an array.
!
!
!  Modified:
!
!    17 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the array on input.
!
!    Input/output, integer A(N+1), the array.  On input, A is assumed
!    to contain N entries.  On output, A actually contains N+1 entries.
!
!    Input, integer POS, the position to be assigned the new entry.
!    1 <= POS <= N+1.
!
!    Input, integer VALUE, the value to be inserted at the given position.
!
  integer n
!
  integer a(n+1)
  integer i
  integer pos
  integer value
!
  if ( pos < 1 .or. pos > n+1 ) then

    write ( *, * ) ' '
    write ( *, * ) 'IVEC_INSERT - Fatal error!'
    write ( *, * ) '  Illegal insertion position = ', pos
    stop

  else

    do i = n+1, pos+1, -1
      a(i) = a(i-1)
    end do

    a(pos) = value

  end if

  return
end
subroutine ivec_max ( n, a, amax )
!
!*******************************************************************************
!
!! IVEC_MAX computes the maximum element of an integer array.
!
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer AMAX, the value of the largest entry.
!
  integer n
!
  integer a(n)
  integer amax
!
  amax = maxval ( a(1:n) )

  return
end
subroutine ivec_mean ( n, a, mean )
!
!*******************************************************************************
!
!! IVEC_MEAN returns the mean of an integer vector.
!
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector whose mean is desired.
!
!    Output, real MEAN, the mean, or average, of the vector entries.
!
  integer n
!
  integer a(n)
  real mean
!
  mean = real ( sum ( a ) ) / real ( n )

  return
end
subroutine ivec_median ( n, a, median )
!
!*******************************************************************************
!
!! IVEC_MEDIAN returns the median of an unsorted integer vector.
!
!
!  Discussion:
!
!    Hoare's algorithm is used.  The values of the vector are
!    rearranged by this routine.
!
!  Modified:
!
!    18 September 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, integer A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.
!
!    Output, integer MEDIAN, the value of the median of A.
!
  integer n
!
  integer a(n)
  integer k
  integer median
!
  k = ( n + 1 ) / 2

  call ivec_frac ( n, a, k, median )

  return
end
subroutine ivec_merge_a ( na, a, nb, b, nc, c )
!
!*******************************************************************************
!
!! IVEC_MERGE_A merges two ascending sorted integer arrays.
!
!
!  Discussion:
!
!    The elements of A and B should be sorted in ascending order.
!
!    The elements in the output array C will also be in ascending order,
!    and unique.
!
!    The output vector C may share storage with A or B.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, integer A(NA), the first sorted array.
!
!    Input, integer NB, the dimension of B.
!
!    Input, integer B(NB), the second sorted array.
!
!    Output, integer NC, the number of elements in the output array.
!    Note that C should usually be dimensioned at least NA+NB in the
!    calling routine.
!
!    Output, integer C(NC), the merged unique sorted array.
!
  integer na
  integer nb
!
  integer a(na)
  integer b(nb)
  integer c(na+nb)
  integer d(na+nb)
  integer j
  integer ja
  integer jb
  integer na2
  integer nb2
  integer nc
  integer order
!
  na2 = na
  nb2 = nb
!
  ja = 0
  jb = 0
  nc = 0

  call ivec_order_type ( na2, a, order )

  if ( order < 0 .or. order > 2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_MERGE - Fatal error!'
    write ( *, * ) '  The input array A is not ascending sorted!'
    stop
  end if

  call ivec_order_type ( nb2, b, order )

  if ( order < 0 .or. order > 2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_MERGE - Fatal error!'
    write ( *, * ) '  The input array B is not ascending sorted!'
    stop
  end if

  do
!
!  If we've used up all the entries of A, stick the rest of B on the end.
!
    if ( ja >= na2 ) then

      do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = b(jb)
        else if ( d(nc) < b(jb) ) then
          nc = nc + 1
          d(nc) = b(jb)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  If we've used up all the entries of B, stick the rest of A on the end.
!
    else if ( jb >= nb2 ) then

      do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = a(ja)
        else if ( d(nc) < a(ja) ) then
          nc = nc + 1
          d(nc) = a(ja)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  Otherwise, if the next entry of A is smaller, that's our candidate.
!
    else if ( a(ja+1) <= b(jb+1) ) then

      ja = ja + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = a(ja)
      else if ( d(nc) < a(ja) ) then
        nc = nc + 1
        d(nc) = a(ja)
      end if
!
!  ...or if the next entry of B is the smaller, consider that.
!
    else

      jb = jb + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = b(jb)
      else if ( d(nc) < b(jb) ) then
        nc = nc + 1
        d(nc) = b(jb)
      end if
    end if

  end do

  return
end
subroutine ivec_min ( n, a, amin )
!
!*******************************************************************************
!
!! IVEC_MIN computes the minimum element of an integer array.
!
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer AMAX, the value of the smallest entry.
!
  integer n
!
  integer a(n)
  integer amin
!
  amin = minval ( a(1:n) )

  return
end
function ivec_nonzero ( n, a )
!
!*******************************************************************************
!
!! IVEC_NONZERO counts the nonzero entries in an integer vector
!
!
!  Modified:
!
!    01 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the input array.
!
!    Input, integer A(N), an array.
!
!    Output, integer IVEC_NONZERO, the number of nonzero entries.
!
  integer n
!
  integer a(n)
  integer i
  integer ivec_nonzero
!
  ivec_nonzero = 0

  do i = 1, n
    if ( a(i) /= 0 ) then
      ivec_nonzero = ivec_nonzero + 1
    end if
  end do

  return
end
subroutine ivec_order_type ( n, a, order )
!
!*******************************************************************************
!
!! IVEC_ORDER_TYPE determines if an integer array is (non)strictly ascending/descending.
!
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the array.
!
!    Input, integer A(N), the array to be checked.
!
!    Output, integer ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  integer n
!
  integer a(n)
  integer i
  integer order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( i > n ) then
      order = 0
      return
    end if

    if ( a(i) > a(1) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
function ivec_pairwise_prime ( n, a )
!
!*******************************************************************************
!
!! IVEC_PAIRWISE_PRIME checks whether a vector of integers is pairwise prime.
!
!
!  Discussion:
!
!    Two positive integers I and J are pairwise prime if they have no common
!    factor greater than 1.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values to check.
!
!    Input, integer A(N), the vector of integers.
!
!    Output, logical IVEC_PAIRWISE_PRIME, is TRUE if the vector of integers
!    is pairwise prime.
!
  integer n
!
  integer a(n)
  integer i
  integer i_gcd
  logical ivec_pairwise_prime
  integer j
!
  ivec_pairwise_prime = .false.

  do i = 1, n
    do j = i+1, n
      if ( i_gcd ( a(i), a(j) ) /= 1 ) then
        return
      end if
    end do
  end do

  ivec_pairwise_prime = .true.

  return
end
subroutine ivec_part ( n, nval, a )
!
!*******************************************************************************
!
!! IVEC_PART partitions an integer NVAL into N nearly equal parts.
!
!
!  Example:
!
!    Input:
!
!      N = 5, NVAL = 17
!
!    Output:
!
!      A = ( 4, 4, 3, 3, 3 ).
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
!    Input, integer N, the number of entries in the array.
!
!    Input, integer NVAL, the integer to be partitioned.
!    NVAL may be positive, zero, or negative.
!
!    Output, integer A(N), the partition of NVAL.  The entries of
!    A add up to NVAL.  The entries of A are either all equal, or
!    differ by at most 1.  The entries of A all have the same sign
!    as NVAL, and the "largest" entries occur first.
!
  integer n
!
  integer a(n)
  integer i
  integer j
  integer nval
!
  a(1:n) = 0

  if ( nval > 0 ) then

    j = 1
    do i = 1, nval
      a(j) = a(j) + 1
      j = j + 1
      if ( j > n ) then
        j = 1
      end if
    end do

  else if ( nval < 0 ) then

    j = 1
    do i = nval, -1
      a(j) = a(j) - 1
      j = j + 1
      if ( j > n ) then
        j = 1
      end if
    end do

  end if

  return
end
subroutine ivec_part_quick_a ( n, a, l, r )
!
!*******************************************************************************
!
!! IVEC_PART_QUICK_A reorders an integer vector as part of a quick sort.
!
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1) as a key,
!    all entries of A that are less than or equal to the key will
!    precede the key which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
!
!    Output:
!
!      L = 3, R = 6
!
!      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
!            -------        -------
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of A.
!
!    Input/output, integer A(N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1).  Then
!    I <= L                 A(I) < KEY;
!         L < I < R         A(I) = KEY;
!                 R <= I    A(I) > KEY.
!
  integer n
!
  integer a(n)
  integer i
  integer key
  integer l
  integer m
  integer r
  integer temp
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_PART_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key = a(1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( a(l+1) > key ) then
      r = r - 1
      call i_swap ( a(r), a(l+1) )
    else if ( a(l+1) == key ) then
      m = m + 1
      call i_swap ( a(m), a(l+1) )
      l = l + 1
    else if ( a(l+1) < key ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(i) = a(i+m)
  end do
!
!  Out of bounds here, occasionally.
!
  l = l - m

  a(l+1:l+m) = key

  return
end
subroutine ivec_permute ( n, a, p )
!
!*******************************************************************************
!
!! IVEC_PERMUTE permutes an integer vector in place.
!
!
!  Note:
!
!    This routine permutes an array of integer "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (   1,   2,   3,   4,   5 )
!
!    Output:
!
!      A    = (   2,   4,   5,   1,   3 ).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input/output, integer A(N), the array to be permuted.
!
!    Input, integer P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  Pmust be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  integer n
!
  integer a(n)
  integer a_temp
  integer i
  integer ierror
  integer iget
  integer iput
  integer istart
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_PERMUTE - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. iget > n ) then
          write ( *, * ) ' '
          write ( *, * ) 'IVEC_PERMUTE - Fatal error!'
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine ivec_permute_random ( n, a )
!
!*******************************************************************************
!
!! IVEC_PERMUTE_RANDOM randomly permutes an integer vector.
!
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input/output, integer A(N), the array to be permuted.
!
  integer n
!
  integer a(n)
  integer p(n)
!
  call perm_random ( n, p )

  call ivec_permute ( n, a, p )

  return
end
subroutine ivec_pop ( n, x, stack1_max, stack1_num, stack1, stack2_max, &
  stack2_num, stack2 )
!
!*******************************************************************************
!
!! IVEC_POP pops an integer vector off of a stack.
!
!
!  Discussion:
!
!    If there are no more objects in the stack, N is returned as -1.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer N, the dimension of the vector.
!
!    Output, integer X(*), the value of the vector.
!
!    Input, integer STACK1_MAX, the maximum size of STACK1.
!
!    Input/output, integer STACK1_NUM, the current size of STACK1.
!
!    Input/output, integer STACK1(STACK1_MAX), the vector dimension stack.
!
!    Input, integer STACK2_MAX, the maximum size of STACK2.
!
!    Input/output, integer STACK2_NUM, the current size of STACK2.
!
!    Input/output, integer STACK2(STACK2_MAX), the vector value stack.
!
  integer n
  integer stack1_max
  integer stack2_max
!
  integer stack1(stack1_max)
  integer stack1_num
  integer stack2(stack2_max)
  integer stack2_num
  integer x(*)
!
  if ( stack1_num < 1 ) then
    n = -1
    return
  end if

  n = stack1(stack1_num)
  stack1_num = stack1_num - 1

  stack2_num = stack2_num - n
  x(1:n) = stack2(stack2_num+1:stack2_num+n)

  return
end
subroutine ivec_print ( n, a, title )
!
!*******************************************************************************
!
!! IVEC_PRINT prints an integer vector.
!
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  integer a(n)
  integer big
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, * ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i6,i11)' ) i, a(i)
    end do
  end if

  return
end
subroutine ivec_print_some ( n, a, max_print )
!
!*******************************************************************************
!
!! IVEC_PRINT_SOME prints "some" of an integer vector.
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
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
  integer n
!
  integer a(n)
  integer i
  integer max_print
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,i10)' ) i, a(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,i10)' ) i, a(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i6,2x,i10)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,i10)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(i6,2x,i10,2x,a)' ) i, a(i), '...more entries...'

  end if

  return
end
subroutine ivec_push ( n, x, stack1_max, stack1_num, stack1, stack2_max, &
  stack2_num, stack2 )
!
!*******************************************************************************
!
!! IVEC_PUSH pushes an integer vector onto a stack.
!
!
!  Discussion:
!
!    STACK1 contains a list of the dimensions of the objects stored.
!    Therefore, STACK1_MAX should be at least as big as the maximum number
!    of objects to be considered.
!
!    STACK2 contains the values of the objects.  Therefore, STACK2_MAX
!    should probably be as big as the maximum total length of the maximum
!    number of objects stored.
!
!    On first call, the user should have set STACK1_NUM and STACK2_NUM to zero.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.  N may be zero.
!
!    Input, integer X(N), the value of the vector.
!
!    Input, integer STACK1_MAX, the maximum size of STACK1.
!
!    Input/output, integer STACK1_NUM, the current size of STACK1.
!
!    Input/output, integer STACK1(STACK1_MAX), the vector dimension stack.
!
!    Input, integer STACK2_MAX, the maximum size of STACK2.
!
!    Input/output, integer STACK2_NUM, the current size of STACK2.
!
!    Input/output, integer STACK2(STACK2_MAX), the vector value stack.
!
  integer n
  integer stack1_max
  integer stack2_max
!
  integer stack1(stack1_max)
  integer stack1_num
  integer stack2(stack2_max)
  integer stack2_num
  integer x(n)
!
  if ( n < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_PUSH - Fatal error!'
    write ( *, * ) '  Input dimension N is negative.'
    stop
  end if

  if ( stack1_num + 1 > stack1_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_PUSH - Fatal error!'
    write ( *, * ) '  Exceeding size of stack #1.'
    stop
  end if

  if ( stack2_num + n > stack2_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_PUSH - Fatal error!'
    write ( *, * ) '  Exceeding size of stack #2.'
    stop
  end if

  stack1_num = stack1_num + 1
  stack1(stack1_num) = n

  stack2(stack2_num+1:stack2_num+n) = x(1:n)
  stack2_num = stack2_num + n

  return
end
subroutine ivec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! IVEC_RANDOM returns a random integer vector in a given range.
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
!    Input, integer ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, integer A(N), the vector of randomly chosen integers.
!
  integer n
!
  integer a(n)
  integer ahi
  integer alo
  integer i
!
  do i = 1, n

    call i_random ( alo, ahi, a(i) )

  end do

  return
end
subroutine ivec_red ( n, a, incx, ifact )
!
!*******************************************************************************
!
!! IVEC_RED divides out common factors in a vector.
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
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer A(*), the vector to be reduced.
!
!    If A is a simple vector, then it has dimension N.
!
!    If A is a row of a matrix, then INCX will not be 1, and
!    the actual dimension of A is at least 1+(N-1)*INCX.
!
!    On output, the entries of A have no common factor
!    greater than 1.
!
!    Input, integer INCX, the distance between successive
!    entries of A that are to be checked.
!
!    If A is a simple vector, then INCX is 1, and we simply
!    check the first N entries of A.
!
!    If A is a row of a matrix, then INCX will be the number
!    of rows declared in the matrix, in order to allow us to
!    "skip" along the row.
!
!    Output, integer IFACT, the common factor that was divided out.
!
  integer a(*)
  integer i
  integer ifact
  integer i_gcd
  integer incx
  integer indx
  integer n
!
!  Find the smallest nonzero value.
!
  ifact = 0
  indx = 1

  do i = 1, n

    if ( a(indx) /= 0 ) then

      if ( ifact == 0 ) then
        ifact = abs ( a(indx) )
      else
        ifact = min ( ifact, abs ( a(indx) ) )
      end if

    end if

    indx = indx + incx

  end do

  if ( ifact == 0 ) then
    return
  end if
!
!  Find the greatest common factor of the entire vector.
!
  indx = 1
  do i = 1, n
    ifact = i_gcd ( a(indx), ifact )
    indx = indx + incx
  end do

  if ( ifact == 1 ) then
    return
  end if
!
!  Divide out the common factor.
!
  indx = 1
  do i = 1, n
    a(indx) = a(indx) / ifact
    indx = indx + incx
  end do

  return
end
subroutine ivec_reverse ( n, a )
!
!*******************************************************************************
!
!! IVEC_REVERSE reverses the elements of an integer vector.
!
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  integer n
!
  integer a(n)
  integer i
!
  do i = 1, n/2
    call i_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine ivec_rotate ( n, m, a )
!
!*******************************************************************************
!
!! IVEC_ROTATE rotates an object in place.
!
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A = ( 1, 2, 3, 4, 5 )
!
!    Output:
!
!      A = ( 4, 5, 1, 2, 3 ).
!
!  Modified:
!
!    07 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input, integer M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, integer A(N), the array to be rotated.
!
  integer n
!
  integer a(n)
  integer i_modp
  integer iget
  integer iput
  integer istart
  integer m
  integer mcopy
  integer nset
  integer temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( istart > n ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy

      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( nset >= n ) then
      exit
    end if

  end do

  return
end
subroutine ivec_search_binary_a ( n, a, b, indx )
!
!*******************************************************************************
!
!! IVEC_SEARCH_BINARY_A searches an ascending sorted vector of integers.
!
!
!  Discussion:
!
!    Binary search is used.
!
!  Reference:
!
!    Algorithm 1.9,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vector.
!
!    Input, integer A(N), the array to be searched.  A must
!    be sorted in ascending order.
!
!    Input, integer B, the value to be searched for.
!
!    Output, integer INDX, the result of the search.
!    0, B does not occur in A.
!    I, A(I) = B.
!
  integer n
!
  integer a(n)
  integer b
  integer high
  integer indx
  integer low
  integer mid
!
  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( a(mid) < b ) then
      low = mid + 1
    else if ( a(mid) > b ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine ivec_search_binary_d ( n, a, b, indx )
!
!*******************************************************************************
!
!! IVEC_SEARCH_BINARY_D searches a descending sorted vector of integers.
!
!
!  Discussion:
!
!    Binary search is used.
!
!  Reference:
!
!    Algorithm 1.9,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vector.
!
!    Input, integer A(N), the array to be searched.  A must
!    be sorted in descending order.
!
!    Input, integer B, the value to be searched for.
!
!    Output, integer INDX, the result of the search.
!    0, B does not occur in A.
!    I, A(I) = B.
!
  integer n
!
  integer a(n)
  integer b
  integer high
  integer indx
  integer low
  integer mid
!
  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( a(mid) > b ) then
      low = mid + 1
    else if ( a(mid) < b ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine ivec_sort_bubble_a ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_BUBBLE_A ascending sorts an integer array using bubble sort.
!
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  integer n
!
  integer a(n)
  integer i
  integer j
!
  do i = 1, n-1
    do j = i+1, n
      if ( a(i) > a(j) ) then
        call i_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine ivec_sort_heap_a ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_HEAP_A ascending sorts an integer array using heap sort.
!
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  integer n
!
  integer a(n)
  integer n1
!
  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call ivec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call ivec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i_swap ( a(1), a(n1) )

  end do

  return
end
subroutine ivec_sort_heap_d ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_HEAP_D descending sorts an integer array using heap sort.
!
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  integer n
!
  integer a(n)
  integer n1
!
  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into ascending heap form.
!
  call ivec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call ivec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    call i_swap ( a(1), a(n1) )

  end do

  return
end
subroutine ivec_sort_heap_index_a ( n, a, indx )
!
!*******************************************************************************
!
!! IVEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an integer vector.
!
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call IVEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), an array to be index-sorted.
!
!    Output, integer INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  integer n
!
  integer a(n)
  integer aval
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l
!
  call ivec_identity ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( l > 1 ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine ivec_sort_heap_index_d ( n, a, indx )
!
!*******************************************************************************
!
!! IVEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an integer vector.
!
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call IVEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), an array to be index-sorted.
!
!    Output, integer INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  integer n
!
  integer a(n)
  integer aval
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l
!
  call ivec_identity ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( l > 1 ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        return
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) > a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval > a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine ivec_sort_insert_a ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_INSERT_A uses an ascending insertion sort on an integer vector.
!
!
!  Reference:
!
!    Algorithm 1.1,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer A(N).
!
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in ascending order.
!
  integer n
!
  integer a(n)
  integer i
  integer j
  integer x
!
  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( j >= 1 )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine ivec_sort_insert_d ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_INSERT_D uses a descending insertion sort on an integer vector.
!
!
!  Reference:
!
!    Algorithm 1.1,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer A(N).
!
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in descending order.
!
  integer n
!
  integer a(n)
  integer i
  integer j
  integer x
!
  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( j >= 1 )

      if ( a(j) >= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine ivec_sort_quick_a ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_QUICK_A ascending sorts an integer vector using quick sort.
!
!
!  Example:
!
!    Input:
!
!      N = 7
!
!      A = ( 6, 7, 3, 2, 9, 1, 8 )
!
!    Output:
!
!      A = ( 1, 2, 3, 6, 7, 8, 9 )
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  integer, parameter :: MAXLEVEL = 25
!
  integer n
!
  integer a(n)
  integer base
  integer l_segment
  integer level
  integer n_segment
  integer rsave(MAXLEVEL)
  integer r_segment
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IVEC_SORT_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call ivec_part_quick_a ( n_segment, a(base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( l_segment > 1 ) then

      if ( level > MAXLEVEL ) then
        write ( *, * ) ' '
        write ( *, * ) 'IVEC_SORT_QUICK_A - Fatal error!'
        write ( *, * ) '  Exceeding recursion maximum of ', MAXLEVEL
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( n_segment > 0 ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine ivec_sort_shell_a ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_SHELL_A ascending sorts an integer array using Shell's sort.
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, an array to be sorted.
!    On output, the sorted array.
!
  integer n
!
  integer a(n)
  integer asave
  integer i
  integer ifree
  integer inc
  integer ipow
  integer j
  integer k
  integer maxpow
!
  if ( n <= 1 ) then
    return
  end if
!
!  Determine the smallest MAXPOW so that
!    N <= ( 3**MAXPOW - 1 ) / 2
!
  maxpow = 1

  do while ( 2 * n + 1 > 3**maxpow )
    maxpow = maxpow + 1
  end do

  if ( maxpow > 1 ) then
    maxpow = maxpow - 1
  end if
!
!  Now sort groups of size ( 3**IPOW - 1 ) / 2.
!
  do ipow = maxpow, 1, -1

    inc = ( 3**ipow - 1 ) / 2
!
!  Sort the values with indices equal to K mod INC.
!
    do k = 1, inc
!
!  Insertion sort of the items with index
!  INC+K, 2*INC+K, 3*INC+K, ...
!
      do i = inc+k, n, inc

        asave = a(i)
        ifree = i
        j = i - inc

        do

          if ( j < 1 ) then
            exit
          end if

          if ( a(j) <= asave ) then
            exit
          end if

          ifree = j
          a(j+inc) = a(j)
          j = j - inc

        end do

        a(ifree) = asave

      end do

    end do

  end do

  return
end
subroutine ivec_split_unsort ( n, a, split, isplit )
!
!*******************************************************************************
!
!! IVEC_SPLIT_UNSORT "splits" an unsorted vector based on a splitting value.
!
!
!  Discussion:
!
!    If the vector is already sorted, it is simpler to do a binary search
!    on the data than to call this routine.
!
!    The vector is not assumed to be sorted before input, and is not
!    sorted during processing.  If sorting is not needed, then it is
!    more efficient to use this routine.
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, integer A(N), the array to split.  On output,
!    all the entries of A that are less than or equal to SPLIT
!    are in A(1:ISPLIT).
!
!    Input, integer SPLIT, the value used to split the vector.
!    It is not necessary that any value of A actually equal SPLIT.
!
!    Output, integer ISPLIT, indicates the position of the last
!    entry of the split vector that is less than or equal to SPLIT.
!
  integer n
!
  integer a(n)
  integer i
  integer i1
  integer i2
  integer i3
  integer isplit
  integer j1
  integer j2
  integer j3
  integer split
!
!  Partition the vector into A1, A2, A3, where
!    A1 = A(I1:J1) holds values <= SPLIT,
!    A2 = A(I2:J2) holds untested values,
!    A3 = A(I3:J3) holds values > SPLIT.
!
  i1 = 1
  j1 = 0

  i2 = 1
  j2 = n

  i3 = n+1
  j3 = n
!
!  Pick the next item from A2, and move it into A1 or A3.
!  Adjust indices appropriately.
!
  do i = 1, n

    if ( a(i2) <= split ) then
      i2 = i2 + 1
      j1 = j1 + 1
    else
      call i_swap ( a(i2), a(i3-1) )
      i3 = i3 - 1
      j2 = j2 - 1
    end if

  end do

  isplit = j1

  return
end
subroutine ivec_swap ( n, a1, a2 )
!
!*******************************************************************************
!
!! IVEC_SWAP swaps the entries of two integer vectors.
!
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the arrays.
!
!    Input/output, integer A1(N), A2(N), the
!    two arrays whose entries are to be swapped.
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer i
!
  do i = 1, n
    call i_swap ( a1(i), a2(i) )
  end do

  return
end
subroutine ivec_uniq ( n, a, nuniq )
!
!*******************************************************************************
!
!! IVEC_UNIQ finds the number of unique elements in a sorted integer array.
!
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in A.
!
!    Input/output, integer A(N).  On input, the sorted
!    integer array.  On output, the unique elements in A.
!
!    Output, integer NUNIQ, the number of unique elements in A.
!
  integer n
!
  integer a(n)
  integer itest
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if

  end do

  return
end
subroutine ivec_value_index ( n, a, value, max_index, n_index, value_index )
!
!*******************************************************************************
!
!! IVEC_VALUE_INDEX indexes integer vector entries equal to a given value.
!
!
!  Example:
!
!    Input:
!
!      N = 10
!      A = (  2, 3, 1, 3, 2, 4, 2, 3, 5, 3 )
!      X_VALUE = 3
!
!    Output:
!
!      N_INDEX = 4
!      VALUE_INDEX = ( 2, 4, 8, 10 ).
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input, integer A(N), the array to be indexed.
!
!    Input, integer VALUE, a value to be searched for.
!
!    Input, integer MAX_INDEX, the maximum number of indices to find.
!
!    Output, integer N_INDEX, the number of entries equal to VALUE.
!
!    Output, integer VALUE_INDEX(MAX_INDEX), the indices of entries
!    equal to VALUE.
!
  integer max_index
  integer n
!
  integer a(n)
  integer i
  integer n_index
  integer value
  integer value_index(max_index)
!
  n_index = 0

  do i = 1, n

    if ( a(i) == value ) then

      if ( n_index >= max_index ) then
        return
      end if

      n_index = n_index + 1
      value_index(n_index) = i

    end if

  end do

  return
end
subroutine ivec_variance ( n, a, variance )
!
!*******************************************************************************
!
!! IVEC_VARIANCE returns the variance of an integer vector.
!
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector whose variance is desired.
!
!    Output, real VARIANCE, the variance of the vector entries.
!
  integer n
!
  integer a(n)
  integer i
  real mean
  real variance
!
  call ivec_mean ( n, a, mean )

  variance = 0.0E+00
  do i = 1, n
    variance = variance + ( real ( a(i) ) - mean )**2
  end do

  if ( n > 1 ) then
    variance = variance / real ( n - 1 )
  else
    variance = 0.0E+00
  end if

  return
end
subroutine l_memory ( action, name, lval )
!
!*******************************************************************************
!
!! L_MEMORY manages a set of runtime logical variables.
!
!
!  Discussion:
!
!    L_MEMORY allows the user to define the name of a logical variable,
!    set it, negate it, push a new value onto a stack of up to
!    five values, pop a value off the stack, or get the current value.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, desired action.
!
!    'GET',  return value of NAME in LVAL.
!    'INIT', reset all values to zero, wipe out all names.
!    'NAME', add a variable of the given name.
!    'NOT',  negate the variable.
!    'POP',  pop the stack, retrieving the last pushed value.
!    'PRINT', print the value of NAME, and return in LVAL.
!        NAME = '*' prints all variables, and returns LVAL = .FALSE.
!    'PUSH', push the previous value down the stack, insert a new one.
!    'SET',  set variable NAME to LVAL.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    NAME may be blank for the 'INIT' command, but never any other
!    time.
!
!    Input/output, logical LVAL.
!
!    For the 'NAME', 'PUSH', and 'SET' commands, LVAL must
!    contain the initial value, push value, or set value
!    of the variable on input.
!
!    For the 'GET', 'NOT', and 'POP' commands, LVAL will contain the
!    current value, negated value, or popped value of the named
!    variable on output.
!
  integer, parameter :: maxnam = 100
  integer, parameter :: maxcol = 5
!
  character ( len = * ) action
  integer i
  integer, save, dimension ( maxnam ) :: ipoint = (/ ( 0, i = 1, maxnam ) /)
  integer j
  logical lval
  logical, save, dimension ( maxnam, maxcol ) :: lvals
  character ( len = * ) name
  character ( len = 20 ), save, dimension ( maxnam ) :: names = &
    (/ ( ' ', i = 1, maxnam ) /)
  integer, save :: numnam = 0
  logical s_eqi
!
  data ( ( lvals(i,j), i = 1, maxnam ), j = 1, maxcol ) / 500 * .false. /
!
  if ( name == ' ' .and. .not. s_eqi ( action, 'INIT' ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  The value of NAME cannot be blank!'
    stop
  end if
!
!  GET: Get the current value of a variable.
!
  if ( s_eqi ( action, 'GET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        lval = lvals(i,ipoint(i))
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to get value of unknown variable.'
    write ( *, '(''  Variable is '', a )' ) trim ( name )
    stop
!
!  INIT: Initialize everything.
!
  else if ( s_eqi ( action, 'INIT' ) ) then

    numnam = 0

    do i = 1, maxnam

      do j = 1, maxcol
        lvals(i,j) = .false.
      end do

      names(i) = ' '
      ipoint(i) = 0

    end do
!
!  NAME: Declare the name of something and set it.
!
  else if ( s_eqi ( action, 'NAME' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        write ( *, * ) ' '
        write ( *, * ) 'L_MEMORY - Warning!'
        write ( *, '(''  There is ALREADY a variable '', a ) ' ) trim ( name )
        write ( *, * ) '  The new value has been stored.'
        ipoint(i) = 1
        lvals(i,ipoint(i)) = lval
        return
      end if

    end do

    if ( numnam < maxnam ) then
      numnam = numnam + 1
      i = numnam
      names(i) = name
      ipoint(i) = 1
      lvals(i,ipoint(i)) = lval
    else
      write ( *, * ) ' '
      write ( *, * ) 'L_MEMORY - Fatal error!'
      write ( *, * ) '  We have reached the name limit of ',maxnam
      stop
    end if
!
!  NOT: Negate the value
!
  else if ( s_eqi ( action, 'NOT' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        lvals(i,ipoint(i)) =  .not. lvals(i,ipoint(i))
        lval = lvals(i,ipoint(i))
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write ( *, '(''  Variable is '', a)' ) trim ( name )
    stop
!
!  POP: "Pop" a value, decrement pointer.
!
  else if ( s_eqi ( action, 'POP' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) > 1 ) then
          ipoint(i) = ipoint(i) - 1
          lval = lvals(i,ipoint(i))
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'L_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to pop the stack to 0.'
          write ( *, '(''  Variable name is '',a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to pop an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  PRINT: "Print" the value, and return in LVAL.
!
  else if ( s_eqi ( action, 'PRINT' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) .or. name == '*' ) then
        lval = lvals(i,ipoint(i))
        write ( *, '(a,a,a,l1)' ) &
          'L_MEMORY - Value of ', trim ( names(i) ), ' is ', lval
      end if

      if ( s_eqi ( name, names(i) ) ) then
        return
      end if

    end do

    if ( name == '*' ) then
      lval = .false.
      return
    end if

    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to print an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  PUSH: "Push" a value, increment the pointer.
!
  else if ( s_eqi ( action, 'PUSH' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) < maxcol ) then
          ipoint(i) = ipoint(i) + 1
          lvals(i,ipoint(i)) = lval
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'L_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to push the stack past ', maxcol
          write ( *, '(''  Variable name is '',a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to push an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  SET: Set something.
!
  else if ( s_eqi ( action, 'SET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        lvals(i,ipoint(i)) = lval
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  Unrecognized action.
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'L_MEMORY - Fatal error!'
    write ( *, * ) '  Unrecognized action:'
    write ( *, '(a)' ) trim ( action )
    stop

  end if

  return
end
subroutine normal_01_sample ( x )
!
!*******************************************************************************
!
!! NORMAL_01_SAMPLE samples the standard Normal PDF.
!
!
!  Discussion:
!
!    The standard normal distribution has mean 0 and standard
!    deviation 1.
!
!  Method:
!
!    The Box-Muller method is used.
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
!    Output, real X, a sample of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
!
  integer, save :: iset = -1
  integer, save :: seed = 0
  real uniform_01_sample
  real v1
  real v2
  real x
  real, save :: xsave = 0.0E+00
!
  if ( iset == -1 ) then
    call random_initialize ( seed )
    iset = 0
  end if

  if ( iset == 0 ) then

!   call random_number ( harvest = v1 )

    v1 = uniform_01_sample ( seed )

    if ( v1 <= 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V1 <= 0.'
      write ( *, * ) '  V1 = ', v1
      stop
    end if

!   call random_number ( harvest = v2 )

    v2 = uniform_01_sample ( seed )

    if ( v2 <= 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V2 <= 0.'
      write ( *, * ) '  V2 = ', v2
      stop
    end if

    x = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * PI * v2 )

    xsave = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
function pause_input ( )
!
!*******************************************************************************
!
!! PAUSE_INPUT waits until an input character is entered.
!
!
!  Modified:
!
!    02 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character PAUSE_INPUT, the character that was entered.
!
  integer ios
  character pause_input
!
  write ( *, * ) 'Press RETURN to continue.'
  read ( *, '(a)', iostat = ios ) pause_input

  return
end
subroutine perm_check ( n, p, ierror )
!
!*******************************************************************************
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries.
!
!    Input, integer P(N), the array to check.
!
!    Output, integer IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  integer n
!
  integer ierror
  integer ifind
  integer iseek
  integer p(n)
!
  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      return
    end if

  end do

  return
end
subroutine perm_cycle ( n, p, isgn, ncycle, iopt )
!
!*******************************************************************************
!
!! PERM_CYCLE analyzes a permutation.
!
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    09 July 2000
!
!  Parameters:
!
!    Input/output, integer P(N).  On input, P describes a
!    permutation, in the sense that entry I is to be moved to P(I).
!    If IOPT = 0, then P will not be changed by this routine.
!    If IOPT = 1, then on output, P will be "tagged".  That is,
!    one element of every cycle in P will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of P
!    which is negative, moving to I2 = ABS(P(I1)), then to
!    P(I2), and so on, until returning to I1.
!
!    Input, integer N, the number of objects being permuted.
!
!    Output, integer ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer NCYCLE, the number of cycles in the permutation.
!
!    Input, integer IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
  integer n
!
  integer i
  integer i1
  integer i2
  integer ierror
  integer iopt
  integer is
  integer isgn
  integer ncycle
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PERM_CYCLE - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
  is = 1
  ncycle = n

  do i = 1, n

    i1 = p(i)

    do while ( i1 > i )
      ncycle = ncycle - 1
      i2 = p(i1)
      p(i1) = - i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - isign ( 1, p(i) )
    end if

    p(i) = isign ( p(i), is )

  end do

  isgn = 1 - 2 * mod ( n - ncycle, 2 )

  return
end
subroutine perm_free ( ipart, npart, ifree, nfree )
!
!*******************************************************************************
!
!! PERM_FREE reports the number of unused items in a partial permutation.
!
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
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
!    Input, integer IPART(NPART), the partial permutation, which should
!    contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer NPART, the number of entries in IPART.  NPART may be 0.
!
!    Output, integer IFREE(NFREE), the integers between 1 and NPART+NFREE
!    that were not used in IPART.
!
!    Input, integer NFREE, the number of integers that have not been
!    used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
  integer nfree
  integer npart
!
  integer i
  integer ifree(nfree)
  integer ipart(npart)
  integer j
  integer k
  integer match
  integer n
!
  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'PERM_FREE - Fatal error!'
    write ( *, * ) '  NPART < 0.'
    stop

  else if ( npart == 0 ) then

    call ivec_identity ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'PERM_FREE - Fatal error!'
    write ( *, * ) '  NFREE < 0.'
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      match = 0

      do j = 1, npart
        if ( ipart(j) == i ) then
          match = j
          exit
        end if
      end do

      if ( match == 0 ) then

        k = k + 1

        if ( k > nfree ) then
          write ( *, * ) ' '
          write ( *, * ) 'PERM_FREE - Fatal error!'
          write ( *, * ) '  The partial permutation is illegal.'
          write ( *, * ) '  It should contain, at most once, some of'
          write ( *, * ) '  the integers between 1 and ', n
          stop
        end if

        ifree(k) = i

      end if

    end do

  end if

  return
end
subroutine perm_next ( n, p, more, even )
!
!*******************************************************************************
!
!! PERM_NEXT computes all of the permutations on N objects, one at a time.
!
!
!  Discussion:
!
!    Note that if this routine is called with MORE = .TRUE., any
!    permutation in P, and EVEN = .TRUE., then the successor of the input
!    permutation will be produced, unless P is the last permutation
!    on N letters, in which case P(1) will be set to 0 on return.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    12 March 2001
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N).
!    On input, P contains the previous permutation.
!    On output, P contains the next permutation.
!
!    Input/output, logical MORE.
!    On input, MORE = FALSE means this is the first call.
!    On output, MORE = FALSE means there are no more permutations.
!
!    Output, logical EVEN, is TRUE if the output permutation is even.
!
  integer n
!
  logical even
  integer i
  integer i1
  integer ia
  integer id
  integer is
  integer j
  integer l
  integer m
  logical more
  integer p(n)
!
  if ( .not. more ) then

    call ivec_identity ( n, p )

    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n-3
      if ( p(i+1) /= p(i) + 1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n-3
        if ( p(i+1) /= p(i) + 1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      is = 0
      more = .false.

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( p(j) > ia ) then
            id = id + 1
          end if
        end do

        is = id + is

        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is+1, 2 ) * ( n + 1 )

    do j = 1, i

      if ( sign ( 1, p(j)-ia ) /= sign ( 1, p(j)-m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end
subroutine perm_print ( n, p )
!
!*******************************************************************************
!
!! PERM_PRINT prints a permutation.
!
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Modified:
!
!    06 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer which is to be partitioned.
!
!    Input, integer P(N), the permutation to be converted.
!
  integer, parameter :: inc = 20
!
  integer n
!
  integer i
  integer ihi
  integer ilo
  integer p(n)
!
  do ilo = 1, n, inc
    ihi = min ( n, ilo + inc - 1 )
    write ( *, * ) ' '
    write ( *, '(20i4)' ) ( i, i = ilo, ihi )
    write ( *, '(20i4)' ) p(ilo:ihi)
  end do

  return
end
subroutine perm_random ( n, p )
!
!*******************************************************************************
!
!! PERM_RANDOM selects a random permutation of N objects.
!
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the number of objects to be permuted.
!
!    Output, integer P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!
  integer n
!
  integer i
  integer j
  integer p(n)
!
  call ivec_identity ( n, p )

  do i = 1, n
    call i_random ( i, n, j )
    call i_swap ( p(i), p(j) )
  end do

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
subroutine points_nearest_point_bins_2d ( n, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, p, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_BINS_2D finds the nearest point to a given point in 2D.
!
!
!  Discussion:
!
!    A set of N points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point P, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if P lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing P, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to P.  We now know that
!       we don't need to search any cell whose points will all be further
!       from P than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       P than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    15 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the set.
!
!    Input, real PSET(2,N), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(N), the index of the next element of the bin
!    containing this element.
!
!    Input, real P(2), the point to be tested.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to P.
!
!    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
!
!    Output, integer COMPARES, the number of point-to-point comparisons.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin)
  integer bin_next(n)
  real c_max(ndim)
  real c_min(ndim)
  integer compares
  real d_min
  real d_min_sq
  real d_sq
  integer i
  integer i_min
  integer ic
  integer il
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real p(ndim)
  real pset(ndim,n)
  real search_radius
!
  compares = 0
!
!  Special cases.
!
  if ( n <= 0 ) then
    d_min = huge ( d_min )
    i_min = 0
    return
  end if

  if ( n == 1 ) then
    d_min = sqrt ( sum ( ( p(1:ndim) - pset(1:ndim,1) )**2 ) )
    compares = 1
    i_min = 1
    return
  end if
!
!  Initialize.
!
  d_min_sq = huge ( d_min_sq )
  i_min = 0
  search_radius = 0.0E+00
  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Determine the bin coordinates of the point P.
!
  call r2_to_bin_even2 ( nbin, bin_min, bin_max, p, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r2_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
  ic = bin(1)
  jc = bin(2)
!
!  Set
!  * the current layer,
!  * the starting bin of the current layer,
!  * the current bin
!
  layer = 0
  il = ic
  jl = jc
  i = il
  j = jl

  do
!
!  Search all legal bins in layer LAYER.
!
    do
!
!  Search BIN I, J.
!
      if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

        node = bin_start(i,j)

        do while ( node > 0 )

          d_sq = sum ( ( p(1:ndim) - pset(1:ndim,node) )**2 )
          compares = compares + 1

          if ( d_sq < d_min_sq ) then
            d_min_sq = d_sq
            i_min = node
          end if

          node = bin_next(node)

        end do

      end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
      more_bins = .true.

      do

        if ( i < ic + layer .and. j == jc - layer ) then
          i = i + 1
        else if ( i == ic + layer .and. j < jc + layer ) then
          j = j + 1
        else if ( ic - layer < i .and. j == jc + layer ) then
          i = i - 1
        else if ( i == ic - layer .and. jc - layer + 1 < j ) then
          j = j - 1
        else
          more_bins = .false.
          exit
        end if

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin ) then
          exit
        end if

      end do

      if ( .not. more_bins ) then
        exit
      end if

    end do
!
!  We've completed this layer.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = min ( &
        minval ( abs ( p(1:ndim) - c_min(1:ndim) ) ), &
        minval ( abs ( p(1:ndim) - c_max(1:ndim) ) ) )
    else
      search_radius = search_radius + layer_width
    end if
!
!  We can stop if:
!
!    * We've found at least one neighbor N;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with P at the center and the nearest N on the circumference.
!
   if ( i_min /= 0 ) then
     d_min = sqrt ( d_min_sq )
     if ( search_radius >= d_min ) then
       exit
     end if
   end if
!
!  Prepare to search the next layer.
!
    layer = layer + 1

    il = ic - layer
    jl = jc - layer

    i = il
    j = jl

  end do

  return
end
subroutine points_nearest_point_naive_2d ( n, pset, p, i_min, d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_NAIVE_2D finds the nearest point to a given point in 2D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the set.
!
!    Input, real PSET(2,N), the coordinates of the points in the set.
!
!    Input, real P(2), the point whose nearest neighbor is sought.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to P.
!
!    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
!
  integer n
  integer, parameter :: ndim = 2
!
  real d
  real d_min
  integer i
  integer i_min
  real p(ndim)
  real pset(ndim,n)
!
  d_min = huge ( d_min )
  i_min = 0

  do i = 1, n
    d = sum ( ( p(1:ndim) - pset(1:ndim,i) )**2 )
    if ( d < d_min ) then
      d_min = d
      i_min = i
    end if
  end do

  d_min = sqrt ( d_min )

  return
end
subroutine points_nearest_point_naive_3d ( n, pset, p, i_min, d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_NAIVE_3D finds the nearest point to a given point in 3D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
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
!    Input, integer N, the number of points in the set.
!
!    Input, real PSET(3,N), the coordinates of the points in the set.
!
!    Input, real P(3), the point whose nearest neighbor is sought.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to P.
!
!    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
!
  integer n
  integer, parameter :: ndim = 3
!
  real d
  real d_min
  integer i
  integer i_min
  real p(ndim)
  real pset(ndim,n)
!
  d_min = huge ( d_min )
  i_min = 0

  do i = 1, n
    d = sum ( ( p(1:ndim) - pset(1:ndim,i) )**2 )
    if ( d < d_min ) then
      d_min = d
      i_min = i
    end if
  end do

  d_min = sqrt ( d_min )

  return
end
subroutine points_nearest_points_bins2_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS2_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by calling
!    a subroutine to compute the next bin index.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    08 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  integer nbin
  integer, parameter :: ndim = 2
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer il
  integer itest
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r2_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r2_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)

    layer = 0
!
!  Search all legal bins in layer LAYER.
!
    do

      more_bins = .false.
      call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )
!
!  In layer LAYER, search each BIN I, J.
!
      do

        if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

          node = bin_start(i,j)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

            node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  I'D LIKE TO DO THIS AUTOMATICALLY...
!  Once you get past the last bin, exit because you're done the layer.
!
        do

          call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        d_min(itest) = sqrt ( d_min_sq )
        if ( search_radius >= d_min(itest) ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_bins3_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS3_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
!    user to specify the number of bins in each dimension.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN(1) by NBIN(2) regular grid of
!    cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = (/ 5, 4 /)
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN(2), the number of cells in the horizontal and
!    vertical directions.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  integer nbin(2)
  integer, parameter :: ndim = 2
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin(1),nbin(2))
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer il
  integer itest
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min(1) )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if
!
!  The efficiency of the code will suffer if the data in the vector
!
!    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) )
!
!  varies significantly.
!
  layer_width = minval ( &
    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) ) )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r2_to_bin_even3 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r2_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)

    layer = 0
!
!  Search all legal bins in layer LAYER.
!
    do

      more_bins = .false.
      call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )
!
!  In layer LAYER, search each BIN I, J.
!
      do

        if ( 1 <= i .and. i <= nbin(1) .and. 1 <= j .and. j <= nbin(2) ) then

          node = bin_start(i,j)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

            node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  I'D LIKE TO DO THIS AUTOMATICALLY...
!  Once you get past the last bin, exit because you're done the layer.
!
        do

          call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin(1) .and. &
               1 <= j .and. j <= nbin(2) ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        d_min(itest) = sqrt ( d_min_sq )
        if ( search_radius >= d_min(itest) ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_bins2_3d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS2_3D finds the nearest point to given points in 3D.
!
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J,K) and BIN_LAST(I,J,K) are given the coordinates
!    I, J, K of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J, K cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    15 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(3,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.  NBIN must be at least 3.
!
!    Input, real BIN_MIN(3), BIN_MAX(3), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(3,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  integer nbin
  integer, parameter :: ndim = 3
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin,nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer itest
  integer j
  integer jc
  integer k
  integer kc
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r3_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r3_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)
    kc = bin(3)

    layer = 0
!
!  Search all legal bins in layer LAYER.
!
    do

      more_bins = .false.

      call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, i, j, k, &
        more_bins )
!
!  In layer LAYER, search each BIN I, J, K.
!
      do

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin .and. &
             1 <= k .and. k <= nbin ) then

          node = bin_start(i,j,k)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

             node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J, K.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
        do

          call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, &
            i, j, k, more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin .and. &
               1 <= k .and. k <= nbin ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
     if ( i_min(itest) /= 0 ) then
       d_min(itest) = sqrt ( d_min_sq )
       if ( search_radius >= d_min(itest) ) then
         exit
       end if
     end if
!
!  Prepare to search the next layer.
!
      layer = layer + 1

    end do

  end do

  return
end
subroutine points_nearest_points_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    08 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  integer nbin
  integer, parameter :: ndim = 2
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer il
  integer itest
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r2_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r2_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)
!
!  Set
!  * the current layer,
!  * the starting bin of the current layer,
!  * the current bin
!
    layer = 0
    il = ic
    jl = jc
    i = il
    j = jl

    do
!
!  Search all legal bins in layer LAYER.
!
      do
!
!  Search BIN I, J.
!
        if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

          node = bin_start(i,j)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

            node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
        more_bins = .true.

        do

          if ( i < ic + layer .and. j == jc - layer ) then
            i = i + 1
          else if ( i == ic + layer .and. j < jc + layer ) then
            j = j + 1
          else if ( ic - layer < i .and. j == jc + layer ) then
            i = i - 1
          else if ( i == ic - layer .and. jc - layer + 1 < j ) then
            j = j - 1
          else
            more_bins = .false.
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed this layer.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
     if ( i_min(itest) /= 0 ) then
       d_min(itest) = sqrt ( d_min_sq )
       if ( search_radius >= d_min(itest) ) then
         exit
       end if
     end if
!
!  Prepare to search the next layer.
!
      layer = layer + 1

      il = ic - layer
      jl = jc - layer

      i = il
      j = jl

    end do

  end do

  return
end
subroutine points_nearest_points_naive_2d ( nset, pset, ntest, ptest, i_min, &
  d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_NAIVE_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
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
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST)
!    and PSET(*,I_MIN).
!
  integer nset
  integer, parameter :: ndim = 2
  integer ntest
!
  real d
  real d_min(ntest)
  integer i
  integer i_min(ntest)
  integer itest
  real pset(ndim,nset)
  real ptest(ndim,ntest)
!
  do itest = 1, ntest

    d_min(itest) = huge ( d_min )
    i_min(itest) = 0

    do i = 1, nset
      d = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,i) )**2 )
      if ( d < d_min(itest) ) then
        d_min(itest) = d
        i_min(itest) = i
      end if
    end do

    d_min(itest) = sqrt ( d_min(itest) )

  end do

  return
end
subroutine points_nearest_points_naive_3d ( nset, pset, ntest, ptest, i_min, &
  d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_NAIVE_3D finds the nearest point to given points in 3D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
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
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(3,NSET), the coordinates of the points in the set.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(3,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST)
!    and PSET(*,I_MIN).
!
  integer nset
  integer, parameter :: ndim = 3
  integer ntest
!
  real d
  real d_min(ntest)
  integer i
  integer i_min(ntest)
  integer itest
  real pset(ndim,nset)
  real ptest(ndim,ntest)
!
  do itest = 1, ntest

    d_min(itest) = huge ( d_min )
    i_min(itest) = 0

    do i = 1, nset
      d = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,i) )**2 )
      if ( d < d_min(itest) ) then
        d_min(itest) = d
        i_min(itest) = i
      end if
    end do

    d_min(itest) = sqrt ( d_min(itest) )

  end do

  return
end
function pounds_to_kilograms ( lb )
!
!*******************************************************************************
!
!! POUNDS_TO_KILOGRAMS converts a measurement in pounds to kilograms.
!
!
!  Modified:
!
!    20 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real LB, the weight in pounds.
!
!    Output, real POUNDS_TO_KILOGRAMS, the corresponding weight in kilograms.
!
  real lb
  real pounds_to_kilograms
!
  pounds_to_kilograms = 0.4535924E+00 * lb

  return
end
function prime ( n )
!
!*******************************************************************************
!
!! PRIME returns any of the first MAXPRIME prime numbers.
!
!
!  Note:
!
!    MAXPRIME is 1400, and the largest prime stored is 11657.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer N, the index of the desired prime number.
!    N = -1 returns MAXPRIME, the index of the largest prime available.
!    N = 0 is legal, returning NPRIME = 1.
!    It should generally be true that 0 <= N <= MAXPRIME.
!
!    Output, integer PRIME, the N-th prime.  If N is out of range, PRIME
!    is returned as 0.
!
  integer, parameter :: maxprime = 1400
!
  integer i
  integer n
  integer, save, dimension ( maxprime ) :: npvec
  integer prime
!
  data ( npvec(i), i = 1, 100 ) / &
          2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
         31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
         73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
        127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
        179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
        233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
        283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
        353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
        419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
        467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /

  data ( npvec(i), i = 101, 200 ) / &
        547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
        607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
        661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
        739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
        811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
        877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
        947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
       1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
       1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
       1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /

  data ( npvec(i), i = 201, 300 ) / &
       1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
       1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
       1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
       1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
       1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
       1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
       1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
       1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
       1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
       1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /

  data ( npvec(i), i = 301, 400 ) / &
       1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
       2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
       2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
       2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
       2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
       2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
       2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
       2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
       2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
       2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /

  data ( npvec(i), i = 401, 500 ) / &
       2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
       2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
       2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
       3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
       3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
       3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
       3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
       3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
       3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
       3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /

  data ( npvec(i), i = 501, 600 ) / &
       3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
       3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
       3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
       3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
       3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
       4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
       4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
       4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
       4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
       4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /

  data ( npvec(i), i = 601, 700 ) / &
       4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
       4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
       4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
       4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
       4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
       4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
       4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
       5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
       5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
       5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /

  data ( npvec(i), i = 701, 800 ) / &
       5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
       5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
       5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
       5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
       5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
       5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
       5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
       5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
       5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
       6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /

  data ( npvec(i), i = 801, 900 ) / &
       6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
       6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
       6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
       6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
       6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
       6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
       6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
       6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
       6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
       6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /

  data ( npvec(i), i = 901, 1000 ) / &
       7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
       7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
       7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
       7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
       7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
       7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
       7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
       7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
       7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
       7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /

  data ( npvec(i), i = 1001, 1100 ) / &
       7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
       8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
       8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
       8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
       8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
       8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
       8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
       8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
       8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
       8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /

  data ( npvec(i), i = 1101, 1200 ) / &
       8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
       8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
       9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
       9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
       9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
       9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
       9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
       9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
       9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
       9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /

  data ( npvec(i), i = 1201, 1300 ) / &
       9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
       9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
       9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
      10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
      10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
      10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
      10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
      10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
      10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
      10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /

  data ( npvec(i), i = 1301, 1400 ) / &
      10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
      10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
      10861,10867,10883,10889,10891,10903,10909,19037,10939,10949, &
      10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
      11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
      11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
      11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
      11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
      11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
      11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /

!
  if ( n == -1 ) then
    prime = maxprime
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= maxprime ) then
    prime = npvec(n)
  else
    prime = 0
    write ( *, * ) ' '
    write ( *, * ) 'PRIMES - Fatal error!'
    write ( *, * ) '  Illegal prime index N = ', n
    write ( *, * ) '  but N must be between 0 and MAXPRIME =',  maxprime
    stop
  end if

  return
end
subroutine primer ( iprime, n )
!
!*******************************************************************************
!
!! PRIMER computes the prime numbers up to a given limit.
!
!
!  Discussion:
!
!    PRIMER returns the results of its computations in the vector
!    IPRIME.  IPRIME(I) is -1 if the number I is not prime, and
!    1 if I is prime.
!
!    The algorithm is a simple-minded sieve of Eristothenes, with
!    no attempt at efficiency.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IPRIME(N), records the results for each integer.
!    IPRIME(I) = -1 if I is not prime, and IPRIME(I) = 1 if I is
!    prime.  By convention, IPRIME(1) will be set to -1.
!
!    Input, integer N, the dimension of IPRIME, and the maximum
!    value that will be considered.
!
  integer n
!
  integer i
  integer iprime(n)
  integer next
!
!  IPRIME(I) = 0 means we don't know if I is prime.
!
  iprime(1:n) = 0
!
!  By convention, 1 is not prime.
!
  iprime(1) = - 1
  next = 1
!
!  Examine the integers in order.
!
  do next = 2, n

    if ( iprime(next) == 0 ) then
      iprime(next) = 1
      do i = 2 * next, n, next
        iprime(i) = - 1
      end do
    end if

  end do

  return
end
subroutine r2_cheby ( n, a, alo, ahi )
!
!*******************************************************************************
!
!! R2_CHEBY sets up the Chebyshev abscissas in a real interval.
!
!
!  Discussion:
!
!    The routine sets up a vector of X values spaced between the values
!    XLO and XHI in a similar way to the spacing of the Chebyshev
!    points of the same order in the interval [-1,1].
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points to compute.
!
!    Output, real A(N), the computed X values.
!
!    Input, real ALO, AHI, the range.
!
  integer n
!
  real a(n)
  real ahi
  real alo
  real arg
  integer i
  real pi
!
  if ( n == 1 ) then

    a(1) = 0.5E+00 * ( alo + ahi )

  else if ( n > 1 ) then

    do i = 1, n

      arg = real ( 2 * i - 1 ) * pi ( ) / real ( 2 * n )

      a(i) = 0.5E+00 * ( ( 1.0E+00 + cos ( arg ) ) * alo &
        + ( 1.0E+00 - cos ( arg ) ) * ahi )

    end do

  end if

  return
end
function r2_eq ( a1, a2 )
!
!*******************************************************************************
!
!! R2_EQ == ( A1 == A2 ) for R2 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 == A2  <=>  A1(1) == A2(1) and A1(2) == A2(2).
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A1(2), A2(2).  Two R2 vectors to be compared.
!
!    Output, logical R2_EQ, is TRUE if and only if A1 == A2.
!
  integer, parameter :: ndim = 2
!
  real a1(ndim)
  real a2(ndim)
  logical r2_eq
!
  if ( all ( a1(1:ndim) == a2(1:ndim) ) ) then
    r2_eq = .true.
  else
    r2_eq = .false.
  end if

  return
end
function r2_ge ( a1, a2 )
!
!*******************************************************************************
!
!! R2_GE == ( A1 >= A2 ) for R2 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 >= A2  <=>  A1(1) > A2(1) or ( A1(1) == A2(1) and A1(2) >= A2(2) ).
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
!    Input, real A1(2), A2(2).  Two R2 vectors to be compared.
!
!    Output, logical R2_GE, is TRUE if and only if A1 >= A2.
!
  integer, parameter :: ndim = 2
!
  real a1(ndim)
  real a2(ndim)
  integer i
  logical r2_ge
!
  r2_ge = .true.

  do i = 1, ndim

    if ( a1(i) > a2(i) ) then
      r2_ge = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r2_ge = .false.
      exit
    end if

  end do

  return
end
function r2_gt ( a1, a2 )
!
!*******************************************************************************
!
!! R2_GT == ( A1 > A2 ) for R2 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>  A1(1) > A2(1) or ( A1(1) == A2(1) and A1(2) > A2(2) ).
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
!    Input, real A1(2), A2(2).  Two R2 vectors to be compared.
!
!    Output, logical R2_GT, is TRUE if and only if A1 > A2.
!
  integer, parameter :: ndim = 2
!
  real a1(ndim)
  real a2(ndim)
  integer i
  logical r2_gt
!
  r2_gt = .false.

  do i = 1, ndim

    if ( a1(i) > a2(i) ) then
      r2_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r2_gt = .false.
      exit
    end if

  end do

  return
end
function r2_le ( a1, a2 )
!
!*******************************************************************************
!
!! R2_LE == ( A1 <= A2 ) for R2 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 <= A2  <=>  A1(1) < A2(1) or ( A1(1) == A2(1) and A1(2) <= A2(2) ).
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
!    Input, real A1(2), A2(2).  Two R2 vectors to be compared.
!
!    Output, logical R2_LE, is TRUE if and only if A1 <= A2.
!
  integer, parameter :: ndim = 2
!
  real a1(ndim)
  real a2(ndim)
  integer i
  logical r2_le
!
  r2_le = .true.

  do i = 1, ndim

    if ( a1(i) < a2(i) ) then
      r2_le = .true.
      exit
    else if ( a1(i) > a2(i) ) then
      r2_le = .false.
      exit
    end if

  end do

  return
end
function r2_lt ( a1, a2 )
!
!*******************************************************************************
!
!! R2_LT == ( A1 < A2 ) for R2 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>  A1(1) < A2(1) or ( A1(1) == A2(1) and A1(2) < A2(2) ).
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
!    Input, real A1(2), A2(2).  Two R2 vectors to be compared.
!
!    Output, logical R2_LT, is TRUE if and only if A1 < A2.
!
  integer, parameter :: ndim = 2
!
  real a1(ndim)
  real a2(ndim)
  integer i
  logical r2_lt
!
  r2_lt = .false.

  do i = 1, ndim

    if ( a1(i) < a2(i) ) then
      r2_lt = .true.
      exit
    else if ( a1(i) > a2(i) ) then
      r2_lt = .false.
      exit
    end if

  end do

  return
end
function r2_ne ( a1, a2 )
!
!*******************************************************************************
!
!! R2_NE == ( A1 /= A2 ) for R2 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 /= A2  <=>  A1(1) /= A2(1) or A1(2) /= A2(2).
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
!    Input, real A1(2), A2(2).  Two R2 vectors to be compared.
!
!    Output, logical R2_NE, is TRUE if and only if A1 /= A2.
!
  integer, parameter :: ndim = 2
!
  real a1(ndim)
  real a2(ndim)
  logical r2_ne
!
  if ( any ( a1(1:ndim) /= a2(1:ndim) ) ) then
    r2_ne = .true.
  else
    r2_ne = .false.
  end if

  return
end
subroutine r2_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R2_RANDOM returns a random R2 value in a given range.
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
!    Input, real RLO(2), RHI(2), the minimum and maximum values.
!
!    Output, real R(2), the randomly chosen value.
!
  integer, parameter :: ndim = 2
!
  integer i
  real r(ndim)
  real rhi(ndim)
  real rlo(ndim)
!
  do i = 1, ndim
    call r_random ( rlo(i), rhi(i), r(i) )
  end do

  return
end
subroutine r2_swap ( x, y )
!
!*******************************************************************************
!
!! R2_SWAP swaps two R2 values.
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
!    Input/output, real X(2), Y(2).  On output, the values of X and
!    Y have been interchanged.
!
  integer, parameter :: ndim = 2
!
  real x(ndim)
  real y(ndim)
  real z(ndim)
!
  z(1:ndim) = x(1:ndim)
  x(1:ndim) = y(1:ndim)
  y(1:ndim) = z(1:ndim)

  return
end
subroutine r2_to_bin_even ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R2_TO_BIN_EVEN determines the appropriate "bin" for an R2 value.
!
!
!  Discussion:
!
!    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN-2
!    equal subintervals or bins.  Boundary bins take care of extreme values.
!
!  Example:
!
!    NBIN = 7, A(1) = 5,  A(2) = 0,
!              B(1) = 15, B(2) = 20.
!
!      C      BIN
!   ------  ------
!    8 -2    3  1
!    0  1    1  2
!    6  9    2  4
!   10 11    4  4
!   14 23    6  7
!   25 13    7  5
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins in each dimension.
!    NBIN is normally at least 3.  If NBIN is 1 or 2, then everything
!    is assigned to bin 1.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(2), a value to be placed in a bin.
!
!    Output, integer BIN(2), the index of the bin to which C is assigned.
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  real bin(ndim)
  real c(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call r_to_bin_even ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r2_to_bin_even2 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R2_TO_BIN_EVEN2 determines the appropriate "bin" for an R2 value.
!
!
!  Discussion:
!
!    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN
!    equal subintervals or bins.  Boundary bins take care of extreme values.
!
!  Example:
!
!    NBIN = 5, A(1) = 5,  A(2) = 0,
!              B(1) = 15, B(2) = 20.
!
!   20 +    +    +    +    +    +
!        15 | 25 | 35 | 45 | 55
!   16 +----+----+----+----+----+
!        14 | 24 | 34 | 44 | 54
!   12 +----+----+----+----+----+
!        13 | 23 | 33 | 43 | 53
!    8 +----+----+----+----+----+
!        12 | 22 | 32 | 42 | 52
!    4 +----+----+----+----+----+
!        11 | 21 | 31 | 41 | 51
!    0 +    +    +    +    +    +
!      5    7    9   11   13   15
!
!      C      BIN
!   ------  ------
!    8 -2    2  1
!    0  1    1  1
!    6  9    1  3
!   10 11    3  3
!   14 23    5  5
!   25 13    5  4
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
!    Input, integer NBIN, the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(2), a value to be placed in a bin.
!
!    Output, integer BIN(2), the index of the bin to which C is assigned.
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  real bin(ndim)
  real c(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call r_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r2_to_bin_even3 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R2_TO_BIN_EVEN3 determines the appropriate "bin" for an R2 value.
!
!
!  Discussion:
!
!    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
!    or bins.
!
!  Example:
!
!    NBIN = (/ 4, 5 /),
!
!      A(1) = 1,  A(2) = 0,
!      B(1) = 17, B(2) = 20.
!
!   20 +    +    +    +    +
!        15 | 25 | 35 | 45
!   16 +----+----+----+----+
!        14 | 24 | 34 | 44
!   12 +----+----+----+----+
!        13 | 23 | 33 | 43
!    8 +----+----+----+----+
!        12 | 22 | 32 | 42
!    4 +----+----+----+----+
!        11 | 21 | 31 | 41
!    0 +    +    +    +    +
!      1    5    9   13   17
!
!      C      BIN
!   ------  ------
!    8 -2    2  1
!    0  1    1  1
!    6  9    2  3
!   10 11    3  3
!   14 23    4  5
!   25 13    4  4
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(2), a value to be placed in a bin.
!
!    Output, integer BIN(2), the index of the bin to which C is assigned.
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  real bin(ndim)
  real c(ndim)
  integer i
  integer nbin(ndim)
!
  do i = 1, ndim
    call r_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r2vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BIN_EVEN bins an R2 array into evenly spaced bins.
!
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in both X and Y directions, making a total
!    of NBIN**2 2D bins.  Each set of 1D bins begins and ends with an
!    "open-ended" bin that catches extreme values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4
!      ----+-----+-----+----
!      1,3 | 2,3 | 3,3 | 4,3
!      ----+-----+-----+----
!      1,2 | 2,2 | 3,2 | 4,2
!      ----+-----+-----+----
!      1,1 | 2,1 | 3,1 | 4,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (4,4).
!
!  Modified:
!
!    09 February 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real A(2,N), the R2 data to be binned.
!
!    Input, integer NBIN, the (square root of) the number of bins.
!    NBIN must be at least 3.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), contains the index
!    of the first and last elements of A that went into each bin, or
!    -1 if there were no entries in the bin.
!
!    Output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer i1
  integer i2
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r2_to_bin_even ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)

    if ( bin_start(i1,i2) == -1 ) then
      bin_start(i1,i2) = j
    else
      k = bin_last(i1,i2)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2) = j

  end do

  return
end
subroutine r2vec_bin_even2 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BIN_EVEN2 bins an R2 array into evenly spaced bins.
!
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in both X and Y directions, making a total
!    of NBIN**2 2D bins.  Each set of 1D bins begins and ends at user
!    specified mininum and maximum values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4
!      ----+-----+-----+----
!      1,3 | 2,3 | 3,3 | 4,3
!      ----+-----+-----+----
!      1,2 | 2,2 | 3,2 | 4,2
!      ----+-----+-----+----
!      1,1 | 2,1 | 3,1 | 4,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (4,4).
!
!  Modified:
!
!    09 February 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real A(2,N), the R2 data to be binned.
!
!    Input, integer NBIN, the (square root of) the number of bins.
!    NBIN must be at least 1.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), contains the
!    index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer i1
  integer i2
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r2_to_bin_even2 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)

    if ( bin_start(i1,i2) == -1 ) then
      bin_start(i1,i2) = j
    else
      k = bin_last(i1,i2)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2) = j

  end do

  return
end
subroutine r2vec_bin_even3 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BIN_EVEN3 bins an R2 array into evenly spaced bins.
!
!
!  Discussion:
!
!    A different number of bins may be used in each dimension.
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, making a 
!    total of NBIN(1) * NBIN(2) 2D bins.  Each set of 1D bins begins and 
!    ends at user specified mininum and maximum values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4 | 5,4
!      ----+-----+-----+-----+-----
!      1,3 | 2,3 | 3,3 | 4,3 | 5,3
!      ----+-----+-----+-----+-----
!      1,2 | 2,2 | 3,2 | 4,2 | 5,2
!      ----+-----+-----+-----+-----
!      1,1 | 2,1 | 3,1 | 4,1 | 5,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (5,4).
!
!  Modified:
!
!    26 March 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real A(2,N), the R2 data to be binned.
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)), contains the
!    index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer, parameter :: ndim = 2
!
  integer n
  integer nbin(ndim)
!
  real a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2))
  real bin_max(ndim)
  real bin_min(ndim)
  integer i1
  integer i2
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin(1),1:nbin(2)) = -1
  bin_start(1:nbin(1),1:nbin(2)) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r2_to_bin_even3 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)

    if ( bin_start(i1,i2) == -1 ) then
      bin_start(i1,i2) = j
    else
      k = bin_last(i1,i2)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2) = j

  end do

  return
end
subroutine r2vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BINNED_REORDER reorders a binned R2 data vector.
!
!
!  Discussion:
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
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
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN, the (square root of the) number of bins.
!
!    Input/output, BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), contains the
!    index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Input/output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  real a2(ndim,n)
  integer bin_last(nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin)
  real bin_max
  real bin_min
  integer i
  integer i1
  integer i2
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin

    do i2 = 1, nbin

      j = bin_start(i1,i2)

      if ( j > 0 ) then
        bin_start(i1,i2) = k + 1
      end if

      do while ( j > 0 )
        k = k + 1
        bin_last(i1,i2) = k
        a2(1:ndim,k) = a(1:ndim,j)
        j = bin_next(j)
      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin
    do i2 = 1, nbin

      k = bin_last(i1,i2)

      if ( k > 0 ) then
        bin_next(k) = 0
      end if

    end do
  end do

  return
end
subroutine r2vec_binned_reorder2 ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BINNED_REORDER2 reorders a binned R2 data vector.
!
!
!  Discussion:
!
!    This routine allows there to be a different number of bins in
!    each dimension.
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
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
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN(2), the number of bins in each direction.
!
!    Input/output, BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)), 
!    contains the index of the first and last element of A that went into 
!    each bin, or -1 if there are no entries in the bin.
!
!    Input/output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer, parameter :: ndim = 2
!
  integer n
  integer nbin(ndim)
!
  real a(ndim,n)
  real a2(ndim,n)
  integer bin_last(nbin(1),nbin(2))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2))
  real bin_max
  real bin_min
  integer i
  integer i1
  integer i2
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      j = bin_start(i1,i2)

      if ( j > 0 ) then
        bin_start(i1,i2) = k + 1
      end if

      do while ( j > 0 )
        k = k + 1
        bin_last(i1,i2) = k
        a2(1:ndim,k) = a(1:ndim,j)
        j = bin_next(j)
      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin(1)
    do i2 = 1, nbin(2)

      k = bin_last(i1,i2)

      if ( k > 0 ) then
        bin_next(k) = 0
      end if

    end do
  end do

  return
end
subroutine r2vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )
!
!*******************************************************************************
!
!! R2VEC_BINNED_SORT_A sorts each bin of an R2 binned data vector.
!
!
!  Discussion:
!
!    Presumably, the data vector was first binned by R2VEC_BIN_EVEN,
!    then reordered by R2VEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R2 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
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
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN, the (square root of the) number of bins.
!
!    Input, BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), contains the index
!    of the first and last element of A that went into each bin, or -1
!    if there are no entries in this bin.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer bin_last(nbin,nbin)
  integer bin_start(nbin,nbin)
  integer i1
  integer i2
  integer j1
  integer j2
  integer n1
!
  do i1 = 1, nbin

    do i2 = 1, nbin

      j1 = bin_start(i1,i2)

      if ( j1 > 0 ) then

        j2 = bin_last(i1,i2)

        n1 = j2 + 1 - j1

        if ( n1 > 1 ) then
          call r2vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
        end if

      end if

    end do

  end do

  return
end
subroutine r2vec_binned_sort_a2 ( n, a, nbin, bin_start, bin_last )
!
!*******************************************************************************
!
!! R2VEC_BINNED_SORT_A2 sorts each bin of an R2 binned data vector.
!
!
!  Discussion:
!
!    This routine allows a difference number of bins in each dimension.
!
!    Presumably, the data vector was first binned by R2VEC_BIN_EVEN3,
!    then reordered by R2VEC_BINNED_REORDER2.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R2 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
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
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!
!    Input, BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)), the index
!    of the first and last element of A that went into each bin, or -1
!    if there are no entries in this bin.
!
  integer, parameter :: ndim = 2
!
  integer n
  integer nbin(ndim)
!
  real a(ndim,n)
  integer bin_last(nbin(1),nbin(2))
  integer bin_start(nbin(1),nbin(2))
  integer i1
  integer i2
  integer j1
  integer j2
  integer n1
!
  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      j1 = bin_start(i1,i2)

      if ( j1 > 0 ) then

        j2 = bin_last(i1,i2)

        n1 = j2 + 1 - j1

        if ( n1 > 1 ) then
          call r2vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
        end if

      end if

    end do

  end do

  return
end
subroutine r2vec_part_quick_a ( n, a, l, r )
!
!*******************************************************************************
!
!! R2VEC_PART_QUICK_A reorders an R2 vector as part of a quick sort.
!
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1:2,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
!             -----------          ----------------------------------
!             LEFT          KEY    RIGHT
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
!    Input, integer N, the number of entries of A.
!
!    Input/output, real A(2,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1:2,1).  Then
!    I <= L                 A(1:2,I) < KEY;
!         L < I < R         A(1:2,I) = KEY;
!                 R <= I    A(1:2,I) > KEY.
!
  integer n
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer i
  real key(ndim)
  integer l
  integer m
  integer r
  logical r2_eq
  logical r2_gt
  logical r2_lt
  real temp
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R2VEC_PART_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:ndim) = a(1:ndim,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r2_gt ( a(1:ndim,l+1), key(1:ndim) ) ) then
      r = r - 1
      call r2_swap ( a(1:ndim,r), a(1:ndim,l+1) )
    else if ( r2_eq ( a(1:ndim,l+1), key(1:ndim) ) ) then
      m = m + 1
      call r2_swap ( a(1:ndim,m), a(1:ndim,l+1) )
      l = l + 1
    else if ( r2_lt ( a(1:ndim,l+1), key(1:ndim) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:ndim,i) = a(1:ndim,i+m)
  end do

  l = l - m

  do i = 1, ndim
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r2vec_print ( n, a, title )
!
!*******************************************************************************
!
!! R2VEC_PRINT prints an R2 vector.
!
!
!  Modified:
!
!    03 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(2,N), A2(N), the R2 vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,(5g14.6))' ) i, a(1:ndim,i)
  end do

  return
end
subroutine r2vec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! R2VEC_RANDOM returns a random R2 vector in a given range.
!
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO(2), AHI(2), the minimum and maximum values allowed
!    for A(1,*) and A(2,*).
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(2,N), the vector of randomly chosen values.
!
  integer n
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  real ahi(ndim)
  real alo(ndim)
  integer i
  integer j
!
  do i = 1, ndim
    do j = 1, n
      call r_random ( alo(i), ahi(i), a(i,j) )
    end do
  end do

  return
end
subroutine r2vec_sort_quick_a ( n, a )
!
!*******************************************************************************
!
!! R2VEC_SORT_QUICK_A ascending sorts an R2 vector using quick sort.
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(2,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  integer, parameter :: MAXLEVEL = 25
!
  integer n
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer base
  integer l_segment
  integer level
  integer n_segment
  integer rsave(MAXLEVEL)
  integer r_segment
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R2VEC_SORT_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r2vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( l_segment > 1 ) then

      if ( level > MAXLEVEL ) then
        write ( *, * ) ' '
        write ( *, * ) 'R2VEC_SORT_QUICK_A - Fatal error!'
        write ( *, * ) '  Exceeding recursion maximum of ', MAXLEVEL
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( n_segment > 0 ) then
          exit
        end if

      end do

    end if

  end do

  return
end
function r3_eq ( a1, a2 )
!
!*******************************************************************************
!
!! R3_EQ == ( A1 == A2 ) for R3 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
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
!    Input, real A1(3), A2(3), two vectors to be compared.
!
!    Output, logical R3_EQ, is TRUE if and only if A1 == A2.
!
  integer, parameter :: ndim = 3
!
  real a1(ndim)
  real a2(ndim)
  logical r3_eq
!
  if ( all ( a1(1:ndim) == a2(1:ndim) ) ) then
    r3_eq = .true.
  else
    r3_eq = .false.
  end if

  return
end
function r3_gt ( a1, a2 )
!
!*******************************************************************************
!
!! R3_GT == ( A1 > A2 ) for R3 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
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
!    Input, real A1(3), A2(3), two vectors to be compared.
!
!    Output, logical R3_GT, is TRUE if and only if A1 > A2.
!
  integer, parameter :: ndim = 3
!
  real a1(ndim)
  real a2(ndim)
  integer i
  logical r3_gt
!
  r3_gt = .false.

  do i = 1, ndim

    if ( a1(i) > a2(i) ) then
      r3_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r3_gt = .false.
      exit
    end if

  end do

  return
end
function r3_lt ( a1, a2 )
!
!*******************************************************************************
!
!! R3_LT == ( A1 < A2 ) for R3 vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
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
!    Input, real A1(3), A2(3), two vectors to be compared.
!
!    Output, logical R3_LT, is TRUE if and only if A1 < A2.
!
  integer, parameter :: ndim = 3
!
  real a1(ndim)
  real a2(ndim)
  integer i
  logical r3_lt
!
  r3_lt = .false.

  do i = 1, ndim

    if ( a1(i) < a2(i) ) then
      r3_lt = .true.
      exit
    else if ( a1(i) > a2(i) ) then
      r3_lt = .false.
      exit
    end if

  end do

  return
end
subroutine r3_swap ( x, y )
!
!*******************************************************************************
!
!! R3_SWAP swaps two R3 values.
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
!    Input/output, real X(3), Y(3).  On output, the values of X and
!    Y have been interchanged.
!
  integer, parameter :: ndim = 3
!
  real x(ndim)
  real y(ndim)
  real z(ndim)
!
  z(1:ndim) = x(1:ndim)
  x(1:ndim) = y(1:ndim)
  y(1:ndim) = z(1:ndim)

  return
end
subroutine r3_to_bin_even2 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R3_TO_BIN_EVEN2 determines the appropriate "bin" for an R3 value.
!
!
!  Discussion:
!
!    The intervals [A(I),B(I)] are each divided into NBIN
!    equal subintervals or bins.  Boundary bins take care of extreme values.
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
!    Input, integer NBIN, the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real A(3), B(3), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(3), a value to be placed in a bin.
!
!    Output, integer BIN(3), the index of the bin to which C is assigned.
!
  integer, parameter :: ndim = 3
!
  real a(ndim)
  real b(ndim)
  real bin(ndim)
  real c(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call r_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r3vec_bin_even2 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! R3VEC_BIN_EVEN2 bins an R3 array into evenly spaced bins.
!
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in each coordinate, making a total
!    of NBIN**NDIM bins.  Each set of 1D bins begins and ends at user
!    specified mininum and maximum values.
!
!    The bins are indexed by the 1D bins that construct them,
!    and ordered lexicographically by these indices:
!
!  Modified:
!
!    09 February 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real A(3,N), the data to be binned.
!
!    Input, integer NBIN, the (cube root of) the number of bins.
!    NBIN must be at least 1.
!
!    Input, real BIN_MIN(3), BIN_MAX(3), the bin limits.
!
!    Output, BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN), contains the
!    index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin,nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer i1
  integer i2
  integer i3
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r3_to_bin_even2 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)
    i3 = bin(3)

    if ( bin_start(i1,i2,i3) == -1 ) then
      bin_start(i1,i2,i3) = j
    else
      k = bin_last(i1,i2,i3)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2,i3) = j

  end do

  return
end
subroutine r3vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! R3VEC_BINNED_REORDER reorders a binned R3 data vector.
!
!
!  Discussion:
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
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
!    Input, integer N, the number of data points.
!
!    Input/output, real A(3,N), the data to be sorted.
!
!    Input, integer NBIN, the (cube root of the) number of bins.
!
!    Input/output, BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN), contains the
!    index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Input/output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  real a2(ndim,n)
  integer bin_last(nbin,nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin,nbin)
  real bin_max
  real bin_min
  integer i
  integer i1
  integer i2
  integer i3
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin

    do i2 = 1, nbin

      do i3 = 1, nbin

        j = bin_start(i1,i2,i3)

        if ( j > 0 ) then
          bin_start(i1,i2,i3) = k + 1
        end if

        do while ( j > 0 )
          k = k + 1
          bin_last(i1,i2,i3) = k
          a2(1:ndim,k) = a(1:ndim,j)
          j = bin_next(j)
        end do

      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin
    do i2 = 1, nbin
      do i3 = 1, nbin
        k = bin_last(i1,i2,i3)

        if ( k > 0 ) then
          bin_next(k) = 0
        end if

      end do
    end do
  end do

  return
end
subroutine r3vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )
!
!*******************************************************************************
!
!! R3VEC_BINNED_SORT_A sorts each bin of an R3 binned data vector.
!
!
!  Discussion:
!
!    Presumably, the data vector was first binned by R3VEC_BIN_EVEN,
!    then reordered by R3VEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R3 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
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
!    Input, integer N, the number of data points.
!
!    Input/output, real A(3,N), the data to be sorted.
!
!    Input, integer NBIN, the (cube root of the) number of bins.
!
!    Input, BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN), contains
!    the index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in this bin.
!
  integer n
  integer nbin
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer bin_last(nbin,nbin,nbin)
  integer bin_start(nbin,nbin,nbin)
  integer i1
  integer i2
  integer i3
  integer j1
  integer j2
  integer n1
!
  do i1 = 1, nbin

    do i2 = 1, nbin

      do i3 = 1, nbin

        j1 = bin_start(i1,i2,i3)

        if ( j1 > 0 ) then

          j2 = bin_last(i1,i2,i3)

          n1 = j2 + 1 - j1

          if ( n1 > 1 ) then
            call r3vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
          end if

        end if

      end do

    end do

  end do

  return
end
subroutine r3vec_part_quick_a ( n, a, l, r )
!
!*******************************************************************************
!
!! R3VEC_PART_QUICK_A reorders an R3 vector as part of a quick sort.
!
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1:3,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
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
!    Input, integer N, the number of entries of A.
!
!    Input/output, real A(3,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1:3,1).  Then
!    I <= L                 A(1:3,I) < KEY;
!         L < I < R         A(1:3,I) = KEY;
!                 R <= I    A(1:3,I) > KEY.
!
  integer n
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer i
  real key(ndim)
  integer l
  integer m
  integer r
  logical r3_eq
  logical r3_gt
  logical r3_lt
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R3VEC_PART_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:ndim) = a(1:ndim,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r3_gt ( a(1:ndim,l+1), key(1:ndim) ) ) then
      r = r - 1
      call r3_swap ( a(1:ndim,r), a(1:ndim,l+1) )
    else if ( r3_eq ( a(1:ndim,l+1), key(1:ndim) ) ) then
      m = m + 1
      call r3_swap ( a(1:ndim,m), a(1:ndim,l+1) )
      l = l + 1
    else if ( r3_lt ( a(1:ndim,l+1), key(1:ndim) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:ndim,i) = a(1:ndim,i+m)
  end do

  l = l - m

  do i = 1, ndim
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r3vec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! R3VEC_RANDOM returns a random R3 vector in a given range.
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
!    Input, real ALO(3), AHI(3), the minimum and maximum values.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(3,N), the vector of randomly chosen values.
!
  integer n
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  real ahi(ndim)
  real alo(ndim)
  integer i
  integer j
!
  do i = 1, ndim
    do j = 1, n
      call r_random ( alo(i), ahi(i), a(i,j) )
    end do
  end do

  return
end
subroutine r3vec_sort_quick_a ( n, a )
!
!*******************************************************************************
!
!! R3VEC_SORT_QUICK_A ascending sorts an R3 vector using quick sort.
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(3,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  integer, parameter :: MAXLEVEL = 25
!
  integer n
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer base
  integer l_segment
  integer level
  integer n_segment
  integer rsave(MAXLEVEL)
  integer r_segment
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R3VEC_SORT_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r3vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( l_segment > 1 ) then

      if ( level > MAXLEVEL ) then
        write ( *, * ) ' '
        write ( *, * ) 'R3VEC_SORT_QUICK_A - Fatal error!'
        write ( *, * ) '  Exceeding recursion maximum of ', MAXLEVEL
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then
  
      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1
 
        if ( n_segment > 0 ) then
          exit
        end if

      end do

    end if

  end do

  return
end
function r_chop ( iplace, x )
!
!*******************************************************************************
!
!! R_CHOP chops a real number to a given number of binary places.
!
!
!  Example:
!
!    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
!
!    The following values would be returned for the 'chopped' value of
!    3.875:
!
!    IPLACE  Value
!
!       1      2
!       2      3     = 2 + 1
!       3      3.5   = 2 + 1 + 1/2
!       4      3.75  = 2 + 1 + 1/2 + 1/4
!       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IPLACE, the number of binary places to preserve.
!
!    IPLACE = 0 means return the integer part of X.
!    IPLACE = 1 means return the value of X, correct to 1/2.
!    IPLACE = 2 means return the value of X, correct to 1/4.
!
!    IPLACE = -1 means return the value of X, correct to 2.
!
!    Input, real X, the number to be chopped.
!
!    Output, real R_CHOP, the chopped number.
!
  real fac
  integer i_log_2
  integer iplace
  integer itemp
  real r_chop
  real x
!
  itemp = i_log_2 ( abs ( x ) )
  fac = 2.0E+00**( itemp - iplace + 1 )
  r_chop = real ( int ( x / fac ) ) * fac

  return
end
function r_cube_root ( x )
!
!*******************************************************************************
!
!! R_CUBE_ROOT returns the cube root of a real number.
!
!
!  Discussion:
!
!    R_CUBE_ROOT is designed to avoid the possible problems that can occur
!    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose cube root is desired.
!
!    Output, real R_CUBE_ROOT, the cube root of X.
!
  real r_cube_root
  real x
!
  if ( x > 0.0E+00 ) then
    r_cube_root = x**(1.0E+00/3.0E+00)
  else if ( x == 0.0E+00 ) then
    r_cube_root = 0.0E+00
  else
    r_cube_root = - ( abs ( x ) )**(1.0E+00/3.0E+00)
  end if

  return
end
function r_diff ( x, y, n )
!
!*******************************************************************************
!
!! R_DIFF computes (X-Y) to a specified accuracy.
!
!
!  Discussion:
!
!    The value returned by DIFF differs from the value X-Y,
!    because the user controls how many binary digits of accuracy
!    are to be used.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the two values whose difference is desired.
!
!    Input, integer N, determines the accuracy of the value of
!    DIFF.  If n = 10, for example, only 11 binary places will
!    be used in the arithmetic.  In general, only N+1 binary
!    places will be used.
!
!    N may be zero.  However, a negative value of N should
!    not be used, since this will cause both X and Y to look like 0.
!
!    Output, real R_DIFF, the value of X-Y.
!
  real cx
  real cy
  integer n
  real pow2
  real r_diff
  real size
  real x
  real y
!
  if ( x == y ) then
    r_diff = 0.0E+00
    return
  end if

  pow2 = 2.0E+00**n
!
!  Compute the magnitude of X and Y, and take the larger of the
!  two.  At least one of the two values is not zero!
!
  size = max ( abs ( x ), abs ( y ) )
!
!  Make normalized copies of X and Y.  One of the two values will
!  actually be equal to 1.
!
  cx = x / size
  cy = y / size
!
!  Here's where rounding comes in.  We know that the larger of the
!  the two values equals 1.  We multiply both values by 2**N,
!  where N+1 is the number of binary digits of accuracy we want
!  to use, truncate the values, and divide back by 2**N.
!
  cx = real ( int ( cx * pow2 + sign ( 0.5E+00, cx ) ) ) / pow2
  cy = real ( int ( cy * pow2 + sign ( 0.5E+00, cy ) ) ) / pow2
!
!  Take the difference now.
!
  r_diff = cx - cy
!
!  Undo the scaling.
!
  r_diff = r_diff * size

  return
end
subroutine r_digit ( x, idigit, digit )
!
!*******************************************************************************
!
!! R_DIGIT returns a particular decimal digit of a real number.
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
!    Input, real X, the number whose NDIG-th decimal digit is desired.
!    Note that if X is zero, all digits will be returned as 0.
!
!    Input, integer IDIGIT, the position of the desired decimal digit.
!    A value of 1 means the leading digit, a value of 2 the second digit
!    and so on.
!
!    Output, integer DIGIT, the value of the IDIGIT-th decimal digit of X.
!
  integer digit
  integer i
  integer idigit
  integer ival
  real x
  real xcopy
!
  if ( x == 0.0E+00 ) then
    digit = 0
    return
  end if

  if ( idigit <= 0 ) then
    digit = 0
    return
  end if
!
!  Set XCOPY = X, and then force XCOPY to lie between 1 and 10.
!
  xcopy = abs ( x )

  do while ( xcopy < 1.0E+00 )
    xcopy = xcopy * 10.0E+00
  end do

  do while ( xcopy >= 10.0E+00 )
    xcopy = xcopy / 10.0E+00
  end do

  do i = 1, idigit
    ival = int ( xcopy )
    xcopy = ( xcopy - ival ) * 10.0E+00
  end do

  digit = ival

  return
end
function r_inf ( )
!
!*******************************************************************************
!
!! R_INF returns a real number set to "Inf".
!
!
!  Modified:
!
!    24 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R_INF, the value of real infinity.
!
  real r_inf
!
  r_inf = 1.0E+00 / 0.0E+00

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
function r_log_2 ( x )
!
!*******************************************************************************
!
!! R_LOG_2 returns the logarithm base 2 of a real number.
!
!
!  Discussion:
!
!    For positive X, the formula
!
!      Log_2 ( X ) = Log ( X ) / Log ( 2.0 )
!
!    could be used.  We will use a more basic, tedious approach.
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose base 2 logarithm is desired.
!    X must be greater than 0.
!
!    Output, real R_LOG_2, the logarithm base 2 of the absolute
!    value of X.  It should be true that X = 2**R_LOG_2.
!
  integer i
  integer j
  integer k
  real r_log_2
  real temp
  real x
  real y
!
  r_log_2 = 0.0E+00

  if ( x <= 0.0E+00 ) then
    return
  end if

  y = x
!
!  Shift X so that 1 <= X < 2.
!
  i = 0

  do while ( y >= 2.0E+00 )
    i = i + 1
    y = y / 2.0E+00
  end do

  do while ( y < 1.0E+00 )
    i = i - 1
    y = y * 2.0E+00
  end do
!
!  Now it is the case that:
!
!    1 <= X < 2
!
!  and
!
!    Xold = X * 2**I
!
  j = 0
  k = 1

  temp = 0.0E+00
!
!  We need to come up with a sensible termination criterion.
!
  do while ( y /= 1.0E+00 .and. y /= 0.0E+00 .and. k < 1000000 )

    y = y * y

    j = 2 * j
    k = 2 * k

    if ( y >= 2.0E+00 ) then
      j = j + 1
      y = y / 2.0E+00
    else

    end if

    temp = real ( i ) + real ( j ) / real ( k )

  end do

  r_log_2 = temp

  return
end
subroutine r_mant ( x, is, r, l )
!
!*******************************************************************************
!
!! R_MANT computes the "mantissa" or "fraction part" of X.
!
!
!  Formula:
!
!    X = IS * R * 2**L
!
!    IS is +1 or -1,
!    R is a real between 1.0 and 2.0,
!    L is an integer.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the real number to be decomposed.
!
!    Output, integer IS, the "sign" of the number.
!    IS will be -1 if X is less than 0, and +1 if X is greater
!    than or equal to zero.
!
!    Output, real R, the mantissa of X.  R will be greater
!    than or equal to 1, and strictly less than 2.  The one
!    exception occurs if X is zero, in which case R will also
!    be zero.
!
!    Output, integer L, the integer part of the logarithm (base 2) of X.
!
  integer is
  integer l
  real r
  real x
!
!  Determine the sign.
!
  if ( x < 0.0E+00 ) then
    is = - 1
  else
    is = 1
  end if
!
!  Set R to the absolute value of X, and L to zero.
!  Then force R to lie between 1 and 2.
!
  if ( x < 0.0E+00 ) then
    r = - x
  else
    r = x
  end if

  l = 0
!
!  Time to bail out if X is zero.
!
  if ( x == 0.0E+00 ) then
    return
  end if

  do while ( r >= 2.0E+00 )
    r = r / 2.0E+00
    l = l + 1
  end do

  do while ( r < 1.0E+00 )
    r = r * 2.0E+00
    l = l - 1
  end do

  return
end
subroutine r_memory ( action, name, rval )
!
!*******************************************************************************
!
!! R_MEMORY manages a set of runtime real variables.
!
!
!  Discussion:
!
!    R_MEMORY allows the user to define the name of an real variable,
!    set it, increment it, push a new value onto a stack of up to
!    five values, pop a value off the stack, or get the current value.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, desired action.
!
!    'GET',  return value of NAME in RVAL.
!    'INC',  increment variable NAME by RVAL.
!    'INIT', reset all values to zero, wipe out all names.
!    'NAME', add a variable of the given name.
!    'POP',  pop the stack, retrieving the last pushed value.
!    'PRINT', print NAME, and return in RVAL.
!        NAME = '*' prints all variables, and returns RVAL = 0.
!    'PUSH', push the previous value down the stack, insert a new one.
!    'SET',  set variable NAME to RVAL.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    NAME may be blank for the 'INIT' command, but never any other
!    time.
!
!    Input/output, real RVAL.
!
!    For the 'INC', 'NAME', 'PUSH', and 'SET' commands, RVAL must
!    contain the increment, initial value, push value, or set value
!    of the variable on input.
!
!    For the 'GET' and 'POP' commands, RVAL will contain the current
!    value or popped value of the named variable on output.
!
  integer, parameter :: maxnam = 100
  integer, parameter :: maxcol = 5
!
  character ( len = * ) action
  integer i
  integer, save, dimension ( maxnam ) :: ipoint = (/ ( 0, i = 1, maxnam ) /)
  integer j
  character ( len = * ) name
  character ( len = 20 ), save, dimension ( maxnam ) ::  names = &
    (/ ( ' ', i = 1, maxnam ) /)
  integer, save :: numnam = 0
  real rval
  real, save, dimension ( maxnam, maxcol ) :: rvals(maxnam,maxcol)
  logical s_eqi
!
  data ( ( rvals(i,j), i = 1, maxnam), j = 1, maxcol ) / 500 * 0.0E+00 /
!
  if ( name == ' ' .and. .not. s_eqi ( action, 'INIT' ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  The value of NAME cannot be blank!'
    stop
  end if
!
!  GET: Get the current value of a variable.
!
  if ( s_eqi ( action, 'GET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        rval = rvals(i,ipoint(i))
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to get value of unknown variable.'
    write ( *, '(''  Variable name is '', a)' ) trim ( name )
    stop
!
!  INC: Increment something.
!
  else if ( s_eqi ( action, 'INC' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        rvals(i,ipoint(i)) = rvals(i,ipoint(i)) + rval
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to increment unknown variable.'
    write ( *, '(''  Variable name is '', a)' ) trim ( name )
    stop
!
!  INIT: Initialize everything.
!
  else if ( s_eqi ( action, 'INIT' ) ) then

    numnam = 0
    rvals(1:maxnam,1:maxcol) = 0.0E+00
    names(1:maxnam) = ' '
    ipoint(1:maxnam) = 0
!
!  NAME: Declare the name of something and set it.
!
  else if ( s_eqi ( action, 'NAME' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        write ( *, * ) ' '
        write ( *, * ) 'R_MEMORY - Warning!'
        write ( *, '(''  There is ALREADY a variable '', a)' ) trim ( name )
        write ( *, * ) '  The new value has been stored.'
        ipoint(i) = 1
        rvals(i,ipoint(i)) = rval
        return
      end if

    end do

    if ( numnam < maxnam ) then
      numnam = numnam + 1
      i = numnam
      names(i) = name
      ipoint(i) = 1
      rvals(i,ipoint(i)) = rval
    else
      write ( *, * ) ' '
      write ( *, * ) 'R_MEMORY - Fatal error!'
      write ( *, * ) '  We have reached the name limit of ', maxnam
      stop
    end if
!
!  POP: "Pop" a value, decrement pointer.
!
  else if ( s_eqi ( action, 'POP' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) > 1 ) then
          ipoint(i) = ipoint(i) - 1
          rval = rvals(i,ipoint(i))
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'R_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to pop the stack to 0.'
          write ( *, '(''  Variable name is '', a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to pop an unknown variable.'
    write ( *, '(''  Variable name is '', a)' ) trim ( name )
    stop
!
!  PRINT: "Print" the value, and return in RVAL.
!
  else if ( s_eqi ( action, 'PRINT' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) .or. name == '*' ) then
        rval = rvals(i,ipoint(i))
        write ( *, '(a,a,a,g14.6)' ) &
          'R_MEMORY - Value of ', trim ( names(i) ), ' is ', rval
      end if

      if ( s_eqi ( name, names(i) ) ) then
        return
      end if

    end do

    if ( name == '*' ) then
      rval = 0.0E+00
      return
    end if

    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to print an unknown variable.'
    write ( *, '(''  Variable name is '', a)' ) trim ( name )
    stop
!
!  PUSH: "Push" a value, increment the pointer.
!
  else if ( s_eqi ( action, 'PUSH' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) < maxcol ) then
          ipoint(i) = ipoint(i) + 1
          rvals(i,ipoint(i)) = rval
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'R_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to push the stack past ', maxcol
          write ( *, '(''  Variable name is '', a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to push an unknown variable.'
    write ( *, '(''  Variable name is '', a )' ) trim ( name )
    stop
!
!  SET: Set something.
!
  else if ( s_eqi ( action, 'SET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        rvals(i,ipoint(i)) = rval
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write ( *, '(''  Variable name is '', a)' ) trim ( name )
    stop
!
!  Unrecognized action.
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'R_MEMORY - Fatal error!'
    write ( *, * ) '  Unrecognized action:'
    write ( *, '(a)' ) trim ( action )
    stop

  end if

  return
end
function r_modp ( x, y )
!
!*******************************************************************************
!
!! R_MODP returns the nonnegative remainder of real division.
!
!
!  Formula:
!
!    If
!      REM = R_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R_MODP(A,360.0) is between 0 and 360, always.
!
!  Examples:
!
!        I         J     MOD  R_MODP   R_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
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
!    Input, real X, the number to be divided.
!
!    Input, real Y, the number that divides X.
!
!    Output, real R_MODP, the nonnegative remainder when X is divided by Y.
!
  real r_modp
  real x
  real y
!
  if ( y == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R_MODP - Fatal error!'
    write ( *, * ) '  R_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r_modp = mod ( x, y )

  if ( r_modp < 0.0E+00 ) then
    r_modp = r_modp + abs ( y )
  end if

  return
end
function r_nan ( )
!
!*******************************************************************************
!
!! R_NAN returns a real number set to "NaN".
!
!
!  Discussion:
!
!    R_NAN is stored as "0xFFFA5A5A" on the SGI.
!
!  Modified:
!
!    24 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R_NAN, the value of NaN.
!
  real r_nan
!
  r_nan = 0.0E+00 / 0.0E+00

  return
end
subroutine r_power ( r, p, rp, mults )
!
!*******************************************************************************
!
!! R_POWER computes the P-th power of R, for real R and integer P.
!
!
!  Discussion:
!
!    Obviously, R**P can be computed using P-1 multiplications.
!
!    However, R**P can also be computed using at most 2*LOG2(P) multiplications.
!    To do the calculation this way, let N = LOG2(P).
!    Compute A, A**2, A**4, ..., A**N by N-1 successive squarings.
!    Start the value of R**P at A, and each time that there is a 1 in
!    the binary expansion of P, multiply by the current result of the squarings.
!
!    This algorithm is not optimal.  For small exponents, and for special
!    cases, the result can be computed even more quickly.
!
!  Modified:
!
!    30 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the base.
!
!    Input, integer P, the power, which may be negative.
!
!    Output, real RP, the value of R**P.
!
!    Output, integer MULTS, the number of multiplications and divisions.
!
  integer mults
  integer p
  integer p_mag
  integer p_sign
  real r
  real r2
  real rp
!
  mults = 0
!
!  Special bases.
!
  if ( r == 1.0E+00 ) then
    rp = 1.0E+00
    return
  end if

  if ( r == - 1.0E+00 ) then

    if ( mod ( p, 2 ) == 1 ) then
      rp = - 1.0E+00
    else
      rp = 1.0E+00
    end if

    return

  end if

  if ( r == 0.0E+00 ) then

    if ( p <= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'R_POWER - Fatal error!'
      write ( *, * ) '  Base is zero, and exponent is negative.'
      stop
    end if

    rp = 0.0E+00
    return

  end if
!
!  Special powers.
!
  if ( p == -1 ) then
    rp = 1.0E+00 / r
    mults = mults + 1
    return
  else if ( p == 0 ) then
    rp = 1.0E+00
    return
  else if ( p == 1 ) then
    rp = r
    return
  end if
!
!  Some work to do.
!
  p_mag = abs ( p )
  p_sign = sign ( 1, p )

  rp = 1.0E+00
  r2 = r

  do while ( p_mag > 0 )

    if ( mod ( p_mag, 2 ) == 1 ) then
      rp = rp * r2
      mults = mults + 1
    end if

    p_mag = p_mag / 2
    r2 = r2 * r2
    mults = mults + 1

  end do

  if ( p_sign == -1 ) then
    rp = 1.0E+00 / rp
    mults = mults + 1
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
!  Discussion:
!
!    Calls to the FORTRAN 90 random number generator should go through
!    this routine, to guarantee that the random number seed has been set.
!
!  Modified:
!
!    05 February 2001
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
  real r
  real rhi
  real rlo
  integer, save :: seed = 0
  logical, save :: seeded = .false.
  real t
  real uniform_01_sample
!
!  Make sure the random number generator has been seeded.
!
  if ( .not. seeded ) then
    call random_initialize ( seed )
    seeded = .true.
  end if
!
!  Pick T, a random number in (0,1).
!
! call random_number ( harvest = t )
!
  t = uniform_01_sample ( seed )
!
!  Set R in ( RLO, RHI ).
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r_round2 ( nplace, x, xround )
!
!*******************************************************************************
!
!! R_ROUND2 rounds a number to a specified number of binary digits.
!
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 2**L
!
!    where S is plus or minus 1, L is an integer, and J is a binary
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.5 and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 2**L
!
!    where S and L are unchanged, and K is a binary mantissa which
!    agrees with J in the first NPLACE binary digits and is zero
!    thereafter.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPLACE, the number of binary digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
!    or 0.75.
!
!    Input, real X, the real number to be decomposed.
!
!    Output, real XROUND, the rounded value of X.
!
  integer iplace
  integer is
  integer l
  integer nplace
  real x
  real xmant
  real xround
  real xtemp
!
  xround = 0.0E+00
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0E+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS
!
  if ( x > 0.0E+00 ) then
    is = 1
    xtemp = x
  else
    is = - 1
    xtemp = - x
  end if
!
!  3: Force XTEMP to lie between 1 and 2, and compute the
!  logarithm L.
!
  l = 0

  do while ( xtemp >= 2.0E+00 )
    xtemp = xtemp / 2.0E+00
    l = l + 1
  end do

  do while ( xtemp < 1.0E+00 )
    xtemp = xtemp * 2.0E+00
    l = l - 1
  end do
!
!  4: Strip out the digits of the mantissa as XMANT, and decrease L.
!
  xmant = 0.0E+00
  iplace = 0

  do

    xmant = 2.0E+00 * xmant

    if ( xtemp >= 1.0E+00 ) then
      xmant = xmant + 1.0E+00
      xtemp = xtemp - 1.0E+00
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0E+00 .or. iplace >= nplace ) then
      xround = is * xmant * 2.0E+00**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * 2.0E+00

  end do

  return
end
subroutine r_roundb ( ibase, nplace, x, xround )
!
!*******************************************************************************
!
!! R_ROUNDB rounds a number to a given number of digits in a given base.
!
!
!  NOTE:
!
!    The code does not seem to do a good job of rounding when
!    the base is negative!
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * IBASE**L
!
!    where S is plus or minus 1, L is an integer, and J is a
!    mantissa base IBASE which is either exactly zero, or greater
!    than or equal to (1/IBASE) and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * IBASE**L
!
!    where S and L are unchanged, and K is a mantissa base IBASE
!    which agrees with J in the first NPLACE digits and is zero
!    thereafter.
!
!    Note that because of rounding, for most bases, most numbers
!    with a fractional quantities cannot be stored exactly in the
!    computer, and hence will have trailing "bogus" digits.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IBASE, the base of the arithmetic.
!    IBASE must not be zero.  Theoretically, IBASE may be negative.
!
!    Input, integer NPLACE, the number of digits base IBASE to
!    preserve.  NPLACE should be 0 or positive.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0,
!    1/IBASE, 2/IBASE, ..., (IBASE-1)/IBASE.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0,
!    IBASE/IBASE**2, (IBASE+1)/IBASE**2, ...,
!    IBASE**2-2/IBASE**2, IBASE**2-1/IBASE**2.
!
!    Input, real X, the real number to be decomposed.
!
!    Output, real XROUND, the rounded value of X.
!
  integer ibase
  integer iplace
  integer is
  integer js
  integer l
  integer nplace
  real x
  real xmant
  real xround
  real xtemp
!
  xround = 0.0E+00
!
!  0: Error checks.
!
  if ( ibase == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R_ROUNDB - Fatal error!'
    write ( *, * ) '  The base IBASE cannot be zero.'
    stop
  end if
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0E+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( x > 0.0E+00 ) then
    is = 1
    xtemp = x
  else
    is = - 1
    xtemp = - x
  end if
!
!  3: Force XTEMP to lie between 1 and ABS(IBASE), and compute the
!  logarithm L.
!
  l = 0

  do while ( abs ( xtemp ) >= abs ( ibase ) )

    xtemp = xtemp / real ( ibase )

    if ( xtemp < 0.0E+00 ) then
      is = - is
      xtemp = - xtemp
    end if

    l = l + 1

  end do

  do while ( abs ( xtemp ) < 1.0E+00 )

    xtemp = xtemp * ibase

    if ( xtemp < 0.0E+00 ) then
      is = - is
      xtemp = - xtemp
    end if

    l = l - 1

  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0E+00
  iplace = 0
  js = is

  do

    xmant = ibase * xmant

    if ( xmant < 0.0E+00 ) then
      js = - js
      xmant = - xmant
    end if

    if ( xtemp >= 1.0E+00 ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0E+00 .or. iplace >= nplace ) then
      xround = js * xmant * real ( ibase )**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * ibase

    if ( xtemp < 0.0E+00 ) then
      is = - is
      xtemp = - xtemp
    end if

  end do

  return
end
subroutine r_roundx ( nplace, x, xround )
!
!*******************************************************************************
!
!! R_ROUNDX rounds a real number to a specified number of decimal digits.
!
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 10**L
!
!    where S is plus or minus 1, L is an integer, and J is a decimal
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.1 and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 10**L
!
!    where S and L are unchanged, and K is a decimal mantissa which
!    agrees with J in the first NPLACE decimal digits and is zero
!    thereafter.
!
!    Note that because of rounding, most decimal fraction quantities
!    cannot be stored exactly in the computer, and hence will have
!    trailing "bogus" digits.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPLACE, the number of decimal digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
!    0.2, ..., or 0.9.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
!    0.03, ..., 0.98, 0.99.
!
!    Input, real X, the real number to be decomposed.
!
!    Output, real XROUND, the rounded value of X.
!
  integer iplace
  integer is
  integer l
  integer nplace
  real x
  real xmant
  real xround
  real xtemp
!
  xround = 0.0E+00
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0E+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( x > 0.0E+00 ) then
    is = 1
    xtemp = x
  else
    is = - 1
    xtemp = - x
  end if
!
!  3: Force XTEMP to lie between 1 and 10, and compute the
!  logarithm L.
!
  l = 0

  do while ( xtemp >= 10.0E+00 )
    xtemp = xtemp / 10.0E+00
    l = l + 1
  end do

  do while ( xtemp < 1.0E+00 )
    xtemp = xtemp * 10.0E+00
    l = l - 1
  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0E+00
  iplace = 0

  do

    xmant = 10.0E+00 * xmant

    if ( xtemp >= 1.0E+00 ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0E+00 .or. iplace >= nplace ) then
      xround = is * xmant * ( 10.0E+00**l )
      exit
    end if

    l = l - 1
    xtemp = xtemp * 10.0E+00

  end do

  return
end
function r_sign ( x )
!
!*******************************************************************************
!
!! R_SIGN evaluates the sign of a real argument.
!
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose sign is desired.
!
!    Output, real R_SIGN, a result based on the sign of X:
!
!    -1.0, if X < 0.0.
!     0.0, if X = 0.0.
!    +1.0, if X > 0.0.
!
  real r_sign
  real x
!
  if ( x < 0.0E+00 ) then
    r_sign = - 1.0E+00
  else if ( x == 0.0E+00 ) then
    r_sign = 0.0E+00
  else
    r_sign = + 1.0E+00
  end if

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP swaps two real values.
!
!
!  Modified:
!
!    01 May 2000
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
subroutine r_swap3 ( x, y, z )
!
!*******************************************************************************
!
!! R_SWAP3 swaps three real items.
!
!
!  Example:
!
!    Input:
!
!      X = 1, Y = 2, Z = 3
!
!    Output:
!
!      X = 2, Y = 3, Z = 1
!
!  Modified:
!
!    08 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y, Z, three values to be swapped.
!
  real w
  real x
  real y
  real z
!
  w = x
  x = y
  y = z
  z = w

  return
end
subroutine r_to_bin_even ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R_TO_BIN_EVEN determines the appropriate "bin" for C in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
!    An initial bin takes everything less than A, and a final bin takes
!    everything greater than B.
!
!  Example:
!
!    NBIN = 7, A = 5, B = 15
!
!    C   BIN
!
!    1    1
!    3    1
!    4.9  1
!    5    2
!    6    2
!    7    3
!    8    3
!    9.5  4
!   13    6
!   14    6
!   15    6
!   15.1  7
!   99    7
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.  NBIN is normally at least 3.
!    If NBIN is 1 or 2, then everything is assigned to bin 1.
!
!    Input, real A, B, the lower and upper limits of the bin interval.
!    While A is expected to be less than B, the code should return useful
!    results if A is actually greater than B.
!
!    Input, real C, a value to be placed in a bin.
!
!    Output, integer BIN, the index of the bin to which C is assigned.
!
  real a
  real a2
  real b
  real b2
  integer bin
  real c
  integer nbin
  logical switch
!
!  Take care of special cases.
!
  if ( nbin < 1 ) then
    bin = 0
    return
  else if ( nbin == 1 .or. nbin == 2 ) then
    bin = 1
    return
  end if

  if ( b == a ) then
    bin = 0
    return
  end if
!
!  If the limits are descending, then we switch them now, and
!  unswitch the results at the end.
!
  if ( a < b ) then
    switch = .false.
    a2 = a
    b2 = b
  else
    switch = .true.
    a2 = b
    b2 = a
  end if
!
!  Compute the bin.
!
  if ( c < a2 ) then
    bin = 1
  else if ( c == a2 ) then
    bin = 2
  else if ( c == b2 ) then
    bin = nbin - 1
  else if ( c > b2 ) then
    bin = nbin
  else
    bin = 2 + int ( real ( nbin - 2 ) * ( c - a2 ) / ( b2 - a2 ) )
    bin = max ( bin, 2 )
    bin = min ( bin, nbin-1 )
  end if
!
!  Reverse the switching.
!
  if ( switch ) then
    bin = nbin + 1 - bin
  end if

  return
end
subroutine r_to_bin_even2 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R_TO_BIN_EVEN2 determines the appropriate "bin" for C in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A = 5, B = 15
!
!    <-1-+-2-+-3-+-4-+-5->
!    5   7   9  11  13  15
!
!
!    C   BIN
!
!    1    1
!    3    1
!    4.9  1
!    5    1
!    6    1
!    7.1  2
!    8    2
!    9.5  3
!   12    4
!   14    5
!   15    5
!   15.1  5
!   99    5
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, real A, B, the lower and upper limits of the bin interval.
!    While A is expected to be less than B, the code should return useful
!    results if A is actually greater than B.
!
!    Input, real C, a value to be placed in a bin.
!
!    Output, integer BIN, the index of the bin to which C is assigned.
!
  real a
  real a2
  real b
  real b2
  integer bin
  real c
  integer nbin
  logical switch
!
!  Take care of special cases.
!
  if ( nbin < 1 ) then
    bin = 0
    return
  end if

  if ( nbin == 1 ) then
    bin = 1
    return
  end if

  if ( b == a ) then
    bin = 1
    return
  end if
!
!  If the limits are descending, then we switch them now, and
!  unswitch the results at the end.
!
  if ( a < b ) then
    switch = .false.
    a2 = a
    b2 = b
  else
    switch = .true.
    a2 = b
    b2 = a
  end if
!
!  Compute the bin.
!
  if ( c <= a2 ) then
    bin = 1
  else if ( c >= b2 ) then
    bin = nbin
  else
    bin = 1 + int ( real ( nbin ) * ( c - a2 ) / ( b2 - a2 ) )
    bin = max ( bin, 1 )
    bin = min ( bin, nbin )
  end if
!
!  Reverse the switching.
!
  if ( switch ) then
    bin = nbin + 1 - bin
  end if

  return
end
subroutine r_to_dhms ( r, d, h, m, s )
!
!*******************************************************************************
!
!! R_TO_DHMS converts a real number of days into days, hours, minutes, seconds.
!
!
!  Modified:
!
!    08 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, a real number representing a time period measured in days.
!
!    Output, integer D, H, M, S, the equivalent number of days, hours,
!    minutes and seconds.
!
  integer d
  integer h
  integer m
  real r
  real r_copy
  integer s
!
  r_copy = abs ( r )

  d = int ( r_copy )

  r_copy = r_copy - d
  r_copy = 24.0E+00 * r_copy
  h = int ( r_copy )

  r_copy = r_copy - h
  r_copy = 60.0E+00 * r_copy
  m = int ( r_copy )

  r_copy = r_copy - m
  r_copy = 60.0E+00 * r_copy
  s = int ( r_copy )

  if ( r < 0.0E+00 ) then
    d = - d
    h = - h
    m = - m
    s = - s
  end if

  return
end
subroutine r_to_r_discrete ( r, rmin, rmax, nr, rd )
!
!*******************************************************************************
!
!! R_TO_R_DISCRETE maps real X to XD in [XDMIN, XDMAX] with NDX possible values.
!
!
!  Formula:
!
!    if ( X < XDMIN ) then
!      XD = XDMIN
!    else if ( X > XDMAX ) then
!      XD = XDMAX
!    else
!      T = nint ( ( NDX - 1 ) * ( X - XDMIN ) / ( XDMAX - XDMIN ) )
!      XD = XDMIN + T * ( XDMAX - XDMIN ) / real ( NDX - 1 )
!
!    In the special case where NDX = 1, when
!
!      XD = 0.5 * ( XDMAX + XDMIN )
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the real number to be converted.
!
!    Input, real XDMAX, XDMIN, the maximum and minimum values for XD.
!
!    Input, integer NR, the number of allowed values for XD.
!    NR should be at least 1.
!
!    Output, real XD, the corresponding discrete value.
!
  integer i
  integer nr
  real t
  real r
  real rd
  real rmax
  real rmin
!
!  Check for errors.
!
  if ( nr < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R_TO_R_DISCRETE - Fatal error!'
    write ( *, * ) '  NR = ', nr
    write ( *, * ) '  but NR must be at least 1.'
    stop
  end if

  if ( nr == 1 ) then
    rd = 0.5E+00 * ( rmin + rmax )
    return
  end if

  if ( rmax == rmin ) then
    rd = rmax
    return
  end if

  if ( rmax < rmin ) then
    write ( *, * ) ' '
    write ( *, * ) 'R_TO_R_DISCRETE - Fatal error!'
    write ( *, * ) '  RMIN = ', rmin
    write ( *, * ) '  RMAX = ', rmax
    write ( *, * ) '  but RMIN should be less than RMAX.'
    stop
  end if
!
!  Set extreme values.
!
  if ( r <= rmin ) then
    rd = rmin
    return
  else if ( r >= rmax ) then
    rd = rmax
    return
  end if
!
!  For RMIN < R < RMAX, T is the corresponding value in [0,NR-1].
!
  t = real ( nr - 1 ) * ( r - rmin ) / ( rmax - rmin )
!
!  I should be the corresponding integer in 0, ..., NR-1.
!
  i = nint ( t )

  rd = rmin + ( real ( i ) / real ( nr - 1 ) ) * ( rmax - rmin )

  return
end
subroutine r_unswap3 ( x, y, z )
!
!*******************************************************************************
!
!! R_UNSWAP3 unswaps three real items.
!
!
!  Example:
!
!    Input:
!
!      X = 2, Y = 3, Z = 1
!
!    Output:
!
!      X = 1, Y = 2, Z = 3
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y, Z, three values to be swapped.
!
  real w
  real x
  real y
  real z
!
  w = z
  z = y
  y = x
  x = w

  return
end
function r_zeta ( p )
!
!*******************************************************************************
!
!! R_ZETA estimates the Riemann Zeta function.
!
!
!  Definition:
!
!    For 1 < P, the Riemann Zeta function is defined as:
!
!      ZETA ( P ) = Sum ( 1 <= N < Infinity ) 1 / N**P
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
!    08 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real P, the power to which the integers are raised.
!    P must be greater than 1.  For integral P up to 20, a
!    precomputed value of R_ZETA is returned; otherwise the infinite
!    sum is approximated.
!
!    Output, real R_ZETA, an approximation to the Riemann Zeta function.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
!
  integer n
  real p
  real r_zeta
  double precision sum
  double precision sum_old
!
  if ( p <= 1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R_ZETA - Fatal error!'
    write ( *, * ) '  Exponent P <= 1.0.'
    stop
  else if ( p == 2.0E+00 ) then
    r_zeta = PI**2 / 6.0E+00
  else if ( p == 3.0E+00 ) then
    r_zeta = 1.2020569032E+00
  else if ( p == 4.0E+00 ) then
    r_zeta = PI**4 / 90.0E+00
  else if ( p == 5.0E+00 ) then
    r_zeta = 1.0369277551E+00
  else if ( p == 6.0E+00 ) then
    r_zeta = PI**6 / 945.0E+00
  else if ( p == 7.0E+00 ) then
    r_zeta = 1.0083492774E+00
  else if ( p == 8.0E+00 ) then
    r_zeta = PI**8 / 9450.0E+00
  else if ( p == 9.0E+00 ) then
    r_zeta = 1.0020083928E+00
  else if ( p == 10.0E+00 ) then
    r_zeta = PI**10 / 93555.0E+00
  else if ( p == 11.0E+00 ) then
    r_zeta = 1.0004941886E+00
  else if ( p == 12.0E+00 ) then
    r_zeta = 1.0002460866E+00
  else if ( p == 13.0E+00 ) then
    r_zeta = 1.0001227133E+00
  else if ( p == 14.0E+00 ) then
    r_zeta = 1.0000612482E+00
  else if ( p == 15.0E+00 ) then
    r_zeta = 1.0000305882E+00
  else if ( p == 16.0E+00 ) then
    r_zeta = 1.0000152823E+00
  else if ( p == 17.0E+00 ) then
    r_zeta = 1.0000076372E+00
  else if ( p == 18.0E+00 ) then
    r_zeta = 1.0000038173E+00
  else if ( p == 19.0E+00 ) then
    r_zeta = 1.0000019082E+00
  else if ( p == 20.0E+00 ) then
    r_zeta = 1.0000009540E+00
  else

    sum = 0.0D+00
    n = 0

    do

      n = n + 1
      sum_old = sum
      sum = sum + 1.0D+00 / ( dble ( n ) )**p

      if ( real ( sum ) <= real ( sum_old ) .or. n >= 1000 ) then
        exit
      end if

    end do

    r_zeta = real ( sum )

  end if

  return
end
subroutine random_initialize ( seed )
!
!*******************************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!    And this is the FORTRAN 90 people's idea of convenience?
!
!    And I still get poorly randomized values, somehow, having to do
!    with a bad seed, or something.  I am about ready to go back to
!    using my own damn routine!
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer SEED, a seed value.
!
  integer date_time(8)
  integer i
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
  real t
  integer value
!
!  Initialize the random number seed.
!
  call random_seed
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )
!
!  Get the current date and time.
!
  call date_and_time ( values = date_time )
!
!  Construct a slightly random value.
!
  seed = 0
  do i = 1, 8
    seed = ieor ( seed, date_time(i) )
  end do
!
!  Make slightly random assignments to SEED_VECTOR.
!
  do i = 1, seed_size
    seed_vector(i) = ieor ( seed, i )
  end do
!
!  Set the random number seed value.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Because EVEN THIS DOESN'T SEEM TO PROPERLY MIX UP THE RANDOM
!  NUMBERS, call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do
!
!  I STILL GET LOUSY RESULTS.  THE HELL WITH IT!
!
  return
end
subroutine rat_factor ( m, n, maxfactor, nfactor, factor, power, mleft, nleft )
!
!*******************************************************************************
!
!! RAT_FACTOR factors a rational value into a product of prime factors.
!
!
!  Formula:
!
!    ( M / N ) = ( MLEFT / NLEFT ) * Product ( I = 1 to NFACTOR )
!      FACTOR(I)**POWER(I).
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the top and bottom of a rational value.
!    The ratio of M and N must be positive.
!
!    Input, integer MAXFACTOR, the maximum number of factors for
!    which storage has been allocated.
!
!    Output, integer NPRIME, the number of prime factors of M/N.
!
!    Output, integer FACTOR(MAXFACTOR), the prime factors of M/N.
!
!    Output, integer POWER(MAXFACTOR).  POWER(I) is the power of
!    the FACTOR(I) in the representation of M/N.
!
!    Output, integer MLEFT, NLEFT, the top and bottom of the factor of
!    M / N that remains.  If ABS ( MLEFT / NLEFT ) is not 1, then
!    the rational value was not completely factored.
!
  integer maxfactor
!
  integer factor(maxfactor)
  integer i
  integer m
  integer maxprime
  integer mleft
  integer n
  integer nleft
  integer nfactor
  integer p
  integer power(maxfactor)
  integer prime
!
  nfactor = 0

  mleft = m
  nleft = n
!
!  NLEFT should be nonnegative.
!
  if ( nleft < 0 ) then
    mleft = - mleft
    nleft = - nleft
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  if ( m == n ) then
    nfactor = 1
    factor(1) = 1
    power(1) = 1
    return
  end if
!
!  Find out how many primes we stored.
!
  maxprime = prime ( -1 )

  do i = 1, maxprime

    p = prime ( i )

    if ( mod ( nleft, p ) == 0 .or. &
         mod ( abs ( mleft ), p ) == 0 ) then

      if ( nfactor < maxfactor ) then

        nfactor = nfactor + 1
        factor(nfactor) = p
        power(nfactor) = 0
!
!  Divide MLEFT by PRIME(I) as often as you can.
!
        if ( mod ( abs ( mleft ), p ) == 0  ) then

          do

            power(nfactor) = power(nfactor) + 1
            mleft = mleft / p

            if ( mod ( abs ( mleft ), p ) /= 0 ) then
              exit
            end if

          end do

        end if
!
!  Divide NLEFT by PRIME(I) as often as you can.
!
        if ( mod ( nleft, p ) == 0  ) then

          do

            power(nfactor) = power(nfactor) - 1
            nleft = nleft / p

            if ( mod ( nleft, p ) /= 0 ) then
              exit
            end if

          end do

        end if

        if ( power(nfactor) == 0 ) then
          nfactor = nfactor - 1
        end if

      end if

    end if

  end do

  return
end
subroutine rcol_compare ( lda, m, n, a, i, j, isgn )
!
!*******************************************************************************
!
!! RCOL_COMPARE compares columns I and J of a real array of column data.
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
!      ISGN = -1
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
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), the M by N array.
!
!    Input, integer I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column I > column J.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer isgn
  integer j
  integer k
  integer m
!
!  Check.
!
  if ( i < 1 .or. i > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'RCOL_COMPARE - Fatal error!'
    write ( *, * ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. j > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'RCOL_COMPARE - Fatal error!'
    write ( *, * ) '  Column index J is out of bounds.'
    stop
  end if
!
  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = - 1
      return
    else if ( a(k,i) > a(k,j) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine rcol_find ( lda, m, n, a, x, icol )
!
!*******************************************************************************
!
!! RCOL_FIND seeks a table column equal to a real vector.
!
!
!  Example:
!
!    Input:
!
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      x = ( 3.,
!            7.,
!           11. )
!
!    Output:
!
!      ICOL = 3
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), a table of numbers, regarded as
!    N columns of vectors of length M.
!
!    Input, real X(M), a vector to be matched with a column of A.
!
!    Output, integer ICOL, the index of the first column of A
!    which exactly matches every entry of X, or 0 if no match
!    could be found.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer icol
  integer j
  real x(m)
!
  icol = 0

  do j = 1, n

    icol = j

    do i = 1, m
      if ( x(i) /= a(i,j) ) then
        icol = 0
        exit
      end if
    end do

    if ( icol /= 0 ) then
      return
    end if

  end do

  return
end
subroutine rcol_insert ( lda, maxcol, m, n, a, x, icol )
!
!*******************************************************************************
!
!! RCOL_INSERT inserts a column into a sorted table, if it is new.
!
!
!  Example:
!
!    Input:
!
!      MAXCOL = 10,
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      X = ( 3., 4., 18. )
!
!    Output:
!
!      N = 5,
!
!      A = (
!        1.  2.  3.  3.  4.
!        5.  6.  4.  7.  8.
!        9. 10. 18. 11. 12. )
!
!      ICOL = 3
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer MAXCOL, the maximum number of columns in A.
!
!    Input, integer M, the number of rows.
!
!    Input/output, integer N, the number of columns.
!    If the new column is inserted into the table, then the output
!    value of N will be increased by 1.
!
!    Input/output, real A(LDA,MAXCOL), a table of numbers, regarded
!    as an array of columns.  The columns must have been sorted
!    lexicographically.
!
!    Input, real X(M), a vector of data which will be inserted
!    into the table if it does not already occur.
!
!    Output, integer ICOL.
!    I, X was inserted into column I.
!    -I, column I was already equal to X.
!    0, N = MAXCOL.
!
  integer lda
  integer m
  integer maxcol
!
  real a(lda,maxcol)
  integer high
  integer i
  integer icol
  integer isgn
  integer j
  integer low
  integer mid
  integer n
  real x(m)
!
!  Refuse to work if N >= MAXCOL.
!
  if ( n >= maxcol ) then
    icol = 0
    return
  end if
!
!  Stick X temporarily in column N+1, just so it's easy to use RCOL_COMPARE.
!
  a(1:m,n+1) = x(1:m)
!
!  Do a binary search.
!
  low = 1
  high = n

  do

    if ( low > high ) then
      icol = low
      exit
    end if

    mid = ( low + high ) / 2

    call rcol_compare ( lda, m, n+1, a, mid, n+1, isgn )

    if ( isgn == 0 ) then
      icol = - mid
      return
    else if ( isgn == -1 ) then
      low = mid + 1
    else if ( isgn == +1 ) then
      high = mid - 1
    end if

  end do
!
!  Shift part of the table up to make room.
!
  do j = n, icol, -1
    a(1:m,j+1) = a(1:m,j)
  end do
!
!  Insert the new column.
!
  a(1:m,icol) = x(1:m)

  n = n + 1

  return
end
subroutine rcol_max ( lda, m, n, a, imax, amax )
!
!*******************************************************************************
!
!! RCOL_MAX returns the maximums of columns of a real array.
!
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, integer IMAX(N); IMAX(I) is the row of A in which
!    the maximum for column I occurs.
!
!    Output, real AMAX(N), the maximums of the columns.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real amax(n)
  integer i
  integer imax(n)
  integer j
!
  do j = 1, n

    imax(j) = 1
    amax(j) = a(1,j)
    do i = 2, m
      if ( a(i,j) > amax(j) ) then
        imax(j) = i
        amax(j) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine rcol_mean ( lda, m, n, a, mean )
!
!*******************************************************************************
!
!! RCOL_MEAN returns the means of columns of a real array.
!
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      1.5  4.0  5.0
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, real MEAN(N), the means, or averages, of the columns.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer j
  real mean(n)
!
  do j = 1, n
    mean(j) = sum ( a(1:m,j) )
  end do

  mean(1:n) = mean(1:n) / real ( m )

  return
end
subroutine rcol_min ( lda, m, n, a, imin, amin )
!
!*******************************************************************************
!
!! RCOL_MIN returns the minimums of columns of a real array.
!
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, integer IMIN(N); IMIN(I) is the row of A in which
!    the minimum for column I occurs.
!
!    Output, real AMIN(N), the minimums of the columns.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real amin(n)
  integer i
  integer imin(n)
  integer j
!
  do j = 1, n

    imin(j) = 1
    amin(j) = a(1,j)
    do i = 2, m
      if ( a(i,j) < amin(j) ) then
        imin(j) = i
        amin(j) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine rcol_sort_a ( lda, m, n, a )
!
!*******************************************************************************
!
!! RCOL_SORT_A ascending sorts a real array of columns.
!
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real A(LDA,N).
!    On input, the real array of N columns of real M-vectors.
!    On output, the columns of A have been sorted in lexicographic order.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer indx
  integer isgn
  integer j
  integer m
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call rcol_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call rcol_compare ( lda, m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine rcol_sortr_a ( lda, m, n, a, key )
!
!*******************************************************************************
!
!! RCOL_SORTR_A ascending sorts one column of a 2D array, adjusting all entries.
!
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array, which
!    must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real A(LDA,N).
!
!    On input, an unsorted M by N array.
!
!    On output, rows of the array have been shifted in such
!    a way that column KEY of the array is in nondecreasing order.
!
!    Input, integer KEY, the column in which the "key" value
!    is stored.  On output, column KEY of the array will be
!    in nondecreasing order.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer indx
  integer isgn
  integer j
  integer key
  integer m
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call rrow_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( a(i,key) < a(j,key) ) then
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
subroutine rcol_sum ( lda, m, n, a, colsum )
!
!*******************************************************************************
!
!! RCOL_SUM sums the entries of columns of a real array.
!
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, real COLSUM(N), the sums of the columns.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real colsum(n)
  integer i
  integer j
!
  do j = 1, n
    colsum(j) = sum ( a(1:m,j) )
  end do

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
subroutine rcol_to_rvec ( lda, m, n, a, x )
!
!*******************************************************************************
!
!! RCOL_TO_RVEC converts a matrix of columns into a vector.
!
!
!  Example:
!
!    M = 3, N = 4
!
!    A =
!      11 12 13 14
!      21 22 23 24
!      31 32 33 34
!
!    X = ( 11, 21, 31, 12, 22, 32, 13, 23, 33, 14, 24, 34 )
!
!  Modified:
!
!    13 July 2000
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
!    Input, real A(LDA,N), the M by N array.
!
!    Output, real X(M*N), a vector containing the N columns of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer m
  real x(m*n)
!
  i = 1
  do j = 1, n
    x(i:i+m-1) = a(1:m,j)
    i = i + m
  end do

  return
end
subroutine rcol_uniq ( lda, m, n, a, nuniq )
!
!*******************************************************************************
!
!! RCOL_UNIQ keeps only the unique elements in a sorted real array of columns.
!
!
!  Modified:
!
!    09 July 2000
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
!    Input/output, real A(LDA,N).
!    On input, the sorted real array of N columns of real M-vectors.
!    On output, a sorted real array of NUNIQ columns of real M-vectors.
!
!    Output, integer NUNIQ, the number of unique columns.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer isgn
  integer itest
  integer m
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    call rcol_compare ( lda, m, n, a, itest, nuniq, isgn )

    if ( isgn /= 0 ) then
      nuniq = nuniq + 1
      a(1:m,nuniq) = a(1:m,itest)
    end if

  end do

  return
end
subroutine rcol_variance ( lda, m, n, a, variance )
!
!*******************************************************************************
!
!! RCOL_VARIANCE returns the variances of the columns of a real array.
!
!
!  Modified:
!
!    29 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real A(LDA,N), the array whose variances are desired.
!
!    Output, real VARIANCE(N), the variances of the rows.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer j
  real mean
  real variance(n)
!
  do j = 1, n

    mean = sum ( a(1:m,j) ) / real ( m )

    variance(j) = 0.0E+00
    do i = 1, m
      variance(j) = variance(j) + ( a(i,j) - mean )**2
    end do

    if ( m > 1 ) then
      variance(j) = variance(j) / real ( m - 1 )
    else
      variance(j) = 0.0E+00
    end if

  end do

  return
end
subroutine rint_to_iint ( rmin, rmax, r, imin, imax, i )
!
!*******************************************************************************
!
!! RINT_TO_IINT maps a real interval to an integer interval.
!
!
!  Formula:
!
!    I := IMIN + ( IMAX - IMIN ) * ( R - RMIN ) / ( RMAX - RMIN )
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RMIN, RMAX, the real range.
!
!    Input, real R, the real number to be converted.
!
!    Input, integer IMAX, IMIN, the integer range.
!
!    Output, integer I, the corresponding value in the range [IMIN,IMAX].
!
  integer i
  integer imax
  integer imin
  real r
  real rmax
  real rmin
!
  if ( rmax == rmin ) then

    i = ( imax + imin ) / 2

  else

    i = nint ( ( ( rmax - r ) * real ( imin ) + ( r - rmin ) * real ( imax ) ) &
      / ( rmax - rmin ) )

  end if

  return
end
subroutine rint_to_rint ( rmin, rmax, r, r2min, r2max, r2 )
!
!*******************************************************************************
!
!! RINT_TO_RINT maps a real interval to another real interval.
!
!
!  Formula:
!
!    R2 := R2MIN + ( R2MAX - R2MIN ) * ( R - RMIN ) / ( RMAX - RMIN )
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RMIN, RMAX, the first real range.
!
!    Input, real R, the real number to be converted.
!
!    Input, real R2MAX, R2MIN, the second real range.
!
!    Output, real R2, the corresponding value in the range [R2MIN,R2MAX].
!
  real r
  real rmax
  real rmin
  real r2
  real r2max
  real r2min
!
  if ( rmax == rmin ) then

    r2 = ( r2max + r2min ) / 2.0E+00

  else

    r2 = ( ( ( rmax - r ) * r2min + ( r - rmin ) * r2max ) / ( rmax - rmin ) )

  end if

  return
end
subroutine rmat_aat ( lda, lda2, m, n, a, aat )
!
!*******************************************************************************
!
!! RMAT_AAT computes the normal matrix A*A'.
!
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the M by N matrix A.
!
!    Input, integer LDA2, the leading dimension of the M by M matrix AAT.
!
!    Input, integer M, N, the number of rows and columns of the matrix A.
!
!    Input, real A(LDA,N), the M by N matrix.
!
!    Output, real AAT(LDA2,M), the M by M normal matrix A*A'.
!
  integer lda
  integer lda2
  integer m
  integer n
!
  real a(lda,n)
  real aat(lda2,m)
  integer i
  integer j
  integer k
  double precision sum
!
  do i = 1, m

    do j = 1, m

      sum = 0.0D+00
      do k = 1, m
        sum = sum + dble ( a(i,k) ) * dble ( a(j,k) )
      end do

      aat(i,j) = sngl ( sum )

    end do
  end do

  return
end
subroutine rmat_ata ( lda, lda2, m, n, a, ata )
!
!*******************************************************************************
!
!! RMAT_ATA computes the normal matrix A'*A.
!
!
!  Modified:
!
!    19 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the M by N matrix A.
!
!    Input, integer LDA2, the leading dimension of the N by N matrix ATA.
!
!    Input, integer M, N, the number of rows and columns of the matrix A.
!
!    Input, real A(LDA,N), the M by N matrix.
!
!    Output, real ATA(LDA2,N), the N by N normal matrix A'*A.
!
  integer lda
  integer lda2

  integer n
!
  real a(lda,n)
  real ata(lda2,n)
  integer i
  integer j
  integer k
  integer m
  double precision sum
!
  do i = 1, n

    do j = 1, n

      sum = 0.0D+00
      do k = 1, m
        sum = sum + dble ( a(k,i) ) * dble ( a(k,j) )
      end do

      ata(i,j) = sngl ( sum )

    end do
  end do

  return
end
subroutine rmat_cholesky_factor ( lda, n, a, c, ierror )
!
!*******************************************************************************
!
!! RMAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is a lower triangular matrix L such that:
!
!      A = L * L'
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the N by N matrix A.
!
!    Input, integer N, the number of rows and columns of the matrix A.
!
!    Input, real A(LDA,N), the N by N matrix.
!
!    Output, real C(LDA,N), the N by N lower triangular Cholesky factor.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, warning, the matrix is positive semidefinite.  The factorization
!    was carried out, but the matrix is singular.
!    2, error, the matrix has at least one negative eigenvalue.  The
!    factorization could not be completed.
!
  integer lda
  integer n
!
  real a(lda,n)
  real c(lda,n)
  integer i
  integer ierror
  integer j
  integer k
  real sum2
!
  ierror = 0

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0E+00

    do i = j, n

      sum2 = c(j,i)
      do k = 1, j-1
        sum2 = sum2 - c(j,k) * c(i,k)
      end do

      if ( i == j ) then
        if ( sum2 < 0.0E+00 ) then
          ierror = 2
          return
        else if ( sum2 == 0.0E+00 ) then
          ierror = 1
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0E+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0E+00
        end if
      end if

    end do

  end do

  return
end
subroutine rmat_cholesky_solve ( lda, n, a, b, x )
!
!*******************************************************************************
!
!! RMAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
!
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the N by N matrix A.
!
!    Input, integer N, the number of rows and columns of the matrix A.
!
!    Input, real A(LDA,N), the N by N Cholesky factor of the system matrix.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution of the linear system.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  real x(n)
!
!  Solve L * y = b.
!
  call rmat_l_solve ( lda, n, a, b, x )
!
!  Solve L' * x = y.
!
  call rmat_lt_solve ( lda, n, a, x, x )

  return
end
function rmat_det_2d ( a )
!
!*******************************************************************************
!
!! RMAT_DET_2D computes the determinant of a 2 by 2 matrix.
!
!
!  Formula:
!
!    The determinant of a 2 by 2 matrix is
!
!      a11 * a22 - a12 * a21.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(2,2), the matrix whose determinant is desired.
!
!    Output, real RMAT_DET_2D, the determinant of the matrix.
!
  real a(2,2)
  real rmat_det_2d
!
  rmat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

  return
end
function rmat_det_3d ( a )
!
!*******************************************************************************
!
!! RMAT_DET_3D computes the determinant of a 3 by 3 matrix.
!
!
!  Formula:
!
!    The determinant of a 3 by 3 matrix is
!
!        a11 * a22 * a33 - a11 * a23 * a32
!      + a12 * a23 * a31 - a12 * a21 * a33
!      + a13 * a21 * a32 - a13 * a22 * a31
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the matrix whose determinant is desired.
!
!    Output, real RMAT_DET_3D, the determinant of the matrix.
!
  real a(3,3)
  real rmat_det_3d
!
  rmat_det_3d = &
         a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
function rmat_det_4d ( a )
!
!*******************************************************************************
!
!! RMAT_DET_4D computes the determinant of a 4 by 4 matrix.
!
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(4,4), the matrix whose determinant is desired.
!
!    Output, real RMAT_DET_4D, the determinant of the matrix.
!
  real a(4,4)
  real rmat_det_4d
!
  rmat_det_4d = &
         a(1,1) * ( &
             a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
             a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
             a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
             a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
           + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
function rmat_det_5d ( a )
!
!*******************************************************************************
!
!! RMAT_DET_5D computes the determinant of a 5 by 5 matrix.
!
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(5,5), the matrix whose determinant is desired.
!
!    Output, real RMAT_DET_5D, the determinant of the matrix.
!
  real a(5,5)
  real b(4,4)
  integer i
  integer inc
  integer j
  integer k
  real rmat_det_4d
  real rmat_det_5d
!
!  Expand the determinant into the sum of the determinants of the
!  five 4 by 4 matrices created by dropping row 1, and column k.
!
  rmat_det_5d = 0.0E+00

  do k = 1, 5

    do i = 1, 4
      do j = 1, 4

        if ( j < k ) then
          inc = 0
        else
          inc = 1
        end if

        b(i,j) = a(i+1,j+inc)

      end do
    end do

    rmat_det_5d = rmat_det_5d + (-1)**( k + 1 ) * a(1,k) * rmat_det_4d ( b )

  end do

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
subroutine rmat_diag_set_vector ( lda, n, a, v )
!
!*******************************************************************************
!
!! RMAT_DIAG_SET_VECTOR sets the diagonal of a matrix to a vector.
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
!    Input/output, real A(LDA,N), the N by N matrix.
!
!    Input, real V(N), the vector to be assigned to the diagonal of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  real v(n)
!
  do i = 1, n
    a(i,i) = v(i)
  end do

  return
end
subroutine rmat_diag_add_vector ( lda, n, a, v )
!
!*******************************************************************************
!
!! RMAT_DIAG_ADD_VECTOR adds a vector to the diagonal of a matrix.
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
!    Input/output, real A(LDA,N), the N by N matrix.
!
!    Input, real V(N), the vector to be added to the diagonal of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  real v(n)
!
  do i = 1, n
    a(i,i) = a(i,i) + v(i)
  end do

  return
end
subroutine rmat_expand_linear ( lda, m, n, a, ldb, mb, nb, b )
!
!*******************************************************************************
!
!! RMAT_EXPAND_LINEAR expands a real array by linear interpolation.
!
!
!  Modified:
!
!    30 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the matrix A.
!
!    Input, integer M, N, the number of rows and columns in A.
!
!    Input, real A(LDA,N), a "small" M by N array.
!
!    Input, integer LDB, the leading dimension of the matrix B.
!
!    Input, integer MB, NB, the number of rows and columns in B.
!
!    Output, real B(LDB,NB), the "big" MB by NB array, which contains an
!    interpolated version of the data in SMALL.
!
  integer lda
  integer ldb
  integer n
  integer nb
!
  real a(lda,n)
  real b(ldb,nb)
  integer i
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer m
  integer mb
  real r
  real r1
  real r2
  real s
  real s1
  real s2
!
  do i = 1, mb

    if ( mb == 1 ) then
      r = 0.5E+00
    else
      r = real ( i - 1 ) / real ( mb - 1 )
    end if

    i1 = 1 + int ( r * ( m - 1 ) )
    i2 = i1 + 1
    if ( i2 > m ) then
      i1 = m - 1
      i2 = m
    end if

    r1 = real ( i1 - 1 ) / real ( m - 1 )
    r2 = real ( i2 - 1 ) / real ( m - 1 )

    do j = 1, nb

      if ( nb == 1 ) then
        s = 0.5E+00
      else
        s = real ( j - 1 ) / real ( nb - 1 )
      end if

      j1 = 1 + int ( s * real ( n - 1 ) )
      j2 = j1 + 1
      if ( j2 > n ) then
        j1 = n - 1
        j2 = n
      end if

      s1 = real ( j1 - 1 ) / real ( n - 1 )
      s2 = real ( j2 - 1 ) / real ( n - 1 )

      b(i,j) = &
        ( ( r2 - r ) * ( s2 - s ) * a(i1,j1) &
        + ( r - r1 ) * ( s2 - s ) * a(i2,j1) &
        + ( r2 - r ) * ( s - s1 ) * a(i1,j2) &
        + ( r - r1 ) * ( s - s1 ) * a(i2,j2) ) &
        / ( ( r2 - r1 ) * ( s2 - s1 ) )

    end do

  end do

  return
end
subroutine rmat_givens_post ( lda, n, a, g, irow, jcol )
!
!*******************************************************************************
!
!! RMAT_GIVENS_POST computes the Givens postmultiplier rotation matrix G(IROW,JCOL).
!
!
!  Discussion:
!
!    G(IROW,JCOL) has the property that the (IROW,JCOL)-th entry of A*G is zero.
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
!    Input, integer LDA, the leading dimension of the arrays.
!
!    Input, integer N, the order of the matrices A and G.
!
!    Input, real A(LDA,N), the matrix to be operated upon.
!
!    Output, real G(LDA,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
!    Input, integer IROW, JCOL, the row and column of the
!    entry of A*G which is to be zeroed out.
!
  integer lda
  integer n
!
  real a(lda,n)
  real g(lda,n)
  integer irow
  integer jcol
  real theta
!
  call rmat_identity ( lda, n, g )

  theta = atan2 ( a(irow,jcol), a(irow,irow) )

  g(irow,irow) =   cos ( theta )
  g(irow,jcol) = - sin ( theta )
  g(jcol,irow) =   sin ( theta )
  g(jcol,jcol) =   cos ( theta )

  return
end
subroutine rmat_givens_pre ( lda, n, a, g, irow, jcol )
!
!*******************************************************************************
!
!! RMAT_GIVENS_PRE computes the Givens premultiplier rotation matrix G(IROW,JCOL).
!
!
!  Discussion:
!
!    G(IROW,JCOL) has the property that the (IROW,JCOL)-th entry of
!    G*A is zero.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the arrays.
!
!    Input, integer N, the order of the matrices A and G.
!
!    Input, real A(LDA,N), the matrix to be operated upon.
!
!    Output, real G(LDA,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
!    Input, integer IROW, JCOL, the row and column of the
!    entry of the G*A which is to be zeroed out.
!
  integer lda
  integer n
!
  real a(lda,n)
  real g(lda,n)
  integer irow
  integer jcol
  real theta
!
  call rmat_identity ( lda, n, g )

  theta = atan2 ( a(irow,jcol), a(jcol,jcol) )

  g(irow,irow) =   cos ( theta )
  g(irow,jcol) = - sin ( theta )
  g(jcol,irow) =   sin ( theta )
  g(jcol,jcol) =   cos ( theta )

  return
end
subroutine rmat_house_axh ( lda, n, a, v, ah )
!
!*******************************************************************************
!
!! RMAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
!
!
!  Discussion:
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Modified:
!
!    26 March 2000
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
!    Output, real A(LDA,N), the matrix to be postmultiplied.
!
!    Input, real V(N), a vector defining a Householder matrix.
!
!    Output, real AH(LDA,N), the product A*H.
!
  integer lda
  integer n
!
  real a(lda,n)
  real ah(lda,n)
  real ah_temp(n,n)
  integer i
  integer j
  integer k
  real v(n)
  real v_normsq
!
  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ah_temp(i,j) = a(i,j)
      do k = 1, n
        ah_temp(i,j) = ah_temp(i,j) - 2.0E+00 * a(i,k) * v(k) * v(j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into AH.
!  Doing it this way means the user can identify the input arguments A and AH.
!
  ah(1:n,1:n) = ah_temp(1:n,1:n)

  return
end
subroutine rmat_house_form ( lda, n, v, h )
!
!*******************************************************************************
!
!! RMAT_HOUSE_FORM constructs a Householder matrix from its compact form.
!
!
!  Discussion:
!
!    H(v) = I - 2 * v * v' / ( v' * v )
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
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real V(N), the vector defining the Householder matrix.
!
!    Output, real H(LDA,N), the Householder matrix.
!
  integer lda
  integer n
!
  real beta
  real h(lda,n)
  integer i
  integer j
  real v(n)
!
!  Compute the L2 norm of V.
!
  beta = sum ( v(1:n)**2 )
!
!  Form the matrix H.
!
  call rmat_identity ( lda, n, h )

  do i = 1, n
    do j = 1, n
      h(i,j) = h(i,j) - 2.0E+00 * v(i) * v(j) / beta
    end do
  end do

  return
end
subroutine rmat_house_hxa ( lda, n, a, v, ha )
!
!*******************************************************************************
!
!! RMAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
!
!
!  Discussion:
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Modified:
!
!    26 March 2000
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
!    Output, real A(LDA,N), the matrix to be premultiplied.
!
!    Input, real V(N), a vector defining a Householder matrix.
!
!    Output, real HA(LDA,N), the product H*A.
!
  integer lda
  integer n
!
  real a(lda,n)
  real ha(lda,n)
  real ha_temp(n,n)
  integer i
  integer j
  integer k
  real v(n)
  real v_normsq
!
  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ha_temp(i,j) = a(i,j)
      do k = 1, n
        ha_temp(i,j) = ha_temp(i,j) - 2.0E+00 * v(i) * v(k) * a(k,j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into HA.
!  Doing it this way means the user can identify the input arguments A and HA.
!
  ha(1:n,1:n) = ha_temp(1:n,1:n)

  return
end
subroutine rmat_house_post ( lda, n, a, h, irow, jcol )
!
!*******************************************************************************
!
!! RMAT_HOUSE_POST computes a Householder post-multiplier matrix.
!
!
!  Discussion:
!
!    H(IROW,JCOL) has the property that the IROW-th column of
!    A*H(IROW,JCOL) is zero from entry JCOL+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    IROW = JCOL.
!
!  Modified:
!
!    25 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the arrays.
!
!    Input, integer N, the order of the matrices.
!
!    Input, real A(LDA,N), the matrix whose Householder matrix is to be
!    computed.
!
!    Output, real H(LDA,N), the Householder matrix.
!
!    Input, integer IROW, JCOL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same row, but higher column, will be zeroed out if
!    A is postmultiplied by H.
!
  integer lda
  integer n
!
  real a(lda,n)
  real h(lda,n)
  integer irow
  integer j
  integer jcol
  real v(n)
  real w(n)
!
!  Set up the vector V.
!
  w(1:jcol-1) = 0.0E+00
  w(jcol:n) = a(irow,jcol:n)

  call rvec_house_column ( n, w, jcol, v )
!
!  Form the matrix H(V).
!
  call rmat_house_form ( lda, n, v, h )

  return
end
subroutine rmat_house_pre ( lda, n, a, h, irow, jcol )
!
!*******************************************************************************
!
!! RMAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
!
!
!  Discussion:
!
!    H(IROW,JCOL) has the property that the JCOL-th column of
!    H(IROW,JCOL)*A is zero from entry IROW+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    IROW = JCOL.
!
!  Modified:
!
!    25 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the arrays.
!
!    Input, integer N, the order of the matrices.
!
!    Input, real A(LDA,N), the matrix whose Householder matrix is to be
!    computed.
!
!    Output, real H(LDA,N), the Householder matrix.
!
!    Input, integer IROW, JCOL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same column, but higher rows, will be zeroed out if A is
!    premultiplied by H.
!
  integer lda
  integer n
!
  real a(lda,n)
  real h(lda,n)
  integer irow
  integer jcol
  real v(n)
  real w(n)
!
!  Set up the vector V.
!
  w(1:irow-1) = 0.0E+00
  w(irow:n) = a(irow:n,jcol)

  call rvec_house_column ( n, w, irow, v )
!
!  Form the matrix H(V).
!
  call rmat_house_form ( lda, n, v, h )

  return
end
subroutine rmat_identity ( lda, n, a )
!
!*******************************************************************************
!
!! RMAT_IDENTITY sets the square matrix A to the identity.
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
!    Output, real A(LDA,N), the N by N identity matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
!
  a(1:n,1:n) = 0.0E+00

  do i = 1, n
    a(i,i) = 1.0E+00
  end do

  return
end
subroutine rmat_imax ( lda, m, n, a, i, j )
!
!*******************************************************************************
!
!! RMAT_IMAX returns the location of the maximum of a real M by N matrix.
!
!
!  Modified:
!
!    17 July 2000
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
!    Input, real A(LDA,N), the M by N matrix.
!
!    Output, integer I, J, the indices of the maximum entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ii
  integer j
  integer jj
  integer m
!
  i = 0
  j = 0

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(ii,jj) > a(i,j) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
subroutine rmat_imin ( lda, m, n, a, i, j )
!
!*******************************************************************************
!
!! RMAT_IMIN returns the location of the minimum of a real M by N matrix.
!
!
!  Modified:
!
!    17 July 2000
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
!    Input, real A(LDA,N), the M by N matrix.
!
!    Output, integer I, J, the indices of the minimum entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer ii
  integer j
  integer jj
  integer m
!
  i = 0
  j = 0

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(ii,jj) < a(i,j) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
subroutine rmat_inverse_2d ( a, b, det )
!
!*******************************************************************************
!
!! RMAT_INVERSE_2D inverts a 2 by 2 real matrix using Cramer's rule.
!
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(2,2), the matrix to be inverted.
!
!    Output, real B(2,2), the inverse of the matrix A.
!
!    Output, real DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  real a(2,2)
  real b(2,2)
  real det
  integer i
  integer j
  real rmat_det_2d
!
!  Compute the determinant of A.
!
  det = rmat_det_2d ( a )

  if ( det == 0.0E+00 ) then

    b(1:2,1:2) = 0.0E+00

  else

    b(1,1) = + a(2,2) / det
    b(1,2) = - a(1,2) / det
    b(2,1) = - a(2,1) / det
    b(2,2) = + a(1,1) / det

  end if

  return
end
subroutine rmat_inverse_3d ( a, b, det )
!
!*******************************************************************************
!
!! RMAT_INVERSE_3D inverts a 3 by 3 real matrix using Cramer's rule.
!
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the matrix to be inverted.
!
!    Output, real B(3,3), the inverse of the matrix A.
!
!    Output, real DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  real a(3,3)
  real b(3,3)
  real det
  real rmat_det_3d
!
!  Compute the determinant of A.
!
  det = rmat_det_3d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0E+00 ) then
    b(1:3,1:3) = 0.0E+00
    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit
!  formula.
!
  b(1,1) = + ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
  b(1,2) = - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
  b(1,3) = + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

  b(2,1) = - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
  b(2,2) = + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
  b(2,3) = - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

  b(3,1) = + ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
  b(3,2) = - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
  b(3,3) = + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

  return
end
subroutine rmat_inverse_4d ( a, b, det )
!
!*******************************************************************************
!
!! RMAT_INVERSE_4D inverts a 4 by 4 real matrix using Cramer's rule.
!
!
!  Modified:
!
!    13 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(4,4), the matrix to be inverted.
!
!    Output, real B(4,4), the inverse of the matrix A.
!
!    Output, real DET, the determinant of the matrix A.
!
  real a(4,4)
  real b(4,4)
  real det
  integer i
  integer j
  real rmat_det_4d
!
!  Compute the determinant of A.
!
  det = rmat_det_4d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0E+00 ) then

    b(1:4,1:4) = 0.0E+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = +( &
        + a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,1) = -( &
        + a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,1) = +( &
        + a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(2,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,1) = -( &
        + a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(2,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,2) = -( &
        + a(1,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(1,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,2) = +( &
        + a(1,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,2) = -( &
        + a(1,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(1,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,2) = +( &
        + a(1,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(1,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(1,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,3) = +( &
        + a(1,2) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,2) - a(2,2) * a(4,4) ) &
        + a(1,4) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        ) / det

  b(2,3) = -( &
        + a(1,1) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,3) - a(2,3) * a(4,1) ) &
        ) / det

  b(3,3) = +( &
        + a(1,1) * ( a(2,2) * a(4,4) - a(2,4) * a(4,2) ) &
        + a(1,2) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(4,3) = -( &
        + a(1,1) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        + a(1,2) * ( a(2,3) * a(4,1) - a(2,1) * a(4,3) ) &
        + a(1,3) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(1,4) = -( &
        + a(1,2) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,2) - a(2,2) * a(3,4) ) &
        + a(1,4) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        ) / det

  b(2,4) = +( &
        + a(1,1) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) &
        ) / det

  b(3,4) = -( &
        + a(1,1) * ( a(2,2) * a(3,4) - a(2,4) * a(3,2) ) &
        + a(1,2) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  b(4,4) = +( &
        + a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
        + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  return
end
subroutine rmat_l1_inverse ( lda, n, a, b )
!
!*******************************************************************************
!
!! RMAT_L1_INVERSE inverts a real unit lower triangular matrix.
!
!
!  Discussion:
!
!    A unit lower triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's above the main diagonal.
!
!    The inverse of a unit lower triangular matrix is also
!    a unit lower triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call rmat_l1_inverse ( lda, n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    23 March 2000
!
!  Parameters:
!
!    Input, integer LDA, the declared first dimension of the arrays.
!    LDA must be at least N.
!
!    Input, integer N, number of rows and columns in the matrix.
!
!    Input, real A(LDA,N), the unit lower triangular matrix.
!
!    Output, real B(LDA,N), the inverse matrix.
!
  integer n
  integer lda
!
  real a(lda,n)
  real b(lda,n)
  integer i
  integer j
!
  do i = 1, n

    do j = 1, n

      if ( j > i ) then
        b(i,j) = 0.0E+00
      else if ( j == i ) then
        b(i,j) = 1.0E+00
      else
        b(i,j) = - dot_product ( a(i,1:i-1), b(1:i-1,j) )
      end if

    end do
  end do

  return
end
subroutine rmat_l_inverse ( lda, n, a, b )
!
!*******************************************************************************
!
!! RMAT_L_INVERSE inverts a real lower triangular matrix.
!
!
!  Discussion:
!
!    A lower triangular matrix is a matrix whose only nonzero entries
!    occur on or below the diagonal.
!
!    The inverse of a lower triangular matrix is a lower triangular matrix.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    23 March 2000
!
!  Parameters:
!
!    Input, integer LDA, the declared first dimension of the arrays.
!    LDA must be at least N.
!
!    Input, integer N, number of rows and columns in the matrix.
!
!    Input, real A(LDA,N), the lower triangular matrix.
!
!    Output, real B(LDA,N), the inverse matrix.
!
  integer n
  integer lda
!
  real a(lda,n)
  real b(lda,n)
  integer i
  integer j
!
  do j = 1, n

    do i = 1, n

      if ( j > i ) then
        b(i,j) = 0.0E+00
      else if ( j == i ) then
        b(i,j) = 1.0E+00 / a(i,j)
      else
        b(i,j) = - dot_product ( a(i,1:i-1), b(1:i-1,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine rmat_l_solve ( lda, n, a, b, x )
!
!*******************************************************************************
!
!! RMAT_L_SOLVE solves a lower triangular linear system.
!
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
!    Input, integer LDA, the leading dimension of the N by N matrix A.
!
!    Input, integer N, the number of rows and columns of the matrix A.
!
!    Input, real A(LDA,N), the N by N lower triangular matrix.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution of the linear system.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  real x(n)
!
!  Solve L * x = b.
!
  do i = 1, n
    x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
  end do

  return
end
subroutine rmat_lt_solve ( lda, n, a, b, x )
!
!*******************************************************************************
!
!! RMAT_LT_SOLVE solves a transposed lower triangular linear system.
!
!
!  Discussion:
!
!    Given the lower triangular matrix A, the linear system to be solved is:
!
!      A' * x = b
!
!  Modified:
!
!    20 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the N by N matrix A.
!
!    Input, integer N, the number of rows and columns of the matrix A.
!
!    Input, real A(LDA,N), the N by N lower triangular matrix.
!
!    Input, real B(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution of the linear system.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer i
  real x(n)
!
!  Solve L'*x = b.
!
  do i = n, 1, -1
    x(i) = ( b(i) - dot_product ( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
  end do

  return
end
subroutine rmat_lu ( lda, m, n, a, l, p, u )
!
!*******************************************************************************
!
!! RMAT_LU computes the LU factorization of a rectangular matrix.
!
!
!  Discussion:
!
!    The routine is given an M by N matrix A, and produces
!
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix, and
!      P, an M by M permutation matrix P,
!
!    so that
!
!      A = TRANS(P)*L*U.
!
!  Modified:
!
!    01 March 1999
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
!    Output, real P(LDA,M), the M by M permutation matrix.
!
!    Output, real U(LDA,N), the M by N upper triangular factor.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer ipiv
  integer j
  integer jj
  real l(lda,m)
  real p(lda,m)
  real pivot
  real u(lda,n)
!
!  Initialize:
!
!    U:=A
!    L:=Identity
!    P:=Identity
!
  u(1:m,1:n) = a(1:m,1:n)

  call rmat_identity ( lda, m, l )

  p(1:m,1:m) = l(1:m,1:m)
!
!  On step J, find the pivot row, IPIV, and the pivot value PIVOT.
!
  do j = 1, min ( m - 1, n )

    pivot = 0.0E+00
    ipiv = 0

    do i = j, m

      if ( abs ( u(i,j) ) > pivot ) then
        pivot = abs ( u(i,j) )
        ipiv = i
      end if

    end do
!
!  Unless IPIV is zero, swap rows J and IPIV.
!
    if ( ipiv /= 0 ) then

      call rrow_swap ( lda, m, n, u, j, ipiv )

      call rrow_swap ( lda, m, m, l, j, ipiv )

      call rrow_swap ( lda, m, m, p, j, ipiv )
!
!  Zero out the entries in column J, from row J+1 to M.
!
      do i = j+1, m

        if ( u(i,j) /= 0.0E+00 ) then
          l(i,j) = u(i,j) / u(j,j)
          u(i,j) = 0.0E+00

          do jj = j+1, n
            u(i,jj) = u(i,jj) - l(i,j) * u(j,jj)
          end do

        end if

      end do

    end if

  end do

  return
end
subroutine rmat_mat_mult ( lda, n, a, b, c )
!
!*******************************************************************************
!
!! RMAT_MAT_MULT computes a simple matrix*matrix product C = A*B.
!
!
!  Formula:
!
!    C(I,J) = SUM ( K = 1, N) A(I,K) * B(K,J)
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
!    Input, integer LDA, the leading dimension of A, B and C.
!
!    Input, integer N, the order of A, B and C.
!
!    Input, real A(LDA,N), B(LDA,N), two N by N matrices to multiply.
!
!    Output, real C(LDA,N), the product of A and B.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
!
  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  return
end
subroutine rmat_mat_mult2 ( lda, ldb, ldc, l, m, n, a, b, c )
!
!*******************************************************************************
!
!! RMAT_MAT_MULT2 computes the general matrix*matrix product C = A*B.
!
!
!  Formula:
!
!    C(I,J) = SUM (K=1 to M) A(I,K) * B(K,J), I=1 to L, J=1 to N
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
!    Input, integer LDA, LDB, LDC, the leading dimensions of A, B and C.
!
!    Input, integer L, the number of rows of A, and of the product C.
!
!    Input, integer M, the number of columns of A, and rows of B.
!
!    Input, integer N, the number of columns of B, and of the product C.
!
!    Input, real A(LDA,M), the L by M matrix.
!
!    Input, real B(LDB,N), the M by N matrix.
!
!    Output, real C(LDC,N), the L by M product of A and B.
!
  integer lda
  integer ldb
  integer ldc
  integer m
  integer n
!
  real a(lda,m)
  real b(ldb,n)
  real c(ldc,n)
  integer l
!
  c(1:l,1:n) = matmul ( a(1:l,1:m), b(1:m,1:n) )

  return
end
function rmat_max ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_MAX returns the maximum entry of a real M by N matrix.
!
!
!  Modified:
!
!    17 January 2000
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
!    Input, real A(LDA,N), the M by N matrix.
!
!    Output, real RMAT_MAX, the maximum entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer m
  real rmat_max
!
  rmat_max = maxval ( a(1:m,1:n) )

  return
end
function rmat_maxcol_minrow ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_MAXCOL_MINROW returns the maximum column minimum row of an M by N matrix.
!
!
!  Definition:
!
!    RMAT_MAXCOL_MINROW = max ( 1 <= I <= N ) ( min ( 1 <= J <= M ) A(I,J) )
!
!  Discussion:
!
!    For a given matrix, RMAT_MAXCOL_MINROW <= RMAT_MINROW_MAXCOL
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
!    Input, real A(LDA,N), the matrix.
!
!    Output, real RMAT_MAXCOL_MINROW, the maximum column minimum row entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer m
  real rmat_maxcol_minrow
  real rmat_minrow
!
  rmat_maxcol_minrow = 0.0E+00

  do i = 1, m

    rmat_minrow = minval ( a(i,1:n) )

    if ( i == 1 ) then
      rmat_maxcol_minrow = rmat_minrow
    else
      rmat_maxcol_minrow = max ( rmat_maxcol_minrow, rmat_minrow )
    end if

  end do

  return
end
function rmat_maxrow_mincol ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_MAXROW_MINCOL returns the maximum row minimum column of an M by N matrix.
!
!
!  Definition:
!
!    RMAT_MAXROW_MINCOL = max ( 1 <= J <= N ) ( min ( 1 <= I <= M ) A(I,J) )
!
!  Discussion:
!
!    For a given matrix, RMAT_MAXROW_MINCOL <= RMAT_MINCOL_MAXROW
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
!    Input, real A(LDA,N), the matrix.
!
!    Output, real RMAT_MAXROW_MINCOL, the maximum row minimum column entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer m
  real rmat_maxrow_mincol
  real rmat_mincol
!
  rmat_maxrow_mincol = 0.0E+00

  do j = 1, n

    rmat_mincol = minval ( a(1:m,j) )

    if ( j == 1 ) then
      rmat_maxrow_mincol = rmat_mincol
    else
      rmat_maxrow_mincol = max ( rmat_maxrow_mincol, rmat_mincol )
    end if

  end do

  return
end
function rmat_min ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_MIN returns the minimum entry of an M by N matrix.
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
!    Input, real A(LDA,N), the matrix.
!
!    Output, real RMAT_MIN, the minimum entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer m
  real rmat_min
!
  rmat_min = minval ( a(1:m,1:n) )

  return
end
function rmat_mincol_maxrow ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_MINCOL_MAXROW returns the minimum column maximum row of an M by N matrix.
!
!
!  Definition:
!
!    RMAT_MINCOL_MAXROW = min ( 1 <= I <= N ) ( max ( 1 <= J <= M ) A(I,J) )
!
!  Discussion:
!
!    For a given matrix, RMAT_MAXROW_MINCOL <= RMAT_MINCOL_MAXROW
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
!    Input, real A(LDA,N), the matrix.
!
!    Output, real RMAT_MINCOL_MAXROW, the minimum column maximum row entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer m
  real rmat_mincol_maxrow
  real rmat_maxrow
!
  rmat_mincol_maxrow = 0.0E+00

  do i = 1, m

    rmat_maxrow = maxval ( a(i,1:n) )

    if ( i == 1 ) then
      rmat_mincol_maxrow = rmat_maxrow
    else
      rmat_mincol_maxrow = min ( rmat_mincol_maxrow, rmat_maxrow )
    end if

  end do

  return
end
function rmat_minrow_maxcol ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_MINROW_MAXCOL returns the minimum row maximum column of an M by N matrix.
!
!
!  Definition:
!
!    RMAT_MINROW_MAXCOL = min ( 1 <= J <= N ) ( max ( 1 <= I <= M ) A(I,J) )
!
!  Discussion:
!
!    For a given matrix, RMAT_MAXCOL_MINROW <= RMAT_MINROW_MAXCOL
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
!    Input, real A(LDA,N), the matrix.
!
!    Output, real RMAT_MINROW_MAXCOL, the minimum row maximum column entry of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer m
  real rmat_minrow_maxcol
  real rmat_maxcol
!
  rmat_minrow_maxcol = 0.0E+00

  do j = 1, n

    rmat_maxcol = maxval ( a(1:m,j) )

    if ( j == 1 ) then
      rmat_minrow_maxcol = rmat_maxcol
    else
      rmat_minrow_maxcol = min ( rmat_minrow_maxcol, rmat_maxcol )
    end if

  end do

  return
end
subroutine rmat_mtv ( lda, m, n, a, x, y )
!
!*******************************************************************************
!
!! RMAT_MTV multiplies a transposed matrix times a vector
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
!    Input, integer LDA, the leading dimension of the M by N matrix.
!
!    Input, integer M, N, the number of rows and columns of the matrix.
!
!    Input/output, real A(LDA,N), the M by N matrix to be modified.
!
!    Input, real X(M), the vector to be multiplied by A.
!
!    Output, real Y(N), the product A'*X.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer j
  real x(m)
  double precision xa
  real y(n)
!
  do i = 1, n
    xa = 0.0D+00
    do j = 1, m
      xa = xa + dble ( x(j) ) * dble ( a(j,i) )
    end do
    y(i) = sngl ( xa )
  end do

  return
end
subroutine rmat_mv ( lda, m, n, a, x, y )
!
!*******************************************************************************
!
!! RMAT_MV multiplies a matrix times a vector.
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
!    Input, integer LDA, the leading dimension of the M by N matrix.
!
!    Input, integer M, N, the number of rows and columns of the matrix.
!
!    Input/output, real A(LDA,N), the M by N matrix to be modified.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real Y(M), the product A*X.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  double precision ax
  integer i
  integer j
  real x(n)
  real y(m)
!
  do i = 1, m
    ax = 0.0D+00
    do j = 1, n
      ax = ax + dble ( a(i,j) ) * dble ( x(j) )
    end do
    y(i) = sngl ( ax )
  end do

  return
end
function rmat_norm1 ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_NORM1 returns the matrix 1-norm of an M by N matrix.
!
!
!  Definition:
!
!    The matrix 1-norm is defined as:
!
!      RMAT_NORM1 =  Max ( 1 <= J <= N )
!        Sum ( 1 <= I <= M ) abs ( A(I,J) ).
!
!  Note:
!
!    The matrix 1-norm is derived from the vector 1-norm, and
!    satisifies:
!
!      Norm1 ( A*x ) <= Norm1 ( A ) * Norm1 ( x ).
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
!    Input, real A(LDA,N), the matrix whose 1-norm is desired.
!
!    Output, real RMAT_NORM1, the 1-norm of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer j
  integer m
  real rmat_norm1
!
  rmat_norm1 = 0.0E+00

  do j = 1, n
    rmat_norm1 = max ( rmat_norm1, sum ( abs ( a(1:m,j) ) ) )
  end do

  return
end
function rmat_norm2 ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_NORM2 returns the matrix 2-norm of an M by N matrix.
!
!
!  Definition:
!
!    The matrix 2-norm is defined as:
!
!      RMAT_NORM2 =  Sqrt ( Max ( 1 <= I <= M ) LAMBDA(I) )
!
!    where LAMBDA contains the eigenvalues of A * Transpose(A).
!
!  Note:
!
!    The matrix 2-norm is derived from the vector 2-norm, and
!    satisifies:
!
!      Norm2 ( A*x ) <= Norm2 ( A ) * Norm2 ( x ).
!
!  Modified:
!
!    19 March 2001
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
!    Input, real A(LDA,N), the matrix whose 2-norm is desired.
!
!    Output, real RMAT_NORM2, the 2-norm of A.  
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(m,m)
  real diag(m)
  real rmat_norm2
!
!  Compute B = A * Transpose(A).
!
  b(1:m,1:m) = matmul ( a(1:m,1:n), transpose ( a(1:m,1:n) ) )
!
!  Diagonalize B.
!
  call rmat_symm_jacobi ( m, m, b )
!
!  Find the maximum eigenvalue, and take its square root.
!
  call rmat_diag_get_vector ( m, m, b, diag )

  rmat_norm2 = sqrt ( maxval ( diag(1:m) ) )

  return
end
function rmat_norme ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_NORME returns the EISPACK norm of an M by N matrix.
!
!
!  Definition:
!
!    The EISPACK norm is defined as:
!
!      RMAT_NORME =  Sum i = 1 to M ( Sum j = 1 to N ( Abs ( A(I,J) ) ) )
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
!    Input, real A(LDA,N), the matrix whose EISPACK norm is desired.
!
!    Output, real RMAT_NORME, the EISPACK norm of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer m
  real rmat_norme
!
  rmat_norme = sum ( abs ( a(1:m,1:n) ) )

  return
end
function rmat_normf ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_NORMF returns the Frobenius norm of an M by N matrix.
!
!
!  Definition:
!
!    The Frobenius norm is defined as
!
!      RMAT_NORMF = Sqrt ( Sum i = 1 to M ( Sum j = 1 to N ( A(I,J)**2) ) )
!
!  Note:
!
!    The matrix Frobenius-norm is not derived from a vector norm, but
!    is compatible with the vector 2-norm, so that:
!
!      Norm2 ( A*x ) <= Normf ( A ) * Norm2 ( x ).
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
!    Input, real A(LDA,N), the matrix whose Frobenius norm is desired.
!
!    Output, real RMAT_NORMF, the Frobenius norm of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer m
  real rmat_normf
!
  rmat_normf = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
function rmat_normi ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_NORMI returns the matrix infinity-norm of an M by N matrix.
!
!
!  Definition:
!
!    The matrix infinity-norm is defined as:
!
!      RMAT_NORMI =  Max ( 1 <= I <= M ) Sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!  Note:
!
!    The matrix infinity-norm is derived from the vector infinity-norm,
!    and satisifies:
!
!      Normi ( A*x ) <= Normi ( A ) * Normi ( x ).
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
!    Input, real A(LDA,N), the matrix whose infinity-norm is desired.
!
!    Output, real RMAT_NORMI, the infinity-norm of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer m
  real rmat_normi
!
  rmat_normi = 0.0E+00

  do i = 1, m
    rmat_normi = max ( rmat_normi, sum ( abs ( a(i,1:n) ) ) )
  end do

  return
end
subroutine rmat_orth_random ( lda, n, a )
!
!*******************************************************************************
!
!! RMAT_ORTH_RANDOM returns a random orthogonal matrix.
!
!
!  Properties:
!
!    The inverse of A is equal to A'.
!
!    A * A'  = A' * A = I.
!
!    Columns and rows of A have unit Euclidean norm.
!
!    Distinct pairs of columns of A are orthogonal.
!
!    Distinct pairs of rows of A are orthogonal.
!
!    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
!
!    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
!
!    The determinant of A is +1 or -1.
!
!    All the eigenvalues of A have modulus 1.
!
!    All singular values of A are 1.
!
!    All entries of A are between -1 and 1.
!
!  Discussion:
!
!    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
!    National Academy of Sciences of Belarus, for convincingly
!    pointing out the severe deficiencies of an earlier version of
!    this routine.
!
!    Essentially, the computation involves saving the Q factor of the
!    QR factorization of a matrix whose entries are normally distributed.
!    However, it is only necessary to generate this matrix a column at
!    a time, since it can be shown that when it comes time to annihilate
!    the subdiagonal elements of column K, these (transformed) elements of
!    column K are still normally distributed random values.  Hence, there
!    is no need to generate them at the beginning of the process and
!    transform them K-1 times.
!
!    For computational efficiency, the individual Householder transformations
!    could be saved, as recommended in the reference, instead of being
!    accumulated into an explicit matrix format.
!
!  Reference:
!
!    G W Stewart,
!    Efficient Generation of Random Orthogonal Matrices With an Application
!    to Condition Estimators,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 3, June 1980, pages 403-409.
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
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of A.
!
!    Output, real A(LDA,N), the orthogonal matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  real v(n)
  real x(n)
!
!  Start with A = the identity matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 1.0E+00
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do
!
!  Now behave as though we were computing the QR factorization of
!  some other random matrix.  Generate the N elements of the first column,
!  compute the Householder matrix H1 that annihilates the subdiagonal elements,
!  and set A := A * H1' = A * H.
!
!  On the second step, generate the lower N-1 elements of the second column,
!  compute the Householder matrix H2 that annihilates them,
!  and set A := A * H2' = A * H2 = H1 * H2.
!
!  On the N-1 step, generate the lower 2 elements of column N-1,
!  compute the Householder matrix HN-1 that annihilates them, and
!  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
!  This is our random orthogonal matrix.
!
  do j = 1, n-1
!
!  Set the vector that represents the J-th column to be annihilated.
!
    x(1:j-1) = 0.0E+00

    do i = j, n
      call normal_01_sample ( x(i) )
    end do
!
!  Compute the vector V that defines a Householder transformation matrix
!  H(V) that annihilates the subdiagonal elements of X.
!
    call rvec_house_column ( n, x, j, v )
!
!  Postmultiply the matrix A by H'(V) = H(V).
!
    call rmat_house_axh ( lda, n, a, v, a )

  end do

  return
end
subroutine rmat_poly_char ( lda, n, a, p )
!
!*******************************************************************************
!
!! RMAT_POLY_CHAR computes the characteristic polynomial of a matrix.
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
!    Input, integer LDA, the leading dimension of the matrix A.
!
!    Input, integer N, the order of the matrix A.
!
!    Input, real A(LDA,N), the N by N matrix.
!
!    Output, real P(0:N), the coefficients of the characteristic
!    polynomial of A.  P(N) contains the coefficient of X**N
!    (which will be 1), P(I) contains the coefficient of X**I,
!    and P(0) contains the constant term.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer iorder
  real p(0:n)
  real rmat_trace
  real trace
  real work1(n,n)
  real work2(n,n)
!
!  Initialize WORK1 to the identity matrix.
!
  call rmat_identity ( n, n, work1 )

  p(n) = 1.0E+00

  do iorder = n-1, 0, -1
!
!  Work2 = A * WORK1.
!
    work2(1:n,1:n) = matmul ( a(1:n,1:n), work1(1:n,1:n) )
!
!  Take the trace.
!
    trace = rmat_trace ( n, n, work2 )
!
!  P(IORDER) = - Trace ( WORK2 ) / ( N - IORDER )
!
    p(iorder) = - trace / real ( n - iorder )
!
!  WORK1 := WORK2 + P(IORDER) * Identity.
!
    work1(1:n,1:n) = work2(1:n,1:n)

    do i = 1, n
      work1(i,i) = work1(i,i) + p(iorder)
    end do

  end do

  return
end
subroutine rmat_power ( lda, n, a, npow, b )
!
!*******************************************************************************
!
!! RMAT_POWER computes A**NPOW, the nonnegative power of a real square matrix.
!
!
!  Discussion:
!
!    The algorithm is:
!
!      B=I
!      Do NPOW times:
!        B=A*B
!      End do
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A and B.
!
!    Input, integer N, the order of A.
!
!    Input, real A(LDA,N), the matrix to be raised to a power.
!
!    Input, integer NPOW, the power to which A is to be raised.
!
!    Output, real B(LDA,N), the value of A**NPOW.
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(lda,n)
  integer ipow
  integer npow
!
  call rmat_identity ( lda, n, b )

  do ipow = 1, npow
    b(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )
  end do

  return
end
subroutine rmat_power_method ( lda, n, a, r, v )
!
!*******************************************************************************
!
!! RMAT_POWER_METHOD applies the power method to a matrix.
!
!
!  Discussion:
!
!    If the power method has not converged, then calling the routine
!    again immediately with the output from the previous call will
!    continue the iteration.
!
!  Modified:
!
!    06 February 2001
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
!    Input, real A(LDA,N), the matrix.
!
!    Output, real R, V(N), the estimated eigenvalue and eigenvector.
!
  real, parameter :: IT_EPS = 0.0001E+00
  integer, parameter :: IT_MAX = 100
  integer, parameter :: IT_MIN = 10
!
  integer lda
  integer n
!
  real a(lda,n)
  real av(n)
  real eps
  integer i
  integer j
  real r
  real r2
  real r_old
  real v(n)
!
  eps = sqrt ( epsilon ( 1.0E+00 ) )

  r = sqrt ( sum ( v(1:n)**2 ) )

  if ( r == 0.0E+00 ) then
    v(1:n) = 1.0E+00
    r = sqrt ( real ( n ) )
  end if

  v(1:n) = v(1:n) / r

  do i = 1, IT_MAX

    av(1:n) = matmul ( a(1:n,1:n), v(1:n) )

    r_old = r
    r = sqrt ( sum ( av(1:n)**2 ) )

    if ( i > IT_MIN ) then
      if ( abs ( r - r_old ) <= IT_EPS * ( 1.0E+00 + abs ( r ) ) ) then
        exit
      end if
    end if

    v(1:n) = av(1:n)

    if ( r /= 0.0E+00 ) then
      v(1:n) = v(1:n) / r
    end if
!
!  Perturb V a bit, to avoid cases where the initial guess is exactly
!  the eigenvector of a smaller eigenvalue.
!
    if ( i < IT_MAX / 2 ) then
      j = 1 + mod ( i-1, n )
      v(j) = v(j) + eps * ( 1.0E+00 + abs ( v(j) ) )
      r2 = sqrt ( sum ( v(1:n)**2 ) )
      v(1:n) = v(1:n) / r2
    end if

  end do

  return
end
subroutine rmat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! RMAT_PRINT prints a real matrix.
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
subroutine rmat_print2 ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_PRINT2 prints out the M by N matrix A.
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
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A.
!
!    Input, real A(LDA,N), the M by N matrix to be printed.
!
  integer lda
  integer n
!
  real a(lda,n)
  real amax
  real amin
  integer i
  character ( len = 10 ) iform
  integer ihi
  integer ilo
  logical integ
  integer j
  integer jhi
  integer jlo
  integer lmax
  integer log10
  integer m
  integer npline
  logical r_is_int
!
!  Check if all entries are integral.
!
  integ = .true.

  do i = 1, m
    do j = 1, n

      if ( integ ) then
        if ( .not. r_is_int ( a(i,j) ) ) then
          integ = .false.
        end if
      end if

    end do
  end do
!
!  Find the maximum and minimum entries.
!
  amax = maxval ( a(1:m,1:n) )
  amin = minval ( a(1:m,1:n) )
!
!  Use the information about the maximum size of an entry to
!  compute an intelligent format for use with integer entries.
!
!  Later, we might also do this for real matrices.
!
  lmax = log10 ( amax )

  if ( integ ) then
    npline = 79 / ( lmax + 3 )
    write ( iform, '(''('',i2,''I'',i2,'')'')' ) npline, lmax+3
  else
    npline = 5
    iform = ' '
  end if
!
!  Print a scalar quantity.
!
  if ( m == 1 .and. n == 1 ) then

    if ( integ ) then
      write ( *, iform ) int ( a(1,1) )
    else
      write ( *, '(g14.6)' ) a(1,1)
    end if
!
!  Column vector of length M,
!
  else if ( n == 1 ) then

    do ilo = 1, m, npline

      ihi = min ( ilo+npline-1, m )

      if ( integ ) then
        write ( *, iform ) ( int ( a(i,1) ), i = ilo, ihi )
      else
        write ( *, '(5g14.6)' ) a(ilo:ihi,1)
      end if

    end do
!
!  Row vector of length N,
!
  else if ( m == 1 ) then

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( integ ) then
        write ( *, iform ) int ( a(1,jlo:jhi) )
      else
        write ( *, '(5g14.6)' ) a(1,jlo:jhi)
      end if

    end do
!
!  M by N Array
!
  else

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( n > npline ) then
        write ( *, * ) ' '
        write ( *, * ) 'Matrix columns ', jlo, ' to ', jhi
        write ( *, * ) ' '
      end if

      do i = 1, m

        if ( integ ) then
          write ( *, iform ) int ( a(i,jlo:jhi) )
        else
          write ( *, '(5g14.6)' ) a(i,jlo:jhi)
        end if

      end do
    end do

  end if

  return
end
subroutine rmat_print_some ( lda, m, n, a, ihi, ilo, jhi, jlo )
!
!*******************************************************************************
!
!! RMAT_PRINT_SOME prints out a portion of a dense matrix.
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
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real A(LDA,NCOL), contains an NROW by NCOL matrix to be printed.
!
!    Input, integer IHI, ILO, the first and last rows to print.
!
!    Input, integer JHI, JLO, the first and last columns to print.
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
  write ( *, * ) ' '

  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, * ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''Columns '',5a14)' ) ctemp(1:inc)
    write ( *, * ) '  Row'
    write ( *, * ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( r_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, * ) ' '

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
!    04 February 2001
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
subroutine rmat_solve ( a, n, nrhs, info )
!
!*******************************************************************************
!
!! RMAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A(N,N+NRHS), contains in rows and columns 1
!    to N the coefficient matrix, and in columns N+1 through
!    N+NRHS, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Input, integer NRHS, the number of right hand sides.  NRHS
!    must be at least 0.
!
!    Output, integer INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  integer n
  integer nrhs
!
  real a(n,n+nrhs)
  real apivot
  real factor
  integer i
  integer info
  integer ipivot
  integer j
  integer k
!
  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j+1, n
      if ( abs ( a(i,j) ) > abs ( apivot ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0E+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    do i = 1, n + nrhs
      call r_swap ( a(ipivot,i), a(j,i) )
    end do
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0E+00
    a(j,j+1:n+nrhs) = a(j,j+1:n+nrhs) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0E+00
        a(i,j+1:n+nrhs) = a(i,j+1:n+nrhs) - factor * a(j,j+1:n+nrhs)

      end if

    end do

  end do

  return
end
subroutine rmat_symm_jacobi ( lda, n, a )
!
!***************************************************************************************
!
!! RMAT_SYMM_JACOBI applies the Jacobi eigenvalue iteration to a symmetric matrix.
!
!
!  Modified:
!
!    21 July 2000
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
!    Input/output, real A(LDA,N), a symmetric N by N matrix.
!    On output, the matrix has been overwritten by an approximately
!    diagonal matrix, with the eigenvalues on the diagonal.
!
  real, parameter :: EPS = 0.00001E+00
  real, parameter :: MAXIT = 100
!
  integer lda
  integer n
!
  real a(lda,n)
  real c
  integer i
  integer it
  integer j
  integer k
  real normf
  real rmat_normf
  real s
  real sum2
  real t
  real t1
  real t2
  real u
!
  normf = rmat_normf ( lda, n, n, a )

  it = 0

  do

    it = it + 1

    do i = 1, n
      do j = 1, i - 1

        if ( a(i,j) /= 0.0E+00 ) then

          u = ( a(j,j) - a(i,i) ) / ( a(i,j) + a(j,i) )

          t = sign ( 1.0E+00, u ) / ( abs ( u ) + sqrt ( u * u + 1.0E+00 ) )
          c = 1.0E+00 / sqrt ( t * t + 1.0E+00 )
          s = t * c
!
!  A -> A * Q.
!
          do k = 1, n
            t1 = a(i,k)
            t2 = a(j,k)
            a(i,k) = t1 * c - t2 * s
            a(j,k) = t1 * s + t2 * c
          end do
!
!  A -> QT * A
!
          do k = 1, n
            t1 = a(k,i)
            t2 = a(k,j)
            a(k,i) = c * t1 - s * t2
            a(k,j) = s * t1 + c * t2
          end do

        end if
      end do
    end do
!
!  Test the size of the off-diagonal elements.
!
    sum2 = 0.0E+00
    do i = 1, n
      do j = 1, i - 1
        sum2 = sum2 + abs ( a(i,j) )
      end do
    end do

    if ( sum2 <= EPS * ( normf + 1.0E+00 ) .or. it >= MAXIT ) then
      exit
    end if

  end do

  return
end
subroutine rmat_symm_random ( lda, n, a, q, x )
!
!***************************************************************************************
!
!! RMAT_SYMM_RANDOM returns a "random" symmetric matrix with given eigenvalues.
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
!    Input, integer LDA, the leading dimension of A, which must be
!    at least N.
!
!    Input, integer N, the order of A.
!
!    Output, real A(LDA,N), a random symmetric N by N matrix with
!    eigenvalues X and eigenvectors the columns of Q.
!
!    Output, real Q(LDA,N), the eigenvector matrix of A.
!
!    Input, real X(N), the desired eigenvalues for the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer k
  real q(lda,n)
  real x(n)
!
!  Get a random orthogonal matrix Q.
!
  call rmat_orth_random ( lda, n, q )
!
!  Set A = Q * Lambda * Transpose ( Q ).
!
  do i = 1, n
    do j = 1, n
      a(i,j) = 0.0E+00
      do k = 1, n
        a(i,j) = a(i,j) + q(i,k) * x(k) * q(j,k)
      end do
    end do
  end do

  return
end
function rmat_trace ( lda, n, a )
!
!*******************************************************************************
!
!! RMAT_TRACE computes the trace of a real matrix.
!
!
!  Definition:
!
!    The trace of a square matrix is the sum of the diagonal elements.
!
!  Modified:
!
!    20 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(LDA,N), the N by N matrix whose trace is desired.
!
!    Input, integer LDA, the leading dimension of the matrix A.
!
!    Input, integer N, the order of the matrix A.
!
!    Output, real RMAT_TRACE, the trace of the matrix.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  real rmat_trace
!
  rmat_trace = 0.0E+00
  do i = 1, n
    rmat_trace = rmat_trace + a(i,i)
  end do

  return
end
subroutine rmat_transpose ( lda, m, n, a )
!
!*******************************************************************************
!
!! RMAT_TRANSPOSE transposes a real matrix.
!
!
!  Modified:
!
!    28 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!    LDA must be greater than or equal to the maximum of M and N.
!
!    Input, integer M, the number of rows in the input matrix A.
!
!    Input, integer N, the number of columns in the input matrix A.
!
!    Input/output, real A(LDA,*);
!    On input, the M by N matrix to be transposed;
!    On output, the N by M tranposed matrix.
!
  integer lda
!
  real a(lda,*)
  integer i
  integer j
  integer m
  integer n
!
!  Flip the square portion.
!
  do i = 1, min ( m, n )
    do j = 1, i -1
      call r_swap ( a(i,j), a(j,i) )
    end do
  end do
!
!  Fill in the right rectangle, if any.
!
  do j = n+1, m
    a(1:n,j) = a(j,1:n)
  end do
!
!  Fill in the lower rectangle, if any.
!
  do i = m+1, n
    a(i,1:m) = a(1:m,i)
  end do

  return
end
subroutine rmat_u1_inverse ( lda, n, a, b )
!
!*******************************************************************************
!
!! RMAT_U1_INVERSE inverts a real unit upper triangular matrix.
!
!
!  Discussion:
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.
!
!    The inverse of a unit upper triangular matrix is also
!    a unit upper triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call rmat_u1_inverse ( lda, n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    01 March 2001
!
!  Parameters:
!
!    Input, integer LDA, the declared first dimension of the arrays.
!    LDA must be at least N.
!
!    Input, integer N, number of rows and columns in the matrix.
!
!    Input, real A(LDA,N), the unit upper triangular matrix.
!
!    Output, real B(LDA,N), the inverse matrix.
!
  integer n
  integer lda
!
  real a(lda,n)
  real b(lda,n)
  integer i
  integer j
!
  do j = n, 1, -1

    do i = n, 1, -1

      if ( i > j ) then
        b(i,j) = 0.0E+00
      else if ( i == j ) then
        b(i,j) = 1.0E+00
      else
        b(i,j) = - dot_product ( a(i,i+1:j), b(i+1:j,j) )
      end if

    end do
  end do

  return
end
subroutine rmat_u_inverse ( lda, n, a, b )
!
!*******************************************************************************
!
!! RMAT_U_INVERSE inverts a real upper triangular matrix.
!
!
!  Discussion:
!
!    An upper triangular matrix is a matrix whose only nonzero entries
!    occur on or above the diagonal.
!
!    The inverse of an upper triangular matrix is an upper triangular matrix.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    23 March 2000
!
!  Parameters:
!
!    Input, integer LDA, the declared first dimension of the arrays.
!    LDA must be at least N.
!
!    Input, integer N, number of rows and columns in the matrix.
!
!    Input, real A(LDA,N), the upper triangular matrix.
!
!    Output, real B(LDA,N), the inverse matrix.
!
  integer n
  integer lda
!
  real a(lda,n)
  real b(lda,n)
  integer i
  integer j
!
  do j = n, 1, -1

    do i = n, 1, -1

      if ( i > j ) then
        b(i,j) = 0.0E+00
      else if ( i == j ) then
        b(i,j) = 1.0E+00 / a(i,j)
      else
        b(i,j) = - dot_product ( a(i,i+1:n), b(i+1:n,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine rmat_vand2 ( lda, n, a, x )
!
!*******************************************************************************
!
!! RMAT_VAND2 returns the N by N row Vandermonde matrix A.
!
!
!  Discussion:
!
!    The row Vandermonde matrix returned by VAND2 reads "across" rather
!    than down.  In particular, each row begins with a 1, followed by
!    some value X, followed by successive powers of X.
!
!  Formula:
!
!    A(I,J) = X(I)**(J-1)
!
!  Properties:
!
!    A is nonsingular if, and only if, the X values are distinct.
!
!    The determinant of A is
!
!      DET(A) = PRODUCT ( I = 2 to N ) (
!        PRODUCT ( J = 1 to I-1 ) ( ( X(I) - X(J) ) ) ).
!
!    The matrix A is generally ill-conditioned.
!
!  Example:
!
!    N = 5, X = (2, 3, 4, 5, 6)
!
!    1 2  4   8   16
!    1 3  9  27   81
!    1 4 16  64  256
!    1 5 25 125  625
!    1 6 36 216 1296
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array A.
!
!    Input, integer N, the order of the matrix desired.
!
!    Output, real A(LDA,N), the N by N row Vandermonde matrix.
!
!    Input, real X(N), the values that define A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  real x(n)
!
  do i = 1, n
    do j = 1, n

      if ( j == 1 .and. x(i) == 0.0E+00 ) then
        a(i,j) = 1.0E+00
      else
        a(i,j) = x(i)**(j-1)
      end if

    end do
  end do

  return
end
subroutine roots_to_ipoly ( n, x, c )
!
!*******************************************************************************
!
!! ROOTS_TO_IPOLY converts polynomial roots to polynomial coefficients.
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
!    Input, integer N, the number of roots specified.
!
!    Input, integer X(N), the roots.
!
!    Output, integer C(0:N), the coefficients of the polynomial.
!
  integer n
!
  integer c(0:n)
  integer i
  integer j
  integer x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0
  c(n) = 1
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n+1-j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine roots_to_rpoly ( n, x, c )
!
!*******************************************************************************
!
!! ROOTS_TO_RPOLY converts polynomial roots to polynomial coefficients.
!
!
!  Modified:
!
!    17 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of roots specified.
!
!    Input, real X(N), the roots.
!
!    Output, real C(0:N), the coefficients of the polynomial.
!
  integer n
!
  real c(0:n)
  integer i
  integer j
  real x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0.0E+00
  c(n) = 1.0E+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n+1-j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine rpoly2_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )
!
!*******************************************************************************
!
!! RPOLY2_EX finds the extremal point of a parabola determined by three points.
!
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
!    on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real X, Y, the X coordinate of the extremal point of the
!    parabola, and the value of the parabola at that point.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  real bot
  integer ierror
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if

  bot = ( x2 - x3 ) * y1 - ( x1 - x3 ) * y2 + ( x1 - x2 ) * y3

  if ( bot == 0.0E+00 ) then
    ierror = 2
    return
  end if

  x = 0.5E+00 * ( &
          x1**2 * ( y3 - y2 ) &
        + x2**2 * ( y1 - y3 ) &
        + x3**2 * ( y2 - y1 ) ) / bot

  y = ( &
         ( x - x2 ) * ( x - x3 ) * ( x2 - x3 ) * y1 &
       - ( x - x1 ) * ( x - x3 ) * ( x1 - x3 ) * y2 &
       + ( x - x1 ) * ( x - x2 ) * ( x1 - x2 ) * y3 ) / &
       ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) )

  return
end
subroutine rpoly2_ex2 ( x1, y1, x2, y2, x3, y3, x, y, a, b, c, ierror )
!
!*******************************************************************************
!
!! RPOLY2_EX2 finds the extremal point of a parabola determined by three points.
!
!
!  Modified:
!
!    29 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
!    on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real X, Y, the X coordinate of the extremal point of the
!    parabola, and the value of the parabola at that point.
!
!    Output, real A, B, C, the coefficients that define the parabola:
!    P(X) = A * X**2 + B * X + C.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  real a
  real b
  real c
  real det
  integer ierror
  real v(3,3)
  real w(3,3)
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if
!
!  Set up the Vandermonde matrix.
!
  v(1,1) = 1.0E+00
  v(1,2) = x1
  v(1,3) = x1 * x1

  v(2,1) = 1.0E+00
  v(2,2) = x2
  v(2,3) = x2 * x2

  v(3,1) = 1.0E+00
  v(3,2) = x3
  v(3,3) = x3 * x3
!
!  Get the inverse.
!
  call rmat_inverse_3d ( v, w, det )
!
!  Compute the parabolic coefficients.
!
  c = w(1,1) * y1 + w(1,2) * y2 + w(1,3) * y3
  b = w(2,1) * y1 + w(2,2) * y2 + w(2,3) * y3
  a = w(3,1) * y1 + w(3,2) * y2 + w(3,3) * y3
!
!  Determine the extremal point.
!
  if ( a == 0.0E+00 ) then
    ierror = 2
    return
  end if

  x = - b / ( 2.0E+00 * a )
  y = a * x * x + b * x + c

  return
end
subroutine rpoly2_root ( a, b, c, r1, r2 )
!
!*******************************************************************************
!
!! RPOLY2_ROOT returns the two roots of a quadratic polynomial.
!
!
!  Modified:
!
!    02 March 1999
!
!  Parameters:
!
!    Input, real A, B, C, the coefficients of the quadratic polynomial
!    A * X**2 + B * X + C = 0 whose roots are desired.  A must
!    not be zero.
!
!    Output, complex R1, R2, the two roots of the polynomial, which
!    might be real and distinct, real and equal, or complex conjugates.
!
  real a
  real b
  real c
  complex disc
  complex q
  complex r1
  complex r2
!
  if ( a == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY2_ROOT - Fatal error!'
    write ( *, * ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0E+00 * a * c
  q = - 0.5E+00 * ( b + sign ( 1.0E+00, b ) * sqrt ( disc ) )
  r1 = q / a
  r2 = c / q

  return
end
subroutine rpoly2_rroot ( a, b, c, r1, r2 )
!
!*******************************************************************************
!
!! RPOLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
!
!
!  Example:
!
!    A    B    C       roots              R1   R2
!   --   --   --     ------------------   --   --
!    1   -4    3     1          3          1    3
!    1    0    4     2*i      - 2*i        0    0
!    2   -6    5     3 +   i    3 -   i    3    3
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
!    Input, real A, B, C, the coefficients of the quadratic polynomial
!    A * X**2 + B * X + C = 0 whose roots are desired.
!    A must not be zero.
!
!    Output, real R1, R2, the real parts of the two roots of the polynomial.
!
  real a
  real b
  real c
  real disc
  real q
  real r1
  real r2
!
  if ( a == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY2_RROOT - Fatal error!'
    write ( *, * ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0E+00 * a * c
  disc = max ( disc, 0.0E+00 )

  q = ( b + sign ( 1.0E+00, b ) * sqrt ( disc ) )
  r1 = - 0.5E+00 * q / a
  r2 = - 2.0E+00 * c / q

  return
end
subroutine rpoly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )
!
!*******************************************************************************
!
!! RPOLY2_VAL evaluates a parabola defined by three data values.
!
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, three pairs of data values.
!    If the X values are distinct, then all the Y values represent
!    actual values of the parabola.
!
!    Three special cases are allowed:
!
!      X1 == X2 /= X3: Y2 is the derivative at X1;
!      X1 /= X2 == X3: Y3 is the derivative at X3;
!      X1 == X2 == X3: Y2 is the derivative at X1, and
!                      Y3 is the second derivative at X1.
!
!    Input, real X, an abscissa at which the parabola is to be
!    evaluated.
!
!    Output, real Y, YP, YPP, the values of the parabola and
!    its first and second derivatives at X.
!
  integer distinct
  real dif1
  real dif2
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
  real yp
  real ypp
!
!  If any X's are equal, put them and the Y data first.
!
  if ( x1 == x2 .and. x2 == x3 ) then
    distinct = 1
  else if ( x1 == x2 ) then
    distinct = 2
  else if ( x1 == x3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY2_VAL - Fatal error!'
    write ( *, * ) '  X1 = X3 =/= X2.'
    stop
  else if ( x2 == x3 ) then
    distinct = 2
    call r_swap ( x1, x2 )
    call r_swap ( x2, x3 )
    call r_swap ( y1, y2 )
    call r_swap ( y2, y3 )
  else
    distinct = 3
  end if
!
!  Set up the coefficients.
!
  if ( distinct == 1 ) then

    dif1 = y2
    dif2 = 0.5E+00 * y3

  else if ( distinct == 2 ) then

    dif1 = y2
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 ) - y2 ) / ( x3 - x2 )

  else if ( distinct == 3 ) then

    dif1 = ( y2 - y1 ) / ( x2 - x1 )
    dif2 =  ( ( y3 - y1 ) / ( x3 - x1 ) &
            - ( y2 - y1 ) / ( x2 - x1 ) ) / ( x3 - x2 )

  end if
!
!  Evaluate.
!
  y = y1 + ( x - x1 ) * dif1 + ( x - x1 ) * ( x - x2 ) * dif2
  yp = dif1 + ( 2.0E+00 * x - x1 - x2 ) * dif2
  ypp = 2.0E+00 * dif2

  return
end
subroutine rpoly2_val2 ( ndim, ndata, tdata, ydata, left, tval, yval )
!
!*******************************************************************************
!
!! RPOLY2_VAL2 evaluates a parabolic interpolant through tabular data.
!
!
!  Discussion:
!
!    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
!    It constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of a single data point.
!    NDIM must be at least 1.
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.  The
!    values in TDATA must be in strictly ascending order.
!
!    Input, real YDATA(NDIM,NDATA), the data points corresponding to
!    the abscissas.
!
!    Input, integer LEFT, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= LEFT <= NDATA - 2.
!
!    Input, real TVAL, the value of T at which the parabolic interpolant
!    is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA), and
!    the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real YVAL(NDIM), the value of the parabolic interpolant
!    at TVAL.
!
  integer ndata
  integer ndim
!
  real dif1
  real dif2
  integer i
  integer left
  real t1
  real t2
  real t3
  real tval
  real tdata(ndata)
  real ydata(ndim,ndata)
  real y1
  real y2
  real y3
  real yval(ndim)
!
!  Check.
!
  if ( left < 1 .or. left > ndata-2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY2_VAL2 - Fatal error!'
    write ( *, * ) '  LEFT < 1 or LEFT > NDATA-2.'
    stop
  end if

  if ( ndim < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY2_VAL2 - Fatal error!'
    write ( *, * ) '  NDIM < 1.'
    stop
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t1 >= t2 .or. t2 >= t3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY2_VAL2 - Fatal error!'
    write ( *, * ) '  T1 >= T2 or T2 >= T3.'
    stop
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  do i = 1, ndim

    y1 = ydata(i,left)
    y2 = ydata(i,left+1)
    y3 = ydata(i,left+2)

    dif1 = ( y2 - y1 ) / ( t2 - t1 )
    dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
           - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

    yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

  end do

  return
end
subroutine rpoly3_root ( a, b, c, d, r1, r2, r3 )
!
!*******************************************************************************
!
!! RPOLY3_ROOT returns the three roots of a cubic polynomial.
!
!
!  Modified:
!
!    02 March 1999
!
!  Parameters:
!
!    Input, real A, B, C, D, the coefficients of the cubic polynomial
!    A * X**3 + B * X**2 + C * X + D = 0 whose roots are desired.
!    A must not be zero.
!
!    Output, complex R1, R2, R3, the roots of the polynomial, which
!    will include at least one real root.
!
  real a
  real b
  real c
  real d
  complex I
  complex one
  real pi
  real q
  real r
  complex r1
  complex r2
  complex r3
  real s1
  real s2
  real temp
  real theta
!
  if ( a == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY3_ROOT - Fatal error!'
    write ( *, * ) '  A must not be zero!'
    stop
  end if

  one = cmplx ( 1.0E+00, 0.0E+00 )
  I = sqrt ( - one )

  q = ( ( b / a )**2 - 3.0E+00 * ( c / a ) ) / 9.0E+00

  r = ( 2.0E+00 * ( b / a )**3 - 9.0E+00 * ( b / a ) * ( c / a ) &
      + 27.0E+00 * ( d / a ) ) / 54.0E+00

  if ( r * r < q * q * q ) then

    theta = acos ( r / sqrt ( q**3 ) )
    r1 = - 2.0E+00 * sqrt(q) * cos (   theta                      / 3.0E+00 )
    r2 = - 2.0E+00 * sqrt(q) * cos ( ( theta + 2.0E+00 * pi ( ) ) / 3.0E+00 )
    r3 = - 2.0E+00 * sqrt(q) * cos ( ( theta + 4.0E+00 * pi ( ) ) / 3.0E+00 )

  else if ( r * r >= q * q * q ) then

    temp = - r + sqrt ( r**2 - q**3 )
    s1 = sign ( 1.0E+00, temp ) * ( abs ( temp ) )**(1.0/3.0)

    temp = - r - sqrt ( r**2 - q**3 )
    s2 = sign ( 1.0E+00, temp ) * ( abs ( temp ) )**(1.0/3.0)

    r1 = s1 + s2
    r2 = - 0.5E+00 * ( s1 + s2 ) + I * 0.5E+00 * sqrt ( 3.0E+00 ) * ( s1 - s2 )
    r3 = - 0.5E+00 * ( s1 + s2 ) - I * 0.5E+00 * sqrt ( 3.0E+00 ) * ( s1 - s2 )

  end if

  r1 = r1 - b / ( 3.0E+00 * a )
  r2 = r2 - b / ( 3.0E+00 * a )
  r3 = r3 - b / ( 3.0E+00 * a )

  return
end
subroutine rpoly4_root ( a, b, c, d, e, r1, r2, r3, r4 )
!
!*******************************************************************************
!
!! RPOLY4_ROOT returns the four roots of a quartic polynomial.
!
!
!  Modified:
!
!    27 August 1999
!
!  Parameters:
!
!    Input, real A, B, C, D, the coefficients of the polynomial
!    A * X**4 + B * X**3 + C * X**2 + D * X + E = 0 whose roots are
!    desired.  A must not be zero.
!
!    Output, complex R1, R2, R3, R4, the roots of the polynomial.
!
  real a
  real a3
  real a4
  real b
  real b3
  real b4
  real c
  real c3
  real c4
  real d
  real d3
  real d4
  real e
  complex p
  complex q
  complex r
  complex r1
  complex r2
  complex r3
  complex r4
  complex zero
!
  zero = cmplx ( 0.0E+00, 0.0E+00 )

  if ( a == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY4_ROOT - Fatal error!'
    write ( *, * ) '  A must not be zero!'
    stop
  end if

  a4 = b / a
  b4 = c / a
  c4 = d / a
  d4 = e / a
!
!  Set the coefficients of the resolvent cubic equation.
!
  a3 = 1.0E+00
  b3 = - b4
  c3 = a4 * c4 - 4.0E+00 * d4
  d3 = - a4 * a4 * d4 + 4.0E+00 * b4 * d4 - c4 * c4
!
!  Find the roots of the resolvent cubic.
!
  call rpoly3_root ( a3, b3, c3, d3, r1, r2, r3 )
!
!  Choose one root of the cubic, here R1.
!
!  Set R = sqrt ( 0.25E+00 * A4**2 - B4 + R1 )
!
  r = sqrt ( 0.25E+00 * a4**2 - b4  + r1 )

  if ( r /= zero ) then

    p = sqrt ( 0.75E+00 * a4**2 - r**2 - 2.0E+00 * b4 &
        + 0.25E+00 * ( 4.0E+00 * a4 * b4 - 8.0E+00 * c4 - a4**3 ) / r )

    q = sqrt ( 0.75E+00 * a4**2 - r**2 - 2.0E+00 * b4 &
        - 0.25E+00 * ( 4.0E+00 * a4 * b4 - 8.0E+00 * c4 - a4**3 ) / r )

  else

    p = sqrt ( 0.75E+00 * a4**2 - 2.0E+00 * b4 &
      + 2.0E+00 * sqrt ( r1**2 - 4.0E+00 * d4 ) )

    q = sqrt ( 0.75E+00 * a4**2 - 2.0E+00 * b4 &
      - 2.0E+00 * sqrt ( r1**2 - 4.0E+00 * d4 ) )

  end if
!
!  Set the roots.
!
  r1 = - 0.25E+00 * a4 + 0.5E+00 * r + 0.5E+00 * p
  r2 = - 0.25E+00 * a4 + 0.5E+00 * r - 0.5E+00 * p
  r3 = - 0.25E+00 * a4 - 0.5E+00 * r + 0.5E+00 * q
  r4 = - 0.25E+00 * a4 - 0.5E+00 * r - 0.5E+00 * q

  return
end
subroutine rpoly_degree ( na, a, degree )
!
!*******************************************************************************
!
!! RPOLY_DEGREE returns the degree of a polynomial.
!
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
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
!    Input, integer NA, the dimension of A.
!
!    Input, real A(0:NA), the coefficients of the polynomials.
!
!    Output, integer DEGREE, the degree of A.
!
  integer na
!
  real a(0:na)
  integer degree
!
  degree = na

  do while ( degree > 0 )

    if ( a(degree) /= 0.0E+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine rpoly_lagrange_coef ( npol, ipol, xpol, pcof )
!
!*******************************************************************************
!
!! RPOLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.
!
!
!  Definition:
!
!    Given distinct abscissas XPOL(1:NPOL), the IPOL-th Lagrange
!    polynomial P(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real XPOL(NPOL), the abscissas of the Lagrange polynomials.
!    The entries in XPOL must be distinct.
!
!    Output, real PCOF(0:NPOL-1), the polynomial coefficients of the
!    IPOL-th Lagrange polynomial.
!
!  P(IPOL)(X) = SUM ( 0 <= I <= NPOL-1 ) PCOF(I) * X**I
!
  integer npol
!
  integer i
  integer indx
  integer ipol
  integer j
  real pcof(0:npol-1)
  logical rvec_distinct
  real xpol(npol)
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. ipol > npol ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY_LAGRANGE_COEF - Fatal error!'
    write ( *, * ) '  1 <= IPOL <= NPOL is required.'
    stop
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. rvec_distinct ( npol, xpol ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY_LAGRANGE_COEF - Fatal error!'
    write ( *, * ) '  Two or more entries of XPOL are equal:'
    stop
  end if

  pcof(0) = 1.0E+00
  pcof(1:npol-1) = 0.0E+00

  indx = 0

  do i = 1, npol

    if ( i /= ipol ) then

      indx = indx + 1

      do j = indx, 0, -1

        pcof(j) = - xpol(i) * pcof(j) / ( xpol(ipol) - xpol(i) )

        if ( j > 0 ) then
          pcof(j) = pcof(j) + pcof(j-1) / ( xpol(ipol) - xpol(i) )
        end if

      end do

    end if

  end do

  return
end
subroutine rpoly_lagrange_factor ( npol, xpol, xval, wval, dwdx )
!
!*******************************************************************************
!
!! RPOLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.
!
!
!  Formula:
!
!    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!  Discussion:
!
!    Suppose F(X) is at least N times continuously differentiable in the
!    interval [A,B].  Pick NPOL distinct points XPOL(I) in [A,B] and compute
!    the interpolating polynomial P(X) of order NPOL ( and degree NPOL-1)
!    which passes through all the points ( XPOL(I), F(XPOL(I)) ).
!    Then in the interval [A,B], the maximum error
!
!      abs ( F(X) - P(X) )
!
!    is bounded by:
!
!      C * FNMAX * W(X)
!
!    where
!
!      C is a constant,
!      FNMAX is the maximum value of the NPOL-th derivative of F in [A,B],
!      W(X) is the Lagrange factor.
!
!    Thus, the value of W(X) is useful as part of an estimated bound
!    for the interpolation error.
!
!    Note that the Chebyshev abscissas have the property that they minimize
!    the value of W(X) over the interval [A,B].  Hence, if the abscissas may
!    be chosen arbitrarily, the Chebyshev abscissas have this advantage over
!    other choices.
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real XPOL(NPOL), the abscissas, which should be distinct.
!
!    Input, real XVAL, the point at which the Lagrange factor is to be evaluated.
!
!    Output, real WVAL, the value of the Lagrange factor at XVAL.
!
!    Output, real DWDX, the derivative of W with respect to XVAL.
!
  integer npol
!
  real dwdx
  integer i
  integer j
  real term
  real wval
  real xpol(npol)
  real xval
!
  wval = product ( xval - xpol(1:npol) )

  dwdx = 0.0E+00

  do i = 1, npol

    term = 1.0E+00

    do j = 1, npol
      if ( i /= j ) then
        term = term * ( xval - xpol(j) )
      end if
    end do

    dwdx = dwdx + term

  end do

  return
end
subroutine rpoly_lagrange_val ( npol, ipol, xpol, xval, pval, dpdx )
!
!*******************************************************************************
!
!! RPOLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.
!
!
!  Definition:
!
!    Given NPOL distinct abscissas, XPOL(*), the IPOL-th Lagrange
!    polynomial P(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real XPOL(NPOL), the abscissas of the Lagrange polynomials.
!    The entries in XPOL must be distinct.
!
!    Input, real XVAL, the point at which the IPOL-th Lagrange polynomial
!    is to be evaluated.
!
!    Output, real PVAL, the value of the IPOL-th Lagrange polynomial at XVAL.
!
!    Output, real DPDX, the derivative of the IPOL-th Lagrange polynomial at XVAL.
!
  integer npol
!
  real dpdx
  integer i
  integer ipol
  integer j
  real p2
  real pval
  logical rvec_distinct
  real xpol(npol)
  real xval
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. ipol > npol ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY_LAGRANGE_VAL - Fatal error!'
    write ( *, * ) '  1 <= IPOL <= NPOL is required.'
    stop
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. rvec_distinct ( npol, xpol ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY_LAGRANGE_VAL - Fatal error!'
    write ( *, * ) '  Two or more entries of XPOL are equal:'
    stop
  end if
!
!  Evaluate the polynomial.
!
  pval = 1.0E+00

  do i = 1, npol

    if ( i /= ipol ) then

      pval = pval * ( xval - xpol(i) ) / ( xpol(ipol) - xpol(i) )

    end if

  end do
!
!  Evaluate the derivative, which can be found by summing up the result
!  of differentiating one factor at a time, successively.
!
  dpdx = 0.0E+00

  do i = 1, npol

    if ( i /= ipol ) then

      p2 = 1.0E+00
      do j = 1, npol

        if ( j == i ) then
          p2 = p2                      / ( xpol(ipol) - xpol(j) )
        else if ( j /= ipol ) then
          p2 = p2 * ( xval - xpol(j) ) / ( xpol(ipol) - xpol(j) )
        end if

      end do

      dpdx = dpdx + p2

    end if

  end do

  return
end
subroutine rpoly_ls_set ( b, c, d, f, npoint, nterms, w, x )
!
!*******************************************************************************
!
!! RPOLY_LS_SET defines a least squares polynomial for given data.
!
!
!  Discussion:
!
!    The polynomial may be evaluated at any point X by calling RPOLY_LS_VAL.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real B(NTERMS), C(NTERMS), D(NTERMS), are quantities
!    defining the least squares polynomial for the input data,
!    which will be needed to evaluate the polynomial.
!
!    Input, real F(NPOINT), the data values at the points X(*).
!
!    Input, integer NTERMS, the number of terms to use in the
!    approximating polynomial.  NTERMS must be at least 1.
!    The degree of the polynomial is NTERMS-1.
!
!    Input, real W(NPOINT), the weights associated with the data points.
!    Each entry of W should be positive.
!
!    Input, real X(NPOINT), the abscissas of the data points.
!    At least NTERMS-1 of the values in X must be distinct.
!
  integer npoint
  integer nterms
!
  real b(nterms)
  real c(nterms)
  real d(nterms)
  real f(npoint)
  integer i
  integer j
  integer nuniq
  real p
  real pj(npoint)
  real pjm1(npoint)
  real s(nterms)
  real w(npoint)
  real x(npoint)
!
!  Make sure at least NTERMS-1 X values are unique.
!
  call rvec_uniq_count ( npoint, x, nuniq )

  if ( nuniq < nterms-1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RPOLY_LS_SET - Fatal error!'
    write ( *, * ) '  The number of distinct X values must be'
    write ( *, * ) '  at least NTERMS-1 = ', nterms - 1
    write ( *, * ) '  but the input data has only ', nuniq
    write ( *, * ) '  distinct entries.'
    return
  end if
!
!  Make sure all W entries are positive.
!
  do i = 1, npoint
    if ( w(i) <= 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RPOLY_LS_SET - Fatal error!'
      write ( *, * ) '  All weights W must be positive,'
      write ( *, * ) '  but weight ', i, ' is ', w(i)
      return
    end if
  end do
!
!  Start inner product summations at zero.
!
  b(1:nterms) = 0.0E+00
  c(1:nterms) = 0.0E+00
  d(1:nterms) = 0.0E+00
  s(1:nterms) = 0.0E+00
!
!  Set the values of P(-1,X) and P(0,X) at all data points.
!
  pjm1(1:npoint) = 0.0E+00
  pj(1:npoint) = 1.0E+00
!
!  Now compute the value of P(J,X(I)) as
!
!    P(J,X(I)) = ( X(I) - B(J) ) * P(J-1,X(I)) - C(J) * P(J-2,X(I))
!
!  where
!
!    S(J) = < P(J,X), P(J,X) >
!    B(J) = < x*P(J,X), P(J,X) > / < P(J,X), P(J,X) >
!    C(J) = S(J) / S(J-1)
!
!  and the least squares coefficients are
!
!    D(J) = < F(X), P(J,X) > / < P(J,X), P(J,X) >
!
  do j = 1, nterms

    d(j) = d(j) + sum ( w(1:npoint) * f(1:npoint) * pj(1:npoint) )
    b(j) = b(j) + sum ( w(1:npoint) * x(1:npoint) * pj(1:npoint)**2 )
    s(j) = s(j) + sum ( w(1:npoint) * pj(1:npoint)**2 )

    d(j) = d(j) / s(j)

    if ( j == nterms ) then
      c(j) = 0.0E+00
      return
    end if

    b(j) = b(j) / s(j)

    if ( j == 1 ) then
      c(j) = 0.0E+00
    else
      c(j) = s(j) / s(j-1)
    end if

    do i = 1, npoint
      p = pj(i)
      pj(i) = ( x(i) - b(j) ) * pj(i) - c(j) * pjm1(i)
      pjm1(i) = p
    end do

  end do

  return
end
subroutine rpoly_ls_val ( b, c, d, nterms, x, px )
!
!*******************************************************************************
!
!! RPOLY_LS_VAL evaluates a least squares polynomial defined by RPOLY_LS_SET.
!
!
!  Discussion:
!
!    The least squares polynomial is assumed to be defined as a sum
!
!      P(X) = SUM ( I = 1 to NTERMS ) D(I) * P(I-1,X)
!
!    where the orthogonal basis polynomials P(I,X) satisfy the following
!    three term recurrence:
!
!      P(-1,X) = 0
!      P(0,X) = 1
!      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
!
!    Therefore, the least squares polynomial can be evaluated as follows:
!
!    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
!
!    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
!    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
!    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
!    can be eliminated from the sum, and its coefficient merged in with
!    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
!    and so on until a single term remains.
!    P(NTERMS,X) of P(NTERMS-1,X)
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real B(NTERMS), C(NTERMS), D(NTERMS), the information
!    computed by RPOLY_LS_SET.
!
!    Input, integer NTERMS, the number of terms in the least squares
!    polynomial.  NTERMS must be at least 1.  The input value of NTERMS
!    may be reduced from the value given to RPOLY_LS_SET.  This will
!    evaluate the least squares polynomial of the lower degree specified.
!
!    Input, real X, the point at which the least squares polynomial
!    is to be evaluated.
!
!    Output, real PX, the value of the least squares polynomial at X.
!
  integer nterms
!
  real b(nterms)
  real c(nterms)
  real d(nterms)
  integer i
  real prev
  real prev2
  real px
  real x
!
  px = d(nterms)
  prev = 0.0E+00

  do i = nterms-1, 1, -1

    prev2 = prev
    prev = px

    if ( i == nterms-1 ) then
      px = d(i) + ( x - b(i) ) * prev
    else
      px = d(i) + ( x - b(i) ) * prev - c(i+1) * prev2
    end if

  end do

  return
end
subroutine rpoly_ls_val2 ( b, c, d, nterms, x, px, pxp )
!
!*******************************************************************************
!
!! RPOLY_LS_VAL2 evaluates a least squares polynomial defined by RPOLY_LS_SET.
!
!
!  Discussion:
!
!    RPOLY_LS_VAL2 also computes the derivative of the polynomial.
!
!    The least squares polynomial is assumed to be defined as a sum
!
!      P(X) = SUM ( I = 1 to NTERMS ) D(I) * P(I-1,X)
!
!    where the orthogonal basis polynomials P(I,X) satisfy the following
!    three term recurrence:
!
!      P(-1,X) = 0
!      P(0,X) = 1
!      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
!
!    Therefore, the least squares polynomial can be evaluated as follows:
!
!    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
!
!    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
!    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
!    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
!    can be eliminated from the sum, and its coefficient merged in with
!    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
!    and so on until a single term remains.
!    P(NTERMS,X) of P(NTERMS-1,X)
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real B(NTERMS), C(NTERMS), D(NTERMS), the information
!    computed by RPOLY_LS_SET.
!
!    Input, integer NTERMS, the number of terms in the least squares
!    polynomial.  NTERMS must be at least 1.  The value of NTERMS
!    may be reduced from the value given to RPOLY_LS_SET.
!    This will cause RPOLY_LS_VAL to evaluate the least squares polynomial
!    of the lower degree specified.
!
!    Input, real X, the point at which the least squares polynomial
!    is to be evaluated.
!
!    Output, real PX, PXP, the value and derivative of the least
!    squares polynomial at X.
!
  integer nterms
!
  real b(nterms)
  real c(nterms)
  real d(nterms)
  integer i
  real px
  real pxm1
  real pxm2
  real pxp
  real pxpm1
  real pxpm2
  real x
!
  px = d(nterms)
  pxp = 0.0E+00
  pxm1 = 0.0E+00
  pxpm1 = 0.0E+00

  do i = nterms-1, 1, -1

    pxm2 = pxm1
    pxpm2 = pxpm1
    pxm1 = px
    pxpm1 = pxp

    if ( i == nterms-1 ) then
      px = d(i) + ( x - b(i) ) * pxm1
      pxp = pxm1 + ( x - b(i) ) * pxpm1
    else
      px = d(i) + ( x - b(i) ) * pxm1 - c(i+1) * pxm2
      pxp = pxm1 + ( x - b(i) ) * pxpm1 - c(i+1) * pxpm2
    end if

  end do

  return
end
subroutine rpoly_print ( n, a, title )
!
!*******************************************************************************
!
!! RPOLY_PRINT prints out a polynomial.
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
!    Input, integer N, the dimension of A.
!
!    Input, real A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, character ( len = * ), an optional title.
!
  integer n
!
  real a(0:n)
  integer i
  real mag
  integer n2
  character plus_minus
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '

  call rpoly_degree ( n, a, n2 )

  if ( a(n2) < 0.0E+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( n2 >= 2 ) then
    write ( *, '( '' p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( '' p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( '' p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0E+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0E+00 ) then

      if ( i >= 2 ) then
        write ( *, ' ( ''        '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''        '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''        '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

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
!    Input, integer N, the dimension of C.
!
!    Input, real C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X**I.
!
!    Input, real X, the point at which the polynomial is to be evaluated.
!
!    Output, real CX, the value of the polynomial at X.
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
subroutine rrow_max ( lda, m, n, a, imax, amax )
!
!*******************************************************************************
!
!! RROW_MAX returns the maximums of rows of a real array.
!
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, integer IMAX(M); IMAX(I) is the column of A in which
!    the maximum for row I occurs.
!
!    Output, real AMAX(M), the maximums of the rows.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real amax(m)
  integer i
  integer imax(m)
  integer j
!
  do i = 1, m

    imax(i) = 1
    amax(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) > amax(i) ) then
        imax(i) = j
        amax(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine rrow_mean ( lda, m, n, a, mean )
!
!*******************************************************************************
!
!! RROW_MEAN returns the means of rows of a real array.
!
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      2
!      5
!
!  Modified:
!
!    24 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, real MEAN(M), the means, or averages, of the rows.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer j
  real mean(m)
!
  do i = 1, m
    mean(i) = sum ( a(i,1:n) ) / real ( n )
  end do

  return
end
subroutine rrow_min ( lda, m, n, a, imin, amin )
!
!*******************************************************************************
!
!! RROW_MIN returns the minimums of rows of a real array.
!
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real A(LDA,N), the array to be examined.
!
!    Output, integer IMIN(M); IMIN(I) is the column of A in which
!    the minimum for row I occurs.
!
!    Output, real AMIN(M), the minimums of the rows.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real amin(m)
  integer i
  integer imin(m)
  integer j
!
  do i = 1, m

    imin(i) = 1
    amin(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) < amin(i) ) then
        imin(i) = j
        amin(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine rrow_sum ( lda, m, n, a, rowsum )
!
!*******************************************************************************
!
!! RROW_SUM returns the sums of the rows of a table.
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
!    Input, real A(LDA,N), the M by N array.
!
!    Output, real ROWSUM(M), the sum of the entries of each row of A.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  real rowsum(m)
!
  do i = 1, m
    rowsum(i) = sum ( a(i,1:n) )
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
subroutine rrow_to_rvec ( lda, m, n, a, x )
!
!*******************************************************************************
!
!! RROW_TO_RVEC converts a matrix of rows into a vector.
!
!
!  Example:
!
!    M = 3, N = 4
!
!    A =
!      11 12 13 14
!      21 22 23 24
!      31 32 33 34
!
!    X = ( 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34 )
!
!  Modified:
!
!    13 July 2000
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
!    Input, real A(LDA,N), the M by N array.
!
!    Output, real X(M*N), a vector containing the M rows of A.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer m
  real x(m*n)
!
  j = 1
  do i = 1, m
    x(j:j+n-1) = a(i,1:n)
    j = j + n
  end do

  return
end
subroutine rrow_variance ( lda, m, n, a, variance )
!
!*******************************************************************************
!
!! RROW_VARIANCE returns the variances of the rows of a real array.
!
!
!  Modified:
!
!    24 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real A(LDA,N), the array whose variances are desired.
!
!    Output, real VARIANCE(M), the variances of the rows.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  integer i
  integer j
  real mean
  real variance(m)
!
  do i = 1, m

    mean = sum ( a(i,1:n) ) / real ( n )

    variance(i) = 0.0E+00
    do j = 1, n
      variance(i) = variance(i) + ( a(i,j) - mean )**2
    end do

    if ( n > 1 ) then
      variance(i) = variance(i) / real ( n - 1 )
    else
      variance(i) = 0.0E+00
    end if

  end do

  return
end
subroutine rvec2_compare ( n, a1, a2, i, j, isgn )
!
!*******************************************************************************
!
!! RVEC2_COMPARE compares pairs of reals stored in two vectors.
!
!
!  Discussion:
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    A1(I) A2(I)   A1(J) A2(J)  ISGN
!    -----------   -----------  ----
!    1.0   5.0  <  1.0   6.0     -1
!    1.0   5.0  <  2.0   8.0     -1
!    1.0   5.0  <  9.0   1.0     -1
!    1.0   5.0  =  1.0   5.0      0
!    1.0   5.0  >  0.0   2.0     +1
!    1.0   5.0  >  0.0   5.0     +1
!    1.0   5.0  >  1.0   3.0     +1
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data items.
!
!    Input, real A1(N), A2(N), contain the two components of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  integer n
!
  real a1(n)
  real a2(n)
  integer i
  integer isgn
  integer j
!
  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(i) > a2(j) ) then
      isgn = +1
    end if

  else if ( a1(i) > a1(j) ) then

    isgn = +1

  end if

  return
end
subroutine rvec2_print ( n, a1, a2, title )
!
!*******************************************************************************
!
!! RVEC2_PRINT prints a pair of real vectors.
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
subroutine rvec2_print_some ( n, x1, x2, max_print )
!
!*******************************************************************************
!
!! RVEC2_PRINT_SOME prints "some" of two real vectors.
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
subroutine rvec2_sort_a ( n, a1, a2 )
!
!*******************************************************************************
!
!! RVEC2_SORT_A ascending sorts a vector of real (X,Y) data.
!
!
!  Discussion:
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, real A1(N), A2(N), the data to be sorted.
!
  integer n
!
  real a1(n)
  real a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call r_swap ( a1(i), a1(j) )
      call r_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call rvec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine rvec2_sort_d ( n, a1, a2 )
!
!*******************************************************************************
!
!! RVEC2_SORT_A descending sorts a vector of real (X,Y) data.
!
!
!  Discussion:
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, real A1(N), A2(N), the data to be sorted.
!
  integer n
!
  real a1(n)
  real a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call r_swap ( a1(i), a1(j) )
      call r_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!  Reverse the value of ISGN to effect a descending sort.
!
    else if ( indx < 0 ) then

      call rvec2_compare ( n, a1, a2, i, j, isgn )

      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine rvec2_sort_heap_index_a ( n, x, y, indx )
!
!*******************************************************************************
!
!! RVEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of (X,Y) data.
!
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    ( X(I), Y(I) ) < ( X(J), Y(J) ) if:
!
!    * X(I) < X(J), or
!
!    * X(I) = X(J), and Y(I) < Y(J).
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      ( X(INDX(I)), Y(INDX(I) ), for I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call RVEC_PERMUTE ( N, X, INDX )
!      call RVEC_PERMUTE ( N, Y, INDX )
!
!    after which ( X(I), Y(I) ), I = 1 to N is sorted.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real X(N),Y(N), pairs of X, Y coordinates of points.
!
!    Output, integer INDX(N), contains the sort index.  The
!    I-th element of the sorted array has coordinates ( X(INDX(I)), Y(INDX(I) ).
!
  integer n
!
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l
  real x(n)
  real xval
  real y(n)
  real yval
!
  call ivec_identity ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( l > 1 ) then

      l = l - 1
      indxt = indx(l)
      xval = x(indxt)
      yval = y(indxt)

    else

      indxt = indx(ir)
      xval = x(indxt)
      yval = y(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then

        if ( x(indx(j)) < x(indx(j+1)) .or. &
          ( x(indx(j)) == x(indx(j+1)) .and. y(indx(j)) < y(indx(j+1)) ) ) then
          j = j + 1
        end if

      end if

      if ( xval < x(indx(j)) .or. &
          ( xval == x(indx(j)) .and. yval < y(indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine rvec2_sum_imax ( n, a, b, imax )
!
!*******************************************************************************
!
!! RVEC2_SUM_IMAX returns the index of the maximum sum of two real vectors.
!
!
!  Modified:
!
!    16 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), B(N), two arrays whose sum is to be examined.
!
!    Output, integer IMAX, the index of the largest entry in A+B.
!
  integer n
!
  real a(n)
  real b(n)
  integer i
  integer imax
  real sum_max
!
  if ( n <= 0 ) then

    imax = 0

  else

    imax = 1
    sum_max = a(1) + b(1)

    do i = 2, n
      if ( a(i) + b(i) > sum_max ) then
        sum_max = a(i) + b(i)
        imax = i
      end if
    end do

  end if

  return
end
subroutine rvec2_uniq ( n, a1, a2, nuniq )
!
!*******************************************************************************
!
!! RVEC2_UNIQ keeps the unique elements in a array of pairs of reals.
!
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items.
!
!    Input/output, real A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer NUNIQ, the number of unique items.
!
  integer n
!
  real a1(n)
  real a2(n)
  integer itest
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a1(itest) /= a1(nuniq) .or. a2(itest) /= a2(nuniq) ) then

      nuniq = nuniq + 1

      a1(nuniq) = a1(itest)
      a2(nuniq) = a2(itest)

    end if

  end do

  return
end
function rvec3_compare ( x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! RVEC3_COMPARE compares two R3 vectors.
!
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, the first vector.
!
!    Input, real X2, Y2, Z2, the second vector.
!
!    Output, character RVEC3_COMPARE: '<', '>' or '=' if the first vector
!    is less, greater or equal to the second.
!
  character c
  character rvec3_compare
  real x1
  real x2
  real y1
  real y2
  real z1
  real z2
!
  if ( x1 < x2 ) then
    c = '<'
  else if ( x1 > x2 ) then
    c = '>'
  else if ( y1 < y2 ) then
    c = '<'
  else if ( y1 > y2 ) then
    c = '>'
  else if ( z1 < z2 ) then
    c = '<'
  else if ( z1 > z2 ) then
    c = '>'
  else
    c = '='
  end if

  rvec3_compare = c

  return
end
subroutine rvec3_index_insert_unique ( maxn, n, x, y, z, indx, &
  xval, yval, zval, ival, ierror )
!
!*******************************************************************************
!
!! RVEC3_INDEX_INSERT_UNIQUE inserts a unique R3 value in an indexed sorted list.
!
!
!  Discussion:
!
!    If the input value does not occur in the current list, it is added,
!    and N, X, Y, Z and INDX are updated.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input/output, integer N, the size of the list.
!
!    Input/output, real X(N), Y(N), Z(N), the list of R3 vectors.
!
!    Input/output, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, YVAL, ZVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL, ZVAL.
!
!    Output, integer IERROR, 0 for no error, 1 if an error occurred.
!
  integer maxn
!
  integer equal
  integer ierror
  integer indx(maxn)
  integer ival
  integer less
  integer more
  integer n
  real x(maxn)
  real xval
  real y(maxn)
  real yval
  real z(maxn)
  real zval
!
  ierror = 0

  if ( n <= 0 ) then

    if ( maxn <= 0 ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'RVEC3_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, * ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = xval
    y(1) = yval
    z(1) = zval
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
!
  call rvec3_index_search ( maxn, n, x, y, z, indx, xval, yval, zval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( n >= maxn ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'RVEC3_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, * ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = xval
    y(n+1) = yval
    z(n+1) = zval
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine rvec3_index_search ( maxn, n, x, y, z, indx, xval, yval, &
  zval, less, equal, more )
!
!*******************************************************************************
!
!! RVEC3_INDEX_SEARCH searches for an R3 value in an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input, integer N, the size of the current list.
!
!    Input, real X(N), Y(N), Z(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, YVAL, ZVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  integer maxn
!
  character c
  integer equal
  integer hi
  integer indx(maxn)
  integer less
  integer lo
  integer mid
  integer more
  integer n
  character rvec3_compare
  real x(maxn)
  real xhi
  real xlo
  real xmid
  real xval
  real y(maxn)
  real yhi
  real ylo
  real ymid
  real yval
  real z(maxn)
  real zhi
  real zlo
  real zmid
  real zval
!
  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))
  zlo = z(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))
  zhi = z(indx(hi))

  c = rvec3_compare ( xval, yval, zval, xlo, ylo, zlo )

  if ( c == '<' ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( c == '=' ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  c = rvec3_compare ( xval, yval, zval, xhi, yhi, zhi )

  if ( c == '>' ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( c == '=' ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))
    zmid = z(indx(mid))

    c = rvec3_compare ( xval, yval, zval, xmid, ymid, zmid )

    if ( c == '=' ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( c == '<' ) then
      hi = mid
    else if ( c == '>' ) then
      lo = mid
    end if

  end do

  return
end
subroutine rvec3_print ( n, a1, a2, a3, title )
!
!*******************************************************************************
!
!! RVEC3_PRINT prints a trio of real vectors.
!
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A1(N), A2(N), A3(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,3g14.6)' ) i, a1(i), a2(i), a3(i)
  end do

  return
end
function rvec3_striple ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! RVEC3_STRIPLE computes the scalar triple product of three vectors.
!
!
!  Definition:
!
!    STRIPLE = V1 dot (V2 x V3).
!
!    STRIPLE is the volume of the parallelogram whose sides are
!    formed by V1, V2 and V3.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the three vectors.
!
!    Output, real RVEC3_STRIPLE, the scalar triple product.
!
  real rvec3_striple
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  rvec3_striple =  x1 * ( y2 * z3 - z2 * y3 ) &
                 + y1 * ( z2 * x3 - x2 * z3 ) &
                 + z1 * ( x2 * y3 - y2 * x3 )

  return
end
subroutine rvec3_vtriple ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z )
!
!*******************************************************************************
!
!! RVEC3_VTRIPLE computes the vector triple product of three vectors.
!
!
!  Definition:
!
!    VTRIPLE = V1 x (V2 x V3)
!
!    VTRIPLE is a vector perpendicular to V1, lying in the plane
!    spanned by V2 and V3.  The norm of VTRIPLE is the product
!    of the norms of V1, V2 and V3.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the three vectors.
!
!    Output, real X, Y, Z, the components of the vector triple product.
!
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
  real z
  real z1
  real z2
  real z3
  real z4
!
  call rvec_cross_3d ( x2, y2, z2, x3, y3, z3, x4, y4, z4 )

  call rvec_cross_3d ( x1, y1, z1, x4, y4, z4, x, y, z )

  return
end
subroutine rvec4_unit_euclidean_4d ( w, x, y, z )
!
!*******************************************************************************
!
!! RVEC4_UNIT_EUCLIDEAN_4D Euclidean normalizes a vector in 4D.
!
!
!  Modified:
!
!    17 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, W, X, Y, Z, the components of the vector.
!
  real norm
  real w
  real x
  real y
  real z
!
  norm = sqrt ( w * w + x * x + y * y + z * z )

  if ( norm /= 0.0E+00 ) then
    w = w / norm
    x = x / norm
    y = y / norm
    z = z / norm
  end if

  return
end
subroutine rvec_01_to_ab ( n, a, amax, amin )
!
!*******************************************************************************
!
!! RVEC_01_TO_AB shifts and rescales data to lie within given bounds.
!
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input/output, real A(N).
!
!    On input, A contains the original data, which is presumed to lie
!    between 0 and 1.  However, it is not necessary that this be so.
!
!    On output, A has been shifted and rescaled so that all entries which
!    on input lay in [0,1] now lie between AMIN and AMAX.  Other entries will
!    be mapped in a corresponding way.
!
!    Input, real AMAX, AMIN, the maximum and minimum values allowed for A.
!
  integer n
!
  real a(n)
  real amax
  real amax2
  real amax3
  real amin
  real amin2
  real amin3
  integer i
!
  if ( amax == amin ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_01_TO_AB - Fatal error!'
    write ( *, * ) '  AMAX = AMIN'
    write ( *, * ) '  AMAX = ', amax
    write ( *, * ) '  AMIN = ', amin
    stop
  end if

  amax2 = max ( amax, amin )
  amin2 = min ( amax, amin )

  call rvec_min ( n, a, amin3 )
  call rvec_max ( n, a, amax3 )

  if ( amax3 /= amin3 ) then

    do i = 1, n
      a(i) = ( ( amax3 - a(i) ) * amin2 + ( a(i) - amin3 ) * amax2 ) &
        / ( amax3 - amin3 )
    end do

  else

    a(1:n) = 0.5E+00 * ( amax2 + amin2 )

  end if

  return
end
subroutine rvec_ab_to_01 ( n, a )
!
!*******************************************************************************
!
!! RVEC_AB_TO_01 shifts and rescales data to lie within [0,1].
!
!
!  Formula:
!
!    A(I) := ( A(I) - AMIN ) / ( AMAX - AMIN )
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input/output, real A(N).
!
!    On input, A contains the original data.  On output, A has been shifted
!    and scaled so that all entries lie between 0 and 1.
!
  integer n
!
  real a(n)
  real amax
  real amin
!
  call rvec_max ( n, a, amax )
  call rvec_min ( n, a, amin )

  if ( amin == amax ) then
    a(1:n) = 0.5E+00
  else
    a(1:n) = ( a(1:n) - amin ) / ( amax - amin )
  end if

  return
end
subroutine rvec_ab_to_cd ( n, a, amax, amin )
!
!*******************************************************************************
!
!! RVEC_AB_TO_CD shifts and rescales data to lie within a given pair of bounds.
!
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input/output, real A(N).
!
!    On input, A contains the original data.  On output,
!    A has been shifted and rescaled so that all entries lie
!    between AMIN and AMAX.
!
!    Input, real AMAX, AMIN, the maximum and minimum values allowed for A.
!
  integer n
!
  real a(n)
  real amax
  real amax2
  real amax3
  real amin
  real amin2
  real amin3
  integer i
!
  if ( amax == amin ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_AB_TO_CD - Fatal error!'
    write ( *, * ) '  AMAX = AMIN:'
    write ( *, * ) '  AMAX = ', amax
    write ( *, * ) '  AMIN = ', amin
    stop
  end if

  amax2 = max ( amax, amin )
  amin2 = min ( amax, amin )

  call rvec_min ( n, a, amin3 )
  call rvec_max ( n, a, amax3 )

  if ( amax3 /= amin3 ) then

    do i = 1, n
      a(i) = ( ( amax3 - a(i) ) * amin2 + ( a(i) - amin3 ) * amax2 &
        ) / ( amax3 - amin3 )
    end do

  else

    a(1:n) = 0.5E+00 * ( amax2 + amin2 )

  end if

  return
end
subroutine rvec_amax ( n, a, amax )
!
!*******************************************************************************
!
!! RVEC_AMAX returns the maximum absolute value in a real vector.
!
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), the array.
!
!    Output, real AMAX, the value of the entry of largest magnitude.
!
  integer n
!
  real a(n)
  real amax
  integer i
!
  amax = maxval ( abs ( a(1:n) ) )

  return
end
subroutine rvec_amin ( n, a, amin )
!
!*******************************************************************************
!
!! RVEC_AMIN returns the minimum absolute value in a real vector.
!
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), the array.
!
!    Output, real AMIN, the value of the entry of smallest magnitude.
!
  integer n
!
  real a(n)
  real amin
  integer i
!
  amin = minval ( abs ( a(1:n) ) )

  return
end
subroutine rvec_bin ( n, a, nbin, bin_min, bin_max, bin, bin_limit )
!
!*******************************************************************************
!
!! RVEC_BIN bins a real vector, returning the population of each bin.
!
!
!  Discussion:
!
!    The user specifies minimum and maximum bin values, BIN_MIN and
!    BIN_MAX, and the number of bins, NBIN.  This determines a
!    "bin width":
!
!      H = ( BIN_MAX - BIN_MIN ) / NBIN
!
!    so that bin I will count all entries X(J) such that
!
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!
!    The array does NOT have to be sorted.
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of A.
!
!    Input, real A(N), an (unsorted) array to be binned.
!
!    Input, integer NBIN, the number of bins.  Two extra bins, #0 and
!    #NBIN+1, count extreme values.
!
!    Input, real BIN_MIN, BIN_MAX, define the range and size of the bins.
!    BIN_MIN and BIN_MAX must be distinct.
!    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
!    this, but proper results will be computed if BIN_MIN > BIN_MAX.
!
!    Output, integer BIN(0:NBIN+1).
!    BIN(0) counts entries of A less than BIN_MIN.
!    BIN(NBIN+1) counts entries greater than or equal to BIN_MAX.
!    For 1 <= I <= NBIN, BIN(I) counts the entries X(J) such that
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!    where H is the bin spacing.
!
!    Output, real BIN_LIMIT(0:NBIN), the "limits" of the bins.
!    BIN(I) counts the number of entries X(J) such that
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!
  integer n
  integer nbin
!
  real a(n)
  integer bin(0:nbin+1)
  real bin_limit(0:nbin)
  real bin_max
  real bin_min
  integer i
  integer j
  real t
!
  if ( bin_max == bin_min ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_BIN - Fatal error!'
    write ( *, * ) '  BIN_MIN = BIN_MAX.'
    stop
  end if

  bin(0:nbin+1) = 0

  do i = 1, n

    t = ( a(i) - bin_min ) / ( bin_max - bin_min )

    if ( t < 0.0E+00 ) then
      j = 0
    else if ( t >= 1.0E+00 ) then
      j = nbin + 1
    else
      j = 1 + int ( real ( nbin ) * t )
    end if

    bin(j) = bin(j) + 1

  end do
!
!  Compute the bin limits.
!
  do i = 0, nbin
    bin_limit(i) = ( real ( nbin - i ) * bin_min &
      + real ( i ) * bin_max ) / real ( nbin )
  end do

  return
end
subroutine rvec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! RVEC_BIN_EVEN bins a real array into evenly spaced bins.
!
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!  Modified:
!
!    08 February 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the set.
!
!    Input, real A(N), the data to be binned.
!
!    Input, integer NBIN, the number of bins.
!
!    Input, real BIN_MIN, BIN_MAX, the bin limits.
!
!    Output, BIN_START(NBIN), BIN_LAST(NBIN), the index of the first and
!    last element of A that went into each bin, or -1 if there are no
!    entries in this bin.
!
!    Output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer n
  integer nbin
!
  real a(n)
  integer bin_last(nbin)
  integer bin_next(n)
  integer bin_start(nbin)
  real bin_max
  real bin_min
  integer i
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin) = -1
  bin_start(1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do i = 1, n

    call r_to_bin_even ( nbin, bin_min, bin_max, a(i), j )

    if ( bin_start(j) == -1 ) then
      bin_start(j) = i
    else
      k = bin_last(j)
      bin_next(k) = i
    end if

    bin_next(i) = 0

    bin_last(j) = i

  end do

  return
end
subroutine rvec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! RVEC_BINNED_REORDER reorders a real binned data vector.
!
!
!  Discussion:
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START, BIN_LAST and BIN_NEXT arrays have also been
!    updated so that they still correspond to the (rearranged) vector A.
!
!  Modified:
!
!    08 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(N), the data to be sorted.
!
!    Input, integer NBIN, the number of bins.
!
!    Input/output, BIN_START(NBIN), BIN_LAST(NBIN), contains the index of
!    the first and last element of A that went into each bin, or -1 if
!    there are no entries in the bin.
!
!    Input/output, BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  integer n
  integer nbin
!
  real a(n)
  real a2(n)
  integer bin_last(nbin)
  integer bin_next(n)
  integer bin_start(nbin)
  real bin_max
  real bin_min
  integer i
  integer j
  integer k
!
  k = 0

  do i = 1, nbin

    j = bin_start(i)

    if ( j > 0 ) then
      bin_start(i) = k + 1
    end if

    do while ( j > 0 )
      k = k + 1
      bin_last(i) = k
      a2(k) = a(j)
      j = bin_next(j)
    end do

  end do

  a(1:n) = a2(1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i = 1, nbin
    k = bin_last(i)
    if ( k > 0 ) then
      bin_next(k) = 0
    end if
  end do

  return
end
subroutine rvec_binned_sort_a ( n, a, nbin, bin_start, bin_last )
!
!*******************************************************************************
!
!! RVEC_BINNED_SORT_A ascending sorts a real binned reordered data vector.
!
!
!  Discussion:
!
!    Presumably, the data vector was first binned by RVEC_BIN_EVEN,
!    then reordered by RVEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time, which results in sorting
!    the entire vector.
!
!  Modified:
!
!    08 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(N), the data to be sorted.
!
!    Input, integer NBIN, the number of bins.
!
!    Input, BIN_START(NBIN), BIN_LAST(NBIN), contains the index of the first
!    and last element of A that went into each bin, or -1 if there were no
!    elements in the bin.
!
  integer n
  integer nbin
!
  real a(n)
  integer bin_last(nbin)
  integer bin_start(nbin)
  integer i
  integer i1
  integer i2
  integer n1
!
  do i = 1, nbin

    i1 = bin_start(i)

    if ( i1 > 0 ) then

      i2 = bin_last(i)

      n1 = i2 + 1 - i1

      call rvec_sort_quick_a ( n1, a(i1:i2) )

    end if

  end do

  return
end
subroutine rvec_blend ( n, t1, x1, t2, x2, x )
!
!*******************************************************************************
!
!! RVEC_BLEND interpolates a vector, given two vectors and weight factors.
!
!
!  Formula:
!
!    x(i) = t * x1(i) + (1-t) * x2(i)
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in each  vector.
!
!    Input, real T1, the weight factor for vector 1.
!
!    Input, real X1(N), the first vector.
!
!    Input, real T2, the weight factor for vector 2.
!
!    Input, real X2(N), the second vector.
!
!    Output, real X(N), the interpolated or extrapolated value.
!
  integer n
!
  real t1
  real t2
  real x(n)
  real x1(n)
  real x2(n)
!
  x(1:n) = t1 * x1(1:n) + t2 * x2(1:n)

  return
end
subroutine rvec_bracket ( n, x, xval, left, right )
!
!*******************************************************************************
!
!! RVEC_BRACKET searches a sorted array for successive brackets of a value.
!
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input, real X(N), an array that has been sorted into ascending order.
!
!    Input, real XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      XVAL > X(N), when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  integer n
!
  integer i
  integer left
  integer right
  real x(n)
  real xval
!
  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
subroutine rvec_bracket2 ( n, x, xval, start, left, right )
!
!*******************************************************************************
!
!! RVEC_BRACKET2 searches a sorted array for successive brackets of a value.
!
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    containing the given value.
!
!    RVEC_BRACKET2 is a variation on RVEC_BRACKET.  It seeks to reduce
!    the search time by allowing the user to suggest an interval that
!    probably contains the value.  The routine will look in that interval
!    and the intervals to the immediate left and right.  If this does
!    not locate the point, a binary search will be carried out on
!    appropriate subportion of the sorted array.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of the input array.
!
!    Input, real X(N), an array that has been sorted into ascending order.
!
!    Input, real XVAL, a value to be bracketed by entries of X.
!
!    Input, integer START, between 1 and N, specifies that XVAL
!    is likely to be in the interval:
!
!      [ X(START), X(START+1) ]
!
!    or, if not in that interval, then either
!
!      [ X(START+1), X(START+2) ]
!    or
!      [ X(START-1), X(START) ].
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
!    and X(LEFT) <= XVAL <= X(RIGHT).
!
!    Special cases:
!      Value is less than all data values:
!    LEFT = -1, RIGHT = 1, and XVAL < X(RIGHT).
!      Value is greater than all data values:
!    LEFT = N, RIGHT = -1, and X(LEFT) < XVAL.
!      Value is equal to a data value:
!    LEFT = RIGHT, and X(LEFT) = X(RIGHT) = XVAL.
!
  integer n
!
  integer high
  integer left
  integer low
  integer right
  integer start
  real x(n)
  real xval
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_BRACKET2 - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if

  if ( start < 1 .or. start > n ) then
    start = ( n + 1 ) / 2
  end if
!
!  XVAL = X(START)?
!
  if ( x(start) == xval ) then

    left = start
    right = start
    return
!
!  X(START) < XVAL?
!
  else if ( x(start) < xval ) then
!
!  X(START) = X(N) < XVAL < Infinity?
!
    if ( start + 1 > n ) then

      left = start
      right = - 1
      return
!
!  XVAL = X(START+1)?
!
    else if ( xval == x(start+1) ) then

      left = start + 1
      right = start + 1
      return
!
!  X(START) < XVAL < X(START+1)?
!
    else if ( xval < x(start+1) ) then

      left = start
      right = start + 1
      return
!
!  X(START+1) = X(N) < XVAL < Infinity?
!
    else if ( start + 2 > n ) then

      left = start + 1
      right = - 1
      return
!
!  XVAL = X(START+2)?
!
    else if ( xval == x(start+2) ) then

      left = start + 2
      right = start + 2
      return
!
!  X(START+1) < XVAL < X(START+2)?
!
    else if ( xval < x(start+2) ) then

      left = start + 1
      right = start + 2
      return
!
!  Binary search for XVAL in [ X(START+2), X(N) ],
!  where XVAL is guaranteed to be greater than X(START+2).
!
    else

      low = start + 2
      high = n
      call rvec_bracket ( high + 1 - low, x(low), xval, left, right )
      left = left + low - 1
      right = right + low - 1

    end if
!
!  - Infinity < XVAL < X(START) = X(1).
!
  else if ( start == 1 ) then

    left = - 1
    right = start
    return
!
!  XVAL = X(START-1)?
!
  else if ( xval == x(start-1) ) then

    left = start - 1
    right = start - 1
    return
!
!  X(START-1) < XVAL < X(START)?
!
  else if ( x(start-1) <= xval ) then

    left = start - 1
    right = start
    return
!
!  Binary search for XVAL in [ X(1), X(START-1) ],
!  where XVAL is guaranteed to be less than X(START-1).
!
  else

    low = 1
    high = start - 1
    call rvec_bracket ( high + 1 - low, x(1), xval, left, right )

  end if

  return
end
subroutine rvec_bracket3 ( n, t, tval, left )
!
!*******************************************************************************
!
!! RVEC_BRACKET3 finds the interval containing or nearest a given value.
!
!
!  Discussion:
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of the input array.
!
!    Input, real T(N), an array that has been sorted into ascending order.
!
!    Input, real TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer LEFT.
!
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  integer n
!
  integer high
  integer left
  integer low
  integer mid
  real t(n)
  real tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_BRACKET3 - Fatal error!'
    write ( *, * ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. left > n - 1 ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( tval >= t(left-1) ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( tval >= t(mid) ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in {T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( tval > t(left+1) ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( tval >= t(n-1) ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( tval >= t(mid) ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

  end if

  return
end
subroutine rvec_compare ( n, a1, a2, isgn )
!
!*******************************************************************************
!
!! RVEC_COMPARE compares two real vectors.
!
!
!  Discussion:
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2.0, 6.0, 2.0 )
!      A2 = ( 2.0, 8.0, 12.0 )
!
!    Output:
!
!      ISGN = -1
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, integer A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  integer n
!
  real a1(n)
  real a2(n)
  integer isgn
  integer k
!
  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = - 1
      return
    else if ( a1(k) > a2(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine rvec_convolve_circ ( n, x, y, z )
!
!*******************************************************************************
!
!! RVEC_CONVOLVE_CIRC returns the discrete circular convolution of two vectors.
!
!
!  Formula:
!
!    z(1+m) = xCCy(m) = sum ( 0 <= k <= n-1 ) x(1+k) * y(1+m-k)
!
!    Here, if the index of Y becomes nonpositive, it is "wrapped around"
!    by having N added to it.
!
!  Discussion:
!
!    The circular convolution is equivalent to multiplication of Y by a
!    circulant matrix formed from the vector X.
!
!  Example:
!
!    Input:
!
!      X = (/ 1, 2, 3, 4 /)
!      Y = (/ 1, 2, 4, 8 /)
!
!    Output:
!
!      Circulant form:
!
!      Z = ( 1 4 3 2 )   ( 1 )
!          ( 2 1 4 3 )   ( 2 )
!          ( 3 2 1 4 ) * ( 4 )
!          ( 4 3 2 1 )   ( 8 )
!
!      The formula:
!
!      Z = (/ 1*1 + 2*8 + 3*4 + 4*2,
!             1*2 + 2*1 + 3*8 + 4*4,
!             1*4 + 2*2 + 3*1 + 4*8,
!             1*8 + 2*4 + 3*2 + 4*1 /)
!
!      Result:
!
!      Z = (/ 37, 44, 43, 26 /)
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
!    Input, integer N, the dimension of the vectors.
!
!    Input, real X(N), Y(N), the vectors to be circularly convolved.
!
!    Output, real Z(N), the circular convolution of X and Y.
!
  integer n
!
  integer m
  real x(n)
  real y(n)
  real z(n)
!
  do m = 1, n
    z(m) = dot_product ( x(1:m), y(m:1:-1) ) &
         + dot_product ( x(m+1:n), y(n:m+1:-1) )
  end do

  return
end
subroutine rvec_cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! RVEC_CROSS_3D computes the cross product of two vectors in 3D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real X3, Y3, Z3, the cross product vector.
!
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  x3 = y1 * z2 - z1 * y2
  y3 = z1 * x2 - x1 * z2
  z3 = x1 * y2 - y1 * x2

  return
end
subroutine rvec_cum ( n, a, a_cum )
!
!*******************************************************************************
!
!! RVEC_CUM computes the cumulutive sum of the entries of a vector.
!
!
!  Example:
!
!    Input:
!
!      A = ( 1.0, 2.0, 3.0, 4.0 )
!
!    Output:
!
!      A_CUM = ( 0.0, 1.0, 3.0, 6.0, 10.0 )
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector to be summed.
!
!    Output, real A_CUM(N+1), the cumulative sum of the entries of A.
!
  integer n
!
  real a(n)
  real a_cum(n+1)
  real a_sum
  integer i
!
  a_sum = 0.0E+00

  do i = 1, n
    a_cum(i) = a_sum
    a_sum = a_sum + a(i)
  end do

  a_cum(n+1) = a_sum

  return
end
subroutine rvec_dif ( cof, h, n )
!
!*******************************************************************************
!
!! RVEC_DIF computes coefficients for estimating the N-th derivative.
!
!
!  Discussion:
!
!    The routine computes the N+1 coefficients for a centered finite difference
!    estimate of the N-th derivative of a function.
!
!    The estimate has the form
!
!      FDIF(N,X) = Sum (I = 0 to N) COF(I) * F ( X(I) )
!
!    To understand the computation of the coefficients, it is enough
!    to realize that the first difference approximation is
!
!      FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)
!
!    and that the second difference approximation can be regarded as
!    the first difference approximation repeated:
!
!      FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
!         = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)
!
!    and so on for higher order differences.
!
!    Thus, the next thing to consider is the integer coefficients of
!    the sampled values of F, which are clearly the Pascal coefficients,
!    but with an alternating negative sign.  In particular, if we
!    consider row I of Pascal's triangle to have entries j = 0 through I,
!    then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
!    and P(0,0) = 1.
!
!       1
!      -1  1
!       1 -2   1
!      -1  3  -3   1
!       1 -4   6  -4   1
!      -1  5 -10  10  -5  1
!       1 -6  15 -20  15 -6 1
!
!    Next, note that the denominator of the approximation for the
!    N-th derivative will be (2*DX)**N.
!
!    And finally, consider the location of the N+1 sampling
!    points for F:
!
!      X-N*DX, X-(N-2)*DX, X-(N-4)*DX, ..., X+(N-4)*DX, X+(N-2*DX), X+N*DX.
!
!    Thus, a formula for evaluating FDIF(N,X) is
!
!      fdif = 0.0
!      do i = 0, n
!        xi = x + (2*i-n) * h
!        fdif = fdif + cof(i) * f(xi)
!      end do
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real COF(0:N), the coefficients needed to approximate
!    the N-th derivative of a function F.
!
!    Input, real H, the half spacing between points.  H must be positive.
!
!    Input, integer N, the order of the derivative to be approximated.
!    N must be 0 or greater.
!
  integer n
!
  real cof(0:n)
  real h
  integer i
  integer j
!
  if ( n < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_DIF - Fatal error!'
    write ( *, * ) '  Derivative order N = ', n
    write ( *, * ) '  but N must be at least 0.'
    stop
  end if

  if ( h <= 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_DIF - Fatal error!'
    write ( *, * ) '  The half sampling spacing is H = ', h
    write ( *, * ) '  but H must be positive.'
    stop
  end if

  do i = 0, n

    cof(i) = 1.0E+00

    do j = i-1, 1, -1
      cof(j) = - cof(j) + cof(j-1)
    end do

    if ( i > 0 ) then
      cof(0) = - cof(0)
    end if

  end do

  cof(0:n) = cof(0:n) / ( 2.0E+00 * h )**n

  return
end
function rvec_distinct ( n, a )
!
!*******************************************************************************
!
!! RVEC_DISTINCT is true if the entries in a real vector are distinct.
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
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector to be checked.
!
!    Output, logical RVEC_DISTINCT is .TRUE. if the elements of A are distinct.
!
  integer n
!
  real a(n)
  integer i
  integer j
  logical rvec_distinct
!
  rvec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1
      if ( a(i) == a(j) ) then
        return
      end if
    end do
  end do

  rvec_distinct = .true.

  return
end
function rvec_dot0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! RVEC_DOT0_3D computes the dot product of (P1-P2) and (P3-P2) in 3D.
!
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
!    the points P1, P2 and P3.
!
!    Output, real RVEC_DOT0_3D, the dot product of (P1-P2) and (P3-P2).
!
  real rvec_dot0_3d
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  rvec_dot0_3d = ( x1 - x2 ) * ( x3 - x2 ) + &
                 ( y1 - y2 ) * ( y3 - y2 ) + &
                 ( z1 - z2 ) * ( z3 - z2 )

  return
end
function rvec_dot_2d ( x1, y1, x2, y2 )
!
!*******************************************************************************
!
!! RVEC_DOT_2D computes the dot product of a pair of vectors in 2D.
!
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the coordinates of the vectors.
!
!    Output, real RVEC_DOT_2D, the dot product.
!
  real rvec_dot_2d
  real x1
  real x2
  real y1
  real y2
!
  rvec_dot_2d = x1 * x2 + y1 * y2

  return
end
function rvec_dot_3d ( x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! RVEC_DOT_3D computes the dot product of a pair of vectors in 3D.
!
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real RVEC_DOT_3D, the dot product.
!
  real rvec_dot_3d
  real x1
  real x2
  real y1
  real y2
  real z1
  real z2
!
  rvec_dot_3d = x1 * x2 + y1 * y2 + z1 * z2

  return
end
function rvec_eq ( n, a1, a2 )
!
!*******************************************************************************
!
!! RVEC_EQ is true if every pair of entries in two vectors is equal.
!
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real A1(N), A2(N), two vectors to compare.
!
!    Output, logical RVEC_EQ.
!    RVEC_EQ is .TRUE. if every pair of elements A1(I) and A2(I) are equal,
!    and .FALSE. otherwise.
!
  integer n
!
  real a1(n)
  real a2(n)
  logical rvec_eq
!
  rvec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
subroutine rvec_even ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  integer n
!
  real a(n)
  real ahi
  real alo
  integer i
!
  if ( n == 1 ) then

    a(1) = 0.5E+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
    end do

  end if

  return
end
subroutine rvec_even2 ( maxval, nfill, nold, nval, xold, xval )
!
!*******************************************************************************
!
!! RVEC_EVEN2 linearly interpolates new numbers into a vector of data.
!
!
!  Discussion:
!
!    The number of values created between two old values can vary from
!    one pair of values to the next.
!
!    The interpolated values are evenly spaced.
!
!    RVEC_EVEN2 is a generalization of RVEC_EVEN.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXVAL, the size of the XVAL array, as declared by the
!    user.  MAXVAL must be large enough to hold the NVAL values computed by
!    this routine.  In other words, MAXVAL must be at least equal to
!    NOLD + SUM(i = 1 to NOLD-1) NFILL(I).
!
!    Input, integer NFILL(NOLD-1), the number of values
!    to be interpolated between XOLD(I) and XOLD(I+1).
!    NFILL(I) does not count the endpoints.  Thus, if
!    NFILL(I) is 1, there will be one new point generated
!    between XOLD(I) and XOLD(I+1).
!
!    NFILL(I) must be nonnegative.
!
!    Input, integer NOLD, the number of values XOLD,
!    between which extra values are to be interpolated.
!
!    Output, integer NVAL, the number of values computed
!    in the XVAL array.
!
!      NVAL = NOLD + SUM(I = 1 to NOLD-1) NFILL(I)
!
!    Input, real XOLD(NOLD), the original vector of numbers
!    between which new values are to be interpolated.
!
!    Output, real XVAL(MAXVAL).  On output, XVAL contains the
!    NOLD values of XOLD, as well as the interpolated
!    values, making a total of NVAL values.
!
  integer maxval
  integer nold
!
  integer i
  integer nadd
  integer nfill(nold-1)
  integer nval
  real xold(nold)
  real xval(maxval)
!
  nval = 1

  do i = 1, nold-1

    if ( nfill(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'RVEC_EVEN2 - Fatal error!'
      write ( *, * ) '  NFILL(I) is negative for I = ', i
      write ( *, * ) '  NFILL(I) = ', nfill(i)
      stop
    end if

    if ( nval + nfill(i) + 1 > maxval ) then
      write ( *, * ) ' '
      write ( *, * ) 'RVEC_EVEN2 - Fatal error!'
      write ( *, * ) '  MAXVAL is not large enough.  '
      write ( *, * ) '  MAXVAL = ', maxval
      write ( *, * ) '  which is exceeded by storage requirements'
      write ( *, * ) '  for interpolating in interval ', i
      stop
    end if

    nadd = nfill(i) + 2

    call rvec_even ( xold(i), xold(i+1), nadd, xval(nval) )

    nval = nval + nfill(i) + 1

  end do

  return
end
subroutine rvec_even3 ( nold, nval, xold, xval )
!
!*******************************************************************************
!
!! RVEC_EVEN3 evenly interpolates new data into a vector.
!
!
!  Discussion:
!
!    RVEC_EVEN3 accepts a short vector of numbers, and returns a longer
!    vector of numbers, created by interpolating new values between
!    the given values.
!
!    Between any two original values, new values are evenly interpolated.
!
!    Over the whole vector, the new numbers are interpolated in
!    such a way as to try to minimize the largest distance interval size.
!
!    The algorithm employed is not "perfect".
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NOLD, the number of values XOLD, between which extra
!    values are to be interpolated.
!
!    Input, integer NVAL, the number of values to be computed
!    in the XVAL array.  NVAL should be at least NOLD.
!
!    Input, real XOLD(NOLD), the original vector of numbers
!    between which new values are to be interpolated.
!
!    Output, real XVAL(NVAL).  On output, XVAL contains the
!    NOLD values of XOLD, as well as interpolated
!    values, making a total of NVAL values.
!
  integer nval
  integer nold
!
  real density
  integer i
  integer ival
  integer nmaybe
  integer npts
  integer ntemp
  integer ntot
  real xlen
  real xleni
  real xlentot
  real xold(nold)
  real xval(nval)
!
  xlen = 0
  do i = 1, nold - 1
    xlen = xlen + abs ( xold(i+1) - xold(i) )
  end do

  ntemp = nval - nold

  density = real ( ntemp ) / xlen

  ival = 1
  ntot = 0
  xlentot = 0.0E+00

  do i = 1, nold-1

    xleni = abs ( xold(i+1) - xold(i) )
    npts = int ( density * xleni )
    ntot = ntot + npts
!
!  Determine if we have enough left-over density that it should
!  be changed into a point.  A better algorithm would agonize
!  more over where that point should go.
!
    xlentot = xlentot + xleni
    nmaybe = nint ( xlentot * density )

    if ( nmaybe > ntot ) then
      npts = npts + nmaybe - ntot
      ntot = nmaybe
    end if

    call rvec_even ( xold(i), xold(i+1), npts+2, xval(ival) )
    ival = ival + npts + 1

  end do

  return
end
subroutine rvec_even_select ( xlo, xhi, n, ival, xval )
!
!*******************************************************************************
!
!! RVEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!
!  Formula:
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XLO, XHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Input, integer IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
!
  integer n
!
  integer ival
  real xhi
  real xlo
  real xval
!
  if ( n == 1 ) then

    xval = 0.5E+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival ) * xlo + real ( ival - 1 ) * xhi ) &
      / real ( n - 1 )

  end if

  return
end
subroutine rvec_frac ( n, a, k, afrac )
!
!*******************************************************************************
!
!! RVEC_FRAC searches for the K-th smallest entry in an N-vector.
!
!
!  Discussion:
!
!    Hoare's algorithm is used.
!
!  Modified:
!
!    17 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, real A(N).
!    On input, A is the array to search.
!    On output, the elements of A have been somewhat rearranged.
!
!    Input, integer K, the fractile to be sought.  If K = 1, the minimum
!    entry is sought.  If K = N, the maximum is sought.  Other values
!    of K search for the entry which is K-th in size.  K must be at
!    least 1, and no greater than N.
!
!    Output, real AFRAC, the value of the K-th fractile of A.
!
  integer n
!
  real a(n)
  real afrac
  integer i
  integer iryt
  integer j
  integer k
  integer left
  real w
  real x
!
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_FRAC - Fatal error!'
    write ( *, * ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_FRAC - Fatal error!'
    write ( *, * ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( k > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_FRAC - Fatal error!'
    write ( *, * ) '  Illegal K > N, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( left >= iryt ) then
      afrac = a(k)
      exit
    end if

    x = a(k)
    i = left
    j = iryt

    do

      if ( i > j ) then
        if ( j < k ) then
          left = i
        end if
        if ( k < i ) then
          iryt = j
        end if
        exit
      end if
!
!  Find I so that X< = A(I)
!
      do while ( a(i) < x )
        i = i + 1
      end do
!
!  Find J so that A(J) < =  X
!
      do while ( x < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call r_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine rvec_heap_a ( n, a )
!
!*******************************************************************************
!
!! RVEC_HEAP_A reorders an array of reals into an ascending heap.
!
!
!  Definition:
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the size of the input array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  integer n
!
  real a(n)
  integer i
  integer ifree
  real key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( m > n ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m+1) > a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine rvec_heap_d ( n, a )
!
!*******************************************************************************
!
!! RVEC_HEAP_D reorders an array of reals into a descending heap.
!
!
!  Definition:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the size of the input array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  integer n
!
  real a(n)
  integer i
  integer ifree
  real key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( m > n ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) >= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, the value KEY
!  moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine rvec_house_column ( n, a, k, v )
!
!*******************************************************************************
!
!! RVEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
!
!
!  Discussion:
!
!    The routine returns a vector V that defines a Householder
!    premultiplier matrix H(V) that zeros out the subdiagonal entries of
!    column K of the matrix A.
!
!      H(V) = I - 2 * v * v'
!
!  Modified:
!
!    25 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix A.
!
!    Input, real A(N), column K of the matrix A.
!
!    Input, integer K, the column of the matrix to be modified.
!
!    Output, real V(N), a vector of unit L2 norm which defines an
!    orthogonal Householder premultiplier matrix H with the property
!    that the K-th column of H*A is zero below the diagonal.
!
  integer n
!
  real a(n)
  real alpha
  real column_norm
  integer i
  integer k
  real rvec_norm2
  real v(n)
!
  v(1:n) = 0.0E+00

  if ( k < 1 .or. k >= n ) then
    return
  end if

  column_norm = rvec_norm2 ( n+1-k, a(k) )

  if ( column_norm == 0.0E+00 ) then
    return
  end if

  alpha = - sign ( 1.0E+00, a(k) ) * column_norm

  v(k) = sqrt ( 0.5E+00 * ( 1.0E+00 - a(k) / alpha ) )

  v(k+1:n) = - 0.5E+00 * a(k+1:n) / ( alpha * v(k) )

  return
end
subroutine rvec_iamax ( n, a, iamax )
!
!*******************************************************************************
!
!! RVEC_IAMAX returns the index of the maximum absolute value in a real vector.
!
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real X(N), the array.
!
!    Output, integer IAMAX, the index of the entry of largest magnitude.
!
  integer n
!
  real a(n)
  real aamax
  integer i
  integer iamax
!
  if ( n <= 0 ) then

    iamax = 0.0E+00

  else

    iamax = 1
    aamax = abs ( a(1) )

    do i = 2, n
      if ( abs ( a(i) ) > aamax ) then
        iamax = i
        aamax = abs ( a(i) )
      end if
    end do

  end if

  return
end
subroutine rvec_iamin ( n, a, iamin )
!
!*******************************************************************************
!
!! RVEC_IAMIN returns the index of the minimum absolute value in a real vector.
!
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real X(N), the array.
!
!    Output, integer IAMIN, the index of the entry of smallest magnitude.
!
  integer n
!
  real a(n)
  real aamin
  integer i
  integer iamin
!
  if ( n <= 0 ) then

    iamin = 0

  else

    iamin = 1
    aamin = abs ( a(1) )

    do i = 2, n
      if ( abs ( a(i) ) < aamin ) then
        iamin = i
        aamin = abs ( a(i) )
      end if
    end do

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
subroutine rvec_imax ( n, a, imax )
!
!*******************************************************************************
!
!! RVEC_IMAX returns the index of the maximum value in a real vector.
!
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), the array.
!
!    Output, integer IMAX, the index of the largest entry.
!
  integer n
!
  real a(n)
  real amax
  integer i
  integer imax
!
  if ( n <= 0 ) then

    imax = 0

  else

    imax = 1
    amax = a(1)

    do i = 2, n
      if ( a(i) > amax ) then
        amax = a(i)
        imax = i
      end if
    end do

  end if

  return
end
subroutine rvec_imin ( n, a, imin )
!
!*******************************************************************************
!
!! RVEC_IMIN returns the index of the minimum value in a real vector.
!
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), the array.
!
!    Output, integer IMIN, the index of the smallest entry.
!
  integer n
!
  real a(n)
  real amin
  integer i
  integer imin
!
  if ( n <= 0 ) then

    imin = 0

  else

    imin = 1
    amin = a(1)

    do i = 2, n
      if ( a(i) < amin ) then
        amin = a(i)
        imin = i
      end if
    end do

  end if

  return
end
subroutine rvec_index_delete_all ( n, x, indx, xval )
!
!*******************************************************************************
!
!! RVEC_INDEX_DELETE_ALL deletes all occurrences of a real value from an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, real X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer equal1
  integer equal2
  integer get
  integer i
  integer indx(*)
  integer j
  integer less
  integer more
  integer put
  real x(*)
  real xval
!
  if ( n < 1 ) then
    n = 0
    return
  end if

  call rvec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal /= 0 ) then

    equal1 = equal

    do

      if ( equal1 <= 1 ) then
        exit
      end if

      if ( x(indx(equal1-1)) /= xval ) then
        exit
      end if

      equal1 = equal1 - 1

    end do

    equal2 = equal

    do

      if ( equal2 >= n ) then
        exit
      end if

      if ( x(indx(equal2+1)) /= xval ) then
        exit
      end if

      equal2 = equal2 + 1

    end do
!
!  Discard certain X values.
!
    put = 0

    do get = 1, n

      if ( x(get) /= xval ) then
        put = put + 1
        x(put) = x(get)
      end if

    end do

    x(put+1:n) = 0.0E+00
!
!  Adjust the INDX values.
!
    do equal = equal1, equal2
      do i = 1, n
        if ( indx(i) > indx(equal) ) then
          indx(i) = indx(i) - 1
        end if
      end do
    end do
!
!  Discard certain INDX values.
!
    indx(equal1:n+equal1-equal2-1) = indx(equal2+1:n)
    indx(n+equal1-equal2:n) = 0
!
!  Adjust N.
!
    n = put

  end if

  return
end
subroutine rvec_index_delete_dupes ( n, x, indx, n2 )
!
!*******************************************************************************
!
!! RVEC_INDEX_DELETE_DUPES deletes duplicate real values from an indexed sorted list.
!
!
!  Discussion:
!
!    If any value occurs more than once in the input list, all duplicate
!    values are removed.
!
!    On output, the list has been sorted, and the index array has
!    been set to the identity.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input/output, real X(N), the list.  On output, the list
!    has only unique entries, and they have been sorted explicitly.
!
!    Input/output, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, the value to be sought.
!
!    Output, integer N2, the number of unique entries in X.
!
  integer n
!
  integer i
  integer indx(*)
  integer n2
  real x(n)
  real y(n)
!
  i = 0
  n2 = 0

  do

    i = i + 1

    if ( i > n ) then
      exit
    end if

    if ( i > 1 ) then
      if ( x(indx(i)) == y(n2) ) then
        cycle
      end if
    end if

    n2 = n2 + 1
    y(n2) = x(indx(i))

  end do

  x(1:n2) = y(1:n2)
  call ivec_identity ( n2, indx )

  return
end
subroutine rvec_index_delete_one ( n, x, indx, xval )
!
!*******************************************************************************
!
!! RVEC_INDEX_DELETE_ONE deletes one copy of a real value from an indexed sorted list.
!
!
!  Discussion:
!
!    If the value occurs in the list more than once, only one copy is deleted.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, real X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer i
  integer indx(*)
  integer j
  integer less
  integer more
  real x(*)
  real xval
!
  if ( n < 1 ) then
    n = 0
    return
  end if

  call rvec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal /= 0 ) then
    j = indx(equal)
    x(j:n-1) = x(j+1:n)
    indx(equal:n-1) = indx(equal+1:n)
    do i = 1, n-1
      if ( indx(i) > j ) then
        indx(i) = indx(i) - 1
      end if
    end do
    n = n - 1
  end if

  return
end
subroutine rvec_index_insert ( n, x, indx, xval )
!
!*******************************************************************************
!
!! RVEC_INDEX_INSERT inserts a real value in an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, real X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer indx(*)
  integer less
  integer more
  real x(*)
  real xval
!
  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if

  call rvec_index_search ( n, x, indx, xval, less, equal, more )

  x(n+1) = xval
  indx(n+1:more+1:-1) = indx(n:more:-1)
  indx(more) = n + 1
  n = n + 1

  return
end
subroutine rvec_index_insert_unique ( n, x, indx, xval )
!
!*******************************************************************************
!
!! RVEC_INDEX_INSERT_UNIQUE inserts a unique real value in an indexed sorted list.
!
!
!  Discussion:
!
!    If the value does not occur in the list, it is included in the list,
!    and N, X and INDX are updated.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N, the size of the current list.
!
!    Input/output, real X(N), the list.
!
!    Input/output, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, the value to be sought.
!
  integer n
!
  integer equal
  integer indx(*)
  integer less
  integer more
  real x(*)
  real xval
!
  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if
!
!  Does XVAL already occur in X?
!
  call rvec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    x(n+1) = xval
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1
  end if

  return
end
subroutine rvec_index_order ( n, x, indx )
!
!*******************************************************************************
!
!! RVEC_INDEX_ORDER sorts an integer vector using an index vector.
!
!
!  Discussion:
!
!    The index vector itself is not modified.  Therefore, the pair
!    (X,INDX) no longer represents an index sorted vector.  If this
!    relationship is to be preserved, then simply set INDX(1:N)=(1:N).
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input/output, real X(N), the list.  On output, the list
!    has been sorted.
!
!    Input, integer INDX(N), the sort index of the list.
!
  integer n
!
  integer i
  integer indx(n)
  real x(n)
  real y(n)
!
  y(1:n) = x(indx(1:n))
  x(1:n) = y(1:n)

  return
end
subroutine rvec_index_search ( n, x, indx, xval, less, equal, more )
!
!*******************************************************************************
!
!! RVEC_INDEX_SEARCH searches for a real value in an indexed sorted list.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input, real X(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, real XVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  integer n
!
  integer equal
  integer hi
  integer indx(n)
  integer less
  integer lo
  integer mid
  integer more
  real x(n)
  real xhi
  real xlo
  real xmid
  real xval
!
  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n
  xlo = x(indx(lo))
  xhi = x(indx(hi))

  if ( xval < xlo ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( xval == xlo ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  if ( xval > xhi ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( xval == xhi ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))

    if ( xval == xmid ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( xval < xmid ) then
      hi = mid
    else if ( xval > xmid ) then
      lo = mid
    end if

  end do

  return
end
subroutine rvec_index_sort_unique ( n, x, indx, n2 )
!
!*******************************************************************************
!
!! RVEC_INDEX_SORT_UNIQUE creates a sort index for a real vector.
!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the current list.
!
!    Input/output, real X(N), the list.  On output, X contains only
!    unique elements.
!
!    Output, integer INDX(N), the sort index of the list.
!
!    Output, integer N2, the number of unique elements in X.
!
  integer n
!
  integer i
  integer indx(n)
  integer n2
  real x(n)
  real y(n)
!
  n2 = 0

  do i = 1, n
    call rvec_index_insert_unique ( n2, y, indx, x(i) )
  end do

  x(1:n2) = y(1:n2)

  x(n2+1:n) = 0.0E+00
  indx(n2+1:n) = 0

  return
end
subroutine rvec_insert ( n, a, pos, value )
!
!*******************************************************************************
!
!! RVEC_INSERT inserts a value into an array.
!
!
!  Modified:
!
!    17 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the array on input.
!
!    Input/output, real A(N+1), the array.  On input, A is assumed to
!    contain only N entries, while on output, A actually contains N+1 entries.
!
!    Input, integer POS, the position to be assigned the new entry.
!    1 <= POS <= N+1.
!
!    Input, real VALUE, the value to be inserted at the given position.
!
  integer n
!
  real a(n+1)
  integer i
  integer pos
  real value
!
  if ( pos < 1 .or. pos > n+1 ) then

    write ( *, * ) ' '
    write ( *, * ) 'RVEC_INSERT - Fatal error!'
    write ( *, * ) '  Illegal insertion position = ', pos
    stop

  else

    do i = n+1, pos+1, -1
      a(i) = a(i-1)
    end do

    a(pos) = value

  end if

  return
end
subroutine rvec_max ( n, a, amax )
!
!*******************************************************************************
!
!! RVEC_MAX returns the maximum value in a real vector.
!
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), the array.
!
!    Output, real AMAX, the value of the largest entry.
!
  integer n
!
  real a(n)
  real amax
!
  amax = maxval ( a(1:n) )

  return
end
subroutine rvec_mean ( n, a, mean )
!
!*******************************************************************************
!
!! RVEC_MEAN returns the mean of a real vector.
!
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector whose mean is desired.
!
!    Output, real MEAN, the mean, or average, of the vector entries.
!
  integer n
!
  real a(n)
  real mean
!
  mean = sum ( a(1:n) ) / real ( n )

  return
end
subroutine rvec_median ( n, a, median )
!
!*******************************************************************************
!
!! RVEC_MEDIAN returns the median of an unsorted real vector.
!
!
!  Discussion:
!
!    Hoare's algorithm is used.  The values of the vector are
!    rearranged by this routine.
!
!  Modified:
!
!    18 September 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, real A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.
!
!    Output, real MEDIAN, the value of the median of A.
!
  integer n
!
  real a(n)
  integer k
  real median
!
  k = ( n + 1 ) / 2

  call rvec_frac ( n, a, k, median )

  return
end
subroutine rvec_merge_a ( na, a, nb, b, nc, c )
!
!*******************************************************************************
!
!! RVEC_MERGE_A merges two ascending sorted real arrays.
!
!
!  Discussion:
!
!    The elements of A and B should be sorted in ascending order.
!
!    The elements in the output array C will also be in ascending order,
!    and unique.
!
!    The output vector C may share storage with A or B.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real A(NA), the first sorted array.
!
!    Input, integer NB, the dimension of B.
!
!    Input, real B(NB), the second sorted array.
!
!    Output, integer NC, the number of elements in the output array.
!    Note that C should usually be dimensioned at least NA+NB in the
!    calling routine.
!
!    Output, real C(NC), the merged unique sorted array.
!
  integer na
  integer nb
!
  real a(na)
  real b(nb)
  real c(na+nb)
  real d(na+nb)
  integer j
  integer ja
  integer jb
  integer na2
  integer nb2
  integer nc
  integer order
!
  na2 = na
  nb2 = nb
!
  ja = 0
  jb = 0
  nc = 0

  call rvec_order_type ( na2, a, order )

  if ( order < 0 .or. order > 2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_MERGE - Fatal error!'
    write ( *, * ) '  The input array A is not ascending sorted!'
    stop
  end if

  call rvec_order_type ( nb2, b, order )

  if ( order < 0 .or. order > 2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_MERGE - Fatal error!'
    write ( *, * ) '  The input array B is not ascending sorted!'
    stop
  end if

  do
!
!  If we've used up all the entries of A, stick the rest of B on the end.
!
    if ( ja >= na2 ) then

      do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = b(jb)
        else if ( d(nc) < b(jb) ) then
          nc = nc + 1
          d(nc) = b(jb)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  If we've used up all the entries of B, stick the rest of A on the end.
!
    else if ( jb >= nb2 ) then

      do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = a(ja)
        else if ( d(nc) < a(ja) ) then
          nc = nc + 1
          d(nc) = a(ja)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  Otherwise, if the next entry of A is smaller, that's our candidate.
!
    else if ( a(ja+1) <= b(jb+1) ) then

      ja = ja + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = a(ja)
      else if ( d(nc) < a(ja) ) then
        nc = nc + 1
        d(nc) = a(ja)
      end if
!
!  ...or if the next entry of B is the smaller, consider that.
!
    else

      jb = jb + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = b(jb)
      else if ( d(nc) < b(jb) ) then
        nc = nc + 1
        d(nc) = b(jb)
      end if
    end if

  end do

  return
end
subroutine rvec_min ( n, a, amin )
!
!*******************************************************************************
!
!! RVEC_MIN returns the minimum value of a real array.
!
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), the array.
!
!    Output, real AMIN, the value of the smallest entry.
!
  integer n
!
  real a(n)
  real amin
!
  amin = minval ( a(1:n) )

  return
end
function rvec_norm1 ( n, a )
!
!*******************************************************************************
!
!! RVEC_NORM1 returns the 1-norm of a vector.
!
!
!  Definition:
!
!    The vector 1-norm is defined as:
!
!      RVEC_NORM1 = Sum ( 1 <= I <= N ) abs ( A(I) ).
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real A(N), the vector whose 1-norm is desired.
!
!    Output, real RVEC_NORM1, the 1-norm of A.
!
  integer n
!
  real a(n)
  real rvec_norm1
!
  rvec_norm1 = sum ( abs ( a(1:n) ) )

  return
end
function rvec_norm2 ( n, a )
!
!*******************************************************************************
!
!! RVEC_NORM2 returns the 2-norm of a vector.
!
!
!  Definition:
!
!    The vector 2-norm is defined as:
!
!      RVEC_NORM2 = Sqrt ( Sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real A(N), the vector whose 2-norm is desired.
!
!    Output, real RVEC_NORM2, the 2-norm of A.
!
  integer n
!
  real a(n)
  real rvec_norm2
!
  rvec_norm2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function rvec_normi ( n, a )
!
!*******************************************************************************
!
!! RVEC_NORMI returns the infinity-norm of a vector.
!
!
!  Definition:
!
!    The vector infinity-norm is defined as:
!
!      RVEC_NORMI = Max ( 1 <= I <= N ) abs ( A(I) ).
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real A(N), the vector whose infinity-norm is desired.
!
!    Output, real RVEC_NORMI, the infinity-norm of A.
!
  integer n
!
  real a(n)
  real rvec_normi
!
  rvec_normi = maxval ( abs ( a(1:n) ) )

  return
end
subroutine rvec_order_type ( n, a, order )
!
!*******************************************************************************
!
!! RVEC_ORDER_TYPE determines if a real array is (non)strictly ascending/descending.
!
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the array.
!
!    Input, real A(N), the array to be checked.
!
!    Output, integer ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  integer n
!
  real a(n)
  integer i
  integer order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( i > n ) then
      order = 0
      return
    end if

    if ( a(i) > a(1) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do

    i = i + 1
    if ( i > n ) then
      exit
    end if

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine rvec_part_quick_a ( n, a, l, r )
!
!*******************************************************************************
!
!! RVEC_PART_QUICK_A reorders a real vector as part of a quick sort.
!
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1) as the key,
!    all entries of A that are less than or equal to the key will
!    precede the key which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
!
!    Output:
!
!      L = 3, R = 6
!
!      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
!            -------        -------
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of A.
!
!    Input/output, real A(N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1).  Then
!    I <= L                 A(I) < KEY;
!         L < I < R         A(I) = KEY;
!                 R <= I    A(I) > KEY.
!
  integer n
!
  real a(n)
  integer i
  real key
  integer l
  integer m
  integer r
  real temp
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_PART_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key = a(1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( a(l+1) > key ) then
      r = r - 1
      temp = a(r)
      a(r) = a(l+1)
      a(l+1) = temp
    else if ( a(l+1) == key ) then
      m = m + 1
      temp = a(m)
      a(m) = a(l+1)
      a(l+1) = temp
      l = l + 1
    else if ( a(l+1) < key ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(i) = a(i+m)
  end do
!
!  Out of bounds here, occasionally
!
  l = l - m

  a(l+1:l+m) = key

  return
end
subroutine rvec_permute ( n, a, p )
!
!*******************************************************************************
!
!! RVEC_PERMUTE permutes a real vector in place.
!
!
!  Note:
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
!
!    Output:
!
!      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input/output, real A(N), the array to be permuted.
!
!    Input, integer P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  integer n
!
  real a(n)
  real a_temp
  integer i
  integer ierror
  integer iget
  integer iput
  integer istart
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_PERMUTE - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. iget > n ) then
          write ( *, * ) ' '
          write ( *, * ) 'RVEC_PERMUTE - Fatal error!'
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine rvec_permute_random ( n, a )
!
!*******************************************************************************
!
!! RVEC_PERMUTE_RANDOM randomly permutes an real vector.
!
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input/output, real A(N), the array to be permuted.
!
  integer n
!
  real a(n)
  integer p(n)
!
  call perm_random ( n, p )

  call rvec_permute ( n, a, p )

  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real vector.
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
subroutine rvec_print2 ( n, a )
!
!*******************************************************************************
!
!! RVEC_PRINT2 prints out the N vector A.
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
!    Input, integer N, the number of entries of A.
!
!    Input, real A(N), the vector to be printed.
!
  integer n
!
  real a(n)
  real amax
  real amin
  integer i
  character ( len = 8 ) iform
  logical integ
  integer lmax
  integer log10
  logical r_is_int
!
!  Check if all entries are integral.
!
  integ = .true.

  do i = 1, n

    if ( .not. r_is_int ( a(i) ) ) then
      integ = .false.
      exit
    end if

  end do
!
!  Find the range of the array.
!
  amax = maxval ( abs ( a(1:n) ) )
  amin = minval ( abs ( a(1:n) ) )
!
!  Use the information about the maximum size of an entry to
!  compute an intelligent format for use with integer entries.
!
!  Later, we might also do this for real vectors.
!
  lmax = log10 ( amax )

  if ( integ ) then
    write ( iform, '( ''(I'', i2, '')'' )' ) lmax + 3
  else
    iform = ' '
  end if

  do i = 1, n

    if ( integ ) then
      write ( *, iform ) int ( a(i) )
    else
      write ( *, '(g14.6)' ) a(i)
    end if

  end do

  return
end
subroutine rvec_print_2d ( x, y )
!
!*******************************************************************************
!
!! RVEC_PRINT_2D prints a 2D vector.
!
!
!  Discussion:
!
!    A format is used which suggests a coordinate pair:
!
!  Example:
!
!    ( 1.23, 7.45 )
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
!    Input, real X, Y, the coordinates of the vector.
!
  real x
  real y
!
  write ( *, '( a1, g14.6, a1, g14.6, a1 )' ) '(', x, ',', y, ')'

  return
end
subroutine rvec_print_3d ( x, y, z )
!
!*******************************************************************************
!
!! RVEC_PRINT_3D prints a 3D vector.
!
!
!  Discussion:
!
!    A format is used which suggests a coordinate triple:
!
!  Example:
!
!    ( 1.23, 7.45, -1.45 )
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
!    Input, real X, Y, Z, the coordinates of the vector.
!
  real x
  real y
  real z
!
  write ( *, '( a1, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
    '(', x, ',', y, ',', z, ')'

  return
end
subroutine rvec_print_some ( n, a, max_print )
!
!*******************************************************************************
!
!! RVEC_PRINT_SOME prints "some" of a real vector.
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
!    Input, real A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
  integer n
!
  real a(n)
  integer i
  integer max_print
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i6,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(i6,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

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
!    04 February 2001
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
subroutine rvec_range ( n, x, xmax, xmin, y, ymax, ymin, none )
!
!*******************************************************************************
!
!! RVEC_RANGE finds the range of Y's within a restricted X range.
!
!
!  Discussion:
!
!    The routine is given a set of pairs of points (X,Y), and a range
!    XMIN to XMAX of valid X values.  Over this range, it seeks
!    YMIN and YMAX, the minimum and maximum values of Y for
!    valid X's.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real X(N), the X array.
!
!    Input, real XMAX, XMIN, the largest and smallest entries valid
!    values of X.
!
!    Input, real Y(N), the Y array.
!
!    Output, real YMAX, YMIN, the largest and smallest entries of
!    Y that occur paired with an entry of X between XMIN and XMAX.
!
!    Output, logical NONE, is TRUE if no values of X were valid.
!    In this case, YMIN = YMAX = 0.0, but these values are meaningless.
!
  integer n
!
  logical none
  integer i
  real x(n)
  real xmax
  real xmin
  real y(n)
  real ymax
  real ymin
!
  none = .true.

  ymin = 0.0E+00
  ymax = 0.0E+00

  do i = 1, n

    if ( xmin <= x(i) .and. x(i) <= xmax ) then

      if ( none ) then
        ymin = y(i)
        ymax = y(i)
        none = .false.
      else
        ymin = min ( ymin, y(i) )
        ymax = max ( ymax, y(i) )
      end if

    end if

  end do

  return
end
subroutine rvec_range_2 ( n, a, amax, amin )
!
!*******************************************************************************
!
!! RVEC_RANGE_2 updates a range to include a new array.
!
!
!  Discussion:
!
!    Given a range AMIN to AMAX, and an array A, the routine will
!    decrease AMIN if necessary, or increase AMAX if necessary, so that
!    every entry of A is between AMIN and AMAX.
!
!    However, AMIN will not be increased, nor AMAX decreased.
!
!    This routine may be used to compute the maximum and minimum of a
!    collection of arrays one at a time.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), the array.
!
!    Input/output, real AMAX, AMIN.  On input, AMAX and AMIN
!    are arbitrary values, which may have been computed as
!    the largest and smallest entries of one or more arrays.
!
!    On output, AMAX has possibly been increased, and XMIN
!    has possibly been decreased, so that every entry of X
!    is between AMIN and AMAX.
!
  integer n
!
  real a(n)
  real amax
  real amin
!
  amax = max ( amax, maxval ( a(1:n) ) )
  amin = min ( amin, minval ( a(1:n) ) )

  return
end
subroutine rvec_reverse ( n, a )
!
!*******************************************************************************
!
!! RVEC_REVERSE reverses the elements of a real vector.
!
!
!  Example:
!
!    Input:
!
!      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N), the array to be reversed.
!
  integer n
!
  real a(n)
  integer i
!
  do i = 1, n/2
    call r_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine rvec_rotate ( n, a, m )
!
!*******************************************************************************
!
!! RVEC_ROTATE rotates an object in place.
!
!
!  Note:
!
!    This routine rotates an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A    = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
!
!    Output:
!
!      A    = ( 4.0, 5.0, 1.0, 2.0, 3.0 ).
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
!    Input, integer N, the number of objects.
!
!    Input, integer M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, real A(N), the array to be rotated.
!
  integer n
!
  real a(n)
  integer i_modp
  integer iget
  integer iput
  integer istart
  integer m
  integer mcopy
  integer nset
  real temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( istart > n ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget
  
      iget = iget - mcopy
      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( nset >= n ) then
      exit
    end if

  end do

  return
end
subroutine rvec_search_binary_a ( n, a, aval, indx )
!
!*******************************************************************************
!
!! RVEC_SEARCH_BINARY_A searches an ascending sorted real vector.
!
!
!  Discussion:
!
!    Binary search is used.
!
!  Reference:
!
!    Algorithm 1.9,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the array.
!
!    Input, real A(N), the array to be searched.  The array must
!    be sorted in ascending order.
!
!    Input, real AVAL, the value to be searched for.
!
!    Output, integer INDX, the result of the search.
!    0, AVAL does not occur in the array.
!    I, A(I) = AVAL.
!
  integer n
!
  real a(n)
  real aval
  integer high
  integer indx
  integer low
  integer mid
!
  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == aval ) then
      indx = mid
      exit
    else if ( a(mid) < aval ) then
      low = mid + 1
    else if ( a(mid) > aval ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine rvec_sort2_a ( n, x, y )
!
!*******************************************************************************
!
!! RVEC_SORT2_A ascending sorts a real array and adjusts an associated real array.
!
!
!  Discussion:
!
!    The routine sorts the elements of X, and whenever
!    an element of X is moved, the corresponding element of
!    Y is moved in the same way.  This action means that after
!    the sorting, every element of X is still paired to the
!    same Y value.
!
!    If you have more than one array associated with X, or
!    an integer array, or some other complication, you may want to
!    look at doing an "indexed sort" instead.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input/output, real X(N).  On input, an unsorted array.
!    On output, X has been sorted.
!
!    Input/output, real Y(N), an array which is to be
!    shifted corresponding to the shifts made in X.
!
  integer n
!
  integer i
  integer indx
  integer isgn
  integer j
  real x(n)
  real y(n)
!
  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( indx > 0 ) then

      call r_swap ( x(i), x(j) )
      call r_swap ( y(i), y(j) )

    else if ( indx < 0 ) then

      if ( x(i) <= x(j) ) then
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
subroutine rvec_sort_bubble_a ( n, a )
!
!*******************************************************************************
!
!! RVEC_SORT_BUBBLE_A ascending sorts a real array using bubble sort.
!
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  integer n
!
  real a(n)
  integer i
  integer j
!
  do i = 1, n-1
    do j = i+1, n
      if ( a(i) > a(j) ) then
        call r_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine rvec_sort_heap_a ( n, a )
!
!*******************************************************************************
!
!! RVEC_SORT_HEAP_A ascending sorts a real array using heap sort.
!
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  integer n
!
  real a(n)
  integer n1
!
  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into heap form.
!
  call rvec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call r_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call rvec_heap_a ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call r_swap ( a(1), a(n1) )

  end do

  return
end
subroutine rvec_sort_heap_index_a ( n, a, indx )
!
!*******************************************************************************
!
!! RVEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of a real vector.
!
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call RVEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), an array to be index-sorted.
!
!    Output, integer INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  integer n
!
  real a(n)
  real aval
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l
!
  call ivec_identity ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( l > 1 ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine rvec_sort_heap_index_d ( n, a, indx )
!
!*******************************************************************************
!
!! RVEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of a real vector.
!
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call RVEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), an array to be index-sorted.
!
!    Output, integer INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  integer n
!
  real a(n)
  real aval
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l
!
  call ivec_identity ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( l > 1 ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) > a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval > a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine rvec_sort_insert_a ( n, a )
!
!*******************************************************************************
!
!! RVEC_SORT_INSERT_A ascending sorts a real vector using an insertion sort.
!
!
!  Reference:
!
!    Algorithm 1.1,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, real A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  integer n
!
  real a(n)
  integer i
  integer j
  real x
!
  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( j >= 1 )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine rvec_sort_insert_index_a ( n, a, indx )
!
!*******************************************************************************
!
!! RVEC_SORT_INSERT_INDEX_A ascending index sorts a real vector using insertion.
!
!
!  Reference:
!
!    Algorithm 1.1,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items in the vector.
!    N must be positive.
!
!    Input, real A(N), the array to be sorted.
!
!    Output, integer INDX(N), the sorted indices.  The array is sorted
!    when listed from A(INDX(1)) through A(INDX(N)).
!
  integer n
!
  real a(n)
  integer i
  integer indx(n)
  integer j
  real x
!
  call ivec_identity ( n, indx )

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( j >= 1 )

      if ( a(indx(j)) <= x ) then
        exit
      end if

      indx(j+1) = indx(j)
      j = j - 1

    end do

    indx(j+1) = i

  end do

  return
end
subroutine rvec_sort_insert_index_d ( n, a, indx )
!
!*******************************************************************************
!
!! RVEC_SORT_INSERT_INDEX_D descending index sorts a real vector using insertion.
!
!
!  Reference:
!
!    Algorithm 1.1,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Modified:
!
!    07 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items in the vector.
!    N must be positive.
!
!    Input, real A(N), the array to be sorted.
!
!    Output, integer INDX(N), the sorted indices.  The array is sorted
!    when listed from A(INDX(1)) through A(INDX(N)).
!
  integer n
!
  real a(n)
  integer i
  integer indx(n)
  integer j
  real x
!
  call ivec_identity ( n, indx )

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( j >= 1 )

      if ( a(indx(j)) >= x ) then
        exit
      end if

      indx(j+1) = indx(j)
      j = j - 1

    end do

    indx(j+1) = i

  end do

  return
end
subroutine rvec_sort_quick_a ( n, a )
!
!*******************************************************************************
!
!! RVEC_SORT_QUICK_A ascending sorts a real vector using quick sort.
!
!
!  Example:
!
!    Input:
!
!      N = 7
!      A = ( 6, 7, 3, 2, 9, 1, 8 )
!
!    Output:
!
!      A = ( 1, 2, 3, 6, 7, 8, 9 )
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
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  integer, parameter :: MAXLEVEL = 25
!
  integer n
!
  real a(n)
  integer base
  integer l_segment
  integer level
  integer n_segment
  integer rsave(MAXLEVEL)
  integer r_segment
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_SORT_QUICK_A - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call rvec_part_quick_a ( n_segment, a(base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( l_segment > 1 ) then

      if ( level > MAXLEVEL ) then
        write ( *, * ) ' '
        write ( *, * ) 'RVEC_SORT_QUICK_A - Fatal error!'
        write ( *, * ) '  Exceeding recursion maximum of ', MAXLEVEL
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( n_segment > 0 ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine rvec_sort_shell_a ( n, a )
!
!*******************************************************************************
!
!! RVEC_SORT_SHELL_A ascending sorts a real array using Shell's sort.
!
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, an array to be sorted.
!    On output, the sorted array.
!
  integer n
!
  real a(n)
  real asave
  integer i
  integer ifree
  integer inc
  integer ipow
  integer j
  integer k
  integer maxpow
!
  if ( n <= 1 ) then
    return
  end if
!
!  Determine the smallest MAXPOW so that
!    N <= ( 3**MAXPOW - 1 ) / 2
!
  maxpow = 1

  do while ( 2 * n + 1 > 3**maxpow )
    maxpow = maxpow + 1
  end do

  if ( maxpow > 1 ) then
    maxpow = maxpow - 1
  end if
!
!  Now sort groups of size ( 3**IPOW - 1 ) / 2.
!
  do ipow = maxpow, 1, -1

    inc = ( 3**ipow - 1 ) / 2
!
!  Sort the values with indices equal to K mod INC.
!
    do k = 1, inc
!
!  Insertion sort of the items with index
!  INC+K, 2*INC+K, 3*INC+K, ...
!
      do i = inc+k, n, inc

        asave = a(i)
        ifree = i
        j = i - inc

        do

          if ( j < 1 ) then
            exit
          end if

          if ( a(j) <= asave ) then
            exit
          end if

          ifree = j
          a(j+inc) = a(j)
          j = j - inc

        end do

        a(ifree) = asave

      end do

    end do

  end do

  return
end
subroutine rvec_split_sort ( n, a, split, i_lt, i_gt )
!
!*******************************************************************************
!
!! RVEC_SPLIT_SORT "splits" a sorted vector, given a splitting value.
!
!
!  Discussion:
!
!    Given SPLIT, the routine seeks indices I_LT and I_GT so that
!
!      A(I_LT) < SPLIT < A(I_GT),
!
!    and if there are intermediate index values between I_LT and
!    I_GT, then those entries of A are exactly equal to SPLIT.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer N, the number of entries in A.
!
!    Input, real A(N), a sorted array.
!
!    Input, real SPLIT, a value to which the entries in A are
!    to be compared.
!
!    Output, integer I_LT:
!    0 if no entries are less than SPLIT;
!    N if all entries are less than SPLIT;
!    otherwise, the index of the last entry in A less than SPLIT.
!
!    Output, integer I_GT:
!    1 if all entries are greater than SPLIT;
!    N+1 if no entries are greater than SPLIT;
!    otherwise the index of the first entry in A greater than SPLIT.
!
  integer n
!
  real a(n)
  integer hi
  integer i
  integer i_gt
  integer i_lt
  integer lo
  integer mid
  real split
!
  if ( a(1) > split ) then
    i_lt = 0
    i_gt = 1
    return
  end if

  if ( a(n) < split ) then
    i_lt = n
    i_gt = n + 1
    return
  end if

  lo = 1
  hi = n

  do

    if ( lo + 1 == hi ) then
      i_lt = lo
      exit
    end if

    mid = ( lo + hi ) / 2

    if ( a(mid) >= split ) then
      hi = mid
    else
      lo = mid
    end if

  end do

  do i = i_lt + 1, n
    if ( a(i) > split ) then
      i_gt = i
      return
    end if
  end do

  i_gt = n + 1

  return
end
subroutine rvec_split_unsort ( n, a, split, isplit )
!
!*******************************************************************************
!
!! RVEC_SPLIT_UNSORT "splits" an unsorted vector based on a splitting value.
!
!
!  Discussion:
!
!    If the vector is already sorted, it is simpler to do a binary search
!    on the data than to call this routine.
!
!    The vector is not assumed to be sorted before input, and is not
!    sorted during processing.  If sorting is not needed, then it is
!    more efficient to use this routine.
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, real A(N), the array to split.  On output,
!    all the entries of A that are less than or equal to SPLIT
!    are in A(1:ISPLIT).
!
!    Input, real SPLIT, the value used to split the vector.
!    It is not necessary that any value of A actually equal SPLIT.
!
!    Output, integer ISPLIT, indicates the position of the last
!    entry of the split vector that is less than or equal to SPLIT.
!
  integer n
!
  real a(n)
  integer i
  integer i1
  integer i2
  integer i3
  integer isplit
  integer j1
  integer j2
  integer j3
  real split
!
!  Partition the vector into A1, A2, A3, where
!    A1 = A(I1:J1) holds values <= SPLIT,
!    A2 = A(I2:J2) holds untested values,
!    A3 = A(I3:J3) holds values > SPLIT.
!
  i1 = 1
  j1 = 0

  i2 = 1
  j2 = n

  i3 = n+1
  j3 = n
!
!  Pick the next item from A2, and move it into A1 or A3.
!  Adjust indices appropriately.
!
  do i = 1, n

    if ( a(i2) <= split ) then
      i2 = i2 + 1
      j1 = j1 + 1
    else
      call i_swap ( a(i2), a(i3-1) )
      i3 = i3 - 1
      j2 = j2 - 1
    end if

  end do

  isplit = j1

  return
end
subroutine rvec_uniq ( n, a, nuniq )
!
!*******************************************************************************
!
!! RVEC_UNIQ keeps only the unique elements in a sorted real array.
!
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, real A(N).
!    On input, the sorted array of N elements;
!    On output, the sorted unique array of NUNIQ elements.
!
!    Output, integer NUNIQ, the number of unique elements of A.
!
  integer n
!
  real a(n)
  integer itest
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1


  do itest = 2, n

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if

  end do

  return
end
subroutine rvec_uniq_count ( n, a, nuniq )
!
!*******************************************************************************
!
!! RVEC_UNIQ_COUNT counts the unique elements in an unsorted real array.
!
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input, real A(N), the array to examine, which does NOT have to
!    be sorted.
!
!    Output, integer NUNIQ, the number of unique elements of A.
!
  integer n
!
  real a(n)
  integer i
  integer j
  integer nuniq
!
  nuniq = 0

  do i = 1, n

    nuniq = nuniq + 1

    do j = 1, i-1

      if ( a(i) == a(j) ) then
        nuniq = nuniq - 1
        exit
      end if

    end do

  end do

  return
end
subroutine rvec_uniq_hist ( n, a, maxuniq, nuniq, auniq, acount )
!
!*******************************************************************************
!
!! RVEC_UNIQ_HIST computes a histogram of the unique elements of a sorted vector.
!
!
!  Modified:
!
!    09 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input, real A(N), the array to examine.  The elements of A
!    should have been sorted.
!
!    Input, integer MAXUNIQ, the maximum number of unique elements
!    that can be handled.  If there are more than MAXUNIQ unique
!    elements in A, the excess will be ignored.
!
!    Output, integer NUNIQ, the number of unique elements of A.
!
!    Output, real AUNIQ(NUNIQ), the unique elements of A.
!
!    Output, integer ACOUNT(NUNIQ), the number of times each element
!    of AUNIQ occurs in A.
!
  integer maxuniq
  integer n
!
  real a(n)
  integer acount(maxuniq)
  real auniq(maxuniq)
  integer i
  integer nuniq
!
!  Start taking statistics.
!
  nuniq = 0

  do i = 1, n

    if ( i == 1 ) then

      nuniq = 1
      auniq(nuniq) = a(1)
      acount(nuniq) = 1

    else if ( a(i) == auniq(nuniq) ) then

      acount(nuniq) = acount(nuniq) + 1

    else if ( nuniq < maxuniq ) then

      nuniq = nuniq + 1
      auniq(nuniq) = a(i)
      acount(nuniq) = 1

    end if

  end do

  return
end
subroutine rvec_unit_euclidean ( n, a )
!
!*******************************************************************************
!
!! RVEC_UNIT_EUCLIDEAN Euclidean normalizes a vector in ND.
!
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input/output, A(N), the vector to be normalized.
!
  integer n
!
  real a(n)
  real norm
!
  norm = sum ( a(1:n)**2 )

  norm = sqrt ( norm )

  if ( norm /= 0.0E+00 ) then
    a(1:n) = a(1:n) / norm
  end if

  return
end
subroutine rvec_unit_euclidean_2d ( x, y )
!
!*******************************************************************************
!
!! RVEC_UNIT_EUCLIDEAN_2D Euclidean normalizes a vector in 2D.
!
!
!  Modified:
!
!    17 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, X, Y, the components of the vector.
!
  real norm
  real x
  real y
!
  norm = sqrt ( x * x + y * y )

  if ( norm /= 0.0E+00 ) then
    x = x / norm
    y = y / norm
  end if

  return
end
subroutine rvec_unit_euclidean_3d ( x, y, z )
!
!*******************************************************************************
!
!! RVEC_UNIT_EUCLIDEAN_3D Euclidean normalizes a vector in 3D.
!
!
!  Modified:
!
!    17 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, X, Y, Z, the components of the vector.
!
  real norm
  real x
  real y
  real z
!
  norm = sqrt ( x * x + y * y + z * z )

  if ( norm /= 0.0E+00 ) then
    x = x / norm
    y = y / norm
    z = z / norm
  end if

  return
end
subroutine rvec_unit_sum ( n, a )
!
!*******************************************************************************
!
!! RVEC_UNIT_SUM normalizes a real vector to have unit sum.
!
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
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, real A(N), the vector to be normalized.  On output,
!    the entries of A should have unit sum.  However, if the input vector
!    has zero sum, the routine halts.
!
  integer n
!
  real a(n)
  real a_sum
!
  a_sum = sum ( a(1:n) )

  if ( a_sum == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_UNIT_SUM - Fatal error!'
    write ( *, * ) '  The vector entries sum to 0.'
    stop
  end if

  a(1:n) = a(1:n) / a_sum

  return
end
subroutine rvec_variance ( n, a, variance )
!
!*******************************************************************************
!
!! RVEC_VARIANCE returns the variance of a real vector.
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
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector whose variance is desired.
!
!    Output, real VARIANCE, the variance of the vector entries.
!
  integer n
!
  real a(n)
  real mean
  real variance
!
  call rvec_mean ( n, a, mean )

  variance = sum ( ( a(1:n) - mean )**2 )

  if ( n > 1 ) then
    variance = variance / real ( n - 1 )
  else
    variance = 0.0E+00
  end if

  return
end
subroutine s_cat ( s1, s2, s3 )
!
!*******************************************************************************
!
!! S_CAT concatenates two strings to make a third string.
!
!
!  Modified:
!
!    11 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3
!
  s3 = trim ( s1 ) // trim ( s2 )

  return
end
function s_eqi ( strng1, strng2 )
!
!*******************************************************************************
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!
!  Example:
!
!    STRING_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  integer i
  integer len1
  integer len2
  integer lenc
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2
  logical s_eqi
!
  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
function s_gti ( strng1, strng2 )
!
!*******************************************************************************
!
!! S_GTI = STRNG1 is lexically greater than STRNG2.
!
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_GTI, the result of the comparison.
!
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_gti
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2
!
  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( lgt ( s1, s2 ) ) then
      s_gti = .true.
      return
    else if ( llt ( s1, s2 ) ) then
      s_gti = .false.
      return
    end if

  end do

  if ( len1 <= len2 ) then
    s_gti = .false.
  else
    s_gti = .true.
  end if

  return
end
function s_lei ( strng1, strng2 )
!
!*******************************************************************************
!
!! S_LEI = STRNG1 is lexically less than or equal to STRNG2.
!
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_LEI, the result of the comparison.
!
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_lei
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2
!
  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( llt ( s1, s2 ) ) then
      s_lei = .true.
      return
    else if ( lgt ( s1, s2 ) ) then
      s_lei = .false.
      return
    end if

  end do

  if ( len1 <= len2 ) then
    s_lei = .true.
  else
    s_lei = .false.
  end if

  return
end
function s_lti ( strng1, strng2 )
!
!*******************************************************************************
!
!! S_LTI = STRNG1 is lexically less than STRNG2.
!
!
!  Discussion:
!
!    The comparison is done in a case-insensitive way.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_LTI, the result of the comparison.
!
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_lti
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2
!
  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( llt ( s1, s2 ) ) then
      s_lti = .true.
      return
    else if ( lgt ( s1, s2 ) ) then
      s_lti = .false.
      return
    end if

  end do

  if ( len1 < len2 ) then
    s_lti = .true.
  else
    s_lti = .false.
  end if

  return
end
subroutine s_memory ( action, name, cval )
!
!*******************************************************************************
!
!! S_MEMORY manages a set of runtime string variables.
!
!
!  Discussion:
!
!    S_MEMORY allows the user to define the name of a string variable,
!    set it, push a new value onto a stack of up to
!    five values, pop a value off the stack, or get the current value.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, desired action.
!
!    'APPEND',  append CVAL to NAME.
!    'GET',     return value of NAME in CVAL.
!    'INIT',    reset all values to zero, wipe out all names.
!    'NAME',    add a variable of the given name.
!    'POP',     pop the stack, retrieving the last pushed value.
!    'PREPEND', prepend CVAL to NAME.
!    'PRINT',   print the value of NAME, and return in CVAL.
!               NAME = '*' prints all variables, and returns ' ' in CVAL.
!    'PUSH',    push the previous value down the stack, insert a new one.
!    'SET',     set variable NAME to CVAL.
!
!    Input, character ( len = * ) NAME, the name of the variable.
!    NAME may be blank for the 'INIT' command, but never any other time.
!
!    Input/output, integer CVAL.
!
!    For the 'APPEND', 'NAME', 'PREPEND', 'PUSH', and 'SET' commands,
!    CVAL must contain the appending value, initial value, prepending
!    value, push value, or set value of the variable on input.
!
!    For the 'GET' and 'POP' commands, CVAL will contain the current
!    value or popped value of the named variable on output.
!
  integer, parameter :: maxnam = 100
  integer, parameter :: maxcol = 5
!
  character ( len = * ) action
  character ( len = * ) cval
  character ( len = 80 ), save, dimension ( maxnam, maxcol ) :: cvals
  integer i
  integer, save, dimension ( maxnam ) :: ipoint = (/ ( 0, i = 1, maxnam ) /)
  integer j
  character ( len = * ) name
  character ( len = 20 ), save, dimension ( maxnam ) :: names = &
    (/ ( ' ', i = 1, maxnam ) /)
  integer, save :: numnam = 0
  logical s_eqi
!
  data ( ( cvals(i,j), i = 1, maxnam ), j = 1, maxcol ) / 500 * ' ' /
!
  if ( name == ' ' .and. .not. s_eqi ( action, 'INIT' ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  The value of NAME cannot be blank!'
    stop
  end if
!
!  APPEND: Append interesting part of CVAL to NAME:
!
  if ( s_eqi (action, 'APPEND' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        call s_cat ( cvals(i,ipoint(i)), cval, cvals(i,ipoint(i)) )
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to append to unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  GET: Get the current value of a variable.
!
  else if ( s_eqi ( action, 'GET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        cval = cvals(i,ipoint(i))
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to get value of unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  INIT: Initialize everything.
!
  else if ( s_eqi ( action, 'INIT' ) ) then

    numnam = 0
    cvals(1:maxnam,1:maxcol) = ' '
    names(1:maxnam) = ' '
    ipoint(1:maxnam) = 0
!
!  NAME: Declare the name of something and set it.
!
  else if ( s_eqi ( action, 'NAME' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        write ( *, * ) ' '
        write ( *, * ) 'S_MEMORY - Warning!'
        write ( *, '(''  There is ALREADY a variable '', a )' ) trim ( name )
        write ( *, * ) '  The new value has been stored.'
        ipoint(i) = 1
        cvals(i,ipoint(i)) = cval
        return
      end if

    end do

    if ( numnam < maxnam ) then
      numnam = numnam + 1
      i = numnam
      names(i) = name
      ipoint(i) = 1
      cvals(i,ipoint(i)) = cval
    else
      write ( *, * ) ' '
      write ( *, * ) 'S_MEMORY - Fatal error!'
      write ( *, * ) '  We have reached the name limit of ', maxnam
      stop
    end if
!
!  POP: "Pop" a value, decrement pointer.
!
  else if ( s_eqi ( action, 'POP' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) > 1 ) then
          ipoint(i) = ipoint(i) - 1
          cval = cvals(i,ipoint(i))
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'S_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to pop the stack to 0.'
          write ( *, '(''  Variable name is '',a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to pop an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  PREPEND: Prepend interesting part of CVAL to NAME:
!
  else if ( s_eqi ( action, 'PREPEND' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        call s_cat ( cval, cvals(i,ipoint(i)), cvals(i,ipoint(i)) )
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to prepend to unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  PRINT: "Print" the value, and return in CVAL.
!
  else if ( s_eqi ( action, 'PRINT' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) .or. name == '*' ) then

        cval = cvals(i,ipoint(i))
        write ( *, '(a,a,a,a)' ) &
          'C_MEMORY - Value of ', trim ( names(i) ), ' is '
        write ( *, '(a)' ) trim ( cval )

      end if

      if ( s_eqi ( name, names(i) ) ) then
        return
      end if

    end do

    if ( s_eqi ( name, '*' ) ) then
      cval = ' '
      return
    end if

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to print an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  PUSH: "Push" a value, increment the pointer.
!
  else if ( s_eqi ( action, 'PUSH' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        if ( ipoint(i) < maxcol ) then
          ipoint(i) = ipoint(i) + 1
          cvals(i,ipoint(i)) = cval
          return
        else
          write ( *, * ) ' '
          write ( *, * ) 'S_MEMORY - Fatal error!'
          write ( *, * ) '  Attempt to push the stack past ', maxcol
          write ( *, '(''  Variable name is '',a)' ) trim ( name )
          stop
        end if
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to push an unknown variable.'
    write ( *, '(''  Variable name is '',a)' ) trim ( name )
    stop
!
!  SET: Set something.
!
  else if ( s_eqi ( action, 'SET' ) ) then

    do i = 1, numnam

      if ( s_eqi ( name, names(i) ) ) then
        cvals(i,ipoint(i)) = cval
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write ( *, '(''  Variable name is '', a)' ) trim ( name )
    stop
!
!  Unrecognized action.
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'S_MEMORY - Fatal error!'
    write ( *, * ) '  Unrecognized action:'
    write ( *, '(a)' ) trim ( action )
    stop

  end if

  return
end
subroutine s_swap ( s1, s2 )
!
!*******************************************************************************
!
!! S_SWAP swaps two strings.
!
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S1, S2.  On output, the values of S1 and
!    S2 have been interchanged.
!
  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = 256 ) s3
!
  s3 = s1
  s1 = s2
  s2 = s3

  return
end
subroutine sge_hess ( fx, lda, n, x, h )
!
!*******************************************************************************
!
!! SGE_HESS approximates a Hessian matrix via finite differences.
!
!
!  Discussion:
!
!    H(I,J) = d2 F / d X(I) d X(J)
!
!    The values returned by this routine will be only approximate.
!    In some cases, they will be so poor that they are useless.
!    However, one of the best applications of this routine is for
!    checking your own Hessian calculations, since as Heraclitus
!    said, you'll never get the same result twice when you differentiate
!    a complicated expression by hand.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FX, the name of the user function routine.
!    FX should have the form:
!
!      subroutine fx ( n, x, f )
!      integer n
!      real f
!      real x(n)
!
!    Input, integer LDA, the leading dimension of H.
!
!    Input, integer N, the number of variables.
!
!    Input, real X(N), the values of the variables.
!
!    Output, real H(LDA,N), the approximated N by N Hessian matrix.
!
  integer lda
  integer n
!
  real eps
  real f00
  real fmm
  real fmp
  real fpm
  real fpp
  real h(lda,n)
  integer i
  integer j
  real s(n)
  real x(n)
  real xi
  real xj
!
  external fx
!
!  Choose the stepsizes.
!
  eps = ( epsilon ( eps ) )**0.33E+00

  do i = 1, n
    s(i) = eps * max ( abs ( x(i) ), 1.0E+00 )
  end do
!
!  Calculate the diagonal elements.
!
  do i = 1, n

    xi = x(i)

    call fx ( n, x, f00 )

    x(i) = xi + s(i)
    call fx ( n, x, fpp )

    x(i) = xi - s(i)
    call fx ( n, x, fmm )

    h(i,i) = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s(i)**2

    x(i) = xi

  end do
!
!  Calculate the off diagonal elements.
!
  do i = 1, n

    xi = x(i)

    do j = i+1, n

      xj = x(j)

      x(i) = xi + s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fpp )

      x(i) = xi + s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fpm )

      x(i) = xi - s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fmp )

      x(i) = xi - s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fmm )

      h(j,i) = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0E+00 * s(i) * s(j) )

      h(i,j) = h(j,i)

      x(j) = xj

    end do

    x(i) = xi

  end do

  return
end
subroutine sge_jac ( eps, fprime, fx, lda, m, n, x )
!
!*******************************************************************************
!
!! SGE_JAC estimates a dense jacobian matrix of the function FX.
!
!
!  Discussion:
!
!    FPRIME(I,J) = d F(I) / d X(J).
!
!    The jacobian is assumed to be dense, and the LINPACK/LAPACK
!    real general matrix storage mode ("SGE") is used.
!
!    Forward differences are used, requiring N+1 function evaluations.
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real EPS, a tolerance to be used for shifting the
!    X values during the finite differencing.  No single value
!    of EPS will be reliable for all vectors X and functions FX.
!
!    Values of EPS have typically been chosen between
!    sqrt ( EPSMCH ) and sqrt ( sqrt ( EPSMCH ) ) where EPSMCH is the
!    machine tolerance.
!
!    If EPS is too small, then F(X+EPS) will be the same as
!    F(X), and the jacobian will be full of zero entries.
!
!    If EPS is too large, the finite difference estimate will
!    be inaccurate.
!
!    Output, real FPRIME(LDA,N), the M by N estimated jacobian matrix.
!
!    Input, external FX, the name of the user written
!    routine which evaluates the function at a given point X.
!
!    FX should have the form:
!
!      subroutine fx ( m, n, x, f )
!      integer m
!      integer n
!      real f(m)
!      real x(n)
!      f(1:m) = ...
!      return
!      end
!
!    Input, integer LDA, the leading dimension of FPRIME.
!    LDA must be at least M.
!
!    Input, integer M, the number of functions.
!
!    Input, integer N, the number of variables.
!
!    Input, real X(N), the point where the jacobian is to be estimated.
!
  integer lda
  integer m
  integer n
!
  real del
  real eps
  real fprime(lda,n)
  integer i
  integer j
  real x(n)
  real xsave
  real work1(m)
  real work2(m)
!
  external fx
!
!  Evaluate the function at the base point, X.
!
  call fx ( m, n, x, work2 )
!
!  Now, one by one, vary each component J of the base point X, and
!  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
!
  do j = 1, n

    xsave = x(j)
    del = eps * ( 1.0E+00 + abs ( x(j) ) )
    x(j) = x(j) + del
    call fx ( m, n, x, work1 )
    x(j) = xsave
    fprime(1:m,j) = ( work1(1:m) - work2(1:m) ) / del

  end do

  return
end
subroutine sge_solve ( a, b, lda, n, x, ierror )
!
!*******************************************************************************
!
!! SGE_SOLVE computes the solution of an N by N linear system.
!
!
!  Discussion:
!
!    The linear system may be represented as
!
!      A*X = B
!
!    If the linear system is singular, but consistent, then the routine will
!    still produce a solution.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A(LDA,N).
!    On input, A is the coefficient matrix to be inverted.
!    On output, A has been overwritten.
!
!    Input/output, real B(N).
!    On input, B is the right hand side of the system.
!    On output, B has been overwritten.
!
!    Input, integer LDA, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer N, the number of equations.
!
!    Output, real X(N), the solution of the linear system.
!
!    Output, integer IERROR.
!    0, no error detected.
!    1, consistent singularity.
!    2, inconsistent singularity.
!    3, the internal array is too small (FATAL).
!
  integer, parameter :: maxn = 100
!
  integer lda
  integer n
!
  real a(lda,n)
  real amax
  real b(n)
  integer i
  integer ierror
  integer imax
  integer ipiv(maxn)
  integer j
  integer k
  real x(n)
!
  ierror = 0
!
  if ( n > maxn ) then
    ierror = 3
    write ( *, * ) ' '
    write ( *, * ) 'SGE_SOLVE - Fatal error!'
    write ( *, * ) '  Number of equations N = ', n
    write ( *, * ) '  exceeds internal limit MAXN = ', maxn
    stop
  end if

  ipiv(1:n) = 0
  x(1:n) = 0.0E+00
!
!  Process the matrix.
!
  do k = 1, n
!
!  In column K:
!    Seek the row IMAX with the properties that:
!  IMAX has not already been used as a pivot;
!  A(IMAX,K) is larger in magnitude than any other candidate.
!
    amax = 0.0E+00
    imax = 0
    do i = 1, n
      if ( ipiv(i) == 0 ) then
        if ( abs ( a(i,k) ) > amax ) then
          imax = i
          amax = abs ( a(i,k) )
        end if
      end if
    end do
!
!  If you found a pivot row IMAX, then,
!    eliminate the K-th entry in all rows that have not been used for pivoting.
!
    if ( imax /= 0 ) then

      ipiv(imax) = k
      a(imax,k+1:n) = a(imax,k+1:n) / a(imax,k)
      b(imax) = b(imax) / a(imax,k)
      a(imax,k) = 1.0E+00

      do i = 1, n

        if ( ipiv(i) == 0 ) then
          a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(imax,k+1:n)
          b(i) = b(i) - a(i,k) * b(imax)
          a(i,k) = 0.0E+00
        end if

      end do

    end if

  end do
!
!  Now, every row with nonzero IPIV begins with a 1, and
!  all other rows are all zero.  Begin solution.
!
  do j = n, 1, -1

    imax = 0
    do k = 1, n
      if ( ipiv(k) == j ) then
        imax = k
      end if
    end do

    if ( imax == 0 ) then

      x(j) = 0.0E+00

      if ( b(j) == 0.0E+00 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'SGE_SOLVE - Warning:'
        write ( *, * ) '  Consistent singularity, equation = ', j
      else
        ierror = 2
        write ( *, * ) ' '
        write ( *, * ) 'SGE_SOLVE - Error:'
        write ( *, * ) '  Inconsistent singularity, equation = ', j
      end if

    else

      x(j) = b(imax)

      do i = 1, n
        if ( i /= imax ) then
          b(i) = b(i) - a(i,j) * x(j)
        end if
      end do

    end if

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
subroutine svec_merge_a ( na, a, nb, b, nc, c )
!
!*******************************************************************************
!
!! SVEC_MERGE_A merges two ascending sorted string arrays.
!
!
!  Discussion:
!
!    The elements of A and B should be sorted in ascending order.
!
!    The elements in the output array C will be in ascending order, and unique.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, character ( len = * ) A(NA), the first sorted array.
!
!    Input, integer NB, the dimension of B.
!
!    Input, character ( len = * ) B(NB), the second sorted array.
!
!    Output, integer NC, the number of elements in the output array.
!    Note that C should usually be dimensioned at least NA+NB in the
!    calling routine.
!
!    Output, character ( len = * ) C(NC), the merged unique sorted array.
!
  integer na
  integer nb
!
  character ( len = * ) a(na)
  character ( len = * ) b(nb)
  character ( len = * ) c(na+nb)
  integer j
  integer ja
  integer jb
  integer na2
  integer nb2
  integer nc
!
  na2 = na
  nb2 = nb
!
  ja = 0
  jb = 0
  nc = 0

  do
!
!  If we've used up all the entries of A, stick the rest of B on the end.
!
    if ( ja >= na2 ) then

      do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 ) then
          nc = nc + 1
          c(nc) = b(jb)
        else if ( llt ( c(nc), b(jb) ) ) then
          nc = nc + 1
          c(nc) = b(jb)
        end if
      end do

      exit
!
!  If we've used up all the entries of B, stick the rest of A on the end.
!
    else if ( jb >= nb2 ) then

      do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 ) then
          nc = nc + 1
          c(nc) = a(ja)
        else if ( llt ( c(nc), a(ja) ) ) then
          nc = nc + 1
          c(nc) = a(ja)
        end if
      end do

      exit
!
!  Otherwise, if the next entry of A is smaller, that's our candidate.
!
    else if ( lle ( a(ja+1), b(jb+1) ) ) then

      ja = ja + 1
      if ( nc == 0 ) then
        nc = nc + 1
        c(nc) = a(ja)
      else if ( llt ( c(nc), a(ja) ) ) then
        nc = nc + 1
        c(nc) = a(ja)
      end if
!
!  ...or if the next entry of B is the smaller, consider that.
!
    else

      jb = jb + 1
      if ( nc == 0 ) then
        nc = nc + 1
        c(nc) = b(jb)
      else if ( llt ( c(nc), b(jb) ) ) then
        nc = nc + 1
        c(nc) = b(jb)
      end if
    end if

  end do

  return
end
subroutine svec_permute ( n, a, p )
!
!*******************************************************************************
!
!! SVEC_PERMUTE permutes a string vector in place.
!
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (  3,     2,     4,       2,      1 )
!      A = ( 'ONE', 'TWO', 'THREE', 'FOUR', 'FIVE' )
!
!    Output:
!
!      A    = ( 'FIVE', 'FOUR', 'ONE', 'THREE', 'TWO' ).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input/output, character ( len = * ) A(N), the array to be permuted.
!
!    Input, integer P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  integer n
!
  character ( len = * ) a(n)
  character ( len = 256 ) a_temp
  integer i
  integer ierror
  integer iget
  integer iput
  integer istart
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SVEC_PERMUTE - Fatal error!'
    write ( *, * ) '  The input array does not represent'
    write ( *, * ) '  a proper permutation.  In particular, the'
    write ( *, * ) '  array is missing the value ', ierror
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. iget > n ) then
          write ( *, * ) ' '
          write ( *, * ) 'CVEC_PERMUTE - Fatal error!'
          write ( *, * ) '  An array element had an illegal value.'
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)

  return
end
subroutine svec_reverse ( n, a )
!
!*******************************************************************************
!
!! SVEC_REVERSE reverses the elements of a string vector.
!
!
!  Example:
!
!    Input:
!
!      N = 4,
!      A = ( 'Bob', 'Carol', 'Ted', 'Alice' ).
!
!    Output:
!
!      A = ( 'Alice', 'Ted', 'Carol', 'Bob' ).
!
!  Modified:
!
!    28 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, character ( len = * ) A(N), the array to be reversed.
!
  integer n
!
  character ( len = * ) a(n)
  character ( len = 256 ) a_temp
  integer i
!
  do i = 1, n/2
    a_temp = a(i)
    a(i) = a(n+1-i)
    a(n+1-i) = a_temp
  end do

  return
end
subroutine svec_search_binary_a ( n, a, b, indx )
!
!*******************************************************************************
!
!! SVEC_SEARCH_BINARY_A searches an ascending sorted string vector.
!
!
!  Discussion:
!
!    Binary search is used.
!
!  Reference:
!
!    Algorithm 1.9,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vector.
!
!    Input, character ( len = * ) A(N), the array to be searched.  A must
!    be sorted in increasing order.
!
!    Input, character ( len = * ) B, the value to be searched for.
!
!    Output, integer INDX, the result of the search.
!    0, B does not occur in A.
!    I, A(I) = B.
!
  integer n
!
  character ( len = * ) a(n)
  character ( len = * ) b
  integer high
  integer indx
  integer low
  integer mid
!
  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( llt ( a(mid), b ) ) then
      low = mid + 1
    else if ( lgt ( a(mid), b ) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine svec_sort_heap_a ( n, a )
!
!*******************************************************************************
!
!! SVEC_SORT_HEAP_A ascending sorts a vector of character strings using heap sort.
!
!
!  Discussion:
!
!    The ASCII collating sequence is used.  This means
!      A < B < C < .... < Y < Z < a < b < .... < z.
!    Numbers and other symbols may also occur, and will be sorted according to
!    the ASCII ordering.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of strings
!
!    Input/output, character ( len = * ) A(N);
!    On input, an array of strings to be sorted.
!    On output, the sorted array.
!
  integer n
!
  character ( len = * ) a(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Do the sorting using the external heap sort routine.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( indx > 0 ) then

      call s_swap ( a(i), a(j) )

    else if ( indx < 0 ) then

      if ( lle ( a(i), a(j) ) ) then
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
subroutine svec_uniq ( n, a, nuniq )
!
!*******************************************************************************
!
!! SVEC_UNIQ finds the number of unique entries in a vector of strings.
!
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the array.
!
!    Input/output, character ( len = * ) A(N).
!    On input, the sorted list of strings.
!    On output, the unique elements, in sorted order.
!
!    Output, integer NUNIQ, the number of unique elements in the array.
!
  integer n
!
  character ( len = * ) a(n)
  integer itest
  integer nuniq
!
  nuniq = 0
  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if

  end do

  return
end
subroutine sveci_search_binary_a ( n, a, b, indx )
!
!*******************************************************************************
!
!! SVECI_SEARCH_BINARY_A searches an ascending sorted vector of implicitly capitalized strings.
!
!
!  Reference:
!
!    Algorithm 1.9,
!    Donald Kreher and Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vector.
!
!    Input, character ( len = * ) A(N), the array to be searched.  A must
!    be sorted in increasing order.
!
!    Input, character ( len = * ) B, the value to be searched for.
!
!    Output, integer INDX, the result of the search.
!    0, B does not occur in A.
!    I, A(I) = B, ignoring capitalization.
!
  integer n
!
  character ( len = * ) a(n)
  character ( len = * ) b
  integer high
  integer indx
  integer low
  integer mid
  logical s_eqi
  logical s_gti
  logical s_lti
!
  indx = 0

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( s_eqi ( a(mid), b ) ) then
      indx = mid
      exit
    else if ( s_lti ( a(mid), b ) ) then
      low = mid + 1
    else if ( s_gti ( a(mid), b ) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine sveci_sort_heap_a ( n, a )
!
!*******************************************************************************
!
!! SVECI_SORT_HEAP_A ascending sorts a vector of implicitly capitalized strings using heap sort.
!
!
!  Discussion:
!
!    The characters in an implicitly capitalized string are treated as
!    though they had been capitalized.  Thus, the letters 'a' and 'A'
!    are considered equal, both 'a' and 'A' precede 'B', and
!    'Fox' and 'fOx' are considered equal.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of strings
!
!    Input/output, character ( len = * ) A(N);
!    On input, an array of strings to be sorted.
!    On output, the sorted array.
!
  integer n
!
  character ( len = * ) a(n)
  integer i
  integer indx
  integer isgn
  integer j
  logical s_lei
!
!  Do the sorting using the external heap sort routine.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( indx > 0 ) then

      call s_swap ( a(i), a(j) )

    else if ( indx < 0 ) then

      if ( s_lei ( a(i), a(j) ) ) then
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
function tand ( angle )
!
!*******************************************************************************
!
!! TAND returns the tangent of an angle given in degrees.
!
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle, in degrees.
!
!    Output, real TAND, the tangent of the angle.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
  real, parameter :: degrees_to_radians = PI / 180.0E+00
!
  real angle
  real tand
!
  tand  =  sin ( degrees_to_radians * angle ) &
         / cos ( degrees_to_radians * angle )

  return
end
subroutine tuple_next2 ( n, xmin, xmax, x, rank )
!
!*******************************************************************************
!
!! TUPLE_NEXT2 computes the next element of an integer tuple space.
!
!
!  Discussion:
!
!    The elements X are N vectors.
!
!    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
!
!    The elements are produced one at a time.
!
!    The first element is
!      (XMIN(1), XMIN(2), ..., XMIN(N)),
!    the second is (probably)
!      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
!    and the last element is
!      (XMAX(1), XMAX(2), ..., XMAX(N))
!
!    Intermediate elements are produced in a lexicographic order, with
!    the first index more important than the last, and the ordering of
!    values at a fixed index implicitly defined by the sign of
!    XMAX(I) - XMIN(I).
!
!  Examples:
!
!    N = 2,
!    XMIN = (/ 1, 10 /)
!    XMAX = (/ 3,  8 /)
!
!    RANK    X
!    ----  -----
!      1   1 10
!      2   1  9
!      3   1  8
!      4   2 10
!      5   2  9
!      6   2  8
!      7   3 10
!      8   3  9
!      9   3  8
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
!    Input, integer N, the number of components.
!
!    Input, integer XMIN(N), XMAX(N), the "minimum" and "maximum" entry values.
!    These values are minimum and maximum only in the sense of the lexicographic
!    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
!    than XMAX(I).
!
!    Input/output, integer X(N), on input the previous tuple.
!    On output, the next tuple.
!
!    Input/output, integer RANK, the rank of the item.  On first call,
!    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
!    there are no more items in the sequence.
!
  integer n
!
  integer i
  integer rank
  integer x(n)
  integer xmin(n)
  integer xmax(n)
!
  if ( rank < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, * ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank > product ( 1 + abs ( xmax(1:n) - xmin(1:n) ) ) ) then
    write ( *, * ) ' '
    write ( *, * ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, * ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank == 0 ) then
    x(1:n) = xmin(1:n)
    rank = 1
    return
  end if

  rank = rank + 1
  i = n

  do

    if ( x(i) /= xmax(i) ) then
      x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
      exit
    end if

    x(i) = xmin(i)

    if ( i == 1 ) then
      rank = 0
      exit
    end if

    i = i - 1

  end do

  return
end
subroutine tvec_even ( nt, t )
!
!*******************************************************************************
!
!! TVEC_EVEN computes an evenly spaced set of angles between 0 and 2*PI.
!
!
!  Discussion:
!
!    The computation realizes that 0 = 2 * PI.
!
!  Example:
!
!    NT = 4
!
!    T = ( 0, PI/2, PI, 3*PI/2 )
!
!  Modified:
!
!    14 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NT, the number of values to compute.
!
!    Output, real TVEC(NT), the evenly spaced angles, in radians.
!
  integer nt
!
  integer i
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real t(nt)
!
  do i = 1, nt
    t(i) = real ( 2 * ( i - 1 ) ) * pi / real ( nt )
  end do

  return
end
subroutine tvec_even_bracket ( theta1, theta2, nt, t )
!
!*******************************************************************************
!
!! TVEC_EVEN_BRACKET computes an evenly spaced set of angles between THETA1 and THETA2.
!
!
!  Example:
!
!    NT = 4
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 30, 50, 70, 90 )
!
!  Modified:
!
!    13 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision THETA1, THETA2, the limiting angles.
!
!    Input, integer NT, the number of values to compute.
!
!    Output, double precision TVEC(NT), the evenly spaced angles.
!
  integer nt
!
  integer i
  double precision t(nt)
  double precision theta1
  double precision theta2
!
  if ( nt == 1 ) then
    t(i) = ( theta1 + theta2 ) / 2.0D+00
  else
    do i = 1, nt
      t(i) = dble ( dble ( nt - i ) * theta1 &
                  + dble ( i - 1 ) * theta2 ) / dble ( nt - 1 )
    end do
  end if

  return
end
subroutine tvec_even_bracket2 ( theta1, theta2, nt, t )
!
!*******************************************************************************
!
!! TVEC_EVEN_BRACKET2 computes an evenly spaced set of angles between THETA1 and THETA2.
!
!
!  Example:
!
!    NT = 5
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 40, 50, 60, 70, 80 )
!
!  Modified:
!
!    13 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision THETA1, THETA2, the limiting angles.
!
!    Input, integer NT, the number of values to compute.
!
!    Output, double precision TVEC(NT), the evenly spaced angles.
!
  integer nt
!
  integer i
  double precision t(nt)
  double precision theta1
  double precision theta2
!
  do i = 1, nt
    t(i) = dble ( dble ( nt + 1 - i ) * theta1 &
                  + dble ( i ) * theta2 ) / dble ( nt + 1 )
  end do

  return
end
subroutine tvec_even_bracket3 ( theta1, theta2, nt, t )
!
!*******************************************************************************
!
!! TVEC_EVEN_BRACKET3 computes an evenly spaced set of angles between THETA1 and THETA2.
!
!
!  Example:
!
!    NT = 3
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 40, 60, 80 )
!
!  Modified:
!
!    13 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision THETA1, THETA2, the limiting angles.
!
!    Input, integer NT, the number of values to compute.
!
!    Output, double precision TVEC(NT), the evenly spaced angles.
!
  integer nt
!
  integer i
  double precision t(nt)
  double precision theta1
  double precision theta2
!
  do i = 1, nt
    t(i) = dble ( dble ( 2 * nt + 1 - 2 * i ) * theta1 &
                  + dble ( 2 * i - 1 ) * theta2 ) / dble ( 2 * nt )
  end do

  return
end
function uniform_01_sample ( iseed )
!
!*******************************************************************************
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!
!  Formula:
!
!    ISEED = ISEED * (7**5) mod (2**31 - 1)
!    RANDOM = ISEED * / ( 2**31 - 1 )
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input/output, integer ISEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real UNIFORM_01_SAMPLE, a random value between 0 and 1.
!

!
!  IA = 7**5
!  IB = 2**15
!  IB16 = 2**16
!  IP = 2**31-1
!
  integer, parameter :: ia = 16807
  integer, parameter :: ib15 = 32768
  integer, parameter :: ib16 = 65536
  integer, parameter :: ip = 2147483647
!
  integer iprhi
  integer iseed
  integer ixhi
  integer k
  integer leftlo
  integer loxa
  real uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real ( iseed ) * 4.656612875E-10

  return
end
