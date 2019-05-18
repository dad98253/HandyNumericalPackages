!  prob.f90  11 July 2000
!
subroutine anglit_cdf ( x, cdf )
!
!*******************************************************************************
!
!! ANGLIT_CDF evaluates the Anglit CDF.
!
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real cdf
  real x
!
  if ( x <  - 0.25 * PI ) then
    cdf = 0.0
  else if ( x < 0.25 * PI ) then
    cdf = 0.5 - 0.5 * cos ( 2.0 * x + PI / 2.0 )
  else
    cdf = 1.0
  end if

  return
end
subroutine anglit_cdf_inv ( cdf, x )
!
!*******************************************************************************
!
!! ANGLIT_CDF_INV inverts the Anglit CDF.
!
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, real X, the corresponding argument.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ANGLIT_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = 0.5 * ( acos ( 1.0 - 2.0 * cdf ) - PI / 2.0 )

  return
end
subroutine anglit_mean ( mean )
!
!*******************************************************************************
!
!! ANGLIT_MEAN returns the mean of the Anglit PDF.
!
! 
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MEAN, the mean of the PDF.
!
  real mean
!
  mean = 0.0

  return
end
subroutine anglit_pdf ( x, pdf )
!
!*******************************************************************************
!
!! ANGLIT_PDF evaluates the Anglit PDF.
!
!
!  Formula:
!
!    PDF(X) = SIN ( 2 * X + PI / 2 ) for -PI/4 <= X <= PI/4
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real pdf
  real x
!
  if ( x <= - 0.25 * PI .or. x >= 0.25 * PI ) then
    pdf = 0.0
  else
    pdf = sin ( 2.0 * x + 0.25 * PI )
  end if

  return
end
subroutine anglit_sample ( iseed, x )
!
!*******************************************************************************
!
!! ANGLIT_SAMPLE samples the Anglit PDF.
!
!
!  Modified:
!
!    28 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call anglit_cdf_inv ( cdf, x )

  return
end
subroutine anglit_variance ( variance )
!
!*******************************************************************************
!
!! ANGLIT_VARIANCE returns the variance of the Anglit PDF.
!
! 
!  Discussion:
!
!    Variance = 
!      Integral ( -PI/4 <= X <= PI/4 ) X**2 * SIN ( 2 * X + PI / 2 ) 
!
!    Antiderivative = 
!      0.5 * X * SIN ( 2 * X + PI / 2 )
!      + ( 0.25 - 0.5 * X**2 ) * COS ( 2 * X + PI / 2 )
!
!  Modified:
!
!    29 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real variance
!
  variance = 0.0625 * PI**2 - 0.5

  return
end
subroutine arcsin_cdf ( x, cdf )
!
!*******************************************************************************
!
!! ARCSIN_CDF evaluates the Arcsin CDF.
!
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real cdf
  real x
!
  if ( x < 0.0 ) then
    cdf = 0.0
  else if ( x < 1.0 ) then
    cdf = 2.0 * asin ( x ) / PI
  else
    cdf = 1.0
  end if

  return
end
subroutine arcsin_cdf_inv ( cdf, x )
!
!*******************************************************************************
!
!! ARCSIN_CDF_INV inverts the Arcsin CDF.
!
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, real X, the corresponding argument.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ARCSIN_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = sin ( PI * cdf / 2.0 )

  return
end
subroutine arcsin_mean ( mean )
!
!*******************************************************************************
!
!! ARCSIN_MEAN returns the mean of the Arcsin PDF.
!
! 
!  Modified:
!
!    17 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MEAN, the mean of the PDF.
!
  real mean
!
  mean = 0.5

  return
end
subroutine arcsin_pdf ( x, pdf )
!
!*******************************************************************************
!
!! ARCSIN_PDF evaluates the Arcsin PDF.
!
!
!  Discussion:
!
!    The LOGISTIC EQUATION has the form:
!
!      X(N+1) = 4.0 * LAMBDA * ( 1.0 - X(N) ).
!
!    where 0 < LAMBDA <= 1.  This nonlinear difference equation maps
!    the unit interval into itself, and is a simple example of a system
!    exhibiting chaotic behavior.  Ulam and von Neumann studied the
!    logistic equation with LAMBDA = 1, and showed that iterates of the
!    function generated a sequence of pseudorandom numbers with 
!    the Arcsin probability density function.
!
!    The derived sequence
!
!      Y(N) = ( 2 / PI ) * Arcsin ( SQRT ( X(N) ) )
!
!    is a pseudorandom sequence with the uniform probability density
!    function on [0,1].  For certain starting values, such as X(0) = 0, 0.75,
!    or 1.0, the sequence degenerates into a constant sequence, and for
!    values very near these, the sequence takes a while before becoming
!    chaotic.
!
!  Formula:
!
!    PDF(X) = 1 / ( PI * Sqrt ( X * ( 1 - X ) ) )
!
!  Reference:
!
!    Daniel Zwillinger and Stephen Kokoska,
!    CRC Standard Probability and Statistics Tables and Formulae,
!    Chapman and Hall/CRC, 2000, pages 114-115.
!
!  Modified:
!
!    24 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 < X < 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real pdf
  real x
!
  if ( x <= 0.0 .or. x >= 1.0 ) then
    pdf = 0.0
  else
    pdf = 1.0 / ( PI * sqrt ( x * ( 1.0 - x ) ) )
  end if

  return
end
subroutine arcsin_sample ( iseed, x )
!
!*******************************************************************************
!
!! ARCSIN_SAMPLE samples the Arcsin PDF.
!
!
!  Modified:
!
!    17 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call arcsin_cdf_inv ( cdf, x )

  return
end
subroutine arcsin_variance ( variance )
!
!*******************************************************************************
!
!! ARCSIN_VARIANCE returns the variance of the Arcsin PDF.
!
! 
!  Modified:
!
!    17 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real variance
!
  variance = 0.125

  return
end
subroutine benford_pdf ( x, pdf )
!
!*******************************************************************************
!
!! BENFORD_PDF returns the Benford probability of one or more significant digits.
!
!
!  Comments:
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
!    PDF(X) = LOG10 ( ( X + 1 ) / X ).
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
!    Input, integer X, the string of significant digits to be checked.
!    If X is 1, then we are asking for the Benford probability that
!    a value will have first digit 1.  If X is 123, we are asking for
!    the probability that the first three digits will be 123, and so on.
!
!    Output, real PDF, the Benford probability that an item taken
!    from a real world distribution will have the initial digits X.
!
  real pdf
  integer x
!
  if ( x <= 0 ) then
    pdf = 0.0
  else
    pdf = log10 ( real ( x + 1 ) / real ( x ) )
  end if

  return
end
subroutine bernoulli_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! BERNOULLI_CDF evaluates the Bernoulli CDF.
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
!    Input, integer X, the number of successes on a single trial.
!    X = 0 or 1.
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real cdf
  integer x
!
  if ( x < 0 ) then
    cdf = 0.0
  else if ( x == 0 ) then
    cdf = 1.0 - a
  else
    cdf = 1.0
  end if

  return
end
subroutine bernoulli_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! BERNOULLI_CDF_INV inverts the Bernoulli CDF.
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
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, the parameter of the PDF.
!    0.0 <= A <= 1.0.
!
!    Output, integer X, the corresponding argument.
!
  real a
  real cdf
  integer x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BERNOULLI_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 1.0 - a ) then
    x = 0
  else
    x = 1
  end if

  return
end
subroutine bernoulli_check ( a )
!
!*******************************************************************************
!
!! BERNOULLI_CHECK checks the parameter of the Bernoulli CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 <= A <= 1.0.
!
  real a
!
  if ( a < 0.0 .or. a > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BERNOULLI_CHECK - Fatal error!'
    write ( *, * ) '  A < 0 or 1 < A.'
    stop
  end if

  return
end
subroutine bernoulli_mean ( a, mean )
!
!*******************************************************************************
!
!! BERNOULLI_MEAN returns the mean of the Bernoulli PDF.
!
! 
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the probability of success.
!    0.0 <= A <= 1.0.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real mean
!
  mean = a

  return
end
subroutine bernoulli_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! BERNOULLI_PDF evaluates the Bernoulli PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = A**X * ( 1.0 - A )**( X - 1 )
! 
!    X = 0 or 1.
!
!  Discussion:
!
!    The Bernoulli PDF describes the simple case in which a single trial 
!    is carried out, with two possible outcomes, called "success" and 
!    "failure"; the probability of success is A.
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
!    Input, integer X, the number of successes on a single trial.
!    X = 0 or 1.
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real pdf
  integer x
!
  if ( x < 0 ) then
    pdf = 0.0
  else if ( x == 0 ) then
    pdf = 1.0 - a
  else if ( x == 1 ) then
    pdf = a
  else
    pdf = 0.0
  end if

  return
end
subroutine bernoulli_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! BERNOULLI_SAMPLE samples the Bernoulli PDF.
!
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
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  real a
  real cdf
  integer iseed
  real uniform_01_sample
  integer x
!
  cdf = uniform_01_sample ( iseed )

  call bernoulli_cdf_inv ( cdf, a, x )

  return
end
subroutine bernoulli_variance ( a, variance )
!
!*******************************************************************************
!
!! BERNOULLI_VARIANCE returns the variance of the Bernoulli PDF.
!
! 
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real variance
!
  variance = a * ( 1.0 - a )

  return
end
function beta ( a, b )
!
!*******************************************************************************
!
!! BETA returns the value of the Beta function.
!
!
!  Formula:
!
!    BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
!              = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
!
!  Modified:
!
!    10 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the function.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real BETA, the value of the function.
!
  real a
  real b
  real beta
  real gamma_log
!
  if ( a <= 0.0 .or. b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA - Fatal error!'
    write ( *, * ) '  Both A and B must be greater than 0.'
    stop
  end if

  beta = exp ( gamma_log ( a ) + gamma_log ( b ) - gamma_log ( a + b ) )

  return
end
subroutine beta_binomial_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! BETA_BINOMIAL_CDF evaluates the Beta Binomial CDF.
!
!
!  Discussion:
!
!    A simple-minded summing approach is used.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the CDF.
!
!    Input, real A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real beta
  integer c
  real cdf
  real pdf
  integer x
  integer y
!
  if ( x < 0 ) then

    cdf = 0.0

  else if ( x < c ) then

    cdf = 0.0
    do y = 0, x
      pdf = beta ( a + real ( y ), b + real ( c - y ) ) / ( real ( c + 1 ) &
        * beta ( real ( y + 1), real ( c - y + 1 ) ) * beta ( a, b ) )
      cdf = cdf + pdf
    end do

  else if ( x >= c ) then

    cdf = 1.0

  end if

  return
end
subroutine beta_binomial_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! BETA_BINOMIAL_CDF_INV inverts the Beta Binomial CDF.
!
!
!  Discussion:
!
!    A simple-minded discrete approach is used.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
!    Output, integer X, the smallest X whose cumulative density function
!    is greater than or equal to CDF.
!
  real a
  real b
  real beta
  integer c
  real cdf
  real cum
  real pdf
  integer x
  integer y
!
  if ( cdf <= 0.0 ) then

    x = 0

  else

    cum = 0.0

    do y = 0, c

      pdf = beta ( a + real ( y ), b + real ( c - y ) ) / ( real ( c + 1 ) &
        * beta ( real ( y + 1), real ( c - y + 1 ) ) * beta ( a, b ) )

      cum = cum + pdf

      if ( cdf <= cum ) then
        x = y
        return
      end if

    end do

    x = c

  end if

  return
end
subroutine beta_binomial_check ( a, b, c )
!
!*******************************************************************************
!
!! BETA_BINOMIAL_CHECK checks the parameters of the Beta Binomial PDF.
!
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
  real a
  real b
  integer c
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( c < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, * ) '  C < 0.'
    stop
  end if

  return
end
subroutine beta_binomial_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! BETA_BINOMIAL_MEAN returns the mean of the Beta Binomial PDF.
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
!    Input, real A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= N.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  integer c
  real mean
!
  mean = real ( c ) * a / ( a + b )

  return
end
subroutine beta_binomial_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! BETA_BINOMIAL_PDF evaluates the Beta Binomial PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = Beta(A+X,B+C-X) 
!      / ( (C+1) * Beta(X+1,C-X+1) * Beta(A,B) )  for 0 <= X <= C.
!
!    This PDF can be reformulated as:
!
!      The beta binomial probability density function for X successes
!      out of N trials is
!
!      PDF2(X)( N, MU, THETA ) =
!        C(N,X) * Product ( 0 <= R <= X - 1 ) ( MU + R * THETA )
!               * Product ( 0 <= R <= N - X - 1 ) ( 1 - MU + R * THETA )
!               / Product ( 0 <= R <= N - 1 )  ( 1 + R * THETA )
!
!      where 
!
!        C(N,X) is the combinatorial coefficient;
!        MU is the expectation of the underlying Beta distribution;
!        THETA is a shape parameter.  
!
!      A THETA value of 0 ( or A+B --> Infinity ) results in the binomial 
!      distribution:
!
!        PDF2(X) ( N, MU, 0 ) = C(N,X) * MU**X * ( 1 - MU )**(N-X)
!
!    Given A, B, C for PDF, then the equivalent PDF2 has:
!
!      N     = C
!      MU    = A / ( A + B )
!      THETA = 1 / ( A + B )
!
!    Given N, MU, THETA for PDF2, the equivalent PDF has:
!
!      A = MU / THETA
!      B = ( 1 - MU ) / THETA
!      C = N
!
!  Discussion:
!
!    BETA_BINOMIAL_PDF(X)(1,1,C) = UNIFORM_DISCRETE_PDF(X)(0,C-1)
!
!  Modified:
!
!     18 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!
!    Input, real A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real beta
  integer c
  real pdf
  integer x
!
  if ( x < 0 ) then

    pdf = 0.0

  else if ( x <= c ) then

    pdf = beta ( a + real ( x ), b + real ( c - x ) ) / ( real ( c + 1 ) &
      * beta ( real ( x + 1 ), real ( c - x + 1 ) ) * beta ( a, b ) )

  else if ( x > c ) then

    pdf = 0.0

  end if

  return
end
subroutine beta_binomial_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! BETA_BINOMIAL_SAMPLE samples the Beta Binomial CDF.
!
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer X, a sample of the PDF.
!
  real a
  real b
  integer c
  real cdf
  integer iseed
  integer x
!
  call r_random ( 0.0, 1.0, iseed, cdf )

  call beta_binomial_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine beta_binomial_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! BETA_BINOMIAL_VARIANCE returns the variance of the Beta Binomial PDF.
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
!    Input, real A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  integer c
  real variance
!
  variance = ( real ( c ) * a * b ) * ( a + b + real ( c ) ) &
    / ( ( a + b )**2 * ( a + b + 1.0 ) )

  return
end
subroutine beta_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! BETA_CDF evaluates the Beta CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real beta_inc
  real cdf
  real x
!
  if ( x <= 0.0 ) then
    cdf = 0.0
  else if ( x <= 1.0 ) then
    cdf = beta_inc ( x, a, b )
  else
    cdf = 1.0
  end if

  return
end
subroutine beta_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! BETA_CDF_INV inverts the Beta CDF.
!
!
!  Reference:
!
!    Algorithm 724,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 4, December 1993, pages 481-483.
!
!  Modified:
!
!    11 April 1999
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real X, the argument of the CDF.
!
  integer, parameter :: MAXK = 20
  real, parameter :: error = 0.0001
  real, parameter :: errapp = 0.01
!
  real a
  real b
  real bcoeff
  real cdf
  real cdf_x
  real d(2:MAXK,0:MAXK-2)
  integer i
  integer j
  integer k
  integer loopct
  real pdf_x
  real q
  real s1
  real s2
  real sum
  real t
  real tail
  real x
  real xold
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if
!
!  Estimate the solution.
!
  x = a / ( a + b )

  xold = 0.0
  loopct = 2

10    continue

  if ( abs ( ( x - xold ) / x ) >= errapp .and. loopct /= 0 ) then

    xold = x
    loopct = loopct - 1
!
!  CDF_X = PROB { BETA(A,B) <= X }.
!  Q = ( CDF - CDF_X ) / PDF_X.
!
    call beta_cdf ( x, a, b, cdf_x )

    call beta_pdf ( x, a, b, pdf_x )

    q = ( cdf - cdf_x ) / pdf_x
!
!  D(N,K) = C(N,K) * Q**(N+K-1) / (N-1)!
!
    t = 1.0 - x
    s1 = q * ( b - 1.0 ) / t
    s2 = q * ( 1.0 - a ) / x
    d(2,0) = s1 + s2
    tail = d(2,0) * q / 2.0
    x = x + q + tail
    k = 3

20      continue

    if ( abs ( tail / x ) > error .and. k <= MAXK ) then
!
!  Find D(2,K-2).
!
      s1 = q * ( real ( k ) - 2.0 ) * s1 / t
      s2 = q * ( 2.0 - real ( k ) ) * s2 / x
      d(2,k-2) = s1 + s2
!
!  Find D(3,K-3), D(4,K-4), D(5,K-5), ... , D(K-1,1).
!
      do i = 3, k-1
        sum = d(2,0) * d(i-1,k-i)
        bcoeff = 1.0
        do j = 1, k-i
          bcoeff = ( bcoeff * real ( k - i - j + 1 ) ) / real ( j )
          sum = sum + bcoeff * d(2,j) * d(i-1,k-i-j)
        end do
        d(i,k-i) = sum + d(i-1,k-i+1) / real ( i - 1 )
      end do
!
!  Compute D(K,0) and use it to expand the series.
!
      d(k,0) = d(2,0) * d(k-1,0) + d(k-1,1) / real ( k - 1 )
      tail = d(k,0) * q / real ( k )
      x = x + tail
!
!  Check for divergence.
!
      if ( x <= 0.0 .or. x >= 1.0 )  then
        write ( *, * ) ' '
        write ( *, * ) 'BETA_CDF_INV - Fatal error!'
        write ( *, * ) '  The series has diverged.'
        write ( *, * ) '  X = ', x
        x = - 1.0
        return
      end if

      k = k + 1
      go to 20

    end if

    go to 10

  end if

  return
end
subroutine beta_check ( a, b )
!
!*******************************************************************************
!
!! BETA_CHECK checks the parameters of the Beta PDF.
!
!
!  Modified:
!
!    08 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
  real a
  real b
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
function beta_inc ( x, a, b )
!
!*******************************************************************************
!
!! BETA_INC returns the value of the incomplete Beta function.
!
!
!  Formula:
!
!    BETA_INC(X,A,B)
!
!      =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
!        / Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT
!
!      =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
!        / BETA(A,B)
!
!  Reference:
!
!    Algorithm AS63,
!    Applied Statistics,
!    1973, volume 22, number 3.
!
!  Modified:
!
!    03 February 1999
!
!  Parameters:
!
!    Input, real X, the argument of the function.
!    Normally, 0.0 <= X <= 1.0.
!
!    Input, real A, B, the parameters of the function.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real BETA_INC, the value of the function.
!
  real, parameter :: TOL = 1.0E-7
!
  real a
  real b
  real beta
  real beta_inc
  real cx
  integer i
  logical indx
  integer ns
  real pp
  real psq
  real qq
  real rx
  real temp
  real term
  real x
  real xx
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_INC - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_INC - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( x <= 0.0 ) then
    beta_inc = 0.0
    return
  else if ( x >= 1.0 ) then
    beta_inc = 1.0
    return
  end if
!
!  Change tail if necessary and determine S.
!
  psq = a + b
  cx = 1.0 - x

  if ( a < ( a + b ) * x ) then
    xx = cx
    cx = x
    pp = b
    qq = a
    indx = .true.
  else
    xx = x
    pp = a
    qq = b
    indx = .false.
  end if

  term = 1.0
  i = 1
  beta_inc = 1.0
  ns = qq + cx * ( a + b )
!
!  Use Soper's reduction formulas.
!
  rx = xx / cx

10    continue

  temp = qq - real ( i )
  if ( ns == 0 ) then
    rx = xx
  end if

20    continue

  term = term * temp * rx / ( pp + real ( i ) )
  beta_inc = beta_inc + term
  temp = abs ( term )

  if ( temp > TOL .or. temp > TOL * beta_inc ) then

    i = i + 1
    ns = ns - 1

    if ( ns >= 0 ) then
      go to 10
    else
      temp = psq
      psq = psq + 1.0
      go to 20
    end if

  end if
!
!  Finish calculation.
!
  beta_inc = beta_inc * exp ( pp * log ( xx ) + ( qq - 1.0 ) * log ( cx ) ) &
    / ( beta ( a, b ) * pp )

  if ( indx ) then
    beta_inc = 1.0 - beta_inc
  end if

  return
end
subroutine beta_mean ( a, b, mean )
!
!*******************************************************************************
!
!! BETA_MEAN returns the mean of the Beta PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a / ( a + b )

  return
end
subroutine beta_pascal_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! BETA_PASCAL_CDF evaluates the Beta Pascal CDF.
!
!
!  Discussion:
!
!    A simple-minded summing approach is used.
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
!    Input, integer X, the argument of the CDF.
!
!    Input, integer A, real B, C, parameters of the PDF.
!    0 < A,
!    0.0 < B < C.
!
!    Output, real CDF, the value of the CDF.
!
  integer a
  real b
  real c
  real cdf
  real pdf
  integer x
  integer y
!
  cdf = 0.0

  do y = a, x

    call beta_pascal_pdf ( y, a, b, c, pdf )

    cdf = cdf + pdf

  end do

  return
end
subroutine beta_pascal_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! BETA_PASCAL_CDF_INV inverts the Beta Pascal CDF.
!
!
!  Discussion:
!
!    A simple-minded discrete approach is used.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, integer A, real B, C, parameters of the PDF.
!    0 < A,
!    0.0 < B < C.
!
!    Output, integer X, the smallest X whose cumulative density function
!    is greater than or equal to CDF.
!
  integer, parameter :: XMAX = 1000
!
  integer a
  real b
  real c
  real cdf
  real cum
  real pdf
  integer x
!
  if ( cdf <= 0.0 ) then

    x = a

  else if ( cdf < 1.0 ) then

    x = a
    call beta_pascal_pdf ( x, a, b, c, pdf )
    cum = pdf

    do while ( cdf > cum .and. x < XMAX )
      x = x + 1
      call beta_pascal_pdf ( x, a, b, c, pdf )
      cum = cum + pdf
    end do

  else

    write ( *, * ) ' '
    write ( *, * ) 'BETA_PASCAL_CDF_INV - Fatal error!'
    write ( *, * ) '  Input CDF >= 1.'
    stop

  end if

  return
end
subroutine beta_pascal_check ( a, b, c )
!
!*******************************************************************************
!
!! BETA_PASCAL_CHECK checks the parameters of the Beta Pascal PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, real B, C, parameters of the PDF.
!    0 < A,
!    0.0 < B < C.
!
  integer a
  real b
  real c
!
  if ( a <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_PASCAL_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_PASCAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( c <= b ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_PASCAL_CHECK - Fatal error!'
    write ( *, * ) '  C <= B.'
    stop
  end if

  return
end
subroutine beta_pascal_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! BETA_PASCAL_MEAN returns the mean of the Beta Pascal PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, real B, C, parameters of the PDF.
!    0 < A,
!    0.0 < B < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  integer a
  real b
  real c
  real mean
!
  if ( b <= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_PASCAL_MEAN - Fatal error!'
    write ( *, * ) '  The mean is undefined for B <= 1.'
    stop
  end if

  mean = ( real ( a ) * ( c - 1.0 ) ) / ( b - 1.0 )

  return
end
subroutine beta_pascal_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! BETA_PASCAL_PDF evaluates the Beta Pascal PDF.
!
!
!  Discussion:
!
!    Something is wrong.  If A = C, the results seem OK, but as
!    A increases, the PDF gets too small.
!    (5,3,4), the CDF seems to get stuck at 0.1428 = 1/7
!    (6,3,4), the CDF gets stuck at 1/56.
!
!  Formula:
!
!    PDF(X)(A,B,C) = Gamma(X) * Gamma(C) * Gamma(B+C) * Gamma(X-A+C-B)
!      / ( Gamma(A) * Gamma(X-A+1) * Gamma(B) * Gamma(C-B) * Gamma(X+C) )
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!    A <= X.
!
!    Input, integer A, real B, C, parameters of the PDF.
!    0 < A,
!    0.0 < B < C.
!
!    Output, real PDF, the value of the PDF.
!
  integer a
  real b
  real c
  real gamma_log
  real pdf
  real pdf_log
  real ra
  real rx
  integer x
!
  if ( x < a ) then
    pdf = 0.0
  else

    ra = real ( a )
    rx = real ( x )

    pdf_log = &
      + gamma_log ( rx ) &
      + gamma_log ( c ) &
      + gamma_log ( b + c ) &
      + gamma_log ( rx - ra + c - b ) &
      - gamma_log ( ra ) &
      - gamma_log ( rx + 1.0 - ra ) &
      - gamma_log ( b ) &
      - gamma_log ( c - b ) &
      - gamma_log ( rx + c )

    pdf = exp ( pdf_log )

  end if

  return
end
subroutine beta_pascal_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! BETA_PASCAL_SAMPLE samples the Beta Pascal CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, real B, C, parameters of the PDF.
!    0 < A,
!    0.0 < B < C.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer X, a sample of the PDF.
!
  integer a
  real b
  real c
  real cdf
  integer iseed
  integer x
!
  call r_random ( 0.0, 1.0, iseed, cdf )

  call beta_pascal_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine beta_pascal_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! BETA_PASCAL_VARIANCE returns the variance of the Beta Pascal PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, real B, C, parameters of the PDF.
!    0 < A,
!    0.0 < B < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer a
  real b
  real bot
  real c
  real top
  real variance
!
  if ( b <= 2.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BETA_PASCAL_VARIANCE - Fatal error!'
    write ( *, * ) '  The mean is undefined for B <= 2.'
    stop
  end if

  top = real ( a ) * ( real ( a ) + b - 1.0 ) * ( c - 1.0 ) * ( c - b )
  bot = ( b - 1.0 )**2 * ( b - 2.0 )

  variance = top / bot

  return
end
subroutine beta_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! BETA_PDF evaluates the Beta PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = X**(A-1) * (1-X)**(B-1) / BETA(A,B).
!
!  Discussion:
!
!    A = B = 1 yields the Uniform distribution on [0,1].
!    A = B = 1/2 yields the Arcsin distribution.
!        B = 1 yields the power function distribution.
!    A = B -> Infinity tends to the Normal distribution.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 <= X <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real beta
  real pdf
  real x
!
  if ( x < 0.0 .or. x > 1.0 ) then
    pdf = 0.0
  else
    pdf = x**( a - 1.0 ) * ( 1.0 - x )**( b - 1.0 ) / beta ( a, b )
  end if

  return
end
subroutine beta_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! BETA_SAMPLE samples the Beta PDF.
!
!
!  Reference:
!
!    Algorithm BN,
!    William Kennedy and James Gentle,
!    Statistical Computing,
!    Dekker, 1980.
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  integer iseed
  real mu
  real stdev
  real test
  real u
  real uniform_01_sample
  real x
  real y
!
  mu = ( a - 1.0 ) / ( a + b - 2.0 )
  stdev = 0.5 / sqrt ( a + b - 2.0 )

  do

    call normal_01_sample ( iseed, y )

    x = mu + stdev * y

    if ( x < 0.0 .or. x > 1.0 ) then
      cycle
    end if

    u = uniform_01_sample ( iseed )

    test =     ( a - 1.0 )     * log (         x   / ( a - 1.0 ) ) &
             + ( b - 1.0 )     * log ( ( 1.0 - x ) / ( b - 1.0 ) ) &
             + ( a + b - 2.0 ) * log ( a + b - 2.0 ) + 0.5 * y**2

    if ( log ( u ) <= test ) then
      exit
    end if

  end do

  return
end
subroutine beta_variance ( a, b, variance )
!
!*******************************************************************************
!
!! BETA_VARIANCE returns the variance of the Beta PDF.
!
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = ( a * b ) / ( ( a + b )**2 * ( 1.0 + a + b ) )

  return
end
subroutine binomial_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! BINOMIAL_CDF evaluates the Binomial CDF.
!
!
!  Definition:
!
!    CDF(X)(A,B) is the probability of at most X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!  Discussion:
!
!    A sequence of trials with fixed probability of success on
!    any trial is known as a sequence of Bernoulli trials.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the desired number of successes.
!    0 <= X <= A.
!
!    Input, integer A, the number of trials.
!    1 <= A.
!
!    Input, real B, the probability of success on one trial.
!    0.0 <= B <= 1.0.
!
!    Output, real CDF, the value of the CDF.
!
  integer a
  real b
  integer cnk
  real cdf
  integer j
  real pr
  integer x
!
  if ( x < 0 ) then

    cdf = 0.0

  else if ( x >= a ) then

    cdf = 1.0      

  else if ( b == 0.0 ) then

    cdf = 1.0

  else if ( b == 1.0 ) then

    cdf = 0.0

  else

    cdf = 0.0

    do j = 0, x

      call binomial_coef ( a, j, cnk )

      pr = real ( cnk ) * b**j * ( 1.0 - b )**( a - j )

      cdf = cdf + pr

    end do

  end if

  return
end
subroutine binomial_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! BINOMIAL_CDF_INV inverts the Binomial CDF.
!
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, integer A, the number of trials.
!    1 <= A.
!
!    Input, real B, the probability of success on one trial.
!    0.0 <= B <= 1.0.
!
!    Output, integer X, the corresponding argument.
!
  integer a
  real b
  real cdf
  real cdf2
  real pdf
  integer x
  integer x2
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BINOMIAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.0

  do x2 = 0, a

    call binomial_pdf ( x2, a, b, pdf )

    cdf2 = cdf2 + pdf

    if ( cdf <= cdf2 ) then
      x = x2
      return
    end if

  end do

  return
end
subroutine binomial_check ( a, b )
!
!*******************************************************************************
!
!! BINOMIAL_CHECK checks the parameter of the Binomial PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of trials.
!    1 <= A.
!
!    Input, real B, the probability of success on one trial.
!    0.0 <= B <= 1.0.
!
  integer a
  real b
!
  if ( a < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BINOMIAL_CHECK - Fatal error!'
    write ( *, * ) '  A < 1.'
    stop
  end if

  if ( b < 0.0 .or. b > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BINOMIAL_CHECK - Fatal error!'
    write ( *, * ) '  B < 0 or 1 < B.'
    stop
  end if

  return
end
subroutine binomial_coef ( n, k, cnk )
!
!*******************************************************************************
!
!! BINOMIAL_COEF computes the Binomial coefficient C(N,K).
!
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!  Formula:
!
!    CNK = C(N,K) = N! / ( K! * (N-K)! )
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
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, integer CNK, the number of combinations of N
!    things taken K at a time.
!
  integer cnk
  integer i
  integer k
  integer mn
  integer mx
  integer n
!
  mn = min ( k, n-k )

  if ( mn < 0 ) then

    cnk = 0

  else if ( mn == 0 ) then

    cnk = 1

  else

    mx = max ( k, n-k )
    cnk = mx + 1

    do i = 2, mn
      cnk = ( cnk * ( mx + i ) ) / i
    end do

  end if

  return
end
subroutine binomial_coef_log ( n, k, cnk_log )
!
!*******************************************************************************
!
!! BINOMIAL_COEF_LOG computes the logarithm of the Binomial coefficient.
!
!
!  Formula:
!
!    CNK_LOG = LOG ( C(N,K) ) = LOG ( N! / ( K! * (N-K)! ) ).
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, real CNK_LOG, the logarithm of C(N,K).
!
  real cnk_log
  real factorial_log
  integer k
  integer n
!
  cnk_log = factorial_log ( n ) - factorial_log ( k ) - factorial_log ( n - k )

  return
end
subroutine binomial_mean ( a, b, mean )
!
!*******************************************************************************
!
!! BINOMIAL_MEAN returns the mean of the Binomial PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of trials.
!    1 <= A.
!
!    Input, real B, the probability of success on one trial.
!    0.0 <= B <= 1.0.
!
!    Output, real MEAN, the expected value of the number of
!    successes in A trials.
!
  integer a
  real b
  real mean
!
  mean = real ( a ) * b

  return
end
subroutine binomial_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! BINOMIAL_PDF evaluates the Binomial PDF.
!
!
!  Definition:
!
!    PDF(X)(A,B) is the probability of exactly X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!  Formula:
!
!    PDF(X)(A,B) = C(N,X) * B**X * ( 1.0 - B )**( A - X )
!
!  Discussion:
!
!    Binomial_PDF(X)(1,B) = Bernoulli_PDF(X)(B).
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the desired number of successes.
!    0 <= X <= A.
!
!    Input, integer A, the number of trials.
!    1 <= A.
!
!    Input, real B, the probability of success on one trial.
!    0.0 <= B <= 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  integer a
  real b
  integer cnk
  real pdf
  integer x
!
  if ( a < 1 ) then

    pdf = 0.0

  else if ( x < 0 .or. x > a ) then

    pdf = 0.0

  else if ( b == 0.0 ) then

    if ( x == 0 ) then
      pdf = 1.0
    else
      pdf = 0.0
    end if

  else if ( b == 1.0 ) then

    if ( x == a ) then
      pdf = 1.0
    else
      pdf = 0.0
    end if
    
  else

    call binomial_coef ( a, x, cnk )

    pdf = real ( cnk ) * b**x * ( 1.0 - b )**( a - x )

  end if

  return
end
subroutine binomial_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! BINOMIAL_SAMPLE samples the Binomial PDF.
!
!
!  Reference:
!
!    Algorithm BU,
!    William Kennedy and James Gentle,
!    Statistical Computing,
!    Dekker, 1980.
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
!    Input, integer A, the number of trials.
!    1 <= A.
!
!    Input, real B, the probability of success on one trial.
!    0.0 <= B <= 1.0.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  integer a
  real b
  integer iseed
  integer i
  real u
  real uniform_01_sample
  integer x
!
  x = 0

  do i = 1, a

    u = uniform_01_sample ( iseed )

    if ( u <= b ) then
      x = x + 1
    end if

  end do

  return
end
subroutine binomial_variance ( a, b, variance )
!
!*******************************************************************************
!
!! BINOMIAL_VARIANCE returns the variance of the Binomial PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of trials.
!    1 <= A.
!
!    Input, real B, the probability of success on one trial.
!    0.0 <= B <= 1.0.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer a
  real b
  real variance
!
  variance = real ( a ) * b * ( 1.0 - b )

  return
end
subroutine bradford_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! BRADFORD_CDF evaluates the Bradford CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    A < B,
!    0.0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real x
!
  if ( x <= a ) then
    cdf = 0.0
  else if ( x <= b ) then
    cdf = log ( 1.0 + c * ( x - a ) / ( b - a ) ) / log ( c + 1.0 )
  else if ( x > b ) then
    cdf = 1.0
  end if

  return
end
subroutine bradford_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! BRADFORD_CDF_INV inverts the Bradford CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, C, the parameters of the PDF.
!    A < B,
!    0.0 < C.
!
!    Output, real X, the corresponding argument of the CDF.
!
  real a
  real b
  real c
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BRADFORD_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.0 ) then
    x = a
  else if ( cdf < 1.0 ) then
    x = a + ( b - a ) * ( ( c + 1.0 )**cdf - 1.0 ) / c
  else if ( cdf >= 1.0 ) then
    x = b
  end if

  return
end
subroutine bradford_check ( a, b, c )
!
!*******************************************************************************
!
!! BRADFORD_CHECK checks the parameters of the Bradford PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    A < B,
!    0.0 < C.
!
  real a
  real b
  real c
!
  if ( a >= b ) then
    write ( *, * ) ' '
    write ( *, * ) 'BRADFORD_CHECK - Fatal error!'
    write ( *, * ) '  A >= B.'
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BRADFORD_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    stop
  end if

  return
end
subroutine bradford_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! BRADFORD_MEAN returns the mean of the Bradford PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    A < B,
!    0.0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real mean
!
  mean = ( c * ( b - a ) + log ( c + 1.0 ) * ( a * ( c + 1.0 ) - b ) ) &
    / ( c * log ( c + 1.0 ) )

  return
end
subroutine bradford_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! BRADFORD_PDF evaluates the Bradford PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = 
!      C / ( ( C * ( X - A ) + B - A ) * log ( C + 1 ) )
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, C, the parameters of the PDF.
!    A < B,
!    0.0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real pdf
  real x
!
  if ( x <= a ) then
    pdf = 0.0
  else if ( x <= b ) then
    pdf = c / ( ( c * ( x - a ) + b - a ) * log ( c + 1.0 ) )
  else if ( x > b ) then
    pdf = 0.0
  end if

  return
end
subroutine bradford_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! BRADFORD_SAMPLE samples the Bradford PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    A < B,
!    0.0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real c
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  x = a + ( b - a ) * ( ( c + 1.0 )**cdf - 1.0 ) / c

  return
end
subroutine bradford_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! BRADFORD_VARIANCE returns the variance of the Bradford PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    A < B,
!    0.0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real c
  real variance
!
  variance = ( b - a )**2 * &
    ( c * ( log ( c + 1.0 ) - 2.0 ) + 2.0 * log ( c + 1.0 ) ) &
    / ( 2.0 * c * ( log ( c + 1.0 ) )**2 )

  return
end
subroutine burr_cdf ( x, a, b, c, d, cdf )
!
!*******************************************************************************
!
!! BURR_CDF evaluates the Burr CDF.
!
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real d
  real x
!
  if ( x <= a ) then

    cdf = 0.0

  else

    cdf = 1.0 / ( 1.0 + ( b / ( x - a ) )**c  )**d

  end if

  return
end
subroutine burr_cdf_inv ( cdf, a, b, c, d, x )
!
!*******************************************************************************
!
!! BURR_CDF_INV inverts the Burr CDF.
!
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real c
  real cdf
  real d
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BURR_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a + b / ( ( 1.0 / cdf )**(1.0 / d ) - 1.0 )**( 1.0 / c )

  return
end
subroutine burr_check ( a, b, c, d )
!
!*******************************************************************************
!
!! BURR_CHECK checks the parameters of the Burr CDF.
!
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
  real a
  real b
  real c
  real d
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BURR_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( c <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'BURR_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    stop
  end if

  return
end
subroutine burr_mean ( a, b, c, d, mean )
!
!*******************************************************************************
!
!! BURR_MEAN returns the mean of the Burr PDF.
!
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real d
  real gamma
  real mean
!
  mean = a + b * gamma ( 1.0 - 1.0 / c ) * gamma ( d + 1.0 / c ) / gamma ( d )

  return
end
subroutine burr_pdf ( x, a, b, c, d, pdf )
!
!*******************************************************************************
!
!! BURR_PDF evaluates the Burr PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C,D) = ( C * D / B ) * ( ( X - A ) / B )**( - C - 1 )
!      * ( 1 + ( ( X - A ) / B )**( - C ) )**( - D - 1 )
!
!  Reference:
!
!    M E Johnson,
!    Multivariate Statistical Simulation,
!    Wiley, New York, 1987.
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real d
  real pdf
  real x
  real y
!
  if ( x <= a ) then
    pdf = 0.0
  else

    y = ( x - a ) / b

    pdf = ( c * d / b ) * y**( - c - 1.0 ) * ( 1.0 + y**( - c ) )**( - d - 1.0 )

  end if

  return
end
subroutine burr_sample ( a, b, c, d, iseed, x )
!
!*******************************************************************************
!
!! BURR_SAMPLE samples the Burr PDF.
!
!
!  Modified:
!
!    01 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real c
  real cdf
  real d
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call burr_cdf_inv ( cdf, a, b, c, d, x )

  return
end
subroutine burr_variance ( a, b, c, d, variance )
!
!*******************************************************************************
!
!! BURR_VARIANCE returns the variance of the Burr PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real c
  real d
  real gamma
  real k
  real variance
!
  if ( c <= 2.0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'BURR_VARIANCE - Warning!'
    write ( *, * ) '  Variance undefined for C <= 2.'
    variance = huge ( variance )

  else

    k = gamma ( d ) * gamma ( 1.0 - 2.0 / c ) * gamma ( d + 2.0 / c ) &
      - ( gamma ( 1.0 - 1.0 / c ) * gamma ( d + 1.0 / c ) )**2
  
    variance = k * b**2 / ( gamma ( d ) )**2

  end if

  return
end
subroutine cauchy_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! CAUCHY_CDF evaluates the Cauchy CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  real x
  real y
!
  y = ( x - a ) / b

  cdf = 0.5 + atan ( y ) / PI

  return
end
subroutine cauchy_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! CAUCHY_CDF_INV inverts the Cauchy CDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  real x
!
  x = a + b * tan ( PI * ( cdf - 0.5 ) )

  return
end
subroutine cauchy_check ( a, b )
!
!*******************************************************************************
!
!! CAUCHY_CHECK checks the parameters of the Cauchy CDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CAUCHY_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine cauchy_mean ( a, b, mean )
!
!*******************************************************************************
!
!! CAUCHY_MEAN returns the mean of the Cauchy PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine cauchy_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! CAUCHY_PDF evaluates the Cauchy PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = 1 / ( PI * B * ( 1 + ( ( X - A ) / B )**2 ) )
!
!  Discussion:
!
!    The Cauchy PDF is also known as the Breit-Wigner PDF.  It
!    has some unusual properties.  In particular, the integrals for the
!    expected value and higher order moments are "singular", in the
!    sense that the limiting values do not exist.  A result can be
!    obtained if the upper and lower limits of integration are set
!    equal to +T and -T, and the limit as T=>INFINITY is taken, but
!    this is a very weak and unreliable sort of limit.
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
  real y
!
  y = ( x - a ) / b

  pdf = 1.0 / ( PI * b * ( 1.0 + y**2 ) )

  return
end
subroutine cauchy_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! CAUCHY_SAMPLE samples the Cauchy PDF.
!
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call cauchy_cdf_inv ( cdf, a, b, x )

  return
end
subroutine cauchy_variance ( a, b, variance )
!
!*******************************************************************************
!
!! CAUCHY_VARIANCE returns the variance of the Cauchy PDF.
!
!
!  Discussion:
!
!    The variance of the Cauchy PDF is not well defined.  This routine
!    is made available for completeness only, and simply returns
!    a "very large" number.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the mean of the PDF.
!
  real a
  real b
  real variance
!
  variance = huge ( variance )

  return
end
subroutine chi_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! CHI_CDF evaluates the Chi CDF.
!
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real gamma_inc
  real p2
  real x
  real x2
  real y
!
  if ( x <= a ) then

    cdf = 0.0

  else

    y = ( x - a ) / b
    x2 = 0.5 * y**2
    p2 = 0.5 * c

    cdf = gamma_inc ( x2, p2 )

  end if

  return
end
subroutine chi_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! CHI_CDF_INV inverts the Chi CDF.
!
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Modified:
!
!    05 January 2000
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real X, the corresponding argument of the CDF.
!
  integer, parameter :: IT_MAX = 10
  real, parameter :: TOL = 0.0001
!
  real a
  real b
  real c
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
!
  if ( cdf <= 0.0 ) then
    x = a
    return
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
    return
  end if

  x1 = a
  cdf1 = 0.0

  x2 = a + 1.0

  do

    call chi_cdf ( x2, a, b, c, cdf2 )

    if ( cdf2 > cdf ) then
      exit
    end if

    x2 = a + 2.0 * ( x2 - a )

  end do
!
!  Now use bisection.
!
  it = 0

  do

    it = it + 1

    x3 = 0.5 * ( x1 + x2 )
    call chi_cdf ( x3, a, b, c, cdf3 )

    if ( abs ( cdf3 - cdf ) < TOL ) then
      x = x3
      return
    end if

    if ( it > IT_MAX ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHI_CDF_INV - Fatal error!'
      write ( *, * ) '  Iteration limit exceeded.'
      stop
    end if

    if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  return
end
subroutine chi_check ( a, b, c )
!
!*******************************************************************************
!
!! CHI_CHECK checks the parameters of the Chi CDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
  real a
  real b
  real c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHI_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.0.'
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHI_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.0.'
    stop
  end if

  return
end
subroutine chi_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! CHI_MEAN returns the mean of the Chi PDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real MEAN, the mean value.
!
  real a
  real b
  real c
  real gamma
  real mean
!
  mean = a + sqrt ( 2.0 ) * b * gamma ( 0.5 * ( c + 1.0 ) ) / gamma ( 0.5 * c ) 

  return
end
subroutine chi_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! CHI_PDF evaluates the Chi PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = EXP ( - 0.5 * ( ( X - A ) / B )**2 ) 
!      * ( ( X - A ) / B )**( C - 1 ) /
!      ( 2**( 0.5 * C - 1 ) * B * GAMMA ( 0.5 * C ) )
!      
!  Discussion:
!
!    CHI(A,B,1) is the Half Normal PDF;
!    CHI(0,B,2) is the Rayleigh PDF;
!    CHI(0,B,3) is the Maxwell PDF.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real gamma
  real pdf
  real x
  real y
!
  if ( x <= a ) then

    pdf = 0.0

  else

    y = ( x - a ) / b

    pdf = exp ( - 0.5 * y**2 ) * y**( c - 1.0 ) / &
      ( 2.0**( 0.5 * c - 1.0 ) * b * gamma ( 0.5 * c ) )

  end if

  return
end
subroutine chi_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! CHI_SAMPLE samples the Chi PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real c
  integer iseed
  real x
!
  call chisquare_central_sample ( c, iseed, x )

  x = a + b * sqrt ( x )

  return
end
subroutine chi_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! CHI_VARIANCE returns the variance of the Chi PDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real c
  real gamma
  real variance
!
  variance = b**2 * ( c - 2.0 * &
    ( gamma ( 0.5 * ( c + 1.0 ) ) / gamma ( 0.5 * c ) )**2 )

  return
end
subroutine chisquare_central_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! CHISQUARE_CENTRAL_CDF evaluates the chi-squared CDF.
!
!
!  Modified:
!
!    02 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the value of the random deviate.
!
!    Input, real A, the parameter of the distribution, usually
!    the number of degrees of freedom.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real a2
  real b2
  real c2
  real cdf
  real x
  real x2
!
  x2 = 0.5 * x

  a2 = 0.0
  b2 = 1.0
  c2 = 0.5 * a

  call gamma_cdf ( x2, a2, b2, c2, cdf )

  return
end
subroutine chisquare_central_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! CHISQUARE_CENTRAL_CDF_INV inverts the chi-squared PDF.
!
!
!  Reference:
!
!    Best and Roberts,
!    The Percentage Points of the Chi-Squared Distribution,
!    Algorithm AS 91,
!    Applied Statistics,
!    Volume 24, Number ?, pages 385-390, 1975.
!
!  Modified:
!
!    30 November 1999
!
!  Parameters:
!
!    Input, real CDF, a value of the chi-squared cumulative probability
!    density function.
!    0.000002 <= CDF <= 0.999998.
!
!    Input, real A, the parameter of the chi-squared probability density
!    function.  A > 0.
!
!    Output, real X, the value of the chi-squared random deviate
!    with the property that the probability that a chi-squared random
!    deviate with parameter A is less than or equal to PPCHI2 is P.
!
  real, parameter :: aa = 0.6931471806
  real, parameter :: c1 = 0.01
  real, parameter :: c2 = 0.222222
  real, parameter :: c3 = 0.32
  real, parameter :: c4 = 0.4
  real, parameter :: c5 = 1.24
  real, parameter :: c6 = 2.2
  real, parameter :: c7 = 4.67
  real, parameter :: c8 = 6.66
  real, parameter :: c9 = 6.73
  real, parameter :: c10 = 13.32
  real, parameter :: c11 = 60.0
  real, parameter :: c12 = 70.0
  real, parameter :: c13 = 84.0
  real, parameter :: c14 = 105.0
  real, parameter :: c15 = 120.0
  real, parameter :: c16 = 127.0
  real, parameter :: c17 = 140.0
  real, parameter :: c18 = 175.0
  real, parameter :: c19 = 210.0
  real, parameter :: c20 = 252.0
  real, parameter :: c21 = 264.0
  real, parameter :: c22 = 294.0
  real, parameter :: c23 = 346.0
  real, parameter :: c24 = 420.0
  real, parameter :: c25 = 462.0
  real, parameter :: c26 = 606.0
  real, parameter :: c27 = 672.0
  real, parameter :: c28 = 707.0
  real, parameter :: c29 = 735.0
  real, parameter :: c30 = 889.0
  real, parameter :: c31 = 932.0
  real, parameter :: c32 = 966.0
  real, parameter :: c33 = 1141.0
  real, parameter :: c34 = 1182.0
  real, parameter :: c35 = 1278.0
  real, parameter :: c36 = 1740.0
  real, parameter :: c37 = 2520.0
  real, parameter :: c38 = 5040.0
  real, parameter :: cdf_max = 0.999998
  real, parameter :: cdf_min = 0.000002
  real, parameter :: e = 0.0000005
  integer, parameter :: maxit = 20
!
  real a
  real a2
  real b
  real c
  real cdf
  real ch
  real g
  real gamma_inc
  integer i
  real gamma_log
  real p1
  real p2
  real q
  real s1
  real s2
  real s3
  real s4
  real s5
  real s6
  real t
  real x
  real x2
  real xx
!
  if ( cdf < cdf_min ) then
    x = - 1.0
    write ( *, * ) ' '
    write ( *, * ) 'CHISQUARE_CENTRAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < CDF_MIN.'
    stop
  end if

  if ( cdf > cdf_max ) then
    x = - 1.0
    write ( *, * ) ' '
    write ( *, * ) 'CHISQUARE_CENTRAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF > CDF_MAX.'
    stop
  end if

  xx = 0.5 * a
  c = xx - 1.0
!
!  Compute Log ( Gamma ( A/2 ) ).
!
  g = gamma_log ( a / 2.0 )
!
!  Starting approximation for small chi-squared.
!
  if ( a < - c5 * log ( cdf ) ) then

    ch = ( cdf * xx * exp ( g + xx * aa ) )**( 1.0 / xx )

    if ( ch < e ) then
      x = ch
      return
    end if
!
!  Starting approximation for A less than or equal to 0.32.
!
  else if ( a <= c3 ) then

    ch = c4
    a2 = log ( 1.0 - cdf )

    do

      q = ch
      p1 = 1.0 + ch * ( c7 + ch )
      p2 = ch * ( c9 + ch * ( c8 + ch ) )

      t = - 0.5 + ( c7 + 2.0 * ch ) / p1 - ( c9 + ch * ( c10 + 3.0 * ch ) ) / p2

      ch = ch - ( 1.0 - exp ( a2 + g + 0.5 * ch + c * aa ) * p2 / p1 ) / t

      if ( abs ( q / ch - 1.0 ) <= c1 ) then
        exit
      end if

    end do
!
!  Call to algorithm AS 111.
!  Note that P has been tested above.
!  AS 241 could be used as an alternative.
!
  else

    call normal_01_cdf_inv ( cdf, x2 )
!
!  Starting approximation using Wilson and Hilferty estimate.
!
    p1 = c2 / a
    ch = a * ( x2 * sqrt ( p1 ) + 1.0 - p1 )**3
!
!  Starting approximation for P tending to 1.
!
    if ( ch > c6 * a + 6.0 ) then
      ch = - 2.0 * ( log ( 1.0 - cdf ) - c * log ( 0.5 * ch ) + g )
    end if

  end if
!
!  Call to algorithm AS 239 and calculation of seven term Taylor series.
!
  do i = 1, maxit

    q = ch
    p1 = 0.5 * ch
    p2 = cdf - gamma_inc ( p1, xx )
    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) )
    b = t / ch
    a2 = 0.5 * t - b * c

    s1 = ( c19 + a2 &
       * ( c17 + a2 &
       * ( c14 + a2 &
       * ( c13 + a2 &
       * ( c12 + a2 &
       *   c11 ) ) ) ) ) / c24

    s2 = ( c24 + a2 &
       * ( c29 + a2 &
       * ( c32 + a2 &
       * ( c33 + a2 &
       *   c35 ) ) ) ) / c37

    s3 = ( c19 + a2 &
       * ( c25 + a2 &
       * ( c28 + a2 &
       *   c31 ) ) ) / c37

    s4 = ( c20 + a2 &
       * ( c27 + a2 &
       *   c34 ) + c &
       * ( c22 + a2 &
       * ( c30 + a2 &
       *   c36 ) ) ) / c38

    s5 = ( c13 + c21 * a2 + c * ( c18 + c26 * a2 ) ) / c37

    s6 = ( c15 + c * ( c23 + c16 * c ) ) / c38

    ch = ch + t * ( 1.0 + 0.5 * t * s1 - b * c &
      * ( s1 - b &
      * ( s2 - b &
      * ( s3 - b &
      * ( s4 - b &
      * ( s5 - b &
      *   s6 ) ) ) ) ) )

    if ( abs ( q / ch - 1.0 ) > e ) then
      x = ch
      return
    end if

  end do

  x = ch
  write ( *, * ) ' '
  write ( *, * ) 'CHISQUARE_CENTRAL_CDF_INV - Warning!'
  write ( *, * ) '  Convergence not reached.'

  return
end
subroutine chisquare_central_check ( a )
!
!*******************************************************************************
!
!! CHISQUARE_CENTRAL_CHECK checks the parameter of the central Chi-Squared PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the distribution.
!    1 <= A.
!
  real a
!
  if ( a < 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHISQUARE_CENTRAL_CHECK - Fatal error!'
    write ( *, * ) '  A < 1.0.'
    stop
  end if

  return
end
subroutine chisquare_central_mean ( a, mean )
!
!*******************************************************************************
!
!! CHISQUARE_CENTRAL_MEAN returns the mean of the central Chi-Squared PDF.
!
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the distribution.
!    1 <= A.
!
!    Output, real MEAN, the mean value.
!
  real a
  real mean
!
  mean = a

  return
end
subroutine chisquare_central_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! CHISQUARE_CENTRAL_PDF evaluates the central Chi-Squared PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = 
!      EXP ( - X / 2 ) * X**((A-2)/2) / ( 2**(A/2) * GAMMA ( A/2 ) )
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 <= X 
!
!    Input, real A, the parameter of the PDF.
!    1 <= A.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real gamma
  real pdf
  real x
!
  if ( x < 0.0 ) then
    pdf = 0.0
  else
    b = a / 2.0
    pdf = exp ( - 0.5 * x ) * x**( b - 1.0 ) / ( 2.0**b * gamma ( b ) )
  end if

  return
end
subroutine chisquare_central_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! CHISQUARE_CENTRAL_SAMPLE samples the central Chisquare PDF.
!
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    1 <= A.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  integer, parameter :: IT_MAX = 100
!
  real a
  real a2
  real b2
  real c2
  integer i
  integer iseed
  integer n
  real x
  real x2
!
  n = int ( a )

  if ( real ( n ) == a .and. n <= IT_MAX ) then

    x = 0.0
    do i = 1, n
      call normal_01_sample ( iseed, x2 )
      x = x + x2**2
    end do

  else

    a2 = 0.0
    b2 = 1.0
    c2 = a / 2.0

    call gamma_sample ( a2, b2, c2, iseed, x )

    x = 2.0 * x

  end if

  return
end
subroutine chisquare_central_variance ( a, variance )
!
!*******************************************************************************
!
!! CHISQUARE_CENTRAL_VARIANCE returns the variance of the central Chisquare PDF.
!
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the distribution.
!    1 <= A.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real variance
!
  variance = 2.0 * a

  return
end
subroutine chisquare_noncentral_check ( a, b )
!
!*******************************************************************************
!
!! CHISQUARE_NONCENTRAL_CHECK checks the parameters of the noncentral Chi-Squared PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the parameter of the PDF.
!    1.0 <= A.
!
!    Input, real B, the noncentrality parameter of the PDF.
!    0.0 <= B.
!
  real a
  real b
!
  if ( a < 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHISQUARE_NONCENTRAL_CHECK - Fatal error!'
    write ( *, * ) '  A < 1.'
    stop
  end if

  if ( b < 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHISQUARE_NONCENTRAL_CHECK - Fatal error!'
    write ( *, * ) '  B < 0.'
    stop
  end if

  return
end
subroutine chisquare_noncentral_mean ( a, b, mean )
!
!*******************************************************************************
!
!! CHISQUARE_NONCENTRAL_MEAN returns the mean of the noncentral Chi-Squared PDF.
!
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the parameter of the PDF.
!    1.0 <= A.
!
!    Input, real B, the noncentrality parameter of the PDF.
!    0.0 <= B.
!
!    Output, real MEAN, the mean value.
!
  real a
  real b
  real mean
!
  mean = a + b

  return
end
subroutine chisquare_noncentral_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! CHISQUARE_NONCENTRAL_SAMPLE samples the noncentral Chisquare PDF.
!
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the parameter of the PDF.
!    1.0 <= A.
!
!    Input, real B, the noncentrality parameter of the PDF.
!    0.0 <= B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real a1
  real a2
  real b
  real b2
  integer iseed
  real x
  real x1
  real x2
!
  a1 = a - 1.0

  call chisquare_central_sample ( a1, iseed, x1 )

  a2 = sqrt ( b )
  b2 = 1.0
  call normal_sample ( a2, b2, iseed, x2 )

  x = x1 + x2**2

  return
end
subroutine chisquare_noncentral_variance ( a, b, variance )
!
!*******************************************************************************
!
!! CHISQUARE_NONCENTRAL_VARIANCE returns the variance of the noncentral Chi-Squared PDF.
!
!
!  Modified:
!
!    30 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    1 <= A.
!
!    Input, real B, the noncentrality parameter of the PDF.
!    0.0 <= B.
!
!    Output, real VARIANCE, the variance value.
!
  real a
  real b
  real variance
!
  variance = 2.0 * ( a + 2.0 * b )

  return
end
subroutine circle_sample ( a, b, c, iseed, x1, x2 )
!
!*******************************************************************************
!
!! CIRCLE_SAMPLE samples points from a circle.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the circle.
!    The circle is centered at (A,B) and has radius C.
!    0.0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X1, X2, a sampled point of the circle.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real angle
  real b
  real c
  integer iseed
  real radius_frac
  real uniform_01_sample
  real x1
  real x2
!
  radius_frac = sqrt ( uniform_01_sample ( iseed ) )

  call r_random ( 0.0, 2.0 * PI, iseed, angle )

  x1 = a + c * radius_frac * cos ( angle )
  x2 = b + c * radius_frac * sin ( angle )

  return
end
subroutine circular_normal_01_mean ( mean )
!
!*******************************************************************************
!
!! CIRCULAR_NORMAL_01_MEAN returns the mean of the Circular Normal 01 PDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MEAN(2), the mean of the PDF.
!
  real mean(2)
!
  mean(1) = 0.0
  mean(2) = 0.0

  return
end
subroutine circular_normal_01_pdf ( x, pdf )
!
!*******************************************************************************
!
!! CIRCULAR_NORMAL_01_PDF evaluates the  Circular Normal 01 PDF.
!
!
!  Formula:
!
!    PDF(X) = EXP ( - 0.5 * ( X(1)**2 + X(2)**2 ) ) / ( 2 * PI )
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X(2), the argument of the PDF.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real pdf
  real x(2)
!
  pdf = exp ( - 0.5 * ( x(1)**2 + x(2)**2 ) ) / ( 2.0 * PI )

  return
end
subroutine circular_normal_01_sample ( iseed, x )
!
!*******************************************************************************
!
!! CIRCULAR_NORMAL_01_SAMPLE samples the Circular Normal 01 PDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X(2), a sample of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  integer iseed
  real uniform_01_sample
  real v1
  real v2
  real x(2)
!
  v1 = uniform_01_sample ( iseed )

  v2 = uniform_01_sample ( iseed )

  x(1) = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * PI * v2 )
  x(2) = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * PI * v2 )

  return
end
subroutine circular_normal_01_variance ( variance )
!
!*******************************************************************************
!
!! CIRCULAR_NORMAL_01_VARIANCE returns the variance of the Circular Normal 01 PDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VARIANCE(2), the variance of the PDF.
!
  real variance(2)
!
  variance(1) = 1.0
  variance(2) = 1.0

  return
end
subroutine cosine_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! COSINE_CDF evaluates the Cosine CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  real x
  real y
!
  if ( x <= a - PI * b ) then

    cdf = 0.0

  else if ( x <= a + PI * b ) then

    y = ( x - a ) / b

    cdf = 0.5 + ( y + sin ( y ) ) / ( 2.0 * PI )

  else if ( x > a + PI * b ) then

    cdf = 1.0

  end if

  return
end
subroutine cosine_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! COSINE_CDF_INV inverts the Cosine CDF.
!
!
!  Discussion:
!
!    A simple bisection method is used on the interval 
!    [ A - PI * B, A + PI * B ].
!
!  Modified:
!
!    30 December 1999
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument of the CDF.
!
  integer, parameter :: IT_MAX = 100
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
  real, parameter :: TOL = 0.0001
!
  real a
  real b
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
!
  if ( cdf <= 0.0 ) then
    x = a - PI * b
    return
  else if ( cdf >= 1.0 ) then
    x = a + PI * b
    return
  end if
!
  x1 = a - PI * b
  cdf1 = 0.0

  x2 = a + PI * b
  cdf2 = 1.0
!
!  Now use bisection.
!
  it = 0

  do it = 1, IT_MAX

    x3 = 0.5 * ( x1 + x2 )
    call cosine_cdf ( x3, a, b, cdf3 )

    if ( abs ( cdf3 - cdf ) < TOL ) then
      x = x3
      return
    end if

    if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
      x1 = x3
      cdf1 = cdf3
    else
      x2 = x3
      cdf2 = cdf3
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'COSINE_CDF_INV - Fatal error!'
  write ( *, * ) '  Iteration limit exceeded.'

  stop
end
subroutine cosine_check ( a, b )
!
!*******************************************************************************
!
!! COSINE_CHECK checks the parameters of the Cosine CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'COSINE_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.0'
    stop
  end if

  return
end
subroutine cosine_mean ( a, b, mean )
!
!*******************************************************************************
!
!! COSINE_MEAN returns the mean of the Cosine PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine cosine_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! COSINE_PDF evaluates the Cosine PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = ( 1 / ( 2 * PI * B ) ) * COS ( ( X - A ) / B )
!    for A - PI * B <= X <= A + PI * B
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
  real y
!
  if ( x < a - PI * b ) then

    pdf = 0.0

  else if ( x <= a + PI * b ) then

    y = ( x - a ) / b

    pdf = 1.0 / ( 2.0 * PI * b ) * cos ( y )

  else if ( x > a + PI * b ) then

    pdf = 0.0

  end if

  return
end
subroutine cosine_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! COSINE_SAMPLE samples the Cosine PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )
  call cosine_cdf_inv ( cdf, a, b, x )

  return
end
subroutine cosine_variance ( a, b, variance )
!
!*******************************************************************************
!
!! COSINE_VARIANCE returns the variance of the Cosine PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real variance
!
  variance = ( PI**2 / 3.0 - 2.0 ) * b**2

  return
end
subroutine coupon_simulate ( n_type, iseed, coupon, n_coupon )
!
!*******************************************************************************
!
!! COUPON_SIMULATE simulates the coupon collector's problem.
!
!
!  Discussion:
!
!    The coupon collector needs to collect one of each of N_TYPE
!    coupons.  The collector may draw one coupon on each trial,
!    and takes as many trials as necessary to complete the task.
!    On each trial, the probability of picking any particular type
!    of coupon is always 1 / N_TYPE.
!
!    The most interesting question is, what is the expected number
!    of drawings necessary to complete the collection?
!    how does this number vary as N_TYPE increases?  What is the
!    distribution of the numbers of each type of coupon in a typical 
!    collection when it is just completed?
!
!    As N increases, the number of coupons necessary to be 
!    collected in order to get a complete set in any simulation 
!    strongly tends to the value N_TYPE * LOG ( N_TYPE ).
!
!    If N_TYPE is 1, the simulation ends with a single drawing.
!
!    If N_TYPE is 2, then we may call the coupon taken on the first drawing 
!    a "Head", say, and the process then is similar to the question of the 
!    length, plus one, of a run of Heads or Tails in coin flipping.
!
!  Modified:
!
!    31 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_TYPE, the number of types of coupons.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer COUPON(N_TYPE), the number of coupons of each type
!    that were collected during the simulation.
!
!    Output, integer N_COUPON, the total number of coupons collected.
!
  integer, parameter :: max_coupon = 2000
!
  integer n_type
!
  integer coupon(n_type)
  integer i
  integer iseed
  integer n_coupon
  integer straight
!
  coupon(1:n_type) = 0

  straight = 0
  n_coupon = 0
!
!  Draw another coupon.
!
  do while ( n_coupon < max_coupon )

    call i_random ( 1, n_type, iseed, i )
!
!  Increment the number of I coupons.
! 
    coupon(i) = coupon(i) + 1
    n_coupon = n_coupon + 1
!
!  If I is the next one we needed, increase STRAIGHT by 1.
!
    if ( i == straight + 1 ) then

      do

        straight = straight + 1
!
!  If STRAIGHT = N_TYPE, we have all of them.
!
        if ( straight >= n_type ) then
          return
        end if
!
!  If the next coupon has not been collected, our straight is over.
!
        if ( coupon(straight+1) <= 0 ) then
          exit
        end if

      end do

    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'COUPON_SIMULATE - Fatal error!'
  write ( *, * ) '  Maximum number of coupons drawn without success.'

  stop
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

  if ( csc == 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CSC - Fatal error!'
    write ( *, * ) '  CSC undefined for THETA = ', theta
    stop
  end if

  csc = 1.0 / csc

  return
end
subroutine deranged_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! DERANGED_CDF evaluates the Deranged CDF.
!
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the maximum number of items in their correct places.
!    0 <= X <= A.
!
!    Input, integer A, the number of items.
!    1 <= A.
!
!    Output, real CDF, the value of the CDF.
!
  integer a
  real cdf
  integer cnk
  integer deranged_enum
  integer dnmk
  integer i_factorial
  integer nfact
  integer sum
  integer x
  integer x2
!
  if ( x < 0 .or. x > a ) then
    cdf = 0.0
  else
    sum = 0
    do x2 = 0, x
      call binomial_coef ( a, x2, cnk )
      dnmk = deranged_enum ( a-x2 )
      sum = sum + cnk * dnmk
    end do
    nfact = i_factorial ( a )
    cdf = real ( sum ) / real ( nfact )
  end if

  return
end
subroutine deranged_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! DERANGED_CDF_INV inverts the Deranged CDF.
!
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, integer A, the number of items.
!    1 <= A.
!
!    Output, integer X, the corresponding argument.
!
  integer a
  real cdf
  real cdf2
  real pdf
  integer x
  integer x2
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DERANGED_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.0

  do x2 = 0, a

    call deranged_pdf ( x2, a, pdf )

    cdf2 = cdf2 + pdf

    if ( cdf <= cdf2 ) then
      x = x2
      return
    end if

  end do

  x = a

  return
end
subroutine deranged_check ( a )
!
!*******************************************************************************
!
!! DERANGED_CHECK checks the parameter of the Deranged PDF.
!
!
!  Modified:
!
!    18 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the total number of items.
!    1 <= A.
!
  integer a
!
  if ( a < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DERANGED_CHECK - Fatal error!'
    write ( *, * ) '  A < 1.'
    stop
  end if

  return
end
function deranged_enum ( n )
!
!*******************************************************************************
!
!! DERANGED_ENUM returns the number of derangements of N objects.
!
!
!  Definition:
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!  Recursion:
!
!      D(0) = 1
!      D(1) = 0
!      D(2) = 1
!      D(N) = (N-1) * ( D(N-1) + D(N-2) )
!
!    or
!
!      D(0) = 1
!      D(1) = 0
!      D(N) = N * D(N-1) + (-1)**N
!
!  Formula:
!
!    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
!
!    Based on the inclusion/exclusion law.
!
!  Comments:
!
!    D(N) is the number of ways of placing N non-attacking rooks on 
!    an N by N chessboard with one diagonal deleted.
!
!    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
!
!    The number of permutations with exactly K items in the right
!    place is COMB(N,K) * D(N-K).
!
!  First values:
!
!     N         D(N)
!     0           1
!     1           0
!     2           1
!     3           2
!     4           9
!     5          44
!     6         265
!     7        1854
!     8       14833
!     9      133496
!    10     1334961
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Output, integer DERANGED_ENUM, the number of derangements of N objects.
!
  integer deranged_enum
  integer dn
  integer dnm1
  integer dnm2
  integer i
  integer n
!
  if ( n < 0 ) then
  
    dn = 0

  else if ( n == 0 ) then

    dn = 1

  else if ( n == 1 ) then

    dn = 0

  else if ( n == 2 ) then

    dn = 1

  else
  
    dnm1 = 0
    dn = 1
    
    do i = 3, n
      dnm2 = dnm1
      dnm1 = dn
      dn = ( i - 1 ) * ( dnm1 + dnm2 )
    end do
            
  end if
  
  deranged_enum = dn

  return
end
subroutine deranged_mean ( a, mean )
!
!*******************************************************************************
!
!! DERANGED_MEAN returns the mean of the Deranged CDF.
!
!
!  Discussion:
!
!    The mean is computed by straightforward summation.
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of items.
!    1 <= A.
!
!    Output, real MEAN, the mean of the PDF.
!
  integer a
  real mean
  real pdf
  integer x
!
  mean = 0.0
  do x = 1, a
    call deranged_pdf ( x, a, pdf )
    mean = mean + pdf * x
  end do

  return
end
subroutine deranged_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! DERANGED_PDF evaluates the Deranged PDF.
!
!
!  Definition:
!
!    PDF(X)(A) is the probability that exactly X items will occur in
!    their proper place after a random permutation of A items.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the number of items in their correct places.
!    0 <= X <= A.
!
!    Input, integer A, the total number of items.
!    1 <= A.
!
!    Output, real PDF, the value of the PDF.
!
  integer a
  integer cnk
  integer deranged_enum
  integer dnmk
  integer i_factorial
  integer nfact
  real pdf
  integer x
!
  if ( x < 0 .or. x > a ) then
    pdf = 0.0
  else
    call binomial_coef ( a, x, cnk )
    dnmk = deranged_enum ( a-x )
    nfact = i_factorial ( a )
    pdf = real ( cnk * dnmk ) / real ( nfact )
  end if

  return
end
subroutine deranged_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! DERANGED_SAMPLE samples the Deranged PDF.
!
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of items.
!    1 <= A.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  integer a
  real cdf
  integer iseed
  real uniform_01_sample
  integer x
!
  cdf = uniform_01_sample ( iseed )

  call deranged_cdf_inv ( cdf, a, x )

  return
end
subroutine deranged_variance ( a, variance )
!
!*******************************************************************************
!
!! DERANGED_VARIANCE returns the variance of the Deranged CDF.
!
!
!  Discussion:
!
!    The variance is computed by straightforward summation.
!
!  Modified:
!
!    06 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of items.
!    1 <= A.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer a
  real mean
  real pdf
  integer x
  real variance
!
  call deranged_mean ( a, mean )

  variance = 0.0
  do x = 1, a
    call deranged_pdf ( x, a, pdf )
    variance = variance + pdf * ( x - mean )**2
  end do

  return
end
function digamma ( x )
!
!*******************************************************************************
!
!! DIGAMMA calculates the digamma or Psi function = d ( LOG ( GAMMA ( X ) ) ) / dX
!
!
!  Reference:
!
!    J Bernardo,
!    Psi ( Digamma ) Function,
!    Algorithm AS 103,
!    Applied Statistics,
!    Volume 25, Number 3, pages 315-317, 1976.
!
!  Modified:
!
!    03 January 2000
!
!  Parameters:
!
!    Input, real X, the argument of the digamma function.
!    0 < X.
!
!    Output, real DIGAMMA, the value of the digamma function at X.
!
  real, parameter :: c = 8.5
  real, parameter :: d1 = -0.5772156649
  real, parameter :: s = 0.00001
  real, parameter :: s3 = 0.08333333333
  real, parameter :: s4 = 0.0083333333333
  real, parameter :: s5 = 0.003968253968
!
  real digamma
  real r
  real x
  real y
!
!  The argument must be positive.
!
  if ( x <= 0.0 ) then

    digamma = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'DIGAMMA - Fatal error!'
    write ( *, * ) '  X <= 0.'
    stop
!
!  Use approximation if argument <= S.
!
  else if ( x <= s ) then

    digamma = d1 - 1.0 / x
!
!  Reduce the argument to DIGAMMA(X + N) where (X + N) >= C.
!
  else

    digamma = 0.0
    y = x

    do while ( y < c )
      digamma = digamma - 1.0 / y
      y = y + 1.0
    end do
!
!  Use Stirling's (actually de Moivre's) expansion if argument > C.
!
    r = 1.0 / y**2
    digamma = digamma + log ( y ) - 0.5 / y - r * ( s3 - r * ( s4 - r * s5 ) )

  end if

  return
end
subroutine dipole_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! DIPOLE_CDF evaluates the Dipole CDF.
!
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
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI 
!      is interesting,
!    and -1.0 <= B <= 1.0.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  real x
!
  cdf = 0.5 + ( 1.0 / PI ) * atan ( x ) + b**2 * ( x * cos ( 2.0 * a ) &
    - sin ( 2.0 * a ) ) / ( PI * ( 1.0 + x**2 ) )

  return
end
subroutine dipole_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! DIPOLE_CDF_INV inverts the Dipole CDF.
!
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Modified:
!
!    04 January 2000
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    -1.0 <= B <= 1.0.
!
!    Output, real X, the corresponding argument of the CDF.
!
  integer, parameter :: IT_MAX = 100
  real, parameter :: TOL = 0.0001
!
  real a
  real b
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
!
!  Take care of horrible input.
!
  if ( cdf <= 0.0 ) then
    x = - huge ( x )
    return
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
    return
  end if
!
!  Seek X1 < X < X2.
!
  x1 = - 1.0

  do

    call dipole_cdf ( x1, a, b, cdf1 ) 

    if ( cdf1 <= cdf ) then
      exit
    end if

    x1 = 2.0 * x1

  end do

  x2 = 1.0

  do

    call dipole_cdf ( x2, a, b, cdf2 )

    if ( cdf2 >= cdf ) then
      exit
    end if

    x2 = 2.0 * x2

  end do
!
!  Now use bisection.
!
  it = 0

30    continue

  it = it + 1

  x3 = 0.5 * ( x1 + x2 )
  call dipole_cdf ( x3, a, b, cdf3 )

  if ( abs ( cdf3 - cdf ) < TOL ) then
    x = x3
    return
  end if

  if ( it > IT_MAX ) then
    write ( *, * ) ' '
    write ( *, * ) 'DIPOLE_CDF_INV - Fatal error!'
    write ( *, * ) '  Iteration limit exceeded.'
    stop
  end if

  if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
    x1 = x3
    cdf1 = cdf3
  else
    x2 = x3
    cdf2 = cdf3
  end if

  go to 30
end
subroutine dipole_check ( a, b )
!
!*******************************************************************************
!
!! DIPOLE_CHECK checks the parameters of the Dipole CDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI 
!      is interesting,
!    and -1.0 <= B <= 1.0.
!
  real a
  real b
!
  if ( b < -1.0 .or. b > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DIPOLE_CHECK - Fatal error!'
    write ( *, * ) '  -1.0 <= B <= 1.0 is required.'
    write ( *, * ) '  The input B = ', b
    stop
  end if

  return
end
subroutine dipole_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! DIPOLE_PDF evaluates the Dipole PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = 
!        1 / ( PI * ( 1 + X**2 ) )
!      + B**2 * ( ( 1 - X**2 ) * cos ( 2 * A ) + 2.0 * X * sin ( 2 * A ) )
!      / ( PI * ( 1 + X )**2 )
!
!  Discussion:
!
!    Densities of this kind commonly occur in the analysis of resonant
!    scattering of elementary particles.
!
!    DIPOLE_PDF(X)(A,0) = CAUCHY_PDF(X)(A)
!    A = 0, B = 1 yields the single channel dipole distribution.
!
!  Reference:
!
!    Robert Knop,
!    Algorithm 441,
!    ACM Transactions on Mathematical Software.
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
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI 
!      is interesting,
!    and -1.0 <= B <= 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
!
  pdf = 1.0 / ( PI * ( 1.0 + x**2 ) ) &
    + b**2 * ( ( 1.0 - x**2 ) * cos ( 2.0 * a ) &
    + 2.0 * x * sin ( 2.0 * x ) ) / ( PI * ( 1.0 + x**2 )**2 )

  return
end
subroutine dipole_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! DIPOLE_SAMPLE samples the Dipole PDF.
!
!
!  Reference:
!
!    Robert Knop,
!    Algorithm 441,
!    ACM Transactions on Mathematical Software.
!
!  Modified:
!
!    04 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI 
!      is interesting,
!    and -1.0 <= B <= 1.0.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real a2
  real b
  real b2
  real c2
  integer iseed
  real x
  real x1
  real x2
!
!  Find (X1,X2) at random in a circle.
!
  a2 = b * sin ( a )
  b2 = b * cos ( a )
  c2 = 1.0

  call circle_sample ( a2, b2, c2, iseed, x1, x2 )
!
!  The dipole variate is the ratio X1 / X2.
!
  x = x1 / x2

  return
end
subroutine dirichlet_check ( n, a )
!
!*******************************************************************************
!
!! DIRICHLET_CHECK checks the parameters of the Dirichlet PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
  integer n
!
  real a(n)
  integer i
  logical positive
!
  positive = .false.

  do i = 1, n

    if ( a(i) < 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'DIRICHLET_CHECK - Fatal error!'
      write ( *, * ) '  A(I) < 0.'
      write ( *, * ) '  For I = ', i
      write ( *, * ) '  A(I) = ', a(i)
      stop
    else if ( a(i) > 0.0 ) then
      positive = .true.
    end if

  end do

  if ( .not. positive ) then
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_CHECK - Fatal error!'
    write ( *, * ) '  All parameters are zero!'
    stop
  end if

  return
end
subroutine dirichlet_mean ( n, a, mean )
!
!*******************************************************************************
!
!! DIRICHLET_MEAN returns the means of the Dirichlet PDF.
!
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real MEAN(N), the means of the PDF.
!
  integer n
!
  real a(n)
  real mean(n)
!
  mean(1:n) = a(1:n)

  call rvec_unit_sum ( n, mean )

  return
end
subroutine dirichlet_mix_check ( comp_num, elem_max, elem_num, a, comp_weight )
!
!*******************************************************************************
!
!! DIRICHLET_MIX_CHECK checks the parameters of a Dirichlet mixture PDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_MAX, the leading dimension of A, which must
!    be at least ELEM_NUM.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(ELEM_MAX,COMP_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
  integer comp_num
  integer elem_max
  integer elem_num
!
  real a(elem_max,comp_num)
  integer comp_i
  real comp_weight(comp_num)
  integer elem_i
  logical positive
!
  do comp_i = 1, comp_num

    do elem_i = 1, elem_num
      if ( a(elem_i,comp_i) < 0.0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'DIRICHLET_MIX_CHECK - Fatal error!'
        write ( *, * ) '  A(ELEM,COMP) < 0.'
        write ( *, * ) '  COMP = ', comp_i
        write ( *, * ) '  ELEM = ', elem_i
        write ( *, * ) '  A(COMP,ELEM) = ', a(elem_i,comp_i)
        stop
      end if
    end do

  end do

  positive = .false.

  do comp_i = 1, comp_num

    if ( comp_weight(comp_i) < 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'DIRICHLET_MIX_CHECK - Fatal error!'
      write ( *, * ) '  COMP_WEIGHT(COMP) < 0.'
      write ( *, * ) '  COMP = ', comp_i
      write ( *, * ) '  COMP_WEIGHT(COMP) = ', comp_weight(comp_i)
      stop
    else if ( comp_weight(comp_i) > 0.0 ) then
      positive = .true.
    end if

  end do

  if ( .not. positive ) then
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_MIX_CHECK - Fatal error!'
    write ( *, * ) '  All component weights are zero.'
    stop
  end if

  return
end
subroutine dirichlet_mix_mean ( comp_num, elem_max, elem_num, a, comp_weight, &
  mean )
!
!*******************************************************************************
!
!! DIRICHLET_MIX_MEAN returns the means of a Dirichlet mixture PDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_MAX, the leading dimension of A, which must
!    be at least ELEM_NUM.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(ELEM_MAX,COMP_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
!    Output, real MEAN(ELEM_NUM), the means for each element.
!
  integer comp_num
  integer elem_max
  integer elem_num
!
  real a(elem_max,comp_num)
  integer comp_i
  real comp_mean(elem_num)
  real comp_weight(comp_num)
  real comp_weight_sum
  integer elem_i
  real mean(elem_num)
!
  comp_weight_sum = sum ( comp_weight )

  mean(1:elem_num) = 0.0

  do comp_i = 1, comp_num
    call dirichlet_mean ( elem_num, a(1,comp_i), comp_mean )
    mean(1:elem_num) = mean(1:elem_num) &
      + comp_weight(comp_i) * comp_mean(1:elem_num)
  end do

  mean(1:elem_num) = mean(1:elem_num) / comp_weight_sum

  return
end
subroutine dirichlet_mix_pdf ( x, comp_num, elem_max, elem_num, a, &
  comp_weight, pdf )
!
!*******************************************************************************
!
!! DIRICHLET_MIX_PDF evaluates a Dirichlet mixture PDF.
!
!
!  Discussion:
!
!    The PDF is a weighted sum of Dirichlet PDF's.  Each PDF is a 
!    "component", with an associated weight.
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X(ELEM_NUM), the argument of the PDF.
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_MAX, the leading dimension of A, which must
!    be at least ELEM_NUM.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(ELEM_MAX,COMP_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
!    Output, real PDF, the value of the PDF.
!
  integer comp_num
  integer elem_max
  integer elem_num
!
  real a(elem_max,comp_num)
  integer comp_i
  real comp_pdf
  real comp_weight(comp_num)
  real comp_weight_sum
  real pdf
  real x(elem_num)
!
  comp_weight_sum = sum ( comp_weight )

  pdf = 0.0
  do comp_i = 1, comp_num

    call dirichlet_pdf ( x, elem_num, a(1,comp_i), comp_pdf )

    pdf = pdf + comp_weight(comp_i) * comp_pdf / comp_weight_sum

  end do

  return
end
subroutine dirichlet_mix_sample ( comp_num, elem_max, elem_num, a, &
  comp_weight, iseed, comp, x )
!
!*******************************************************************************
!
!! DIRICHLET_MIX_SAMPLE samples a Dirichlet mixture PDF.
!
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COMP_NUM, the number of components in the Dirichlet
!    mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    Input, integer ELEM_MAX, the leading dimension of A, which must
!    be at least ELEM_NUM.
!
!    Input, integer ELEM_NUM, the number of elements of an observation.
!
!    Input, real A(ELEM_MAX,COMP_NUM), the probabilities for element ELEM_NUM
!    in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    Input, integer COMP_WEIGHT(COMP_NUM), the mixture weights of the densities.
!    These do not need to be normalized.  The weight of a given component is
!    the relative probability that that component will be used to generate
!    the sample.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer COMP, the index of the component of the Dirichlet
!    mixture that was chosen to generate the sample.
!
!    Output, real X(ELEM_NUM), a sample of the PDF.
!
  integer comp_num
  integer elem_max
  integer elem_num
!
  real a(elem_max,comp_num)
  integer comp
  real comp_weight(comp_num)
  integer iseed
  real x(elem_num)
!
!  Choose a particular density component COMP.
!
  call discrete_sample ( comp_num, comp_weight, iseed, comp )
!
!  Sample the density number COMP.
!
  call dirichlet_sample ( elem_num, a(1,comp), iseed, x )

  return
end
subroutine dirichlet_moment2 ( n, a, m2 )
!
!*******************************************************************************
!
!! DIRICHLET_MOMENT2 returns the second moments of the Dirichlet PDF.
!
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real M2(N,N), the second moments of the PDF.
!
  integer n
!
  real a(n)
  real a_sum
  real m2(n,n)
  integer i
  integer j
!
  a_sum = sum ( a )

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        m2(i,j) = a(i) * ( a(i) + 1.0 ) / ( a_sum * ( a_sum + 1.0 ) )
      else
        m2(i,j) = a(i) * a(j) / ( a_sum * ( a_sum + 1.0 ) )
      end if
    end do
  end do

  return
end
subroutine dirichlet_multinomial_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! DIRICHLET_MULTINOMIAL_PDF evaluates a Dirichlet Multinomial PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = Comb(A,B,X) * ( Gamma(C_Sum) / Gamma(C_Sum+A) )
!      Product ( 1 <= I <= B ) Gamma(C(I)+X(I)) / Gamma(C(I))
!
!    where:
!
!      Comb(A,B,X) is the multinomial coefficient C( A; X(1), X(2), ..., X(B) ),
!      C_Sum = Sum ( 1 <= I <= B ) C(I)
!
!  Reference:
!
!    Kenneth Lange,
!    Mathematical and Statistical Methods for Genetic Analysis,
!    Springer, 1997, page 45.
!
!  Modified:
!
!    17 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X(B); X(I) counts the number of occurrences of
!    outcome I, out of the total of A trials.
!
!    Input, integer A, the total number of trials.
!
!    Input, integer B, the number of different possible outcomes on
!    one trial.
!
!    Input, integer C(B); C(I) is the Dirichlet parameter associated
!    with outcome I.
!
!    Output, real PDF, the value of the Dirichlet multinomial PDF.
!
  integer b
!
  integer a
  real c(b)
  real c_sum
  real gamma_log
  integer i
  real pdf
  real pdf_log
  integer x(b)
!
  c_sum = sum ( c )

  pdf_log = - gamma_log ( c_sum + real ( a ) ) + gamma_log ( c_sum ) &
            + gamma_log ( real ( a + 1 ) )

  do i = 1, b
    pdf_log = pdf_log + gamma_log ( c(i) + real ( x(i) ) ) &
      - gamma_log ( c(i) ) - gamma_log ( real ( x(i) + 1 ) )
  end do

  pdf = exp ( pdf_log )

  return
end
subroutine dirichlet_pdf ( x, n, a, pdf )
!
!*******************************************************************************
!
!! DIRICHLET_PDF evaluates the Dirichlet PDF.
!
!
!  Definition:
!
!    PDF(X)(N,A) = Product ( 1 <= I <= N ) X(I)**( A(I) - 1 ) 
!      * Gamma ( A_SUM ) / A_PROD
!
!    where 
!
!      0 <= A(I) for all I;
!      0 <= X(I) for all I;
!      Sum ( 1 <= I <= N ) X(I) = 1;
!      A_SUM = Sum ( 1 <= I <= N ) A(I).
!      A_PROD = Product ( 1 <= I <= N ) Gamma ( A(I) )
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
!    Input, real X(N), the argument of the PDF.  Each X(I) should
!    be greater than 0.0, and the X(I)'s must add up to 1.0.
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be
!    positive.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: TOL = 0.0001
!
  integer n
!
  real a(n)
  real a_prod
  real a_sum
  real gamma
  integer i
  real pdf
  real x(n)
  real x_sum
!
  do i = 1, n
    if ( x(i) <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'DIRICHLET_PDF - Fatal error!'
      write ( *, * ) '  X(I) <= 0.'
    end if
  end do

  x_sum = sum ( x )

  if ( abs ( x_sum - 1.0 ) > TOL ) then
    write ( *, * ) ' '
    write ( *, * ) 'DIRICHLET_PDF - Fatal error!'
    write ( *, * ) '  SUM X(I) =/= 1.'
  end if

  a_sum = sum ( a )

  a_prod = 1.0
  do i = 1, n
    a_prod = a_prod * gamma ( a(i) )
  end do

  pdf = gamma ( a_sum ) / a_prod
  do i = 1, n
    pdf = pdf * x(i)**( a(i) - 1.0 )
  end do

  return
end
subroutine dirichlet_sample ( n, a, iseed, x )
!
!*******************************************************************************
!
!! DIRICHLET_SAMPLE samples the Dirichlet PDF.
!
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 169.
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be
!    positive.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X(N), a sample of the PDF.  The entries of X should
!    sum to 1.
!
  integer n
!
  real a(n)
  real a2
  real b2
  real c2
  integer i
  integer iseed
  real x(n)
!
  a2 = 0.0
  b2 = 1.0

  do i = 1, n
    c2 = a(i)
    call gamma_sample ( a2, b2, c2, iseed, x(i) )
  end do
!
!  Rescale the vector to have unit sum.
!
  call rvec_unit_sum ( n, x )

  return
end
subroutine dirichlet_variance ( n, a, variance )
!
!*******************************************************************************
!
!! DIRICHLET_VARIANCE returns the variances of the Dirichlet PDF.
!
!
!  Modified:
!
!    03 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, real A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!    Output, real VARIANCE(N), the variances of the PDF.
!
  integer n
!
  real a(n)
  real a_sum
  integer i
  real variance(n)
!
  a_sum = sum ( a )

  do i = 1, n
    variance(i) = a(i) * ( a_sum - a(i) ) / ( a_sum**2 * ( a_sum + 1.0 ) )
  end do

  return
end
subroutine discrete_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! DISCRETE_CDF evaluates the Discrete CDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the item whose probability is desired.
!
!    Input, integer A, the number of probabilities assigned.
!
!    Input, real B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
!    Output, real CDF, the value of the CDF.
!
  integer a
!
  real b(a)
  real cdf
  integer x
!
  if ( x < 1 ) then
    cdf = 0.0
  else if ( x < a ) then
    cdf = sum ( b(1:x) ) / sum ( b(1:a) )
  else if ( x >= a ) then
    cdf = 1.0
  end if

  return
end
subroutine discrete_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! DISCRETE_CDF_INV inverts the Discrete CDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, integer A, the number of probabilities assigned.
!
!    Input, real B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
!    Output, integer X, the corresponding argument for which
!    CDF(X-1) < CDF <= CDF(X)
!
  integer a
!
  real b(a)
  real b_sum
  real cdf
  real cum
  integer j
  integer x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DISCRETE_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  b_sum = sum ( b )

  cum = 0.0

  do j = 1, a

    cum = cum + b(j) / b_sum

    if ( cdf <= cum ) then
      x = j
      return
    end if

  end do

  x = a
  
  return
end
subroutine discrete_check ( a, b )
!
!*******************************************************************************
!
!! DISCRETE_CHECK checks the parameters of the Discrete CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of probabilities assigned.
!
!    Input, real B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
  integer a
!
  real b(a)
  real b_sum
  integer j
!
  do j = 1, a
    if ( b(j) < 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'DISCRETE_CHECK - Fatal error!'
      write ( *, * ) '  Negative probabilities not allowed.'
      stop
    end if
  end do

  b_sum = sum ( b )

  if ( b_sum == 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DISCRETE_CHECK - Fatal error!'
    write ( *, * ) '  Total probablity is zero.'
    stop
  end if
  
  return
end
subroutine discrete_mean ( a, b, mean )
!
!*******************************************************************************
!
!! DISCRETE_MEAN evaluates the mean of the Discrete PDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of probabilities assigned.
!
!    Input, real B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
!    Output, real MEAN, the mean of the PDF.
!
  integer a
!
  real b(a)
  real b_sum
  integer j
  real mean
!
  b_sum = sum ( b )

  mean = 0.0
  do j = 1, a
    mean = mean + real ( j ) * b(j)
  end do

  mean = mean / b_sum

  return
end
subroutine discrete_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! DISCRETE_PDF evaluates the Discrete PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = B(X) if 1 <= X <= A
!                = 0    otherwise
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the item whose probability is desired.
!
!    Input, integer A, the number of probabilities assigned.
!
!    Input, real B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
!    Output, real PDF, the value of the PDF.
!
  integer a
!
  real b(a)
  real b_sum
  real pdf
  integer x
!
  b_sum = sum ( b )

  if ( 1 <= x .and. x <= a ) then
    pdf = b(x) / b_sum
  else
    pdf = 0.0
  end if

  return
end
subroutine discrete_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! DISCRETE_SAMPLE samples the Discrete PDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of probabilities assigned.
!
!    Input, real B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  integer a
!
  real b(a)
  real b_sum
  real cdf
  integer iseed
  real uniform_01_sample
  integer x
!
  b_sum = sum ( b )

  cdf = uniform_01_sample ( iseed )

  call discrete_cdf_inv ( cdf, a, b, x )

  return
end
subroutine discrete_variance ( a, b, variance )
!
!*******************************************************************************
!
!! DISCRETE_VARIANCE evaluates the variance of the Discrete PDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of probabilities assigned.
!
!    Input, real B(A), the relative probabilities of outcomes 1 through A.
!    Each entry must be nonnegative.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer a
!
  real b(a)
  real b_sum
  integer j
  real mean
  real variance
!
  b_sum = sum ( b )

  mean = 0.0
  do j = 1, a
    mean = mean + real ( j ) * b(j)
  end do

  mean = mean / b_sum

  variance = 0.0
  do j = 1, a
    variance = variance + b(j) * ( j - mean )**2 
  end do

  variance = variance / b_sum

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
!   "E" was named in honor of Euler.
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
!    Output, real E_CONSTANT, the base of the natural logarithm system.
!
  real e_constant
!
  e_constant = 2.71828182845904523536028747135266249775724709369995

  return
end
subroutine erlang_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! ERLANG_CDF evaluates the Erlang CDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, integer C, the parameters of the PDF.
!    0.0 < B.
!    0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  integer c
  real cdf
  real gamma_inc
  real p2
  real x
  real x2
!
  if ( x < a ) then

    cdf = 0.0

  else

    x2 = ( x - a ) / b
    p2 = real ( c )

    cdf = gamma_inc ( x2, p2 )

  end if

  return
end
subroutine erlang_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! ERLANG_CDF_INV inverts the Erlang CDF.
!
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Modified:
!
!    05 January 2000
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, integer C, the parameters of the PDF.
!    0.0 < B.
!    0 < C.
!
!    Output, real X, the corresponding argument of the CDF.
!
  integer, parameter :: IT_MAX = 10
  real, parameter :: TOL = 0.0001
!
  real a
  real b
  integer c
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
!
  if ( cdf <= 0.0 ) then
    x = a
    return
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
    return
  end if

  x1 = a
  cdf1 = 0.0

  x2 = a + 1.0

  do

    call erlang_cdf ( x2, a, b, c, cdf2 )

    if ( cdf2 > cdf ) then
      exit
    end if

    x2 = a + 2.0 * ( x2 - a )

  end do
!
!  Now use bisection.
!
  it = 0

20    continue

  it = it + 1

  x3 = 0.5 * ( x1 + x2 )
  call erlang_cdf ( x3, a, b, c, cdf3 )

  if ( abs ( cdf3 - cdf ) < TOL ) then
    x = x3
    return
  end if

  if ( it > IT_MAX ) then
    write ( *, * ) ' '
    write ( *, * ) 'ERLANG_CDF_INV - Fatal error!'
    write ( *, * ) '  Iteration limit exceeded.'
    stop
  end if

  if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
    x1 = x3
    cdf1 = cdf3
  else
    x2 = x3
    cdf2 = cdf3
  end if

  go to 20
end
subroutine erlang_check ( a, b, c )
!
!*******************************************************************************
!
!! ERLANG_CHECK checks the parameters of the Erlang PDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, integer C, the parameters of the PDF.
!    0.0 < B.
!    0 < C.
!
  real a
  real b
  integer c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ERLANG_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.0'
    stop
  end if

  if ( c <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ERLANG_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    stop
  end if

  return
end
subroutine erlang_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! ERLANG_MEAN returns the mean of the Erlang PDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, integer C, the parameters of the PDF.
!    0.0 < B.
!    0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  integer c
  real mean
!
  mean =  a + b * real ( c )

  return
end
subroutine erlang_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! ERLANG_PDF evaluates the Erlang PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = ( ( X - A ) / B )**( C - 1 ) 
!      / ( B * Gamma ( C ) * EXP ( ( X - A ) / B ) )
!
!    for B > 0, C integer > 0, X >= A.
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, integer C, the parameters of the PDF.
!    0.0 < B.
!    0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  integer c
  real pdf
  real r_factorial
  real x
  real y
!
  if ( x <= a ) then

    pdf = 0.0

  else

    y = ( x - a ) / b

    pdf = y**( c - 1 ) / ( b * r_factorial ( c - 1 ) * exp ( y ) )

  end if

  return
end
subroutine erlang_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! ERLANG_SAMPLE samples the Erlang PDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, integer C, the parameters of the PDF.
!    0.0 < B.
!    0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real a2
  real b
  real b2
  integer c
  integer i
  integer iseed
  real x
  real x2
!
  a2 = 0.0
  b2 = b
  x = a
  do i = 1, c
    call exponential_sample ( a2, b2, iseed, x2 )
    x = x + x2
  end do

  return
end
subroutine erlang_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! ERLANG_VARIANCE returns the variance of the Erlang PDF.
!
!
!  Modified:
!
!    31 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, integer C, the parameters of the PDF.
!    0.0 < B.
!    0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  integer c
  real variance
!
  variance =  b**2 * real ( c )

  return
end
function error_function ( x )
!
!*******************************************************************************
!
!! ERROR_FUNCTION evaluates the error function ERF(X).
!
!
!  Definition:
!
!    ERF(X) = ( 2 / SQRT ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T**2 ) dT.
!
!  Properties:
!
!    Limit ( X -> -Infinity ) ERF(X) =          -1.0;
!                             ERF(0) =           0.0;
!                             ERF(0.476936...) = 0.5;
!    Limit ( X -> +Infinity ) ERF(X) =          +1.0.
!
!    0.5 * ( ERF(X/SQRT(2)) + 1 ) = Normal_CDF(X)
!
!  Reference:
!
!    W J Cody,
!    "Rational Chebyshev Approximations for the Error Function",
!    Math. Comp., 1969, pages 631-638.
!
!  Modified:
!
!    06 December 1999
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
!    Input, real X, the argument of the error function.
!
!    Output, real ERROR_FUNCTION, the value of the error function.
!
  double precision, parameter :: SQRPI = 0.56418958354775628695D+00
  double precision, parameter :: THRESH = 0.46875D+00
  double precision, parameter :: XBIG = 26.543D+00
  double precision, parameter :: XSMALL = 1.11D-16
!
  double precision, parameter, dimension ( 5 ) :: a = (/ &
    3.16112374387056560D+00, &
    1.13864154151050156D+02, &
    3.77485237685302021D+02, &
    3.20937758913846947D+03, &
    1.85777706184603153D-01 /)
  double precision, parameter, dimension ( 4 ) :: b = (/ &
    2.36012909523441209D+01, &
    2.44024637934444173D+02, &
    1.28261652607737228D+03, &
    2.84423683343917062D+03 /)
  double precision, parameter, dimension ( 9 ) :: c = (/ &
    5.64188496988670089D-01, &
    8.88314979438837594D+00, &
    6.61191906371416295D+01, &
    2.98635138197400131D+02, &
    8.81952221241769090D+02, &
    1.71204761263407058D+03, &
    2.05107837782607147D+03, &
    1.23033935479799725D+03, &
    2.15311535474403846D-08 /)
  double precision, parameter, dimension ( 8 ) :: d = (/ &
    1.57449261107098347D+01, &
    1.17693950891312499D+02, &
    5.37181101862009858D+02, &
    1.62138957456669019D+03, &
    3.29079923573345963D+03, &
    4.36261909014324716D+03, &
    3.43936767414372164D+03, &
    1.23033935480374942D+03 /)
  double precision del
  real error_function
  double precision erfxd
  integer i
  double precision, parameter, dimension ( 6 ) :: p = (/ &
    3.05326634961232344D-01, &
    3.60344899949804439D-01, &
    1.25781726111229246D-01, &
    1.60837851487422766D-02, &
    6.58749161529837803D-04, &
    1.63153871373020978D-02 /)
  double precision, parameter, dimension ( 5 ) :: q = (/ &
    2.56852019228982242D+00, &
    1.87295284992346047D+00, &
    5.27905102951428412D-01, &
    6.05183413124413191D-02, &
    2.33520497626869185D-03 /)
  real x
  double precision xabs
  double precision xden
  double precision xnum
  double precision xsq
!
  xabs = abs ( dble ( x ) )
!
!  Evaluate ERF(X) for |X| <= 0.46875.
!
  if ( xabs <= THRESH ) then

    if ( xabs > XSMALL ) then
      xsq = xabs**2
    else
      xsq = 0.0
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do

    erfxd = dble ( x ) * ( xnum + a(4) ) / ( xden + b(4) )
!
!  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
!
  else if ( xabs <= 4.0 ) then

    xnum = c(9) * xabs
    xden = xabs
    do i = 1, 7
      xnum = ( xnum + c(i) ) * xabs
      xden = ( xden + d(i) ) * xabs
    end do

    erfxd = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( xabs * 16.0 ) / 16.0
    del = ( xabs - xsq ) * ( xabs + xsq )
    erfxd = exp ( - xsq**2 ) * exp ( - del ) * erfxd

    erfxd = ( 0.5 - erfxd ) + 0.5

    if ( x < 0.0 ) then
      erfxd = - erfxd
    end if
!
!  Evaluate ERFC(X) for 4.0 < |X|.
!
  else

    if ( xabs >= XBIG ) then

      if ( x > 0.0 ) then
        erfxd = 1.0
      else
        erfxd = - 1.0
      end if

    else

      xsq = 1.0 / xabs**2

      xnum = p(6) * xsq
      xden = xsq
      do i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
      end do

      erfxd = xsq * ( xnum + p(5) ) / ( xden + q(5) )
      erfxd = ( SQRPI - erfxd ) / xabs
      xsq = aint ( xabs * 16.0 ) / 16.0
      del = ( xabs - xsq ) * ( xabs + xsq )
      erfxd = exp ( - xsq**2 ) * exp ( - del ) * erfxd

      erfxd = ( 0.5 - erfxd ) + 0.5

      if ( x < 0.0 ) then
        erfxd = - erfxd
      end if

    end if

  end if

  error_function = real ( erfxd )

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
!    Output, real EULER_CONSTANT, the value of the Euler-Mascheroni 
!    constant.
!
  real euler_constant
!
  euler_constant = 0.577215664901532860606512090082402431042

  return
end
subroutine exponential_01_cdf ( x, cdf )
!
!*******************************************************************************
!
!! EXPONENTIAL_01_CDF evaluates the Exponential 01 CDF.
!
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
!    Input, real X, the argument of the PDF.
!
!    Output, real CDF, the value of the CDF.
!
  real cdf
  real x
!
  if ( x <= 0.0 ) then
    cdf = 0.0
  else
    cdf = 1.0 - exp ( - x ) 
  end if

  return
end
subroutine exponential_01_cdf_inv ( cdf, x )
!
!*******************************************************************************
!
!! EXPONENTIAL_01_CDF_INV inverts the Exponential 01 CDF.
!
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
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, real X, the corresponding argument.
!
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'EXPONENTIAL_01_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = - log ( 1.0 - cdf )

  return
end
subroutine exponential_01_mean ( mean )
!
!*******************************************************************************
!
!! EXPONENTIAL_01_MEAN returns the mean of the Exponential 01 PDF.
!
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
!    Output, real MEAN, the mean of the PDF.
!
  real mean
!
  mean = 1.0

  return
end
subroutine exponential_01_pdf ( x, pdf )
!
!*******************************************************************************
!
!! EXPONENTIAL_01_PDF evaluates the Exponential 01 PDF.
!
!
!  Formula:
!
!    PDF(X) = EXP ( - X )
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
!    Input, real X, the argument of the PDF.
!    0.0 <= X
!
!    Output, real PDF, the value of the PDF.
!
  real pdf
  real x
!
  if ( x < 0.0 ) then
    pdf = 0.0
  else
    pdf = exp ( - x )
  end if

  return
end
subroutine exponential_01_sample ( iseed, x )
!
!*******************************************************************************
!
!! EXPONENTIAL_01_SAMPLE samples the Exponential PDF with parameter 1.
!
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
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call exponential_01_cdf_inv ( cdf, x )

  return
end
subroutine exponential_01_variance ( variance )
!
!*******************************************************************************
!
!! EXPONENTIAL_01_VARIANCE returns the variance of the Exponential 01 PDF.
!
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
!    Output, real VARIANCE, the variance of the PDF.
!
  real variance
!
  variance = 1.0

  return
end
subroutine exponential_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! EXPONENTIAL_CDF evaluates the Exponential CDF.
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
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( x <= a ) then
    cdf = 0.0
  else
    cdf = 1.0 - exp ( ( a - x ) / b ) 
  end if

  return
end
subroutine exponential_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! EXPONENTIAL_CDF_INV inverts the Exponential CDF.
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
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'EXPONENTIAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( 1.0 - cdf )

  return
end
subroutine exponential_check ( a, b )
!
!*******************************************************************************
!
!! EXPONENTIAL_CHECK checks the parameters of the Exponential CDF.
!
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
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'EXPONENTIAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.0'
    stop
  end if

  return
end
subroutine exponential_mean ( a, b, mean )
!
!*******************************************************************************
!
!! EXPONENTIAL_MEAN returns the mean of the Exponential PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a + b

  return
end
subroutine exponential_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! EXPONENTIAL_PDF evaluates the Exponential PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = ( 1 / B ) * EXP ( ( A - X ) / B )
!
!  Discussion:
!
!    The time interval between two Poisson events is a random 
!    variable with the Exponential PDF.  The parameter B is the
!    average interval between events.
!
!    In another context, the Exponential PDF is related to
!    the Boltzmann distribution, which describes the relative 
!    probability of finding a system, which is in thermal equilibrium 
!    at absolute temperature T, in a given state having energy E.  
!    The relative probability is
!
!      Boltzmann_Relative_Probability(E,T) = exp ( - E / ( k * T ) ),
!
!    where k is the Boltzmann constant, 
!
!      k = 1.38 * 10**(-23) joules / degree Kelvin
!
!    and normalization requires a determination of the possible
!    energy states of the system.
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
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  if ( x < a ) then
    pdf = 0.0
  else
    pdf = ( 1.0 / b ) * exp ( ( a - x ) / b )
  end if

  return
end
subroutine exponential_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! EXPONENTIAL_SAMPLE samples the Exponential PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call exponential_cdf_inv ( cdf, a, b, x )

  return
end
subroutine exponential_variance ( a, b, variance )
!
!*******************************************************************************
!
!! EXPONENTIAL_VARIANCE returns the variance of the Exponential PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = b**2

  return
end
subroutine extreme_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! EXTREME_CDF evaluates the Extreme Value CDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
  real y
!
  y = ( x - a ) / b

  cdf = exp ( - exp ( - y ) )

  return
end
subroutine extreme_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! EXTREME_CDF_INV inverts the Extreme Value CDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'EXTREME_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( - log ( cdf ) )

  return
end
subroutine extreme_check ( a, b )
!
!*******************************************************************************
!
!! EXTREME_CHECK checks the parameters of the Extreme Value CDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'EXTREME_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine extreme_mean ( a, b, mean )
!
!*******************************************************************************
!
!! EXTREME_MEAN returns the mean of the Extreme Value PDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real euler_constant
  real mean
!
  mean = a + b * euler_constant ( )

  return
end
subroutine extreme_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! EXTREME_PDF evaluates the Extreme Value PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = 
!      ( 1 / B ) * 
!      EXP ( 
!        ( A - X ) / B - EXP ( ( A - X ) / B  ) 
!      ).
!
!  Discussion:
!
!    The Extreme Value PDF is also known as the Fisher-Tippet PDF
!    and the Log-Weibull PDF.
!
!    The special case A = 0 and B = 1 is the Gumbel PDF.
!
!    The Extreme Value PDF is the limiting distribution for the
!    smallest or largest value in a large sample drawn from
!    any of a great variety of distributions.
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998.
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  pdf = ( 1.0 / b ) * exp ( ( a - x ) / b - exp ( ( a - x ) / b ) )

  return
end
subroutine extreme_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! EXTREME_SAMPLE samples the Extreme Value PDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call extreme_cdf_inv ( cdf, a, b, x )

  return
end
subroutine extreme_variance ( a, b, variance )
!
!*******************************************************************************
!
!! EXTREME_VARIANCE returns the variance of the Extreme Value PDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real variance
!
  variance = PI**2 * b**2 / 6.0

  return
end
subroutine f_central_cdf ( x, m, n, cdf )
!
!*******************************************************************************
!
!! F_CENTRAL_CDF evaluates the F CDF.
!
!
!  Reference:
!
!    Formula 26.5.28
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, integer M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real beta_inc
  real cdf
  integer m
  integer n
  real x
  real xx
!
  if ( m < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_CDF - Fatal error!'
    write ( *, * ) '  M < 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_CDF - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if
!
  if ( x <= 0.0 ) then

    cdf = 0.0

  else

    a = 0.5 * real ( m )
    b = 0.5 * real ( n )
    xx = b / ( b + a * x )

    cdf = 1.0 - beta_inc ( xx, a, b )

  end if

  return
end
subroutine f_central_mean ( m, n, mean )
!
!*******************************************************************************
!
!! F_CENTRAL_MEAN returns the mean of the F PDF.
!
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the mean is not defined unless 3 <= N.
!
!    Output, real MEAN, the mean of the PDF.
!
  integer m
  real mean
  integer n
!
  if ( m < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_MEAN - Fatal error!'
    write ( *, * ) '  M < 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_MEAN - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if

  if ( n < 3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_MEAN - Fatal error!'
    write ( *, * ) '  The mean is not defined for N < 3.'
    stop
  end if
!
  mean = real ( n ) / real ( n - 2 )

  return
end
subroutine f_central_pdf ( x, m, n, pdf )
!
!*******************************************************************************
!
!! F_CENTRAL_PDF evaluates the F PDF.
!
!
!  Formula:
!
!    PDF(X)(M,N) = M**(M/2) * X**((M-2)/2)
!      / ( Beta(M/2,N/2) * N**(M/2) * ( 1 + (M/N) * X )**((M+N)/2)
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 <= X
!
!    Input, integer M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real beta
  real bot1
  real bot2
  integer m
  integer n
  real pdf
  real top
  real x
!
  if ( m < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_PDF - Fatal error!'
    write ( *, * ) '  M < 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_PDF - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if

  if ( x < 0.0 ) then

    pdf = 0.0

  else

    a = real ( m )
    b = real ( n )

    top = sqrt ( a**m * b**n * x**( m - 2 ) )
    bot1 = beta ( a / 2.0, b / 2.0 ) 
    bot2 =  sqrt ( ( b + a * x )**( m + n ) )

    pdf = top / ( bot1 * bot2 )

  end if

  return
end
subroutine f_central_sample ( m, n, iseed, x )
!
!*******************************************************************************
!
!! F_CENTRAL_SAMPLE samples the F PDF.
!
!
!  Modified:
!
!    18 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  integer iseed
  integer m
  integer n
  real x
  real xm
  real xn
!
  if ( m < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_SAMPLE - Fatal error!'
    write ( *, * ) '  M < 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_SAMPLE - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if
!
  a = real ( m )
  call chisquare_central_sample ( a, iseed, xm )

  a = real ( n )
  call chisquare_central_sample ( a, iseed, xn )

  x = real ( n ) * xm / ( real ( m ) * xn )

  return
end
subroutine f_central_variance ( m, n, variance )
!
!*******************************************************************************
!
!! F_CENTRAL_VARIANCE returns the variance of the F PDF.
!
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the variance is not defined unless 5 <= N.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer m
  integer n
  real variance
!
  if ( m < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  M < 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if

  if ( n < 5 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_CENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  The variance is not defined for N < 5.'
    stop
  end if
!
  variance = real ( 2 * n**2 * ( m + n - 2 ) ) / &
    real ( m * ( n - 2 )**2 * ( n - 4 ) )

  return
end
subroutine f_noncentral_mean ( a, m, n, mean )
!
!*******************************************************************************
!
!! F_NONCENTRAL_MEAN returns the mean of the F PDF.
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
!    Input, real A, a parameter of the PDF.
!
!    Input, integer M, N, parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the mean is not defined unless 3 <= N.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  integer m
  real mean
  integer n
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_MEAN - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( m < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_MEAN - Fatal error!'
    write ( *, * ) '  M < 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_MEAN - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if

  if ( n < 3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_MEAN - Fatal error!'
    write ( *, * ) '  The mean is not defined for N < 3.'
    stop
  end if
!
  mean = ( real ( m ) + a ) * real ( n ) / ( real ( m ) * real ( n - 2 ) )

  return
end
subroutine f_noncentral_variance ( a, m, n, variance )
!
!*******************************************************************************
!
!! F_NONCENTRAL_VARIANCE returns the variance of the F PDF.
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
!    Input, real A, a parameter of the PDF.
!
!    Input, integer M, N, parameters of the PDF.
!    1 <= M,
!    1 <= N.
!    Note, however, that the variance is not defined unless 5 <= N.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  integer m
  integer n
  real variance
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( m < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  M < 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  N < 1.'
    stop
  end if

  if ( n < 5 ) then
    write ( *, * ) ' '
    write ( *, * ) 'F_NONCENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  The variance is not defined for N < 5.'
    stop
  end if
!
  variance = ( ( real ( m ) + a )**2 + 2.0 * ( real ( m ) + a ) &
    * real ( n )**2 ) / real ( ( n - 2 ) * ( n - 4 ) * m**2 ) - &
    ( real ( m ) + a )**2 * real ( n )**2 / real ( ( n - 2 )**2 * m**2 )

  return
end
function factorial_log ( n )
!
!*******************************************************************************
!
!! FACTORIAL_LOG returns the logarithm of N!.
!
!
!  Definition:
!
!    N! = Product ( 1 <= I <= N ) I
!
!  Method:
!
!    N! = Gamma(N+1).
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the function.
!    0 <= N.
!
!    Output, real FACTORIAL_LOG, the logarithm of N!.
!
  real factorial_log
  integer i
  integer n
!
  if ( n < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FACTORIAL_LOG - Fatal error!'
    write ( *, * ) '  N < 0.'
    stop
  end if

  factorial_log = 0.0

  do i = 2, n
    factorial_log = factorial_log + log ( real ( i ) )
  end do

  return
end
function factorial_stirling ( n )
!
!*******************************************************************************
!
!! FACTORIAL_STIRLING computes Stirling's approximation to N!.
!
!
!  Definition:
!
!    N! = Product ( 1 <= I <= N ) I
!
!    Stirling ( N ) = sqrt ( 2 * PI * N ) * ( N / E )**N * E**(1/(12*N) )
!
!  Discussion:
!
!    This routine returns the raw approximation for all nonnegative
!    values of N.  If N is less than 0, the value is returned as 0,
!    and if N is 0, the value of 1 is returned.  In all other cases,
!    Stirling's formula is used.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the function.
!
!    Output, real FACTORIAL_STIRLING, an approximation to N!.
!
  real, parameter :: E_NATURAL = 2.71828182845904523536028747135266249775724709369995
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real factorial_stirling
  integer n
  real value
!
  if ( n < 0 ) then

    value = 0.0

  else if ( n == 0 ) then

    value = 1.0

  else

    value = sqrt ( 2.0 * PI * real ( n ) ) * ( real ( n ) / E_NATURAL )**n &
      * exp ( 1.0 / real ( 12 * n ) )

  end if

  factorial_stirling = value

  return
end
subroutine fisk_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! FISK_CDF evaluates the Fisk CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real x
!
  if ( x <= a ) then
    cdf = 0.0
  else
    cdf = 1.0 / ( 1.0 + ( b / ( x - a ) )**c )
  end if

  return
end
subroutine fisk_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! FISK_CDF_INV inverts the Fisk CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real X, the corresponding argument of the CDF.
!
  real a
  real b
  real c
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FISK_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.0 ) then
    x = a
  else if ( cdf < 1.0 ) then
    x = a + b * ( cdf / ( 1.0 - cdf ) )**( 1.0 / c )
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
  end if

  return
end
subroutine fisk_check ( a, b, c )
!
!*******************************************************************************
!
!! FISK_CHECK checks the parameters of the Fisk PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
  real a
  real b
  real c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FISK_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FISK_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    stop
  end if

  return
end
subroutine fisk_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! FISK_MEAN returns the mean of the Fisk PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real c
  real csc
  real mean
!
  if ( c <= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FISK_MEAN - Fatal error!'
    write ( *, * ) '  No mean defined for C <= 1.0'
    stop
  end if

  mean = a + PI * ( b / c ) * csc ( PI / c )

  return
end
subroutine fisk_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! FISK_PDF evaluates the Fisk PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = 
!      ( C / B ) * ( ( X - A ) / B )**( C - 1 ) /
!      ( 1 + ( ( X - A ) / B )**C )**2
!
!  Discussion:
!
!    The Fisk PDF is also known as the Log Logistic PDF.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real pdf
  real x
  real y
!
  if ( x <= a ) then

    pdf = 0.0

  else

    y = ( x - a ) / b

    pdf = ( c / b ) * y**( c - 1.0 ) / ( 1.0 + y**c )**2

  end if

  return
end
subroutine fisk_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! FISK_SAMPLE samples the Fisk PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real c
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call fisk_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine fisk_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! FISK_VARIANCE returns the variance of the Fisk PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real c
  real csc
  real g
  real variance
!
  if ( c <= 2.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FISK_VARIANCE - Fatal error!'
    write ( *, * ) '  No variance defined for C <= 2.0'
    stop
  end if

  g = PI / c

  variance = b**2 * ( 2.0 * g * csc ( 2.0 * g ) - ( g * csc ( g ) )**2 )

  return
end
subroutine folded_normal_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! FOLDED_NORMAL_CDF evaluates the Folded Normal CDF.
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
!    Input, real X, the argument of the CDF.
!    0.0 <= X.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A,
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real cdf1
  real cdf2
  real x
  real x1
  real x2
!
  if ( x < 0.0 ) then
    cdf = 0.0
  else
    x1 = ( x - a ) / b
    call normal_01_cdf ( x1, cdf1 )
    x2 = ( - x - a ) / b
    call normal_01_cdf ( x2, cdf2 ) 
    cdf = cdf1 - cdf2
  end if

  return
end
subroutine folded_normal_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! FOLDED_NORMAL_CDF_INV inverts the Folded Normal CDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A,
!    0.0 < B.
!
!    Output, real X, the argument of the CDF.
!    0.0 <= X.
!
  integer, parameter :: IT_MAX = 100
  real, parameter :: TOL = 0.0001
!
  real a
  real b
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
  real xa
  real xb
!
  if ( cdf <= 0.0 ) then
    x = 0.0
    return
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
    return
  end if
!
!  Find X1, for which the value of CDF will be too small.
!
  if ( a >= 0.0 ) then
    call normal_cdf_inv ( cdf, a, b, x1 )
  else
    call normal_cdf_inv ( cdf, -a, b, x1 )
  end if
  x1 = max ( x1, 0.0 )
  call folded_normal_cdf ( x1, a, b, cdf1 )
!
!  Find X2, for which the value of CDF will be too big.
!
  cdf2 = ( 1.0 - cdf ) / 2.0

  call normal_cdf_inv ( cdf2, a, b, xa )
  call normal_cdf_inv ( cdf2, -a, b, xb )
  x2 = max ( abs ( xa ), abs ( xb ) )
  call folded_normal_cdf ( x2, a, b, cdf2 )
!
!  Now use bisection.
!
  it = 0

10    continue

  it = it + 1

  x3 = 0.5 * ( x1 + x2 )
  call folded_normal_cdf ( x3, a, b, cdf3 )

  if ( abs ( cdf3 - cdf ) < TOL ) then
    x = x3
    return
  end if

  if ( it > IT_MAX ) then
    write ( *, * ) ' '
    write ( *, * ) 'FOLDED_NORMAL_CDF_INV - Fatal error!'
    write ( *, * ) '  Iteration limit exceeded.'
    stop
  end if

  if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
    x1 = x3
    cdf1 = cdf3
  else
    x2 = x3
    cdf2 = cdf3
  end if

  go to 10
end
subroutine folded_normal_check ( a, b )
!
!*******************************************************************************
!
!! FOLDED_NORMAL_CHECK checks the parameters of the Folded Normal CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A,
!    0.0 < B.
!
  real a
  real b
!
  if ( a < 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FOLDED_NORMAL_CHECK - Fatal error!'
    write ( *, * ) '  A < 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FOLDED_NORMAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine folded_normal_mean ( a, b, mean )
!
!*******************************************************************************
!
!! FOLDED_NORMAL_MEAN returns the mean of the Folded Normal PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A,
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real a2
  real b
  real cdf
  real mean
!
  a2 = a / b

  call normal_01_cdf ( a2, cdf )

  mean = b * sqrt ( 2.0 / PI ) * exp ( - 0.5 * a2**2 ) - a * ( 1.0 - 2.0 * cdf )

  return
end
subroutine folded_normal_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! FOLDED_NORMAL_PDF evaluates the Folded Normal PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = SQRT ( 2 / PI ) * ( 1 / B ) * COSH ( A * X / B**2 )
!      * EXP ( - 0.5 * ( X**2 + A**2 ) / B**2 )
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
!    Input, real X, the argument of the PDF.
!    0.0 <= X
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A,
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
!
  if ( x < 0.0 ) then
    pdf = 0.0
  else
    pdf = sqrt ( 2.0 / PI ) * ( 1.0 / b ) * cosh ( a * x / b**2 ) &
      * exp ( - 0.5 * ( x**2 + a**2 ) / b**2 )
  end if

  return
end
subroutine folded_normal_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! FOLDED_NORMAL_SAMPLE samples the Folded Normal PDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ISEED, a seed for the random number generator.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A,
!    0.0 < B.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )
  
  call folded_normal_cdf_inv ( cdf, a, b, x )

  return
end
subroutine folded_normal_variance ( a, b, variance )
!
!*******************************************************************************
!
!! FOLDED_NORMAL_VARIANCE returns the variance of the Folded Normal PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A,
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real mean
  real variance
!
  call folded_normal_mean ( a, b, mean )

  variance = a**2 + b**2 - mean**2

  return
end
function gamma ( x )
!
!*******************************************************************************
!
!! GAMMA calculates the Gamma function for a real argument X.
!
!
!  Definition:
!
!    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
!
!  Recursion:
!
!    GAMMA(X+1) = X * GAMMA(X)
!
!  Special values:
!
!    GAMMA(0.5) = SQRT(PI)
!    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the GAMMA
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for X .GE. 12 are from reference 2.
!    The accuracy achieved depends on the arithmetic system, the
!    compiler, the intrinsic functions, and proper selection of the
!    machine-dependent constants.
!
!  Machine-dependent constants:
!
!    BETA: radix for the floating-point representation.
!    MAXEXP: the smallest positive power of BETA that overflows.
!    XBIG: the largest argument for which GAMMA(X) is representable
!      in the machine, i.e., the solution to the equation
!      GAMMA(XBIG) = BETA**MAXEXP.
!    XINF: the largest machine representable floating-point number;
!      approximately BETA**MAXEXP.
!    EPS: the smallest positive floating-point number such that
!      1.0+EPS .GT. 1.0.
!    XMININ: the smallest positive floating-point number such that
!      1/XMININ is machine representable.
!
!    Approximate values for some important machines are:
!
!                               BETA       MAXEXP        XBIG
!
!    CRAY-1         (S.P.)        2         8191        966.961
!    Cyber 180/855
!      under NOS    (S.P.)        2         1070        177.803
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)        2          128        35.040
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)        2         1024        171.624
!    IBM 3033       (D.P.)       16           63        57.574
!    VAX D-Format   (D.P.)        2          127        34.844
!    VAX G-Format   (D.P.)        2         1023        171.489
!
!                               XINF         EPS        XMININ
!
!    CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
!    Cyber 180/855
!      under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
!    IEEE (IBM/XT,
!      SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
!    IEEE (IBM/XT,
!      SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
!    IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
!    VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
!    VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!  Reference: 
!
!    W J Cody,
!    "An Overview of Software Development for Special Functions", 
!    Lecture Notes in Mathematics, 506, 
!    Numerical Analysis Dundee, 1975, 
!    G. A. Watson (ed.),
!    Springer Verlag, Berlin, 1976.
!
!    Hart et al,
!    Computer Approximations, 
!    Wiley and sons, New York, 1968.
!
!  Author: 
!
!    W. J. Cody and L. Stoltz,
!    Applied Mathematics Division,
!    Argonne National Laboratory,
!    Argonne, Illinois, 60439.
!
!  Parameters:
!
!    Input, real X, the argument of the function.
!
!    Output, real GAMMA, the value of the function.  The program 
!    returns the value XINF for singularities or when overflow would occur.  
!    The computation is believed to be free of underflow and overflow.
!
  real EPS
  real PI
  real SQRTPI
  real XBIG
  real XINF
  real XMININ

  real, parameter :: EPS = 1.19e-7
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
  real, parameter :: SQRTPI = 0.9189385332046727417803297
  real, parameter :: XBIG = 35.040
  real, parameter :: XINF = 3.4e38
  real, parameter :: XMININ = 1.18e-38
!
  real, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728E-03, &
     8.4171387781295E-04, &
    -5.952379913043012E-04, &
     7.93650793500350248E-04, &
    -2.777777777777681622553E-03, &
     8.333333333333333331554247E-02, &
     5.7083835261E-03 /)
  real fact
  real gamma
  integer i
  integer n
  real, parameter, dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811e+0, &
     2.47656508055759199108314e+1, &
    -3.79804256470945635097577e+2, &
     6.29331155312818442661052e+2, &
     8.66966202790413211295064e+2, &
    -3.14512729688483675254357e+4, &
    -3.61444134186911729807069e+4, &
     6.64561438202405440627855e+4 /)
  logical parity
  real, parameter, dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353e+1, &
     3.15350626979604161529144e+2, &
    -1.01515636749021914166146e+3, &
    -3.10777167157231109440444e+3, &
     2.25381184209801510330112e+4, &
     4.75584627752788110767815e+3, &
    -1.34659959864969306392456e+5, &
    -1.15132259675553483497211e+5 /)
  real sum
  real x
  real xden
  real xnum
  real y
  real y1
  real ysq
  real z
!
  parity = .false.
  fact = 1.0
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0 ) then

    y = - x
    y1 = aint ( y )
    gamma = y - y1

    if ( gamma /= 0.0 ) then

      if ( y1 /= aint ( y1 * 0.5 ) * 2.0 ) then
        parity = .true.
      end if

      fact = - PI / sin ( PI * gamma )
      y = y + 1.0

    else

      gamma = XINF
      return

    end if

  end if
!
!  Argument .LT. EPS
!
  if ( y < EPS ) then

    if (y >= XMININ ) then
      gamma = 1.0 / y
    else
      gamma = XINF
      return
    end if

  else if ( y < 12.0 ) then

    y1 = y
!
!  0.0 < argument < 1.0
!
    if ( y < 1.0 ) then
      z = y
      y = y + 1.0
!
!  1.0 < argument < 12.0, reduce argument if necessary.
!
    else
      n = int ( y ) - 1
      y = y - real ( n )
      z = y - 1.0
    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0
    xden = 1.0
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    gamma = xnum / xden + 1.0
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then
      gamma = gamma / y1
!
!  Adjust result for case  2.0 < argument < 12.0.
!
    else if ( y1 > y ) then

      do i = 1, n
        gamma = gamma * y
        y = y + 1.0
      end do

    end if
!
!  Evaluate for 12 <= argument.
!
  else

    if ( y <= XBIG ) then

      ysq = y**2
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + SQRTPI
      sum = sum + ( y - 0.5 ) * log ( y )
      gamma = exp ( sum )

    else

      gamma = XINF
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    gamma = - gamma
  end if

  if ( fact /= 1.0 ) then
    gamma = fact / gamma
  end if

  return
end
subroutine gamma_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! GAMMA_CDF evaluates the Gamma CDF.
!

!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B, 
!    0.0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real gamma_inc
  real p2
  real x
  real x2
!
  x2 = ( x - a ) / b
  p2 = c

  cdf = gamma_inc ( x2, p2 )

  return
end
subroutine gamma_check ( a, b, c )
!
!*******************************************************************************
!
!! GAMMA_CHECK checks the parameters of the Gamma PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
  real a
  real b
  real c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GAMMA_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    write ( *, * ) '  B = ', b
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GAMMA_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    write ( *, * ) '  C = ', c
    stop
  end if

  return
end
function gamma_inc ( x, p )
!
!*******************************************************************************
!
!! GAMMA_INC computes the incomplete Gamma function.
!
!
!  Definition:
!
!    GAMMA_INC(X,P) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
!
!  Discussion:
!
!    GAMMA_INC(       0,P) = 0, 
!    GAMMA_INC(Infinity,P) = 1.
!
!  Reference:
!
!    B L Shea,
!    Chi-squared and Incomplete Gamma Integral,
!    Algorithm AS239,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Modified:
!
!    05 January 2000
!
!  Parameters:
!
!    Input, real X, the integral limit parameter.
!    If X is less than or equal to 0, GAMMA_INC is returned as 0.
!
!    Input, real P, the exponent parameter.
!    0.0 < P.
!
!    Output, real GAMMA_INC, the value of the function.
!
  real, parameter :: EXP_ARG_MIN = -88.0
  real, parameter :: OVERFLOW = 1.0E+37
  real, parameter :: PLIMIT = 1000.0
  real, parameter :: TOL = 1.0e-7
  real, parameter :: XBIG = 1.0E+08
!
  real a
  real arg
  real b
  real c
  real cdf
  real gamma_inc
  real gamma_log
  real p
  real pn1
  real pn2
  real pn3
  real pn4
  real pn5
  real pn6
  real rn
  real x
!
  gamma_inc = 0.0

  if ( p <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GAMMA_INC - Fatal error!'
    write ( *, * ) '  Parameter P <= 0.'
    stop
  end if

  if ( x <= 0.0 ) then
    gamma_inc = 0.0
  end if
!
!  Use a normal approximation if PLIMIT < P.
!
  if ( p > PLIMIT ) then
    pn1 = 3.0 * sqrt ( p ) * ( ( x / p ) ** ( 1.0 / 3.0 ) + 1.0 / ( 9.0 * p ) &
      - 1.0 )
    call normal_01_cdf ( pn1, cdf )
    gamma_inc = cdf
    return
  end if
!
!  Is X extremely large compared to P?
!
  if ( x > XBIG ) then
    gamma_inc = 1.0
    return
  end if
!
!  Use Pearson's series expansion.
!  (P is not large enough to force overflow in the log of Gamma.
!
  if ( x <= 1.0 .or. x < p ) then

    arg = p * log ( x ) - x - gamma_log ( p + 1.0 )
    c = 1.0
    gamma_inc = 1.0
    a = p

10      continue

    a = a + 1.0
    c = c * x / a
    gamma_inc = gamma_inc + c

    if ( c > TOL ) then
      go to 10
    end if

    arg = arg + log ( gamma_inc )

    if ( arg >= EXP_ARG_MIN ) then
      gamma_inc = exp ( arg )
    else
      gamma_inc = 0.0
    end if

  else
!
!  Use a continued fraction expansion.
!
    arg = p * log ( x ) - x - gamma_log ( p )
    a = 1.0 - p
    b = a + x + 1.0
    c = 0.0
    pn1 = 1.0
    pn2 = x
    pn3 = x + 1.0
    pn4 = x * b
    gamma_inc = pn3 / pn4

20      continue

    a = a + 1.0
    b = b + 2.0
    c = c + 1.0
    pn5 = b * pn3 - a * c * pn1
    pn6 = b * pn4 - a * c * pn2

    if ( abs ( pn6 ) > 0.0 ) then

      rn = pn5 / pn6

      if ( abs ( gamma_inc - rn ) <= min ( TOL, TOL * rn ) ) then

        arg = arg + log ( gamma_inc )

        if ( arg >= EXP_ARG_MIN ) then
          gamma_inc = 1.0 - exp ( arg )
        else
          gamma_inc = 1.0
        end if

        return

      end if

      gamma_inc = rn

    end if

    pn1 = pn3
    pn2 = pn4
    pn3 = pn5
    pn4 = pn6
!
!  Re-scale terms in continued fraction if terms are large.
!
    if ( abs ( pn5 ) >= OVERFLOW ) then
      pn1 = pn1 / OVERFLOW
      pn2 = pn2 / OVERFLOW
      pn3 = pn3 / OVERFLOW
      pn4 = pn4 / OVERFLOW
    end if

    go to 20

  end if

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
  real, parameter :: d1 = - 5.772156649015328605195174e-1
  real, parameter :: d2 =   4.227843350984671393993777e-1
  real, parameter :: d4 =   1.791759469228055000094023e0
  real, parameter :: EPS = 1.19e-7
  real, parameter :: FRTBIG = 1.42e9
  real, parameter :: PNT68 = 0.6796875
  real, parameter :: SQRTPI = 0.9189385332046727417803297
  real, parameter :: XBIG = 4.08e36
  real, parameter :: XINF = 3.401e38
!
  real, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728e-03, &
     8.4171387781295e-04, &
    -5.952379913043012e-04, &
     7.93650793500350248e-04, &
    -2.777777777777681622553e-03, &
     8.333333333333333331554247e-02, &
     5.7083835261e-03 /)
  real corr
  integer i
  real gamma_log
  real, parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888e0, &
    2.018112620856775083915565e2, &
    2.290838373831346393026739e3, &
    1.131967205903380828685045e4, &
    2.855724635671635335736389e4, &
    3.848496228443793359990269e4, &
    2.637748787624195437963534e4, &
    7.225813979700288197698961e3 /)
  real, parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064e0, &
    5.424138599891070494101986e2, &
    1.550693864978364947665077e4, &
    1.847932904445632425417223e5, &
    1.088204769468828767498470e6, &
    3.338152967987029735917223e6, &
    5.106661678927352456275255e6, &
    3.074109054850539556250927e6 /)
  real, parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062e4, &
    2.426813369486704502836312e6, &
    1.214755574045093227939592e8, &
    2.663432449630976949898078e9, &
    2.940378956634553899906876e10, &
    1.702665737765398868392998e11, &
    4.926125793377430887588120e11, &
    5.606251856223951465078242e11 /)
  real, parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036e1, &
    1.113332393857199323513008e3, &
    7.738757056935398733233834e3, &
    2.763987074403340708898585e4, &
    5.499310206226157329794414e4, &
    6.161122180066002127833352e4, &
    3.635127591501940507276287e4, &
    8.785536302431013170870835e3 /)
  real, parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942e2, &
    7.765049321445005871323047e3, &
    1.331903827966074194402448e5, &
    1.136705821321969608938755e6, &
    5.267964117437946917577538e6, &
    1.346701454311101692290052e7, &
    1.782736530353274213975932e7, &
    9.533095591844353613395747e6 /)
  real, parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843e3, &
    6.393885654300092398984238e5, &
    4.135599930241388052042842e7, &
    1.120872109616147941376570e9, &
    1.488613728678813811542398e10, &
    1.016803586272438228077304e11, &
    3.417476345507377132798597e11, &
    4.463158187419713286462081e11 /)
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
  if ( x <= 0.0 .or. x > XBIG ) then
    gamma_log = XINF
    return
  end if

  if ( x <= EPS ) then

    res = - log ( x )

  else if ( x <= 1.5 ) then

    if ( x < PNT68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0
      xm1 = ( x - 0.5 ) - 0.5
    end if

    if ( x <= 0.5 .or. x >= PNT68 ) then

      xden = 1.0
      xnum = 0.0

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5 ) - 0.5
      xden = 1.0
      xnum = 0.0
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0 ) then

    xm2 = x - 2.0
    xden = 1.0
    xnum = 0.0
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0 ) then

    xm4 = x - 4.0
    xden = - 1.0
    xnum = 0.0
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0

    if ( x <= FRTBIG ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + SQRTPI - 0.5 * corr
    res = res + x * ( corr - 1.0 )

  end if

  gamma_log = res

  return
end
function gamma_log_int ( n )
!
!*******************************************************************************
!
!! GAMMA_LOG_INT computes the logarithm of Gamma of an integer N.
!
!
!  Modified:
!
!    08 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the logarithm of the Gamma function.  
!    0 < N.
!
!    Output, real GAMMA_LOG_INT, the logarithm of the Gamma function of N.
!
  real gamma_log
  real gamma_log_int
  integer n
!
  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GAMMA_LOG_INT - Fatal error!'
    write ( *, * ) '  Illegal input value of N = ', n
    write ( *, * ) '  But N must be strictly positive.'
    stop
  end if

  gamma_log_int = gamma_log ( real ( n ) )

  return
end
subroutine gamma_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! GAMMA_MEAN returns the mean of the Gamma PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real mean
!
  mean = a + b * c

  return
end
subroutine gamma_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! GAMMA_PDF evaluates the Gamma PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = EXP ( ( A - X ) / B ) * ( ( X - A ) / B )**(C-1) 
!      / ( B * GAMMA ( C ) )
!
!  Discussion:
!
!    GAMMA_PDF(A,B,C), where C is an integer, is the Erlang PDF.
!    GAMMA_PDF(A,B,1) is the Exponential PDF.
!    GAMMA_PDF(0,2,C/2) is the Chi Squared PDF with C degrees of freedom.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B, 
!    0.0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real gamma
  real pdf
  real x
  real y
!
  if ( x <= a ) then

    pdf = 0.0

  else

    y = ( x - a ) / b

    pdf = y**( c - 1.0 ) / ( b * gamma ( c ) * exp ( y ) )

  end if

  return
end
subroutine gamma_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! GAMMA_SAMPLE samples the Gamma PDF.
!
!
!  References:
!
!    J H Ahrens and U Dieter,
!    Generating Gamma Variates by a Modified Rejection Technique,
!    Communications of the ACM, 
!    Volume 25, Number 1, January 1982, pages 47 - 54.
!
!    J H Ahrens and U Dieter,
!    Computer Methods for Sampling from Gamma, Beta, Poisson and
!      Binomial Distributions.
!    Computing, Volume 12, 1974, pages 223 - 246.
!
!  Modified:
!
!    02 January 2000
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B, 
!    0.0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real, parameter :: a1 =   0.3333333
  real, parameter :: a2 = - 0.2500030
  real, parameter :: a3 =   0.2000062
  real, parameter :: a4 = - 0.1662921
  real, parameter :: a5 =   0.1423657
  real, parameter :: a6 = - 0.1367177
  real, parameter :: a7 =   0.1233795
  real b
  real bcoef
  real c
  real co
  real d
  real e
  real, parameter :: e1 = 1.0
  real, parameter :: e2 = 0.4999897
  real, parameter :: e3 = 0.1668290
  real, parameter :: e4 = 0.0407753
  real, parameter :: e5 = 0.0102930
  integer iseed
  real p
  real q
  real q0
  real, parameter :: q1 =   0.04166669
  real, parameter :: q2 =   0.02083148
  real, parameter :: q3 =   0.00801191
  real, parameter :: q4 =   0.00144121
  real, parameter :: q5 = - 0.00007388
  real, parameter :: q6 =   0.00024511
  real, parameter :: q7 =   0.00024240
  real r
  real s
  real si
  real s2
  real t
  real u
  real uniform_01_sample
  real v
  real w
  real x
!
!  C < 1.
!
  if ( c < 1.0 ) then

10      continue

    p = ( 1.0 + 0.3678794 * c ) * uniform_01_sample ( iseed )

    call exponential_01_sample ( iseed, e )

    if ( p >= 1.0 ) then

      x = - log ( ( 1.0 + 0.3678794 * a - p ) / c )

      if ( e >= ( 1.0 - c ) * log ( x ) ) then
        x = a + b * x
        return
      end if

    else

      x = exp ( log ( p ) / c )

      if ( e >= x ) then
        x = a + b * x
        return
      end if

    end if

    go to 10
!
!  1 <= C.
!
  else

    s2 = c - 0.5
    s = sqrt ( c - 0.5 )
    d = sqrt ( 32.0 ) - 12.0 * sqrt ( c - 0.5 )

    call normal_01_sample ( iseed, t )
    x = ( sqrt ( c - 0.5 ) + 0.5 * t )**2

    if ( t >= 0.0 ) then
      x = a + b * x
      return
    end if

    u = uniform_01_sample ( iseed )

    if ( d * u <= t**3 ) then
      x = a + b * x
      return
    end if

    r = 1.0 / c

    q0 = ( ( ( ( ( ( &
           q7   * r &
         + q6 ) * r &
         + q5 ) * r &
         + q4 ) * r &
         + q3 ) * r &
         + q2 ) * r &
         + q1 ) * r

    if ( c <= 3.686 ) then
      bcoef = 0.463 + s - 0.178 * s2
      si = 1.235
      co = 0.195 / s - 0.079 + 0.016 * s
    else if ( c <= 13.022 ) then
      bcoef = 1.654 + 0.0076 * s2
      si = 1.68 / s + 0.275
      co = 0.062 / s + 0.024
    else
      bcoef = 1.77
      si = 0.75
      co = 0.1515 / s
    end if

    if ( sqrt ( c - 0.5 ) + 0.5 * t > 0.0 ) then

      v = 0.5 * t / s

      if ( abs ( v ) > 0.25 ) then
        q = q0 - s * t + 0.25 * t**2 + 2.0 * s2 * log ( 1.0 + v )
      else
        q = q0 + 0.5 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
      end if

      if ( log ( 1.0 - u ) <= q ) then
        x = a + b * x
        return
      end if

    end if

20      continue

    call exponential_01_sample ( iseed, e )

    u = 2.0 * uniform_01_sample ( iseed ) - 1.0
    t = bcoef + sign ( si * e, u )

    if ( t >= - 0.7187449 ) then

      v = 0.5 * t / s

      if ( abs ( v ) > 0.25 ) then
        q = q0 - s * t + 0.25 * t**2 + 2.0 * s2 * log ( 1.0 + v )
      else
        q = q0 + 0.5 * t**2 * ( ( ( ( ( ( &
             a7   * v &
           + a6 ) * v &
           + a5 ) * v &
           + a4 ) * v &
           + a3 ) * v &
           + a2 ) * v &
           + a1 ) * v
      end if

      if ( q > 0.0 ) then

        if ( q > 0.5 ) then
          w = exp ( q ) - 1.0
        else
          w = ( ( ( ( &
                  e5   * q &
                + e4 ) * q &
                + e3 ) * q &
                + e2 ) * q &
                + e1 ) * q
        end if

        if ( co * abs ( u ) <= w * exp ( e - 0.5 * t**2 ) ) then
          x = a + b * ( s + 0.5 * t )**2
          return
        end if

      end if

    end if

    go to 20

  end if

end
subroutine gamma_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! GAMMA_VARIANCE returns the variance of the Gamma PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real c
  real variance
!
  variance = b**2 * c

  return
end
subroutine genlogistic_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! GENLOGISTIC_CDF evaluates the Generalized Logistic CDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real x
  real y
!
  y = ( x - a ) / b

  cdf = 1.0 / ( 1.0 + exp ( - y ) )**c

  return
end
subroutine genlogistic_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! GENLOGISTIC_CDF_INV inverts the Generalized Logistic CDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real c
  real cdf
  real x
!
  if ( cdf <= 0.0 ) then
    x = - huge ( x )
  else if ( cdf < 1.0 ) then
    x = a - b * log ( cdf**( - 1.0 / c ) - 1.0 )
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
  end if

  return
end
subroutine genlogistic_check ( a, b, c )
!
!*******************************************************************************
!
!! GENLOGISTIC_CHECK checks the parameters of the Generalized Logistic CDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
  real a
  real b
  real c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GENLOGISTIC_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GENLOGISTIC_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    stop
  end if

  return
end
subroutine genlogistic_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! GENLOGISTIC_MEAN returns the mean of the Generalized Logistic PDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real digamma
  real euler_constant
  real mean
!
  mean = a + b * ( euler_constant ( ) + digamma ( c ) )

  return
end
subroutine genlogistic_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! GENLOGISTIC_PDF evaluates the Generalized Logistic PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = ( C / B ) * EXP ( ( A - X ) / B ) /
!      ( ( 1 + EXP ( ( A - X ) / B ) )**(C+1) )
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real pdf
  real x
  real y
!
  y = ( x - a ) / b

  pdf = ( c / b ) * exp ( - y ) / ( 1.0 + exp ( - y ) )**( c + 1.0 )

  return
end
subroutine genlogistic_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! GENLOGISTIC_SAMPLE samples the Generalized Logistic PDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real c
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call genlogistic_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine genlogistic_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! GENLOGISTIC_VARIANCE returns the variance of the Generalized Logistic PDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real c
  real trigamma
  real variance
!
  variance = b**2 * ( PI**2.0 / 6.0 + trigamma ( c ) )

  return
end
subroutine geometric_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! GEOMETRIC_CDF evaluates the Geometric CDF.
!
!
!  Definition:
!
!    CDF(X,P) is the probability that there will be at least one
!    successful trial in the first X Bernoulli trials, given that
!    the probability of success in a single trial is P.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the maximum number of trials.
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real cdf
  integer x
!
  if ( x <= 0 ) then
    cdf = 0.0
  else if ( a == 0.0 ) then
    cdf = 0.0
  else if ( a == 1.0 ) then
    cdf = 1.0
  else
    cdf = 1.0 - ( 1.0 - a )**x
  end if

  return
end
subroutine geometric_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! GEOMETRIC_CDF_INV inverts the Geometric CDF.
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
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, integer X, the corresponding value of X.
!
  integer, parameter :: INTMAX = 2147483647
!
  real a
  real cdf
  integer x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GEOMETRIC_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( a == 1.0 ) then
    x = 1
  else if ( a == 0.0 ) then
    x = INTMAX
  else
    x = 1 + int ( log ( 1.0 - cdf ) / log ( 1.0 - a ) )
  end if

  return
end
subroutine geometric_check ( a )
!
!*******************************************************************************
!
!! GEOMETRIC_CHECK checks the parameter of the Geometric CDF.
!
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
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
  real a
!
  if ( a < 0.0 .or. a > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GEOMETRIC_CHECK - Fatal error!'
    write ( *, * ) '  A < 0 or 1 < A.'
    stop
  end if

  return
end
subroutine geometric_mean ( a, mean )
!
!*******************************************************************************
!
!! GEOMETRIC_MEAN returns the mean of the Geometric PDF.
!
!
!  Discussion:
!
!    MEAN is the expected value of the number of trials required
!    to obtain a single success.
! 
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real mean
!
  mean = 1.0 / a

  return
end
subroutine geometric_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! GEOMETRIC_PDF evaluates the Geometric PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = A * ( 1 - A )**(X-1)
!
!  Definition:
!
!    PDF(X)(A) is the probability that exactly X Bernoulli trials, each
!    with probability of success A, will be required to achieve 
!    a single success.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the number of trials.
!    0 < X
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real pdf
  integer x
!
!  Special cases.
!
  if ( x < 1 ) then

    pdf = 0.0

  else if ( a == 0.0 ) then

    pdf = 0.0

  else if ( a == 1.0 ) then

    if ( x == 1 ) then
      pdf = 1.0
    else
      pdf = 0.0
    end if

  else

    pdf = a * ( 1.0 - a )**( x - 1 )

  end if

  return
end
subroutine geometric_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! GEOMETRIC_SAMPLE samples the Geometric PDF.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  real a
  real cdf
  integer iseed
  real uniform_01_sample
  integer x
!
  cdf = uniform_01_sample ( iseed )

  call geometric_cdf_inv ( cdf, a, x )

  return
end
subroutine geometric_variance ( a, variance )
!
!*******************************************************************************
!
!! GEOMETRIC_VARIANCE returns the variance of the Geometric PDF.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the probability of success on one trial.
!    0.0 <= A <= 1.0.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real variance
!
  variance = ( 1.0 - a ) / a**2

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
!    If a FORTRAN77 compiler is used, then a different set of routines
!    must be invoked in order to get the date and time.
!
!  Modified:
!
!    23 March 2000
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
!  FORTRAN90 declarations
!
  character ( len = 10 ) time90
  character ( len = 8 ) today90
  character ( len = 5 ) zone
!
!  FORTRAN77 declarations
!
!     integer day
!     integer hour
!     integer minute
!     integer month
!     integer second
!     character*8 time77
!     character*9 today77
!     integer year
!
!  FORTRAN90
!
  call date_and_time ( today90, time90, zone, values )
!
!  FORTRAN77
!
!     call date ( today77 )
!     call s_to_ymd ( today77, 'DD-NNN-YY', year, month, day )
!     call time ( time77 )
!     call s_to_hms ( time77, 'hh:mm:ss', hour, minute, second )
!
!     if ( year == 99 ) then
!       year = 1900 + year
!     else
!       year = 2000 + year
!     end if
!
!     values(1) = year
!     values(2) = month
!     values(3) = day
!     values(4) = 0
!     values(5) = hour
!     values(6) = minute
!     values(7) = second
!     values(8) = 0
!
  temp = 0.0

  temp = temp + dble ( values(2) - 1 ) / 11.0
  temp = temp + dble ( values(3) - 1 ) / 30.0
  temp = temp + dble ( values(5) ) / 23.0
  temp = temp + dble ( values(6) ) / 59.0
  temp = temp + dble ( values(7) ) / 59.0
  temp = temp + dble ( values(8) ) / 999.0
  temp = temp / 6.0

  if ( temp <= 0.0 ) then
    temp = 1.0 / 3.0
  else if ( temp >= 1.0 ) then
    temp = 2.0 / 3.0
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
subroutine gompertz_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! GOMPERTZ_CDF evaluates the Gompertz CDF.
!
!
!  Reference:
!
!    Johnson, Kotz, and Balakrishnan,
!    Continuous Univariate Distributions, Volume 2, second edition,
!    Wiley, 1994, pages 25-26.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A > 1, B > 0.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( x <= 0.0 ) then
    cdf = 0.0
  else
    cdf = 1.0 - exp ( - b * ( a**x - 1.0 ) / log ( a ) )
  end if

  return
end
subroutine gompertz_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! GOMPERTZ_CDF_INV inverts the Gompertz CDF.
!
!
!  Reference:
!
!    Johnson, Kotz, and Balakrishnan,
!    Continuous Univariate Distributions, Volume 2, second edition,
!    Wiley, 1994, pages 25-26.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A > 1, B > 0.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 ) then
    x = 0.0
  else if ( cdf < 1.0 ) then
    x = log ( 1.0 - log ( 1.0 - cdf ) * log ( a ) / b  ) / log ( a )
  else
    x = huge ( x )
  end if

  return
end
subroutine gompertz_check ( a, b )
!
!*******************************************************************************
!
!! GOMPERTZ_CHECK checks the parameters of the Gompertz PDF.
!
!
!  Reference:
!
!    Johnson, Kotz, and Balakrishnan,
!    Continuous Univariate Distributions, Volume 2, second edition,
!    Wiley, 1994, pages 25-26.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A > 1, B > 0.
!
  real a
  real b
!
  if ( a <= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GOMPERTZ_CHECK - Fatal error!'
    write ( *, * ) '  A <= 1.0!'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GOMPERTZ_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.0!'
    stop
  end if

  return
end
subroutine gompertz_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! GOMPERTZ_PDF evaluates the Gompertz PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = B * A**X / exp ( B * ( A**X - 1 ) / log ( A ) )     
!
!    for
!
!      X >= 0.0, 
!      A >  1.0,
!      B >  0.0.
!
!  Reference:
!
!    Johnson, Kotz, and Balakrishnan,
!    Continuous Univariate Distributions, Volume 2, second edition,
!    Wiley, 1994, pages 25-26.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A > 1, B > 0.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x

  if ( x < 0.0 ) then

    pdf = 0.0

  else if ( a > 1.0 ) then

    pdf = exp ( log ( b ) + x * log ( a ) - ( b / log ( a ) ) * ( a**x - 1.0 ) )

  end if

  return
end
subroutine gompertz_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! GOMPERTZ_SAMPLE samples the Gompertz PDF.
!
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A > 1, B > 0.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call gompertz_cdf_inv ( cdf, a, b, x )

  return
end
subroutine gumbel_cdf ( x, cdf )
!
!*******************************************************************************
!
!! GUMBEL_CDF evaluates the Gumbel CDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Output, real CDF, the value of the CDF.
!
  real cdf
  real x
!
  cdf = exp ( - exp ( - x ) )

  return
end
subroutine gumbel_cdf_inv ( cdf, x )
!
!*******************************************************************************
!
!! GUMBEL_CDF_INV inverts the Gumbel CDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, real X, the corresponding argument of the CDF.
!
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GUMBEL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x =  - log ( - log ( cdf ) )

  return
end
subroutine gumbel_mean ( mean )
!
!*******************************************************************************
!
!! GUMBEL_MEAN returns the mean of the Gumbel PDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MEAN, the mean of the PDF.
!
  real euler_constant
  real mean
!
  mean = euler_constant ( )

  return
end
subroutine gumbel_pdf ( x, pdf )
!
!*******************************************************************************
!
!! GUMBEL_PDF evaluates the Gumbel PDF.
!
!
!  Formula:
!
!    PDF(X) = EXP ( - X - EXP ( - X  ) ).
!
!  Discussion:
!
!    GUMBEL_PDF(X) = EXTREME_PDF(X)(0,1)
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998.
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Output, real PDF, the value of the PDF.
!
  real pdf
  real x
!
  pdf = exp ( - x - exp ( - x ) )

  return
end
subroutine gumbel_sample ( iseed, x )
!
!*******************************************************************************
!
!! GUMBEL_SAMPLE samples the Gumbel PDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call gumbel_cdf_inv ( cdf, x )

  return
end
subroutine gumbel_variance ( variance )
!
!*******************************************************************************
!
!! GUMBEL_VARIANCE returns the variance of the Gumbel PDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real variance
!
  variance = PI**2 / 6.0

  return
end
subroutine half_normal_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! HALF_NORMAL_CDF evaluates the Half Normal CDF.
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
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real cdf2
  real x
!
  if ( x <= a ) then
    cdf = 0.0
  else
    call normal_cdf ( x, a, b, cdf2 ) 
    cdf = 2.0 * cdf2 - 1.0
  end if

  return
end
subroutine half_normal_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! HALF_NORMAL_CDF_INV inverts the Half Normal CDF.
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
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real cdf2
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'HALF_NORMAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.5 * ( cdf + 1.0 )

  call normal_cdf_inv ( cdf2, a, b, x ) 

  return
end
subroutine half_normal_check ( a, b )
!
!*******************************************************************************
!
!! HALF_NORMAL_CHECK checks the parameters of the Half Normal PDF.
!
! 
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'HALF_NORMAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine half_normal_mean ( a, b, mean )
!
!*******************************************************************************
!
!! HALF_NORMAL_MEAN returns the mean of the Half Normal PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real mean
!
  mean = a + b * sqrt ( 2.0 / PI )

  return
end
subroutine half_normal_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! HALF_NORMAL_PDF evaluates the Half Normal PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = 
!      SQRT ( 2 / PI ) * ( 1 / B ) * EXP ( - 0.5 * ( ( X - A ) / B )**2 )
!
!    for A <= X
!
!  Discussion:
!
!    The Half Normal PDF is a special case of both the Chi PDF and the
!    Folded Normal PDF.
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
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
  real y
!
  if ( x <= a ) then

    pdf = 0.0

  else

    y = ( x - a ) / b

    pdf = sqrt ( 2.0 / PI ) * ( 1.0 / b ) * exp ( - 0.5 * y**2 )

  end if

  return
end
subroutine half_normal_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! HALF_NORMAL_SAMPLE samples the Half Normal PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )
  
  call half_normal_cdf_inv ( cdf, a, b, x )
  
  return
end
subroutine half_normal_variance ( a, b, variance )
!
!*******************************************************************************
!
!! HALF_NORMAL_VARIANCE returns the variance of the Half Normal PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real variance
!
  variance = b**2 * ( 1.0 - 2.0 / PI )

  return
end
subroutine hypergeometric_cdf ( x, n, m, l, cdf )
!
!*******************************************************************************
!
!! HYPERGEOMETRIC_CDF evaluates the Hypergeometric CDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the CDF.
!
!    Input, integer N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer L, the number of balls to select from.
!    0 <= L.
!
!    Output, real CDF, the value of the CDF.
!
  real cdf
  real c1_log
  real c2_log
  integer l
  integer m
  integer n
  real pdf
  integer x
  integer x2
!
  call binomial_coef_log ( l - m, n, c1_log )
  call binomial_coef_log ( l, n, c2_log )

  pdf = exp ( c1_log - c2_log )
  cdf = pdf

  do x2 = 0, x - 1

    pdf = pdf * real ( ( m - x2 ) * ( n - x2 ) ) &
      / real ( ( x2 + 1 ) * ( l - m - n + x2 + 1 ) )

    cdf = cdf + pdf

  end do

  return
end
subroutine hypergeometric_check ( n, m, l )
!
!*******************************************************************************
!
!! HYPERGEOMETRIC_CHECK checks the parameters of the Hypergeometric CDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer L, the number of balls to select from.
!    0 <= L.
!
  integer l
  integer m
  integer n
!
  if ( n < 0 .or. n > l ) then
    write ( *, * ) ' '
    write ( *, * ) 'HYPERGEOMETRIC_CHECK - Fatal error!'
    write ( *, * ) '  Input N is out of range.'
    stop
  end if

  if ( m < 0 .or. m > l ) then
    write ( *, * ) ' '
    write ( *, * ) 'HYPERGEOMETRIC_CHECK - Fatal error!'
    write ( *, * ) '  Input M is out of range.'
    stop
  end if

  if ( l < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'HYPERGEOMETRIC_CHECK - Fatal error!'
    write ( *, * ) '  Input L is out of range.'
    stop
  end if

  return
end
subroutine hypergeometric_mean ( n, m, l, mean )
!
!*******************************************************************************
!
!! HYPERGEOMETRIC_MEAN returns the mean of the Hypergeometric PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer L, the number of balls to select from.
!    0 <= L.
!
!    Output, real MEAN, the mean of the PDF.
!
  integer l
  integer m
  real mean
  integer n
!
  mean = real ( n * m ) / real ( l )

  return
end
subroutine hypergeometric_pdf ( x, n, m, l, pdf )
!
!*******************************************************************************
!
!! HYPERGEOMETRIC_PDF evaluates the Hypergeometric PDF.
!
!
!  Formula:
!
!    PDF(X)(N,M,L) = C(M,X) * C(L-M,N-X) / C(L,N).
!
!  Definition:
!
!    PDF(X)(N,M,L) is the probability of drawing X white balls in a 
!    single random sample of size N from a population containing 
!    M white balls and a total of L balls.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the desired number of white balls.
!    0 <= X <= N, usually, although any value of X can be given.
!
!    Input, integer N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer L, the number of balls to select from.
!    0 <= L.
!
!    Output, real PDF, the probability of exactly K white balls.
!
  real c1
  real c2
  real c3
  integer l
  integer m
  integer n
  real pdf
  real pdf_log
  integer x
!
!  Special cases.
!
  if ( x < 0 ) then

    pdf = 1.0

  else if ( x > n ) then

    pdf = 0.0

  else if ( x > m ) then

    pdf = 0.0

  else if ( x > l ) then

    pdf = 0.0

  else if ( n == 0 ) then

    if ( x == 0 ) then
      pdf = 1.0
    else
      pdf = 0.0
    end if

  else

    call binomial_coef_log ( m, x, c1 ) 
    call binomial_coef_log ( l-m, n-x, c2 )
    call binomial_coef_log ( l, n, c3 )

    pdf_log = c1 + c2 - c3

    pdf = exp ( pdf_log )

  end if

  return
end
subroutine hypergeometric_sample ( n, m, l, iseed, x )
!
!*******************************************************************************
!
!! HYPERGEOMETRIC_SAMPLE samples the Hypergeometric PDF.
!
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 165.
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
!    Input, integer N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer L, the number of balls to select from.
!    0 <= L.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  real a
  real b
  real c1_log
  real c2_log
  integer iseed
  integer l
  integer m
  integer n
  real u
  real uniform_01_sample
  integer x
!
  call binomial_coef_log ( l - m, n, c1_log )
  call binomial_coef_log ( l, n, c2_log )

  a = exp ( c1_log - c2_log )
  b = a

  u = uniform_01_sample ( iseed )

  x = 0

  do while ( u > a )

    b = b * real ( ( m - x ) * ( n - x ) ) &
      / real ( ( x + 1 ) * ( l - m - n + x + 1 ) )

    a = a + b

    x = x + 1

  end do

  return
end
subroutine hypergeometric_variance ( n, m, l, variance )
!
!*******************************************************************************
!
!! HYPERGEOMETRIC_VARIANCE returns the variance of the Hypergeometric PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of balls selected.
!    0 <= N <= L.
!
!    Input, integer M, the number of white balls in the population.
!    0 <= M <= L.
!
!    Input, integer L, the number of balls to select from.
!    0 <= L.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer l
  integer m
  integer n
  real variance
!
  variance = real ( n * m * ( l - m ) * ( l - n ) )  / real ( l**2 * ( l - 1 ) )

  return
end
function i_factorial ( n )
!
!*******************************************************************************
!
!! I_FACTORIAL returns N!.
!
!
!  Definition:
!
!    N! = Product ( 1 <= I <= N ) I
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
!    Input, integer N, the argument of the function.
!    0 <= N.
!
!    Output, integer I_FACTORIAL, the factorial of N.
!
  integer i_factorial
  integer i
  integer n
!
  if ( n < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_FACTORIAL - Fatal error!'
    write ( *, * ) '  N < 0.'
    stop
  end if

  i_factorial = 1

  do i = 2, n
    i_factorial = i_factorial * i
  end do

  return
end
subroutine i_random ( ilo, ihi, iseed, i )
!
!*******************************************************************************
!
!! I_RANDOM returns a random integer in a given range.
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
!    Input, integer ILO, IHI, the minimum and maximum acceptable values.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer I, the randomly chosen integer.
!
  integer i
  integer ihi
  integer ilo
  integer iseed
  real r
  real rhi
  real rlo
  real uniform_01_sample
!
!  Pick a random number in (0,1).
!
  r = uniform_01_sample ( iseed )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo ) - 0.5
  rhi = real ( ihi ) + 0.5
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0 - r ) * rlo + r * rhi )
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
subroutine i0 ( arg, result )
!
!*******************************************************************************
!
!! I0 evaluates the modified Bessel function of the first kind and order zero.  
!
!
!  Discussion:
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.  This transportable program is patterned after
!    the machine-dependent FUNPACK packet NATSI0, but cannot match
!    that version for efficiency or accuracy.  This version uses
!    rational functions that theoretically approximate I-SUB-0(X)
!    to at least 18 significant decimal digits.  
!
!  Machine-dependent constants:
!
!    beta   = Radix for the floating-point system
!    maxexp = Smallest power of beta that overflows
!    XSMALL = Positive argument such that 1.0 - X = 1.0 to
!             machine precision for all ABS(X) .LE. XSMALL.
!    XINF =   Largest positive machine number; approximately
!             beta**maxexp
!    XMAX =   Largest argument acceptable to BESI0;  Solution to
!             equation:
!               W(X) * (1+1/(8*X)+9/(128*X**2) = beta**maxexp
!             where  W(X) = EXP(X)/SQRT(2*PI*X)
!
!    Approximate values for some important machines are:
!
!                             beta       maxexp       XSMALL
!
!    CRAY-1        (S.P.)       2         8191       3.55E-15
!    Cyber 180/855
!      under NOS   (S.P.)       2         1070       3.55E-15
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)       2          128       2.98E-8
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)       2         1024       5.55D-17
!    IBM 3033      (D.P.)      16           63       6.95D-18
!    VAX           (S.P.)       2          127       2.98E-8
!    VAX D-Format  (D.P.)       2          127       6.95D-18
!    VAX G-Format  (D.P.)       2         1023       5.55D-17
!
!
!                                  XINF          XMAX
!
!    CRAY-1        (S.P.)       5.45E+2465     5682.810
!    Cyber 180/855
!      under NOS   (S.P.)       1.26E+322       745.893
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)       3.40E+38         91.900
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)       1.79D+308       713.986
!    IBM 3033      (D.P.)       7.23D+75        178.182
!    VAX           (S.P.)       1.70D+38         91.203
!    VAX D-Format  (D.P.)       1.70D+38         91.203
!    VAX G-Format  (D.P.)       8.98D+307       713.293
!
!  Authors: 
!
!    W. J. Cody and L. Stoltz,
!    Mathematics and Computer Science Division,
!    Argonne National Laboratory,
!    Argonne, IL 60439.
!
!  Parameters:
!
!    Input, real ARG, the argument.
!
!    Output, real RESULT, the value of the modified Bessel function
!    of the first kind.
!
  real, parameter :: EXP40 = 2.353852668370199854E17
  real, parameter :: REC15 = 6.6666666666666666666E-2
  real, parameter :: XSMALL = 2.98E-8
  real, parameter :: XINF = 3.40E38
  real, parameter :: XMAX = 91.9E0
!
  real a
  real arg
  real b
  integer i
  real, parameter, dimension ( 15 ) :: p = (/ &
    -5.2487866627945699800E-18, &
    -1.5982226675653184646E-14, &
    -2.6843448573468483278E-11, &
    -3.0517226450451067446E-08, &
    -2.5172644670688975051E-05, &
    -1.5453977791786851041E-02, &
    -7.0935347449210549190E+00, &
    -2.4125195876041896775E+03, &
    -5.9545626019847898221E+05, &
    -1.0313066708737980747E+08, &
    -1.1912746104985237192E+10, &
    -8.4925101247114157499E+11, &
    -3.2940087627407749166E+13, &
    -5.5050369673018427753E+14, &
    -2.2335582639474375249E+15 /)
  real, parameter, dimension ( 8 ) :: pp = (/ &
    -3.9843750000000000000E-01, &
     2.9205384596336793945E+00, &
    -2.4708469169133954315E+00, &
     4.7914889422856814203E-01, &
    -3.7384991926068969150E-03, &
    -2.6801520353328635310E-03, &
     9.9168777670983678974E-05, &
    -2.1877128189032726730E-06 /)
  real, parameter, dimension ( 5 ) :: q = (/ &
    -3.7277560179962773046E+03, &
     6.5158506418655165707E+06, &
    -6.5626560740833869295E+09, &
     3.7604188704092954661E+12, &
    -9.7087946179594019126E+14 /)
  real, parameter, dimension ( 7 ) :: qq = (/ &
    -3.1446690275135491500E+01, &
     8.5539563258012929600E+01, &
    -6.0228002066743340583E+01, &
     1.3982595353892851542E+01, &
    -1.1151759188741312645E+00, &
     3.2547697594819615062E-02, &
    -5.5194330231005480228E-04 /)
  real result
  real sump
  real sumq
  real x
  real xx
!
  x = abs ( arg )

  if ( x < XSMALL ) then
    result = 1.0
  else if ( x < 15.0 ) then
!
!  XSMALL <= ABS(ARG) < 15.0
!
    xx = x**2
    sump = p(1)
    do i = 2, 15
      sump = sump * xx + p(i)
    end do

    xx = xx - 225.0
    sumq = (((( &
           xx + q(1) ) &
         * xx + q(2) ) &
         * xx + q(3) ) &
         * xx + q(4) ) &
         * xx + q(5)

    result = sump / sumq

  else if ( x >= 15.0 ) then

    if ( x > XMAX ) then
      result = XINF
    else
!
!  15.0 <= ABS(ARG)
!
      xx = 1.0 / x - REC15

      sump = ((((((  &
                  pp(1) &
           * xx + pp(2) ) &
           * xx + pp(3) ) &
           * xx + pp(4) ) &
           * xx + pp(5) ) &
           * xx + pp(6) ) &
           * xx + pp(7) ) &
           * xx + pp(8)

      sumq = (((((( &
             xx + qq(1) ) &
           * xx + qq(2) ) &
           * xx + qq(3) ) &
           * xx + qq(4) ) &
           * xx + qq(5) ) &
           * xx + qq(6) ) &
           * xx + qq(7)

      result = sump / sumq
!
!  Calculation reformulated to avoid premature overflow.
!
      if ( x <= XMAX - 15.0 ) then
        a = exp ( x )
        b = 1.0
      else
        a = exp ( x - 40.0 )
        b = EXP40
      end if

      result = ( ( result * a - pp(1) * a ) / sqrt ( x ) ) * b
    
    end if

  end if

  return
end
subroutine inverse_gaussian_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! INVERSE_GAUSSIAN_CDF evaluates the Inverse Gaussian CDF.
!
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
!    Input, real X, the argument of the CDF.
!    0.0 < X.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real cdf1
  real cdf2
  real x
  real x1
  real x2
!
  if ( x <= 0.0 ) then

    cdf = 0.0

  else

    x1 = sqrt ( b / x ) * ( x - a ) / a
    call normal_01_cdf ( x1, cdf1 )

    x2 = - sqrt ( b / x ) * ( x + a ) / a
    call normal_01_cdf ( x2, cdf2 )

    cdf = cdf1 + exp ( 2.0 * b / a ) * cdf2

  end if

  return
end
subroutine inverse_gaussian_check ( a, b )
!
!*******************************************************************************
!
!! INVERSE_GAUSSIAN_CHECK checks the parameters of the Inverse Gaussian CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
  real a
  real b
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'INVERSE_GAUSSIAN_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'INVERSE_GAUSSIAN_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine inverse_gaussian_mean ( a, b, mean )
!
!*******************************************************************************
!
!! INVERSE_GAUSSIAN_MEAN returns the mean of the Inverse Gaussian PDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine inverse_gaussian_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! INVERSE_GAUSSIAN_PDF evaluates the Inverse Gaussian PDF.
!
!
!  Discussion:
!
!    The Inverse Gaussian PDF is also known as the Wald PDF
!    and the Inverse Normal PDF.
!
!  Formula:
!
!    PDF(X)(A,B) 
!      = SQRT ( B / ( 2 * PI * X**3 ) ) 
!        * EXP ( - B * ( X - A )**2 / ( 2.0 * A**2 * X ) )
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
!    Input, real X, the argument of the PDF.
!    0.0 < X
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
!
  if ( x <= 0.0 ) then
    pdf = 0.0
  else
    pdf = sqrt ( b / ( 2.0 * PI * x**3 ) ) * &
      exp ( - b * ( x - a )**2 / ( 2.0 * a**2 * x ) )
  end if

  return
end
subroutine inverse_gaussian_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! INVERSE_GAUSSIAN_SAMPLE samples the Inverse Gaussian PDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input/output, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  integer iseed
  real phi
  real t
  real u
  real uniform_01_sample
  real x
  real y
  real z
!
  phi = b / a
  call normal_01_sample ( iseed, z )
  y = z**2

  t = 1.0 + 0.5 * ( y - sqrt ( 4.0 * phi * y + y**2 ) ) / phi
  u = uniform_01_sample ( iseed )

  if ( u * ( 1.0 + t ) <= 1.0 ) then
    x = a * t
  else
    x = a / t
  end if

  return
end
subroutine inverse_gaussian_variance ( a, b, variance )
!
!*******************************************************************************
!
!! INVERSE_GAUSSIAN_VARIANCE returns the variance of the Inverse Gaussian PDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = a**3 / b

  return
end
subroutine irow_max ( lda, m, n, x, ixmax, xmax )
!
!*******************************************************************************
!
!! IROW_MAX returns the maximums of rows of an integer array.
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
!    Input, integer LDA, the leading dimension of X, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, integer X(LDA,N), the array to be examined.
!
!    Output, integer IXMAX(M); IXMAX(I) is the column of X in which
!    the maximum for row I occurs.
!
!    Output, integer XMAX(M), the maximums of the rows of X.
!
  integer lda
  integer m
  integer n
!
  integer i
  integer ixmax(m)
  integer j
  integer x(lda,n)
  integer xmax(m)
!
  do i = 1, m

    ixmax(i) = 1
    xmax(i) = x(i,1)
    do j = 2, n
      if ( x(i,j) > xmax(i) ) then
        ixmax(i) = j
        xmax(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine irow_mean ( maxrow, nrow, ncol, table, mean )
!
!*******************************************************************************
!
!! IROW_MEAN returns the means of the rows of a table.
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
!    Input, integer MAXROW, the first dimension of TABLE.
!
!    Input, integer NROW, NCOL, the number of rows and columns of data in TABLE.
!
!    Input, integer TABLE(MAXROW,NCOL), an array of NROW by NCOL data.
!
!    Output, real MEAN(NROW), the mean or average value of each row of TABLE.
!
  integer maxrow
  integer ncol
  integer nrow
!
  integer i
  integer j
  real mean(nrow)
  integer table(maxrow,ncol)
!
  do i = 1, nrow
    mean(i) = 0.0
    do j = 1, ncol
      mean(i) = mean(i) + table(i,j)
    end do
    mean(i) = mean(i) / real ( ncol )
  end do

  return
end
subroutine irow_min ( lda, m, n, x, ixmin, xmin )
!
!*******************************************************************************
!
!! IROW_MIN returns the minimums of rows of an integer array.
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
!    Input, integer LDA, the leading dimension of X, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, integer X(LDA,N), the array to be examined.
!
!    Output, integer IXMIN(M); IXMIN(I) is the column of X in which
!    the minimum for row I occurs.
!
!    Output, integer XMIN(M), the minimums of the rows of X.
!
  integer lda
  integer m
  integer n
!
  integer i
  integer ixmin(m)
  integer j
  integer x(lda,n)
  integer xmin(m)
!
  do i = 1, m

    ixmin(i) = 1
    xmin(i) = x(i,1)
    do j = 2, n
      if ( x(i,j) < xmin(i) ) then
        ixmin(i) = j
        xmin(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine irow_variance ( maxrow, nrow, ncol, table, variance )
!
!*******************************************************************************
!
!! IROW_VARIANCE returns the variances of the rows of a table.
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
!    Input, integer MAXROW, the first dimension of TABLE.
!
!    Input, integer NROW, NCOL, the number of rows and columns of data in TABLE.
!
!    Input, integer TABLE(MAXROW,NCOL), an array of NROW by NCOL data.
!
!    Output, real VARIANCE(NROW), the variance of each row of TABLE.
!
  integer maxrow
  integer ncol
  integer nrow
!
  integer i
  integer j
  real mean
  integer table(maxrow,ncol)
  real variance(nrow)
!
  call irow_mean ( maxrow, nrow, ncol, table, variance )

  do i = 1, nrow

    mean = variance(i)

    variance(i) = 0.0
    do j = 1, ncol
      variance(i) = variance(i) + ( real ( table(i,j) ) - mean )**2
    end do

    if ( ncol > 1 ) then
      variance(i) = variance(i) / real ( ncol - 1 )
    else
      variance(i) = 0.0
    end if

  end do

  return
end
subroutine ivec_max ( n, iarray, index, imax )
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
!    Input, integer IARRAY(N), the array.
!
!    Output, integer INDEX, the index of the largest entry.
!
!    Output, integer IMAX, the value of the largest entry.
!
  integer n
!
  integer i
  integer iarray(n)
  integer index
  integer imax
!
  if ( n <= 0 ) then

    imax = 0
    index = 0

  else

    imax = iarray(1)
    index = 1
    do i = 2, n

      if ( iarray(i) > imax ) then
        imax = iarray(i)
        index = i
      end if

    end do

  end if

  return
end
subroutine ivec_mean ( n, x, mean )
!
!*******************************************************************************
!
!! IVEC_MEAN returns the mean of an integer vector.
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
!    Input, integer X(N), the vector whose mean is desired.
!
!    Output, real MEAN, the mean, or average, of the vector entries.
!
  integer n
!
  integer i
  real mean
  integer x(n)
!
  mean = 0.0
  do i = 1, n
    mean = mean + real ( x(i) )
  end do

  mean = mean / real ( n )

  return
end
subroutine ivec_min ( n, iarray, index, imin )
!
!*******************************************************************************
!
!! IVEC_MIN computes the minimum element of an integer array.
!
!
!  Modified:
!
!    09 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer IARRAY(N), the array.
!
!    Output, integer INDEX, the index of the smallest entry.
!
!    Output, integer IMIN, the value of the smallest entry.
!
  integer n
!
  integer i
  integer iarray(n)
  integer imin
  integer index
!
  if ( n <= 0 ) then

    imin = 0
    index = 0

  else

    imin = iarray(1)
    index = 1
    do i = 2, n
      if ( iarray(i) < imin ) then
        imin = iarray(i)
        index = i
      end if
    end do

  end if

  return
end
subroutine ivec_variance ( n, x, variance )
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
!    Input, integer X(N), the vector whose variance is desired.
!
!    Output, real VARIANCE, the variance of the vector entries.
!
  integer n
!
  integer i
  real mean
  real variance
  integer x(n)
!
  call ivec_mean ( n, x, mean )

  variance = 0.0
  do i = 1, n
    variance = variance + ( real ( x(i) ) - mean )**2
  end do

  if ( n > 1 ) then
    variance = variance / real ( n - 1 )
  else
    variance = 0.0
  end if

  return
end
subroutine laplace_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! LAPLACE_CDF evaluates the Laplace CDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the PDF.
!
  real a
  real b
  real cdf
  real x
  real y
!
  y = ( x - a ) / b

  if ( x <= a ) then
    cdf = 0.5 * exp ( y )
  else 
    cdf = 1.0 - 0.5 * exp ( - y )
  end if

  return
end
subroutine laplace_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! LAPLACE_CDF_INV inverts the Laplace CDF.
!
!
!  Modified:
!
!    17 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LAPLACE_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.5 ) then
    x = a + b * log ( 2.0 * cdf )
  else
    x = a - b * log ( 2.0 * ( 1.0 - cdf ) )
  end if

  return
end
subroutine laplace_check ( a, b )
!
!*******************************************************************************
!
!! LAPLACE_CHECK checks the parameters of the Laplace PDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LAPLACE_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine laplace_mean ( a, b, mean )
!
!*******************************************************************************
!
!! LAPLACE_MEAN returns the mean of the Laplace PDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine laplace_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! LAPLACE_PDF evaluates the Laplace PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = exp ( - abs ( X - A ) / B ) / ( 2 * B )
!
!  Discussion:
!
!    The Laplace PDF is also known as the Double Exponential PDF.
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  pdf = exp ( - abs ( x - a ) / b ) / ( 2.0 * b )

  return
end
subroutine laplace_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! LAPLACE_SAMPLE samples the Laplace PDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call laplace_cdf_inv ( cdf, a, b, x )

  return
end
subroutine laplace_variance ( a, b, variance )
!
!*******************************************************************************
!
!! LAPLACE_VARIANCE returns the variance of the Laplace PDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = 2.0 * b**2

  return
end
function lerch ( a, b, c )
!
!*******************************************************************************
!
!! LERCH estimates the Lerch transcendent function.
!
!
!  Definition:
!
!    The Lerch transcendent function is defined as:
!
!      LERCH ( A, B, C ) = Sum ( 0 <= K < Infinity ) A**K / ( C + K )**B
!
!    excluding any term with ( C + K ) = 0.
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998.
!
!  Modified:
!
!    17 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Thanks:
!
!    Oscar van Vlijmen
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the function. 
!
!    Output, real LERCH, an approximation to the Lerch transcendent function.
!
  real a
  double precision a_k
  real b
  real c
  integer k
  real lerch
  double precision sum
  double precision sum_old
!
  sum = 0.0
  k = 0
  a_k = 1.0

  do

    sum_old = sum

    if ( dble ( c ) + dble ( k ) == 0.0 ) then
      k = k + 1
      a_k = a_k * dble ( a )
      cycle
    end if
    
    sum = sum + a_k / ( dble ( c ) + dble ( k ) )**b

    if ( real ( sum ) <= real ( sum_old ) ) then
      exit
    end if

    k = k + 1
    a_k = a_k * dble ( a )

  end do

  lerch = real ( sum )

  return
end
subroutine logistic_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! LOGISTIC_CDF evaluates the Logistic CDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  cdf = 1.0 / ( 1.0 + exp ( ( a - x ) / b ) )

  return
end
subroutine logistic_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! LOGISTIC_CDF_INV inverts the Logistic CDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LOGISTIC_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( ( 1.0 - cdf ) / cdf )

  return
end
subroutine logistic_check ( a, b )
!
!*******************************************************************************
!
!! LOGISTIC_CHECK checks the parameters of the Logistic CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LOGISTIC_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine logistic_mean ( a, b, mean )
!
!*******************************************************************************
!
!! LOGISTIC_MEAN returns the mean of the Logistic PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine logistic_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! LOGISTIC_PDF evaluates the Logistic PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = EXP ( ( A - X ) / B ) /
!      ( B * ( 1 + EXP ( ( A - X ) / B ) )**2 )
!
!  Discussion:
!
!    The Logistic PDF is also known as the Sech-Squared PDF.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real temp
  real x
!
  temp = exp ( ( a - x ) / b )

  pdf = temp / ( b * ( 1.0 + temp )**2 )

  return
end
subroutine logistic_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! LOGISTIC_SAMPLE samples the Logistic PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call logistic_cdf_inv ( cdf, a, b, x )

  return
end
subroutine logistic_variance ( a, b, variance )
!
!*******************************************************************************
!
!! LOGISTIC_VARIANCE returns the variance of the Logistic PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real variance
!
  variance = ( PI * b )**2.0 / 3.0

  return
end
subroutine lognormal_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! LOGNORMAL_CDF evaluates the Lognormal CDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 < X.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real logx
  real x
!
  if ( x <= 0.0 ) then

    cdf = 0.0

  else

    logx = log ( x )
  
    call normal_cdf ( logx, a, b, cdf )

  end if

  return
end
subroutine lognormal_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! LOGNORMAL_CDF_INV inverts the Lognormal CDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real logx
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LOGNORMAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  call normal_cdf_inv ( cdf, a, b, logx )

  x = exp ( logx )

  return
end
subroutine lognormal_check ( a, b )
!
!*******************************************************************************
!
!! LOGNORMAL_CHECK checks the parameters of the Lognormal PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LOGNORMAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine lognormal_mean ( a, b, mean )
!
!*******************************************************************************
!
!! LOGNORMAL_MEAN returns the mean of the Lognormal PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = exp ( a + 0.5 * b**2 )

  return
end
subroutine lognormal_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! LOGNORMAL_PDF evaluates the Lognormal PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) 
!      = EXP ( - 0.5 * ( ( LOG ( X ) - A ) / B )**2 ) 
!        / ( B * X * SQRT ( 2 * PI ) )
!
!  Discussion:
!
!    The Lognormal PDF is also known as the Cobb-Douglas PDF,
!    and as the Antilognormal PDF.
!
!    The Lognormal PDF describes a variable X whose logarithm
!    is normally distributed.
!
!    The special case A = 0, B = 1 is known as Gilbrat's PDF.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 < X
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
!
  if ( x <= 0.0 ) then
    pdf = 0.0
  else
    pdf = exp ( - 0.5 * ( ( log ( x ) - a ) / b )**2 ) &
      / ( b * x * sqrt ( 2.0 * PI ) )
  end if

  return
end
subroutine lognormal_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! LOGNORMAL_SAMPLE samples the Lognormal PDF.
!
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call lognormal_cdf_inv ( cdf, a, b, x )

  return
end
subroutine lognormal_variance ( a, b, variance )
!
!*******************************************************************************
!
!! LOGNORMAL_VARIANCE returns the variance of the Lognormal PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = exp ( 2.0 * a + b**2 ) * ( exp ( b**2 ) - 1.0 )

  return
end
subroutine logseries_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! LOGSERIES_CDF evaluates the Logarithmic Series CDF.
!
!
!  Discussion:
!
!    Simple summation is used, with a recursion to generate successive
!    values of the PDF.
!
!  Modified:
!
!    18 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Thanks:
!
!    Oscar van Vlijmen
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!    0 < X 
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A < 1.0.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real cdf
  real pdf
  integer x
  integer x2
!
  cdf = 0.0

  do x2 = 1, x

    if ( x2 == 1 ) then
      pdf = - a / log ( 1.0 - a )
    else
      pdf = real ( x2 - 1 ) * a * pdf / real ( x2 )
    end if

    cdf = cdf + pdf

  end do

  return
end
subroutine logseries_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! LOGSERIES_CDF_INV inverts the Logarithmic Series CDF.
!
!
!  Discussion:
!
!    Simple summation is used.  The only protection against an
!    infinite loop caused by roundoff is that X cannot be larger 
!    than 1000.
!
!  Modified:
!
!    18 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A < 1.0.
!
!    Output, real X, the argument of the CDF for which
!    CDF(X-1) <= CDF <= CDF(X).
!
  integer, parameter :: XMAX = 1000
!
  real a
  real cdf
  real cdf2
  real pdf
  integer x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LOGSERIES_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  cdf2 = 0.0
  x = 1

  do while ( cdf > cdf2 .and. x < XMAX )

    if ( x == 1 ) then
      pdf = - a / log ( 1.0 - a )
    else
      pdf = real ( x - 1 ) * a * pdf / real ( x )
    end if

    cdf2 = cdf2 + pdf

    x = x + 1

  end do

  return
end
subroutine logseries_check ( a )
!
!*******************************************************************************
!
!! LOGSERIES_CHECK checks the parameter of the Logarithmic Series PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A < 1.0.
!
  real a
!
  if ( a <= 0.0 .or. a >= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LOGSERIES_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.0 or 1.0 <= A'
    stop
  end if

  return
end
subroutine logseries_mean ( a, mean )
!
!*******************************************************************************
!
!! LOGSERIES_MEAN returns the mean of the Logarithmic Series PDF.
!
!
!  Modified:
!
!    20 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A < 1.0.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real mean
!
  mean = - a / ( ( 1.0 - a ) * log ( 1.0 - a ) )

  return
end
subroutine logseries_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! LOGSERIES_PDF evaluates the Logarithmic Series PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = - A**X / ( X * log ( 1 - A ) )
!
!  Modified:
!
!    20 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!    0 < X
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A < 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real pdf
  integer x
!
  if ( x <= 0 ) then
    pdf = 0.0
  else
    pdf = - a**x / ( real ( x ) * log ( 1.0 - a ) )
  end if

  return
end
subroutine logseries_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! LOGSERIES_SAMPLE samples the Logarithmic Series PDF.
!
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer-Verlag, New York, 1986, page 547.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A < 1.0.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer K, a sample of the PDF.
!
  real a
  integer iseed
  real u
  real uniform_01_sample
  real v
  integer x
!
  u = uniform_01_sample ( iseed )
  v = uniform_01_sample ( iseed )

  x = int ( 1.0 + log ( v ) / ( log ( 1.0 - ( 1.0 - a )**u ) ) )

  return
end
subroutine logseries_variance ( a, variance )
!
!*******************************************************************************
!
!! LOGSERIES_VARIANCE returns the variance of the Logarithmic Series PDF.
!
!
!  Modified:
!
!    20 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A < 1.0.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real alpha
  real variance
!
  alpha = - 1.0 / log ( 1.0 - a )

  variance = a * alpha * ( 1.0 - alpha * a ) / ( 1.0 - a )**2

  return
end
subroutine lorentz_cdf ( x, cdf )
!
!*******************************************************************************
!
!! LORENTZ_CDF evaluates the Lorentz CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real cdf
  real x
!
  cdf = 0.5 + atan ( x ) / PI

  return
end
subroutine lorentz_cdf_inv ( cdf, x )
!
!*******************************************************************************
!
!! LORENTZ_CDF_INV inverts the Lorentz CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, real X, the corresponding argument.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'LORENTZ_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = tan ( PI * ( cdf - 0.5 ) )

  return
end
subroutine lorentz_mean ( mean )
!
!*******************************************************************************
!
!! LORENTZ_MEAN returns the mean of the Lorentz PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MEAN, the mean of the PDF.
!
  real mean
!
  mean = 0.0

  return
end
subroutine lorentz_pdf ( x, pdf )
!
!*******************************************************************************
!
!! LORENTZ_PDF evaluates the Lorentz PDF.
!
!
!  Formula:
!
!    PDF(X) = 1 / ( PI * ( 1 + X**2 ) )
!
!  Discussion:
!
!    The chief interest of the Lorentz PDF is that it is easily
!    inverted, and can be used to dominate other PDF's in an
!    acceptance/rejection method.
!
!    LORENTZ_PDF(X) = CAUCHY_PDF(X)(0,1)
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real pdf
  real x
!
  pdf = 1.0 / ( PI * ( 1.0 + x**2 ) )

  return
end
subroutine lorentz_sample ( iseed, x )
!
!*******************************************************************************
!
!! LORENTZ_SAMPLE samples the Lorentz PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call lorentz_cdf_inv ( cdf, x )

  return
end
subroutine lorentz_variance ( variance )
!
!*******************************************************************************
!
!! LORENTZ_VARIANCE returns the variance of the Lorentz PDF.
!
!
!  Discussion:
!
!    The variance of the Lorentz PDF is not well defined.  This routine
!    is made available for completeness only, and simply returns
!    a "very large" number.
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VARIANCE, the mean of the PDF.
!
  real variance
!
  variance = huge ( variance )

  return
end
subroutine maxwell_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! MAXWELL_CDF evaluates the Maxwell CDF.
!
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 <= X
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real cdf
  real gamma_inc
  real p2
  real x
  real x2
!
  if ( x <= 0.0 ) then

    cdf = 0.0

  else

    x2 = x / a
    p2 = 1.5

    cdf = gamma_inc ( x2, p2 )

  end if

  return
end
subroutine maxwell_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! MAXWELL_CDF_INV inverts the Maxwell CDF.
!
!
!  Discussion:
!
!    A simple bisection method is used.
!
!  Modified:
!
!    05 January 2000
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
!    Output, real X, the corresponding argument of the CDF.
!
  integer, parameter :: IT_MAX = 10
  real, parameter :: TOL = 0.0001
!
  real a
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
!
  if ( cdf <= 0.0 ) then
    x = 0.0
    return
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
    return
  end if

  x1 = 0.0
  cdf1 = 0.0

  x2 = 1.0

10    continue

  call maxwell_cdf ( x2, a, cdf2 )

  if ( cdf2 <= cdf ) then
    x2 = 2.0 * x2
    if ( x2 > 1000000.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'MAXWELL_CDF_INV - Fatal error!'
      write ( *, * ) '  Initial bracketing effort fails.'
      stop
    end if
    go to 10
  end if
!
!  Now use bisection.
!
  it = 0

20    continue

  it = it + 1

  x3 = 0.5 * ( x1 + x2 )
  call maxwell_cdf ( x3, a, cdf3 )

  if ( abs ( cdf3 - cdf ) < TOL ) then
    x = x3
    return
  end if

  if ( it > IT_MAX ) then
    write ( *, * ) ' '
    write ( *, * ) 'MAXWELL_CDF_INV - Fatal error!'
    write ( *, * ) '  Iteration limit exceeded.'
    stop
  end if

  if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
    x1 = x3
    cdf1 = cdf3
  else
    x2 = x3
    cdf2 = cdf3
  end if

  go to 20
end
subroutine maxwell_check ( a )
!
!*******************************************************************************
!
!! MAXWELL_CHECK checks the parameters of the Maxwell CDF.
!
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
  real a
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'MAXWELL_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.0.'
    stop
  end if

  return
end
subroutine maxwell_mean ( a, mean )
!
!*******************************************************************************
!
!! MAXWELL_MEAN returns the mean of the Maxwell PDF.
!
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
!    Output, real MEAN, the mean value.
!
  real a
  real gamma
  real mean
!
  mean = sqrt ( 2.0 ) * a * gamma ( 2.0 ) / gamma ( 1.5 ) 

  return
end
subroutine maxwell_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! MAXWELL_PDF evaluates the Maxwell PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = EXP ( - 0.5 * ( X / A )**2 ) * ( X / A )**2 /
!      ( SQRT ( 2 ) * A * GAMMA ( 1.5 ) )
!      
!  Discussion:
!
!    MAXWELL_PDF(X)(A) = CHI_PDF(0,A,3)
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0 < X
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real gamma
  real pdf
  real x
  real y
!
  if ( x <= 0.0 ) then

    pdf = 0.0

  else

    y = x / a

    pdf = exp ( - 0.5 * y**2 ) * y**2 / ( sqrt ( 2.0 ) * a * gamma ( 1.5 ) )

  end if

  return
end
subroutine maxwell_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! MAXWELL_SAMPLE samples the Maxwell PDF.
!
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real a2
  integer iseed
  real x
!
  a2 = 3.0
  call chisquare_central_sample ( a2, iseed, x )

  x = a * sqrt ( x )

  return
end
subroutine maxwell_variance ( a, variance )
!
!*******************************************************************************
!
!! MAXWELL_VARIANCE returns the variance of the Maxwell PDF.
!
!
!  Modified:
!
!    05 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real gamma
  real variance
!
  variance = a**2 * ( 3.0 - 2.0 * ( gamma ( 2.0 ) / gamma ( 1.5 ) )**2 )

  return
end
subroutine multicoef_check ( nfactor, factor )
!
!*******************************************************************************
!
!! MULTICOEF_CHECK checks the parameters of the multinomial coefficient.
!
!
!  Modified:
!
!    11 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NFACTOR, the number of factors.
!    1 <= NFACTOR.
!
!    Input, integer FACTOR(NFACTOR), contains the factors.
!    0.0 <= FACTOR(I).
!
  integer nfactor
!
  integer factor(nfactor)
  integer i
!
  if ( nfactor < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'MULTICOEF_CHECK - Fatal error!'
    write ( *, * ) '  NFACTOR < 1.'
    stop
  end if

  do i = 1, nfactor

    if ( factor(i) < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'MULTICOEF_CHECK - Fatal error'
      write ( *, * ) '  Factor ', I, ' = ', factor(i)
      write ( *, * ) '  But this value must be nonnegative.'
      stop
    end if

  end do

  return
end
subroutine multinomial_coef1 ( nfactor, factor, ncomb )
!
!*******************************************************************************
!
!! MULTINOMIAL_COEF1 computes a Multinomial coefficient.
!
!
!  Definition:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!  Formula:
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!  Method:
!
!    The log of the gamma function is used, to avoid overflow.
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
!    Input, integer NFACTOR, the number of factors.
!    1 <= NFACTOR.
!
!    Input, integer FACTOR(NFACTOR), contains the factors.
!    0.0 <= FACTOR(I).
!
!    Output, integer NCOMB, the value of the multinomial coefficient.
!
  integer nfactor
!
  real facn
  integer factor(nfactor)
  real factorial_log
  integer i
  integer n
  integer ncomb
!
  call multicoef_check ( nfactor, factor )
!
!  The factors sum to N.
!
  n = 0
  do i = 1, nfactor
    n = n + factor(i)
  end do

  facn = factorial_log ( n )

  do i = 1, nfactor

    facn = facn - factorial_log ( factor(i) )

  end do

  ncomb = nint ( exp ( facn ) )

  return
end
subroutine multinomial_coef2 ( nfactor, factor, ncomb )
!
!*******************************************************************************
!
!! MULTINOMIAL_COEF2 computes a Multinomial coefficient.
!
!
!  Definition:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!  Formula:
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!  Method:
!
!    A direct method is used, which should be exact.  However, there
!    is a possibility of intermediate overflow of the result.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NFACTOR, the number of factors.
!    1 <= NFACTOR.
!
!    Input, integer FACTOR(NFACTOR), contains the factors.
!    0.0 <= FACTOR(I).
!
!    Output, integer NCOMB, the value of the multinomial coefficient.
!
  integer nfactor
!
  integer factor(nfactor)
  integer i
  integer j
  integer k
  integer ncomb
!
  call multicoef_check ( nfactor, factor )

  ncomb = 1
  k = 0

  do i = 1, nfactor

    do j = 1, factor(i)
      k = k + 1
      ncomb = ( ncomb * k ) / j
    end do

  end do

  return
end
subroutine multinomial_check ( a, b, c )
!
!*******************************************************************************
!
!! MULTINOMIAL_CHECK checks the parameters of the Multinomial PDF.
!
!
!  Modified:
!
!    11 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of trials.
!
!    Input, integer B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0 <= C(I) <= 1.0,
!    SUM ( 1 <= I <= B ) C(I) = 1.0.
!
  integer b
!
  integer a
  real c(b)
  real c_sum
  integer i
!
  if ( b < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'MULTINOMIAL_CHECK - Fatal error!'
    write ( *, * ) '  B < 1.'
    stop
  end if

  do i = 1, b

    if ( c(i) < 0.0 .or. c(i) > 1.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'MULTINOMIAL_CHECK - Fatal error!'
      write ( *, * ) '  Input C(I) is out of range.'
      stop
    end if

  end do

  c_sum = sum ( c )

  if ( abs ( 1.0 - c_sum ) > 0.0001 ) then
    write ( *, * ) ' '
    write ( *, * ) 'MULTINOMIAL_CHECK - Fatal error!'
    write ( *, * ) '  The probabilities do not sum to 1.'
    stop
  end if

  return
end
subroutine multinomial_covariance ( a, b, c, covariance )
!
!*******************************************************************************
!
!! MULTINOMIAL_COVARIANCE returns the covariances of the Multinomial PDF.
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
!    Input, integer A, the number of trials.
!
!    Input, integer B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0 <= C(I) <= 1.0,
!    SUM ( 1 <= I <= B) C(I) = 1.0.
!
!    Output, real COVARIANCE(B,B), the covariance matrix.
!
  integer b
!
  integer a
  real c(b)
  real covariance(b,b)
  integer i
  integer j
!
  do i = 1, b
    do j = 1, b

      if ( i == j ) then
        covariance(i,j) = real ( a ) * c(i) * ( 1.0 - c(i) )
      else
        covariance(i,j) = - real ( a ) * c(i) * c(j)
      end if

    end do
  end do

  return
end
subroutine multinomial_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! MULTINOMIAL_MEAN returns the means of the Multinomial PDF.
!
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of trials.
!
!    Input, integer B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0 <= C(I) <= 1.0,
!    SUM ( 1 <= I <= B) C(I) = 1.0.
!
!    Output, real MEAN(B), MEAN(I) is the expected value of the 
!    number of outcome I in N trials.
!
  integer b
!
  integer a
  real c(b)
  integer i
  real mean(b)
!
  do i = 1, b
    mean(i) = real ( a ) * c(i)
  end do

  return
end
subroutine multinomial_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! MULTINOMIAL_PDF computes a Multinomial PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = Comb(A,B,X) * Product ( 1 <= I <= B ) C(I)**X(I)
!
!    where Comb(A,B,X) is the multinomial coefficient
!      C( A; X(1), X(2), ..., X(B) )
!
!  Discussion:
!
!    PDF(X)(A,B,C) is the probability that in A trials there
!    will be exactly X(I) occurrences of event I, whose probability
!    on one trial is C(I), for I from 1 to B.
!
!    As soon as A or B gets large, the number of possible X's explodes,
!    and the probability of any particular X can become extremely small.
!
!  Modified:
!
!    14 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X(B); X(I) counts the number of occurrences of
!    outcome I, out of the total of A trials.
!
!    Input, integer A, the total number of trials.
!
!    Input, integer B, the number of different possible outcomes on 
!    one trial.
!
!    Input, integer C(B); C(I) is the probability of outcome I on 
!    any one trial.
!
!    Output, real PDF, the value of the multinomial PDF.
!
  integer b
!
  integer a
  real c(b)
  real gamma_log
  integer i
  real pdf
  real pdf_log
  integer x(b)
!
!  To try to avoid overflow, do the calculation in terms of logarithms.
!  Note that Gamma(A+1) = A factorial.
!
  pdf_log = gamma_log ( real ( a + 1 ) )

  do i = 1, b
    pdf_log = pdf_log + x(i) * log ( c(i) ) - gamma_log ( real ( x(i) + 1 ) )
  end do

  pdf = exp ( pdf_log )

  return
end
subroutine multinomial_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! MULTINOMIAL_SAMPLE samples the Multinomial PDF.
!
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer-Verlag, New York, 1986, page 559.
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
!    Input, integer A, the total number of trials.
!    0 <= A.
!
!    Input, integer B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0 <= C(I) <= 1.0,
!    SUM ( 1 <= I <= B) C(I) = 1.0.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X(B); X(I) is the number of
!    occurrences of event I during the N trials.
!
  integer b
!
  integer a
  real c(b)
  integer i
  integer ifactor
  integer iseed
  integer ntot
  real prob
  real sum
  integer x(b)
!
  ntot = a

  sum = 1.0

  x(1:b) = 0

  do ifactor = 1, b - 1

    prob = c(ifactor) / sum
!
!  Generate a binomial random deviate for NTOT trials with 
!  single trial success probability PROB.
!
    call binomial_sample ( ntot, prob, iseed, x(ifactor) )

    ntot = ntot - x(ifactor)
    if ( ntot <= 0 ) then
      return
    end if

    sum = sum - c(ifactor)

  end do
!
!  The last factor gets what's left.
!
  x(b) = ntot

  return
end
subroutine multinomial_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! MULTINOMIAL_VARIANCE returns the variances of the Multinomial PDF.
!
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, the number of trials.
!
!    Input, integer B, the number of outcomes possible on one trial.
!    1 <= B.
!
!    Input, real C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0 <= C(I) <= 1.0,
!    SUM ( 1 <= I <= B ) C(I) = 1.0.
!
!    Output, real VARIANCE(B), VARIANCE(I) is the variance of the 
!    total number of events of type I.
!
  integer b
!
  integer a
  real c(b)
  integer i
  real variance(b)
!
  do i = 1, b
    variance(i) = real ( a ) * c(i) * ( 1.0 - c(i) )
  end do

  return
end
subroutine multivariate_normal_sample ( n, mean, covar_factor, iseed, x )
!
!*******************************************************************************
!
!! MULTIVARIATE_NORMAL_SAMPLE samples the Multivariate Normal PDF.
!
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 167.
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
!    Input, integer N, the dimension of the sphere.
!
!    Input, real MEAN(N), the mean values of the variates.
! 
!    Input, real COVAR_FACTOR(N,N), the lower triangular Cholesky
!    factor of the covariance matrix.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X(N), a point on the unit N sphere, chosen
!    with a uniform probability.
!
  integer n
!
  real covar_factor(n,n)
  integer i
  integer iseed
  integer j
  real mean(n)
  real x(n)
  real z
!
  do i = 1, n
    x(i) = mean(i)
    call normal_01_sample ( iseed, z )
    do j = 1, i
      x(i) = x(i) + covar_factor(i,j) * z
    end do
  end do

  return
end
subroutine nakagami_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! NAKAGAMI_CDF evaluates the Nakagami CDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B
!    0.0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real gamma_inc
  real p2
  real x
  real x2
  real y
!
  if ( x <= 0.0 ) then

    cdf = 0.0

  else if ( x > 0.0 ) then

    y = ( x - a ) / b 
    x2 = c * y**2
    p2 = c

    cdf = gamma_inc ( x2, p2 )

  end if

  return
end
subroutine nakagami_check ( a, b, c )
!
!*******************************************************************************
!
!! NAKAGAMI_CHECK checks the parameters of the Nakagami PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
  real a
  real b
  real c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NAKAGAMI_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NAKAGAMI_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    stop
  end if

  return
end
subroutine nakagami_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! NAKAGAMI_MEAN returns the mean of the Nakagami PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B
!    0.0 < C
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real gamma
  real mean
!
  mean = a + b * gamma ( c + 0.5 ) / ( sqrt ( c ) * gamma ( c ) )

  return
end
subroutine nakagami_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! NAKAGAMI_PDF evaluates the Nakagami PDF.
!
!
!  Modified:
!
!    03 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B
!    0.0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real gamma
  real pdf
  real x
  real y
!
  if ( x <= 0.0 ) then

    pdf = 0.0

  else if ( x > 0.0 ) then

    y = ( x - a ) / b 

    pdf = 2.0 * c**c / ( b * gamma ( c ) ) * y**( 2.0 * c - 1.0 ) &
      * exp ( - c * y**2 )

  end if

  return
end
subroutine nakagami_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! NAKAGAMI_VARIANCE returns the variance of the Nakagami PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B
!    0.0 < C
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real c
  real gamma
  real t1
  real t2
  real variance
!
  t1 = gamma ( c + 0.5 ) 
  t2 = gamma ( c )

  variance = b**2 * ( 1.0 - t1**2 / ( c * t2**2 ) )

  return
end
subroutine normal_01_cdf ( x, cdf )
!
!*******************************************************************************
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!
!  Reference: 
!
!    A G Adams,
!    Areas Under the Normal Curve,
!    Algorithm 39, 
!    Computer j., 
!    Volume 12, pages 197-198, 1969.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: A1 = 0.398942280444
  real, parameter :: A2 = 0.399903438504 
  real, parameter :: A3 = 5.75885480458
  real, parameter :: A4 = 29.8213557808
  real, parameter :: A5 = 2.62433121679
  real, parameter :: A6 = 48.6959930692
  real, parameter :: A7 = 5.92885724438

  real, parameter :: B0 = 0.398942280385
  real, parameter :: B1 = 3.8052E-08
  real, parameter :: B2 = 1.00000615302 
  real, parameter :: B3 = 3.98064794E-04
  real, parameter :: B4 = 1.98615381364
  real, parameter :: B5 = 0.151679116635 
  real, parameter :: B6 = 5.29330324926
  real, parameter :: B7 = 4.8385912808
  real, parameter :: B8 = 15.1508972451
  real, parameter :: B9 = 0.742380924027
  real, parameter :: B10 = 30.789933034
  real, parameter :: B11 = 3.99019417011
!
  real cdf
  real q
  real x
  real y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28 ) then

    y = 0.5 * x**2

    q = 0.5 - abs ( x ) * ( A1 - A2 * y / ( y + A3 - A4 / ( y + A5 &
      + A6 / ( y + A7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7 ) then

    y = 0.5 * x**2

    q = exp ( - y ) * B0 / ( abs ( x ) - B1 &
      + B2 / ( abs ( x ) + B3 &
      + B4 / ( abs ( x ) - B5 &
      + B6 / ( abs ( x ) + B7 &
      - B8 / ( abs ( x ) + B9 &
      + B10 / ( abs ( x ) + B11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0 ) then
    cdf = q
  else
    cdf = 1.0 - q
  end if

  return
end
subroutine normal_01_cdf_inv ( cdf, x )
!
!*******************************************************************************
!
!! NORMAL_01_CDF_INV inverts the Normal 01 CDF.
!
!
!  Reference:
!
!    Algorithm AS 111, 
!    Applied Statistics, 
!    Volume 26, pages 118-121, 1977.
!
!  Modified:
!
!    23 February 1999
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real X, the corresponding argument.
!
  real, parameter :: A0 =   2.50662823884
  real, parameter :: A1 = -18.61500062529
  real, parameter :: A2 =  41.39119773534
  real, parameter :: A3 = -25.44106049637
  real, parameter :: B1 =  -8.47351093090
  real, parameter :: B2 =  23.08336743743
  real, parameter :: B3 = -21.06224101826
  real, parameter :: B4 =   3.13082909833
  real, parameter :: C0 =  -2.78718931138
  real, parameter :: C1 =  -2.29796479134
  real, parameter :: C2 =   4.85014127135
  real, parameter :: C3 =   2.32121276858
  real, parameter :: D1 =   3.54388924762
  real, parameter :: D2 =   1.63706781897
!
  real cdf
  real q
  real r
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NORMAL_01_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  q = cdf - 0.5

  q = min ( q, 0.5 )
  q = max ( q, -0.5 )
!
!  0.08 < CDF < 0.92
!
  if ( abs ( q ) <= 0.42 ) then

    r = q**2

    x = q * ( ( ( &
           A3   * r &
         + A2 ) * r &
         + A1 ) * r &
         + A0 ) / ( ( ( ( &
           B4   * r &
         + B3 ) * r &
         + B2 ) * r &
         + B1 ) * r + 1.0 )
!
!  CDF < 0.08 or 0.92 < CDF.
!
  else

    r = min ( cdf, 1.0 - cdf )
    r = max ( r, 0.0 )
    r = min ( r, 1.0 )

    r = sqrt ( - log ( r ) )

    x = ( ( ( &
           C3   * r &
         + C2 ) * r &
         + C1 ) * r &
         + C0 ) / ( ( &
           D2   * r &
         + D1 ) * r + 1.0 )

    if ( q < 0.0 ) then
      x = - x
    end if

  end if

  return
end
subroutine normal_01_mean ( mean )
!
!*******************************************************************************
!
!! NORMAL_01_MEAN returns the mean of the Normal 01 PDF.
!
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MEAN, the mean of the PDF.
!
  real mean
!
  mean = 0.0

  return
end
subroutine normal_01_pdf ( x, pdf )
!
!*******************************************************************************
!
!! NORMAL_01_PDF evaluates the Normal 01 PDF.
!
!
!  Discussion:
!
!    The Normal 01 PDF is also called the "Standard Normal" PDF, or
!    the Normal PDF with 0 mean and variance 1.
!
!  Formula:
!
!    PDF(X) = EXP ( - 0.5 * X**2 ) / SQRT ( 2 * PI )
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real pdf
  real x
!
  pdf = exp ( - 0.5 * x**2 ) / sqrt ( 2.0 * PI )

  return
end
subroutine normal_01_sample ( iseed, x )
!
!*******************************************************************************
!
!! NORMAL_01_SAMPLE samples the Normal 01 PDF.
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
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  integer iseed
  integer, save :: iset = 0
  real uniform_01_sample
  real v1
  real v2
  real x
  real, save :: xsave = 0.0
!
  if ( iset == 0 ) then

    v1 = uniform_01_sample ( iseed )

    if ( v1 <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V1 <= 0.'
      write ( *, * ) '  V1 = ', v1
      write ( *, * ) '  ISEED = ', iseed
      stop
    end if
    
    v2 = uniform_01_sample ( iseed )

    if ( v2 <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V2 <= 0.'
      write ( *, * ) '  V2 = ', v2
      write ( *, * ) '  ISEED = ', iseed
      stop
    end if

    x = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * PI * v2 )

    xsave = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
subroutine normal_01_variance ( variance )
!
!*******************************************************************************
!
!! NORMAL_01_VARIANCE returns the variance of the Normal 01 PDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real variance
!
  variance = 1.0

  return
end
subroutine normal_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! NORMAL_CDF evaluates the Normal CDF.
!
!
!  Discussion:
!
!    The Normal CDF is related to the Error Function ERF(X) by:
!
!      ERF ( X ) = 2 * NORMAL_CDF ( SQRT ( 2 ) * X ) - 1.0.
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
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
  real y
!
  y = ( x - a ) / b

  call normal_01_cdf ( y, cdf )

  return
end
subroutine normal_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! NORMAL_CDF_INV inverts the Normal CDF.
!
!
!  Reference:
!
!    Algorithm AS 111, 
!    Applied Statistics, 
!    Volume 26, pages 118-121, 1977.
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
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
  real x2
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NORMAL_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  call normal_01_cdf_inv ( cdf, x2 )

  x = a + b * x2

  return
end
subroutine normal_check ( a, b )
!
!*******************************************************************************
!
!! NORMAL_CHECK checks the parameters of the Normal PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NORMAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine normal_mean ( a, b, mean )
!
!*******************************************************************************
!
!! NORMAL_MEAN returns the mean of the Normal PDF.
!
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine normal_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! NORMAL_PDF evaluates the Normal PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) 
!      = EXP ( - 0.5 * ( ( X - A ) / B )**2 )
!      / ( B * SQRT ( 2 * PI ) )
!
!  Discussion:
!
!    The normal PDF is also known as the Gaussian PDF.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
  real y
!
  y = ( x - a ) / b

  pdf = exp ( - 0.5 * y**2 )  / ( b * sqrt ( 2.0 * PI ) )

  return
end
subroutine normal_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! NORMAL_SAMPLE samples the Normal PDF.
!
!
!  Method:
!
!    The Box-Muller method is used.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  integer iseed
  integer, save :: iset = 0
  real uniform_01_sample
  real v1
  real v2
  real x
  real, save :: xsave = 0.0
!
  if ( iset == 0 ) then

    v1 = uniform_01_sample ( iseed )

    if ( v1 <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_SAMPLE - Fatal error!'
      write ( *, * ) '  V1 <= 0.'
      write ( *, * ) '  V1 = ', v1
      write ( *, * ) '  ISEED = ', iseed
      stop
    end if
    
    v2 = uniform_01_sample ( iseed )

    if ( v2 <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_SAMPLE - Fatal error!'
      write ( *, * ) '  V2 <= 0.'
      write ( *, * ) '  V2 = ', v2
      write ( *, * ) '  ISEED = ', iseed
      stop
    end if

    x = a + b * sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * PI * v2 )

    xsave = a + b * sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
subroutine normal_variance ( a, b, variance )
!
!*******************************************************************************
!
!! NORMAL_VARIANCE returns the variance of the Normal PDF.
!
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = b**2

  return
end
subroutine pareto_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! PARETO_CDF evaluates the Pareto CDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( x < a ) then
    cdf = 0.0
  else
    cdf = 1.0 - ( a / x )**b
  end if

  return
end
subroutine pareto_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! PARETO_CDF_INV inverts the Pareto CDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PARETO_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a / ( 1.0 - cdf )**( 1.0 / b )

  return
end
subroutine pareto_check ( a, b )
!
!*******************************************************************************
!
!! PARETO_CHECK checks the parameters of the Pareto CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
  real a
  real b
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PARETO_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PARETO_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine pareto_mean ( a, b, mean )
!
!*******************************************************************************
!
!! PARETO_MEAN returns the mean of the Pareto PDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  if ( b <= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PARETO_MEAN - Fatal error!'
    write ( *, * ) '  For B <= 1, the mean does not exist.'
    mean = 0.0
    return
  end if

  mean = b * a / ( b - 1.0 )

  return
end
subroutine pareto_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! PARETO_PDF evaluates the Pareto PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = B * A**B / X**(B+1).
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  if ( x < a ) then
    pdf = 0.0
  else
    pdf = b * a**b / x**( b + 1.0 )
  end if

  return
end
subroutine pareto_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! PARETO_SAMPLE samples the Pareto PDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call pareto_cdf_inv ( cdf, a, b, x )

  return
end
subroutine pareto_variance ( a, b, variance )
!
!*******************************************************************************
!
!! PARETO_VARIANCE returns the variance of the Pareto PDF.
!
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  if ( b <= 2.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PARETO_VARIANCE - Warning!'
    write ( *, * ) '  For B <= 2, the variance does not exist.'
    variance = 0.0
    return
  end if

  variance = a**2 * b / ( ( b - 1.0 )**2 * ( b - 2.0 ) )

  return
end
subroutine pascal_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! PASCAL_CDF evaluates the Pascal CDF.
!
!
!  Discussion:
!
!    A simple-minded summing approach is used.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the CDF.
!
!    Input, integer A, real B, parameters of the PDF.
!    0 <= A,
!    0 < B <= 1.
!
!    Output, real CDF, the value of the CDF.
!
  integer a
  real b
  real cdf
  integer cnk
  real pdf
  integer x
  integer y
!
  cdf = 0.0

  do y = a, x

    call binomial_coef ( y-1, a-1, cnk )

    pdf = real ( cnk ) * b**a * ( 1.0 - b )**( y - a )

    cdf = cdf + pdf

  end do

  return
end
subroutine pascal_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! PASCAL_CDF_INV inverts the Pascal CDF.
!
!
!  Discussion:
!
!    A simple-minded discrete approach is used.
!
!  Modified:
!
!    06 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, integer A, real B, parameters of the PDF.
!    0 <= A,
!    0 < B <= 1.
!
!    Output, integer X, the smallest X whose cumulative density function
!    is greater than or equal to CDF.
!
  integer, parameter :: x_max = 1000
!
  integer a
  real b
  real cdf
  real cum
  real pdf
  integer x
!
  if ( cdf <= 0.0 ) then

    x = a

  else

    cum = 0.0

    x = a

10      continue

    call pascal_pdf ( x, a, b, pdf )

    cum = cum + pdf

    if ( cum < cdf .and. x < x_max ) then
      x = x + 1
      go to 10
    end if

  end if

  return
end
subroutine pascal_check ( a, b )
!
!*******************************************************************************
!
!! PASCAL_CHECK checks the parameters of the Pascal PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, real B, parameters of the PDF.
!    0 <= A,
!    0 < B <= 1.
!
  integer a
  real b
!
  if ( a < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PASCAL_CHECK - Fatal error!'
    write ( *, * ) '  A < 0.'
    stop
  end if

  if ( b <= 0.0 .or. b > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PASCAL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0 or 1 < B.'
    stop
  end if

  return
end
subroutine pascal_mean ( a, b, mean )
!
!*******************************************************************************
!
!! PASCAL_MEAN returns the mean of the Pascal PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, real B, parameters of the PDF.
!    0 <= A,
!    0 < B <= 1.
!
!    Output, real MEAN, the mean of the PDF.
!
  integer a
  real b
  real mean
!
  mean = real ( a ) / b

  return
end
subroutine pascal_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! PASCAL_PDF evaluates the Pascal PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = C(X-1,A-1) * B**A * ( 1 - B )**(X-A)
!
!  Discussion:
!
!    PDF(X)(A,B) is the probability that the A-th success will
!    occur on the X-th trial, given that the probability
!    of a success on a single trial is B.
!
!    The Pascal PDF is also known as the "negative binomial" PDF or 
!    the "Polya" PDF.
!
!    PASCAL_PDF(X)(1,B) = GEOMETRIC_PDF(X)(B)
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the number of trials.
!    A <= X.
!
!    Input, integer A, the number of successes required.
!    0 <= A <= X, normally.
!
!    Input, real B, the probability of a success on a single trial.
!    0.0 < B <= 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  integer a
  real b
  integer cnk
  real pdf
  integer x
!
  if ( x < a ) then

    pdf = 0.0
    
  else

    call binomial_coef ( x-1, a-1, cnk )

    pdf = real ( cnk ) * b**a * ( 1.0 - b )**( x - a )

  end if

  return
end
subroutine pascal_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! PASCAL_SAMPLE samples the Pascal PDF.
!
!
!  Modified:
!
!    28 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, parameters of the PDF.
!    0 <= A,
!    0 < B <= 1.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer X, a sample of the PDF.
!
  integer, parameter :: INTMAX = 2147483647
!
  integer a
  real b
  integer iseed
  integer num_success
  real r
  real uniform_01_sample
  integer x
!
  if ( b == 1.0 ) then
    x = a
    return
  else if ( b == 0.0 ) then
    x = INTMAX
    return
  end if

  x = 0
  num_success = 0

  do while ( num_success < a )

    x = x + 1
    r = uniform_01_sample ( iseed )

    if ( r <= b ) then
      num_success = num_success + 1
    end if

  end do

  return
end
subroutine pascal_variance ( a, b, variance )
!
!*******************************************************************************
!
!! PASCAL_VARIANCE returns the variance of the Pascal PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, parameters of the PDF.
!    0 <= A,
!    0 < B <= 1.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer a
  real b
  real variance
!
  variance = real ( a ) * ( 1.0 - b ) / b**2

  return
end
subroutine pearson_05_check ( a, b, c )
!
!*******************************************************************************
!
!! PEARSON_05_CHECK checks the parameters of the Pearson 5 PDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
  real a
  real b
  real c
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PEARSON_05_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PEARSON_05_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine pearson_05_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! PEARSON_05_MEAN evaluates the mean of the Pearson 5 PDF.
!
!
!  Discussion:
!
!    The mean is undefined for B <= 1.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real mean
!     
  if ( b <= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PEARSON_05_MEAN - Warning!'
    write ( *, * ) '  MEAN undefined for B <= 1.'
    mean = c
    return
  end if

  mean = c + a / ( b - 1.0 )
  
  return
end
subroutine pearson_05_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! PEARSON_05_PDF evaluates the Pearson 5 PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = A**B * ( X - C )**(-B-1) 
!      * exp ( - A / ( X - C ) ) / Gamma ( B )
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    C < X
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real gamma
  real pdf
  real x
!
  if ( x <= c ) then
    pdf = 0.0
  else
    pdf = a**b * ( x - c )**(-b-1.0) * exp ( - a / ( x - c ) ) / gamma ( b )
  end if

  return
end
subroutine pearson_05_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! PEARSON_05_SAMPLE samples the Pearson 5 PDF.
!
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real a2
  real b
  real b2
  real c
  real c2
  integer iseed
  real x
  real x2
!
  a2 = 0.0
  b2 = b
  c2 = 1.0 / a

  call gamma_sample ( a2, b2, c2, iseed, x2 )

  x = c + 1.0 / x2
  
  return
end
function pi ( )
!
!*******************************************************************************
!
!! PI returns the value of Pi.
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
  pi = 3.14159265358979323846264338327950288419716939937510

  return
end
subroutine planck_mean ( mean )
!
!*******************************************************************************
!
!! PLANCK_MEAN returns the mean of the Planck PDF.
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
!    Output, real MEAN, the mean of the PDF.
!
  real mean
!
  mean = 3.8322

  return
end
subroutine planck_pdf ( x, pdf )
!
!*******************************************************************************
!
!! PLANCK_PDF evaluates the Planck PDF.
!
!
!  Discussion:
!
!    The Planck PDF describes the distribution of frequencies in
!    blackbody radiation at a given temperature T.
!
!    A generalization of the Planck PDF is:
!
!      PDF(X)(A,B) = B**( A + 1 ) * X**A / 
!        ( Gamma ( A + 1 ) * Zeta ( A + 1 ) * ( EXP ( B * X ) - 1 ) )
!
!    The standard PDF has A = 3, B = 1.
!
!  Formula:
!
!    PDF(X) = ( 15 / PI**4 ) * X**3 / ( EXP ( X ) - 1 )
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
!    Input, real X, the argument of the PDF.
!    0.0 <= X
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real pdf
  real x
!
  if ( x < 0.0 ) then
    pdf = 0.0
  else
    pdf = ( 15.0 / PI**4 ) * x**3 / ( exp ( x ) - 1.0 )
  end if

  return
end
subroutine planck_sample ( iseed, x )
!
!*******************************************************************************
!
!! PLANCK_SAMPLE samples the Planck PDF.
!
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer Verlag, 1986, pages 552.
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
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real a2
  real b
  real b2
  real c2
  real g
  integer iseed
  real x
  integer z
!
!  This sampling technique should work for the general Planck PDF,
!  described by the two parameters A and B.
!
  a = 3.0
  b = 1.0

  a2 = 0.0
  b2 = 1.0
  c2 = a + 1.0

  call gamma_sample ( a2, b2, c2, iseed, g )

  call zipf_sample ( c2, iseed, z )

  x = g / ( b * real ( z ) )

  return
end
subroutine planck_variance ( variance )
!
!*******************************************************************************
!
!! PLANCK_VARIANCE returns the variance of the Planck PDF.
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
!    Output, real VARIANCE, the variance of the PDF.
!
  real variance
!
  variance = ( 2.0281 )**2

  return
end
subroutine poisson_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! POISSON_CDF evaluates the Poisson CDF.
!
!
!  Definition:
!
!    CDF(X,A) is the probability that the number of events observed
!    in a unit time period will be no greater than X, given that the 
!    expected number of events in a unit time period is A.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the CDF.
!    X >= 0.
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real cdf
  integer i
  real last
  real new
  real sum
  integer x
!
  if ( x < 0 ) then

    cdf = 0.0

  else

    new = exp ( - a )
    sum = new

    do i = 1, x
      last = new
      new = last * a / real ( i )
      sum = sum + new
    end do

    cdf = sum

  end if

  return
end
subroutine poisson_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! POISSON_CDF_INV inverts the Poisson CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, a value of the CDF.
!    0 <= CDF < 1.
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, integer X, the corresponding argument.
!
  integer, parameter :: XMAX = 100
!
  real a
  real cdf
  integer i
  real last
  real new
  real sum
  real sumold
  integer x
!
!  Now simply start at X = 0, and find the first value for which
!  CDF(X-1) <= CDF <= CDF(X).
!
  sum = 0.0

  do i = 0, XMAX

    sumold = sum

    if ( i == 0 ) then
      new = exp ( - a )
      sum = new 
    else
      last = new
      new = last * a / real ( i )
      sum = sum + new
    end if

    if ( sumold <= cdf .and. cdf <= sum ) then
      x = i
      return
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'POISSON_SAMPLE - Warning!'
  write ( *, * ) '  Exceeded XMAX = ', XMAX

  x = XMAX

  return
end
subroutine poisson_check ( a )
!
!*******************************************************************************
!
!! POISSON_CHECK checks the parameter of the Poisson PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
  real a
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'POISSON_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  return
end
subroutine poisson_mean ( a, mean )
!
!*******************************************************************************
!
!! POISSON_MEAN returns the mean of the Poisson PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real mean
!     
  mean = a

  return
end
subroutine poisson_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! POISSON_PDF evaluates the Poisson PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = EXP ( - A ) * A**X / X!
!
!  Discussion:
!
!    PDF(X)(A) is the probability that the number of events observed
!    in a unit time period will be X, given the expected number 
!    of events in a unit time.
!
!    The parameter A is the expected number of events per unit time.
!
!    The Poisson PDF is a discrete version of the Exponential PDF.
!
!    The time interval between two Poisson events is a random 
!    variable with the Exponential PDF.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!    0 <= X
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real pdf
  real r_factorial
  integer x
!
  if ( x < 0 ) then
    pdf = 0.0
  else
    pdf = exp ( - a ) * a**x / r_factorial ( x )
  end if

  return
end
subroutine poisson_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! POISSON_SAMPLE samples the Poisson PDF.
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
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  real a
  real cdf
  integer iseed
  real uniform_01_sample
  integer x
!
  cdf = uniform_01_sample ( iseed )

  call poisson_cdf_inv ( cdf, a, x )

  return
end
subroutine poisson_variance ( a, variance )
!
!*******************************************************************************
!
!! POISSON_VARIANCE returns the variance of the Poisson PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real variance
!
  variance = a

  return
end
subroutine poly_val ( n, array, x, val )
!
!*******************************************************************************
!
!! POLY_VAL evaluates a polynomial using Horner's method.
!
!
!  Formula:
!
!    p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!         = a(0) + x * (
!           a(1) + x * (
!           ...
!           a(n-1) + x * (
!           a(n) ) ... ) )
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the polynomial.
!
!    Input, real ARRAY(0:N), the coefficients of the polynomial.
!    ARRAY(0) is the constant term and A(N) the coefficient of X**N.
!
!    Input, real X, the evaluation point.
!
!    Output, real VAL, the value of the polynomial at X.
!
  integer n
!
  real array(0:n)
  integer i
  real val
  real x
!
  val = array(n)
  do i = n - 1, 0, -1
    val = array(i) + x * val
  end do

  return
end
subroutine power_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! POWER_CDF evaluates the Power CDF.
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
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B,
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( x <= 0.0 ) then
    cdf = 0.0
  else if ( x <= b ) then
    cdf = ( x / b )**a
  else
    cdf = 1.0
  end if

  return
end
subroutine power_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! POWER_CDF_INV inverts the Power CDF.
!
!
!  Modified:
!
!    11 July 2000
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real X, the argument of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf <= 0.0 ) then
    x = 0.0
  else if ( cdf < 1.0 ) then
    x = b * exp ( log ( cdf ) / a )
  else
    x = b
  end if

  return
end
subroutine power_check ( a, b )
!
!*******************************************************************************
!
!! POWER_CHECK checks the parameter of the Power PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
  real a
  real b
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'POWER_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'POWER_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine power_mean ( a, b, mean )
!
!*******************************************************************************
!
!! POWER_MEAN returns the mean of the Power PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a * b / ( a + 1.0 )

  return
end
subroutine power_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! POWER_PDF evaluates the Power PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = (A/B) * (X/B)**(A-1)
!
!  Reference:
!
!    Daniel Zwillinger and Stephen Kokoska,
!    CRC Standard Probability and Statistics Tables and Formulae,
!    Chapman and Hall/CRC, 2000, pages 152-153.
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
!    Input, real X, the argument of the PDF.
!    0.0 <= X <= B.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  if ( x < 0.0 .or. x > b ) then
    pdf = 0.0
  else
    pdf = ( a / b ) * ( x / b )**( a - 1.0 )
  end if

  return
end
subroutine power_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! POWER_SAMPLE samples the Power PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )        

  call power_cdf_inv ( cdf, a, b, x )

  return
end
subroutine power_variance ( a, b, variance )
!
!*******************************************************************************
!
!! POWER_VARIANCE returns the variance of the Power PDF.
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A, 0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = b**2 * a / ( ( a + 1.0 )**2 * ( a + 2.0 ) )

  return
end
function r_factorial ( n )
!
!*******************************************************************************
!
!! R_FACTORIAL returns N!.
!
!
!  Definition:
!
!    N! = Product ( 1 <= I <= N ) I
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
!    Input, integer N, the argument of the function.
!    0 <= N.
!
!    Output, real R_FACTORIAL, the factorial of N.
!
  integer i
  integer n
  real r_factorial
!
  if ( n < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'R_FACTORIAL - Fatal error!'
    write ( *, * ) '  N < 0.'
    stop
  end if

  r_factorial = 1.0

  do i = 2, n
    r_factorial = r_factorial * real ( i )
  end do

  return
end
subroutine r_random ( rlo, rhi, iseed, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, real R, the randomly chosen value.
!
  integer iseed
  real r
  real rhi
  real rlo
  real t
  real uniform_01_sample
!
!  Pick a random number in (0,1).
!
  t = uniform_01_sample ( iseed )
!
!  Set R.
!
  r = ( 1.0 - t ) * rlo + t * rhi

  return
end
subroutine rayleigh_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! RAYLEIGH_CDF evaluates the Rayleigh CDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!    0.0 <= X.
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real cdf
  real x
!
  if ( x < 0.0 ) then
    cdf = 0.0
  else
    cdf = 1.0 - exp ( - x**2 / ( 2.0 * a**2 ) )
  end if

  return
end
subroutine rayleigh_cdf_inv ( cdf, a, x )
!
!*******************************************************************************
!
!! RAYLEIGH_CDF_INV inverts the Rayleigh CDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, real X, the corresponding argument.
!
  real a
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RAYLEIGH_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = sqrt ( - 2.0 * a**2 * log ( 1.0 - cdf ) ) 

  return
end
subroutine rayleigh_check ( a )
!
!*******************************************************************************
!
!! RAYLEIGH_CHECK checks the parameter of the Rayleigh PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
  real a
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RAYLEIGH_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.'
    stop
  end if

  return
end
subroutine rayleigh_mean ( a, mean )
!
!*******************************************************************************
!
!! RAYLEIGH_MEAN returns the mean of the Rayleigh PDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Output, real MEAN, the mean of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real mean
!
  mean = a * sqrt ( 0.5 * PI )

  return
end
subroutine rayleigh_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! RAYLEIGH_PDF evaluates the Rayleigh PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = ( X / A**2 ) * EXP ( - X**2 / ( 2 * A**2 ) )
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 <= X
!
!    Input, real A, the parameter of the PDF.
!    0 < A.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real pdf
  real x
!
  if ( x < 0.0 ) then
    pdf = 0.0
  else
    pdf = ( x / a**2 ) * exp ( - x**2 / ( 2.0 * a**2 ) )
  end if

  return
end
subroutine rayleigh_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! RAYLEIGH_SAMPLE samples the Rayleigh PDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    0.0 < A.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call rayleigh_cdf_inv ( cdf, a, x )

  return
end
subroutine rayleigh_variance ( a, variance )
!
!*******************************************************************************
!
!! RAYLEIGH_VARIANCE returns the variance of the Rayleigh PDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameters of the PDF.
!    0.0 < A.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real variance
!
  variance = 2.0 * a**2 * ( 1.0 - 0.25 * PI )

  return
end
subroutine reciprocal_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! RECIPROCAL_CDF evaluates the Reciprocal CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A <= B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( x <= 0.0 ) then

    cdf = 0.0

  else if ( x > 0.0 ) then

    cdf = log ( a / x ) / log ( a / b )

  end if

  return
end
subroutine reciprocal_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! RECIPROCAL_CDF_INV inverts the Reciprocal CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A <= B.
!
!    Output, real X, the corresponding argument of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf <= 0.0 ) then
    x = 0.0
  else if ( cdf > 0.0 ) then
    x = b**cdf / a**( cdf - 1.0 )
  end if

  return
end
subroutine reciprocal_check ( a, b )
!
!*******************************************************************************
!
!! RECIPROCAL_CHECK checks the parameters of the Reciprocal CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A <= B.
!
  real a
  real b
!
  if ( a <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RECIPROCAL_CHECK - Fatal error!'
    write ( *, * ) '  A <= 0.0'
    stop
  end if

  if ( b < a ) then
    write ( *, * ) ' '
    write ( *, * ) 'RECIPROCAL_CHECK - Fatal error!'
    write ( *, * ) '  B < A'
    stop
  end if

  return
end
subroutine reciprocal_mean ( a, b, mean )
!
!*******************************************************************************
!
!! RECIPROCAL_MEAN returns the mean of the Reciprocal PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A <= B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = ( a - b ) / log ( a / b )

  return
end
subroutine reciprocal_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! RECIPROCAL_PDF evaluates the Reciprocal PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = 1.0 / ( X * LOG ( B / A ) )
!    for 0.0 <= X
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A <= B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  if ( x <= 0.0 ) then
    pdf = 0.0
  else if ( x > 0.0 ) then
    pdf = 1.0 / ( x * log ( b / a ) )
  end if

  return
end
subroutine reciprocal_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! RECIPROCAL_SAMPLE samples the Reciprocal PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A <= B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  x = b**cdf / a**( cdf - 1.0 )

  return
end
subroutine reciprocal_variance ( a, b, variance )
!
!*******************************************************************************
!
!! RECIPROCAL_VARIANCE returns the variance of the Reciprocal PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < A <= B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real d
  real variance
!
  d = log ( a / b )

  variance = ( a - b )* ( a * ( d - 2.0 ) + b * ( d + 2.0 ) ) / ( 2.0 * d**2 )

  return
end
subroutine rrow_max ( lda, m, n, x, ixmax, xmax )
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
!    Input, integer LDA, the leading dimension of X, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real X(LDA,N), the array to be examined.
!
!    Output, integer IXMAX(M); IXMAX(I) is the column of X in which
!    the maximum for row I occurs.
!
!    Output, real XMAX(M), the maximums of the rows of X.
!
  integer lda
  integer m
  integer n
!
  integer i
  integer ixmax(m)
  integer j
  real x(lda,n)
  real xmax(m)
!
  do i = 1, m

    ixmax(i) = 1
    xmax(i) = x(i,1)
    do j = 2, n
      if ( x(i,j) > xmax(i) ) then
        ixmax(i) = j
        xmax(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine rrow_mean ( lda, m, n, x, mean )
!
!*******************************************************************************
!
!! RROW_MEAN returns the means of rows of a real array.
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
!    Input, integer LDA, the leading dimension of X, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real X(LDA,N), the array whose row means are desired.
!
!    Output, real MEAN(M), the means, or averages, of the rows of X.
!
  integer lda
  integer m
  integer n
!
  integer i
  integer j
  real mean(m)
  real x(lda,n)
!
  do i = 1, m

    mean(i) = 0.0
    do j = 1, n
      mean(i) = mean(i) + x(i,j)
    end do

    mean(i) = mean(i) / real ( n )

  end do

  return
end
subroutine rrow_min ( lda, m, n, x, ixmin, xmin )
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
!    Input, integer LDA, the leading dimension of X, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real X(LDA,N), the array to be examined.
!
!    Output, integer IXMIN(M); IXMIN(I) is the column of X in which
!    the minimum for row I occurs.
!
!    Output, real XMIN(M), the minimums of the rows of X.
!
  integer lda
  integer m
  integer n
!
  integer i
  integer ixmin(m)
  integer j
  real x(lda,n)
  real xmin(m)
!
  do i = 1, m

    ixmin(i) = 1
    xmin(i) = x(i,1)
    do j = 2, n
      if ( x(i,j) < xmin(i) ) then
        ixmin(i) = j
        xmin(i) = x(i,j)
      end if
    end do

  end do

  return
end
subroutine rrow_variance ( lda, m, n, x, variance )
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
!    Input, integer LDA, the leading dimension of X, which should
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, real X(LDA,N), the array whose row means are desired.
!
!    Output, real VARIANCE(M), the variances of the rows of X.
!
  integer lda
  integer m
  integer n
!
  integer i
  integer j
  real mean
  real variance(m)
  real x(lda,n)
!
  do i = 1, m

    mean = 0.0
    do j = 1, n
      mean = mean + x(i,j)
    end do
    mean = mean / real ( n )

    variance(i) = 0.0
    do j = 1, n
      variance(i) = variance(i) + ( x(i,j) - mean )**2
    end do

    if ( n > 1 ) then
      variance(i) = variance(i) / real ( n - 1 )
    else
      variance(i) = 0.0
    end if

  end do

  return
end
subroutine rvec_max ( n, x, index, xmax )
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
!    Input, real X(N), the array.
!
!    Output, integer INDEX, the index of the largest entry.
!
!    Output, real XMAX, the value of the largest entry.
!
  integer n
!
  integer i
  integer index
  real x(n)
  real xmax
!
  if ( n <= 0 ) then

    index = 0
    xmax = 0.0

  else

    index = 1
    xmax = x(1)

    do i = 2, n
      if ( x(i) > xmax ) then
        xmax = x(i)
        index = i
      end if
    end do

  end if

  return
end
subroutine rvec_mean ( n, x, mean )
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
!    Input, real X(N), the vector whose mean is desired.
!
!    Output, real MEAN, the mean, or average, of the vector entries.
!
  integer n
!
  integer i
  real mean
  real x(n)
!
  mean = sum ( x ) / real ( n )

  return
end
subroutine rvec_min ( n, x, index, xmin )
!
!*******************************************************************************
!
!! RVEC_MIN returns the minimum value of a real array.
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
!    Output, integer INDEX, the index of the smallest entry.
!
!    Output, real XMIN, the value of the smallest entry.
!
  integer n
!
  integer i
  integer index
  real x(n)
  real xmin
!
  if ( n <= 0 ) then

    index = 0
    xmin = 0.0

  else

    xmin = x(1)
    index = 1
    do i = 2, n
      if ( x(i) < xmin ) then
        xmin = x(i)
        index = i
      end if
    end do

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
!    08 July 2000
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
  integer i
!
  a_sum = sum ( a )

  if ( a_sum == 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RVEC_UNIT_SUM - Fatal error!'
    write ( *, * ) '  The vector entries sum to 0.'
    stop
  end if

  a(1:n) = a(1:n) / a_sum

  return
end
subroutine rvec_variance ( n, x, variance )
!
!*******************************************************************************
!
!! RVEC_VARIANCE returns the variance of a real vector.
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
!    Input, real X(N), the vector whose variance is desired.
!
!    Output, real VARIANCE, the variance of the vector entries.
!
  integer n
!
  integer i
  real mean
  real variance
  real x(n)
!
  call rvec_mean ( n, x, mean )

  variance = 0.0
  do i = 1, n
    variance = variance + ( x(i) - mean )**2
  end do

  if ( n > 1 ) then
    variance = variance / real ( n - 1 )
  else
    variance = 0.0
  end if

  return
end
function sech ( X )
!
!*******************************************************************************
!
!! SECH returns the hyperbolic secant.
!
!
!  Definition:
!
!    SECH ( X ) = 1.0 / COSH ( X ) = 2.0 / ( EXP ( X ) + EXP ( - X ) )
!
!  Discussion:
!
!    SECH is not a built-in function in FORTRAN, and occasionally it
!    is handier, or more concise, to be able to refer to it directly
!    rather than through its definition in terms of the sine function.
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument.
!
!    Output, real SECH, the hyperbolic secant of X.
!
  real sech
  real x
!
  sech = 1.0 / cosh ( x )

  return
end
subroutine sech_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! SECH_CDF evaluates the Hyperbolic Secant CDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  real x
  real y
!
  y = ( x - a ) / b

  cdf = 2.0 * atan ( exp ( y ) ) / PI

  return
end
subroutine sech_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! SECH_CDF_INV inverts the Hyperbolic Secant CDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf <= 0.0 ) then
    x = - huge ( x )
  else if ( cdf < 1.0 ) then
    x = a + b * log ( tan ( 0.5 * PI * cdf ) )
  else if ( cdf >= 1.0 ) then
    x = huge ( x )
  end if

  return
end
subroutine sech_check ( a, b )
!
!*******************************************************************************
!
!! SECH_CHECK checks the parameters of the Hyperbolic Secant CDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SECH_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.0'
    stop
  end if

  return
end
subroutine sech_mean ( a, b, mean )
!
!*******************************************************************************
!
!! SECH_MEAN returns the mean of the Hyperbolic Secant PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine sech_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! SECH_PDF evaluates the Hypebolic Secant PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = sech ( ( X - A ) / B ) / ( PI * B )
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real sech
  real x
  real y
!
  y = ( x - a ) / b

  pdf = sech ( y ) / ( PI * b )

  return
end
subroutine sech_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! SECH_SAMPLE samples the Hyperbolic Secant PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  x = a + b * log ( tan ( 0.5 * PI * cdf ) )

  return
end
subroutine sech_variance ( a, b, variance )
!
!*******************************************************************************
!
!! SECH_VARIANCE returns the variance of the Hyperbolic Secant PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real variance
!
  variance = 0.25 * ( PI * b )**2

  return
end
subroutine semicircular_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! SEMICIRCULAR_CDF evaluates the Semicircular CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real cdf
  real x
  real y
!
  if ( x <= a - b ) then

    cdf = 0.0

  else if ( x <= a + b ) then

    y = ( x - a ) / b

    cdf = 0.5 + ( y * sqrt ( 1.0 - y**2 ) + asin ( y ) ) / PI

  else if ( x > a + b ) then

    cdf = 1.0

  end if

  return
end
subroutine semicircular_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! SEMICIRCULAR_CDF_INV inverts the Semicircular CDF.
!
!
!  Discussion:
!
!    A simple bisection method is used on the interval [ A - B, A + B ].
!
!  Modified:
!
!    30 December 1999
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real X, the corresponding argument of the CDF.
!
  integer, parameter :: IT_MAX = 100
  real, parameter :: TOL = 0.0001
!
  real a
  real b
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
!
  if ( cdf <= 0.0 ) then
    x = a - b
    return
  else if ( cdf >= 1.0 ) then
    x = a + b
    return
  end if
!
  x1 = a - b
  cdf1 = 0.0

  x2 = a + b
  cdf2 = 1.0
!
!  Now use bisection.
!
  it = 0

10    continue

  it = it + 1

  x3 = 0.5 * ( x1 + x2 )
  call semicircular_cdf ( x3, a, b, cdf3 )

  if ( abs ( cdf3 - cdf ) < TOL ) then
    x = x3
    return
  end if

  if ( it > IT_MAX ) then
    write ( *, * ) ' '
    write ( *, * ) 'SEMICIRCULAR_CDF_INV - Fatal error!'
    write ( *, * ) '  Iteration limit exceeded.'
    stop
  end if

  if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
    x1 = x3
    cdf1 = cdf3
  else
    x2 = x3
    cdf2 = cdf3
  end if

  go to 10
end
subroutine semicircular_check ( a, b )
!
!*******************************************************************************
!
!! SEMICIRCULAR_CHECK checks the parameters of the Semicircular CDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameter of the PDF.
!    0.0 < B.
!
  real a
  real b
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SEMICIRCULAR_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.0'
    stop
  end if

  return
end
subroutine semicircular_mean ( a, b, mean )
!
!*******************************************************************************
!
!! SEMICIRCULAR_MEAN returns the mean of the Semicircular PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine semicircular_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! SEMICIRCULAR_PDF evaluates the Semicircular PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = ( 2 / ( B * PI ) ) * SQRT ( 1 - ( ( X - A ) / B )**2 )
!    for A - B <= X <= A + B
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real pdf
  real x
  real y
!
  if ( x < a - b ) then

    pdf = 0.0

  else if ( x <= a + b ) then

    y = ( x - a ) / b

    pdf = 2.0 / ( b * PI ) * sqrt ( 1.0 - y**2 )

  else if ( x > a + b ) then

    pdf = 0.0

  end if

  return
end
subroutine semicircular_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! SEMICIRCULAR_SAMPLE samples the Semicircular PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real angle
  real b
  integer iseed
  real radius
  real uniform_01_sample
  real x
!
  radius = b * sqrt ( uniform_01_sample ( iseed ) )
  call r_random ( 0.0, PI, iseed, angle )
  x = a + radius * cos ( angle )

  return
end
subroutine semicircular_variance ( a, b, variance )
!
!*******************************************************************************
!
!! SEMICIRCULAR_VARIANCE returns the variance of the Semicircular PDF.
!
!
!  Modified:
!
!    30 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = b**2 / 4.0

  return
end
subroutine student_central_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! STUDENT_CENTRAL_CDF evaluates the central Student T CDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real a2
  real b
  real b2
  real beta_inc
  real c
  real c2
  real cdf
  real x
  real y
!
  y = ( x - a ) / b

  a2 = c / ( c + y**2 )
  b2 = 0.5 * c
  c2 = 0.5

  if ( y <= 0.0 ) then
    cdf = 0.5 * beta_inc ( a2, b2, c2 )
  else
    cdf = 1.0 - 0.5 * beta_inc ( a2, b2, c2 )
  end if

  return
end
subroutine student_central_check ( a, b, c )
!
!*******************************************************************************
!
!! STUDENT_CENTRAL_CHECK checks the parameter of the central Student T CDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
  real a
  real b
  real c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'STUDENT_CENTRAL_CHECK - Fatal error!'
    write ( *, * ) '  B must be greater than 0.'
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'STUDENT_CENTRAL_CHECK - Fatal error!'
    write ( *, * ) '  C must be greater than 0.'
    stop
  end if
  
  return
end
subroutine student_central_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! STUDENT_CENTRAL_MEAN returns the mean of the central Student T PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, C, the parameter of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real mean
!
  mean = a
  
  return
end
subroutine student_central_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! STUDENT_CENTRAL_PDF evaluates the central Student T PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = Gamma ( (C+1)/2 ) /
!      ( Gamma ( C / 2 ) * Sqrt ( PI * C ) 
!      * ( 1 + ((X-A)/B)**2/C )**(C + 1/2 ) )
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real c
  real gamma
  real pdf
  real x
  real y
!
  y = ( x - a ) / b

  pdf = gamma ( 0.5 * ( c + 1.0 ) ) / ( sqrt ( PI * c ) * gamma ( 0.5 * c ) &
    * sqrt ( ( 1.0 + y**2 / c )**( 2 * c + 1.0 ) ) )
  
  return
end
subroutine student_central_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! STUDENT_CENTRAL_SAMPLE samples the central Student T PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, C, the parameters of the PDF.
!    0 < B,
!    0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real a2
  real a3
  real b
  real b2
  real c
  integer iseed
  real x
  real x2
  real x3
!
  if ( c < 3.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'STUDENT_CENTRAL_SAMPLE - Fatal error!'
    write ( *, * ) '  Sampling fails for C < 3.'
    return
  end if

  a2 = 0.0
  b2 = real ( c ) / real ( c - 2 )

  call normal_sample ( a2, b2, iseed, x2 )

  a3 = c
  call chisquare_central_sample ( a3, iseed, x3 )
  x3 = x3 * c / ( c - 2.0 )

  x = a + b * x2 * sqrt ( c ) / x3

  return
end
subroutine student_central_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! STUDENT_CENTRAL_VARIANCE returns the variance of the central Student T PDF.
!
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, C, the parameter of the PDF.
!    0 < B,
!    0 < C.
!    However, the variance is not defined unless 2 < C
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real c
  real variance
!
  if ( c < 3.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'STUDENT_CENTRAL_VARIANCE - Fatal error!'
    write ( *, * ) '  Variance not defined for C <= 2.'
    stop
  end if

  variance = b**2 * c / ( c - 2.0 )
  
  return
end
subroutine student_noncentral_cdf ( x, idf, d, cdf )
!
!*******************************************************************************
!
!! STUDENT_NONCENTRAL_CDF evaluates the noncentral Student T CDF.
!
!
!  Reference:
!
!    Algorithm AS 5,
!    Applied Statistics,
!    Volume 17, 1968, page 193.
!
!  Modified:
!
!    07 March 1999
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, integer IDF, the number of degrees of freedom.
!
!    Input, real D, the noncentrality parameter.
!
!    Output, real CDF, the value of the CDF.
!
  integer, parameter :: A_MAX = 100
  real, parameter :: EMIN = 12.5
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real ak
  real b
  real cdf
  real cdf2
  real d
  real drb
  real f
  real fk
  real fmkm1
  real fmkm2
  real gamma_log
  integer idf
  integer k
  real sum
  real temp
  real tfn
  real x
!
  f = idf

  if ( idf == 1 ) then

    a = x / sqrt ( f )
    b = f / ( f + x**2 )
    drb = d * sqrt ( b )

    call normal_01_cdf ( drb, cdf2 )
    cdf = 1.0 - cdf2 + 2.0 * tfn ( drb, a )

  else if ( idf <= A_MAX ) then

    a = x / sqrt ( f )
    b = f / ( f + x**2 )
    drb = D * sqrt ( b )
    sum = 0.0

    fmkm2 = 0.0
    if ( abs ( drb ) < EMIN ) then
      call normal_01_cdf ( a * drb, cdf2 )
      fmkm2 = a * sqrt ( b ) * exp ( - 0.5 * drb**2 ) * cdf2 / sqrt ( 2.0 * PI )
    end if

    fmkm1 = b * d * a * fmkm2
    if ( abs ( d ) < EMIN ) then
      fmkm1 = fmkm1 + 0.5 * b * a * exp ( - 0.5 * d**2 ) / PI
    end if

    if ( mod ( idf, 2 ) == 0 ) then
      sum = fmkm2
    else
      sum = fmkm1
    end if

    ak = 1.0

    do k = 2, idf - 2, 2

      fk = real ( k )

      fmkm2 = b * ( d * a * ak * fmkm1 + fmkm2 ) * ( fk - 1.0 ) / fk

      ak = 1.0 / ( ak * ( fk - 1.0 ) )
      fmkm1 = b * ( d * a * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0 )

      if ( mod ( idf, 2 ) == 0 ) then
        sum = sum + fmkm2
      else
        sum = sum + fmkm1
      end if

      ak = 1.0 / ( ak * fk )

    end do

    if ( mod ( idf, 2 ) == 0 ) then
      call normal_01_cdf ( d, cdf2 )
      cdf = 1.0 - cdf2 + sum * sqrt ( 2.0 * PI )
    else
      call normal_01_cdf ( drb, cdf2 )
      cdf = 1.0 - cdf2 + 2.0 * ( sum + tfn ( drb, a ) )
    end if
!
!  Normal approximation.
!
  else

    a = sqrt ( 0.5 * f ) * exp ( gamma_log ( 0.5 * ( f - 1.0 ) ) &
      - gamma_log ( 0.5 * f ) ) * d

    temp = ( x - a ) / sqrt ( f * ( 1.0 + d**2 ) / ( f - 2.0 ) - a**2 )

    call normal_01_cdf ( temp, cdf2 )
    cdf = cdf2

  end if

  return
end
function tfn ( x, fx )
!
!*******************************************************************************
!
!! TFN calculates the T function of Owen.
!
!
!  Reference:
!
!    Algorithm AS 76,
!    Applied Statistics,
!    Volume 23, Number 3, 1974.
!
!  Parameters:
!
!    Input, real X, FX, the arguments of the T function.
!
!    Output, real TFN, the value of the T function.
!
  integer, parameter :: NGAUSS = 5
  real, parameter :: TP = 0.159155
  real, parameter :: TV1 = 1.0E-35
  real, parameter :: TV2 = 15.0
  real, parameter :: TV3 = 15.0
  real, parameter :: TV4 = 1.0E-05
!
  real fx
  real fxs
  integer i
  real, parameter, dimension ( NGAUSS ) :: r = (/ &
    0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357 /)
  real r1
  real r2
  real rt
  real tfn
  real, parameter, dimension ( NGAUSS ) :: u = (/ &
    0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533 /)
  real x
  real x1
  real x2
  real xs
!
!  Test for X near zero.
!
  if ( abs ( x ) < TV1 ) then
    tfn = TP * atan ( fx )
!
!  Test for large values of abs(X).
!
  else if ( abs ( X ) > TV2 ) then
    tfn = 0.0
!
!  Test for FX near zero.
!
  else if ( abs ( fx ) < TV1 ) then
    tfn = 0.0
!
!  Test whether abs(FX) is so large that it must be truncated.
!
  else 

    xs = - 0.5 * x**2
    x2 = fx
    fxs = fx**2
!
!  Computation of truncation point by Newton iteration.
!
    if ( log ( 1.0 + fxs ) - xs * fxs >= TV3 ) then

      x1 = 0.5 * fx
      fxs = 0.25 * fxs

   20     continue

      rt = fxs + 1.0
      x2 = x1 + ( xs * fxs + TV3 - log ( rt ) ) / ( 2.0 * x1 * ( 1.0 / rt - xs ) )
      fxs = x2**2

      if ( abs ( x2 - x1 ) >= TV4 ) then
        x1 = x2
        go to 20
      end if

    end if
!
!  Gaussian quadrature.
!
    rt = 0.0
    do i = 1, NGAUSS
      r1 = 1.0 + fxs * ( 0.5 + u(i) )**2
      r2 = 1.0 + fxs * ( 0.5 - u(i) )**2
      rt = rt + r(i) * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 )
    end do

    tfn = rt * x2 * tp

  end if

  return
end
subroutine triangular_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! TRIANGULAR_CDF evaluates the Triangular CDF.
!
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( x <= a ) then
    cdf = 0.0
  else if ( x <= 0.5 * ( a + b ) ) then
    cdf = 2.0 * ( x**2 - 2.0 * a * x + a**2 ) / ( b - a )**2
  else if ( x <= b ) then
    cdf = 0.5 + ( - 2.0 * x**2 + 4.0 * b * x + 0.5 * a**2 &
      - a * b - 1.5 * b**2 ) / ( b - a )**2
  else
    cdf = 1.0
  end if

  return
end
subroutine triangular_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! TRIANGULAR_CDF_INV inverts the Triangular CDF.
!
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TRIANGULAR_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  if ( cdf <= 0.5 ) then
    x = a + 0.5 * ( b - a ) * sqrt ( 2.0 * cdf )
  else
    x = b - 0.5 * ( b - a ) * sqrt ( 2.0 * ( 1.0 - cdf ) )
  end if

  return
end
subroutine triangular_check ( a, b )
!
!*******************************************************************************
!
!! TRIANGULAR_CHECK checks the parameters of the Triangular CDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
  real a
  real b
!
  if ( a >= b ) then
    write ( *, * ) ' '
    write ( *, * ) 'TRIANGULAR_CHECK - Fatal error!'
    write ( *, * ) '  B <= A.'
    stop
  end if

  return
end
subroutine triangular_mean ( a, b, mean )
!
!*******************************************************************************
!
!! TRIANGULAR_MEAN returns the mean of the Triangular PDF.
!
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real MEAN, the mean of the discrete uniform PDF.
!
  real a
  real b
  real mean
!
  mean = 0.5 * ( a + b )

  return
end
subroutine triangular_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! TRIANGULAR_PDF evaluates the Triangular PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = 4 * ( X - A ) / ( B - A )**2 for A <= X <= (A+B)/2
!               = 4 * ( B - X ) / ( B - A )**2 for (A+B)/2 <= X <= B.
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  if ( x <= a ) then
    pdf = 0.0
  else if ( x <= 0.5 * ( a + b ) ) then
    pdf = 4.0 * ( x - a ) / ( b - a )**2
  else if ( x <= b ) then
    pdf = 4.0 * ( b - x ) / ( b - a )**2
  else
    pdf = 0.0
  end if

  return
end
subroutine triangular_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! TRIANGULAR_SAMPLE samples the Triangular PDF.
!
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call triangular_cdf_inv ( cdf, a, b, x )

  return
end
subroutine triangular_variance ( a, b, variance )
!
!*******************************************************************************
!
!! TRIANGULAR_VARIANCE returns the variance of the Triangular PDF.
!
!
!  Modified:
!
!    21 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = ( b - a )**2 / 24.0

  return
end
function trigamma ( x )
!
!*******************************************************************************
!
!! TRIGAMMA calculates trigamma(x) = d**2 log(Gamma(x)) / dx**2.
!
!
!  Reference:
!
!    B Schneider,
!    Trigamma Function,
!    Algorithm AS 121,
!    Applied Statistics, 
!    Volume 27, Number 1, page 97-99, 1978.
!
!  Modified:
!
!    03 January 2000
!
!  Parameters:
!
!    Input, real X, the argument of the trigamma function.
!    0 < X.
!
!    Output, real TRIGAMMA, the value of the trigamma function at X.
!
  real, parameter :: a = 0.0001
  real, parameter :: b = 5.0
  real, parameter :: b2 =   1.0 / 6.0
  real, parameter :: b4 = - 1.0 / 30.0
  real, parameter :: b6 =   1.0 / 42.0
  real, parameter :: b8 = - 1.0 / 30.0
!
  real trigamma
  real x
  real y
  real z
!
!  1): If X is not positive, fail.
!
  if ( x <= 0.0 ) then

    trigamma = 0.0
    write ( *, * ) ' '
    write ( *, * ) 'TRIGAMNA - Fatal error!'
    write ( *, * ) '  X <= 0.'
    stop
!
!  2): If X is smaller than A, use a small value approximation.
!
  else if ( x <= a ) then

    trigamma = 1.0 / x**2
!
!  3): Otherwise, increase the argument to ( X + I ) >= B...
!
  else

    z = x
    trigamma = 0.0

    do while ( z < b )
      trigamma = trigamma + 1.0 / z**2
      z = z + 1.0
    end do
!
!  ...and then apply an asymptotic formula.
!
    y = 1.0 / z**2

    trigamma = trigamma + 0.5 * &
            y + ( 1.0 &
          + y * ( b2 &
          + y * ( b4 &
          + y * ( b6 &
          + y *   b8 )))) / z

  end if

  return
end
subroutine uniform_01_cdf ( x, cdf )
!
!*******************************************************************************
!
!! UNIFORM_01_CDF evaluates the Uniform 01 CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Output, real CDF, the value of the CDF.
!
  real cdf
  real x
!
  if ( x < 0.0 ) then
    cdf = 0.0
  else if ( x > 1.0 ) then
    cdf = 1.0
  else
    cdf = x
  end if

  return
end
subroutine uniform_01_cdf_inv ( cdf, x )
!
!*******************************************************************************
!
!! UNIFORM_01_CDF_INV inverts the Uniform 01 CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Output, real X, the corresponding argument.
!
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'UNIFORM_01_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = cdf

  return
end
subroutine uniform_01_mean ( mean )
!
!*******************************************************************************
!
!! UNIFORM_01_MEAN returns the mean of the Uniform 01 PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MEAN, the mean of the discrete uniform PDF.
!
  real mean
!
  mean = 0.5

  return
end
subroutine uniform_01_order_sample ( n, iseed, x )
!
!*******************************************************************************
!
!! UNIFORM_01_ORDER_SAMPLE samples the Uniform 01 Order PDF.
!
!
!  Discussion:
!
!    In effect, this routine simply generates N samples of the
!    Uniform 01 PDF; but it generates them in order.  (Actually,
!    it generates them in descending order, but stores them in
!    the array in ascending order).  This saves the work of
!    sorting the results.  Moreover, if the order statistics
!    for another PDF are desired, and the inverse CDF is available,
!    then the desired values may be generated, presorted, by
!    calling this routine and using the results as input to the
!    inverse CDF routine.
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 168.
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
!    Input, integer N, the number of elements in the sample.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X(N), N samples of the Uniform 01 PDF, in
!    ascending order.
!
  integer n
!
  integer i
  integer iseed
  real u
  real uniform_01_sample
  real v
  real x(n)
!
  v = 1.0
  do i = n, 1, -1
    u = uniform_01_sample ( iseed )
    v = v * u**( 1.0 / real ( i ) )
    x(i) = v
  end do

  return
end
subroutine uniform_01_pdf ( x, pdf )
!
!*******************************************************************************
!
!! UNIFORM_01_PDF evaluates the Uniform 01 PDF.
!
!
!  Formula:
!
!    PDF(X) = 1 for 0 <= X <= 1
!           = 0 otherwise
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    0.0 <= X <= 1.0.
!
!    Output, real PDF, the value of the PDF.
!
  real pdf
  real x
!
  if ( x < 0.0 .or. x > 1.0 ) then
    pdf = 0.0
  else
    pdf = 1.0
  end if

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
!    UNIFORM_01_SAMPLE = ISEED * / ( 2**31 - 1 )
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
  uniform_01_sample = real ( iseed ) * 4.656612875e-10

  return
end
subroutine uniform_01_variance ( variance )
!
!*******************************************************************************
!
!! UNIFORM_01_VARIANCE returns the variance of the Uniform 01 PDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real variance
!
  variance = 1.0 / 12.0

  return
end
subroutine uniform_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! UNIFORM_CDF evaluates the Uniform CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  real x
!
  if ( x < a ) then
    cdf = 0.0
  else if ( x > b ) then
    cdf = 1.0
  else
    cdf = ( x - a ) / ( b - a )
  end if

  return
end
subroutine uniform_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! UNIFORM_CDF_INV inverts the Uniform CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real X, the corresponding argument.
!
  real a
  real b
  real cdf
  real x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'UNIFORM_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a + ( b - a ) * cdf

  return
end
subroutine uniform_check ( a, b )
!
!*******************************************************************************
!
!! UNIFORM_CHECK checks the parameters of the Uniform CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
  real a
  real b
!
  if ( a >= b ) then
    write ( *, * ) ' '
    write ( *, * ) 'UNIFORM_CHECK - Fatal error!'
    write ( *, * ) '  B <= A.'
    stop
  end if

  return
end
subroutine uniform_mean ( a, b, mean )
!
!*******************************************************************************
!
!! UNIFORM_MEAN returns the mean of the Uniform PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real MEAN, the mean of the discrete uniform PDF.
!
  real a
  real b
  real mean
!
  mean = 0.5 * ( a + b )

  return
end
subroutine uniform_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! UNIFORM_PDF evaluates the Uniform PDF.
!
!
!  Discussion:
!
!    The Uniform PDF is also known as the "Rectangular" PDF.
!
!  Formula:
!
!    PDF(X)(A,B) = 1 / ( B - A ) for A <= X <= B
!               = 0 otherwise
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  real x
!
  if ( x < a .or. x > b ) then
    pdf = 0.0
  else
    pdf = 1.0 / ( b - a )
  end if

  return
end
subroutine uniform_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! UNIFORM_SAMPLE samples the Uniform PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call uniform_cdf_inv ( cdf, a, b, x )

  return
end
subroutine uniform_variance ( a, b, variance )
!
!*******************************************************************************
!
!! UNIFORM_VARIANCE returns the variance of the Uniform PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    A < B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real variance
!
  variance = ( b - a )**2 / 12.0

  return
end
subroutine uniform_discrete_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! UNIFORM_DISCRETE_CDF evaluates the Uniform Discrete CDF.
!
!
!  Modified:
!
!    29 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the CDF.
!
!    Input, integer A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real CDF, the value of the CDF.
!
  integer a
  integer b
  real cdf
  integer x
!
  if ( x < a ) then
    cdf = 0.0
  else if ( x > b ) then
    cdf = 1.0
  else
    cdf = real ( x + 1 - a ) / real ( b + 1 - a )
  end if

  return
end
subroutine uniform_discrete_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! UNIFORM_DISCRETE_CDF_INV inverts the Uniform Discrete CDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, integer A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, integer X, the smallest argument whose CDF is greater
!    than or equal to CDF.
!
  integer a
  real a2
  integer b
  real b2
  real cdf
  integer x
  real x2
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'UNIFORM_DISCRETE_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  a2 = real ( a ) - 0.5
  b2 = real ( b ) + 0.5
  x2 = a + cdf * ( b2 - a2 )

  x = nint ( x2 )

  x = max ( x, a )
  x = min ( x, b )

  return
end
subroutine uniform_discrete_check ( a, b )
!
!*******************************************************************************
!
!! UNIFORM_DISCRETE_CHECK checks the parameters of the Uniform discrete CDF.
!
!
!  Modified:
!
!    13 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, the parameters of the PDF.
!    A <= B.
!
  integer a
  integer b
!
  if ( a > b ) then
    write ( *, * ) ' '
    write ( *, * ) 'UNIFORM_DISCRETE_CHECK - Fatal error!'
    write ( *, * ) '  A > B.'
    stop
  end if

  return
end
subroutine uniform_discrete_mean ( a, b, mean )
!
!*******************************************************************************
!
!! UNIFORM_DISCRETE_MEAN returns the mean of the Uniform discrete PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real MEAN, the mean of the PDF.
!
  integer a
  integer b
  real mean
!
  mean = 0.5 * real ( a + b )

  return
end
subroutine uniform_discrete_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! UNIFORM_DISCRETE_PDF evaluates the Uniform discrete PDF.
!
!
!  Discussion:
!
!    The Uniform Discrete PDF is also known as the "Rectangular"
!    Discrete PDF.
!
!  Formula:
!
!    PDF(X)(A,B) = 1 / ( B + 1 - A ) for A <= X <= B. 
!
!    The parameters define the interval of integers
!    for which the PDF is nonzero.
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!
!    Input, integer A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real PDF, the value of the PDF.
!
  integer a
  integer b
  real pdf
  integer x
!
  if ( x < a .or. x > b ) then
    pdf = 0.0
  else
    pdf = 1.0 / real ( b + 1 - a )
  end if

  return
end
subroutine uniform_discrete_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! UNIFORM_DISCRETE_SAMPLE samples the Uniform discrete PDF.
!
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, the parameters of the PDF.
!    A <= B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, integer X, a sample of the PDF.
!
  integer a
  integer b
  real cdf
  integer iseed
  real uniform_01_sample
  integer x
!
  cdf = uniform_01_sample ( iseed )

  call uniform_discrete_cdf_inv ( cdf, a, b, x )

  return
end
subroutine uniform_discrete_variance ( a, b, variance )
!
!*******************************************************************************
!
!! UNIFORM_DISCRETE_VARIANCE returns the variance of the Uniform discrete PDF.
!
!
!  Modified:
!
!    28 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, the parameters of the PDF.
!    A <= B.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  integer a
  integer b
  real variance
!
  variance = real ( ( b + 1 - a )**2 - 1 ) / 12.0

  return
end
subroutine uniform_nsphere_sample ( n, iseed, x )
!
!*******************************************************************************
!
!! UNIFORM_NSPHERE_SAMPLE samples the Uniform Unit Sphere PDF.
!
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 168.
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
!    Input, integer N, the dimension of the sphere.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X(N), a point on the unit N sphere, chosen
!    with a uniform probability.
!
  integer n
!
  integer i
  integer iseed
  real norm
  real x(n)
!
  do i = 1, n
    call normal_01_sample ( iseed, x(i) )
  end do

  norm = 0.0
  do i = 1, n
    norm = norm + x(i)**2
  end do
  norm = sqrt ( norm )

  x(1:n) = x(1:n) / norm

  return
end
subroutine von_mises_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! VON_MISES_CDF evaluates the Von Mises CDF.
!
!
!  Reference:
!
!    Geoffrey Hill,
!    ACM TOMS Algorithm 518,
!    Incomplete Bessel Function I0: The Von Mises Distribution,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 3, September 1977, pages 279-284.
!
!  Modified:
!
!    05 April 1999
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!    A - PI <= X <= A + PI.
!
!    Input, real A, B, the parameters of the PDF.
!    -PI <= A <= PI,
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real, parameter :: a1 = 12.0
  real, parameter :: a2 = 0.8
  real, parameter :: a3 = 8.0
  real, parameter :: a4 = 1.0
  real, parameter :: ck = 10.5
  real, parameter :: c1 = 56.0
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real arg
  real b
  real c
  real cdf
  real cn
  real erfx
  real error_function
  integer ip
  integer n
  real p
  real r
  real s
  real sn
  real u
  real v
  real x
  real y
  real z
!
!  We expect -PI <= X - A <= PI.
!
  if ( x - a <= -PI ) then
    cdf = 0.0
    return
  end if

  if ( x - a >= PI ) then
    cdf = 1.0
    return
  end if
!
!  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ). 
!
  z = b

  u = mod ( x - a + PI, 2.0 * PI )

  if ( u < 0.0 ) then
    u = u + 2.0 * PI
  end if

  y = u - PI
!
!  For small B, sum IP terms by backwards recursion.
!
  if ( z <= ck ) then

    v = 0.0

    if ( z > 0.0 ) then

      ip = z * a2 - a3 / ( z + a4 ) + a1
      p = real ( ip )
      s = sin ( y )
      c = cos ( y )
      y = p * y
      sn = sin ( y )
      cn = cos ( y )
      r = 0.0
      z = 2.0 / z

      do n = 2, ip
        p = p - 1.0
        y = sn
        sn = sn * c - cn * s
        cn = cn * c + y * s
        r = 1.0 / ( p * z + r )
        v = ( sn / p + v ) * r
      end do

    end if

    cdf = ( u * 0.5 + v ) / PI
!
!  For large B, compute the normal approximation and left tail.
!
  else

    c = 24.0 * z
    v = c - c1
    r = sqrt ( ( 54.0 / ( 347.0 / v + 26.0 - c ) - 6.0 + c ) / 12.0 )
    z = sin ( 0.5 * y ) * r
    s = 2.0 * z**2
    v = v - s + 3.0
    y = ( c - s - s- 16.0 ) / 3.0
    y = ( ( s + 1.75 ) * s + 83.5 ) / v - y
    arg = z * ( 1.0 - s / y**2 )
    erfx = error_function ( arg )
    cdf = 0.5 * erfx + 0.5

  end if

  cdf = max ( cdf, 0.0 )
  cdf = min ( cdf, 1.0 )

  return
end
subroutine von_mises_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! VON_MISES_CDF_INV inverts the Von Mises CDF.
!
!
!  Discussion:
!
!    A simple bisection method is used on the interval [ A - PI, A + PI ].
!
!  Modified:
!
!    19 December 1999
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!
!    Input, real A, B, the parameters of the PDF.
!    -PI <= A <= PI,
!    0.0 < B.
!
!    Output, real X, the corresponding argument of the CDF.
!    A - PI <= X <= A + PI.
!
  integer, parameter :: IT_MAX = 100
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
  real, parameter :: TOL = 0.0001
!
  real a
  real b
  real cdf
  real cdf1
  real cdf2
  real cdf3
  integer it
  real x
  real x1
  real x2
  real x3
!
  if ( cdf <= 0.0 ) then
    x = a - PI
    return
  else if ( cdf >= 1.0 ) then
    x = a + PI
    return
  end if
!
  x1 = a - PI
  cdf1 = 0.0

  x2 = a + PI
  cdf2 = 1.0
!
!  Now use bisection.
!
  it = 0

10    continue

  it = it + 1

  x3 = 0.5 * ( x1 + x2 )
  call von_mises_cdf ( x3, a, b, cdf3 )

  if ( abs ( cdf3 - cdf ) < TOL ) then
    x = x3
    return
  end if

  if ( it > IT_MAX ) then
    write ( *, * ) ' '
    write ( *, * ) 'VON_MISES_CDF_INV - Fatal error!'
    write ( *, * ) '  Iteration limit exceeded.'
    stop
  end if

  if ( sign ( 1.0, cdf3 - cdf ) == sign ( 1.0, cdf1 - cdf ) ) then
    x1 = x3
    cdf1 = cdf3
  else
    x2 = x3
    cdf2 = cdf3
  end if

  go to 10
end
subroutine von_mises_check ( a, b )
!
!*******************************************************************************
!
!! VON_MISES_CHECK checks the parameters of the Von Mises PDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    -PI <= A <= PI,
!    0.0 < B.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
!
  if ( a < -PI .or. a > PI ) then
    write ( *, * ) ' '
    write ( *, * ) 'VON_MISES_CHECK - Fatal error!'
    write ( *, * ) '  A < -PI or PI < A.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'VON_MISES_MEAN - Fatal error!'
    write ( *, * ) '  B <= 0.0'
    stop
  end if

  return
end
subroutine von_mises_mean ( a, b, mean )
!
!*******************************************************************************
!
!! VON_MISES_MEAN returns the mean of the Von Mises PDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    -PI <= A <= PI,
!    0.0 < B.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real mean
!
  mean = a

  return
end
subroutine von_mises_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! VON_MISES_PDF evaluates the Von Mises PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = EXP ( B * COS ( X - A ) ) / ( 2 * PI * I0(B) )
!
!    where:
! 
!      I0(*) is the modified Bessel function of the first
!      kind of order 0.
!
!    The von Mises distribution for points on the unit circle is
!    analogous to the normal distribution of points on a line.
!    The variable X is interpreted as a deviation from the angle A,
!    with B controlling the amount of dispersion.
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 160.
!
!    D J Best and N I Fisher,
!    Efficient Simulation of the von Mises Distribution,
!    Applied Statistics,
!    Volume 28, Number 2, pages 152-157.
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
!    Input, real X, the argument of the PDF.
!    A - PI <= X <= A + PI.
!
!    Input, real A, B, the parameters of the PDF.
!    -PI <= A <= PI,
!    0.0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real i0b
  real pdf
  real x
!
  if ( x < a - PI ) then
    pdf = 0.0
  else if ( x <= a + PI ) then
    call i0 ( b, i0b )
    pdf = exp ( b * cos ( x - a ) ) / ( 2.0 * PI * i0b )
  else
    pdf = 0.0
  end if

  return
end
subroutine von_mises_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! VON_MISES_SAMPLE samples the Von Mises PDF.
!
!
!  Reference:
!
!    D J Best and N I Fisher,
!    Efficient Simulation of the von Mises Distribution,
!    Applied Statistics,
!    Volume 28, Number 2, pages 152-157.
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
!    Input, real A, B, the parameters of the PDF.
!    -PI <= A <= PI,
!    0.0 < B.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  real a
  real b
  real c
  real f
  integer iseed
  real r
  real rho
  real tau
  real u1
  real u2
  real u3
  real uniform_01_sample
  real x
  real z
!
  tau = 1.0 + sqrt ( 1.0 + 4.0 * b**2 )
  rho = ( tau - sqrt ( 2.0 * tau ) ) / ( 2.0 * b )
  r = ( 1.0 + rho**2 ) / ( 2.0 * rho )

10    continue

  u1 = uniform_01_sample ( iseed )
  z = cos ( PI * u1 )
  f = ( 1.0 + r * z ) / ( r + z )
  c = b * ( r - f )

  u2 = uniform_01_sample ( iseed )
  if ( c * ( 2.0 - c ) <= u2 ) then
    if ( log ( c / u2 ) + 1.0 < c ) then
      go to 10
    end if
  end if

  u3 = uniform_01_sample ( iseed )

  x = a + sign ( 1.0, u3 - 0.5 ) * acos ( f )

  return
end
subroutine weibull_cdf ( x, a, b, c, cdf )
!
!*******************************************************************************
!
!! WEIBULL_CDF evaluates the Weibull CDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the CDF.
!    A <= X.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real c
  real cdf
  real x
  real y
!
  if ( x < a ) then
    cdf = 0.0
  else
    y = ( x - a ) / b
    cdf = 1.0 - 1.0 / exp ( y**c )
  end if

  return
end
subroutine weibull_cdf_inv ( cdf, a, b, c, x )
!
!*******************************************************************************
!
!! WEIBULL_CDF_INV inverts the Weibull CDF.
!
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real CDF, the value of the CDF.
!    0.0 < CDF < 1.0.
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real X, the corresponding argument of the CDF.
!
  real a
  real b
  real c
  real cdf
  real x
!
  x = a + b * ( - log ( 1.0 - cdf ) )**( 1.0 / c )

  return
end
subroutine weibull_check ( a, b, c )
!
!*******************************************************************************
!
!! WEIBULL_CHECK checks the parameters of the Weibull CDF.
!
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
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
  real a
  real b
  real c
!
  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'WEIBULL_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  if ( c <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'WEIBULL_CHECK - Fatal error!'
    write ( *, * ) '  C <= 0.'
    stop
  end if

  return
end
subroutine weibull_mean ( a, b, c, mean )
!
!*******************************************************************************
!
!! WEIBULL_MEAN returns the mean of the Weibull PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real MEAN, the mean of the PDF.
!
  real a
  real b
  real c
  real gamma
  real mean
!
  mean = b * gamma ( ( c + 1.0 ) / c ) + a

  return
end
subroutine weibull_pdf ( x, a, b, c, pdf )
!
!*******************************************************************************
!
!! WEIBULL_PDF evaluates the Weibull PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B,C) = ( C / B ) * ( ( X - A ) / B )**( C - 1 ) 
!     * EXP ( - ( ( X - A ) / B )**C ).
!
!  Discussion:
!
!    The Weibull PDF is also known as the Frechet PDF.
!
!    WEIBULL_PDF(X)(A,B,1) is the Exponential PDF.
!
!    WEIBULL_PDF(X)(0,1,2) is the Rayleigh PDF.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the PDF.
!    A <= X
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real c
  real pdf
  real x
  real y
!
  if ( x < a ) then

    pdf = 0.0

  else

    y = ( x - a ) / b

    pdf = ( c / b ) * y**( c - 1.0 )  / exp ( y**c )

  end if

  return
end
subroutine weibull_sample ( a, b, c, iseed, x )
!
!*******************************************************************************
!
!! WEIBULL_SAMPLE samples the Weibull PDF.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Input/output, integer ISEED, the random number generator seed.
!
!    Output, real X, a sample of the PDF.
!
  real a
  real b
  real c
  real cdf
  integer iseed
  real uniform_01_sample
  real x
!
  cdf = uniform_01_sample ( iseed )

  call weibull_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine weibull_variance ( a, b, c, variance )
!
!*******************************************************************************
!
!! WEIBULL_VARIANCE returns the variance of the Weibull PDF.
!
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the parameters of the PDF.
!    0.0 < B,
!    0.0 < C.
!
!    Output, real VARIANCE, the variance of the PDF.
!
  real a
  real b
  real c
  real gamma
  real g1
  real g2
  real variance
!
  g1 = gamma ( ( c + 2.0 ) / c )
  g2 = gamma ( ( c + 1.0 ) / c )

  variance = b**2 * ( g1 - g2**2 )

  return
end
subroutine weibull_discrete_cdf ( x, a, b, cdf )
!
!*******************************************************************************
!
!! WEIBULL_DISCRETE_CDF evaluates the Discrete Weibull CDF.
!
!
!  Modified:
!
!    17 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the CDF.
!    0 <= X.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A <= 1.0,
!    0.0 < B.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real b
  real cdf
  integer x
!
  if ( x < 0 ) then
    cdf = 0.0
  else
    cdf = 1.0 - ( 1.0 - a )**((x+1)**b)
  end if

  return
end
subroutine weibull_discrete_cdf_inv ( cdf, a, b, x )
!
!*******************************************************************************
!
!! WEIBULL_DISCRETE_CDF_INV inverts the Discrete Weibull CDF.
!
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
!    Input, real CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A <= 1.0,
!    0.0 < B.
!
!    Output, integer X, the corresponding argument.
!
  real a
  real b
  real cdf
  integer i_roundup
  integer x
!
  if ( cdf < 0.0 .or. cdf > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'WEIBULL_DISCRETE_CDF_INV - Fatal error!'
    write ( *, * ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = i_roundup ( ( log ( 1.0 - cdf ) / log ( 1.0 - a ) )**( 1.0 / b ) - 1.0 )

  return
end
subroutine weibull_discrete_check ( a, b )
!
!*******************************************************************************
!
!! WEIBULL_DISCRETE_CHECK checks the parameters of the discrete Weibull CDF.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A <= 1.0,
!    0.0 < B.
!
  real a
  real b
!
  if ( a < 0.0 .or. a > 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'WEIBULL_DISCRETE_CHECK - Fatal error!'
    write ( *, * ) '  A < 0 or 1 < A.'
    stop
  end if

  if ( b <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'WEIBULL_DISCRETE_CHECK - Fatal error!'
    write ( *, * ) '  B <= 0.'
    stop
  end if

  return
end
subroutine weibull_discrete_pdf ( x, a, b, pdf )
!
!*******************************************************************************
!
!! WEIBULL_DISCRETE_PDF evaluates the discrete Weibull PDF.
!
!
!  Formula:
!
!    PDF(X)(A,B) = ( 1 - A )**X**B - ( 1 - A )**(X+1)**B.
!
!  Discussion:
!
!    WEIBULL_DISCRETE_PDF(X)(A,1) = GEOMETRIC_PDF(X)(A)
!
!  Modified:
!
!    15 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!    0 <= X
!
!    Input, real A, B, the parameters that define the PDF.
!    0 <= A <= 1,
!    0 < B.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real b
  real pdf
  integer x
!
  if ( x < 0 ) then
    pdf = 0.0
  else
    pdf = ( 1.0 - a )**(x**b) - ( 1.0 - a )**((x+1)**b)
  end if

  return
end
subroutine weibull_discrete_sample ( a, b, iseed, x )
!
!*******************************************************************************
!
!! WEIBULL_DISCRETE_SAMPLE samples the discrete Weibull PDF.
!
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
!    Input, real A, B, the parameters of the PDF.
!    0.0 <= A <= 1.0,
!    0.0 < B.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer X, a sample of the PDF.
!
  real a
  real b
  real cdf
  integer iseed
  real uniform_01_sample
  integer x
!
  cdf = uniform_01_sample ( iseed )

  call weibull_discrete_cdf_inv ( cdf, a, b, x )

  return
end
function zeta ( p )
!
!*******************************************************************************
!
!! ZETA estimates the Riemann Zeta function.
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
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real P, the power to which the integers are raised.
!    P must be greater than 1.  For integral P up to 20, a
!    precomputed value of ZETA is returned; otherwise the infinite
!    sum is approximated.
!
!    Output, real ZETA, an approximation to the Riemann Zeta function.
!
  real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
!
  integer n
  real p
  double precision sum
  double precision sum_old
  real zeta
!
  if ( p <= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ZETA - Fatal error!'
    write ( *, * ) '  Exponent P <= 1.0.'
    stop
  else if ( p == 2.0 ) then
    zeta = PI**2 / 6.0
  else if ( p == 3.0 ) then
    zeta = 1.2020569032
  else if ( p == 4.0 ) then
    zeta = PI**4 / 90.0
  else if ( p == 5.0 ) then
    zeta = 1.0369277551
  else if ( p == 6.0 ) then
    zeta = PI**6 / 945.0
  else if ( p == 7.0 ) then
    zeta = 1.0083492774
  else if ( p == 8.0 ) then
    zeta = PI**8 / 9450.0
  else if ( p == 9.0 ) then
    zeta = 1.0020083928
  else if ( p == 10.0 ) then
    zeta = PI**10 / 93555.0
  else if ( p == 11.0 ) then
    zeta = 1.0004941886
  else if ( p == 12.0 ) then
    zeta = 1.0002460866
  else if ( p == 13.0 ) then
    zeta = 1.0001227133
  else if ( p == 14.0 ) then
    zeta = 1.0000612482
  else if ( p == 15.0 ) then
    zeta = 1.0000305882
  else if ( p == 16.0 ) then
    zeta = 1.0000152823
  else if ( p == 17.0 ) then
    zeta = 1.0000076372
  else if ( p == 18.0 ) then
    zeta = 1.0000038173
  else if ( p == 19.0 ) then
    zeta = 1.0000019082
  else if ( p == 20.0 ) then
    zeta = 1.0000009540
  else
    sum = 0.0
    n = 0
10      continue
    n = n + 1
    sum_old = sum
    sum = sum + 1.0 / ( dble ( n ) )**p
    if ( real ( sum ) > real ( sum_old ) ) then
      go to 10
    end if

    zeta = real ( sum )

  end if

  return
end
subroutine zipf_cdf ( x, a, cdf )
!
!*******************************************************************************
!
!! ZIPF_CDF evaluates the Zipf CDF.
!
!
!  Discussion:
!
!    Simple summation is used.
!
!  Modified:
!
!    19 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X, the argument of the PDF.
!    1 <= N
!
!    Input, real A, the parameter of the PDF.
!    1.0 < A.
!
!    Output, real CDF, the value of the CDF.
!
  real a
  real, save :: asave = 0.0
  real, save :: c = 0.0
  real cdf
  real pdf
  integer x
  integer y
  real zeta
!
  if ( x < 1 ) then

    cdf = 0.0

  else

    if ( a /= asave ) then

      c = zeta ( a )
      asave = a

    end if

    cdf = 0.0
    do y = 1, x
      pdf = ( 1.0 / real ( y )**a ) / c
      cdf = cdf + pdf
    end do

  end if

  return
end
subroutine zipf_check ( a )
!
!*******************************************************************************
!
!! ZIPF_CHECK checks the parameter of the Zipf PDF.
!
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
!    Input, real A, the parameter of the PDF.
!    1.0 < A.
!
  real a
!
  if ( a <= 1.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ZIPF_CHECK - Fatal error!'
    write ( *, * ) '  A <= 1.'
    stop
  end if

  return
end
subroutine zipf_mean ( a, mean )
!
!*******************************************************************************
!
!! ZIPF_MEAN returns the mean of the Zipf PDF.
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
!    Input, real A, the parameter of the PDF.
!    1.0 < A.
!
!    Output, real MEAN, the mean of the PDF.
!    The mean is only defined for A > 2.
!
  real a
  real mean
  real zeta
!
  if ( a <= 2.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ZIPF_MEAN - Fatal error!'
    write ( *, * ) '  No mean defined for A <= 2.'
    stop
  end if
!
  mean = zeta ( a - 1.0 ) / zeta ( a )

  return
end
subroutine zipf_pdf ( x, a, pdf )
!
!*******************************************************************************
!
!! ZIPF_PDF evaluates the Zipf PDF.
!
!
!  Formula:
!
!    PDF(X)(A) = ( 1 / X**A ) / C
!
!    where the normalizing constant is chosen so that
!
!    C = Sum ( 1 <= I < Infinity ) 1 / I**A.
!
!  Discussion:
!
!    From observation, the frequency of different words in long
!    sequences of text seems to follow the Zipf PDF, with
!    parameter A slightly greater than 1.  The Zipf PDF is sometimes
!    known as the "discrete Pareto" PDF.
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
!    Input, integer X, the argument of the PDF.
!    1 <= N
!
!    Input, real A, the parameter of the PDF.
!    1.0 < A.
!
!    Output, real PDF, the value of the PDF.
!
  real a
  real, save :: asave = 0.0
  real, save :: c = 0.0
  real pdf
  integer x
  real zeta
!
  if ( x < 1 ) then

    pdf = 0.0

  else

    if ( a /= asave ) then

      c = zeta ( a )
      asave = a

    end if

    pdf = ( 1.0 / real ( x )**a ) / c

  end if

  return
end
subroutine zipf_sample ( a, iseed, x )
!
!*******************************************************************************
!
!! ZIPF_SAMPLE samples the Zipf PDF.
!
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer Verlag, 1986, pages 550-551.
!
!  Modified:
!
!    06 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the parameter of the PDF.
!    1.0 < A.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer X, a sample of the PDF.
!
  real a
  real b
  integer iseed
  real t
  real u
  real uniform_01_sample
  real v
  real w
  integer x
!
  b = 2.0**( a - 1.0 )

10    continue

  u = uniform_01_sample ( iseed )
  v = uniform_01_sample ( iseed )
  w = int ( 1.0 / u**( 1.0 / ( a - 1.0 ) ) )
  
  t = ( ( w + 1.0 ) / w )**( a - 1.0 )
  if ( v * w * ( t - 1.0 ) * b > t * ( b - 1.0 ) ) then
    go to 10
  end if

  x = int ( w )

  return
end
subroutine zipf_variance ( a, variance )
!
!*******************************************************************************
!
!! ZIPF_VARIANCE returns the variance of the Zipf PDF.
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
!    Input, real A, the parameter of the PDF.
!    1.0 < A.
!
!    Output, real VARIANCE, the variance of the PDF.
!    The variance is only defined for A > 3.
!
  real a
  real mean
  real variance
  real zeta
!
  if ( a <= 3.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ZIPF_VARIANCE - Fatal error!'
    write ( *, * ) '  No variance defined for A <= 3.0.'
    stop
  end if

  call zipf_mean ( a, mean )

  variance = zeta ( a - 2.0 ) / zeta ( a ) - mean**2

  return
end
